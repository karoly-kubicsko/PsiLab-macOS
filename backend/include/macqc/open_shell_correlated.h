// ΨLab 15.0: analytic gradients and Hessians of the spin-orbital engine's methods (MP2-MP5, CISD ...
// CISDTQ, QCISD ... QCISDTQ, CCSD ... CCSDTQ) and of Full-CI on open-shell references, UHF and ROHF
// (open_shell_correlated.cpp).
//
// The energy is the Lagrangian
//   L = <Gamma(t, lambda), H(C)> + z . g(C),
// <Gamma, H> = sum_s D_s h_s + sum G_aa (pq|rs)_aa + G_bb (pq|rs)_bb + G_ab (pq|rs)_ab (the method's
// densities from so_densities, the reference included), g the stationarity conditions of the SCF
// (antisymmetric part of eps_s = sum_s' F^s'_s N^s'_s on the non-redundant rotations: UHF the
// occupied-virtual pairs of each spin, ROHF the closed-open, closed-virtual and open-virtual pairs
// of the common orbitals, in non-canonical form) and z their multipliers (the Z-vector). g is
// linear in H, z . g = sum_s' <R_s', F_s'>, so its density Gamma_z (one-particle R, two-particle
// R x P with Coulomb and exchange parts) joins Gamma: Gamma_tot = Gamma + Gamma_z.
//
// Gradient (orbitals C(x) = C T(x), T = (C' S(x) C)^(-1/2)):
//   dE/dx = <Gamma_tot, H~^x> - 1/2 sum_s <Y_s, S~^x_s> + V_nn^x
// with H~^x the derivative integrals in the reference MOs and Y the generalized Fock matrix
// (orbital Lagrangian) of Gamma_tot: Y^s_tp = 2 sum h_tq D_qp + 4 sum (tq|rs) G_pqrs + 2 sum (tq|rs)_ab
// G_ab,pqrs (summed over the spins of an orbital set). z solves A z = -(Y - Y')_pairs with A the
// orbital Hessian of the reference (UHF or ROHF), so Y_tot is symmetric (and the correlated energy
// is invariant to the redundant rotations).
//
// Hessian (2n+2 rule: first-order responses of the parameters t (amplitudes) and kappa (orbital
// rotations), zeroth-order multipliers lambda and z):
//   E^xy = <Gamma_tot, H^[xy]> + <dD[t^x], H^[y]> + <dD[t^y], H^[x]> + t^y L_tt t^x + V_nn^xy,
// H^[x] = H~^x + H(U^x) the first-order MO Hamiltonian with the orbital response U^x = -S~^x/2 +
// kappa^x (coupled-perturbed UHF or ROHF), and
//   <Gamma_tot, H^[xy]> = <Gamma_tot, H^xy (AO second derivatives)> + <Y[H^[x]], U^y> + <Y[H~^y], U^x>
//                         - <Y, U^x' U^y + S~^x U^y + S~^y U^x + U^x U^y> - <W, S^xy>,
// W = 1/2 sum C Y C' (the second-order orbital terms follow from the orthonormality of the orbitals
// because Y_tot is symmetric). Translations and rotations complete the Hessian from the internal
// directions (vibrational_directions) by invariance; electric-field components give the dipole
// derivatives in the same way.
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/correlated_hessian.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"
#include "macqc/so_methods.h"

namespace macqc {

struct OpenShellGradient {
    std::vector<double> gradient;      // 3 n_atoms (E_h/bohr)
    std::vector<double> nuclear, one_electron, two_electron, overlap;
    Matrix relaxed_density_ao;         // alpha + beta, with the Z-vector term
    std::array<double, 3> relaxed_dipole{};
    double redundant_asymmetry = 0.0;  // largest |Y_pq - Y_qp| on the redundant rotations (invariance)
    double stationarity = 0.0;         // largest |Y_pq - Y_qp| on the non-redundant ones after the Z-vector
    bool zvector_converged = false;
    int zvector_iterations = 0;
    double zvector_residual = 0.0;
    double density_energy = 0.0;       // <Gamma_tot, H> + V_nn (= the energy)
    double seconds = 0.0;
};

// `aa`, `bb`, `ab`: the MO integrals of the reference (UHF: alpha, beta and mixed; ROHF: the same
// integrals three times); `d`: the method's densities (so_densities).
OpenShellGradient open_shell_correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                 const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                 const SODensities& d, int threads = 0, double memory_mb = 1024.0,
                                                 const CPHFCallback& progress = CPHFCallback());

// The analytic Hessian (and with `dipoles` the dipole derivatives); `res` from run_so_method on
// the same reference (its state is used for the Lambda equations and the amplitude responses).
CorrelatedHessianResult open_shell_correlated_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                      const MOIntegrals& aa, const MOIntegrals& bb,
                                                      const MOIntegrals& ab, SOResult& res, const SOOptions& options,
                                                      int threads = 0, double memory_mb = 1024.0, bool dipoles = true,
                                                      const SOProgress& progress = SOProgress());

// The same for the SCF energy alone (no correlated method: the densities of the reference); a
// test of the orbital part against scf_gradient and scf_hessian.
OpenShellGradient open_shell_scf_gradient_check(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                int threads = 0);
CorrelatedHessianResult open_shell_scf_hessian_check(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                     const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                     int threads = 0, bool dipoles = true);

// Full-CI on an open-shell reference (n_alpha != n_beta) in the orbitals of one set, scf.C (ROHF; UHF:
// the alpha orbitals), with the integrals `mo` in these orbitals and the (spin-adapted) CI vector c.
// The energy is invariant to every rotation of the orbitals (all of them are correlated), so the
// generalized Fock matrix of the Full-CI densities is symmetric by itself: no Z-vector, and the
// orbitals follow the geometry by symmetric orthonormalisation alone (U^x = -S~^x/2, no CPHF).
// The Hessian adds the CI vector responses c^x (fci_vector_responses, in the space of the
// requested spin): E^xy = <Gamma, H^[xy]> + 2 <c^x|H^[y]|c> (variational).
OpenShellGradient open_shell_fci_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                          const MOIntegrals& mo, const std::vector<double>& c, int threads = 0,
                                          double memory_mb = 1024.0);
CorrelatedHessianResult open_shell_fci_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                               const MOIntegrals& mo, const std::vector<double>& c,
                                               double electronic_energy, int threads = 0, double memory_mb = 1024.0,
                                               double fci_memory_mb = 4096.0, bool dipoles = true);
double open_shell_fci_gradient_bytes(int n_alpha, int n_beta, int n_mo, int n_basis);
double open_shell_fci_hessian_bytes(int n_alpha, int n_beta, int n_mo, int n_basis, int n_atoms);

// Memory (bytes) of open_shell_correlated_hessian: n orbitals, N basis functions.
double open_shell_hessian_bytes(const std::string& method, const int nocc[2], int n_mo, int n_basis, int n_atoms);

}  // namespace macqc
