// Analytic second derivatives (Hessians) of correlated energies (ΨLab 10.0).
//
// A method supplies, as for its analytic gradient (correlated_gradient.h), the one- and
// two-particle densities D, G of its energy Lagrangian in the canonical RHF orbitals,
//   L = sum D h + sum G (pq|rs) + V_nn,
// and in addition their first-order responses D^(y), G^(y) to every nuclear coordinate y
// (from the response of its own parameters, e.g. the Full CI vector or the coupled-cluster
// amplitudes and multipliers, to the perturbed Hamiltonian H^(y) built here). The Hessian is
// the derivative of the analytic gradient,
//   E^xy = sum D h^(xy) + sum G g^(xy) + sum D^(y) h^(x) + sum G^(y) g^(x) + V_nn^xy,
// with h^(x), g^(x) the MO integrals differentiated with the orbitals, C(x) = C U(x):
//   h^(x) = h~^x + U^x' h + h U^x,  g^(x) = g~^x + one-index transformations of g by U^x,
// h~^x, g~^x the derivative AO integrals in the reference MOs ("skeleton"). U^x is the RHF
// orbital response: U_ai from the coupled-perturbed HF equations of all 3 n_atoms coordinates,
// U_ia = -S~_ia - U_ai, and U_ij = -S~_ij/2, U_ab = -S~_ab/2 (the energy is invariant to
// occupied-occupied and virtual-virtual rotations). The second-order response U^xy appears only
// through the orbital Lagrangian X of the gradient: its symmetric part follows from the
// orthonormality of the orbitals, its occupied-virtual part from the second-order Brillouin
// condition, which the Z-vector of the gradient contracts without solving second-order CPHF
// equations (Handy & Schaefer, J. Chem. Phys. 81, 5031 (1984)). The second-derivative
// integrals are contracted in the AO basis with the relaxed densities of the gradient (the
// relaxed one-particle density, the two-particle density with its Z-vector term, and the
// energy-weighted density).
//
// With the RHF densities and D^(y) = G^(y) = 0 the result is the analytic RHF Hessian, which
// tests the orbital part independently of any correlated method.
#pragma once

#include <cstddef>
#include <vector>

#include "macqc/basis.h"
#include "macqc/correlated_gradient.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

// The first-order orbital response of the RHF reference to every nuclear coordinate and the
// skeleton derivative integrals in the reference MOs.
struct OrbitalResponse {
    std::size_t n3 = 0;                 // 3 n_atoms
    std::vector<Matrix> S;              // S~^x = C' S^x C (n x n)
    std::vector<Matrix> U;              // U^x (n x n), columns = perturbed orbitals
    std::vector<Matrix> h;              // h~^x = C' h^x C
    std::vector<Packed4> g;             // g~^x = (pq|rs)^x in the reference MOs
    std::vector<Matrix> F;              // F~^x: total derivative of the AO Fock matrix in the
                                        // reference MOs, explicit + density response
    std::vector<Matrix> P;              // P^x (AO): derivative of the RHF density with fixed AO
                                        // integrals, 2 (C U^x_o Co' + Co U^x_o' C')
    bool cphf_converged = false;
    int cphf_iterations = 0;
    double cphf_residual = 0.0;
    double integral_seconds = 0.0, cphf_seconds = 0.0;
};

// Memory (bytes) of rhf_orbital_response for n orbitals, N basis functions, n_atoms atoms.
double orbital_response_bytes(std::size_t n_mo, std::size_t n_basis, int n_atoms);

OrbitalResponse rhf_orbital_response(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads = 0,
                                     double memory_mb = 1024.0, double cphf_tol = 1.0e-9,
                                     const CPHFCallback& progress = CPHFCallback());

// The perturbed MO Hamiltonian H^(x) (h^(x) and g^(x), with the orbital response) of one
// coordinate, as MO integrals (the constant V_nn^x is not included).
MOIntegrals perturbed_hamiltonian(const MOIntegrals& mo, const OrbitalResponse& r, std::size_t x, int threads = 0);

// The same for a perturbation along a direction of the nuclear coordinates (3 n_atoms
// components): H^(d) = sum_x d_x H^(x).
MOIntegrals perturbed_hamiltonian(const MOIntegrals& mo, const OrbitalResponse& r, const std::vector<double>& direction,
                                  int threads = 0);

// First-order responses of the Lagrangian densities.
//  - `internal` false: one per Cartesian coordinate y, D^(y), G^(y).
//  - `internal` true: one per direction q of `directions`, the orthonormal Cartesian directions
//    orthogonal to the translations and rotations (vibrational_directions). correlated_hessian
//    then computes the internal block q' H q and completes the Hessian exactly from the
//    translational and rotational invariance (H u = 0 for a translation, the rotated gradient
//    for a rotation; the gradient is <D, H^(x)> + V_nn^x, from the same densities).
//  - `symmetric` (with `internal`; coupled cluster, 2n+2 rule): D, G are only the amplitude
//    parts dD[t^q] = dD(t + s t^q, l)/ds of the density responses and `pair` holds
//    t^q' L_tt t^q (second derivatives of the Lagrangian with the fixed multipliers); the
//    response term of the internal block is <dD[q], H^(q')> + <dD[q'], H^(q)> + pair(q, q'),
//    without the multiplier responses (by the derivative of the Lambda equations).
struct DensityResponse {
    bool internal = false;
    std::vector<std::vector<double>> directions;
    std::vector<Matrix> D;   // D^(y) (n x n)
    std::vector<Packed4> G;  // G^(y)
    bool symmetric = false;
    Matrix pair;             // m x m
};

struct CorrelatedHessianResult {
    Matrix hessian;                   // 3 n_atoms x 3 n_atoms (E_h/bohr^2)
    bool zvector_converged = false;
    int zvector_iterations = 0;
    double zvector_residual = 0.0;
    double asymmetry = 0.0;           // largest |H_xy - H_yx| before symmetrisation
    std::vector<double> gradient;     // <D, H^(x)> + V_nn^x (E_h/bohr)
    std::size_t perturbations = 0;    // directions (or coordinates) with a method response
    double seconds = 0.0, explicit_seconds = 0.0, response_seconds = 0.0;
    // the drivers below: the orbital response and the method's own response equations
    bool cphf_converged = true;
    int cphf_iterations = 0;
    double cphf_residual = 0.0;
    bool method_converged = true;     // e.g. the Full CI vector responses
    int method_iterations = 0;
    double method_residual = 0.0;
    double orbital_seconds = 0.0, method_seconds = 0.0;
};

// `mo`: MO integrals of all orbitals at the reference; `dens`: the Lagrangian densities;
// `dresp`: their responses (empty vectors: zero, e.g. RHF).
CorrelatedHessianResult correlated_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                           const MOIntegrals& mo, const CorrelatedDensity& dens,
                                           const OrbitalResponse& resp, const DensityResponse& dresp,
                                           int threads = 0, double memory_mb = 1024.0);

// Full CI (stage 10.2): the Hessian of the Full CI energy with the normalised CI vector `c` of
// the RHF orbitals (`fci_electronic_energy` without V_nn). The CI vector responses solve
// (H - E) c^y = -(H^(y) - E^(y)) c; the density responses are the bilinear densities
// D(c, c^y) + D(c^y, c) = [D(c + c^y) - D(c - c^y)] / 2 (exact: the densities are quadratic).
// Responses along the 3 n_atoms - 6 (- 5) internal directions only.
CorrelatedHessianResult fci_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                    const MOIntegrals& mo, const std::vector<double>& c, double fci_electronic_energy,
                                    int threads = 0, double memory_mb = 1024.0, double fci_memory_mb = 4096.0);

// CCSD (stage 10.3): the Hessian of the closed-shell CCSD energy with the converged amplitudes
// t1 [i][a], t2 [i][j][a][b] and Lambda multipliers l1, l2 (cc_gradient.h). Only the amplitude
// responses J t^q = -Omega(t; H^(q)) along the internal directions are solved (2n+2 rule, see
// DensityResponse and cc_hessian.cpp); `tol` is the largest residual element of these equations.
CorrelatedHessianResult ccsd_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                     const MOIntegrals& mo, const std::vector<double>& t1, const std::vector<double>& t2,
                                     const std::vector<double>& l1, const std::vector<double>& l2, int threads = 0,
                                     double memory_mb = 1024.0, double tol = 1.0e-9);

// Packed 4-index transformation out_pqrs = sum M_ap M_bq M_cr M_ds in_abcd of an 8-fold
// symmetric array in the layout pair(pair(a,b), pair(c,d)) over ns indices; M is ns x nd.
std::vector<double> transform_packed(const std::vector<double>& in, std::size_t ns, const Matrix& M, int threads = 0);

}  // namespace macqc
