// ΨLab 15.0: analytic derivatives of the SCF energy for every reference (RHF, UHF, ROHF).
//
// Gradients (scf_open_derivatives.cpp): the RHF gradient of rhf_derivatives.h for closed shells;
// for UHF and ROHF the Coulomb and exchange two-particle density of the alpha and beta densities
// and the energy-weighted density W = C eps C', eps = sum_s F_s(MO) n_s.
//
// Hessians (scf_open_hessian.cpp): second-derivative integrals contracted with the same densities
// plus the orbital response. UHF: the coupled-perturbed UHF equations of all 3 n_atoms coordinates
// (alpha and beta occupied-virtual blocks, coupled through the Coulomb interaction). ROHF: the
// response of the one orbital set in its non-redundant rotations (closed-open, closed-virtual,
// open-virtual), solved with the exact ROHF orbital Hessian applied by Fock builds.
#pragma once

#include <array>
#include <functional>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_derivatives.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"

namespace macqc {

// The AO two-particle density of the SCF energy (RHF: P P/2 - P P/4; open shells: P P/2 -
// (P_a P_a + P_b P_b)/2 in exchange form) and its energy-weighted density.
TwoParticleDensity scf_two_particle_density(const RHFResult& scf);
Matrix scf_energy_weighted_density(const RHFResult& scf);

// Analytic gradient of the SCF energy of any reference (`lindep_tol`: the SCF's threshold).
RHFGradientResult scf_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf, int threads = 0,
                               double lindep_tol = 1.0e-7);

// Analytic Hessian and dipole derivatives of the SCF energy of any reference (RHF: rhf_hessian).
RHFHessianResult scf_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf, int threads = 0,
                             double memory_mb = 1024.0, double cphf_tol = 1.0e-9,
                             const CPHFCallback& progress = CPHFCallback());

// Krylov (Galerkin) solver for several symmetric linear systems A x = b at once (the open-shell
// response and Z-vector equations): `apply` maps a batch of vectors to A times them (one call per
// iteration, so the Fock builds of all right-hand sides share one pass over the integrals);
// `diagonal` preconditions.
struct LinearBatchSolution {
    std::vector<std::vector<double>> x;
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
};
LinearBatchSolution solve_symmetric_batch(
    const std::function<std::vector<std::vector<double>>(const std::vector<std::vector<double>>&)>& apply,
    const std::vector<double>& diagonal, const std::vector<std::vector<double>>& rhs, double tol, int max_iter,
    const CPHFCallback& progress = CPHFCallback());

}  // namespace macqc
