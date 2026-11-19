// ΨLab 15.0: analytic gradients and Hessians of the spin-orbital engine's methods
// (so_derivatives.cpp for closed shells through the RHF framework of correlated_gradient.h and
// correlated_hessian.h; uhf_correlated.cpp for open shells).
#pragma once

#include "macqc/basis.h"
#include "macqc/correlated_gradient.h"
#include "macqc/correlated_hessian.h"
#include "macqc/molecule.h"
#include "macqc/mo_integrals.h"
#include "macqc/rhf.h"
#include "macqc/so_methods.h"

namespace macqc {

// The Hamiltonian of a closed-shell (RHF) reference for the engine.
SpatialHamiltonian closed_shell_hamiltonian(const MOIntegrals& mo, int nocc, double nuclear_repulsion);
// A perturbed closed-shell Hamiltonian (h^(x), g^(x) as MO integrals) for so_second_order; `Hq` must
// outlive the result.
HamiltonianDerivative closed_shell_derivative(const MOIntegrals& Hq);

// The analytic gradient of an engine method on an RHF reference (Lambda equations, relaxed
// densities, Z-vector). `res` from run_so_method on closed_shell_hamiltonian(mo, ...).
CorrelatedGradientResult so_closed_shell_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                                  const MOIntegrals& mo, SOResult& res, const SOOptions& options,
                                                  SODensities* densities = nullptr, int threads = 0,
                                                  double memory_mb = 1024.0,
                                                  const CPHFCallback& cphf = CPHFCallback(),
                                                  const SOProgress& progress = SOProgress());

// The analytic Hessian (and, with `dipoles`, the dipole derivatives) of an engine method on an RHF
// reference.
CorrelatedHessianResult so_closed_shell_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                                const MOIntegrals& mo, SOResult& res, const SOOptions& options,
                                                int threads = 0, double memory_mb = 1024.0, bool dipoles = true,
                                                const SOProgress& progress = SOProgress());

}  // namespace macqc
