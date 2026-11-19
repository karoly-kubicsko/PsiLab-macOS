// ΨLab 15.0: analytic gradients of correlated methods on a UHF reference (uhf_correlated.cpp).
#pragma once

#include <array>
#include <vector>

#include "macqc/basis.h"
#include "macqc/correlated_gradient.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"
#include "macqc/so_methods.h"

namespace macqc {

struct UHFCorrelatedGradient {
    std::vector<double> gradient;      // 3 n_atoms (E_h/bohr)
    std::vector<double> nuclear, one_electron, two_electron, overlap;
    Matrix relaxed_density_ao;         // total (alpha + beta) relaxed one-particle density
    std::array<double, 3> relaxed_dipole{};
    double orbital_lagrangian_oo_vv = 0.0;  // largest |X_pq - X_qp| in the oo and vv blocks of either spin
    bool zvector_converged = false;
    int zvector_iterations = 0;
    double zvector_residual = 0.0;
    double density_energy = 0.0;
    double seconds = 0.0, zvector_seconds = 0.0, backtransform_seconds = 0.0, integral_seconds = 0.0;
};

// `aa`, `bb`: the MO integrals of the alpha and the beta orbitals (with h), `ab` the mixed ones
// (transform_to_mo_mixed(C_a, C_b)); `d` the Lagrangian densities of the method (so_densities).
UHFCorrelatedGradient uhf_correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                             const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                             const SODensities& d, int threads = 0, double memory_mb = 1024.0,
                                             const CPHFCallback& progress = CPHFCallback());

}  // namespace macqc
