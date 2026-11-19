// ΨLab 15.0: correlated methods of the spin-orbital engine (so_engine.h) for RHF, UHF and ROHF
// references: Moller-Plesset perturbation theory to fifth order, truncated CI (CISD, CISDT,
// CISDTQ), quadratic CI (QCISD, QCISDT, QCISDTQ), coupled cluster (CCSD, CCSDT, CCSDTQ), and the
// perturbative triples of CISD(T), QCISD(T) and CCSD(T).
//
// Every method is a set of amplitude equations R(t; H) = 0 and an energy E(t; H), polynomials in the
// amplitudes and linear in the Hamiltonian (h, <pq||rs>), assembled from the verified terms of
// backend/tools/so_equations.py (connected coupled-cluster terms up to quadruples with the full Fock
// matrix, and the disconnected terms of linear CI):
//   CC-n     all connected terms with clusters up to rank n;
//   QCI-n    (Pople, Head-Gordon, Raghavachari 1987; He and Cremer 1994 for QCISDT, QCISDTQ) the
//            connected terms linear in the clusters (CI-n made size-extensive), plus in the rank-m
//            equation the quadratic coupling (T2 T_m)_C, i.e. T1 T2 (singles), T2^2/2 (doubles),
//            T2 T3 (triples), T2 T4 and T2^2/2 (quadruples); E = <0|F T1 + V T2|0>;
//   CI-n     <mu|H_N (1 + C)|0> = E_c c_mu with the connected terms linear in C and the
//            disconnected ones, E_c = <0|H_N C|0> (intermediate normalisation);
//   MPn      the order-by-order expansion of the coupled-cluster equations with H0 = f_oo + f_vv
//            (Moller-Plesset partitioning; RHF/UHF: MPn, ROHF: the ROMPn/RMPn of Knowles et al. in
//            invariant form): R^(k)(T^(1..k)) = 0, E^(k) = <0|(V e^T)_C|0>^(k); E(MPn) = E_ref +
//            E^(2) + ... + E^(n). MP4 needs T3^(2), MP5 T4^(3) (connected quadruples);
//   (T)      t3 from <T|[F_oo + F_vv, T3]|0> + <T|(V T2)_C|0> = 0 and
//            E_(T) = <t3|(V T2)_C> + f <t3|(V T1 + F_ov T2)_DC> with f = 1 (CCSD(T)) or 2 (QCISD(T),
//            CISD(T): Psi4's fnocc and Gaussian's QCISD(T) use 2 for the singles-triples term).
// Because the Fock matrix enters in full, the equations hold in any orbitals (invariance to
// occupied-occupied and virtual-virtual rotations) and their solutions are differentiable: the
// Lagrangian L = E + sum lambda R gives the relaxed one- and two-particle densities (derivatives
// with respect to h and (pq|rs), so_derivatives.cpp) and, with the amplitude and multiplier
// responses, the analytic Hessians.
#pragma once

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "macqc/correlated_gradient.h"
#include "macqc/linalg.h"
#include "macqc/mo_integrals.h"
#include "macqc/so_engine.h"

namespace macqc {

// The canonical names of the methods of the engine ("MP2" ... "QCISDTQ", "CCSD", "CCSD(T)", ...).
bool so_engine_method(const std::string& canonical);
// Highest excitation rank stored by a method (MP4 3, MP5 4, CISDT 3, ...).
int so_method_max_rank(const std::string& canonical);

// The MO Hamiltonian of a reference in spatial orbitals of each spin.
struct SpatialHamiltonian {
    int n[2] = {0, 0};        // orbitals per spin (alpha, beta)
    int nocc[2] = {0, 0};     // occupied orbitals per spin
    Matrix h[2];              // MO core Hamiltonian per spin
    // (pq|rs): aa (alpha orbitals), bb (beta), ab (p, q alpha, r, s beta; a mixed transformation, or
    // for a closed-shell reference the alpha integrals themselves)
    const MOIntegrals* aa = nullptr;
    const MOIntegrals* bb = nullptr;
    const MOIntegrals* ab = nullptr;
    double nuclear_repulsion = 0.0;
    bool restricted = true;   // alpha and beta orbitals identical (RHF, ROHF)
    double g(int sp, int sq, std::size_t p, std::size_t q, std::size_t r, std::size_t s) const;
};

// A perturbed Hamiltonian H^(x) (derivative MO integrals) in the same layout: h[2] and the
// two-electron derivative integrals through a callback g(sp, sq, p, q, r, s).
struct HamiltonianDerivative {
    Matrix h[2];
    std::function<double(int, int, std::size_t, std::size_t, std::size_t, std::size_t)> g;  // empty: no 2e part
};

struct SOOptions {
    int max_iter = 150;
    double amplitude_tol = 1.0e-9;   // largest |delta t| (also Lambda, responses)
    double energy_tol = 1.0e-10;
    int diis_size = 8;
    int threads = 0;
    double memory_mb = 4096.0;       // DIIS vectors are limited to a quarter of this
};

struct SOIteration {
    int iteration = 0;
    double energy = 0.0;      // total energy
    double delta_e = 0.0;
    double residual = 0.0;    // largest |delta t|
    std::string stage;        // e.g. "amplitudes", "(T) amplitudes", "MP order 3"
};

using SOProgress = std::function<void(const SOIteration&)>;

struct SOEnergyTerm {
    std::string label;        // "E(2)", "E(3)", "(T)", "E_c", ...
    double value = 0.0;
};

struct SOAmplitudeInfo {
    std::string name;         // "T2", "T2(1)", "C3", "t3 (T)"
    int rank = 0;
    double norm = 0.0;        // over unique amplitudes
    double count = 0.0;       // unique amplitudes
};

class SOMethodState;  // graph, engine and amplitudes (kept for the derivatives)

struct SOResult {
    std::string method;
    std::string reference;    // "RHF", "UHF", "ROHF"
    bool converged = false;
    int iterations = 0;
    double reference_energy = 0.0;
    double correlation_energy = 0.0;   // total - reference
    double total_energy = 0.0;
    std::vector<SOEnergyTerm> terms;    // components (MP orders, CI E_c, (T))
    std::vector<SOAmplitudeInfo> amplitudes;
    double t1_diagnostic = 0.0;
    std::vector<SOIteration> history;
    double seconds = 0.0;
    std::shared_ptr<SOMethodState> state;  // for densities and responses
};

// Size of a method: stored amplitudes (dense spin blocks) and memory estimate (bytes).
struct SOSize {
    double amplitudes = 0.0;
    double memory_bytes = 0.0;
    double flops_per_iteration = 0.0;
};
SOSize estimate_so_method(const std::string& method, const int nocc[2], const int nvir[2]);

SOResult run_so_method(const std::string& method, const SpatialHamiltonian& H, const std::string& reference,
                       const SOOptions& options, const SOProgress& progress = SOProgress());

// Lambda equations of the method (stage by stage) and the relaxed densities of its Lagrangian in the
// spatial orbitals of each spin: E = sum_s D_s h_s + sum G_ss' (pq|rs)_ss' + V_nn.
struct SODensities {
    Matrix D[2];                          // alpha, beta (n x n)
    // two-particle densities in chemists' order over all indices: G_aa[p][q][r][s] etc. (full arrays)
    std::vector<double> Gaa, Gbb, Gab;    // Gab: p, q alpha; r, s beta
    int n[2] = {0, 0};
    bool lambda_converged = false;
    int lambda_iterations = 0;
    double lambda_residual = 0.0;
    double energy = 0.0;                  // the Lagrangian (= the energy)
};
SODensities so_densities(SOResult& result, const SOOptions& options, const SOProgress& progress = SOProgress());

// Closed shell (RHF reference): the spin-summed densities of the RHF framework (correlated_gradient.h).
CorrelatedDensity so_closed_shell_density(const SODensities& d);

// Second-order quantities of the analytic Hessians (2n+2 rule; the Lagrangian is linear in the
// Hamiltonian): for every perturbation x the amplitude response t^x from J t^x = -R_H H^x (stage by
// stage), the amplitude part of the density response dD[t^x] = d/ds dL/dH (t + s t^x, lambda) and the
// pair terms pair(x, y) = t^y L_tt t^x. The response part of the Hessian is then
//   <dD[t^x], H^y> + <dD[t^y], H^x> + pair(x, y)   (no multiplier responses).
// so_densities() must have been called (multipliers).
struct SOSecondOrder {
    std::vector<SODensities> dD;
    Matrix pair;
    bool converged = true;
    int iterations = 0;
    double residual = 0.0;
};
SOSecondOrder so_second_order(SOResult& result, const std::vector<HamiltonianDerivative>& perturbations,
                              const SOOptions& options, const SOProgress& progress = SOProgress());

}  // namespace macqc
