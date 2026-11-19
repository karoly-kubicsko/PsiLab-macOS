// ΨLab 15.0: analytic gradients and Hessians of the methods of the spin-orbital engine.
//
// Closed shells (RHF reference): the spin-summed Lagrangian densities of the engine (so_densities)
// enter the relaxed-density gradient of correlated_gradient.h (Z-vector, derivative integrals); the
// Hessian uses the general framework of correlated_hessian.h in its symmetric (2n+2) form: the
// amplitude responses to the perturbed Hamiltonians H^(q) of the internal directions (and of the
// field components, for the dipole derivatives), the amplitude parts of the density responses and
// the pair terms t^q' L_tt t^q from one forward-over-reverse pass per perturbation (so_second_order).
// Every method of the engine has a Lagrangian linear in the Hamiltonian, so no second derivative
// with respect to H appears.
#include "macqc/so_derivatives.h"

#include <chrono>
#include <cmath>
#include <stdexcept>

#include "macqc/gradient.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"

namespace macqc {

namespace {
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }
}  // namespace

SpatialHamiltonian closed_shell_hamiltonian(const MOIntegrals& mo, int nocc, double nuclear_repulsion) {
    SpatialHamiltonian H;
    H.n[0] = H.n[1] = mo.n;
    H.nocc[0] = H.nocc[1] = nocc;
    H.h[0] = H.h[1] = mo.h;
    H.aa = H.bb = H.ab = &mo;
    H.nuclear_repulsion = nuclear_repulsion;
    H.restricted = true;
    return H;
}

HamiltonianDerivative closed_shell_derivative(const MOIntegrals& Hq) {
    HamiltonianDerivative d;
    d.h[0] = d.h[1] = Hq.h;
    const MOIntegrals* p = &Hq;
    d.g = [p](int, int, std::size_t a, std::size_t b, std::size_t c, std::size_t e) { return (*p)(a, b, c, e); };
    return d;
}

CorrelatedGradientResult so_closed_shell_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                                  const MOIntegrals& mo, SOResult& res, const SOOptions& opt,
                                                  SODensities* densities, int threads, double memory_mb,
                                                  const CPHFCallback& cphf, const SOProgress& progress) {
    SODensities d = so_densities(res, opt, progress);
    const CorrelatedDensity cd = so_closed_shell_density(d);
    CorrelatedGradientResult g = correlated_gradient(mol, basis, rhf, mo, cd, threads, memory_mb, cphf);
    if (densities) *densities = std::move(d);
    return g;
}

CorrelatedHessianResult so_closed_shell_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                                const MOIntegrals& mo, SOResult& res, const SOOptions& opt,
                                                int threads, double memory_mb, bool dipoles,
                                                const SOProgress& progress) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const OrbitalResponse resp = rhf_orbital_response(mol, basis, rhf, nt, memory_mb, 1.0e-10, CPHFCallback(), dipoles);
    const double orbital_seconds = seconds_since(t0);
    const auto t1 = Clock::now();
    SODensities d0 = so_densities(res, opt, progress);
    const CorrelatedDensity D0 = so_closed_shell_density(d0);
    DensityResponse dr;
    dr.internal = true;
    dr.symmetric = true;
    dr.directions = vibrational_directions(mol);
    const std::size_t m = dr.directions.size();
    const std::size_t mr = m + resp.n_fields;
    std::vector<MOIntegrals> Hq(mr);
    std::vector<HamiltonianDerivative> pert(mr);
    for (std::size_t k = 0; k < mr; ++k) {
        Hq[k] = k < m ? perturbed_hamiltonian(mo, resp, dr.directions[k], nt)
                      : perturbed_hamiltonian(mo, resp, resp.n3 + (k - m), nt);
        pert[k] = closed_shell_derivative(Hq[k]);
    }
    const SOSecondOrder so2 = so_second_order(res, pert, opt, progress);
    dr.D.resize(mr);
    dr.G.resize(mr);
    for (std::size_t k = 0; k < mr; ++k) {
        const CorrelatedDensity c = so_closed_shell_density(so2.dD[k]);
        dr.D[k] = c.D;
        dr.G[k] = c.G;
    }
    dr.pair = so2.pair;
    const double method_seconds = seconds_since(t1);
    CorrelatedHessianResult r = correlated_hessian(mol, basis, rhf, mo, D0, resp, dr, nt, memory_mb);
    r.cphf_converged = resp.cphf_converged;
    r.cphf_iterations = resp.cphf_iterations;
    r.cphf_residual = resp.cphf_residual;
    r.method_converged = so2.converged && d0.lambda_converged;
    r.method_iterations = so2.iterations;
    r.method_residual = so2.residual;
    r.orbital_seconds = orbital_seconds;
    r.method_seconds = method_seconds;
    r.seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
