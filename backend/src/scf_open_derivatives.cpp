// ΨLab 15.0: analytic gradients of the open-shell SCF energies (UHF, ROHF).
//
// Both determinants have the energy
//   E = sum_s tr(P_s h) + 1/2 sum (mn|ls) [P_mn P_ls - P^a_ml P^a_ns - P^b_ml P^b_ns] + V_nn
// (P = P_a + P_b), stationary under the orbital rotations each method allows. With the Lagrange
// multipliers of the orthonormality constraint, eps = sum_s F_s(MO) n_s (n_s the occupations of
// spin s; UHF: one orbital set per spin, ROHF: one set for both), the gradient is
//   dE/dx = sum P h^x + sum Gamma (mn|ls)^x - sum W S^x + dV_nn/dx,   W = C eps C' (symmetrised),
// which for UHF is the sum over the spins of the energy-weighted densities and for ROHF holds
// because the ROHF solution makes eps symmetric (its c-o, c-v and o-v rotations are stationary).
// With removed near-linear dependencies the rotation of the kept orbital space adds the term of
// rhf_derivatives.h with F P S replaced by sum_s F_s P_s S.
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/rhf_derivatives.h"
#include "macqc/scf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

Matrix mm(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.rows(), m = b.cols(), kk = a.cols();
    Matrix c(n, m);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t k = 0; k < kk; ++k) {
            const double aik = a(i, k);
            if (aik == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aik * b(k, j);
        }
    return c;
}

Matrix mtm(const Matrix& a, const Matrix& b) {
    const std::size_t n = a.cols(), m = b.cols(), kk = a.rows();
    Matrix c(n, m);
    for (std::size_t k = 0; k < kk; ++k)
        for (std::size_t i = 0; i < n; ++i) {
            const double aki = a(k, i);
            if (aki == 0.0) continue;
            for (std::size_t j = 0; j < m; ++j) c(i, j) += aki * b(k, j);
        }
    return c;
}

// C eps C' with eps = F_MO n (columns scaled by the occupations), symmetrised.
Matrix lagrangian_density(const Matrix& C, const Matrix& F, int nocc, Matrix& W) {
    const std::size_t n = C.rows(), nmo = C.cols();
    const Matrix Fmo = mtm(C, mm(F, C));
    Matrix eps(nmo, nmo);
    for (std::size_t p = 0; p < nmo; ++p)
        for (int q = 0; q < nocc; ++q) eps(p, static_cast<std::size_t>(q)) = Fmo(p, static_cast<std::size_t>(q));
    const Matrix Wr = mm(C, mm(eps, C.transposed()));
    if (W.rows() != n) W = Matrix(n, n);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) W(i, j) += 0.5 * (Wr(i, j) + Wr(j, i));
    return Fmo;
}

}  // namespace

TwoParticleDensity scf_two_particle_density(const RHFResult& scf) {
    TwoParticleDensity g;
    if (!scf.open_shell()) {
        SeparableDensity s;
        s.A = s.B = scf.P;
        s.coulomb = 0.5;
        s.exchange = -0.25;
        g.separable.push_back(s);
        return g;
    }
    SeparableDensity c;
    c.A = c.B = scf.P;
    c.coulomb = 0.5;
    c.exchange = 0.0;
    g.separable.push_back(c);
    for (const Matrix* Ps : {&scf.Pa, &scf.Pb}) {
        SeparableDensity x;
        x.A = x.B = *Ps;
        x.coulomb = 0.0;
        x.exchange = -0.5;
        g.separable.push_back(x);
    }
    return g;
}

Matrix scf_energy_weighted_density(const RHFResult& scf) {
    Matrix W;
    if (!scf.open_shell()) {
        // RHF: F = F_a = F_b, two electrons per occupied orbital
        const std::size_t n = scf.C.rows();
        W = Matrix(n, n);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                for (int k = 0; k < scf.n_occupied; ++k)
                    W(i, j) += 2.0 * scf.orbital_energies[static_cast<std::size_t>(k)] *
                               scf.C(i, static_cast<std::size_t>(k)) * scf.C(j, static_cast<std::size_t>(k));
        return W;
    }
    lagrangian_density(scf.C, scf.Fa, scf.n_alpha, W);
    lagrangian_density(scf.beta_orbitals(), scf.Fb, scf.n_beta, W);
    return W;
}

RHFGradientResult scf_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf, int threads,
                               double lindep_tol) {
    if (!scf.open_shell()) return rhf_gradient(mol, basis, scf, threads, lindep_tol);
    const auto t0 = Clock::now();
    const int nat = static_cast<int>(mol.atoms.size());
    const int nt = resolve_threads(threads);
    const Matrix W = scf_energy_weighted_density(scf);
    RHFGradientResult r;
    r.nuclear = nuclear_repulsion_gradient(mol);
    report_progress(JobStage::Gradient, 1);
    r.one_electron = core_hamiltonian_gradient(basis, mol, scf.P, nt);
    report_progress(JobStage::Gradient, 2);
    r.two_electron = eri_gradient(basis, nat, scf_two_particle_density(scf), nt);
    r.overlap = overlap_gradient(basis, nat, W, nt);
    for (double& v : r.overlap) v = -v;
    r.lindep.assign(static_cast<std::size_t>(3 * nat), 0.0);
    if (scf.removed_functions > 0) {
        const EigenResult se = symmetric_eigen(scf.S);
        const std::size_t N = scf.S.rows();
        std::vector<std::size_t> kept, removed;
        for (std::size_t k = 0; k < N; ++k) (se.values[k] >= lindep_tol ? kept : removed).push_back(k);
        Matrix FPS = mm(mm(scf.Fa, scf.Pa), scf.S);
        const Matrix FPSb = mm(mm(scf.Fb, scf.Pb), scf.S);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j < N; ++j) FPS(i, j) += FPSb(i, j);
        const Matrix& E = se.vectors;
        const Matrix X = mtm(E, mm(FPS, E));
        Matrix Y(N, N);
        for (std::size_t l : removed)
            for (std::size_t k : kept) {
                const double c = 2.0 * X(l, k) / (se.values[k] * (se.values[k] - se.values[l]));
                for (std::size_t m = 0; m < N; ++m)
                    for (std::size_t n = 0; n < N; ++n) Y(m, n) += c * E(m, l) * E(n, k);
            }
        Matrix Ys(N, N);
        for (std::size_t m = 0; m < N; ++m)
            for (std::size_t n = 0; n < N; ++n) Ys(m, n) = 0.5 * (Y(m, n) + Y(n, m));
        r.lindep = overlap_gradient(basis, nat, Ys, nt);
    }
    r.gradient.assign(static_cast<std::size_t>(3 * nat), 0.0);
    for (std::size_t x = 0; x < r.gradient.size(); ++x)
        r.gradient[x] = r.nuclear[x] + r.one_electron[x] + r.two_electron[x] + r.overlap[x] + r.lindep[x];
    r.seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
