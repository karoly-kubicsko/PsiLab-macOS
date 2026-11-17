#include "macqc/correlated_hessian.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

#include "macqc/ci.h"
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/job.h"
#include "macqc/parallel.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

Matrix mul(const Matrix& A, bool ta, const Matrix& B, bool tb, int threads) {
    const std::size_t m = ta ? A.cols() : A.rows();
    const std::size_t k = ta ? A.rows() : A.cols();
    const std::size_t n = tb ? B.rows() : B.cols();
    Matrix C(m, n);
    if (m && n && k) gemm(ta, tb, m, n, k, 1.0, A.data(), A.cols(), B.data(), B.cols(), 0.0, C.data(), n, threads);
    return C;
}

Matrix columns(const Matrix& C, std::size_t c0, std::size_t nc) {
    Matrix out(C.rows(), nc);
    for (std::size_t i = 0; i < C.rows(); ++i)
        for (std::size_t j = 0; j < nc; ++j) out(i, j) = C(i, c0 + j);
    return out;
}

Matrix symmetrized(const Matrix& A) {
    Matrix S(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) S(i, j) = 0.5 * (A(i, j) + A(j, i));
    return S;
}

void add(Matrix& y, double a, const Matrix& x) {
    for (std::size_t i = 0; i < y.rows(); ++i)
        for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) += a * x(i, j);
}

double dot(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    const double* pa = a.data();
    const double* pb = b.data();
    const std::size_t n = a.rows() * a.cols();
    for (std::size_t i = 0; i < n; ++i) s += pa[i] * pb[i];
    return s;
}

// X_tp = 2 sum_q h_tq D_qp + 4 sum_qrs g(t,q,r,s) G_pqrs: the orbital Lagrangian of the
// one-electron matrix h and the two-electron array g (accessor) with the densities D, G.
template <class Integrals>
Matrix orbital_lagrangian(const Matrix& h, const Integrals& g, const Matrix& D, const Packed4& G, int nt) {
    const std::size_t n = D.rows();
    Matrix X = mul(h, false, D, false, nt);
    for (std::size_t t = 0; t < n; ++t)
        for (std::size_t p = 0; p < n; ++p) X(t, p) *= 2.0;
    const std::size_t np = n * (n + 1) / 2;
    std::vector<double> A(n * np), B(n * np), Xq(n * n);
    for (std::size_t q = 0; q < n; ++q) {
        check_cancelled();
        for (std::size_t t = 0; t < n; ++t)
            for (std::size_t r = 0, rs = 0; r < n; ++r)
                for (std::size_t s = 0; s <= r; ++s, ++rs) {
                    const double w = r == s ? 1.0 : 2.0;
                    A[t * np + rs] = g(t, q, r, s);
                    B[t * np + rs] = w * G.get(t, q, r, s);
                }
        gemm(false, true, n, n, np, 4.0, A.data(), np, B.data(), np, 0.0, Xq.data(), n, nt);
        for (std::size_t t = 0; t < n; ++t)
            for (std::size_t p = 0; p < n; ++p) X(t, p) += Xq[t * n + p];
    }
    return X;
}

// sum_pqrs G_pqrs g_pqrs over all index values (G, g 8-fold symmetric).
template <class Integrals>
double contract_full(const Packed4& G, const Integrals& g) {
    const std::size_t n = G.n;
    double e = 0.0;
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q <= p; ++q)
            for (std::size_t r = 0; r <= p; ++r)
                for (std::size_t s = 0; s <= (r == p ? q : r); ++s)
                    e += Packed4::degeneracy(p, q, r, s) * G.get(p, q, r, s) * g(p, q, r, s);
    return e;
}

// out[y] += sum_b w_b A_b U[y] B_b' over nb blocks of n x n matrices (12.0: the products of
// all blocks at once, in chunks of at most ~16 M doubles): A_b[t][u] = fa(b, t, u),
// B_b[p][v] = fb(b, p, v). With the blocks stacked, As[t][(b,u)] and Bs[(b,p)][v], for each y
// D[(b,p)][u] = sum_v Bs[(b,p)][v] U[u][v], E[(b,u)][p] = w_b D[(b,p)][u] and out[y] += As E.
template <class FA, class FB, class FW>
void add_sandwiches(std::size_t nb, std::size_t n, FA fa, FB fb, FW w, const std::vector<Matrix>& U,
                    std::vector<Matrix>& out, int nt) {
    if (nb == 0 || U.empty()) return;
    const std::size_t chunk = std::max<std::size_t>(1, std::min<std::size_t>(nb, (std::size_t{16} << 20) / (n * n)));
    std::vector<double> As(n * chunk * n), Bs(chunk * n * n), D(chunk * n * n), E(chunk * n * n);
    for (std::size_t b0 = 0; b0 < nb; b0 += chunk) {
        check_cancelled();
        const std::size_t nbc = std::min(chunk, nb - b0), wcols = nbc * n;
        for (std::size_t b = 0; b < nbc; ++b)
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t u = 0; u < n; ++u) {
                    As[t * wcols + b * n + u] = fa(b0 + b, t, u);
                    Bs[(b * n + t) * n + u] = fb(b0 + b, t, u);
                }
        for (std::size_t y = 0; y < U.size(); ++y) {
            gemm(false, true, nbc * n, n, n, 1.0, Bs.data(), n, U[y].data(), n, 0.0, D.data(), n, nt);
            for (std::size_t b = 0; b < nbc; ++b) {
                const double wb = w(b0 + b);
                for (std::size_t p = 0; p < n; ++p)
                    for (std::size_t u = 0; u < n; ++u) E[(b * n + u) * n + p] = wb * D[(b * n + p) * n + u];
            }
            gemm(false, false, n, n, wcols, 1.0, As.data(), wcols, E.data(), n, 1.0, out[y].data(), n, nt);
        }
    }
}

}  // namespace

std::vector<double> transform_packed(const std::vector<double>& in, std::size_t ns, const Matrix& M, int threads) {
    const int nt = resolve_threads(threads);
    const std::size_t nd = M.cols();
    const std::size_t sp = ns * (ns + 1) / 2, dp = nd * (nd + 1) / 2;
    auto pair = [](std::size_t a, std::size_t b) { return a >= b ? a * (a + 1) / 2 + b : b * (b + 1) / 2 + a; };
    auto get = [&](std::size_t a, std::size_t b, std::size_t c, std::size_t d) { return in[pair(pair(a, b), pair(c, d))]; };
    // half[(pq)][(cd)] = sum_ab M_ap M_bq in_abcd
    std::vector<double> half(dp * sp, 0.0);
    {
        Matrix S(ns, ns);
        for (std::size_t c = 0, cd = 0; c < ns; ++c)
            for (std::size_t d = 0; d <= c; ++d, ++cd) {
                check_cancelled();
                for (std::size_t a = 0; a < ns; ++a)
                    for (std::size_t b = 0; b <= a; ++b) S(a, b) = S(b, a) = get(a, b, c, d);
                const Matrix A = mul(M, true, mul(S, false, M, false, nt), false, nt);
                for (std::size_t p = 0, pq = 0; p < nd; ++p)
                    for (std::size_t q = 0; q <= p; ++q, ++pq) half[pq * sp + cd] = A(p, q);
            }
    }
    std::vector<double> out(dp * (dp + 1) / 2, 0.0);
    Matrix S(ns, ns);
    for (std::size_t p = 0, pq = 0; p < nd; ++p)
        for (std::size_t q = 0; q <= p; ++q, ++pq) {
            check_cancelled();
            const double* h = half.data() + pq * sp;
            for (std::size_t c = 0, cd = 0; c < ns; ++c)
                for (std::size_t d = 0; d <= c; ++d, ++cd) S(c, d) = S(d, c) = h[cd];
            const Matrix A = mul(M, true, mul(S, false, M, false, nt), false, nt);
            for (std::size_t r = 0, rs = 0; r < nd && rs <= pq; ++r)
                for (std::size_t s = 0; s <= r && rs <= pq; ++s, ++rs) out[pq * (pq + 1) / 2 + rs] = A(r, s);
        }
    return out;
}

double orbital_response_bytes(std::size_t n, std::size_t N, int n_atoms) {
    const double n3 = 3.0 * n_atoms;
    const double np = static_cast<double>(n) * (n + 1) / 2.0, NP = static_cast<double>(N) * (N + 1) / 2.0;
    const double mo = n3 * np * (np + 1) / 2.0;                 // g~^x in MOs
    const double ao = n3 * NP * (NP + 1) / 2.0;                 // AO derivative integrals (transient)
    const double half = NP * np;                                // transformation
    const double mats = n3 * (5.0 * static_cast<double>(n) * n + static_cast<double>(N) * N);
    return 8.0 * (mo + std::max(ao, half) + mats);
}

OrbitalResponse rhf_orbital_response(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf, int threads,
                                     double memory_mb, double cphf_tol, const CPHFCallback& progress, bool fields) {
    if (rhf.removed_functions > 0)
        throw std::runtime_error("analytic Hessian with removed near-linear dependencies");
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t n3 = static_cast<std::size_t>(3 * nat);
    const Matrix& C = rhf.C;
    const std::size_t N = C.rows(), n = C.cols(), o = static_cast<std::size_t>(rhf.n_occupied), v = n - o;
    const std::vector<double>& eps = rhf.orbital_energies;
    const Matrix Co = columns(C, 0, o), Cv = columns(C, o, v);
    Matrix P = mul(Co, false, Co, true, nt);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) P(i, j) *= 2.0;

    OrbitalResponse r;
    r.n3 = n3;
    r.n_fields = fields ? 3 : 0;
    const std::size_t nt3 = n3 + r.n_fields;  // perturbations: coordinates, then the field components
    report_progress(JobStage::Hessian, 1);
    const std::vector<Matrix> Sx = overlap_derivative_matrices(basis, nat, nt);
    const std::vector<Matrix> hx = core_hamiltonian_derivative_matrices(basis, mol, nt);
    const std::vector<Matrix> GxP = eri_derivative_fock(basis, nat, P, nt);
    std::vector<Matrix> Fexpl(nt3);
    r.S.resize(nt3);
    r.h.resize(nt3);
    for (std::size_t x = 0; x < n3; ++x) {
        r.S[x] = transform_matrix(C, Sx[x], C, nt);
        r.h[x] = transform_matrix(C, hx[x], C, nt);
        Matrix F = hx[x];
        add(F, 1.0, GxP[x]);
        Fexpl[x] = transform_matrix(C, F, C, nt);
    }
    if (fields) {
        // a uniform field F adds F . r to the one-electron Hamiltonian (origin 0): h~^F = r (MO),
        // no overlap or two-electron derivative
        const std::vector<Matrix> rd = dipole_matrices(basis, {0.0, 0.0, 0.0});
        for (std::size_t f = 0; f < 3; ++f) {
            r.S[n3 + f] = Matrix(n, n);
            r.h[n3 + f] = transform_matrix(C, rd[f], C, nt);
            Fexpl[n3 + f] = r.h[n3 + f];
        }
    }
    {
        std::vector<std::vector<double>> gao = eri_derivative_packed(basis, nat, nt);
        r.g.resize(nt3);  // the field components: empty (zero)
        for (std::size_t x = 0; x < n3; ++x) {
            report_progress(JobStage::Hessian, 1, 0.0, static_cast<double>(x + 1));
            Packed4 g;
            g.n = n;
            g.v = transform_packed(gao[x], N, C, nt);
            r.g[x] = std::move(g);
            std::vector<double>().swap(gao[x]);
        }
    }
    r.integral_seconds = seconds_since(t0);

    // ---- coupled-perturbed HF for every coordinate (as rhf_hessian)
    const auto t1 = Clock::now();
    FockBuilder fock(basis, nt, memory_mb);
    std::vector<Matrix> rhs(nt3);
    for (std::size_t x = 0; x < nt3; ++x) {
        check_cancelled();
        if (x >= n3) {
            Matrix B(v, o);
            for (std::size_t a = 0; a < v; ++a)
                for (std::size_t i = 0; i < o; ++i) B(a, i) = -Fexpl[x](o + a, i);
            rhs[x] = B;
            continue;
        }
        Matrix Soo(o, o);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j) Soo(i, j) = r.S[x](i, j);
        const Matrix GS = fock.build(mul(mul(Co, false, Soo, false, nt), false, Co, true, nt));
        const Matrix GSvo = transform_matrix(Cv, GS, Co, nt);
        Matrix B(v, o);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i)
                B(a, i) = -Fexpl[x](o + a, i) + eps[i] * r.S[x](o + a, i) + 2.0 * GSvo(a, i);
        rhs[x] = B;
    }
    report_progress(JobStage::CPHF, 0);
    const CPHFSolution sol = solve_cphf(fock, C, eps, static_cast<int>(o), rhs, cphf_tol, 80, progress, nt);
    r.cphf_converged = sol.converged;
    r.cphf_iterations = sol.iterations;
    r.cphf_residual = sol.residual;

    r.U.resize(nt3);
    r.P.resize(nt3);
    r.F.resize(nt3);
    for (std::size_t x = 0; x < nt3; ++x) {
        check_cancelled();
        const Matrix& S = r.S[x];
        Matrix U(n, n);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q)
                if ((p < o) == (q < o)) U(p, q) = -0.5 * S(p, q);
        for (std::size_t a = 0; a < v; ++a)
            for (std::size_t i = 0; i < o; ++i) {
                U(o + a, i) = sol.U[x](a, i);
                U(i, o + a) = -S(i, o + a) - sol.U[x](a, i);
            }
        // P^x = 2 (C U_o Co' + Co U_o' C'), U_o the occupied columns of U
        const Matrix CUo = mul(C, false, columns(U, 0, o), false, nt);
        Matrix Px = mul(CUo, false, Co, true, nt);
        for (std::size_t i = 0; i < N; ++i)
            for (std::size_t j = 0; j <= i; ++j) Px(i, j) = Px(j, i) = 2.0 * (Px(i, j) + Px(j, i));
        Matrix F = Fexpl[x];
        add(F, 1.0, transform_matrix(C, fock.build(Px), C, nt));
        r.U[x] = U;
        r.P[x] = Px;
        r.F[x] = symmetrized(F);
    }
    r.cphf_seconds = seconds_since(t1);
    return r;
}

namespace {

// H^(x) from the orbital response U, the skeleton one-electron integrals hs and two-electron
// integrals gs of one perturbation.
MOIntegrals perturbed_from(const MOIntegrals& mo, const Matrix& U, const Matrix& hs, const Packed4& gx, int nt) {
    const std::size_t n = static_cast<std::size_t>(mo.n);
    MOIntegrals out;
    out.n = mo.n;
    out.h = hs;
    add(out.h, 1.0, mul(U, true, mo.h, false, nt));
    add(out.h, 1.0, mul(mo.h, false, U, false, nt));
    out.h = symmetrized(out.h);
    out.npair = n * (n + 1) / 2;
    out.stride = out.npair;
    // M[(pq)][(rs)] = sum_t (U_tp g_tqrs + U_tq g_ptrs) for every rs; g^(x) = g~ + M + M'
    // (12.0: B_rs = U' G_rs for a chunk of rs pairs in one product, GG[t][(rs, q)] = (tq|rs))
    std::vector<double> M(out.npair * out.npair, 0.0);
    {
        const std::size_t np = out.npair;
        const std::size_t chunk = std::max<std::size_t>(1, std::min<std::size_t>(np, (std::size_t{8} << 20) / (n * n)));
        std::vector<double> GG(n * chunk * n), BB(n * chunk * n);
        std::vector<std::size_t> rr(np), ss(np);
        for (std::size_t r_ = 0, rs = 0; r_ < n; ++r_)
            for (std::size_t s = 0; s <= r_; ++s, ++rs) {
                rr[rs] = r_;
                ss[rs] = s;
            }
        for (std::size_t b0 = 0; b0 < np; b0 += chunk) {
            const std::size_t nbc = std::min(chunk, np - b0), wcols = nbc * n;
            for (std::size_t b = 0; b < nbc; ++b)
                for (std::size_t t = 0; t < n; ++t)
                    for (std::size_t q = 0; q < n; ++q) GG[t * wcols + b * n + q] = mo(t, q, rr[b0 + b], ss[b0 + b]);
            gemm(true, false, n, wcols, n, 1.0, U.data(), n, GG.data(), wcols, 0.0, BB.data(), wcols, nt);
            for (std::size_t b = 0; b < nbc; ++b)
                for (std::size_t p = 0, pq = 0; p < n; ++p)
                    for (std::size_t q = 0; q <= p; ++q, ++pq)
                        M[pq * np + b0 + b] = BB[p * wcols + b * n + q] + BB[q * wcols + b * n + p];
        }
    }
    out.eri.assign(out.npair * out.npair, 0.0);
    for (std::size_t p = 0, pq = 0; p < n; ++p)
        for (std::size_t q = 0; q <= p; ++q, ++pq)
            for (std::size_t r_ = 0, rs = 0; r_ < n; ++r_)
                for (std::size_t s = 0; s <= r_; ++s, ++rs)
                    out.eri[pq * out.npair + rs] =
                        (gx.v.empty() ? 0.0 : gx.get(p, q, r_, s)) + M[pq * out.npair + rs] + M[rs * out.npair + pq];
    return out;
}

}  // namespace

MOIntegrals perturbed_hamiltonian(const MOIntegrals& mo, const OrbitalResponse& r, std::size_t x, int threads) {
    return perturbed_from(mo, r.U[x], r.h[x], r.g[x], resolve_threads(threads));
}

MOIntegrals perturbed_hamiltonian(const MOIntegrals& mo, const OrbitalResponse& r, const std::vector<double>& d,
                                  int threads) {
    if (d.size() != r.n3) throw std::invalid_argument("perturbed_hamiltonian: direction of another molecule");
    const std::size_t n = static_cast<std::size_t>(mo.n);
    Matrix U(n, n), hs(n, n);
    Packed4 gs(n);
    for (std::size_t x = 0; x < r.n3; ++x) {
        if (d[x] == 0.0) continue;
        add(U, d[x], r.U[x]);
        add(hs, d[x], r.h[x]);
        const std::vector<double>& g = r.g[x].v;
        for (std::size_t k = 0; k < g.size(); ++k) gs.v[k] += d[x] * g[k];
    }
    return perturbed_from(mo, U, hs, gs, resolve_threads(threads));
}

CorrelatedHessianResult correlated_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                           const MOIntegrals& mo, const CorrelatedDensity& dens,
                                           const OrbitalResponse& resp, const DensityResponse& dresp, int threads,
                                           double memory_mb) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t n3 = static_cast<std::size_t>(3 * nat);
    const Matrix& C = rhf.C;
    const std::size_t N = C.rows(), n = C.cols(), o = static_cast<std::size_t>(rhf.n_occupied), v = n - o;
    const std::vector<double>& eps = rhf.orbital_energies;
    if (resp.n3 != n3 || resp.U.size() != n3 + resp.n_fields)
        throw std::invalid_argument("correlated_hessian: orbital response of another molecule");
    const bool internal = dresp.internal;
    const bool has_response = internal || !dresp.D.empty();
    const std::size_t m = internal ? dresp.directions.size() : (has_response ? n3 : 0);
    // the field components (dipole derivatives): internal form only; in the 2n+2 form the method
    // supplies their responses after those of the directions
    const std::size_t nf = internal ? resp.n_fields : 0;
    const std::size_t mr = dresp.symmetric ? m + nf : m;  // responses supplied
    if (has_response && (dresp.D.size() != mr || dresp.G.size() != mr))
        throw std::invalid_argument("correlated_hessian: density responses for every coordinate or direction needed");
    if (dresp.symmetric && (!internal || dresp.pair.rows() != mr || dresp.pair.cols() != mr))
        throw std::invalid_argument("correlated_hessian: the 2n+2 form needs internal directions and the pair terms");
    for (const auto& d : dresp.directions)
        if (d.size() != n3) throw std::invalid_argument("correlated_hessian: direction of another molecule");
    const Matrix Co = columns(C, 0, o);
    Matrix P = mul(Co, false, Co, true, nt);
    for (std::size_t i = 0; i < N; ++i)
        for (std::size_t j = 0; j < N; ++j) P(i, j) *= 2.0;
    CorrelatedHessianResult res;
    const Matrix Ds = symmetrized(dens.D);
    const Packed4& G = dens.G;

    // ---- orbital Lagrangian, Z-vector, relaxed and energy-weighted densities (as the gradient)
    report_progress(JobStage::Hessian, 2);
    const Matrix X = orbital_lagrangian(mo.h, mo, Ds, G, nt);
    FockBuilder fock(basis, nt, memory_mb);
    Matrix Yv(v, o);
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t i = 0; i < o; ++i) Yv(a, i) = X(o + a, i) - X(i, o + a);
    const CPHFSolution zs = solve_cphf(fock, C, eps, static_cast<int>(o), {Yv}, 1.0e-10, 80, CPHFCallback(), nt);
    res.zvector_converged = zs.converged;
    res.zvector_iterations = zs.iterations;
    res.zvector_residual = zs.residual;
    const Matrix& Z = zs.U[0];
    Matrix Zf(n, n), DZmo(n, n);  // Zf_ai = Z_ai (a virtual row, i occupied column)
    for (std::size_t a = 0; a < v; ++a)
        for (std::size_t i = 0; i < o; ++i) {
            Zf(o + a, i) = Z(a, i);
            DZmo(o + a, i) = DZmo(i, o + a) = 0.5 * Z(a, i);
        }
    const Matrix DZ = symmetrized(mul(mul(C, false, DZmo, false, nt), false, C, true, nt));  // AO
    const Matrix Gz = symmetrized(transform_matrix(C, fock.build(DZ), C, nt));
    // W (MO, symmetric): coefficient of -S~^xy_pq (and of the products of the second-order
    // orthonormality condition)
    Matrix W(n, n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q) {
            const bool po = p < o, qo = q < o;
            if (po && qo) W(p, q) = 0.25 * (X(p, q) + X(q, p)) - 2.0 * Gz(p, q);
            else if (!po && !qo) W(p, q) = 0.25 * (X(p, q) + X(q, p));
            else {
                const std::size_t i = po ? p : q, a = po ? q : p;  // a: MO index of the virtual
                W(p, q) = 0.5 * X(i, a) - 0.5 * Z(a - o, i) * eps[i];
            }
        }

    // ---- explicit second derivatives with the relaxed densities
    std::vector<std::array<double, 3>> dip_expl(n3, {0.0, 0.0, 0.0});  // d2/dx dF of the explicit terms
    {
        const auto te = Clock::now();
        Matrix D1 = mul(mul(C, false, Ds, false, nt), false, C, true, nt);
        add(D1, -1.0, DZ);
        D1 = symmetrized(D1);
        const std::vector<double> gao = transform_packed(G.v, n, C.transposed(), nt);
        TwoParticleDensity gamma;
        gamma.packed = &gao;
        SeparableDensity zterm;
        zterm.A = DZ;
        zterm.B = P;
        zterm.coulomb = -1.0;
        zterm.exchange = 0.5;
        gamma.separable.push_back(zterm);
        report_progress(JobStage::Hessian, 3);
        res.hessian = core_hamiltonian_hessian(basis, mol, D1, nt);
        add(res.hessian, 1.0, eri_hessian(basis, nat, gamma, nt));
        const Matrix Wao = symmetrized(mul(mul(C, false, W, false, nt), false, C, true, nt));
        add(res.hessian, -1.0, overlap_hessian(basis, nat, Wao, nt));
        add(res.hessian, 1.0, nuclear_repulsion_hessian(mol));
        if (nf > 0) {
            // field: h(F) = h + F . r, V_nn(F) = V_nn - F . sum Z R
            const std::vector<std::vector<Matrix>> rx = dipole_derivative_matrices(basis, nat, {0.0, 0.0, 0.0}, nt);
            for (std::size_t x = 0; x < n3; ++x)
                for (std::size_t f = 0; f < 3; ++f) {
                    dip_expl[x][f] = dot(D1, rx[x][f]);
                    if (x % 3 == f) dip_expl[x][f] -= mol.atoms[x / 3].Z;
                }
        }
        res.explicit_seconds = seconds_since(te);
    }
    Matrix& H = res.hessian;

    // ---- response terms, bilinear in first-order quantities of two perturbations p, p': the
    // Cartesian coordinates, or (internal) the directions q = sum_x Q_xq x, with the
    // first-order quantities combined accordingly
    const auto tr = Clock::now();
    const std::size_t ng = internal ? m : n3;  // geometric perturbations, then the field components
    const std::size_t np = ng + nf;
    std::vector<std::vector<double>> dir(np, std::vector<double>(n3, 0.0));  // zero for a field
    for (std::size_t p = 0; p < ng; ++p) {
        if (internal) dir[p] = dresp.directions[p];
        else dir[p][p] = 1.0;
    }
    auto combo = [&](const std::vector<Matrix>& v, std::size_t p) {
        if (p >= ng) return v[n3 + (p - ng)];
        if (!internal) return v[p];
        Matrix out(v[0].rows(), v[0].cols());
        for (std::size_t x = 0; x < n3; ++x)
            if (dir[p][x] != 0.0) add(out, dir[p][x], v[x]);
        return out;
    };
    std::vector<Matrix> U(np), S(np), F(np), Pp(np), hs(np);
    for (std::size_t p = 0; p < np; ++p) {
        U[p] = combo(resp.U, p);
        S[p] = combo(resp.S, p);
        F[p] = combo(resp.F, p);
        Pp[p] = combo(resp.P, p);
        hs[p] = combo(resp.h, p);
    }
    // (i) sum U^p' X~^[p] + U^p X~^[p'], X~^[p] the orbital Lagrangian of the skeleton integrals
    std::vector<Matrix> Xs(np);
    for (std::size_t p = 0; p < np; ++p) {
        report_progress(JobStage::Hessian, 4, 0.0, static_cast<double>(p + 1));
        Packed4 combined;  // empty for a field (no two-electron part)
        if (internal && p < ng) {
            combined = Packed4(n);
            for (std::size_t x = 0; x < n3; ++x) {
                if (dir[p][x] == 0.0) continue;
                const std::vector<double>& g = resp.g[x].v;
                for (std::size_t k = 0; k < g.size(); ++k) combined.v[k] += dir[p][x] * g[k];
            }
        }
        const Packed4& gx = (internal || p >= ng) ? combined : resp.g[p];
        Xs[p] = orbital_lagrangian(
            hs[p], [&](std::size_t a, std::size_t b, std::size_t c, std::size_t d) { return gx.v.empty() ? 0.0 : gx.get(a, b, c, d); },
            Ds, G, nt);
    }
    // (ii) M^p_tp = 4 sum U^p_uq g_turs G_pqrs + 8 sum U^p_ur g_tqus G_pqrs (double one-index
    // transformations of the two-electron term) and 2 sum D (U^p' h U^p')
    std::vector<Matrix> My(np, Matrix(n, n));
    {
        // 12.0: all blocks in a few large products (add_sandwiches). First sum over r >= s
        // (block b = pair index rs), A_b[t][u] = (tu|rs), B_b[p][v] = G_pvrs, weight 4 (8 for r != s);
        // then over all (q, s), A_b[t][u] = (tq|us), B_b[p][v] = G_pqvs, weight 8.
        std::vector<std::size_t> rr, ss;
        for (std::size_t r_ = 0; r_ < n; ++r_)
            for (std::size_t s = 0; s <= r_; ++s) {
                rr.push_back(r_);
                ss.push_back(s);
            }
        add_sandwiches(
            rr.size(), n, [&](std::size_t b, std::size_t t, std::size_t u) { return mo(t, u, rr[b], ss[b]); },
            [&](std::size_t b, std::size_t p, std::size_t q) { return G.get(p, q, rr[b], ss[b]); },
            [&](std::size_t b) { return (rr[b] == ss[b] ? 1.0 : 2.0) * 4.0; }, U, My, nt);
        add_sandwiches(
            n * n, n, [&](std::size_t b, std::size_t t, std::size_t u) { return mo(t, b / n, u, b % n); },
            [&](std::size_t b, std::size_t p, std::size_t v) { return G.get(p, b / n, v, b % n); },
            [&](std::size_t) { return 8.0; }, U, My, nt);
    }
    std::vector<Matrix> hU(np), DUx(np);  // h U^p, U^p D
    for (std::size_t x = 0; x < np; ++x) {
        hU[x] = mul(mo.h, false, U[x], false, nt);
        DUx[x] = mul(U[x], false, Ds, false, nt);
    }
    // (iii) products of the orthonormality condition: -2 [sum (S^y W) o U^x + (S^x W) o U^y + (U^x W) o U^y]
    std::vector<Matrix> SW(np), UW(np);
    for (std::size_t x = 0; x < np; ++x) {
        SW[x] = mul(S[x], false, W, false, nt);
        UW[x] = mul(U[x], false, W, false, nt);
    }
    // (iv) Z-vector terms of the first-order products: FZ[y] = F^y Zf' + F^y Zf + eps (U^y Zf')
    std::vector<Matrix> FZ(np);
    for (std::size_t y = 0; y < np; ++y) {
        Matrix Am = mul(F[y], false, Zf, true, nt);
        add(Am, 1.0, mul(F[y], false, Zf, false, nt));
        const Matrix UZ = mul(U[y], false, Zf, true, nt);
        for (std::size_t t = 0; t < n; ++t)
            for (std::size_t p = 0; p < n; ++p) Am(t, p) += eps[t] * UZ(t, p);
        FZ[y] = Am;
    }
    // (v) sum DZ G~^x[P^y] for every x and y (one derivative-integral pass per y), projected
    std::vector<std::vector<double>> VZ(np, std::vector<double>(np, 0.0));
    {
        // 12.0: one pass over the derivative integrals for all perturbations
        report_progress(JobStage::Hessian, 5, 0.0, static_cast<double>(np));
        std::vector<TwoParticleDensity> gzs(np);
        for (std::size_t y = 0; y < np; ++y) {
            SeparableDensity t;
            t.A = DZ;
            t.B = Pp[y];
            t.coulomb = 1.0;
            t.exchange = -0.5;
            gzs[y].separable.push_back(t);
        }
        const std::vector<std::vector<double>> vzs = eri_gradients(basis, nat, gzs, nt);
        for (std::size_t y = 0; y < np; ++y)
            for (std::size_t x = 0; x < np; ++x)
                for (std::size_t c = 0; c < n3; ++c) VZ[y][x] += dir[x][c] * vzs[y][c];
    }
    // (vi) -4 sum_k (U^x' Gz U^y)_kk
    std::vector<Matrix> GU(np);
    for (std::size_t y = 0; y < np; ++y) GU[y] = mul(Gz, false, U[y], false, nt);

    Matrix Hb(np, np);
    for (std::size_t x = 0; x < np; ++x)
        for (std::size_t y = 0; y < np; ++y) {
            double s = dot(U[y], Xs[x]) + dot(U[x], Xs[y]);                 // (i)
            s += dot(U[x], My[y]);                                            // (ii) two-electron
            s += 2.0 * dot(DUx[x], hU[y]);                                    // (ii) 2 sum D U^x' h U^y
            s -= 2.0 * (dot(SW[y], U[x]) + dot(SW[x], U[y]) + dot(UW[x], U[y]));  // (iii)
            s -= dot(U[x], FZ[y]) + dot(U[y], FZ[x]);                         // (iv)
            s -= VZ[y][x] + VZ[x][y];                                         // (v)
            double k = 0.0;                                                   // (vi)
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t c = 0; c < o; ++c) k += U[x](t, c) * GU[y](t, c);
            s -= 4.0 * k;
            Hb(x, y) = s;
        }

    // the gradient <D, H^(x)> + V_nn^x = V_nn^x + <D, h~^x> + <G, g~^x> + U^x . X (X: the
    // orbital Lagrangian, the derivative of <D, H> with respect to the orbital rotation)
    res.gradient = nuclear_repulsion_gradient(mol);
    for (std::size_t x = 0; x < n3; ++x) {
        const Packed4& gx = resp.g[x];
        res.gradient[x] += dot(Ds, resp.h[x]) + dot(resp.U[x], X) +
                           contract_full(G, [&](std::size_t a, std::size_t b, std::size_t c, std::size_t d) { return gx.get(a, b, c, d); });
    }
    // the relaxed dipole moment -dE/dF = sum Z R - <D, H^(F)> (origin 0)
    if (nf > 0) {
        for (std::size_t f = 0; f < 3; ++f) {
            double mu = -(dot(Ds, resp.h[n3 + f]) + dot(resp.U[n3 + f], X));
            for (const Atom& a : mol.atoms) mu += a.Z * a.position[f];
            res.dipole[f] = mu;
        }
    }
    // (vii) the response of the method's densities, A(p, k) = <D^(k), H^(p)>
    Matrix A(np, mr);
    {
        std::vector<Matrix> Dk(mr);
        for (std::size_t k = 0; k < mr; ++k) Dk[k] = symmetrized(dresp.D[k]);
        if (mr > 0)
            for (std::size_t p = 0; p < np; ++p) {
                report_progress(JobStage::Hessian, 6, 0.0, static_cast<double>(p + 1));
                const MOIntegrals Hp = p >= ng ? perturbed_hamiltonian(mo, resp, n3 + (p - ng), nt)
                                       : internal ? perturbed_hamiltonian(mo, resp, dir[p], nt)
                                                  : perturbed_hamiltonian(mo, resp, p, nt);
                for (std::size_t k = 0; k < mr; ++k) A(p, k) = dot(Dk[k], Hp.h) + contract_full(dresp.G[k], Hp);
            }
    }
    res.perturbations = m;
    res.response_seconds = seconds_since(tr);
    if (!internal) {
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t y = 0; y < n3; ++y) H(x, y) += Hb(x, y) + (m > 0 ? A(x, y) : 0.0);
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t y = 0; y < x; ++y) {
                res.asymmetry = std::max(res.asymmetry, std::fabs(H(x, y) - H(y, x)));
                H(x, y) = H(y, x) = 0.5 * (H(x, y) + H(y, x));
            }
        res.seconds = seconds_since(t0);
        return res;
    }
    // internal block B = Q' H Q: the explicit terms projected, the bilinear terms and the response
    const std::vector<std::vector<double>>& Q = dresp.directions;
    Matrix B(m, m);
    for (std::size_t k = 0; k < m; ++k)
        for (std::size_t l = 0; l < m; ++l) {
            double b = 0.0;
            for (std::size_t x = 0; x < n3; ++x) {
                if (Q[k][x] == 0.0) continue;
                for (std::size_t y = 0; y < n3; ++y) b += Q[k][x] * 0.5 * (H(x, y) + H(y, x)) * Q[l][y];
            }
            B(k, l) = b + Hb(k, l) + (dresp.symmetric ? A(k, l) + A(l, k) + dresp.pair(k, l) : A(k, l));
        }
    for (std::size_t k = 0; k < m; ++k)
        for (std::size_t l = 0; l < k; ++l) {
            res.asymmetry = std::max(res.asymmetry, std::fabs(B(k, l) - B(l, k)));
            B(k, l) = B(l, k) = 0.5 * (B(k, l) + B(l, k));
        }
    // H = sum_k (H q_k) q_k' + sum_r (H r) r' over the complete orthonormal basis {q_k, r},
    // H q_k = sum_l q_l B_lk + sum_r r (H r . q_k), H r from the invariance
    const RigidMotion rm = rigid_motions(mol, {0.0, 0.0, 0.0}, &res.gradient);
    Matrix Hf(n3, n3);
    for (std::size_t k = 0; k < m; ++k) {
        std::vector<double> hq(n3, 0.0);
        for (std::size_t l = 0; l < m; ++l)
            for (std::size_t x = 0; x < n3; ++x) hq[x] += Q[l][x] * B(l, k);
        for (std::size_t r = 0; r < rm.directions.size(); ++r) {
            double c = 0.0;
            for (std::size_t x = 0; x < n3; ++x) c += rm.hessian_times[r][x] * Q[k][x];
            for (std::size_t x = 0; x < n3; ++x) hq[x] += c * rm.directions[r][x];
        }
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t y = 0; y < n3; ++y) Hf(x, y) += hq[x] * Q[k][y];
    }
    for (std::size_t r = 0; r < rm.directions.size(); ++r)
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t y = 0; y < n3; ++y) Hf(x, y) += rm.hessian_times[r][x] * rm.directions[r][y];
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t y = 0; y < x; ++y) Hf(x, y) = Hf(y, x) = 0.5 * (Hf(x, y) + Hf(y, x));
    res.hessian = Hf;
    // atomic polar tensor: d mu / dq = -d2E/dq dF for the internal directions, completed with the
    // translations and rotations (from the dipole moment)
    if (nf > 0) {
        std::vector<std::array<double, 3>> dq(m, {0.0, 0.0, 0.0});
        for (std::size_t k = 0; k < m; ++k)
            for (std::size_t f = 0; f < 3; ++f) {
                const std::size_t pf = ng + f;
                double e = 0.5 * (Hb(k, pf) + Hb(pf, k));
                for (std::size_t x = 0; x < n3; ++x) e += Q[k][x] * dip_expl[x][f];
                e += dresp.symmetric ? A(k, pf) + A(pf, k) + dresp.pair(k, pf) : A(pf, k);
                dq[k][f] = -e;
            }
        const RigidMotion rd = rigid_motions(mol, res.dipole, &res.gradient);
        res.dipole_derivatives.assign(n3, {0.0, 0.0, 0.0});
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t f = 0; f < 3; ++f) {
                double s = 0.0;
                for (std::size_t k = 0; k < m; ++k) s += Q[k][x] * dq[k][f];
                for (std::size_t r = 0; r < rd.directions.size(); ++r) s += rd.directions[r][x] * rd.dipole_derivatives[r][f];
                res.dipole_derivatives[x][f] = s;
            }
        res.has_dipole_derivatives = true;
    }
    res.seconds = seconds_since(t0);
    return res;
}

CorrelatedHessianResult fci_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                    const MOIntegrals& mo, const std::vector<double>& c, double E, int threads,
                                    double memory_mb, double fci_memory_mb, bool dipoles) {
    const int nt = resolve_threads(threads);
    const int nel = mol.electrons();
    if (nel % 2 != 0) throw std::invalid_argument("fci_hessian: closed shell only");
    const int na = nel / 2;
    const std::size_t n = static_cast<std::size_t>(mo.n);
    const auto t0 = Clock::now();
    const OrbitalResponse resp = rhf_orbital_response(mol, basis, rhf, nt, memory_mb, 1.0e-10, CPHFCallback(), dipoles);
    const double orbital_seconds = seconds_since(t0);
    const auto t1 = Clock::now();
    DensityResponse dr;
    dr.internal = true;
    dr.directions = vibrational_directions(mol);
    const std::size_t m = dr.directions.size();
    const FCIResponse fr = fci_vector_responses(
        mo, na, na, c, E, m, [&](std::size_t k) { return perturbed_hamiltonian(mo, resp, dr.directions[k], nt); }, nt,
        fci_memory_mb, 1.0e-10);
    const CorrelatedDensity dens = fci_densities(mo.n, na, na, c, nt);
    dr.D.resize(m);
    dr.G.resize(m);
    std::vector<double> cp(c.size()), cm(c.size());
    for (std::size_t x = 0; x < m; ++x) {
        check_cancelled();
        const std::vector<double>& y = fr.vectors[x];
        for (std::size_t i = 0; i < c.size(); ++i) {
            cp[i] = c[i] + y[i];
            cm[i] = c[i] - y[i];
        }
        const CorrelatedDensity dp = fci_densities(mo.n, na, na, cp, nt);
        const CorrelatedDensity dm = fci_densities(mo.n, na, na, cm, nt);
        Matrix D(n, n);
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) D(p, q) = 0.5 * (dp.D(p, q) - dm.D(p, q));
        Packed4 G(n);
        for (std::size_t k = 0; k < G.v.size(); ++k) G.v[k] = 0.5 * (dp.G.v[k] - dm.G.v[k]);
        dr.D[x] = std::move(D);
        dr.G[x] = std::move(G);
    }
    const double method_seconds = seconds_since(t1);
    CorrelatedHessianResult r = correlated_hessian(mol, basis, rhf, mo, dens, resp, dr, nt, memory_mb);
    r.cphf_converged = resp.cphf_converged;
    r.cphf_iterations = resp.cphf_iterations;
    r.cphf_residual = resp.cphf_residual;
    r.method_converged = fr.converged;
    r.method_iterations = fr.iterations;
    r.method_residual = fr.residual;
    r.orbital_seconds = orbital_seconds;
    r.method_seconds = method_seconds;
    r.seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
