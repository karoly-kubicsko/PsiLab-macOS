// ΨLab 15.0: analytic gradients and Hessians of the spin-orbital engine's methods on open-shell
// references (UHF, ROHF); the theory is summarised in open_shell_correlated.h.
//
// Orbital sets (as scf_open_hessian.cpp): UHF has two sets (alpha orbitals with the alpha
// occupations, beta orbitals with the beta ones), ROHF one set with both occupations (n^a = n^b = 1
// closed, n^a = 1 open); the non-redundant rotations are the pairs of orbitals of a set whose
// occupations differ. The orbital Hessian is applied with Fock builds (non-canonical orbitals:
// the full Fock matrices of both spins), the CPHF and Z-vector equations are solved together for
// all right-hand sides (solve_symmetric_batch).
//
// MO quantities of the two spins: one-particle matrices per spin, the two-electron integrals and
// densities in the blocks aa, bb (8-fold symmetric) and ab (alpha pair first, 4-fold symmetric).
// The perturbed Hamiltonians are stored as pair matrices [pq][rs] (ROHF: one block for all three,
// the alpha and beta orbitals being the same and rotating together).
#include "macqc/open_shell_correlated.h"

#include "macqc/ci.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <stdexcept>

#include "macqc/correlated_gradient.h"
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/scf_derivatives.h"

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

void add(Matrix& y, double a, const Matrix& x) {
    for (std::size_t i = 0; i < y.rows(); ++i)
        for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) += a * x(i, j);
}

// sum_ij a_ij b_ij
double dot(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    const std::size_t n = a.rows() * a.cols();
    for (std::size_t i = 0; i < n; ++i) s += a.data()[i] * b.data()[i];
    return s;
}

Matrix symmetrized(const Matrix& A) {
    Matrix S(A.rows(), A.cols());
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) S(i, j) = 0.5 * (A(i, j) + A(j, i));
    return S;
}

Matrix scale_columns(const Matrix& M, const std::vector<double>& n) {
    Matrix out = M;
    for (std::size_t i = 0; i < M.rows(); ++i)
        for (std::size_t j = 0; j < M.cols(); ++j) out(i, j) *= n[j];
    return out;
}

inline std::size_t pr(std::size_t p, std::size_t q) { return p >= q ? p * (p + 1) / 2 + q : q * (q + 1) / 2 + p; }

// ---------------------------------------------------------------- the orbital sets of the reference

struct OrbitalModel {
    struct Set {
        Matrix C;                                                  // N x n
        std::vector<double> occ[2];                                // alpha, beta occupation of every orbital
        std::vector<std::pair<std::size_t, std::size_t>> pairs;   // non-redundant (p > q)
        Matrix F[2];                                               // Fock matrices of both spins in these MOs
        Matrix eps;                                                // sum_s' F^s' N^s'
    };
    std::vector<Set> sets;
    int set_of_spin[2] = {0, 0};
    std::size_t N = 0, dim = 0;
    int nt = 1;
    FockBuilder fock;

    // invariant: one orbital set (scf.C) without non-redundant rotations, for a method whose energy
    // is invariant to every orbital rotation (Full-CI): no Z-vector, no CPHF
    bool invariant = false;

    OrbitalModel(const RHFResult& scf, const BasisSet& basis, int threads, double memory_mb, bool invariant_method = false)
        : N(scf.C.rows()), nt(threads), fock(basis, threads, memory_mb), invariant(invariant_method) {
        Matrix Fs[2] = {scf.Fa, scf.Fb};
        if (Fs[0].rows() == 0) {  // closed shell (RHF): F = h + G(P) for both spins
            Matrix F = scf.H;
            add(F, 1.0, fock.build(scf.P));
            Fs[0] = Fs[1] = F;
        }
        auto make_set = [&](const Matrix& C, int na, int nb) {
            Set s;
            s.C = C;
            const std::size_t n = C.cols();
            s.occ[0].assign(n, 0.0);
            s.occ[1].assign(n, 0.0);
            for (int i = 0; i < na; ++i) s.occ[0][static_cast<std::size_t>(i)] = 1.0;
            for (int i = 0; i < nb; ++i) s.occ[1][static_cast<std::size_t>(i)] = 1.0;
            if (!invariant)
                for (std::size_t p = 0; p < n; ++p)
                    for (std::size_t q = 0; q < p; ++q)
                        if (s.occ[0][p] != s.occ[0][q] || s.occ[1][p] != s.occ[1][q]) s.pairs.emplace_back(p, q);
            s.eps = Matrix(n, n);
            for (int sp = 0; sp < 2; ++sp) {
                s.F[sp] = mul(mul(C, true, Fs[sp], false, nt), false, C, false, nt);
                add(s.eps, 1.0, scale_columns(s.F[sp], s.occ[sp]));
            }
            sets.push_back(std::move(s));
        };
        if (scf.reference == "UHF" && !invariant) {
            make_set(scf.C, scf.n_alpha, 0);
            make_set(scf.Cb, 0, scf.n_beta);
            set_of_spin[0] = 0;
            set_of_spin[1] = 1;
        } else {
            const int na = scf.reference == "RHF" ? scf.n_occupied : scf.n_alpha;
            const int nb = scf.reference == "RHF" ? scf.n_occupied : scf.n_beta;
            make_set(scf.C, na, nb);
        }
        for (const Set& s : sets) dim += s.pairs.size();
    }

    // Spin densities of orbital responses U_s (one matrix per set): P_s' = sum_s C (U N + N U') C'.
    std::array<Matrix, 2> response_densities(const std::vector<Matrix>& U) const {
        std::array<Matrix, 2> D{Matrix(N, N), Matrix(N, N)};
        for (std::size_t s = 0; s < sets.size(); ++s)
            for (int sp = 0; sp < 2; ++sp) {
                const Matrix CU = mul(sets[s].C, false, scale_columns(U[s], sets[s].occ[sp]), false, nt);
                const Matrix A = mul(CU, false, sets[s].C, true, nt);
                Matrix& d = D[static_cast<std::size_t>(sp)];
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t j = 0; j < N; ++j) d(i, j) += A(i, j) + A(j, i);
            }
        return D;
    }
    // G_s'(D) = J(D_a + D_b) - K(D_s') of several density pairs
    std::vector<std::array<Matrix, 2>> fock_responses(const std::vector<std::array<Matrix, 2>>& D) const {
        std::vector<Matrix> in;
        for (const auto& d : D) {
            in.push_back(d[0]);
            in.push_back(d[1]);
        }
        std::vector<Matrix> J, K;
        if (!in.empty()) fock.build_jk(in, J, K);
        std::vector<std::array<Matrix, 2>> G(D.size());
        for (std::size_t k = 0; k < D.size(); ++k) {
            Matrix Jt = J[2 * k];
            add(Jt, 1.0, J[2 * k + 1]);
            for (int sp = 0; sp < 2; ++sp) {
                G[k][static_cast<std::size_t>(sp)] = Jt;
                add(G[k][static_cast<std::size_t>(sp)], -1.0, K[2 * k + static_cast<std::size_t>(sp)]);
            }
        }
        return G;
    }
    // first-order eps_s of a response U with the total Fock derivative Ft (AO, per spin)
    Matrix eps_response(std::size_t s, const Matrix& U, const std::array<Matrix, 2>& Ft) const {
        const Set& o = sets[s];
        Matrix e(o.C.cols(), o.C.cols());
        for (int sp = 0; sp < 2; ++sp) {
            Matrix M = mul(U, true, o.F[sp], false, nt);
            add(M, 1.0, mul(o.F[sp], false, U, false, nt));
            add(M, 1.0, mul(mul(o.C, true, Ft[static_cast<std::size_t>(sp)], false, nt), false, o.C, false, nt));
            add(e, 1.0, scale_columns(M, o.occ[sp]));
        }
        return e;
    }
    std::vector<Matrix> unpack(const std::vector<double>& k) const {  // antisymmetric matrix per set
        std::vector<Matrix> K;
        std::size_t pos = 0;
        for (const Set& o : sets) {
            Matrix m(o.C.cols(), o.C.cols());
            for (const auto& pq : o.pairs) {
                m(pq.first, pq.second) = k[pos];
                m(pq.second, pq.first) = -k[pos];
                ++pos;
            }
            K.push_back(m);
        }
        return K;
    }
    std::vector<double> antisym(const std::vector<Matrix>& e) const {  // e_pq - e_qp on the pairs
        std::vector<double> v;
        v.reserve(dim);
        for (std::size_t s = 0; s < sets.size(); ++s)
            for (const auto& pq : sets[s].pairs) v.push_back(e[s](pq.first, pq.second) - e[s](pq.second, pq.first));
        return v;
    }
    // the orbital Hessian: kappa -> antisym(eps^(1)(kappa)) for a batch
    std::vector<std::vector<double>> apply(const std::vector<std::vector<double>>& batch) const {
        std::vector<std::vector<Matrix>> Ks;
        std::vector<std::array<Matrix, 2>> Ds;
        for (const auto& b : batch) {
            Ks.push_back(unpack(b));
            Ds.push_back(response_densities(Ks.back()));
        }
        const std::vector<std::array<Matrix, 2>> G = fock_responses(Ds);
        std::vector<std::vector<double>> out;
        for (std::size_t k = 0; k < batch.size(); ++k) {
            std::vector<Matrix> e;
            for (std::size_t s = 0; s < sets.size(); ++s) e.push_back(eps_response(s, Ks[k][s], G[k]));
            out.push_back(antisym(e));
        }
        return out;
    }
    std::vector<double> diagonal() const {
        std::vector<double> d;
        for (const Set& o : sets)
            for (const auto& pq : o.pairs) {
                double v = 0.0;
                for (int sp = 0; sp < 2; ++sp)
                    v += (o.F[sp](pq.first, pq.first) - o.F[sp](pq.second, pq.second)) *
                         (o.occ[sp][pq.second] - o.occ[sp][pq.first]);
                d.push_back(v);
            }
        return d;
    }
    // occupation of spin sp in the orbitals of that spin
    const std::vector<double>& occupation(int sp) const {
        return sets[static_cast<std::size_t>(set_of_spin[sp])].occ[sp];
    }
    const Matrix& orbitals(int sp) const { return sets[static_cast<std::size_t>(set_of_spin[sp])].C; }
};

// ---------------------------------------------------------------- densities and Hamiltonians

// 8-fold symmetrisation of a full n^4 array into a packed one
Packed4 pack8(const std::vector<double>& G, std::size_t n) {
    Packed4 P(n);
    for (std::size_t p = 0; p < n; ++p)
        for (std::size_t q = 0; q < n; ++q)
            for (std::size_t r = 0; r < n; ++r)
                for (std::size_t s = 0; s < n; ++s) {
                    const double v = G[((p * n + q) * n + r) * n + s];
                    if (v != 0.0) P.add_raw(p, q, r, s, v);
                }
    return P;
}

// 4-fold symmetrisation (p <-> q, r <-> s) of an alpha-beta density, full layout
std::vector<double> sym4(const std::vector<double>& G, std::size_t na, std::size_t nb) {
    std::vector<double> out(G.size(), 0.0);
    auto at = [&](std::size_t a, std::size_t b, std::size_t c, std::size_t d) { return G[((a * na + b) * nb + c) * nb + d]; };
    for (std::size_t p = 0; p < na; ++p)
        for (std::size_t q = 0; q < na; ++q)
            for (std::size_t r = 0; r < nb; ++r)
                for (std::size_t s = 0; s < nb; ++s)
                    out[((p * na + q) * nb + r) * nb + s] =
                        0.25 * (at(p, q, r, s) + at(q, p, r, s) + at(p, q, s, r) + at(q, p, s, r));
    return out;
}

// The Lagrangian densities in the MOs of each spin: E = sum D h + sum Gaa (aa) + Gbb (bb) + Gab (ab).
struct Gamma {
    std::size_t n[2] = {0, 0};
    Matrix D[2];
    Packed4 Gaa, Gbb;
    std::vector<double> Gab;  // full, 4-fold symmetric (p, q alpha; r, s beta)
    double ab(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const {
        return Gab[((p * n[0] + q) * n[1] + r) * n[1] + s];
    }
    double& ab_ref(std::size_t p, std::size_t q, std::size_t r, std::size_t s) {
        return Gab[((p * n[0] + q) * n[1] + r) * n[1] + s];
    }
};

// A Hamiltonian in the MOs of both spins with the two-electron blocks as pair matrices
// [pq][rs] (aa, bb, ab with the alpha pair first).
struct PairHam {
    std::size_t n[2] = {0, 0}, P[2] = {0, 0};
    Matrix h[2];
    std::vector<double> store[3];
    const double* blk[3] = {nullptr, nullptr, nullptr};
    std::size_t ld[3] = {0, 0, 0};
    bool shared = false;  // bb and ab are the aa block (ROHF)

    PairHam() = default;
    PairHam(const PairHam&) = delete;
    PairHam& operator=(const PairHam&) = delete;
    PairHam(PairHam&&) = default;
    PairHam& operator=(PairHam&&) = default;

    double aa(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const { return blk[0][pr(p, q) * ld[0] + pr(r, s)]; }
    double bb(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const { return blk[1][pr(p, q) * ld[1] + pr(r, s)]; }
    double ab(std::size_t p, std::size_t q, std::size_t r, std::size_t s) const { return blk[2][pr(p, q) * ld[2] + pr(r, s)]; }
    double g(int sp, int sq, std::size_t p, std::size_t q, std::size_t r, std::size_t s) const {
        if (sp == 0 && sq == 0) return aa(p, q, r, s);
        if (sp == 1 && sq == 1) return bb(p, q, r, s);
        if (sp == 0) return ab(p, q, r, s);
        return ab(r, s, p, q);
    }
    void own(int k, std::vector<double> v, std::size_t row) {
        store[k] = std::move(v);
        blk[k] = store[k].data();
        ld[k] = row;
    }
};

PairHam reference_hamiltonian(const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab, bool uhf) {
    PairHam H;
    H.n[0] = static_cast<std::size_t>(aa.n);
    H.n[1] = static_cast<std::size_t>(uhf ? bb.n : aa.n);
    for (int s = 0; s < 2; ++s) H.P[s] = H.n[s] * (H.n[s] + 1) / 2;
    H.h[0] = aa.h;
    H.h[1] = uhf ? bb.h : aa.h;
    H.blk[0] = aa.eri.data();
    H.ld[0] = aa.stride;
    if (!uhf) {
        H.shared = true;
        H.blk[1] = H.blk[2] = aa.eri.data();
        H.ld[1] = H.ld[2] = aa.stride;
        return H;
    }
    H.blk[1] = bb.eri.data();
    H.ld[1] = bb.stride;
    // the mixed integrals are stored with the beta (ket) pair first: transpose
    std::vector<double> t(H.P[0] * H.P[1]);
    for (std::size_t pq = 0; pq < H.P[0]; ++pq)
        for (std::size_t rs = 0; rs < H.P[1]; ++rs) t[pq * H.P[1] + rs] = ab.mixed ? ab.eri[rs * ab.stride + pq] : ab.eri[pq * ab.stride + rs];
    H.own(2, std::move(t), H.P[1]);
    return H;
}

// Y^s: the generalized Fock matrices (orbital Lagrangians) of the densities G with the Hamiltonian K,
// Y_tp = d<G, K>/dU_tp for a transformation of the MOs of spin s, C -> C (1 + U).
std::array<Matrix, 2> lagrangian(const Gamma& G, const PairHam& K, int nt) {
    std::array<Matrix, 2> Y;
    for (int s = 0; s < 2; ++s) {
        Y[s] = mul(K.h[s], false, G.D[s], false, nt);
        for (std::size_t i = 0; i < Y[s].rows(); ++i)
            for (std::size_t j = 0; j < Y[s].cols(); ++j) Y[s](i, j) *= 2.0;
    }
    // same spin: 4 sum_qrs K(t,q,r,s) G(p,q,r,s)
    for (int s = 0; s < 2; ++s) {
        const std::size_t n = K.n[s], np = K.P[s];
        const Packed4& Gs = s ? G.Gbb : G.Gaa;
        std::vector<double> A(n * np), B(n * np), Xq(n * n);
        for (std::size_t q = 0; q < n; ++q) {
            check_cancelled();
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t r = 0, rs = 0; r < n; ++r)
                    for (std::size_t u = 0; u <= r; ++u, ++rs) {
                        A[t * np + rs] = s ? K.bb(t, q, r, u) : K.aa(t, q, r, u);
                        B[t * np + rs] = (r == u ? 1.0 : 2.0) * Gs.get(t, q, r, u);
                    }
            gemm(false, true, n, n, np, 4.0, A.data(), np, B.data(), np, 0.0, Xq.data(), n, nt);
            for (std::size_t t = 0; t < n; ++t)
                for (std::size_t p = 0; p < n; ++p) Y[s](t, p) += Xq[t * n + p];
        }
    }
    const std::size_t n0 = K.n[0], n1 = K.n[1], P0 = K.P[0], P1 = K.P[1];
    // alpha side of ab: 2 sum_qrs K_ab(t,q,r,s) G_ab(p,q,r,s)
    {
        std::vector<double> A(n0 * P1), B(n0 * P1), Xq(n0 * n0);
        for (std::size_t q = 0; q < n0; ++q) {
            check_cancelled();
            for (std::size_t t = 0; t < n0; ++t)
                for (std::size_t r = 0, rs = 0; r < n1; ++r)
                    for (std::size_t u = 0; u <= r; ++u, ++rs) {
                        A[t * P1 + rs] = K.ab(t, q, r, u);
                        B[t * P1 + rs] = (r == u ? 1.0 : 2.0) * G.ab(t, q, r, u);
                    }
            gemm(false, true, n0, n0, P1, 2.0, A.data(), P1, B.data(), P1, 0.0, Xq.data(), n0, nt);
            for (std::size_t t = 0; t < n0; ++t)
                for (std::size_t p = 0; p < n0; ++p) Y[0](t, p) += Xq[t * n0 + p];
        }
    }
    // beta side of ab: 2 sum_pqs K_ab(p,q,t,s) G_ab(p,q,r,s)
    {
        std::vector<double> A(n1 * P0), B(n1 * P0), Xs(n1 * n1);
        for (std::size_t u = 0; u < n1; ++u) {
            check_cancelled();
            for (std::size_t t = 0; t < n1; ++t)
                for (std::size_t p = 0, pq = 0; p < n0; ++p)
                    for (std::size_t q = 0; q <= p; ++q, ++pq) {
                        A[t * P0 + pq] = K.ab(p, q, t, u);
                        B[t * P0 + pq] = (p == q ? 1.0 : 2.0) * G.ab(p, q, t, u);
                    }
            gemm(false, true, n1, n1, P0, 2.0, A.data(), P0, B.data(), P0, 0.0, Xs.data(), n1, nt);
            for (std::size_t t = 0; t < n1; ++t)
                for (std::size_t r = 0; r < n1; ++r) Y[1](t, r) += Xs[t * n1 + r];
        }
    }
    return Y;
}

// <G, K> over all index values (Gamma, symmetrised) and for the full arrays of so_densities.
double contract(const Gamma& G, const PairHam& K) {
    double e = 0.0;
    for (int s = 0; s < 2; ++s) e += dot(G.D[s], K.h[s]);
    for (int s = 0; s < 2; ++s) {
        const std::size_t n = K.n[s];
        const Packed4& Gs = s ? G.Gbb : G.Gaa;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q <= p; ++q)
                for (std::size_t r = 0; r <= p; ++r)
                    for (std::size_t u = 0; u <= (r == p ? q : r); ++u)
                        e += Packed4::degeneracy(p, q, r, u) * Gs.get(p, q, r, u) * (s ? K.bb(p, q, r, u) : K.aa(p, q, r, u));
    }
    const std::size_t n0 = K.n[0], n1 = K.n[1];
    for (std::size_t p = 0; p < n0; ++p)
        for (std::size_t q = 0; q < n0; ++q)
            for (std::size_t r = 0; r < n1; ++r)
                for (std::size_t u = 0; u < n1; ++u) e += G.ab(p, q, r, u) * K.ab(p, q, r, u);
    return e;
}

double contract(const SODensities& d, const PairHam& K) {
    double e = 0.0;
    for (int s = 0; s < 2; ++s) e += dot(d.D[s], K.h[s]);
    const std::size_t n0 = K.n[0], n1 = K.n[1];
    for (std::size_t p = 0; p < n0; ++p)
        for (std::size_t q = 0; q < n0; ++q)
            for (std::size_t r = 0; r < n0; ++r) {
                const double* g = d.Gaa.data() + ((p * n0 + q) * n0 + r) * n0;
                for (std::size_t u = 0; u < n0; ++u)
                    if (g[u] != 0.0) e += g[u] * K.aa(p, q, r, u);
            }
    for (std::size_t p = 0; p < n1; ++p)
        for (std::size_t q = 0; q < n1; ++q)
            for (std::size_t r = 0; r < n1; ++r) {
                const double* g = d.Gbb.data() + ((p * n1 + q) * n1 + r) * n1;
                for (std::size_t u = 0; u < n1; ++u)
                    if (g[u] != 0.0) e += g[u] * K.bb(p, q, r, u);
            }
    for (std::size_t p = 0; p < n0; ++p)
        for (std::size_t q = 0; q < n0; ++q)
            for (std::size_t r = 0; r < n1; ++r) {
                const double* g = d.Gab.data() + ((p * n0 + q) * n1 + r) * n1;
                for (std::size_t u = 0; u < n1; ++u)
                    if (g[u] != 0.0) e += g[u] * K.ab(p, q, r, u);
            }
    return e;
}

// packed 8-fold layout (pq >= rs) -> full pair matrix [pq][rs]
std::vector<double> pair_matrix(const std::vector<double>& packed, std::size_t P) {
    std::vector<double> out(P * P);
    for (std::size_t pq = 0; pq < P; ++pq)
        for (std::size_t rs = 0; rs <= pq; ++rs) out[pq * P + rs] = out[rs * P + pq] = packed[pq * (pq + 1) / 2 + rs];
    return out;
}

// (pq|rs) with p, q orbitals of Ca and r, s of Cb from packed AO integrals: [pq][rs]
std::vector<double> transform_packed_mixed(const std::vector<double>& in, std::size_t N, const Matrix& Ca,
                                           const Matrix& Cb, int nt) {
    const std::size_t na = Ca.cols(), nb = Cb.cols(), NP = N * (N + 1) / 2;
    const std::size_t Pa = na * (na + 1) / 2, Pb = nb * (nb + 1) / 2;
    std::vector<double> half(Pa * NP, 0.0);
    Matrix S(N, N);
    for (std::size_t l = 0, ls = 0; l < N; ++l)
        for (std::size_t s = 0; s <= l; ++s, ++ls) {
            check_cancelled();
            for (std::size_t m = 0; m < N; ++m)
                for (std::size_t k = 0; k <= m; ++k) S(m, k) = S(k, m) = in[pr(pr(m, k), ls)];
            const Matrix A = mul(Ca, true, mul(S, false, Ca, false, nt), false, nt);
            for (std::size_t p = 0, pq = 0; p < na; ++p)
                for (std::size_t q = 0; q <= p; ++q, ++pq) half[pq * NP + ls] = A(p, q);
        }
    std::vector<double> out(Pa * Pb, 0.0);
    for (std::size_t pq = 0; pq < Pa; ++pq) {
        check_cancelled();
        const double* h = half.data() + pq * NP;
        for (std::size_t l = 0, ls = 0; l < N; ++l)
            for (std::size_t s = 0; s <= l; ++s, ++ls) S(l, s) = S(s, l) = h[ls];
        const Matrix A = mul(Cb, true, mul(S, false, Cb, false, nt), false, nt);
        for (std::size_t r = 0, rs = 0; r < nb; ++r)
            for (std::size_t s = 0; s <= r; ++s, ++rs) out[pq * Pb + rs] = A(r, s);
    }
    return out;
}

// R[pq][c] = sum_t (U_tp g(t,q,c) + U_tq g(p,t,c)) over the pairs pq of n orbitals and the columns c
template <class G>
std::vector<double> first_pair_rotation(std::size_t n, std::size_t ncol, const Matrix& U, G g, int nt) {
    const std::size_t P = n * (n + 1) / 2;
    std::vector<double> R(P * ncol, 0.0);
    if (ncol == 0 || n == 0) return R;
    const std::size_t chunk = std::max<std::size_t>(1, std::min<std::size_t>(ncol, (std::size_t{8} << 20) / (n * n)));
    std::vector<double> GG(n * chunk * n), BB(n * chunk * n);
    for (std::size_t c0 = 0; c0 < ncol; c0 += chunk) {
        check_cancelled();
        const std::size_t nbc = std::min(chunk, ncol - c0), w = nbc * n;
        for (std::size_t t = 0; t < n; ++t)
            for (std::size_t b = 0; b < nbc; ++b)
                for (std::size_t q = 0; q < n; ++q) GG[t * w + b * n + q] = g(t, q, c0 + b);
        gemm(true, false, n, w, n, 1.0, U.data(), n, GG.data(), w, 0.0, BB.data(), w, nt);
        for (std::size_t b = 0; b < nbc; ++b)
            for (std::size_t p = 0; p < n; ++p)
                for (std::size_t q = 0; q <= p; ++q) R[pr(p, q) * ncol + c0 + b] = BB[p * w + b * n + q] + BB[q * w + b * n + p];
    }
    return R;
}

// out += the MO Hamiltonian H transformed to first order by the orbital responses U[spin]:
// h^(U) = U' h + h U, (pq|rs)^(U) with U on every index. `out` must own its blocks.
void add_rotated(const PairHam& H, const std::array<Matrix, 2>& U, PairHam& out, int nt) {
    for (int s = 0; s < (H.shared ? 1 : 2); ++s) {
        add(out.h[s], 1.0, mul(U[s], true, H.h[s], false, nt));
        add(out.h[s], 1.0, mul(H.h[s], false, U[s], false, nt));
    }
    if (H.shared) out.h[1] = out.h[0];
    std::vector<std::size_t> rr[2], ss[2];
    for (int s = 0; s < 2; ++s)
        for (std::size_t r = 0; r < H.n[s]; ++r)
            for (std::size_t u = 0; u <= r; ++u) {
                rr[s].push_back(r);
                ss[s].push_back(u);
            }
    for (int s = 0; s < (H.shared ? 1 : 2); ++s) {
        const std::size_t n = H.n[s], P = H.P[s];
        const std::vector<double> R = first_pair_rotation(
            n, P, U[s], [&](std::size_t t, std::size_t q, std::size_t c) { return s ? H.bb(t, q, rr[s][c], ss[s][c]) : H.aa(t, q, rr[s][c], ss[s][c]); },
            nt);
        std::vector<double>& o = out.store[s];
        for (std::size_t pq = 0; pq < P; ++pq)
            for (std::size_t rs = 0; rs < P; ++rs) o[pq * P + rs] += R[pq * P + rs] + R[rs * P + pq];
    }
    if (H.shared) return;
    const std::size_t n0 = H.n[0], n1 = H.n[1], P0 = H.P[0], P1 = H.P[1];
    const std::vector<double> Ra = first_pair_rotation(
        n0, P1, U[0], [&](std::size_t t, std::size_t q, std::size_t c) { return H.ab(t, q, rr[1][c], ss[1][c]); }, nt);
    const std::vector<double> Rb = first_pair_rotation(
        n1, P0, U[1], [&](std::size_t t, std::size_t q, std::size_t c) { return H.ab(rr[0][c], ss[0][c], t, q); }, nt);
    std::vector<double>& o = out.store[2];
    for (std::size_t pq = 0; pq < P0; ++pq)
        for (std::size_t rs = 0; rs < P1; ++rs) o[pq * P1 + rs] += Ra[pq * P1 + rs] + Rb[rs * P0 + pq];
}

// AO two-particle density of the alpha-beta part (bra Ca, ket Cb), symmetrised over the exchange of
// the pairs, added to the packed layout of FockBuilder::packed() (as uhf_correlated.cpp).
void back_transform_mixed(const std::vector<double>& Gab, const Matrix& Ca, const Matrix& Cb, std::vector<double>& out,
                          int nt) {
    const std::size_t N = Ca.rows(), na = Ca.cols(), nb = Cb.cols();
    const std::size_t NP = N * (N + 1) / 2, npb = nb * (nb + 1) / 2;
    std::vector<double> half(NP * npb, 0.0);
    Matrix M(na, na);
    for (std::size_t r = 0, rs = 0; r < nb; ++r)
        for (std::size_t s = 0; s <= r; ++s, ++rs) {
            check_cancelled();
            const double w = r == s ? 1.0 : 2.0;
            for (std::size_t p = 0; p < na; ++p)
                for (std::size_t q = 0; q < na; ++q) M(p, q) = w * Gab[((p * na + q) * nb + r) * nb + s];
            const Matrix A = mul(Ca, false, mul(M, false, Ca, true, nt), false, nt);
            for (std::size_t m = 0, mn = 0; m < N; ++m)
                for (std::size_t k = 0; k <= m; ++k, ++mn) half[mn * npb + rs] = A(m, k);
        }
    std::vector<double> full(NP * NP, 0.0);
    Matrix M2(nb, nb);
    for (std::size_t mn = 0; mn < NP; ++mn) {
        check_cancelled();
        const double* h = half.data() + mn * npb;
        for (std::size_t r = 0, rs = 0; r < nb; ++r)
            for (std::size_t s = 0; s <= r; ++s, ++rs) {
                const double x = (r == s ? 1.0 : 0.5) * h[rs];
                M2(r, s) = M2(s, r) = x;
            }
        const Matrix B = mul(Cb, false, mul(M2, false, Cb, true, nt), false, nt);
        for (std::size_t l = 0, ls = 0; l < N; ++l)
            for (std::size_t s = 0; s <= l; ++s, ++ls) full[mn * NP + ls] = B(l, s);
    }
    for (std::size_t mn = 0; mn < NP; ++mn)
        for (std::size_t ls = 0; ls <= mn; ++ls) out[mn * (mn + 1) / 2 + ls] += 0.5 * (full[mn * NP + ls] + full[ls * NP + mn]);
}

// The densities of the SCF energy itself (the orbital-part check): E = sum_s tr(N_s h_s) +
// 1/2 sum_ss' sum_ij (ii|jj) - 1/2 sum_s sum_ij (ij|ji) over the occupied orbitals.
SODensities scf_densities(const OrbitalModel& m) {
    SODensities d;
    for (int s = 0; s < 2; ++s) d.n[s] = static_cast<int>(m.orbitals(s).cols());
    const std::size_t n0 = static_cast<std::size_t>(d.n[0]), n1 = static_cast<std::size_t>(d.n[1]);
    const std::vector<double>& na = m.occupation(0);
    const std::vector<double>& nb = m.occupation(1);
    for (int s = 0; s < 2; ++s) {
        const std::size_t n = s ? n1 : n0;
        const std::vector<double>& o = s ? nb : na;
        d.D[s] = Matrix(n, n);
        for (std::size_t i = 0; i < n; ++i) d.D[s](i, i) = o[i];
    }
    d.Gaa.assign(n0 * n0 * n0 * n0, 0.0);
    d.Gbb.assign(n1 * n1 * n1 * n1, 0.0);
    d.Gab.assign(n0 * n0 * n1 * n1, 0.0);
    for (std::size_t i = 0; i < n0; ++i)
        for (std::size_t j = 0; j < n0; ++j) {
            const double w = na[i] * na[j];
            if (w == 0.0) continue;
            d.Gaa[((i * n0 + i) * n0 + j) * n0 + j] += 0.5 * w;
            d.Gaa[((i * n0 + j) * n0 + j) * n0 + i] -= 0.5 * w;
        }
    for (std::size_t i = 0; i < n1; ++i)
        for (std::size_t j = 0; j < n1; ++j) {
            const double w = nb[i] * nb[j];
            if (w == 0.0) continue;
            d.Gbb[((i * n1 + i) * n1 + j) * n1 + j] += 0.5 * w;
            d.Gbb[((i * n1 + j) * n1 + j) * n1 + i] -= 0.5 * w;
        }
    for (std::size_t i = 0; i < n0; ++i)
        for (std::size_t j = 0; j < n1; ++j)
            if (na[i] * nb[j] != 0.0) d.Gab[((i * n0 + i) * n1 + j) * n1 + j] += 1.0;
    d.lambda_converged = true;
    return d;
}

// ---------------------------------------------------------------- the Lagrangian of the reference

struct Prepared {
    std::unique_ptr<OrbitalModel> model;
    PairHam H0;
    Gamma gam;                  // Gamma_tot = Gamma + Gamma_z
    std::array<Matrix, 2> Y;    // generalized Fock of Gamma_tot with H0, per spin
    std::vector<Matrix> Yset;   // per orbital set: the symmetric part of the sum over its spins
    Matrix Dao, Wao;            // relaxed one-particle density, energy-weighted density (AO)
    std::vector<double> gao;    // AO two-particle density (packed)
    OpenShellGradient info;
};

Prepared prepare(const Molecule& mol, const BasisSet& basis, const RHFResult& scf, const MOIntegrals& aa,
                 const MOIntegrals& bb, const MOIntegrals& ab, const SODensities& d, int nt, double memory_mb,
                 const CPHFCallback& progress, bool invariant = false) {
    if (scf.removed_functions > 0) throw std::runtime_error("analytic correlated derivatives with removed near-linear dependencies");
    Prepared P;
    P.model = std::make_unique<OrbitalModel>(scf, basis, nt, memory_mb, invariant);
    OrbitalModel& M = *P.model;
    const bool uhf = M.sets.size() == 2;
    P.H0 = reference_hamiltonian(aa, bb, ab, uhf);
    Gamma& G = P.gam;
    for (int s = 0; s < 2; ++s) G.n[s] = P.H0.n[s];
    for (int s = 0; s < 2; ++s) G.D[s] = symmetrized(d.D[s]);
    G.Gaa = pack8(d.Gaa, G.n[0]);
    G.Gbb = pack8(d.Gbb, G.n[1]);
    G.Gab = sym4(d.Gab, G.n[0], G.n[1]);

    // ---- Z-vector: A z = -(Y - Y')_pairs with the method's own Lagrangian
    report_progress(JobStage::Gradient, 1);
    const std::array<Matrix, 2> Yc = lagrangian(G, P.H0, nt);
    std::vector<Matrix> Ys(M.sets.size());
    for (std::size_t s = 0; s < M.sets.size(); ++s) Ys[s] = Matrix(M.sets[s].C.cols(), M.sets[s].C.cols());
    for (int sp = 0; sp < 2; ++sp) add(Ys[static_cast<std::size_t>(M.set_of_spin[sp])], 1.0, Yc[sp]);
    for (std::size_t s = 0; s < M.sets.size(); ++s) {
        const OrbitalModel::Set& o = M.sets[s];
        for (std::size_t p = 0; p < o.C.cols(); ++p)
            for (std::size_t q = 0; q < p; ++q)
                if (M.invariant || (o.occ[0][p] == o.occ[0][q] && o.occ[1][p] == o.occ[1][q]))
                    P.info.redundant_asymmetry = std::max(P.info.redundant_asymmetry, std::fabs(Ys[s](p, q) - Ys[s](q, p)));
    }
    std::vector<double> zsol(M.dim, 0.0);
    P.info.zvector_converged = true;
    if (M.dim > 0) {
        std::vector<double> rhs = M.antisym(Ys);
        for (double& v : rhs) v = -v;
        const LinearBatchSolution sol =
            solve_symmetric_batch([&](const std::vector<std::vector<double>>& b) { return M.apply(b); }, M.diagonal(), {rhs},
                                  1.0e-10, 100, progress);
        P.info.zvector_converged = sol.converged;
        P.info.zvector_iterations = sol.iterations;
        P.info.zvector_residual = sol.residual;
        zsol = sol.x[0];
    }
    const std::vector<Matrix> Zhat = M.unpack(zsol);

    // ---- Gamma_z: z . g = sum_s' <Rbar_s', F_s'>, Rbar_pq = Zhat_pq (n_q - n_p) / 2
    Matrix R[2];
    for (int sp = 0; sp < 2; ++sp) {
        const Matrix& Z = Zhat[static_cast<std::size_t>(M.set_of_spin[sp])];
        const std::vector<double>& n = M.occupation(sp);
        R[sp] = Matrix(Z.rows(), Z.cols());
        for (std::size_t p = 0; p < Z.rows(); ++p)
            for (std::size_t q = 0; q < Z.cols(); ++q) R[sp](p, q) = 0.5 * Z(p, q) * (n[q] - n[p]);
        add(G.D[sp], 1.0, R[sp]);
    }
    {
        const std::vector<double>& na = M.occupation(0);
        const std::vector<double>& nb = M.occupation(1);
        const std::size_t n0 = G.n[0], n1 = G.n[1];
        for (std::size_t p = 0; p < n0; ++p)
            for (std::size_t q = 0; q < n0; ++q) {
                const double r = R[0](p, q);
                if (r == 0.0) continue;
                for (std::size_t k = 0; k < n0; ++k)
                    if (na[k] != 0.0) {
                        G.Gaa.add_raw(p, q, k, k, r * na[k]);   // Coulomb, alpha with alpha
                        G.Gaa.add_raw(p, k, k, q, -r * na[k]);  // exchange
                    }
                for (std::size_t k = 0; k < n1; ++k)
                    if (nb[k] != 0.0) G.ab_ref(p, q, k, k) += r * nb[k];  // Coulomb, alpha with beta
            }
        for (std::size_t p = 0; p < n1; ++p)
            for (std::size_t q = 0; q < n1; ++q) {
                const double r = R[1](p, q);
                if (r == 0.0) continue;
                for (std::size_t k = 0; k < n1; ++k)
                    if (nb[k] != 0.0) {
                        G.Gbb.add_raw(p, q, k, k, r * nb[k]);
                        G.Gbb.add_raw(p, k, k, q, -r * nb[k]);
                    }
                for (std::size_t k = 0; k < n0; ++k)
                    if (na[k] != 0.0) G.ab_ref(k, k, p, q) += r * na[k];  // (pq|kk)_ba = (kk|pq)_ab
            }
    }

    // ---- the Lagrangian of Gamma_tot: stationary (symmetric per set), and the AO densities
    P.Y = lagrangian(G, P.H0, nt);
    P.Yset.resize(M.sets.size());
    for (std::size_t s = 0; s < M.sets.size(); ++s) P.Yset[s] = Matrix(M.sets[s].C.cols(), M.sets[s].C.cols());
    for (int sp = 0; sp < 2; ++sp) add(P.Yset[static_cast<std::size_t>(M.set_of_spin[sp])], 1.0, P.Y[sp]);
    for (std::size_t s = 0; s < M.sets.size(); ++s) {
        for (const auto& pq : M.sets[s].pairs)
            P.info.stationarity = std::max(P.info.stationarity, std::fabs(P.Yset[s](pq.first, pq.second) - P.Yset[s](pq.second, pq.first)));
        P.Yset[s] = symmetrized(P.Yset[s]);
    }
    const std::size_t N = scf.C.rows();
    P.Dao = Matrix(N, N);
    P.Wao = Matrix(N, N);
    for (int sp = 0; sp < 2; ++sp) {
        const Matrix& C = M.orbitals(sp);
        add(P.Dao, 1.0, mul(mul(C, false, G.D[sp], false, nt), false, C, true, nt));
        Matrix y = symmetrized(P.Y[sp]);
        for (std::size_t i = 0; i < y.rows(); ++i)
            for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) *= 0.5;
        add(P.Wao, 1.0, mul(mul(C, false, y, false, nt), false, C, true, nt));
    }
    P.Dao = symmetrized(P.Dao);
    P.Wao = symmetrized(P.Wao);
    report_progress(JobStage::Gradient, 2);
    P.gao = back_transform_packed(G.Gaa, M.orbitals(0), nt);
    {
        const std::vector<double> gb = back_transform_packed(G.Gbb, M.orbitals(1), nt);
        for (std::size_t k = 0; k < P.gao.size(); ++k) P.gao[k] += gb[k];
    }
    back_transform_mixed(G.Gab, M.orbitals(0), M.orbitals(1), P.gao, nt);
    P.info.density_energy = scf.nuclear_repulsion + contract(G, P.H0);
    (void)mol;
    return P;
}

OpenShellGradient gradient_of(const Molecule& mol, const BasisSet& basis, const Prepared& P, int nt) {
    OpenShellGradient res = P.info;
    const int nat = static_cast<int>(mol.atoms.size());
    report_progress(JobStage::Gradient, 3);
    res.nuclear = nuclear_repulsion_gradient(mol);
    res.one_electron = core_hamiltonian_gradient(basis, mol, P.Dao, nt);
    TwoParticleDensity gamma;
    gamma.packed = &P.gao;
    res.two_electron = eri_gradient(basis, nat, gamma, nt);
    res.overlap = overlap_gradient(basis, nat, P.Wao, nt);
    res.gradient.assign(static_cast<std::size_t>(3 * nat), 0.0);
    for (std::size_t x = 0; x < res.gradient.size(); ++x) {
        res.overlap[x] = -res.overlap[x];
        res.gradient[x] = res.nuclear[x] + res.one_electron[x] + res.two_electron[x] + res.overlap[x];
    }
    res.relaxed_density_ao = P.Dao;
    res.relaxed_dipole = dipole_moment(mol, basis, P.Dao);
    return res;
}

// ---------------------------------------------------------------- the Hessian

// The method's second-order terms from the first-order Hamiltonians H^[k] (amplitude responses):
// A(p, k) = <dGamma[t^k], H^[p]> and the pair terms, with their convergence in `res`.
using SecondOrder = std::function<void(const std::vector<PairHam>& Hk, Matrix& A, Matrix& pairs, CorrelatedHessianResult& res)>;

CorrelatedHessianResult hessian_core(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                     const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                     const SODensities& d0, const SecondOrder& second, int nt, double memory_mb,
                                     bool dipoles, bool invariant = false) {
    const auto t0 = Clock::now();
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t n3 = static_cast<std::size_t>(3 * nat);
    CorrelatedHessianResult res;
    Prepared P = prepare(mol, basis, scf, aa, bb, ab, d0, nt, memory_mb, CPHFCallback(), invariant);
    OrbitalModel& M = *P.model;
    const bool uhf = M.sets.size() == 2;
    res.zvector_converged = P.info.zvector_converged;
    res.zvector_iterations = P.info.zvector_iterations;
    res.zvector_residual = P.info.zvector_residual;
    const OpenShellGradient grad = gradient_of(mol, basis, P, nt);
    res.gradient = grad.gradient;
    res.dipole = grad.relaxed_dipole;
    const std::size_t N = scf.C.rows();
    const std::size_t nsets = M.sets.size();

    // ---- perturbations: the internal directions, then the field components
    const std::vector<std::vector<double>> Q = vibrational_directions(mol);
    const std::size_t m = Q.size(), nf = dipoles ? 3 : 0, np = m + nf;
    res.perturbations = m;

    // ---- first-order AO quantities and the orbital responses (CPHF of the reference)
    report_progress(JobStage::Hessian, 1);
    const auto ti = Clock::now();
    const std::vector<Matrix> Sx = overlap_derivative_matrices(basis, nat, nt);
    const std::vector<Matrix> hx = core_hamiltonian_derivative_matrices(basis, mol, nt);
    std::vector<std::vector<Matrix>> Gx;  // the skeleton Fock derivatives: only for the CPHF
    if (M.dim > 0) Gx = eri_derivative_fock_multi(basis, nat, {scf.P, scf.Pa, scf.Pb}, {1.0, 0.0, 0.0}, {0.0, -1.0, -1.0}, nt);
    std::vector<Matrix> rf;
    if (nf) rf = dipole_matrices(basis, {0.0, 0.0, 0.0});
    auto combine = [&](const std::vector<Matrix>& v, const std::vector<double>& dir) {
        Matrix out(N, N);
        for (std::size_t x = 0; x < n3; ++x)
            if (dir[x] != 0.0) add(out, dir[x], v[x]);
        return out;
    };
    std::vector<Matrix> Sao(np, Matrix(N, N)), hao(np);
    std::vector<std::array<Matrix, 2>> Fskel(np);
    for (std::size_t k = 0; k < np; ++k) {
        if (k < m) {
            Sao[k] = combine(Sx, Q[k]);
            hao[k] = combine(hx, Q[k]);
            if (M.dim > 0) {
                const Matrix J = combine(Gx[0], Q[k]);
                for (int sp = 0; sp < 2; ++sp) {
                    Fskel[k][static_cast<std::size_t>(sp)] = hao[k];
                    add(Fskel[k][static_cast<std::size_t>(sp)], 1.0, J);
                    add(Fskel[k][static_cast<std::size_t>(sp)], 1.0, combine(Gx[1 + static_cast<std::size_t>(sp)], Q[k]));
                }
            }
        } else {
            hao[k] = rf[k - m];  // H(F) = H + F . r
            Fskel[k] = {hao[k], hao[k]};
        }
    }
    // S~ per set, U0 = -S~/2, right-hand sides
    std::vector<std::vector<Matrix>> Smo(np), U(np);
    std::vector<std::array<Matrix, 2>> D0(np);
    for (std::size_t k = 0; k < np; ++k) {
        for (std::size_t s = 0; s < nsets; ++s) {
            Smo[k].push_back(mul(mul(M.sets[s].C, true, Sao[k], false, nt), false, M.sets[s].C, false, nt));
            Matrix u = Smo[k][s];
            for (std::size_t i = 0; i < u.rows(); ++i)
                for (std::size_t j = 0; j < u.cols(); ++j) u(i, j) *= -0.5;
            U[k].push_back(u);
        }
        D0[k] = M.response_densities(U[k]);
    }
    res.cphf_converged = true;
    if (M.dim > 0) {
        const std::vector<std::array<Matrix, 2>> G0 = M.fock_responses(D0);
        std::vector<std::vector<double>> rhs(np);
        for (std::size_t k = 0; k < np; ++k) {
            std::array<Matrix, 2> Ft = G0[k];
            for (int sp = 0; sp < 2; ++sp) add(Ft[static_cast<std::size_t>(sp)], 1.0, Fskel[k][static_cast<std::size_t>(sp)]);
            std::vector<Matrix> e;
            for (std::size_t s = 0; s < nsets; ++s) e.push_back(M.eps_response(s, U[k][s], Ft));
            rhs[k] = M.antisym(e);
            for (double& v : rhs[k]) v = -v;
        }
        report_progress(JobStage::CPHF, 0);
        const LinearBatchSolution sol = solve_symmetric_batch(
            [&](const std::vector<std::vector<double>>& b) { return M.apply(b); }, M.diagonal(), rhs, 1.0e-10, 100);
        res.cphf_converged = sol.converged;
        res.cphf_iterations = sol.iterations;
        res.cphf_residual = sol.residual;
        for (std::size_t k = 0; k < np; ++k) {
            const std::vector<Matrix> K = M.unpack(sol.x[k]);
            for (std::size_t s = 0; s < nsets; ++s) add(U[k][s], 1.0, K[s]);
        }
    }
    res.orbital_seconds = seconds_since(ti);

    // ---- the first-order MO Hamiltonians H~^k (derivative integrals) and H^[k] = H~^k + H(U^k), and
    // the generalized Fock matrices of Gamma_tot with both
    const std::size_t n0 = P.H0.n[0], n1 = P.H0.n[1], P0 = P.H0.P[0], P1 = P.H0.P[1];
    std::vector<std::vector<double>> gx;
    if (m) gx = eri_derivative_packed(basis, nat, nt);
    std::vector<PairHam> Hk(np);
    std::vector<std::vector<Matrix>> Yk(np), Yder(np);  // per set
    auto per_set = [&](const std::array<Matrix, 2>& Y) {
        std::vector<Matrix> out(nsets);
        for (std::size_t s = 0; s < nsets; ++s) out[s] = Matrix(M.sets[s].C.cols(), M.sets[s].C.cols());
        for (int sp = 0; sp < 2; ++sp) add(out[static_cast<std::size_t>(M.set_of_spin[sp])], 1.0, Y[sp]);
        return out;
    };
    for (std::size_t k = 0; k < np; ++k) {
        report_progress(JobStage::Hessian, 2, 0.0, static_cast<double>(k + 1));
        PairHam Hd;
        for (int s = 0; s < 2; ++s) {
            Hd.n[s] = P.H0.n[s];
            Hd.P[s] = P.H0.P[s];
            Hd.h[s] = mul(mul(M.orbitals(s), true, hao[k], false, nt), false, M.orbitals(s), false, nt);
        }
        if (k < m) {
            std::vector<double> gq(gx.empty() ? 0 : gx[0].size(), 0.0);
            for (std::size_t x = 0; x < n3; ++x) {
                const double w = Q[k][x];
                if (w == 0.0) continue;
                const std::vector<double>& g = gx[x];
                for (std::size_t i = 0; i < gq.size(); ++i) gq[i] += w * g[i];
            }
            Hd.own(0, pair_matrix(transform_packed(gq, N, M.orbitals(0), nt), P0), P0);
            if (uhf) {
                Hd.own(1, pair_matrix(transform_packed(gq, N, M.orbitals(1), nt), P1), P1);
                Hd.own(2, transform_packed_mixed(gq, N, M.orbitals(0), M.orbitals(1), nt), P1);
            }
        } else {
            Hd.own(0, std::vector<double>(P0 * P0, 0.0), P0);
            if (uhf) {
                Hd.own(1, std::vector<double>(P1 * P1, 0.0), P1);
                Hd.own(2, std::vector<double>(P0 * P1, 0.0), P1);
            }
        }
        if (!uhf) {
            Hd.shared = true;
            Hd.blk[1] = Hd.blk[2] = Hd.blk[0];
            Hd.ld[1] = Hd.ld[2] = Hd.ld[0];
        }
        Yder[k] = per_set(lagrangian(P.gam, Hd, nt));
        // H^[k] = H~^k + H(U^k)
        PairHam& H = Hk[k];
        for (int s = 0; s < 2; ++s) {
            H.n[s] = Hd.n[s];
            H.P[s] = Hd.P[s];
            H.h[s] = Hd.h[s];
        }
        H.own(0, std::move(Hd.store[0]), P0);
        if (uhf) {
            H.own(1, std::move(Hd.store[1]), P1);
            H.own(2, std::move(Hd.store[2]), P1);
        }
        std::array<Matrix, 2> Us;
        for (int sp = 0; sp < 2; ++sp) Us[static_cast<std::size_t>(sp)] = U[k][static_cast<std::size_t>(M.set_of_spin[sp])];
        add_rotated(P.H0, Us, H, nt);
        if (!uhf) {
            H.shared = true;
            H.blk[1] = H.blk[2] = H.blk[0];
            H.ld[1] = H.ld[2] = H.ld[0];
            H.h[1] = H.h[0];
        }
        Yk[k] = per_set(lagrangian(P.gam, H, nt));
    }
    (void)n0;
    (void)n1;

    // ---- the method's responses: amplitude responses, density responses and pair terms
    const auto tm = Clock::now();
    Matrix A(np, np), pairs(np, np);
    if (second) second(Hk, A, pairs, res);
    res.method_seconds = seconds_since(tm);

    // ---- explicit second derivatives with Gamma_tot (Cartesian)
    const auto te = Clock::now();
    report_progress(JobStage::Hessian, 3);
    TwoParticleDensity gamma;
    gamma.packed = &P.gao;
    Matrix Hx = core_hamiltonian_hessian(basis, mol, P.Dao, nt);
    add(Hx, 1.0, eri_hessian(basis, nat, gamma, nt));
    add(Hx, -1.0, overlap_hessian(basis, nat, P.Wao, nt));
    add(Hx, 1.0, nuclear_repulsion_hessian(mol));
    std::vector<std::array<double, 3>> dip_expl(n3, {0.0, 0.0, 0.0});
    if (nf) {
        const std::vector<std::vector<Matrix>> rx = dipole_derivative_matrices(basis, nat, {0.0, 0.0, 0.0}, nt);
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t f = 0; f < 3; ++f) {
                dip_expl[x][f] = dot(P.Dao, rx[x][f]);
                if (x % 3 == f) dip_expl[x][f] -= mol.atoms[x / 3].nuclear_charge();
            }
    }
    res.explicit_seconds = seconds_since(te);

    // ---- the bilinear orbital terms and the whole block over the perturbations
    const auto tr = Clock::now();
    Matrix B(np, np);
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t q = 0; q < np; ++q) {
            if (p >= m && q >= m) continue;  // field-field (polarizability): not needed
            double v = A(p, q) + A(q, p) + pairs(p, q);
            for (std::size_t s = 0; s < nsets; ++s) {
                v += dot(Yk[p][s], U[q][s]) + dot(Yder[q][s], U[p][s]);
                Matrix T = mul(U[p][s], true, U[q][s], false, nt);  // U^p' U^q
                add(T, 1.0, mul(Smo[p][s], false, U[q][s], false, nt));
                add(T, 1.0, mul(Smo[q][s], false, U[p][s], false, nt));
                add(T, 1.0, mul(U[p][s], false, U[q][s], false, nt));
                v -= dot(P.Yset[s], T);
            }
            if (p < m && q < m) {
                for (std::size_t x = 0; x < n3; ++x) {
                    if (Q[p][x] == 0.0) continue;
                    for (std::size_t y = 0; y < n3; ++y) v += Q[p][x] * 0.5 * (Hx(x, y) + Hx(y, x)) * Q[q][y];
                }
            } else {
                const std::size_t k = p < m ? p : q, f = (p < m ? q : p) - m;
                for (std::size_t x = 0; x < n3; ++x) v += Q[k][x] * dip_expl[x][f];
            }
            B(p, q) = v;
        }
    for (std::size_t p = 0; p < np; ++p)
        for (std::size_t q = 0; q < p; ++q) {
            if (p < m && q < m) res.asymmetry = std::max(res.asymmetry, std::fabs(B(p, q) - B(q, p)));
            B(p, q) = B(q, p) = 0.5 * (B(p, q) + B(q, p));
        }
    res.response_seconds = seconds_since(tr);

    // ---- the Cartesian Hessian from the internal block and the invariance (translations, rotations)
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
    if (nf) {
        // atomic polar tensor: d mu / dq = -d2E/dq dF, completed with the translations and rotations
        const RigidMotion rd = rigid_motions(mol, res.dipole, &res.gradient);
        res.dipole_derivatives.assign(n3, {0.0, 0.0, 0.0});
        for (std::size_t x = 0; x < n3; ++x)
            for (std::size_t f = 0; f < 3; ++f) {
                double s = 0.0;
                for (std::size_t k = 0; k < m; ++k) s -= Q[k][x] * B(k, m + f);
                for (std::size_t r = 0; r < rd.directions.size(); ++r) s += rd.directions[r][x] * rd.dipole_derivatives[r][f];
                res.dipole_derivatives[x][f] = s;
            }
        res.has_dipole_derivatives = true;
    }
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace

OpenShellGradient open_shell_correlated_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                 const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                 const SODensities& d, int threads, double memory_mb,
                                                 const CPHFCallback& progress) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const Prepared P = prepare(mol, basis, scf, aa, bb, ab, d, nt, memory_mb, progress);
    OpenShellGradient g = gradient_of(mol, basis, P, nt);
    g.seconds = seconds_since(t0);
    return g;
}

CorrelatedHessianResult open_shell_correlated_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                      const MOIntegrals& aa, const MOIntegrals& bb,
                                                      const MOIntegrals& ab, SOResult& res, const SOOptions& options,
                                                      int threads, double memory_mb, bool dipoles,
                                                      const SOProgress& progress) {
    const int nt = resolve_threads(threads);
    const SODensities d = so_densities(res, options, progress);
    const SecondOrder second = [&](const std::vector<PairHam>& Hk, Matrix& A, Matrix& pairs, CorrelatedHessianResult& out) {
        const std::size_t np = Hk.size();
        std::vector<HamiltonianDerivative> pert(np);
        for (std::size_t k = 0; k < np; ++k) {
            pert[k].h[0] = Hk[k].h[0];
            pert[k].h[1] = Hk[k].h[1];
            const PairHam* hp = &Hk[k];
            pert[k].g = [hp](int sp, int sq, std::size_t a, std::size_t b, std::size_t c, std::size_t e) { return hp->g(sp, sq, a, b, c, e); };
        }
        const SOSecondOrder so2 = so_second_order(res, pert, options, progress);
        out.method_converged = so2.converged && d.lambda_converged;
        out.method_iterations = so2.iterations;
        out.method_residual = so2.residual;
        pairs = so2.pair;
        for (std::size_t p = 0; p < np; ++p) {
            report_progress(JobStage::Hessian, 6, 0.0, static_cast<double>(p + 1));
            for (std::size_t k = 0; k < np; ++k) A(p, k) = contract(so2.dD[k], Hk[p]);
        }
    };
    return hessian_core(mol, basis, scf, aa, bb, ab, d, second, nt, memory_mb, dipoles);
}

OpenShellGradient open_shell_fci_gradient(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                          const MOIntegrals& mo, const std::vector<double>& c, int threads,
                                          double memory_mb) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    report_progress(JobStage::Gradient, 0);
    const SODensities d = fci_spin_densities(mo.n, scf.n_alpha, scf.n_beta, c, c, nt);
    const Prepared P = prepare(mol, basis, scf, mo, mo, mo, d, nt, memory_mb, CPHFCallback(), true);
    OpenShellGradient g = gradient_of(mol, basis, P, nt);
    g.seconds = seconds_since(t0);
    return g;
}

CorrelatedHessianResult open_shell_fci_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                               const MOIntegrals& mo, const std::vector<double>& c,
                                               double electronic_energy, int threads, double memory_mb,
                                               double fci_memory_mb, bool dipoles) {
    const int nt = resolve_threads(threads);
    const int na = scf.n_alpha, nb = scf.n_beta;
    const SODensities d = fci_spin_densities(mo.n, na, nb, c, c, nt);
    const SecondOrder second = [&](const std::vector<PairHam>& Hk, Matrix& A, Matrix& pairs, CorrelatedHessianResult& out) {
        const std::size_t np = Hk.size();
        // the first-order Hamiltonians in the one orbital set as MO integrals
        auto ints = [&](std::size_t k) {
            MOIntegrals m;
            m.n = mo.n;
            m.h = symmetrized(Hk[k].h[0]);
            m.npair = Hk[k].P[0];
            m.stride = Hk[k].P[0];
            m.eri.assign(Hk[k].blk[0], Hk[k].blk[0] + m.npair * m.npair);
            return m;
        };
        const FCIResponse fr = fci_vector_responses(mo, na, nb, c, electronic_energy, np, ints, nt, fci_memory_mb, 1.0e-9, 300);
        out.method_converged = fr.converged;
        out.method_iterations = fr.iterations;
        out.method_residual = fr.residual;
        // variational: E^pq = ... + 2 <c^q|H^[p]|c>, A(p, q) = <c^q|H^[p]|c> (symmetrised), no pair term
        for (std::size_t q = 0; q < np; ++q) {
            report_progress(JobStage::Hessian, 6, 0.0, static_cast<double>(q + 1));
            const SODensities dq = fci_spin_densities(mo.n, na, nb, c, fr.vectors[q], nt);
            for (std::size_t p = 0; p < np; ++p) A(p, q) = contract(dq, Hk[p]);
        }
        pairs = Matrix(np, np);
    };
    return hessian_core(mol, basis, scf, mo, mo, mo, d, second, nt, memory_mb, dipoles, true);
}

double open_shell_fci_gradient_bytes(int n_alpha, int n_beta, int n_mo, int n_basis) {
    const double n = n_mo, N = n_basis;
    const double P = n * (n + 1) / 2.0, NP = N * (N + 1) / 2.0;
    auto binom = [](double a, int k) {
        double b = 1.0;
        for (int i = 0; i < k; ++i) b *= (a - i) / (i + 1.0);
        return b;
    };
    const double det = binom(n, n_alpha) * binom(n, n_beta);
    double b = 8.0 * 4.0 * det;                                  // CI vector (the Davidson space is the energy's)
    b += 8.0 * 2.0 * n * n * std::max(binom(n, n_alpha - 1) * binom(n, n_beta - 1),
                                       std::max(binom(n, n_alpha - 2) * binom(n, n_beta), binom(n, n_alpha) * binom(n, n_beta - 2)));
    b += 8.0 * 6.0 * n * n * n * n + 8.0 * 2.0 * P * P;         // densities, MO integrals
    b += 8.0 * NP * NP;                                          // AO two-particle density
    return b;
}

double open_shell_fci_hessian_bytes(int n_alpha, int n_beta, int n_mo, int n_basis, int n_atoms) {
    const double n = n_mo, N = n_basis;
    const double npert = std::max(0, 3 * n_atoms - 6) + 3.0;
    const double P = n * (n + 1) / 2.0, NP = N * (N + 1) / 2.0;
    auto binom = [](double a, int k) {
        double b = 1.0;
        for (int i = 0; i < k; ++i) b *= (a - i) / (i + 1.0);
        return b;
    };
    const double det = binom(n, n_alpha) * binom(n, n_beta);
    double b = 8.0 * (npert + 8.0) * det;                       // CI vector, responses, CG vectors
    b += 8.0 * 2.0 * n * n * std::max(binom(n, n_alpha - 1) * binom(n, n_beta - 1),
                                       std::max(binom(n, n_alpha - 2) * binom(n, n_beta), binom(n, n_alpha) * binom(n, n_beta - 2)));
    b += 8.0 * 6.0 * n * n * n * n;                             // densities (full and packed)
    b += 8.0 * (npert + 2.0) * P * P;                           // perturbed Hamiltonians (pair matrices)
    b += 8.0 * (3.0 * n_atoms * NP * (NP + 1) / 2.0 + NP * NP);  // derivative integrals, back-transformation
    return b;
}

OpenShellGradient open_shell_scf_gradient_check(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                int threads) {
    const int nt = resolve_threads(threads);
    const OrbitalModel model(scf, basis, nt, 1024.0);
    const SODensities d = scf_densities(model);
    const Prepared P = prepare(mol, basis, scf, aa, bb, ab, d, nt, 1024.0, CPHFCallback());
    return gradient_of(mol, basis, P, nt);
}

CorrelatedHessianResult open_shell_scf_hessian_check(const Molecule& mol, const BasisSet& basis, const RHFResult& scf,
                                                     const MOIntegrals& aa, const MOIntegrals& bb, const MOIntegrals& ab,
                                                     int threads, bool dipoles) {
    const int nt = resolve_threads(threads);
    const OrbitalModel model(scf, basis, nt, 1024.0);
    const SODensities d = scf_densities(model);
    return hessian_core(mol, basis, scf, aa, bb, ab, d, SecondOrder(), nt, 1024.0, dipoles);
}

double open_shell_hessian_bytes(const std::string& method, const int nocc[2], int n_mo, int n_basis, int n_atoms) {
    const double n = n_mo, N = n_basis;
    const double npert = std::max(0, 3 * n_atoms - 6) + 3.0;
    const double P = n * (n + 1) / 2.0, NP = N * (N + 1) / 2.0;
    const int nvir[2] = {n_mo - nocc[0], n_mo - nocc[1]};
    const SOSize so = estimate_so_method(method, nocc, nvir);
    double b = 2.0 * so.memory_bytes;                          // the method's graph, amplitudes and multipliers
    b += 8.0 * npert * so.amplitudes * 2.0;                    // amplitude responses
    b += 8.0 * npert * 3.0 * n * n * n * n;                    // density responses (full arrays)
    b += 8.0 * (npert + 1.0) * 3.0 * P * P;                    // perturbed Hamiltonians (pair matrices)
    b += 8.0 * 6.0 * n * n * n * n;                            // the densities (full and packed)
    b += 8.0 * (3.0 * n_atoms * NP * (NP + 1) / 2.0 + NP * NP);  // derivative integrals, back-transformation
    return b;
}

}  // namespace macqc
