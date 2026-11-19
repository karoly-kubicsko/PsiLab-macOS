// ΨLab 15.0: analytic Hessians of the open-shell SCF energies (UHF, ROHF), and dipole derivatives.
//
// The orbitals are grouped in sets s with occupation numbers n^a_s, n^b_s (0 or 1) of the alpha and
// beta electrons in each orbital: UHF has two sets (the alpha orbitals with n^a = 1 for the
// occupied ones, the beta orbitals with n^b = 1), ROHF one set (n^a = n^b = 1 doubly occupied,
// n^a = 1 singly occupied). With P_s' = sum_s C_s N^s'_s C_s' the spin densities, F_s' the Fock
// matrices and the Lagrangian eps_s = sum_s' F^s'_s N^s'_s (F^s'_s = C_s' F_s' C_s), the gradient
//   dE/dx = sum P h^x + sum Gamma (mn|ls)^x - sum W S^x + V_nn^x,  W = sum_s C_s eps_s C_s',
// holds at every geometry. Its derivative along y, with the orbitals C_s(y) = C_s (1 + U^y_s):
//   E^xy = [explicit second derivatives] + sum_s sum_s' 2 tr(U_s N^s'_s F~^s'(x)_s)
//          - sum_s [tr(U_s eps_s S~^x_s) + tr(eps_s U_s' S~^x_s) + tr(eps^y_s S~^x_s)],
//   eps^y_s = sum_s' (U_s' F^s'_s + F^s'_s U_s + C_s' F_s'^(y,total) C_s) N^s'_s,
// F_s'^(x) = h^x + J^x(P) - K^x(P_s') the skeleton Fock derivative, F^(y,total) = F^(y) + J(P^y)
// - K(P^y_s'), S~^x_s = C_s' S^x C_s. U_s = K_s - S~^y_s / 2 with K_s antisymmetric and non-zero only
// for the pairs whose occupations differ (the non-redundant rotations); K solves the linearised
// stationarity condition antisym(eps^y_s) = 0 on those pairs: UHF the coupled-perturbed UHF
// equations, ROHF the Roothaan conditions F_b(c,o) = F_a(o,v) = (F_a + F_b)(c,v) = 0 to first order.
// For RHF the expression is that of rhf_derivatives.h. The energy is invariant to the redundant
// rotations, so their choice (-S~/2) does not change E^xy.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/scf_derivatives.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t0) { return std::chrono::duration<double>(Clock::now() - t0).count(); }

Matrix mul(const Matrix& A, bool ta, const Matrix& B, bool tb, int threads = 1) {
    const std::size_t m = ta ? A.cols() : A.rows();
    const std::size_t k = ta ? A.rows() : A.cols();
    const std::size_t n = tb ? B.rows() : B.cols();
    Matrix C(m, n);
    if (m == 0 || n == 0 || k == 0) return C;
    gemm(ta, tb, m, n, k, 1.0, A.data(), A.cols(), B.data(), B.cols(), 0.0, C.data(), n, threads);
    return C;
}

void axpy(Matrix& y, double a, const Matrix& x) {
    for (std::size_t i = 0; i < y.rows(); ++i)
        for (std::size_t j = 0; j < y.cols(); ++j) y(i, j) += a * x(i, j);
}

double dot(const Matrix& a, const Matrix& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.rows(); ++i)
        for (std::size_t j = 0; j < a.cols(); ++j) s += a(i, j) * b(i, j);
    return s;
}

// tr(A B) = sum_ij A_ij B_ji
double trace_product(const Matrix& A, const Matrix& B) {
    double s = 0.0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        for (std::size_t j = 0; j < A.cols(); ++j) s += A(i, j) * B(j, i);
    return s;
}

// M with its columns scaled by n
Matrix scale_columns(const Matrix& M, const std::vector<double>& n) {
    Matrix out = M;
    for (std::size_t i = 0; i < M.rows(); ++i)
        for (std::size_t j = 0; j < M.cols(); ++j) out(i, j) *= n[j];
    return out;
}

bool solve_dense(std::vector<double> a, std::vector<double> b, std::size_t n, std::vector<double>& x) {
    for (std::size_t col = 0; col < n; ++col) {
        std::size_t piv = col;
        for (std::size_t r = col + 1; r < n; ++r)
            if (std::fabs(a[r * n + col]) > std::fabs(a[piv * n + col])) piv = r;
        if (std::fabs(a[piv * n + col]) < 1.0e-300) return false;
        if (piv != col) {
            for (std::size_t c = 0; c < n; ++c) std::swap(a[piv * n + c], a[col * n + c]);
            std::swap(b[piv], b[col]);
        }
        for (std::size_t r = col + 1; r < n; ++r) {
            const double f = a[r * n + col] / a[col * n + col];
            if (f == 0.0) continue;
            for (std::size_t c = col; c < n; ++c) a[r * n + c] -= f * a[col * n + c];
            b[r] -= f * b[col];
        }
    }
    x.assign(n, 0.0);
    for (std::size_t r = n; r-- > 0;) {
        double s = b[r];
        for (std::size_t c = r + 1; c < n; ++c) s -= a[r * n + c] * x[c];
        x[r] = s / a[r * n + r];
    }
    return true;
}

struct OrbitalSet {
    Matrix C;                         // N x n
    std::vector<double> occ[2];       // alpha, beta occupation of every orbital
    std::vector<std::pair<std::size_t, std::size_t>> pairs;  // non-redundant (p > q)
    Matrix F[2];                      // F_s' in this set's MOs
    Matrix eps;                       // sum_s' F^s' N^s'
};

}  // namespace

// Krylov (Galerkin) solver for several symmetric linear systems A x = b at once, shared by the
// open-shell response equations: `apply` maps a batch of vectors to A times them (one call per
// iteration, so the Fock builds of all perturbations can share one pass over the integrals).
struct LinearBatchSolution {
    std::vector<std::vector<double>> x;
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
};

LinearBatchSolution solve_symmetric_batch(
    const std::function<std::vector<std::vector<double>>(const std::vector<std::vector<double>>&)>& apply,
    const std::vector<double>& diagonal, const std::vector<std::vector<double>>& rhs, double tol, int max_iter,
    const CPHFCallback& progress) {
    const std::size_t nx = rhs.size(), dim = diagonal.size();
    auto vdot = [](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
        return s;
    };
    auto precondition = [&](const std::vector<double>& r) {
        std::vector<double> z(dim);
        for (std::size_t i = 0; i < dim; ++i) z[i] = r[i] / (std::fabs(diagonal[i]) > 1.0e-3 ? diagonal[i] : 1.0e-3);
        return z;
    };
    LinearBatchSolution sol;
    sol.x.assign(nx, std::vector<double>(dim, 0.0));
    std::vector<std::vector<std::vector<double>>> V(nx), AV(nx);
    std::vector<std::vector<double>> next(nx);
    std::vector<bool> done(nx, false);
    for (std::size_t x = 0; x < nx; ++x) {
        if (std::sqrt(vdot(rhs[x], rhs[x])) < 1.0e-14) done[x] = true;
        else next[x] = precondition(rhs[x]);
    }
    double worst = 0.0;
    for (int it = 1; it <= max_iter; ++it) {
        check_cancelled();
        std::vector<std::size_t> active;
        std::vector<std::vector<double>> batch;
        for (std::size_t x = 0; x < nx; ++x) {
            if (done[x]) continue;
            std::vector<double> w = next[x];
            for (int pass = 0; pass < 2; ++pass)
                for (const auto& q : V[x]) {
                    const double c = vdot(q, w);
                    for (std::size_t i = 0; i < dim; ++i) w[i] -= c * q[i];
                }
            const double nw = std::sqrt(vdot(w, w));
            if (nw < 1.0e-14) {
                done[x] = true;
                continue;
            }
            for (double& e : w) e /= nw;
            V[x].push_back(w);
            active.push_back(x);
            batch.push_back(std::move(w));
        }
        if (active.empty()) break;
        const std::vector<std::vector<double>> Ab = apply(batch);
        for (std::size_t k = 0; k < active.size(); ++k) AV[active[k]].push_back(Ab[k]);
        worst = 0.0;
        for (std::size_t x = 0; x < nx; ++x) {
            if (V[x].empty()) continue;
            const std::size_t m = V[x].size();
            std::vector<double> M(m * m), b(m), c;
            for (std::size_t i = 0; i < m; ++i) {
                b[i] = vdot(V[x][i], rhs[x]);
                for (std::size_t j = 0; j < m; ++j)
                    M[i * m + j] = 0.5 * (vdot(V[x][i], AV[x][j]) + vdot(V[x][j], AV[x][i]));
            }
            if (!solve_dense(M, b, m, c)) throw std::runtime_error("orbital response: singular subspace matrix");
            std::vector<double> xs(dim, 0.0), r = rhs[x];
            for (std::size_t i = 0; i < m; ++i)
                for (std::size_t e = 0; e < dim; ++e) {
                    xs[e] += c[i] * V[x][i][e];
                    r[e] -= c[i] * AV[x][i][e];
                }
            sol.x[x] = xs;
            const double rn = std::sqrt(vdot(r, r));
            worst = std::max(worst, rn);
            if (rn < tol || done[x]) done[x] = true;
            else next[x] = precondition(r);
        }
        sol.iterations = it;
        sol.residual = worst;
        report_progress(JobStage::CPHF, it, 0.0, worst);
        if (progress) progress(CPHFProgress{it, worst});
        if (std::all_of(done.begin(), done.end(), [](bool d) { return d; })) break;
    }
    sol.converged = worst < tol;
    return sol;
}

RHFHessianResult scf_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& scf, int threads,
                             double memory_mb, double cphf_tol, const CPHFCallback& progress) {
    if (!scf.open_shell()) return rhf_hessian(mol, basis, scf, threads, memory_mb, cphf_tol, progress);
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int nat = static_cast<int>(mol.atoms.size());
    const std::size_t n3 = static_cast<std::size_t>(3 * nat);
    const std::size_t N = static_cast<std::size_t>(basis.n_functions);
    const Matrix Fs[2] = {scf.Fa, scf.Fb};
    const Matrix Ps[2] = {scf.Pa, scf.Pb};

    // ---- orbital sets
    std::vector<OrbitalSet> sets;
    auto make_set = [&](const Matrix& C, int na, int nb) {
        OrbitalSet s;
        s.C = C;
        const std::size_t n = C.cols();
        s.occ[0].assign(n, 0.0);
        s.occ[1].assign(n, 0.0);
        for (int i = 0; i < na; ++i) s.occ[0][static_cast<std::size_t>(i)] = 1.0;
        for (int i = 0; i < nb; ++i) s.occ[1][static_cast<std::size_t>(i)] = 1.0;
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < p; ++q)
                if (s.occ[0][p] != s.occ[0][q] || s.occ[1][p] != s.occ[1][q]) s.pairs.emplace_back(p, q);
        s.eps = Matrix(n, n);
        for (int sp = 0; sp < 2; ++sp) {
            s.F[sp] = mul(mul(C, true, Fs[sp], false, nt), false, C, false, nt);
            axpy(s.eps, 1.0, scale_columns(s.F[sp], s.occ[sp]));
        }
        sets.push_back(std::move(s));
    };
    if (scf.reference == "UHF") {
        make_set(scf.C, scf.n_alpha, 0);
        make_set(scf.Cb, 0, scf.n_beta);
    } else {
        make_set(scf.C, scf.n_alpha, scf.n_beta);
    }
    std::size_t dim = 0;
    for (const OrbitalSet& s : sets) dim += s.pairs.size();

    RHFHessianResult res;
    // ---- explicit second derivatives and skeleton first derivatives
    report_progress(JobStage::Hessian, 1);
    const std::vector<Matrix> Sx = overlap_derivative_matrices(basis, nat, nt);
    const std::vector<Matrix> hx = core_hamiltonian_derivative_matrices(basis, mol, nt);
    const std::vector<std::vector<Matrix>> Gx =
        eri_derivative_fock_multi(basis, nat, {scf.P, scf.Pa, scf.Pb}, {1.0, 0.0, 0.0}, {0.0, -1.0, -1.0}, nt);
    report_progress(JobStage::Hessian, 2);
    const Matrix W = scf_energy_weighted_density(scf);
    Matrix H = core_hamiltonian_hessian(basis, mol, scf.P, nt);
    axpy(H, 1.0, eri_hessian(basis, nat, scf_two_particle_density(scf), nt));
    axpy(H, -1.0, overlap_hessian(basis, nat, W, nt));
    axpy(H, 1.0, nuclear_repulsion_hessian(mol));
    res.integral_seconds = seconds_since(t0);
    const auto t1 = Clock::now();

    // skeleton Fock derivatives (AO) per spin and per coordinate, and S~ in each set's MOs
    std::vector<std::array<Matrix, 2>> Fskel(n3);
    std::vector<std::vector<Matrix>> Smo(sets.size(), std::vector<Matrix>(n3)), Fmo_x(sets.size() * 2, std::vector<Matrix>(n3));
    for (std::size_t x = 0; x < n3; ++x) {
        for (int sp = 0; sp < 2; ++sp) {
            Fskel[x][static_cast<std::size_t>(sp)] = hx[x];
            axpy(Fskel[x][static_cast<std::size_t>(sp)], 1.0, Gx[0][x]);
            axpy(Fskel[x][static_cast<std::size_t>(sp)], 1.0, Gx[1 + static_cast<std::size_t>(sp)][x]);
        }
        for (std::size_t s = 0; s < sets.size(); ++s) {
            const Matrix& C = sets[s].C;
            Smo[s][x] = mul(mul(C, true, Sx[x], false, nt), false, C, false, nt);
            for (int sp = 0; sp < 2; ++sp)
                Fmo_x[2 * s + static_cast<std::size_t>(sp)][x] =
                    mul(mul(C, true, Fskel[x][static_cast<std::size_t>(sp)], false, nt), false, C, false, nt);
        }
    }

    FockBuilder fock(basis, nt, memory_mb);
    // Spin densities of orbital responses U_s (one matrix per set): P_s' = sum_s C (U N + N U') C'.
    auto response_densities = [&](const std::vector<Matrix>& U) {
        std::array<Matrix, 2> D{Matrix(N, N), Matrix(N, N)};
        for (std::size_t s = 0; s < sets.size(); ++s)
            for (int sp = 0; sp < 2; ++sp) {
                const Matrix CU = mul(sets[s].C, false, scale_columns(U[s], sets[s].occ[sp]), false, nt);  // C U N
                const Matrix A = mul(CU, false, sets[s].C, true, nt);
                Matrix& d = D[static_cast<std::size_t>(sp)];
                for (std::size_t i = 0; i < N; ++i)
                    for (std::size_t j = 0; j < N; ++j) d(i, j) += A(i, j) + A(j, i);
            }
        return D;
    };
    // The two-electron Fock responses G_s'(D) = J(D_a + D_b) - K(D_s') of several density pairs.
    auto fock_responses = [&](const std::vector<std::array<Matrix, 2>>& D) {
        std::vector<Matrix> in;
        for (const auto& d : D) {
            in.push_back(d[0]);
            in.push_back(d[1]);
        }
        std::vector<Matrix> J, K;
        fock.build_jk(in, J, K);
        std::vector<std::array<Matrix, 2>> G(D.size());
        for (std::size_t k = 0; k < D.size(); ++k) {
            Matrix Jt = J[2 * k];
            axpy(Jt, 1.0, J[2 * k + 1]);
            for (int sp = 0; sp < 2; ++sp) {
                G[k][static_cast<std::size_t>(sp)] = Jt;
                axpy(G[k][static_cast<std::size_t>(sp)], -1.0, K[2 * k + static_cast<std::size_t>(sp)]);
            }
        }
        return G;
    };
    // eps^(1)_s of a response U with total Fock derivative Ft (AO, per spin)
    auto eps_response = [&](std::size_t s, const Matrix& U, const std::array<Matrix, 2>& Ft) {
        const OrbitalSet& o = sets[s];
        Matrix e(o.C.cols(), o.C.cols());
        for (int sp = 0; sp < 2; ++sp) {
            Matrix M = mul(U, true, o.F[sp], false, nt);
            axpy(M, 1.0, mul(o.F[sp], false, U, false, nt));
            axpy(M, 1.0, mul(mul(o.C, true, Ft[static_cast<std::size_t>(sp)], false, nt), false, o.C, false, nt));
            axpy(e, 1.0, scale_columns(M, o.occ[sp]));
        }
        return e;
    };
    auto unpack = [&](const std::vector<double>& k) {  // antisymmetric K per set
        std::vector<Matrix> K;
        std::size_t pos = 0;
        for (const OrbitalSet& o : sets) {
            Matrix m(o.C.cols(), o.C.cols());
            for (const auto& pq : o.pairs) {
                m(pq.first, pq.second) = k[pos];
                m(pq.second, pq.first) = -k[pos];
                ++pos;
            }
            K.push_back(m);
        }
        return K;
    };
    auto antisym = [&](const std::vector<Matrix>& e) {  // (e_pq - e_qp) on the pairs
        std::vector<double> v;
        v.reserve(dim);
        for (std::size_t s = 0; s < sets.size(); ++s)
            for (const auto& pq : sets[s].pairs) v.push_back(e[s](pq.first, pq.second) - e[s](pq.second, pq.first));
        return v;
    };
    // The linear map K -> antisym(eps^(1)(K)) for a batch (no skeleton term).
    auto apply = [&](const std::vector<std::vector<double>>& batch) {
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
    };
    std::vector<double> diagonal;
    for (const OrbitalSet& o : sets)
        for (const auto& pq : o.pairs) {
            double d = 0.0;
            for (int sp = 0; sp < 2; ++sp)
                d += (o.F[sp](pq.first, pq.first) - o.F[sp](pq.second, pq.second)) *
                     (o.occ[sp][pq.second] - o.occ[sp][pq.first]);
            diagonal.push_back(d);
        }

    // right-hand sides: the response to U0 = -S~/2 plus the skeleton Fock derivative
    std::vector<std::vector<Matrix>> U0(n3);
    std::vector<std::array<Matrix, 2>> D0(n3);
    for (std::size_t x = 0; x < n3; ++x) {
        for (std::size_t s = 0; s < sets.size(); ++s) {
            Matrix u = Smo[s][x];
            for (std::size_t i = 0; i < u.rows(); ++i)
                for (std::size_t j = 0; j < u.cols(); ++j) u(i, j) *= -0.5;
            U0[x].push_back(u);
        }
        D0[x] = response_densities(U0[x]);
    }
    const std::vector<std::array<Matrix, 2>> G0 = fock_responses(D0);
    std::vector<std::vector<double>> rhs(n3);
    for (std::size_t x = 0; x < n3; ++x) {
        std::array<Matrix, 2> Ft = G0[x];
        for (int sp = 0; sp < 2; ++sp) axpy(Ft[static_cast<std::size_t>(sp)], 1.0, Fskel[x][static_cast<std::size_t>(sp)]);
        std::vector<Matrix> e;
        for (std::size_t s = 0; s < sets.size(); ++s) e.push_back(eps_response(s, U0[x][s], Ft));
        rhs[x] = antisym(e);
        for (double& v : rhs[x]) v = -v;
    }
    report_progress(JobStage::CPHF, 0);
    const LinearBatchSolution sol = solve_symmetric_batch(apply, diagonal, rhs, cphf_tol, 80, progress);
    res.cphf_converged = sol.converged;
    res.cphf_iterations = sol.iterations;
    res.cphf_residual = sol.residual;

    // ---- full responses U^y, density responses and eps^y
    std::vector<std::vector<Matrix>> U(n3);
    std::vector<std::array<Matrix, 2>> Py(n3);
    for (std::size_t y = 0; y < n3; ++y) {
        const std::vector<Matrix> K = unpack(sol.x[y]);
        for (std::size_t s = 0; s < sets.size(); ++s) {
            Matrix u = U0[y][s];
            axpy(u, 1.0, K[s]);
            U[y].push_back(u);
        }
        Py[y] = response_densities(U[y]);
    }
    const std::vector<std::array<Matrix, 2>> Gy = fock_responses(Py);
    for (std::size_t y = 0; y < n3; ++y) {
        check_cancelled();
        std::array<Matrix, 2> Ft = Gy[y];
        for (int sp = 0; sp < 2; ++sp) axpy(Ft[static_cast<std::size_t>(sp)], 1.0, Fskel[y][static_cast<std::size_t>(sp)]);
        for (std::size_t s = 0; s < sets.size(); ++s) {
            const OrbitalSet& o = sets[s];
            const Matrix& u = U[y][s];
            Matrix M = eps_response(s, u, Ft);              // eps^y
            axpy(M, 1.0, mul(u, false, o.eps, false, nt));  // + U eps
            axpy(M, 1.0, mul(o.eps, false, u, true, nt));   // + eps U'
            const Matrix uN[2] = {scale_columns(u, o.occ[0]), scale_columns(u, o.occ[1])};
            for (std::size_t x = 0; x < n3; ++x) {
                double v = 0.0;
                for (int sp = 0; sp < 2; ++sp)
                    v += 2.0 * trace_product(uN[sp], Fmo_x[2 * s + static_cast<std::size_t>(sp)][x]);
                v -= trace_product(M, Smo[s][x]);
                H(x, y) += v;
            }
        }
    }
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t y = 0; y < x; ++y) H(x, y) = H(y, x) = 0.5 * (H(x, y) + H(y, x));
    res.hessian = H;

    // ---- dipole moment and its derivatives
    const std::array<double, 3> origin{0.0, 0.0, 0.0};
    const std::vector<Matrix> r = dipole_matrices(basis, origin);
    res.dipole = dipole_moment(mol, basis, scf.P, origin);
    const std::vector<std::vector<Matrix>> rx = dipole_derivative_matrices(basis, nat, origin, nt);
    res.dipole_derivatives.assign(n3, {0.0, 0.0, 0.0});
    for (std::size_t x = 0; x < n3; ++x)
        for (int m = 0; m < 3; ++m) {
            const Matrix& R = r[static_cast<std::size_t>(m)];
            double val = -dot(Py[x][0], R) - dot(Py[x][1], R) - dot(scf.P, rx[x][static_cast<std::size_t>(m)]);
            if (static_cast<int>(x % 3) == m) val += mol.atoms[x / 3].nuclear_charge();
            res.dipole_derivatives[x][static_cast<std::size_t>(m)] = val;
        }
    (void)Ps;
    res.cphf_seconds = seconds_since(t1);
    res.seconds = seconds_since(t0);
    return res;
}

}  // namespace macqc
