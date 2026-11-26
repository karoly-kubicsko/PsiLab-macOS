// ΨLab 20.0/21.0: nuclear motion (see nuclear_motion.h). The Hamiltonian, its DVR and the solvers
// follow the program cpp_DOPI_expval_array.cpp (reference/nuclear_motion); the line references
// "(original: ...)" name its functions.
#include "macqc/nuclear_motion.h"

#include <algorithm>
#include <atomic>
#include <cstdarg>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "macqc/job.h"
#include "macqc/linalg.h"
#include "macqc/parallel.h"

namespace macqc {
namespace {

// Conversion factors of the original program (masses in u, m_e per u) and CODATA 2018 for the
// output units.
constexpr double kAmuToMe = 1.8228883e3;            // original: amtoau
constexpr double kHartreeToCm = 219474.6313632;      // CODATA 2018 (original: 219474.630)
constexpr double kBohrToAngstrom = 0.529177210903;   // CODATA 2018 (original: 0.529177249)
constexpr double kPi = 3.14159265358979323846;

// Nuclear masses (u) of the original program.
constexpr double kMassH = 1.00727647;
constexpr double kMassD = 2.013553214;
constexpr double kMass16O = 15.990526;
constexpr double kMass17O = 16.9947425;
constexpr double kMass18O = 17.999161;

struct IsotopologueData {
    const char* name;
    double hydrogen;  // u
    double oxygen;    // u
    double vmin;      // E_h: the minimum of the CVRQD surface (symmetric, C2v)
    double rmin;      // bohr
    double amin;      // degrees
};

// The minima were found by minimising cvrqd_potential (r1 = r2 and the angle; checked by the
// test suite). The original program subtracted 0.00358028360318797 (H2(16)O),
// 0.00343991010054165 (H2(17)O) and 0.00337543999746866 (D2(16)O): at most 3.2e-8 E_h
// (0.007 cm-1) from these minima.
IsotopologueData isotopologue_data(CvrqdIsotopologue iso) {
    switch (iso) {
    case CvrqdIsotopologue::H2_16O:
        return {"H2(16)O", kMassH, kMass16O, 0.00358030713628730351, 1.81008106, 104.499962};
    case CvrqdIsotopologue::H2_17O:
        return {"H2(17)O", kMassH, kMass17O, 0.00343994205004556847, 1.81008133, 104.500211};
    case CvrqdIsotopologue::D2_16O:
        return {"D2(16)O", kMassD, kMass16O, 0.00337544075604669785, 1.81004459, 104.490277};
    case CvrqdIsotopologue::H2_18O:
        return {"H2(18)O", kMassH, kMass18O, 0.00331573097906425179, 1.81008158, 104.500431};
    }
    throw std::invalid_argument("unknown isotopologue");
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::string fmt(const char* f, ...) __attribute__((format(printf, 1, 2)));
std::string fmt(const char* f, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, f);
    std::vsnprintf(buf, sizeof buf, f, ap);
    va_end(ap);
    return buf;
}

// ---------------------------------------------------------------------------------------------
// Dense vector helpers

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    const std::size_t n = a.size();
    for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
    return s;
}

double norm(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

void axpy(double alpha, const std::vector<double>& x, std::vector<double>& y) {
    const std::size_t n = x.size();
    for (std::size_t i = 0; i < n; ++i) y[i] += alpha * x[i];
}

void scale(std::vector<double>& x, double s) {
    for (double& v : x) v *= s;
}

// ---------------------------------------------------------------------------------------------
// Sparse symmetric matrix (CSR, all entries of every row, columns ascending)

struct Csr {
    std::size_t n = 0;
    std::vector<std::size_t> row;       // n + 1
    std::vector<int> col;
    std::vector<double> val;
    std::vector<std::size_t> diag_pos;  // position of the diagonal entry of each row
    std::vector<double> diag;
    int threads = 1;

    void multiply(const std::vector<double>& x, std::vector<double>& y) const {
        y.assign(n, 0.0);
        const int t = (n > 20000 && threads > 1) ? threads : 1;
        run_parallel(t, [&](int k) {
            const std::size_t a = n * static_cast<std::size_t>(k) / t, b = n * static_cast<std::size_t>(k + 1) / t;
            for (std::size_t i = a; i < b; ++i) {
                double s = 0.0;
                for (std::size_t p = row[i]; p < row[i + 1]; ++p) s += val[p] * x[static_cast<std::size_t>(col[p])];
                y[i] = s;
            }
        });
    }
};

// ---------------------------------------------------------------------------------------------
// Preconditioners (original: CGM with precondtype 1, 2, 3; ILUprecond; 21.0 also
// ILUprecondfordavidson and ILUprecondforjacobidavidson, i.e. the same factors of H - theta)

enum class Precond { None, Jacobi, Ilu };

Precond parse_precond(const std::string& s) {
    if (s == "none") return Precond::None;
    if (s == "jacobi" || s == "diagonal") return Precond::Jacobi;
    if (s == "ilu" || s == "sgs") return Precond::Ilu;
    throw std::invalid_argument("preconditioner must be none, jacobi or ilu");
}

// M ~ H - shift. ILU: the factors of the original, M = (D' + L) D'^-1 (D' + U) with the strict
// triangles L, U of H and the diagonal D' = D - shift (symmetric Gauss-Seidel); Jacobi: M = D'.
// absolute: D' = |D - shift| (floored), which keeps M symmetric positive definite, as conjugate
// gradients, MINRES and the Chebyshev iteration require; for shift = 0 both are D (H is positive
// definite), the preconditioner of the original.
struct Preconditioner {
    const Csr* a = nullptr;
    Precond kind = Precond::None;
    std::vector<double> d;

    void setup(const Csr& m, Precond k, double shift, bool absolute);
    void apply(const std::vector<double>& r, std::vector<double>& z) const;
};

void Preconditioner::setup(const Csr& m, Precond k, double shift, bool absolute) {
    a = &m;
    kind = k;
    if (k == Precond::None) return;
    d.resize(m.n);
    for (std::size_t i = 0; i < m.n; ++i) {
        double v = m.diag[i] - shift;
        if (absolute) v = std::fabs(v);
        if (std::fabs(v) < 1e-8) v = v < 0 ? -1e-8 : 1e-8;
        d[i] = v;
    }
}

void Preconditioner::apply(const std::vector<double>& r, std::vector<double>& z) const {
    switch (kind) {
    case Precond::None:
        z = r;
        return;
    case Precond::Jacobi: {
        const std::size_t n = r.size();
        z.resize(n);
        for (std::size_t i = 0; i < n; ++i) z[i] = r[i] / d[i];
        return;
    }
    case Precond::Ilu: {
        const Csr& m = *a;
        const std::size_t n = m.n;
        z.resize(n);
        std::vector<double> w(n);
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t q = m.row[i]; q < m.diag_pos[i]; ++q) s += m.val[q] * w[static_cast<std::size_t>(m.col[q])];
            w[i] = (r[i] - s) / d[i];
        }
        for (std::size_t ii = n; ii-- > 0;) {
            double s = 0.0;
            for (std::size_t q = m.diag_pos[ii] + 1; q < m.row[ii + 1]; ++q) s += m.val[q] * z[static_cast<std::size_t>(m.col[q])];
            z[ii] = w[ii] - s / d[ii];
        }
        return;
    }
    }
}

// ---------------------------------------------------------------------------------------------
// Linear solvers for A x = b, A symmetric (H, H - theta, or that projected onto the complement of
// converged eigenvectors), x0 = 0 (original: CGM, minres, minresfordavidson,
// minresforjacobidavidson, CGMfordavidson, chebyshev, BiCGM). `op(v, y)`: y = A v; `prec(r, z)`:
// z = M^-1 r. A trace receives ||r|| (MINRES: its estimate of the preconditioned residual) after
// every iteration, from the start.

enum class Linear { Pcg, Minres, Chebyshev, Bicg };

Linear parse_linear(const std::string& s) {
    if (s == "pcg" || s == "cg") return Linear::Pcg;
    if (s == "minres") return Linear::Minres;
    if (s == "chebyshev") return Linear::Chebyshev;
    if (s == "bicg" || s == "bicgm") return Linear::Bicg;
    throw std::invalid_argument("the linear solver must be pcg, minres, chebyshev or bicg");
}

using VecOp = std::function<void(const std::vector<double>&, std::vector<double>&)>;

struct LinearStats {
    int iterations = 0;
    double residual = 0.0;
    bool converged = false;
    bool breakdown = false;  // conjugate gradients met p^T A p <= 0 (A not positive definite) or BiCG a zero
};

// Preconditioned conjugate gradients; stops when ||b - A x|| <= tol. strict: A must be positive
// definite (H); otherwise (a shifted A) a breakdown ends the solve with the current iterate.
LinearStats pcg(std::size_t n, const VecOp& op, const VecOp& prec, const std::vector<double>& b, std::vector<double>& x,
                double tol, int maxit, bool strict, std::vector<double>* trace) {
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r = b, z, p, q(n);
    prec(r, z);
    p = z;
    double rho = dot(r, z);
    st.residual = norm(r);
    if (trace) trace->push_back(st.residual);
    if (st.residual <= tol) {
        st.converged = true;
        return st;
    }
    for (int it = 1; it <= maxit; ++it) {
        if ((it & 63) == 0) check_cancelled();
        op(p, q);
        const double pq = dot(p, q);
        if (!(pq > 0.0)) {
            if (strict) throw std::runtime_error("conjugate gradients: the Hamiltonian is not positive definite");
            st.breakdown = true;
            return st;
        }
        const double alpha = rho / pq;
        axpy(alpha, p, x);
        axpy(-alpha, q, r);
        st.iterations = it;
        st.residual = norm(r);
        if (trace) trace->push_back(st.residual);
        if (st.residual <= tol) {
            st.converged = true;
            return st;
        }
        prec(r, z);
        const double rho_new = dot(r, z);
        if (!(std::fabs(rho) > 0.0)) {
            st.breakdown = true;
            return st;
        }
        const double beta = rho_new / rho;
        rho = rho_new;
        for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
    return st;
}

// MINRES (Paige and Saunders, SIAM J. Numer. Anal. 12, 617 (1975)) for a symmetric operator,
// indefinite allowed, with a symmetric positive definite preconditioner. Stops when the
// preconditioned residual estimate falls below tol.
LinearStats minres(std::size_t n, const VecOp& op, const VecOp& prec, const std::vector<double>& b, std::vector<double>& x,
                   double tol, int maxit, std::vector<double>* trace) {
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r1 = b, y, r2, v, w(n, 0.0), w1(n, 0.0), w2(n, 0.0);
    prec(r1, y);
    double beta1 = dot(r1, y);
    if (beta1 < 0.0) throw std::runtime_error("MINRES: the preconditioner is not positive definite");
    beta1 = std::sqrt(beta1);
    st.residual = beta1;
    if (trace) trace->push_back(beta1);
    if (beta1 <= tol || beta1 == 0.0) {
        st.converged = true;
        return st;
    }
    r2 = r1;
    double oldb = 0.0, beta = beta1, dbar = 0.0, epsln = 0.0, phibar = beta1, cs = -1.0, sn = 0.0;
    for (int it = 1; it <= maxit; ++it) {
        if ((it & 63) == 0) check_cancelled();
        const double s = 1.0 / beta;
        v = y;
        scale(v, s);
        op(v, y);
        if (it >= 2) axpy(-beta / oldb, r1, y);
        const double alfa = dot(v, y);
        axpy(-alfa / beta, r2, y);
        r1 = r2;
        r2 = y;
        prec(r2, y);
        oldb = beta;
        double b2 = dot(r2, y);
        if (b2 < 0.0) throw std::runtime_error("MINRES: the preconditioner is not positive definite");
        beta = std::sqrt(b2);
        const double oldeps = epsln;
        const double delta = cs * dbar + sn * alfa;
        const double gbar = sn * dbar - cs * alfa;
        epsln = sn * beta;
        dbar = -cs * beta;
        const double gamma = std::max(std::sqrt(gbar * gbar + beta * beta), std::numeric_limits<double>::min());
        cs = gbar / gamma;
        sn = beta / gamma;
        const double phi = cs * phibar;
        phibar = sn * phibar;
        const double denom = 1.0 / gamma;
        w1 = w2;
        w2 = w;
        for (std::size_t i = 0; i < n; ++i) w[i] = (v[i] - oldeps * w1[i] - delta * w2[i]) * denom;
        axpy(phi, w, x);
        st.iterations = it;
        st.residual = phibar;
        if (trace) trace->push_back(phibar);
        if (phibar <= tol) {
            st.converged = true;
            return st;
        }
        if (beta == 0.0) break;
    }
    return st;
}

// Preconditioned Chebyshev iteration (Saad, Iterative Methods for Sparse Linear Systems, 2nd ed.,
// Algorithm 12.1) for the eigenvalues of M^-1 A in [lo, hi], 0 < lo < hi: no inner products but the
// residual norm. The original's chebyshev without the misprint of its first step (alpha = 2/d).
LinearStats chebyshev(std::size_t n, const VecOp& op, const VecOp& prec, const std::vector<double>& b, std::vector<double>& x,
                      double tol, int maxit, double lo, double hi, std::vector<double>* trace) {
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r = b, z, d, q;
    st.residual = norm(r);
    if (trace) trace->push_back(st.residual);
    if (st.residual <= tol) {
        st.converged = true;
        return st;
    }
    const double theta = 0.5 * (hi + lo), delta = 0.5 * (hi - lo);
    const double sigma1 = theta / delta;
    double rho = 1.0 / sigma1;
    prec(r, z);
    d = z;
    scale(d, 1.0 / theta);
    for (int it = 1; it <= maxit; ++it) {
        if ((it & 63) == 0) check_cancelled();
        axpy(1.0, d, x);
        op(d, q);
        axpy(-1.0, q, r);
        st.iterations = it;
        st.residual = norm(r);
        if (trace) trace->push_back(st.residual);
        if (!std::isfinite(st.residual)) {
            st.breakdown = true;
            return st;
        }
        if (st.residual <= tol) {
            st.converged = true;
            return st;
        }
        const double rho_new = 1.0 / (2.0 * sigma1 - rho);
        prec(r, z);
        for (std::size_t i = 0; i < n; ++i) d[i] = rho_new * rho * d[i] + (2.0 * rho_new / delta) * z[i];
        rho = rho_new;
    }
    return st;
}

// Preconditioned biconjugate gradients (the original's BiCGM, there without its preconditioner),
// shadow residual r~0 = r0. For a symmetric A and M this reproduces the iterates of conjugate
// gradients at twice the products A v (r~ = r throughout), but it does not need A to be positive
// definite.
LinearStats bicg(std::size_t n, const VecOp& op, const VecOp& prec, const std::vector<double>& b, std::vector<double>& x,
                 double tol, int maxit, std::vector<double>* trace) {
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r = b, rt = b, z, zt, p, pt, q, qt;
    st.residual = norm(r);
    if (trace) trace->push_back(st.residual);
    if (st.residual <= tol) {
        st.converged = true;
        return st;
    }
    double rho_old = 0.0;
    for (int it = 1; it <= maxit; ++it) {
        if ((it & 63) == 0) check_cancelled();
        prec(r, z);
        prec(rt, zt);  // M^-T = M^-1
        const double rho = dot(z, rt);
        if (!(std::fabs(rho) > 0.0)) {
            st.breakdown = true;
            return st;
        }
        if (it == 1) {
            p = z;
            pt = zt;
        } else {
            const double beta = rho / rho_old;
            for (std::size_t i = 0; i < n; ++i) {
                p[i] = z[i] + beta * p[i];
                pt[i] = zt[i] + beta * pt[i];
            }
        }
        op(p, q);
        op(pt, qt);  // A^T = A
        const double ptq = dot(pt, q);
        if (!(std::fabs(ptq) > 0.0)) {
            st.breakdown = true;
            return st;
        }
        const double alpha = rho / ptq;
        axpy(alpha, p, x);
        axpy(-alpha, q, r);
        axpy(-alpha, qt, rt);
        rho_old = rho;
        st.iterations = it;
        st.residual = norm(r);
        if (trace) trace->push_back(st.residual);
        if (st.residual <= tol) {
            st.converged = true;
            return st;
        }
    }
    return st;
}

// ---------------------------------------------------------------------------------------------
// Symmetric tridiagonal eigenproblem: implicit QL (the algorithm of tqli in the original);
// eigenvalues ascending, optionally the eigenvectors (columns of z, z = identity on entry).

void tridiagonal_eigen(std::vector<double> d, std::vector<double> e, std::vector<double>& values, Matrix* z) {
    const int n = static_cast<int>(d.size());
    e.resize(static_cast<std::size_t>(n), 0.0);
    if (z) {
        *z = Matrix(static_cast<std::size_t>(n), static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) (*z)(static_cast<std::size_t>(i), static_cast<std::size_t>(i)) = 1.0;
    }
    // e[i] couples d[i] and d[i+1]
    for (int l = 0; l < n; ++l) {
        int iter = 0;
        int m;
        do {
            for (m = l; m < n - 1; ++m) {
                const double dd = std::fabs(d[m]) + std::fabs(d[m + 1]);
                if (std::fabs(e[m]) <= std::numeric_limits<double>::epsilon() * dd) break;
            }
            if (m != l) {
                if (++iter > 60) throw std::runtime_error("tridiagonal eigenproblem: no convergence");
                double g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                double r = std::hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0 ? std::fabs(r) : -std::fabs(r)));
                double s = 1.0, c = 1.0, p = 0.0;
                int i;
                for (i = m - 1; i >= l; --i) {
                    double f = s * e[i];
                    const double b = c * e[i];
                    r = std::hypot(f, g);
                    e[i + 1] = r;
                    if (r == 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }
                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;
                    if (z) {
                        for (int k = 0; k < n; ++k) {
                            const std::size_t kk = static_cast<std::size_t>(k);
                            f = (*z)(kk, static_cast<std::size_t>(i + 1));
                            (*z)(kk, static_cast<std::size_t>(i + 1)) = s * (*z)(kk, static_cast<std::size_t>(i)) + c * f;
                            (*z)(kk, static_cast<std::size_t>(i)) = c * (*z)(kk, static_cast<std::size_t>(i)) - s * f;
                        }
                    }
                }
                if (r == 0.0 && i >= l) continue;
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while (m != l);
    }
    // sort ascending
    std::vector<int> order(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort(order.begin(), order.end(), [&](int a, int b) { return d[static_cast<std::size_t>(a)] < d[static_cast<std::size_t>(b)]; });
    values.resize(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) values[static_cast<std::size_t>(i)] = d[static_cast<std::size_t>(order[static_cast<std::size_t>(i)])];
    if (z) {
        Matrix s(static_cast<std::size_t>(n), static_cast<std::size_t>(n));
        for (int k = 0; k < n; ++k)
            for (int i = 0; i < n; ++i)
                s(static_cast<std::size_t>(k), static_cast<std::size_t>(i)) = (*z)(static_cast<std::size_t>(k), static_cast<std::size_t>(order[static_cast<std::size_t>(i)]));
        *z = s;
    }
}

// ---------------------------------------------------------------------------------------------
// DVR (original: main, legmat, legmat2, jacbolbelso2)

struct LegendreDvr {
    std::vector<double> q;  // cos(theta) points, ascending
    Matrix t;               // FBR -> DVR: t(i, k) = component of P_{K+i}^K in DVR function k
    Matrix l2;              // j^2 in the DVR
    std::vector<double> w;  // 21.0: quadrature weights (in cos theta) of the DVR functions
    std::vector<double> s;  // 21.0: their signs at their own points (+1 or -1)
};

LegendreDvr legendre_dvr(int np, int K) {
    Matrix c(static_cast<std::size_t>(np), static_cast<std::size_t>(np), 0.0);
    for (int i = 0; i + 1 < np; ++i) {
        const double a = i + K;
        const double x = std::sqrt((a + 1.0 - K) * (a + 1.0 + K) / (4.0 * (a + 1.0) * (a + 1.0) - 1.0));
        c(static_cast<std::size_t>(i), static_cast<std::size_t>(i + 1)) = x;
        c(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(i)) = x;
    }
    EigenResult e = symmetric_eigen(c);
    LegendreDvr d;
    d.q = e.values;
    d.t = e.vectors;
    d.l2 = Matrix(static_cast<std::size_t>(np), static_cast<std::size_t>(np), 0.0);
    for (int k = 0; k < np; ++k)
        for (int l = 0; l < np; ++l) {
            double s = 0.0;
            for (int i = 0; i < np; ++i) {
                const double a = i + K;
                s += d.t(static_cast<std::size_t>(i), static_cast<std::size_t>(k)) * (a + 1.0) * a *
                     d.t(static_cast<std::size_t>(i), static_cast<std::size_t>(l));
            }
            d.l2(static_cast<std::size_t>(k), static_cast<std::size_t>(l)) = s;
        }
    // 21.0: t(0, k) = s_k w_k^1/2 phi_0(q_k), phi_0 the normalised P_K^K = (1 - x^2)^(K/2) / N_K^1/2,
    // N_K = 2^(2K+1) (K!)^2 / (2K+1)!, positive inside (-1, 1)
    const double log_norm = (2.0 * K + 1.0) * std::log(2.0) + 2.0 * std::lgamma(K + 1.0) - std::lgamma(2.0 * K + 2.0);
    d.w.resize(static_cast<std::size_t>(np));
    d.s.resize(static_cast<std::size_t>(np));
    for (int k = 0; k < np; ++k) {
        const double x = d.q[static_cast<std::size_t>(k)];
        const double phi0 = std::exp(0.5 * K * std::log(std::max(1e-300, 1.0 - x * x)) - 0.5 * log_norm);
        const double t0 = d.t(0, static_cast<std::size_t>(k));
        d.s[static_cast<std::size_t>(k)] = t0 >= 0.0 ? 1.0 : -1.0;
        d.w[static_cast<std::size_t>(k)] = (t0 / phi0) * (t0 / phi0);
    }
    return d;
}

// Coriolis block between K = k (rows) and K = k + 1 (columns), without the radial factor
// (original: legmat2): T_k^T e T_{k+1} sqrt(J(J+1) - k(k+1)), e the matrix of j_+ in the
// associated Legendre functions, with the factor sqrt(2) of the Wang functions for k = 0.
Matrix coriolis_block(const LegendreDvr& dk, const LegendreDvr& dk1, int np, int k, int J) {
    Matrix e(static_cast<std::size_t>(np), static_cast<std::size_t>(np), 0.0);
    const double coef = (k == 0) ? -std::sqrt(2.0) : -1.0;
    for (int i = 1; i < np; ++i) {
        const double a = i + 1 + k;
        e(static_cast<std::size_t>(i), static_cast<std::size_t>(i - 1)) = coef * std::sqrt((a - 1.0) * a - k * (k + 1.0));
    }
    Matrix tmp = matmul(dk.t.transposed(), e);
    Matrix h = matmul(tmp, dk1.t);
    const double f = std::sqrt(J * (J + 1.0) - k * (k + 1.0));
    for (std::size_t i = 0; i < h.rows(); ++i)
        for (std::size_t j = 0; j < h.cols(); ++j) h(i, j) *= f;
    return h;
}

// The coordinates and their bond lengths and angle (original: jacbolbelso2). The atoms in the
// notation of the thesis (Fig. 3.2, Table 3.1): A = O, B = H1 (the end of R1), C = H2 (the end
// of R2); g1 = r(A-P)/r(A-B), g2 = r(A-S)/r(A-C) of Sutcliffe and Tennyson.
//   Jacobi: R1 = r(O-H1), R2 from the centre of mass of O-H1 to H2: g1 = m(H1)/(m(O) + m(H1)), g2 = 0.
//   Radau: R1 = r(K-H1), R2 = r(K-H2), K = alpha O + (1 - alpha)(m(H1) H1 + m(H2) H2)/(m(H1) + m(H2)),
//     alpha = (m(O)/M)^1/2, beta = m(H1)/(m(H1) + m(H2)): g1 = 1 - alpha/(alpha + beta - alpha beta),
//     g2 = 1 - alpha/(1 - beta + alpha beta). The vectors of Sutcliffe and Tennyson are
//     H1 - S = (H1 - K)/k1 and H2 - P = (H2 - K)/k2 with k1 = 1 - beta + alpha beta,
//     k2 = alpha + beta - alpha beta: the same angle, distances scaled by constants. The distances
//     from K have the reduced masses m(H1) and m(H2) (and no kinetic coupling).
struct Coordinates {
    bool radau = false;
    double g1 = 0.0, g2 = 0.0, k1 = 1.0, k2 = 1.0;

    // r1 = r(O-H1), r2 = r(O-H2), cb = cos(H1-O-H2)
    void to_bond(double R1, double R2, double c, double& r1, double& r2, double& cb) const {
        if (!radau) {
            r1 = R1;
            r2 = std::sqrt(g1 * g1 * R1 * R1 + R2 * R2 + 2.0 * g1 * R1 * R2 * c);
            cb = R2 / r2 * c + g1 * R1 / r2;
            return;
        }
        const double a = R1 / k1, b = R2 / k2, den = 1.0 - g1 * g2;
        const double bb = (a * a + g2 * g2 * b * b + 2.0 * g2 * a * b * c) / (den * den);
        const double cc = (b * b + g1 * g1 * a * a + 2.0 * g1 * a * b * c) / (den * den);
        const double bc = (g1 * a * a + g2 * b * b + (1.0 + g1 * g2) * a * b * c) / (den * den);
        r1 = std::sqrt(bb);
        r2 = std::sqrt(cc);
        cb = bc / (r1 * r2);
    }
};

struct Problem {
    int ns = 0, np = 0, nb = 0;
    std::size_t nx = 0, n = 0;
    std::vector<int> K;
    std::vector<double> R;              // radial points (bohr)
    std::vector<double> wr;             // 21.0: their quadrature weights (bohr)
    std::vector<double> sr;             // 21.0: the signs of the radial DVR functions, (-1)^i
    Matrix t1, t2;                      // radial kinetic energy
    std::vector<double> r1, r2;         // 1 / (2 mu R^2)
    std::vector<LegendreDvr> leg;       // per block
    std::vector<Matrix> coriolis;       // between blocks b and b + 1
    Coordinates coords;                 // 21.0: Jacobi or Radau
    std::size_t index(int b, int i, int j, int l) const {
        return ((static_cast<std::size_t>(b) * static_cast<std::size_t>(ns) + static_cast<std::size_t>(i)) *
                    static_cast<std::size_t>(ns) + static_cast<std::size_t>(j)) * static_cast<std::size_t>(np) +
               static_cast<std::size_t>(l);
    }
};

std::vector<int> k_values(int J, int parity_block) {
    std::vector<int> k;
    if (J == 0) return {0};
    for (int K = (parity_block == 0 ? 0 : 1); K <= J; ++K) k.push_back(K);
    return k;
}

void validate(const NuclearMotionInput& in) {
    if (in.J < 0 || in.J > 30) throw std::invalid_argument("J must be from 0 to 30");
    if (in.parity_block != 0 && in.parity_block != 1) throw std::invalid_argument("the parity block must be 0 (K = 0 ... J) or 1 (K = 1 ... J)");
    if (in.J == 0 && in.parity_block == 1) throw std::invalid_argument("J = 0 has only the block K = 0");
    if (in.radial_points < 2 || in.radial_points > 200) throw std::invalid_argument("the radial DVR needs 2 to 200 points");
    if (in.angular_points < 2 || in.angular_points > 400) throw std::invalid_argument("the angular DVR needs 2 to 400 points");
    if (!(in.r_half_width > 0.0) || !(in.r_centre > in.r_half_width))
        throw std::invalid_argument("the radial grid must lie at positive distances: 0 < Rmax < Re");
    if (in.levels < 1 || in.levels > 500) throw std::invalid_argument("the number of levels must be from 1 to 500");
    if (in.max_iterations < 1) throw std::invalid_argument("the maximum number of iterations must be positive");
    if (!(in.eigenvalue_tolerance > 0.0) || !(in.residual_tolerance > 0.0) || !(in.linear_tolerance > 0.0))
        throw std::invalid_argument("the tolerances must be positive");
    const std::string& s = in.eigensolver;
    if (s != "lanczos" && s != "jacobi-davidson" && s != "davidson" && s != "inverse-power" && s != "rayleigh-quotient" &&
        s != "dense")
        throw std::invalid_argument(
            "eigensolver must be lanczos, jacobi-davidson, davidson, inverse-power, rayleigh-quotient or dense");
    if (in.coordinates != "jacobi" && in.coordinates != "radau")
        throw std::invalid_argument("the coordinates must be jacobi or radau");
    if (in.jd_inner_iterations < 0 || in.jd_inner_iterations > 1000 || in.davidson_inner_iterations < 0 ||
        in.davidson_inner_iterations > 1000)
        throw std::invalid_argument("the inner iterations must be from 0 to 1000");
    if (in.max_linear_iterations < 1) throw std::invalid_argument("the maximum number of linear-solver iterations must be positive");
    parse_linear(in.linear_solver);
    parse_precond(in.preconditioner);
}

int thread_count(int requested) {
    if (requested > 0) return requested;
    const unsigned h = std::thread::hardware_concurrency();
    return h == 0 ? 1 : static_cast<int>(h);
}

Problem build_problem(const NuclearMotionInput& in, const IsotopologueData& iso, NuclearMotionResult& res) {
    Problem p;
    p.ns = in.radial_points;
    p.np = in.angular_points;
    p.K = k_values(in.J, in.parity_block);
    p.nb = static_cast<int>(p.K.size());
    p.nx = static_cast<std::size_t>(p.ns) * static_cast<std::size_t>(p.ns) * static_cast<std::size_t>(p.np);
    p.n = p.nx * static_cast<std::size_t>(p.nb);

    // masses: xmass = {H, H, O} as in the original: atom C = xmass[0] (H2), B = xmass[1] (H1),
    // A = xmass[2] (O)
    const double m0 = iso.hydrogen, m1 = iso.hydrogen, m2 = iso.oxygen;
    double tred1, tred2;
    Coordinates& co = p.coords;
    co.radau = in.coordinates == "radau";
    if (!co.radau) {
        // R1 joins atoms 2 and 3 (O-H), R2 joins atom 1 to their centre of mass (original: nkoord = 1)
        co.g1 = m1 / (m1 + m2);
        co.g2 = 0.0;
        tred1 = kAmuToMe / (1.0 / m1 + 1.0 / m2);
        tred2 = kAmuToMe / (1.0 / m0 + co.g1 * co.g1 / m1 + (1.0 - co.g1) * (1.0 - co.g1) / m2);
    } else {
        // original: nkoord = 2 (a, b, g1, g2); the distances from the canonical point K
        const double alpha = std::sqrt(m2 / (m0 + m1 + m2)), beta = m1 / (m0 + m1);
        co.g1 = 1.0 - alpha / (alpha + beta - alpha * beta);
        co.g2 = 1.0 - alpha / (1.0 - beta + alpha * beta);
        co.k1 = 1.0 - beta + alpha * beta;
        co.k2 = alpha + beta - alpha * beta;
        tred1 = kAmuToMe * m1;
        tred2 = kAmuToMe * m0;
    }
    res.g1 = co.g1;
    res.g2 = co.g2;
    res.masses[0] = m0;
    res.masses[1] = m1;
    res.masses[2] = m2;
    res.mu1 = tred1;
    res.mu2 = tred2;

    // Hermite DVR: the roots of H_ns (eigenvalues of the coordinate matrix)
    const int ns = p.ns;
    Matrix hx(static_cast<std::size_t>(ns), static_cast<std::size_t>(ns), 0.0);
    for (int i = 0; i + 1 < ns; ++i) {
        const double x = std::sqrt((i + 1.0) / 2.0);
        hx(static_cast<std::size_t>(i), static_cast<std::size_t>(i + 1)) = x;
        hx(static_cast<std::size_t>(i + 1), static_cast<std::size_t>(i)) = x;
    }
    const EigenResult hx_eigen = symmetric_eigen(hx);
    const std::vector<double>& x = hx_eigen.values;
    const double allk = x.back() / in.r_half_width;
    p.R.resize(static_cast<std::size_t>(ns));
    for (int i = 0; i < ns; ++i) p.R[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)] / allk + in.r_centre;
    // 21.0: weights and signs of the radial DVR functions. T(0, i) = w_i^1/2 h_0(x_i) for the
    // functions positive at their own points (h_0 = pi^-1/4 exp(-x^2/2)); the kinetic-energy
    // formula below (original: hr1) belongs to the functions of sign (-1)^i.
    p.wr.resize(static_cast<std::size_t>(ns));
    p.sr.resize(static_cast<std::size_t>(ns));
    for (int i = 0; i < ns; ++i) {
        const double xi = x[static_cast<std::size_t>(i)];
        const double h0 = std::pow(kPi, -0.25) * std::exp(-0.5 * xi * xi);
        const double t0 = hx_eigen.vectors(0, static_cast<std::size_t>(i)) / h0;
        p.wr[static_cast<std::size_t>(i)] = t0 * t0 / allk;
        p.sr[static_cast<std::size_t>(i)] = (i % 2 == 0) ? 1.0 : -1.0;
    }
    auto kinetic = [&](double tred) {
        Matrix t(static_cast<std::size_t>(ns), static_cast<std::size_t>(ns), 0.0);
        const double coef = allk * allk / (2.0 * tred);
        for (int j = 0; j < ns; ++j)
            for (int k = 0; k < ns; ++k) {
                const double xj = x[static_cast<std::size_t>(j)], xk = x[static_cast<std::size_t>(k)];
                t(static_cast<std::size_t>(j), static_cast<std::size_t>(k)) =
                    j == k ? coef * (2.0 * (ns - 1.0) / 3.0 + 0.5 - xj * xj / 3.0)
                           : (-0.5 + 2.0 / ((xj - xk) * (xj - xk))) * coef;
            }
        return t;
    };
    p.t1 = kinetic(tred1);
    p.t2 = kinetic(tred2);
    p.r1.resize(static_cast<std::size_t>(ns));
    p.r2.resize(static_cast<std::size_t>(ns));
    for (int i = 0; i < ns; ++i) {
        const double r = p.R[static_cast<std::size_t>(i)];
        p.r1[static_cast<std::size_t>(i)] = (1.0 / (2.0 * tred1)) / (r * r);
        p.r2[static_cast<std::size_t>(i)] = (1.0 / (2.0 * tred2)) / (r * r);
    }
    for (int K : p.K) p.leg.push_back(legendre_dvr(p.np, K));
    for (int b = 0; b + 1 < p.nb; ++b)
        p.coriolis.push_back(coriolis_block(p.leg[static_cast<std::size_t>(b)], p.leg[static_cast<std::size_t>(b + 1)], p.np, p.K[static_cast<std::size_t>(b)], in.J));

    res.r_grid = p.R;
    res.r_weights = p.wr;
    res.k_values = p.K;
    for (const LegendreDvr& d : p.leg) res.cos_grid.push_back(d.q);
    return p;
}

Csr build_hamiltonian(const Problem& p, const NuclearMotionInput& in, const IsotopologueData& iso, int threads) {
    const int ns = p.ns, np = p.np, nb = p.nb;
    // the potential at every DVR point (measured from the minimum of the surface)
    std::vector<double> V(p.n);
    std::atomic<long long> done{0};
    const long long blocks = static_cast<long long>(nb) * ns;
    run_parallel(threads, [&](int t) {
        for (long long bi = t; bi < blocks; bi += threads) {
            check_cancelled();
            const int b = static_cast<int>(bi / ns), i = static_cast<int>(bi % ns);
            for (int j = 0; j < ns; ++j)
                for (int l = 0; l < np; ++l) {
                    double r1, r2, cb;
                    p.coords.to_bond(p.R[static_cast<std::size_t>(i)], p.R[static_cast<std::size_t>(j)],
                                     p.leg[static_cast<std::size_t>(b)].q[static_cast<std::size_t>(l)], r1, r2, cb);
                    cb = std::max(-1.0, std::min(1.0, cb));
                    V[p.index(b, i, j, l)] = cvrqd_potential(in.isotopologue, r1, r2, cb) - iso.vmin;
                }
            const long long d = ++done;
            if (t == 0) report_progress_text(fmt("Potential energy at the DVR points: %lld of %lld", d, blocks));
        }
    });

    Csr a;
    a.n = p.n;
    a.threads = threads;
    const std::size_t per_row = static_cast<std::size_t>(2 * ns + np - 2 + (nb > 1 ? 2 * np : 0));
    a.row.reserve(p.n + 1);
    a.col.reserve(p.n * per_row);
    a.val.reserve(p.n * per_row);
    a.diag_pos.resize(p.n);
    a.diag.resize(p.n);
    a.row.push_back(0);
    const double JJ = in.J * (in.J + 1.0);
    for (int b = 0; b < nb; ++b) {
        check_cancelled();
        const int K = p.K[static_cast<std::size_t>(b)];
        const Matrix& l2 = p.leg[static_cast<std::size_t>(b)].l2;
        for (int i = 0; i < ns; ++i)
            for (int j = 0; j < ns; ++j)
                for (int l = 0; l < np; ++l) {
                    const std::size_t rowi = p.index(b, i, j, l);
                    const double ri = p.r1[static_cast<std::size_t>(i)];
                    const double rsum = ri + p.r2[static_cast<std::size_t>(j)];
                    auto push = [&](std::size_t c, double v) {
                        a.col.push_back(static_cast<int>(c));
                        a.val.push_back(v);
                    };
                    if (b > 0) {  // K - 1: the transpose of the block (K - 1, K)
                        const Matrix& c = p.coriolis[static_cast<std::size_t>(b - 1)];
                        for (int lv = 0; lv < np; ++lv)
                            push(p.index(b - 1, i, j, lv), c(static_cast<std::size_t>(lv), static_cast<std::size_t>(l)) * ri);
                    }
                    for (int iv = 0; iv < i; ++iv)
                        push(p.index(b, iv, j, l), p.t1(static_cast<std::size_t>(i), static_cast<std::size_t>(iv)));
                    for (int jv = 0; jv < j; ++jv)
                        push(p.index(b, i, jv, l), p.t2(static_cast<std::size_t>(j), static_cast<std::size_t>(jv)));
                    for (int lv = 0; lv < np; ++lv) {
                        double v = l2(static_cast<std::size_t>(l), static_cast<std::size_t>(lv)) * rsum;
                        if (lv == l) {
                            v += p.t1(static_cast<std::size_t>(i), static_cast<std::size_t>(i)) +
                                 p.t2(static_cast<std::size_t>(j), static_cast<std::size_t>(j)) + V[rowi];
                            if (in.J > 0) v += (JJ - 2.0 * K * K) * ri;
                            a.diag_pos[rowi] = a.col.size();
                            a.diag[rowi] = v;
                        }
                        push(p.index(b, i, j, lv), v);
                    }
                    for (int jv = j + 1; jv < ns; ++jv)
                        push(p.index(b, i, jv, l), p.t2(static_cast<std::size_t>(j), static_cast<std::size_t>(jv)));
                    for (int iv = i + 1; iv < ns; ++iv)
                        push(p.index(b, iv, j, l), p.t1(static_cast<std::size_t>(i), static_cast<std::size_t>(iv)));
                    if (b + 1 < nb) {  // K + 1
                        const Matrix& c = p.coriolis[static_cast<std::size_t>(b)];
                        for (int lv = 0; lv < np; ++lv)
                            push(p.index(b + 1, i, j, lv), c(static_cast<std::size_t>(l), static_cast<std::size_t>(lv)) * ri);
                    }
                    a.row.push_back(a.col.size());
                }
    }
    for (double d : a.diag)
        if (!(d > 0.0)) throw std::runtime_error("the Hamiltonian has a diagonal element that is not positive");
    return a;
}

// The reproducible start vector of the iterative solvers: pseudo-random numbers (or, 21.0, the
// Gaussian of Algorithm 5, see gaussian_start_vector).
std::vector<double> start_vector(std::size_t n, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    std::vector<double> v(n);
    for (double& x : v) x = u(rng);
    scale(v, 1.0 / norm(v));
    return v;
}

// Gram-Schmidt of w against the columns `basis` (twice, for orthogonality to machine precision).
void orthogonalize(std::vector<double>& w, const std::vector<std::vector<double>>& basis) {
    for (int pass = 0; pass < 2; ++pass)
        for (const auto& v : basis) axpy(-dot(w, v), v, w);
}

struct Eigenpairs {
    std::vector<double> values;
    std::vector<std::vector<double>> vectors;  // empty: eigenvalues only (Lanczos)
    std::vector<double> residuals;
};

void emit(const NuclearOutputSink& sink, std::ostringstream& os, const std::string& text) {
    os << text;
    if (sink) sink(text);
}

// The state of an iterative eigensolver: the matrix, the input, the result being filled, the
// linear solver and its preconditioner.
struct Solver {
    const Csr& a;
    const NuclearMotionInput& in;
    NuclearMotionResult& res;
    const NuclearOutputSink& sink;
    std::ostringstream& os;
    Precond pc;
    Linear lin;
    std::vector<double> guess;  // 21.0: the Gaussian start vector (empty: none)
    std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    Preconditioner unshifted;   // M ~ H
    bool have_interval = false; // Chebyshev: the interval of M^-1 H
    double lo = 0.0, hi = 0.0;
    bool noted_indefinite = false;

    Solver(const Csr& a_, const NuclearMotionInput& in_, NuclearMotionResult& res_, const NuclearOutputSink& sink_,
           std::ostringstream& os_)
        : a(a_), in(in_), res(res_), sink(sink_), os(os_), pc(parse_precond(in_.preconditioner)),
          lin(parse_linear(in_.linear_solver)) {
        unshifted.setup(a, pc, 0.0, true);
    }

    std::size_t n() const { return a.n; }

    void H(const std::vector<double>& x, std::vector<double>& y) {
        a.multiply(x, y);
        ++res.matrix_vector_products;
    }

    double elapsed() const { return seconds_since(t0); }

    void record(int step, std::vector<double>&& trace) {
        if (res.first_solve.step == 0) res.first_solve = LinearSolveTrace{step, trace};
        res.last_solve = LinearSolveTrace{step, std::move(trace)};
    }

    // The interval of the eigenvalues of M^-1 A for the Chebyshev iteration: Lanczos in the
    // M inner product (`steps` steps from a fixed random vector); [smallest Ritz value, 1.05 x
    // the largest]. The smallest Ritz value lies above the smallest eigenvalue, which slows the
    // iteration for those components but keeps it stable (their factor stays below 1). For an
    // indefinite A (a shifted matrix) the iteration only filters: [0.01 hi, hi].
    void interval(const VecOp& op, const VecOp& prec, int steps, double& lo_out, double& hi_out) {
        const std::size_t nn = n();
        steps = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(steps), nn));
        std::vector<double> r = start_vector(nn, 7u), z, v, vprev(nn, 0.0), w, y;
        prec(r, z);
        double beta = std::sqrt(std::max(0.0, dot(r, z)));
        std::vector<double> alpha_t, beta_t;
        for (int j = 0; j < steps && beta > 1e-300; ++j) {
            v = r;
            scale(v, 1.0 / beta);
            w = z;
            scale(w, 1.0 / beta);
            op(w, y);
            if (j > 0) axpy(-beta, vprev, y);
            const double al = dot(w, y);
            axpy(-al, v, y);
            alpha_t.push_back(al);
            vprev = v;
            r = y;
            prec(r, z);
            beta = std::sqrt(std::max(0.0, dot(r, z)));
            if (j + 1 < steps) beta_t.push_back(beta);
        }
        beta_t.resize(alpha_t.empty() ? 0 : alpha_t.size() - 1);
        std::vector<double> theta;
        tridiagonal_eigen(alpha_t, beta_t, theta, nullptr);
        hi_out = 1.05 * theta.back();
        lo_out = theta.front();
        if (!(lo_out > 0.0) || lo_out > 0.5 * hi_out) {
            if (!(lo_out > 0.0) && !noted_indefinite) {
                noted_indefinite = true;
                res.notes.push_back("Chebyshev iteration on an indefinite shifted matrix: it acts as a polynomial filter on [0.01 hi, hi]");
            }
            lo_out = 0.01 * hi_out;
        }
    }

    // One linear solve of A x = b with the chosen method. strict: A is H (positive definite).
    LinearStats solve(const VecOp& op, const VecOp& prec, const std::vector<double>& b, std::vector<double>& x, double tol,
                      int maxit, bool strict, bool cached_interval, int step) {
        std::vector<double> trace;
        LinearStats st;
        switch (lin) {
        case Linear::Pcg:
            st = pcg(n(), op, prec, b, x, tol, maxit, strict, &trace);
            break;
        case Linear::Minres:
            st = minres(n(), op, prec, b, x, tol, maxit, &trace);
            break;
        case Linear::Bicg:
            st = bicg(n(), op, prec, b, x, tol, maxit, &trace);
            break;
        case Linear::Chebyshev: {
            double l = 0.0, h = 0.0;
            if (cached_interval) {
                if (!have_interval) {
                    interval(op, prec, 40, lo, hi);
                    have_interval = true;
                    res.chebyshev_min = lo;
                    res.chebyshev_max = hi;
                }
                l = lo;
                h = hi;
            } else {
                interval(op, prec, 10, l, h);
                if (res.chebyshev_max == 0.0) {
                    res.chebyshev_min = l;
                    res.chebyshev_max = h;
                }
            }
            st = chebyshev(n(), op, prec, b, x, tol, maxit, l, h, &trace);
            break;
        }
        }
        res.linear_iterations += st.iterations;
        record(step, std::move(trace));
        return st;
    }

    // H x = b (Lanczos, inverse iteration)
    LinearStats solve_h(const std::vector<double>& b, std::vector<double>& x, int step) {
        VecOp op = [this](const std::vector<double>& v, std::vector<double>& y) { H(v, y); };
        VecOp prec = [this](const std::vector<double>& r, std::vector<double>& z) { unshifted.apply(r, z); };
        return solve(op, prec, b, x, in.linear_tolerance, in.max_linear_iterations, true, true, step);
    }

    // (I - QQ^T)(H - theta)(I - QQ^T) x = b on the complement of the orthonormal columns Q
    // (Jacobi-Davidson, Rayleigh-quotient iteration; Q empty: Davidson); the preconditioner
    // (I - QQ^T) M^-1 (I - QQ^T), M ~ H - theta (absolute diagonal: positive definite).
    LinearStats solve_shifted(double theta, const std::vector<const std::vector<double>*>& Q, const std::vector<double>& b,
                              std::vector<double>& x, double tol, int maxit, int step) {
        Preconditioner m;
        m.setup(a, pc, theta, true);
        auto project = [&Q](std::vector<double>& v) {
            for (const std::vector<double>* q : Q) axpy(-dot(v, *q), *q, v);
        };
        VecOp op = [&, theta](const std::vector<double>& v, std::vector<double>& y) {
            std::vector<double> pv = v;
            project(pv);
            H(pv, y);
            axpy(-theta, pv, y);
            project(y);
        };
        VecOp prec = [&](const std::vector<double>& r, std::vector<double>& z) {
            if (Q.empty()) {
                m.apply(r, z);
                return;
            }
            std::vector<double> pr = r;
            project(pr);
            m.apply(pr, z);
            project(z);
        };
        LinearStats st = solve(op, prec, b, x, tol, maxit, false, false, step);
        project(x);
        return st;
    }

    void step_line(const NuclearSolverStep& s, int third) {
        emit(sink, os, fmt("  %5d  %18.10f  %12.3e  %8d  %5d\n", s.iteration, s.lowest * kHartreeToCm, s.change, third, s.converged));
    }
};

std::vector<double> residuals_of(Solver& sv, const std::vector<double>& values, const std::vector<std::vector<double>>& vectors) {
    std::vector<double> res(values.size());
    std::vector<double> hy;
    for (std::size_t k = 0; k < values.size(); ++k) {
        sv.H(vectors[k], hy);
        axpy(-values[k], vectors[k], hy);
        res[k] = norm(hy);
    }
    return res;
}

// The start vector: the Gaussian (21.0, with a little of the random vector when more than one level
// is wanted: the Gaussian alone has hardly any component of the antisymmetric levels) or random.
std::vector<double> first_vector(const Solver& sv, int want, double admixture) {
    if (sv.guess.empty()) return start_vector(sv.n(), 20u);
    std::vector<double> v = sv.guess;
    if (want > 1 && admixture > 0.0) {
        axpy(admixture, start_vector(sv.n(), 20u), v);
        scale(v, 1.0 / norm(v));
    }
    return v;
}

// Shift-invert Lanczos (original: lanczos): the Krylov space of H^-1, every step a linear solve,
// full reorthogonalisation; the eigenvalues E = 1/theta of the largest Ritz values theta. With the
// eigenvalues only (the original's bsvektor = 0), converged when the wanted eigenvalues change
// by less than the tolerance; otherwise also every Ritz vector's residual below its tolerance.
Eigenpairs lanczos(Solver& sv) {
    const NuclearMotionInput& in = sv.in;
    NuclearMotionResult& res = sv.res;
    const std::size_t n = sv.n();
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    const int mmax = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.max_iterations), n));
    std::vector<std::vector<double>> V;
    V.reserve(static_cast<std::size_t>(mmax) + 1);
    V.push_back(first_vector(sv, want, 1e-2));
    std::vector<double> alpha, beta, prev, w;
    Eigenpairs out;
    for (int j = 0; j < mmax; ++j) {
        check_cancelled();
        const LinearStats ls = sv.solve_h(V[static_cast<std::size_t>(j)], w, j + 1);
        if (!ls.converged)
            res.notes.push_back(fmt("Lanczos step %d: the linear solve stopped after %d iterations (residual %.2e)", j + 1,
                                    ls.iterations, ls.residual));
        const double al = dot(w, V[static_cast<std::size_t>(j)]);
        alpha.push_back(al);
        axpy(-al, V[static_cast<std::size_t>(j)], w);
        if (j > 0) axpy(-beta.back(), V[static_cast<std::size_t>(j - 1)], w);
        orthogonalize(w, V);
        const double be = norm(w);
        res.iterations = j + 1;

        const int m = j + 1;
        NuclearSolverStep step;
        step.iteration = m;
        step.linear_iterations = ls.iterations;
        step.linear_residual = ls.residual;
        bool done = false;
        if (m >= want || be < 1e-14) {
            std::vector<double> theta;
            std::vector<double> e(beta);
            tridiagonal_eigen(alpha, e, theta, nullptr);
            std::vector<double> E;
            for (int k = 0; k < std::min(want, m); ++k) E.push_back(1.0 / theta[static_cast<std::size_t>(m - 1 - k)]);
            double change = std::numeric_limits<double>::infinity();
            if (prev.size() == E.size()) {
                change = 0.0;
                for (std::size_t k = 0; k < E.size(); ++k) change = std::max(change, std::fabs(E[k] - prev[k]));
            }
            prev = E;
            step.lowest = E.front();
            step.change = change;
            const bool values_converged = change < in.eigenvalue_tolerance && static_cast<int>(E.size()) == want;
            if (!in.eigenvectors) {
                int conv = 0;
                for (std::size_t k = 0; k < E.size(); ++k)
                    if (std::isfinite(change) && change < in.eigenvalue_tolerance) ++conv;
                step.converged = conv;
                out.values = E;
                out.vectors.clear();
                out.residuals.assign(E.size(), std::numeric_limits<double>::quiet_NaN());
                if (values_converged || be < 1e-14 || m == mmax) done = true;
                res.converged = values_converged;
            } else if (values_converged || be < 1e-14 || m == mmax) {
                // Ritz vectors and their residuals
                Matrix s;
                std::vector<double> th;
                std::vector<double> e2(beta);
                tridiagonal_eigen(alpha, e2, th, &s);
                out.values.clear();
                out.vectors.clear();
                for (int k = 0; k < std::min(want, m); ++k) {
                    const std::size_t col = static_cast<std::size_t>(m - 1 - k);
                    std::vector<double> y(n, 0.0);
                    for (int q = 0; q < m; ++q) axpy(s(static_cast<std::size_t>(q), col), V[static_cast<std::size_t>(q)], y);
                    scale(y, 1.0 / norm(y));
                    out.values.push_back(1.0 / th[col]);
                    out.vectors.push_back(std::move(y));
                }
                out.residuals = residuals_of(sv, out.values, out.vectors);
                int conv = 0;
                for (double r : out.residuals)
                    if (r < in.residual_tolerance) ++conv;
                step.converged = conv;
                if (conv == static_cast<int>(out.values.size()) && static_cast<int>(out.values.size()) == want) done = true;
                if (be < 1e-14 || m == mmax) done = true;
                res.converged = conv == want && static_cast<int>(out.values.size()) == want;
            }
            report_progress(JobStage::NuclearEigen, m, E.front(), std::isfinite(change) ? change : 0.0);
            report_progress_text(fmt("Lanczos step %d · E1 = %.6f cm-1 · change %.1e E_h · %d linear-solver iterations", m,
                                     E.front() * kHartreeToCm, std::isfinite(change) ? change : 0.0, ls.iterations));
        } else {
            report_progress(JobStage::NuclearEigen, m, 0.0, 0.0);
            report_progress_text(fmt("Lanczos step %d of at least %d · %d linear-solver iterations", m, want, ls.iterations));
        }
        step.seconds = sv.elapsed();
        res.history.push_back(step);
        sv.step_line(step, ls.iterations);
        if (done) break;
        beta.push_back(be);
        scale(w, 1.0 / be);
        V.push_back(w);
    }
    return out;
}

// Davidson and Jacobi-Davidson for the lowest eigenpairs (original: davidson, jacobidavidson).
//   Davidson: the correction t = -M^-1 r (M ~ H - theta: the preconditioner; "none": t = -r),
//     or (21.0, davidson_inner_iterations > 0) a few steps of the linear solver on
//     (H - theta) t = -r with that preconditioner (an exact solve would give t = -u, nothing new).
//   Jacobi-Davidson: t = -(M^-1 r - mu M^-1 u), mu such that t is orthogonal to u (original:
//     jdcorrectionequation_tuti), refined by jd_inner_iterations steps of the linear solver on the
//     projected correction equation (I - QQ^T)(H - sigma)(I - QQ^T) d = -r - (...) t, Q = the
//     converged eigenvectors and u, sigma = theta near convergence (see below).
Eigenpairs davidson(Solver& sv, bool jacobi_davidson) {
    const NuclearMotionInput& in = sv.in;
    NuclearMotionResult& res = sv.res;
    const Csr& a = sv.a;
    const std::size_t n = sv.n();
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    const int maxdim = static_cast<int>(std::min<std::size_t>(
        static_cast<std::size_t>(in.subspace_size > 0 ? in.subspace_size : std::max(2 * want + 20, 40)), n));
    const int keep = std::min(maxdim - 1, std::max(want + 5, want * 3 / 2));
    std::vector<std::vector<double>> V, AV;
    // start: the reproducible (or Gaussian) vector and unit vectors at the smallest diagonal elements
    {
        std::vector<std::size_t> order(n);
        for (std::size_t i = 0; i < n; ++i) order[i] = i;
        const std::size_t m0 = std::min<std::size_t>(static_cast<std::size_t>(want), n);
        std::partial_sort(order.begin(), order.begin() + static_cast<long>(m0), order.end(),
                          [&](std::size_t x, std::size_t y) { return a.diag[x] < a.diag[y]; });
        std::vector<std::vector<double>> init;
        init.push_back(sv.guess.empty() ? start_vector(n, 20u) : sv.guess);
        for (std::size_t k = 0; k < m0 && init.size() < static_cast<std::size_t>(want); ++k) {
            std::vector<double> e(n, 0.0);
            e[order[k]] = 1.0;
            init.push_back(e);
        }
        for (auto& v : init) {
            orthogonalize(v, V);
            const double nv = norm(v);
            if (nv < 1e-10) continue;
            scale(v, 1.0 / nv);
            V.push_back(v);
            std::vector<double> hv;
            sv.H(v, hv);
            AV.push_back(hv);
        }
    }
    Eigenpairs out;
    std::vector<double> u, au, r, t;
    for (int it = 1; it <= in.max_iterations; ++it) {
        check_cancelled();
        const int k = static_cast<int>(V.size());
        Matrix G(static_cast<std::size_t>(k), static_cast<std::size_t>(k));
        for (int i = 0; i < k; ++i)
            for (int j = i; j < k; ++j) {
                const double s = dot(V[static_cast<std::size_t>(i)], AV[static_cast<std::size_t>(j)]);
                G(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = s;
                G(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = s;
            }
        const EigenResult ev = symmetric_eigen(G);
        auto ritz = [&](int idx, std::vector<double>& x, std::vector<double>& ax) {
            x.assign(n, 0.0);
            ax.assign(n, 0.0);
            for (int q = 0; q < k; ++q) {
                const double c = ev.vectors(static_cast<std::size_t>(q), static_cast<std::size_t>(idx));
                axpy(c, V[static_cast<std::size_t>(q)], x);
                axpy(c, AV[static_cast<std::size_t>(q)], ax);
            }
        };
        // the leading converged Ritz pairs and the first unconverged one
        int conv = 0;
        double rnorm = 0.0, theta = 0.0;
        std::vector<std::vector<double>> locked;
        out.values.clear();
        out.vectors.clear();
        out.residuals.clear();
        for (int idx = 0; idx < std::min(want, k); ++idx) {
            ritz(idx, u, au);
            theta = ev.values[static_cast<std::size_t>(idx)];
            r = au;
            axpy(-theta, u, r);
            rnorm = norm(r);
            out.values.push_back(theta);
            out.vectors.push_back(u);
            out.residuals.push_back(rnorm);
            if (rnorm < in.residual_tolerance) {
                ++conv;
                locked.push_back(u);
                continue;
            }
            break;
        }
        res.iterations = it;
        NuclearSolverStep step;
        step.iteration = it;
        step.lowest = ev.values.front();
        step.change = conv < want ? rnorm : 0.0;
        step.converged = conv;
        report_progress(JobStage::NuclearEigen, it, ev.values.front(), step.change);
        report_progress_text(fmt("%s iteration %d · %d of %d levels converged · residual %.1e E_h",
                                 jacobi_davidson ? "Jacobi-Davidson" : "Davidson", it, conv, want, step.change));
        auto finish_step = [&](const LinearStats* ls) {
            if (ls) {
                step.linear_iterations = ls->iterations;
                step.linear_residual = ls->residual;
            }
            step.seconds = sv.elapsed();
            res.history.push_back(step);
            sv.step_line(step, k);
        };
        if (conv >= want) {
            res.converged = true;
            finish_step(nullptr);
            break;
        }
        if (k >= static_cast<int>(n)) {
            finish_step(nullptr);
            break;
        }
        // correction vector for the target pair (theta, u, r); the direct preconditioner with the
        // signed diagonal D - theta, as the original's transform2 and ILUprecondfordavidson
        Preconditioner m;
        m.setup(a, sv.pc, theta, false);
        t.assign(n, 0.0);
        LinearStats inner;
        bool have_inner = false;
        // the shift of the inner solves: the target 0 (below the spectrum of H, so that they steer
        // towards the lowest remaining level) while the residual is large, the Ritz value theta
        // once ||r|| < 0.01 theta (Fokkema, Sleijpen and van der Vorst, SIAM J. Sci. Comput. 20,
        // 94 (1998)); from a random start with theta far up in the spectrum, accurate solves with
        // theta converge to an excited level instead
        const double sigma = rnorm < 1e-2 * std::fabs(theta) ? theta : 0.0;
        if (!jacobi_davidson) {
            if (in.davidson_inner_iterations > 0) {
                std::vector<double> rhs = r;
                scale(rhs, -1.0);
                inner = sv.solve_shifted(sigma, {}, rhs, t, in.linear_tolerance, in.davidson_inner_iterations, it);
                have_inner = true;
                if (!(norm(t) > 0.0)) {
                    m.apply(r, t);
                    scale(t, -1.0);
                }
            } else {
                m.apply(r, t);
                scale(t, -1.0);
            }
        } else {
            std::vector<double> y, rh;
            m.apply(u, y);
            m.apply(r, rh);
            const double uy = dot(u, y);
            const double mu = std::fabs(uy) > 0.0 ? dot(u, rh) / uy : 0.0;
            for (std::size_t i = 0; i < n; ++i) t[i] = -(rh[i] - mu * y[i]);
            if (in.jd_inner_iterations > 0) {
                std::vector<const std::vector<double>*> Q;
                for (const auto& q : locked) Q.push_back(&q);
                Q.push_back(&u);
                auto project = [&Q](std::vector<double>& x) {
                    for (const std::vector<double>* q : Q) axpy(-dot(x, *q), *q, x);
                };
                project(t);
                // rhs = -r - A t with A = (I - QQ^T)(H - sigma)(I - QQ^T)
                std::vector<double> rhs = r, at;
                scale(rhs, -1.0);
                sv.H(t, at);
                axpy(-sigma, t, at);
                project(at);
                axpy(-1.0, at, rhs);
                project(rhs);
                std::vector<double> d;
                inner = sv.solve_shifted(sigma, Q, rhs, d, in.linear_tolerance, in.jd_inner_iterations, it);
                have_inner = true;
                axpy(1.0, d, t);
                project(t);
            }
        }
        finish_step(have_inner ? &inner : nullptr);
        // restart: keep the lowest Ritz vectors
        if (k + 1 > maxdim) {
            std::vector<std::vector<double>> V2, AV2;
            for (int idx = 0; idx < keep; ++idx) {
                std::vector<double> x, ax;
                ritz(idx, x, ax);
                V2.push_back(std::move(x));
                AV2.push_back(std::move(ax));
            }
            V.swap(V2);
            AV.swap(AV2);
        }
        orthogonalize(t, V);
        double nt = norm(t);
        if (!(nt >= 1e-12)) {
            // stagnation: add the residual instead
            t = r;
            orthogonalize(t, V);
            nt = norm(t);
            if (nt < 1e-14) break;
        }
        scale(t, 1.0 / nt);
        std::vector<double> ht;
        sv.H(t, ht);
        V.push_back(t);
        AV.push_back(ht);
    }
    // the residuals of the final vectors (all wanted pairs)
    if (static_cast<int>(out.values.size()) < want) {
        const int k = static_cast<int>(V.size());
        Matrix G(static_cast<std::size_t>(k), static_cast<std::size_t>(k));
        for (int i = 0; i < k; ++i)
            for (int j = 0; j < k; ++j) G(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = dot(V[static_cast<std::size_t>(i)], AV[static_cast<std::size_t>(j)]);
        const EigenResult ev = symmetric_eigen(G);
        out.values.clear();
        out.vectors.clear();
        for (int idx = 0; idx < std::min(want, k); ++idx) {
            std::vector<double> x(n, 0.0);
            for (int q = 0; q < k; ++q) axpy(ev.vectors(static_cast<std::size_t>(q), static_cast<std::size_t>(idx)), V[static_cast<std::size_t>(q)], x);
            out.values.push_back(ev.values[static_cast<std::size_t>(idx)]);
            out.vectors.push_back(std::move(x));
        }
        out.residuals = residuals_of(sv, out.values, out.vectors);
    }
    return out;
}

// Orthonormalise the columns in order (modified Gram-Schmidt, twice); a column that vanishes is
// replaced by a random vector.
void orthonormalize(std::vector<std::vector<double>>& X, unsigned seed) {
    for (std::size_t k = 0; k < X.size(); ++k) {
        for (int pass = 0; pass < 2; ++pass)
            for (std::size_t q = 0; q < k; ++q) axpy(-dot(X[k], X[q]), X[q], X[k]);
        double nk = norm(X[k]);
        if (!(nk > 1e-12)) {
            X[k] = start_vector(X[k].size(), seed + static_cast<unsigned>(k));
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t q = 0; q < k; ++q) axpy(-dot(X[k], X[q]), X[q], X[k]);
            nk = norm(X[k]);
        }
        scale(X[k], 1.0 / nk);
    }
}

// 21.0: inverse iteration (original: invpower): x <- H^-1 x / ||H^-1 x||, the Rayleigh quotient
// theta = x^T H x, until ||H x - theta x|| is below the tolerance; converges to the lowest level
// with the factor E1/E2 per step. For several levels, inverse subspace iteration: a block of
// p = levels + guard vectors (guard: max(2, levels/2)), each solved with H, orthonormalised and
// rotated to the Ritz vectors of H in their span (Rayleigh-Ritz); level k converges with the
// factor E_k/E_(p+1). Converged leading vectors are no longer solved for.
Eigenpairs inverse_power(Solver& sv) {
    const NuclearMotionInput& in = sv.in;
    NuclearMotionResult& res = sv.res;
    const std::size_t n = sv.n();
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    const int p = want == 1 ? 1 : static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(want + std::max(2, want / 2)), n));
    std::vector<std::vector<double>> X(static_cast<std::size_t>(p)), HX(static_cast<std::size_t>(p));
    X[0] = first_vector(sv, want, 0.0);
    for (int k = 1; k < p; ++k) X[static_cast<std::size_t>(k)] = start_vector(n, 100u + static_cast<unsigned>(k));
    orthonormalize(X, 300u);
    int nl = 0;  // converged (locked) leading vectors
    std::vector<double> theta(static_cast<std::size_t>(p), 0.0), resid(static_cast<std::size_t>(p), 0.0);
    std::vector<double> prev_theta;
    Eigenpairs out;
    for (int it = 1; it <= in.max_iterations; ++it) {
        check_cancelled();
        std::vector<std::vector<double>> Y(static_cast<std::size_t>(p));
        int lin_its = 0;
        double lin_res = 0.0;
        for (int k = 0; k < p; ++k) {
            if (k < nl) {
                Y[static_cast<std::size_t>(k)] = X[static_cast<std::size_t>(k)];
                continue;
            }
            const LinearStats ls = sv.solve_h(X[static_cast<std::size_t>(k)], Y[static_cast<std::size_t>(k)], it);
            lin_its += ls.iterations;
            lin_res = std::max(lin_res, ls.residual);
        }
        orthonormalize(Y, 500u + static_cast<unsigned>(it));
        std::vector<std::vector<double>> HY(static_cast<std::size_t>(p));
        for (int k = 0; k < p; ++k) sv.H(Y[static_cast<std::size_t>(k)], HY[static_cast<std::size_t>(k)]);
        if (p == 1) {
            X[0] = Y[0];
            HX[0] = HY[0];
            theta[0] = dot(X[0], HX[0]);
        } else {
            Matrix G(static_cast<std::size_t>(p), static_cast<std::size_t>(p));
            for (int i = 0; i < p; ++i)
                for (int j = i; j < p; ++j) {
                    const double s = dot(Y[static_cast<std::size_t>(i)], HY[static_cast<std::size_t>(j)]);
                    G(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = s;
                    G(static_cast<std::size_t>(j), static_cast<std::size_t>(i)) = s;
                }
            const EigenResult ev = symmetric_eigen(G);
            for (int k = 0; k < p; ++k) {
                std::vector<double> x(n, 0.0), hx(n, 0.0);
                for (int q = 0; q < p; ++q) {
                    const double c = ev.vectors(static_cast<std::size_t>(q), static_cast<std::size_t>(k));
                    axpy(c, Y[static_cast<std::size_t>(q)], x);
                    axpy(c, HY[static_cast<std::size_t>(q)], hx);
                }
                X[static_cast<std::size_t>(k)] = std::move(x);
                HX[static_cast<std::size_t>(k)] = std::move(hx);
                theta[static_cast<std::size_t>(k)] = ev.values[static_cast<std::size_t>(k)];
            }
        }
        for (int k = 0; k < want; ++k) {
            std::vector<double> r = HX[static_cast<std::size_t>(k)];
            axpy(-theta[static_cast<std::size_t>(k)], X[static_cast<std::size_t>(k)], r);
            resid[static_cast<std::size_t>(k)] = norm(r);
        }
        while (nl < want && resid[static_cast<std::size_t>(nl)] < in.residual_tolerance) ++nl;
        res.iterations = it;
        NuclearSolverStep step;
        step.iteration = it;
        step.lowest = theta[0];
        step.change = nl < want ? resid[static_cast<std::size_t>(nl)] : 0.0;
        step.converged = nl;
        step.level = want == 1 ? 1 : 0;
        step.linear_iterations = lin_its;
        step.linear_residual = lin_res;
        step.seconds = sv.elapsed();
        res.history.push_back(step);
        sv.step_line(step, lin_its);
        report_progress(JobStage::NuclearEigen, it, theta[0], step.change);
        report_progress_text(fmt("Inverse iteration %d · %d of %d levels converged · residual %.1e E_h · %d linear-solver iterations",
                                 it, nl, want, step.change, lin_its));
        if (nl >= want) {
            res.converged = true;
            break;
        }
    }
    for (int k = 0; k < want; ++k) {
        out.values.push_back(theta[static_cast<std::size_t>(k)]);
        out.vectors.push_back(X[static_cast<std::size_t>(k)]);
        out.residuals.push_back(resid[static_cast<std::size_t>(k)]);
    }
    return out;
}

// 21.0: Rayleigh-quotient iteration (original: RQiteration): (H - theta) y = x with the Rayleigh
// quotient theta = x^T H x as the shift, x <- y/||y||; cubic convergence to the level nearest to
// theta. The levels one after the other, on the complement of the converged ones (deflation);
// each from inverse-iteration steps (H y = x) until ||H x - theta x|| < 0.01 theta, so that the
// shift is near the lowest remaining level (from the Gaussian start vector, the first level
// directly). The levels are sorted at the end. With conjugate gradients or the Chebyshev
// iteration the shift is theta - ||r|| (see below).
Eigenpairs rayleigh_quotient(Solver& sv) {
    const NuclearMotionInput& in = sv.in;
    NuclearMotionResult& res = sv.res;
    const std::size_t n = sv.n();
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    std::vector<std::vector<double>> locked;
    std::vector<double> values, residuals;
    int it_total = 0;
    bool all_converged = true;
    // conjugate gradients and the Chebyshev iteration need a positive definite matrix: the shift
    // theta - ||r||, below the level when ||r|| is below its gap (Kato-Temple: theta - E <=
    // ||r||^2 / gap), so that (H - shift) is positive definite on the complement of the converged
    // levels; quadratic instead of cubic convergence. MINRES and BiCG: the shift theta.
    const bool spd = sv.lin == Linear::Pcg || sv.lin == Linear::Chebyshev;
    for (int lev = 0; lev < want; ++lev) {
        std::vector<double> x = (lev == 0 && !sv.guess.empty()) ? sv.guess : start_vector(n, 20u + static_cast<unsigned>(lev));
        orthogonalize(x, locked);
        scale(x, 1.0 / norm(x));
        std::vector<double> hx, r, y;
        sv.H(x, hx);
        double theta = dot(x, hx);
        r = hx;
        axpy(-theta, x, r);
        double rho = norm(r);
        bool rqi = lev == 0 && !sv.guess.empty();
        int warm = 0, rejected = 0;
        bool converged = rho < in.residual_tolerance;
        auto inverse_step = [&](const std::vector<double>& from, std::vector<double>& out) {
            const LinearStats st = sv.solve_h(from, out, it_total);
            orthogonalize(out, locked);
            return st;
        };
        for (int it = 1; it <= in.max_iterations && !converged; ++it) {
            check_cancelled();
            ++it_total;
            LinearStats ls;
            const bool shifted = rqi && rejected < 2;
            if (!shifted) {
                ls = inverse_step(x, y);
                ++warm;
            } else {
                std::vector<const std::vector<double>*> Q;
                for (const auto& q : locked) Q.push_back(&q);
                ls = sv.solve_shifted(spd ? theta - rho : theta, Q, x, y, in.linear_tolerance, in.max_linear_iterations, it_total);
            }
            double ny = norm(y);
            if (!(ny > 0.0) || !std::isfinite(ny)) {
                // the shifted solve gave nothing (a breakdown at its first step): an inverse-iteration step
                ls = inverse_step(x, y);
                ny = norm(y);
            }
            std::vector<double> xn = y, hxn;
            scale(xn, 1.0 / ny);
            sv.H(xn, hxn);
            double theta_n = dot(xn, hxn);
            std::vector<double> rn = hxn;
            axpy(-theta_n, xn, rn);
            double rho_n = norm(rn);
            if (shifted && !(rho_n < rho)) {
                // safeguard: a shifted step that does not lower the residual (an inaccurate solve of
                // the nearly singular system, e.g. by the Chebyshev iteration) is replaced by an
                // inverse-iteration step from the previous vector; after two, the level is finished
                // by inverse iteration
                ++rejected;
                ls = inverse_step(x, y);
                xn = y;
                scale(xn, 1.0 / norm(xn));
                sv.H(xn, hxn);
                theta_n = dot(xn, hxn);
                rn = hxn;
                axpy(-theta_n, xn, rn);
                rho_n = norm(rn);
                if (rejected == 2)
                    res.notes.push_back(fmt("Rayleigh-quotient iteration, level %d: the shifted solves stopped lowering the residual; "
                                            "finished by inverse iteration",
                                            lev + 1));
            }
            x = std::move(xn);
            hx = std::move(hxn);
            theta = theta_n;
            r = std::move(rn);
            rho = rho_n;
            converged = rho < in.residual_tolerance;
            if (!rqi && (rho < 1e-2 * theta || warm >= 50)) rqi = true;
            NuclearSolverStep step;
            step.iteration = it_total;
            step.lowest = theta;
            step.change = rho;
            step.converged = lev + (converged ? 1 : 0);
            step.level = lev + 1;
            step.linear_iterations = ls.iterations;
            step.linear_residual = ls.residual;
            step.seconds = sv.elapsed();
            res.history.push_back(step);
            sv.step_line(step, ls.iterations);
            report_progress(JobStage::NuclearEigen, it_total, theta, rho);
            report_progress_text(fmt("Rayleigh-quotient iteration %d · level %d of %d · %s · E = %.6f cm-1 · residual %.1e E_h",
                                     it_total, lev + 1, want, shifted ? "shifted" : "inverse step", theta * kHartreeToCm, rho));
        }
        if (!converged) all_converged = false;
        locked.push_back(x);
        values.push_back(theta);
        residuals.push_back(rho);
    }
    res.iterations = it_total;
    res.converged = all_converged;
    // sorted by energy
    std::vector<std::size_t> order(values.size());
    for (std::size_t k = 0; k < order.size(); ++k) order[k] = k;
    std::sort(order.begin(), order.end(), [&](std::size_t p, std::size_t q) { return values[p] < values[q]; });
    bool reordered = false;
    Eigenpairs out;
    for (std::size_t k = 0; k < order.size(); ++k) {
        if (order[k] != k) reordered = true;
        out.values.push_back(values[order[k]]);
        out.vectors.push_back(std::move(locked[order[k]]));
        out.residuals.push_back(residuals[order[k]]);
    }
    if (reordered) res.notes.push_back("Rayleigh-quotient iteration found the levels out of order; they are sorted");
    return out;
}

Eigenpairs dense(Solver& sv) {
    const Csr& a = sv.a;
    const NuclearMotionInput& in = sv.in;
    NuclearMotionResult& res = sv.res;
    const std::size_t n = a.n;
    if (n > 4000) throw std::invalid_argument(fmt("dense diagonalisation is limited to 4000 basis functions (here %zu)", n));
    Matrix h(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t p = a.row[i]; p < a.row[i + 1]; ++p) h(i, static_cast<std::size_t>(a.col[p])) = a.val[p];
    report_progress_text(fmt("Dense diagonalisation of the %zu x %zu Hamiltonian", n, n));
    const EigenResult e = symmetric_eigen(h);
    Eigenpairs out;
    const std::size_t want = std::min<std::size_t>(static_cast<std::size_t>(in.levels), n);
    for (std::size_t k = 0; k < want; ++k) {
        out.values.push_back(e.values[k]);
        std::vector<double> v(n);
        for (std::size_t i = 0; i < n; ++i) v[i] = e.vectors(i, k);
        out.vectors.push_back(std::move(v));
    }
    out.residuals = residuals_of(sv, out.values, out.vectors);
    res.iterations = 1;
    res.converged = true;
    return out;
}

// 21.0: the Gaussian start vector of the thesis's Algorithm 5 (the original's psikiind): at every
// grid point the bond lengths r1, r2 and the angle a of the bond coordinates (step 1), the product
// f(r1) f(r2) f(a) of Gaussians (steps 2, 3) and the normalisation (step 4).
//   sigma (Eq. 6.36): the standard deviation of the harmonic ground state of a bond and of the
//     angle, sigma^2 = (G/k)^1/2 / 2, with the curvature k of the surface at its minimum and
//     Wilson's G (G_r = 1/m(H) + 1/m(O); G_a = 2/(m(H) r^2) + (2 - 2 cos a)/(m(O) r^2)).
//   mu: the vibrationally averaged values, as the thesis prescribes, estimated from the minimum
//     by the first-order shift of a cubic oscillator, <x - x_e> = -f3 sigma^2/(2 k) (f3 the third
//     derivative); for H2(16)O 1.844 bohr, the original program's averaged 1.8437 bohr.
//   f: the amplitude whose square has the standard deviation sigma, exp(-(x - mu)^2/(4 sigma^2)),
//     i.e. the normal distribution of Eq. 6.33 with the width 2^1/2 sigma (with the width sigma
//     itself the vector is too narrow: overlap 0.89 instead of 0.98-0.99 with the ground state).
// The DVR coefficient of the function is c = s w^1/2 f at the point (s the sign of the DVR
// function, w its weight); the original alternated the sign from point to point instead. For
// J > 0 the K blocks are combined by the lowest eigenvector of H in their span (Rayleigh-Ritz).
std::vector<double> gaussian_start_vector(const Problem& p, const IsotopologueData& iso, const NuclearMotionInput& in,
                                          Solver& sv, NuclearMotionResult& res) {
    const double rq = iso.rmin, aq = iso.amin * kPi / 180.0;
    auto V = [&](double r1, double r2, double a) { return cvrqd_potential(in.isotopologue, r1, r2, std::cos(a)); };
    // second derivatives (h = 1e-3) and third derivatives (h = 1e-2) along a bond and the angle
    const double h = 1e-3, h3 = 1e-2;
    const double v0 = V(rq, rq, aq);
    const double kr = (V(rq + h, rq, aq) - 2.0 * v0 + V(rq - h, rq, aq)) / (h * h);
    const double ka = (V(rq, rq, aq + h) - 2.0 * v0 + V(rq, rq, aq - h)) / (h * h);
    const double fr = (V(rq + 2 * h3, rq, aq) - 2.0 * V(rq + h3, rq, aq) + 2.0 * V(rq - h3, rq, aq) - V(rq - 2 * h3, rq, aq)) /
                      (2.0 * h3 * h3 * h3);
    const double fa = (V(rq, rq, aq + 2 * h3) - 2.0 * V(rq, rq, aq + h3) + 2.0 * V(rq, rq, aq - h3) - V(rq, rq, aq - 2 * h3)) /
                      (2.0 * h3 * h3 * h3);
    const double mh = iso.hydrogen * kAmuToMe, mo = iso.oxygen * kAmuToMe;
    const double gr = 1.0 / mh + 1.0 / mo;
    const double ga = 2.0 / (mh * rq * rq) + (2.0 - 2.0 * std::cos(aq)) / (mo * rq * rq);
    const double sr = std::sqrt(0.5 * std::sqrt(gr / kr)), sa = std::sqrt(0.5 * std::sqrt(ga / ka));
    const double re = rq - fr * sr * sr / (2.0 * kr), ae = aq - fa * sa * sa / (2.0 * ka);
    res.start_mu[0] = res.start_mu[1] = re;
    res.start_mu[2] = ae;
    res.start_sigma[0] = res.start_sigma[1] = sr;
    res.start_sigma[2] = sa;
    std::vector<std::vector<double>> blocks;
    for (int b = 0; b < p.nb; ++b) {
        std::vector<double> g(p.n, 0.0);
        const LegendreDvr& leg = p.leg[static_cast<std::size_t>(b)];
        for (int i = 0; i < p.ns; ++i)
            for (int j = 0; j < p.ns; ++j)
                for (int l = 0; l < p.np; ++l) {
                    double r1, r2, cb;
                    p.coords.to_bond(p.R[static_cast<std::size_t>(i)], p.R[static_cast<std::size_t>(j)], leg.q[static_cast<std::size_t>(l)], r1, r2, cb);
                    const double a = std::acos(std::max(-1.0, std::min(1.0, cb)));
                    const double f = std::exp(-0.25 * ((r1 - re) * (r1 - re) + (r2 - re) * (r2 - re)) / (sr * sr)) *
                                     std::exp(-0.25 * (a - ae) * (a - ae) / (sa * sa));
                    const double s = p.sr[static_cast<std::size_t>(i)] * p.sr[static_cast<std::size_t>(j)] * leg.s[static_cast<std::size_t>(l)];
                    const double w = p.wr[static_cast<std::size_t>(i)] * p.wr[static_cast<std::size_t>(j)] * leg.w[static_cast<std::size_t>(l)];
                    g[p.index(b, i, j, l)] = s * std::sqrt(w) * f;
                }
        const double ng = norm(g);
        if (ng > 0.0) {
            scale(g, 1.0 / ng);
            blocks.push_back(std::move(g));
        }
    }
    if (blocks.empty()) throw std::runtime_error("the Gaussian start vector vanishes on the grid");
    std::vector<double> q;
    if (blocks.size() == 1) {
        q = std::move(blocks[0]);
    } else {
        const std::size_t m = blocks.size();
        std::vector<std::vector<double>> hb(m);
        for (std::size_t k = 0; k < m; ++k) sv.H(blocks[k], hb[k]);
        Matrix G(m, m);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) G(i, j) = 0.5 * (dot(blocks[i], hb[j]) + dot(blocks[j], hb[i]));
        const EigenResult ev = symmetric_eigen(G);
        q.assign(p.n, 0.0);
        for (std::size_t k = 0; k < m; ++k) axpy(ev.vectors(k, 0), blocks[k], q);
        scale(q, 1.0 / norm(q));
    }
    std::vector<double> hq;
    sv.H(q, hq);
    res.start_rayleigh = dot(q, hq);
    return q;
}

// 21.0: the cuts and densities of an eigenvector (WavefunctionPlot).
WavefunctionPlot wavefunction_plot(const Problem& p, const std::vector<double>& x) {
    WavefunctionPlot w;
    const int ns = p.ns, np = p.np;
    auto psi = [&](int b, int i, int j, int l) {
        const LegendreDvr& leg = p.leg[static_cast<std::size_t>(b)];
        const double s = p.sr[static_cast<std::size_t>(i)] * p.sr[static_cast<std::size_t>(j)] * leg.s[static_cast<std::size_t>(l)];
        const double wt = p.wr[static_cast<std::size_t>(i)] * p.wr[static_cast<std::size_t>(j)] * leg.w[static_cast<std::size_t>(l)];
        return x[p.index(b, i, j, l)] * s / std::sqrt(wt);
    };
    double best = -1.0, top = 1.0;
    for (int b = 0; b < p.nb; ++b)
        for (int i = 0; i < ns; ++i)
            for (int j = 0; j < ns; ++j)
                for (int l = 0; l < np; ++l) {
                    const double v = psi(b, i, j, l);
                    if (std::fabs(v) > best) {
                        best = std::fabs(v);
                        top = v;
                        w.block = b;
                        w.i = i;
                        w.j = j;
                        w.l = l;
                    }
                }
    if (!(best > 0.0)) return w;
    const std::size_t uns = static_cast<std::size_t>(ns), unp = static_cast<std::size_t>(np);
    w.psi_r1r2.resize(uns * uns);
    w.psi_r1t.resize(uns * unp);
    w.psi_r2t.resize(uns * unp);
    w.rho_r1r2.assign(uns * uns, 0.0f);
    w.rho_r1t.assign(uns * unp, 0.0f);
    w.rho_r2t.assign(uns * unp, 0.0f);
    for (int i = 0; i < ns; ++i)
        for (int j = 0; j < ns; ++j) w.psi_r1r2[static_cast<std::size_t>(i) * uns + static_cast<std::size_t>(j)] = static_cast<float>(psi(w.block, i, j, w.l) / top);
    for (int i = 0; i < ns; ++i)
        for (int l = 0; l < np; ++l) {
            w.psi_r1t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] = static_cast<float>(psi(w.block, i, w.j, l) / top);
            w.psi_r2t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] = static_cast<float>(psi(w.block, w.i, i, l) / top);
        }
    // densities: |c|^2 / (weights); theta weights w(cos) / sin(theta)
    std::vector<double> d12(uns * uns, 0.0), d1t(uns * unp, 0.0), d2t(uns * unp, 0.0);
    for (int b = 0; b < p.nb; ++b) {
        const LegendreDvr& leg = p.leg[static_cast<std::size_t>(b)];
        for (int i = 0; i < ns; ++i)
            for (int j = 0; j < ns; ++j)
                for (int l = 0; l < np; ++l) {
                    const double c = x[p.index(b, i, j, l)];
                    const double c2 = c * c;
                    d12[static_cast<std::size_t>(i) * uns + static_cast<std::size_t>(j)] += c2;
                    if (b != w.block) continue;
                    const double ql = leg.q[static_cast<std::size_t>(l)];
                    const double wt = leg.w[static_cast<std::size_t>(l)] / std::sqrt(std::max(1e-300, 1.0 - ql * ql));
                    d1t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] += c2 / wt;
                    d2t[static_cast<std::size_t>(j) * unp + static_cast<std::size_t>(l)] += c2 / wt;
                }
    }
    for (int i = 0; i < ns; ++i) {
        const double wi = p.wr[static_cast<std::size_t>(i)];
        for (int j = 0; j < ns; ++j)
            w.rho_r1r2[static_cast<std::size_t>(i) * uns + static_cast<std::size_t>(j)] =
                static_cast<float>(d12[static_cast<std::size_t>(i) * uns + static_cast<std::size_t>(j)] / (wi * p.wr[static_cast<std::size_t>(j)]));
        for (int l = 0; l < np; ++l) {
            w.rho_r1t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] = static_cast<float>(d1t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] / wi);
            w.rho_r2t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] = static_cast<float>(d2t[static_cast<std::size_t>(i) * unp + static_cast<std::size_t>(l)] / wi);
        }
    }
    return w;
}

void expectation_values(const Problem& p, const Eigenpairs& e, NuclearMotionResult& res, bool moments) {
    res.levels.clear();
    const bool vectors = !e.vectors.empty();
    for (std::size_t k = 0; k < e.values.size(); ++k) {
        NuclearLevel lv;
        lv.energy = e.values[k];
        lv.residual = k < e.residuals.size() ? e.residuals[k] : std::numeric_limits<double>::quiet_NaN();
        if (!vectors) {
            res.levels.push_back(lv);
            continue;
        }
        const std::vector<double>& x = e.vectors[k];
        double nrm = 0.0;
        for (double c : x) nrm += c * c;
        lv.k_weights.assign(static_cast<std::size_t>(p.nb), 0.0);
        double s1[6] = {0, 0, 0, 0, 0, 0}, s2[6] = {0, 0, 0, 0, 0, 0}, th = 0.0, ct = 0.0;
        for (int b = 0; b < p.nb; ++b)
            for (int i = 0; i < p.ns; ++i)
                for (int j = 0; j < p.ns; ++j)
                    for (int l = 0; l < p.np; ++l) {
                        const double w = x[p.index(b, i, j, l)] * x[p.index(b, i, j, l)] / nrm;
                        lv.k_weights[static_cast<std::size_t>(b)] += w;
                        if (!moments) continue;
                        double r1, r2, cb;
                        p.coords.to_bond(p.R[static_cast<std::size_t>(i)], p.R[static_cast<std::size_t>(j)],
                                         p.leg[static_cast<std::size_t>(b)].q[static_cast<std::size_t>(l)], r1, r2, cb);
                        cb = std::max(-1.0, std::min(1.0, cb));
                        const double a1[6] = {r1, r1 * r1, r1 * r1 * r1, 1.0 / r1, 1.0 / (r1 * r1), 1.0 / (r1 * r1 * r1)};
                        const double a2[6] = {r2, r2 * r2, r2 * r2 * r2, 1.0 / r2, 1.0 / (r2 * r2), 1.0 / (r2 * r2 * r2)};
                        for (int q = 0; q < 6; ++q) {
                            s1[q] += w * a1[q];
                            s2[q] += w * a2[q];
                        }
                        th += w * std::acos(cb);
                        ct += w * cb;
                    }
        if (moments) {
            auto fill = [](const double* s, DistanceMoments& m) {
                m.mean = s[0];
                m.rms = std::sqrt(s[1]);
                m.cubic = std::cbrt(s[2]);
                m.inverse = 1.0 / s[3];
                m.inverse2 = 1.0 / std::sqrt(s[4]);
                m.inverse3 = 1.0 / std::cbrt(s[5]);
            };
            fill(s1, lv.r1);
            fill(s2, lv.r2);
            lv.theta = th;
            lv.cos_theta = ct;
        }
        if (k < 60) lv.plot = wavefunction_plot(p, x);
        res.levels.push_back(lv);
    }
}

std::string solver_description(const NuclearMotionInput& in) {
    const std::string pc = in.preconditioner == "ilu" ? "ILU (symmetric Gauss-Seidel factors)"
                           : in.preconditioner == "jacobi" ? "Jacobi (diagonal)"
                                                           : "none";
    const std::string lin = in.linear_solver == "minres" ? "MINRES"
                            : in.linear_solver == "chebyshev" ? "Chebyshev iteration"
                            : in.linear_solver == "bicg" ? "BiCG"
                                                         : "preconditioned conjugate gradients";
    if (in.eigensolver == "lanczos")
        return "shift-invert Lanczos with full reorthogonalisation; H x = q by " + lin + ", preconditioner " + pc;
    if (in.eigensolver == "jacobi-davidson")
        return fmt("Jacobi-Davidson; correction equation: the preconditioner (%s), then %d steps of %s on the projected equation",
                   pc.c_str(), in.jd_inner_iterations, lin.c_str());
    if (in.eigensolver == "davidson") {
        if (in.davidson_inner_iterations > 0)
            return fmt("Davidson; correction: %d steps of %s on (H - theta) t = -r, preconditioner %s", in.davidson_inner_iterations,
                       lin.c_str(), pc.c_str());
        return "Davidson; correction: the preconditioner " + pc + " of H - theta";
    }
    if (in.eigensolver == "inverse-power")
        return std::string(in.levels == 1 ? "inverse iteration" : "inverse subspace iteration with Rayleigh-Ritz") + "; H y = x by " +
               lin + ", preconditioner " + pc;
    if (in.eigensolver == "rayleigh-quotient")
        return "Rayleigh-quotient iteration with deflation; (H - theta) y = x by " + lin + ", preconditioner " + pc +
               " (of H - theta); inverse-iteration steps first";
    return "dense diagonalisation (Householder tridiagonalisation and implicit QL)";
}

// 20.0: the picture of the sparsity of rows and columns first ... first + span - 1 (see
// MatrixPicture); bins = min(span, max_bins).
MatrixPicture matrix_picture(const Csr& a, std::size_t first, std::size_t span, int max_bins) {
    MatrixPicture m;
    if (span == 0) return m;
    const int bins = static_cast<int>(std::min<std::size_t>(span, static_cast<std::size_t>(max_bins)));
    m.bins = bins;
    m.first = static_cast<long long>(first);
    m.span = static_cast<long long>(span);
    const std::size_t nb = static_cast<std::size_t>(bins);
    std::vector<unsigned> count(nb * nb, 0);
    std::vector<float> largest(nb * nb, 0.0f);
    auto bin = [&](std::size_t k) { return k * nb / span; };
    for (std::size_t r = first; r < first + span; ++r) {
        const std::size_t br = bin(r - first);
        for (std::size_t p = a.row[r]; p < a.row[r + 1]; ++p) {
            const std::size_t c = static_cast<std::size_t>(a.col[p]);
            if (c < first || c >= first + span || a.val[p] == 0.0) continue;
            const std::size_t q = br * nb + bin(c - first);
            ++count[q];
            largest[q] = std::max(largest[q], static_cast<float>(std::fabs(a.val[p])));
        }
    }
    // elements of each bin row/column: the k with bin(k) = b
    std::vector<std::size_t> width(nb, 0);
    for (std::size_t k = 0; k < span; ++k) ++width[bin(k)];
    double lo = std::numeric_limits<double>::infinity(), hi = -lo;
    for (std::size_t q = 0; q < nb * nb; ++q)
        if (count[q] > 0) {
            const double l = std::log10(static_cast<double>(largest[q]));
            lo = std::min(lo, l);
            hi = std::max(hi, l);
        }
    if (!(hi >= lo)) lo = hi = 0.0;
    if (hi - lo < 1e-6) hi = lo + 1.0;
    m.log_min = lo;
    m.log_max = hi;
    m.density.assign(nb * nb, 0);
    m.magnitude.assign(nb * nb, 0);
    for (std::size_t i = 0; i < nb; ++i)
        for (std::size_t j = 0; j < nb; ++j) {
            const std::size_t q = i * nb + j;
            if (count[q] == 0) continue;
            const double frac = static_cast<double>(count[q]) / static_cast<double>(width[i] * width[j]);
            m.density[q] = static_cast<unsigned char>(std::max(1.0, std::min(255.0, std::round(255.0 * frac))));
            const double t = (std::log10(static_cast<double>(largest[q])) - lo) / (hi - lo);
            m.magnitude[q] = static_cast<unsigned char>(1 + std::lround(254.0 * std::max(0.0, std::min(1.0, t))));
        }
    return m;
}

}  // namespace

std::string isotopologue_name(CvrqdIsotopologue iso) { return isotopologue_data(iso).name; }

CvrqdIsotopologue isotopologue_from_name(const std::string& s) {
    std::string t;
    for (char c : s)
        if (c != '(' && c != ')' && c != ' ' && c != '_' && c != '-') t += static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (t == "H216O" || t == "H2O" || t == "1") return CvrqdIsotopologue::H2_16O;
    if (t == "H217O" || t == "2") return CvrqdIsotopologue::H2_17O;
    if (t == "D216O" || t == "D2O" || t == "3") return CvrqdIsotopologue::D2_16O;
    if (t == "H218O" || t == "4") return CvrqdIsotopologue::H2_18O;
    throw std::invalid_argument("isotopologue must be H2(16)O, H2(17)O, H2(18)O or D2(16)O");
}

NuclearMotionEstimate estimate_nuclear_motion(const NuclearMotionInput& in) {
    NuclearMotionEstimate e;
    try {
        validate(in);
    } catch (const std::exception& ex) {
        e.refusal = ex.what();
        return e;
    }
    const long long nb = static_cast<long long>(k_values(in.J, in.parity_block).size());
    const long long ns = in.radial_points, np = in.angular_points;
    e.dimension = nb * ns * ns * np;
    e.nonzeros = e.dimension * (2 * ns + np - 2 + (nb > 1 ? 2 * np : 0));
    double vectors = 0.0;
    if (in.eigensolver == "lanczos") vectors = std::min<double>(in.max_iterations, static_cast<double>(e.dimension)) + 8;
    else if (in.eigensolver == "dense") vectors = static_cast<double>(e.dimension);
    else if (in.eigensolver == "inverse-power")
        vectors = 3.0 * (in.levels == 1 ? 1 : in.levels + std::max(2, in.levels / 2)) + 12;
    else if (in.eigensolver == "rayleigh-quotient") vectors = in.levels + 14.0;
    else vectors = 2.0 * (in.subspace_size > 0 ? in.subspace_size : std::max(2 * in.levels + 20, 40)) + 10;
    e.memory_mb = (static_cast<double>(e.nonzeros) * 12.0 + static_cast<double>(e.dimension) * (8.0 * vectors + 40.0)) / 1048576.0;
    if (in.eigensolver == "dense" && e.dimension > 4000)
        e.refusal = fmt("dense diagonalisation is limited to 4000 basis functions (here %lld); choose Lanczos or (Jacobi-)Davidson",
                        e.dimension);
    if (e.dimension > 20000000) e.refusal = fmt("the basis has %lld functions; the limit is 20 million", e.dimension);
    return e;
}

NuclearMotionResult run_nuclear_motion(const NuclearMotionInput& in, const NuclearOutputSink& sink) {
    validate(in);
    const NuclearMotionEstimate est = estimate_nuclear_motion(in);
    if (!est.refusal.empty()) throw std::invalid_argument(est.refusal);
    const auto t0 = std::chrono::steady_clock::now();
    NuclearMotionResult res;
    res.input = in;
    const IsotopologueData iso = isotopologue_data(in.isotopologue);
    res.isotopologue_name = iso.name;
    res.potential_minimum = iso.vmin;
    res.minimum_r = iso.rmin;
    res.minimum_angle = iso.amin * kPi / 180.0;
    res.solver_description = solver_description(in);
    const int threads = thread_count(in.threads);
    std::ostringstream os;
    const bool radau = in.coordinates == "radau";

    emit(sink, os, "========================================================================================\n");
    emit(sink, os, "                      ΨLab 21.0 · nuclear motion (rovibrational levels)\n");
    emit(sink, os, "========================================================================================\n\n");
    emit(sink, os, fmt("Molecule          %s on the CVRQD potential energy surface\n", iso.name));
    emit(sink, os, "                  (Barletta, Shirin, Zobov, Polyansky, Tennyson, Valeev and Csaszar,\n"
                   "                  J. Chem. Phys. 125, 204307 (2006)); energies from its minimum\n");
    emit(sink, os, fmt("Minimum           V = %.12f E_h at r(OH) = %.6f bohr (%.6f A), angle HOH = %.4f deg\n", iso.vmin,
                       iso.rmin, iso.rmin * kBohrToAngstrom, iso.amin));
    report_progress(JobStage::NuclearHamiltonian);
    report_progress_text("Discrete variable representation");
    Problem p = build_problem(in, iso, res);
    const char* hname = iso.hydrogen > 1.5 ? "D" : "H";
    emit(sink, os, fmt("Nuclear masses    %.9f u (%s), %.9f u (%s), %.7f u (O)\n", res.masses[0], hname, res.masses[1], hname,
                       res.masses[2]));
    if (!radau) {
        emit(sink, os, fmt("Coordinates       Jacobi: R1 = r(O-H1), R2 = H2 to the centre of mass of O-H1, theta;\n"
                           "                  g1 = %.8f, g2 = 0 (Sutcliffe-Tennyson);\n"
                           "                  reduced masses mu1 = %.6f, mu2 = %.6f m_e (1 u = %.4f m_e)\n",
                           res.g1, res.mu1, res.mu2, kAmuToMe));
    } else {
        const double alpha = std::sqrt(res.masses[2] / (res.masses[0] + res.masses[1] + res.masses[2]));
        emit(sink, os, fmt("Coordinates       Radau: R1 = r(K-H1), R2 = r(K-H2), theta = angle H1-K-H2, K the canonical point\n"
                           "                  K = alpha O + (1 - alpha) (centre of mass of H1, H2), alpha = (m(O)/M)^1/2 = %.8f;\n"
                           "                  g1 = %.8f, g2 = %.8f (Sutcliffe-Tennyson, thesis Table 3.1);\n"
                           "                  reduced masses mu1 = m(H1) = %.6f, mu2 = m(H2) = %.6f m_e (1 u = %.4f m_e)\n",
                           alpha, res.g1, res.g2, res.mu1, res.mu2, kAmuToMe));
    }
    emit(sink, os, fmt("Rotation          J = %d, K = %s (%s)\n", in.J,
                       in.J == 0 ? "0" : (in.parity_block == 0 ? fmt("0 ... %d", in.J).c_str() : fmt("1 ... %d", in.J).c_str()),
                       in.J == 0 ? "no rotation" : "parity-adapted Wang functions; body-fixed z axis along R1"));
    emit(sink, os, fmt("Radial DVR        %d Hermite points in R1 and in R2, %.4f ... %.4f bohr (Re = %.4f, Rmax = %.4f bohr)\n",
                       in.radial_points, p.R.front(), p.R.back(), in.r_centre, in.r_half_width));
    emit(sink, os, fmt("Angular DVR       %d (associated) Legendre points in cos(theta)\n", in.angular_points));
    res.dimension = static_cast<long long>(p.n);
    const auto th0 = std::chrono::steady_clock::now();
    Csr a = build_hamiltonian(p, in, iso, threads);
    res.hamiltonian_seconds = seconds_since(th0);
    res.nonzeros = static_cast<long long>(a.val.size());
    res.matrix_whole = matrix_picture(a, 0, p.n, 1024);
    res.matrix_detail = matrix_picture(a, 0, std::min<std::size_t>(p.n, 1024), 1024);
    emit(sink, os, fmt("Hamiltonian       %lld x %lld, %lld non-zero elements (%.2f per row), built in %.2f s\n",
                       res.dimension, res.dimension, res.nonzeros, static_cast<double>(res.nonzeros) / static_cast<double>(res.dimension),
                       res.hamiltonian_seconds));
    emit(sink, os, fmt("Eigensolver       %s\n", res.solver_description.c_str()));
    emit(sink, os, fmt("Eigenvectors      %s\n", in.eigenvectors ? "computed (expectation values, Eigenvectors tab)"
                                                  : (in.eigensolver == "lanczos"
                                                         ? "not computed: the eigenvalues only (no Ritz vectors; no residual check)"
                                                         : (in.eigensolver == "dense" ? "not kept: the eigenvalues only"
                                                                                       : "not kept: the eigenvalues only (the method needs them internally)"))));
    emit(sink, os, fmt("Convergence       eigenvalues %.1e E_h, eigenvector residual ||Hx - Ex|| %.1e E_h", in.eigenvalue_tolerance,
                       in.residual_tolerance));
    if (in.eigensolver != "dense") emit(sink, os, fmt(", linear solves %.1e", in.linear_tolerance));
    emit(sink, os, fmt("\nThreads           %d\n", threads));

    Solver sv(a, in, res, sink, os);
    if (in.gaussian_start && in.eigensolver != "dense") {
        report_progress_text("Gaussian start vector (Algorithm 5)");
        sv.guess = gaussian_start_vector(p, iso, in, sv, res);
        emit(sink, os, fmt("Start vector      Gaussian of Algorithm 5 (psikiind): mu = %.5f bohr, %.5f bohr, %.3f deg (averaged:\n"
                           "                  minimum + cubic anharmonic shift); sigma = %.5f bohr, %.5f bohr, %.3f deg (harmonic\n"
                           "                  ground state; amplitude exp(-(x - mu)^2/(4 sigma^2))); Rayleigh quotient %.4f cm-1\n",
                           res.start_mu[0], res.start_mu[1], res.start_mu[2] * 180.0 / kPi, res.start_sigma[0], res.start_sigma[1],
                           res.start_sigma[2] * 180.0 / kPi, res.start_rayleigh * kHartreeToCm));
    } else if (in.eigensolver != "dense") {
        emit(sink, os, "Start vector      pseudo-random (reproducible)\n");
    }
    emit(sink, os, "\n");

    const auto ts0 = std::chrono::steady_clock::now();
    sv.t0 = ts0;
    Eigenpairs e;
    if (in.eigensolver == "dense") {
        e = dense(sv);
    } else {
        if (in.eigensolver == "lanczos")
            emit(sink, os, "  Lanczos  lowest level (cm-1)   change (E_h)  lin.it.  conv.\n");
        else if (in.eigensolver == "inverse-power" || in.eigensolver == "rayleigh-quotient")
            emit(sink, os, fmt("     Iter  %s (cm-1)  residual (E_h)  lin.it.  conv.\n",
                               in.eigensolver == "inverse-power" ? "lowest level" : " this level "));
        else
            emit(sink, os, "     Iter  lowest level (cm-1)  residual (E_h)  subspace  conv.\n");
        if (in.eigensolver == "lanczos") e = lanczos(sv);
        else if (in.eigensolver == "inverse-power") e = inverse_power(sv);
        else if (in.eigensolver == "rayleigh-quotient") e = rayleigh_quotient(sv);
        else e = davidson(sv, in.eigensolver == "jacobi-davidson");
        emit(sink, os, "\n");
    }
    res.solver_seconds = seconds_since(ts0);
    if (!sv.guess.empty() && !e.vectors.empty()) res.start_overlap = std::fabs(dot(sv.guess, e.vectors.front()));
    if (!in.eigenvectors) e.vectors.clear();
    res.has_vectors = !e.vectors.empty();
    report_progress_text(res.has_vectors ? "Expectation values" : "Energy levels");
    expectation_values(p, e, res, in.expectation_values);

    const bool lanczos_like = in.eigensolver == "lanczos";
    const bool linear = in.eigensolver != "dense";
    emit(sink, os, fmt("%s after %d %s, %lld products H x%s, %.2f s\n", res.converged ? "Converged" : "NOT CONVERGED",
                       res.iterations, lanczos_like ? "Lanczos steps" : "iterations", res.matrix_vector_products,
                       linear ? fmt(" (%lld linear-solver iterations)", res.linear_iterations).c_str() : "",
                       res.solver_seconds));
    if (in.linear_solver == "chebyshev" && res.chebyshev_max > 0.0)
        emit(sink, os, fmt("Chebyshev         eigenvalues of M^-1 A estimated in [%.4e, %.4e]\n", res.chebyshev_min, res.chebyshev_max));
    if (res.start_overlap >= 0.0)
        emit(sink, os, fmt("Start vector      overlap with the lowest eigenvector |<q1|x1>| = %.6f\n", res.start_overlap));
    emit(sink, os, "\n");
    // levels
    emit(sink, os, "Energy levels (from the minimum of the potential; E - E1 from the lowest level of this run)\n\n");
    emit(sink, os, "     n          E (E_h)            E (cm-1)       E - E1 (cm-1)   ||Hx-Ex|| (E_h)\n");
    const double e1 = res.levels.empty() ? 0.0 : res.levels.front().energy;
    for (std::size_t k = 0; k < res.levels.size(); ++k) {
        const NuclearLevel& lv = res.levels[k];
        emit(sink, os, fmt("  %4zu  %16.12f  %18.6f  %16.6f   %12s\n", k + 1, lv.energy, lv.energy * kHartreeToCm,
                           (lv.energy - e1) * kHartreeToCm, std::isfinite(lv.residual) ? fmt("%.2e", lv.residual).c_str() : "-"));
    }
    if (res.has_vectors && in.expectation_values && !res.levels.empty()) {
        emit(sink, os, "\nExpectation values (angstrom, degrees): r1 = r(O-H1), r2 = r(O-H2)\n\n");
        emit(sink, os, "     n  <r1>     <r1^2>^1/2  <r1^3>^1/3  <r1^-1>^-1  <r1^-2>^-1/2  <r1^-3>^-1/3   <HOH>\n");
        for (std::size_t k = 0; k < res.levels.size(); ++k) {
            const NuclearLevel& lv = res.levels[k];
            auto row = [&](const char* label, const DistanceMoments& m, bool angle) {
                emit(sink, os, fmt("  %4s  %8.5f  %10.5f  %10.5f  %10.5f  %12.5f  %12.5f", label, m.mean * kBohrToAngstrom,
                                   m.rms * kBohrToAngstrom, m.cubic * kBohrToAngstrom, m.inverse * kBohrToAngstrom,
                                   m.inverse2 * kBohrToAngstrom, m.inverse3 * kBohrToAngstrom));
                emit(sink, os, angle ? fmt("  %8.3f\n", lv.theta * 180.0 / kPi) : std::string("\n"));
            };
            row(fmt("%zu r1", k + 1).c_str(), lv.r1, true);
            row("  r2", lv.r2, false);
        }
    }
    if (res.has_vectors && in.J > 0 && !res.levels.empty()) {
        emit(sink, os, "\nWeights of the K blocks\n\n     n");
        for (int K : p.K) emit(sink, os, fmt("     K=%-3d", K));
        emit(sink, os, "\n");
        for (std::size_t k = 0; k < res.levels.size(); ++k) {
            emit(sink, os, fmt("  %4zu", k + 1));
            for (double w : res.levels[k].k_weights) emit(sink, os, fmt("  %8.5f", w));
            emit(sink, os, "\n");
        }
    }
    for (const std::string& note : res.notes) emit(sink, os, "\nNote: " + note);
    if (!res.notes.empty()) emit(sink, os, "\n");
    emit(sink, os, fmt("\n%s termination of ΨLab 21.0 (nuclear motion), %.2f s\n", res.converged ? "Normal" : "Unconverged",
                       seconds_since(t0)));
    res.report = os.str();
    report_progress(JobStage::Done);
    return res;
}

}  // namespace macqc
