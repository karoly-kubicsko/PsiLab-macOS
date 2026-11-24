// ΨLab 20.0: nuclear motion (see nuclear_motion.h). The Hamiltonian, its DVR and the solvers
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
// Preconditioners (original: CGM with precondtype 1, 2, 3; ILUprecond)

enum class Precond { None, Jacobi, Ilu };

Precond parse_precond(const std::string& s) {
    if (s == "none") return Precond::None;
    if (s == "jacobi" || s == "diagonal") return Precond::Jacobi;
    if (s == "ilu" || s == "sgs") return Precond::Ilu;
    throw std::invalid_argument("preconditioner must be none, jacobi or ilu");
}

// z = M^-1 r. ILU: the factors of the original, M = (D + L) D^-1 (D + U) with the diagonal D and
// the strict triangles L, U of H (symmetric Gauss-Seidel); positive definite when H is.
void apply_precond(const Csr& a, Precond p, const std::vector<double>& r, std::vector<double>& z) {
    const std::size_t n = a.n;
    z.resize(n);
    switch (p) {
    case Precond::None:
        z = r;
        return;
    case Precond::Jacobi:
        for (std::size_t i = 0; i < n; ++i) z[i] = r[i] / a.diag[i];
        return;
    case Precond::Ilu: {
        std::vector<double> w(n);
        for (std::size_t i = 0; i < n; ++i) {
            double s = 0.0;
            for (std::size_t q = a.row[i]; q < a.diag_pos[i]; ++q) s += a.val[q] * w[static_cast<std::size_t>(a.col[q])];
            w[i] = (r[i] - s) / a.diag[i];
        }
        for (std::size_t ii = n; ii-- > 0;) {
            double s = 0.0;
            for (std::size_t q = a.diag_pos[ii] + 1; q < a.row[ii + 1]; ++q) s += a.val[q] * z[static_cast<std::size_t>(a.col[q])];
            z[ii] = w[ii] - s / a.diag[ii];
        }
        return;
    }
    }
}

// ---------------------------------------------------------------------------------------------
// Linear solvers for H x = b (original: CGM, minres)

struct LinearStats {
    int iterations = 0;
    double residual = 0.0;
    bool converged = false;
};

// Preconditioned conjugate gradients from x = 0; stops when ||b - H x|| <= tol.
LinearStats pcg(const Csr& a, Precond pc, const std::vector<double>& b, std::vector<double>& x, double tol,
                int maxit, long long& matvecs) {
    const std::size_t n = a.n;
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r = b, z, p, q(n);
    apply_precond(a, pc, r, z);
    p = z;
    double rho = dot(r, z);
    st.residual = norm(r);
    if (st.residual <= tol) {
        st.converged = true;
        return st;
    }
    for (int it = 1; it <= maxit; ++it) {
        if ((it & 63) == 0) check_cancelled();
        a.multiply(p, q);
        ++matvecs;
        const double pq = dot(p, q);
        if (!(pq > 0.0)) throw std::runtime_error("conjugate gradients: the Hamiltonian is not positive definite");
        const double alpha = rho / pq;
        axpy(alpha, p, x);
        axpy(-alpha, q, r);
        st.iterations = it;
        st.residual = norm(r);
        if (st.residual <= tol) {
            st.converged = true;
            return st;
        }
        apply_precond(a, pc, r, z);
        const double rho_new = dot(r, z);
        const double beta = rho_new / rho;
        rho = rho_new;
        for (std::size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
    }
    return st;
}

// MINRES (Paige and Saunders, SIAM J. Numer. Anal. 12, 617 (1975)) for a symmetric operator,
// with a symmetric positive definite preconditioner; x0 = 0. `op(v, out)` applies the matrix.
// Stops when the preconditioned residual estimate falls below tol.
template <class Op, class Prec>
LinearStats minres(std::size_t n, const Op& op, const Prec& prec, const std::vector<double>& b, std::vector<double>& x,
                   double tol, int maxit) {
    LinearStats st;
    x.assign(n, 0.0);
    std::vector<double> r1 = b, y, r2, v, w(n, 0.0), w1(n, 0.0), w2(n, 0.0);
    prec(r1, y);
    double beta1 = dot(r1, y);
    if (beta1 < 0.0) throw std::runtime_error("MINRES: the preconditioner is not positive definite");
    beta1 = std::sqrt(beta1);
    st.residual = beta1;
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
        if (phibar <= tol) {
            st.converged = true;
            return st;
        }
        if (beta == 0.0) break;
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

// Bond coordinates from the Jacobi coordinates (original: jacbolbelso2): r(O-H1) = R1,
// r(O-H2) and cos(H-O-H); g = m(H)/(m(O) + m(H)).
inline void jacobi_to_bond(double g, double R1, double R2, double c, double& r1, double& r2, double& cb) {
    r1 = R1;
    r2 = std::sqrt(g * g * R1 * R1 + R2 * R2 + 2.0 * g * R1 * R2 * c);
    cb = R2 / r2 * c + g * R1 / r2;
}

struct Problem {
    int ns = 0, np = 0, nb = 0;
    std::size_t nx = 0, n = 0;
    std::vector<int> K;
    std::vector<double> R;              // radial points (bohr)
    Matrix t1, t2;                      // radial kinetic energy
    std::vector<double> r1, r2;         // 1 / (2 mu R^2)
    std::vector<LegendreDvr> leg;       // per block
    std::vector<Matrix> coriolis;       // between blocks b and b + 1
    double g = 0;                       // m(H)/(m(O)+m(H)) of the R1 diatom
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
    if (s != "lanczos" && s != "jacobi-davidson" && s != "davidson" && s != "dense")
        throw std::invalid_argument("eigensolver must be lanczos, jacobi-davidson, davidson or dense");
    if (in.linear_solver != "pcg" && in.linear_solver != "minres")
        throw std::invalid_argument("the linear solver must be pcg or minres");
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

    // masses: xmass = {H, H, O} as in the original; R1 joins atoms 2 and 3 (O-H), R2 joins atom 1
    // to their centre of mass
    const double m0 = iso.hydrogen, m1 = iso.hydrogen, m2 = iso.oxygen;
    const double g1 = m1 / (m1 + m2);
    const double tred1 = kAmuToMe / (1.0 / m1 + 1.0 / m2);
    const double tred2 = kAmuToMe / (1.0 / m0 + g1 * g1 / m1 + (1.0 - g1) * (1.0 - g1) / m2);
    p.g = m1 / (m2 + m1);
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
    const std::vector<double> x = symmetric_eigen(hx).values;
    const double allk = x.back() / in.r_half_width;
    p.R.resize(static_cast<std::size_t>(ns));
    for (int i = 0; i < ns; ++i) p.R[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(i)] / allk + in.r_centre;
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
                    jacobi_to_bond(p.g, p.R[static_cast<std::size_t>(i)], p.R[static_cast<std::size_t>(j)],
                                   p.leg[static_cast<std::size_t>(b)].q[static_cast<std::size_t>(l)], r1, r2, cb);
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

// The start vector of the iterative solvers: reproducible pseudo-random numbers (the original
// read a Gaussian guess from thing.dat, written by psikiind).
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
    std::vector<std::vector<double>> vectors;
    std::vector<double> residuals;
};

void emit(const NuclearOutputSink& sink, std::ostringstream& os, const std::string& text) {
    os << text;
    if (sink) sink(text);
}

std::vector<double> residuals_of(const Csr& a, const std::vector<double>& values,
                                 const std::vector<std::vector<double>>& vectors, long long& matvecs) {
    std::vector<double> res(values.size());
    std::vector<double> hy;
    for (std::size_t k = 0; k < values.size(); ++k) {
        a.multiply(vectors[k], hy);
        ++matvecs;
        axpy(-values[k], vectors[k], hy);
        res[k] = norm(hy);
    }
    return res;
}

// Shift-invert Lanczos (original: lanczos): the Krylov space of H^-1, every step a linear solve,
// full reorthogonalisation; the eigenvalues E = 1/theta of the largest Ritz values theta.
Eigenpairs lanczos(const Csr& a, const NuclearMotionInput& in, NuclearMotionResult& res, const NuclearOutputSink& sink,
                   std::ostringstream& os) {
    const std::size_t n = a.n;
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    const int mmax = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.max_iterations), n));
    const Precond pc = parse_precond(in.preconditioner);
    const bool use_minres = in.linear_solver == "minres";
    std::vector<std::vector<double>> V;
    V.reserve(static_cast<std::size_t>(mmax) + 1);
    V.push_back(start_vector(n, 20u));
    std::vector<double> alpha, beta, prev, w;
    Eigenpairs out;
    for (int j = 0; j < mmax; ++j) {
        check_cancelled();
        LinearStats ls;
        if (use_minres) {
            auto op = [&](const std::vector<double>& v, std::vector<double>& y) {
                a.multiply(v, y);
                ++res.matrix_vector_products;
            };
            auto pre = [&](const std::vector<double>& r, std::vector<double>& z) { apply_precond(a, pc, r, z); };
            ls = minres(n, op, pre, V[static_cast<std::size_t>(j)], w, in.linear_tolerance, in.max_linear_iterations);
        } else {
            ls = pcg(a, pc, V[static_cast<std::size_t>(j)], w, in.linear_tolerance, in.max_linear_iterations, res.matrix_vector_products);
        }
        res.linear_iterations += ls.iterations;
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
            if ((change < in.eigenvalue_tolerance && static_cast<int>(E.size()) == want) || be < 1e-14 || m == mmax) {
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
                out.residuals = residuals_of(a, out.values, out.vectors, res.matrix_vector_products);
                int conv = 0;
                for (double r : out.residuals)
                    if (r < in.residual_tolerance) ++conv;
                step.converged = conv;
                if (conv == static_cast<int>(out.values.size()) && static_cast<int>(out.values.size()) == want) done = true;
                if (be < 1e-14 || m == mmax) done = true;
                res.converged = conv == want && static_cast<int>(out.values.size()) == want;
            }
            report_progress(JobStage::NuclearEigen, m, E.front(), std::isfinite(change) ? change : 0.0);
            report_progress_text(fmt("Lanczos step %d · E1 = %.6f cm-1 · change %.1e E_h · %d CG iterations", m,
                                     E.front() * kHartreeToCm, std::isfinite(change) ? change : 0.0, ls.iterations));
        } else {
            report_progress(JobStage::NuclearEigen, m, 0.0, 0.0);
            report_progress_text(fmt("Lanczos step %d of at least %d · %d linear-solver iterations", m, want, ls.iterations));
        }
        res.history.push_back(step);
        emit(sink, os,
             fmt("  %5d  %18.10f  %12.3e  %8d  %5d\n", m, step.lowest * kHartreeToCm, step.change, ls.iterations,
                 step.converged));
        if (done) break;
        beta.push_back(be);
        scale(w, 1.0 / be);
        V.push_back(w);
    }
    return out;
}

// Davidson and Jacobi-Davidson for the lowest eigenpairs (original: davidson, jacobidavidson).
Eigenpairs davidson(const Csr& a, const NuclearMotionInput& in, bool jacobi_davidson, NuclearMotionResult& res,
                    const NuclearOutputSink& sink, std::ostringstream& os) {
    const std::size_t n = a.n;
    const int want = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(in.levels), n));
    const int maxdim = static_cast<int>(std::min<std::size_t>(
        static_cast<std::size_t>(in.subspace_size > 0 ? in.subspace_size : std::max(2 * want + 20, 40)), n));
    const int keep = std::min(maxdim - 1, std::max(want + 5, want * 3 / 2));
    std::vector<std::vector<double>> V, AV;
    // start: the reproducible vector and unit vectors at the smallest diagonal elements
    {
        std::vector<std::size_t> order(n);
        for (std::size_t i = 0; i < n; ++i) order[i] = i;
        const std::size_t m0 = std::min<std::size_t>(static_cast<std::size_t>(want), n);
        std::partial_sort(order.begin(), order.begin() + static_cast<long>(m0), order.end(),
                          [&](std::size_t x, std::size_t y) { return a.diag[x] < a.diag[y]; });
        std::vector<std::vector<double>> init;
        init.push_back(start_vector(n, 20u));
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
            a.multiply(v, hv);
            ++res.matrix_vector_products;
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
        res.history.push_back(step);
        report_progress(JobStage::NuclearEigen, it, ev.values.front(), step.change);
        report_progress_text(fmt("%s iteration %d · %d of %d levels converged · residual %.1e E_h",
                                 jacobi_davidson ? "Jacobi-Davidson" : "Davidson", it, conv, want, step.change));
        emit(sink, os, fmt("  %5d  %18.10f  %12.3e  %8d  %5d\n", it, step.lowest * kHartreeToCm, step.change, k, conv));
        if (conv >= want) {
            res.converged = true;
            break;
        }
        if (k >= static_cast<int>(n)) break;
        // correction vector for the target pair (theta, u, r)
        std::vector<double> dinv(n);
        for (std::size_t i = 0; i < n; ++i) {
            double d = a.diag[i] - theta;
            if (std::fabs(d) < 1e-8) d = d < 0 ? -1e-8 : 1e-8;
            dinv[i] = 1.0 / d;
        }
        t.assign(n, 0.0);
        if (!jacobi_davidson) {
            for (std::size_t i = 0; i < n; ++i) t[i] = -r[i] * dinv[i];
        } else {
            // one-step approximation (original: jdcorrectionequation_tuti): t = -(K^-1 r - mu K^-1 u)
            // with mu such that t is orthogonal to u
            std::vector<double> y(n), rh(n);
            for (std::size_t i = 0; i < n; ++i) {
                y[i] = dinv[i] * u[i];
                rh[i] = dinv[i] * r[i];
            }
            const double mu = dot(u, rh) / dot(u, y);
            for (std::size_t i = 0; i < n; ++i) t[i] = -(rh[i] - mu * y[i]);
            if (in.jd_inner_iterations > 0) {
                // refine: (I - QQ^T)(H - theta)(I - QQ^T) d = -r - (...) t, Q = locked pairs and u
                std::vector<std::vector<double>> Q = locked;
                Q.push_back(u);
                auto project = [&](std::vector<double>& x) {
                    for (const auto& q : Q) axpy(-dot(x, q), q, x);
                };
                auto op = [&](const std::vector<double>& v, std::vector<double>& y2) {
                    std::vector<double> pv = v;
                    project(pv);
                    a.multiply(pv, y2);
                    ++res.matrix_vector_products;
                    axpy(-theta, pv, y2);
                    project(y2);
                };
                auto ident = [](const std::vector<double>& x, std::vector<double>& z) { z = x; };
                project(t);
                std::vector<double> rhs = r, at;
                scale(rhs, -1.0);
                op(t, at);
                axpy(-1.0, at, rhs);
                project(rhs);
                std::vector<double> d;
                const LinearStats ls = minres(n, op, ident, rhs, d, in.linear_tolerance, in.jd_inner_iterations);
                res.linear_iterations += ls.iterations;
                axpy(1.0, d, t);
                project(t);
            }
        }
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
        if (nt < 1e-12) {
            // stagnation: add the residual instead
            t = r;
            orthogonalize(t, V);
            nt = norm(t);
            if (nt < 1e-14) break;
        }
        scale(t, 1.0 / nt);
        std::vector<double> ht;
        a.multiply(t, ht);
        ++res.matrix_vector_products;
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
        out.residuals = residuals_of(a, out.values, out.vectors, res.matrix_vector_products);
    }
    return out;
}

Eigenpairs dense(const Csr& a, const NuclearMotionInput& in, NuclearMotionResult& res) {
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
    out.residuals = residuals_of(a, out.values, out.vectors, res.matrix_vector_products);
    res.iterations = 1;
    res.converged = true;
    return out;
}

void expectation_values(const Problem& p, const Eigenpairs& e, NuclearMotionResult& res, bool moments) {
    res.levels.clear();
    for (std::size_t k = 0; k < e.values.size(); ++k) {
        NuclearLevel lv;
        lv.energy = e.values[k];
        lv.residual = k < e.residuals.size() ? e.residuals[k] : 0.0;
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
                        jacobi_to_bond(p.g, p.R[static_cast<std::size_t>(i)], p.R[static_cast<std::size_t>(j)],
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
        res.levels.push_back(lv);
    }
}

std::string solver_description(const NuclearMotionInput& in) {
    auto pc = [&]() -> std::string {
        if (in.preconditioner == "ilu") return "ILU (symmetric Gauss-Seidel factors of H)";
        if (in.preconditioner == "jacobi") return "Jacobi (diagonal)";
        return "none";
    };
    if (in.eigensolver == "lanczos")
        return std::string("shift-invert Lanczos with full reorthogonalisation; H x = q by ") +
               (in.linear_solver == "minres" ? "MINRES" : "preconditioned conjugate gradients") + ", preconditioner " + pc();
    if (in.eigensolver == "jacobi-davidson")
        return fmt("Jacobi-Davidson; correction equation: diagonal preconditioner and %d MINRES steps", in.jd_inner_iterations);
    if (in.eigensolver == "davidson") return "Davidson with the diagonal preconditioner";
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

    emit(sink, os, "========================================================================================\n");
    emit(sink, os, "                      ΨLab 20.0 · nuclear motion (rovibrational levels)\n");
    emit(sink, os, "========================================================================================\n\n");
    emit(sink, os, fmt("Molecule          %s on the CVRQD potential energy surface\n", iso.name));
    emit(sink, os, "                  (Barletta, Shirin, Zobov, Polyansky, Tennyson, Valeev and Csaszar,\n"
                   "                  J. Chem. Phys. 125, 204307 (2006)); energies from its minimum\n");
    emit(sink, os, fmt("Minimum           V = %.12f E_h at r(OH) = %.6f bohr (%.6f A), angle HOH = %.4f deg\n", iso.vmin,
                       iso.rmin, iso.rmin * kBohrToAngstrom, iso.amin));
    report_progress(JobStage::NuclearHamiltonian);
    report_progress_text("Discrete variable representation");
    Problem p = build_problem(in, iso, res);
    emit(sink, os, fmt("Nuclear masses    %.9f u (%s), %.9f u (%s), %.7f u (O)\n", res.masses[0],
                       iso.hydrogen > 1.5 ? "D" : "H", res.masses[1], iso.hydrogen > 1.5 ? "D" : "H", res.masses[2]));
    emit(sink, os, fmt("Coordinates       Jacobi: R1 = r(O-H1), R2 = H2 to the centre of mass of O-H1, theta;\n"
                       "                  reduced masses mu1 = %.6f, mu2 = %.6f m_e (1 u = %.4f m_e)\n",
                       res.mu1, res.mu2, kAmuToMe));
    emit(sink, os, fmt("Rotation          J = %d, K = %s (%s)\n", in.J,
                       in.J == 0 ? "0" : (in.parity_block == 0 ? fmt("0 ... %d", in.J).c_str() : fmt("1 ... %d", in.J).c_str()),
                       in.J == 0 ? "no rotation" : "parity-adapted Wang functions; body-fixed z axis along R1"));
    emit(sink, os, fmt("Radial DVR        %d Hermite points in R1 and in R2, %.4f ... %.4f bohr (Re = %.4f, Rmax = %.4f)\n",
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
    emit(sink, os, fmt("Convergence       eigenvalues %.1e E_h, eigenvector residual ||Hx - Ex|| %.1e E_h", in.eigenvalue_tolerance,
                       in.residual_tolerance));
    if (in.eigensolver == "lanczos") emit(sink, os, fmt(", linear solves %.1e", in.linear_tolerance));
    emit(sink, os, fmt("\nThreads           %d\n\n", threads));

    const auto ts0 = std::chrono::steady_clock::now();
    Eigenpairs e;
    if (in.eigensolver == "dense") {
        e = dense(a, in, res);
    } else {
        if (in.eigensolver == "lanczos")
            emit(sink, os, "  Lanczos  lowest level (cm-1)   change (E_h)  lin.it.  conv.\n");
        else
            emit(sink, os, "     Iter  lowest level (cm-1)  residual (E_h)  subspace  conv.\n");
        e = in.eigensolver == "lanczos" ? lanczos(a, in, res, sink, os)
                                        : davidson(a, in, in.eigensolver == "jacobi-davidson", res, sink, os);
        emit(sink, os, "\n");
    }
    res.solver_seconds = seconds_since(ts0);
    report_progress_text("Expectation values");
    expectation_values(p, e, res, in.expectation_values);

    emit(sink, os, fmt("%s after %d %s, %lld products H x%s, %.2f s\n\n", res.converged ? "Converged" : "NOT CONVERGED",
                       res.iterations, in.eigensolver == "lanczos" ? "Lanczos steps" : "iterations", res.matrix_vector_products,
                       in.eigensolver == "lanczos" ? fmt(" (%lld linear-solver iterations)", res.linear_iterations).c_str() : "",
                       res.solver_seconds));
    // levels
    emit(sink, os, "Energy levels (from the minimum of the potential; E - E1 from the lowest level of this run)\n\n");
    emit(sink, os, "     n          E (E_h)            E (cm-1)       E - E1 (cm-1)   ||Hx-Ex|| (E_h)\n");
    const double e1 = res.levels.empty() ? 0.0 : res.levels.front().energy;
    for (std::size_t k = 0; k < res.levels.size(); ++k) {
        const NuclearLevel& lv = res.levels[k];
        emit(sink, os, fmt("  %4zu  %16.12f  %18.6f  %16.6f   %12.2e\n", k + 1, lv.energy, lv.energy * kHartreeToCm,
                           (lv.energy - e1) * kHartreeToCm, lv.residual));
    }
    if (in.expectation_values && !res.levels.empty()) {
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
    if (in.J > 0 && !res.levels.empty()) {
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
    emit(sink, os, fmt("\n%s termination of ΨLab 20.0 (nuclear motion), %.2f s\n", res.converged ? "Normal" : "Unconverged",
                       seconds_since(t0)));
    res.report = os.str();
    report_progress(JobStage::Done);
    return res;
}

}  // namespace macqc
