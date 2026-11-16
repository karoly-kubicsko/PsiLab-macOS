// Analytic CCSD Hessian (ΨLab 10.0, stage 10.3).
//
// The CCSD amplitudes t and multipliers l satisfy Omega(t; H) = 0 and L_t(t, l; H) = 0, with
// Omega the residual of the closed-shell amplitude equations (ccsd_amplitude_residual) and
// L(t, l; H) = <D(t, l), H> the CCSD Lagrangian (Lagrangian densities D, linear in H, linear in l,
// a polynomial of degree 4 in t). The response of the densities enters the Hessian as
// <D^(y), H^(x)> (correlated_hessian.h). With J t^x = -Omega(t; H^(x)) and the derivative of the
// Lambda equations, L_t(H^(y)) + L_tt t^y + L_tl l^y = 0, the multiplier responses drop out
// (2n+2 rule):
//   <D^(y), H^(x)> = L_t(H^(x)) t^y + L_t(H^(y)) t^x + t^x L_tt t^y,
// L_t(H) w = <dD(t + s w, l)/ds, H>. Only the amplitude responses t^q along the internal
// directions are solved, by the Jacobi iteration with the Fock-diagonal denominators (J ~ -D)
// and DIIS. The Jacobian products and the directional derivatives of the densities are exact:
// for a polynomial of degree <= 4 the symmetric five-point formula
//   f'(0) = [8 (f(s) - f(-s)) - (f(2s) - f(-2s))] / (12 s)
// has no truncation error (only rounding, with s d of order 0.1), and the pair terms
// t^k L_tt t^l follow from the five-point second derivatives of L along t^k, t^l and t^k + t^l.
#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "macqc/cc.h"
#include "macqc/cc_gradient.h"
#include "macqc/correlated_hessian.h"
#include "macqc/gradient.h"
#include "macqc/gto_integrals.h"
#include "macqc/job.h"
#include "macqc/parallel.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

using Vec = std::vector<double>;

// A pair of single and double blocks (t1-like [i][a], t2-like [i][j][a][b]).
struct Amp {
    Vec a, b;
};

double max_abs(const Amp& x) {
    double m = 0.0;
    for (double v : x.a) m = std::max(m, std::fabs(v));
    for (double v : x.b) m = std::max(m, std::fabs(v));
    return m;
}

Amp combine(const Amp& x, double fx, const Amp& y, double fy) {
    Amp z;
    z.a.resize(x.a.size());
    z.b.resize(x.b.size());
    for (std::size_t k = 0; k < x.a.size(); ++k) z.a[k] = fx * x.a[k] + fy * y.a[k];
    for (std::size_t k = 0; k < x.b.size(); ++k) z.b[k] = fx * x.b[k] + fy * y.b[k];
    return z;
}

// Exact derivative d/ds f(x0 + s d) at s = 0 for f polynomial of degree <= 4 (five-point formula).
template <class Value, class F, class Combine>
Value poly_derivative(F f, const Amp& x0, const Amp& d, Combine comb) {
    const double dm = max_abs(d);
    const double s = dm > 0.0 ? 0.1 / dm : 1.0;
    const Value p1 = f(combine(x0, 1.0, d, s)), m1 = f(combine(x0, 1.0, d, -s));
    const Value p2 = f(combine(x0, 1.0, d, 2.0 * s)), m2 = f(combine(x0, 1.0, d, -2.0 * s));
    // [8 (p1 - m1) - (p2 - m2)] / (12 s)
    return comb({&p1, &m1, &p2, &m2}, {8.0 / (12.0 * s), -8.0 / (12.0 * s), -1.0 / (12.0 * s), 1.0 / (12.0 * s)});
}

Amp amp_comb(std::initializer_list<const Amp*> xs, std::initializer_list<double> fs) {
    Amp z;
    z.a.assign((*xs.begin())->a.size(), 0.0);
    z.b.assign((*xs.begin())->b.size(), 0.0);
    auto f = fs.begin();
    for (const Amp* x : xs) {
        for (std::size_t k = 0; k < z.a.size(); ++k) z.a[k] += *f * x->a[k];
        for (std::size_t k = 0; k < z.b.size(); ++k) z.b[k] += *f * x->b[k];
        ++f;
    }
    return z;
}

CorrelatedDensity density_comb(std::initializer_list<const CorrelatedDensity*> xs, std::initializer_list<double> fs) {
    const CorrelatedDensity& first = **xs.begin();
    const std::size_t n = first.D.rows();
    CorrelatedDensity z;
    z.D = Matrix(n, n);
    z.G = Packed4(n);
    auto f = fs.begin();
    for (const CorrelatedDensity* x : xs) {
        for (std::size_t p = 0; p < n; ++p)
            for (std::size_t q = 0; q < n; ++q) z.D(p, q) += *f * x->D(p, q);
        for (std::size_t k = 0; k < z.G.v.size(); ++k) z.G.v[k] += *f * x->G.v[k];
        ++f;
    }
    return z;
}

// Pulay DIIS on (vector, error) pairs.
class Diis {
public:
    explicit Diis(std::size_t size) : size_(size) {}
    Vec extrapolate(const Vec& x, const Vec& err) {
        xs_.push_back(x);
        es_.push_back(err);
        if (xs_.size() > size_) {
            xs_.erase(xs_.begin());
            es_.erase(es_.begin());
        }
        const std::size_t m = xs_.size();
        if (m < 2) return x;
        Matrix B(m + 1, m + 1);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < m; ++j) {
                double s = 0.0;
                for (std::size_t k = 0; k < es_[i].size(); ++k) s += es_[i][k] * es_[j][k];
                B(i, j) = s;
            }
        for (std::size_t i = 0; i < m; ++i) B(i, m) = B(m, i) = -1.0;
        // solve B c = (0, ..., 0, -1) by Gaussian elimination
        std::vector<double> A((m + 1) * (m + 1)), rhs(m + 1, 0.0);
        for (std::size_t i = 0; i <= m; ++i)
            for (std::size_t j = 0; j <= m; ++j) A[i * (m + 1) + j] = B(i, j);
        rhs[m] = -1.0;
        const std::size_t n = m + 1;
        for (std::size_t c = 0; c < n; ++c) {
            std::size_t piv = c;
            for (std::size_t r = c + 1; r < n; ++r)
                if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) piv = r;
            if (std::fabs(A[piv * n + c]) < 1.0e-300) return x;
            if (piv != c) {
                for (std::size_t k = 0; k < n; ++k) std::swap(A[piv * n + k], A[c * n + k]);
                std::swap(rhs[piv], rhs[c]);
            }
            for (std::size_t r = c + 1; r < n; ++r) {
                const double f = A[r * n + c] / A[c * n + c];
                for (std::size_t k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
                rhs[r] -= f * rhs[c];
            }
        }
        std::vector<double> c(n);
        for (std::size_t r = n; r-- > 0;) {
            double s = rhs[r];
            for (std::size_t k = r + 1; k < n; ++k) s -= A[r * n + k] * c[k];
            c[r] = s / A[r * n + r];
        }
        Vec out(x.size(), 0.0);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t k = 0; k < x.size(); ++k) out[k] += c[i] * xs_[i][k];
        return out;
    }

private:
    std::size_t size_;
    std::vector<Vec> xs_, es_;
};

// The ladder term Z(tau) = sum_cd (ac|bd) tau_ijcd of the amplitude equations is linear in tau,
// and along t + s d, tau(s) = tau + s a + s^2 b with a_ijab = d2_ijab + d1_ia t1_jb + t1_ia d1_jb.
// In the odd-part difference formulas only the coefficient of s^1 survives, so Z(tau(s)) may be
// replaced by s Z(a): one ladder application per Jacobian product instead of one per residual.
class NullLadder final : public LadderOperator {
public:
    void apply(const std::vector<double>& tau, std::vector<double>& z) override { z.assign(tau.size(), 0.0); }
    std::string description() const override { return "none"; }
};

class DirectionalLadder final : public LadderOperator {
public:
    DirectionalLadder(LadderOperator& base, std::size_t o, std::size_t v) : base_(base), o_(o), v_(v) {}
    void set_direction(const Amp& t, const Amp& d) {
        const std::size_t o = o_, v = v_;
        Vec a(d.b);
        for (std::size_t i = 0; i < o; ++i)
            for (std::size_t j = 0; j < o; ++j)
                for (std::size_t p = 0; p < v; ++p)
                    for (std::size_t q = 0; q < v; ++q)
                        a[((i * o + j) * v + p) * v + q] += d.a[i * v + p] * t.a[j * v + q] + t.a[i * v + p] * d.a[j * v + q];
        base_.apply(a, za_);
    }
    void set_step(double s) { s_ = s; }
    void apply(const std::vector<double>& tau, std::vector<double>& z) override {
        z.resize(tau.size());
        for (std::size_t k = 0; k < z.size(); ++k) z[k] = s_ * za_[k];
    }
    std::string description() const override { return base_.description(); }

private:
    LadderOperator& base_;
    std::size_t o_, v_;
    Vec za_;
    double s_ = 0.0;
};

struct LinearSolve {
    Amp x;
    bool converged = false;
    int iterations = 0;
    double residual = 0.0;
};

// Solves A x = -b for A ~ -D (diagonal of the Fock denominators) by x <- x + (A x + b) / D with DIIS.
template <class Apply>
LinearSolve solve_preconditioned(Apply apply, const Amp& b, const Amp& D, double tol, int max_iter,
                                 const Amp* start = nullptr) {
    LinearSolve s;
    Amp x;
    if (start) {
        x = *start;
    } else {
        x.a.resize(b.a.size());
        x.b.resize(b.b.size());
        for (std::size_t k = 0; k < b.a.size(); ++k) x.a[k] = b.a[k] / D.a[k];
        for (std::size_t k = 0; k < b.b.size(); ++k) x.b[k] = b.b[k] / D.b[k];
    }
    Diis diis(10);
    for (int it = 1; it <= max_iter; ++it) {
        check_cancelled();
        const Amp Ax = apply(x);
        Amp r = combine(Ax, 1.0, b, 1.0);
        s.iterations = it;
        s.residual = max_abs(r);
        if (s.residual < tol) {
            s.converged = true;
            break;
        }
        Vec xv, ev;
        xv.reserve(x.a.size() + x.b.size());
        ev.reserve(xv.capacity());
        for (std::size_t k = 0; k < x.a.size(); ++k) {
            const double e = r.a[k] / D.a[k];
            xv.push_back(x.a[k] + e);
            ev.push_back(e);
        }
        for (std::size_t k = 0; k < x.b.size(); ++k) {
            const double e = r.b[k] / D.b[k];
            xv.push_back(x.b[k] + e);
            ev.push_back(e);
        }
        const Vec y = diis.extrapolate(xv, ev);
        std::copy(y.begin(), y.begin() + static_cast<std::ptrdiff_t>(x.a.size()), x.a.begin());
        std::copy(y.begin() + static_cast<std::ptrdiff_t>(x.a.size()), y.end(), x.b.begin());
    }
    s.x = std::move(x);
    return s;
}

// Solves J x = -b with the exact operator `exact` (four evaluations of the residual) by defect
// correction: the inner iterations use the cheap approximation `approx` (central difference, two
// evaluations, relative error ~1e-6) plus the correction c = J x_ref - approx(x_ref) of the
// last outer point, so that the inner equations approach the exact ones as x converges.
template <class Exact, class Approx>
LinearSolve solve_defect_corrected(Exact exact, Approx approx, const Amp& b, const Amp& D, double tol, int max_iter) {
    LinearSolve total;
    Amp c;
    c.a.assign(b.a.size(), 0.0);
    c.b.assign(b.b.size(), 0.0);
    const Amp* start = nullptr;
    Amp x;
    for (int pass = 0; pass < 6; ++pass) {
        const double inner_tol = pass == 0 ? std::max(0.5 * tol, 1.0e-6 * max_abs(b)) : 0.5 * tol;
        auto corrected = [&](const Amp& d) { return combine(approx(d), 1.0, c, 1.0); };
        LinearSolve in = solve_preconditioned(corrected, b, D, inner_tol, max_iter - total.iterations, start);
        total.iterations += in.iterations;
        x = std::move(in.x);
        const Amp Jx = exact(x);
        total.residual = max_abs(combine(Jx, 1.0, b, 1.0));
        if (total.residual < tol) {
            total.converged = true;
            break;
        }
        if (total.iterations >= max_iter) break;
        c = combine(Jx, 1.0, approx(x), -1.0);
        start = &x;
    }
    total.x = std::move(x);
    return total;
}

}  // namespace

CorrelatedHessianResult ccsd_hessian(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                                     const MOIntegrals& mo, const std::vector<double>& t1, const std::vector<double>& t2,
                                     const std::vector<double>& l1, const std::vector<double>& l2, int threads,
                                     double memory_mb, double tol) {
    const auto t0 = Clock::now();
    const int nt = resolve_threads(threads);
    const int o = rhf.n_occupied;
    const int v = static_cast<int>(rhf.C.cols()) - o;
    const OrbitalResponse resp = rhf_orbital_response(mol, basis, rhf, nt, memory_mb, 1.0e-10);
    const double orbital_seconds = seconds_since(t0);
    const auto t1s = Clock::now();

    const CCIntegrals ints0 = cc_integrals_from_mo(mo, o);
    const std::unique_ptr<LadderOperator> ladder0 = make_mo_ladder(mo, o);
    const Amp T{t1, t2}, L{l1, l2};
    // Fock-diagonal denominators of the reference
    Amp D;
    {
        const std::size_t oo = static_cast<std::size_t>(o), vv = static_cast<std::size_t>(v);
        D.a.resize(oo * vv);
        D.b.resize(oo * oo * vv * vv);
        for (std::size_t i = 0; i < oo; ++i)
            for (std::size_t a = 0; a < vv; ++a) D.a[i * vv + a] = ints0.fock(i, i) - ints0.fock(oo + a, oo + a);
        for (std::size_t i = 0; i < oo; ++i)
            for (std::size_t j = 0; j < oo; ++j)
                for (std::size_t a = 0; a < vv; ++a)
                    for (std::size_t b = 0; b < vv; ++b)
                        D.b[((i * oo + j) * vv + a) * vv + b] = ints0.fock(i, i) + ints0.fock(j, j) -
                                                                ints0.fock(oo + a, oo + a) - ints0.fock(oo + b, oo + b);
    }
    auto omega = [&](const CCIntegrals& I, LadderOperator& lad, const Amp& t) {
        Amp r;
        ccsd_amplitude_residual(I, lad, t.a, t.b, r.a, r.b, nt);
        return r;
    };
    auto dens = [&](const Amp& t) { return ccsd_lagrangian_densities(o, v, t.a, t.b, L.a, L.b, nt); };
    // The Lagrangian with the reference Hamiltonian and the fixed multipliers, without the
    // reference energy: <D(t, l), H> - E_ref = E_corr(t) + 2 l1 . Omega1 + lt . Omega2, with
    // lt_ijab = 2 l2_ijab - l2_ijba. The ladder part of lt . Omega2 is Z(lt) . tau(t) ((ac|bd) =
    // (ca|db)), so the residual is evaluated without the ladder and Z(lt) once.
    NullLadder no_ladder;
    Vec lt(l2.size()), zlt;
    {
        const std::size_t oo = static_cast<std::size_t>(o), vv = static_cast<std::size_t>(v);
        for (std::size_t i = 0; i < oo; ++i)
            for (std::size_t j = 0; j < oo; ++j)
                for (std::size_t a = 0; a < vv; ++a)
                    for (std::size_t b = 0; b < vv; ++b)
                        lt[((i * oo + j) * vv + a) * vv + b] =
                            2.0 * l2[((i * oo + j) * vv + a) * vv + b] - l2[((i * oo + j) * vv + b) * vv + a];
        ladder0->apply(lt, zlt);
    }
    auto lagrangian = [&](const Amp& t) {
        const Amp r = omega(ints0, no_ladder, t);
        double e = ccsd_correlation_energy(ints0, *ladder0, t.a, t.b);
        for (std::size_t k = 0; k < r.a.size(); ++k) e += 2.0 * l1[k] * r.a[k];
        const std::size_t oo = static_cast<std::size_t>(o), vv = static_cast<std::size_t>(v);
        for (std::size_t i = 0; i < oo; ++i)
            for (std::size_t j = 0; j < oo; ++j)
                for (std::size_t a = 0; a < vv; ++a)
                    for (std::size_t b = 0; b < vv; ++b) {
                        const std::size_t k = ((i * oo + j) * vv + a) * vv + b;
                        e += lt[k] * r.b[k] + zlt[k] * (t.b[k] + t.a[i * vv + a] * t.a[j * vv + b]);
                    }
        return e;
    };
    const CorrelatedDensity D0 = dens(T);
    const double L0 = lagrangian(T);

    DensityResponse dr;
    dr.internal = true;
    dr.symmetric = true;
    dr.directions = vibrational_directions(mol);
    const std::size_t m = dr.directions.size();
    dr.D.resize(m);
    dr.G.resize(m);
    // Jacobian products along d: exact (five-point formula, four residuals) and approximate
    // (central difference with a small step, two residuals, relative error ~1e-6), both with the
    // ladder term applied once along the direction
    DirectionalLadder dl(*ladder0, static_cast<std::size_t>(o), static_cast<std::size_t>(v));
    auto jacobian = [&](const Amp& d) {
        dl.set_direction(T, d);
        const double dm = max_abs(d);
        const double h = dm > 0.0 ? 0.1 / dm : 1.0;
        auto at = [&](double s) {
            dl.set_step(s);
            return omega(ints0, dl, combine(T, 1.0, d, s));
        };
        const Amp p1 = at(h), m1 = at(-h), p2 = at(2.0 * h), m2 = at(-2.0 * h);
        const double w = 1.0 / (12.0 * h);
        return amp_comb({&p1, &m1, &p2, &m2}, {8.0 * w, -8.0 * w, -w, w});
    };
    auto central = [&](const Amp& d) {
        dl.set_direction(T, d);
        const double dm = max_abs(d);
        const double h = dm > 0.0 ? 1.0e-3 / dm : 1.0;
        dl.set_step(h);
        const Amp p = omega(ints0, dl, combine(T, 1.0, d, h));
        dl.set_step(-h);
        const Amp q = omega(ints0, dl, combine(T, 1.0, d, -h));
        return combine(p, 0.5 / h, q, -0.5 / h);
    };
    std::vector<Amp> tq(m);
    bool all_converged = true;
    int max_it = 0;
    double worst = 0.0;
    for (std::size_t k = 0; k < m; ++k) {
        report_progress(JobStage::Hessian, 7, 0.0, static_cast<double>(k + 1));
        const MOIntegrals Hq = perturbed_hamiltonian(mo, resp, dr.directions[k], nt);
        const CCIntegrals intsQ = cc_integrals_from_mo(Hq, o);
        const std::unique_ptr<LadderOperator> ladderQ = make_mo_ladder(Hq, o);
        // amplitude response: J t^q = -Omega(t; H^(q))
        const Amp b = omega(intsQ, *ladderQ, T);
        const LinearSolve ts = solve_defect_corrected(jacobian, central, b, D, tol, 200);
        tq[k] = ts.x;
        // amplitude part of the density response, dD(t + s t^q)/ds
        const CorrelatedDensity dT = poly_derivative<CorrelatedDensity>(dens, T, ts.x, density_comb);
        dr.D[k] = dT.D;
        dr.G[k] = dT.G;
        all_converged = all_converged && ts.converged;
        max_it = std::max(max_it, ts.iterations);
        worst = std::max(worst, ts.residual);
    }
    // pair terms t^k L_tt t^l: second directional derivatives of the Lagrangian (a polynomial of
    // degree 4 in t: the five-point formula is exact), by polarisation
    auto second = [&](const Amp& w) {
        const double dm = max_abs(w);
        if (dm == 0.0) return 0.0;
        const double h = 0.1 / dm;
        const double p1 = lagrangian(combine(T, 1.0, w, h)), m1 = lagrangian(combine(T, 1.0, w, -h));
        const double p2 = lagrangian(combine(T, 1.0, w, 2.0 * h)), m2 = lagrangian(combine(T, 1.0, w, -2.0 * h));
        return (-(p2 + m2) + 16.0 * (p1 + m1) - 30.0 * L0) / (12.0 * h * h);
    };
    dr.pair = Matrix(m, m);
    std::vector<double> diag(m);
    for (std::size_t k = 0; k < m; ++k) {
        check_cancelled();
        diag[k] = second(tq[k]);
        dr.pair(k, k) = diag[k];
    }
    for (std::size_t k = 0; k < m; ++k)
        for (std::size_t l = 0; l < k; ++l) {
            check_cancelled();
            const double c = 0.5 * (second(combine(tq[k], 1.0, tq[l], 1.0)) - diag[k] - diag[l]);
            dr.pair(k, l) = dr.pair(l, k) = c;
        }
    const double method_seconds = seconds_since(t1s);
    CorrelatedHessianResult r = correlated_hessian(mol, basis, rhf, mo, D0, resp, dr, nt, memory_mb);
    r.cphf_converged = resp.cphf_converged;
    r.cphf_iterations = resp.cphf_iterations;
    r.cphf_residual = resp.cphf_residual;
    r.method_converged = all_converged;
    r.method_iterations = max_it;
    r.method_residual = worst;
    r.orbital_seconds = orbital_seconds;
    r.method_seconds = method_seconds;
    r.seconds = seconds_since(t0);
    return r;
}

}  // namespace macqc
