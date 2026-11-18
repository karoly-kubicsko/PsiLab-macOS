// ΨLab 14.0: integrals of semi-local effective core potentials (see ecp.h).
#include "macqc/ecp.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <unordered_map>

namespace macqc {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kMaxProjector = 6;                                  // projectors l < 6 (local part up to l = 6)
constexpr int kMaxLambda = kMaxProjector - 1 + kEcpMaxShellL;     // Bessel orders and harmonics
constexpr int kMaxDegree = kMaxLambda + kMaxProjector + kEcpMaxShellL + 2;  // monomial degrees on the sphere
constexpr int kDefaultPoints = 64;
constexpr double kWindow = 9.0;  // the radial quadrature covers r0 +- 9 / sqrt(P)

inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }
// Index of x^i y^j z^k (i + j + k = l) in the order (l,0,0), (l-1,1,0), (l-1,0,1), (l-2,2,0), ...
inline int cart_index(int /*x*/, int y, int z) { return ((y + z) * (y + z + 1)) / 2 + z; }

double factorial(int n) {
    double r = 1.0;
    for (int k = 2; k <= n; ++k) r *= k;
    return r;
}

double binomial(int n, int k) {
    if (k < 0 || k > n) return 0.0;
    double r = 1.0;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;
    return r;
}

inline double ipow(double x, int n) {
    double r = 1.0;
    for (int k = 0; k < n; ++k) r *= x;
    return r;
}

struct Tables {
    std::vector<std::array<int, 3>> comp[kMaxDegree + 1];
    // Real spherical harmonics normalised on the unit sphere, Y_lm = sqrt((2l+1)/4pi) S_lm with the
    // solid harmonics S_lm of Helgaker, Jorgensen, Olsen (Eqs. 6.4.47-6.4.50), rows m = -l .. l,
    // sparse over the monomials of degree l.
    std::vector<std::vector<std::pair<int, double>>> ylm[kMaxLambda + 1];
    double dfact[2 * kMaxDegree + 4];  // dfact[k + 1] = k!! for k >= -1

    Tables() {
        for (int d = 0; d <= kMaxDegree; ++d)
            for (int x = d; x >= 0; --x)
                for (int y = d - x; y >= 0; --y) comp[d].push_back({x, y, d - x - y});
        dfact[0] = 1.0;  // (-1)!!
        dfact[1] = 1.0;  // 0!!
        for (int k = 1; k < 2 * kMaxDegree + 3; ++k) dfact[k + 1] = k * (k >= 2 ? dfact[k - 1] : 1.0);
        for (int l = 0; l <= kMaxLambda; ++l) {
            const int nc = ncart(l);
            ylm[l].resize(static_cast<std::size_t>(2 * l + 1));
            const double scale = std::sqrt((2.0 * l + 1.0) / (4.0 * kPi));
            for (int m = -l; m <= l; ++m) {
                std::vector<double> row(static_cast<std::size_t>(nc), 0.0);
                const int am = std::abs(m);
                const double norm = 1.0 / (std::pow(2.0, am) * factorial(l)) *
                                    std::sqrt(2.0 * factorial(l + am) * factorial(l - am) / (m == 0 ? 2.0 : 1.0));
                const int w0 = m < 0 ? 1 : 0;
                for (int t = 0; t <= (l - am) / 2; ++t)
                    for (int u = 0; u <= t; ++u)
                        for (int w = w0; w <= am; w += 2) {
                            const int sign_exp = t + (w - w0) / 2;
                            const double c = ((sign_exp % 2) ? -1.0 : 1.0) * std::pow(0.25, t) * binomial(l, t) *
                                             binomial(l - t, am + t) * binomial(t, u) * binomial(am, w);
                            const int x = 2 * t + am - 2 * u - w, y = 2 * u + w, z = l - 2 * t - am;
                            if (x < 0) continue;
                            row[static_cast<std::size_t>(cart_index(x, y, z))] += norm * c;
                        }
                for (int c = 0; c < nc; ++c)
                    if (row[static_cast<std::size_t>(c)] != 0.0)
                        ylm[l][static_cast<std::size_t>(m + l)].emplace_back(c, scale * row[static_cast<std::size_t>(c)]);
            }
        }
    }

    // integral of x^i y^j z^k over the unit sphere
    double sphere(int i, int j, int k) const {
        if ((i | j | k) & 1) return 0.0;
        return 4.0 * kPi * dfact[i] * dfact[j] * dfact[k] / dfact[i + j + k + 2];
    }
};

const Tables& tables() {
    static const Tables t;
    return t;
}

// Gauss-Legendre nodes and weights on [-1, 1].
struct GaussLegendre {
    std::vector<double> x, w;
    explicit GaussLegendre(int n) : x(static_cast<std::size_t>(n)), w(static_cast<std::size_t>(n)) {
        for (int i = 0; i < (n + 1) / 2; ++i) {
            double z = std::cos(kPi * (i + 0.75) / (n + 0.5)), dp = 0.0;
            for (int it = 0; it < 100; ++it) {
                double p0 = 1.0, p1 = z;
                for (int k = 2; k <= n; ++k) {
                    const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = p2;
                }
                if (n == 1) p0 = 1.0;
                dp = n * (z * p1 - p0) / (z * z - 1.0);
                const double dz = p1 / dp;
                z -= dz;
                if (std::fabs(dz) < 1e-16) break;
            }
            x[static_cast<std::size_t>(i)] = -z;
            x[static_cast<std::size_t>(n - 1 - i)] = z;
            w[static_cast<std::size_t>(i)] = w[static_cast<std::size_t>(n - 1 - i)] = 2.0 / ((1.0 - z * z) * dp * dp);
        }
    }
};

const GaussLegendre& gauss_legendre(int n) {
    if (n == kDefaultPoints) {
        static const GaussLegendre g(kDefaultPoints);
        return g;
    }
    thread_local std::unordered_map<int, std::unique_ptr<GaussLegendre>> other;
    auto& p = other[n];
    if (!p) p = std::make_unique<GaussLegendre>(n);
    return *p;
}

// The angular factors of a shell of angular momentum la at a = A - C for the projector l:
// T[c][m][n][lambda] = sum_{pqs, p+q+s = n} B_c(pqs) sum_mu Y_lambda,mu(a^) <Y_lambda,mu Y_lm x^p y^q z^s>,
// B_c(pqs) = C(cx,p) C(cy,q) C(cz,s) (-a_x)^(cx-p) (-a_y)^(cy-q) (-a_z)^(cz-s).
struct Angular {
    int la = 0, l = 0, nlam = 1;  // lambda = 0 .. nlam - 1
    double amag = 0.0;
    double max_abs = 0.0;
    std::vector<double> T;
    std::size_t at(int c, int m, int n, int lam) const {
        return ((static_cast<std::size_t>(c) * (2 * l + 1) + static_cast<std::size_t>(m)) * (la + 1) +
                static_cast<std::size_t>(n)) * nlam + static_cast<std::size_t>(lam);
    }
};

std::shared_ptr<const Angular> make_angular(int la, int l, const double a[3]) {
    const Tables& tb = tables();
    auto r = std::make_shared<Angular>();
    Angular& g = *r;
    g.la = la;
    g.l = l;
    g.amag = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    const bool at_centre = g.amag < 1e-12;
    // at the centre exp(0) = 1 = 4 pi Y00 Y00: only lambda = 0 (any direction a^)
    const int lam_max = at_centre ? 0 : l + la;
    g.nlam = lam_max + 1;
    const double ahat[3] = {at_centre ? 0.0 : a[0] / g.amag, at_centre ? 0.0 : a[1] / g.amag,
                            at_centre ? 1.0 : a[2] / g.amag};
    const int nc = ncart(la), nm = 2 * l + 1;
    g.T.assign(static_cast<std::size_t>(nc) * nm * (la + 1) * g.nlam, 0.0);
    // powers of -a
    double ma[3][kEcpMaxShellL + 1];
    for (int d = 0; d < 3; ++d)
        for (int k = 0; k <= la; ++k) ma[d][k] = ipow(at_centre ? 0.0 : -a[d], k);
    std::vector<double> Q, R, H;
    for (int lam = 0; lam <= lam_max; ++lam) {
        // Q_lambda(r^) = sum_mu Y_lambda,mu(a^) Y_lambda,mu(r^), a polynomial of degree lambda
        Q.assign(static_cast<std::size_t>(ncart(lam)), 0.0);
        for (const auto& row : tb.ylm[lam]) {
            double y = 0.0;
            for (const auto& [idx, c] : row) {
                const auto& e = tb.comp[lam][static_cast<std::size_t>(idx)];
                y += c * ipow(ahat[0], e[0]) * ipow(ahat[1], e[1]) * ipow(ahat[2], e[2]);
            }
            if (y == 0.0) continue;
            for (const auto& [idx, c] : row) Q[static_cast<std::size_t>(idx)] += y * c;
        }
        const int dR = lam + l;
        for (int m = 0; m < nm; ++m) {
            // R = Q_lambda Y_lm, degree lambda + l
            R.assign(static_cast<std::size_t>(ncart(dR)), 0.0);
            for (std::size_t t = 0; t < Q.size(); ++t) {
                if (Q[t] == 0.0) continue;
                const auto& et = tb.comp[lam][t];
                for (const auto& [u, c] : tb.ylm[l][static_cast<std::size_t>(m)]) {
                    const auto& eu = tb.comp[l][static_cast<std::size_t>(u)];
                    R[static_cast<std::size_t>(cart_index(et[0] + eu[0], et[1] + eu[1], et[2] + eu[2]))] += Q[t] * c;
                }
            }
            // H[n][pqs] = int R x^p y^q z^s over the sphere, n = p + q + s <= la (parity lambda + l + n even)
            for (int n = 0; n <= la; ++n) {
                if (((lam + l + n) & 1) || lam > l + n) continue;
                const int npqs = ncart(n);
                H.assign(static_cast<std::size_t>(npqs), 0.0);
                bool any = false;
                for (int q = 0; q < npqs; ++q) {
                    const auto& e = tb.comp[n][static_cast<std::size_t>(q)];
                    double h = 0.0;
                    for (std::size_t v = 0; v < R.size(); ++v) {
                        if (R[v] == 0.0) continue;
                        const auto& ev = tb.comp[dR][v];
                        h += R[v] * tb.sphere(ev[0] + e[0], ev[1] + e[1], ev[2] + e[2]);
                    }
                    H[static_cast<std::size_t>(q)] = h;
                    any = any || h != 0.0;
                }
                if (!any) continue;
                // distribute to the components c of the shell: pqs <= c componentwise
                for (int c = 0; c < nc; ++c) {
                    const auto& ec = tb.comp[la][static_cast<std::size_t>(c)];
                    for (int p = 0; p <= std::min(ec[0], n); ++p)
                        for (int q = 0; q <= std::min(ec[1], n - p); ++q) {
                            const int s = n - p - q;
                            if (s > ec[2]) continue;
                            const double coef = binomial(ec[0], p) * binomial(ec[1], q) * binomial(ec[2], s) *
                                                ma[0][ec[0] - p] * ma[1][ec[1] - q] * ma[2][ec[2] - s];
                            if (coef == 0.0) continue;
                            g.T[g.at(c, m, n, lam)] += coef * H[static_cast<std::size_t>(cart_index(p, q, s))];
                        }
                }
            }
        }
    }
    for (double v : g.T) g.max_abs = std::max(g.max_abs, std::fabs(v));
    return r;
}

// Thread-local cache of the angular factors: the same (l_a, l, a) recurs for every partner shell
// and in the derivative passes (component shells).
struct AngularKey {
    int la, l;
    double a[3];
    bool operator==(const AngularKey& o) const {
        return la == o.la && l == o.l && a[0] == o.a[0] && a[1] == o.a[1] && a[2] == o.a[2];
    }
};
struct AngularKeyHash {
    std::size_t operator()(const AngularKey& k) const {
        std::size_t h = std::hash<int>()(k.la * 64 + k.l);
        for (double v : k.a) h = h * 1000003u ^ std::hash<double>()(v);
        return h;
    }
};

std::shared_ptr<const Angular> angular(int la, int l, const double a[3]) {
    thread_local std::unordered_map<AngularKey, std::shared_ptr<const Angular>, AngularKeyHash> cache;
    AngularKey k{la, l, {a[0], a[1], a[2]}};
    auto it = cache.find(k);
    if (it != cache.end()) return it->second;
    if (cache.size() > 4096) cache.clear();
    auto g = make_angular(la, l, a);
    cache.emplace(k, g);
    return g;
}

// 1D overlap recurrence S(i, j), i <= imax, j <= jmax (without the prefactor).
void overlap_1d(double p, double PA, double PB, int imax, int jmax, double* S) {
    const int nj = jmax + 1;
    const double o2p = 0.5 / p;
    S[0] = 1.0;
    for (int i = 0; i < imax; ++i) S[(i + 1) * nj] = PA * S[i * nj] + (i > 0 ? i * o2p * S[(i - 1) * nj] : 0.0);
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i <= imax; ++i) {
            double v = PB * S[i * nj + j];
            if (i > 0) v += i * o2p * S[(i - 1) * nj + j];
            if (j > 0) v += j * o2p * S[i * nj + j - 1];
            S[i * nj + j + 1] = v;
        }
}

// The local part: sum_k d_k <a| exp(-z_k |r - C|^2) |b>, three-centre overlaps.
void add_local(const Shell& A, const Shell& B, const std::array<double, 3>& C, const std::vector<EcpTerm>& terms,
               std::vector<double>& cart) {
    const Tables& tb = tables();
    const int la = A.l, lb = B.l, nb = ncart(lb), nj = lb + 1;
    std::vector<double> S[3];
    for (auto& s : S) s.resize(static_cast<std::size_t>((la + 1) * (lb + 1)));
    double AB2 = 0.0, AC2 = 0.0, BC2 = 0.0;
    for (int d = 0; d < 3; ++d) {
        AB2 += (A.center[d] - B.center[d]) * (A.center[d] - B.center[d]);
        AC2 += (A.center[d] - C[d]) * (A.center[d] - C[d]);
        BC2 += (B.center[d] - C[d]) * (B.center[d] - C[d]);
    }
    for (std::size_t i = 0; i < A.exponents.size(); ++i)
        for (std::size_t j = 0; j < B.exponents.size(); ++j)
            for (const EcpTerm& t : terms) {
                const double al = A.exponents[i], be = B.exponents[j], ze = t.exponent, p = al + be + ze;
                const double K = std::exp(-(al * be * AB2 + al * ze * AC2 + be * ze * BC2) / p);
                const double pref = A.coefficients[i] * B.coefficients[j] * t.coefficient * K * std::pow(kPi / p, 1.5);
                if (pref == 0.0) continue;
                for (int d = 0; d < 3; ++d) {
                    const double P = (al * A.center[d] + be * B.center[d] + ze * C[d]) / p;
                    overlap_1d(p, P - A.center[d], P - B.center[d], la, lb, S[d].data());
                }
                for (int ca = 0; ca < ncart(la); ++ca) {
                    const auto& ea = tb.comp[la][static_cast<std::size_t>(ca)];
                    for (int cb = 0; cb < nb; ++cb) {
                        const auto& eb = tb.comp[lb][static_cast<std::size_t>(cb)];
                        cart[static_cast<std::size_t>(ca) * nb + cb] +=
                            pref * S[0][ea[0] * nj + eb[0]] * S[1][ea[1] * nj + eb[1]] * S[2][ea[2] * nj + eb[2]];
                    }
                }
            }
}

// The semi-local part of the projector l: int r^2 dr [U_l - U_L](r) sum_m F^A_cm(r) F^B_dm(r).
void add_semilocal(const Shell& A, const Shell& B, const std::array<double, 3>& C, int l,
                   const std::vector<EcpTerm>& terms, int points, std::vector<double>& cart) {
    const int la = A.l, lb = B.l, na = ncart(la), nb = ncart(lb), nm = 2 * l + 1;
    const double a[3] = {A.center[0] - C[0], A.center[1] - C[1], A.center[2] - C[2]};
    const double b[3] = {B.center[0] - C[0], B.center[1] - C[1], B.center[2] - C[2]};
    const std::shared_ptr<const Angular> TA = angular(la, l, a), TB = angular(lb, l, b);
    if (TA->max_abs == 0.0 || TB->max_abs == 0.0) return;
    const GaussLegendre& gl = gauss_legendre(points);
    const double am = TA->amag, bm = TB->amag;
    std::vector<double> SA(static_cast<std::size_t>(na) * nm), SB(static_cast<std::size_t>(nb) * nm);
    double besA[2 * kEcpMaxShellL + 2], besB[2 * kEcpMaxShellL + 2], rp[kEcpMaxShellL + 3];
    const double four_pi_sq = 16.0 * kPi * kPi;
    for (std::size_t i = 0; i < A.exponents.size(); ++i)
        for (std::size_t j = 0; j < B.exponents.size(); ++j)
            for (const EcpTerm& t : terms) {
                const double al = A.exponents[i], be = B.exponents[j];
                const double P = al + be + t.exponent;
                const double r0 = (al * am + be * bm) / P;
                const double logE = P * r0 * r0 - al * am * am - be * bm * bm;  // <= 0
                const double cc = A.coefficients[i] * B.coefficients[j] * t.coefficient * four_pi_sq;
                const double w = kWindow / std::sqrt(P);
                const double lo = std::max(0.0, r0 - w), hi = r0 + w;
                // bound of the contribution; skip negligible products
                const double bound = std::fabs(cc) * std::exp(logE) * TA->max_abs * TB->max_abs * nm *
                                     ipow(std::max(1.0, hi), la + lb + 2) * std::sqrt(kPi / P);
                if (bound < 1e-17) continue;
                const double mid = 0.5 * (hi + lo), half = 0.5 * (hi - lo);
                const double ka = 2.0 * al * am, kb = 2.0 * be * bm;
                for (std::size_t g = 0; g < gl.x.size(); ++g) {
                    const double r = mid + half * gl.x[g];
                    const double wt = half * gl.w[g] * std::exp(logE - P * (r - r0) * (r - r0)) * r * r * cc;
                    if (wt == 0.0) continue;
                    scaled_bessel_i(TA->nlam - 1, ka * r, besA);
                    scaled_bessel_i(TB->nlam - 1, kb * r, besB);
                    rp[0] = 1.0;
                    for (int k = 1; k <= std::max(la, lb); ++k) rp[k] = rp[k - 1] * r;
                    for (int c = 0; c < na; ++c)
                        for (int m = 0; m < nm; ++m) {
                            double s = 0.0;
                            for (int n = 0; n <= la; ++n) {
                                const double* T = &TA->T[TA->at(c, m, n, 0)];
                                double u = 0.0;
                                for (int lam = (l + n) & 1; lam < TA->nlam; lam += 2) u += T[lam] * besA[lam];
                                s += rp[n] * u;
                            }
                            SA[static_cast<std::size_t>(c) * nm + m] = s;
                        }
                    for (int d = 0; d < nb; ++d)
                        for (int m = 0; m < nm; ++m) {
                            double s = 0.0;
                            for (int n = 0; n <= lb; ++n) {
                                const double* T = &TB->T[TB->at(d, m, n, 0)];
                                double u = 0.0;
                                for (int lam = (l + n) & 1; lam < TB->nlam; lam += 2) u += T[lam] * besB[lam];
                                s += rp[n] * u;
                            }
                            SB[static_cast<std::size_t>(d) * nm + m] = wt * s;
                        }
                    for (int c = 0; c < na; ++c)
                        for (int d = 0; d < nb; ++d) {
                            double s = 0.0;
                            for (int m = 0; m < nm; ++m)
                                s += SA[static_cast<std::size_t>(c) * nm + m] * SB[static_cast<std::size_t>(d) * nm + m];
                            cart[static_cast<std::size_t>(c) * nb + d] += s;
                        }
                }
            }
}

}  // namespace

void scaled_bessel_i(int nmax, double x, double* out) {
    if (nmax < 0) return;
    if (x < 0.0) throw std::invalid_argument("scaled_bessel_i: negative argument");
    if (x < 1e-300) {
        out[0] = 1.0;
        for (int n = 1; n <= nmax; ++n) out[n] = 0.0;
        return;
    }
    if (x <= 0.5) {
        // series: i_n(x) = x^n / (2n+1)!! sum_k (x^2/2)^k / (k! (2n+3)(2n+5)...(2n+2k+1))
        const double e = std::exp(-x), h = 0.5 * x * x;
        double lead = 1.0;  // x^n / (2n+1)!!
        for (int n = 0; n <= nmax; ++n) {
            if (n > 0) lead *= x / (2.0 * n + 1.0);
            double term = 1.0, sum = 1.0;
            for (int k = 1; k < 60; ++k) {
                term *= h / (k * (2.0 * n + 2.0 * k + 1.0));
                sum += term;
                if (term < 1e-17 * sum) break;
            }
            out[n] = e * lead * sum;
        }
        return;
    }
    const double e2 = std::exp(-2.0 * x);
    const double i0 = (1.0 - e2) / (2.0 * x);
    if (x > 2.0 * nmax + 10.0) {
        // upward recurrence (stable for n < x): i_(n+1) = i_(n-1) - (2n+1)/x i_n
        out[0] = i0;
        if (nmax >= 1) out[1] = ((x - 1.0) + (x + 1.0) * e2) / (2.0 * x * x);
        for (int n = 1; n < nmax; ++n) out[n + 1] = out[n - 1] - (2.0 * n + 1.0) / x * out[n];
        return;
    }
    // Miller's downward recurrence i_(n-1) = i_(n+1) + (2n+1)/x i_n from n = start (f = 0 above),
    // normalised with i_0 (the values grow by less than e^80 over the range used here)
    const int start = nmax + 40 + static_cast<int>(x);
    double vals[2 * kEcpMaxShellL + 2] = {};
    double above = 0.0, f = 1e-30;  // f_(n+1), f_n
    for (int n = start; n >= 1; --n) {
        const double below = above + (2.0 * n + 1.0) / x * f;  // f_(n-1)
        above = f;
        f = below;
        if (n - 1 <= nmax) vals[n - 1] = f;
    }
    const double scale = i0 / vals[0];
    for (int n = 0; n <= nmax; ++n) out[n] = vals[n] * scale;
}

void add_ecp_cartesian(const Shell& A, const Shell& B, const Atom& center, std::vector<double>& cart, int points) {
    if (!center.ecp) return;
    if (A.l > kEcpMaxShellL || B.l > kEcpMaxShellL)
        throw std::invalid_argument("ECP integrals: angular momentum beyond the supported maximum");
    const AtomEcp& U = *center.ecp;
    if (static_cast<int>(U.projected.size()) > kMaxProjector)
        throw std::invalid_argument("ECP integrals: projector angular momentum beyond the supported maximum");
    const std::size_t need = static_cast<std::size_t>(ncart(A.l)) * ncart(B.l);
    if (cart.size() < need) cart.resize(need, 0.0);
    if (!U.local.empty()) add_local(A, B, center.position, U.local, cart);
    for (std::size_t l = 0; l < U.projected.size(); ++l)
        if (!U.projected[l].empty())
            add_semilocal(A, B, center.position, static_cast<int>(l), U.projected[l], points > 0 ? points : kDefaultPoints,
                          cart);
}

}  // namespace macqc
