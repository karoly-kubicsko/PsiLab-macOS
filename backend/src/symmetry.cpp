#include "macqc/symmetry.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <map>

#include "macqc/constants.h"

namespace macqc {

namespace {

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<std::array<double, 3>, 3>;
constexpr double kPi = 3.14159265358979323846;

double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 rotate_vector(const Mat3& R, const Vec3& v) {
    return {R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2], R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
            R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2]};
}
Mat3 multiply(const Mat3& A, const Mat3& B) {
    Mat3 C{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
    return C;
}
// rotation by angle t about the unit vector a (Rodrigues)
Mat3 rotation(const Vec3& a, double t) {
    const double c = std::cos(t), s = std::sin(t), u = 1.0 - c;
    return {{{c + a[0] * a[0] * u, a[0] * a[1] * u - a[2] * s, a[0] * a[2] * u + a[1] * s},
             {a[1] * a[0] * u + a[2] * s, c + a[1] * a[1] * u, a[1] * a[2] * u - a[0] * s},
             {a[2] * a[0] * u - a[1] * s, a[2] * a[1] * u + a[0] * s, c + a[2] * a[2] * u}}};
}
Mat3 reflection(const Vec3& n) {
    Mat3 R{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) R[i][j] = (i == j ? 1.0 : 0.0) - 2.0 * n[i] * n[j];
    return R;
}

// Character tables: classes in order, their sizes, irreps with characters per class.
struct Table {
    std::vector<std::string> classes;
    std::vector<int> sizes;
    std::vector<std::pair<std::string, std::vector<double>>> irreps;
};

const std::map<std::string, Table>& tables() {
    static const std::map<std::string, Table> t = {
        {"C1", {{"E"}, {1}, {{"A", {1}}}}},
        {"Cs", {{"E", "sigma_h"}, {1, 1}, {{"A'", {1, 1}}, {"A''", {1, -1}}}}},
        {"Ci", {{"E", "i"}, {1, 1}, {{"Ag", {1, 1}}, {"Au", {1, -1}}}}},
        {"C2", {{"E", "C2"}, {1, 1}, {{"A", {1, 1}}, {"B", {1, -1}}}}},
        {"C2v",
         {{"E", "C2", "sigma_v", "sigma_v'"},
          {1, 1, 1, 1},
          {{"A1", {1, 1, 1, 1}}, {"A2", {1, 1, -1, -1}}, {"B1", {1, -1, 1, -1}}, {"B2", {1, -1, -1, 1}}}}},
        {"C2h",
         {{"E", "C2", "i", "sigma_h"},
          {1, 1, 1, 1},
          {{"Ag", {1, 1, 1, 1}}, {"Bg", {1, -1, 1, -1}}, {"Au", {1, 1, -1, -1}}, {"Bu", {1, -1, -1, 1}}}}},
        {"C3", {{"E", "C3"}, {1, 2}, {{"A", {1, 1}}, {"E", {2, -1}}}}},
        {"C3v", {{"E", "C3", "sigma_v"}, {1, 2, 3}, {{"A1", {1, 1, 1}}, {"A2", {1, 1, -1}}, {"E", {2, -1, 0}}}}},
        {"D3h",
         {{"E", "C3", "C2", "sigma_h", "S3", "sigma_v"},
          {1, 2, 3, 1, 2, 3},
          {{"A1'", {1, 1, 1, 1, 1, 1}},
           {"A2'", {1, 1, -1, 1, 1, -1}},
           {"E'", {2, -1, 0, 2, -1, 0}},
           {"A1''", {1, 1, 1, -1, -1, -1}},
           {"A2''", {1, 1, -1, -1, -1, 1}},
           {"E''", {2, -1, 0, -2, 1, 0}}}}},
        // 9.0: further groups of custom molecules
        {"D2",
         {{"E", "C2(z)", "C2(y)", "C2(x)"},
          {1, 1, 1, 1},
          {{"A", {1, 1, 1, 1}}, {"B1", {1, 1, -1, -1}}, {"B2", {1, -1, 1, -1}}, {"B3", {1, -1, -1, 1}}}}},
        {"D2h",
         {{"E", "C2(z)", "C2(y)", "C2(x)", "i", "sigma(xy)", "sigma(xz)", "sigma(yz)"},
          {1, 1, 1, 1, 1, 1, 1, 1},
          {{"Ag", {1, 1, 1, 1, 1, 1, 1, 1}},
           {"B1g", {1, 1, -1, -1, 1, 1, -1, -1}},
           {"B2g", {1, -1, 1, -1, 1, -1, 1, -1}},
           {"B3g", {1, -1, -1, 1, 1, -1, -1, 1}},
           {"Au", {1, 1, 1, 1, -1, -1, -1, -1}},
           {"B1u", {1, 1, -1, -1, -1, -1, 1, 1}},
           {"B2u", {1, -1, 1, -1, -1, 1, -1, 1}},
           {"B3u", {1, -1, -1, 1, -1, 1, 1, -1}}}}},
        {"D3", {{"E", "C3", "C2"}, {1, 2, 3}, {{"A1", {1, 1, 1}}, {"A2", {1, 1, -1}}, {"E", {2, -1, 0}}}}},
        {"D3d",
         {{"E", "C3", "C2", "i", "S6", "sigma_d"},
          {1, 2, 3, 1, 2, 3},
          {{"A1g", {1, 1, 1, 1, 1, 1}},
           {"A2g", {1, 1, -1, 1, 1, -1}},
           {"Eg", {2, -1, 0, 2, -1, 0}},
           {"A1u", {1, 1, 1, -1, -1, -1}},
           {"A2u", {1, 1, -1, -1, -1, 1}},
           {"Eu", {2, -1, 0, -2, 1, 0}}}}},
        {"Td",
         {{"E", "C3", "C2", "S4", "sigma_d"},
          {1, 8, 3, 6, 6},
          {{"A1", {1, 1, 1, 1, 1}},
           {"A2", {1, 1, 1, -1, -1}},
           {"E", {2, -1, 2, 0, 0}},
           {"T1", {3, 0, -1, 1, -1}},
           {"T2", {3, 0, -1, -1, 1}}}}},
    };
    return t;
}

}  // namespace

namespace {

// Axis of a proper rotation (eigenvector with eigenvalue 1) or of an improper one (eigenvalue
// -1); the zero vector for E and i.
Vec3 operation_axis(const Mat3& R, bool proper) {
    // R' = R (proper) or -R (improper; -sigma_h C(theta) is the rotation C(theta + pi) about the
    // same axis) is a rotation by theta' about a: R' + R'^T + (1 - tr R') I = 2 (1 - cos theta') a a^T
    const double s = proper ? 1.0 : -1.0;
    Vec3 best{0.0, 0.0, 0.0};
    Mat3 Rs{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) Rs[i][j] = s * R[i][j];
    const double tr = Rs[0][0] + Rs[1][1] + Rs[2][2];
    for (int k = 0; k < 3; ++k) {
        const Vec3 col{Rs[0][k] + Rs[k][0] + (k == 0 ? 1.0 - tr : 0.0), Rs[1][k] + Rs[k][1] + (k == 1 ? 1.0 - tr : 0.0),
                       Rs[2][k] + Rs[k][2] + (k == 2 ? 1.0 - tr : 0.0)};
        if (norm(col) > norm(best)) best = col;
    }
    const double l = norm(best);
    if (l < 1.0e-8) return {0.0, 0.0, 0.0};
    for (double& c : best) c /= l;
    return best;
}

int expected_order(const std::string& name) {
    if (name == "C1") return 1;
    if (name == "Cs" || name == "Ci") return 2;
    if (name == "T") return 12;
    if (name == "Td" || name == "Th" || name == "O") return 24;
    if (name == "Oh") return 48;
    if (name == "I") return 60;
    if (name == "Ih") return 120;
    // Cn, Cnv, Cnh, Sn, Dn, Dnh, Dnd
    std::size_t i = 1;
    int n = 0;
    while (i < name.size() && std::isdigit(static_cast<unsigned char>(name[i]))) n = 10 * n + (name[i++] - '0');
    const std::string suffix = name.substr(i);
    if (name[0] == 'C') return suffix.empty() ? n : 2 * n;
    if (name[0] == 'S') return n;
    if (name[0] == 'D') return suffix.empty() ? 2 * n : 4 * n;
    return 0;
}

}  // namespace

PointGroup point_group(const Molecule& mol, double tol) {
    if (tol <= 0.0) {
        PointGroup g;
        for (double t : {2.0e-3, 5.0e-4, 1.0e-4, 2.0e-5}) {
            g = point_group(mol, t);
            if (g.consistent) return g;
        }
        return g;
    }
    PointGroup g;
    const std::size_t n = mol.atoms.size();
    // centre of mass
    Vec3 com{0.0, 0.0, 0.0};
    double M = 0.0;
    for (const Atom& a : mol.atoms) {
        const double m = isotope_mass(a.Z);
        for (int k = 0; k < 3; ++k) com[static_cast<std::size_t>(k)] += m * a.position[static_cast<std::size_t>(k)];
        M += m;
    }
    for (double& c : com) c /= M;
    std::vector<Vec3> x(n);
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t k = 0; k < 3; ++k) x[a][k] = mol.atoms[a].position[k] - com[k];
    const Mat3 I{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    const Mat3 inversion{{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}};
    auto test = [&](const Mat3& R, const std::string& kind) {
        std::vector<int> perm(n, -1);
        for (std::size_t a = 0; a < n; ++a) {
            const Vec3 y = rotate_vector(R, x[a]);
            for (std::size_t b = 0; b < n; ++b) {
                if (mol.atoms[b].Z != mol.atoms[a].Z) continue;
                const Vec3 d{y[0] - x[b][0], y[1] - x[b][1], y[2] - x[b][2]};
                if (norm(d) < tol) {
                    perm[a] = static_cast<int>(b);
                    break;
                }
            }
            if (perm[a] < 0) return false;
        }
        // the same kind of operation with the same permutation of the atoms is the same symmetry
        // element (candidate axes that differ by the residual asymmetry give nearly equal matrices)
        for (const auto& op : g.operations)
            if (op.kind == kind && op.permutation == perm) return true;
        g.operations.push_back({kind, "", R, perm});
        return true;
    };

    // 9.0: linear molecules, C_infinity_v and D_infinity_h (the rotations about the axis are not
    // listed; E and, for D_infinity_h, the inversion are, for the g/u species of the modes)
    {
        std::size_t a0 = 0, a1 = 0;
        double far = -1.0;
        for (std::size_t a = 0; a < n; ++a)
            for (std::size_t b = a + 1; b < n; ++b) {
                const Vec3 d{x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]};
                if (norm(d) > far) {
                    far = norm(d);
                    a0 = a;
                    a1 = b;
                }
            }
        if (n >= 2) {
            Vec3 e{x[a1][0] - x[a0][0], x[a1][1] - x[a0][1], x[a1][2] - x[a0][2]};
            for (double& c : e) c /= far;
            bool linear = true;
            for (std::size_t a = 0; a < n && linear; ++a) {
                const Vec3 d{x[a][0] - x[a0][0], x[a][1] - x[a0][1], x[a][2] - x[a0][2]};
                const double along = dot(d, e);
                const Vec3 perp{d[0] - along * e[0], d[1] - along * e[1], d[2] - along * e[2]};
                if (norm(perp) > tol) linear = false;
            }
            if (linear) {
                test(I, "E");
                const bool centro = test(inversion, "i");
                g.name = centro ? "D\xE2\x88\x9Eh" : "C\xE2\x88\x9Ev";
                g.linear = true;
                g.axis = e;
                g.symmetry_number = centro ? 2 : 1;
                g.has_table = false;
                g.consistent = true;
                for (auto& op : g.operations) op.symmetry_class = op.kind;
                return g;
            }
        }
    }

    // candidate axes / plane normals: Cartesian and principal axes, atoms, and for pairs of
    // symmetry-equivalent atoms (same element, same distance from the centre of mass) their
    // midpoint, difference and cross product; normals of triangles of equivalent atoms (C_n axes
    // without atoms on them)
    std::vector<Vec3> cand;
    std::vector<std::array<long long, 3>> keys;
    auto add = [&](Vec3 v) {
        const double l = norm(v);
        if (l < 1.0e-6) return;
        for (double& c : v) c /= l;
        // canonical hemisphere, then a coarse key for duplicates
        const double s = (v[2] > 1.0e-9 || (std::fabs(v[2]) <= 1.0e-9 && (v[1] > 1.0e-9 || (std::fabs(v[1]) <= 1.0e-9 && v[0] > 0.0)))) ? 1.0 : -1.0;
        for (double& c : v) c *= s;
        const std::array<long long, 3> key{std::llround(v[0] * 1.0e7), std::llround(v[1] * 1.0e7), std::llround(v[2] * 1.0e7)};
        if (std::find(keys.begin(), keys.end(), key) != keys.end()) return;
        for (const auto& w : cand)
            if (std::fabs(std::fabs(dot(v, w)) - 1.0) < 1.0e-10) return;
        keys.push_back(key);
        cand.push_back(v);
    };
    for (int k = 0; k < 3; ++k) {
        Vec3 e{0.0, 0.0, 0.0};
        e[static_cast<std::size_t>(k)] = 1.0;
        add(e);
    }
    {
        // principal axes of the (mass-weighted) inertia tensor, by Jacobi rotations
        double T[3][3] = {};
        for (std::size_t a = 0; a < n; ++a) {
            const double m = isotope_mass(mol.atoms[a].Z), r2 = dot(x[a], x[a]);
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) T[i][j] += m * ((i == j ? r2 : 0.0) - x[a][static_cast<std::size_t>(i)] * x[a][static_cast<std::size_t>(j)]);
        }
        double V[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        for (int sweep = 0; sweep < 50; ++sweep) {
            double off = std::fabs(T[0][1]) + std::fabs(T[0][2]) + std::fabs(T[1][2]);
            if (off < 1.0e-14) break;
            for (int p = 0; p < 2; ++p)
                for (int q = p + 1; q < 3; ++q) {
                    if (std::fabs(T[p][q]) < 1.0e-300) continue;
                    const double th = 0.5 * std::atan2(2.0 * T[p][q], T[q][q] - T[p][p]);
                    const double c = std::cos(th), s = std::sin(th);
                    for (int k = 0; k < 3; ++k) {
                        const double tkp = T[k][p], tkq = T[k][q];
                        T[k][p] = c * tkp - s * tkq;
                        T[k][q] = s * tkp + c * tkq;
                    }
                    for (int k = 0; k < 3; ++k) {
                        const double tpk = T[p][k], tqk = T[q][k];
                        T[p][k] = c * tpk - s * tqk;
                        T[q][k] = s * tpk + c * tqk;
                    }
                    for (int k = 0; k < 3; ++k) {
                        const double vkp = V[k][p], vkq = V[k][q];
                        V[k][p] = c * vkp - s * vkq;
                        V[k][q] = s * vkp + c * vkq;
                    }
                }
        }
        for (int k = 0; k < 3; ++k) add({V[0][k], V[1][k], V[2][k]});
    }
    auto equivalent = [&](std::size_t a, std::size_t b) {
        return mol.atoms[a].Z == mol.atoms[b].Z && std::fabs(norm(x[a]) - norm(x[b])) < 2.0 * tol;
    };
    for (std::size_t a = 0; a < n; ++a) {
        add(x[a]);
        std::vector<std::size_t> set;
        for (std::size_t b = a + 1; b < n; ++b) {
            if (!equivalent(a, b)) continue;
            set.push_back(b);
            add({x[a][0] + x[b][0], x[a][1] + x[b][1], x[a][2] + x[b][2]});
            add({x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]});
            add(cross(x[a], x[b]));
        }
        if (set.size() >= 2 && set.size() <= 12)
            for (std::size_t i = 0; i < set.size(); ++i)
                for (std::size_t j = i + 1; j < set.size(); ++j) {
                    const std::size_t b = set[i], c = set[j];
                    add(cross({x[b][0] - x[a][0], x[b][1] - x[a][1], x[b][2] - x[a][2]},
                              {x[c][0] - x[a][0], x[c][1] - x[a][1], x[c][2] - x[a][2]}));
                }
    }
    test(I, "E");
    test(inversion, "i");
    constexpr int kMaxN = 8;  // C2 ... C8 and S3 ... S8 (and their powers)
    for (const auto& a : cand) {
        for (int nn = 2; nn <= kMaxN; ++nn) {
            // C_n^1 first; its powers only if it is a symmetry operation (C_n^k with gcd(n, k) = 1
            // is one exactly when C_n is)
            if (!test(rotation(a, 2.0 * kPi / nn), "C" + std::to_string(nn))) continue;
            for (int k = 2; k < nn; ++k) {
                int gg = nn, kk = k;
                while (kk) {
                    const int r = gg % kk;
                    gg = kk;
                    kk = r;
                }
                const int m = nn / gg;
                test(rotation(a, 2.0 * kPi * k / nn), "C" + std::to_string(m));
            }
        }
        test(reflection(a), "sigma");
        for (int nn = 3; nn <= kMaxN; ++nn) {
            if (!test(multiply(reflection(a), rotation(a, 2.0 * kPi / nn)), "S" + std::to_string(nn))) continue;
            for (int k = 2; k < nn; ++k) {
                int gg = nn, kk = k;
                while (kk) {
                    const int r = gg % kk;
                    gg = kk;
                    kk = r;
                }
                const int m = nn / gg;
                const Mat3 S = multiply(reflection(a), rotation(a, 2.0 * kPi * k / nn));
                if (m == 1) continue;  // sigma_h, found above
                if (m == 2) {
                    test(S, "i");
                    continue;
                }
                test(S, "S" + std::to_string(m));
            }
        }
    }
    // the group from its elements
    int proper = 0, nsigma = 0, ninv = 0, nmax = 1;
    std::map<std::string, int> count;
    for (const auto& op : g.operations) {
        ++count[op.kind];
        const bool improper = op.kind == "sigma" || op.kind == "i" || op.kind[0] == 'S';
        if (!improper) {
            ++proper;
            if (op.kind != "E") nmax = std::max(nmax, std::stoi(op.kind.substr(1)));
        }
        if (op.kind == "sigma") ++nsigma;
        if (op.kind == "i") ++ninv;
    }
    g.symmetry_number = proper;
    const std::size_t order = g.operations.size();
    auto axes_of = [&](const std::string& kind) {
        std::vector<Vec3> axes;
        for (const auto& op : g.operations) {
            if (op.kind != kind) continue;
            const Vec3 ax = operation_axis(op.R, kind[0] == 'C');
            bool seen = false;
            for (const auto& b : axes)
                if (std::fabs(std::fabs(dot(ax, b)) - 1.0) < 1.0e-3) seen = true;
            if (!seen) axes.push_back(ax);
        }
        return axes;
    };
    std::string name;
    const int nC3axes = static_cast<int>(axes_of("C3").size());
    if (nC3axes > 1) {
        // cubic and icosahedral groups
        if (count["C5"] > 0) name = ninv ? "Ih" : "I";
        else if (count["C4"] > 0) name = ninv ? "Oh" : "O";
        else if (ninv) name = "Th";
        else if (nsigma) name = "Td";
        else name = "T";
    } else if (nmax >= 2) {
        const std::vector<Vec3> principal_axes = axes_of("C" + std::to_string(nmax));
        Vec3 z = principal_axes.empty() ? Vec3{0.0, 0.0, 1.0} : principal_axes[0];
        if (nmax == 2 && principal_axes.size() == 3) {
            // D2, D2h, D2d share three C2 axes: for D2d the principal axis is the S4 axis
            const std::vector<Vec3> s4 = axes_of("S4");
            if (!s4.empty()) z = s4[0];
        }
        int perp_c2 = 0;
        for (const auto& ax : axes_of("C2"))
            if (std::fabs(dot(ax, z)) < 1.0e-3) ++perp_c2;
        bool sigma_h = false;
        int sigma_v = 0;
        for (const auto& op : g.operations) {
            if (op.kind != "sigma") continue;
            const Vec3 nrm = operation_axis(op.R, false);
            if (std::fabs(std::fabs(dot(nrm, z)) - 1.0) < 1.0e-3) sigma_h = true;
            else if (std::fabs(dot(nrm, z)) < 1.0e-3) ++sigma_v;
        }
        const std::string N = std::to_string(nmax);
        if (perp_c2 >= nmax) name = sigma_h ? "D" + N + "h" : (sigma_v >= nmax ? "D" + N + "d" : "D" + N);
        else if (sigma_h) name = "C" + N + "h";
        else if (sigma_v >= nmax) name = "C" + N + "v";
        else if (count["S" + std::to_string(2 * nmax)] > 0) name = "S" + std::to_string(2 * nmax);
        else name = "C" + N;
    } else {
        name = nsigma ? "Cs" : (ninv ? "Ci" : "C1");
    }
    g.name = name;
    g.consistent = expected_order(name) == static_cast<int>(order);
    g.has_table = g.consistent && tables().count(name) > 0;
    if (!g.consistent) g.name = name + " (order " + std::to_string(order) + ")";
    // classes
    Vec3 c3axis{0.0, 0.0, 0.0};
    for (const auto& op : g.operations)
        if (op.kind == "C3") {
            c3axis = operation_axis(op.R, true);
            break;
        }
    // C2v: sigma_v'(yz) is the plane with more atoms in it (the molecular plane of a planar molecule)
    int best_atoms = -1, best = -1;
    for (std::size_t o = 0; o < g.operations.size(); ++o) {
        if (g.operations[o].kind != "sigma") continue;
        int in = 0;
        for (std::size_t a = 0; a < n; ++a)
            if (g.operations[o].permutation[a] == static_cast<int>(a)) ++in;
        if (in > best_atoms) {
            best_atoms = in;
            best = static_cast<int>(o);
        }
    }
    // D2, D2h (9.0): the axes x, y, z (Mulliken, J. Chem. Phys. 23, 1997 (1955)): a planar molecule
    // has x perpendicular to its plane; z is the C2 axis through the most atoms (then y)
    std::vector<Vec3> d2axes;
    if (g.name == "D2" || g.name == "D2h") {
        std::vector<Vec3> c2 = axes_of("C2");
        auto on_axis = [&](const Vec3& ax) {
            int k = 0;
            for (std::size_t a = 0; a < n; ++a) {
                const double along = dot(x[a], ax);
                const Vec3 perp{x[a][0] - along * ax[0], x[a][1] - along * ax[1], x[a][2] - along * ax[2]};
                if (norm(perp) < tol) ++k;
            }
            return k;
        };
        auto in_plane = [&](const Vec3& nrm) {
            int k = 0;
            for (std::size_t a = 0; a < n; ++a)
                if (std::fabs(dot(x[a], nrm)) < tol) ++k;
            return k;
        };
        // x: the normal of the plane containing all atoms, if any
        int xi = -1;
        for (std::size_t k = 0; k < c2.size(); ++k)
            if (in_plane(c2[k]) == static_cast<int>(n)) xi = static_cast<int>(k);
        std::vector<int> rest;
        for (int k = 0; k < static_cast<int>(c2.size()); ++k)
            if (k != xi) rest.push_back(k);
        std::stable_sort(rest.begin(), rest.end(), [&](int a, int b) { return on_axis(c2[static_cast<std::size_t>(a)]) > on_axis(c2[static_cast<std::size_t>(b)]); });
        if (c2.size() == 3) {
            // order z, y, x
            if (xi >= 0) d2axes = {c2[static_cast<std::size_t>(rest[0])], c2[static_cast<std::size_t>(rest[1])], c2[static_cast<std::size_t>(xi)]};
            else d2axes = {c2[static_cast<std::size_t>(rest[0])], c2[static_cast<std::size_t>(rest[1])], c2[static_cast<std::size_t>(rest[2])]};
        }
    }
    const char* const d2names[3] = {"C2(z)", "C2(y)", "C2(x)"};
    const char* const d2planes[3] = {"sigma(xy)", "sigma(xz)", "sigma(yz)"};
    for (std::size_t o = 0; o < g.operations.size(); ++o) {
        SymmetryOperation& op = g.operations[o];
        if (op.kind == "sigma") {
            if (g.name == "Cs") op.symmetry_class = "sigma_h";
            else if (g.name == "C2h") op.symmetry_class = "sigma_h";
            else if (g.name == "C3v") op.symmetry_class = "sigma_v";
            else if (g.name == "C2v") op.symmetry_class = static_cast<int>(o) == best ? "sigma_v'" : "sigma_v";
            else if (g.name == "D3h") {
                const Vec3 nn = operation_axis(op.R, false);
                op.symmetry_class = std::fabs(dot(nn, c3axis)) > 0.9 ? "sigma_h" : "sigma_v";
            } else if (g.name == "D3d" || g.name == "Td") {
                op.symmetry_class = "sigma_d";
            } else if (g.name == "D2h" && d2axes.size() == 3) {
                const Vec3 nn = operation_axis(op.R, false);
                for (std::size_t k = 0; k < 3; ++k)
                    if (std::fabs(dot(nn, d2axes[k])) > 0.9) op.symmetry_class = d2planes[k];
            } else {
                op.symmetry_class = op.kind;
            }
        } else if (op.kind == "C2" && (g.name == "D2" || g.name == "D2h") && d2axes.size() == 3) {
            const Vec3 ax = operation_axis(op.R, true);
            for (std::size_t k = 0; k < 3; ++k)
                if (std::fabs(dot(ax, d2axes[k])) > 0.9) op.symmetry_class = d2names[k];
        } else {
            op.symmetry_class = op.kind;
        }
    }
    return g;
}

bool has_degenerate_species(const PointGroup& g) {
    if (g.linear) return true;
    const auto it = tables().find(g.name);
    if (it == tables().end()) return false;
    for (const auto& ir : it->second.irreps)
        if (ir.second[0] > 1.5) return true;
    return false;
}

std::string mode_species(const PointGroup& g, const std::vector<std::vector<double>>& modes) {
    if (g.linear) {
        // 9.0: Sigma (displacements along the axis) or Pi (perpendicular); g/u from the inversion
        double along = 0.0, total = 0.0, chi_i = 0.0;
        const SymmetryOperation* inv = nullptr;
        for (const auto& op : g.operations)
            if (op.kind == "i") inv = &op;
        for (const auto& l : modes) {
            const std::size_t nat = l.size() / 3;
            for (std::size_t a = 0; a < nat; ++a) {
                const Vec3 v{l[3 * a], l[3 * a + 1], l[3 * a + 2]};
                along += dot(v, g.axis) * dot(v, g.axis);
                total += dot(v, v);
                if (inv) {
                    const std::size_t b = static_cast<std::size_t>(inv->permutation[a]);
                    chi_i -= v[0] * l[3 * b] + v[1] * l[3 * b + 1] + v[2] * l[3 * b + 2];
                }
            }
        }
        if (total <= 0.0) return "";
        const double f = along / total;
        std::string base;
        if (f > 0.99 && modes.size() == 1) base = "Sigma";
        else if (f < 0.01 && modes.size() == 2) base = "Pi";
        else return "";
        if (!inv) return base == "Sigma" ? "Sigma+" : "Pi";
        const std::string gu = chi_i > 0.5 ? "g" : (chi_i < -0.5 ? "u" : "");
        if (gu.empty()) return "";
        return base == "Sigma" ? "Sigma_" + gu + "+" : "Pi_" + gu;
    }
    if (!g.has_table) return "";
    const auto it = tables().find(g.name);
    if (it == tables().end()) return "";
    const Table& t = it->second;
    const double h = static_cast<double>(g.operations.size());
    // character of each operation in the space of the modes
    std::vector<double> chi(g.operations.size(), 0.0);
    for (std::size_t o = 0; o < g.operations.size(); ++o) {
        const SymmetryOperation& op = g.operations[o];
        double c = 0.0;
        for (const auto& l : modes) {
            const std::size_t nat = l.size() / 3;
            for (std::size_t a = 0; a < nat; ++a) {
                const Vec3 v{l[3 * a], l[3 * a + 1], l[3 * a + 2]};
                const Vec3 rv = rotate_vector(op.R, v);
                const std::size_t b = static_cast<std::size_t>(op.permutation[a]);
                c += rv[0] * l[3 * b] + rv[1] * l[3 * b + 1] + rv[2] * l[3 * b + 2];
            }
        }
        chi[o] = c;
    }
    std::string label;
    double total_dim = 0.0;
    for (const auto& ir : t.irreps) {
        double s = 0.0;
        for (std::size_t o = 0; o < g.operations.size(); ++o) {
            const auto ci = std::find(t.classes.begin(), t.classes.end(), g.operations[o].symmetry_class);
            if (ci == t.classes.end()) return "";
            s += chi[o] * ir.second[static_cast<std::size_t>(ci - t.classes.begin())];
        }
        const double m = s / h;
        const long k = std::lround(m);
        if (std::fabs(m - static_cast<double>(k)) > 0.2) return "";  // not a representation: numerical noise
        for (long r = 0; r < k; ++r) label += (label.empty() ? "" : "+") + ir.first;
        total_dim += static_cast<double>(k) * ir.second[0];
    }
    if (std::fabs(total_dim - static_cast<double>(modes.size())) > 0.5) return "";
    return label;
}

}  // namespace macqc
