#include "macqc/symmetry.h"

#include <algorithm>
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
    };
    return t;
}

}  // namespace

PointGroup point_group(const Molecule& mol, double tol) {
    if (tol <= 0.0) {
        PointGroup g;
        for (double t : {2.0e-3, 5.0e-4, 1.0e-4, 2.0e-5}) {
            g = point_group(mol, t);
            if (g.has_table) return g;
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
    // candidate axes / plane normals: atoms, midpoints and normals of atom pairs, Cartesian axes
    std::vector<Vec3> cand;
    auto add = [&](Vec3 v) {
        const double l = norm(v);
        if (l < 1.0e-6) return;
        for (double& c : v) c /= l;
        for (const auto& w : cand)
            if (std::fabs(std::fabs(dot(v, w)) - 1.0) < 1.0e-10) return;
        cand.push_back(v);
    };
    for (int k = 0; k < 3; ++k) {
        Vec3 e{0.0, 0.0, 0.0};
        e[static_cast<std::size_t>(k)] = 1.0;
        add(e);
    }
    for (std::size_t a = 0; a < n; ++a) {
        add(x[a]);
        for (std::size_t b = a + 1; b < n; ++b) {
            add({x[a][0] + x[b][0], x[a][1] + x[b][1], x[a][2] + x[b][2]});
            add({x[a][0] - x[b][0], x[a][1] - x[b][1], x[a][2] - x[b][2]});
            add(cross(x[a], x[b]));
        }
    }
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
            if (perm[a] < 0) return;
        }
        // the same kind of operation with the same permutation of the atoms is the same symmetry
        // element (candidate axes that differ by the residual asymmetry give nearly equal matrices)
        for (const auto& op : g.operations)
            if (op.kind == kind && op.permutation == perm) return;
        g.operations.push_back({kind, "", R, perm});
    };
    const Mat3 I{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    test(I, "E");
    test({{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}, "i");
    for (const auto& a : cand) {
        test(rotation(a, kPi), "C2");
        test(rotation(a, 2.0 * kPi / 3.0), "C3");
        test(rotation(a, 4.0 * kPi / 3.0), "C3");
        test(reflection(a), "sigma");
        test(multiply(reflection(a), rotation(a, 2.0 * kPi / 3.0)), "S3");
        test(multiply(reflection(a), rotation(a, 4.0 * kPi / 3.0)), "S3");
    }
    int nC2 = 0, nC3 = 0, ns = 0, ni = 0, nS3 = 0;
    for (const auto& op : g.operations) {
        if (op.kind == "C2") ++nC2;
        if (op.kind == "C3") ++nC3;
        if (op.kind == "sigma") ++ns;
        if (op.kind == "i") ++ni;
        if (op.kind == "S3") ++nS3;
    }
    g.symmetry_number = 1 + nC2 + nC3;
    const std::size_t order = g.operations.size();
    if (nC3 == 2 && nC2 == 3 && ns == 4 && nS3 == 2 && order == 12) g.name = "D3h";
    else if (nC3 == 2 && ns == 3 && order == 6) g.name = "C3v";
    else if (nC3 == 2 && order == 3) g.name = "C3";
    else if (nC2 == 1 && ns == 2 && order == 4) g.name = "C2v";
    else if (nC2 == 1 && ns == 1 && ni == 1 && order == 4) g.name = "C2h";
    else if (nC2 == 1 && order == 2) g.name = "C2";
    else if (ns == 1 && order == 2) g.name = "Cs";
    else if (ni == 1 && order == 2) g.name = "Ci";
    else if (order == 1) g.name = "C1";
    else {
        g.name = "other (order " + std::to_string(order) + ")";
        g.has_table = false;
    }
    // classes
    Vec3 c3axis{0.0, 0.0, 0.0};
    for (const auto& op : g.operations)
        if (op.kind == "C3") {
            // axis: eigenvector with eigenvalue 1, from R - R^T
            c3axis = {op.R[2][1] - op.R[1][2], op.R[0][2] - op.R[2][0], op.R[1][0] - op.R[0][1]};
            const double l = norm(c3axis);
            for (double& c : c3axis) c /= l;
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
    for (std::size_t o = 0; o < g.operations.size(); ++o) {
        SymmetryOperation& op = g.operations[o];
        if (op.kind == "sigma") {
            if (g.name == "Cs") op.symmetry_class = "sigma_h";
            else if (g.name == "C2h") op.symmetry_class = "sigma_h";
            else if (g.name == "C3v") op.symmetry_class = "sigma_v";
            else if (g.name == "C2v") op.symmetry_class = static_cast<int>(o) == best ? "sigma_v'" : "sigma_v";
            else if (g.name == "D3h") {
                // normal of the plane: eigenvector with eigenvalue -1, I - R = 2 n n^T
                Vec3 nn{0.0, 0.0, 0.0};
                for (int k = 0; k < 3; ++k) {
                    const Vec3 col{(k == 0 ? 1.0 : 0.0) - op.R[0][k], (k == 1 ? 1.0 : 0.0) - op.R[1][k],
                                   (k == 2 ? 1.0 : 0.0) - op.R[2][k]};
                    if (norm(col) > norm(nn)) nn = col;
                }
                const double l = norm(nn);
                for (double& c : nn) c /= l;
                op.symmetry_class = std::fabs(dot(nn, c3axis)) > 0.9 ? "sigma_h" : "sigma_v";
            }
        } else {
            op.symmetry_class = op.kind;
        }
    }
    return g;
}

std::string mode_species(const PointGroup& g, const std::vector<std::vector<double>>& modes) {
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
