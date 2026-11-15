// Tests of the molecular orbitals of ΨLab 7.0, stage 7.5 (CTest):
//  - basis functions against PySCF (backend/tests/data/orbitals_reference.json from
//    reference/pyscf_orbitals_reference.py): for every contracted shell of O and Li (aug-cc-pV5Z),
//    N (aug-cc-pCV5Z), He and H (aug-cc-pV5Z), l = 0 ... 5, the kernel sum_m phi_lm(r) phi_lm(r'),
//    which does not depend on how the real harmonics are ordered or signed;
//  - RHF canonical orbitals of H2O, NH3 and H2O2 (aug-cc-pVDZ) at points around the molecules:
//    energies, every occupied orbital and the lowest four virtual ones (up to the sign; the
//    degenerate e orbitals of NH3 through the sum of their squares) and the electron density;
//  - Full CI natural orbitals of LiH (aug-cc-pVDZ): the occupations, the Full CI density
//    sum_k n_k chi_k(r)^2 and the orthonormality of the natural orbitals (C^T S C = 1);
//  - the two Szabo-Ostlund orbitals of HeH+ (STO-3G): normalised and orthogonal on a grid;
//  - isosurfaces (marching tetrahedra): a Gaussian sphere (radius, normals, a closed and
//    consistently oriented surface, enclosed volume), the lobes of the HOMO of H2O (|psi| at the
//    vertices, the two lobes mirror images of each other);
//  - the "orbitals" object of the output (energy, Full CI, a job, Szabo-Ostlund, estimate) and the
//    C API (macqc_orbital_surface, macqc_orbital_plane, errors).
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_integrals.h"
#include "macqc/macqc_api.h"
#include "macqc/molecular.h"
#include "macqc/molecule.h"
#include "macqc/orbitals.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-84s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
}

#define CHECK(cond)                                                     \
    do {                                                                \
        ++g_checks;                                                     \
        if (!(cond)) {                                                  \
            ++g_failures;                                               \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        }                                                               \
    } while (0)

using Clock = std::chrono::steady_clock;
static double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static std::string basis_dir() {
    const char* d = std::getenv("PSILAB_BASIS_DIR");
    return d ? std::string(d) : default_basis_directory();
}

static json run(const json& input) {
    char* r = macqc_run(input.dump().c_str());
    const json j = json::parse(r);
    macqc_free(r);
    return j;
}

static std::array<double, 3> point(const json& p) { return {p[0].get<double>(), p[1].get<double>(), p[2].get<double>()}; }

// ---------------------------------------------------------------- basis functions

static void test_shell_kernels(const json& ref) {
    for (const json& c : ref["shells"]) {
        const std::string sym = c["element"], basis = c["basis"];
        Molecule mol;
        mol.name = "atom";
        Atom a;
        a.symbol = sym;
        a.Z = sym == "H" ? 1 : sym == "He" ? 2 : sym == "Li" ? 3 : sym == "N" ? 7 : 8;
        a.position = point(c["center"]);
        mol.atoms = {a};
        const BasisSet bs = build_basis(mol, basis, basis_dir());
        std::vector<std::array<double, 3>> pts;
        for (const json& p : c["points"]) pts.push_back(point(p));
        std::vector<std::vector<double>> phi(pts.size(), std::vector<double>(static_cast<std::size_t>(bs.n_functions)));
        for (std::size_t i = 0; i < pts.size(); ++i) basis_function_values(bs.shells, bs.n_functions, pts[i], phi[i].data());
        std::vector<bool> used(bs.shells.size(), false);
        double worst = 0.0;
        int matched = 0, max_l = 0;
        for (const json& rs : c["shells"]) {
            const int l = rs["l"];
            const std::vector<double> ex = rs["exponents"].get<std::vector<double>>();
            const std::vector<double> co = rs["coefficients"].get<std::vector<double>>();
            int found = -1;
            for (std::size_t k = 0; k < bs.shells.size() && found < 0; ++k) {
                const Shell& s = bs.shells[k];
                if (used[k] || s.l != l || s.exponents.size() != ex.size()) continue;
                bool same = true;
                double dot = 0.0, n1 = 0.0, n2 = 0.0;
                for (std::size_t p = 0; p < ex.size() && same; ++p) {
                    same = std::fabs(s.exponents[p] - ex[p]) <= 1e-12 * ex[p];
                    // ΨLab's coefficients include the primitive normalisation (proportional to
                    // a^((2l+3)/4)), PySCF's bas_ctr_coeff does not
                    const double mine = s.coefficients[p] / std::pow(s.exponents[p], (2.0 * l + 3.0) / 4.0);
                    dot += mine * co[p];
                    n1 += mine * mine;
                    n2 += co[p] * co[p];
                }
                if (same && std::fabs(dot) >= (1.0 - 1e-10) * std::sqrt(n1 * n2)) found = static_cast<int>(k);
            }
            if (found < 0) continue;
            used[static_cast<std::size_t>(found)] = true;
            ++matched;
            max_l = std::max(max_l, l);
            const Shell& s = bs.shells[static_cast<std::size_t>(found)];
            double scale = 0.0, diff = 0.0;
            for (std::size_t i = 0; i < pts.size(); ++i)
                for (std::size_t j = 0; j < pts.size(); ++j) {
                    double k = 0.0;
                    for (int m = 0; m < s.size(); ++m)
                        k += phi[i][static_cast<std::size_t>(s.first_function + m)] * phi[j][static_cast<std::size_t>(s.first_function + m)];
                    const double r = rs["kernel"][i][j];
                    scale = std::max(scale, std::fabs(r));
                    diff = std::max(diff, std::fabs(k - r));
                }
            worst = std::max(worst, diff / std::max(scale, 1e-300));
        }
        CHECK(matched == static_cast<int>(c["shells"].size()));
        CHECK(matched == static_cast<int>(bs.shells.size()));
        char label[160];
        std::snprintf(label, sizeof label, "%s %s: %d shells (l <= %d), kernel sum_m phi(r) phi(r') vs PySCF", sym.c_str(),
                      basis.c_str(), matched, max_l);
        check_max(label, worst, 1e-10);
    }
}

// ---------------------------------------------------------------- canonical orbitals

static void test_canonical(const json& ref) {
    for (const json& c : ref["canonical"]) {
        const std::string name = c["molecule"];
        const auto t0 = Clock::now();
        const json out = run({{"molecule", name}, {"basis", c["basis"]}, {"method", "RHF"}, {"basisDirectory", basis_dir()}});
        CHECK(out["error"].is_null());
        CHECK(out["orbitals"].is_object());
        if (!out["orbitals"].is_object()) continue;
        const OrbitalData d = orbitals_from_json(out["orbitals"].dump());
        check_max(name + " RHF/aug-cc-pVDZ energy vs PySCF (E_h)",
                  std::fabs(out["rhf"]["totalEnergy"].get<double>() - c["energy"].get<double>()), 1e-8);
        CHECK(d.sets.size() == 1);
        const OrbitalSet& set = d.sets[0];
        const int nocc = c["occupied"];
        const int n = static_cast<int>(c["orbitalEnergies"].size());
        CHECK(static_cast<int>(set.orbitals.size()) == nocc + 10);
        CHECK(set.orbitals[static_cast<std::size_t>(nocc - 1)].label == "HOMO");
        CHECK(set.orbitals[static_cast<std::size_t>(nocc)].label == "LUMO");
        CHECK(set.orbitals[0].occupation == 2.0 && set.orbitals[static_cast<std::size_t>(nocc)].occupation == 0.0);
        double de = 0.0;
        for (int k = 0; k < n; ++k)
            de = std::max(de, std::fabs(set.orbitals[static_cast<std::size_t>(k)].energy - c["orbitalEnergies"][k].get<double>()));
        check_max(name + ": orbital energies vs PySCF (E_h)", de, 1e-6);
        std::vector<std::array<double, 3>> pts;
        for (const json& p : c["points"]) pts.push_back(point(p));
        // values of the orbitals at the points
        std::vector<std::vector<double>> v(static_cast<std::size_t>(n), std::vector<double>(pts.size()));
        for (int k = 0; k < n; ++k)
            for (std::size_t i = 0; i < pts.size(); ++i)
                v[static_cast<std::size_t>(k)][i] = orbital_value(d, set.coefficients[static_cast<std::size_t>(k)], pts[i]);
        double worst_single = 0.0, worst_group = 0.0, worst_density = 0.0;
        int degenerate_groups = 0;
        for (int k = 0; k < n;) {
            int e = k + 1;
            while (e < n && std::fabs(c["orbitalEnergies"][e].get<double>() - c["orbitalEnergies"][k].get<double>()) < 1e-5) ++e;
            double scale = 0.0;
            for (int q = k; q < e; ++q)
                for (std::size_t i = 0; i < pts.size(); ++i) scale = std::max(scale, std::fabs(c["values"][q][i].get<double>()));
            if (e - k == 1) {
                double plus = 0.0, minus = 0.0;
                for (std::size_t i = 0; i < pts.size(); ++i) {
                    const double r = c["values"][k][i];
                    plus = std::max(plus, std::fabs(v[static_cast<std::size_t>(k)][i] - r));
                    minus = std::max(minus, std::fabs(v[static_cast<std::size_t>(k)][i] + r));
                }
                worst_single = std::max(worst_single, std::min(plus, minus) / scale);
            } else {
                ++degenerate_groups;
                for (std::size_t i = 0; i < pts.size(); ++i) {
                    double a = 0.0, b = 0.0;
                    for (int q = k; q < e; ++q) {
                        a += v[static_cast<std::size_t>(q)][i] * v[static_cast<std::size_t>(q)][i];
                        b += c["values"][q][i].get<double>() * c["values"][q][i].get<double>();
                    }
                    worst_group = std::max(worst_group, std::fabs(a - b) / (scale * scale));
                }
            }
            k = e;
        }
        double rho_scale = 0.0;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            double a = 0.0, b = 0.0;
            for (int q = 0; q < nocc; ++q) {
                a += 2.0 * v[static_cast<std::size_t>(q)][i] * v[static_cast<std::size_t>(q)][i];
                b += 2.0 * c["values"][q][i].get<double>() * c["values"][q][i].get<double>();
            }
            rho_scale = std::max(rho_scale, b);
            worst_density = std::max(worst_density, std::fabs(a - b));
        }
        check_max(name + ": orbitals at 24 points vs PySCF (up to the sign; relative)", worst_single, 2e-5);
        if (degenerate_groups > 0)
            check_max(name + ": degenerate orbitals, sum of squares vs PySCF (relative)", worst_group, 2e-5);
        check_max(name + ": electron density 2 sum_occ psi^2 vs PySCF (relative)", worst_density / rho_scale, 1e-6);
        CHECK(name != "NH3" || degenerate_groups >= 2);
        std::printf("     %s: %.1f s\n", name.c_str(), since(t0));
    }
}

// ---------------------------------------------------------------- natural orbitals

static void test_natural(const json& ref) {
    const json& c = ref["natural"];
    MolecularInput in;
    in.molecule = "LiH";
    in.basis = c["basis"];
    in.method = "FCI";
    in.basis_directory = basis_dir();
    in.print_level = 0;
    MolecularResult m = prepare_molecular(in);
    run_molecular(m);
    CHECK(m.rhf && m.fci);
    if (!m.rhf || !m.fci) return;
    check_max("LiH FCI/aug-cc-pVDZ energy vs PySCF (E_h)", std::fabs(m.fci->total_energy - c["energy"].get<double>()), 1e-8);
    const OrbitalData all = make_orbital_data(m.molecule, m.basis, *m.rhf, &*m.fci, 1000);
    CHECK(all.sets.size() == 2);
    if (all.sets.size() != 2) return;
    const OrbitalSet& nat = all.sets[1];
    CHECK(nat.kind == "natural");
    CHECK(nat.orbitals.size() == static_cast<std::size_t>(m.rhf->n_mo));
    double docc = 0.0, sum = 0.0;
    for (std::size_t k = 0; k < nat.orbitals.size(); ++k) {
        docc = std::max(docc, std::fabs(nat.orbitals[k].occupation - c["occupations"][k].get<double>()));
        sum += nat.orbitals[k].occupation;
    }
    check_max("LiH: natural occupations vs PySCF", docc, 1e-7);
    check_max("LiH: sum of the natural occupations - 4", std::fabs(sum - 4.0), 1e-10);
    // C^T S C = 1 for the natural orbitals
    const Matrix S = overlap_matrix(m.basis);
    double ortho = 0.0;
    for (std::size_t a = 0; a < nat.coefficients.size(); ++a)
        for (std::size_t b = 0; b <= a; ++b) {
            double s = 0.0;
            for (std::size_t mu = 0; mu < S.rows(); ++mu)
                for (std::size_t nu = 0; nu < S.cols(); ++nu) s += nat.coefficients[a][mu] * S(mu, nu) * nat.coefficients[b][nu];
            ortho = std::max(ortho, std::fabs(s - (a == b ? 1.0 : 0.0)));
        }
    check_max("LiH: natural orbitals orthonormal (max |C^T S C - 1|)", ortho, 1e-10);
    double worst = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < c["points"].size(); ++i) {
        const std::array<double, 3> p = point(c["points"][i]);
        double rho = 0.0;
        for (std::size_t k = 0; k < nat.coefficients.size(); ++k) {
            const double x = orbital_value(all, nat.coefficients[k], p);
            rho += nat.orbitals[k].occupation * x * x;
        }
        const double r = c["density"][i];
        scale = std::max(scale, r);
        worst = std::max(worst, std::fabs(rho - r));
    }
    check_max("LiH: Full CI density sum_k n_k chi_k^2 at 18 points vs PySCF (relative)", worst / scale, 1e-7);
    // the output keeps the occupied and ten more natural orbitals
    const OrbitalData shown = make_orbital_data(m.molecule, m.basis, *m.rhf, &*m.fci);
    CHECK(shown.sets.size() == 2 && shown.sets[1].orbitals.size() == 12);
}

// ---------------------------------------------------------------- Szabo-Ostlund HeH+

static void test_szabo_ostlund() {
    const json out = run({{"molecule", "HeH+"}, {"basis", "STO-3G"}});
    CHECK(out["error"].is_null());
    CHECK(out["orbitals"].is_object());
    if (!out["orbitals"].is_object()) return;
    const OrbitalData d = orbitals_from_json(out["orbitals"].dump());
    CHECK(d.basis_name == "STO-3G" && d.n_basis == 2 && d.atoms.size() == 2);
    CHECK(d.sets.size() == 1 && d.sets[0].orbitals.size() == 2);
    CHECK(d.sets[0].orbitals[0].label == "HOMO" && d.sets[0].orbitals[1].label == "LUMO");
    CHECK(std::fabs(d.atoms[1].position[2] - 1.4632) < 1e-12);
    CHECK(std::fabs(d.sets[0].orbitals[0].energy - out["hf"]["orbitalEnergies"][0].get<double>()) < 1e-14);
    // overlaps on a grid (the trapezoidal rule is exact to ~1e-12 for these Gaussians)
    GridSpec g;
    const double h = 0.07;
    for (int k = 0; k < 3; ++k) {
        g.origin[k] = (k == 2 ? 0.7316 : 0.0) - 8.0;
        g.n[k] = static_cast<int>(16.0 / h) + 1;
        g.axes[k] = {0.0, 0.0, 0.0};
        g.axes[k][k] = h;
    }
    const std::vector<double> a = orbital_grid(d, d.sets[0].coefficients[0], g);
    const std::vector<double> b = orbital_grid(d, d.sets[0].coefficients[1], g);
    double aa = 0.0, bb = 0.0, ab = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        aa += a[i] * a[i];
        bb += b[i] * b[i];
        ab += a[i] * b[i];
    }
    const double dv = h * h * h;
    // The Szabo-Ostlund program takes S11 = S22 = 1 for its STO-3G functions, which are normalised
    // only to about 1e-5 (the fit coefficients have six digits).
    check_max("HeH+ STO-3G: integral of psi1^2 - 1 on a grid", std::fabs(aa * dv - 1.0), 2e-4);
    check_max("HeH+ STO-3G: integral of psi2^2 - 1 on a grid", std::fabs(bb * dv - 1.0), 2e-4);
    check_max("HeH+ STO-3G: integral of psi1 psi2 on a grid", std::fabs(ab * dv), 2e-4);
}

// ---------------------------------------------------------------- isosurfaces

static double enclosed_volume(const Mesh& m) {
    double v = 0.0;
    for (std::size_t t = 0; t < m.triangle_count(); ++t) {
        const float* p[3];
        for (int c = 0; c < 3; ++c) p[c] = &m.positions[3 * static_cast<std::size_t>(m.indices[3 * t + c])];
        v += (p[0][0] * (p[1][1] * p[2][2] - p[1][2] * p[2][1]) - p[0][1] * (p[1][0] * p[2][2] - p[1][2] * p[2][0]) +
              p[0][2] * (p[1][0] * p[2][1] - p[1][1] * p[2][0])) / 6.0;
    }
    return v;
}

// Every directed edge once and its reverse once: closed and consistently oriented.
static bool closed_and_oriented(const Mesh& m) {
    std::map<std::pair<int, int>, int> directed;
    for (std::size_t t = 0; t < m.triangle_count(); ++t)
        for (int c = 0; c < 3; ++c) ++directed[{m.indices[3 * t + c], m.indices[3 * t + (c + 1) % 3]}];
    for (const auto& [e, n] : directed) {
        if (n != 1) return false;
        auto r = directed.find({e.second, e.first});
        if (r == directed.end() || r->second != 1) return false;
    }
    return true;
}

static void test_sphere() {
    GridSpec g;
    const double h = 0.05;
    for (int k = 0; k < 3; ++k) {
        g.origin[k] = -1.5 + 0.013 * k;  // not symmetric about the centre
        g.n[k] = 61;
        g.axes[k] = {0.0, 0.0, 0.0};
        g.axes[k][k] = h;
    }
    const std::array<double, 3> c0{0.11, -0.07, 0.05};
    std::vector<double> f(g.size());
    for (int k = 0; k < g.n[2]; ++k)
        for (int j = 0; j < g.n[1]; ++j)
            for (int i = 0; i < g.n[0]; ++i) {
                const auto p = g.point(i, j, k);
                double r2 = 0.0;
                for (int d = 0; d < 3; ++d) r2 += (p[d] - c0[d]) * (p[d] - c0[d]);
                f[g.index(i, j, k)] = std::exp(-r2);
            }
    const Mesh m = isosurface(f, g, 0.5);
    const double R = std::sqrt(std::log(2.0));
    double dr = 0.0, dn = 0.0;
    for (std::size_t v = 0; v < m.vertex_count(); ++v) {
        double r[3], rr = 0.0;
        for (int d = 0; d < 3; ++d) {
            r[d] = m.positions[3 * v + d] - c0[d];
            rr += r[d] * r[d];
        }
        rr = std::sqrt(rr);
        dr = std::max(dr, std::fabs(rr - R));
        double dot = 0.0;
        for (int d = 0; d < 3; ++d) dot += m.normals[3 * v + d] * r[d] / rr;
        dn = std::max(dn, 1.0 - dot);
    }
    std::printf("     sphere: %zu vertices, %zu triangles\n", m.vertex_count(), m.triangle_count());
    CHECK(m.triangle_count() > 1000);
    check_max("isosurface of exp(-r^2) = 1/2: vertex radius - sqrt(ln 2) (bohr; spacing 0.05)", dr, 2e-3);
    check_max("isosurface: 1 - normal . radial direction", dn, 2e-3);
    CHECK(closed_and_oriented(m));
    const double vol = enclosed_volume(m) - 0.0;  // the centre is inside
    check_max("isosurface: enclosed volume / (4/3 pi r^3) - 1 (outward orientation)",
              std::fabs(vol / (4.0 / 3.0 * 3.14159265358979323846 * R * R * R) - 1.0), 1e-2);
    // an empty and a full region give no surface
    CHECK(isosurface(f, g, 2.0).triangle_count() == 0);
    CHECK(isosurface(f, g, -1.0).triangle_count() == 0);
}

static void test_water_homo() {
    const json out = run({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}, {"basisDirectory", basis_dir()}});
    if (!out["orbitals"].is_object()) {
        CHECK(false);
        return;
    }
    const std::string text = out["orbitals"].dump();
    const OrbitalData d = orbitals_from_json(text);
    const int homo = 4;
    const std::vector<double>& c = d.sets[0].coefficients[homo];
    const double iso = 0.05;
    const auto t0 = Clock::now();
    const OrbitalSurfaces s = orbital_surfaces(d, c, iso, 200000);
    std::printf("     H2O HOMO surfaces: %zu + %zu triangles, grid %d x %d x %d, %.2f s\n", s.positive.triangle_count(),
                s.negative.triangle_count(), s.grid.n[0], s.grid.n[1], s.grid.n[2], since(t0));
    CHECK(s.positive.triangle_count() > 500 && s.negative.triangle_count() > 500);
    CHECK(s.grid.size() <= 200000);
    CHECK(s.max_abs > iso);
    double dev = 0.0;
    for (int lobe = 0; lobe < 2; ++lobe) {
        const Mesh& m = lobe == 0 ? s.positive : s.negative;
        for (std::size_t v = 0; v < m.vertex_count(); v += 7) {
            const std::array<double, 3> p{m.positions[3 * v], m.positions[3 * v + 1], m.positions[3 * v + 2]};
            const double x = orbital_value(d, c, p) * (lobe == 0 ? 1.0 : -1.0);
            dev = std::max(dev, std::fabs(x - iso) / iso);
        }
    }
    check_max("H2O HOMO: |psi| at the surface vertices / isovalue - 1 (after the Newton steps)", dev, 1e-3);
    const double vp = enclosed_volume(s.positive), vn = enclosed_volume(s.negative);
    check_max("H2O HOMO (1b1): the two lobes enclose equal volumes (relative difference)", std::fabs(vp - vn) / vp, 0.02);
    CHECK(vp > 0.0 && vn > 0.0);

    // C API: the same surfaces, the plane, errors
    macqc_mesh pos{}, neg{};
    char* err = macqc_orbital_surface(text.c_str(), 0, homo, iso, 200000, &pos, &neg);
    CHECK(err == nullptr);
    CHECK(pos.triangle_count == static_cast<int>(s.positive.triangle_count()));
    CHECK(neg.vertex_count == static_cast<int>(s.negative.vertex_count()));
    CHECK(pos.positions && pos.normals && pos.indices);
    bool same = pos.positions && pos.triangle_count > 0 && pos.positions[0] == s.positive.positions[0];
    CHECK(same);
    macqc_mesh_free(&pos);
    macqc_mesh_free(&neg);
    CHECK(pos.positions == nullptr && pos.vertex_count == 0);
    const double origin[3] = {-3.0, 0.0, -3.0}, u[3] = {0.1, 0.0, 0.0}, w[3] = {0.0, 0.0, 0.1};
    std::vector<float> plane(61 * 61);
    err = macqc_orbital_plane(text.c_str(), 0, 2, origin, u, w, 61, 61, plane.data());
    CHECK(err == nullptr);
    double dp = 0.0;
    for (int j = 0; j < 61; j += 5)
        for (int i = 0; i < 61; i += 5) {
            const std::array<double, 3> p{origin[0] + i * u[0], 0.0, origin[2] + j * w[2]};
            dp = std::max(dp, std::fabs(plane[static_cast<std::size_t>(j * 61 + i)] - orbital_value(d, d.sets[0].coefficients[2], p)));
        }
    check_max("macqc_orbital_plane = orbital_value (float precision)", dp, 1e-6);
    err = macqc_orbital_surface(text.c_str(), 3, 0, iso, 0, &pos, &neg);
    CHECK(err != nullptr && std::string(err).find("no such orbital set") != std::string::npos);
    macqc_free(err);
    err = macqc_orbital_surface(text.c_str(), 0, 999, iso, 0, &pos, &neg);
    CHECK(err != nullptr);
    macqc_free(err);
    err = macqc_orbital_surface(text.c_str(), 0, 0, -1.0, 0, &pos, &neg);
    CHECK(err != nullptr && std::string(err).find("isovalue") != std::string::npos);
    macqc_free(err);
    err = macqc_orbital_surface("{\"broken\":", 0, 0, 0.05, 0, &pos, &neg);
    CHECK(err != nullptr && pos.positions == nullptr && neg.positions == nullptr);
    macqc_free(err);
    err = macqc_orbital_plane(text.c_str(), 0, 0, origin, u, w, 0, 10, plane.data());
    CHECK(err != nullptr);
    macqc_free(err);
}

// ---------------------------------------------------------------- the output JSON

static void test_output() {
    const std::string dir = basis_dir();
    // a job: the orbitals belong to the final (optimised) geometry
    const json opt = run({{"molecule", "H2"}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}, {"job", "opt"}, {"basisDirectory", dir}});
    CHECK(opt["error"].is_null() && opt["orbitals"].is_object());
    if (opt["orbitals"].is_object()) {
        const json& a = opt["orbitals"]["atoms"];
        const json& f = opt["optimization"]["finalGeometry"]["atoms"];
        double d = 0.0;
        for (int i = 0; i < 2; ++i)
            for (int k = 0; k < 3; ++k) d = std::max(d, std::fabs(a[i]["position"][k].get<double>() - f[i]["position"][k].get<double>()));
        check_max("H2 opt: orbitals at the optimised geometry", d, 1e-12);
        CHECK(opt["orbitals"]["sets"].size() == 1);
        CHECK(opt["orbitals"]["sets"][0]["orbitals"].size() == 11);  // 1 occupied + 10 virtual
    }
    // CCSD: canonical orbitals only; Full CI: also natural orbitals
    const json cc = run({{"molecule", "H2"}, {"basis", "aug-cc-pVDZ"}, {"method", "CCSD"}, {"basisDirectory", dir}});
    CHECK(cc["orbitals"].is_object() && cc["orbitals"]["sets"].size() == 1);
    const json f = run({{"molecule", "HeH+"}, {"basis", "aug-cc-pVDZ"}, {"method", "FCI"}, {"basisDirectory", dir}});
    CHECK(f["orbitals"].is_object() && f["orbitals"]["sets"].size() == 2);
    if (f["orbitals"].is_object() && f["orbitals"]["sets"].size() == 2) {
        const json& n = f["orbitals"]["sets"][1];
        CHECK(n["kind"] == "natural" && n["orbitals"][0]["energy"].is_null());
        CHECK(std::fabs(n["orbitals"][0]["occupation"].get<double>() - f["fullCI"]["naturalOccupations"][0].get<double>()) < 1e-14);
    }
    // an estimate has no orbitals; a refused calculation neither
    const json e = run({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}, {"estimateOnly", true}, {"basisDirectory", dir}});
    CHECK(e["orbitals"].is_null());
    const json r = run({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "FCI"}, {"basisDirectory", dir}});
    CHECK(r["orbitals"].is_null());
    // round trip of the JSON
    if (opt["orbitals"].is_object()) {
        const OrbitalData d = orbitals_from_json(opt["orbitals"].dump());
        CHECK(json::parse(orbitals_to_json(d)) == opt["orbitals"]);
    }
}

int main() {
    const auto t0 = Clock::now();
    try {
        const json ref = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/orbitals_reference.json"));
        test_shell_kernels(ref);
        test_canonical(ref);
        test_natural(ref);
        test_szabo_ostlund();
        test_sphere();
        test_water_homo();
        test_output();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t0));
    return g_failures == 0 ? 0 : 1;
}
