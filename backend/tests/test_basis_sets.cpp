// Tests of the basis sets added in ΨLab 14.0 (CTest): the Pople 6-31G family (Cartesian d
// functions, 6D; spherical f, 7F), the 6-311G family and the Karlsruhe def2 sets (spherical).
//  - the library: 33 sets in five families, names and Gaussian's spellings (6-31G(d), def2SVP),
//    the files (BSE header, elements, angular momenta);
//  - Cartesian shells: normalisation and overlaps of the six d functions, the function counts;
//  - against PySCF (backend/tests/data/basis_sets_reference.json from
//    reference/pyscf_basis_sets_reference.py): RHF energies of H2O in all 25 new sets and of
//    molecules with Li, C, N, Cl, K, Ca, Zn (6D 7F), Br; CCSD(T) energies; analytic RHF and CCSD
//    gradients; the RHF Hessian, dipole moment and orbital values of H2O in 6-31G*;
//  - integral-direct SCF with Cartesian shells, the analytic correlated Hessians (CCSD, CCSD(T),
//    Full CI) with Cartesian d functions against differences of analytic gradients;
//  - the orbitals JSON (Cartesian shells survive the round trip; 13.0 data without "pure"), the
//    output file (6D in the basis block, labels such as "3dxy"), the terminal input and refusals.
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/ecp.h"
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/hessian.h"
#include "macqc/macqc_api.h"
#include "macqc/molecule.h"
#include "macqc/orbitals.h"
#include "macqc/rhf.h"
#include "macqc/structure_file.h"
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
static bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static json run(const json& input) {
    char* r = macqc_run(input.dump().c_str());
    const json j = json::parse(r);
    macqc_free(r);
    return j;
}

static Molecule molecule_from_xyz(const std::string& xyz, int charge = 0) {
    const StructureFile f = parse_structure_file(xyz, "xyz");
    return make_custom_molecule(f.atoms, charge, f.title);
}

static json atoms_json(const std::string& xyz) {
    char* r = macqc_parse_structure(xyz.c_str(), "xyz");
    const json p = json::parse(r);
    macqc_free(r);
    return p["atoms"];
}

static ElectronicOptions tight() {
    ElectronicOptions eo;
    eo.scf_energy_tol = 1e-12;
    eo.scf_gradient_tol = 1e-10;
    eo.cc_amplitude_tol = 1e-10;
    eo.cc_energy_tol = 1e-12;
    return eo;
}

// ---------------------------------------------------------------------------------------------

static void test_library() {
    std::printf("-- the basis-set library\n");
    const auto& lib = basis_set_library();
    CHECK(lib.size() == 44 && supported_basis_sets().size() == 44);  // 15.0: with STO-1G/2G/3G, cc-pVXZ, cc-pCVXZ
    std::vector<std::string> families;
    for (const auto& b : lib)
        if (families.empty() || families.back() != b.family) families.push_back(b.family);
    CHECK((families == std::vector<std::string>{"Pople STO-nG", "Pople 6-31G", "Pople 6-311G", "Dunning cc-pVXZ",
                                                "Dunning aug-cc-pVXZ", "Dunning cc-pCVXZ", "Dunning aug-cc-pCVXZ",
                                                "Karlsruhe def2"}));
    int cartesian = 0;
    for (const auto& b : lib) cartesian += b.cartesian_d ? 1 : 0;
    CHECK(cartesian == 7);
    const char* requested[] = {"6-31G", "6-31G*", "6-31+G", "6-31+G*", "6-31+G**", "6-31++G*", "6-31++G**",
                               "6-311G", "6-311G*", "6-311+G", "6-311+G*", "6-311+G**", "6-311++G*", "6-311++G**",
                               "6-311++G(2d,2p)", "6-311++G(3df,3pd)", "def2-SVP", "def2-mTZVP", "def2-mTZVPP",
                               "def2-TZVP", "def2-TZVPP", "def2-TZVPPD", "def2-QZVP", "def2-QZVPP", "def2-QZVPPD"};
    for (const char* n : requested) CHECK(canonical_basis_name(n) == n);
    // spellings: case, blanks, Gaussian's (d) and (d,p), def2 without the hyphen
    CHECK(canonical_basis_name("6-31g(d)") == "6-31G*");
    CHECK(canonical_basis_name("6-31+G(d,p)") == "6-31+G**");
    CHECK(canonical_basis_name("6-31++g(d,p)") == "6-31++G**");
    CHECK(canonical_basis_name("6-311++G(d,p)") == "6-311++G**");
    CHECK(canonical_basis_name("6-311G(d)") == "6-311G*");
    CHECK(canonical_basis_name(" 6-311++G(2D,2P) ") == "6-311++G(2d,2p)");
    CHECK(canonical_basis_name("Def2SVP") == "def2-SVP");
    CHECK(canonical_basis_name("DEF2-QZVPPD") == "def2-QZVPPD");
    CHECK(canonical_basis_name("def2mtzvpp") == "def2-mTZVPP");
    CHECK(canonical_basis_name("6-31G**").empty());     // not in the library
    CHECK(canonical_basis_name("def2-TZVPD").empty());  // not in the library
    CHECK(basis_file_name("6-311++G(3df,3pd)") == "6-311ppg_3df_3pd_.gbs");
    CHECK(basis_file_name("6-31+G**") == "6-31pgss.gbs");
    CHECK(contains(basis_function_types("6-31G*"), "6D"));
    CHECK(contains(basis_function_types("6-311G*"), "5D"));
    CHECK(contains(basis_function_types("def2-TZVP"), "5D"));

    // every file: BSE 0.12 header with the name, elements up to Kr, l <= 4 (g)
    struct Coverage {
        const char* basis;
        const char* elements;  // expected element_ranges text, checked through a refusal below
        int count;
    };
    // 14.0, all elements the BSE defines: I in 6-311G(*); the def2 sets up to Rn with their ECPs
    // (Rb-Rn; the diffuse ones without Ce-Lu), def2-mTZVP(P) also Th-Lr
    const Coverage cov[] = {{"6-31G", "H–Kr", 36},       {"6-31+G*", "H–Ar", 18},   {"6-311G*", "H–Ca, Ga–Kr, I", 27},
                            {"6-311+G", "H–Ca", 20},     {"6-311++G**", "H, Li–Ca", 19},
                            {"6-311++G(3df,3pd)", "H, Li–Ar", 17},                  {"def2-SVP", "H–Rn", 86},
                            {"def2-QZVPPD", "H–La, Hf–Rn", 72},                     {"def2-mTZVPP", "H–Rn, Th–Lr", 100}};
    for (const char* n : requested) {
        const BasisFile f = parse_gaussian94(read_file(default_basis_directory() + "/" + basis_file_name(n)));
        CHECK(contains(f.header, "Basis Set Exchange") && contains(f.header, "Version 0.12"));
        CHECK(contains(f.header, std::string("Basis set: ") + n));
        int maxl = 0;
        for (const auto& [Z, e] : f.elements) {
            CHECK(Z >= 1 && Z <= 103);
            // an ECP exactly for the elements beyond Kr of the def2 sets (6-311G: I all-electron)
            CHECK((e.ecp != nullptr) == (Z > 36 && std::string(n).rfind("def2", 0) == 0));
            for (const auto& sh : e.shells) maxl = std::max(maxl, sh.l);
        }
        CHECK(maxl <= 5);  // h functions: def2-QZVPP of the heavier elements
        for (const auto& c : cov)
            if (std::string(c.basis) == n) CHECK(static_cast<int>(f.elements.size()) == c.count);
    }
}

static void test_cartesian_shells() {
    std::printf("-- Cartesian d shells\n");
    const Molecule h2o = make_molecule("H2O", default_bond_length("H2O"), kDefaultWaterAngle);
    // function counts: 6D in the 6-31G family, 5D elsewhere
    const BasisSet b6 = build_basis(h2o, "6-31G*", "");
    const BasisSet b5 = build_basis(h2o, "6-311G*", "");
    CHECK(b6.n_functions == 19 && b6.cartesian_d && b6.atoms[0].functions == 15 && b6.atoms[0].contraction == "[3s2p1d]");
    CHECK(b5.n_functions == 24 && !b5.cartesian_d && b5.atoms[0].functions == 18);
    const Shell* d = nullptr;
    for (const Shell& s : b6.shells)
        if (s.l == 2) d = &s;
    CHECK(d && !d->pure && d->size() == 6);
    for (const Shell& s : b6.shells) CHECK(s.pure == (s.l != 2));
    for (const Shell& s : b5.shells) CHECK(s.pure);
    CHECK(cartesian_transform(2).size() == 6 && cartesian_transform(3).size() == 10 && shell_transform(*d).size() == 6);

    // normalised components: <xx|xx> = <xy|xy> = 1, <xx|yy> = 1/3, <xx|xy> = 0 (same exponents)
    if (d) {
        const Matrix S = overlap_matrix(b6);
        const int f = d->first_function;  // xx, xy, xz, yy, yz, zz
        double worst = 0.0;
        for (int i = 0; i < 6; ++i) worst = std::max(worst, std::fabs(S(f + i, f + i) - 1.0));
        check_max("Cartesian d of O (6-31G*): unit diagonal of S", worst, 1e-13);
        check_max("<xx|yy> = 1/3", std::fabs(S(f + 0, f + 3) - 1.0 / 3.0), 1e-13);
        check_max("<xx|zz> = 1/3, <yy|zz> = 1/3",
                  std::max(std::fabs(S(f + 0, f + 5) - 1.0 / 3.0), std::fabs(S(f + 3, f + 5) - 1.0 / 3.0)), 1e-13);
        check_max("<xx|xy> = <xy|xz> = 0", std::max(std::fabs(S(f + 0, f + 1)), std::fabs(S(f + 1, f + 2))), 1e-15);
    }
    // unit diagonal of S in all 25 new sets (H2O), and for heavier elements
    double worst = 0.0;
    for (const auto& info : basis_set_library()) {
        if (std::string(info.family).rfind("Dunning", 0) == 0) continue;
        const BasisSet b = build_basis(h2o, info.name, "");
        const Matrix S = overlap_matrix(b);
        for (int i = 0; i < b.n_functions; ++i) worst = std::max(worst, std::fabs(S(i, i) - 1.0));
    }
    check_max("unit diagonal of S, H2O in the 25 new sets", worst, 1e-12);
    const Molecule hbr = molecule_from_xyz("2\n\nBr 0 0 0\nH 0 0 1.4145\n");
    for (const char* n : {"6-31G*", "6-311G*", "def2-QZVPP"}) {
        const BasisSet b = build_basis(hbr, n, "");
        const Matrix S = overlap_matrix(b);
        double w = 0.0;
        for (int i = 0; i < b.n_functions; ++i) w = std::max(w, std::fabs(S(i, i) - 1.0));
        check_max(std::string("unit diagonal of S, HBr in ") + n, w, 1e-12);
    }
}

static void test_energies(const json& ref) {
    std::printf("-- energies against PySCF\n");
    double worst_rhf = 0.0, worst_cc = 0.0;
    for (const auto& r : ref["energies"]) {
        const std::string name = r["molecule"], basis = r["basis"], method = r["method"];
        const json& s = ref["structures"][name];
        const auto t0 = Clock::now();
        const json out = run({{"molecule", "custom"}, {"atoms", atoms_json(s["xyz"])}, {"title", s["title"]},
                              {"charge", s["charge"]}, {"basis", basis}, {"method", method}, {"printLevel", 0}});
        if (!out["error"].is_null()) {
            ++g_checks;
            ++g_failures;
            std::printf("FAIL %s %s %s: %s\n", name.c_str(), basis.c_str(), method.c_str(), out["error"].get<std::string>().c_str());
            continue;
        }
        const double e = method == "RHF" ? out["rhf"]["totalEnergy"].get<double>() : out["coupledCluster"]["totalEnergy"].get<double>();
        const double diff = std::fabs(e - r["energy"].get<double>());
        CHECK(out["basis"]["functions"] == r["basisFunctions"]);
        std::printf("     %-8s %-18s %-8s %3d functions %20.10f  diff %.1e  %5.1f s\n", name.c_str(), basis.c_str(),
                    method.c_str(), out["basis"]["functions"].get<int>(), e, diff, since(t0));
        double& w = method == "RHF" ? worst_rhf : worst_cc;
        w = std::max(w, diff);
    }
    check_max("RHF energies vs PySCF (25 sets for H2O; Li, C, N, Cl, K, Ca, Zn 6D 7F, Br; f and g)", worst_rhf, 1e-8);
    check_max("CCSD(T) energies vs PySCF (6-31G*, def2-SVP, 6-311++G**)", worst_cc, 1e-8);

    // the built-in H2O is the same molecule as the custom one of the reference
    const json bi = run({{"molecule", "H2O"}, {"basis", "6-31G*"}, {"method", "RHF"}, {"printLevel", 0}});
    double want = 0.0;
    for (const auto& r : ref["energies"])
        if (r["molecule"] == "H2O" && r["basis"] == "6-31G*" && r["method"] == "RHF") want = r["energy"];
    check_max("built-in H2O, RHF/6-31G* (Gaussian's 6D) vs PySCF (cart=True)", std::fabs(bi["rhf"]["totalEnergy"].get<double>() - want), 1e-8);

    // integral-direct SCF with Cartesian shells (no integrals kept in memory)
    const Molecule m = molecule_from_xyz(ref["structures"]["H2O-distorted"]["xyz"]);
    const BasisSet b = build_basis(m, "6-31++G**", "");
    ElectronicOptions stored = tight(), direct = tight();
    direct.memory_mb = 0.0;
    const SinglePoint ps = run_single_point(m, b, "RHF", stored);
    const SinglePoint pd = run_single_point(m, b, "RHF", direct);
    CHECK(ps.rhf.eri_stored && !pd.rhf.eri_stored);
    check_max("integral-direct vs stored-integral RHF, 6-31++G** (6D)", std::fabs(ps.energy - pd.energy), 1e-10);
}

static void test_gradients(const json& ref) {
    std::printf("-- analytic gradients against PySCF\n");
    for (const auto& r : ref["gradients"]) {
        const std::string name = r["molecule"], basis = r["basis"], method = r["method"];
        const Molecule m = molecule_from_xyz(ref["structures"][name]["xyz"]);
        const BasisSet b = build_basis(m, basis, "");
        const auto t0 = Clock::now();
        const GradientResult g = run_gradient(m, b, method, tight());
        CHECK(g.analytic);
        double worst = 0.0;
        for (std::size_t x = 0; x < g.gradient.size(); ++x)
            worst = std::max(worst, std::fabs(g.gradient[x] - r["gradient"][x].get<double>()));
        check_max(name + " " + method + "/" + basis + " energy vs PySCF", std::fabs(g.point.energy - r["energy"].get<double>()),
                  method == "RHF" ? 1e-9 : 1e-8);
        // PySCF's CCSD gradient solves its response equations to about 1e-8
        check_max(name + " " + method + "/" + basis + " analytic gradient vs PySCF (" +
                      std::to_string(static_cast<int>(since(t0) * 10) / 10.0).substr(0, 4) + " s)",
                  worst, method == "RHF" ? 2e-8 : 2e-7);
    }
}

static void test_h2o_631gs(const json& ref) {
    std::printf("-- H2O in 6-31G*: RHF Hessian, dipole moment, orbitals\n");
    const json& r = ref["h2o_631gs"];
    const Molecule m = molecule_from_xyz(ref["structures"]["H2O"]["xyz"]);
    const BasisSet b = build_basis(m, "6-31G*", "");
    HessianOptions ho;
    const HessianResult h = run_hessian(m, b, "RHF", tight(), ho);
    CHECK(h.analytic);
    double worst = 0.0;
    for (std::size_t i = 0; i < 9; ++i)
        for (std::size_t j = 0; j < 9; ++j) worst = std::max(worst, std::fabs(h.hessian(i, j) - r["hessian"][i * 9 + j].get<double>()));
    check_max("RHF/6-31G* analytic Hessian vs PySCF (E_h/bohr^2)", worst, 2e-7);
    double wd = 0.0;
    for (int k = 0; k < 3; ++k) wd = std::max(wd, std::fabs(h.dipole[static_cast<std::size_t>(k)] - r["dipole"][k].get<double>()));
    check_max("RHF/6-31G* dipole moment vs PySCF (e bohr)", wd, 1e-8);

    // orbital values at points (Cartesian d functions in basis_function_values); the sign of an
    // orbital is arbitrary: fixed by the first point
    const SinglePoint p = run_single_point(m, b, "RHF", tight());
    const OrbitalData d = make_orbital_data(m, b, p.rhf);
    const int homo = p.rhf.n_occupied - 1;
    double wv = 0.0;
    for (const char* which : {"HOMO", "LUMO"}) {
        int idx = -1;
        for (std::size_t k = 0; k < d.sets[0].orbitals.size(); ++k)
            if (d.sets[0].orbitals[k].label == which) idx = static_cast<int>(k);
        CHECK(idx >= 0);
        if (idx < 0) continue;
        const auto& c = d.sets[0].coefficients[static_cast<std::size_t>(idx)];
        const auto& pts = r["points"];
        const auto& want = r["orbitals"][which];
        double sign = 0.0;
        for (std::size_t k = 0; k < pts.size(); ++k) {
            const std::array<double, 3> x{pts[k][0].get<double>(), pts[k][1].get<double>(), pts[k][2].get<double>()};
            const double v = orbital_value(d, c, x);
            if (k == 0) sign = (v * want[0].get<double>() >= 0.0) ? 1.0 : -1.0;
            wv = std::max(wv, std::fabs(sign * v - want[k].get<double>()));
        }
    }
    check_max("RHF/6-31G* HOMO and LUMO values at 5 points vs PySCF", wv, 1e-8);
    double we = 0.0;
    for (int k = 0; k <= homo + 1; ++k)
        we = std::max(we, std::fabs(p.rhf.orbital_energies[static_cast<std::size_t>(k)] - r["orbitalEnergies"][k].get<double>()));
    check_max("RHF/6-31G* orbital energies vs PySCF", we, 1e-8);

    // the orbitals JSON keeps the Cartesian shells; 13.0 data (no "pure") means spherical
    const OrbitalData back = orbitals_from_json(orbitals_to_json(d));
    bool same = back.shells.size() == d.shells.size();
    for (std::size_t k = 0; same && k < d.shells.size(); ++k) same = back.shells[k].pure == d.shells[k].pure;
    CHECK(same);
    const std::array<double, 3> x{0.4, -0.3, 0.9};
    check_max("orbital value after the JSON round trip", std::fabs(orbital_value(back, d.sets[0].coefficients[homo], x) -
                                                                   orbital_value(d, d.sets[0].coefficients[homo], x)), 1e-15);
    json j = json::parse(orbitals_to_json(d));
    CHECK(j["shells"][0]["pure"] == true);
    for (auto& s : j["shells"]) s.erase("pure");
    const OrbitalData old = orbitals_from_json(j.dump());
    bool all_pure = true;
    for (const Shell& s : old.shells) all_pure = all_pure && s.pure;
    CHECK(all_pure);
}

// Hessian and dipole derivatives by 4-point central differences of analytic gradients and their
// dipole moments (step h bohr), as in test_correlated_hessian.cpp.
struct FDReference {
    Matrix H;
    std::vector<std::array<double, 3>> apt;
};
static FDReference fd_hessian(const Molecule& m, const BasisSet& b, const std::string& method, double h) {
    const std::size_t n3 = 3 * m.atoms.size();
    FDReference r;
    r.H = Matrix(n3, n3);
    r.apt.assign(n3, {0.0, 0.0, 0.0});
    ElectronicOptions eo = tight();
    eo.scf_energy_tol = 1e-13;
    eo.scf_gradient_tol = 1e-11;
    eo.fci_residual_tol = 1e-10;
    eo.cc_amplitude_tol = 1e-11;
    eo.cc_energy_tol = 1e-13;
    for (std::size_t x = 0; x < n3; ++x) {
        std::vector<double> g[4];
        std::array<double, 3> mu[4];
        const double steps[4] = {h, -h, 2 * h, -2 * h};
        for (int k = 0; k < 4; ++k) {
            Molecule d = m;
            d.atoms[x / 3].position[x % 3] += steps[k];
            const GradientResult gr = run_gradient(d, basis_at_geometry(b, d), method, eo);
            g[k] = gr.gradient;
            mu[k] = gr.dipole;
        }
        for (std::size_t y = 0; y < n3; ++y) r.H(x, y) = (8.0 * (g[0][y] - g[1][y]) - (g[2][y] - g[3][y])) / (12.0 * h);
        for (std::size_t c = 0; c < 3; ++c) r.apt[x][c] = (8.0 * (mu[0][c] - mu[1][c]) - (mu[2][c] - mu[3][c])) / (12.0 * h);
    }
    return r;
}

static void test_correlated_hessians() {
    std::printf("-- analytic Hessians with Cartesian d functions vs 4-point differences of analytic gradients\n");
    struct Case {
        const char* xyz;
        const char* basis;
        const char* method;
    };
    const Case cases[] = {
        {"3\n\nO 0 0 0.02\nH 0.80 0.03 0.55\nH -0.73 0 0.61\n", "6-31G*", "RHF"},
        {"3\n\nO 0 0 0.02\nH 0.80 0.03 0.55\nH -0.73 0 0.61\n", "6-31+G*", "CCSD"},
        {"3\n\nO 0 0 0.02\nH 0.80 0.03 0.55\nH -0.73 0 0.61\n", "6-31G*", "CCSD(T)"},
        {"2\n\nLi 0 0 0\nH 0 0 1.62\n", "6-31G*", "FCI"},
    };
    for (const auto& c : cases) {
        const Molecule m = molecule_from_xyz(c.xyz);
        const BasisSet b = build_basis(m, c.basis, "");
        const auto t0 = Clock::now();
        HessianOptions ho;
        ho.response_tol = 1e-11;
        const HessianResult ha = run_hessian(m, b, c.method, tight(), ho);
        CHECK(ha.analytic && ha.has_dipole_derivatives);
        const FDReference fr = fd_hessian(m, b, c.method, 0.005);
        double worst = 0.0, wa = 0.0;
        const std::size_t n = ha.hessian.rows();
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) worst = std::max(worst, std::fabs(ha.hessian(i, j) - fr.H(i, j)));
        for (std::size_t x = 0; x < n && x < ha.dipole_derivatives.size(); ++x)
            for (std::size_t k = 0; k < 3; ++k) wa = std::max(wa, std::fabs(ha.dipole_derivatives[x][k] - fr.apt[x][k]));
        const std::string name = std::string(c.method) + "/" + c.basis + " (" + std::to_string(m.atoms.size()) + " atoms, " +
                                 std::to_string(static_cast<int>(since(t0))) + " s)";
        check_max(name + ": analytic Hessian vs 4-point differences (E_h/bohr^2)", worst, 1e-8);
        check_max(name + ": dipole derivatives vs 4-point differences (e)", wa, 1e-7);
    }
}

static void test_output_and_input() {
    std::printf("-- output file, terminal input, refusals\n");
    const json out = run({{"molecule", "H2O"}, {"basis", "6-31G(d)"}, {"method", "RHF"}, {"printLevel", 1}});
    CHECK(out["error"].is_null());
    CHECK(out["input"]["basis"] == "6-31G*");
    const std::string text = out["report"].get<std::string>() + (out["output"].is_string() ? out["output"].get<std::string>() : "");
    CHECK(contains(text, "6D"));
    // The 13.0 output file (appendix) labels Cartesian functions by their component
    const json job = run({{"molecule", "H2O"}, {"basis", "6-31G*"}, {"method", "RHF"}, {"job", "energy"}, {"printLevel", 1}});
    const std::string o = job["output"].is_string() ? job["output"].get<std::string>() : "";
    CHECK(contains(o, "Cartesian d functions (6D)"));
    CHECK(contains(o, "3dxy") && contains(o, "3dzz"));
    CHECK(contains(o, "Cartesian d: xx, xy, xz, yy, yz, zz"));
    const json sph = run({{"molecule", "H2O"}, {"basis", "def2-TZVP"}, {"method", "RHF"}, {"job", "energy"}, {"printLevel", 1}});
    const std::string os = sph["output"].is_string() ? sph["output"].get<std::string>() : "";
    CHECK(contains(os, "spherical harmonics (5D, 7F)") && contains(os, "3d+2") && !contains(os, "3dxy"));

    // refusals with the coverage of the set
    auto error_of = [](const json& in) {
        const json r = run(in);
        return r["error"].is_string() ? r["error"].get<std::string>() : std::string();
    };
    const std::string he = error_of({{"molecule", "HeH+"}, {"basis", "6-311++G**"}, {"method", "RHF"}});
    CHECK(contains(he, "6-311++G** is not available for He") && contains(he, "H, Li–Ca"));
    const std::string k = error_of({{"molecule", "custom"}, {"atoms", atoms_json("2\n\nK 0 0 0\nH 0 0 2.24\n")}, {"basis", "6-31+G*"}});
    CHECK(contains(k, "not available for K") && contains(k, "H–Ar"));
    // 14.0: the def2 sets end at Rn (def2-mTZVP(P): Th-Lr); the diffuse ones skip Ce-Lu; Rb needs an ECP
    const std::string u = error_of({{"molecule", "custom"}, {"atoms", atoms_json("2\n\nU 0 0 0\nO 0 0 1.8\n")}, {"basis", "def2-SVP"}});
    CHECK(contains(u, "not available for U") && contains(u, "H–Rn") && contains(u, "Th–Lr"));
    const std::string ce = error_of({{"molecule", "custom"}, {"atoms", atoms_json("2\n\nCe 0 0 0\nO 0 0 1.8\n")}, {"basis", "def2-TZVPPD"}});
    CHECK(contains(ce, "not available for Ce") && contains(ce, "lanthanides"));
    const std::string rb = error_of({{"molecule", "custom"}, {"atoms", atoms_json("2\n\nRb 0 0 0\nH 0 0 2.37\n")}, {"basis", "aug-cc-pVDZ"}});
    CHECK(contains(rb, "not available for Rb") && contains(rb, "effective core potentials") && contains(rb, "def2"));
    const std::string unknown = error_of({{"molecule", "H2O"}, {"basis", "6-31G**"}});
    CHECK(contains(unknown, "6-31G**") && contains(unknown, "6-311++G(3df,3pd)"));
    // HeH+ in a Pople set is the general molecule (all methods), not the Szabo-Ostlund calculation
    const json heh = run({{"molecule", "HeH+"}, {"basis", "6-311G"}, {"method", "FCI"}, {"printLevel", 0}});
    CHECK(heh["error"].is_null() && heh["fullCI"].is_object() && heh["basis"]["functions"] == 6);
}

// ---- 14.0: elements beyond Kr: the def2 sets with their effective core potentials ---------------

static void test_ecp(const json& ref) {
    std::printf("-- effective core potentials (def2 sets beyond Kr) and all-electron iodine\n");
    // modified spherical Bessel functions e^-x i_n(x) against closed forms and the recurrence, in
    // every regime of scaled_bessel_i (series x <= 0.5, Miller's recurrence, upward recurrence)
    {
        double worst = 0.0, worst_rec = 0.0;
        for (double x : {0.6, 1.0, 2.5, 7.0, 18.0, 29.9, 30.1, 45.0, 120.0, 4000.0}) {
            double v[16];
            scaled_bessel_i(14, x, v);
            const double e2 = std::exp(-2.0 * x);
            const double sh = 0.5 * (1.0 - e2), ch = 0.5 * (1.0 + e2);  // e^-x sinh x, e^-x cosh x
            const double i0 = sh / x, i1 = ch / x - sh / (x * x), i2 = (3.0 / (x * x) + 1.0) * sh / x - 3.0 * ch / (x * x);
            worst = std::max({worst, std::fabs(v[0] - i0) / i0, std::fabs(v[1] - i1) / i1, std::fabs(v[2] - i2) / i2});
            for (int n = 1; n < 14; ++n)
                worst_rec = std::max(worst_rec, std::fabs(v[n - 1] - v[n + 1] - (2.0 * n + 1.0) / x * v[n]) / v[n - 1]);
        }
        check_max("e^-x i_n(x), n = 0, 1, 2, against closed forms (relative)", worst, 1e-13);
        check_max("e^-x i_n(x), n <= 14: recurrence i_(n-1) - i_(n+1) = (2n+1)/x i_n (relative)", worst_rec, 1e-12);
        double a[16], b[16], cont = 0.0;
        for (double x : {0.5, 30.0}) {  // the regime boundaries for nmax = 10 (2 nmax + 10 = 30)
            scaled_bessel_i(10, x, a);
            scaled_bessel_i(10, x * (1.0 + 1e-12), b);
            for (int n = 0; n <= 10; ++n) cont = std::max(cont, std::fabs(a[n] - b[n]) / a[n]);
        }
        check_max("e^-x i_n(x): continuity across the regimes (relative)", cont, 1e-10);
    }
    const json& S = ref["structures"];
    auto mol_of = [&](const std::string& name) { return molecule_from_xyz(S[name]["xyz"].get<std::string>(), S[name]["charge"].get<int>()); };
    // the radial quadrature: 64 points (default) against 160 for every shell pair of HI in def2-QZVPP
    {
        const Molecule m0 = mol_of("HI");
        const BasisSet b = build_basis(m0, "def2-QZVPP", "");
        const Molecule m = attach_ecps(m0, b);
        CHECK(m.atoms[1].ecp && m.atoms[1].ecp->core_electrons == 28 && !m.atoms[0].ecp);
        double worst = 0.0, largest = 0.0;
        std::vector<double> c1, c2;
        for (const Shell& A : b.shells)
            for (const Shell& B : b.shells) {
                c1.assign(static_cast<std::size_t>(A.cartesian_size()) * B.cartesian_size(), 0.0);
                c2 = c1;
                add_ecp_cartesian(A, B, m.atoms[1], c1);
                add_ecp_cartesian(A, B, m.atoms[1], c2, 160);
                for (std::size_t k = 0; k < c1.size(); ++k) {
                    worst = std::max(worst, std::fabs(c1[k] - c2[k]));
                    largest = std::max(largest, std::fabs(c1[k]));
                }
            }
        std::printf("     largest ECP integral of HI/def2-QZVPP: %.3f E_h\n", largest);
        check_max("ECP integrals HI/def2-QZVPP: 64 against 160 radial points (relative to the largest)", worst / largest,
                  1e-13);
    }
    // energies against PySCF
    double worst_rhf = 0.0, worst_cc = 0.0;
    for (const auto& r : ref["energies"]) {
        const std::string name = r["molecule"], basis = r["basis"], method = r["method"];
        const auto t0 = Clock::now();
        const json out = run({{"molecule", "custom"}, {"atoms", atoms_json(S[name]["xyz"])}, {"title", S[name]["title"]},
                              {"charge", S[name]["charge"]}, {"basis", basis}, {"method", method}, {"printLevel", 0}});
        if (!out["error"].is_null()) {
            ++g_checks;
            ++g_failures;
            std::printf("FAIL %s %s %s: %s\n", name.c_str(), basis.c_str(), method.c_str(), out["error"].get<std::string>().c_str());
            continue;
        }
        const double e = method == "RHF" ? out["rhf"]["totalEnergy"].get<double>() : out["coupledCluster"]["totalEnergy"].get<double>();
        const double diff = std::fabs(e - r["energy"].get<double>());
        CHECK(out["basis"]["functions"] == r["basisFunctions"]);
        CHECK(out["molecule"]["electrons"] == r["electrons"]);
        CHECK(out["molecule"]["coreElectrons"] == r["coreElectrons"]);
        std::printf("     %-13s %-12s %-8s %3d functions %3d+%2d electrons %20.10f  diff %.1e  %5.1f s\n", name.c_str(),
                    basis.c_str(), method.c_str(), out["basis"]["functions"].get<int>(), out["molecule"]["electrons"].get<int>(),
                    out["molecule"]["coreElectrons"].get<int>(), e, diff, since(t0));
        double& w = method == "RHF" ? worst_rhf : worst_cc;
        w = std::max(w, diff);
    }
    check_max("RHF energies with ECP28/46/60, lanthanide ECP28 (Lu), actinide ECP60MWB (Th); I all-electron, vs PySCF",
              worst_rhf, 1e-8);
    check_max("CCSD(T)/def2-SVP of HI (ECP28) vs PySCF", worst_cc, 1e-8);
    // analytic gradients (the ECP derivatives through translational invariance) and Hessians
    for (const auto& r : ref["gradients"]) {
        const std::string name = r["molecule"], basis = r["basis"], method = r["method"];
        const Molecule m0 = mol_of(name);
        const BasisSet b = build_basis(m0, basis, "");
        const Molecule m = attach_ecps(m0, b);
        const auto t0 = Clock::now();
        const GradientResult g = run_gradient(m, b, method, tight());
        CHECK(g.analytic);
        double worst = 0.0;
        for (std::size_t x = 0; x < g.gradient.size(); ++x)
            worst = std::max(worst, std::fabs(g.gradient[x] - r["gradient"][x].get<double>()));
        // PySCF's CCSD gradient solves its response equations to about 1e-8
        check_max(name + " " + method + "/" + basis + " analytic gradient vs PySCF (" +
                      std::to_string(static_cast<int>(since(t0))) + " s)",
                  worst, method == "RHF" ? 2e-8 : 2e-7);
    }
    for (const auto& r : ref["hessians"]) {
        const std::string name = r["molecule"], basis = r["basis"];
        const Molecule m0 = mol_of(name);
        const BasisSet b = build_basis(m0, basis, "");
        const Molecule m = attach_ecps(m0, b);
        const HessianResult h = run_hessian(m, b, "RHF", tight(), HessianOptions());
        CHECK(h.analytic);
        const std::size_t n = h.hessian.rows();
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) worst = std::max(worst, std::fabs(h.hessian(i, j) - r["hessian"][i * n + j].get<double>()));
        check_max(name + " RHF/" + basis + " analytic Hessian vs PySCF (E_h/bohr^2)", worst, 2e-7);
    }
    // the analytic CCSD(T) Hessian with an ECP against 4-point differences of analytic gradients (Full CI
    // is impractical here: the lightest ECP atoms keep 9 or more valence electrons)
    {
        const Molecule m0 = mol_of("HI-stretched");
        const BasisSet b = build_basis(m0, "def2-SVP", "");
        const Molecule m = attach_ecps(m0, b);
        const auto t0 = Clock::now();
        HessianOptions ho;
        ho.response_tol = 1e-11;
        const HessianResult ha = run_hessian(m, b, "CCSD(T)", tight(), ho);
        CHECK(ha.analytic && ha.has_dipole_derivatives);
        const FDReference fr = fd_hessian(m, b, "CCSD(T)", 0.005);
        double worst = 0.0, wa = 0.0;
        for (std::size_t i = 0; i < ha.hessian.rows(); ++i)
            for (std::size_t j = 0; j < ha.hessian.cols(); ++j) worst = std::max(worst, std::fabs(ha.hessian(i, j) - fr.H(i, j)));
        for (std::size_t x = 0; x < ha.dipole_derivatives.size(); ++x)
            for (std::size_t k = 0; k < 3; ++k) wa = std::max(wa, std::fabs(ha.dipole_derivatives[x][k] - fr.apt[x][k]));
        check_max("CCSD(T)/def2-SVP HI (ECP28), analytic Hessian vs 4-point differences, " +
                      std::to_string(static_cast<int>(since(t0))) + " s (E_h/bohr^2)",
                  worst, 1e-8);
        check_max("CCSD(T)/def2-SVP HI (ECP28), dipole derivatives vs 4-point differences (e)", wa, 1e-7);
    }
    // rotational invariance (the projectors about the ECP centre): LuF3 and AuH turned arbitrarily
    {
        const double c1 = std::cos(0.7), s1 = std::sin(0.7), c2 = std::cos(-1.9), s2 = std::sin(-1.9);
        auto turn = [&](Molecule m) {
            for (Atom& a : m.atoms) {
                const auto p = a.position;
                const double x = c1 * p[0] - s1 * p[1], y = s1 * p[0] + c1 * p[1], z = p[2];
                a.position = {x, c2 * y - s2 * z, s2 * y + c2 * z};
                a.position[0] += 0.3;
            }
            return m;
        };
        for (const auto& [name, basis] : {std::pair<std::string, std::string>{"AuH", "def2-TZVPP"}, {"PbO", "def2-SVP"}}) {
            const Molecule m0 = mol_of(name), m1 = turn(m0);
            const BasisSet b0 = build_basis(m0, basis, ""), b1 = build_basis(m1, basis, "");
            const double e0 = run_single_point(attach_ecps(m0, b0), b0, "RHF", tight()).energy;
            const double e1 = run_single_point(attach_ecps(m1, b1), b1, "RHF", tight()).energy;
            check_max(name + "/" + basis + ": energy of the turned and shifted molecule (E_h)", std::fabs(e0 - e1), 1e-9);
        }
    }
    // the output: electrons, core electrons and the ECP in the basis block and the appendix
    const json out = run({{"molecule", "custom"}, {"atoms", atoms_json(S["HI"]["xyz"])}, {"basis", "def2-SVP"},
                          {"method", "RHF"}, {"job", "energy"}, {"printLevel", 1}});
    const std::string o = out["output"].is_string() ? out["output"].get<std::string>() : "";
    CHECK(contains(o, "26 electrons (and 28 core electrons in ECPs)"));
    CHECK(contains(o, "26 functions, ECP28"));
    CHECK(contains(o, "ECP28: 28 core electrons"));
    CHECK(out["basis"]["atoms"][1]["ecp"] == "ECP28" && out["basis"]["atoms"][1]["ecpElectrons"] == 28);
    CHECK(out["basis"]["atoms"][0]["ecp"].is_null());
}

int main() {
    const std::string data = MACQC_TEST_DATA_DIR;
    try {
        test_library();
        test_cartesian_shells();
        const json ref = json::parse(read_file(data + "/basis_sets_reference.json"));
        test_energies(ref);
        test_gradients(ref);
        test_h2o_631gs(ref);
        test_correlated_hessians();
        test_output_and_input();
        test_ecp(json::parse(read_file(data + "/ecp_reference.json")));
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
