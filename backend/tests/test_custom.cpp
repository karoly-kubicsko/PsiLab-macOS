// Tests of the custom molecules of ΨLab 9.0 (CTest): any closed-shell molecule of the elements
// H-Ar and Ga-Kr from Cartesian coordinates.
//  - structure files: XYZ (title line, element symbols or atomic numbers, labels such as C1,
//    several frames) and Gaussian input (Link 0, route, title, charge and multiplicity, freeze
//    codes, Units=AU), the errors (Z-matrix, wrong atom count, ghost and dummy atoms, unknown
//    elements) and the C API macqc_parse_structure;
//  - the molecular formula (Hill order) and the input checks (odd electron count, multiplicity,
//    elements without basis functions);
//  - point groups and rotational symmetry numbers of molecules of many groups (C1 ... Oh, C∞v, D∞h),
//    and the species of the normal modes (D2h, Td, D3d, D∞h);
//  - against PySCF (backend/tests/data/custom_reference.json from
//    reference/pyscf_custom_reference.py): RHF energies of molecules with the elements added in
//    9.0 (aug-cc-pVDZ and aug-cc-pCVDZ; Cu and Zn), ions, single atoms, CCSD(T) energies, the RHF gradient of HOCl;
//  - the order of the atoms does not matter (NH3 and H2O from files against the built-in molecules);
//  - the Cartesian numerical gradients (energies) and Hessians (gradients) of custom molecules
//    against the analytic ones, IR intensities included;
//  - a job through the C API (optimisation + frequencies, relaxed scan) and the refusals.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/frequencies.h"
#include "macqc/gradient.h"
#include "macqc/hessian.h"
#include "macqc/macqc_api.h"
#include "macqc/molecule.h"
#include "macqc/structure_file.h"
#include "macqc/symmetry.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_max(const std::string& name, double worst, double tol) {
    ++g_checks;
    const bool ok = worst <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-80s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
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

static json parse_c_api(const std::string& text, const char* hint) {
    char* r = macqc_parse_structure(text.c_str(), hint);
    const json j = json::parse(r);
    macqc_free(r);
    return j;
}

// The parse error of a structure text ("" if it parses).
static std::string parse_error(const std::string& text, const std::string& hint = "") {
    try {
        parse_structure_file(text, hint);
        return "";
    } catch (const std::exception& e) {
        return e.what();
    }
}

static Molecule molecule_from_xyz(const std::string& xyz, int charge = 0) {
    const StructureFile f = parse_structure_file(xyz, "xyz");
    return make_custom_molecule(f.atoms, charge, f.title);
}

// Atoms of a structure as the JSON input of the C API ("atoms" in angstrom).
static json input_for(const std::string& xyz, int charge, const std::string& basis, const std::string& method) {
    const json p = parse_c_api(xyz, "xyz");
    return {{"molecule", "custom"}, {"atoms", p["atoms"]}, {"title", p["title"]}, {"charge", charge},
            {"basis", basis},       {"method", method},    {"printLevel", 0}};
}

static void test_structure_files() {
    std::printf("-- structure files\n");
    {
        const StructureFile f = parse_structure_file("3\nCO$_2$ (carbon dioxide)\nC 0 0 0\nO 0 0 1.16\r\nO 0 0 -1.16\n", "co2.xyz");
        CHECK(f.format == "xyz");
        CHECK(f.title == "CO$_2$ (carbon dioxide)");
        CHECK(f.atoms.size() == 3 && f.atoms[0].Z == 6 && f.atoms[1].symbol == "O");
        check_max("XYZ: angstrom -> bohr", std::fabs(f.atoms[1].position[2] - 1.16 / kBohrInAngstrom), 1e-14);
        CHECK(!f.has_charge && f.charge == 0 && f.multiplicity == 1 && f.warnings.empty());
    }
    {
        // atomic numbers, labels with numbers, lower case, tabs, extra columns, D exponents
        const StructureFile f = parse_structure_file("4\n\n6\t0.0 0.0 0.0  extra\nCl1 1.77D0 0 0\nh 0 1.09 0\nBR -1.94 0 0\n");
        CHECK(f.atoms.size() == 4 && f.atoms[0].Z == 6 && f.atoms[1].Z == 17 && f.atoms[2].Z == 1 && f.atoms[3].Z == 35);
        CHECK(f.title.empty());
        check_max("XYZ: Fortran D exponent", std::fabs(f.atoms[1].position[0] * kBohrInAngstrom - 1.77), 1e-12);
    }
    {
        const StructureFile f = parse_structure_file("2\nframe 1\nH 0 0 0\nH 0 0 0.74\n2\nframe 2\nH 0 0 0\nH 0 0 0.80\n", "xyz");
        CHECK(f.atoms.size() == 2 && f.title == "frame 1" && f.warnings.size() == 1);
    }
    CHECK(contains(parse_error("3\nx\nO 0 0 0\nH 0 0 1\n", "xyz"), "ends after 2 of 3 atoms"));
    CHECK(contains(parse_error("2\nx\nO 0 0 0\nH 0 0 1\nH 1 0 0\n", "xyz"), "more lines than the 2 atoms"));
    CHECK(contains(parse_error("1\nx\nQq 0 0 0\n", "xyz"), "unknown element 'Qq'"));
    CHECK(contains(parse_error("1\nx\nXx 0 0 0\n", "xyz"), "dummy atoms"));
    CHECK(contains(parse_error("1\nx\nK 0 0 0\n", "xyz"), ""));  // K is known (no basis set: refused later)
    CHECK(parse_error("1\nx\nRb 0 0 0\n", "xyz").empty());   // any element is read (Z = 1 ... 118)
    CHECK(parse_error("1\nx\nOg 0 0 0\n", "xyz").empty());
    CHECK(parse_error("1\nx\n118 0 0 0\n", "xyz").empty());
    CHECK(contains(parse_error("1\nx\n119 0 0 0\n", "xyz"), "Z = 1 to 118"));
    CHECK(contains(parse_error("1\nx\nO 0 zero 0\n", "xyz"), "must be numbers"));
    CHECK(contains(parse_error("O 0 0 0\n", "xyz"), "starts with the number of atoms"));
    CHECK(contains(parse_error("hello\n"), "not a structure file"));
    CHECK(parse_error("2\nx\nO 0 0 0\nH 0 0 1\n", "molecule.txt").empty());  // other extensions: recognised from the text

    const std::string gjf =
        "%chk=water.chk\n%mem=1GB\n#p HF/aug-cc-pVDZ Opt\n\nwater, H$_2$O\nsecond title line\n\n0 1\n"
        "O     0.000000    0.000000    0.117300\nH(Fragment=1)  0  0.000000    0.757200   -0.469200\n"
        "H–H_  -1  0.000000   -0.757200   -0.469200\n\nB 1 2 S 10 0.05\n\n";
    {
        const StructureFile f = parse_structure_file(gjf, "water.gjf");
        CHECK(f.format == "gaussian");
        CHECK(f.title == "water, H$_2$O second title line");
        CHECK(f.has_charge && f.charge == 0 && f.multiplicity == 1);
        CHECK(f.atoms.size() == 3 && f.atoms[0].Z == 8 && f.atoms[1].Z == 1 && f.atoms[2].Z == 1);
        check_max("Gaussian input: freeze code skipped", std::fabs(f.atoms[2].position[1] * kBohrInAngstrom + 0.7572), 1e-12);
    }
    {
        const StructureFile f = parse_structure_file("# HF/STO-3G units=au\n\nOH-\n\n-1 1\nO 0 0 0\nH 0 0 1.8\n");
        CHECK(f.format == "gaussian" && f.charge == -1 && f.atoms.size() == 2);
        check_max("Gaussian input: Units=AU (bohr)", std::fabs(f.atoms[1].position[2] - 1.8), 1e-14);
    }
    {
        const StructureFile f = parse_structure_file("#p B3LYP\n\noxygen\n\n0 3\nO 0 0 0\nO 0 0 1.21\n", "o2.com");
        CHECK(f.multiplicity == 3);
    }
    {
        // 13.0: Gaussian's comments ("!" anywhere on a line; comment-only lines anywhere, also
        // before the route and inside a section), as in the Gaussian input that ΨLab writes,
        // with the general basis input after the coordinates
        const std::string commented =
            "! CO2: geometry optimised at CCSD(T)/aug-cc-pCVDZ with PsiLab, E = -188.1 E_h\n"
            "! all electrons correlated, as in PsiLab (Full: no frozen core)\n"
            "#p CCSD(T,Full)/Gen 5D 7F\n\nCO2   ! the title\n! a comment line inside the title section\n\n0 1\n"
            "C   0.0 0.0 0.0\n! a comment line between atoms\nO   0.0 0.0 1.16   ! oxygen\nO 0.0 0.0 -1.16\n\n"
            "C 0\nS   2 1.00\n  1.0E+03  5.0E-01\n  1.0E+02  5.0E-01\n****\n\n";
        const StructureFile f = parse_structure_file(commented, "co2.gjf");
        CHECK(f.format == "gaussian" && f.title == "CO2" && f.charge == 0 && f.multiplicity == 1);
        CHECK(f.atoms.size() == 3 && f.atoms[0].Z == 6 && f.atoms[1].Z == 8 && f.atoms[2].Z == 8);
        CHECK(parse_structure_file(commented).format == "gaussian");  // recognised from the text
        // messages give the line of the file (line 8: the Z-matrix line after two comment lines)
        CHECK(contains(parse_error("! one\n! two\n# HF\n\nz\n\n0 1\nO\n", "x.gjf"), "line 8:"));
    }
    CHECK(contains(parse_error("# HF\n\nzmat\n\n0 1\nO\nH 1 0.96\nH 1 0.96 2 104.5\n", "x.gjf"), "Z-matrix"));
    CHECK(contains(parse_error("# HF\n\nghost\n\n0 1\nO 0 0 0\nBq 0 0 1\n", "x.gjf"), "ghost atoms"));
    CHECK(contains(parse_error("# HF\n\ndummy\n\n0 1\nX 0 0 0\nO 0 0 1\n", "x.gjf"), "dummy atoms"));
    CHECK(contains(parse_error("# HF geom=check\n\nx\n\n0 1\n", "x.gjf"), "checkpoint"));
    CHECK(contains(parse_error("%chk=a\n\ntitle\n", "x.gjf"), "route section"));

    // the C API
    const json p = parse_c_api(gjf, "water.gjf");
    CHECK(p["error"].is_null() && p["format"] == "gaussian" && p["atoms"].size() == 3 && p["hasCharge"] == true);
    CHECK(p["atoms"][0]["symbol"] == "O" && std::fabs(p["atoms"][1]["y"].get<double>() - 0.7572) < 1e-12);
    const json bad = parse_c_api("2\nx\nO 0 0 0\n", nullptr);
    CHECK(bad["error"].is_string() && bad["atoms"].empty());
}

static void test_formula_and_checks() {
    std::printf("-- formula and input checks\n");
    CHECK(molecule_from_xyz("3\n\nO 0 0 0\nC 0 0 1.16\nO 0 0 2.32\n").formula() == "CO2");
    CHECK(molecule_from_xyz("5\n\nH 0 0 1\nH 0 1 0\nC 0 0 0\nCl 1 0 0\nBr -1 0 0\n").formula() == "CH2BrCl");
    CHECK(molecule_from_xyz("3\n\nH 0 0 0.96\nO 0 0 0\nH 0.96 0 0\n").formula() == "H2O");
    CHECK(molecule_from_xyz("5\n\nN 0 0 0\nH 0.59 0.59 0.59\nH -0.59 -0.59 0.59\nH -0.59 0.59 -0.59\nH 0.59 -0.59 -0.59\n", 1)
              .formula() == "H4N+");
    CHECK(molecule_from_xyz("5\n\nS 0 0 0\nO 1.5 0 0\nO -1.5 0 0\nO 0 1.5 0\nO 0 -1.5 0\n", -2).formula() == "O4S 2-");
    bool threw = false;
    try {
        molecule_from_xyz("2\n\nH 0 0 0\nH 0 0 0.1\n");
    } catch (const std::invalid_argument& e) {
        threw = contains(e.what(), "closer than 0.3 bohr");
    }
    CHECK(threw);
    auto estimate = [](const json& extra) {
        json in = {{"molecule", "custom"}, {"estimateOnly", true}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}};
        for (auto it = extra.begin(); it != extra.end(); ++it) in[it.key()] = it.value();
        const json j = run(in);
        if (!j["error"].is_null()) return j["error"].get<std::string>();
        if (j.contains("job") && j["job"].is_object() && !j["job"]["estimate"]["refusal"].is_null())
            return j["job"]["estimate"]["refusal"].get<std::string>();
        return std::string();
    };
    const json oh = json::array({{{"symbol", "O"}, {"x", 0}, {"y", 0}, {"z", 0}}, {{"symbol", "H"}, {"x", 0}, {"y", 0}, {"z", 0.97}}});
    CHECK(contains(estimate({{"atoms", oh}}), "an odd number"));
    CHECK(estimate({{"atoms", oh}, {"charge", -1}}).empty());
    // 15.0: open shells with UHF or ROHF; RHF itself needs a closed-shell singlet
    CHECK(contains(estimate({{"atoms", oh}, {"charge", -1}, {"multiplicity", 3}}), "needs a closed-shell singlet"));
    CHECK(estimate({{"atoms", oh}, {"charge", -1}, {"multiplicity", 3}, {"method", "HF"}}).empty());
    CHECK(estimate({{"atoms", oh}, {"multiplicity", 2}, {"method", "UMP2"}}).empty());
    CHECK(contains(estimate({}), "needs \"atoms\""));
    const json kr = json::array({{{"symbol", "Kr"}, {"x", 0}, {"y", 0}, {"z", 0}}});
    CHECK(estimate({{"atoms", kr}}).empty());
    CHECK(contains(estimate({{"atoms", kr}, {"basis", "aug-cc-pCVDZ"}}), "not available for Kr"));
    CHECK(contains(estimate({{"atoms", kr}, {"job", "opt"}}), "single atom"));
    const json k = json::array({{{"symbol", "K"}, {"x", 0}, {"y", 0}, {"z", 0}}, {{"symbol", "H"}, {"x", 0}, {"y", 0}, {"z", 2.2}}});
    CHECK(contains(estimate({{"atoms", k}}), "no aug-cc-pVXZ sets for K and Ca"));
    const json rb = json::array({{{"symbol", "Rb"}, {"x", 0}, {"y", 0}, {"z", 0}}, {{"symbol", "H"}, {"x", 0}, {"y", 0}, {"z", 2.4}}});
    CHECK(contains(estimate({{"atoms", rb}}), "effective core potentials"));
    const json zn = json::array({{{"symbol", "Zn"}, {"x", 0}, {"y", 0}, {"z", 0}}});
    CHECK(estimate({{"atoms", zn}, {"basis", "aug-cc-pVQZ"}}).empty());
    CHECK(contains(estimate({{"atoms", zn}, {"basis", "aug-cc-pV5Z"}}), "i functions"));
    const json h2o = json::array({{{"symbol", "H"}, {"x", 0.757}, {"y", 0}, {"z", 0.587}},
                                  {{"symbol", "O"}, {"x", 0}, {"y", 0}, {"z", 0}},
                                  {{"symbol", "H"}, {"x", -0.757}, {"y", 0}, {"z", 0.587}}});
    CHECK(contains(estimate({{"atoms", h2o}, {"method", "CCSDT"}, {"job", "freq"}}), "Hessian from energies is not available"));
    CHECK(estimate({{"atoms", h2o}, {"method", "CCSDT"}, {"job", "opt"}}).empty());  // numerical gradients: 6 energies
    const json e = run({{"molecule", "custom"}, {"atoms", h2o}, {"estimateOnly", true}, {"job", "opt+freq"}, {"method", "CCSD"}});
    CHECK(e["error"].is_null() && e["molecule"]["formula"] == "H2O" && e["molecule"]["custom"] == true);
    CHECK(e["input"]["atoms"].size() == 3 && e["input"]["multiplicity"] == 1);
    // 10.0: the analytic CCSD Hessian (3 internal directions); with "hessian": "numerical" the
    // central differences of 6 analytic gradients along the 3 Cartesian directions
    auto has_note = [](const json& j, const std::string& part) {
        for (const auto& n : j["job"]["estimate"]["notes"])
            if (contains(n.get<std::string>(), part)) return true;
        return false;
    };
    CHECK(has_note(e, "amplitude responses along the 3 internal directions"));
    const json en = run({{"molecule", "custom"}, {"atoms", h2o}, {"estimateOnly", true}, {"job", "opt+freq"}, {"method", "CCSD"},
                         {"hessian", "numerical"}});
    CHECK(has_note(en, "6 analytic gradients along the 3 Cartesian directions"));
}

struct GroupCase {
    const char* name;
    const char* xyz;
    const char* group;
    int sigma;
};

static void test_point_groups() {
    std::printf("-- point groups\n");
    const GroupCase cases[] = {
        {"CHFClBr", "5\n\nC 0 0 0\nH 0 0 1.09\nF 1.3 0 -0.45\nCl -0.85 1.45 -0.6\nBr -0.9 -1.55 -0.65\n", "C1", 1},
        {"HOCl", "3\n\nO 0 0 0\nH 0.96 0 0.1\nCl -0.4 0 1.65\n", "Cs", 1},
        {"H2O (H, O, H)", "3\n\nH 0.757 0 0.587\nO 0 0 0\nH -0.757 0 0.587\n", "C2v", 2},
        {"H2O2", "4\n\nO -0.7278 0 0\nO 0.7278 0 0\nH -0.9 0.8 0.4\nH 0.9 -0.8 0.4\n", "C2", 2},
        {"NH3 (H first)", "4\n\nH 0.9377 0 -0.381\nN 0 0 0\nH -0.46885 0.81207 -0.381\nH -0.46885 -0.81207 -0.381\n", "C3v", 3},
        {"BF3", "4\n\nB 0 0 0\nF 1.31 0 0\nF -0.655 1.13449 0\nF -0.655 -1.13449 0\n", "D3h", 6},
        {"CH4", "5\n\nC 0 0 0\nH 0.629 0.629 0.629\nH -0.629 -0.629 0.629\nH -0.629 0.629 -0.629\nH 0.629 -0.629 -0.629\n", "Td", 12},
        {"C2H4", "6\n\nC 0 0 0.6695\nC 0 0 -0.6695\nH 0 0.9289 1.2321\nH 0 -0.9289 1.2321\nH 0 0.9289 -1.2321\nH 0 -0.9289 -1.2321\n", "D2h", 4},
        {"allene", "7\n\nC 0 0 0\nC 0 0 1.31\nC 0 0 -1.31\nH 0.93 0 1.87\nH -0.93 0 1.87\nH 0 0.93 -1.87\nH 0 -0.93 -1.87\n", "D2d", 4},
        {"C2H6 staggered",
         "8\n\nC 0 0 0.765\nC 0 0 -0.765\nH 1.0 0 1.16\nH -0.5 0.866025 1.16\nH -0.5 -0.866025 1.16\nH -1.0 0 -1.16\n"
         "H 0.5 0.866025 -1.16\nH 0.5 -0.866025 -1.16\n",
         "D3d", 6},
        {"C2H6 eclipsed",
         "8\n\nC 0 0 0.765\nC 0 0 -0.765\nH 1.0 0 1.16\nH -0.5 0.866025 1.16\nH -0.5 -0.866025 1.16\nH 1.0 0 -1.16\n"
         "H -0.5 0.866025 -1.16\nH -0.5 -0.866025 -1.16\n",
         "D3h", 6},
        {"B(OH)3",
         "7\n\nB 0 0 0\nO 1.36 0 0\nO -0.68 1.177795 0\nO -0.68 -1.177795 0\nH 1.6 0.9 0\nH -1.579423 0.935641 0\n"
         "H -0.020577 -1.835641 0\n",
         "C3h", 3},
        {"C4H4 (square)",
         "8\n\nC 0.72 0.72 0\nC -0.72 0.72 0\nC -0.72 -0.72 0\nC 0.72 -0.72 0\nH 1.48 1.48 0\nH -1.48 1.48 0\nH -1.48 -1.48 0\n"
         "H 1.48 -1.48 0\n",
         "D4h", 8},
        {"C6H6",
         "12\n\nC 1.397 0 0\nC 0.6985 1.209838 0\nC -0.6985 1.209838 0\nC -1.397 0 0\nC -0.6985 -1.209838 0\n"
         "C 0.6985 -1.209838 0\nH 2.481 0 0\nH 1.2405 2.148609 0\nH -1.2405 2.148609 0\nH -2.481 0 0\n"
         "H -1.2405 -2.148609 0\nH 1.2405 -2.148609 0\n",
         "D6h", 12},
        {"SF6", "7\n\nS 0 0 0\nF 1.56 0 0\nF -1.56 0 0\nF 0 1.56 0\nF 0 -1.56 0\nF 0 0 1.56\nF 0 0 -1.56\n", "Oh", 24},
        {"CO2", "3\n\nO 0 0 1.16\nC 0 0 0\nO 0 0 -1.16\n", "D\xE2\x88\x9Eh", 2},
        {"HCN", "3\n\nH 0 0 -1.0655\nC 0 0 0\nN 0 0 1.1532\n", "C\xE2\x88\x9Ev", 1},
        {"C2H2 (tilted axis)", "4\n\nC 0.35 0.35 0.35\nC -0.35 -0.35 -0.35\nH 0.96 0.96 0.96\nH -0.96 -0.96 -0.96\n", "D\xE2\x88\x9Eh", 2},
    };
    for (const auto& c : cases) {
        const Molecule m = molecule_from_xyz(c.xyz);
        const PointGroup g = point_group(m);
        const bool ok = g.name == c.group && g.symmetry_number == c.sigma && symmetry_number(m) == c.sigma;
        ++g_checks;
        if (!ok) ++g_failures;
        std::printf("%s point group of %-20s %-6s (expected %s), symmetry number %d (expected %d)\n", ok ? "ok  " : "FAIL",
                    c.name, g.name.c_str(), c.group, g.symmetry_number, c.sigma);
    }
}

static void test_energies_against_pyscf(const json& ref) {
    std::printf("-- energies against PySCF\n");
    double worst_rhf = 0.0, worst_cc = 0.0;
    for (const auto& r : ref["energies"]) {
        const std::string name = r["molecule"], basis = r["basis"], method = r["method"];
        const json& s = ref["structures"][name];
        const auto t0 = Clock::now();
        const json out = run(input_for(s["xyz"], s["charge"], basis, method));
        if (!out["error"].is_null()) {
            ++g_checks;
            ++g_failures;
            std::printf("FAIL %s %s %s: %s\n", name.c_str(), basis.c_str(), method.c_str(), out["error"].get<std::string>().c_str());
            continue;
        }
        const double e = method == "RHF" ? out["rhf"]["totalEnergy"].get<double>() : out["coupledCluster"]["totalEnergy"].get<double>();
        const double d = std::fabs(e - r["energy"].get<double>());
        CHECK(out["basis"]["functions"] == r["basisFunctions"]);
        std::printf("     %-6s %-13s %-8s %18.10f  diff %.1e  %5.1f s\n", name.c_str(), basis.c_str(), method.c_str(), e, d, since(t0));
        (method == "RHF" ? worst_rhf : worst_cc) = std::max(method == "RHF" ? worst_rhf : worst_cc, d);
    }
    check_max("RHF energies of custom molecules (H–Kr, ions, one atom) vs PySCF", worst_rhf, 1e-8);
    check_max("CCSD(T) energies (HCl, CH4) vs PySCF", worst_cc, 1e-8);

    // HOCl: the analytic RHF gradient, and the Cartesian numerical gradient (energies)
    const json& s = ref["structures"]["HOCl"];
    const Molecule m = molecule_from_xyz(s["xyz"].get<std::string>());
    const BasisSet b = build_basis(m, "aug-cc-pVDZ", "");
    ElectronicOptions eo;
    eo.scf_energy_tol = 1e-12;
    eo.scf_gradient_tol = 1e-10;
    const GradientResult ga = run_gradient(m, b, "RHF", eo);
    double worst = 0.0;
    const auto& gr = ref["gradients"][0]["gradient"];
    for (std::size_t x = 0; x < ga.gradient.size(); ++x) worst = std::max(worst, std::fabs(ga.gradient[x] - gr[x].get<double>()));
    check_max("HOCl RHF/aug-cc-pVDZ analytic gradient vs PySCF", worst, 2e-8);
    const GradientResult gn = run_gradient(m, b, "RHF", eo, nullptr, ElectronicCallbacks(), true);
    worst = 0.0;
    for (std::size_t x = 0; x < ga.gradient.size(); ++x) worst = std::max(worst, std::fabs(ga.gradient[x] - gn.gradient[x]));
    CHECK(gn.energy_evaluations == 6 && contains(gn.description, "3 Cartesian directions"));
    // central differences, h = 0.005 bohr: truncation error h^2/6 times the third derivative, ~1e-5
    check_max("HOCl numerical gradient (6 energies, Cartesian directions) vs analytic", worst, 3e-5);

    // the numerical Hessian (Cartesian directions, analytic gradients) against the analytic one
    HessianOptions ho;
    const HessianResult ha = run_hessian(m, b, "RHF", eo, ho);
    ho.force_numerical = true;
    const HessianResult hn = run_hessian(m, b, "RHF", eo, ho);
    CHECK(ha.analytic && !hn.analytic && hn.gradient_evaluations == 6);
    const HarmonicAnalysis aa = harmonic_analysis(m, ha.hessian, &ha.dipole_derivatives);
    const HarmonicAnalysis an = harmonic_analysis(m, hn.hessian, &hn.dipole_derivatives);
    double wn = 0.0, wi = 0.0;
    for (std::size_t k = 0; k < aa.modes.size(); ++k) {
        wn = std::max(wn, std::fabs(aa.modes[k].wavenumber - an.modes[k].wavenumber));
        wi = std::max(wi, std::fabs(aa.modes[k].ir_intensity - an.modes[k].ir_intensity) / std::max(1.0, aa.modes[k].ir_intensity));
    }
    check_max("HOCl numerical Hessian (Cartesian directions, not at a minimum): wavenumbers vs analytic (cm-1)", wn, 0.2);
    check_max("HOCl numerical Hessian: IR intensities vs analytic (relative)", wi, 1e-3);
}

static void test_atom_order() {
    std::printf("-- the order of the atoms\n");
    // NH3 and H2O at the built-in geometries, from files with the hydrogens first
    const json nh3 = run({{"molecule", "NH3"}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}, {"printLevel", 0}});
    json atoms = json::array();
    const auto& a = nh3["molecule"]["atoms"];
    for (int i : {1, 2, 0, 3}) {
        const auto& p = a[i]["position"];
        atoms.push_back({{"symbol", a[i]["symbol"]}, {"x", p[0].get<double>() * kBohrInAngstrom},
                         {"y", p[1].get<double>() * kBohrInAngstrom}, {"z", p[2].get<double>() * kBohrInAngstrom}});
    }
    const json c = run({{"molecule", "custom"}, {"atoms", atoms}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"},
                        {"printLevel", 0}, {"job", "freq"}});
    const json b = run({{"molecule", "NH3"}, {"basis", "aug-cc-pVDZ"}, {"method", "RHF"}, {"printLevel", 0}, {"job", "freq"}});
    CHECK(c["error"].is_null() && b["error"].is_null());
    check_max("NH3 (H, H, N, H) RHF energy vs the built-in NH3",
              std::fabs(c["rhf"]["totalEnergy"].get<double>() - nh3["rhf"]["totalEnergy"].get<double>()), 1e-10);
    double w = 0.0;
    std::string sc, sb;
    for (std::size_t k = 0; k < 6; ++k) {
        w = std::max(w, std::fabs(c["frequencies"]["modes"][k]["wavenumber"].get<double>() -
                                  b["frequencies"]["modes"][k]["wavenumber"].get<double>()));
        sc += c["frequencies"]["modes"][k]["symmetry"].get<std::string>() + " ";
        sb += b["frequencies"]["modes"][k]["symmetry"].get<std::string>() + " ";
    }
    check_max("NH3 (H, H, N, H) analytic RHF wavenumbers vs the built-in NH3 (cm-1)", w, 1e-6);
    CHECK(sc == sb && c["frequencies"]["pointGroup"] == "C3v");
    CHECK(std::fabs(c["thermochemistry"]["totals"]["gibbs"].get<double>() - b["thermochemistry"]["totals"]["gibbs"].get<double>()) < 1e-9);
}

static void test_species() {
    std::printf("-- species of the normal modes\n");
    struct S {
        const char* xyz;
        const char* expected;
    };
    const S cases[] = {
        // C2H4 in Mulliken's axes (x perpendicular to the plane, z along C=C): 12 modes
        {"6\n\nC 0 0 0.6695\nC 0 0 -0.6695\nH 0 0.9289 1.2321\nH 0 -0.9289 1.2321\nH 0 0.9289 -1.2321\nH 0 -0.9289 -1.2321\n", ""},
        {"3\n\nO 0 0 1.16\nC 0 0 0\nO 0 0 -1.16\n", "Pi_u Pi_u Sigma_g+ Sigma_u+"},
    };
    for (const auto& c : cases) {
        const Molecule m = molecule_from_xyz(c.xyz);
        const BasisSet b = build_basis(m, "aug-cc-pVDZ", "");
        ElectronicOptions eo;
        const HessianResult h = run_hessian(m, b, "RHF", eo, HessianOptions());
        const HarmonicAnalysis an = harmonic_analysis(m, h.hessian, &h.dipole_derivatives);
        std::string s;
        bool all = true;
        for (const auto& md : an.modes) {
            s += (s.empty() ? "" : " ") + md.symmetry;
            if (md.symmetry.empty()) all = false;
        }
        std::printf("     %-8s %s\n", an.point_group.c_str(), s.c_str());
        CHECK(all);
        if (*c.expected) CHECK(s == c.expected);
        if (an.point_group == "D2h") {
            // ethylene: 3 Ag, 1 Au, 2 B1u, 1 B1g, 2 B2g, 1 B2u, ... the IR-active ones are the u species
            // except Au; the g species have zero intensity
            int ag = 0;
            for (const auto& md : an.modes) {
                if (md.symmetry == "Ag") ++ag;
                if (md.symmetry.back() == 'g') CHECK(md.ir_intensity < 1e-6);
            }
            CHECK(ag == 3);
        }
    }
}

static void test_job() {
    std::printf("-- a job through the C API\n");
    const std::string xyz = "3\nH$_2$S\nH 0.9616 0 0.9269\nS 0 0 0\nH -0.9616 0 0.9269\n";
    json in = input_for(xyz, 0, "aug-cc-pVDZ", "RHF");
    in["job"] = "opt+freq";
    in["printLevel"] = 2;
    const auto t0 = Clock::now();
    const json j = run(in);
    std::printf("     H2S RHF/aug-cc-pVDZ opt+freq: %.1f s\n", since(t0));
    CHECK(j["error"].is_null());
    CHECK(j["molecule"]["formula"] == "H2S" && j["molecule"]["title"] == "H$_2$S");
    CHECK(j["optimization"]["converged"] == true);
    CHECK(j["frequencies"]["pointGroup"] == "C2v" && j["thermochemistry"]["symmetryNumber"] == 2);
    const std::string out = j["output"];
    CHECK(contains(out, "H2S, 18 electrons, charge 0, closed-shell singlet"));
    CHECK(contains(out, "H$_2$S (3 atoms, Cartesian coordinates)"));
    CHECK(contains(out, "Normal termination"));
    CHECK(j["input"]["atoms"].size() == 3 && j["input"]["title"] == "H$_2$S");
    // the optimised atoms of the custom molecule are in "optimization.finalGeometry.atoms"
    CHECK(j["optimization"]["finalGeometry"]["atoms"].size() == 3);

    // relaxed scan of the H-S-H angle (atoms 1 2 3), three points
    json sc = input_for(xyz, 0, "aug-cc-pVDZ", "RHF");
    sc["job"] = "scan";
    sc["scan"] = "A 1 2 3 S 2 5.0";
    sc["scanStart"] = 90.0;
    const json s = run(sc);
    CHECK(s["error"].is_null());
    CHECK(s["scan"]["points"].size() == 3);
}

int main() {
    const std::string data = MACQC_TEST_DATA_DIR;
    try {
        test_structure_files();
        test_formula_and_checks();
        test_point_groups();
        const json ref = json::parse(read_file(data + "/custom_reference.json"));
        test_energies_against_pyscf(ref);
        test_atom_order();
        test_species();
        test_job();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
