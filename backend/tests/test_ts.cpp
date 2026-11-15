// Tests of the transition-state optimisation of ΨLab 7.0, stage 7.3 (CTest):
//  - Opt=TS searches against geomeTRIC's transition-state optimiser (transition=True, Hessian by
//    finite differences at the first step) with PySCF gradients (backend/tests/data/ts_reference.json
//    from reference/pyscf_ts_reference.py): trans and cis H2O2 and planar NH3 with RHF (CalcFC and
//    CalcAll), NH3 with CCSD(T) (numerical CalcFC Hessian from analytic gradients);
//  - ts+freq: exactly one imaginary frequency, its species (Au, A2, A2'') and the wavenumbers against
//    PySCF's analytic Hessian at geomeTRIC's transition state;
//  - an asymmetric start (R(O1-H3) != R(O2-H4)) reaches the symmetric trans transition state;
//  - the eigenvalue test: a search started at the minimum stops (0 negative eigenvalues), with the
//    test off P-RFO climbs from the minimum to the trans transition state; linear H2O (a second-order
//    saddle point) fails the test, and without it ts+freq ends with an error (two imaginary modes);
//  - ts+freq+sp: the single point (CCSD(T)/aug-cc-pVTZ) at the transition state equals an independent
//    calculation at that geometry, and the final energies are the single point plus the thermal
//    corrections of the frequency step;
//  - refusals (diatomics, CCSDT/CCSDTQ of four atoms), input errors, the estimate, the output file.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gradient.h"
#include "macqc/macqc_api.h"
#include "macqc/molecule.h"
#include "macqc/optimizer.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = 180.0 / kPi;

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

static json run_logged(const std::string& label, const json& input) {
    const auto t0 = Clock::now();
    const json j = run(input);
    std::printf("     %-66s %6.1f s%s\n", label.c_str(), since(t0),
                j["error"].is_null() ? "" : (" error: " + j["error"].get<std::string>()).c_str());
    return j;
}

static const json& ref_case(const json& ref, const std::string& label, const std::string& method) {
    for (const auto& c : ref["cases"])
        if (c["label"] == label && c["method"] == method) return c;
    throw std::runtime_error("no reference for " + label + " " + method);
}

// Largest differences of bond lengths (angstrom) and angles (degrees; |dihedral|) from the reference.
static std::pair<double, double> geometry_diff(const json& g, const json& rc) {
    double dr = 0.0, da = 0.0;
    for (const char* k : {"ROO", "R1", "R2", "R3"})
        if (g.contains(k) && rc.contains(k)) dr = std::max(dr, std::fabs(g[k].get<double>() * kBohrInAngstrom - rc[k].get<double>()));
    for (const char* k : {"angle1", "angle2", "angle3"})
        if (g.contains(k) && rc.contains(k)) da = std::max(da, std::fabs(g[k].get<double>() - rc[k].get<double>()));
    if (g.contains("dihedral") && rc.contains("dihedral") && rc["molecule"] == "H2O2")
        da = std::max(da, std::fabs(std::fabs(g["dihedral"].get<double>()) - std::fabs(rc["dihedral"].get<double>())));
    return {dr, da};
}

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    const auto t_all = Clock::now();
    const std::string data = MACQC_TEST_DATA_DIR;
    const json ref = json::parse(read_file(data + "/ts_reference.json"));

    // ---- transition states against geomeTRIC (transition=True)
    std::printf("Transition states against geomeTRIC (transition=True) and PySCF frequencies\n");
    struct TsRun {
        std::string label, molecule, method, hessian;
        json start;
        std::string species, group;
    };
    const std::vector<TsRun> runs = {
        {"H2O2 trans", "H2O2", "RHF", "calcfc", {{"dihedral", 160.0}}, "Au", "C2h"},
        {"H2O2 trans", "H2O2", "RHF", "calcall", {{"dihedral", 160.0}}, "Au", "C2h"},
        {"H2O2 cis", "H2O2", "RHF", "calcfc", {{"dihedral", 20.0}}, "A2", "C2v"},
        {"NH3 planar", "NH3", "RHF", "calcfc", {{"angle", 119.0}}, "A2''", "D3h"},
        {"NH3 planar", "NH3", "CCSD(T)", "calcfc", {{"angle", 119.0}}, "", ""},
    };
    int steps_fc = 0, steps_all = 0;
    for (const auto& tr : runs) {
        const json rc = ref_case(ref, tr.label, tr.method);
        const bool freq = tr.method == "RHF";
        json in = {{"molecule", tr.molecule}, {"basis", "aug-cc-pVDZ"}, {"method", tr.method},
                   {"job", freq ? "ts+freq" : "ts"}, {"tsHessian", tr.hessian}, {"optConvergence", "tight"}};
        for (auto it = tr.start.begin(); it != tr.start.end(); ++it) in[it.key()] = it.value();
        const std::string tag = "  " + tr.label + " " + tr.method + " (" + tr.hessian + ")";
        const json r = run_logged(tag.substr(2) + (freq ? ", ts+freq" : ", ts"), in);
        CHECK(r["error"].is_null());
        if (!r["error"].is_null()) continue;
        const json& o = r["optimization"];
        CHECK(o["converged"].get<bool>() && o["transitionState"].get<bool>());
        CHECK(o["negativeEigenvalues"].get<int>() == 1);
        std::printf("     %d steps, %d computed Hessian%s\n", o["steps"].get<int>(), o["hessianEvaluations"].get<int>(),
                    o["hessianEvaluations"].get<int>() == 1 ? "" : "s");
        if (tr.label == "H2O2 trans" && tr.hessian == "calcfc") steps_fc = o["steps"].get<int>();
        if (tr.label == "H2O2 trans" && tr.hessian == "calcall") {
            steps_all = o["steps"].get<int>();
            CHECK(o["hessianEvaluations"].get<int>() == steps_all);
        }
        if (tr.hessian == "calcfc") CHECK(o["hessianEvaluations"].get<int>() == 1);
        const auto [dr, da] = geometry_diff(o["finalGeometry"], rc);
        check_max(tag + ": energy (E_h)", std::fabs(o["finalEnergy"].get<double>() - rc["energy"].get<double>()), 2e-9);
        check_max(tag + ": bond lengths (angstrom)", dr, 4e-5);
        check_max(tag + ": angles, |dihedral| (degrees)", da, tr.molecule == "H2O2" ? 1e-2 : 5e-3);
        if (!freq) continue;
        const json& f = r["frequencies"];
        int imaginary = 0;
        for (const auto& m : f["modes"]) imaginary += m["imaginary"].get<bool>() ? 1 : 0;
        CHECK(imaginary == 1);
        CHECK(f["pointGroup"] == tr.group);
        CHECK(f["modes"][0]["symmetry"] == tr.species);
        double dw = 0.0;
        for (std::size_t k = 0; k < 6; ++k)
            dw = std::max(dw, std::fabs(f["modes"][k]["wavenumber"].get<double>() - rc["frequencies"][k].get<double>()));
        check_max(tag + ": wavenumbers vs PySCF at geomeTRIC's transition state (cm-1)", dw, 0.5);
        std::printf("     %s, %.2fi cm-1 (%s)\n", f["pointGroup"].get<std::string>().c_str(),
                    std::fabs(f["modes"][0]["wavenumber"].get<double>()), f["modes"][0]["symmetry"].get<std::string>().c_str());
        const std::string out = r["output"].get<std::string>();
        CHECK(contains(out, "---- Transition-state optimisation: HF/aug-cc-pVDZ"));
        CHECK(contains(out, "TS Optimisation step 1"));
        CHECK(contains(out, "P-RFO step: lambda_p"));
        CHECK(contains(out, "Transition state found"));
        CHECK(contains(out, "One imaginary frequency"));
        CHECK(contains(out, tr.hessian == "calcall" ? "computed Hessian (CalcAll)" : "Bofill update"));
    }
    CHECK(steps_all > 0 && steps_all < steps_fc);
    std::printf("     trans H2O2: CalcFC %d steps, CalcAll %d steps\n", steps_fc, steps_all);

    // ---- asymmetric start, through the optimiser
    std::printf("Asymmetric start\n");
    {
        Molecule m = make_molecule("H2O2", default_bond_length("H2O2"), default_angle("H2O2"), 0.0, 155.0);
        for (std::size_t k = 0; k < 3; ++k) m.atoms[2].position[k] += 0.04 * (m.atoms[2].position[k] - m.atoms[0].position[k]);
        OptimizationOptions oo;
        oo.tight = true;
        oo.transition_state = true;
        ElectronicOptions eo;
        eo.scf_energy_tol = 1e-12;
        eo.scf_gradient_tol = 1e-10;
        const auto t0 = Clock::now();
        const OptimizationResult r = optimize_geometry(m, build_basis(m, "aug-cc-pVDZ", default_basis_directory()), "RHF", eo, oo);
        std::printf("     start R(O1-H3) - R(O2-H4) = %.3f bohr, dihedral 155 degrees: %d steps, %.1f s\n",
                    simple_internals(m).values[1] - simple_internals(m).values[2], r.steps, since(t0));
        CHECK(r.converged && r.final_negative_eigenvalues == 1);
        const SimpleInternals si = simple_internals(r.final_geometry);
        const json rc = ref_case(ref, "H2O2 trans", "RHF");
        check_max("  energy vs the trans transition state (E_h)", std::fabs(r.final_energy - rc["energy"].get<double>()), 2e-9);
        check_max("  symmetric at the end: |R(O1-H3) - R(O2-H4)| (angstrom)",
                  std::fabs(si.values[1] - si.values[2]) * kBohrInAngstrom, 4e-5);
        check_max("  |dihedral| - 180 (degrees)", std::fabs(std::fabs(si.values[5] * kDeg) - 180.0), 1e-2);
    }

    // ---- the eigenvalue test
    std::printf("Eigenvalue test\n");
    {
        const json a = run_logged("H2O2 HF ts from the minimum region (dihedral 113.7)", {{"molecule", "H2O2"}, {"job", "ts"}});
        CHECK(!a["error"].is_null() && contains(a["error"].get<std::string>(), "eigenvalue test failed (0 negative"));
        CHECK(a["optimization"]["eigenTestFailed"].get<bool>() && a["optimization"]["steps"].get<int>() == 1);
        const json b = run_logged("  the same with eigenTest false (NoEigenTest), ts+freq",
                                  {{"molecule", "H2O2"}, {"job", "ts+freq"}, {"eigenTest", false}, {"maxSteps", 60}});
        CHECK(b["error"].is_null());
        if (b["error"].is_null()) {
            CHECK(b["frequencies"]["pointGroup"] == "C2h" && b["frequencies"]["modes"][0]["symmetry"] == "Au");
            check_max("  climbed to the trans transition state (E_h)",
                      std::fabs(b["optimization"]["finalEnergy"].get<double>() - ref_case(ref, "H2O2 trans", "RHF")["energy"].get<double>()),
                      1e-7);
        }
        const json c = run_logged("H2O at 180 degrees (second-order saddle point), ts", {{"molecule", "H2O"}, {"job", "ts"}, {"angle", 180.0}});
        CHECK(!c["error"].is_null() && contains(c["error"].get<std::string>(), "(2 negative eigenvalues"));
        const json d = run_logged("  the same with eigenTest false, ts+freq",
                                  {{"molecule", "H2O"}, {"job", "ts+freq"}, {"angle", 180.0}, {"eigenTest", false}});
        CHECK(!d["error"].is_null() && contains(d["error"].get<std::string>(), "a saddle point of order 2, not a transition state"));
        CHECK(d["job"]["termination"] == "error");
    }

    // ---- ts+freq+sp
    std::printf("Transition state + frequencies + single-point energy\n");
    {
        const json r = run_logged("NH3 HF/aug-cc-pVDZ ts+freq, CCSD(T)/aug-cc-pVTZ single point",
                                  {{"molecule", "NH3"}, {"method", "HF"}, {"job", "ts+freq+sp"}, {"angle", 119.0},
                                   {"spMethod", "CCSD(T)"}, {"spBasis", "aug-cc-pVTZ"}, {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        if (r["error"].is_null()) {
            const json& j = r["job"];
            CHECK(j["type"] == "ts+freq+sp" && j["termination"] == "normal");
            CHECK(j["estimate"]["steps"].size() == 3 && j["estimate"]["steps"][2]["name"] == "Single-point energy");
            CHECK(j["estimate"]["steps"][2]["how"] == "single point at the transition state");
            CHECK(j["finalEnergyLabel"] == "CCSD(T)/aug-cc-pVTZ");
            CHECK(r["frequencies"]["pointGroup"] == "D3h" && r["frequencies"]["modes"][0]["symmetry"] == "A2''");
            // the single point again, independently, at the transition state
            Molecule m = make_molecule("NH3", 1.9, 120.0);
            const json& atoms = r["optimization"]["finalGeometry"]["atoms"];
            for (std::size_t a = 0; a < 4; ++a)
                for (std::size_t k = 0; k < 3; ++k) m.atoms[a].position[k] = atoms[a]["position"][k].get<double>();
            ElectronicOptions eo;
            eo.scf_energy_tol = 1e-12;
            eo.scf_gradient_tol = 1e-10;
            eo.cc_amplitude_tol = 1e-10;
            eo.cc_energy_tol = 1e-12;
            const SinglePoint sp = run_single_point(m, build_basis(m, "aug-cc-pVTZ", default_basis_directory()), "CCSD(T)", eo);
            const double e_sp = r["singlePoint"]["energy"].get<double>();
            check_max("  single point at the transition state vs an independent CCSD(T)/aug-cc-pVTZ energy (E_h)",
                      std::fabs(e_sp - sp.energy), 2e-9);
            check_max("  final energy = the single point (E_h)", std::fabs(j["finalEnergy"].get<double>() - e_sp), 1e-12);
            const json& t = r["thermochemistry"];
            const double kj = 1.0 / (4.3597447222071e-18 * 6.02214076e23);  // J/mol -> E_h (CODATA 2018)
            check_max("  G = single point + thermal correction to G (E_h)",
                      std::fabs(t["totals"]["gibbs"].get<double>() - (e_sp + t["corrections"]["gibbs"].get<double>() * kj)), 1e-10);
            check_max("  H = single point + thermal correction to H (E_h)",
                      std::fabs(t["totals"]["enthalpy"].get<double>() - (e_sp + t["corrections"]["enthalpy"].get<double>() * kj)),
                      1e-10);
            const std::string out = r["output"].get<std::string>();
            CHECK(contains(out, "---- Single-point energy at the transition state: CCSD(T)/aug-cc-pVTZ"));
            CHECK(contains(out, "One imaginary frequency"));
            CHECK(out.find("One imaginary frequency") < out.find("---- Single-point energy at the transition state"));
        }
        // a failed transition-state check stops before the single point
        const json d = run_logged("linear H2O, eigenTest false, ts+freq+sp",
                                  {{"molecule", "H2O"}, {"job", "ts+freq+sp"}, {"angle", 180.0}, {"eigenTest", false},
                                   {"spMethod", "CCSD"}});
        CHECK(!d["error"].is_null() && contains(d["error"].get<std::string>(), "a saddle point of order 2"));
        CHECK(d["singlePoint"].is_null());
    }

    // ---- refusals, input errors, estimate
    std::printf("Refusals, input errors, estimate\n");
    {
        auto err = [&](const json& in) {
            const json j = run(in);
            return j["error"].is_null() ? std::string() : j["error"].get<std::string>();
        };
        CHECK(contains(err({{"molecule", "LiH"}, {"job", "ts"}}), "A diatomic molecule has no transition state"));
        CHECK(contains(err({{"molecule", "NH3"}, {"method", "CCSDT"}, {"job", "ts"}}), "CCSDT transition-state searches for NH3"));
        CHECK(contains(err({{"molecule", "H2O"}, {"job", "ts"}, {"tsHessian", "readfc"}}), "tsHessian must be calcfc or calcall"));
        const json e = run({{"molecule", "NH3"}, {"method", "CCSD(T)"}, {"job", "ts+freq"}, {"estimateOnly", true}});
        CHECK(e["error"].is_null());
        const json& st = e["job"]["estimate"]["steps"];
        CHECK(st[0]["name"] == "Transition-state optimisation" &&
              contains(st[0]["how"].get<std::string>(), "Hessian computed at the first step (CalcFC)"));
        CHECK(st[1]["name"] == "Hessian");
        const json ea = run({{"molecule", "NH3"}, {"method", "CCSD(T)"}, {"job", "ts"}, {"tsHessian", "calcall"}, {"estimateOnly", true}});
        CHECK(contains(ea["job"]["estimate"]["steps"][0]["how"].get<std::string>(), "at every step (CalcAll)"));
        CHECK(ea["job"]["estimate"]["totalSeconds"].get<double>() > e["job"]["estimate"]["steps"][0]["seconds"].get<double>());
    }

    std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t_all));
    return g_failures == 0 ? 0 : 1;
}
