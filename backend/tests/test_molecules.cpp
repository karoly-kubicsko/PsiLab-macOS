// Tests of the molecules of ΨLab 7.0, stage 7.2: H2O2 and NH3 (CTest).
//  - geometry: the parameters of the built molecules (IUPAC sign of the H2O2 dihedral, planar NH3
//    at 120 degrees), the natural internal coordinates (Z-matrix round trip, B rows against finite
//    differences), the improper torsions of NH3 (complete set also at the planar geometry);
//  - against PySCF (backend/tests/data/molecules_reference.json from
//    reference/pyscf_molecules_reference.py): RHF energies and analytic gradients (aug-cc-pVDZ,
//    aug-cc-pVTZ, aug-cc-pCVDZ), CCSD(T)/aug-cc-pVDZ energies and the NH3 CCSD(T) gradient, tight
//    optimisations (geomeTRIC), RHF wavenumbers at the minima;
//  - symmetry: point groups, species and symmetry numbers at the minima (C2 H2O2, C3v NH3) and at
//    the planar structures from backend/tests/data/scan_reference.json (trans C2h and cis C2v H2O2,
//    D3h NH3: one imaginary mode each, against PySCF's wavenumbers);
//  - numerical gradients (energies) and Hessians (gradients) against the analytic ones;
//  - relaxed scans: the H2O2 torsion (RHF 0-180 degrees; CCSD(T) frozen at 180) and the NH3 umbrella
//    (improper dihedral H2-N1-H3-H4 to -180 degrees, planar) against constrained minima from
//    PySCF + geomeTRIC; refusals (Full CI, CCSDT frequencies of four atoms) and input errors.
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
#include "macqc/internal_coordinates.h"
#include "macqc/macqc_api.h"
#include "macqc/molecule.h"
#include "macqc/optimizer.h"
#include "macqc/symmetry.h"
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
    std::printf("     %-62s %6.1f s%s\n", label.c_str(), since(t0),
                j["error"].is_null() ? "" : (" error: " + j["error"].get<std::string>()).c_str());
    return j;
}

static ElectronicOptions tight() {
    ElectronicOptions o;
    o.scf_energy_tol = 1e-12;
    o.scf_gradient_tol = 1e-10;
    o.cc_amplitude_tol = 1e-10;
    o.cc_energy_tol = 1e-12;
    return o;
}

static Molecule default_molecule(const std::string& name) {
    return make_molecule(name, default_bond_length(name), default_angle(name));
}

static Molecule from_json_geometry(const std::string& name, const json& g) {
    Molecule m = default_molecule(name);
    for (std::size_t a = 0; a < m.atoms.size(); ++a)
        for (std::size_t k = 0; k < 3; ++k) m.atoms[a].position[k] = g[a][k].get<double>();
    return m;
}

static double grad_diff(const std::vector<double>& g, const json& ref) {
    double d = 0.0;
    for (std::size_t a = 0; a < ref.size(); ++a)
        for (std::size_t k = 0; k < 3; ++k) d = std::max(d, std::fabs(g[3 * a + k] - ref[a][k].get<double>()));
    return d;
}

static std::string species(const HarmonicAnalysis& h) {
    std::string s;
    for (const auto& m : h.modes) s += (s.empty() ? "" : " ") + m.symmetry;
    return s;
}

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    const auto t_all = Clock::now();
    const std::string data = MACQC_TEST_DATA_DIR;
    const json ref = json::parse(read_file(data + "/molecules_reference.json"));
    const json scanref = json::parse(read_file(data + "/scan_reference.json"));
    const std::string dir = default_basis_directory();

    // ---- geometry and internal coordinates
    std::printf("Geometry and internal coordinates\n");
    {
        const Molecule p = default_molecule("H2O2");
        const SimpleInternals si = simple_internals(p);
        CHECK(si.names.size() == 6 && si.names[5] == "D(H3–O1–O2–H4)");
        check_max("H2O2 default: R(O–O) 1.4556, R(O–H) 0.9670 angstrom",
                  std::max(std::fabs(si.values[0] * kBohrInAngstrom - 1.4556),
                           std::max(std::fabs(si.values[1] * kBohrInAngstrom - 0.9670),
                                    std::fabs(si.values[2] * kBohrInAngstrom - 0.9670))),
                  1e-12);
        check_max("  angles O–O–H 102.32, dihedral H–O–O–H +113.70 degrees (IUPAC sign)",
                  std::max(std::fabs(si.values[3] * kDeg - 102.32),
                           std::max(std::fabs(si.values[4] * kDeg - 102.32), std::fabs(si.values[5] * kDeg - 113.70))),
                  1e-10);
        const Molecule a = default_molecule("NH3");
        const SimpleInternals sa = simple_internals(a);
        CHECK(sa.names[5] == "D(H4–N1–H2–H3)");
        check_max("NH3 default: R(N–H) 1.012 angstrom, H–N–H 106.7 degrees",
                  std::max(std::fabs(sa.values[0] * kBohrInAngstrom - 1.012), std::fabs(sa.values[3] * kDeg - 106.7)), 1e-10);
        const Molecule flat = make_molecule("NH3", default_bond_length("NH3"), 120.0);
        double zmax = 0.0;
        for (const auto& at : flat.atoms) zmax = std::max(zmax, std::fabs(at.position[2]));
        check_max("NH3 with H–N–H 120 degrees is planar (|z|, bohr)", zmax, 1e-12);
        // Z-matrix round trips (also distorted, and planar NH3)
        double worst = 0.0;
        for (const Molecule& m0 : {p, a, flat, make_molecule("H2O2", 1.9, 95.0, 2.9, -150.0)}) {
            Molecule m = m0;
            m.atoms[2].position[0] += 0.07;
            m.atoms[3].position[2] -= 0.05;
            const SimpleInternals s0 = simple_internals(m);
            const Molecule r = geometry_from_internals(m, s0.values);
            for (std::size_t i = 0; i < m.atoms.size(); ++i)
                for (std::size_t k = 0; k < 3; ++k) worst = std::max(worst, std::fabs(r.atoms[i].position[k] - m.atoms[i].position[k]));
        }
        check_max("geometry_from_internals(simple_internals(x)) = x (H2O2, NH3, planar NH3; bohr)", worst, 1e-12);
        // B rows of the natural coordinates against central differences
        double wb = 0.0;
        for (const Molecule& m : {p, a, flat}) {
            const SimpleInternals s0 = simple_internals(m);
            for (std::size_t x = 0; x < 12; ++x) {
                Molecule mp = m, mm = m;
                mp.atoms[x / 3].position[x % 3] += 1e-5;
                mm.atoms[x / 3].position[x % 3] -= 1e-5;
                const SimpleInternals sp = simple_internals(mp), sm = simple_internals(mm);
                for (std::size_t i = 0; i < 6; ++i)
                    wb = std::max(wb, std::fabs(wrap_angle(sp.values[i] - sm.values[i]) / 2e-5 - s0.B[i][x]));
            }
        }
        check_max("B rows of the natural coordinates vs central differences", wb, 1e-8);
        // redundant set of NH3: 3 stretches, 3 bends, 3 improper torsions; rank 6 also when planar
        const InternalCoordinates ic(flat);
        CHECK(ic.size() == 9);
        int impropers = 0;
        for (const auto& c : ic.coordinates())
            if (c.type == InternalType::Dihedral) ++impropers;
        CHECK(impropers == 3);
        const Matrix Bm = ic.wilson_b(flat);
        Matrix G(ic.size(), ic.size());
        for (std::size_t i = 0; i < ic.size(); ++i)
            for (std::size_t j = 0; j < ic.size(); ++j)
                for (std::size_t x = 0; x < 12; ++x) G(i, j) += Bm(i, x) * Bm(j, x);
        int rank = 0;
        generalized_inverse(G, 1e-8, &rank);
        CHECK(rank == 6);
        const InternalCoordinates bends_only = InternalCoordinates(make_molecule("H2O", 1.8, 104.5));
        CHECK(bends_only.size() == 3);  // no torsions for H2O
    }

    // ---- energies and gradients against PySCF
    std::printf("RHF and CCSD(T) energies and gradients against PySCF\n");
    for (const auto& p : ref["points"]) {
        const std::string name = p["molecule"], basis = p["basis"], method = p["method"];
        const Molecule m = default_molecule(name);
        const BasisSet b = build_basis(m, basis, dir);
        CHECK(b.n_functions == p["nbf"].get<int>());
        const auto t0 = Clock::now();
        const bool grad = !p["gradient"].is_null();
        if (grad) {
            const GradientResult g = run_gradient(m, b, method, tight());
            CHECK(g.analytic);
            const std::string tag = name + " " + method + "/" + basis;
            check_max(tag + ": energy (E_h)", std::fabs(g.point.energy - p["energy"].get<double>()), 1e-9);
            check_max(tag + ": analytic gradient (E_h/bohr)", grad_diff(g.gradient, p["gradient"]), method == "RHF" ? 2e-9 : 2e-8);
        } else {
            const SinglePoint sp = run_single_point(m, b, method, tight());
            check_max(name + " " + method + "/" + basis + ": energy (E_h)", std::fabs(sp.energy - p["energy"].get<double>()), 1e-9);
        }
        std::printf("     (%.1f s)\n", since(t0));
    }

    // ---- optimisations and frequencies at the minima
    std::printf("Optimisations (tight) against PySCF + geomeTRIC; frequencies at the RHF minima\n");
    for (const auto& o : ref["optimizations"]) {
        const std::string name = o["molecule"], method = o["method"];
        const json r = run_logged(name + " " + method + "/aug-cc-pVDZ opt" + (method == "RHF" ? "+freq" : ""),
                                  {{"molecule", name}, {"basis", "aug-cc-pVDZ"}, {"method", method},
                                   {"job", method == "RHF" ? "opt+freq" : "opt"}, {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        if (!r["error"].is_null()) continue;
        const json& g = r["optimization"]["finalGeometry"];
        double dr = 0.0, da = 0.0;
        const std::string tag = "  " + name + " " + method;
        if (name == "H2O2") {
            dr = std::max({std::fabs(g["ROO"].get<double>() * kBohrInAngstrom - o["ROO"].get<double>()),
                           std::fabs(g["R1"].get<double>() * kBohrInAngstrom - o["R1"].get<double>()),
                           std::fabs(g["R2"].get<double>() * kBohrInAngstrom - o["R2"].get<double>())});
            da = std::max({std::fabs(g["angle1"].get<double>() - o["angle1"].get<double>()),
                           std::fabs(g["angle2"].get<double>() - o["angle2"].get<double>()),
                           std::fabs(std::fabs(g["dihedral"].get<double>()) - std::fabs(o["dihedral"].get<double>()))});
        } else {
            for (const char* k : {"R1", "R2", "R3"})
                dr = std::max(dr, std::fabs(g[k].get<double>() * kBohrInAngstrom - o[k].get<double>()));
            for (const char* k : {"angle1", "angle2", "angle3"}) da = std::max(da, std::fabs(g[k].get<double>() - o[k].get<double>()));
        }
        check_max(tag + ": bond lengths (angstrom)", dr, 2e-5);
        check_max(tag + ": angles and dihedral (degrees)", da, 3e-3);
        check_max(tag + ": energy (E_h)", std::fabs(r["optimization"]["finalEnergy"].get<double>() - o["energy"].get<double>()), 1e-9);
        if (method != "RHF") continue;
        const json& f = r["frequencies"];
        double dw = 0.0;
        for (std::size_t k = 0; k < 6; ++k)
            dw = std::max(dw, std::fabs(f["modes"][k]["wavenumber"].get<double>() - o["frequencies"][k].get<double>()));
        check_max(tag + ": wavenumbers at ΨLab's minimum vs PySCF's (cm-1)", dw, 0.1);
        std::string sp;
        for (const auto& md : f["modes"]) sp += (sp.empty() ? "" : " ") + md["symmetry"].get<std::string>();
        std::printf("     point group %s, species %s, symmetry number %d\n", f["pointGroup"].get<std::string>().c_str(),
                    sp.c_str(), r["thermochemistry"]["symmetryNumber"].get<int>());
        if (name == "H2O2") {
            // the two O-H stretches lie 0.4 cm-1 apart, B below A
            CHECK(f["pointGroup"] == "C2" && sp == "A A B A B A");
            CHECK(r["thermochemistry"]["symmetryNumber"].get<int>() == 2);
        } else {
            CHECK(f["pointGroup"] == "C3v" && sp == "A1 E E A1 E E");
            CHECK(r["thermochemistry"]["symmetryNumber"].get<int>() == 3);
        }
        int imaginary = 0;
        for (const auto& md : f["modes"]) imaginary += md["imaginary"].get<bool>() ? 1 : 0;
        CHECK(imaginary == 0);
    }

    // ---- planar transition structures: one imaginary mode each
    std::printf("Planar structures (transition states by symmetry): frequencies, species, point groups\n");
    for (const auto& c : scanref["cases"]) {
        if (c["method"] != "RHF") continue;
        for (const auto& pt : c["points"]) {
            if (!pt.contains("frequencies")) continue;
            const std::string name = c["molecule"];
            const Molecule m = from_json_geometry(name, pt["geometry"]);
            const BasisSet b = build_basis(m, "aug-cc-pVDZ", dir);
            const HessianResult h = run_hessian(m, b, "RHF", tight());
            const HarmonicAnalysis a = harmonic_analysis(m, h.hessian, &h.dipole_derivatives);
            double dw = 0.0;
            int imaginary = 0;
            for (std::size_t k = 0; k < 6; ++k) {
                dw = std::max(dw, std::fabs(a.modes[k].wavenumber - pt["frequencies"][k].get<double>()));
                imaginary += a.modes[k].imaginary ? 1 : 0;
            }
            const double v = pt["value"].get<double>();
            const std::string tag = "  " + name + " " + c["coordinate"].get<std::string>() + " = " + std::to_string(static_cast<int>(v));
            check_max(tag + ": wavenumbers vs PySCF (cm-1)", dw, 0.05);
            CHECK(imaginary == 1);
            std::printf("     %s: point group %s, %.2fi cm-1 (%s), species %s, symmetry number %d\n", tag.c_str(),
                        a.point_group.c_str(), std::fabs(a.modes[0].wavenumber), a.modes[0].symmetry.c_str(),
                        species(a).c_str(), symmetry_number(m));
            if (name == "H2O2" && v == 180.0) CHECK(a.point_group == "C2h" && species(a) == "Au Ag Bu Ag Ag Bu");
            if (name == "H2O2" && v == 0.0) CHECK(a.point_group == "C2v" && a.modes[0].symmetry == "A2");
            if (name == "NH3") {
                CHECK(a.point_group == "D3h" && species(a) == "A2'' E' E' A1' E' E'");
                CHECK(symmetry_number(m) == 6);
            }
        }
    }

    // ---- numerical gradients and Hessians against the analytic ones
    std::printf("Numerical gradients (energies) and Hessians (gradients) against the analytic ones\n");
    {
        for (const std::string name : {"H2O2", "NH3"}) {
            const Molecule m = default_molecule(name);
            const BasisSet b = build_basis(m, "aug-cc-pVDZ", dir);
            const GradientResult ga = run_gradient(m, b, "RHF", tight());
            const GradientResult gn = run_gradient(m, b, "RHF", tight(), nullptr, ElectronicCallbacks(), true);
            double d = 0.0;
            for (std::size_t x = 0; x < 12; ++x) d = std::max(d, std::fabs(ga.gradient[x] - gn.gradient[x]));
            check_max(name + " RHF: numerical gradient (" + std::to_string(gn.energy_evaluations) + " energies) vs analytic", d, 1e-5);
            CHECK(gn.energy_evaluations == (name == "H2O2" ? 8 : 12));
            HessianOptions hn;
            hn.force_numerical = true;
            const HessianResult ha = run_hessian(m, b, "RHF", tight());
            const HessianResult hg = run_hessian(m, b, "RHF", tight(), hn);
            CHECK(hg.gradient_evaluations == 12);
            const HarmonicAnalysis A = harmonic_analysis(m, ha.hessian, &ha.dipole_derivatives);
            const HarmonicAnalysis N = harmonic_analysis(m, hg.hessian, &hg.dipole_derivatives);
            double dw = 0.0, di = 0.0;
            for (std::size_t k = 0; k < 6; ++k) {
                dw = std::max(dw, std::fabs(A.modes[k].wavenumber - N.modes[k].wavenumber));
                di = std::max(di, std::fabs(A.modes[k].ir_intensity - N.modes[k].ir_intensity));
            }
            check_max(name + " RHF: numerical Hessian (12 gradients) vs analytic, wavenumbers (cm-1)", dw, 0.1);
            check_max(name + " RHF:   IR intensities (km/mol)", di, 0.01);
            CHECK(species(A) == species(N));
        }
    }

    // ---- relaxed scans against constrained minima
    std::printf("Relaxed scans against constrained minima from PySCF + geomeTRIC\n");
    struct ScanRun {
        std::string molecule, method, line, coordinate;
        double start;
    };
    const std::vector<ScanRun> runs = {{"H2O2", "RHF", "3 1 2 4 S 6 30", "D 3 1 2 4", 0.0},
                                       {"NH3", "RHF", "2 1 3 4 S 4 -10", "D 2 1 3 4", -140.0}};
    for (const auto& sr : runs) {
        const json* rc = nullptr;
        for (const auto& c : scanref["cases"])
            if (c["molecule"] == sr.molecule && c["method"] == sr.method && c["coordinate"] == sr.coordinate) rc = &c;
        CHECK(rc != nullptr);
        if (!rc) continue;
        const json r = run_logged(sr.molecule + " " + sr.method + "/aug-cc-pVDZ scan " + sr.line,
                                  {{"molecule", sr.molecule}, {"basis", "aug-cc-pVDZ"}, {"method", sr.method}, {"job", "scan"},
                                   {"scan", sr.line}, {"scanStart", sr.start}, {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        if (!r["error"].is_null()) continue;
        CHECK(r["scan"]["converged"].get<bool>());
        for (const auto& rp : (*rc)["points"]) {
            const double v = rp["value"].get<double>();
            const json* p = nullptr;
            for (const auto& q : r["scan"]["points"])
                if (std::fabs(q["value"].get<double>() - v) < 1e-8) p = &q;
            if (!p) continue;
            const std::string tag = "  " + sr.molecule + " " + sr.method + " " + sr.coordinate + " = " + std::to_string(static_cast<int>(v));
            check_max(tag + ": energy (E_h)", std::fabs((*p)["energy"].get<double>() - rp["energy"].get<double>()), 2e-9);
            const json& g = (*p)["geometry"];
            double dr = 0.0, da = 0.0;
            for (const char* k : {"R1", "R2", "R3", "ROO"})
                if (g.contains(k) && rp.contains(k)) dr = std::max(dr, std::fabs(g[k].get<double>() * kBohrInAngstrom - rp[k].get<double>()));
            for (const char* k : {"angle1", "angle2", "angle3"})
                if (g.contains(k) && rp.contains(k)) da = std::max(da, std::fabs(g[k].get<double>() - rp[k].get<double>()));
            check_max(tag + ": bond lengths (angstrom)", dr, 4e-5);
            // the constrained valley of NH3 at -140 degrees is flat (8e-10 E_h for 5e-3 degrees)
            check_max(tag + ": angles (degrees)", da, sr.molecule == "NH3" ? 1e-2 : 5e-3);
            check_max(tag + ": frozen coordinate at its value", std::fabs((*p)["measured"].get<double>() - v), 1e-6);
        }
        if (sr.molecule == "NH3") {
            const json& last = r["scan"]["points"][4];
            check_max("  NH3 at -180 degrees is planar: 360 - sum of the H–N–H angles (degrees)",
                      std::fabs(last["geometry"]["pyramidal"].get<double>()), 1e-6);
        }
    }

    // CCSD(T): one constrained optimisation (trans H2O2) through the optimiser itself
    {
        const json* rp = nullptr;
        for (const auto& c : scanref["cases"])
            if (c["molecule"] == "H2O2" && c["method"] == "CCSD(T)") rp = &c["points"][0];
        CHECK(rp != nullptr);
        if (rp) {
            const auto t0 = Clock::now();
            const Molecule m = make_molecule("H2O2", default_bond_length("H2O2"), default_angle("H2O2"), 0.0, 180.0);
            OptimizationOptions oo;
            oo.tight = true;
            CoordinateConstraint con;
            con.coordinate.type = InternalType::Dihedral;
            con.coordinate.atoms = {2, 0, 1, 3};
            con.target = kPi;
            oo.constraints = {con};
            const OptimizationResult r = optimize_geometry(m, build_basis(m, "aug-cc-pVDZ", dir), "CCSD(T)", tight(), oo);
            std::printf("     H2O2 CCSD(T)/aug-cc-pVDZ, D(H3–O1–O2–H4) frozen at 180 degrees: %d steps, %.1f s\n", r.steps, since(t0));
            CHECK(r.converged);
            const SimpleInternals si = simple_internals(r.final_geometry);
            check_max("  H2O2 CCSD(T) D = 180: energy (E_h)", std::fabs(r.final_energy - (*rp)["energy"].get<double>()), 2e-9);
            check_max("  H2O2 CCSD(T) D = 180: R(O–O), R(O–H) (angstrom)",
                      std::max({std::fabs(si.values[0] * kBohrInAngstrom - (*rp)["ROO"].get<double>()),
                                std::fabs(si.values[1] * kBohrInAngstrom - (*rp)["R1"].get<double>()),
                                std::fabs(si.values[2] * kBohrInAngstrom - (*rp)["R2"].get<double>())}),
                      4e-5);
            check_max("  H2O2 CCSD(T) D = 180: angles O–O–H (degrees)",
                      std::max(std::fabs(si.values[3] * kDeg - (*rp)["angle1"].get<double>()),
                               std::fabs(si.values[4] * kDeg - (*rp)["angle2"].get<double>())),
                      5e-3);
        }
    }

    // ---- HeH+ with correlation-consistent basis sets (7.4): a general molecule, charge +1
    std::printf("HeH+ with correlation-consistent basis sets against PySCF\n");
    if (ref.contains("HeH+")) {
        const json& h = ref["HeH+"];
        for (const auto& p : h["points"]) {
            const std::string basis = p["basis"], method = p["method"];
            const Molecule m = make_molecule("HeH+", 1.4632, 0.0);
            CHECK(m.charge == 1 && m.electrons() == 2);
            const BasisSet b = build_basis(m, basis, dir);
            CHECK(b.n_functions == p["nbf"].get<int>());
            const std::string tag = "  HeH+ " + method + "/" + basis;
            if (method == "RHF") {
                const GradientResult g = run_gradient(m, b, "RHF", tight());
                check_max(tag + ": energy (E_h)", std::fabs(g.point.energy - p["energy"].get<double>()), 1e-9);
                check_max(tag + ": analytic gradient (E_h/bohr)", grad_diff(g.gradient, p["gradient"]), 2e-9);
            } else {
                const SinglePoint sp = run_single_point(m, b, method, tight());
                check_max(tag + ": energy (E_h)", std::fabs(sp.energy - p["energy"].get<double>()), 1e-9);
            }
        }
        const json& o = h["optimization"];
        for (const std::string method : {"CCSD", "FCI"}) {
            const json r = run_logged("HeH+ " + method + "/aug-cc-pVTZ opt+freq (two electrons: exact)",
                                      {{"molecule", "HeH+"}, {"basis", "aug-cc-pVTZ"}, {"method", method},
                                       {"job", "opt+freq"}, {"optConvergence", "tight"}});
            CHECK(r["error"].is_null());
            if (!r["error"].is_null()) continue;
            CHECK(r["molecule"]["charge"].get<int>() == 1);
            check_max("  HeH+ " + method + ": r_e vs PySCF + geomeTRIC CCSD (angstrom)",
                      std::fabs(r["optimization"]["finalGeometry"]["R"].get<double>() * kBohrInAngstrom - o["R"].get<double>()), 2e-5);
            check_max("  HeH+ " + method + ": energy at r_e (E_h)",
                      std::fabs(r["optimization"]["finalEnergy"].get<double>() - o["energy"].get<double>()), 1e-9);
            check_max("  HeH+ " + method + ": harmonic wavenumber vs PySCF energies' second derivative (cm-1)",
                      std::fabs(r["frequencies"]["modes"][0]["wavenumber"].get<double>() - o["wavenumber"].get<double>()), 0.5);
            CHECK(r["frequencies"]["modes"][0]["symmetry"] == "Sigma+" && r["thermochemistry"]["symmetryNumber"].get<int>() == 1);
        }
        // the Szabo-Ostlund path is unchanged
        const json legacy = run({{"molecule", "HeH+"}});
        CHECK(legacy["error"].is_null() && !legacy["hf"].is_null());
        check_max("  HeH+ without a basis set: the Szabo–Ostlund STO-3G energy (E_h)",
                  std::fabs(legacy["hf"]["totalEnergy"].get<double>() - (-2.860661992)), 1e-9);
        CHECK(run({{"molecule", "HeH+"}, {"basis", "STO-1G"}})["input"]["nGaussians"].get<int>() == 1);
        CHECK(contains(run({{"molecule", "HeH+"}, {"job", "opt"}})["error"].get<std::string>(), "choose a Gaussian basis set"));
        CHECK(contains(run({{"molecule", "HeH+"}, {"basis", "6-31G**"}})["error"].get<std::string>(), "STO-1G, STO-2G, STO-3G"));
        CHECK(contains(run({{"molecule", "HeH+"}, {"basis", "aug-cc-pVDZ"}, {"job", "ts"}})["error"].get<std::string>(),
                       "no transition state"));
    }

    // ---- refusals, input errors
    std::printf("Refusals and input errors\n");
    {
        auto err = [&](const json& in) {
            const json j = run(in);
            return j["error"].is_null() ? std::string() : j["error"].get<std::string>();
        };
        CHECK(contains(err({{"molecule", "NH3"}, {"method", "FCI"}}), "Full CI for NH3 is impossible"));
        CHECK(contains(err({{"molecule", "H2O2"}, {"method", "FCI"}}), "Full CI for H2O2 is impossible"));
        CHECK(contains(err({{"molecule", "NH3"}, {"method", "CCSDT"}, {"job", "freq"}}),
                       "CCSDT frequencies for NH3 are not possible"));
        CHECK(contains(err({{"molecule", "NH3"}, {"angle", 125.0}}), "between 60 and 120"));
        CHECK(contains(err({{"molecule", "H2O2"}, {"angle", 178.0}}), "between 30 and 175"));
        CHECK(contains(err({{"molecule", "H2O2"}, {"ROO", 0.1}}), "R(O–O)"));
        CHECK(contains(err({{"molecule", "H2O2"}, {"job", "scan"}, {"scan", "3 1 2 5 S 2 10"}}), "atom 5 does not exist"));
        const json e = run({{"molecule", "H2O2"}, {"method", "CCSD(T)"}, {"job", "freq"}, {"estimateOnly", true}});
        CHECK(e["error"].is_null());
        // 10.0: analytic Hessian (6 internal directions and 3 field components); numerical on request
        CHECK(contains(e["job"]["estimate"]["steps"][0]["how"].get<std::string>(),
                       "analytic Hessian: amplitude responses along the 6 internal directions"));
        CHECK(e["job"]["estimate"]["steps"][0]["evaluations"].get<int>() == 9);
        const json en6 = run({{"molecule", "H2O2"}, {"method", "CCSD(T)"}, {"job", "freq"}, {"estimateOnly", true},
                              {"hessian", "numerical"}});
        CHECK(contains(en6["job"]["estimate"]["steps"][0]["how"].get<std::string>(), "central differences of 12 analytic gradients"));
        const json en = run({{"molecule", "NH3"}, {"basis", "aug-cc-pVTZ"}, {"method", "HF"}});
        CHECK(en["error"].is_null());
        CHECK(en["molecule"]["angle"].get<double>() == 106.7);
        CHECK(contains(en["output"].get<std::string>(), "NH3, 10 electrons"));
    }

    std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t_all));
    return g_failures == 0 ? 0 : 1;
}
