// Tests of the relaxed potential-energy-surface scan of ΨLab 7.0, stage 7.1 (CTest):
//  - torsions: value (IUPAC sign) and Wilson B row against finite differences (also at 0 and 180
//    degrees), differences across +-180 degrees, a back-transformed torsion step through 180
//    degrees; coordinates added by hand (a non-bonded distance, a bend kept above 175 degrees);
//  - the frozen coordinate: apply_constraints changes only the frozen coordinate (also to 180
//    degrees through the linear-bend components);
//  - scan points against independent constrained minima from PySCF + geomeTRIC
//    (backend/tests/data/scan_reference.json from reference/pyscf_scan_reference.py): H2O RHF scans
//    of a bond and of the angle (up to 180 degrees), a CCSD(T) bond scan;
//  - diatomic scans (single points equal to energy jobs), the ModRedundant line (with and without
//    the type letter), the summary table, the estimate, input errors and refusals, the streamed
//    output and cancellation during a scan.
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "macqc/internal_coordinates.h"
#include "macqc/jobs.h"
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
    std::printf("%s %-78s %.2e (tol %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), worst, tol);
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

static json run(const json& input, macqc_job* job = nullptr) {
    char* r = macqc_run_job(input.dump().c_str(), job);
    const json j = json::parse(r);
    macqc_free(r);
    return j;
}

static json run_logged(const std::string& label, const json& input) {
    const auto t0 = Clock::now();
    const json j = run(input);
    std::printf("     %-60s %6.1f s%s\n", label.c_str(), since(t0),
                j["error"].is_null() ? "" : (" error: " + j["error"].get<std::string>()).c_str());
    return j;
}

static std::string error_of(const json& input) {
    const json j = run(input);
    return j["error"].is_null() ? "" : j["error"].get<std::string>();
}

static Atom atom(int Z, double x, double y, double z) {
    Atom a;
    a.Z = Z;
    a.symbol = element_symbol(Z);
    a.position = {x, y, z};
    return a;
}

// A four-atom chain O1-O2 with H3 on O1 and H4 on O2 (the torsion H3-O1-O2-H4), built directly
// (H2O2 itself is a molecule of stage 7.2).
static Molecule chain(double dihedral_deg) {
    const double roo = 2.75, r = 1.83, a = 102.0 / kDeg, h = 0.5 * dihedral_deg / kDeg;
    Molecule m;
    m.name = "chain";
    m.atoms = {atom(8, -0.5 * roo, 0.0, 0.0), atom(8, 0.5 * roo, 0.0, 0.0),
               atom(1, -0.5 * roo + r * std::cos(a), r * std::sin(a) * std::sin(h), r * std::sin(a) * std::cos(h)),
               atom(1, 0.5 * roo - r * std::cos(a), -r * std::sin(a) * std::sin(h), r * std::sin(a) * std::cos(h))};
    return m;
}

static InternalCoordinate coord(InternalType t, std::vector<int> atoms) {
    InternalCoordinate c;
    c.type = t;
    c.atoms = std::move(atoms);
    return c;
}

// The first point of a scan result with the given value (angstrom or degrees).
static const json* point_at(const json& scan, double value) {
    for (const auto& p : scan["points"])
        if (std::fabs(p["value"].get<double>() - value) < 1e-8) return &p;
    return nullptr;
}

int main() {
    std::setvbuf(stdout, nullptr, _IOLBF, 0);
    const auto t_all = Clock::now();
    const std::string data = MACQC_TEST_DATA_DIR;
    const json ref = json::parse(read_file(data + "/scan_reference.json"));

    // ---- torsions
    std::printf("Torsions: value, Wilson B row, periodicity, back-transformation\n");
    {
        // IUPAC sign: the chain above is built with the given dihedral angle
        double worst = 0.0;
        for (double d : {-170.0, -113.7, -30.0, 0.0, 45.0, 113.7, 179.0, 180.0}) {
            const double v = internal_value(chain(d), coord(InternalType::Dihedral, {2, 0, 1, 3})) * kDeg;
            worst = std::max(worst, std::fabs(wrap_angle((v - d) / kDeg)) * kDeg);
        }
        check_max("dihedral H3-O1-O2-H4 of the built chain (IUPAC sign), degrees", worst, 1e-10);
        // B row against central differences, random and special geometries
        std::mt19937 rng(7);
        std::normal_distribution<double> nd(0.0, 0.25);
        double worst_b = 0.0;
        for (int trial = 0; trial < 12; ++trial) {
            Molecule m = chain(trial < 4 ? std::vector<double>{0.0, 180.0, 90.0, -179.9}[static_cast<std::size_t>(trial)]
                                         : 37.0 * trial);
            if (trial >= 4)
                for (auto& a : m.atoms)
                    for (double& x : a.position) x += nd(rng);
            const InternalCoordinate c = coord(InternalType::Dihedral, {2, 0, 1, 3});
            const std::vector<double> row = wilson_b_row(m, c);
            for (std::size_t x = 0; x < 12; ++x) {
                Molecule mp = m, mm = m;
                const double h = 1e-5;
                mp.atoms[x / 3].position[x % 3] += h;
                mm.atoms[x / 3].position[x % 3] -= h;
                const double fd = wrap_angle(internal_value(mp, c) - internal_value(mm, c)) / (2.0 * h);
                worst_b = std::max(worst_b, std::fabs(fd - row[x]));
            }
        }
        check_max("Wilson B row of a torsion vs central differences (also at 0, 90, 180 deg)", worst_b, 1e-8);
        // the set of a four-atom chain: 3 bonds, 2 bends, 1 torsion
        const Molecule m170 = chain(170.0);
        const InternalCoordinates ic(m170);
        CHECK(ic.size() == 6);
        int dih = -1;
        for (std::size_t i = 0; i < ic.size(); ++i)
            if (ic.coordinates()[i].type == InternalType::Dihedral) dih = static_cast<int>(i);
        CHECK(dih >= 0);
        CHECK(ic.coordinates()[static_cast<std::size_t>(dih)].definition(m170) == "D(H3-O1-O2-H4)");
        // differences across +-180 degrees
        std::vector<double> q1 = ic.values(m170), q0 = q1;
        q1[static_cast<std::size_t>(dih)] = -175.0 / kDeg;
        q0[static_cast<std::size_t>(dih)] = 175.0 / kDeg;
        check_max("difference of dihedrals -175 and 175 degrees = 10 degrees",
                  std::fabs(ic.difference(q1, q0)[static_cast<std::size_t>(dih)] * kDeg - 10.0), 1e-12);
        // a torsion step of 20 degrees from 170 degrees, through 180 (to -170)
        std::vector<double> dq(ic.size(), 0.0);
        dq[static_cast<std::size_t>(dih)] = 20.0 / kDeg;
        const BackTransform bt = back_transform(ic, m170, dq);
        const std::vector<double> d = ic.difference(ic.values(bt.molecule), ic.values(m170));
        double other = 0.0;
        for (std::size_t i = 0; i < d.size(); ++i)
            if (static_cast<int>(i) != dih) other = std::max(other, std::fabs(d[i]));
        CHECK(bt.converged);
        check_max("back-transformed torsion step 170 -> 190 = -170 degrees",
                  std::fabs(internal_value(bt.molecule, ic.coordinates()[static_cast<std::size_t>(dih)]) * kDeg + 170.0), 1e-8);
        check_max("  the other internal coordinates unchanged (bohr, rad)", other, 1e-9);
    }

    // ---- coordinates added by hand, the frozen coordinate
    std::printf("Coordinates added by hand; apply_constraints\n");
    {
        const Molecule w = make_molecule("H2O", 0.9572 / kBohrInAngstrom, 104.52);
        const InternalCoordinates ic(w, {coord(InternalType::Stretch, {1, 2})});
        CHECK(ic.size() == 4);
        CHECK(ic.find(coord(InternalType::Stretch, {2, 1})) == 3);
        const Molecule w178 = make_molecule("H2O", 0.9572 / kBohrInAngstrom, 178.0);
        const InternalCoordinates plain(w178), kept(w178, {coord(InternalType::Bend, {1, 0, 2})});
        CHECK(plain.size() == 4 && plain.coordinates()[2].type == InternalType::LinearBend);
        CHECK(kept.size() == 3 && kept.coordinates()[2].type == InternalType::Bend);
        CHECK(!kept.needs_rebuild(w178));

        CoordinateConstraint c;
        c.coordinate = coord(InternalType::Stretch, {0, 1});
        c.target = 1.2 / kBohrInAngstrom;
        double res = 1.0;
        const Molecule w2 = apply_constraints(w, {c}, &res);
        const InternalCoordinates i0(w);
        const std::vector<double> d = i0.difference(i0.values(w2), i0.values(w));
        check_max("apply_constraints: R(O1-H2) set to 1.2 angstrom", std::fabs(d[0] - (c.target - 0.9572 / kBohrInAngstrom)) + res,
                  1e-9);
        check_max("  R(O1-H3) and the angle unchanged", std::max(std::fabs(d[1]), std::fabs(d[2])), 1e-9);
        c.coordinate = coord(InternalType::Bend, {1, 0, 2});
        c.target = kPi;
        const Molecule lin = apply_constraints(w, {c}, &res);
        const InternalCoordinates il(lin, {coord(InternalType::LinearBend, {1, 0, 2})});
        double r1 = internal_value(lin, coord(InternalType::Stretch, {0, 1})),
               r2 = internal_value(lin, coord(InternalType::Stretch, {0, 2}));
        check_max("apply_constraints: angle set to 180 degrees (linear-bend components 0)",
                  std::fabs(internal_value(lin, coord(InternalType::Bend, {1, 0, 2})) * kDeg - 180.0), 1e-6);
        check_max("  bond lengths nearly unchanged (bohr; bend to 179 degrees, then the linear-bend components)",
                  std::max(std::fabs(r1 - 0.9572 / kBohrInAngstrom), std::fabs(r2 - 0.9572 / kBohrInAngstrom)), 1e-5);
        check_max("  both bonds equal (bohr)", std::fabs(r1 - r2), 1e-10);
    }

    // ---- scan points against PySCF + geomeTRIC
    std::printf("Scan points against constrained minima from PySCF + geomeTRIC\n");
    struct ScanRun {
        std::string method, line, coordinate;
        double start;
    };
    const std::vector<ScanRun> runs = {{"RHF", "1 2 S 1 0.35", "B 1 2", 0.85},
                                       {"RHF", "A 2 1 3 S 9 10.0", "A 2 1 3", 90.0},
                                       {"CCSD(T)", "B 1 2 S 1 0.1", "B 1 2", 1.1}};
    for (const auto& sr : runs) {
        const json* rc = nullptr;
        for (const auto& c : ref["cases"])
            if (c["molecule"] == "H2O" && c["method"] == sr.method && c["coordinate"] == sr.coordinate) rc = &c;
        CHECK(rc != nullptr);
        if (!rc) continue;
        const json r = run_logged("H2O " + sr.method + "/aug-cc-pVDZ scan " + sr.line,
                                  {{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", sr.method}, {"job", "scan"},
                                   {"scan", sr.line}, {"scanStart", sr.start}, {"optConvergence", "tight"}});
        CHECK(r["error"].is_null());
        if (!r["error"].is_null()) continue;
        const json& scan = r["scan"];
        CHECK(scan["converged"].get<bool>());
        CHECK(scan["relaxed"].get<bool>());
        for (const auto& rp : (*rc)["points"]) {
            const double v = rp["value"].get<double>();
            const json* p = point_at(scan, v);
            CHECK(p != nullptr);
            if (!p) continue;
            const json& g = (*p)["geometry"];
            const double r1 = g["R1"].get<double>() * kBohrInAngstrom, r2 = g["R2"].get<double>() * kBohrInAngstrom;
            const double ang = g["angle"].get<double>();
            const std::string tag = "  " + sr.method + " " + sr.coordinate + " = " + std::to_string(v).substr(0, 6);
            check_max(tag + ": energy (E_h)", std::fabs((*p)["energy"].get<double>() - rp["energy"].get<double>()), 2e-9);
            check_max(tag + ": R(O1-H2), R(O1-H3) (angstrom)",
                      std::max(std::fabs(r1 - rp["R1"].get<double>()), std::fabs(r2 - rp["R2"].get<double>())), 4e-5);
            check_max(tag + ": angle (degrees)", std::fabs(ang - rp["angle"].get<double>()), 5e-3);
            check_max(tag + ": frozen coordinate at its value", std::fabs((*p)["measured"].get<double>() - v), 1e-6);
        }
        if (sr.method == "RHF" && sr.coordinate == "A 2 1 3") {
            const std::string out = r["output"].get<std::string>();
            CHECK(contains(out, "Summary of the relaxed scan: HF/aug-cc-pVDZ"));
            CHECK(contains(out, "Scan point 10 of 10: A(H2-O1-H3) = 180.000000 deg"));
            CHECK(contains(out, "initial Hessian: Hessian of the previous scan point"));
            CHECK(contains(out, "frozen"));
            CHECK(scan["points"].size() == 10);
            CHECK(scan["lowest"].get<int>() == 2);   // 110 degrees, the nearest to the RHF minimum (105.5)
            CHECK(scan["highest"].get<int>() == 9);  // linear
            int steps = 0;
            for (const auto& p : scan["points"]) steps += p["steps"].get<int>();
            std::printf("     %d optimisation steps for 10 points\n", steps);
        }
    }

    // ---- diatomics: single points; the line without the type letter
    std::printf("Diatomic scans (single points)\n");
    {
        const json r = run_logged("LiH HF/aug-cc-pVDZ scan 1 2 S 2 0.1 from 1.5 angstrom",
                                  {{"molecule", "LiH"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "scan"},
                                   {"scan", "1 2 S 2 0.1"}, {"scanStart", 1.5}});
        CHECK(r["error"].is_null());
        const json& scan = r["scan"];
        CHECK(!scan["relaxed"].get<bool>());
        CHECK(scan["coordinate"]["type"] == "bond" && scan["coordinate"]["unit"] == "angstrom");
        CHECK(scan["points"].size() == 3);
        double worst = 0.0;
        for (const auto& p : scan["points"]) {
            CHECK(p["steps"].get<int>() == 0);
            const double R = p["value"].get<double>() / kBohrInAngstrom;
            const json e = run({{"molecule", "LiH"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"R", R}});
            worst = std::max(worst, std::fabs(e["rhf"]["totalEnergy"].get<double>() - p["energy"].get<double>()));
        }
        check_max("LiH scan points = energy jobs at the same R (E_h)", worst, 1e-10);
        CHECK(contains(r["output"].get<std::string>(), "a series of single-point energies"));
        const json h = run_logged("H2 FCI/aug-cc-pVDZ scan B 1 2 S 1 -0.05",
                                  {{"molecule", "H2"}, {"basis", "aug-cc-pVDZ"}, {"method", "FCI"}, {"job", "scan"},
                                   {"scan", "B 1 2 S 1 -0.05"}});
        CHECK(h["error"].is_null());
        CHECK(h["scan"]["points"].size() == 2);
        check_max("H2 scan starts at the input bond length (angstrom)",
                  std::fabs(h["scan"]["points"][0]["value"].get<double>() - 0.74144), 1e-12);
        check_max("  and steps by -0.05 angstrom",
                  std::fabs(h["scan"]["points"][1]["value"].get<double>() - 0.69144), 1e-12);
    }

    // ---- input errors, refusals, estimate
    std::printf("Input errors, refusals, estimate\n");
    {
        const json base = {{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "scan"}};
        auto with = [&](const json& extra) {
            json j = base;
            for (auto it = extra.begin(); it != extra.end(); ++it) j[it.key()] = it.value();
            return j;
        };
        CHECK(contains(error_of(base), "needs \"scan\""));
        CHECK(contains(error_of(with({{"scan", "3 1 2 4 S 2 10"}})), "atom 4 does not exist"));
        CHECK(contains(error_of(with({{"scan", "1 1 S 2 0.1"}})), "appears twice"));
        CHECK(contains(error_of(with({{"scan", "1 2 S 0 0.1"}})), "number of scan steps"));
        CHECK(contains(error_of(with({{"scan", "1 2 S 2 0"}})), "step size"));
        CHECK(contains(error_of(with({{"scan", "X 1 2 S 2 0.1"}})), "unknown scan coordinate type"));
        CHECK(contains(error_of(with({{"scan", "1 2 0.1"}})), "ModRedundant"));
        CHECK(contains(error_of(with({{"scan", "2 1 3 S 4 30"}, {"scanStart", 100.0}})), "angles must stay between"));
        CHECK(contains(error_of(with({{"scan", "1 2 S 4 -0.2"}})), "bonds must stay between"));
        CHECK(contains(error_of({{"molecule", "LiH"}, {"job", "scan"}, {"scan", "2 1 3 S 2 10"}}), "does not exist"));
        CHECK(contains(error_of(with({{"method", "CCSDTQ"}, {"scan", "1 2 S 2 0.05"}})), "not possible"));
        const json e = run(with({{"scan", "2 1 3 S 8 5.0"}, {"estimateOnly", true}}));
        CHECK(e["error"].is_null());
        const json& st = e["job"]["estimate"]["steps"][0];
        CHECK(st["name"] == "Relaxed scan" && st["evaluationKind"] == "steps");
        CHECK(st["evaluations"].get<int>() == 6 + 4 * 8);
        CHECK(contains(st["how"].get<std::string>(), "9 points"));
    }

    // ---- streamed output, cancellation
    std::printf("Streamed output and cancellation\n");
    {
        macqc_job* job = macqc_job_create();
        std::string result;
        std::atomic<bool> done{false};
        std::thread worker([&] {
            char* r = macqc_run_job(json({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "scan"},
                                          {"scan", "1 2 S 3 0.1"}})
                                        .dump()
                                        .c_str(),
                                    job);
            result = r;
            macqc_free(r);
            done = true;
        });
        std::string streamed;
        long long offset = 0;
        bool scan_progress = false;
        while (true) {
            const bool finished = done.load();
            long long total = 0;
            char* t = macqc_job_output(job, offset, &total);
            const std::string piece = t ? t : "";
            macqc_free(t);
            streamed += piece;
            offset += static_cast<long long>(piece.size());
            char* p = macqc_job_progress_text(job);
            if (p && std::string(p).rfind("Scan point", 0) == 0) scan_progress = true;
            macqc_free(p);
            if (finished) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
        }
        worker.join();
        const json r = json::parse(result);
        CHECK(r["error"].is_null());
        CHECK(streamed == r["output"].get<std::string>());
        CHECK(scan_progress);
        macqc_job_free(job);
    }
    {
        macqc_job* job = macqc_job_create();
        std::string result;
        std::thread worker([&] {
            char* r = macqc_run_job(json({{"molecule", "H2O"}, {"basis", "aug-cc-pVDZ"}, {"method", "HF"}, {"job", "scan"},
                                          {"scan", "2 1 3 S 20 3.0"}})
                                        .dump()
                                        .c_str(),
                                    job);
            result = r;
            macqc_free(r);
        });
        // cancel once the second point has started
        for (int i = 0; i < 20000; ++i) {
            char* p = macqc_job_progress_text(job);
            const std::string s = p ? p : "";
            macqc_free(p);
            if (s.rfind("Scan point 2 of 21", 0) == 0) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        const auto tc = Clock::now();
        macqc_job_cancel(job);
        worker.join();
        const json r = json::parse(result);
        CHECK(r["error"] == "calculation cancelled");
        char* t = macqc_job_output(job, 0, nullptr);
        CHECK(t && contains(t, "Calculation cancelled"));
        CHECK(t && !contains(t, "Scan point 4 of 21"));
        macqc_free(t);
        std::printf("     cancelled during scan point 2 after %.3f s\n", since(tc));
        macqc_job_free(job);
    }

    std::printf("\n%d checks, %d failures (%.1f s)\n", g_checks, g_failures, since(t_all));
    return g_failures == 0 ? 0 : 1;
}
