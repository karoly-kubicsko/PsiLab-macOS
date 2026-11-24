// ΨLab 20.0: tests of the nuclear-motion module (CTest).
//  - the CVRQD surfaces (backend/src/cvrqd.cpp, translated from Fortran) against values of the
//    gfortran build of the original routines (reference/nuclear_motion/pot_*_cvrqd.f);
//  - the minima used as the zero of the energies;
//  - J = 0 and J = 2 levels and expectation values against the original program
//    cpp_DOPI_expval_array.cpp (10 Hermite points in R1 and R2, 12 Legendre points; its energies
//    are shifted by the difference of its potential constant and the minimum);
//  - every eigensolver and linear solver gives the same levels; J = 1 in both parity blocks;
//  - the C API (JSON) and refusals.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "macqc/cvrqd.h"
#include "macqc/macqc_api.h"
#include "macqc/nuclear_motion.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check(bool ok, const std::string& what) {
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

static void check_near(const std::string& what, double got, double want, double tol) {
    const bool ok = std::fabs(got - want) <= tol;
    ++g_checks;
    if (!ok) {
        ++g_failures;
        std::printf("FAIL: %s: got %.15g, want %.15g (difference %.3g, tolerance %.3g)\n", what.c_str(), got, want,
                    got - want, tol);
    }
}

constexpr double kCm = 219474.6313632;     // E_h -> cm-1 (CODATA 2018), as the module
constexpr double kCmOriginal = 219474.630;  // the original program's factor
constexpr double kPi = 3.14159265358979323846;

int main() {
    // ---- 1. CVRQD surfaces: values of the gfortran build (r1, r2 bohr; cos HOH)
    struct Point {
        double r1, r2, c, v;
    };
    const Point h216o[] = {{1.810, 1.810, -0.250, 0.0035803250219873335},
                           {1.600, 2.300, -0.500, 0.065328905721448463},
                           {3.900, 1.700, 0.200, 0.1624535878602768},
                           {2.200, 2.600, -0.900, 0.1562271665165739}};
    const Point h217o[] = {{1.810, 1.810, -0.250, 0.0034399602641121681},
                           {1.600, 2.300, -0.500, 0.065188288187538551},
                           {3.900, 1.700, 0.200, 0.16231312253466287},
                           {2.200, 2.600, -0.900, 0.15608645470691482}};
    const Point h218o[] = {{1.810, 1.810, -0.250, 0.003315749486441442},
                           {1.600, 2.300, -0.500, 0.065063853724746795},
                           {3.900, 1.700, 0.200, 0.16218882276316468},
                           {2.200, 2.600, -0.900, 0.15596193681816703}};
    const Point d216o[] = {{1.810, 1.810, -0.250, 0.0033754464469109126},
                           {1.600, 2.300, -0.500, 0.065124304453378565},
                           {3.900, 1.700, 0.200, 0.16219298396120171},
                           {2.200, 2.600, -0.900, 0.1559901909625665}};
    struct Set {
        CvrqdIsotopologue iso;
        const Point* p;
        const char* name;
    };
    const Set sets[] = {{CvrqdIsotopologue::H2_16O, h216o, "H2(16)O"},
                        {CvrqdIsotopologue::H2_17O, h217o, "H2(17)O"},
                        {CvrqdIsotopologue::H2_18O, h218o, "H2(18)O"},
                        {CvrqdIsotopologue::D2_16O, d216o, "D2(16)O"}};
    for (const Set& s : sets)
        for (int k = 0; k < 4; ++k) {
            const Point& q = s.p[k];
            // the translation keeps the arithmetic of gfortran: identical up to the last bit or two
            check_near(std::string("CVRQD ") + s.name + " point " + std::to_string(k + 1),
                       cvrqd_potential(s.iso, q.r1, q.r2, q.c), q.v, 4e-17 * std::max(1.0, std::fabs(q.v) * 10));
        }
    // symmetric in r1 and r2
    check_near("CVRQD symmetric", cvrqd_potential(CvrqdIsotopologue::H2_16O, 1.7, 2.1, -0.3),
               cvrqd_potential(CvrqdIsotopologue::H2_16O, 2.1, 1.7, -0.3), 1e-15);

    // ---- 2. the minima (the zero of the energies): the reported minimum is the lowest value
    //         of the surface around it, and equals the value at the reported geometry
    for (const Set& s : sets) {
        NuclearMotionInput in;
        in.isotopologue = s.iso;
        in.radial_points = 3;
        in.angular_points = 3;
        in.levels = 1;
        in.eigensolver = "dense";
        const NuclearMotionResult r = run_nuclear_motion(in);
        const double vmin = r.potential_minimum, rm = r.minimum_r, am = r.minimum_angle;
        check_near(std::string("minimum value ") + s.name, cvrqd_potential(s.iso, rm, rm, std::cos(am)), vmin, 2e-15);
        bool lowest = true;
        for (double dr1 : {-1e-3, 0.0, 1e-3})
            for (double dr2 : {-1e-3, 0.0, 1e-3})
                for (double da : {-1e-3, 0.0, 1e-3})
                    if (cvrqd_potential(s.iso, rm + dr1, rm + dr2, std::cos(am + da)) < vmin - 1e-15) lowest = false;
        check(lowest, std::string("minimum is lowest nearby ") + s.name);
    }

    // ---- 3. J = 0, H2(16)O against the original program (ns1 = 10, np = 12, 8 levels)
    NuclearMotionInput base;
    base.radial_points = 10;
    base.angular_points = 12;
    base.levels = 8;
    base.residual_tolerance = 1e-9;
    const double c_original_h216o = 0.00358028360318796825754850132739;  // the original program's shift
    {
        const NuclearMotionResult r = run_nuclear_motion(base);
        check(r.converged, "J=0 Lanczos converged");
        check(r.dimension == 1200, "J=0 dimension 10*10*12");
        const double ref[8] = {4639.62426, 6211.51736, 7726.42692, 8297.79594,
                               8397.37576, 9441.46382, 9845.55659, 9931.52455};
        for (int k = 0; k < 8; ++k)
            check_near("J=0 level " + std::to_string(k + 1) + " vs original (cm-1)",
                       (r.levels[static_cast<std::size_t>(k)].energy + r.potential_minimum - c_original_h216o) * kCmOriginal, ref[k],
                       2e-5);
        // expectation values (original: expectval_rot; angstrom with 0.529177249, degrees)
        const double a = 0.529177249;
        const NuclearLevel& l1 = r.levels[0];
        check_near("<r1> level 1", l1.r1.mean * a, 0.975602, 1e-6);
        check_near("<r1^2>^1/2 level 1", l1.r1.rms * a, 0.978032, 1e-6);
        check_near("<r1^-3>^-1/3 level 1", l1.r1.inverse3 * a, 0.965883, 1e-6);
        check_near("<r2> level 1", l1.r2.mean * a, 0.975534, 1e-6);
        check_near("<HOH> level 1", l1.theta * 180.0 / kPi, 104.5802, 1e-4);
        const NuclearLevel& l8 = r.levels[7];
        check_near("<r1> level 8", l8.r1.mean * a, 0.997129, 1e-6);
        check_near("<r1^3>^1/3 level 8", l8.r1.cubic * a, 1.00727, 1e-5);
        check_near("<r2^-1>^-1 level 8", l8.r2.inverse * a, 0.986514, 1e-6);
        check_near("<HOH> level 8", l8.theta * 180.0 / kPi, 102.8045, 1e-4);
        for (const NuclearLevel& l : r.levels) check(l.residual < 1e-9, "J=0 residual below tolerance");

        // every solver: the same levels
        const char* solvers[][3] = {{"dense", "pcg", "ilu"},           {"davidson", "pcg", "ilu"},
                                    {"jacobi-davidson", "pcg", "ilu"}, {"lanczos", "minres", "ilu"},
                                    {"lanczos", "pcg", "jacobi"},      {"lanczos", "pcg", "none"},
                                    {"lanczos", "minres", "none"}};
        for (const auto& s : solvers) {
            NuclearMotionInput in = base;
            in.eigensolver = s[0];
            in.linear_solver = s[1];
            in.preconditioner = s[2];
            const NuclearMotionResult q = run_nuclear_motion(in);
            const std::string name = std::string(s[0]) + "/" + s[1] + "/" + s[2];
            check(q.converged, name + " converged");
            for (int k = 0; k < 8; ++k)
                check_near(name + " level " + std::to_string(k + 1), q.levels[static_cast<std::size_t>(k)].energy,
                           r.levels[static_cast<std::size_t>(k)].energy, 1e-9);
        }
    }

    // ---- 4. other isotopologues against the original program (lowest levels)
    {
        NuclearMotionInput in = base;
        in.levels = 3;
        in.isotopologue = CvrqdIsotopologue::D2_16O;
        NuclearMotionResult r = run_nuclear_motion(in);
        const double cd = 0.00337543999746865977740843511867;
        const double refd[3] = {3369.64229, 4678.99994, 5362.74861};
        for (int k = 0; k < 3; ++k)
            check_near("D2(16)O level " + std::to_string(k + 1), (r.levels[static_cast<std::size_t>(k)].energy + r.potential_minimum - cd) * kCmOriginal,
                       refd[k], 2e-5);
        in.isotopologue = CvrqdIsotopologue::H2_17O;
        r = run_nuclear_motion(in);
        const double c17 = 0.00343991010054164597375307543814;
        const double ref17[3] = {4631.96363, 6196.71260, 7727.00398};
        for (int k = 0; k < 3; ++k)
            check_near("H2(17)O level " + std::to_string(k + 1), (r.levels[static_cast<std::size_t>(k)].energy + r.potential_minimum - c17) * kCmOriginal,
                       ref17[k], 2e-5);
        in.isotopologue = CvrqdIsotopologue::H2_18O;
        r = run_nuclear_motion(in);
        // the heavier oxygen lowers the zero-point energy (no reference: the original program did not
        // evaluate the H2(18)O surface)
        check(r.levels[0].energy * kCm > 4600.0 && r.levels[0].energy * kCm < 4639.6, "H2(18)O zero-point energy below H2(16)O");
    }

    // ---- 5. J = 2 in both parity blocks against the original program; K weights
    {
        NuclearMotionInput in = base;
        in.levels = 3;
        in.J = 2;
        in.parity_block = 0;
        NuclearMotionResult r = run_nuclear_motion(in);
        check(r.dimension == 3600, "J=2 K=0..2 dimension");
        check_near("J=2 K=0..2 lowest vs original", (r.levels[0].energy + r.potential_minimum - c_original_h216o) * kCmOriginal,
                   4706.99801, 2e-5);
        double wsum = 0.0;
        for (double w : r.levels[0].k_weights) wsum += w;
        check_near("K weights sum to 1", wsum, 1.0, 1e-12);
        NuclearMotionInput jd = in;
        jd.eigensolver = "jacobi-davidson";
        const NuclearMotionResult rj = run_nuclear_motion(jd);
        for (int k = 0; k < 3; ++k)
            check_near("J=2 Jacobi-Davidson level " + std::to_string(k + 1), rj.levels[static_cast<std::size_t>(k)].energy,
                       r.levels[static_cast<std::size_t>(k)].energy, 1e-9);
        in.parity_block = 1;
        r = run_nuclear_motion(in);
        check(r.dimension == 2400, "J=2 K=1..2 dimension");
        check_near("J=2 K=1..2 lowest vs original", (r.levels[0].energy + r.potential_minimum - c_original_h216o) * kCmOriginal,
                   4732.63910, 2e-5);
    }

    // ---- 6. J = 1: the rotational levels 1_01, 1_11 (K = 0..1) and 1_10 (K = 1) above the
    //         zero-point level, in order (this small basis: 22.4, 36.8 and 40.6 cm-1)
    {
        NuclearMotionInput in = base;
        in.levels = 1;
        const double zpe = run_nuclear_motion(in).levels[0].energy;
        in.J = 1;
        in.levels = 2;
        in.parity_block = 0;
        const NuclearMotionResult a = run_nuclear_motion(in);
        in.parity_block = 1;
        in.levels = 1;
        const NuclearMotionResult b = run_nuclear_motion(in);
        const double e101 = (a.levels[0].energy - zpe) * kCm, e111 = (a.levels[1].energy - zpe) * kCm,
                     e110 = (b.levels[0].energy - zpe) * kCm;
        check(e101 > 20 && e101 < 25, "J=1 1_01 above the zero-point level");
        check(e111 > e101 && e111 < e110, "J=1 1_11 between 1_01 and 1_10");
        check(e110 > 38 && e110 < 44, "J=1 1_10 above the zero-point level");
    }

    // ---- 7. C API: JSON round trip, estimate, refusals
    {
        const char* input = R"({"isotopologue":"H2(16)O","radialPoints":8,"angularPoints":10,"levels":4,
                                "eigensolver":"jacobi-davidson","threads":2})";
        char* raw = macqc_nuclear_run_job(input, nullptr);
        check(raw != nullptr, "macqc_nuclear_run_job returns a result");
        if (raw) {
            const json j = json::parse(raw);
            macqc_free(raw);
            check(j["error"].is_null(), "JSON result has no error");
            check(j["levels"].size() == 4, "JSON result has 4 levels");
            check(j["dimension"].get<long long>() == 640, "JSON dimension");
            check(j["input"]["eigensolver"] == "jacobi-davidson", "JSON echoes the input");
            check(j["output"].get<std::string>().find("Normal termination") != std::string::npos, "JSON output file");
            check(j["levels"][0]["r1"]["mean"].get<double>() > 1.7, "JSON expectation values");
        }
        char* est = macqc_nuclear_estimate(R"({"radialPoints":30,"angularPoints":40,"J":3})");
        const json e = json::parse(est);
        macqc_free(est);
        check(e["dimension"].get<long long>() == 4LL * 30 * 30 * 40, "estimate: dimension of J=3, K=0..3");
        check(e["refusal"].is_null(), "estimate: no refusal");
        char* bad = macqc_nuclear_run_job(R"({"radialPoints":1})", nullptr);
        const json b = json::parse(bad);
        macqc_free(bad);
        check(b["error"].is_string(), "refusal of one radial point");
        char* bad2 = macqc_nuclear_estimate(R"({"eigensolver":"dense","radialPoints":20,"angularPoints":20})");
        const json b2 = json::parse(bad2);
        macqc_free(bad2);
        check(b2["refusal"].is_string(), "refusal of a dense diagonalisation of 8000 functions");
        bool threw = false;
        try {
            isotopologue_from_name("HDO");
        } catch (const std::exception&) {
            threw = true;
        }
        check(threw, "unknown isotopologue refused");
        check(isotopologue_from_name("D2O") == CvrqdIsotopologue::D2_16O, "D2O is D2(16)O");
    }

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
