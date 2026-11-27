// ΨLab 20.0/21.0: tests of the nuclear-motion module (CTest).
//  - the CVRQD surfaces (backend/src/cvrqd.cpp, translated from Fortran) against values of the
//    gfortran build of the original routines (reference/nuclear_motion/pot_*_cvrqd.f);
//  - the minima used as the zero of the energies;
//  - J = 0 and J = 2 levels and expectation values against the original program
//    cpp_DOPI_expval_array.cpp (10 Hermite points in R1 and R2, 12 Legendre points; its energies
//    are shifted by the difference of its potential constant and the minimum);
//  - every eigensolver and linear solver gives the same levels; J = 1 in both parity blocks;
//  - the C API (JSON) and refusals;
//  - 21.0: Radau coordinates (the same levels as Jacobi), inverse power and Rayleigh-quotient
//    iteration, Chebyshev and BiCG, the linear solvers in Davidson and Jacobi-Davidson, the
//    eigenvalues only, the Gaussian start vector of Algorithm 5, the eigenvector pictures.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "macqc/cvrqd.h"
#include "macqc/macqc_api.h"
#include "macqc/nuclear_motion.h"
#include "nlohmann/json.hpp"

namespace macqc {
std::string archive_from_output_text(const std::string& text);  // calculation.h (22.0 tests)
}

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

    // ---- 6a. (21.0) J = 0 with the parity block of an earlier J > 0 (K = 1 ... J): the block K = 0,
    //          as for parity block 0 (20.0 refused it: "J = 0 has only the block K = 0")
    {
        NuclearMotionInput in = base;
        in.levels = 2;
        in.parity_block = 1;
        const NuclearMotionEstimate est = estimate_nuclear_motion(in);
        check(est.refusal.empty() && est.dimension == 1200, "J = 0 with parity block 1: no refusal, one K block");
        const NuclearMotionResult r1 = run_nuclear_motion(in);
        in.parity_block = 0;
        const NuclearMotionResult r0 = run_nuclear_motion(in);
        check(r1.input.parity_block == 0 && r1.k_values.size() == 1 && r1.k_values[0] == 0, "J = 0: the block K = 0");
        check_near("J = 0 with parity block 1 = parity block 0", r1.levels[0].energy, r0.levels[0].energy, 1e-14);
    }

    // ---- 6b. the sparsity picture: binned for a large matrix (2400 functions in 1024 pixels: every
    //          pixel row covers 2 or 3 rows), exact for a small one (one pixel per element)
    {
        NuclearMotionInput in;
        in.radial_points = 10;
        in.angular_points = 12;
        in.J = 1;
        in.levels = 1;
        in.eigensolver = "jacobi-davidson";
        const NuclearMotionResult r = run_nuclear_motion(in);
        check(r.matrix_whole.bins == 1024 && r.matrix_whole.span == 2400, "picture of 2400 functions: 1024 pixels");
        check(r.matrix_detail.bins == 1024 && r.matrix_detail.span == 1024, "detail: rows 1 ... 1024, one pixel per element");
        // the Coriolis block (K = 0 rows, K = 1 columns) is not empty, the corner far from both is
        bool coriolis = false;
        for (int j = 512; j < 1024 && !coriolis; ++j)
            if (r.matrix_whole.density[static_cast<std::size_t>(0 * 1024 + j)]) coriolis = true;
        check(coriolis, "picture: Coriolis coupling of K = 0 and K = 1");
        check(r.matrix_whole.log_max > r.matrix_whole.log_min, "picture: a range of |H_ij|");
    }
    // ---- 6b. the sparsity picture: exact for a small matrix (one pixel per element)
    {
        NuclearMotionInput in;
        in.radial_points = 4;
        in.angular_points = 5;
        in.levels = 2;
        in.eigensolver = "dense";
        const NuclearMotionResult r = run_nuclear_motion(in);
        const MatrixPicture& m = r.matrix_whole;
        check(m.bins == 80 && m.span == 80, "picture of an 80 x 80 matrix: one pixel per element");
        long long nz = 0;
        bool symmetric = true;
        for (int i = 0; i < m.bins; ++i)
            for (int j = 0; j < m.bins; ++j) {
                const unsigned char d = m.density[static_cast<std::size_t>(i * m.bins + j)];
                if (d) ++nz;
                if ((d != 0) != (m.density[static_cast<std::size_t>(j * m.bins + i)] != 0)) symmetric = false;
                if (i == j && d != 255) symmetric = false;
            }
        check(nz == r.nonzeros, "picture: as many pixels as non-zero elements");
        check(symmetric, "picture: symmetric pattern with the full diagonal");
        // ((i*4 + j)*5 + l): T(R2) couples j and j' (offset 5), T(R1) i and i' (offset 20)
        check(m.density[static_cast<std::size_t>(0 * 80 + 5)] == 255 && m.density[static_cast<std::size_t>(0 * 80 + 20)] == 255, "picture: kinetic-energy diagonals");
        check(m.density[static_cast<std::size_t>(0 * 80 + 6)] == 0, "picture: no element between them");
    }


    // ---- 8. (21.0) Radau coordinates: the same converged levels and expectation values as Jacobi;
    //         the reduced masses of the distances from the canonical point are the atomic masses
    {
        NuclearMotionInput in;
        in.radial_points = 16;
        in.angular_points = 20;
        in.levels = 5;
        const NuclearMotionResult j = run_nuclear_motion(in);
        in.coordinates = "radau";
        const NuclearMotionResult r = run_nuclear_motion(in);
        check(r.converged, "Radau converged");
        check_near("Radau mu1 = m(H)", r.mu1, 1.00727647 * 1.8228883e3, 1e-9);
        check_near("Radau mu2 = m(H)", r.mu2, 1.00727647 * 1.8228883e3, 1e-9);
        check(r.g1 > 0.0 && r.g1 < 0.05 && std::fabs(r.g1 - r.g2) < 1e-15, "Radau g1 = g2 small for H2O");
        for (int k = 0; k < 5; ++k)
            check_near("Radau = Jacobi level " + std::to_string(k + 1) + " (cm-1)", r.levels[static_cast<std::size_t>(k)].energy * kCm,
                       j.levels[static_cast<std::size_t>(k)].energy * kCm, 0.01);
        check_near("Radau <r1> = Jacobi", r.levels[0].r1.mean, j.levels[0].r1.mean, 1e-4);
        check_near("Radau <r2> = Jacobi", r.levels[0].r2.mean, j.levels[0].r2.mean, 1e-4);
        check_near("Radau <HOH> = Jacobi", r.levels[0].theta, j.levels[0].theta, 1e-4);
        check(r.report.find("Radau: R1 = r(K-H1)") != std::string::npos, "Radau in the output file");
        // D2O (other masses) runs as well
        in.isotopologue = CvrqdIsotopologue::D2_16O;
        in.levels = 2;
        const NuclearMotionResult d = run_nuclear_motion(in);
        in.coordinates = "jacobi";
        const NuclearMotionResult dj = run_nuclear_motion(in);
        check_near("D2O Radau = Jacobi zero-point level (cm-1)", d.levels[0].energy * kCm, dj.levels[0].energy * kCm, 0.01);
    }

    // ---- 9. (21.0) the new eigensolvers and linear solvers give the levels of the dense
    //         diagonalisation (8 levels, 1200 functions)
    {
        NuclearMotionInput d = base;
        d.eigensolver = "dense";
        const NuclearMotionResult ref = run_nuclear_motion(d);
        struct Case {
            const char *e, *l, *p;
            int inner;
        };
        const Case cases[] = {{"inverse-power", "pcg", "ilu", 0},       {"inverse-power", "minres", "jacobi", 0},
                              {"inverse-power", "chebyshev", "ilu", 0}, {"inverse-power", "bicg", "none", 0},
                              {"rayleigh-quotient", "minres", "ilu", 0}, {"rayleigh-quotient", "pcg", "ilu", 0},
                              {"rayleigh-quotient", "bicg", "jacobi", 0}, {"lanczos", "chebyshev", "jacobi", 0},
                              {"lanczos", "bicg", "ilu", 0},             {"davidson", "pcg", "ilu", 5},
                              {"davidson", "minres", "jacobi", 5},       {"davidson", "chebyshev", "none", 5},
                              {"davidson", "bicg", "ilu", 5},            {"davidson", "pcg", "none", 0},
                              {"jacobi-davidson", "minres", "ilu", 15},  {"jacobi-davidson", "chebyshev", "ilu", 15},
                              {"jacobi-davidson", "bicg", "jacobi", 15}, {"jacobi-davidson", "pcg", "none", 15}};
        for (const Case& c : cases) {
            NuclearMotionInput in = base;
            in.eigensolver = c.e;
            in.linear_solver = c.l;
            in.preconditioner = c.p;
            in.davidson_inner_iterations = c.inner;
            in.jd_inner_iterations = c.inner;
            for (int g = 0; g < 2; ++g) {
                in.gaussian_start = g == 1;
                const NuclearMotionResult q = run_nuclear_motion(in);
                const std::string name = std::string(c.e) + "/" + c.l + "/" + c.p + (g ? "/gauss" : "");
                check(q.converged, name + " converged");
                for (int k = 0; k < 8; ++k)
                    check_near(name + " level " + std::to_string(k + 1), q.levels[static_cast<std::size_t>(k)].energy,
                               ref.levels[static_cast<std::size_t>(k)].energy, 1e-9);
                for (const NuclearLevel& l : q.levels) check(l.residual < 1e-9, name + " residual");
            }
        }
        // one level: the inverse power iteration of the original (a single vector)
        NuclearMotionInput one = base;
        one.levels = 1;
        one.eigensolver = "inverse-power";
        const NuclearMotionResult ip = run_nuclear_motion(one);
        check_near("inverse power: the lowest level", ip.levels[0].energy, ref.levels[0].energy, 1e-12);
        check(ip.history.size() > 10 && ip.history.front().level == 1, "inverse power: history of level 1");
    }

    // ---- 10. (21.0) BiCG with r~0 = r0 reproduces conjugate gradients for the symmetric H; the
    //          linear solver and the preconditioner change the inner iterations, not the Lanczos levels
    {
        NuclearMotionInput in = base;
        const NuclearMotionResult cg = run_nuclear_motion(in);
        in.linear_solver = "bicg";
        const NuclearMotionResult bi = run_nuclear_motion(in);
        check(bi.linear_iterations == cg.linear_iterations, "BiCG: the iterations of CG");
        check(bi.matrix_vector_products > cg.matrix_vector_products * 3 / 2, "BiCG: about twice the products H x");
        in.linear_solver = "minres";
        in.preconditioner = "none";
        const NuclearMotionResult mr = run_nuclear_motion(in);
        check(mr.linear_iterations > 2 * cg.linear_iterations, "MINRES without a preconditioner: more inner iterations than ILU-CG");
        check(cg.first_solve.step == 1 && cg.first_solve.residuals.size() > 2, "linear-solve trace of step 1");
        check(cg.last_solve.step == cg.iterations && !cg.last_solve.residuals.empty(), "linear-solve trace of the last step");
        check(cg.first_solve.residuals.size() != mr.first_solve.residuals.size(), "the traces differ between the solvers");
        for (std::size_t k = 0; k < 8; ++k) check_near("Lanczos levels independent of the solver", mr.levels[k].energy, cg.levels[k].energy, 1e-10);
        in.linear_solver = "chebyshev";
        in.preconditioner = "jacobi";
        const NuclearMotionResult ch = run_nuclear_motion(in);
        check(ch.chebyshev_max > ch.chebyshev_min && ch.chebyshev_min > 0.0, "Chebyshev interval estimated");
    }

    // ---- 11. (21.0) the eigenvalues only (the original's sv = 0)
    {
        NuclearMotionInput in = base;
        in.eigenvectors = false;
        const NuclearMotionResult r = run_nuclear_motion(in);
        const NuclearMotionResult v = run_nuclear_motion(base);
        check(r.converged && !r.has_vectors, "eigenvalues only: converged, no vectors");
        for (std::size_t k = 0; k < 8; ++k) {
            check_near("eigenvalues only: level " + std::to_string(k + 1), r.levels[k].energy, v.levels[k].energy, 1e-9);
            check(std::isnan(r.levels[k].residual), "eigenvalues only: no residual");
            check(r.levels[k].plot.psi_r1r2.empty() && r.levels[k].k_weights.empty(), "eigenvalues only: no picture, no weights");
        }
        check(r.matrix_vector_products < v.matrix_vector_products, "eigenvalues only: no residual products");
        check(r.report.find("Expectation values") == std::string::npos, "eigenvalues only: no expectation values in the output");
        in.eigensolver = "jacobi-davidson";
        const NuclearMotionResult jd = run_nuclear_motion(in);
        check(!jd.has_vectors && jd.levels[0].residual < 1e-9, "eigenvalues only (Jacobi-Davidson): vectors not kept, residuals known");
    }

    // ---- 12. (21.0) the Gaussian start vector of Algorithm 5: close to the ground state, fewer iterations
    {
        NuclearMotionInput in = base;
        in.levels = 1;
        in.eigensolver = "jacobi-davidson";
        in.preconditioner = "ilu";
        const NuclearMotionResult plain = run_nuclear_motion(in);
        in.gaussian_start = true;
        const NuclearMotionResult g = run_nuclear_motion(in);
        check(g.start_overlap > 0.99, "Gaussian start vector: overlap with the ground state above 0.99");
        check(g.start_rayleigh > g.levels[0].energy && (g.start_rayleigh - g.levels[0].energy) * kCm < 100.0,
              "Gaussian start vector: Rayleigh quotient within 100 cm-1 of the zero-point level");
        check(g.iterations < plain.iterations, "Gaussian start vector: fewer Jacobi-Davidson iterations");
        check_near("Gaussian start: sigma(r) of the harmonic ground state (bohr)", g.start_sigma[0], 0.1278, 0.002);
        check_near("Gaussian start: averaged r (bohr; the original's 1.8437)", g.start_mu[0], 1.8437, 0.01);
        check(g.report.find("Algorithm 5") != std::string::npos, "Gaussian start vector in the output file");
        // J > 0: the K blocks combined by Rayleigh-Ritz
        in.J = 2;
        in.levels = 2;
        const NuclearMotionResult j2 = run_nuclear_motion(in);
        in.gaussian_start = false;
        in.eigensolver = "lanczos";
        const NuclearMotionResult j2l = run_nuclear_motion(in);
        check(j2.converged && j2.start_overlap > 0.9, "Gaussian start vector for J = 2");
        check_near("Gaussian start J = 2 level", j2.levels[0].energy, j2l.levels[0].energy, 1e-9);
    }

    // ---- 13. (21.0) the eigenvector pictures: sizes, the normalisation of the density, the
    //          nodeless ground state, the node of the bending fundamental along theta
    {
        NuclearMotionInput in = base;
        in.levels = 2;
        const NuclearMotionResult r = run_nuclear_motion(in);
        const WavefunctionPlot& p0 = r.levels[0].plot;
        check(p0.psi_r1r2.size() == 100 && p0.psi_r1t.size() == 120 && p0.rho_r2t.size() == 120, "picture sizes");
        check_near("psi = 1 at its largest", p0.psi_r1r2[static_cast<std::size_t>(p0.i * 10 + p0.j)], 1.0, 1e-6);
        double norm = 0.0, cross = 0.0;
        for (int i = 0; i < 10; ++i)
            for (int j = 0; j < 10; ++j) {
                norm += p0.rho_r1r2[static_cast<std::size_t>(i * 10 + j)] * r.r_weights[static_cast<std::size_t>(i)] * r.r_weights[static_cast<std::size_t>(j)];
                cross += p0.rho_r1r2[static_cast<std::size_t>(i * 10 + j)] - p0.rho_r1r2[static_cast<std::size_t>(j * 10 + i)];
            }
        check_near("density integrates to 1", norm, 1.0, 1e-5);
        check(std::fabs(cross) < 1e-3, "ground-state density symmetric in R1, R2 (Radau-like near symmetry of Jacobi)");
        float most_negative = 0.0f;
        for (float v : p0.psi_r1r2) most_negative = std::min(most_negative, v);
        for (float v : p0.psi_r1t) most_negative = std::min(most_negative, v);
        check(most_negative > -0.02f, "ground state nodeless (the DVR signs are right)");
        // the bending fundamental (level 2): a change of sign along theta through the maximum
        const WavefunctionPlot& p1 = r.levels[1].plot;
        float lo = 0.0f, hi = 0.0f;
        for (int l = 0; l < 12; ++l) {
            const float v = p1.psi_r1t[static_cast<std::size_t>(p1.i * 12 + l)];
            lo = std::min(lo, v);
            hi = std::max(hi, v);
        }
        check(hi > 0.99f && lo < -0.3f, "bending fundamental: a node along theta");
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
            // the sparsity pictures: 640 functions, one pixel per element (at most 1024 x 1024)
            check(j["matrixWhole"]["bins"].get<int>() == 640 && j["matrixWhole"]["span"].get<int>() == 640, "matrix picture: bins");
            check(j["matrixWhole"]["density"].get<std::string>().size() == (640 * 640 + 2) / 3 * 4, "matrix picture: base64 size");
            check(j["matrixDetail"]["span"].get<int>() == 640, "matrix detail: all 640 rows");
        }
        // 21.0 keys
        char* raw2 = macqc_nuclear_run_job(R"({"radialPoints":8,"angularPoints":10,"levels":3,"coordinates":"radau",
                                               "eigensolver":"rayleigh-quotient","linearSolver":"minres","gaussianStart":true})",
                                           nullptr);
        const json j2 = json::parse(raw2);
        macqc_free(raw2);
        check(j2["error"].is_null(), "21.0 JSON: no error");
        check(j2["input"]["coordinates"] == "radau" && j2["input"]["gaussianStart"] == true, "21.0 JSON echoes the input");
        check(j2["hasVectors"] == true && j2["startVector"]["overlap"].get<double>() > 0.5, "21.0 JSON: vectors, start vector");
        check(j2["levels"][0]["wavefunction"]["psiR1R2"].get<std::string>().size() == (8 * 8 * 4 + 2) / 3 * 4, "21.0 JSON: float32 picture");
        check(j2["firstSolve"]["residuals"].size() > 1 && j2["history"][0]["level"] == 1, "21.0 JSON: traces and levels of the steps");
        char* raw3 = macqc_nuclear_run_job(R"({"radialPoints":6,"angularPoints":8,"levels":2,"eigenvectors":false})", nullptr);
        const json j3 = json::parse(raw3);
        macqc_free(raw3);
        check(j3["hasVectors"] == false && j3["levels"][0]["residual"].is_null() && j3["levels"][0]["wavefunction"].is_null(),
              "21.0 JSON: eigenvalues only");
        // 21.0: the export of the eigenvalues, the eigenvectors and the Gaussian start vector
        {
            char* rawx = macqc_nuclear_run_job(R"({"radialPoints":6,"angularPoints":8,"levels":3,"gaussianStart":true,
                                                   "eigensolver":"jacobi-davidson","J":1})", nullptr);
            const json jx = json::parse(rawx);
            macqc_free(rawx);
            const std::string id = jx["exportId"].get<std::string>();
            const std::string path = "nuclear_export_test.txt";
            auto request = [&](const std::string& what, const std::string& eid) {
                char* e = macqc_nuclear_export(("{\"exportId\":\"" + eid + "\",\"path\":\"" + path + "\"," + what + "}").c_str());
                const json r = json::parse(e);
                macqc_free(e);
                return r;
            };
            check(request(R"("eigenvalues":true,"eigenvectors":true,"startVector":true)", id)["error"].is_null(), "export: all three");
            // parse: the level table and the vector rows
            std::FILE* f = std::fopen(path.c_str(), "r");
            check(f != nullptr, "export: file written");
            std::vector<double> energies, sum2(3, 0.0), overlap(1, 0.0);
            double q2 = 0.0;
            long long rows = 0;
            bool table_header = false;
            char line[4096];
            while (f && std::fgets(line, sizeof line, f)) {
                const std::string t = line;
                if (t.find("||Hx-Ex|| (E_h)") != std::string::npos) table_header = true;
                if (t.empty() || t[0] == '#' || t[0] == '\n' || t.find("Energy levels") == 0 || t.find("     n ") == 0) continue;
                std::vector<double> x;
                char* p = line;
                char* end = nullptr;
                for (double v = std::strtod(p, &end); end != p; v = std::strtod(p, &end)) {
                    x.push_back(v);
                    p = end;
                }
                if (x.size() == 5) energies.push_back(x[1]);
                if (x.size() == 14) {  // n K i j l R1 R2 theta s w q1 c1 c2 c3
                    ++rows;
                    q2 += x[10] * x[10];
                    for (int k = 0; k < 3; ++k) sum2[static_cast<std::size_t>(k)] += x[static_cast<std::size_t>(11 + k)] * x[static_cast<std::size_t>(11 + k)];
                    overlap[0] += x[10] * x[11];
                }
            }
            if (f) std::fclose(f);
            check(table_header && energies.size() == 3, "export: the level table");
            for (std::size_t k = 0; k < energies.size(); ++k)
                check_near("export: E (E_h) of level " + std::to_string(k + 1), energies[k], jx["levels"][k]["energy"].get<double>(), 1e-11);
            check(rows == jx["dimension"].get<long long>(), "export: one row per basis function (J = 1: two K blocks)");
            for (int k = 0; k < 3; ++k) check_near("export: eigenvector normalised", sum2[static_cast<std::size_t>(k)], 1.0, 1e-9);
            check_near("export: start vector normalised", q2, 1.0, 1e-9);
            check_near("export: |<q1|x1>| as reported", std::fabs(overlap[0]), jx["startVector"]["overlap"].get<double>(), 1e-9);
            check(request(R"("eigenvalues":true)", "nuclear-0-0")["error"].is_string(), "export: unknown result refused");
            check(request(R"("eigenvalues":false)", id)["error"].is_string(), "export: nothing chosen refused");
            std::remove(path.c_str());
            // eigenvalues only: no eigenvectors; no Gaussian start vector without it
            char* rawy = macqc_nuclear_run_job(R"({"radialPoints":5,"angularPoints":6,"levels":2,"eigenvectors":false})", nullptr);
            const json jy = json::parse(rawy);
            macqc_free(rawy);
            const std::string idy = jy["exportId"].get<std::string>();
            check(request(R"("eigenvectors":true)", idy)["error"].is_string(), "export: no eigenvectors with the eigenvalues only");
            check(request(R"("startVector":true)", idy)["error"].is_string(), "export: no start vector without it");
            check(request(R"("eigenvalues":true)", idy)["error"].is_null(), "export: eigenvalues only");
            check(request(R"("eigenvalues":true)", id)["error"].is_string(), "export: a later run replaces the earlier result");
            std::remove(path.c_str());
        }
        char* bad3 = macqc_nuclear_estimate(R"({"coordinates":"valence"})");
        const json b3 = json::parse(bad3);
        macqc_free(bad3);
        check(b3["refusal"].is_string(), "refusal of unknown coordinates");
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

    // ---- 22.0: the input files of psilab, the layout of the output file, its archive, the
    // pictures of every level
    {
        using nlohmann::json;
        const std::string text = "# the input of the app's screenshot\nmolecule = H2(16)O\ncoordinates = Radau\nns1 = 50\nnp  = 50\n"
                                 "Re (bohr) = 2\n";
        bool threw = false;
        try {
            parse_nuclear_input_text(text);  // "Re (bohr)" is not a key
        } catch (const std::runtime_error& e) {
            threw = std::string(e.what()).rfind("line 6:", 0) == 0;
        }
        check(threw, "22.0 input: an unknown key names its line");
        const NuclearMotionInput in = parse_nuclear_input_text(
            "calculation = nuclear\nmolecule = H2(16)O\ncoordinates = Radau\nns1 = 50\nnp = 50  # angular\nre = 2\nrmax = 0,8\n"
            "jmax = 0\nlevels = 100\neigenvectors = yes\nexpectation_values = yes\neigensolver = Lánczos\nlinear_solver = CG\n"
            "preconditioner = ILU\ngaussian_start = ON\n");
        check(in.coordinates == "radau" && in.radial_points == 50 && in.angular_points == 50 && in.levels == 100 &&
                  std::fabs(in.r_half_width - 0.8) < 1e-15 && in.eigensolver == "lanczos" && in.linear_solver == "pcg" &&
                  in.preconditioner == "ilu" && in.gaussian_start && in.eigenvectors,
              "22.0 input: the settings of the app's screenshot");
        const NuclearMotionInput back = parse_nuclear_input_text(nuclear_input_text(in));
        check(back.coordinates == in.coordinates && back.levels == in.levels && back.r_centre == in.r_centre &&
                  back.r_half_width == in.r_half_width && back.linear_solver == in.linear_solver &&
                  back.eigenvalue_tolerance == in.eigenvalue_tolerance && back.gaussian_start == in.gaussian_start,
              "22.0 input: the Input section of the output file reads back");
        check(std::fabs(parse_nuclear_input_text("re = 1.0 angstrom\n").r_centre - 1.0 / 0.529177210903) < 1e-12,
              "22.0 input: angstrom");
        check(is_nuclear_input_text("calculation = nuclear\n") && is_nuclear_input_text("molecule = D2(16)O\n") &&
                  is_nuclear_input_text("molecule = H2(18)O # heavy\n") && !is_nuclear_input_text("molecule = H2O\nbasis = STO-3G\n"),
              "22.0 input: nuclear-motion inputs are recognised");
        threw = false;
        try {
            parse_nuclear_input_text("levels = 3\nlevels = 4\n");
        } catch (const std::runtime_error& e) {
            threw = std::string(e.what()).find("given twice") != std::string::npos;
        }
        check(threw, "22.0 input: a key given twice");
        check(parse_nuclear_input_text(nuclear_example_input()).coordinates == "radau", "22.0 input: the example reads");

        // 70 levels: pictures of all of them (21.0: the first 60), the level numbers right-aligned
        char* raw = macqc_nuclear_run_job(R"({"radialPoints":6,"angularPoints":8,"levels":105,"eigensolver":"dense"})", nullptr);
        const json j = json::parse(raw);
        macqc_free(raw);
        check(j["error"].is_null() && j["levels"].size() == 105, "22.0: 105 levels");
        check(!j["levels"][104]["wavefunction"].is_null(), "22.0: a picture of level 105 (21.0: up to 60)");
        const std::string out = j["output"].get<std::string>();
        auto has = [&](const std::string& t) { return out.find(t) != std::string::npos; };
        check(has("\n     n  <r1>     <r1^2>^1/2") && has("\n  9 r1   ") && has("\n 10 r1   ") && has("\n100 r1   ") &&
                  has("\n    r2   ") && !has("  10 r1"),
              "22.0: the level numbers right-aligned (1, 10, 100 in one column)");
        check(has("\n" + std::string(39, ' ') + "\xCE\xA8Lab 2"), "22.0: the program name centred");
        check(has("\n---- Input ---") && has("\n  calculation               = nuclear\n"), "22.0: the Input section");
        check(has("\n========================================================================================\n  Normal termination of "
                  "\xCE\xA8Lab ") &&
                  out.size() > 89 && out.compare(out.size() - 89, 89, std::string(88, '=') + "\n") == 0,
              "22.0: the termination banner ends the output file");
        // the archive: the result as the app receives it, the matrix pictures with shortened zero runs
        const std::string a = archive_from_output_text(out);
        check(!a.empty(), "22.0: the output file has an archive");
        const json ja = json::parse(a.empty() ? std::string("{}") : a);
        check(ja["archive"]["kind"] == "nuclear-motion" && ja["levels"].size() == 105 && !ja.contains("output") &&
                  !ja.contains("exportId"),
              "22.0: the archive is the nuclear-motion result");
        check(ja["levels"][104]["wavefunction"] == j["levels"][104]["wavefunction"] && ja["history"] == j["history"] &&
                  ja["input"] == j["input"] && ja["rGrid"] == j["rGrid"],
              "22.0: the archive has the pictures, the history, the input and the grids");
        auto expand = [](const std::string& s) {
            if (s.rfind("rle:", 0) != 0) return s;
            std::string o;
            for (std::size_t i = 4; i < s.size();) {
                if (s[i] == '~') {
                    const std::size_t dot = s.find('.', i);
                    o.append(static_cast<std::size_t>(std::stoul(s.substr(i + 1, dot - i - 1))), 'A');
                    i = dot + 1;
                } else {
                    o += s[i++];
                }
            }
            return o;
        };
        const std::string rle = ja["matrixWhole"]["density"].get<std::string>();
        check(rle.rfind("rle:", 0) == 0 && rle.size() < j["matrixWhole"]["density"].get<std::string>().size() / 2 &&
                  expand(rle) == j["matrixWhole"]["density"].get<std::string>() &&
                  expand(ja["matrixDetail"]["magnitude"].get<std::string>()) == j["matrixDetail"]["magnitude"].get<std::string>(),
              "22.0: the matrix pictures of the archive expand to those of the result");
        check(nuclear_termination_text("error", 1.0, "boom").find("  Error termination of \xCE\xA8Lab ") != std::string::npos,
              "22.0: the error banner");
    }

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
