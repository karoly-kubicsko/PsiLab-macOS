// Tests of the analytic correlated Hessians of ΨLab 10.0 (CTest). Arguments: "ccsd", "ccsdt" run
// only the CCSD or CCSD(T) tests; "ccsd-hof" (HOF without symmetry) and "ccsdt-extra" (HOF and
// NH3 with degenerate orbitals) are longer checks run by hand (10-30 minutes).
//  - the first-derivative electron-repulsion integrals against central differences of the
//    integrals themselves (transformed to MOs with fixed orbitals);
//  - the perturbed MO Hamiltonian H^(x) against central differences of the MO integrals with
//    the RHF orbitals of the displaced geometries (made to follow the reference orbitals);
//  - the general Hessian with the RHF densities against the analytic RHF Hessian (the orbital
//    part of the formalism, independent of any correlated method);
//  - the analytic Full CI Hessian against 4-point central differences of analytic Full CI
//    gradients (H2, HeH+ and H3+ at unsymmetric geometries, LiH);
//  - the analytic CCSD Hessian against differences of analytic CCSD gradients (LiH, H2O, HOF):
//    also the first test of the Z-vector terms (zero for RHF and Full CI).
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/cc.h"
#include "macqc/cc_gradient.h"
#include "macqc/ci.h"
#include "macqc/correlated_hessian.h"
#include "macqc/gradient.h"
#include "macqc/gto_derivatives.h"
#include "macqc/molecular.h"
#include "macqc/molecule.h"
#include "macqc/mo_integrals.h"
#include "macqc/rhf.h"
#include "macqc/rhf_derivatives.h"
#include "macqc/structure_file.h"

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

static Molecule molecule_from_xyz(const std::string& xyz, int charge = 0) {
    const StructureFile f = parse_structure_file(xyz, "xyz");
    return make_custom_molecule(f.atoms, charge, f.title);
}

static RHFOptions tight() {
    RHFOptions o;
    o.energy_tol = 1.0e-13;
    o.gradient_tol = 1.0e-11;
    o.max_iter = 200;
    return o;
}

// Molecules for the tests: water with the atoms in the built-in order, a non-symmetric
// geometry of HOF (no symmetry: every Hessian element independent), LiH.
static const char* kWater = "3\nH2O\nO 0 0 0\nH 0.7570 0 0.5859\nH -0.7570 0 0.5859\n";
static const char* kHOF = "3\nHOF\nO 0 0 0\nH 0.95 0.08 0.15\nF -0.32 0.05 1.36\n";
static const char* kLiH = "2\nLiH\nLi 0 0 0\nH 0.05 0.02 1.60\n";

static void test_rhf_through_general(const char* xyz, const std::string& basis_name, const std::string& label) {
    const Molecule m = molecule_from_xyz(xyz);
    const BasisSet b = build_basis(m, basis_name, "");
    const RHFResult r = run_rhf(m, b, tight());
    const auto t0 = Clock::now();
    const RHFHessianResult ref = rhf_hessian(m, b, r, 0, 1024.0, 1.0e-11);
    const double tref = since(t0);
    Matrix h_ao = r.H;
    const MOIntegrals mo = transform_to_mo(b, h_ao, r.C, 0);
    const auto t1 = Clock::now();
    const OrbitalResponse resp = rhf_orbital_response(m, b, r, 0, 1024.0, 1.0e-11);
    const CorrelatedDensity d = rhf_density(static_cast<int>(r.C.cols()), r.n_occupied);
    const CorrelatedHessianResult h = correlated_hessian(m, b, r, mo, d, resp, DensityResponse(), 0, 1024.0);
    const double tgen = since(t1);
    double worst = 0.0, big = 0.0;
    for (std::size_t x = 0; x < h.hessian.rows(); ++x)
        for (std::size_t y = 0; y < h.hessian.cols(); ++y) {
            worst = std::max(worst, std::fabs(h.hessian(x, y) - ref.hessian(x, y)));
            big = std::max(big, std::fabs(ref.hessian(x, y)));
        }
    std::printf("     %s %s: analytic RHF Hessian %.2f s, general framework %.2f s, largest element %.4f, asymmetry %.1e\n",
                label.c_str(), basis_name.c_str(), tref, tgen, big, h.asymmetry);
    check_max(label + " " + basis_name + ": general Hessian with the RHF densities vs rhf_hessian (E_h/bohr^2)", worst, 1.0e-8);
}

static ElectronicOptions tight_electronic() {
    ElectronicOptions eo;
    eo.scf_energy_tol = 1.0e-13;
    eo.scf_gradient_tol = 1.0e-11;
    eo.scf_max_iter = 200;
    eo.fci_residual_tol = 1.0e-10;
    eo.cc_amplitude_tol = 1.0e-11;
    eo.cc_energy_tol = 1.0e-13;
    eo.cc_max_iter = 300;
    return eo;
}

// Hessian and dipole derivatives by 4-point central differences of analytic gradients and their
// relaxed dipole moments (step h bohr).
struct FDReference {
    Matrix H;
    std::vector<std::array<double, 3>> apt;
    std::array<double, 3> dipole{};
};
static FDReference fd_hessian(const Molecule& m, const BasisSet& b, const std::string& method, double h) {
    const std::size_t n3 = 3 * m.atoms.size();
    FDReference r;
    r.H = Matrix(n3, n3);
    r.apt.assign(n3, {0.0, 0.0, 0.0});
    const ElectronicOptions eo = tight_electronic();
    r.dipole = run_gradient(m, b, method, eo).dipole;
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
        for (std::size_t y = 0; y < n3; ++y)
            r.H(x, y) = (8.0 * (g[0][y] - g[1][y]) - (g[2][y] - g[3][y])) / (12.0 * h);
        for (std::size_t c = 0; c < 3; ++c)
            r.apt[x][c] = (8.0 * (mu[0][c] - mu[1][c]) - (mu[2][c] - mu[3][c])) / (12.0 * h);
    }
    return r;
}

static void check_dipoles(const std::string& name, const CorrelatedHessianResult& h, const FDReference& ref) {
    CHECK(h.has_dipole_derivatives);
    if (!h.has_dipole_derivatives) return;
    double wd = 0.0, wa = 0.0;
    for (std::size_t c = 0; c < 3; ++c) wd = std::max(wd, std::fabs(h.dipole[c] - ref.dipole[c]));
    for (std::size_t x = 0; x < ref.apt.size(); ++x)
        for (std::size_t c = 0; c < 3; ++c) wa = std::max(wa, std::fabs(h.dipole_derivatives[x][c] - ref.apt[x][c]));
    check_max(name + ": relaxed dipole moment (e bohr) vs the analytic gradient", wd, 1.0e-8);
    check_max(name + ": dipole derivatives (e) vs differences of relaxed dipoles", wa, 2.0e-7);
}

static void test_fci(const char* xyz, int charge, const std::string& basis_name, const std::string& label, double tol) {
    const Molecule m = molecule_from_xyz(xyz, charge);
    const BasisSet b = build_basis(m, basis_name, "");
    const RHFResult r = run_rhf(m, b, tight());
    const MOIntegrals mo = transform_to_mo(b, r.H, r.C, 0);
    FCIOptions fo;
    fo.keep_vector = true;
    fo.residual_tol = 1.0e-11;
    fo.energy_tol = 1.0e-14;
    fo.max_iter = 300;
    const int na = m.electrons() / 2;
    const FCIResultGeneral f = run_fci_determinants(mo, na, na, r.nuclear_repulsion, r.total_energy, fo);
    const auto t0 = Clock::now();
    const CorrelatedHessianResult h = fci_hessian(m, b, r, mo, f.vector, f.electronic_energy, 0, 1024.0, 4096.0);
    const double ta = since(t0);
    const auto t1 = Clock::now();
    const FDReference fr = fd_hessian(m, b, "FCI", 0.005);
    const Matrix& ref = fr.H;
    const double tf = since(t1);
    double worst = 0.0, big = 0.0;
    for (std::size_t x = 0; x < ref.rows(); ++x)
        for (std::size_t y = 0; y < ref.cols(); ++y) {
            worst = std::max(worst, std::fabs(h.hessian(x, y) - ref(x, y)));
            big = std::max(big, std::fabs(ref(x, y)));
        }
    std::printf("     %s %s: %zu determinants; analytic %.2f s (CI response: %d iterations, residual %.1e; Z-vector %d), "
                "finite differences %.1f s; largest element %.4f, asymmetry %.1e\n",
                label.c_str(), basis_name.c_str(), f.determinants, ta, h.method_iterations, h.method_residual,
                h.zvector_iterations, tf, big, h.asymmetry);
    CHECK(h.method_converged && h.cphf_converged && h.zvector_converged);
    check_max(label + " " + basis_name + ": analytic Full CI Hessian vs differences of analytic gradients", worst, tol);
    check_dipoles(label + " " + basis_name + " Full CI", h, fr);
}

static void test_ccsd(const char* xyz, int charge, const std::string& basis_name, const std::string& label, double tol,
                      bool triples = false) {
    const std::string method = triples ? "CCSD(T)" : "CCSD";
    const Molecule m = molecule_from_xyz(xyz, charge);
    const BasisSet b = build_basis(m, basis_name, "");
    const RHFResult r = run_rhf(m, b, tight());
    const MOIntegrals mo = transform_to_mo(b, r.H, r.C, 0);
    const int o = r.n_occupied;
    const CCIntegrals ints = cc_integrals_from_mo(mo, o);
    const auto ladder = make_mo_ladder(mo, o);
    CCOptions co;
    co.amplitude_tol = 1.0e-11;
    co.energy_tol = 1.0e-13;
    co.max_iter = 300;
    const CCResult cc = run_ccsd_closed_shell(ints, *ladder, r.total_energy, triples, co);
    const CCDensityResult dl = ccsd_densities(ints, *ladder, cc.t1, cc.t2, triples, co);
    const auto t0 = Clock::now();
    const CorrelatedHessianResult h = ccsd_hessian(m, b, r, mo, cc.t1, cc.t2, dl.l1, dl.l2, triples, 0, 1024.0, 1.0e-10);
    const double ta = since(t0);
    const auto t1 = Clock::now();
    const FDReference fr = fd_hessian(m, b, method, 0.005);
    const Matrix& ref = fr.H;
    const double tf = since(t1);
    double worst = 0.0, big = 0.0;
    for (std::size_t x = 0; x < ref.rows(); ++x)
        for (std::size_t y = 0; y < ref.cols(); ++y) {
            worst = std::max(worst, std::fabs(h.hessian(x, y) - ref(x, y)));
            big = std::max(big, std::fabs(ref(x, y)));
        }
    std::printf("     %s %s %s: o=%d v=%d; analytic %.2f s (orbital %.2f s, amplitude/multiplier response %.2f s: %d iterations, "
                "residual %.1e; Z-vector %d; explicit %.2f s, response terms %.2f s), finite differences %.1f s; largest element %.4f, "
                "asymmetry %.1e\n",
                method.c_str(), label.c_str(), basis_name.c_str(), o, cc.virtuals, ta, h.orbital_seconds, h.method_seconds,
                h.method_iterations, h.method_residual, h.zvector_iterations, h.explicit_seconds, h.response_seconds, tf, big,
                h.asymmetry);
    CHECK(cc.converged && dl.lambda_converged);
    CHECK(h.method_converged && h.cphf_converged && h.zvector_converged);
    check_max(label + " " + basis_name + ": analytic " + method + " Hessian vs differences of analytic gradients", worst, tol);
    check_dipoles(label + " " + basis_name + " " + method, h, fr);
}

int main(int argc, char** argv) {
    const std::string only = argc > 1 ? argv[1] : "";
    if (only == "ccsd") {
        try {
            test_ccsd(kLiH, 0, "aug-cc-pVDZ", "LiH", 2.0e-7);
            test_ccsd(kWater, 0, "aug-cc-pVDZ", "H2O", 2.0e-7);
        } catch (const std::exception& e) {
            std::printf("FAIL exception: %s\n", e.what());
            ++g_failures;
        }
        std::printf("%d checks, %d failures\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    if (only == "ccsdt") {
        try {
            test_ccsd(kLiH, 0, "aug-cc-pVDZ", "LiH", 2.0e-7, true);
            test_ccsd(kWater, 0, "aug-cc-pVDZ", "H2O", 2.0e-7, true);
        } catch (const std::exception& e) {
            std::printf("FAIL exception: %s\n", e.what());
            ++g_failures;
        }
        std::printf("%d checks, %d failures\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    if (only == "ccsdt-extra") {
        try {
            test_ccsd(kHOF, 0, "aug-cc-pVDZ", "HOF (C1)", 2.0e-7, true);
            test_ccsd("4\nNH3\nN 0 0 0.1160\nH 0 0.9387 -0.2707\nH 0.8129 -0.4693 -0.2707\nH -0.8129 -0.4693 -0.2707\n", 0,
                      "aug-cc-pVDZ", "NH3 (C3v, degenerate orbitals)", 2.0e-7, true);
        } catch (const std::exception& e) {
            std::printf("FAIL exception: %s\n", e.what());
            ++g_failures;
        }
        std::printf("%d checks, %d failures\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    if (only == "ccsd-hof") {
        try {
            test_ccsd(kHOF, 0, "aug-cc-pVDZ", "HOF (C1)", 2.0e-7);
        } catch (const std::exception& e) {
            std::printf("FAIL exception: %s\n", e.what());
            ++g_failures;
        }
        std::printf("%d checks, %d failures\n", g_checks, g_failures);
        return g_failures == 0 ? 0 : 1;
    }
    try {
        test_rhf_through_general(kLiH, "aug-cc-pVDZ", "LiH");
        test_rhf_through_general(kWater, "aug-cc-pVDZ", "H2O");
        test_rhf_through_general(kHOF, "aug-cc-pVDZ", "HOF (C1)");
        test_fci("2\nH2\nH 0 0 0\nH 0.03 -0.02 0.76\n", 0, "aug-cc-pVDZ", "H2", 2.0e-7);
        test_fci("2\nHeH+\nHe 0 0 0\nH 0.02 0.01 0.80\n", 1, "aug-cc-pVDZ", "HeH+", 2.0e-7);
        test_fci("3\nH3+\nH 0 0 0\nH 0.90 0 0\nH 0.42 0.78 0.05\n", 1, "aug-cc-pVDZ", "H3+ (C1)", 2.0e-7);
        test_fci(kLiH, 0, "aug-cc-pVDZ", "LiH", 2.0e-7);
        test_ccsd(kLiH, 0, "aug-cc-pVDZ", "LiH", 2.0e-7);
        test_ccsd(kWater, 0, "aug-cc-pVDZ", "H2O", 2.0e-7);
        test_ccsd(kLiH, 0, "aug-cc-pVDZ", "LiH", 2.0e-7, true);
        test_ccsd(kWater, 0, "aug-cc-pVDZ", "H2O", 2.0e-7, true);
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        ++g_failures;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
