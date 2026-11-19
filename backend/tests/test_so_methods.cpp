// Tests of the correlated methods of ΨLab 15.0 (CTest):
//  - the energies of MP2-MP5, CISD, CISD(T), CISDT, CISDTQ, QCISD, QCISD(T), QCISDT, QCISDTQ and
//    CCSD ... CCSDTQ of the spin-orbital engine, on RHF and UHF references, against brute-force
//    references in the determinant space (backend/tests/data/so_reference.json, written by
//    reference/so_brute_force.py from the definitions of the methods, independently of the
//    engine's generated equations);
//  - closed-shell analytic gradients against differences of energies and an analytic Hessian
//    against differences of analytic gradients;
//  - UHF-based analytic gradients against differences of energies;
//  - (open_shell_correlated.h) the orbital part of the open-shell derivatives against the UHF and
//    ROHF SCF gradients and Hessians, ROHF-based analytic gradients against differences of
//    energies, and UHF- and ROHF-based analytic Hessians and dipole derivatives against
//    differences of analytic gradients and relaxed dipole moments (OH, H2O+);
//  - open-shell Full CI (ROHF and UHF orbitals): the energy of the spin-resolved densities, the
//    gradient against differences of energies, the Hessian and dipole derivatives against
//    differences of analytic gradients and dipole moments (OH, H2O+).
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/ci.h"
#include "macqc/gto_derivatives.h"
#include "macqc/molecule.h"
#include "macqc/mo_integrals.h"
#include "macqc/open_shell_correlated.h"
#include "macqc/rhf.h"
#include "macqc/scf_derivatives.h"
#include "macqc/so_derivatives.h"
#include "macqc/so_methods.h"
#include "macqc/uhf_correlated.h"
#include "nlohmann/json.hpp"

using nlohmann::json;
using namespace macqc;

static int g_checks = 0, g_failures = 0;

static void check_near(const std::string& name, double got, double want, double tol) {
    ++g_checks;
    const bool ok = std::fabs(got - want) <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-44s got %.12f expected %.12f (diff %.1e)\n", ok ? "ok  " : "FAIL", name.c_str(), got, want,
                got - want);
}

static void check_small(const std::string& name, double value, double tol) {
    ++g_checks;
    const bool ok = std::isfinite(value) && std::fabs(value) <= tol;
    if (!ok) ++g_failures;
    std::printf("%s %-44s %.2e (tolerance %.0e)\n", ok ? "ok  " : "FAIL", name.c_str(), value, tol);
}

static std::string read_file(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("cannot open " + path);
    return std::string(std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>());
}

static std::string basis_name(const std::string& file) {
    if (file == "sto-3g.gbs") return "STO-3G";
    if (file == "6-31g.gbs") return "6-31G";
    throw std::runtime_error("unexpected basis file " + file);
}

// One SCF reference with its MO integrals and the engine's Hamiltonian.
struct Reference {
    RHFResult scf;
    MOIntegrals aa, bb, ab;
    SpatialHamiltonian H;
};

static void build_reference(Reference& r, const Molecule& mol, const BasisSet& b0, const std::string& ref, int mult,
                            const RHFResult* guess) {
    const BasisSet b = basis_at_geometry(b0, mol);
    RHFOptions o;
    o.reference = ref;
    o.multiplicity = mult;
    o.energy_tol = 1e-13;
    o.gradient_tol = 1e-11;
    o.max_iter = 400;
    if (guess) {
        o.guess_orbitals = &guess->C;
        if (ref == "UHF") o.guess_orbitals_b = &guess->Cb;
    }
    r.scf = run_scf(mol, b, o);
    r.aa = transform_to_mo(b, r.scf.H, r.scf.C, 0);
    if (ref == "UHF") {
        r.bb = transform_to_mo(b, r.scf.H, r.scf.Cb, 0);
        r.ab = transform_to_mo_mixed(b, r.scf.C, r.scf.Cb, 0);
        SpatialHamiltonian& H = r.H;
        H.n[0] = H.n[1] = r.aa.n;
        H.nocc[0] = r.scf.n_alpha;
        H.nocc[1] = r.scf.n_beta;
        H.h[0] = r.aa.h;
        H.h[1] = r.bb.h;
        H.aa = &r.aa;
        H.bb = &r.bb;
        H.ab = &r.ab;
        H.nuclear_repulsion = r.scf.nuclear_repulsion;
        H.restricted = false;
    } else {
        r.H = closed_shell_hamiltonian(r.aa, r.scf.n_occupied, r.scf.nuclear_repulsion);
        if (ref == "ROHF") {  // one set of orbitals, alpha and beta occupations
            r.H.nocc[0] = r.scf.n_alpha;
            r.H.nocc[1] = r.scf.n_beta;
        }
    }
}

static SOOptions tight() {
    SOOptions so;
    so.amplitude_tol = 1e-11;
    so.energy_tol = 1e-13;
    so.max_iter = 300;
    return so;
}

static Molecule molecule_from(const json& atoms, int charge) {
    std::vector<Atom> list;
    for (const auto& a : atoms) {
        Atom at;
        at.Z = atomic_number(a[0].get<std::string>());
        for (int k = 0; k < 3; ++k) at.position[static_cast<std::size_t>(k)] = a[1][static_cast<std::size_t>(k)].get<double>();
        list.push_back(at);
    }
    return make_custom_molecule(list, charge, "test");
}

static void test_energies() {
    const json data = json::parse(read_file(std::string(MACQC_TEST_DATA_DIR) + "/so_reference.json"));
    for (const auto& c : data["cases"]) {
        const std::string name = c["name"].get<std::string>();
        const std::string ref = c["reference"].get<std::string>();
        const int mult = c["multiplicity"].get<int>();
        Molecule mol = molecule_from(c["atoms"], c["charge"].get<int>());
        BasisSet b = build_basis(mol, basis_name(c["basis_file"].get<std::string>()), "");
        mol = attach_ecps(mol, b);
        Reference r;
        build_reference(r, mol, b, ref, mult, nullptr);
        check_near(name + " " + ref + " SCF", r.scf.total_energy, c["scf"].get<double>(), 1e-9);
        for (auto it = c["energies"].begin(); it != c["energies"].end(); ++it) {
            // the quadruples of LiH/6-31G take minutes in the dense spin-orbital engine (o^4 v^6 per
            // iteration); H2O/STO-3G, OH and CH2 test them
            const std::string key = it.key();
            if (name == "LiH" && key.size() > 2 && key.compare(key.size() - 2, 2, "TQ") == 0) {
                std::printf("skip %s %s (the quadruples are tested on H2O, OH, CH2)\n", name.c_str(), key.c_str());
                continue;
            }
            SOResult res = run_so_method(it.key(), r.H, ref, tight());
            check_near(name + " " + ref + " " + it.key(), res.total_energy, it.value().get<double>(), 2e-8);
        }
    }
}

static Molecule water(double shift) {
    // a distorted water molecule (no symmetry), positions in angstrom
    const double xyz[3][3] = {{0.0, 0.0, 0.02 + shift}, {0.80, 0.03, 0.55}, {-0.73, 0.0, 0.61}};
    const int Z[3] = {8, 1, 1};
    std::vector<Atom> atoms;
    for (int i = 0; i < 3; ++i) {
        Atom a;
        a.Z = Z[i];
        for (int k = 0; k < 3; ++k) a.position[static_cast<std::size_t>(k)] = xyz[i][k] / kBohrInAngstrom;
        atoms.push_back(a);
    }
    return make_custom_molecule(atoms, 0, "water");
}

// Four-point central difference of f along every Cartesian coordinate.
template <class F>
static std::vector<double> differences(const Molecule& mol, F f) {
    const double h = 0.005;
    std::vector<double> d;
    for (std::size_t x = 0; x < 3 * mol.atoms.size(); ++x) {
        double e[4];
        const double s[4] = {-2, -1, 1, 2};
        for (int k = 0; k < 4; ++k) {
            Molecule m = mol;
            m.atoms[x / 3].position[x % 3] += s[k] * h;
            e[k] = f(m);
        }
        d.push_back((e[0] - 8 * e[1] + 8 * e[2] - e[3]) / (12 * h));
    }
    return d;
}

static void test_closed_shell_derivatives() {
    Molecule mol = water(0.0);
    BasisSet b = build_basis(mol, "STO-3G", "");
    for (const char* method : {"MP3", "CISD(T)", "QCISD(T)", "QCISDT"}) {
        Reference r;
        build_reference(r, mol, b, "RHF", 1, nullptr);
        SOResult res = run_so_method(method, r.H, "RHF", tight());
        const CorrelatedGradientResult g = so_closed_shell_gradient(mol, b, r.scf, r.aa, res, tight());
        const std::vector<double> d = differences(mol, [&](const Molecule& m) {
            Reference rm;
            build_reference(rm, m, b, "RHF", 1, &r.scf);
            return run_so_method(method, rm.H, "RHF", tight()).total_energy;
        });
        double worst = 0.0;
        for (std::size_t x = 0; x < d.size(); ++x) worst = std::max(worst, std::fabs(d[x] - g.gradient[x]));
        check_small(std::string("H2O/STO-3G ") + method + " gradient - differences", worst, 1e-8);
        check_near(std::string("H2O/STO-3G ") + method + " relaxed-density energy", g.density_energy,
                   res.total_energy, 1e-9);
        if (std::string(method) == "QCISD(T)") {
            const CorrelatedHessianResult H = so_closed_shell_hessian(mol, b, r.scf, r.aa, res, tight());
            double wh = 0.0;
            const std::size_t n3 = d.size();
            for (std::size_t y = 0; y < n3; ++y) {
                const std::vector<double> col = differences(mol, [&](const Molecule& m) {
                    Reference rm;
                    build_reference(rm, m, b, "RHF", 1, &r.scf);
                    SOResult rs = run_so_method(method, rm.H, "RHF", tight());
                    return so_closed_shell_gradient(m, basis_at_geometry(b, m), rm.scf, rm.aa, rs, tight()).gradient[y];
                });
                for (std::size_t x = 0; x < n3; ++x) wh = std::max(wh, std::fabs(col[x] - H.hessian(y, x)));
            }
            check_small("H2O/STO-3G QCISD(T) Hessian - differences of gradients", wh, 1e-7);
        }
    }
}

static void test_uhf_derivatives() {
    // OH radical (doublet), distorted from the axis
    std::vector<Atom> atoms(2);
    atoms[0].Z = 8;
    atoms[1].Z = 1;
    atoms[1].position = {0.1, 0.05, 0.98 / kBohrInAngstrom};
    Molecule mol = make_custom_molecule(atoms, 0, "OH");
    BasisSet b = build_basis(mol, "STO-3G", "");
    for (const char* method : {"MP3", "QCISD", "CISD(T)"}) {
        Reference r;
        build_reference(r, mol, b, "UHF", 2, nullptr);
        SOResult res = run_so_method(method, r.H, "UHF", tight());
        const SODensities dens = so_densities(res, tight());
        const UHFCorrelatedGradient g = uhf_correlated_gradient(mol, b, r.scf, r.aa, r.bb, r.ab, dens);
        const std::vector<double> d = differences(mol, [&](const Molecule& m) {
            Reference rm;
            build_reference(rm, m, b, "UHF", 2, &r.scf);
            return run_so_method(method, rm.H, "UHF", tight()).total_energy;
        });
        double worst = 0.0;
        for (std::size_t x = 0; x < d.size(); ++x) worst = std::max(worst, std::fabs(d[x] - g.gradient[x]));
        check_small(std::string("OH/STO-3G U") + method + " gradient - differences", worst, 1e-8);
    }
}

// ---- 15.0: open-shell derivatives (UHF and ROHF references, open_shell_correlated.h)
static const MOIntegrals& beta_ints(const Reference& r) { return r.scf.reference == "UHF" ? r.bb : r.aa; }
static const MOIntegrals& mixed_ints(const Reference& r) { return r.scf.reference == "UHF" ? r.ab : r.aa; }

static double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    double w = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) w = std::max(w, std::fabs(a[i] - b[i]));
    return w;
}

static OpenShellGradient open_gradient(const Molecule& m, const BasisSet& b, const std::string& ref, int mult,
                                       const std::string& method, const RHFResult* guess) {
    Reference r;
    build_reference(r, m, b, ref, mult, guess);
    SOResult res = run_so_method(method, r.H, ref, tight());
    const SODensities d = so_densities(res, tight());
    return open_shell_correlated_gradient(m, basis_at_geometry(b, m), r.scf, r.aa, beta_ints(r), mixed_ints(r), d);
}

static void test_open_shell_derivatives() {
    std::vector<Atom> oh(2);
    oh[0].Z = 8;
    oh[1].Z = 1;
    oh[1].position = {0.1, 0.05, 0.98 / kBohrInAngstrom};
    const Molecule moh = make_custom_molecule(oh, 0, "OH");
    const BasisSet boh = build_basis(moh, "STO-3G", "");
    for (const char* ref : {"UHF", "ROHF"}) {
        const std::string R = ref;
        Reference r;
        build_reference(r, moh, boh, R, 2, nullptr);
        // the orbital part alone: the SCF energy as the method
        const RHFGradientResult g0 = scf_gradient(moh, boh, r.scf);
        const OpenShellGradient g1 = open_shell_scf_gradient_check(moh, boh, r.scf, r.aa, beta_ints(r), mixed_ints(r));
        check_small("OH/STO-3G " + R + " SCF: general gradient - scf_gradient", max_abs_diff(g0.gradient, g1.gradient), 1e-9);
        const RHFHessianResult h0 = scf_hessian(moh, boh, r.scf);
        const CorrelatedHessianResult h1 = open_shell_scf_hessian_check(moh, boh, r.scf, r.aa, beta_ints(r), mixed_ints(r));
        double wh = 0.0, wd = 0.0;
        for (std::size_t x = 0; x < h0.hessian.rows(); ++x)
            for (std::size_t y = 0; y < h0.hessian.cols(); ++y) wh = std::max(wh, std::fabs(h0.hessian(x, y) - h1.hessian(x, y)));
        for (std::size_t x = 0; x < h0.dipole_derivatives.size(); ++x)
            for (std::size_t k = 0; k < 3; ++k)
                wd = std::max(wd, std::fabs(h0.dipole_derivatives[x][k] - h1.dipole_derivatives[x][k]));
        check_small("OH/STO-3G " + R + " SCF: general Hessian - scf_hessian", wh, 1e-8);
        check_small("OH/STO-3G " + R + " SCF: dipole derivatives", wd, 1e-8);
    }
    // ROHF-based analytic gradients against differences of energies
    for (const char* method : {"MP2", "CCSD(T)"}) {
        Reference r;
        build_reference(r, moh, boh, "ROHF", 2, nullptr);
        const OpenShellGradient g = open_gradient(moh, boh, "ROHF", 2, method, nullptr);
        const std::vector<double> d = differences(moh, [&](const Molecule& m) {
            Reference rm;
            build_reference(rm, m, boh, "ROHF", 2, &r.scf);
            return run_so_method(method, rm.H, "ROHF", tight()).total_energy;
        });
        check_small(std::string("OH/STO-3G RO") + method + " gradient - differences", max_abs_diff(d, g.gradient), 1e-8);
        check_near(std::string("OH/STO-3G RO") + method + " relaxed-density energy", g.density_energy,
                   run_so_method(method, r.H, "ROHF", tight()).total_energy, 1e-9);
    }
    // UHF- and ROHF-based analytic Hessians against differences of analytic gradients
    std::vector<Atom> w(3);
    w[0].Z = 8;
    w[1].Z = 1;
    w[2].Z = 1;
    w[0].position = {0.0, 0.0, 0.02 / kBohrInAngstrom};
    w[1].position = {0.80 / kBohrInAngstrom, 0.03 / kBohrInAngstrom, 0.55 / kBohrInAngstrom};
    w[2].position = {-0.73 / kBohrInAngstrom, 0.0, 0.61 / kBohrInAngstrom};
    const Molecule mw = make_custom_molecule(w, 1, "H2O+");
    const BasisSet bw = build_basis(mw, "STO-3G", "");
    struct Case {
        const Molecule* mol;
        const BasisSet* basis;
        const char* name;
        const char* ref;
        const char* method;
    };
    const Case cases[] = {{&moh, &boh, "OH/STO-3G", "UHF", "CCSD"},
                          {&moh, &boh, "OH/STO-3G", "ROHF", "CISD"},
                          {&mw, &bw, "H2O+/STO-3G", "UHF", "MP3"},
                          {&mw, &bw, "H2O+/STO-3G", "ROHF", "QCISD(T)"}};
    for (const Case& c : cases) {
        const std::string label = std::string(c.name) + " " + (std::string(c.ref) == "UHF" ? "U" : "RO") + c.method;
        Reference r;
        build_reference(r, *c.mol, *c.basis, c.ref, 2, nullptr);
        SOResult res = run_so_method(c.method, r.H, c.ref, tight());
        const CorrelatedHessianResult H =
            open_shell_correlated_hessian(*c.mol, *c.basis, r.scf, r.aa, beta_ints(r), mixed_ints(r), res, tight());
        const std::size_t n3 = 3 * c.mol->atoms.size();
        const double h = 0.005;
        double wh = 0.0, wd = 0.0;
        for (std::size_t x = 0; x < n3; ++x) {
            std::vector<std::vector<double>> gs;
            std::vector<std::array<double, 3>> ds;
            const double f[4] = {-2, -1, 1, 2};
            for (int k = 0; k < 4; ++k) {
                Molecule m = *c.mol;
                m.atoms[x / 3].position[x % 3] += f[k] * h;
                const OpenShellGradient g = open_gradient(m, *c.basis, c.ref, 2, c.method, &r.scf);
                gs.push_back(g.gradient);
                ds.push_back(g.relaxed_dipole);
            }
            for (std::size_t y = 0; y < n3; ++y)
                wh = std::max(wh, std::fabs((gs[0][y] - 8 * gs[1][y] + 8 * gs[2][y] - gs[3][y]) / (12 * h) - H.hessian(y, x)));
            for (std::size_t k = 0; k < 3; ++k)
                wd = std::max(wd, std::fabs((ds[0][k] - 8 * ds[1][k] + 8 * ds[2][k] - ds[3][k]) / (12 * h) - H.dipole_derivatives[x][k]));
        }
        check_small(label + " Hessian - differences of gradients", wh, 1e-7);
        check_small(label + " dipole derivatives - differences", wd, 1e-7);
    }
}

// ---- 15.0: Full CI on open-shell references (orbital-invariant: no Z-vector)
struct FCIPoint {
    Reference r;
    FCIResultGeneral fci;
};

static void fci_point(FCIPoint& p, const Molecule& m, const BasisSet& b, const std::string& ref, const RHFResult* guess) {
    build_reference(p.r, m, b, ref, 2, guess);
    FCIOptions fo;
    fo.keep_vector = true;
    fo.residual_tol = 1e-10;
    fo.energy_tol = 1e-14;
    fo.max_iter = 300;
    p.fci = run_fci_determinants(p.r.aa, p.r.scf.n_alpha, p.r.scf.n_beta, p.r.scf.nuclear_repulsion, p.r.scf.total_energy, fo);
}

static void test_open_shell_fci() {
    std::vector<Atom> oh(2);
    oh[0].Z = 8;
    oh[1].Z = 1;
    oh[1].position = {0.1, 0.05, 0.98 / kBohrInAngstrom};
    const Molecule moh = make_custom_molecule(oh, 0, "OH");
    const BasisSet boh = build_basis(moh, "STO-3G", "");
    std::vector<Atom> w(3);
    w[0].Z = 8;
    w[1].Z = 1;
    w[2].Z = 1;
    w[0].position = {0.0, 0.0, 0.02 / kBohrInAngstrom};
    w[1].position = {0.80 / kBohrInAngstrom, 0.03 / kBohrInAngstrom, 0.55 / kBohrInAngstrom};
    w[2].position = {-0.73 / kBohrInAngstrom, 0.0, 0.61 / kBohrInAngstrom};
    const Molecule mw = make_custom_molecule(w, 1, "H2O+");
    const BasisSet bw = build_basis(mw, "STO-3G", "");
    // the gradient against differences of the Full CI energies
    {
        FCIPoint p;
        fci_point(p, moh, boh, "ROHF", nullptr);
        const OpenShellGradient g = open_shell_fci_gradient(moh, boh, p.r.scf, p.r.aa, p.fci.vector);
        check_near("OH/STO-3G FCI (ROHF orbitals) density energy", g.density_energy, p.fci.total_energy, 1e-10);
        check_small("OH/STO-3G FCI generalized Brillouin asymmetry", g.redundant_asymmetry, 1e-8);
        const std::vector<double> d = differences(moh, [&](const Molecule& m) {
            FCIPoint q;
            fci_point(q, m, boh, "ROHF", &p.r.scf);
            return q.fci.total_energy;
        });
        check_small("OH/STO-3G FCI gradient - differences", max_abs_diff(d, g.gradient), 1e-8);
    }
    // the Hessian and dipole derivatives against differences of analytic gradients
    struct Case {
        const Molecule* mol;
        const BasisSet* basis;
        const char* name;
        const char* ref;
    };
    const Case cases[] = {{&moh, &boh, "OH/STO-3G", "UHF"}, {&mw, &bw, "H2O+/STO-3G", "ROHF"}};
    for (const Case& c : cases) {
        const std::string label = std::string(c.name) + " FCI (" + c.ref + " orbitals)";
        FCIPoint p;
        fci_point(p, *c.mol, *c.basis, c.ref, nullptr);
        const CorrelatedHessianResult H =
            open_shell_fci_hessian(*c.mol, *c.basis, p.r.scf, p.r.aa, p.fci.vector, p.fci.electronic_energy);
        const std::size_t n3 = 3 * c.mol->atoms.size();
        const double h = 0.005;
        double wh = 0.0, wd = 0.0;
        for (std::size_t x = 0; x < n3; ++x) {
            std::vector<std::vector<double>> gs;
            std::vector<std::array<double, 3>> ds;
            const double f[4] = {-2, -1, 1, 2};
            for (int k = 0; k < 4; ++k) {
                Molecule m = *c.mol;
                m.atoms[x / 3].position[x % 3] += f[k] * h;
                FCIPoint q;
                fci_point(q, m, *c.basis, c.ref, &p.r.scf);
                const OpenShellGradient g = open_shell_fci_gradient(m, basis_at_geometry(*c.basis, m), q.r.scf, q.r.aa, q.fci.vector);
                gs.push_back(g.gradient);
                ds.push_back(g.relaxed_dipole);
            }
            for (std::size_t y = 0; y < n3; ++y)
                wh = std::max(wh, std::fabs((gs[0][y] - 8 * gs[1][y] + 8 * gs[2][y] - gs[3][y]) / (12 * h) - H.hessian(y, x)));
            for (std::size_t k = 0; k < 3; ++k)
                wd = std::max(wd, std::fabs((ds[0][k] - 8 * ds[1][k] + 8 * ds[2][k] - ds[3][k]) / (12 * h) - H.dipole_derivatives[x][k]));
        }
        check_small(label + " Hessian - differences of gradients", wh, 1e-7);
        check_small(label + " dipole derivatives - differences", wd, 1e-7);
        ++g_checks;
        const bool ok = H.method_converged;
        if (!ok) ++g_failures;
        std::printf("%s %-44s %d iterations\n", ok ? "ok  " : "FAIL", (label + " CI responses converged").c_str(), H.method_iterations);
    }
}

// Arguments (optional): the groups to run (energies, closed, uhf, open, fci); none: all.
int main(int argc, char** argv) {
    auto run = [&](const char* group) {
        if (argc < 2) return true;
        for (int i = 1; i < argc; ++i)
            if (std::string(argv[i]) == group) return true;
        return false;
    };
    try {
        if (run("energies")) test_energies();
        if (run("closed")) test_closed_shell_derivatives();
        if (run("uhf")) test_uhf_derivatives();
        if (run("open")) test_open_shell_derivatives();
        if (run("fci")) test_open_shell_fci();
    } catch (const std::exception& e) {
        std::printf("FAIL exception: %s\n", e.what());
        return 1;
    }
    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
