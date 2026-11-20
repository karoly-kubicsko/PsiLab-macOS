#include "macqc/gradient.h"
#include "macqc/open_shell_correlated.h"

#include <chrono>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <stdexcept>

#include "macqc/gto_derivatives.h"
#include "macqc/gto_integrals.h"
#include "macqc/internal_coordinates.h"
#include "macqc/job.h"
#include "macqc/mo_integrals.h"
#include "macqc/rhf_derivatives.h"
#include "macqc/scf_derivatives.h"
#include "macqc/so_derivatives.h"
#include "macqc/uhf_correlated.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string format(const char* fmt, ...) {
    char buf[1024];
    va_list args;
    va_start(args, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, args);
    va_end(args);
    return buf;
}

RHFOptions rhf_options(const ElectronicOptions& o, const Matrix* guess) {
    RHFOptions r;
    r.max_iter = o.scf_max_iter;
    r.energy_tol = o.scf_energy_tol;
    r.gradient_tol = o.scf_gradient_tol;
    r.lindep_tol = o.lindep_tol;
    r.threads = o.threads;
    r.memory_limit_mb = o.memory_mb;
    r.guess_orbitals = guess;
    r.field = o.field;
    r.reference = o.reference.empty() ? "RHF" : o.reference;
    r.multiplicity = o.multiplicity;
    return r;
}

SOOptions so_options(const ElectronicOptions& o) {
    SOOptions s;
    s.threads = o.threads;
    s.max_iter = o.cc_max_iter;
    s.amplitude_tol = o.cc_amplitude_tol;
    s.energy_tol = o.cc_energy_tol;
    s.memory_mb = o.fci_memory_mb;
    return s;
}

// The MO Hamiltonian of an SCF reference for the spin-orbital engine (kept together: H points into
// the integrals).
struct SOSetup {
    MOIntegrals aa, bb, ab;
    SpatialHamiltonian H;
    SOSetup(const RHFResult& scf, const BasisSet& basis, int threads) {
        aa = transform_to_mo(basis, scf.H, scf.C, threads);
        H.n[0] = H.n[1] = aa.n;
        H.nocc[0] = scf.n_alpha;
        H.nocc[1] = scf.n_beta;
        H.h[0] = aa.h;
        H.nuclear_repulsion = scf.nuclear_repulsion;
        if (scf.reference == "UHF") {
            bb = transform_to_mo(basis, scf.H, scf.Cb, threads);
            ab = transform_to_mo_mixed(basis, scf.C, scf.Cb, threads);
            H.h[1] = bb.h;
            H.aa = &aa;
            H.bb = &bb;
            H.ab = &ab;
            H.restricted = false;
        } else {
            H.h[1] = aa.h;
            H.aa = H.bb = H.ab = &aa;
            H.restricted = true;
        }
    }
    SOSetup(const SOSetup&) = delete;
    SOSetup& operator=(const SOSetup&) = delete;
};

// The method label with the reference (Gaussian's prefixes for open shells: UMP2, ROCCSD(T)).
std::string labelled(const std::string& method, const std::string& reference) {
    if (method == "RHF") return reference.empty() ? "RHF" : reference;
    if (reference == "UHF") return "U" + method;
    if (reference == "ROHF") return "RO" + method;
    return method;
}

CCOptions cc_options(const ElectronicOptions& o) {
    CCOptions c;
    c.threads = o.threads;
    c.memory_limit_mb = o.fci_memory_mb;
    c.eri_memory_mb = o.memory_mb;
    c.max_iter = o.cc_max_iter;
    c.energy_tol = o.cc_energy_tol;
    c.amplitude_tol = o.cc_amplitude_tol;
    return c;
}

FCIOptions fci_options(const ElectronicOptions& o, bool keep) {
    FCIOptions f;
    f.threads = o.threads;
    f.memory_limit_mb = o.fci_memory_mb;
    f.residual_tol = o.fci_residual_tol;
    f.energy_tol = 1.0e-12;
    f.max_iter = 200;
    f.keep_vector = keep;
    return f;
}

double norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

void check_rhf(SinglePoint& sp) {
    if (!sp.rhf.converged)
        sp.warnings.push_back(format("The SCF did not converge in %d iterations.", sp.rhf.iterations));
}

// Correlated energy after RHF (MO integrals given when available).
void correlated_energy(SinglePoint& sp, const BasisSet& basis, const ElectronicOptions& o,
                       const ElectronicCallbacks& cb, const MOIntegrals* mo_in) {
    const std::string& m = sp.method;
    const int nocc = sp.rhf.n_occupied;
    if (m == "RHF") {
        sp.energy = sp.rhf.total_energy;
        sp.converged = sp.rhf.converged;
        return;
    }
    if (!sp.rhf.converged) {
        sp.energy = sp.rhf.total_energy;
        sp.converged = false;
        return;
    }
    if (uses_so_engine(m, sp.rhf.reference)) {
        report_progress(JobStage::Transformation);
        const SOSetup setup(sp.rhf, basis, o.threads);
        SOProgress prog;
        if (cb.cc)
            prog = [&](const SOIteration& it) {
                CCIteration c;
                c.iteration = it.iteration;
                c.energy = it.energy;
                c.correlation = it.energy - sp.rhf.total_energy;
                c.delta_e = it.delta_e;
                c.residual = it.residual;
                cb.cc(c);
            };
        sp.so = run_so_method(m, setup.H, sp.rhf.reference, so_options(o), prog);
        sp.energy = sp.so->total_energy;
        sp.converged = sp.so->converged;
        if (!sp.so->converged)
            sp.warnings.push_back(format("The %s amplitude equations did not converge in %d iterations.",
                                         labelled(m, sp.rhf.reference).c_str(), sp.so->iterations));
        return;
    }
    if (m == "FCI" && sp.rhf.n_alpha != sp.rhf.n_beta) {
        // open shell: Full-CI in the (ROHF) orbitals of one set, alpha and beta strings
        MOIntegrals local;
        if (!mo_in) local = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
        const MOIntegrals& mo = mo_in ? *mo_in : local;
        sp.fci = run_fci_determinants(mo, sp.rhf.n_alpha, sp.rhf.n_beta, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                      fci_options(o, false), cb.fci);
        sp.energy = sp.fci->total_energy;
        sp.converged = sp.fci->converged;
        return;
    }
    if (m != "FCI" && sp.rhf.open_shell())
        throw std::runtime_error(m + " on an open-shell reference is computed by the spin-orbital engine");
    if (m == "FCI") {
        MOIntegrals local;
        if (!mo_in) local = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
        const MOIntegrals& mo = mo_in ? *mo_in : local;
        sp.fci = run_fci_determinants(mo, nocc, nocc, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                      fci_options(o, false), cb.fci);
        sp.energy = sp.fci->total_energy;
        sp.converged = sp.fci->converged;
        if (!sp.fci->converged)
            sp.warnings.push_back(format("The Full-CI (Davidson) did not converge in %d iterations.", sp.fci->iterations));
        return;
    }
    CCMethod ccm = CCMethod::CCSD;
    if (m == "CCSD(T)") ccm = CCMethod::CCSD_T;
    else if (m == "CCSDT") ccm = CCMethod::CCSDT;
    else if (m == "CCSDTQ") ccm = CCMethod::CCSDTQ;
    else if (m != "CCSD") throw std::invalid_argument("unknown method " + m);
    const CCOptions co = cc_options(o);
    if (cc_max_rank(ccm) == 2) {
        const CCIntegrals ints = mo_in ? cc_integrals_from_mo(*mo_in, nocc)
                                       : cc_integrals_from_ao(basis, sp.rhf.H, sp.rhf.C, nocc, o.threads);
        auto ladder = make_ao_ladder(basis, sp.rhf.C, nocc, o.threads, o.memory_mb);
        sp.cc = run_ccsd_closed_shell(ints, *ladder, sp.rhf.total_energy, ccm == CCMethod::CCSD_T, co, cb.cc);
    } else {
        MOIntegrals local;
        if (!mo_in) local = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
        const MOIntegrals& mo = mo_in ? *mo_in : local;
        sp.cc = run_cc_spin_orbital(mo, nocc, cc_max_rank(ccm), sp.rhf.nuclear_repulsion, co, cb.cc);
        sp.cc->reference_energy = sp.rhf.total_energy;
        sp.cc->cc_energy = sp.cc->reference_energy + sp.cc->correlation_energy;
        sp.cc->total_energy = sp.cc->cc_energy;
    }
    sp.energy = sp.cc->total_energy;
    sp.converged = sp.cc->converged;
    if (!sp.cc->converged)
        sp.warnings.push_back(format("The %s iterations did not converge in %d iterations.", sp.cc->method.c_str(),
                                     sp.cc->iterations));
}

}  // namespace

bool uses_so_engine(const std::string& m, const std::string& reference) {
    if (!so_engine_method(m)) return false;
    const bool cc = m == "CCSD" || m == "CCSD(T)" || m == "CCSDT" || m == "CCSDTQ";
    if (cc && (reference.empty() || reference == "RHF")) return false;  // closed shell: the CC engines of 4.0
    return true;
}

bool analytic_gradient_method(const std::string& m, const std::string& reference) {
    const bool closed = reference.empty() || reference == "RHF";
    if (m == "RHF") return true;  // the SCF of every reference
    if (closed) return m == "CCSD" || m == "CCSD(T)" || m == "FCI" || uses_so_engine(m, "RHF");
    // 15.0: UHF (coupled-perturbed UHF) and ROHF (ROHF orbital Hessian) Z-vectors, open_shell_correlated.h;
    // Full-CI (orbital-invariant: no Z-vector)
    if (reference == "UHF" || reference == "ROHF") return uses_so_engine(m, reference) || m == "FCI";
    return false;
}

double analytic_gradient_bytes(const std::string& m, int nocc, int nmo, int nbasis) {
    if (m == "RHF") return 0.0;
    const double base = correlated_gradient_bytes(static_cast<std::size_t>(nmo), static_cast<std::size_t>(nbasis));
    if (m == "CCSD" || m == "CCSD(T)") return base + ccsd_density_bytes(nocc, nmo - nocc, m == "CCSD(T)");
    if (m == "FCI") return base + fci_density_bytes(nmo, nocc);
    return 0.0;
}

SinglePoint run_single_point(const Molecule& mol, const BasisSet& basis, const std::string& method,
                             const ElectronicOptions& o, const Matrix* guess, const ElectronicCallbacks& cb) {
    const auto t0 = Clock::now();
    SinglePoint sp;
    sp.molecule = mol;
    sp.method = method;
    sp.reference = o.reference.empty() ? "RHF" : o.reference;
    sp.rhf = run_scf(mol, basis, rhf_options(o, guess), cb.scf);
    check_rhf(sp);
    correlated_energy(sp, basis, o, cb, nullptr);
    if (sp.so) sp.so->state.reset();  // energy only: no derivatives follow
    sp.seconds = seconds_since(t0);
    return sp;
}

SimpleInternals simple_internals(const Molecule& mol) {
    SimpleInternals s;
    const std::size_t n = mol.atoms.size();
    auto stretch = [&](std::size_t a, std::size_t b, const std::string& name) {
        double d[3];
        for (int k = 0; k < 3; ++k) d[k] = mol.atoms[b].position[static_cast<std::size_t>(k)] - mol.atoms[a].position[static_cast<std::size_t>(k)];
        const double r = norm3(d);
        std::vector<double> row(3 * n, 0.0);
        for (int k = 0; k < 3; ++k) {
            row[3 * a + static_cast<std::size_t>(k)] = -d[k] / r;
            row[3 * b + static_cast<std::size_t>(k)] = d[k] / r;
        }
        s.names.push_back(name);
        s.values.push_back(r);
        s.B.push_back(row);
    };
    if (n == 2) {
        stretch(0, 1, "R(" + mol.atoms[0].symbol + "1-" + mol.atoms[1].symbol + "2)");
        return s;
    }
    if (n == 4) {
        auto add = [&](InternalType t, std::vector<int> atoms) {
            InternalCoordinate c;
            c.type = t;
            c.atoms = std::move(atoms);
            s.names.push_back(c.definition(mol));
            s.values.push_back(internal_value(mol, c));
            s.B.push_back(wilson_b_row(mol, c));
        };
        if (mol.atoms[0].Z == 7) {
            // NH3: N1, H2, H3, H4
            add(InternalType::Stretch, {0, 1});
            add(InternalType::Stretch, {0, 2});
            add(InternalType::Stretch, {0, 3});
            add(InternalType::Bend, {1, 0, 2});
            add(InternalType::Bend, {1, 0, 3});
            add(InternalType::Dihedral, {3, 0, 1, 2});
        } else {
            // H2O2: O1, O2, H3 (on O1), H4 (on O2)
            add(InternalType::Stretch, {0, 1});
            add(InternalType::Stretch, {0, 2});
            add(InternalType::Stretch, {1, 3});
            add(InternalType::Bend, {1, 0, 2});
            add(InternalType::Bend, {0, 1, 3});
            add(InternalType::Dihedral, {2, 0, 1, 3});
        }
        return s;
    }
    if (n != 3) throw std::invalid_argument("simple_internals: two, three or four atoms");
    stretch(0, 1, "R(O1–H2)");
    stretch(0, 2, "R(O1–H3)");
    double u[3], v[3];
    for (int k = 0; k < 3; ++k) {
        u[k] = mol.atoms[1].position[static_cast<std::size_t>(k)] - mol.atoms[0].position[static_cast<std::size_t>(k)];
        v[k] = mol.atoms[2].position[static_cast<std::size_t>(k)] - mol.atoms[0].position[static_cast<std::size_t>(k)];
    }
    const double ru = norm3(u), rv = norm3(v);
    for (int k = 0; k < 3; ++k) {
        u[k] /= ru;
        v[k] /= rv;
    }
    const double c = std::max(-1.0, std::min(1.0, u[0] * v[0] + u[1] * v[1] + u[2] * v[2]));
    const double th = std::acos(c), sn = std::sin(th);
    std::vector<double> row(9, 0.0);
    if (sn > 1.0e-8)
        for (int k = 0; k < 3; ++k) {
            const double da = (c * u[k] - v[k]) / (ru * sn), db = (c * v[k] - u[k]) / (rv * sn);
            row[3 + static_cast<std::size_t>(k)] = da;
            row[6 + static_cast<std::size_t>(k)] = db;
            row[static_cast<std::size_t>(k)] = -(da + db);
        }
    s.names.push_back("A(H2–O1–H3)");
    s.values.push_back(th);
    s.B.push_back(row);
    return s;
}

bool uses_cartesian_displacements(const Molecule& mol) { return mol.custom && mol.atoms.size() >= 3; }

RigidMotion rigid_motions(const Molecule& mol, const std::array<double, 3>& mu, const std::vector<double>* grad) {
    RigidMotion rm;
    const std::size_t n = mol.atoms.size(), n3 = 3 * n;
    double c[3] = {0.0, 0.0, 0.0};
    for (const Atom& a : mol.atoms)
        for (std::size_t k = 0; k < 3; ++k) c[k] += a.position[k] / static_cast<double>(n);
    const double Q = static_cast<double>(mol.charge);
    // translations: (1/sqrt N) e_k on every atom; d mu / ds = Q / sqrt N e_k
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<double> t(n3, 0.0);
        for (std::size_t a = 0; a < n; ++a) t[3 * a + k] = 1.0 / std::sqrt(static_cast<double>(n));
        std::array<double, 3> d{0.0, 0.0, 0.0};
        d[k] = Q / std::sqrt(static_cast<double>(n));
        rm.directions.push_back(t);
        rm.dipole_derivatives.push_back(d);
        if (grad) rm.hessian_times.push_back(std::vector<double>(n3, 0.0));
    }
    // rotations about the centroid (orthogonal to the translations), Gram-Schmidt; the dipole
    // derivatives follow the same linear combinations
    double largest = 0.0;
    std::vector<std::vector<double>> raw, raw_h;
    std::vector<std::array<double, 3>> raw_mu;
    for (std::size_t k = 0; k < 3; ++k) {
        double e[3] = {0.0, 0.0, 0.0};
        e[k] = 1.0;
        std::vector<double> r(n3, 0.0);
        for (std::size_t a = 0; a < n; ++a) {
            double p[3];
            for (std::size_t i = 0; i < 3; ++i) p[i] = mol.atoms[a].position[i] - c[i];
            r[3 * a + 0] = e[1] * p[2] - e[2] * p[1];
            r[3 * a + 1] = e[2] * p[0] - e[0] * p[2];
            r[3 * a + 2] = e[0] * p[1] - e[1] * p[0];
        }
        const double m[3] = {mu[0] - Q * c[0], mu[1] - Q * c[1], mu[2] - Q * c[2]};
        raw.push_back(r);
        raw_mu.push_back({e[1] * m[2] - e[2] * m[1], e[2] * m[0] - e[0] * m[2], e[0] * m[1] - e[1] * m[0]});
        std::vector<double> hr(n3, 0.0);
        if (grad)
            for (std::size_t a = 0; a < n; ++a) {
                const double* g = &(*grad)[3 * a];
                hr[3 * a + 0] = e[1] * g[2] - e[2] * g[1];
                hr[3 * a + 1] = e[2] * g[0] - e[0] * g[2];
                hr[3 * a + 2] = e[0] * g[1] - e[1] * g[0];
            }
        raw_h.push_back(hr);
        double s = 0.0;
        for (double v : r) s += v * v;
        largest = std::max(largest, std::sqrt(s));
    }
    std::vector<std::vector<double>> rot;
    for (std::size_t k = 0; k < 3; ++k) {
        std::vector<double> r = raw[k], hr = raw_h[k];
        std::array<double, 3> d = raw_mu[k];
        for (std::size_t p = 0; p < rot.size(); ++p) {
            double s = 0.0;
            for (std::size_t i = 0; i < n3; ++i) s += rot[p][i] * r[i];
            for (std::size_t i = 0; i < n3; ++i) r[i] -= s * rot[p][i];
            for (std::size_t i = 0; i < 3; ++i) d[i] -= s * rm.dipole_derivatives[3 + p][i];
            if (grad)
                for (std::size_t i = 0; i < n3; ++i) hr[i] -= s * rm.hessian_times[3 + p][i];
        }
        double s = 0.0;
        for (double v : r) s += v * v;
        s = std::sqrt(s);
        if (s < 1.0e-6 * std::max(1.0, largest)) continue;  // linear molecule: rotation about the axis
        for (double& v : r) v /= s;
        for (double& v : d) v /= s;
        for (double& v : hr) v /= s;
        rot.push_back(r);
        rm.directions.push_back(r);
        rm.dipole_derivatives.push_back(d);
        if (grad) rm.hessian_times.push_back(hr);
    }
    return rm;
}

std::vector<std::vector<double>> vibrational_directions(const Molecule& mol, bool* linear) {
    const std::size_t n3 = 3 * mol.atoms.size();
    const RigidMotion rm = rigid_motions(mol, {0.0, 0.0, 0.0});
    if (linear) *linear = rm.directions.size() == 5;
    std::vector<std::vector<double>> basis = rm.directions, out;
    const std::size_t want = n3 > rm.directions.size() ? n3 - rm.directions.size() : 0;
    for (std::size_t x = 0; x < n3 && out.size() < want; ++x) {
        std::vector<double> v(n3, 0.0);
        v[x] = 1.0;
        for (int pass = 0; pass < 2; ++pass)  // twice for numerical orthogonality
            for (const auto& b : basis) {
                double s = 0.0;
                for (std::size_t i = 0; i < n3; ++i) s += b[i] * v[i];
                for (std::size_t i = 0; i < n3; ++i) v[i] -= s * b[i];
            }
        double s = 0.0;
        for (double w : v) s += w * w;
        s = std::sqrt(s);
        if (s < 1.0e-3) continue;
        for (double& w : v) w /= s;
        basis.push_back(v);
        out.push_back(v);
    }
    return out;
}

Molecule geometry_from_internals(const Molecule& ref, const std::vector<double>& q) {
    Molecule m = ref;
    const std::size_t n = ref.atoms.size();
    double e1[3];
    for (int k = 0; k < 3; ++k) e1[k] = ref.atoms[1].position[static_cast<std::size_t>(k)] - ref.atoms[0].position[static_cast<std::size_t>(k)];
    const double r = norm3(e1);
    for (double& x : e1) x /= r;
    const auto& O = ref.atoms[0].position;
    if (n == 2) {
        for (int k = 0; k < 3; ++k) m.atoms[1].position[static_cast<std::size_t>(k)] = O[static_cast<std::size_t>(k)] + q[0] * e1[k];
        return m;
    }
    if (n == 4) {
        // A Z-matrix in the frame of the reference: atom 0 fixed, atom 1 along the reference
        // 0-1 direction e1; p is the unit vector perpendicular to e1 towards the reference atom 2,
        // t = e1 x p. A right-handed rotation about e1 by tau is the IUPAC dihedral angle seen
        // along 0 -> 1.
        double p[3], w[3];
        for (int k = 0; k < 3; ++k) w[k] = ref.atoms[2].position[static_cast<std::size_t>(k)] - O[static_cast<std::size_t>(k)];
        const double proj = w[0] * e1[0] + w[1] * e1[1] + w[2] * e1[2];
        for (int k = 0; k < 3; ++k) p[k] = w[k] - proj * e1[k];
        double pn = norm3(p);
        if (pn < 1.0e-10) {
            const double a[3] = {std::fabs(e1[0]) < 0.9 ? 1.0 : 0.0, std::fabs(e1[0]) < 0.9 ? 0.0 : 1.0, 0.0};
            const double d = a[0] * e1[0] + a[1] * e1[1] + a[2] * e1[2];
            for (int k = 0; k < 3; ++k) p[k] = a[k] - d * e1[k];
            pn = norm3(p);
        }
        for (double& x : p) x /= pn;
        const double t[3] = {e1[1] * p[2] - e1[2] * p[1], e1[2] * p[0] - e1[0] * p[2], e1[0] * p[1] - e1[1] * p[0]};
        const bool ammonia = ref.atoms[0].Z == 7;
        for (int k = 0; k < 3; ++k) {
            const std::size_t kk = static_cast<std::size_t>(k);
            const double a1 = O[kk] + q[0] * e1[k];
            m.atoms[1].position[kk] = a1;
            m.atoms[2].position[kk] = O[kk] + q[1] * (std::cos(q[3]) * e1[k] + std::sin(q[3]) * p[k]);
            if (ammonia) {
                // H4 on N1: angle H2-N1-H4 = q4, D(H4-N1-H2-H3) = q5 (H3 is H4 rotated by +q5)
                const double perp = std::cos(q[5]) * p[k] - std::sin(q[5]) * t[k];
                m.atoms[3].position[kk] = O[kk] + q[2] * (std::cos(q[4]) * e1[k] + std::sin(q[4]) * perp);
            } else {
                // H4 on O2: angle O1-O2-H4 = q4, D(H3-O1-O2-H4) = q5
                const double perp = std::cos(q[5]) * p[k] + std::sin(q[5]) * t[k];
                m.atoms[3].position[kk] = a1 + q[2] * (-std::cos(q[4]) * e1[k] + std::sin(q[4]) * perp);
            }
        }
        return m;
    }
    double w[3];
    for (int k = 0; k < 3; ++k) w[k] = ref.atoms[2].position[static_cast<std::size_t>(k)] - O[static_cast<std::size_t>(k)];
    // in-plane unit vector perpendicular to e1, towards atom 3
    const double proj = w[0] * e1[0] + w[1] * e1[1] + w[2] * e1[2];
    double p[3];
    for (int k = 0; k < 3; ++k) p[k] = w[k] - proj * e1[k];
    double pn = norm3(p);
    if (pn < 1.0e-10) {  // linear: any perpendicular direction
        const double a[3] = {std::fabs(e1[0]) < 0.9 ? 1.0 : 0.0, std::fabs(e1[0]) < 0.9 ? 0.0 : 1.0, 0.0};
        const double d = a[0] * e1[0] + a[1] * e1[1] + a[2] * e1[2];
        for (int k = 0; k < 3; ++k) p[k] = a[k] - d * e1[k];
        pn = norm3(p);
    }
    for (double& x : p) x /= pn;
    for (int k = 0; k < 3; ++k) {
        m.atoms[1].position[static_cast<std::size_t>(k)] = O[static_cast<std::size_t>(k)] + q[0] * e1[k];
        m.atoms[2].position[static_cast<std::size_t>(k)] =
            O[static_cast<std::size_t>(k)] + q[1] * (std::cos(q[2]) * e1[k] + std::sin(q[2]) * p[k]);
    }
    return m;
}

GradientResult run_gradient(const Molecule& mol, const BasisSet& basis, const std::string& method,
                            const ElectronicOptions& o, const Matrix* guess, const ElectronicCallbacks& cb,
                            bool force_numerical) {
    const auto t0 = Clock::now();
    GradientResult g;
    const std::size_t n3 = 3 * mol.atoms.size();
    g.gradient.assign(n3, 0.0);
    SinglePoint& sp = g.point;
    sp.molecule = mol;
    sp.method = method;
    sp.reference = o.reference.empty() ? "RHF" : o.reference;
    sp.rhf = run_scf(mol, basis, rhf_options(o, guess), cb.scf);
    check_rhf(sp);
    const int nocc = sp.rhf.n_occupied, nmo = sp.rhf.n_mo, nbf = basis.n_functions;
    const std::string& ref = sp.rhf.reference;
    const bool so = uses_so_engine(method, ref);

    // ---- analytic or numerical?
    std::string reason;
    if (force_numerical) reason = "requested";
    else if (!analytic_gradient_method(method, ref))
        reason = "no analytic gradient for " + labelled(method, ref) + " in this version";
    else if (method != "RHF" && sp.rhf.removed_functions > 0)
        reason = format("%d near-linear dependencies were removed; the analytic %s gradient does not contain the "
                        "rotation of the kept orbital space",
                        sp.rhf.removed_functions, method.c_str());
    else if (method == "FCI" && sp.rhf.open_shell()) {
        const double need = open_shell_fci_gradient_bytes(sp.rhf.n_alpha, sp.rhf.n_beta, nmo, nbf) / (1024.0 * 1024.0);
        if (need > o.fci_memory_mb)
            reason = format("the analytic gradient would need about %.1f GB (limit %.1f GB)", need / 1024.0,
                            o.fci_memory_mb / 1024.0);
    } else if (method != "RHF" && so) {
        const int on[2] = {sp.rhf.n_alpha, sp.rhf.n_beta};
        const int vn[2] = {nmo - sp.rhf.n_alpha, nmo - sp.rhf.n_beta};
        const double need = (1.6 * estimate_so_method(method, on, vn).memory_bytes +
                             correlated_gradient_bytes(static_cast<std::size_t>(nmo), static_cast<std::size_t>(nbf)) * (ref == "UHF" ? 3.0 : 1.0)) /
                            (1024.0 * 1024.0);
        if (need > o.fci_memory_mb)
            reason = format("the analytic gradient would need about %.1f GB (limit %.1f GB)", need / 1024.0,
                            o.fci_memory_mb / 1024.0);
    } else if (method != "RHF") {
        const double need = analytic_gradient_bytes(method, nocc, nmo, nbf) / (1024.0 * 1024.0);
        if (need > o.fci_memory_mb)
            reason = format("the analytic gradient would need about %.1f GB (limit %.1f GB)", need / 1024.0,
                            o.fci_memory_mb / 1024.0);
    }

    if (reason.empty() && sp.rhf.converged) {
        g.analytic = true;
        if (method == "RHF") {
            sp.energy = sp.rhf.total_energy;
            sp.converged = true;
            g.rhf_gradient = scf_gradient(mol, basis, sp.rhf, o.threads, o.lindep_tol);
            g.gradient = g.rhf_gradient->gradient;
            g.description = "analytic (" + ref + ": derivative integrals, energy-weighted density" +
                            std::string(sp.rhf.removed_functions ? ", rotation of the kept orbital space" : "") + ")";
            g.dipole = dipole_moment(mol, basis, sp.rhf.P);
            g.has_dipole = true;
        } else if (method == "FCI" && sp.rhf.open_shell()) {
            // 15.0: Full-CI on an open-shell reference (orbitals of one set: ROHF, or the alpha orbitals of
            // UHF); the energy is invariant to the orbital rotations, so no Z-vector
            report_progress(JobStage::Transformation);
            const MOIntegrals mo = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
            FCIOptions fo = fci_options(o, true);
            fo.residual_tol = std::min(fo.residual_tol, 1.0e-9);
            sp.fci = run_fci_determinants(mo, sp.rhf.n_alpha, sp.rhf.n_beta, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                          fo, cb.fci);
            sp.energy = sp.fci->total_energy;
            sp.converged = sp.fci->converged;
            if (!sp.fci->converged)
                sp.warnings.push_back(format("The Full-CI (Davidson) did not converge in %d iterations.", sp.fci->iterations));
            if (sp.converged) {
                if (cb.note) cb.note("Full-CI spin densities from the CI vector");
                const OpenShellGradient og = open_shell_fci_gradient(mol, basis, sp.rhf, mo, sp.fci->vector, o.threads, o.memory_mb);
                CorrelatedGradientResult cg;
                cg.gradient = og.gradient;
                cg.nuclear = og.nuclear;
                cg.one_electron = og.one_electron;
                cg.two_electron = og.two_electron;
                cg.overlap = og.overlap;
                cg.relaxed_density_ao = og.relaxed_density_ao;
                cg.relaxed_dipole = og.relaxed_dipole;
                cg.orbital_lagrangian_oo_vv = og.redundant_asymmetry;
                cg.zvector_converged = true;
                cg.density_energy = og.density_energy;
                cg.seconds = og.seconds;
                g.correlated = cg;
                g.gradient = og.gradient;
                g.dipole = og.relaxed_dipole;
                g.has_dipole = true;
            }
            sp.fci->vector.clear();
            sp.fci->vector.shrink_to_fit();
            g.description = "analytic (Full-CI on " + ref +
                            " orbitals: spin-resolved one- and two-particle densities of the CI vector; the energy is "
                            "invariant to the orbital rotations, so there is no Z-vector)";
        } else if (so) {
            report_progress(JobStage::Transformation);
            const SOSetup setup(sp.rhf, basis, o.threads);
            const SOOptions sopt = so_options(o);
            SOProgress prog;
            if (cb.cc)
                prog = [&](const SOIteration& it) {
                    if (it.stage == "Lambda") {
                        if (cb.lambda) {
                            LambdaIteration l;
                            l.iteration = it.iteration;
                            l.residual = it.residual;
                            cb.lambda(l);
                        }
                        return;
                    }
                    CCIteration c;
                    c.iteration = it.iteration;
                    c.energy = it.energy;
                    c.correlation = it.energy - sp.rhf.total_energy;
                    c.delta_e = it.delta_e;
                    c.residual = it.residual;
                    cb.cc(c);
                };
            sp.so = run_so_method(method, setup.H, ref, sopt, prog);
            sp.energy = sp.so->total_energy;
            sp.converged = sp.so->converged;
            if (sp.converged) {
                if (cb.note) cb.note(labelled(method, ref) + " Lambda equations and relaxed densities (spin-orbital engine)");
                if (ref == "UHF") {
                    const SODensities d = so_densities(*sp.so, sopt, prog);
                    g.lambda_converged = d.lambda_converged;
                    g.lambda_iterations = d.lambda_iterations;
                    g.lambda_residual = d.lambda_residual;
                    const UHFCorrelatedGradient ug =
                        uhf_correlated_gradient(mol, basis, sp.rhf, setup.aa, setup.bb, setup.ab, d, o.threads, o.memory_mb, cb.cphf);
                    CorrelatedGradientResult cg;
                    cg.gradient = ug.gradient;
                    cg.nuclear = ug.nuclear;
                    cg.one_electron = ug.one_electron;
                    cg.two_electron = ug.two_electron;
                    cg.overlap = ug.overlap;
                    cg.relaxed_density_ao = ug.relaxed_density_ao;
                    cg.relaxed_dipole = ug.relaxed_dipole;
                    cg.orbital_lagrangian_oo_vv = ug.orbital_lagrangian_oo_vv;
                    cg.zvector_converged = ug.zvector_converged;
                    cg.zvector_iterations = ug.zvector_iterations;
                    cg.zvector_residual = ug.zvector_residual;
                    cg.density_energy = ug.density_energy;
                    cg.seconds = ug.seconds;
                    g.correlated = cg;
                } else if (ref == "ROHF") {
                    // 15.0: the ROHF Z-vector (non-redundant closed-open, closed-virtual and open-virtual
                    // rotations of the common orbitals)
                    const SODensities d = so_densities(*sp.so, sopt, prog);
                    g.lambda_converged = d.lambda_converged;
                    g.lambda_iterations = d.lambda_iterations;
                    g.lambda_residual = d.lambda_residual;
                    const OpenShellGradient og = open_shell_correlated_gradient(mol, basis, sp.rhf, setup.aa, setup.aa, setup.aa, d,
                                                                                o.threads, o.memory_mb, cb.cphf);
                    CorrelatedGradientResult cg;
                    cg.gradient = og.gradient;
                    cg.nuclear = og.nuclear;
                    cg.one_electron = og.one_electron;
                    cg.two_electron = og.two_electron;
                    cg.overlap = og.overlap;
                    cg.relaxed_density_ao = og.relaxed_density_ao;
                    cg.relaxed_dipole = og.relaxed_dipole;
                    cg.orbital_lagrangian_oo_vv = og.redundant_asymmetry;
                    cg.zvector_converged = og.zvector_converged;
                    cg.zvector_iterations = og.zvector_iterations;
                    cg.zvector_residual = og.zvector_residual;
                    cg.density_energy = og.density_energy;
                    cg.seconds = og.seconds;
                    g.correlated = cg;
                } else {
                    SODensities d;
                    g.correlated = so_closed_shell_gradient(mol, basis, sp.rhf, setup.aa, *sp.so, sopt, &d, o.threads,
                                                            o.memory_mb, cb.cphf, prog);
                    g.lambda_converged = d.lambda_converged;
                    g.lambda_iterations = d.lambda_iterations;
                    g.lambda_residual = d.lambda_residual;
                }
                if (!g.lambda_converged)
                    sp.warnings.push_back(format("The %s Lambda equations did not converge in %d iterations.",
                                                 labelled(method, ref).c_str(), g.lambda_iterations));
                g.gradient = g.correlated->gradient;
                g.dipole = g.correlated->relaxed_dipole;
                g.has_dipole = true;
                if (!g.correlated->zvector_converged) sp.warnings.push_back("The Z-vector equations did not converge.");
            }
            g.description = "analytic (" + labelled(method, ref) +
                            ": Lagrangian of the spin-orbital engine differentiated automatically, relaxed one- and "
                            "two-particle densities, orbital response by Z-vector" +
                            std::string(ref == "UHF" ? " (coupled-perturbed UHF)" : ref == "ROHF" ? " (ROHF orbital Hessian)" : "") + ")";
            g.seconds = seconds_since(t0);
            return g;
        } else {
            report_progress(JobStage::Transformation);
            const MOIntegrals mo = transform_to_mo(basis, sp.rhf.H, sp.rhf.C, o.threads);
            CorrelatedDensity dens;
            if (method == "FCI") {
                sp.fci = run_fci_determinants(mo, nocc, nocc, sp.rhf.nuclear_repulsion, sp.rhf.total_energy,
                                              fci_options(o, true), cb.fci);
                sp.energy = sp.fci->total_energy;
                sp.converged = sp.fci->converged;
                if (cb.note) cb.note("Full-CI densities from the CI vector");
                dens = fci_densities(nmo, nocc, nocc, sp.fci->vector, o.threads);
                sp.fci->vector.clear();
                sp.fci->vector.shrink_to_fit();
                g.description = "analytic (Full-CI one- and two-particle densities, orbital response by Z-vector)";
            } else {
                const bool triples = method == "CCSD(T)";
                const CCIntegrals ints = cc_integrals_from_mo(mo, nocc);
                auto ladder = make_ao_ladder(basis, sp.rhf.C, nocc, o.threads, o.memory_mb);
                const CCOptions co = cc_options(o);
                sp.cc = run_ccsd_closed_shell(ints, *ladder, sp.rhf.total_energy, triples, co, cb.cc);
                sp.energy = sp.cc->total_energy;
                sp.converged = sp.cc->converged;
                if (sp.cc->converged) {
                    CCDensityResult d = ccsd_densities(ints, *ladder, sp.cc->t1, sp.cc->t2, triples, co, cb.lambda);
                    g.lambda_converged = d.lambda_converged;
                    g.lambda_iterations = d.lambda_iterations;
                    g.lambda_residual = d.lambda_residual;
                    g.l1 = std::move(d.l1);
                    g.l2 = std::move(d.l2);
                    if (!d.lambda_converged)
                        sp.warnings.push_back(format("The %s Lambda equations did not converge in %d iterations.",
                                                     method.c_str(), d.lambda_iterations));
                    dens = std::move(d.density);
                }
                g.description = std::string("analytic (") + method +
                                " Lambda equations, relaxed one- and two-particle densities, orbital response by "
                                "Z-vector)";
            }
            if (sp.converged) {
                g.correlated = correlated_gradient(mol, basis, sp.rhf, mo, dens, o.threads, o.memory_mb, cb.cphf);
                g.gradient = g.correlated->gradient;
                g.dipole = g.correlated->relaxed_dipole;
                g.has_dipole = true;
                if (!g.correlated->zvector_converged)
                    sp.warnings.push_back("The Z-vector equations did not converge.");
            }
        }
        g.seconds = seconds_since(t0);
        return g;
    }

    // ---- numerical: central differences of the energy along internal coordinates
    correlated_energy(sp, basis, o, cb, nullptr);
    g.numerical_reason = reason;
    if (!sp.converged) {
        g.seconds = seconds_since(t0);
        return g;
    }
    const double h = o.numerical_step;
    SimpleInternals si;
    bool linear = false;
    if (!uses_cartesian_displacements(mol)) {
        si = simple_internals(mol);
        linear = mol.atoms.size() == 3 && std::fabs(std::sin(si.values[2])) < 1.0e-6;
    }
    auto energy_at = [&](const Molecule& m) {
        const SinglePoint p = run_single_point(m, basis_at_geometry(basis, m), method, o, &sp.rhf.C, ElectronicCallbacks());
        ++g.energy_evaluations;
        if (!p.converged) throw std::runtime_error("a displaced calculation of the numerical gradient did not converge");
        return p.energy;
    };
    if (uses_cartesian_displacements(mol)) {
        // 9.0, custom molecules: central differences along the vibrational Cartesian directions
        const std::vector<std::vector<double>> dirs = vibrational_directions(mol);
        for (std::size_t k = 0; k < dirs.size(); ++k) {
            Molecule mp = mol, mm = mol;
            for (std::size_t x = 0; x < n3; ++x) {
                mp.atoms[x / 3].position[x % 3] += h * dirs[k][x];
                mm.atoms[x / 3].position[x % 3] -= h * dirs[k][x];
            }
            const double ep = energy_at(mp), em = energy_at(mm);
            const double dEds = (ep - em) / (2.0 * h);
            for (std::size_t x = 0; x < n3; ++x) g.gradient[x] += dEds * dirs[k][x];
            if (cb.note) cb.note(format("  numerical gradient: direction %zu of %zu, dE/ds = %+.8f", k + 1, dirs.size(), dEds));
        }
        g.description = format("numerical (central differences of %s energies along %zu Cartesian directions orthogonal "
                               "to the translations and rotations, %d energies, step %.4f bohr)",
                               method.c_str(), dirs.size(), g.energy_evaluations, h);
        g.seconds = seconds_since(t0);
        return g;
    }
    if (linear) {
        // Cartesian central differences (the bend of a linear triatomic is degenerate).
        for (std::size_t x = 0; x < n3; ++x) {
            Molecule mp = mol, mm = mol;
            mp.atoms[x / 3].position[x % 3] += h;
            mm.atoms[x / 3].position[x % 3] -= h;
            const double ep = energy_at(mp), em = energy_at(mm);
            g.gradient[x] = (ep - em) / (2.0 * h);
            if (cb.note) cb.note(format("  numerical gradient: Cartesian coordinate %zu of %zu", x + 1, n3));
        }
        g.description = format("numerical (central differences of %s energies, %d Cartesian displacements of %.4f bohr)",
                               method.c_str(), g.energy_evaluations, h);
        g.seconds = seconds_since(t0);
        return g;
    }
    const std::size_t nq = si.values.size();
    std::vector<double> gq(nq, 0.0);
    const bool peroxide = nq == 6 && mol.atoms[0].Z == 8;
    const bool symmetric = (nq == 3 && std::fabs(si.values[0] - si.values[1]) < 1.0e-9) ||
                           (peroxide && std::fabs(si.values[1] - si.values[2]) < 1.0e-9 &&
                            std::fabs(si.values[3] - si.values[4]) < 1.0e-9);
    struct Disp {
        std::string name;
        std::vector<double> dir;  // displacement of (q1, q2, ...) per unit step
    };
    std::vector<Disp> disps;
    if (nq == 1) {
        disps.push_back({si.names[0], {1.0}});
    } else if (symmetric && nq == 3) {
        const double s = 1.0 / std::sqrt(2.0);
        disps.push_back({"S1 = [R(O1–H2) + R(O1–H3)]/sqrt(2) (symmetric stretch)", {s, s, 0.0}});
        disps.push_back({"S2 = A(H2–O1–H3) (bend)", {0.0, 0.0, 1.0}});
    } else if (symmetric) {
        // C2-symmetric H2O2: the antisymmetric (B) derivatives vanish by symmetry
        const double s = 1.0 / std::sqrt(2.0);
        disps.push_back({"S1 = R(O1–O2) (O–O stretch)", {1, 0, 0, 0, 0, 0}});
        disps.push_back({"S2 = [R(O1–H3) + R(O2–H4)]/sqrt(2) (symmetric O–H stretch)", {0, s, s, 0, 0, 0}});
        disps.push_back({"S3 = [A(O2–O1–H3) + A(O1–O2–H4)]/sqrt(2) (symmetric bend)", {0, 0, 0, s, s, 0}});
        disps.push_back({"S4 = D(H3–O1–O2–H4) (torsion)", {0, 0, 0, 0, 0, 1}});
    } else {
        for (std::size_t k = 0; k < nq; ++k) {
            std::vector<double> d(nq, 0.0);
            d[k] = 1.0;
            disps.push_back({si.names[k], d});
        }
    }
    for (std::size_t k = 0; k < disps.size(); ++k) {
        std::vector<double> qp = si.values, qm = si.values;
        for (std::size_t j = 0; j < nq; ++j) {
            qp[j] += h * disps[k].dir[j];
            qm[j] -= h * disps[k].dir[j];
        }
        const double ep = energy_at(geometry_from_internals(mol, qp));
        const double em = energy_at(geometry_from_internals(mol, qm));
        const double dEdS = (ep - em) / (2.0 * h);
        // dE/dq_j = sum_k dE/dS_k dS_k/dq_j with orthonormal directions: dS_k/dq = dir_k
        for (std::size_t j = 0; j < nq; ++j) gq[j] += dEdS * disps[k].dir[j];
        g.internal_coordinates.push_back(disps[k].name);
        if (cb.note) cb.note(format("  numerical gradient: %s, dE/dq = %+.8f", disps[k].name.c_str(), dEdS));
    }
    for (std::size_t x = 0; x < n3; ++x)
        for (std::size_t j = 0; j < nq; ++j) g.gradient[x] += si.B[j][x] * gq[j];
    g.description = format("numerical (central differences of %s energies along %zu %sinternal coordinate%s, %d "
                           "energies, step %.4f bohr / rad)",
                           method.c_str(), disps.size(), symmetric ? "symmetry-adapted " : "",
                           disps.size() == 1 ? "" : "s", g.energy_evaluations, h);
    g.seconds = seconds_since(t0);
    return g;
}

}  // namespace macqc
