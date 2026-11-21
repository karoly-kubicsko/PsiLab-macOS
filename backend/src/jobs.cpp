#include "macqc/jobs.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <sstream>

#include "macqc/cc.h"
#include "macqc/constants.h"
#include "macqc/gto_derivatives.h"
#include "macqc/job.h"
#include "macqc/output.h"

namespace macqc {

namespace {

using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string human_time(double s) { return estimate_duration(s); }

std::string method_label(const std::string& m) { return m == "RHF" ? "HF" : (m == "FCI" ? "Full CI" : m); }
// 15.0: with the reference: HF (closed shell), UHF, ROHF; UMP2, ROCCSD(T); Full CI
std::string job_label(const std::string& m, const std::string& ref) {
    if (m == "FCI") return "Full CI";
    if (ref.empty() || ref == "RHF") return method_label(m);
    return macqc::method_label(m, ref);
}

ElectronicOptions electronic_options(const MolecularResult& p) {
    const MolecularInput& in = p.input;
    ElectronicOptions o;
    o.reference = p.reference;
    o.multiplicity = in.multiplicity;
    o.threads = in.threads;
    o.memory_mb = in.memory_mb;
    o.fci_memory_mb = in.fci_memory_mb;
    o.scf_max_iter = in.max_iter;
    o.scf_energy_tol = std::min(in.energy_tol, 1.0e-10);
    o.scf_gradient_tol = std::min(in.gradient_tol, 1.0e-9);
    return o;
}

// Cost of one gradient (seconds); sets *how. Calibrated against measured times (4 threads,
// x86-64): an RHF gradient costs about 4 RHF energies (3.2-4.0 measured for H2O aug-cc-pVDZ to
// aug-cc-pVQZ); correlated analytic gradients 1.1-2 times the sum of their parts, so that sum
// is scaled by 1.3 plus 0.2 s of set-up; numerical gradients (n + 1) energies, each with a fixed
// set-up time (about 0.25 s for the general-order CCSDT/CCSDTQ engine).
double gradient_seconds(const MolecularResult& p, bool* analytic, std::string* how, int* energies) {
    const std::string& m = p.input.method;
    const CostEstimate& e = p.estimate;
    const int nbf = p.basis.n_functions, nocc = p.molecule.electrons() / 2;
    const bool c2v_or_diatomic = true;  // the sidebar geometries are C2v (H2O2: C2)
    const std::size_t nat = p.molecule.atoms.size();
    const int numerical_energies =
        uses_cartesian_displacements(p.molecule)
            ? static_cast<int>(2 * vibrational_directions(p.molecule).size())  // 9.0: custom molecules
            : (nat == 2 ? 2 : (nat == 4 ? (p.molecule.name == "H2O2" ? 8 : 12) : (c2v_or_diatomic ? 4 : 6)));
    *energies = 0;
    if (e.so && analytic_gradient_method(m, p.reference)) {
        // 15.0: the spin-orbital engine: amplitudes, multipliers (about as many iterations), densities, Z-vector
        *analytic = true;
        *how = "analytic gradient (" + e.so->method + ": Lambda equations of the automatically differentiated "
               "Lagrangian, relaxed densities, Z-vector)";
        const double n5 = 6.0 * std::pow(static_cast<double>(nbf), 5) / (4.0e9 * std::max(1, e.threads));
        return 2.5 * e.total_seconds + 3.0 * e.rhf_seconds + n5 + 0.2;
    }
    if (analytic_gradient_method(m, p.reference)) {
        const double need = analytic_gradient_bytes(m, nocc, nbf, nbf) / (1024.0 * 1024.0);
        if (need <= p.input.fci_memory_mb) {
            *analytic = true;
            if (m == "RHF") {
                *how = "analytic gradient (derivative integrals)";
                return 4.0 * e.rhf_seconds;
            }
            const double n = nbf;
            const double n5 = 6.0 * std::pow(n, 5) / (4.0e9 * std::max(1, e.threads));
            if (m == "FCI") {
                *how = p.reference.empty() || p.reference == "RHF"
                           ? "analytic gradient (Full CI densities, Z-vector)"
                           : "analytic gradient (Full CI spin densities; orbital-invariant energy, no Z-vector)";
                return 1.3 * (e.total_seconds * 1.3 + 2.5 * e.rhf_seconds + n5) + 0.2;
            }
            *how = "analytic gradient (" + m + " Lambda equations, relaxed densities, Z-vector)";
            double t = e.total_seconds + 2.5 * e.rhf_seconds + n5;
            if (e.cc) t += e.cc->iterations * e.cc->iteration_seconds + 6.0 * e.cc->triples_seconds;
            return 1.3 * t + 0.2;
        }
        *how = strprintf("numerical gradient (%d %s energies; the analytic gradient would need %s)",
                         numerical_energies, m.c_str(), estimate_memory_mb(need).c_str());
    } else {
        *how = strprintf("numerical gradient (%d %s energies per gradient, central differences)", numerical_energies,
                         m.c_str());
    }
    *analytic = false;
    *energies = numerical_energies;
    const double setup = (m == "CCSDT" || m == "CCSDTQ") ? 0.25 : 0.05;
    return (numerical_energies + 1) * (e.total_seconds + setup);
}

}  // namespace

std::string canonical_job_type(const std::string& type) {
    std::string t;
    for (char c : type)
        if (!std::isspace(static_cast<unsigned char>(c)) && c != '_' && c != '-')
            t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "energy" || t == "sp" || t == "singlepoint" || t == "1") return "energy";
    if (t == "opt" || t == "optimization" || t == "optimisation" || t == "geometryoptimization" ||
        t == "geometryoptimisation" || t == "2")
        return "opt";
    if (t == "freq" || t == "frequencies" || t == "frequency" || t == "3") return "freq";
    if (t == "opt+freq" || t == "optfreq" || t == "4") return "opt+freq";
    if (t == "opt+freq+sp" || t == "optfreqsp" || t == "5") return "opt+freq+sp";
    if (t == "scan" || t == "relaxedscan" || t == "pesscan" || t == "relaxedpesscan" || t == "modredundant" ||
        t == "optmodredundant" || t == "6")
        return "scan";
    if (t == "ts" || t == "opt=ts" || t == "optts" || t == "transitionstate" || t == "tsopt" || t == "7") return "ts";
    if (t == "ts+freq" || t == "tsfreq" || t == "opt=tsfreq" || t == "optts+freq" || t == "8") return "ts+freq";
    if (t == "ts+freq+sp" || t == "tsfreqsp" || t == "opt=tsfreqsp" || t == "optts+freq+sp" || t == "9") return "ts+freq+sp";
    return "";
}

std::string job_title(const std::string& t) {
    if (t == "opt") return "Geometry optimisation";
    if (t == "freq") return "Frequencies";
    if (t == "opt+freq") return "Optimisation + frequencies";
    if (t == "opt+freq+sp") return "Optimisation + frequencies + single-point energy";
    if (t == "scan") return "Relaxed potential energy surface scan";
    if (t == "ts") return "Transition-state optimisation";
    if (t == "ts+freq") return "Transition-state optimisation + frequencies";
    if (t == "ts+freq+sp") return "Transition-state optimisation + frequencies + single-point energy";
    return "Energy";
}

bool job_has_optimization(const std::string& t) {
    return t == "opt" || t == "opt+freq" || t == "opt+freq+sp" || t == "ts" || t == "ts+freq" || t == "ts+freq+sp";
}
bool job_has_frequencies(const std::string& t) {
    return t == "freq" || t == "opt+freq" || t == "opt+freq+sp" || t == "ts+freq" || t == "ts+freq+sp";
}
bool job_has_ts(const std::string& t) { return t == "ts" || t == "ts+freq" || t == "ts+freq+sp"; }
bool job_has_single_point(const std::string& t) { return t == "opt+freq+sp" || t == "ts+freq+sp"; }
bool job_has_scan(const std::string& t) { return t == "scan"; }
bool job_uses_optimizer(const std::string& t) { return job_has_optimization(t) || job_has_scan(t); }  // incl. ts

std::string parse_scan_coordinate(const std::string& spec, int n_atoms, ScanCoordinate* out) {
    std::string text;
    for (char c : spec) text += (c == ',' || c == ';' || c == '(' || c == ')') ? ' ' : c;
    std::vector<std::string> tok;
    {
        std::string cur;
        for (char c : text) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!cur.empty()) tok.push_back(cur);
                cur.clear();
            } else {
                cur += c;
            }
        }
        if (!cur.empty()) tok.push_back(cur);
    }
    if (tok.empty()) return "the scan coordinate is empty (e.g. 1 2, 2 1 3 or B 1 2)";
    std::string kind;
    for (char c : tok[0]) kind += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (std::isdigit(static_cast<unsigned char>(tok[0][0]))) {
        // no type letter: the number of atoms decides, as in Gaussian
        const std::size_t n = tok.size();
        if (n < 2 || n > 4)
            return strprintf("a scan coordinate has 2 (bond), 3 (angle) or 4 (dihedral) atom numbers, got %zu", n);
        tok.insert(tok.begin(), n == 2 ? "B" : (n == 3 ? "A" : "D"));
        kind = n == 2 ? "b" : (n == 3 ? "a" : "d");
    }
    ScanCoordinate sc;
    std::size_t need = 0;
    if (kind == "b" || kind == "r" || kind == "bond" || kind == "stretch" || kind == "distance") {
        sc.coordinate.type = InternalType::Stretch;
        sc.letter = "B";
        sc.kind = "bond";
        need = 2;
    } else if (kind == "a" || kind == "angle" || kind == "bend") {
        sc.coordinate.type = InternalType::Bend;
        sc.letter = "A";
        sc.kind = "angle";
        need = 3;
    } else if (kind == "d" || kind == "dihedral" || kind == "torsion") {
        sc.coordinate.type = InternalType::Dihedral;
        sc.letter = "D";
        sc.kind = "dihedral";
        need = 4;
    } else {
        return "unknown scan coordinate type '" + tok[0] + "' (B = bond, A = angle, D = dihedral)";
    }
    if (tok.size() - 1 != need)
        return strprintf("a %s needs %zu atom numbers, got %zu (e.g. %s)", sc.kind.c_str(), need, tok.size() - 1,
                         need == 2 ? "B 1 2" : (need == 3 ? "A 2 1 3" : "D 3 1 2 4"));
    for (std::size_t k = 1; k < tok.size(); ++k) {
        char* end = nullptr;
        const long v = std::strtol(tok[k].c_str(), &end, 10);
        if (!end || *end != '\0') return "atom numbers must be integers, got '" + tok[k] + "'";
        if (v < 1 || v > n_atoms) return strprintf("atom %ld does not exist (the molecule has %d atoms)", v, n_atoms);
        const int a = static_cast<int>(v - 1);
        if (std::find(sc.coordinate.atoms.begin(), sc.coordinate.atoms.end(), a) != sc.coordinate.atoms.end())
            return strprintf("atom %ld appears twice in the scan coordinate", v);
        sc.coordinate.atoms.push_back(a);
    }
    if (out) *out = sc;
    return "";
}

std::string parse_scan_line(const std::string& line, int n_atoms, JobSettings* settings) {
    std::istringstream in(line);
    std::vector<std::string> tok;
    for (std::string t; in >> t;) tok.push_back(t);
    std::size_t sp = tok.size();
    for (std::size_t i = 0; i < tok.size(); ++i)
        if (tok[i] == "S" || tok[i] == "s") sp = i;
    if (sp == tok.size() || tok.size() != sp + 3)
        return "expected '[B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize' (Gaussian ModRedundant), e.g. '1 2 S 10 0.05' or "
               "'3 1 2 4 S 18 10.0'";
    std::string coord;
    for (std::size_t i = 0; i < sp; ++i) coord += (coord.empty() ? "" : " ") + tok[i];
    ScanCoordinate sc;
    const std::string e = parse_scan_coordinate(coord, n_atoms, &sc);
    if (!e.empty()) return e;
    char* end = nullptr;
    const double n = std::strtod(tok[sp + 1].c_str(), &end);
    if (!end || *end != '\0' || n != std::floor(n) || n < 1 || n > 360)
        return "the number of scan steps must be an integer from 1 to 360, got '" + tok[sp + 1] + "'";
    const double step = std::strtod(tok[sp + 2].c_str(), &end);
    if (!end || *end != '\0' || !std::isfinite(step) || step == 0.0)
        return "the scan step size must be a non-zero number, got '" + tok[sp + 2] + "'";
    if (settings) {
        settings->scan_coordinate = coord;
        settings->scan_steps = static_cast<int>(n);
        settings->scan_step = sc.letter == "B" ? step / kBohrInAngstrom : step;
    }
    return "";
}

std::string plan_scan(const Molecule& mol, const JobSettings& s, ScanPlan* plan) {
    ScanPlan pl;
    const std::string e = parse_scan_coordinate(s.scan_coordinate, static_cast<int>(mol.atoms.size()), &pl.coordinate);
    if (!e.empty()) return e;
    pl.coordinate.definition = pl.coordinate.coordinate.definition(mol);
    if (s.scan_steps < 1 || s.scan_steps > 360) return "the number of scan steps must be between 1 and 360";
    if (!std::isfinite(s.scan_step) || s.scan_step == 0.0) return "the scan step must be a non-zero number";
    const double deg = 180.0 / 3.14159265358979323846;
    const InternalCoordinate& c = pl.coordinate.coordinate;
    if (c.type == InternalType::Dihedral) {
        InternalCoordinate b1, b2;
        b1.type = b2.type = InternalType::Bend;
        b1.atoms = {c.atoms[0], c.atoms[1], c.atoms[2]};
        b2.atoms = {c.atoms[1], c.atoms[2], c.atoms[3]};
        const double a1 = internal_value(mol, b1) * deg, a2 = internal_value(mol, b2) * deg;
        if (a1 > 175.0 || a2 > 175.0 || a1 < 5.0 || a2 < 5.0)
            return strprintf("the dihedral angle %s is not defined at this geometry (angles %.1f and %.1f degrees; both "
                             "must be between 5 and 175)",
                             pl.coordinate.definition.c_str(), a1, a2);
    }
    const double v0 = internal_value(mol, c) * (c.type == InternalType::Stretch ? 1.0 : deg);
    pl.input_value = v0;
    pl.start = s.scan_start ? *s.scan_start : v0;
    pl.step = s.scan_step;
    pl.steps = s.scan_steps;
    if (!std::isfinite(pl.start)) return "the scan start must be a finite number";
    if (c.type == InternalType::Stretch && std::fabs(pl.step) > 2.0) return "the scan step of a bond must be at most 2 bohr";
    if (c.type != InternalType::Stretch && std::fabs(pl.step) > 90.0)
        return "the scan step of an angle must be at most 90 degrees";
    for (int k = 0; k <= pl.steps; ++k) {
        const double v = pl.start + k * pl.step;
        if (c.type == InternalType::Stretch) {
            if (v < 0.3 - 1e-12 || v > 20.0 + 1e-12)
                return strprintf("the scan reaches %s = %.4f bohr (%.4f Å); bonds must stay between 0.3 and 20 bohr",
                                 pl.coordinate.definition.c_str(), v, v * kBohrInAngstrom);
            pl.values.push_back(v);
            pl.targets.push_back(v);
        } else if (c.type == InternalType::Bend) {
            if (v < 10.0 - 1e-9 || v > 180.0 + 1e-9)
                return strprintf("the scan reaches %s = %.4f degrees; angles must stay between 10 and 180 degrees",
                                 pl.coordinate.definition.c_str(), v);
            const double vv = std::min(v, 180.0);
            pl.values.push_back(vv);
            pl.targets.push_back(std::fabs(vv - 180.0) < 1e-9 ? 3.14159265358979323846 : vv / deg);
        } else {
            pl.values.push_back(v);
            pl.targets.push_back(wrap_angle(v / deg));
        }
    }
    if (plan) *plan = pl;
    return "";
}

static bool fewer_than_two_atoms(const MolecularResult& p) { return p.molecule.atoms.size() < 2; }

JobEstimate estimate_job(const MolecularResult& p, const JobSettings& s, const MolecularResult* sp) {
    JobEstimate est;
    const std::string& m = p.input.method;
    const std::string& t = s.type;
    const CostEstimate& e = p.estimate;
    if (!p.refusal.empty()) est.refusal = p.refusal;
    const int nbf = p.basis.n_functions, nocc = p.molecule.electrons() / 2;
    // 15.0: alpha and beta electrons (open-shell Hessians)
    int n_alpha = nocc, n_beta = nocc;
    try {
        const auto ab = alpha_beta_electrons(p.molecule.electrons(), p.input.multiplicity);
        n_alpha = ab.first;
        n_beta = ab.second;
    } catch (const std::exception&) {
    }
    const double hessian_mb = analytic_hessian_bytes(m, p.reference, n_alpha, n_beta, nbf, nbf,
                                                     static_cast<int>(p.molecule.atoms.size())) / 1048576.0;
    if (t == "energy") {
        JobStepEstimate st{"Energy", m, p.input.basis, "single point", true, 1, "energies", e.total_seconds, 0.0};
        est.steps.push_back(st);
        est.total_seconds = e.total_seconds;
        return est;
    }
    bool analytic = false;
    std::string how;
    int grad_energies = 0;
    const double tg = gradient_seconds(p, &analytic, &how, &grad_energies);
    const double grad_mb = analytic && m != "RHF" ? analytic_gradient_bytes(m, nocc, nbf, nbf) / (1024.0 * 1024.0) : 0.0;
    est.memory_mb = grad_mb;
    if (e.fci) est.memory_mb = std::max(est.memory_mb, e.fci->memory_bytes / (1024.0 * 1024.0));
    if (e.cc) est.memory_mb = std::max(est.memory_mb, e.cc->memory_bytes / (1024.0 * 1024.0));
    const int typical_steps = 6;
    // 10.0: the analytic Full CI, CCSD or CCSD(T) Hessian at the first point of a minimisation
    // ("auto" or "analytic", where it fits into memory; no dipole derivatives): extra seconds
    const bool correlated_hessian_fits =
        analytic && m != "RHF" && analytic_hessian_method(m, p.reference) && p.molecule.atoms.size() >= 2 &&
        hessian_mb <= p.input.fci_memory_mb;
    const bool analytic_initial = correlated_hessian_fits && s.initial_hessian != "lindh";
    const double initial_hessian_seconds =
        analytic_initial
            ? ((m == "FCI" ? 0.6 : 1.5) * static_cast<double>(vibrational_directions(p.molecule).size()) + 3.0) * tg
            : 0.0;
    auto hessian_estimate = [&]() {
        JobStepEstimate st;
        st.name = "Hessian";
        st.method = m;
        st.basis = p.input.basis;
        const bool diatomic = p.molecule.atoms.size() == 2;
        if (m == "RHF") {
            st.analytic = true;
            st.how = p.reference == "RHF" ? "analytic Hessian (second-derivative integrals, CPHF), analytic dipole derivatives"
                                          : "analytic " + p.reference + " Hessian (second-derivative integrals, coupled-perturbed " +
                                                p.reference + "), analytic dipole derivatives";
            st.evaluations = 1;
            st.evaluation_kind = "Hessians";
            st.seconds = 7.0 * tg;
        } else if (analytic && m != "RHF" && analytic_hessian_method(m, p.reference) && s.hessian != "numerical" &&
                   !fewer_than_two_atoms(p) &&
                   hessian_mb <= p.input.fci_memory_mb) {
            // 10.0: analytic correlated Hessian, responses along the internal directions and the
            // field components (dipole derivatives)
            st.analytic = true;
            const int dirs = static_cast<int>(vibrational_directions(p.molecule).size());
            const int responses = dirs + 3;
            st.evaluations = responses;
            st.evaluation_kind = "responses";
            const bool open = !(p.reference.empty() || p.reference == "RHF");
            if (m == "FCI" && open)
                st.how = strprintf("analytic Hessian: CI-vector responses along the %d internal direction%s and the 3 "
                                   "electric-field components (IR intensities); orbital-invariant energy, no CPHF or Z-vector",
                                   dirs, dirs == 1 ? "" : "s");
            else
                st.how = strprintf("analytic Hessian: %s responses along the %d internal direction%s and the 3 electric-field "
                                   "components (IR intensities), %s orbital response of all %zu coordinates, Z-vector",
                                   m == "FCI" ? "CI-vector" : "amplitude", dirs, dirs == 1 ? "" : "s",
                                   open ? ("coupled-perturbed " + p.reference).c_str() : "CPHF", 3 * p.molecule.atoms.size());
            st.seconds = ((m == "FCI" ? 0.6 : 1.5) * responses + 3.0) * tg;
            st.memory_mb = hessian_mb;
        } else if (analytic && uses_cartesian_displacements(p.molecule)) {
            // 9.0: custom molecules, Cartesian directions orthogonal to translations and rotations
            st.analytic = false;
            st.evaluations = static_cast<int>(2 * vibrational_directions(p.molecule).size());
            st.evaluation_kind = "gradients";
            st.how = strprintf("numerical Hessian: central differences of %d analytic gradients along the %d Cartesian "
                               "directions orthogonal to the translations and rotations; IR intensities from the "
                               "relaxed dipole moments",
                               st.evaluations, st.evaluations / 2);
            st.seconds = (st.evaluations + 1) * tg;
        } else if (analytic) {
            st.analytic = false;
            const bool four = p.molecule.atoms.size() == 4;
            st.evaluations = diatomic ? 2 : (four ? 12 : 6);
            st.evaluation_kind = "gradients";
            st.how = strprintf("numerical Hessian: central differences of %d analytic gradients along %s; IR intensities "
                               "from the relaxed dipole moments",
                               st.evaluations,
                               diatomic ? "the bond"
                                        : (four ? (p.molecule.name == "H2O2" ? "the symmetry-adapted coordinates S1 ... S6"
                                                                              : "the six internal coordinates")
                                                : "the symmetry-adapted coordinates S1, S2, S3"));
            st.seconds = (st.evaluations + 1) * tg;
        } else {
            st.analytic = false;
            const NumericalHessianCost c = numerical_hessian_cost(p.molecule, false, true);
            st.evaluations = c.energies + c.field_energies;
            st.evaluation_kind = "energies";
            st.how = c.field_energies > 0
                         ? strprintf("numerical Hessian: %d %s energies (second differences) + %d finite-field energies for "
                                     "the dipole derivatives",
                                     c.energies, m.c_str(), c.field_energies)
                         : strprintf("numerical Hessian: %d %s energies (second differences); no IR intensities needed "
                                     "(no dipole moment change by symmetry)",
                                     c.energies, m.c_str());
            st.seconds = (st.evaluations + 1) * (e.total_seconds + ((m == "CCSDT" || m == "CCSDTQ") ? 0.25 : 0.05));
        }
        return st;
    };
    if (job_has_optimization(t)) {
        JobStepEstimate st;
        st.name = job_has_ts(t) ? "Transition-state optimisation" : "Geometry optimisation";
        st.method = m;
        st.basis = p.input.basis;
        st.analytic = analytic;
        st.how = how + strprintf("; typically %d steps (at most %d)", typical_steps, s.max_steps);
        st.evaluations = typical_steps;
        st.evaluation_kind = "steps";
        st.seconds = typical_steps * tg;
        const std::string hmode = job_has_ts(t) ? s.ts_hessian : s.opt_hessian;
        if (!hmode.empty()) {
            // the Hessian computed at the first point (CalcFC), at every point (CalcAll) or at the
            // first and every N-th point (15.0: RecalcFC=N)
            const JobStepEstimate h = hessian_estimate();
            const bool all = hmode == "calcall";
            const int computed = all ? typical_steps
                                     : hmode == "recalcfc" ? 1 + (typical_steps - 1) / std::max(1, s.recalc_fc) : 1;
            std::string hh = h.how;
            const auto ir = hh.find("; IR");
            if (ir != std::string::npos) hh = hh.substr(0, ir);
            const auto dip = hh.find(", analytic dipole");
            if (dip != std::string::npos) hh = hh.substr(0, dip);
            st.how += std::string(job_has_ts(t) ? "; P-RFO steps, Bofill updates" : "; RFO steps, BFGS updates") +
                      "; Hessian computed " +
                      (all ? std::string("at every step (CalcAll)")
                       : hmode == "recalcfc" ? strprintf("at step 1 and every %d steps (RecalcFC=%d)", s.recalc_fc, s.recalc_fc)
                                             : std::string("at the first step (CalcFC)")) +
                      ": " + hh;
            st.seconds += computed * h.seconds;
            st.analytic = st.analytic && h.analytic;
        } else if (m == "RHF" && nbf <= 100) {
            st.seconds += 6.0 * tg;  // analytic initial Hessian
        } else if (analytic_initial) {
            st.seconds += initial_hessian_seconds;  // 10.0: analytic correlated initial Hessian
            st.how += "; analytic initial Hessian";
        }
        st.memory_mb = grad_mb;
        est.steps.push_back(st);
    }
    if (job_has_scan(t)) {
        JobStepEstimate st;
        st.name = "Relaxed scan";
        st.method = m;
        st.basis = p.input.basis;
        st.analytic = analytic;
        const int points = std::max(1, s.scan_steps + 1);
        if (p.molecule.atoms.size() == 2) {
            st.how = strprintf("%d single-point energies (a diatomic has no other coordinate to relax)", points);
            st.evaluations = points;
            st.evaluation_kind = "energies";
            st.seconds = points * e.total_seconds;
            st.analytic = true;
        } else {
            // the first point from the input Hessian guess, later points start from the previous
            // point's geometry and Hessian
            const int first = typical_steps, later = 4;
            const int total = first + later * (points - 1);
            st.how = how + strprintf("; %d points, typically %d optimisation steps at the first and %d at each further "
                                     "point (at most %d per point)",
                                     points, first, later, s.max_steps);
            st.evaluations = total;
            st.evaluation_kind = "steps";
            st.seconds = total * tg;
            if (m == "RHF" && nbf <= 100) st.seconds += 6.0 * tg;  // analytic initial Hessian
            if (analytic_initial) {
                st.seconds += initial_hessian_seconds;  // 10.0: at the first point
                st.how += "; analytic initial Hessian at the first point";
            }
        }
        st.memory_mb = grad_mb;
        est.steps.push_back(st);
    }
    if (job_has_frequencies(t)) est.steps.push_back(hessian_estimate());
    if (job_has_single_point(t) && sp) {
        JobStepEstimate st;
        st.name = "Single-point energy";
        st.method = sp->input.method;
        st.basis = sp->input.basis;
        st.how = job_has_ts(t) ? "single point at the transition state" : "single point at the optimised geometry";
        st.evaluations = 1;
        st.evaluation_kind = "energies";
        st.seconds = sp->estimate.total_seconds;
        if (sp->estimate.fci) st.memory_mb = sp->estimate.fci->memory_bytes / (1024.0 * 1024.0);
        if (sp->estimate.cc) st.memory_mb = std::max(st.memory_mb, sp->estimate.cc->memory_bytes / (1024.0 * 1024.0));
        est.memory_mb = std::max(est.memory_mb, st.memory_mb);
        est.steps.push_back(st);
        if (!sp->refusal.empty() && est.refusal.empty())
            est.refusal = "The single-point step is not possible: " + sp->refusal;
    }
    for (const auto& st : est.steps) est.total_seconds += st.seconds;
    // refusals
    if (est.refusal.empty() && job_has_scan(t)) {
        ScanPlan plan;
        const std::string why = plan_scan(p.molecule, s, &plan);
        if (!why.empty()) est.refusal = "The scan is not possible: " + why + ".";
    }
    if (est.refusal.empty() && p.molecule.atoms.size() == 1)
        est.refusal = "A single atom has no geometry to optimise, scan or vibrate: only the energy job is possible.";
    if (est.refusal.empty() && (m == "CCSDT" || m == "CCSDTQ") && uses_cartesian_displacements(p.molecule) &&
        (job_has_frequencies(t) || job_has_ts(t))) {
        const NumericalHessianCost c = numerical_hessian_cost(p.molecule, false, true);
        est.refusal = strprintf(
            "%s %s for a custom molecule of %zu atoms are not possible: %s has no analytic gradient, and a Hessian "
            "from energies is not available (it would need about %d %s energies of %s each). Use CCSD(T) for the "
            "Hessian, and %s in the single-point step of job 5.",
            m.c_str(), job_has_ts(t) ? "transition-state searches" : "frequencies", p.molecule.atoms.size(), m.c_str(),
            c.energies, m.c_str(), human_time(e.total_seconds).c_str(), m.c_str());
    }
    if (est.refusal.empty() && job_has_ts(t) && p.molecule.atoms.size() == 2)
        est.refusal = "A diatomic molecule has no transition state: its energy depends on the bond length only. Use a "
                      "relaxed scan (job scan) for the potential curve.";
    if (est.refusal.empty() && (m == "CCSDT" || m == "CCSDTQ") && p.molecule.atoms.size() == 4 &&
        (job_has_frequencies(t) || job_has_ts(t))) {
        est.refusal = strprintf(
            "%s %s for %s are not possible: %s has no analytic gradient, and a Hessian from energies is "
            "not available for four atoms (it would need about 140 %s energies of %s each). Use CCSD(T) for the "
            "Hessian.",
            m.c_str(), job_has_ts(t) ? "transition-state searches" : "frequencies", p.molecule.name.c_str(), m.c_str(),
            m.c_str(), human_time(e.total_seconds).c_str());
    }
    if (est.refusal.empty() && m == "CCSDTQ" && p.molecule.name == "H2O") {
        const NumericalHessianCost c = numerical_hessian_cost(p.molecule, false, true);
        est.refusal = strprintf(
            "CCSDTQ %s for H2O is not possible: CCSDTQ has no analytic gradient, so %s would need %s CCSDTQ energies "
            "of %s each (%.1f GB of memory). Use CCSD(T) or CCSDT for the geometry and the frequencies, and CCSDTQ in "
            "the single-point step of job 5.",
            job_has_frequencies(t) ? "frequencies" : (job_has_scan(t) ? "scans" : "optimisation"),
            job_has_frequencies(t) ? "the Hessian" : "every optimisation step",
            job_has_frequencies(t) ? std::to_string(c.energies + c.field_energies).c_str() : "4",
            human_time(e.total_seconds).c_str(), e.cc ? e.cc->memory_bytes / 1073741824.0 : 0.0);
    }
    if (est.refusal.empty() && est.total_seconds > kMaxCCSeconds)
        est.refusal = strprintf("The %s would take %s (limit: three days). Use a smaller basis set or a cheaper method.",
                                job_title(t).c_str(), human_time(est.total_seconds).c_str());
    if (est.refusal.empty() && est.total_seconds > kLongCalculationSeconds)
        est.warnings.push_back(strprintf("Long calculation: %s with %s/%s, estimated %s on %d threads.", job_title(t).c_str(),
                                         method_label(m).c_str(), p.input.basis.c_str(), human_time(est.total_seconds).c_str(),
                                         e.threads));
    for (const auto& st : est.steps)
        est.notes.push_back(st.name + ": " + st.how + (st.seconds > 0 ? " — " + human_time(st.seconds) : ""));
    return est;
}

namespace {

std::string scan_value_text(const ScanPlan& pl, double v) {
    return pl.angular() ? strprintf("%.6f deg", v) : strprintf("%.8f Å (%.8f bohr)", v * kBohrInAngstrom, v);
}

// Text shared by both job paths.
void write_header(OutputWriter& w, const JobSettings& s, const std::string& echo, int threads) {
    const std::string rule(88, '=');
    w.line(rule);
    w.line("                                ΨLab " + std::string(kProgramVersion));
    w.line("  Hartree–Fock (RHF, UHF, ROHF), Møller–Plesset, CI, QCI, coupled cluster and Full CI");
    w.line(rule);
    w.line("  Program      ΨLab " + std::string(kProgramVersion) + " (C++17 backend; integrals, SCF, CC, CI, derivatives)");
    w.line("  Started      " + local_time_string());
    w.line("  Job          " + job_title(s.type) + " (" + s.type + ")");
    w.line(strprintf("  Threads      %d", threads));
    w.line(std::string("  Matrices     ") + gemm_backend_name());
    w.section("Input");
    w.write(echo);
    if (!echo.empty() && echo.back() != '\n') w.line();
}

void write_molecule(OutputWriter& w, const MolecularResult& p, const std::string& what) {
    const Molecule& m = p.molecule;
    w.section("Molecule");
    const int mult = p.input.multiplicity;
    const char* spin_names[] = {"singlet", "doublet", "triplet", "quartet", "quintet", "sextet", "septet"};
    const std::string spin = mult == 1 ? "closed-shell singlet"
                                       : (mult <= 7 ? spin_names[mult - 1] : strprintf("multiplicity %d", mult)) +
                                             std::string(", reference ") + p.reference;
    w.line(strprintf("    %s, %d electrons%s, charge %d, %s", m.display_name().c_str(), m.electrons(),
                     m.core_electrons() ? strprintf(" (and %d core electrons in ECPs)", m.core_electrons()).c_str() : "",
                     m.charge, spin.c_str()));
    if (m.custom) w.line(strprintf("    %s (%zu atoms, Cartesian coordinates)", m.title.empty() ? "custom molecule" : m.title.c_str(), m.atoms.size()));
    w.line("    " + what);
    w.write(geometry_table(m));
    w.write(parameter_lines(m));
    w.line(strprintf("    Nuclear repulsion energy %20.10f E_h", m.nuclear_repulsion()));
    w.section("Basis set");
    w.write(basis_block(p.basis, p.notes));
    const CostEstimate& e = p.estimate;
    w.line(strprintf("    two-electron integrals: %zu unique, %s", e.unique_eris,
                     e.eri_stored ? strprintf("kept in memory (%.1f MB)", e.eri_memory_mb).c_str()
                                  : strprintf("integral-direct (storing them would need %s)", estimate_memory_mb(e.eri_memory_mb).c_str()).c_str()));
}

void write_estimate(OutputWriter& w, const MolecularResult& p, const JobEstimate& est) {
    w.section("Size, memory and time (estimates)");
    const CostEstimate& e = p.estimate;
    auto count = [](double n) {
        if (!std::isfinite(n) || n > 1.0e300) return std::string("more than 10^300");
        return n < 1.0e9 ? strprintf("%.0f", n) : strprintf("%.3g", n);
    };
    auto bytes = [](double b) { return estimate_memory_mb(b / 1048576.0); };
    if (e.fci)
        w.line(strprintf("    Full CI: %s determinants (%d electrons in %d orbitals), about %s", count(e.fci->determinants).c_str(),
                         e.fci->alpha + e.fci->beta, e.fci->orbitals, bytes(e.fci->memory_bytes).c_str()));
    if (e.cc)
        w.line(strprintf("    %s: %d occupied, %d virtual orbitals, %s amplitudes, about %s", e.cc->method.c_str(),
                         e.cc->occupied, e.cc->virtuals, count(e.cc->amplitudes).c_str(), bytes(e.cc->memory_bytes).c_str()));
    if (e.so)
        w.line(strprintf("    %s: %s amplitudes in dense spin blocks (spin-orbital engine), about %s", e.so->method.c_str(),
                         count(e.so->amplitudes).c_str(), bytes(e.so->memory_bytes).c_str()));
    for (const auto& n : est.notes) w.line("    " + n);
    if (est.memory_mb > 0) w.line(strprintf("    memory of the derivative steps: about %s", estimate_memory_mb(est.memory_mb).c_str()));
    w.line("    total: " + human_time(est.total_seconds) + strprintf(" on %d threads", e.threads));
    const bool job_level = !(est.steps.size() == 1 && est.steps[0].name == "Energy");
    for (const auto& x : p.warnings)  // the job-level estimate replaces the single-energy "long calculation" warning
        if (!(job_level && x.rfind("Long calculation", 0) == 0)) w.line("    Warning: " + x);
    for (const auto& x : est.warnings) w.line("    Warning: " + x);
    if (!est.refusal.empty()) w.line("    NOT POSSIBLE: " + est.refusal);
}

void write_correlated(OutputWriter& w, const SinglePoint& sp) {
    if (sp.so) w.write(so_report_section(*sp.so, 1));
    if (sp.cc) w.write(cc_report_section(*sp.cc, 1));
    if (sp.fci) w.write(fci_report_section(*sp.fci, sp.rhf.total_energy, 1));
}

void write_termination(OutputWriter& w, JobResult& r, double seconds, const AppendixWriter& appendix) {
    // 13.0: the basis set in full, the molecular orbitals and the archive of the results
    if (appendix && r.termination != "cancelled") {
        std::string text;
        try {
            text = appendix(r);
        } catch (const std::exception& e) {
            text = std::string("\n    (the appendix could not be written: ") + e.what() + ")\n";
        }
        w.write(text);
    }
    w.line();
    const std::string rule(88, '=');
    w.line(rule);
    if (r.termination == "normal")
        w.line("  Normal termination of ΨLab " + std::string(kProgramVersion) + " at " + local_time_string() +
               "; total time " + duration_string(seconds) + ".");
    else if (r.termination == "cancelled")
        w.line("  Calculation cancelled at " + local_time_string() + " after " + duration_string(seconds) + ".");
    else
        w.line("  Error termination of ΨLab " + std::string(kProgramVersion) + " at " + local_time_string() + ": " + r.message);
    w.line(rule);
}

struct Streams {
    OutputWriter* w = nullptr;
    std::string context;  // progress-text prefix, e.g. "Optimisation step 3 · "
    std::string method;
    std::string reference = "RHF";
    bool scf_started = false, cc_started = false, fci_started = false, lambda_started = false, cphf_started = false;
    double scf_energy = 0.0;
    int scf_iterations = 0;
    void reset() { scf_started = cc_started = fci_started = lambda_started = cphf_started = false; }
    void scf_done_line() {
        if (scf_started)
            w->line(strprintf("    SCF converged: E(%s) = %.10f E_h after %d iterations", reference.c_str(), scf_energy, scf_iterations));
    }
};

ElectronicCallbacks stream_callbacks(Streams& s) {
    ElectronicCallbacks cb;
    cb.scf = [&s](const RHFIteration& it) {
        if (!s.scf_started) {
            s.w->line(s.reference == "UHF"    ? "\n    SCF (unrestricted Hartree–Fock, DIIS)"
                      : s.reference == "ROHF" ? "\n    SCF (restricted open-shell Hartree–Fock, Roothaan effective Fock matrix, DIIS)"
                                              : "\n    SCF (restricted Hartree–Fock, DIIS)");
            s.w->write(scf_header());
            s.scf_started = true;
        }
        s.w->write(scf_iteration_line(it));
        s.scf_energy = it.energy;
        s.scf_iterations = it.iteration;
        report_progress_text(s.context + strprintf("SCF iteration %d · E = %.10f E_h", it.iteration, it.energy));
    };
    cb.cc = [&s](const CCIteration& it) {
        if (!s.cc_started) {
            s.scf_done_line();
            s.w->line();
            s.w->write(cc_header(s.method == "CCSD(T)" && s.reference == "RHF" ? "CCSD" : job_label(s.method, s.reference)));
            s.cc_started = true;
        }
        s.w->write(cc_iteration_line(it));
        report_progress_text(s.context + strprintf("%s iteration %d · E = %.10f E_h · max |dt| %.1e",
                                                   s.method == "CCSD(T)" ? "CCSD" : s.method.c_str(), it.iteration,
                                                   it.energy, it.residual));
    };
    cb.fci = [&s](const FCIIteration& it) {
        if (!s.fci_started) {
            s.scf_done_line();
            s.w->line();
            s.w->write(fci_header());
            s.fci_started = true;
        }
        s.w->write(fci_iteration_line(it));
        report_progress_text(s.context + strprintf("Davidson iteration %d · E = %.10f E_h · residual %.1e", it.iteration,
                                                   it.energy, it.residual));
    };
    cb.lambda = [&s](const LambdaIteration& it) {
        if (!s.lambda_started) {
            s.w->line("\n    " + s.method + " Lambda equations (Jacobi updates with DIIS)\n    iter   max |d lambda|");
            s.lambda_started = true;
        }
        s.w->line(strprintf("    %4d   %12.3e", it.iteration, it.residual));
        report_progress_text(s.context + strprintf("Lambda iteration %d · max change %.1e", it.iteration, it.residual));
    };
    cb.cphf = [&s](const CPHFProgress& it) {
        if (!s.cphf_started) {
            s.w->line("\n    Orbital response (coupled-perturbed HF / Z-vector, Krylov subspace)\n    iter   largest residual");
            s.cphf_started = true;
        }
        s.w->line(strprintf("    %4d   %12.3e", it.iteration, it.residual));
        report_progress_text(s.context + strprintf("Orbital response iteration %d · residual %.1e", it.iteration, it.residual));
    };
    cb.note = [&s](const std::string& text) {
        s.w->line("    " + text);
        report_progress_text(s.context + text);
    };
    return cb;
}

void write_point_details(OutputWriter& w, const SinglePoint& sp, const BasisSet& b, bool full_orbitals) {
    w.line();
    w.write(scf_result_block(sp.rhf, sp.molecule, b, full_orbitals));
    write_correlated(w, sp);
    for (const auto& x : sp.warnings) w.line("    Warning: " + x);
}

std::string final_summary(const JobResult& r, const MolecularResult& final) {
    std::string s;
    const Molecule& m = r.final_geometry;
    s += strprintf("    Job: %s\n", job_title(r.settings.type).c_str());
    s += "    Final geometry:\n" + parameter_lines(m, "      ");
    s += strprintf("    Final energy (%s): %.10f E_h\n", r.final_energy_label.c_str(), r.final_energy);
    if (r.optimization)
        s += strprintf("    Optimisation: %s, %d steps\n", r.optimization->converged ? "converged" : "NOT converged",
                       r.optimization->steps);
    if (r.scan && !r.scan->points.empty()) {
        const ScanResult& sc = *r.scan;
        const ScanPoint& lo = sc.points[static_cast<std::size_t>(sc.lowest)];
        const ScanPoint& hi = sc.points[static_cast<std::size_t>(sc.highest)];
        s += strprintf("    Relaxed scan of %s: %zu points (the final geometry is the last one); %s\n",
                       sc.plan.coordinate.definition.c_str(), sc.points.size(), sc.message.c_str());
        s += "      lowest point " + std::to_string(lo.index + 1) + " (" + scan_value_text(sc.plan, lo.value) +
             strprintf("), E = %.10f E_h\n", lo.energy);
        s += "      highest point " + std::to_string(hi.index + 1) + " (" + scan_value_text(sc.plan, hi.value) +
             strprintf("), E = %.10f E_h, %.4f kcal/mol above the lowest\n", hi.energy,
                       (hi.energy - lo.energy) * codata2018::kHartreeKcalPerMol);
    }
    if (r.harmonic) {
        s += "    Harmonic wavenumbers (cm-1):";
        for (const auto& md : r.harmonic->modes)
            s += md.imaginary ? strprintf(" %.2fi", std::fabs(md.wavenumber)) : strprintf(" %.2f", md.wavenumber);
        s += "\n";
    }
    if (r.thermo) {
        const ThermoResult& t = *r.thermo;
        const double base = r.final_energy;
        s += strprintf("    Thermochemistry at %.2f K and %.0f Pa (corrections from %s):\n", t.temperature, t.pressure,
                       (method_label(r.point ? r.point->method : final.input.method) + "/" + final.input.basis).c_str());
        s += strprintf("      E + ZPE %20.10f E_h\n", base + joule_per_mol_to_hartree(t.zero_point_energy));
        s += strprintf("      U       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.thermal_energy));
        s += strprintf("      H       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.enthalpy));
        s += strprintf("      G       %20.10f E_h\n", base + joule_per_mol_to_hartree(t.gibbs));
        s += strprintf("      S       %20.4f J/(mol K)\n", t.total.entropy);
    }
    return s;
}

std::string summary_lines(const MolecularResult& m) {
    // The short summary of versions 3.0 and 4.0 (also at the top of psilab's output file).
    std::string s;
    if (m.rhf) {
        s += strprintf("    Hartree–Fock (%s) total energy:  %.10f E_h (%d iterations)\n", m.reference.c_str(),
                       m.rhf->total_energy, m.rhf->iterations);
        if (m.rhf->open_shell()) s += strprintf("    <S^2> = %.6f (S(S+1) = %.6f)\n", m.rhf->s2, m.rhf->s2_exact);
        if (m.so)
            for (const auto& t : m.so->terms)
                if (t.label.rfind("E(", 0) == 0) s += strprintf("    %s correction:  %.10f E_h\n", t.label.c_str(), t.value);
        if (m.so) s += strprintf("    %s total energy:  %.10f E_h\n", job_label(m.so->method, m.reference).c_str(), m.so->total_energy);
        if (m.cc) s += strprintf("    %s total energy:  %.10f E_h\n", m.cc->method.c_str(), m.cc->total_energy);
        if (m.fci) s += strprintf("    Full CI total energy:  %.10f E_h\n", m.fci->total_energy);
    }
    return s;
}

}  // namespace

std::string job_estimate_text(const MolecularResult& p, const JobEstimate& est) {
    std::string s;
    OutputWriter w([&s](const std::string& t) { s += t; });
    write_estimate(w, p, est);
    return s;
}

MolecularResult final_point_result(const JobResult& job, const MolecularInput& input) {
    MolecularResult r;
    r.input = input;
    const SinglePoint* sp = job.single_point ? &*job.single_point : (job.point ? &*job.point : nullptr);
    r.input.method = sp ? sp->method : input.method;
    if (job.single_point) r.input.basis = job.single_point_basis;
    r.molecule = sp ? sp->molecule : job.final_geometry;
    if (r.molecule.name.empty()) r.molecule.name = canonical_molecule_name(input.molecule);
    // 9.0: a custom molecule keeps its atoms (the final geometry); R is the first interatomic distance
    if (r.molecule.custom) {
        r.input.atoms = r.molecule.atoms;
        if (r.molecule.atoms.size() >= 2) {
            double s = 0.0;
            for (std::size_t k = 0; k < 3; ++k)
                s += std::pow(r.molecule.atoms[1].position[k] - r.molecule.atoms[0].position[k], 2);
            r.input.R = std::sqrt(s);
        }
    }
    // geometry parameters of the final point (R of the first bond, the H-O-H angle)
    const SimpleInternals si = r.molecule.custom ? SimpleInternals() : simple_internals(r.molecule);
    const double deg = 180.0 / 3.14159265358979323846;
    if (!si.values.empty()) r.input.R = si.values[0];
    if (si.values.size() == 3) r.input.angle = si.values[2] * deg;
    if (si.values.size() == 6 && r.molecule.name == "H2O2") {
        r.input.r_oo = si.values[0];
        r.input.R = si.values[1];
        r.input.angle = si.values[3] * deg;
        r.input.dihedral = si.values[5] * deg;
    } else if (si.values.size() == 6) {
        r.input.angle = si.values[3] * deg;  // NH3: R(N1-H2), A(H2-N1-H3)
    }
    r.basis = build_basis(r.molecule, r.input.basis, input.basis_directory);
    r.molecule = attach_ecps(r.molecule, r.basis);  // 14.0: effective core potentials
    for (const auto& a : r.basis.atoms)
        if (!a.note.empty() && std::find(r.notes.begin(), r.notes.end(), a.note) == r.notes.end()) r.notes.push_back(a.note);
    if (sp) {
        r.rhf = sp->rhf;
        r.cc = sp->cc;
        r.fci = sp->fci;
        r.so = sp->so;
        r.reference = sp->rhf.reference;
        r.warnings = sp->warnings;
    }
    r.report = molecular_report(r);
    return r;
}

JobResult run_energy_job(MolecularResult& p, const JobSettings& s, const OutputSink& sink, const std::string& echo,
                         const AppendixWriter& appendix) {
    const auto t0 = Clock::now();
    JobResult r;
    r.settings = s;
    OutputWriter w(sink);
    write_header(w, s, echo, p.estimate.threads);
    write_molecule(w, p, "Geometry");
    r.estimate = estimate_job(p, s, nullptr);
    write_estimate(w, p, r.estimate);
    r.initial_geometry = r.final_geometry = p.molecule;
    if (!p.refusal.empty()) {
        r.termination = "refused";
        r.message = p.refusal;
        write_termination(w, r, seconds_since(t0), appendix);
        r.output = w.text();
        return r;
    }
    Streams st;
    st.w = &w;
    st.method = p.input.method;
    st.reference = p.reference;
    w.section("Energy: " + job_label(p.input.method, p.reference) + "/" + p.input.basis);
    const ElectronicCallbacks cb = stream_callbacks(st);
    try {
        run_molecular(p, cb.scf, cb.fci, cb.cc);
        SinglePoint sp;
        sp.molecule = p.molecule;
        sp.method = p.input.method;
        sp.rhf = *p.rhf;
        sp.cc = p.cc;
        sp.fci = p.fci;
        sp.so = p.so;
        sp.reference = p.reference;
        sp.energy = p.so ? p.so->total_energy
                         : (p.cc ? p.cc->total_energy : (p.fci ? p.fci->total_energy : p.rhf->total_energy));
        sp.converged = p.rhf->converged && (!p.cc || p.cc->converged) && (!p.fci || p.fci->converged) &&
                       (!p.so || p.so->converged);
        sp.warnings = p.warnings;
        w.subsection("Results");
        write_point_details(w, sp, p.basis, true);
        r.point = sp;
        r.final_energy = sp.energy;
        r.final_energy_label = job_label(sp.method, p.reference) + "/" + p.input.basis;
        w.section("Summary");
        w.write(summary_lines(p));
        w.write(final_summary(r, p));
        r.completed = sp.converged;
        r.termination = sp.converged ? "normal" : "error";
        if (!sp.converged) r.message = "the calculation did not converge";
    } catch (const Cancelled&) {
        r.termination = "cancelled";
        write_termination(w, r, seconds_since(t0), appendix);
        r.output = w.text();
        throw;
    } catch (const std::exception& e) {
        r.termination = "error";
        r.message = e.what();
    }
    r.seconds = seconds_since(t0);
    write_termination(w, r, r.seconds, appendix);
    r.output = w.text();
    return r;
}

namespace {

// Local minima and maxima inside the scan, from the parabola through each point and its neighbours.
std::vector<ScanExtremum> scan_extrema(const std::vector<ScanPoint>& pts) {
    std::vector<ScanExtremum> out;
    for (std::size_t i = 1; i + 1 < pts.size(); ++i) {
        const double em = pts[i - 1].energy, e0 = pts[i].energy, ep = pts[i + 1].energy;
        const bool mn = e0 < em && e0 < ep, mx = e0 > em && e0 > ep;
        if (!mn && !mx) continue;
        const double h = pts[i + 1].value - pts[i].value;
        const double c = ep - 2.0 * e0 + em;
        ScanExtremum x;
        x.kind = mn ? "minimum" : "maximum";
        x.index = static_cast<int>(i);
        x.value = pts[i].value;
        x.energy = e0;
        if (std::fabs(c) > 1.0e-14) {
            x.value = pts[i].value - 0.5 * h * (ep - em) / c;
            x.energy = e0 - (ep - em) * (ep - em) / (8.0 * c);
        }
        out.push_back(x);
    }
    return out;
}

std::string scan_summary_block(const ScanResult& sr) {
    const ScanPlan& pl = sr.plan;
    std::string s;
    if (sr.points.empty()) return s;
    const double emin = sr.points[static_cast<std::size_t>(sr.lowest)].energy;
    const std::string unit = pl.angular() ? "deg" : "Å";
    s += "    Point  " + pad_text(pl.coordinate.definition + "/" + unit, 16, true) +
         strprintf("  %20s  %12s  %11s  %6s  %s\n", "Energy (E_h)", "dE (kcal/mol)", "dE (kJ/mol)", "Steps", "Converged");
    for (const auto& pt : sr.points) {
        const double de = pt.energy - emin;
        s += strprintf("    %5d  %16.6f  %20.10f  %12.4f  %11.4f  %6s  %s\n", pt.index + 1,
                       pl.angular() ? pt.value : pt.value * kBohrInAngstrom, pt.energy,
                       de * codata2018::kHartreeKcalPerMol, de * codata2018::kHartreeKJPerMol,
                       sr.relaxed ? std::to_string(pt.steps).c_str() : "--",
                       sr.relaxed ? (pt.converged ? "yes" : "NO") : "(single point)");
    }
    const auto& lo = sr.points[static_cast<std::size_t>(sr.lowest)];
    const auto& hi = sr.points[static_cast<std::size_t>(sr.highest)];
    s += strprintf("    dE relative to the lowest point, point %d (%s); highest point %d (%s), %.4f kcal/mol "
                   "(%.4f kJ/mol) above it\n",
                   lo.index + 1, scan_value_text(pl, lo.value).c_str(), hi.index + 1, scan_value_text(pl, hi.value).c_str(),
                   (hi.energy - lo.energy) * codata2018::kHartreeKcalPerMol,
                   (hi.energy - lo.energy) * codata2018::kHartreeKJPerMol);
    // the other parameters along the scan
    if (sr.relaxed && !sr.points.front().parameters.empty()) {
        const auto& names = sr.points.front().parameters;
        s += "\n    Geometry along the scan (Å and degrees)\n    Point";
        for (const auto& p : names) s += "  " + pad_text(p.name, 16, true);
        s += "\n";
        for (const auto& pt : sr.points) {
            s += strprintf("    %5d", pt.index + 1);
            for (const auto& n : names) {
                auto it = std::find_if(pt.parameters.begin(), pt.parameters.end(),
                                       [&](const GeometricParameter& q) { return q.name == n.name; });
                if (it == pt.parameters.end()) s += strprintf("  %16s", "--");
                else s += strprintf("  %16.6f", it->angle ? (std::fabs(it->value) < 5.0e-7 ? 0.0 : it->value)
                                                          : it->value * kBohrInAngstrom);
            }
            s += "\n";
        }
    }
    if (!sr.extrema.empty()) {
        s += "\n    Stationary points along the scan (parabola through the point and its two neighbours)\n";
        for (const auto& x : sr.extrema)
            s += strprintf("    %-8s near %s = %s, E = %.10f E_h (%+.4f kcal/mol relative to the lowest point), at point %d\n",
                           x.kind.c_str(), pl.coordinate.definition.c_str(), scan_value_text(pl, x.value).c_str(), x.energy,
                           (x.energy - emin) * codata2018::kHartreeKcalPerMol, x.index + 1);
        s += "    (estimates from the scan; a transition-state optimisation locates a maximum exactly)\n";
    }
    return s;
}

}  // namespace

static void run_scan(OutputWriter& w, Streams& st, ElectronicCallbacks& cb, const MolecularResult& p, const JobSettings& s,
                     const ElectronicOptions& eo, JobResult& r) {
    const auto t0 = Clock::now();
    const std::string method = p.input.method;
    ScanResult sr;
    {
        const std::string why = plan_scan(p.molecule, s, &sr.plan);
        if (!why.empty()) throw std::runtime_error("the scan is not possible: " + why);
    }
    const ScanPlan& pl = sr.plan;
    const int npts = pl.steps + 1;
    sr.relaxed = p.molecule.atoms.size() > 2;
    const ConvergenceCriteria c = convergence_criteria(s.tight);
    w.section("Relaxed potential energy surface scan: " + job_label(method, p.reference) + "/" + p.input.basis);
    std::string gaussian = pl.coordinate.letter;
    for (int a : pl.coordinate.coordinate.atoms) gaussian += " " + std::to_string(a + 1);
    w.line(strprintf("    Scanned coordinate   %s  (Gaussian: %s S %d %g%s)", pl.coordinate.definition.c_str(),
                     gaussian.c_str(), pl.steps, pl.angular() ? pl.step : pl.step * kBohrInAngstrom,
                     pl.angular() ? "" : " in Å"));
    w.line("    At the input         " + scan_value_text(pl, pl.input_value));
    w.line(strprintf("    Points               %d, from %s to %s in steps of %s", npts, scan_value_text(pl, pl.values.front()).c_str(),
                     scan_value_text(pl, pl.values.back()).c_str(),
                     (pl.angular() ? strprintf("%.6f deg", pl.step) : strprintf("%.8f Å", pl.step * kBohrInAngstrom)).c_str()));
    if (sr.relaxed) {
        w.line("    At each point the scanned coordinate is frozen at its value and all other coordinates are optimised");
        w.line("    (redundant internal coordinates; projector P~ = P - P C (C' P C)^-1 C' P removes the frozen coordinate");
        w.line("    from the RFO step and the convergence test). Each point starts from the optimised geometry of the");
        w.line("    previous one, with the scanned coordinate set to its new value, and from its updated Hessian.");
        w.line(strprintf("    %s convergence, at most %d steps per point; thresholds: max force %.1e, RMS force %.1e; "
                         "max step %.1e, RMS step %.1e; energy change %.0e E_h; constraint %.1e",
                         s.tight ? "Tight" : "Normal", s.max_steps, c.max_force, c.rms_force, c.max_displacement,
                         c.rms_displacement, c.energy, 0.01 * c.max_displacement));
    } else {
        w.line("    A diatomic has no other coordinate to relax: the relaxed scan is a series of single-point energies.");
    }
    Molecule geom = p.molecule;
    Matrix orbitals, hessian;
    std::vector<InternalCoordinate> hessian_coords;
    std::vector<int> failed;
    for (int k = 0; k < npts; ++k) {
        check_cancelled();
        const auto tp = Clock::now();
        ScanPoint pt;
        pt.index = k;
        pt.value = pl.values[static_cast<std::size_t>(k)];
        w.subsection(strprintf("Scan point %d of %d: %s = %s", k + 1, npts, pl.coordinate.definition.c_str(),
                               scan_value_text(pl, pt.value).c_str()));
        CoordinateConstraint con;
        con.coordinate = pl.coordinate.coordinate;
        con.target = pl.targets[static_cast<std::size_t>(k)];
        double residual = 0.0;
        geom = apply_constraints(geom, {con}, &residual);
        if (residual > 1.0e-6)
            pt.warnings.push_back(strprintf("The coordinate could be set to its value only within %.1e.", residual));
        const std::string ctx = strprintf("Scan point %d of %d · ", k + 1, npts);
        report_progress(JobStage::Scan, k + 1, k > 0 ? sr.points.back().energy : 0.0, npts);
        if (!sr.relaxed) {
            st.reset();
            st.context = ctx;
            w.write(geometry_table(geom));
            report_progress_text(ctx + "energy");
            const SinglePoint sp = run_single_point(geom, basis_at_geometry(p.basis, geom), method, eo,
                                                    orbitals.rows() ? &orbitals : nullptr, cb);
            if (!sp.converged) throw std::runtime_error(strprintf("the energy at scan point %d did not converge", k + 1));
            st.scf_done_line();
            orbitals = sp.rhf.C;
            pt.energy = sp.energy;
            pt.converged = true;
            pt.geometry = geom;
            pt.measured = internal_value(geom, pl.coordinate.coordinate);
            w.line(strprintf("\n    Energy %.10f E_h", sp.energy));
            r.point = sp;
        } else {
            OptimizationOptions oo;
            oo.tight = s.tight;
            oo.max_steps = s.max_steps;
            oo.initial_hessian = s.initial_hessian;
            oo.constraints = {con};
            oo.guess_orbitals = orbitals.rows() ? &orbitals : nullptr;
            if (hessian.rows()) {
                oo.initial_hessian_matrix = &hessian;
                oo.hessian_coordinates = hessian_coords;
            }
            oo.on_step_start = [&](int step, const Molecule& g) {
                st.reset();
                st.context = ctx + strprintf("step %d · ", step);
                w.line(strprintf("\n  Scan point %d, optimisation step %d", k + 1, step));
                w.write(geometry_table(g));
                report_progress_text(st.context + "energy and gradient");
            };
            const OptimizationResult opt = optimize_geometry(geom, p.basis, method, eo, oo, cb, [&](const OptimizationStep& stp) {
                w.line(stp.step == 1 ? strprintf("\n    Energy %.10f E_h", stp.energy)
                                     : strprintf("\n    Energy %.10f E_h (change %+.3e E_h)", stp.energy, stp.delta_e));
                w.write(gradient_block(stp.geometry, stp.cartesian_gradient, stp.gradient_description));
                w.write(optimization_step_block(stp, c));
                report_progress_text(ctx + strprintf("step %d · E = %.10f E_h · max force %.2e", stp.step, stp.energy,
                                                     stp.convergence[0].value));
            });
            if (opt.history.empty() || !opt.final_point.point.converged)
                throw std::runtime_error(strprintf("the %s calculation at scan point %d did not converge", method.c_str(), k + 1));
            orbitals = opt.final_point.point.rhf.C;
            hessian = opt.final_hessian;
            hessian_coords = opt.final_coordinates;
            pt.energy = opt.final_energy;
            pt.converged = opt.converged;
            pt.steps = opt.steps;
            pt.geometry = opt.final_geometry;
            pt.max_force = opt.history.back().convergence[0].value;
            for (const auto& h : opt.history) {
                pt.step_energies.push_back(h.energy);
                pt.step_forces.push_back(h.convergence[0].value);
                for (const auto& x2 : h.warnings) pt.warnings.push_back(strprintf("step %d: ", h.step) + x2);
            }
            pt.measured = internal_value(opt.final_geometry, pl.coordinate.coordinate);
            geom = opt.final_geometry;
            r.point = opt.final_point.point;
            if (!opt.converged) {
                failed.push_back(k + 1);
                pt.warnings.push_back(opt.message);
            }
        }
        if (pl.coordinate.coordinate.type != InternalType::Stretch) {
            pt.measured *= 180.0 / 3.14159265358979323846;
            // report the measured dihedral on the branch of the nominal value
            if (pl.coordinate.coordinate.type == InternalType::Dihedral)
                pt.measured = pt.value + wrap_angle((pt.measured - pt.value) * 3.14159265358979323846 / 180.0) * 180.0 /
                                             3.14159265358979323846;
        }
        pt.parameters = geometric_parameters(pt.geometry);
        pt.seconds = seconds_since(tp);
        w.line(strprintf("\n    Point %d: E = %.10f E_h, %s = %s (deviation %.1e)%s", k + 1, pt.energy,
                         pl.coordinate.definition.c_str(), scan_value_text(pl, pt.measured).c_str(),
                         std::fabs(pt.measured - pt.value),
                         sr.relaxed ? strprintf(", %d optimisation step%s, %s", pt.steps, pt.steps == 1 ? "" : "s",
                                                pt.converged ? "converged" : "NOT converged").c_str()
                                    : ""));
        if (sr.relaxed) w.write(parameter_lines(pt.geometry, "      "));
        for (const auto& x2 : pt.warnings) w.line("    Warning: " + x2);
        sr.points.push_back(pt);
    }
    sr.lowest = sr.highest = 0;
    for (std::size_t i = 0; i < sr.points.size(); ++i) {
        if (sr.points[i].energy < sr.points[static_cast<std::size_t>(sr.lowest)].energy) sr.lowest = static_cast<int>(i);
        if (sr.points[i].energy > sr.points[static_cast<std::size_t>(sr.highest)].energy) sr.highest = static_cast<int>(i);
    }
    sr.extrema = scan_extrema(sr.points);
    sr.converged = failed.empty();
    if (!failed.empty()) {
        std::string list;
        for (int f : failed) list += (list.empty() ? "" : ", ") + std::to_string(f);
        sr.message = strprintf("The constrained optimisation did not converge within %d steps at point%s %s; the other "
                               "points are complete.",
                               s.max_steps, failed.size() == 1 ? "" : "s", list.c_str());
    } else {
        sr.message = strprintf("All %d points %s.", npts, sr.relaxed ? "converged" : "computed");
    }
    sr.seconds = seconds_since(t0);
    w.section("Summary of the relaxed scan: " + job_label(method, p.reference) + "/" + p.input.basis);
    w.write(scan_summary_block(sr));
    w.line("    " + sr.message + " Total time " + duration_string(sr.seconds) + ".");
    const ScanPoint& last = sr.points.back();
    r.final_geometry = last.geometry;
    r.final_energy = last.energy;
    r.final_energy_label = job_label(method, p.reference) + "/" + p.input.basis;
    if (r.point) {
        st.reset();
        w.subsection("Electronic structure at the last scan point");
        write_point_details(w, *r.point, basis_at_geometry(p.basis, last.geometry), true);
    }
    r.scan = sr;
    if (!sr.converged) {
        r.termination = "error";
        r.message = sr.message;
        throw std::runtime_error(sr.message);
    }
}

JobResult run_job(const MolecularInput& input, const JobSettings& s, const OutputSink& sink, const std::string& echo,
                  const AppendixWriter& appendix) {
    const auto t0 = Clock::now();
    JobResult r;
    r.settings = s;
    MolecularResult p = prepare_molecular(input);
    std::optional<MolecularResult> psp;
    if (job_has_single_point(s.type)) {
        MolecularInput in2 = input;
        in2.method = s.sp_method.empty() ? input.method : s.sp_method;
        in2.basis = s.sp_basis.empty() ? input.basis : s.sp_basis;
        psp = prepare_molecular(in2);
        r.single_point_method = psp->input.method;
        r.single_point_basis = psp->input.basis;
    }
    OutputWriter w(sink);
    write_header(w, s, echo, p.estimate.threads);
    write_molecule(w, p, job_uses_optimizer(s.type) ? "Starting geometry" : "Geometry (input)");
    r.estimate = estimate_job(p, s, psp ? &*psp : nullptr);
    write_estimate(w, p, r.estimate);
    r.initial_geometry = r.final_geometry = p.molecule;
    if (!r.estimate.refusal.empty()) {
        r.termination = "refused";
        r.message = r.estimate.refusal;
        write_termination(w, r, seconds_since(t0), appendix);
        r.output = w.text();
        return r;
    }
    const std::string method = p.input.method;
    const ElectronicOptions eo = electronic_options(p);
    Streams st;
    st.w = &w;
    st.method = method;
    st.reference = p.reference;
    ElectronicCallbacks cb = stream_callbacks(st);
    Matrix orbitals;
    const Matrix* guess = nullptr;
    try {
        Molecule x = p.molecule;
        // ---- relaxed scan (7.0)
        if (job_has_scan(s.type)) run_scan(w, st, cb, p, s, eo, r);
        // ---- optimisation
        if (job_has_optimization(s.type)) {
            const bool ts = job_has_ts(s.type);
            const std::string label = ts ? "TS optimisation step" : "Optimisation step";
            w.section((ts ? "Transition-state optimisation: " : "Geometry optimisation: ") + job_label(method, p.reference) + "/" +
                      p.input.basis);
            const ConvergenceCriteria c = convergence_criteria(s.tight);
            if (ts) {
                w.line("    Berny-type search for a first-order saddle point (Gaussian: Opt=TS) in redundant internal");
                w.line(std::string("    coordinates: Hessian computed ") +
                       (s.ts_hessian == "calcall" ? std::string("at every step (CalcAll)")
                        : s.ts_hessian == "recalcfc"
                            ? strprintf("at step 1 and every %d steps (RecalcFC=%d), Bofill updates between", s.recalc_fc, s.recalc_fc)
                            : std::string("at the first step (CalcFC), then Bofill updates")) +
                       "; partitioned RFO steps (Baker 1986) that");
                w.line("    maximise along the followed Hessian eigenvector (the lowest at the start, then the one most");
                w.line(strprintf("    like the previous one) and minimise along the others; trust radius %.2f bohr/rad (fixed);",
                                 OptimizationOptions().trust_radius));
                w.line(strprintf("    eigenvalue test %s; %s convergence, at most %d steps",
                                 s.eigen_test ? "on (exactly one negative eigenvalue)" : "off (NoEigenTest)",
                                 s.tight ? "tight" : "normal", s.max_steps));
            } else {
                w.line(strprintf("    Redundant internal coordinates, %s steps with a trust radius, BFGS updates; %s "
                                 "convergence, at most %d steps",
                                 s.opt_algorithm == "gediis" ? "hybrid GEDIIS/RFO (Li and Frisch 2006)" : "RFO",
                                 s.tight ? "tight" : "normal", s.max_steps));
                if (!s.opt_hessian.empty())
                    w.line(std::string("    Hessian computed ") +
                           (s.opt_hessian == "calcall" ? std::string("at every step (CalcAll)")
                            : s.opt_hessian == "recalcfc"
                                ? strprintf("at step 1 and every %d steps (RecalcFC=%d), BFGS updates between", s.recalc_fc,
                                            s.recalc_fc)
                                : std::string("at the first step (CalcFC), then BFGS updates")) +
                           " (negative eigenvalues replaced by their magnitude)");
            }
            w.line(strprintf("    thresholds: max force %.1e, RMS force %.1e E_h/bohr; max step %.1e, RMS step %.1e bohr; "
                             "energy change %.0e E_h",
                             c.max_force, c.rms_force, c.max_displacement, c.rms_displacement, c.energy));
            OptimizationOptions oo;
            oo.tight = s.tight;
            oo.max_steps = s.max_steps;
            oo.transition_state = ts;
            oo.initial_hessian = s.initial_hessian;
            oo.ts_hessian = s.ts_hessian;
            oo.computed_hessian = ts ? "" : s.opt_hessian;
            oo.recalc_every = s.recalc_fc;
            oo.algorithm = s.opt_algorithm;
            oo.numerical_hessian = s.hessian == "numerical";
            oo.eigen_test = s.eigen_test;
            oo.on_step_start = [&](int step, const Molecule& geom) {
                st.reset();
                st.context = strprintf("%s %d · ", label.c_str(), step);
                w.subsection(strprintf("%s %d", label.c_str(), step));
                w.write(geometry_table(geom));
                report_progress_text(st.context + "energy and gradient");
            };
            const OptimizationResult opt = optimize_geometry(x, p.basis, method, eo, oo, cb, [&](const OptimizationStep& stp) {
                w.line(stp.step == 1 ? strprintf("\n    Energy %.10f E_h", stp.energy)
                                     : strprintf("\n    Energy %.10f E_h (change %+.3e E_h)", stp.energy, stp.delta_e));
                w.write(gradient_block(stp.geometry, stp.cartesian_gradient, stp.gradient_description));
                w.write(optimization_step_block(stp, c));
                report_progress_text(strprintf("%s %d · E = %.10f E_h · max force %.2e", label.c_str(), stp.step, stp.energy,
                                               stp.convergence[0].value));
            });
            w.write(optimization_result_block(opt));
            r.optimization = opt;
            x = opt.final_geometry;
            r.final_geometry = x;
            r.final_energy = opt.final_energy;
            r.final_energy_label = job_label(method, p.reference) + "/" + p.input.basis;
            orbitals = opt.final_point.point.rhf.C;
            guess = &orbitals;
            {
                SinglePoint sp = opt.final_point.point;
                r.point = sp;
                st.reset();
                w.subsection("Electronic structure at the " +
                             std::string(opt.converged ? (ts ? "transition state" : "optimised geometry") : "last geometry"));
                write_point_details(w, sp, basis_at_geometry(p.basis, x), true);
                if (opt.final_point.has_dipole && method != "RHF") {
                    const auto& d = opt.final_point.dipole;
                    const double nd = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                    w.line("\n    Dipole moment (relaxed " + method + " density, origin at 0)");
                    w.line(strprintf("      x %12.8f   y %12.8f   z %12.8f   |mu| %12.8f e bohr   (%.6f D)", d[0], d[1], d[2],
                                     nd, nd / codata2018::kDebyeInAtomicUnits));
                }
            }
            if (!opt.converged) {
                r.termination = "error";
                r.message = opt.message + (s.type == "ts+freq" ? " The frequency step was not run: it needs a transition state."
                                           : s.type == "ts+freq+sp" ? " The frequency and single-point steps were not run: they "
                                                                      "need a transition state."
                                           : (s.type != "opt" && s.type != "ts" ? " The frequency and single-point steps were "
                                                                                    "not run: they need a stationary point."
                                                                                  : ""));
                throw std::runtime_error(r.message);
            }
        }
        // ---- frequencies
        if (job_has_frequencies(s.type)) {
            const bool at_input = s.type == "freq";
            w.section("Harmonic frequencies: " + job_label(method, p.reference) + "/" + p.input.basis);
            w.line(std::string("    at the ") +
                   (at_input ? "input geometry" : (job_has_ts(s.type) ? "transition state" : "optimised geometry")));
            w.write(geometry_table(x));
            st.reset();
            st.context = "Hessian · ";
            HessianOptions ho;
            ho.force_numerical = s.hessian == "numerical";
            if (r.optimization && r.optimization->converged) {
                ho.reference = &r.optimization->final_point;  // energy and gradient of the last step
                w.line("    energy, gradient and orbitals at this geometry: from the last optimisation step");
            }
            int last_index = 0;
            const HessianResult h = run_hessian(x, basis_at_geometry(p.basis, x), method, eo, ho, cb, guess,
                                                [&](int i, int n, const std::string& what) {
                                                    if (i == 1) w.line("\n    Displaced calculations");
                                                    last_index = i;
                                                    w.line(strprintf("    displacement %d of %d: %s", i, n, what.c_str()));
                                                    report_progress_text(strprintf("Hessian: displacement %d of %d (%s)", i, n, what.c_str()));
                                                    st.reset();
                                                });
            (void)last_index;
            r.hessian = h;
            if (at_input) {
                r.point = h.reference.point;
                r.final_energy = h.reference.point.energy;
                r.final_energy_label = job_label(method, p.reference) + "/" + p.input.basis;
                w.subsection("Electronic structure at the input geometry");
                write_point_details(w, h.reference.point, basis_at_geometry(p.basis, x), true);
                if (h.reference.has_dipole && method != "RHF") {
                    const auto& d = h.reference.dipole;
                    const double nd = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                    w.line("\n    Dipole moment (relaxed " + method + " density, origin at 0)");
                    w.line(strprintf("      x %12.8f   y %12.8f   z %12.8f   |mu| %12.8f e bohr   (%.6f D)", d[0], d[1], d[2],
                                     nd, nd / codata2018::kDebyeInAtomicUnits));
                }
                w.line();
                w.write(gradient_block(x, h.reference.gradient, h.reference.description));
                double mx = 0.0;
                for (double g : h.reference.gradient) mx = std::max(mx, std::fabs(g));
                if (mx > convergence_criteria(false).max_force) {
                    const std::string warn = strprintf(
                        "The input geometry is not a stationary point (largest gradient component %.2e E_h/bohr): the "
                        "frequencies and the thermochemistry refer to a non-stationary geometry. Optimise first "
                        "(job opt+freq).",
                        mx);
                    w.line("    Warning: " + warn);
                    r.warnings.push_back(warn);
                }
            }
            HarmonicAnalysis a = harmonic_analysis(x, h.hessian, h.has_dipole_derivatives ? &h.dipole_derivatives : nullptr);
            w.line();
            w.write(frequency_block(x, h, a));
            report_progress_text("Thermochemistry");
            ThermoResult t = thermochemistry(x, a, r.final_energy, s.temperature, s.pressure, symmetry_number(x));
            w.section("Thermochemistry (ideal gas, rigid rotor, harmonic oscillators)");
            w.write(thermochemistry_block(t));
            for (const auto& x2 : t.warnings) r.warnings.push_back(x2);
            r.harmonic = a;
            r.thermo = t;
            if (job_has_ts(s.type)) {
                // a transition state has exactly one imaginary frequency
                int imaginary = 0;
                for (const auto& md : a.modes) imaginary += md.imaginary ? 1 : 0;
                w.section("Transition-state check");
                if (imaginary == 1) {
                    const NormalMode& md = a.modes.front();
                    w.line(strprintf("    One imaginary frequency, %.2fi cm-1%s: a first-order saddle point (transition state).",
                                     std::fabs(md.wavenumber), md.symmetry.empty() ? "" : (" (" + md.symmetry + ")").c_str()));
                    w.line("    Its normal mode (normalised Cartesian displacements) is the reaction coordinate at the "
                           "transition state:");
                    for (std::size_t i = 0; i < x.atoms.size(); ++i)
                        w.line(strprintf("      %-2s%-3zu %12.6f %12.6f %12.6f", x.atoms[i].symbol.c_str(), i + 1,
                                         md.displacement[i][0], md.displacement[i][1], md.displacement[i][2]));
                } else {
                    r.termination = "error";
                    r.message = strprintf(
                        "The frequencies give %d imaginary frequencies, not one: the stationary point is %s, not a "
                        "transition state.",
                        imaginary, imaginary == 0 ? "a minimum" : strprintf("a saddle point of order %d", imaginary).c_str());
                    w.line("    " + r.message);
                    throw std::runtime_error(r.message);
                }
            }
        }
        // ---- single point
        if (job_has_single_point(s.type)) {
            const std::string m2 = psp->input.method, b2 = psp->input.basis;
            w.section(std::string("Single-point energy at the ") +
                      (job_has_ts(s.type) ? "transition state: " : "optimised geometry: ") + job_label(m2, psp->reference) + "/" + b2);
            w.write(geometry_table(x));
            const BasisSet bs = build_basis(x, b2, p.input.basis_directory);
            const Molecule xs = attach_ecps(x, bs);  // 14.0: the ECPs of the single point's basis set
            w.write(basis_block(bs, psp->notes));
            st.reset();
            st.method = m2;
            st.reference = psp->reference;
            st.context = "Single point · ";
            ElectronicCallbacks cb2 = stream_callbacks(st);
            ElectronicOptions eo2 = eo;
            eo2.reference = psp->reference;
            const SinglePoint sp = run_single_point(xs, bs, m2, eo2, b2 == p.input.basis && eo2.reference == eo.reference ? guess : nullptr, cb2);
            write_point_details(w, sp, bs, true);
            r.single_point = sp;
            r.final_energy = sp.energy;
            r.final_energy_label = job_label(m2, psp->reference) + "/" + b2;
            if (!sp.converged) throw std::runtime_error("the single-point calculation did not converge");
            if (r.thermo) {
                // the final energies: single point + thermal corrections of the frequency step
                ThermoResult& t = *r.thermo;
                t.electronic_energy = sp.energy;
                t.total_zpe = sp.energy + joule_per_mol_to_hartree(t.zero_point_energy);
                t.total_u = sp.energy + joule_per_mol_to_hartree(t.thermal_energy);
                t.total_h = sp.energy + joule_per_mol_to_hartree(t.enthalpy);
                t.total_g = sp.energy + joule_per_mol_to_hartree(t.gibbs);
                t.total_a = sp.energy + joule_per_mol_to_hartree(t.helmholtz);
                w.line();
                w.line("    Final energies: " + r.final_energy_label + " single point + thermal corrections of " +
                       job_label(method, p.reference) + "/" + p.input.basis + strprintf(" (%.2f K, %.0f Pa)", t.temperature, t.pressure));
                w.line(strprintf("      E(%s)%*s %20.10f E_h", r.final_energy_label.c_str(),
                                 std::max(0, 24 - static_cast<int>(r.final_energy_label.size())), "", sp.energy));
                w.line(strprintf("      E + ZPE                      %20.10f E_h", t.total_zpe));
                w.line(strprintf("      U = E + thermal energy       %20.10f E_h", t.total_u));
                w.line(strprintf("      H = E + thermal enthalpy     %20.10f E_h", t.total_h));
                w.line(strprintf("      A = E + thermal Helmholtz    %20.10f E_h", t.total_a));
                w.line(strprintf("      G = E + thermal Gibbs        %20.10f E_h", t.total_g));
            }
        }
        r.completed = true;
        r.termination = "normal";
    } catch (const Cancelled&) {
        r.termination = "cancelled";
        r.seconds = seconds_since(t0);
        write_termination(w, r, r.seconds, appendix);
        r.output = w.text();
        throw;
    } catch (const std::exception& e) {
        r.termination = "error";
        if (r.message.empty()) r.message = e.what();
    }
    // ---- summary
    w.section("Summary");
    const MolecularResult fp = final_point_result(r, p.input);
    w.write(summary_lines(fp));
    w.write(final_summary(r, p));
    for (const auto& x : r.warnings) w.line("    Warning: " + x);
    r.seconds = seconds_since(t0);
    write_termination(w, r, r.seconds, appendix);
    r.output = w.text();
    return r;
}

}  // namespace macqc
