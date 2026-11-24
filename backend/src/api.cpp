#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <exception>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>

#include "macqc/calculation.h"
#include "macqc/job.h"
#include "macqc/jobs.h"
#include "macqc/macqc_api.h"
#include "macqc/nuclear_motion.h"
#include "macqc/orbitals.h"
#include "macqc/output.h"
#include "macqc/structure_file.h"
#include "nlohmann/json.hpp"

namespace macqc {

using nlohmann::json;

namespace {

json matrix_json(const Matrix& m) {
    json rows = json::array();
    for (std::size_t i = 0; i < m.rows(); ++i) {
        json row = json::array();
        for (std::size_t j = 0; j < m.cols(); ++j) row.push_back(m(i, j));
        rows.push_back(row);
    }
    return rows;
}

json input_json_general(const Input& in);

// 9.0: atoms of a custom molecule as in the input, {"symbol", "Z", "x", "y", "z"} in angstrom.
json custom_atoms_json(const std::vector<Atom>& atoms) {
    json a = json::array();
    for (const Atom& at : atoms)
        a.push_back({{"symbol", at.symbol.empty() ? element_symbol(at.Z) : at.symbol},
                     {"Z", at.Z},
                     {"x", at.position[0] * kBohrInAngstrom},
                     {"y", at.position[1] * kBohrInAngstrom},
                     {"z", at.position[2] * kBohrInAngstrom}});
    return a;
}

json input_json(const Input& in) {
    if (!is_heh_input(in) && (!canonical_molecule_name(in.molecule).empty() || find_molecule_preset(in.molecule))) {
        json j = input_json_general(in);
        if (canonical_molecule_name(in.molecule) == "custom" || !in.atoms.empty()) {
            j["atoms"] = custom_atoms_json(in.atoms);
            j["title"] = in.title;
        }
        if (canonical_molecule_name(in.molecule) == "custom" || in.chargeGiven) j["charge"] = in.charge;
        if (canonical_molecule_name(in.molecule) == "custom" || in.multiplicityGiven) j["multiplicity"] = in.multiplicity;
        j["reference"] = in.reference.empty() ? json(nullptr) : json(in.reference);
        if (canonical_molecule_name(in.molecule) == "HeH+") j["szaboOstlund"] = false;
        return j;
    }
    return {{"R", in.R},           {"zeta1", in.zeta1},         {"zeta2", in.zeta2},
            {"ZA", in.ZA},         {"ZB", in.ZB},               {"nGaussians", in.nGaussians},
            {"exactErf", in.exactErf}, {"maxIter", in.maxIter}, {"convergence", in.convergence},
            {"printLevel", in.printLevel}};
}

json input_json_general(const Input& in) {
    {
        return {{"molecule", in.molecule},
                {"basis", in.basis},
                {"R", in.R},
                {"angle", in.angle},
                {"ROO", canonical_molecule_name(in.molecule) == "H2O2" ? json(in.ROO) : json(nullptr)},
                {"dihedral", canonical_molecule_name(in.molecule) == "H2O2" ? json(in.dihedral) : json(nullptr)},
                {"method", in.method},
                {"basisDirectory", in.basisDirectory},
                {"threads", in.threads},
                {"memoryMB", in.memoryMB},
                {"fciMemoryMB", in.fciMemoryMB},
                {"scfMaxIter", in.scfMaxIter},
                {"energyConvergence", in.energyConvergence},
                {"gradientConvergence", in.gradientConvergence},
                {"estimateOnly", in.estimateOnly},
                {"printLevel", in.printLevel},
                {"job", in.job},
                {"spMethod", in.spMethod},
                {"spBasis", in.spBasis},
                {"temperature", in.temperature},
                {"pressure", in.pressure},
                {"optConvergence", in.optConvergence},
                {"maxSteps", in.maxSteps},
                {"scan", in.scan.empty() ? json(nullptr) : json(in.scan)},
                {"scanStart", in.scanStart ? json(*in.scanStart) : json(nullptr)},
                {"tsHessian", in.tsHessian},
                {"eigenTest", in.eigenTest},
                {"hessian", in.hessian},
                {"initialHessian", in.initialHessian},
                {"optHessian", in.optHessian},
                {"recalcFC", in.recalcFC},
                {"optAlgorithm", in.optAlgorithm}};
    }
}

double get_number(const json& j, const char* key, double fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_number()) throw std::runtime_error(std::string("\"") + key + "\" must be a number");
    return it->get<double>();
}

int get_int(const json& j, const char* key, int fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_number()) throw std::runtime_error(std::string("\"") + key + "\" must be an integer");
    const double v = it->get<double>();
    if (v != std::floor(v) || std::fabs(v) > 1.0e9)
        throw std::runtime_error(std::string("\"") + key + "\" must be an integer");
    return static_cast<int>(v);
}

std::string get_string(const json& j, const char* key, const std::string& fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_string()) throw std::runtime_error(std::string("\"") + key + "\" must be a string");
    return it->get<std::string>();
}

json molecular_json(const MolecularResult& m) {
    json mol;
    mol["name"] = m.molecule.name;
    mol["formula"] = m.molecule.formula();
    mol["custom"] = m.molecule.custom;
    mol["title"] = m.molecule.custom ? json(m.molecule.title) : json(nullptr);
    mol["electrons"] = m.molecule.electrons();          // 14.0: without the core electrons of ECPs
    mol["coreElectrons"] = m.molecule.core_electrons();  // 14.0: replaced by ECPs
    mol["charge"] = m.molecule.charge;
    mol["R"] = m.input.R;
    const bool has_angle = m.molecule.name == "H2O" || m.molecule.name == "H2O2" || m.molecule.name == "NH3";
    mol["angle"] = has_angle ? json(m.input.angle) : json(nullptr);
    mol["ROO"] = m.molecule.name == "H2O2" ? json(m.input.r_oo) : json(nullptr);
    mol["dihedral"] = m.molecule.name == "H2O2" ? json(m.input.dihedral) : json(nullptr);
    mol["nuclearRepulsion"] = m.molecule.nuclear_repulsion();
    mol["multiplicity"] = m.input.multiplicity;  // 15.0
    mol["reference"] = m.reference;
    mol["preset"] = m.input.preset.empty() ? json(nullptr) : json(m.input.preset);
    json atoms = json::array();
    for (const Atom& a : m.molecule.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    mol["atoms"] = atoms;

    json basis;
    basis["name"] = m.basis.name;
    basis["functions"] = m.basis.n_functions;
    basis["shells"] = m.basis.shells.size();
    basis["maxL"] = m.basis.max_l;
    basis["cartesianD"] = m.basis.cartesian_d;                    // 14.0: 6D (the 6-31G family)
    basis["functionTypes"] = basis_function_types(m.basis.name);  // 14.0: e.g. "spherical harmonics (5D, 7F)"
    json ba = json::array();
    for (const auto& a : m.basis.atoms)
        ba.push_back({{"symbol", a.symbol},
                      {"basis", a.basis},
                      {"contraction", a.contraction},
                      {"primitives", a.primitives},
                      {"functions", a.functions},
                      {"ecp", a.ecp.empty() ? json(nullptr) : json(a.ecp)},  // 14.0: e.g. "ECP28"
                      {"ecpElectrons", a.ecp_electrons}});
    basis["atoms"] = ba;

    const CostEstimate& e = m.estimate;
    json est = {{"basisFunctions", e.basis_functions},
                {"uniqueIntegrals", e.unique_eris},
                {"integralMemoryMB", e.eri_memory_mb},
                {"integralsStored", e.eri_stored},
                {"threads", e.threads},
                {"rhfSeconds", e.rhf_seconds},
                {"totalSeconds", e.total_seconds},
                {"longCalculation", e.total_seconds > kLongCalculationSeconds}};
    if (e.fci) {
        const FCISize& f = *e.fci;
        est["fci"] = {{"orbitals", f.orbitals},
                      {"alphaElectrons", f.alpha},
                      {"betaElectrons", f.beta},
                      {"alphaStrings", f.alpha_strings},
                      {"betaStrings", f.beta_strings},
                      {"determinants", f.determinants},
                      {"memoryMB", f.memory_bytes / (1024.0 * 1024.0)},
                      {"seconds", f.seconds},
                      {"feasible", f.feasible},
                      {"reason", f.reason.empty() ? json(nullptr) : json(f.reason)}};
    } else {
        est["fci"] = nullptr;
    }
    if (e.cc) {
        const CCSize& c = *e.cc;
        est["coupledCluster"] = {{"method", c.method},
                                 {"engine", c.engine},
                                 {"occupied", c.occupied},
                                 {"virtuals", c.virtuals},
                                 {"amplitudes", c.amplitudes},
                                 {"memoryMB", c.memory_bytes / (1024.0 * 1024.0)},
                                 {"ladderIntegralsStored", c.ladder_stored},
                                 {"iterations", c.iterations},
                                 {"iterationSeconds", c.iteration_seconds},
                                 {"triplesSeconds", c.triples_seconds},
                                 {"seconds", c.seconds},
                                 {"feasible", c.feasible},
                                 {"reason", c.reason.empty() ? json(nullptr) : json(c.reason)}};
    } else {
        est["coupledCluster"] = nullptr;
    }
    if (e.so) {
        const CostEstimate::SOCost& c = *e.so;
        est["correlated"] = {{"method", c.method},
                             {"engine", "spin-orbital, dense spin blocks"},
                             {"amplitudes", c.amplitudes},
                             {"memoryMB", c.memory_bytes / (1024.0 * 1024.0)},
                             {"seconds", c.seconds},
                             {"feasible", c.feasible},
                             {"reason", c.reason.empty() ? json(nullptr) : json(c.reason)}};
    } else {
        est["correlated"] = nullptr;
    }

    json out;
    out["molecule"] = mol;
    out["basis"] = basis;
    out["estimate"] = est;
    out["notes"] = m.notes;
    out["warnings"] = m.warnings;
    out["refusal"] = m.refusal.empty() ? json(nullptr) : json(m.refusal);
    if (m.fci) {
        const FCIResultGeneral& f = *m.fci;
        json fc;
        fc["converged"] = f.converged;
        fc["iterations"] = f.iterations;
        fc["totalEnergy"] = f.total_energy;
        fc["electronicEnergy"] = f.electronic_energy;
        fc["correlationEnergy"] = f.correlation_energy;
        fc["S2"] = f.s2;
        fc["hfWeight"] = f.hf_weight;
        fc["orbitals"] = f.orbitals;
        fc["alphaElectrons"] = f.alpha;
        fc["betaElectrons"] = f.beta;
        fc["alphaStrings"] = f.alpha_strings;
        fc["betaStrings"] = f.beta_strings;
        fc["determinants"] = f.determinants;
        if (f.spin_adapted) {  // 15.0
            fc["csfs"] = f.csfs;
            fc["spin"] = f.spin;
        }
        json lead = json::array();
        for (const auto& d : f.leading)
            lead.push_back({{"coefficient", d.coefficient}, {"alpha", d.alpha}, {"beta", d.beta},
                            {"excitation", d.excitation}});
        fc["leadingDeterminants"] = lead;
        fc["naturalOccupations"] = f.natural_occupations;
        json hist = json::array();
        for (const auto& it : f.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"residual", it.residual},
                            {"subspace", it.subspace}});
        fc["history"] = hist;
        fc["transformSeconds"] = f.transform_seconds;
        fc["ciSeconds"] = f.ci_seconds;
        out["fullCI"] = fc;
    } else {
        out["fullCI"] = nullptr;
    }
    if (m.cc) {
        const CCResult& c = *m.cc;
        json cj;
        cj["method"] = c.method;
        cj["engine"] = c.engine;
        cj["converged"] = c.converged;
        cj["iterations"] = c.iterations;
        cj["occupied"] = c.occupied;
        cj["virtuals"] = c.virtuals;
        cj["referenceEnergy"] = c.reference_energy;
        cj["correlationEnergy"] = c.correlation_energy;
        cj["ccEnergy"] = c.cc_energy;
        cj["triplesCorrection"] = c.has_triples ? json(c.triples_correction) : json(nullptr);
        cj["totalEnergy"] = c.total_energy;
        cj["t1Diagnostic"] = c.t1_diagnostic;
        cj["amplitudeCounts"] = c.amplitude_counts;
        cj["amplitudeNorms"] = c.rank_norms;
        json lead = json::array();
        for (const auto& a : c.leading)
            lead.push_back({{"value", a.value},
                            {"rank", a.rank},
                            {"occupied", a.occupied},
                            {"virtuals", a.virtuals},
                            {"spins", a.spins}});
        cj["leadingAmplitudes"] = lead;
        json hist = json::array();
        for (const auto& it : c.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"correlationEnergy", it.correlation},
                            {"deltaE", it.delta_e},
                            {"residual", it.residual}});
        cj["history"] = hist;
        cj["ladder"] = c.ladder.empty() ? json(nullptr) : json(c.ladder);
        cj["transformSeconds"] = c.transform_seconds;
        cj["ccSeconds"] = c.cc_seconds;
        cj["triplesSeconds"] = c.has_triples ? json(c.triples_seconds) : json(nullptr);
        out["coupledCluster"] = cj;
    } else {
        out["coupledCluster"] = nullptr;
    }
    if (m.rhf) {
        const RHFResult& h = *m.rhf;
        json r;
        r["converged"] = h.converged;
        r["iterations"] = h.iterations;
        r["totalEnergy"] = h.total_energy;
        r["electronicEnergy"] = h.electronic_energy;
        r["nuclearRepulsion"] = h.nuclear_repulsion;
        r["kineticEnergy"] = h.kinetic_energy;
        r["nuclearAttractionEnergy"] = h.nuclear_attraction_energy;
        r["electronRepulsionEnergy"] = h.two_electron_energy;
        r["virialRatio"] = h.virial_ratio;
        r["orbitalEnergies"] = h.orbital_energies;
        r["occupiedOrbitals"] = h.n_occupied;
        const int no = h.n_occupied, nmo = static_cast<int>(h.orbital_energies.size());
        r["homo"] = no >= 1 ? json(h.orbital_energies[no - 1]) : json(nullptr);
        r["lumo"] = no < nmo ? json(h.orbital_energies[no]) : json(nullptr);
        r["smallestOverlapEigenvalue"] = h.smallest_overlap_eigenvalue;
        r["removedFunctions"] = h.removed_functions;
        r["integralsStored"] = h.eri_stored;
        r["uniqueIntegrals"] = h.unique_integrals;
        r["threads"] = h.threads;
        r["integralSeconds"] = h.integral_seconds;
        r["scfSeconds"] = h.scf_seconds;
        json hist = json::array();
        for (const auto& it : h.history)
            hist.push_back({{"iteration", it.iteration},
                            {"totalEnergy", it.energy},
                            {"deltaE", it.delta_e},
                            {"gradient", it.gradient},
                            {"diis", it.diis}});
        r["scfHistory"] = hist;
        // 15.0: the reference and the spin
        r["reference"] = h.reference;
        r["multiplicity"] = h.multiplicity;
        r["alphaElectrons"] = h.n_alpha;
        r["betaElectrons"] = h.n_beta;
        r["S2"] = h.s2;
        r["S2Exact"] = h.s2_exact;
        r["orbitalEnergiesBeta"] = h.reference == "UHF" ? json(h.orbital_energies_b) : json(nullptr);
        if (h.reference == "UHF") {
            const int nb = h.n_beta, nmb = static_cast<int>(h.orbital_energies_b.size());
            r["homoBeta"] = nb >= 1 ? json(h.orbital_energies_b[nb - 1]) : json(nullptr);
            r["lumoBeta"] = nb < nmb ? json(h.orbital_energies_b[nb]) : json(nullptr);
        }
        out["rhf"] = r;
    } else {
        out["rhf"] = nullptr;
    }
    out["reference"] = m.reference;
    if (m.so) {
        const SOResult& c = *m.so;
        json cj;
        cj["method"] = c.method;
        cj["label"] = method_label(c.method, c.reference);
        cj["reference"] = c.reference;
        cj["converged"] = c.converged;
        cj["iterations"] = c.iterations;
        cj["referenceEnergy"] = c.reference_energy;
        cj["correlationEnergy"] = c.correlation_energy;
        cj["totalEnergy"] = c.total_energy;
        cj["t1Diagnostic"] = c.t1_diagnostic;
        json terms = json::array();
        for (const auto& t : c.terms) terms.push_back({{"label", t.label}, {"value", t.value}});
        cj["terms"] = terms;
        json amps = json::array();
        for (const auto& a : c.amplitudes)
            amps.push_back({{"name", a.name}, {"rank", a.rank}, {"count", a.count}, {"norm", a.norm}});
        cj["amplitudes"] = amps;
        json hist = json::array();
        for (const auto& it : c.history)
            hist.push_back({{"iteration", it.iteration}, {"totalEnergy", it.energy}, {"deltaE", it.delta_e},
                            {"residual", it.residual}, {"stage", it.stage}});
        cj["history"] = hist;
        cj["seconds"] = c.seconds;
        out["correlated"] = cj;
    } else {
        out["correlated"] = nullptr;
    }
    return out;
}

json atoms_json(const Molecule& m) {
    json atoms = json::array();
    for (const Atom& a : m.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    return atoms;
}

double distance(const Atom& a, const Atom& b) {
    double s = 0.0;
    for (int k = 0; k < 3; ++k) s += (a.position[k] - b.position[k]) * (a.position[k] - b.position[k]);
    return std::sqrt(s);
}

// The parameters of the inputs: R (bohr; H2O: both O-H distances and their mean) and the
// H-O-H angle (degrees); H2O2: R (mean O-H), R1, R2, ROO, angle (mean O-O-H), angle1, angle2,
// dihedral; NH3: R (mean N-H), R1, R2, R3, angle (mean H-N-H), angle1, angle2, angle3 and
// "pyramidal" (360 minus the sum of the H-N-H angles, 0 for the planar molecule).
json geometry_parameters_json(const Molecule& m) {
    json g;
    g["atoms"] = atoms_json(m);
    if (m.custom && m.atoms.size() != 2) {
        // 9.0: a custom molecule is described by its Cartesian coordinates only
        g["R"] = nullptr;
        g["angle"] = nullptr;
        return g;
    }
    if (m.atoms.size() == 2) {
        g["R"] = distance(m.atoms[0], m.atoms[1]);
        g["angle"] = nullptr;
        return g;
    }
    if (m.atoms.size() == 4) {
        const double deg = 180.0 / 3.14159265358979323846;
        auto angle = [&](int a, int b, int c) {
            double u[3], v[3], uu = 0, vv = 0, uv = 0;
            for (int k = 0; k < 3; ++k) {
                u[k] = m.atoms[a].position[k] - m.atoms[b].position[k];
                v[k] = m.atoms[c].position[k] - m.atoms[b].position[k];
                uu += u[k] * u[k];
                vv += v[k] * v[k];
                uv += u[k] * v[k];
            }
            return std::acos(std::max(-1.0, std::min(1.0, uv / std::sqrt(uu * vv)))) * deg;
        };
        if (m.atoms[0].Z == 7) {
            const double r1 = distance(m.atoms[0], m.atoms[1]), r2 = distance(m.atoms[0], m.atoms[2]),
                         r3 = distance(m.atoms[0], m.atoms[3]);
            const double a1 = angle(1, 0, 2), a2 = angle(1, 0, 3), a3 = angle(2, 0, 3);
            g["R1"] = r1;
            g["R2"] = r2;
            g["R3"] = r3;
            g["R"] = (r1 + r2 + r3) / 3.0;
            g["angle1"] = a1;
            g["angle2"] = a2;
            g["angle3"] = a3;
            g["angle"] = (a1 + a2 + a3) / 3.0;
            g["pyramidal"] = 360.0 - (a1 + a2 + a3);
            return g;
        }
        const SimpleInternals si = simple_internals(m);
        g["ROO"] = si.values[0];
        g["R1"] = si.values[1];
        g["R2"] = si.values[2];
        g["R"] = 0.5 * (si.values[1] + si.values[2]);
        g["angle1"] = si.values[3] * deg;
        g["angle2"] = si.values[4] * deg;
        g["angle"] = 0.5 * (si.values[3] + si.values[4]) * deg;
        g["dihedral"] = si.values[5] * deg;
        return g;
    }
    int o = 0;
    for (std::size_t i = 0; i < m.atoms.size(); ++i)
        if (m.atoms[i].Z == 8) o = static_cast<int>(i);
    std::vector<int> h;
    for (std::size_t i = 0; i < m.atoms.size(); ++i)
        if (static_cast<int>(i) != o) h.push_back(static_cast<int>(i));
    if (h.size() != 2) {
        g["R"] = nullptr;
        g["angle"] = nullptr;
        return g;
    }
    const Atom &O = m.atoms[o], &H1 = m.atoms[h[0]], &H2 = m.atoms[h[1]];
    const double r1 = distance(O, H1), r2 = distance(O, H2);
    double dot = 0.0;
    for (int k = 0; k < 3; ++k) dot += (H1.position[k] - O.position[k]) * (H2.position[k] - O.position[k]);
    const double c = std::max(-1.0, std::min(1.0, dot / (r1 * r2)));
    g["R"] = 0.5 * (r1 + r2);
    g["R1"] = r1;
    g["R2"] = r2;
    g["angle"] = std::acos(c) * 180.0 / 3.14159265358979323846;
    return g;
}

json estimate_json(const JobEstimate& e) {
    json steps = json::array();
    for (const auto& s : e.steps)
        steps.push_back({{"name", s.name},
                         {"method", s.method},
                         {"basis", s.basis},
                         {"how", s.how},
                         {"analytic", s.analytic},
                         {"evaluations", s.evaluations},
                         {"evaluationKind", s.evaluation_kind},
                         {"seconds", s.seconds},
                         {"memoryMB", s.memory_mb}});
    return {{"steps", steps},
            {"totalSeconds", e.total_seconds},
            {"memoryMB", e.memory_mb},
            {"longCalculation", e.total_seconds > kLongCalculationSeconds},
            {"refusal", e.refusal.empty() ? json(nullptr) : json(e.refusal)},
            {"warnings", e.warnings},
            {"notes", e.notes}};
}

json point_json(const SinglePoint& p) {
    json j;
    j["method"] = p.method;
    j["energy"] = p.energy;
    j["converged"] = p.converged;
    j["rhfEnergy"] = p.rhf.total_energy;
    j["correlationEnergy"] = p.energy - p.rhf.total_energy;
    j["geometry"] = geometry_parameters_json(p.molecule);
    j["warnings"] = p.warnings;
    j["seconds"] = p.seconds;
    return j;
}

json optimization_json(const OptimizationResult& o) {
    json j;
    j["converged"] = o.converged;
    j["steps"] = o.steps;
    j["message"] = o.message;
    j["initialHessian"] = o.initial_hessian;
    j["transitionState"] = o.transition_state;
    j["tsHessian"] = o.transition_state ? json(o.ts_hessian) : json(nullptr);
    j["hessianEvaluations"] = o.hessian_evaluations;
    j["negativeEigenvalues"] = o.final_negative_eigenvalues >= 0 ? json(o.final_negative_eigenvalues) : json(nullptr);
    j["eigenTestFailed"] = o.eigen_test_failed;
    j["computedHessian"] = o.computed_hessian.empty() ? json(nullptr) : json(o.computed_hessian);
    j["recalcFC"] = o.recalc_every;
    j["algorithm"] = o.algorithm;
    j["gediisSteps"] = o.gediis_steps;
    j["criteria"] = {{"maxForce", o.criteria.max_force},
                     {"rmsForce", o.criteria.rms_force},
                     {"maxDisplacement", o.criteria.max_displacement},
                     {"rmsDisplacement", o.criteria.rms_displacement},
                     {"energy", o.criteria.energy}};
    j["initialEnergy"] = o.initial_energy;
    j["finalEnergy"] = o.final_energy;
    j["initialGeometry"] = geometry_parameters_json(o.initial);
    j["finalGeometry"] = geometry_parameters_json(o.final_geometry);
    auto params = [](const std::vector<GeometricParameter>& v) {
        json a = json::array();
        for (const auto& p : v) a.push_back({{"name", p.name}, {"value", p.value}, {"angle", p.angle}});
        return a;
    };
    j["initialParameters"] = params(o.initial_parameters);
    j["finalParameters"] = params(o.final_parameters);
    j["symmetricStart"] = o.symmetric_start;
    j["symmetricEnd"] = o.symmetric_end;
    j["symmetryDeviation"] = o.symmetry_deviation;
    j["gradient"] = o.history.empty() ? json(nullptr) : json(o.history.front().gradient_description);
    json hist = json::array();
    for (const auto& s : o.history) {
        json h;
        h["step"] = s.step;
        h["energy"] = s.energy;
        h["deltaE"] = s.step == 1 ? json(nullptr) : json(s.delta_e);
        h["geometry"] = geometry_parameters_json(s.geometry);
        json coords = json::array();
        for (std::size_t k = 0; k < s.names.size(); ++k)
            coords.push_back({{"name", s.names[k]},
                              {"definition", k < s.definitions.size() ? s.definitions[k] : ""},
                              {"value", k < s.q.size() ? s.q[k] : 0.0},
                              {"gradient", k < s.gq.size() ? s.gq[k] : 0.0},
                              {"step", k < s.dq.size() ? json(s.dq[k]) : json(nullptr)}});
        h["coordinates"] = coords;
        json conv = json::array();
        for (const auto& c : s.convergence)
            conv.push_back({{"name", c.name},
                            {"value", c.applicable ? json(c.value) : json(nullptr)},
                            {"threshold", c.threshold},
                            {"met", c.met},
                            {"applicable", c.applicable}});
        h["convergence"] = conv;
        h["converged"] = s.converged;
        h["hessianUpdate"] = s.hessian_update;
        h["trustRadius"] = s.trust_radius;
        h["stepNorm"] = s.step_norm;
        h["stepScaled"] = s.step_scaled;
        h["predictedChange"] = s.predicted_change;
        h["ratio"] = s.has_ratio ? json(s.ratio) : json(nullptr);
        h["warnings"] = s.warnings;
        h["seconds"] = s.seconds;
        if (s.ts_step) {
            h["ts"] = {{"hessianEigenvalues", s.hessian_eigenvalues},
                       {"negativeEigenvalues", s.negative_eigenvalues},
                       {"followedEigenvalue", s.followed_eigenvalue},
                       {"followedGradient", s.followed_gradient},
                       {"followedOverlap", s.followed_overlap},
                       {"followedMode", s.followed_mode},
                       {"lambdaP", s.lambda_p},
                       {"lambdaN", s.lambda_n},
                       {"hessianComputed", s.hessian_computed},
                       {"stepType", s.step_type},
                       {"gediisCoefficients", s.gediis_coefficients},
                       {"gediisSteps", s.gediis_steps}};
        } else {
            h["ts"] = nullptr;
        }
        hist.push_back(h);
    }
    j["history"] = hist;
    j["seconds"] = o.seconds;
    return j;
}

json frequencies_json(const HessianResult& hr, const HarmonicAnalysis& h) {
    json j;
    j["analytic"] = hr.analytic;
    j["description"] = hr.description;
    j["numericalReason"] = hr.numerical_reason.empty() ? json(nullptr) : json(hr.numerical_reason);
    j["gradientEvaluations"] = hr.gradient_evaluations;
    j["energyEvaluations"] = hr.energy_evaluations;
    j["coordinates"] = hr.coordinates;
    j["intensities"] = h.has_intensities;
    j["intensityDescription"] = hr.dipole_description;
    j["dipole"] = hr.has_dipole ? json({hr.dipole[0], hr.dipole[1], hr.dipole[2]}) : json(nullptr);
    j["linear"] = h.linear;
    j["pointGroup"] = h.point_group.empty() ? json(nullptr) : json(h.point_group);
    j["vibrations"] = h.n_vib;
    j["masses"] = h.masses;
    j["isotopes"] = h.isotopes;
    j["geometry"] = geometry_parameters_json(hr.reference.point.molecule);
    j["energy"] = hr.reference.point.energy;
    j["lowFrequencies"] = h.low_frequencies;
    json modes = json::array();
    for (const auto& m : h.modes) {
        json d = json::array();
        for (const auto& v : m.displacement) d.push_back({v[0], v[1], v[2]});
        modes.push_back({{"wavenumber", m.wavenumber},
                         {"imaginary", m.imaginary},
                         {"symmetry", m.symmetry},
                         {"reducedMass", m.reduced_mass},
                         {"forceConstant", m.force_constant},
                         {"irIntensity", m.ir_intensity >= 0.0 ? json(m.ir_intensity) : json(nullptr)},
                         {"displacement", d}});
    }
    j["modes"] = modes;
    j["seconds"] = hr.seconds;
    return j;
}

json thermo_json(const ThermoResult& t) {
    auto term = [](const ThermoTerm& x) {
        return json{{"energy", x.energy}, {"entropy", x.entropy}, {"cv", x.cv}, {"q", x.q}};
    };
    auto arr3 = [](const std::array<double, 3>& a) { return json{a[0], a[1], a[2]}; };
    json j;
    j["temperature"] = t.temperature;
    j["pressure"] = t.pressure;
    j["symmetryNumber"] = t.symmetry_number;
    j["linear"] = t.linear;
    j["mass"] = t.mass;
    j["moments"] = arr3(t.moments);
    j["rotationalConstantsGHz"] = arr3(t.rotational_constants_ghz);
    j["rotationalConstantsCm"] = arr3(t.rotational_constants_cm);
    j["rotationalTemperatures"] = arr3(t.rotational_temperatures);
    j["vibrationalTemperatures"] = t.vibrational_temperatures;
    j["excludedModes"] = t.excluded_modes;
    j["warnings"] = t.warnings;
    // J/mol, J/(mol K): energies of the terms include the zero-point energy (vibration).
    j["terms"] = {{"translation", term(t.translation)},
                  {"rotation", term(t.rotation)},
                  {"vibration", term(t.vibration)},
                  {"electronic", term(t.electronic)},
                  {"total", term(t.total)}};
    j["zeroPointEnergy"] = t.zero_point_energy;  // J/mol
    j["corrections"] = {{"zeroPoint", t.zero_point_energy},
                        {"thermalEnergy", t.thermal_energy},
                        {"enthalpy", t.enthalpy},
                        {"gibbs", t.gibbs},
                        {"helmholtz", t.helmholtz}};  // J/mol
    j["electronicEnergy"] = t.electronic_energy;      // E_h
    j["totals"] = {{"zeroPoint", t.total_zpe},
                   {"thermalEnergy", t.total_u},
                   {"enthalpy", t.total_h},
                   {"gibbs", t.total_g},
                   {"helmholtz", t.total_a}};  // E_h
    return j;
}

json job_json(const JobResult& r, const std::string& type) {
    json j;
    j["type"] = type;
    j["title"] = job_title(type);
    j["termination"] = r.termination;
    j["completed"] = r.completed;
    j["message"] = r.message.empty() ? json(nullptr) : json(r.message);
    j["estimate"] = estimate_json(r.estimate);
    j["warnings"] = r.warnings;
    j["seconds"] = r.seconds;
    if (r.termination == "estimate") {
        j["initialGeometry"] = nullptr;
        j["finalGeometry"] = nullptr;
        j["finalEnergy"] = nullptr;
        j["finalEnergyLabel"] = nullptr;
    } else {
        j["initialGeometry"] = geometry_parameters_json(r.initial_geometry);
        j["finalGeometry"] = geometry_parameters_json(r.final_geometry);
        j["finalEnergy"] = r.final_energy_label.empty() ? json(nullptr) : json(r.final_energy);
        j["finalEnergyLabel"] = r.final_energy_label.empty() ? json(nullptr) : json(r.final_energy_label);
    }
    j["singlePointMethod"] = r.single_point_method.empty() ? json(nullptr) : json(r.single_point_method);
    j["singlePointBasis"] = r.single_point_basis.empty() ? json(nullptr) : json(r.single_point_basis);
    return j;
}

json params_json(const std::vector<GeometricParameter>& v) {
    json a = json::array();
    for (const auto& p : v) a.push_back({{"name", p.name}, {"value", p.value}, {"angle", p.angle}});
    return a;
}

json scan_json(const ScanResult& sr) {
    const ScanPlan& pl = sr.plan;
    json j;
    json atoms = json::array();
    for (int a : pl.coordinate.coordinate.atoms) atoms.push_back(a + 1);
    // values in angstrom (bond) or degrees, as in the ModRedundant line
    const double f = pl.angular() ? 1.0 : kBohrInAngstrom;
    j["coordinate"] = {{"type", pl.coordinate.kind},
                       {"letter", pl.coordinate.letter},
                       {"atoms", atoms},
                       {"definition", pl.coordinate.definition},
                       {"unit", pl.angular() ? "degrees" : "angstrom"}};
    j["start"] = pl.start * f;
    j["stepSize"] = pl.step * f;
    j["steps"] = pl.steps;
    j["inputValue"] = pl.input_value * f;
    j["relaxed"] = sr.relaxed;
    j["converged"] = sr.converged;
    j["message"] = sr.message;
    j["lowest"] = sr.lowest;
    j["highest"] = sr.highest;
    const double emin = sr.lowest >= 0 ? sr.points[static_cast<std::size_t>(sr.lowest)].energy : 0.0;
    json pts = json::array();
    for (const auto& p : sr.points) {
        json steps = json::array();
        for (std::size_t k = 0; k < p.step_energies.size(); ++k)
            steps.push_back({{"step", k + 1}, {"energy", p.step_energies[k]}, {"maxForce", p.step_forces[k]}});
        pts.push_back({{"index", p.index},
                       {"value", p.value * f},
                       {"measured", p.measured * f},
                       {"energy", p.energy},
                       {"relativeEnergy", p.energy - emin},
                       {"converged", p.converged},
                       {"steps", p.steps},
                       {"maxForce", p.max_force},
                       {"geometry", geometry_parameters_json(p.geometry)},
                       {"parameters", params_json(p.parameters)},
                       {"optimizationSteps", steps},
                       {"warnings", p.warnings},
                       {"seconds", p.seconds}});
    }
    j["points"] = pts;
    json ex = json::array();
    for (const auto& x : sr.extrema)
        ex.push_back({{"kind", x.kind}, {"index", x.index}, {"value", x.value * f}, {"energy", x.energy}});
    j["extrema"] = ex;
    j["seconds"] = sr.seconds;
    return j;
}

// "job", "optimization", "frequencies", "thermochemistry", "singlePoint", "scan", "output".
void add_job_json(json& j, const Output& out) {
    if (!out.job) {
        j["job"] = nullptr;
        j["optimization"] = nullptr;
        j["frequencies"] = nullptr;
        j["thermochemistry"] = nullptr;
        j["singlePoint"] = nullptr;
        j["scan"] = nullptr;
        j["output"] = nullptr;
        return;
    }
    const JobResult& r = *out.job;
    const std::string type = canonical_job_type(out.input.job).empty() ? "energy" : canonical_job_type(out.input.job);
    j["job"] = job_json(r, type);
    j["optimization"] = r.optimization ? optimization_json(*r.optimization) : json(nullptr);
    j["frequencies"] = (r.hessian && r.harmonic) ? frequencies_json(*r.hessian, *r.harmonic) : json(nullptr);
    j["thermochemistry"] = r.thermo ? thermo_json(*r.thermo) : json(nullptr);
    j["singlePoint"] = r.single_point ? point_json(*r.single_point) : json(nullptr);
    j["scan"] = r.scan ? scan_json(*r.scan) : json(nullptr);
    j["output"] = r.termination == "estimate" ? json(nullptr) : json(r.output);
}

bool get_bool(const json& j, const char* key, bool fallback) {
    auto it = j.find(key);
    if (it == j.end() || it->is_null()) return fallback;
    if (!it->is_boolean()) throw std::runtime_error(std::string("\"") + key + "\" must be a boolean");
    return it->get<bool>();
}

// ΨLab 8.0: the "orbitals" object: atoms, basis shells (normalised contraction coefficients,
// spherical functions m = -l .. l from "first" on; 14.0: atoms with an "ecp" object (core
// electrons, local part, projected parts), "pure": false for a Cartesian shell,
// its (l+1)(l+2)/2 functions xx, xy, xz, yy, yz, zz ...; absent means spherical) and the orbital
// sets with their AO coefficients. The app passes it back to macqc_orbital_surface and
// macqc_orbital_plane.
json orbitals_json(const OrbitalData& d) {
    json j;
    j["basis"] = d.basis_name;
    j["functions"] = d.n_basis;
    json atoms = json::array();
    for (const Atom& a : d.atoms) {
        json at = {{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}};
        if (a.ecp) {
            // 14.0: the effective core potential (for Gaussian input; not needed to draw orbitals):
            // core electrons, L of the local part, its terms and those of U_l - U_L, [exponent, coefficient]
            auto terms = [](const std::vector<EcpTerm>& t) {
                json out = json::array();
                for (const EcpTerm& x : t) out.push_back({x.exponent, x.coefficient});
                return out;
            };
            json projected = json::array();
            for (const auto& p : a.ecp->projected) projected.push_back(terms(p));
            at["ecp"] = {{"core", a.ecp->core_electrons}, {"L", a.ecp->local_l}, {"local", terms(a.ecp->local)},
                         {"projected", projected}};
        }
        atoms.push_back(at);
    }
    j["atoms"] = atoms;
    json shells = json::array();
    for (const Shell& s : d.shells)
        shells.push_back({{"l", s.l},
                          {"atom", s.atom},
                          {"center", {s.center[0], s.center[1], s.center[2]}},
                          {"exponents", s.exponents},
                          {"coefficients", s.coefficients},
                          {"first", s.first_function},
                          {"pure", s.pure}});
    j["shells"] = shells;
    json sets = json::array();
    for (const OrbitalSet& set : d.sets) {
        json orbitals = json::array();
        for (const OrbitalInfo& o : set.orbitals)
            orbitals.push_back({{"number", o.number},
                                {"energy", o.has_energy ? json(o.energy) : json(nullptr)},
                                {"occupation", o.occupation},
                                {"label", o.label}});
        sets.push_back({{"kind", set.kind}, {"title", set.title}, {"orbitals", orbitals}, {"coefficients", set.coefficients}});
    }
    j["sets"] = sets;
    return j;
}

OrbitalData parse_orbitals_json(const std::string& text) {
    const json j = json::parse(text);
    OrbitalData d;
    d.basis_name = j.value("basis", std::string());
    d.n_basis = j.at("functions").get<int>();
    for (const auto& a : j.at("atoms")) {
        Atom atom;
        atom.symbol = a.at("symbol").get<std::string>();
        atom.Z = a.at("Z").get<int>();
        for (int k = 0; k < 3; ++k) atom.position[k] = a.at("position").at(k).get<double>();
        d.atoms.push_back(atom);
    }
    for (const auto& s : j.at("shells")) {
        Shell shell;
        shell.l = s.at("l").get<int>();
        if (shell.l < 0 || shell.l > kMaxAngularMomentum) throw std::runtime_error("orbitals: angular momentum out of range");
        shell.atom = s.value("atom", 0);
        for (int k = 0; k < 3; ++k) shell.center[k] = s.at("center").at(k).get<double>();
        shell.exponents = s.at("exponents").get<std::vector<double>>();
        shell.coefficients = s.at("coefficients").get<std::vector<double>>();
        shell.first_function = s.at("first").get<int>();
        shell.pure = s.value("pure", true);  // archives of 13.0 have spherical shells only
        if (shell.exponents.size() != shell.coefficients.size() || shell.first_function < 0 ||
            shell.first_function + shell.size() > d.n_basis)
            throw std::runtime_error("orbitals: inconsistent shell");
        d.shells.push_back(shell);
    }
    for (const auto& s : j.at("sets")) {
        OrbitalSet set;
        set.kind = s.value("kind", std::string());
        set.title = s.value("title", std::string());
        for (const auto& o : s.at("orbitals")) {
            OrbitalInfo info;
            info.number = o.value("number", 0);
            info.has_energy = o.contains("energy") && o.at("energy").is_number();
            info.energy = info.has_energy ? o.at("energy").get<double>() : 0.0;
            info.occupation = o.value("occupation", 0.0);
            info.label = o.value("label", std::string());
            set.orbitals.push_back(info);
        }
        set.coefficients = s.at("coefficients").get<std::vector<std::vector<double>>>();
        for (const auto& c : set.coefficients)
            if (c.size() != static_cast<std::size_t>(d.n_basis))
                throw std::runtime_error("orbitals: the coefficients do not match the basis");
        d.sets.push_back(set);
    }
    return d;
}

// The orbital chosen by set and orbital index, or an exception with a message.
const std::vector<double>& chosen_orbital(const OrbitalData& d, int set, int orbital) {
    if (set < 0 || set >= static_cast<int>(d.sets.size())) throw std::runtime_error("orbitals: no such orbital set");
    const OrbitalSet& s = d.sets[static_cast<std::size_t>(set)];
    if (orbital < 0 || orbital >= static_cast<int>(s.coefficients.size()))
        throw std::runtime_error("orbitals: no such orbital");
    return s.coefficients[static_cast<std::size_t>(orbital)];
}

}  // namespace

std::string orbitals_to_json(const OrbitalData& d) { return orbitals_json(d).dump(); }

OrbitalData orbitals_from_json(const std::string& text) { return parse_orbitals_json(text); }

Input parse_input_json(const std::string& text) {
    json j;
    try {
        j = json::parse(text);
    } catch (const std::exception& e) {
        throw std::runtime_error(std::string("invalid JSON input: ") + e.what());
    }
    if (!j.is_object()) throw std::runtime_error("input JSON must be an object");
    Input in;
    in.molecule = get_string(j, "molecule", in.molecule);
    std::string canonical = canonical_molecule_name(in.molecule);
    const MoleculePreset* preset = canonical.empty() ? find_molecule_preset(in.molecule) : nullptr;
    if (canonical.empty() && !preset)
        throw std::runtime_error("\"molecule\" must be HeH+, H2, LiH, H2O, H2O2, NH3, an atom H–Ar, H2+, H3+, H5+, N2, O2, "
                                 "F2, Cl2, NO, NO2 or custom");
    in.molecule = preset ? preset->key : canonical;
    // 15.0: charge and multiplicity of every molecule, the reference
    if (auto it = j.find("charge"); it != j.end() && !it->is_null()) {
        in.charge = get_int(j, "charge", 0);
        in.chargeGiven = true;
    }
    if (auto it = j.find("multiplicity"); it != j.end() && !it->is_null()) {
        in.multiplicity = get_int(j, "multiplicity", 1);
        in.multiplicityGiven = true;
    }
    in.reference = get_string(j, "reference", "");
    if (!in.reference.empty() && canonical_reference_name(in.reference).empty())
        throw std::runtime_error("\"reference\" must be RHF, UHF or ROHF");
    in.szaboOstlund = get_bool(j, "szaboOstlund", true);
    if (canonical != "HeH+") in.R = 0.0;  // default bond length of the molecule unless given
    if (canonical == "HeH+") {
        // 7.0: a basis set of the library (correlation-consistent; 14.0 also Pople and def2)
        // selects the general path; STO-1G/2G/3G (or no basis) the Szabo-Ostlund calculation
        const std::string b = get_string(j, "basis", "");
        const std::string cb = canonical_basis_name(b);
        const bool sto = cb == "STO-1G" || cb == "STO-2G" || cb == "STO-3G";
        if (!b.empty() && !cb.empty() && (!sto || !in.szaboOstlund)) {
            in.hehGeneral = true;  // 15.0: STO-nG for HeH+ also in the general calculation (szaboOstlund false)
        } else if (!b.empty()) {
            std::string t;
            for (char c : b)
                if (c != '-' && c != ' ') t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (t == "sto1g" || t == "sto2g" || t == "sto3g") in.nGaussians = t[3] - '0';
            else throw std::runtime_error("\"basis\" for HeH+ must be STO-1G, STO-2G, STO-3G or a basis set of the library "
                                          "(6-31G ..., 6-311G ..., (aug-)cc-pVXZ, (aug-)cc-pCVXZ, def2-...)");
        }
    }
    if (canonical != "HeH+" && canonical != "custom") in.angle = default_angle(canonical);
    if (canonical == "custom" || (preset && j.find("atoms") != j.end())) {
        // 9.0: "atoms": [{"symbol": "C", "x": 0, "y": 0, "z": 1.16}, ...] in angstrom (or with
        // "position": [x, y, z]; "Z" instead of "symbol"), "charge", "multiplicity", "title"
        auto it = j.find("atoms");
        if (it == j.end() || !it->is_array() || it->empty())
            throw std::runtime_error("a custom molecule needs \"atoms\": [{\"symbol\", \"x\", \"y\", \"z\"}, ...] in angstrom");
        int k = 0;
        for (const json& a : *it) {
            ++k;
            const std::string where = "atom " + std::to_string(k) + ": ";
            if (!a.is_object()) throw std::runtime_error(where + "must be an object");
            Atom at;
            const std::string sym = get_string(a, "symbol", "");
            at.Z = sym.empty() ? get_int(a, "Z", 0) : atomic_number(sym);
            if (at.Z < 1 || at.Z > kMaxAtomicNumber)
                throw std::runtime_error(where + "unknown element" + (sym.empty() ? std::string() : " '" + sym + "'") +
                                         " (Z = 1 to 118)");
            at.symbol = element_symbol(at.Z);
            std::array<double, 3> x{};
            if (auto p = a.find("position"); p != a.end() && !p->is_null()) {
                if (!p->is_array() || p->size() != 3) throw std::runtime_error(where + "\"position\" must be [x, y, z]");
                for (std::size_t c = 0; c < 3; ++c) {
                    if (!(*p)[c].is_number()) throw std::runtime_error(where + "the coordinates must be numbers");
                    x[c] = (*p)[c].get<double>();
                }
            } else {
                const char* keys[3] = {"x", "y", "z"};
                for (std::size_t c = 0; c < 3; ++c) {
                    if (a.find(keys[c]) == a.end()) throw std::runtime_error(where + "\"x\", \"y\" and \"z\" are needed (angstrom)");
                    x[c] = get_number(a, keys[c], 0.0);
                }
            }
            for (std::size_t c = 0; c < 3; ++c) {
                if (!std::isfinite(x[c])) throw std::runtime_error(where + "the coordinates must be finite numbers");
                at.position[c] = x[c] / kBohrInAngstrom;
            }
            in.atoms.push_back(at);
        }
        in.charge = get_int(j, "charge", preset ? preset->charge : 0);
        in.multiplicity = get_int(j, "multiplicity", preset ? preset->multiplicity : 1);
        in.title = get_string(j, "title", "");
        in.R = 0.0;
    }
    if (preset) in.R = 0.0;
    in.R = get_number(j, "R", in.R);
    in.zeta1 = get_number(j, "zeta1", in.zeta1);
    in.zeta2 = get_number(j, "zeta2", in.zeta2);
    in.ZA = get_number(j, "ZA", in.ZA);
    in.ZB = get_number(j, "ZB", in.ZB);
    in.nGaussians = get_int(j, "nGaussians", in.nGaussians);
    in.exactErf = get_bool(j, "exactErf", in.exactErf);
    in.maxIter = get_int(j, "maxIter", in.maxIter);
    in.convergence = get_number(j, "convergence", in.convergence);
    in.printLevel = get_int(j, "printLevel", in.printLevel);
    if (!is_heh_input(in)) in.basis = get_string(j, "basis", in.basis);
    in.angle = get_number(j, "angle", in.angle);
    in.ROO = get_number(j, "ROO", in.ROO);
    in.dihedral = get_number(j, "dihedral", in.dihedral);
    in.method = get_string(j, "method", in.method);
    in.basisDirectory = get_string(j, "basisDirectory", in.basisDirectory);
    in.threads = get_int(j, "threads", in.threads);
    in.memoryMB = get_number(j, "memoryMB", in.memoryMB);
    in.fciMemoryMB = get_number(j, "fciMemoryMB", in.fciMemoryMB);
    in.scfMaxIter = get_int(j, "scfMaxIter", in.scfMaxIter);
    in.energyConvergence = get_number(j, "energyConvergence", in.energyConvergence);
    in.gradientConvergence = get_number(j, "gradientConvergence", in.gradientConvergence);
    in.estimateOnly = get_bool(j, "estimateOnly", in.estimateOnly);
    in.job = get_string(j, "job", in.job);
    in.spMethod = get_string(j, "spMethod", in.spMethod);
    in.spBasis = get_string(j, "spBasis", in.spBasis);
    in.temperature = get_number(j, "temperature", in.temperature);
    in.pressure = get_number(j, "pressure", in.pressure);
    in.optConvergence = get_string(j, "optConvergence", in.optConvergence);
    in.maxSteps = get_int(j, "maxSteps", in.maxSteps);
    in.scan = get_string(j, "scan", in.scan);
    in.tsHessian = get_string(j, "tsHessian", in.tsHessian);
    for (char& c : in.tsHessian) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    in.eigenTest = get_bool(j, "eigenTest", in.eigenTest);
    in.hessian = get_string(j, "hessian", in.hessian);
    in.initialHessian = get_string(j, "initialHessian", in.initialHessian);
    for (char& c : in.initialHessian) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    in.optHessian = get_string(j, "optHessian", in.optHessian);
    for (char& c : in.optHessian) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    in.optAlgorithm = get_string(j, "optAlgorithm", in.optAlgorithm);
    for (char& c : in.optAlgorithm) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (j.contains("recalcFC") && j["recalcFC"].is_number_integer()) in.recalcFC = j["recalcFC"].get<int>();
    for (char& c : in.hessian) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (auto it = j.find("scanStart"); it != j.end() && !it->is_null()) in.scanStart = get_number(j, "scanStart", 0.0);
    return in;
}

namespace {

json output_json(const Output& out, const OrbitalData* orbitals = nullptr) {
    json j;
    j["input"] = input_json(out.input);
    j["error"] = out.error.empty() ? json(nullptr) : json(out.error);
    if (!is_heh_input(out.input) && general_molecule_name(out.input.molecule)) {
        // General molecules: {"input","error","molecule","basis","estimate","notes",
        // "warnings","refusal","rhf","fullCI","coupledCluster","job","optimization",
        // "frequencies","thermochemistry","singlePoint","output","hf":null,"fci":null,"report"}.
        // rhf/fullCI/coupledCluster/report describe the final energy calculation of the job.
        if (out.molecular) {
            const json m = molecular_json(*out.molecular);
            for (auto it = m.begin(); it != m.end(); ++it) j[it.key()] = it.value();
        } else {
            j["molecule"] = nullptr;
            j["basis"] = nullptr;
            j["estimate"] = nullptr;
            j["notes"] = json::array();
            j["warnings"] = json::array();
            j["refusal"] = nullptr;
            j["rhf"] = nullptr;
            j["fullCI"] = nullptr;
            j["coupledCluster"] = nullptr;
        }
        add_job_json(j, out);
        j["orbitals"] = nullptr;
        if (orbitals) {
            j["orbitals"] = orbitals_json(*orbitals);  // 13.0: built once for the appendix
        } else if (out.molecular && out.molecular->rhf && out.molecular->rhf->C.rows() > 0) {
            try {
                const MolecularResult& m = *out.molecular;
                const OrbitalData d = make_orbital_data(m.molecule, m.basis, *m.rhf, m.fci ? &*m.fci : nullptr);
                j["orbitals"] = orbitals_json(d);
            } catch (const std::exception&) {
                j["orbitals"] = nullptr;
            }
        }
        j["hf"] = nullptr;
        j["fci"] = nullptr;
        j["report"] = out.report;
        return j;
    }

    if (out.hf) {
        const HFResult& h = *out.hf;
        json hf;
        hf["converged"] = h.converged;
        hf["iterations"] = h.iterations;
        hf["electronicEnergy"] = h.electronic_energy;
        hf["nuclearRepulsion"] = h.nuclear_repulsion;
        hf["totalEnergy"] = h.total_energy;
        hf["orbitalEnergies"] = h.orbital_energies;
        hf["C"] = matrix_json(h.C);
        hf["P"] = matrix_json(h.P);
        hf["S"] = matrix_json(h.S);
        hf["X"] = matrix_json(h.X);
        hf["H"] = matrix_json(h.H);
        hf["F"] = matrix_json(h.F);
        hf["G"] = matrix_json(h.G);
        hf["PS"] = matrix_json(h.PS);
        const Integrals& v = h.integrals;
        hf["integrals"] = {{"S12", v.S12},     {"T11", v.T11},     {"T12", v.T12},   {"T22", v.T22},
                           {"V11A", v.V11A},   {"V12A", v.V12A},   {"V22A", v.V22A}, {"V11B", v.V11B},
                           {"V12B", v.V12B},   {"V22B", v.V22B},   {"V1111", v.V1111}, {"V2111", v.V2111},
                           {"V2121", v.V2121}, {"V2211", v.V2211}, {"V2221", v.V2221}, {"V2222", v.V2222}};
        json hist = json::array();
        for (const auto& s : h.history)
            hist.push_back({{"iteration", s.iteration},
                            {"electronicEnergy", s.electronic_energy},
                            {"deltaP", s.delta_p}});
        hf["scfHistory"] = hist;
        j["hf"] = hf;
        // 8.0: the two orbitals, He (atom A) at the origin and H (atom B) at z = R
        try {
            const Input& in = out.input;
            const std::vector<double> C = {h.C(0, 0), h.C(0, 1), h.C(1, 0), h.C(1, 1)};
            j["orbitals"] = orbitals_json(
                make_sto_orbital_data(in.nGaussians, in.R, in.zeta1, in.zeta2, in.ZA, in.ZB, C, h.orbital_energies));
        } catch (const std::exception&) {
            j["orbitals"] = nullptr;
        }
    } else {
        j["hf"] = nullptr;
        j["orbitals"] = nullptr;
    }

    if (out.fci) {
        const FCIResult& f = *out.fci;
        json fci;
        fci["totalEnergy"] = f.total_energy;
        fci["electronicEnergy"] = f.electronic_energy;
        fci["correlationEnergy"] = f.correlation_energy;
        json states = json::array();
        for (std::size_t k = 0; k < f.states.size(); ++k)
            states.push_back({{"energy", f.states[k].energy},
                              {"S2", f.states[k].S2},
                              {"spin", f.states[k].spin},
                              {"vector", f.vectors[k]}});
        fci["states"] = states;
        fci["groundStateVector"] = f.ground_state_vector;
        fci["determinants"] = f.determinants;
        json w = json::array();
        for (const auto& s : f.singlet_weights)
            w.push_back({{"configuration", s.configuration}, {"weight", s.weight}});
        fci["singletWeights"] = w;
        j["fci"] = fci;
    } else {
        j["fci"] = nullptr;
    }

    j["report"] = out.report;
    return j;
}


// ---- 13.0: the appendix of the output file ------------------------------------------------

// Normalisation of the primitive x^l exp(-a r^2), as the basis set applies it (a contraction of
// one primitive with coefficient 1 is that primitive, normalised).
double primitive_norm(int l, double a) { return normalized_coefficients(l, {a}, {1.0})[0]; }

std::string atom_label(const OrbitalData& d, int atom) {
    const Atom& a = d.atoms[static_cast<std::size_t>(atom)];
    return (a.symbol.empty() ? element_symbol(a.Z) : a.symbol) + std::to_string(atom + 1);
}

// x^i y^j z^k as "xxy" (14.0, labels of Cartesian functions); the components of l in basis order.
std::vector<std::string> cartesian_component_names(int l) {
    std::vector<std::string> names;
    for (int x = l; x >= 0; --x)
        for (int y = l - x; y >= 0; --y)
            names.push_back(std::string(static_cast<std::size_t>(x), 'x') + std::string(static_cast<std::size_t>(y), 'y') +
                            std::string(static_cast<std::size_t>(l - x - y), 'z'));
    return names;
}

// Labels of the basis functions in basis order: atom, shell ("2p": the second p shell of the
// atom) and m (real solid harmonics m = -l ... l), e.g. "O1 2p-1"; for a Cartesian shell
// (14.0) the component, e.g. "O1 3dxy".
std::vector<std::string> basis_function_labels(const OrbitalData& d) {
    std::vector<std::string> labels(static_cast<std::size_t>(d.n_basis));
    std::vector<std::vector<int>> count(d.atoms.size(), std::vector<int>(kMaxAngularMomentum + 1, 0));
    for (const Shell& sh : d.shells) {
        const int n = sh.l + 1 + count[static_cast<std::size_t>(sh.atom)][static_cast<std::size_t>(sh.l)]++;
        const std::string head = strprintf("%-5s %d%c", atom_label(d, sh.atom).c_str(), n, angular_letter(sh.l));
        if (!sh.pure) {
            const auto names = cartesian_component_names(sh.l);
            for (std::size_t c = 0; c < names.size(); ++c) {
                const std::size_t f = static_cast<std::size_t>(sh.first_function) + c;
                if (f < labels.size()) labels[f] = head + names[c];
            }
            continue;
        }
        for (int m = -sh.l; m <= sh.l; ++m) {
            const std::size_t f = static_cast<std::size_t>(sh.first_function + m + sh.l);
            if (f < labels.size())
                labels[f] = head + (sh.l == 0 ? "" : (m == 0 ? "0" : strprintf("%+d", m)));
        }
    }
    return labels;
}

// The basis set in full: exponents (bohr^-2) and the contraction coefficients of normalised
// primitives (every contracted function normalised); atoms with the same functions as an
// earlier atom of their element are listed once.
std::string basis_set_text(const OrbitalData& d) {
    OutputWriter w;
    w.section("Basis set in full");
    bool cartesian = false;
    for (const Shell& sh : d.shells) cartesian = cartesian || !sh.pure;
    if (!cartesian)
        w.line(strprintf("    %s: %d spherical basis functions in %zu shells (real solid harmonics, m = -l ... l)",
                         d.basis_name.c_str(), d.n_basis, d.shells.size()));
    else
        w.line(strprintf("    %s: %d basis functions in %zu shells (Cartesian d: xx, xy, xz, yy, yz, zz, each normalised;\n"
                         "    the other shells real solid harmonics, m = -l ... l)",
                         d.basis_name.c_str(), d.n_basis, d.shells.size()));
    w.line("    exponents in bohr^-2; contraction coefficients of normalised primitives, every contracted");
    w.line("    function normalised; the functions are numbered as in the molecular orbitals below");
    std::vector<std::vector<const Shell*>> by_atom(d.atoms.size());
    for (const Shell& sh : d.shells) by_atom[static_cast<std::size_t>(sh.atom)].push_back(&sh);
    for (std::size_t a = 0; a < d.atoms.size(); ++a) {
        // same element and the same shells (exponents, l) as an earlier atom: name it
        int same = -1;
        for (std::size_t b = 0; b < a && same < 0; ++b) {
            if (d.atoms[b].Z != d.atoms[a].Z || by_atom[b].size() != by_atom[a].size()) continue;
            bool equal = true;
            for (std::size_t k = 0; k < by_atom[a].size() && equal; ++k)
                equal = by_atom[a][k]->l == by_atom[b][k]->l && by_atom[a][k]->pure == by_atom[b][k]->pure &&
                        by_atom[a][k]->exponents == by_atom[b][k]->exponents &&
                        by_atom[a][k]->coefficients == by_atom[b][k]->coefficients;
            if (equal) same = static_cast<int>(b);
        }
        const int first = by_atom[a].empty() ? 0 : by_atom[a].front()->first_function + 1;
        int last = first - 1;
        for (const Shell* sh : by_atom[a]) last = sh->first_function + sh->size();
        if (same >= 0) {
            w.line(strprintf("\n    %-5s functions %d-%d: the same as %s", atom_label(d, static_cast<int>(a)).c_str(),
                             first, last, atom_label(d, same).c_str()));
            continue;
        }
        w.line(strprintf("\n    %-5s functions %d-%d", atom_label(d, static_cast<int>(a)).c_str(), first, last));
        w.line("      shell  functions          exponent        coefficient");
        std::vector<int> count(kMaxAngularMomentum + 1, 0);
        for (const Shell* sh : by_atom[a]) {
            const int n = sh->l + 1 + count[static_cast<std::size_t>(sh->l)]++;
            const std::string name = strprintf("%d%c", n, angular_letter(sh->l));
            const std::string range = sh->size() == 1 ? strprintf("%d", sh->first_function + 1)
                                                      : strprintf("%d-%d", sh->first_function + 1,
                                                                  sh->first_function + sh->size());
            for (std::size_t k = 0; k < sh->exponents.size(); ++k) {
                const double c = sh->coefficients[k] / primitive_norm(sh->l, sh->exponents[k]);
                w.line(strprintf("      %-5s  %-9s %18.10f %18.10f", k == 0 ? name.c_str() : "",
                                 k == 0 ? range.c_str() : "", sh->exponents[k], c));
            }
        }
        // 14.0: the effective core potential of the atom
        if (const auto& ecp = d.atoms[a].ecp) {
            const char L = angular_letter(ecp->local_l);
            w.line(strprintf("      ECP%d: %d core electrons; U(r) = U_%c(r) + sum_l |lm> [U_l(r) - U_%c(r)] <lm|, "
                             "terms c exp(-a r^2)",
                             ecp->core_electrons, ecp->core_electrons, L, L));
            w.line("      part                     exponent        coefficient");
            auto part = [&](const std::string& label, const std::vector<EcpTerm>& terms) {
                for (std::size_t k = 0; k < terms.size(); ++k)
                    w.line(strprintf("      %-12s %18.10f %18.10f", k == 0 ? label.c_str() : "", terms[k].exponent,
                                     terms[k].coefficient));
            };
            part(std::string("U_") + L, ecp->local);
            for (std::size_t l = 0; l < ecp->projected.size(); ++l)
                part(std::string("U_") + angular_letter(static_cast<int>(l)) + " - U_" + L, ecp->projected[l]);
        }
    }
    return w.text();
}

// The orbitals of the Orbitals tab with their AO coefficients, five per block.
std::string orbitals_text(const OrbitalData& d) {
    OutputWriter w;
    w.section("Molecular orbitals (the orbitals of the Orbitals tab)");
    w.line("    AO coefficients of the basis functions of the basis set above");
    const std::vector<std::string> labels = basis_function_labels(d);
    for (const OrbitalSet& set : d.sets) {
        w.line(strprintf("\n  %s (%zu orbitals)", set.title.c_str(), set.orbitals.size()));
        for (std::size_t first = 0; first < set.orbitals.size(); first += 5) {
            const std::size_t last = std::min(first + 5, set.orbitals.size());
            std::string number = "\n                     ", label = "                     ",
                        energy = "    Energy (E_h)     ", occupation = "    Occupation       ";
            bool any_energy = false;
            for (std::size_t k = first; k < last; ++k) {
                const OrbitalInfo& o = set.orbitals[k];
                number += strprintf("%12d", o.number);
                label += strprintf("%12s", o.label.c_str());
                energy += o.has_energy ? strprintf("%12.6f", o.energy) : std::string(12, ' ');
                occupation += strprintf("%12.6f", o.occupation);
                any_energy = any_energy || o.has_energy;
            }
            w.line(number);
            if (label.find_first_not_of(' ') != std::string::npos) w.line(label);
            if (any_energy) w.line(energy);
            w.line(occupation);
            for (int f = 0; f < d.n_basis; ++f) {
                std::string row = strprintf("    %4d  %-13s", f + 1, labels[static_cast<std::size_t>(f)].c_str());
                for (std::size_t k = first; k < last; ++k) {
                    const std::vector<double>& c = set.coefficients[k];
                    row += strprintf("%12.6f", static_cast<std::size_t>(f) < c.size() ? c[static_cast<std::size_t>(f)] : 0.0);
                }
                w.line(row);
            }
        }
    }
    return w.text();
}

constexpr const char* kArchiveBegin = "@@PSILAB-ARCHIVE-BEGIN";
constexpr const char* kArchiveEnd = "@@PSILAB-ARCHIVE-END";

// The archive: the results as ASCII JSON (non-ASCII characters escaped), 100 columns per line.
std::string archive_text(json j) {
    j.erase("output");
    j["archive"] = {{"format", 1},
                    {"program", std::string("\xCE\xA8Lab ") + kProgramVersion},
                    {"written", local_time_string()}};
    const std::string data = j.dump(-1, ' ', true, json::error_handler_t::replace);
    OutputWriter w;
    w.section("\xCE\xA8Lab archive: the complete results, read by File > Open Output File\xE2\x80\xA6 of the app");
    w.line("    JSON (format 1); the lines between the markers are joined without line breaks");
    w.line(std::string(kArchiveBegin) + " 1");
    for (std::size_t k = 0; k < data.size(); k += 100) w.line(data.substr(k, 100));
    w.line(kArchiveEnd);
    return w.text();
}

}  // namespace

namespace {
// 16.0: JSON has no infinity: non-finite numbers (the estimates of impossible calculations, e.g.
// Full CI determinant counts beyond 1e308) become the largest finite double instead of null.
void finite_numbers(json& j) {
    if (j.is_number_float()) {
        const double v = j.get<double>();
        if (!std::isfinite(v)) j = std::isnan(v) ? 0.0 : (v > 0 ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max());
    } else if (j.is_structured()) {
        for (auto& x : j) finite_numbers(x);
    }
}
}  // namespace

std::string output_to_json(const Output& out, int indent) {
    json j = output_json(out);
    finite_numbers(j);
    return j.dump(indent, ' ', false, json::error_handler_t::replace);
}

std::string output_appendix(const Output& out) {
    std::string text;
    std::optional<OrbitalData> d;
    if (out.molecular && out.molecular->rhf && out.molecular->rhf->C.rows() > 0) {
        try {
            const MolecularResult& m = *out.molecular;
            d = make_orbital_data(m.molecule, m.basis, *m.rhf, m.fci ? &*m.fci : nullptr);
            text += basis_set_text(*d);
            text += orbitals_text(*d);
        } catch (const std::exception& e) {
            text += std::string("\n    (the molecular orbitals could not be written: ") + e.what() + ")\n";
        }
    }
    text += archive_text(output_json(out, d ? &*d : nullptr));
    return text;
}

std::string archive_from_output_text(const std::string& text) {
    const std::size_t begin = text.rfind(kArchiveBegin);
    if (begin == std::string::npos) return "";
    std::size_t pos = text.find('\n', begin);
    std::string data;
    while (pos != std::string::npos && pos + 1 < text.size()) {
        const std::size_t next = text.find('\n', pos + 1);
        std::string line = text.substr(pos + 1, next == std::string::npos ? std::string::npos : next - pos - 1);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.rfind(kArchiveEnd, 0) == 0) return data;
        data += line;
        pos = next;
    }
    return "";  // no end marker: incomplete
}

}  // namespace macqc

struct macqc_job {
    macqc::Job job;
};

extern "C" {

char* macqc_run(const char* input_json) { return macqc_run_job(input_json, nullptr); }

macqc_job* macqc_job_create(void) {
    try {
        return new macqc_job();
    } catch (...) {
        return nullptr;
    }
}

void macqc_job_cancel(macqc_job* job) {
    if (job) job->job.cancelled.store(true);
}

void macqc_job_progress(const macqc_job* job, int* stage, int* iteration, double* energy, double* residual) {
    if (!job) return;
    if (stage) *stage = job->job.stage.load(std::memory_order_relaxed);
    if (iteration) *iteration = job->job.iteration.load(std::memory_order_relaxed);
    if (energy) *energy = job->job.energy.load(std::memory_order_relaxed);
    if (residual) *residual = job->job.residual.load(std::memory_order_relaxed);
}

void macqc_job_free(macqc_job* job) { delete job; }

namespace {
char* copy_string(const std::string& s) {
    char* buf = static_cast<char*>(std::malloc(s.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, s.c_str(), s.size() + 1);
    return buf;
}
}  // namespace

char* macqc_job_output(const macqc_job* job, long long offset, long long* total) {
    if (total) *total = 0;
    if (!job) return copy_string("");
    try {
        std::size_t t = 0;
        const std::string s = job->job.output_since(offset < 0 ? 0 : static_cast<std::size_t>(offset), &t);
        if (total) *total = static_cast<long long>(t);
        return copy_string(s);
    } catch (...) {
        return nullptr;
    }
}

char* macqc_job_progress_text(const macqc_job* job) {
    if (!job) return copy_string("");
    try {
        return copy_string(job->job.progress_text());
    } catch (...) {
        return nullptr;
    }
}

}  // extern "C"

// ---- ΨLab 20.0: nuclear motion (JSON in and out)

namespace macqc {

NuclearMotionInput parse_nuclear_input_json(const std::string& text) {
    using nlohmann::json;
    json j = json::parse(text.empty() ? std::string("{}") : text);
    if (!j.is_object()) throw std::invalid_argument("the input must be a JSON object");
    NuclearMotionInput in;
    if (j.contains("isotopologue") && j["isotopologue"].is_string())
        in.isotopologue = isotopologue_from_name(j["isotopologue"].get<std::string>());
    in.J = get_int(j, "J", in.J);
    in.parity_block = get_int(j, "parityBlock", in.parity_block);
    in.radial_points = get_int(j, "radialPoints", in.radial_points);
    in.angular_points = get_int(j, "angularPoints", in.angular_points);
    in.r_centre = get_number(j, "rCentre", in.r_centre);
    in.r_half_width = get_number(j, "rHalfWidth", in.r_half_width);
    in.levels = get_int(j, "levels", in.levels);
    in.eigensolver = get_string(j, "eigensolver", in.eigensolver);
    in.linear_solver = get_string(j, "linearSolver", in.linear_solver);
    in.preconditioner = get_string(j, "preconditioner", in.preconditioner);
    in.eigenvalue_tolerance = get_number(j, "eigenvalueTolerance", in.eigenvalue_tolerance);
    in.residual_tolerance = get_number(j, "residualTolerance", in.residual_tolerance);
    in.linear_tolerance = get_number(j, "linearTolerance", in.linear_tolerance);
    in.max_iterations = get_int(j, "maxIterations", in.max_iterations);
    in.max_linear_iterations = get_int(j, "maxLinearIterations", in.max_linear_iterations);
    in.jd_inner_iterations = get_int(j, "jdInnerIterations", in.jd_inner_iterations);
    in.subspace_size = get_int(j, "subspaceSize", in.subspace_size);
    in.expectation_values = get_bool(j, "expectationValues", in.expectation_values);
    in.threads = get_int(j, "threads", in.threads);
    return in;
}

nlohmann::json nuclear_input_json(const NuclearMotionInput& in) {
    using nlohmann::json;
    return json{{"isotopologue", isotopologue_name(in.isotopologue)},
                {"J", in.J},
                {"parityBlock", in.parity_block},
                {"radialPoints", in.radial_points},
                {"angularPoints", in.angular_points},
                {"rCentre", in.r_centre},
                {"rHalfWidth", in.r_half_width},
                {"levels", in.levels},
                {"eigensolver", in.eigensolver},
                {"linearSolver", in.linear_solver},
                {"preconditioner", in.preconditioner},
                {"eigenvalueTolerance", in.eigenvalue_tolerance},
                {"residualTolerance", in.residual_tolerance},
                {"linearTolerance", in.linear_tolerance},
                {"maxIterations", in.max_iterations},
                {"maxLinearIterations", in.max_linear_iterations},
                {"jdInnerIterations", in.jd_inner_iterations},
                {"subspaceSize", in.subspace_size},
                {"expectationValues", in.expectation_values},
                {"threads", in.threads}};
}

nlohmann::json nuclear_result_json(const NuclearMotionResult& r) {
    using nlohmann::json;
    auto moments = [](const DistanceMoments& m) {
        return json{{"mean", m.mean}, {"rms", m.rms}, {"cubic", m.cubic},
                    {"inverse", m.inverse}, {"inverse2", m.inverse2}, {"inverse3", m.inverse3}};
    };
    json j;
    j["error"] = nullptr;
    j["input"] = nuclear_input_json(r.input);
    j["isotopologue"] = r.isotopologue_name;
    j["masses"] = {r.masses[0], r.masses[1], r.masses[2]};
    j["mu1"] = r.mu1;
    j["mu2"] = r.mu2;
    j["potentialMinimum"] = r.potential_minimum;
    j["minimumR"] = r.minimum_r;
    j["minimumAngle"] = r.minimum_angle * 180.0 / 3.14159265358979323846;
    j["rGrid"] = r.r_grid;
    j["kValues"] = r.k_values;
    j["cosGrid"] = r.cos_grid;
    j["dimension"] = r.dimension;
    j["nonzeros"] = r.nonzeros;
    j["hamiltonianSeconds"] = r.hamiltonian_seconds;
    j["solverSeconds"] = r.solver_seconds;
    j["iterations"] = r.iterations;
    j["matrixVectorProducts"] = r.matrix_vector_products;
    j["linearIterations"] = r.linear_iterations;
    j["converged"] = r.converged;
    j["solver"] = r.solver_description;
    json hist = json::array();
    for (const NuclearSolverStep& s : r.history)
        hist.push_back(json{{"iteration", s.iteration}, {"lowest", s.lowest}, {"change", std::isfinite(s.change) ? json(s.change) : json(nullptr)},
                            {"converged", s.converged}, {"linearIterations", s.linear_iterations}});
    j["history"] = hist;
    json levels = json::array();
    for (const NuclearLevel& l : r.levels)
        levels.push_back(json{{"energy", l.energy}, {"residual", l.residual}, {"r1", moments(l.r1)}, {"r2", moments(l.r2)},
                              {"theta", l.theta * 180.0 / 3.14159265358979323846}, {"cosTheta", l.cos_theta},
                              {"kWeights", l.k_weights}});
    j["levels"] = levels;
    j["notes"] = r.notes;
    j["output"] = r.report;
    return j;
}

}  // namespace macqc

namespace {
char* nuclear_copy(const std::string& s) {
    char* buf = static_cast<char*>(std::malloc(s.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, s.c_str(), s.size() + 1);
    return buf;
}
}  // namespace

extern "C" {

char* macqc_nuclear_run_job(const char* input_json, macqc_job* job) {
    macqc::JobScope scope(job ? &job->job : nullptr);
    using nlohmann::json;
    std::string result;
    try {
        json j;
        try {
            const macqc::NuclearMotionInput in = macqc::parse_nuclear_input_json(input_json ? input_json : "");
            macqc::NuclearOutputSink sink;
            if (job) sink = [job](const std::string& text) { job->job.append_output(text); };
            j = macqc::nuclear_result_json(macqc::run_nuclear_motion(in, sink));
        } catch (const std::exception& e) {
            j = json{{"error", e.what()}, {"levels", json::array()}};
            if (job) {
                std::size_t total = 0;
                j["output"] = job->job.output_since(0, &total);
            }
        }
        result = j.dump(-1, ' ', false, json::error_handler_t::replace);
        macqc::report_progress(macqc::JobStage::Done);
    } catch (...) {
        result = "{\"error\":\"internal error\",\"levels\":[]}";
    }
    return nuclear_copy(result);
}

char* macqc_nuclear_estimate(const char* input_json) {
    using nlohmann::json;
    json j;
    try {
        const macqc::NuclearMotionInput in = macqc::parse_nuclear_input_json(input_json ? input_json : "");
        const macqc::NuclearMotionEstimate e = macqc::estimate_nuclear_motion(in);
        j = json{{"error", nullptr}, {"dimension", e.dimension}, {"nonzeros", e.nonzeros}, {"memoryMB", e.memory_mb},
                 {"refusal", e.refusal.empty() ? json(nullptr) : json(e.refusal)}};
    } catch (const std::exception& e) {
        j = json{{"error", e.what()}, {"refusal", e.what()}};
    } catch (...) {
        j = json{{"error", "internal error"}, {"refusal", "internal error"}};
    }
    return nuclear_copy(j.dump(-1, ' ', false, json::error_handler_t::replace));
}

char* macqc_run_job(const char* input_json, macqc_job* job) {
    macqc::JobScope scope(job ? &job->job : nullptr);
    std::string result;
    try {
        macqc::Output out;
        try {
            const macqc::Input in = macqc::parse_input_json(input_json ? input_json : "");
            macqc::OutputSink sink;
            if (job) sink = [job](const std::string& text) { job->job.append_output(text); };
            out = macqc::run_calculation(in, sink);
        } catch (const std::exception& e) {
            out = macqc::Output();
            out.error = e.what();
        }
        result = macqc::output_to_json(out);
        macqc::report_progress(macqc::JobStage::Done);
    } catch (...) {
        result = "{\"error\":\"internal error\",\"hf\":null,\"fci\":null,\"report\":\"\"}";
    }
    char* buf = static_cast<char*>(std::malloc(result.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, result.c_str(), result.size() + 1);
    return buf;
}

void macqc_free(char* p) { std::free(p); }

}  // extern "C"

namespace {

template <class T>
T* copy_array(const std::vector<T>& v) {
    if (v.empty()) return nullptr;
    T* p = static_cast<T*>(std::malloc(v.size() * sizeof(T)));
    if (p) std::memcpy(p, v.data(), v.size() * sizeof(T));
    return p;
}

void clear_mesh(macqc_mesh* m) {
    if (!m) return;
    m->vertex_count = 0;
    m->triangle_count = 0;
    m->positions = nullptr;
    m->normals = nullptr;
    m->indices = nullptr;
}

void fill_mesh(const macqc::Mesh& src, macqc_mesh* m) {
    clear_mesh(m);
    if (!m || src.vertex_count() == 0 || src.triangle_count() == 0) return;
    m->positions = copy_array(src.positions);
    m->normals = copy_array(src.normals);
    m->indices = copy_array(src.indices);
    if (!m->positions || !m->normals || !m->indices) {
        macqc_mesh_free(m);
        throw std::bad_alloc();
    }
    m->vertex_count = static_cast<int>(src.vertex_count());
    m->triangle_count = static_cast<int>(src.triangle_count());
}

}  // namespace

extern "C" {

char* macqc_orbital_surface(const char* orbitals_json, int set, int orbital, double isovalue, int max_points,
                            macqc_mesh* positive, macqc_mesh* negative) {
    clear_mesh(positive);
    clear_mesh(negative);
    try {
        if (!orbitals_json) throw std::runtime_error("orbitals: no input");
        const macqc::OrbitalData d = macqc::orbitals_from_json(orbitals_json);
        const std::vector<double>& c = macqc::chosen_orbital(d, set, orbital);
        const macqc::OrbitalSurfaces s =
            macqc::orbital_surfaces(d, c, isovalue, max_points > 0 ? max_points : 400000, 0);
        fill_mesh(s.positive, positive);
        fill_mesh(s.negative, negative);
        return nullptr;
    } catch (const std::exception& e) {
        macqc_mesh_free(positive);
        macqc_mesh_free(negative);
        return copy_string(e.what());
    }
}

char* macqc_orbital_plane(const char* orbitals_json, int set, int orbital, const double* origin, const double* u,
                          const double* v, int nu, int nv, float* values) {
    try {
        if (!orbitals_json || !origin || !u || !v || !values) throw std::runtime_error("orbitals: missing argument");
        if (nu < 1 || nv < 1 || static_cast<long long>(nu) * nv > 16000000LL)
            throw std::runtime_error("orbitals: the plane must have 1 to 16 million points");
        const macqc::OrbitalData d = macqc::orbitals_from_json(orbitals_json);
        const std::vector<double>& c = macqc::chosen_orbital(d, set, orbital);
        macqc::GridSpec g;
        for (int k = 0; k < 3; ++k) {
            g.origin[k] = origin[k];
            g.axes[0][k] = u[k];
            g.axes[1][k] = v[k];
            g.axes[2][k] = 0.0;
        }
        g.n = {nu, nv, 1};
        const std::vector<double> w = macqc::orbital_grid(d, c, g, 0);
        for (std::size_t k = 0; k < w.size(); ++k) values[k] = static_cast<float>(w[k]);
        return nullptr;
    } catch (const std::exception& e) {
        return copy_string(e.what());
    }
}

void macqc_mesh_free(macqc_mesh* mesh) {
    if (!mesh) return;
    std::free(mesh->positions);
    std::free(mesh->normals);
    std::free(mesh->indices);
    clear_mesh(mesh);
}

char* macqc_parse_structure(const char* text, const char* hint) {
    using nlohmann::json;
    json j;
    try {
        const macqc::StructureFile f = macqc::parse_structure_file(text ? text : "", hint ? hint : "");
        j["error"] = nullptr;
        j["format"] = f.format;
        j["title"] = f.title;
        j["charge"] = f.charge;
        j["multiplicity"] = f.multiplicity;
        j["hasCharge"] = f.has_charge;
        j["warnings"] = f.warnings;
        j["atoms"] = macqc::custom_atoms_json(f.atoms);
    } catch (const std::exception& e) {
        j = {{"error", e.what()}, {"atoms", json::array()}, {"warnings", json::array()}};
    } catch (...) {
        j = {{"error", "internal error"}, {"atoms", json::array()}, {"warnings", json::array()}};
    }
    const std::string s = j.dump(-1, ' ', false, json::error_handler_t::replace);
    char* buf = static_cast<char*>(std::malloc(s.size() + 1));
    if (!buf) return nullptr;
    std::memcpy(buf, s.c_str(), s.size() + 1);
    return buf;
}

}  // extern "C"
