#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <new>
#include <stdexcept>
#include <string>

#include "macqc/calculation.h"
#include "macqc/job.h"
#include "macqc/jobs.h"
#include "macqc/macqc_api.h"
#include "macqc/orbitals.h"
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
    if (!is_heh_input(in) && !canonical_molecule_name(in.molecule).empty()) {
        json j = input_json_general(in);
        if (canonical_molecule_name(in.molecule) == "custom") {
            j["atoms"] = custom_atoms_json(in.atoms);
            j["charge"] = in.charge;
            j["multiplicity"] = in.multiplicity;
            j["title"] = in.title;
        }
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
                {"initialHessian", in.initialHessian}};
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
    mol["electrons"] = m.molecule.electrons();
    mol["charge"] = m.molecule.charge;
    mol["R"] = m.input.R;
    const bool has_angle = m.molecule.name == "H2O" || m.molecule.name == "H2O2" || m.molecule.name == "NH3";
    mol["angle"] = has_angle ? json(m.input.angle) : json(nullptr);
    mol["ROO"] = m.molecule.name == "H2O2" ? json(m.input.r_oo) : json(nullptr);
    mol["dihedral"] = m.molecule.name == "H2O2" ? json(m.input.dihedral) : json(nullptr);
    mol["nuclearRepulsion"] = m.molecule.nuclear_repulsion();
    json atoms = json::array();
    for (const Atom& a : m.molecule.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    mol["atoms"] = atoms;

    json basis;
    basis["name"] = m.basis.name;
    basis["functions"] = m.basis.n_functions;
    basis["shells"] = m.basis.shells.size();
    basis["maxL"] = m.basis.max_l;
    json ba = json::array();
    for (const auto& a : m.basis.atoms)
        ba.push_back({{"symbol", a.symbol},
                      {"basis", a.basis},
                      {"contraction", a.contraction},
                      {"primitives", a.primitives},
                      {"functions", a.functions}});
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
        out["rhf"] = r;
    } else {
        out["rhf"] = nullptr;
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
                       {"hessianComputed", s.hessian_computed}};
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
// spherical functions m = -l .. l from "first" on) and the orbital sets with their AO
// coefficients. The app passes it back to macqc_orbital_surface and macqc_orbital_plane.
json orbitals_json(const OrbitalData& d) {
    json j;
    j["basis"] = d.basis_name;
    j["functions"] = d.n_basis;
    json atoms = json::array();
    for (const Atom& a : d.atoms)
        atoms.push_back({{"symbol", a.symbol}, {"Z", a.Z}, {"position", {a.position[0], a.position[1], a.position[2]}}});
    j["atoms"] = atoms;
    json shells = json::array();
    for (const Shell& s : d.shells)
        shells.push_back({{"l", s.l},
                          {"atom", s.atom},
                          {"center", {s.center[0], s.center[1], s.center[2]}},
                          {"exponents", s.exponents},
                          {"coefficients", s.coefficients},
                          {"first", s.first_function}});
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
    const std::string canonical = canonical_molecule_name(in.molecule);
    if (canonical.empty()) throw std::runtime_error("\"molecule\" must be HeH+, H2, LiH, H2O, H2O2 or NH3");
    in.molecule = canonical;
    if (canonical != "HeH+") in.R = 0.0;  // default bond length of the molecule unless given
    if (canonical == "HeH+") {
        // 7.0: a correlation-consistent basis set selects the general path; STO-1G/2G/3G (or no
        // basis) the Szabo-Ostlund calculation
        const std::string b = get_string(j, "basis", "");
        if (!b.empty() && !canonical_basis_name(b).empty()) {
            in.hehGeneral = true;
        } else if (!b.empty()) {
            std::string t;
            for (char c : b)
                if (c != '-' && c != ' ') t += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
            if (t == "sto1g" || t == "sto2g" || t == "sto3g") in.nGaussians = t[3] - '0';
            else throw std::runtime_error("\"basis\" for HeH+ must be STO-1G, STO-2G, STO-3G or a correlation-consistent set");
        }
    }
    if (canonical != "HeH+" && canonical != "custom") in.angle = default_angle(canonical);
    if (canonical == "custom") {
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
        in.charge = get_int(j, "charge", 0);
        in.multiplicity = get_int(j, "multiplicity", 1);
        in.title = get_string(j, "title", "");
        in.R = 0.0;
    }
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
    for (char& c : in.hessian) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (auto it = j.find("scanStart"); it != j.end() && !it->is_null()) in.scanStart = get_number(j, "scanStart", 0.0);
    return in;
}

std::string output_to_json(const Output& out, int indent) {
    json j;
    j["input"] = input_json(out.input);
    j["error"] = out.error.empty() ? json(nullptr) : json(out.error);
    if (!is_heh_input(out.input) && !canonical_molecule_name(out.input.molecule).empty()) {
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
        if (out.molecular && out.molecular->rhf && out.molecular->rhf->C.rows() > 0) {
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
        return j.dump(indent, ' ', false, json::error_handler_t::replace);
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
    return j.dump(indent, ' ', false, json::error_handler_t::replace);
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
