// Job types of ΨLab 5.0 for H2, LiH and H2O, with the complete output file:
//   energy        single-point energy (as in 4.0)
//   opt           geometry optimisation
//   freq          harmonic frequencies and thermochemistry at the input geometry
//   opt+freq      optimisation, then frequencies at the optimised geometry
//   opt+freq+sp   optimisation, frequencies, then a single-point energy with its own method
//                 and basis set at the optimised geometry; the final energy is reported with the
//                 thermal corrections of the frequency step
// If the optimisation does not converge, jobs 4 and 5 stop after it with a message.
// ΨLab 7.0 adds:
//   scan          relaxed potential-energy-surface scan of one bond, angle or dihedral angle
//                 (Gaussian: Opt=ModRedundant with the line "N1 N2 [N3 [N4]] S nsteps stepsize"):
//                 at each of the nsteps + 1 points the coordinate is frozen and all others are
//                 optimised; the job ends with a table of the energy against the coordinate
//   ts            transition-state optimisation (Gaussian: Opt=TS, with CalcFC or CalcAll):
//                 P-RFO steps, Bofill updates, eigenvalue test (see optimizer.h)
//   ts+freq       then harmonic frequencies at the transition state; the job ends with an error
//                 unless there is exactly one imaginary frequency
//   ts+freq+sp    then, at the confirmed transition state, a single-point energy with its own
//                 method and basis set; final U, H, G, A = single point + thermal corrections
//                 (the transition-state counterpart of opt+freq+sp)
#pragma once

#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "macqc/frequencies.h"
#include "macqc/gradient.h"
#include "macqc/hessian.h"
#include "macqc/molecular.h"
#include "macqc/optimizer.h"

namespace macqc {

struct JobSettings {
    std::string type = "energy";      // energy | opt | freq | opt+freq | opt+freq+sp | scan | ts | ts+freq | ts+freq+sp
    std::string sp_method, sp_basis;  // job 5 (default: the method and basis set of the job)
    double temperature = 298.15;      // K
    double pressure = 101325.0;       // Pa
    bool tight = false;               // optimisation convergence
    int max_steps = 50;
    // 7.0 relaxed scan (see parse_scan_line): coordinate with 1-based atom numbers ("B 1 2",
    // "A 2 1 3", "D 3 1 2 4"), number of steps (points = steps + 1), step size (bohr for bonds,
    // degrees for angles) and the value of the first point (empty: the input geometry's)
    std::string scan_coordinate;
    int scan_steps = 10;
    double scan_step = 0.0;
    std::optional<double> scan_start;
    // 7.0 transition states: the Hessian computed at the first point ("calcfc", Gaussian's
    // default for Opt=(TS,CalcFC)) or at every point ("calcall"); the eigenvalue test
    std::string ts_hessian = "calcfc";
    bool eigen_test = true;
    // 10.0: "analytic" or "numerical" Hessians of Full CI, CCSD and CCSD(T)
    std::string hessian = "analytic";
    std::string initial_hessian = "auto";  // minimisations: auto | analytic | lindh
    // 15.0: computed Hessians of minimisations ("" | calcfc | calcall | recalcfc; ts_hessian also
    // takes recalcfc), N of RecalcFC=N, the step of minimisations (rfo | gediis)
    std::string opt_hessian;
    int recalc_fc = 0;
    std::string opt_algorithm = "rfo";
};

// The scanned coordinate: type (Stretch, Bend, Dihedral) and 0-based atoms.
struct ScanCoordinate {
    InternalCoordinate coordinate;
    std::string letter;                 // "B", "A", "D"
    std::string kind;                   // "bond", "angle", "dihedral"
    std::string definition;             // e.g. "D(H3–O1–O2–H4)"
};
// Parses a coordinate "N1 N2 [N3 [N4]]" (2 atoms: bond, 3: angle, 4: dihedral), optionally
// preceded by its type as in Gaussian's ModRedundant input: "B 1 2" (also "R", "bond"),
// "A 2 1 3" ("angle"), "D 3 1 2 4" ("dihedral"); commas allowed. Returns "" or an error message
// (unknown type, wrong number of atoms, atom out of range or repeated).
std::string parse_scan_coordinate(const std::string& spec, int n_atoms, ScanCoordinate* out);

// Parses a Gaussian ModRedundant scan line "[B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize" (step size
// in angstrom for a bond, degrees for an angle or a dihedral angle) into settings->scan_coordinate,
// scan_steps and scan_step (bohr or degrees). Returns "" or an error message.
std::string parse_scan_line(const std::string& line, int n_atoms, JobSettings* settings);

// The points of a scan: nominal values (bohr or degrees, start + k step) and the targets of the
// constraint (bohr or rad; dihedrals reduced to (-pi, pi]). Returns "" or why the scan is not
// possible at this geometry (bond lengths outside 0.3-20 bohr, angles outside 10-180 degrees, an
// undefined dihedral).
struct ScanPlan {
    ScanCoordinate coordinate;
    double start = 0.0, step = 0.0;     // bohr or degrees
    int steps = 0;
    std::vector<double> values;         // bohr or degrees
    std::vector<double> targets;        // bohr or rad
    double input_value = 0.0;           // at the input geometry (bohr or degrees)
    bool angular() const { return coordinate.letter != "B"; }
};
std::string plan_scan(const Molecule& mol, const JobSettings& settings, ScanPlan* plan);

struct ScanPoint {
    int index = 0;                      // 0-based
    double value = 0.0;                 // nominal (bohr or degrees)
    double measured = 0.0;              // at the optimised geometry (bohr or degrees)
    double energy = 0.0;
    bool converged = false;
    int steps = 0;                      // optimisation steps (0: single point, diatomic)
    double max_force = 0.0;             // projected, at the last step
    Molecule geometry;
    std::vector<GeometricParameter> parameters;
    std::vector<double> step_energies, step_forces;  // per optimisation step
    std::vector<std::string> warnings;
    double seconds = 0.0;
};

struct ScanExtremum {
    std::string kind;                   // "minimum" or "maximum"
    int index = 0;                      // the scan point
    double value = 0.0, energy = 0.0;   // parabola through the point and its neighbours
};

struct ScanResult {
    ScanPlan plan;
    bool relaxed = true;                // false for a diatomic (nothing to relax: single points)
    std::vector<ScanPoint> points;
    int lowest = -1, highest = -1;
    std::vector<ScanExtremum> extrema;
    bool converged = true;              // every point
    std::string message;
    double seconds = 0.0;
};

// Canonical job type for any spelling ("opt+freq", "OPT FREQ", "optfreq", "4"), or "".
std::string canonical_job_type(const std::string& type);
// "Energy", "Geometry optimisation", "Frequencies", "Optimisation + frequencies",
// "Optimisation + frequencies + single-point energy".
std::string job_title(const std::string& type);
bool job_has_optimization(const std::string& type);
bool job_has_frequencies(const std::string& type);
bool job_has_single_point(const std::string& type);
bool job_has_scan(const std::string& type);
bool job_has_ts(const std::string& type);  // ts, ts+freq, ts+freq+sp (job_has_optimization is true for them too)
// Jobs that run the optimiser (and take opt_convergence and max_steps): opt, opt+freq,
// opt+freq+sp, scan, ts, ts+freq.
bool job_uses_optimizer(const std::string& type);

struct JobStepEstimate {
    std::string name;          // "Geometry optimisation", "Hessian", "Single-point energy", "Energy"
    std::string method, basis;
    std::string how;           // e.g. "analytic gradients (Lambda equations, relaxed densities, Z-vector)"
    bool analytic = true;
    int evaluations = 0;       // energies or gradients (per optimisation step for "opt")
    std::string evaluation_kind;  // "gradients", "energies", "steps"
    double seconds = 0.0;
    double memory_mb = 0.0;
};

struct JobEstimate {
    std::vector<JobStepEstimate> steps;
    double total_seconds = 0.0;
    double memory_mb = 0.0;
    std::string refusal;                 // non-empty: the job is not possible
    std::vector<std::string> warnings;   // long calculation etc.
    std::vector<std::string> notes;      // what is analytic / numerical, ...
};

// Uses prepare_molecular() for the cost of single energies (sp_prepared: the method and basis
// of the single-point step of job 5, or nullptr).
JobEstimate estimate_job(const MolecularResult& prepared, const JobSettings& settings,
                         const MolecularResult* sp_prepared);

struct JobResult {
    JobSettings settings;
    JobEstimate estimate;
    bool completed = false;
    std::string termination;           // "normal", "error", "refused", "cancelled"
    std::string message;               // reason for an error or a stop
    Molecule initial_geometry, final_geometry;
    std::optional<SinglePoint> point;           // energy job / single point at the final geometry
    std::optional<OptimizationResult> optimization;
    std::optional<HessianResult> hessian;
    std::optional<HarmonicAnalysis> harmonic;
    std::optional<ThermoResult> thermo;
    std::optional<SinglePoint> single_point;    // job 5
    std::optional<ScanResult> scan;             // 7.0 relaxed scan
    std::string single_point_method, single_point_basis;
    double final_energy = 0.0;
    std::string final_energy_label;             // e.g. "CCSD(T)/aug-cc-pVTZ"
    std::vector<std::string> warnings;
    double seconds = 0.0;
    std::string output;                         // the complete output file
};

using OutputSink = std::function<void(const std::string&)>;

// 13.0: writes the appendix of the output file (the basis set in full, the molecular orbitals and
// the archive of the complete results) for a finished job; called before the termination banner
// (not for a cancelled job).
using AppendixWriter = std::function<std::string(const JobResult&)>;

// Runs a job (not "energy": that path is run_energy_job). `input` gives molecule, starting
// geometry, method, basis and limits. Never throws for calculation problems: they end the
// job with termination "error" and the reason in the output; Cancelled is rethrown after the
// output says so.
JobResult run_job(const MolecularInput& input, const JobSettings& settings, const OutputSink& sink,
                  const std::string& input_echo, const AppendixWriter& appendix = AppendixWriter());

// The energy job with the 5.0 output file: runs run_molecular (4.0 results unchanged) and
// writes the output around it.
JobResult run_energy_job(MolecularResult& prepared, const JobSettings& settings, const OutputSink& sink,
                         const std::string& input_echo, const AppendixWriter& appendix = AppendixWriter());

// The "Size, memory and time (estimates)" section of the output file.
std::string job_estimate_text(const MolecularResult& prepared, const JobEstimate& estimate);

// A MolecularResult (4.0 structure) for the final energy calculation of a job.
MolecularResult final_point_result(const JobResult& job, const MolecularInput& input);

// The program version.
constexpr const char* kProgramVersion = "16.0";

}  // namespace macqc
