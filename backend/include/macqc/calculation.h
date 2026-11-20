// High-level driver: validated input -> HF -> FCI -> text report / JSON.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "macqc/fci.h"
#include "macqc/jobs.h"
#include "macqc/molecular.h"
#include "macqc/scf.h"

namespace macqc {

struct Input {
    // "HeH+" (default): the Szabo-Ostlund STO-nG calculation (HF + FCI), unchanged
    // since version 1.0, which uses the keys R .. printLevel below. With a basis set of the library
    // (hehGeneral, 7.0; 14.0 also Pople and def2) HeH+ is a general molecule like the others (charge +1).
    // "H2", "LiH", "H2O": general molecules (version 3.0), see the keys after printLevel;
    // "H2O2", "NH3" (7.0).
    std::string molecule = "HeH+";
    double R = 1.4632;  // bohr; for H2O the O-H distance
    double zeta1 = 2.0925;
    double zeta2 = 1.24;
    double ZA = 2.0;
    double ZB = 1.0;
    int nGaussians = 3;
    bool exactErf = false;
    int maxIter = 25;
    double convergence = 1.0e-4;
    int printLevel = 2;

    // General molecules (molecule != "HeH+", or HeH+ with hehGeneral).
    bool hehGeneral = false;            // HeH+ with a basis set of the library (all methods and jobs)
    std::string basis = "aug-cc-pVDZ";
    double angle = kDefaultWaterAngle;  // degrees: H-O-H (H2O), O-O-H (H2O2), H-N-H (NH3)
    double ROO = 0.0;                   // H2O2: R(O-O) in bohr (<= 0: the default 1.4556 Å)
    double dihedral = kDefaultPeroxideDihedral;  // H2O2: H-O-O-H dihedral angle (degrees)
    // 9.0, molecule "custom": the atoms (positions in bohr; the JSON key "atoms" gives them in
    // angstrom), the charge, the spin multiplicity (1) and the title from the structure file
    std::vector<Atom> atoms;
    int charge = 0;
    int multiplicity = 1;
    std::string title;
    std::string geometryFile;           // terminal version: the structure file the atoms were read from
    std::string method = "HF";          // 15.0 default "HF" (14.0: "RHF"); "RHF" (HF), "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ" or "FCI"; 15.0: "UHF",
                                        // "ROHF", "MP2" ... "MP5", "CISD" ... "QCISDTQ", prefixes U/RO
    // 15.0: the reference ("RHF", "UHF", "ROHF"; "": RHF for closed-shell singlets, UHF otherwise);
    // "charge" and "multiplicity" apply to every molecule (given: chargeGiven, multiplicityGiven);
    // a preset of the molecule menu ("O2", "NO2", "H5+", an atom "N") runs as a custom molecule;
    // HeH+ with STO-nG is the Szabo-Ostlund calculation unless szaboOstlund is false.
    std::string reference;
    bool chargeGiven = false;
    bool multiplicityGiven = false;
    bool szaboOstlund = true;
    std::string basisDirectory;         // "" -> default_basis_directory()
    int threads = 0;                    // 0 -> all hardware threads
    double memoryMB = 1024.0;           // limit for keeping the two-electron integrals in memory
    double fciMemoryMB = 4096.0;        // hard memory limit of Full-CI and coupled cluster
    int scfMaxIter = 100;
    double energyConvergence = 1.0e-10;
    double gradientConvergence = 1.0e-7;
    bool estimateOnly = false;          // only molecule, basis, size and time estimate, warnings

    // ΨLab 5.0 job types (general molecules).
    std::string job = "energy";         // energy | opt | freq | opt+freq | opt+freq+sp | scan | ts | ts+freq | ts+freq+sp
    std::string spMethod;               // opt+freq+sp, ts+freq+sp: method of the single point ("": method)
    std::string spBasis;                // opt+freq+sp, ts+freq+sp: basis set of the single point ("": basis)
    double temperature = 298.15;        // K (frequency jobs)
    double pressure = 101325.0;         // Pa (frequency jobs)
    std::string optConvergence = "normal";  // normal | tight
    int maxSteps = 50;                  // optimisation steps (per point of a scan)

    // ΨLab 7.0 relaxed scan (job "scan"): the Gaussian ModRedundant line with 1-based atoms,
    // "[B|A|D] N1 N2 [N3 [N4]] S nsteps stepsize", step in angstrom (bond) or degrees (angle,
    // dihedral), e.g. "1 2 S 10 0.05" or "A 2 1 3 S 8 10.0"; nsteps + 1 points
    std::string scan;
    std::optional<double> scanStart;    // value of the first point (angstrom or degrees); default: the input geometry
    // ΨLab 7.0 transition states (jobs "ts", "ts+freq")
    std::string tsHessian = "calcfc";   // calcfc (Hessian at the first point) | calcall (at every point) | recalcfc (15.0)
    bool eigenTest = true;              // exactly one negative Hessian eigenvalue (Gaussian: EigenTest)
    // ΨLab 10.0: the Hessian of frequency jobs and transition-state searches with Full-CI, CCSD and
    // CCSD(T): "analytic" (where it fits in memory) or "numerical" (differences of analytic gradients)
    std::string hessian = "analytic";
    // 10.0: the initial Hessian of minimisations: "auto" (analytic for HF, Lindh's model Hessian
    // otherwise), "analytic" (also for Full-CI, CCSD, CCSD(T)) or "model" (Lindh)
    std::string initialHessian = "auto";
    // 15.0 (Gaussian Opt=CalcFC, CalcAll, RecalcFC=N): the computed Hessians of minimisations,
    // "auto" (the initial Hessian of initialHessian, then BFGS updates) | "calcfc" | "calcall" |
    // "recalcfc" (at the first step and every recalcFC-th step thereafter, BFGS updates between);
    // tsHessian also takes "recalcfc" (Bofill updates between)
    std::string optHessian = "auto";
    int recalcFC = 0;
    // 15.0: the step of minimisations, "rfo" (default) or "gediis" (the hybrid GEDIIS/RFO step of
    // Li and Frisch 2006, Gaussian's default)
    std::string optAlgorithm = "rfo";
    std::string echo;                   // input text for the output file ("": generated from the keys)
};

// True if the input selects the HeH+ STO-nG calculation of versions 1 and 2.
bool is_heh_input(const Input& in);

struct Output {
    Input input;
    std::string error;  // empty: no error
    std::optional<HFResult> hf;
    std::optional<FCIResult> fci;
    std::string hf_report;   // HF part only (as printed by the reference program)
    std::string report;      // hf_report + FCI section (general molecules: their report)
    std::optional<MolecularResult> molecular;  // general molecules only (the final energy calculation)
    std::optional<JobResult> job;              // general molecules (5.0): the job and its output file
};

// Returns an empty string if the input is valid, else a message.
std::string validate(const Input& in);

// Never throws; problems are reported through Output::error. General molecules: the output
// file is passed to `sink` piece by piece while the job runs (it is also in Output::job->output).
Output run_calculation(const Input& in, const OutputSink& sink = OutputSink());

// The input as "key = value" lines (for the output file).
std::string input_echo(const Input& in);

// Parse input JSON text (missing keys take defaults). Throws std::runtime_error
// with a descriptive message on malformed input.
Input parse_input_json(const std::string& text);

// Serialise to the documented output JSON (never throws).
std::string output_to_json(const Output& out, int indent = -1);

// 13.0: the appendix of an output file, written before its termination banner: for the general
// molecules the basis set in full and the molecular orbitals of the app's Orbitals tab (AO
// coefficients), then the archive: the complete results (the JSON of output_to_json without
// "output") between the lines "@@PSILAB-ARCHIVE-BEGIN 1" and "@@PSILAB-ARCHIVE-END", ASCII,
// wrapped at 100 columns. The app reads it back (File > Open Output File…).
std::string output_appendix(const Output& out);

// The archive JSON in the text of an output file (the lines between the markers joined); ""
// if there is none.
std::string archive_from_output_text(const std::string& text);

}  // namespace macqc
