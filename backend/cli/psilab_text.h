// ΨLab terminal version: plain-text input file -> plain-text output file.
//
// Input: one "key = value" per line, '#' starts a comment. Keys (case-insensitive):
//   molecule   HeH+ (default) | H2 | LiH | H2O | H2O2 | NH3
//   basis      aug-cc-pVXZ, aug-cc-pCVXZ (X = D, T, Q, 5); HeH+ also STO-1G | STO-2G | STO-3G, the
//              Szabo-Ostlund calculation (HeH+ with a correlation-consistent set: all methods and jobs)
//   R          bond length (H2O, H2O2: O-H; NH3: N-H), a number optionally followed by a unit: bohr (default)
//              or angstrom; default for H2, LiH, H2O: the experimental equilibrium value
//   angle      degrees: H-O-H (H2O; default 104.52), O-O-H (H2O2; 102.32), H-N-H (NH3; 106.7, 120 planar)
//   R_OO       O-O bond length (H2O2; default 1.4556 angstrom)   dihedral  H-O-O-H (H2O2; 113.70)
//   method     HF | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI (H2, LiH, H2O; refused where impossible)
//   job        energy (default) | opt | freq | opt+freq | opt+freq+sp | scan | ts | ts+freq | ts+freq+sp
//   scan       job = scan: the Gaussian ModRedundant line N1 N2 [N3 [N4]] S nsteps stepsize
//              (1-based atoms, optionally preceded by B, A or D; step in angstrom for a bond,
//              degrees for an angle or dihedral), e.g. 1 2 S 10 0.05, 2 1 3 S 8 10.0 or
//              3 1 2 4 S 18 10.0 (H2O2 torsion)
//   scan_start value of the first scan point (angstrom or degrees; default: the input geometry)
//   ts_hessian job = ts, ts+freq, ts+freq+sp (Gaussian Opt=TS): calcfc (Hessian at the first step, default) | calcall
//   eigen_test yes (default: exactly one negative Hessian eigenvalue) | no (Gaussian NoEigenTest)
//   sp_method, sp_basis   method and basis set of the single point of opt+freq+sp
//   temperature  K (default 298.15)       pressure  Pa (default 101325; units Pa, kPa, bar, atm)
//   opt_convergence  normal | tight       max_steps  optimisation steps (default 50)
//   threads, memory, fci_memory (MB or GB), basis_dir      (H2, LiH, H2O)
//   zeta_He, zeta_H, exact_erf, scan, scan_step (HeH+ only):
//   zeta_He    Slater exponent of He          zeta_H     Slater exponent of H
//   exact_erf  yes | no                       scan       yes | no
//   scan_step  step of the potential-curve scan (bohr), R = 0.8 ... 5.0 bohr
// The output text matches the macOS app's File > Export Results as Text. For H2, LiH and H2O
// it is the output file of ΨLab 5.0, written while the job runs.
#pragma once

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "macqc/calculation.h"

namespace psilab {

// ΨLab version (same as the macOS app).
constexpr const char* kVersion = "8.0";

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrToAngstrom = 0.529177210544;

constexpr double kScanStart = 0.8;  // bohr
constexpr double kScanEnd = 5.0;    // bohr

struct Settings {
    macqc::Input input;     // R always in bohr
    bool scan = false;
    double scanStep = 0.1;  // bohr
    bool rGiven = false;    // R was in the input (else the molecule's default)
    std::string sourceText; // the input file (echoed in the output file)
};

// The "Input" section of the output file: the keys as understood, then the file as given.
std::string input_echo_text(const Settings& settings);

// Parses the text input. Throws std::runtime_error("line N: ...") on any problem,
// including values the calculation would reject.
Settings parse_input_text(const std::string& text);

// A commented example input file (the Szabo-Ostlund default case).
std::string example_input();

struct CurvePoint {
    double R, hf, fci;  // bohr, E_h, E_h
};

struct CurveMinimum {
    std::string method;
    double R = 0.0, energy = 0.0;
    bool refined = false;  // parabola through the lowest point and its neighbours
};

std::vector<double> scan_grid(double step);
std::optional<CurveMinimum> curve_minimum(const std::vector<std::pair<double, double>>& samples,
                                          const std::string& method);

struct RunResult {
    std::string text;   // the complete output file (empty only if nothing could be computed)
    std::string error;  // empty on success
    std::optional<macqc::Output> output;  // the single-point calculation
};

// Runs the calculation (and the scan if requested) and formats the output. H2, LiH, H2O: the
// output file is also passed to `sink` piece by piece while the job runs.
RunResult run_to_text(const Settings& settings, const macqc::OutputSink& sink = macqc::OutputSink());

// H2, LiH, H2O: size of the calculation, time estimate and warnings, without computing
// (empty for HeH+).
std::string estimate_text(const Settings& settings);
// The warnings of the estimate (e.g. a long calculation), for the terminal.
std::vector<std::string> estimate_warnings(const Settings& settings);

}  // namespace psilab
