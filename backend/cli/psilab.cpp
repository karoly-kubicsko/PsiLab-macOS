// ΨLab for the terminal (Linux, macOS): reads a plain-text input file, writes a
// plain-text output file.  Usage: psilab INPUT > OUTPUT (15.0) or psilab INPUT OUTPUT
// (see psilab --help)
#include <algorithm>
#include <cctype>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <chrono>

#include "psilab_text.h"
#include "macqc/nuclear_motion.h"
#include "macqc/output.h"

#if defined(__APPLE__)
#include <mach-o/dyld.h>
#endif
#include <sys/stat.h>
#include <unistd.h>

namespace {

void usage(const char* prog) {
    std::fprintf(stdout,
                 "\xCE\xA8Lab - Hartree–Fock, perturbation theory, CI, QCI, coupled cluster and Full CI in the terminal\n\n"
                 "Usage:\n"
                 "  %s INPUT > OUTPUT   run INPUT, write the output file to the standard output (15.0;\n"
                 "                      the summary and messages go to the standard error)\n"
                 "  %s INPUT OUTPUT     run INPUT, write the output file to OUTPUT\n"
                 "  %s INPUT -o         run INPUT, write INPUT with the extension replaced by .out\n"
                 "  %s INPUT key=value ...  settings on the command line (19.0), e.g. method=MP2 job=opt\n"
                 "                      (quote values with parentheses: method='CCSD(T)')\n"
                 "  %s INPUT1 INPUT2 ... several jobs, one after the other, each into its INPUT.out (19.0)\n"
                 "\n"
                 "  INPUT is an input file (.inp, keys below) or a structure file (19.0): .xyz (Å) or a\n"
                 "  Gaussian .gjf/.com/.gau, run as molecule = custom (default HF/aug-cc-pVDZ, job = energy):\n"
                 "    psilab water.xyz                                     HF energy, output on the terminal\n"
                 "    psilab water.xyz method='CCSD(T)' basis=aug-cc-pVTZ job=opt+freq -o   -> water.out\n"
                 "    psilab *.xyz method=MP2 job=opt                      every .xyz file of the folder\n"
                 "    psilab radical.xyz multiplicity=2 method=UHF         open shell (charge=, multiplicity=)\n"
                 "\n"
                 "  %s --estimate INPUT [key=value ...]  size of the calculation and time estimate, no calculation\n"
                 "  %s --example        print a commented example input file (HeH+)\n"
                 "  %s --example nuclear  print an example nuclear-motion input file (22.0)\n"
                 "  %s --version        print the version\n"
                 "  %s --help           this help\n\n"
                 "Input: one 'key = value' per line, '#' starts a comment. Keys:\n"
                 "  molecule = HeH+ | H2 | LiH | H2O | H2O2 | NH3 | custom   (default HeH+)\n"
                 "             15.0: the atoms H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar (ground-state\n"
                 "             multiplicity) and H2+ H3+ H5+ N2 O2 F2 Cl2 NO NO2 CO CO2 N2O H2CO CH4 (fixed geometries)\n"
                 "  geometry = <file>   molecule = custom: any structure, an .xyz file (Å; the comment line is\n"
                 "             the title) or a Gaussian .gjf/.com (charge and multiplicity from the file),\n"
                 "             relative to the input file; any element H–Lr, open shells with multiplicity\n"
                 "  title = <text>       name of a custom molecule (default: the file's title)\n"
                 "  * xyz <charge> <multiplicity>     19.0, as in ORCA: the atoms in the input file, one per line\n"
                 "    O  0.0000  0.0000  0.1173       (element x y z in Å), closed by a line '*'; means\n"
                 "    H  0.0000  0.7572 -0.4692       molecule = custom with this charge and multiplicity\n"
                 "    H  0.0000 -0.7572 -0.4692\n"
                 "    *\n"
                 "  * xyzfile <charge> <multiplicity> <file.xyz>   19.0, as in ORCA: the atoms from a structure file\n"
                 "  R = <number> [bohr | angstrom]       bond length (H2O, H2O2: O–H; NH3: N–H)\n"
                 "  charge = <n>   multiplicity = <2S+1>   (15.0: every molecule; default: the molecule's)\n"
                 "  reference = RHF | UHF | ROHF   (15.0; default RHF for closed-shell singlets, else UHF)\n"
                 "HeH+ (Hartree–Fock and Full CI, Szabo–Ostlund):\n"
                 "  basis = STO-1G | STO-2G | STO-3G      zeta_He = <number>   zeta_H = <number>\n"
                 "  exact_erf = yes | no   scan = yes | no   scan_step = <bohr>  (R = 0.8 ... 5.0 bohr)\n"
                 "  szabo_ostlund = no   HeH+ in STO-nG as a general molecule (all methods and jobs, 15.0)\n"
                 "General molecules (every method and job; HeH+ with any other basis set):\n"
                 "  basis = STO-1G | STO-2G | STO-3G   (15.0: all molecules; STO-1G H–Ca, STO-2G/3G H–Xe)\n"
                 "          6-31G | 6-31G* | 6-31+G | 6-31+G* | 6-31+G** | 6-31++G* | 6-31++G**   (6D)\n"
                 "          6-311G | 6-311G* | 6-311+G | 6-311+G* | 6-311+G** | 6-311++G* | 6-311++G** |\n"
                 "          6-311++G(2d,2p) | 6-311++G(3df,3pd) |\n"
                 "          def2-SVP | def2-mTZVP | def2-mTZVPP | def2-TZVP | def2-TZVPP | def2-TZVPPD |\n"
                 "          def2-QZVP | def2-QZVPP | def2-QZVPPD |\n"
                 "          aug-cc-pVDZ | aug-cc-pVTZ | aug-cc-pVQZ | aug-cc-pV5Z |\n"
                 "          aug-cc-pCVDZ | aug-cc-pCVTZ | aug-cc-pCVQZ | aug-cc-pCV5Z\n"
                 "          (also Gaussian's spellings 6-31G(d), 6-31+G(d,p), def2SVP; the 6-31G family\n"
                 "          has Cartesian d functions, 6D, as in Gaussian; all other sets are spherical)\n"
                 "          custom molecules: the def2 sets cover H–Rn (mTZVP, mTZVPP also Th–Lr) with\n"
                 "          effective core potentials from Rb on; 6-311G, 6-311G* include I\n"
                 "  angle = <degrees> (H2O: H–O–H, H2O2: O–O–H, NH3: H–N–H, 120 = planar)   threads = <n>\n"
                 "  R_OO = <number> [bohr | angstrom], dihedral = <degrees>   (H2O2: O–O bond, H–O–O–H)\n"
                 "  method = HF | RHF | UHF | ROHF | MP2 | MP3 | MP4 | MP5 | CISD | CISD(T) | CISDT | CISDTQ |\n"
                 "           QCISD | QCISD(T) | QCISDT | QCISDTQ | CCSD | CCSD(T) | CCSDT | CCSDTQ | FCI\n"
                 "           (all electrons correlated; prefixes U and RO select the reference: UMP2, ROCCSD(T))\n"
                 "  memory = <MB> (two-electron integrals kept in memory up to this size)\n"
                 "  fci_memory = <MB> (memory limit of Full CI and coupled cluster, default 4096)\n"
                 "  basis_dir = <folder with the .gbs files> (or environment variable PSILAB_BASIS_DIR)\n"
                 "  job = energy | opt | freq | opt+freq | opt+freq+sp | scan | ts | ts+freq | ts+freq+sp\n"
                 "        (default energy)\n"
                 "        opt: geometry optimisation (redundant internal coordinates, RFO, BFGS)\n"
                 "        freq: harmonic frequencies, IR intensities and thermochemistry\n"
                 "        opt+freq+sp, ts+freq+sp: then a single-point energy with sp_method / sp_basis\n"
                 "  sp_method = <method>   sp_basis = <basis>    (opt+freq+sp, ts+freq+sp; default: method, basis)\n"
                 "  temperature = <K> (298.15)   pressure = <number> [Pa | kPa | bar | atm] (101325 Pa)\n"
                 "  opt_convergence = normal | tight   max_steps = <n> (50, per point of a scan)\n"
                 "  scan = N1 N2 [N3 [N4]] S <steps> <step>   relaxed scan (job = scan), Gaussian ModRedundant line:\n"
                 "        2 atoms: bond (step in angstrom), 3: angle, 4: dihedral (degrees); nsteps + 1 points\n"
                 "  scan_start = <angstrom | degrees>     first point (default: the input geometry)\n"
                 "  ts_hessian = calcfc | calcall | recalcfc=N   eigen_test = yes | no   (ts, ts+freq: Gaussian Opt=TS)\n"
                 "  opt_hessian = auto | calcfc | calcall | recalcfc=N   (opt jobs, 15.0: Gaussian Opt=CalcFC,\n"
                 "        CalcAll, RecalcFC=N: the Hessian computed at step 1 and every N-th step, updated between)\n"
                 "  recalc_fc = N   (the same as opt_hessian / ts_hessian = recalcfc=N)\n"
                 "  opt_algorithm = rfo | gediis   (opt jobs, 15.0: RFO, or the hybrid GEDIIS/RFO of Gaussian)\n"
                 "  hessian = analytic | numerical   (frequencies and TS searches with FCI, CCSD, CCSD(T); 10.0)\n"
                 "  initial_hessian = auto | analytic | model   (minimisations: auto = analytic for HF, FCI, CCSD, CCSD(T)\n"
                 "                    where possible, else Lindh's model; model = Lindh; 10.0)\n\n"
                 "Nuclear motion (22.0): the rovibrational levels of H2(16)O, H2(17)O, H2(18)O, D2(16)O on the\n"
                 "CVRQD surfaces, as in the app (psilab --example nuclear prints a commented input):\n"
                 "  calculation = nuclear   (or molecule = H2(16)O | H2(17)O | H2(18)O | D2(16)O alone)\n"
                 "  molecule = H2(16)O | H2(17)O | H2(18)O | D2(16)O   coordinates = jacobi | radau\n"
                 "  ns1 = <n> (radial DVR points)   np = <n> (angular)   re = <bohr>   rmax = <bohr> (also angstrom)\n"
                 "  jmax = <J> (0: vibrational levels)   parity = 0 | 1 (J > 0: K from 0 or 1 to J)   levels = <n>\n"
                 "  eigenvectors = yes | no (no: the eigenvalues only)   expectation_values = yes | no\n"
                 "  eigensolver = lanczos | jacobi-davidson | davidson | inverse-power | rayleigh-quotient | dense\n"
                 "  linear_solver = cg | minres | chebyshev | bicg   preconditioner = ilu | jacobi | none\n"
                 "  gaussian_start = yes | no   eigenvalue_tolerance, residual_tolerance, linear_tolerance (E_h)\n"
                 "  max_iterations, max_linear_iterations, jd_inner_iterations, davidson_inner_iterations,\n"
                 "  subspace_size, threads\n"
                 "  The output file ends with the archive of the results: File > Open Output File... of the app\n"
                 "  shows it in every tab (Levels, Geometry, Matrix, Eigenvectors, Convergence, Output).\n\n"
                 "The output file is written while the job runs (tail -f OUTPUT follows it).\n",
                 prog, prog, prog, prog, prog, prog, prog, prog, prog, prog);
}

bool is_directory(const std::string& path) {
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

// <folder of the executable>/../share/psilab/basis, if it exists (installed psilab).
std::string installed_basis_directory() {
    std::string exe;
#if defined(__APPLE__)
    char buf[4096];
    uint32_t size = sizeof buf;
    if (_NSGetExecutablePath(buf, &size) == 0) exe = buf;
#elif defined(__linux__)
    char buf[4096];
    const ssize_t n = readlink("/proc/self/exe", buf, sizeof buf - 1);
    if (n > 0) exe.assign(buf, static_cast<std::size_t>(n));
#endif
    const auto slash = exe.find_last_of('/');
    if (slash == std::string::npos) return "";
    const std::string dir = exe.substr(0, slash) + "/../share/psilab/basis";
    return is_directory(dir) ? dir : "";
}

// 19.0: a structure file (.xyz, Gaussian .gjf/.com/.gau) is run directly, as molecule = custom.
std::string lower_extension(const std::string& path) {
    const auto slash = path.find_last_of('/');
    const auto dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
    std::string e = path.substr(dot);
    for (auto& c : e) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return e;
}

bool is_structure_file(const std::string& path) {
    const std::string e = lower_extension(path);
    return e == ".xyz" || e == ".gjf" || e == ".com" || e == ".gau";
}

bool is_input_file(const std::string& path) {
    return is_structure_file(path) || lower_extension(path) == ".inp";
}

// The input text of a run: the input file, or for a structure file "molecule = custom" and
// "geometry = <file>"; then the settings given on the command line (key=value), one per line.
// Errors in the added lines name the command-line option (or the structure file).
std::string read_input(const std::string& path, const std::vector<std::string>& options, psilab::Settings& settings,
                       int& status) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "error: cannot open %s\n", path.c_str());
        status = 2;
        return "";
    }
    const auto slash = path.find_last_of('/');
    const std::string folder = slash == std::string::npos ? std::string(".") : path.substr(0, slash);
    const bool structure = is_structure_file(path);
    std::string text;
    if (structure) {
        text = "molecule = custom\ngeometry = " + (slash == std::string::npos ? path : path.substr(slash + 1)) + "\n";
    } else {
        text.assign((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        if (!text.empty() && text.back() != '\n') text += '\n';
    }
    const int fileLines = static_cast<int>(std::count(text.begin(), text.end(), '\n'));
    for (const auto& o : options) text += o + "\n";
    try {
        settings = psilab::parse_input_text(text, folder);
    } catch (const std::exception& e) {
        std::string message = e.what();
        int line = 0;
        if (std::sscanf(message.c_str(), "line %d:", &line) == 1) {
            const auto colon = message.find(':');
            const std::string rest = message.substr(colon + 2);
            if (line > fileLines && line - fileLines <= static_cast<int>(options.size())) {
                message = "command-line option '" + options[static_cast<std::size_t>(line - fileLines - 1)] + "': " + rest;
                if (rest.find("is given twice") != std::string::npos)
                    message += structure ? " (molecule = custom and geometry come from the structure file)"
                                         : " (it is set in the input file, or in another option: an option adds a setting,"
                                           " it does not replace one)";
            }
            else if (structure)
                message = rest + " (molecule = custom and geometry come from the structure file; do not give them as options)";
        }
        std::fprintf(stderr, "error in %s: %s\n", path.c_str(), message.c_str());
        status = 2;
        return "";
    }
    for (const auto& w : settings.structureWarnings) std::fprintf(stderr, "note: %s\n", w.c_str());
    // Basis-set folder: input file > PSILAB_BASIS_DIR > installed copy > source tree.
    const char* env = std::getenv("PSILAB_BASIS_DIR");
    if (settings.input.basisDirectory.empty() && !(env && *env))
        settings.input.basisDirectory = installed_basis_directory();
    status = 0;
    return text;
}

std::string default_output_name(const std::string& input) {
    const auto slash = input.find_last_of('/');
    const auto dot = input.find_last_of('.');
    if (dot != std::string::npos && (slash == std::string::npos || dot > slash) && dot != 0 &&
        input.substr(dot) != ".out")
        return input.substr(0, dot) + ".out";
    return input + ".out";
}

// The arguments: input files (.inp) and structure files (.xyz, .gjf, .com, .gau), settings as
// key=value (also --key=value), -o, and at most one output file (with a single input).
struct Arguments {
    std::vector<std::string> inputs;
    std::vector<std::string> options;
    std::string output;
    bool defaultOutput = false;
};

bool parse_arguments(int first, int argc, char** argv, Arguments& a, const char* prog) {
    std::vector<std::string> positional;
    for (int i = first; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-o") {
            a.defaultOutput = true;
        } else if (arg.find('=') != std::string::npos && !is_input_file(arg)) {
            if (arg.rfind("--", 0) == 0) arg = arg.substr(2);
            a.options.push_back(arg);
        } else if (arg.size() > 1 && arg[0] == '-') {
            std::fprintf(stderr, "error: unknown option %s (settings are given as key=value; see %s --help)\n", arg.c_str(), prog);
            return false;
        } else {
            positional.push_back(arg);
        }
    }
    for (std::size_t i = 0; i < positional.size(); ++i) {
        // the first argument is always an input; later ones are inputs if they look like one
        if (i == 0 || is_input_file(positional[i])) {
            a.inputs.push_back(positional[i]);
        } else if (a.output.empty()) {
            a.output = positional[i];
        } else {
            std::fprintf(stderr, "error: too many arguments (see %s --help)\n", prog);
            return false;
        }
    }
    if (a.inputs.empty()) {
        std::fprintf(stderr, "error: no input file or structure file (see %s --help)\n", prog);
        return false;
    }
    if (!a.output.empty() && (a.inputs.size() > 1 || a.defaultOutput)) {
        std::fprintf(stderr, "error: an output file name goes with a single input and without -o (see %s --help)\n", prog);
        return false;
    }
    return true;
}

// 22.0: the input text of a nuclear-motion run (the file and the command-line settings), empty
// if the input is not a nuclear-motion input (or cannot be read: read_input reports that).
std::string nuclear_input(const std::string& path, const std::vector<std::string>& options) {
    if (is_structure_file(path)) return "";
    std::ifstream file(path, std::ios::binary);
    if (!file) return "";
    std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    if (!text.empty() && text.back() != '\n') text += '\n';
    for (const auto& o : options) text += o + "\n";
    return macqc::is_nuclear_input_text(text) ? text : "";
}

bool parse_nuclear(const std::string& path, const std::string& text, std::size_t options, macqc::NuclearMotionInput& in) {
    try {
        in = macqc::parse_nuclear_input_text(text);
    } catch (const std::exception& e) {
        std::string message = e.what();
        int line = 0;
        const int fileLines = static_cast<int>(std::count(text.begin(), text.end(), '\n')) - static_cast<int>(options);
        if (std::sscanf(message.c_str(), "line %d:", &line) == 1 && line > fileLines)
            message = "command-line option " + std::to_string(line - fileLines) + ": " + message.substr(message.find(':') + 2);
        std::fprintf(stderr, "error in %s: %s\n", path.c_str(), message.c_str());
        return false;
    }
    return true;
}

std::string nuclear_estimate_text(const macqc::NuclearMotionInput& in) {
    const macqc::NuclearMotionEstimate e = macqc::estimate_nuclear_motion(in);
    if (!e.refusal.empty()) return "Error: " + e.refusal + "\n";
    char buf[512];
    std::snprintf(buf, sizeof buf,
                  "Nuclear motion: %s, %s coordinates, J = %d, ns1 = %d, np = %d\n"
                  "  basis functions (dimension of H)  %lld\n"
                  "  non-zero elements of H            %lld\n"
                  "  memory                            %.1f MB\n",
                  macqc::isotopologue_name(in.isotopologue).c_str(), in.coordinates.c_str(), in.J, in.radial_points,
                  in.angular_points, e.dimension, e.nonzeros, e.memory_mb);
    return buf;
}

// Runs a nuclear-motion input (22.0): the output file as the app writes it, archive included.
int run_nuclear(const std::string& inputPath, const std::string& text, std::size_t options, const std::string& outputPath) {
    macqc::NuclearMotionInput in;
    if (!parse_nuclear(inputPath, text, options, in)) return 2;
    const macqc::NuclearMotionEstimate est = macqc::estimate_nuclear_motion(in);
    if (!est.refusal.empty()) {
        std::fprintf(stderr, "error in %s: %s\n", inputPath.c_str(), est.refusal.c_str());
        return 2;
    }
    const bool to_stdout = outputPath.empty();
    std::ofstream out;
    bool write_failed = false, written = false;
    auto sink = [&](const std::string& piece) {
        if (write_failed) return;
        written = true;
        if (to_stdout) {
            std::fwrite(piece.data(), 1, piece.size(), stdout);
            std::fflush(stdout);
            if (std::ferror(stdout)) write_failed = true;
            return;
        }
        if (!out.is_open()) {
            out.open(outputPath, std::ios::binary | std::ios::trunc);
            if (!out) {
                write_failed = true;
                std::fprintf(stderr, "error: cannot write %s\n", outputPath.c_str());
                return;
            }
        }
        out << piece;
        out.flush();
        if (!out) write_failed = true;
    };
    const auto t0 = std::chrono::steady_clock::now();
    macqc::NuclearMotionResult r;
    std::string error;
    try {
        r = macqc::run_nuclear_motion(in, sink, macqc::nuclear_archive_text);
    } catch (const std::exception& e) {
        error = e.what();
        if (written)
            sink(macqc::nuclear_termination_text(
                "error", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(), error));
    }
    if (out.is_open()) out.close();
    if (write_failed) {
        std::fprintf(stderr, "error: cannot write %s\n", to_stdout ? "the standard output" : outputPath.c_str());
        return 2;
    }
    if (!error.empty()) {
        std::fprintf(stderr, "error: %s\n", error.c_str());
        return 1;
    }
    FILE* summary = to_stdout ? stderr : stdout;
    if (!to_stdout) std::fprintf(summary, "Wrote %s\n", outputPath.c_str());
    std::fprintf(summary, "  %s, J = %d: %zu levels, %s; %lld basis functions\n", r.isotopologue_name.c_str(), in.J,
                 r.levels.size(), r.converged ? "converged" : "NOT converged", r.dimension);
    if (!r.levels.empty()) {
        const double e1 = r.levels.front().energy;
        std::fprintf(summary, "  E1 = %.10f E_h (%.4f cm^-1 above the minimum); E - E1 (cm^-1):", e1, e1 * 219474.6313632);
        for (std::size_t k = 1; k < r.levels.size() && k < 8; ++k)
            std::fprintf(summary, " %.4f", (r.levels[k].energy - e1) * 219474.6313632);
        std::fprintf(summary, r.levels.size() > 8 ? " ...\n" : "\n");
    }
    return r.converged ? 0 : 1;
}

// Runs one input; the output goes to `outputPath`, or to the standard output if it is empty.
int run_one(const std::string& inputPath, const std::vector<std::string>& options, const std::string& outputPath) {
    // 22.0: a nuclear-motion input (calculation = nuclear, or an isotopologue of water)
    const std::string nuclear = nuclear_input(inputPath, options);
    if (!nuclear.empty()) return run_nuclear(inputPath, nuclear, options.size(), outputPath);
    psilab::Settings settings;
    int status = 0;
    read_input(inputPath, options, settings, status);
    if (status) return status;
    for (const auto& w : psilab::estimate_warnings(settings)) std::fprintf(stderr, "warning: %s\n", w.c_str());

    // The output file is opened at the first piece of output and flushed after every piece, so
    // `tail -f` follows the job; HeH+ (Szabo–Ostlund) writes it at the end.
    const bool to_stdout = outputPath.empty();
    std::ofstream out;
    bool write_failed = false;
    auto sink = [&](const std::string& text) {
        if (write_failed) return;
        if (to_stdout) {
            std::fwrite(text.data(), 1, text.size(), stdout);
            std::fflush(stdout);
            if (std::ferror(stdout)) write_failed = true;
            return;
        }
        if (!out.is_open()) {
            out.open(outputPath, std::ios::binary | std::ios::trunc);
            if (!out) {
                write_failed = true;
                std::fprintf(stderr, "error: cannot write %s\n", outputPath.c_str());
                return;
            }
        }
        out << text;
        out.flush();
        if (!out) write_failed = true;
    };
    const bool streamed = !macqc::is_heh_input(settings.input);
    const psilab::RunResult result = psilab::run_to_text(settings, streamed ? macqc::OutputSink(sink) : macqc::OutputSink());
    if (result.text.empty()) {
        std::fprintf(stderr, "error: %s\n", result.error.c_str());
        return 1;
    }
    if (!streamed) sink(result.text);
    if (write_failed) {
        std::fprintf(stderr, "error: cannot write %s\n", to_stdout ? "the standard output" : outputPath.c_str());
        return 2;
    }
    if (out.is_open()) out.close();

    const auto& o = *result.output;
    // the summary: on the terminal (standard output) when the output went to a file, on the
    // standard error when the output file is the standard output
    FILE* summary = to_stdout ? stderr : stdout;
    if (!to_stdout) std::fprintf(summary, "Wrote %s\n", outputPath.c_str());
    if (o.job && o.input.job != "energy") {
        const auto& j = *o.job;
        if (j.optimization)
            std::fprintf(summary, "  %s: %s in %d steps, E = %.10f E_h\n",
                        j.optimization->transition_state ? "transition-state optimisation" : "optimisation",
                        j.optimization->converged ? "converged" : "NOT converged", j.optimization->steps,
                        j.optimization->final_energy);
        if (j.scan) {
            const auto& sc = *j.scan;
            const bool ang = sc.plan.angular();
            std::fprintf(summary, "  relaxed scan of %s: %zu points, %s\n", sc.plan.coordinate.definition.c_str(), sc.points.size(),
                        sc.message.c_str());
            if (!sc.points.empty()) {
                const auto& lo = sc.points[static_cast<std::size_t>(sc.lowest)];
                const auto& hi = sc.points[static_cast<std::size_t>(sc.highest)];
                std::fprintf(summary, "  lowest point %d (%.6f %s): E = %.10f E_h; highest point %d (%.6f %s): %.4f kcal/mol above\n",
                            lo.index + 1, ang ? lo.value : lo.value * psilab::kBohrToAngstrom, ang ? "deg" : "angstrom",
                            lo.energy, hi.index + 1, ang ? hi.value : hi.value * psilab::kBohrToAngstrom,
                            ang ? "deg" : "angstrom", (hi.energy - lo.energy) * 627.5094740631);
            }
        }
        if (j.harmonic) {
            std::fprintf(summary, "  frequencies (cm^-1):");
            for (const auto& m : j.harmonic->modes) std::fprintf(summary, " %.2f%s", std::fabs(m.wavenumber), m.imaginary ? "i" : "");
            std::fprintf(summary, "\n");
        }
        if (j.thermo)
            std::fprintf(summary, "  G(%.2f K) = %.10f E_h, H = %.10f E_h, ZPE = %.10f E_h\n", j.thermo->temperature,
                        j.thermo->total_g, j.thermo->total_h, j.thermo->total_zpe - j.thermo->electronic_energy);
        if (!j.final_energy_label.empty())
            std::fprintf(summary, "  final energy %s: %.10f E_h\n", j.final_energy_label.c_str(), j.final_energy);
    } else if (o.molecular && o.molecular->rhf) {
        std::fprintf(summary, "  %s total energy: %.10f E_h\n", o.molecular->reference.c_str(), o.molecular->rhf->total_energy);
        if (o.molecular->so)
            std::fprintf(summary, "  %s total energy: %.10f E_h\n", o.molecular->so->method.c_str(), o.molecular->so->total_energy);
        if (o.molecular->cc) std::fprintf(summary, "  %s total energy: %.10f E_h\n", o.molecular->cc->method.c_str(), o.molecular->cc->total_energy);
        if (o.molecular->fci) std::fprintf(summary, "  FCI total energy: %.10f E_h\n", o.molecular->fci->total_energy);
    } else if (o.hf) {
        std::fprintf(summary, "  HF  total energy: %.8f E_h\n", o.hf->total_energy);
        if (o.fci) std::fprintf(summary, "  FCI total energy: %.8f E_h\n", o.fci->total_energy);
    }
    if (!result.error.empty()) {
        const bool failed = o.job && o.job->termination != "normal";
        std::fprintf(stderr, "%s: %s\n", failed ? "error" : "warning", result.error.c_str());
        return 1;
    }
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) {
        usage(argv[0]);
        return argc < 2 ? 2 : 0;
    }
    if (!std::strcmp(argv[1], "--version")) {
        std::printf("\xCE\xA8Lab %s\n", psilab::kVersion);
        return 0;
    }
    if (!std::strcmp(argv[1], "--example")) {
        // 22.0: psilab --example nuclear (a nuclear-motion input)
        if (argc > 2 && (!std::strcmp(argv[2], "nuclear") || !std::strcmp(argv[2], "nuclear-motion"))) {
            std::cout << macqc::nuclear_example_input();
            return 0;
        }
        std::cout << psilab::example_input();
        return 0;
    }
    const bool estimate = !std::strcmp(argv[1], "--estimate");
    Arguments args;
    if (!parse_arguments(estimate ? 2 : 1, argc, argv, args, argv[0])) return 2;

    if (estimate) {
        int worst = 0;
        for (const auto& input : args.inputs) {
            if (args.inputs.size() > 1) std::printf("== %s\n", input.c_str());
            const std::string nuclear = nuclear_input(input, args.options);
            if (!nuclear.empty()) {
                macqc::NuclearMotionInput in;
                if (!parse_nuclear(input, nuclear, args.options.size(), in)) {
                    worst = std::max(worst, 2);
                    continue;
                }
                const std::string t = nuclear_estimate_text(in);
                std::fputs(t.c_str(), stdout);
                if (t.rfind("Error: ", 0) == 0) worst = std::max(worst, 1);
                continue;
            }
            psilab::Settings settings;
            int status = 0;
            read_input(input, args.options, settings, status);
            if (status) {
                worst = std::max(worst, status);
                continue;
            }
            const std::string t = psilab::estimate_text(settings);
            if (t.empty()) {
                std::printf("HeH+ (STO-nG): a two-function calculation, it takes milliseconds.\n");
                continue;
            }
            std::fputs(t.c_str(), stdout);
            if (t.rfind("Error: ", 0) == 0) worst = std::max(worst, 1);
        }
        return worst;
    }

    // One input: the standard output (15.0), OUTPUT or with -o INPUT.out. 19.0: several inputs
    // (e.g. psilab *.xyz method=MP2) run one after the other, each into its own INPUT.out.
    if (args.inputs.size() == 1) {
        const std::string& input = args.inputs.front();
        return run_one(input, args.options, args.defaultOutput ? default_output_name(input) : args.output);
    }
    int worst = 0;
    for (std::size_t i = 0; i < args.inputs.size(); ++i) {
        const std::string& input = args.inputs[i];
        std::printf("== %s (%zu of %zu)\n", input.c_str(), i + 1, args.inputs.size());
        std::fflush(stdout);
        worst = std::max(worst, run_one(input, args.options, default_output_name(input)));
    }
    return worst;
}
