// ΨLab for the terminal (Linux, macOS): reads a plain-text input file, writes a
// plain-text output file.  Usage: psilab INPUT > OUTPUT (15.0) or psilab INPUT OUTPUT
// (see psilab --help)
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

#include "psilab_text.h"

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
                 "  %s INPUT > OUTPUT   run INPUT, write the output file to the standard output (15.0)\n"
                 "  %s INPUT OUTPUT     run INPUT, write the output file to OUTPUT\n"
                 "  %s INPUT -o         run INPUT, write INPUT with the extension replaced by .out\n"
                 "                      (the summary and messages go to the standard error)\n"
                 "  %s --estimate INPUT size of the calculation and time estimate, no calculation\n"
                 "  %s --example        print a commented example input file (HeH+)\n"
                 "  %s --version        print the version\n"
                 "  %s --help           this help\n\n"
                 "Input: one 'key = value' per line, '#' starts a comment. Keys:\n"
                 "  molecule = HeH+ | H2 | LiH | H2O | H2O2 | NH3 | custom   (default HeH+)\n"
                 "             15.0: the atoms H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar (ground-state\n"
                 "             multiplicity) and H2+ H3+ H5+ N2 O2 F2 Cl2 NO NO2 CO CO2 N2O H2CO CH4 (fixed geometries)\n"
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
                 "The output file is written while the job runs (tail -f OUTPUT follows it).\n",
                 prog, prog, prog, prog, prog, prog, prog);
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

std::string read_input(const std::string& path, psilab::Settings& settings, int& status) {
    std::ifstream file(path, std::ios::binary);
    if (!file) {
        std::fprintf(stderr, "error: cannot open %s\n", path.c_str());
        status = 2;
        return "";
    }
    const std::string text((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    try {
        const auto slash = path.find_last_of('/');
        settings = psilab::parse_input_text(text, slash == std::string::npos ? std::string(".") : path.substr(0, slash));
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error in %s: %s\n", path.c_str(), e.what());
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
        std::cout << psilab::example_input();
        return 0;
    }
    if (!std::strcmp(argv[1], "--estimate")) {
        if (argc != 3) {
            std::fprintf(stderr, "usage: %s --estimate INPUT\n", argv[0]);
            return 2;
        }
        psilab::Settings settings;
        int status = 0;
        read_input(argv[2], settings, status);
        if (status) return status;
        const std::string t = psilab::estimate_text(settings);
        if (t.empty()) {
            std::printf("HeH+ (STO-nG): a two-function calculation, it takes milliseconds.\n");
            return 0;
        }
        std::fputs(t.c_str(), stdout);
        return t.rfind("Error: ", 0) == 0 ? 1 : 0;
    }
    if (argc > 3) {
        std::fprintf(stderr, "error: too many arguments (see %s --help)\n", argv[0]);
        return 2;
    }

    const std::string inputPath = argv[1];
    psilab::Settings settings;
    int status = 0;
    read_input(inputPath, settings, status);
    if (status) return status;
    for (const auto& w : psilab::estimate_warnings(settings)) std::fprintf(stderr, "warning: %s\n", w.c_str());

    // The output file is opened at the first piece of output and flushed after every piece, so
    // `tail -f` follows the job (H2, LiH, H2O); HeH+ writes it at the end.
    // 15.0: without an output file the output goes to the standard output (psilab INPUT > OUTPUT);
    // "-o": the default name INPUT.out
    const bool to_stdout = argc == 2;
    const std::string outputPath =
        argc == 3 ? (std::strcmp(argv[2], "-o") == 0 ? default_output_name(inputPath) : std::string(argv[2])) : "";
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
