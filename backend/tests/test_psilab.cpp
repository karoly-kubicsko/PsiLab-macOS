// Tests for the ΨLab terminal version (text input -> text output).
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "psilab_text.h"

static int g_checks = 0, g_failures = 0;

#define CHECK(cond)                                                            \
    do {                                                                       \
        ++g_checks;                                                            \
        if (!(cond)) {                                                         \
            ++g_failures;                                                      \
            std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond);        \
        }                                                                      \
    } while (0)

static bool contains(const std::string& s, const std::string& part) { return s.find(part) != std::string::npos; }

static bool rejects(const std::string& text, const std::string& expected_part) {
    try {
        psilab::parse_input_text(text);
    } catch (const std::exception& e) {
        if (!contains(e.what(), expected_part)) std::printf("  message was: %s\n", e.what());
        return contains(e.what(), expected_part);
    }
    return false;
}

int main() {
    // The example input is the Szabo-Ostlund default case.
    const psilab::Settings ex = psilab::parse_input_text(psilab::example_input());
    CHECK(ex.input.nGaussians == 3);
    CHECK(std::fabs(ex.input.R - 1.4632) < 1e-12);
    CHECK(std::fabs(ex.input.zeta1 - 2.0925) < 1e-12 && std::fabs(ex.input.zeta2 - 1.24) < 1e-12);
    CHECK(!ex.input.exactErf && ex.scan && std::fabs(ex.scanStep - 0.1) < 1e-12);

    // Units, case, aliases, comments.
    const psilab::Settings a = psilab::parse_input_text(
        "  BASIS = sto-2g  # comment\nr = 0.7743 Angstrom\nZeta_He: 2\nzeta-h = 1.2\nExact_Erf = YES\n");
    CHECK(a.input.nGaussians == 2);
    CHECK(std::fabs(a.input.R - 0.7743 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(a.input.zeta1 == 2.0 && a.input.zeta2 == 1.2 && a.input.exactErf && !a.scan);
    CHECK(std::fabs(psilab::parse_input_text("R = 1.0 \xC3\x85\n").input.R - 1.0 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(psilab::parse_input_text("R = 2\n").input.R == 2.0);  // default unit bohr

    // Errors carry the line number and say what is wrong.
    CHECK(rejects("basis = STO-3G\nbasis = STO-2G\n", "line 2"));
    CHECK(rejects("basis = 6-31G\n", "STO-1G, STO-2G or STO-3G"));
    CHECK(rejects("R = abc\n", "needs a number"));
    CHECK(rejects("R = 1.4 parsec\n", "unknown length unit"));
    CHECK(rejects("colour = blue\n", "unknown key"));
    CHECK(rejects("just some text\n", "expected 'key = value'"));
    CHECK(rejects("exact_erf = maybe\n", "yes or no"));
    CHECK(rejects("scan_step = 0\n", "scan_step"));
    CHECK(rejects("R = -1\n", "R must be"));
    CHECK(rejects("zeta_H = 0\n", "zeta2"));

    // Scan grid: 0.8 ... 5.0 inclusive.
    const auto grid = psilab::scan_grid(0.1);
    CHECK(grid.size() == 43 && std::fabs(grid.front() - 0.8) < 1e-12 && std::fabs(grid.back() - 5.0) < 1e-9);
    CHECK(std::fabs(psilab::scan_grid(0.25).back() - 5.0) < 1e-12);  // end point appended

    // Full run of the example: energies, report, curve and minima.
    const psilab::RunResult r = psilab::run_to_text(ex);
    CHECK(r.error.empty());
    CHECK(contains(r.text, "\xCE\xA8Lab \xE2\x80\x93 HeH+ (STO-3G)"));
    CHECK(contains(r.text, "R = 0.7743 \xC3\x85 (1.4632 bohr)"));
    CHECK(contains(r.text, "total energy:         -2.86066199 E_h"));
    CHECK(contains(r.text, "total energy:         -2.88071185 E_h"));
    CHECK(contains(r.text, "correlation energy:   -0.02004986 E_h"));
    CHECK(contains(r.text, "state 2: -2.07384042 E_h  <S^2> = 2.0000  triplet"));
    CHECK(contains(r.text, "STO-3G FOR ATOMIC NUMBERS  2.00 AND  1.00"));  // FORTRAN-style report
    CHECK(contains(r.text, "CORRELATION ENERGY ="));
    CHECK(contains(r.text, "POTENTIAL CURVE"));
    CHECK(contains(r.text, "HF minimum (parabolic fit): R = 0.7323 \xC3\x85 (1.3838 bohr), E = -2.86277416 E_h"));
    CHECK(contains(r.text, "FCI minimum (parabolic fit): R = 0.7397 \xC3\x85 (1.3978 bohr), E = -2.88223178 E_h"));
    CHECK(!contains(r.text, "failed point"));

    // General molecules: molecule, basis set, geometry defaults, keys per molecule.
    const psilab::Settings w = psilab::parse_input_text("molecule = water\nbasis = AUG-CC-PVTZ\n");
    CHECK(w.input.molecule == "H2O" && w.input.basis == "aug-cc-pVTZ");
    CHECK(std::fabs(w.input.R - 0.9572 / psilab::kBohrToAngstrom) < 1e-12 && w.input.angle == 104.52);
    const psilab::Settings w2 = psilab::parse_input_text(
        "molecule = H2O\nR = 1.0 angstrom\nangle = 100\nmethod = rhf\nthreads = 2\nmemory = 2 GB\n");
    CHECK(std::fabs(w2.input.R - 1.0 / psilab::kBohrToAngstrom) < 1e-12 && w2.input.angle == 100.0);
    CHECK(w2.input.threads == 2 && w2.input.memoryMB == 2048.0 && w2.input.basis == "aug-cc-pVDZ");
    CHECK(std::fabs(psilab::parse_input_text("molecule = LiH\n").input.R - 1.5949 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(rejects("molecule = CH4\n", "molecule must be HeH+, H2, LiH, H2O, H2O2, NH3 or custom"));
    CHECK(rejects("molecule = H2\nzeta_H = 1.2\n", "line 2: 'zeta_H' applies to HeH+ in an STO-nG basis (Szabo-Ostlund) only"));
    CHECK(rejects("molecule = H2\nangle = 100\n", "applies to H2O, H2O2 and NH3 only"));
    CHECK(rejects("angle = 100\n", "applies to H2O, H2O2 and NH3 only"));
    CHECK(rejects("method = RHF\n", "HeH+ always runs Hartree-Fock and Full CI"));
    CHECK(rejects("molecule = H2O\nbasis = STO-3G\n", "basis must be one of aug-cc-pVDZ"));
    CHECK(rejects("molecule = H2\nmethod = MP2\n", "method must be HF (RHF), CCSD, CCSD(T), CCSDT, CCSDTQ or FCI"));
    CHECK(psilab::parse_input_text("molecule = H2\nmethod = ccsd(t)\n").input.method == "CCSD(T)");
    CHECK(psilab::parse_input_text("molecule = H2\nmethod = HF\n").input.method == "RHF");
    CHECK(psilab::parse_input_text("molecule = LiH\nmethod = CCSDTQ\n").input.method == "CCSDTQ");
    CHECK(rejects("molecule = H2\nfci_memory = 1 PB\n", "fci_memory unit"));
    const psilab::Settings fs = psilab::parse_input_text("molecule = H2\nmethod = Full-CI\nfci_memory = 2 GB\n");
    CHECK(fs.input.method == "FCI" && fs.input.fciMemoryMB == 2048.0);
    CHECK(rejects("molecule = H2O\nangle = 10\n", "angle must be between 30 and 180"));
    CHECK(rejects("molecule = H2\nthreads = 1.5\n", "threads must be an integer"));
    CHECK(rejects("molecule = H2\nmemory = 3 TB\n", "memory unit"));

    // H2, LiH, H2O: the output file of ΨLab 5.0 (the same text the app exports), streamed.
    std::string streamed;
    int pieces = 0;
    const psilab::RunResult h2 = psilab::run_to_text(psilab::parse_input_text("molecule = H2\nbasis = aug-cc-pVDZ\n"),
                                                     [&](const std::string& t) {
                                                         streamed += t;
                                                         ++pieces;
                                                     });
    CHECK(h2.error.empty());
    CHECK(streamed == h2.text && pieces > 10);
    CHECK(contains(h2.text, "\xCE\xA8Lab 13.0"));
    CHECK(contains(h2.text, "Job          Energy (energy)"));
    CHECK(contains(h2.text, "    | molecule = H2"));  // the input file as given
    CHECK(contains(h2.text, "Hartree-Fock (RHF) total energy:  -1.1287937367 E_h"));
    CHECK(contains(h2.text, "total energy (RHF)               -1.1287937367"));
    CHECK(contains(h2.text, "18 basis functions"));
    CHECK(contains(h2.text, "H   aug-cc-pVDZ   (5s2p)           -> [3s2p]"));
    CHECK(contains(h2.text, "Normal termination of \xCE\xA8Lab 13.0"));
    CHECK(contains(psilab::estimate_text(psilab::parse_input_text("molecule = H2O\nbasis = aug-cc-pCVDZ\n")),
                   "aug-cc-pCVDZ is not defined for H"));
    // Full CI: H2 runs, H2O is refused with an explanation (in the output file and the error).
    const psilab::RunResult h2fci =
        psilab::run_to_text(psilab::parse_input_text("molecule = H2\nbasis = aug-cc-pVDZ\nmethod = FCI\n"));
    CHECK(h2fci.error.empty());
    CHECK(contains(h2fci.text, "Energy: Full CI/aug-cc-pVDZ"));
    CHECK(contains(h2fci.text, "Full CI total energy:  -1.1646244795 E_h"));
    CHECK(contains(h2fci.text, "324 determinants"));
    CHECK(contains(h2fci.text, "Natural orbital occupation numbers"));
    // Coupled cluster: CCSD(T) for H2O (PySCF: -76.2760053269 E_h).
    const psilab::RunResult h2occ =
        psilab::run_to_text(psilab::parse_input_text("molecule = H2O\nbasis = aug-cc-pVDZ\nmethod = CCSD(T)\n"));
    CHECK(h2occ.error.empty());
    CHECK(contains(h2occ.text, "Energy: CCSD(T)/aug-cc-pVDZ"));
    CHECK(contains(h2occ.text, "CCSD(T) total energy:  -76.2760053269 E_h"));
    CHECK(contains(h2occ.text, "(T) correction"));
    const psilab::RunResult h2ofci = psilab::run_to_text(psilab::parse_input_text("molecule = H2O\nmethod = FCI\n"));
    CHECK(contains(h2ofci.error, "impossible") && contains(h2ofci.text, "NOT POSSIBLE") &&
          contains(h2ofci.text, "Error termination"));
    CHECK(contains(psilab::estimate_text(psilab::parse_input_text("molecule = LiH\nbasis = aug-cc-pVQZ\nmethod = FCI\n")),
                   "NOT POSSIBLE"));
    const psilab::RunResult lih =
        psilab::run_to_text(psilab::parse_input_text("molecule = LiH\nbasis = aug-cc-pCV5Z\n"));
    CHECK(lih.text.empty() && contains(lih.error, "not available for Li"));

    // ΨLab 5.0 keys: job types and their options.
    const psilab::Settings j5 = psilab::parse_input_text(
        "molecule = H2O\njob = Opt+Freq+SP\nmethod = HF\nsp_method = ccsd(t)\nsp_basis = AUG-CC-PVTZ\n"
        "temperature = 350 K\npressure = 1 bar\nopt_convergence = tight\nmax_steps = 30\n");
    CHECK(j5.input.job == "opt+freq+sp" && j5.input.spMethod == "CCSD(T)" && j5.input.spBasis == "aug-cc-pVTZ");
    CHECK(j5.input.temperature == 350.0 && j5.input.pressure == 1.0e5);
    CHECK(j5.input.optConvergence == "tight" && j5.input.maxSteps == 30);
    CHECK(psilab::parse_input_text("molecule = H2\njob = freq\npressure = 1 atm\n").input.pressure == 101325.0);
    CHECK(psilab::parse_input_text("molecule = H2\njob = freq\npressure = 100 kPa\n").input.pressure == 1.0e5);
    CHECK(rejects("molecule = H2\njob = sweep\n", "job must be energy, opt, freq, opt+freq, opt+freq+sp, scan, ts, ts+freq or ts+freq+sp"));
    CHECK(rejects("job = opt\n", "applies to H2, LiH, H2O, H2O2, NH3 and HeH+ in a correlation-consistent basis set only"));
    CHECK(rejects("molecule = H2\njob = opt\nsp_method = CCSD\n", "'sp_method' applies to job = opt+freq+sp and ts+freq+sp only"));
    CHECK(rejects("molecule = H2\njob = opt\ntemperature = 300\n",
                  "applies to job = freq, opt+freq, opt+freq+sp, ts+freq and ts+freq+sp only"));
    CHECK(rejects("molecule = H2\njob = freq\nmax_steps = 10\n",
                  "applies to job = opt, opt+freq, opt+freq+sp, scan, ts, ts+freq and ts+freq+sp only"));
    CHECK(rejects("molecule = H2\njob = opt\nopt_convergence = loose\n", "opt_convergence must be normal or tight"));
    CHECK(rejects("molecule = H2\njob = opt\nmax_steps = 0\n", "max_steps must be an integer from 1 to 1000"));
    CHECK(rejects("molecule = H2\njob = freq\ntemperature = -3\n", "temperature must be between"));
    CHECK(rejects("molecule = H2\njob = freq\npressure = 3 psi\n", "pressure unit"));
    CHECK(rejects("molecule = H2\njob = opt+freq+sp\nsp_basis = 6-31G\n", "unknown basis set"));
    // estimates of jobs: steps, analytic or numerical, refusals
    const std::string est = psilab::estimate_text(psilab::parse_input_text("molecule = H2O\nmethod = CCSD(T)\njob = opt+freq\n"));
    CHECK(contains(est, "Geometry optimisation: analytic gradient"));
    CHECK(contains(est, "analytic Hessian: amplitude responses along the 3 internal directions"));
    // 10.0: hessian = numerical keeps the differences of analytic gradients; initial_hessian
    const std::string estn = psilab::estimate_text(
        psilab::parse_input_text("molecule = H2O\nmethod = CCSD(T)\njob = opt+freq\nhessian = numerical\n"));
    CHECK(contains(estn, "numerical Hessian: central differences of 6 analytic gradients"));
    const std::string esti = psilab::estimate_text(
        psilab::parse_input_text("molecule = H2O\nmethod = CCSD\njob = opt\ninitial_hessian = analytic\n"));
    CHECK(contains(esti, "analytic initial Hessian"));
    // the default: the analytic initial Hessian also for the correlated methods; "model": Lindh
    CHECK(contains(psilab::estimate_text(psilab::parse_input_text("molecule = H2O\nmethod = CCSD\njob = opt\n")),
                   "analytic initial Hessian"));
    CHECK(!contains(psilab::estimate_text(
                        psilab::parse_input_text("molecule = H2O\nmethod = CCSD\njob = opt\ninitial_hessian = model\n")),
                    "analytic initial Hessian"));
    CHECK(rejects("molecule = H2\njob = freq\nhessian = exact\n", "hessian must be analytic or numerical"));
    CHECK(rejects("molecule = H2\njob = energy\nhessian = analytic\n", "applies to job = freq"));
    CHECK(rejects("molecule = H2\njob = ts\ninitial_hessian = analytic\n", "applies to job = opt, opt+freq, opt+freq+sp and scan only"));
    const std::string refused = psilab::estimate_text(psilab::parse_input_text("molecule = H2O\nmethod = CCSDTQ\njob = freq\nfci_memory = 8 GB\n"));
    CHECK(refused.rfind("Error: ", 0) == 0 && contains(refused, "CCSDTQ frequencies for H2O is not possible"));
    // a short optimisation job with psilab
    const psilab::RunResult opt =
        psilab::run_to_text(psilab::parse_input_text("molecule = LiH\nmethod = HF\njob = opt\nR = 1.5 angstrom\n"));
    CHECK(opt.error.empty());
    CHECK(contains(opt.text, "Optimisation step 1") && contains(opt.text, "Stationary point found"));
    CHECK(contains(opt.text, "R(Li1-H2)"));
    CHECK(opt.output && opt.output->job && opt.output->job->optimization && opt.output->job->optimization->converged);

    // ΨLab 7.0: relaxed scan, the Gaussian ModRedundant line
    const psilab::Settings sc = psilab::parse_input_text(
        "molecule = H2O\njob = scan\nscan_start = 100\nscan = 2 1 3 S 8 10.0   # H-O-H\nmax_steps = 20\n");
    CHECK(sc.input.job == "scan" && sc.input.scan == "2 1 3 S 8 10.0");
    CHECK(sc.input.scanStart && *sc.input.scanStart == 100.0 && sc.input.maxSteps == 20);
    const psilab::Settings sb = psilab::parse_input_text("molecule = H2O\njob = scan\nscan = B 1 2 S 4 0.05\nscan_start = 1.8 bohr\n");
    CHECK(sb.input.scanStart && std::fabs(*sb.input.scanStart - 1.8 * psilab::kBohrToAngstrom) < 1e-12);
    CHECK(psilab::parse_input_text("molecule = H2O\njob = scan\nscan = 1 2 S 4 0.05\nscan_start = 1.0 angstrom\n").input.scanStart == 1.0);
    CHECK(rejects("molecule = H2O\njob = scan\nscan = 3 1 2 4 S 2 10\n", "atom 4 does not exist"));
    CHECK(rejects("molecule = H2O\njob = scan\nscan = 1 2 S 2\n", "ModRedundant"));
    CHECK(rejects("molecule = H2O\njob = scan\n", "job = scan needs 'scan ="));
    CHECK(rejects("molecule = H2O\njob = opt\nscan = 1 2 S 2 0.1\n", "'scan' applies to job = scan only"));
    CHECK(rejects("molecule = H2O\njob = scan\nscan = 1 2 S 2 0.1\nscan_start = 90 degrees\n", "degrees apply to angles"));
    CHECK(rejects("molecule = H2O\njob = scan\nscan = 2 1 3 S 2 10\nscan_start = 1 bohr\n", "bohr applies to a bond"));
    CHECK(rejects("scan_start = 1\n", "applies to H2, LiH, H2O, H2O2, NH3 and HeH+ in a correlation-consistent basis set only"));
    const std::string sest = psilab::estimate_text(psilab::parse_input_text("molecule = H2O\njob = scan\nscan = 2 1 3 S 8 5\n"));
    CHECK(contains(sest, "Relaxed scan: analytic gradient (derivative integrals); 9 points"));
    const psilab::RunResult srun = psilab::run_to_text(
        psilab::parse_input_text("molecule = LiH\nmethod = HF\njob = scan\nscan = 1 2 S 2 0.1\nscan_start = 1.5\n"));
    CHECK(srun.error.empty());
    CHECK(contains(srun.text, "Summary of the relaxed scan: HF/aug-cc-pVDZ"));
    CHECK(contains(srun.text, "scan             = 1 2 S 2 0.1   (bond, 3 points; step in angstrom)"));
    CHECK(srun.output && srun.output->job && srun.output->job->scan && srun.output->job->scan->points.size() == 3);

    // ΨLab 7.0 stage 7.2: H2O2 and NH3
    const psilab::Settings hp = psilab::parse_input_text(
        "molecule = HOOH\nR_OO = 1.45 angstrom\nR = 0.96 angstrom\nangle = 100\ndihedral = 180\n");
    CHECK(hp.input.molecule == "H2O2" && std::fabs(hp.input.ROO - 1.45 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(std::fabs(hp.input.R - 0.96 / psilab::kBohrToAngstrom) < 1e-12 && hp.input.angle == 100.0 && hp.input.dihedral == 180.0);
    const psilab::Settings hd = psilab::parse_input_text("molecule = H2O2\n");
    CHECK(hd.input.angle == 102.32 && hd.input.dihedral == 113.70 && std::fabs(hd.input.R - 0.9670 / psilab::kBohrToAngstrom) < 1e-12);
    const psilab::Settings am = psilab::parse_input_text("molecule = ammonia\nangle = 120\n");
    CHECK(am.input.molecule == "NH3" && am.input.angle == 120.0 && std::fabs(am.input.R - 1.012 / psilab::kBohrToAngstrom) < 1e-12);
    CHECK(psilab::parse_input_text("molecule = NH3\n").input.angle == 106.7);
    CHECK(rejects("molecule = H2O\nR_OO = 1.4 angstrom\n", "'R_OO' applies to H2O2 only"));
    CHECK(rejects("molecule = NH3\ndihedral = 100\n", "'dihedral' applies to H2O2 only"));
    CHECK(rejects("molecule = NH3\nangle = 130\n", "between 60 and 120"));
    const psilab::Settings hs = psilab::parse_input_text("molecule = H2O2\njob = scan\nscan = 3 1 2 4 S 18 10.0\nscan_start = 0\n");
    CHECK(hs.input.scan == "3 1 2 4 S 18 10.0" && hs.input.scanStart && *hs.input.scanStart == 0.0);
    CHECK(rejects("molecule = NH3\njob = scan\nscan = 2 1 3 5 S 2 10\n", "atom 5 does not exist"));

    // ΨLab 7.0 stage 7.3: transition states
    const psilab::Settings ts = psilab::parse_input_text(
        "molecule = NH3\njob = ts+freq\nangle = 116\nts_hessian = CalcAll\neigen_test = no\n");
    CHECK(ts.input.job == "ts+freq" && ts.input.tsHessian == "calcall" && !ts.input.eigenTest);
    const psilab::Settings tsd = psilab::parse_input_text("molecule = H2O2\njob = ts\ndihedral = 170\n");
    CHECK(tsd.input.job == "ts" && tsd.input.tsHessian == "calcfc" && tsd.input.eigenTest);
    CHECK(psilab::parse_input_text("molecule = H2O2\njob = Opt=TS\n").input.job == "ts");
    CHECK(rejects("molecule = NH3\njob = ts\nts_hessian = readfc\n", "ts_hessian must be calcfc or calcall"));
    CHECK(rejects("molecule = NH3\njob = opt\nts_hessian = calcfc\n", "'ts_hessian' applies to job = ts, ts+freq and ts+freq+sp only"));
    CHECK(rejects("molecule = NH3\njob = opt\neigen_test = no\n", "'eigen_test' applies to job = ts, ts+freq and ts+freq+sp only"));
    const psilab::Settings tsp = psilab::parse_input_text(
        "molecule = NH3\njob = ts+freq+sp\nangle = 119\nsp_method = CCSD(T)\nsp_basis = aug-cc-pVTZ\n");
    CHECK(tsp.input.job == "ts+freq+sp" && tsp.input.spMethod == "CCSD(T)" && tsp.input.spBasis == "aug-cc-pVTZ");
    const std::string tsest = psilab::estimate_text(psilab::parse_input_text("molecule = NH3\njob = ts\n"));
    CHECK(contains(tsest, "Transition-state optimisation: analytic gradient"));
    CHECK(contains(tsest, "Hessian computed at the first step (CalcFC)"));

    // ΨLab 7.0 stage 7.4: HeH+ with a correlation-consistent basis set is a general molecule
    const psilab::Settings hg = psilab::parse_input_text("molecule = HeH+\nbasis = aug-cc-pVTZ\nmethod = FCI\njob = opt+freq\n");
    CHECK(hg.input.hehGeneral && hg.input.basis == "aug-cc-pVTZ" && hg.input.method == "FCI" && hg.input.job == "opt+freq");
    CHECK(std::fabs(hg.input.R - 1.4632) < 1e-12);  // the Szabo-Ostlund bond length (0.7743 angstrom)
    CHECK(!psilab::parse_input_text("molecule = HeH+\nbasis = STO-2G\n").input.hehGeneral);
    CHECK(rejects("molecule = HeH+\nbasis = aug-cc-pVDZ\nzeta_He = 2\n", "'zeta_He' applies to HeH+ in an STO-nG basis"));
    CHECK(rejects("molecule = HeH+\nmethod = CCSD\n", "applies to H2, LiH, H2O, H2O2, NH3 and HeH+ in a correlation-consistent"));
    CHECK(rejects("molecule = HeH+\nbasis = 6-31G\n", "or a correlation-consistent set"));
    const psilab::RunResult hrun = psilab::run_to_text(
        psilab::parse_input_text("molecule = HeH+\nbasis = aug-cc-pVDZ\nmethod = CCSD\njob = scan\nscan = 1 2 S 2 0.1\n"));
    CHECK(hrun.error.empty());
    CHECK(contains(hrun.text, "HeH+, 2 electrons, charge 1"));
    CHECK(contains(hrun.text, "Summary of the relaxed scan: CCSD/aug-cc-pVDZ"));

    // Without scan: no curve section.
    psilab::Settings single = ex;
    single.scan = false;
    CHECK(!contains(psilab::run_to_text(single).text, "POTENTIAL CURVE"));

    std::printf("%d checks, %d failures\n", g_checks, g_failures);
    return g_failures == 0 ? 0 : 1;
}
