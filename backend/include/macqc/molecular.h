// ΨLab 3.0 driver for general molecules (H2, LiH, H2O; 7.0: H2O2, NH3): basis set, size and
// time estimate, warnings and refusals, restricted Hartree-Fock, text report.
#pragma once

#include <optional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/cc.h"
#include "macqc/ci.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"
#include "macqc/so_methods.h"

namespace macqc {

struct MolecularInput {
    std::string molecule = "H2O";  // ... or (9.0) "custom": the atoms below
    double R = 0.0;              // bond length in bohr (H2O, H2O2: O-H; NH3: N-H); <= 0 means the default
    double angle = kDefaultWaterAngle;  // degrees: H-O-H (H2O), O-O-H (H2O2), H-N-H (NH3)
    double r_oo = 0.0;           // H2O2: O-O bond length in bohr; <= 0 means the default
    double dihedral = kDefaultPeroxideDihedral;  // H2O2: H-O-O-H dihedral angle in degrees
    // 9.0, molecule "custom": the atoms (Z, positions in bohr), the charge, the spin multiplicity
    // (RHF: 1) and the title (name of the molecule; may contain LaTeX, e.g. "CO$_2$")
    std::vector<Atom> atoms;
    int charge = 0;
    int multiplicity = 1;
    std::string title;
    // 15.0: charge and multiplicity of every molecule (built-in molecules and presets take their own
    // unless these are set), the reference ("RHF", "UHF", "ROHF"; "" = RHF for closed-shell singlets,
    // UHF otherwise, ROHF orbitals for an open-shell Full CI) and the preset of the molecule menu
    // (atoms H-Ar, H2+, H3+, H5+, N2, O2, F2, Cl2, NO, NO2: find_molecule_preset).
    bool charge_given = false;
    bool multiplicity_given = false;
    std::string reference;
    std::string preset;
    std::string basis = "aug-cc-pVDZ";
    // "RHF" (shown as HF), "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ" or "FCI"; the correlated
    // methods follow RHF and correlate all electrons.
    std::string method = "RHF";
    std::string basis_directory;  // "" -> default_basis_directory()
    int threads = 0;              // 0 -> all hardware threads
    double memory_mb = 1024.0;    // limit for keeping the two-electron integrals in memory
    double fci_memory_mb = 4096.0;  // hard memory limit of a Full CI or coupled-cluster calculation
    int max_iter = 100;
    double energy_tol = 1.0e-10;
    double gradient_tol = 1.0e-7;
    int print_level = 2;          // 0 silent, 1 summary, 2 full
};

// Size and cost of a calculation, known before it runs.
struct CostEstimate {
    int basis_functions = 0;
    int shells = 0;
    int max_l = 0;
    std::size_t unique_eris = 0;
    double eri_memory_mb = 0.0;
    bool eri_stored = false;
    int threads = 1;
    double rhf_seconds = 0.0;    // rough estimate of the wall time
    std::optional<FCISize> fci;  // method FCI only
    std::optional<CCSize> cc;    // coupled-cluster methods only
    // 15.0: the spin-orbital engine (MPn, CI, QCI, open-shell coupled cluster)
    struct SOCost {
        std::string method;      // label with the reference, e.g. "UMP4", "QCISD(T)"
        double amplitudes = 0.0; // stored amplitudes (dense spin blocks)
        double memory_bytes = 0.0;
        double seconds = 0.0;
        bool feasible = true;
        std::string reason;
    };
    std::optional<SOCost> so;
    double total_seconds = 0.0;  // RHF (+ Full CI or coupled cluster)
};

struct MolecularResult {
    MolecularInput input;
    Molecule molecule;
    BasisSet basis;
    CostEstimate estimate;
    std::vector<std::string> notes;     // e.g. which basis set is used for H
    std::vector<std::string> warnings;  // e.g. long calculation, removed functions
    std::string refusal;                // non-empty: the calculation is not possible (reason)
    std::optional<RHFResult> rhf;
    std::optional<FCIResultGeneral> fci;
    std::optional<CCResult> cc;
    std::optional<SOResult> so;         // 15.0
    std::string reference = "RHF";      // 15.0: the reference used (effective)
    std::string report;
};

// Canonical method name ("RHF", "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ", "FCI"; 15.0: "MP2" ... "MP5",
// "CISD", "CISD(T)", "CISDT", "CISDTQ", "QCISD", "QCISD(T)", "QCISDT", "QCISDTQ") of a name
// given in any case, with or without spaces and hyphens ("HF" and "SCF" mean RHF,
// "Full CI" means FCI; 15.0: "UHF" and "ROHF" also mean the SCF, "UMP2" MP2 and so on, the prefix
// giving the reference, see parse_method_spec); "" if unknown. "RHF" stands for the SCF method.
std::string canonical_method_name(const std::string& method);

// 15.0: a method name with Gaussian's reference prefixes: "UHF" -> {RHF (the SCF), UHF},
// "ROMP2" -> {MP2, ROHF}, "UCCSD(T)" -> {CCSD(T), UHF}, "HF" -> {RHF, ""}, "RHF" -> {RHF, RHF}.
struct MethodSpec {
    std::string method;     // canonical, "" if unknown
    std::string reference;  // "", "RHF", "UHF" or "ROHF"
};
MethodSpec parse_method_spec(const std::string& text);
// The reference a calculation uses: the requested one, else RHF for a singlet, UHF for open shells
// (ROHF orbitals for Full CI). Throws std::runtime_error if RHF is requested for an open shell.
std::string effective_reference(const std::string& canonical_method, const std::string& requested, int multiplicity);
// "MP2", "UMP2", "ROMP2", "UHF", "ROHF", "RHF": the method with the prefix of its reference.
std::string method_label(const std::string& canonical_method, const std::string& reference);
// True for the coupled-cluster methods; sets *method to the CC method.
bool coupled_cluster_method(const std::string& canonical, CCMethod* method = nullptr);

// Validates the input; returns "" or a message (unknown molecule/basis/method, bad geometry).
std::string validate_molecular(const MolecularInput& in);

// Builds molecule and basis and estimates the cost without computing anything.
// Throws std::runtime_error with an explanation if the input cannot be set up (e.g. basis
// set not defined for an element). A Full CI that is too large is not an exception here:
// result.refusal explains why it is not possible (run_molecular then refuses to run).
MolecularResult prepare_molecular(const MolecularInput& in);

// Wall-time estimate (seconds) for an RHF calculation in this basis: a cost model of the
// integral work per shell quartet, times the expected number of SCF iterations when the
// integrals are recomputed (integral-direct). Typically within a factor of 2.
double estimate_rhf_seconds(const BasisSet& basis, bool eri_stored, int threads, int iterations);

// Typical number of SCF iterations for the molecule near its equilibrium geometry.
int expected_scf_iterations(const std::string& molecule);

// Runs the calculation (after prepare_molecular). Throws std::runtime_error on errors,
// including a refused calculation (result.refusal).
void run_molecular(MolecularResult& result, const RHFProgress& progress = RHFProgress(),
                   const FCIProgress& fci_progress = FCIProgress(), const CCProgress& cc_progress = CCProgress());

// The Full CI and coupled-cluster sections of the report (also used by the job output).
std::string fci_report_section(const FCIResultGeneral& f, double rhf_energy, int level);
std::string cc_report_section(const CCResult& c, int level);
// 15.0: the section of a method of the spin-orbital engine (energies by order or component,
// amplitude norms, T1 diagnostic, convergence).
std::string so_report_section(const SOResult& r, int level);

// Report text (also stored in result.report by run_molecular).
std::string molecular_report(const MolecularResult& result);

// Estimated-time threshold above which a warning is issued (seconds).
constexpr double kLongCalculationSeconds = 30.0;

}  // namespace macqc
