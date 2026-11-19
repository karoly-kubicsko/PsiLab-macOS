// Closed-shell restricted Hartree-Fock for general molecules in contracted
// spherical Gaussian basis sets (H2, LiH, H2O in ΨLab 3.0).
#pragma once

#include <array>
#include <cstddef>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gto_integrals.h"
#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

// Two-electron part of the Fock matrix, G(P) = J(P) - K(P)/2 for a total density P.
// The unique integrals (ij|kl) are kept in memory when they fit into the memory
// limit ("stored"); otherwise they are recomputed for every Fock build ("direct",
// with Schwarz and density screening).
class FockBuilder {
public:
    FockBuilder(const BasisSet& basis, int threads, double memory_limit_mb);

    Matrix build(const Matrix& P) const;

    // 15.0 (open shells): the Coulomb and exchange matrices of several symmetric densities in one
    // pass over the integrals, J(D)_pq = sum_rs (pq|rs) D_rs, K(D)_pr = sum_qs (pq|rs) D_qs.
    void build_jk(const std::vector<Matrix>& D, std::vector<Matrix>& J, std::vector<Matrix>& K) const;

    bool stored() const { return stored_; }
    std::size_t unique_integrals() const { return unique_; }
    std::size_t significant_shell_quartets() const { return quartets_; }
    double setup_seconds() const { return setup_seconds_; }
    int threads() const { return threads_; }

    // Packed unique (ij|kl), i >= j, k >= l, ij >= kl (stored mode only).
    const std::vector<double>& packed() const { return eri_; }
    static std::size_t pair_index(std::size_t i, std::size_t j) {
        return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i;
    }

private:
    Matrix build_stored(const Matrix& P) const;
    Matrix build_direct(const Matrix& P) const;
    void build_jk_stored(const std::vector<Matrix>& D, std::vector<Matrix>& J, std::vector<Matrix>& K) const;
    void build_jk_direct(const std::vector<Matrix>& D, std::vector<Matrix>& J, std::vector<Matrix>& K) const;

    const BasisSet& basis_;
    ShellPairList pairs_;
    int threads_ = 1;
    bool stored_ = false;
    std::size_t unique_ = 0, quartets_ = 0;
    double setup_seconds_ = 0.0;
    std::vector<double> eri_;
};

// Number of unique two-electron integrals for n basis functions, n^4/8 asymptotically.
std::size_t unique_eri_count(std::size_t n);

struct RHFOptions {
    int max_iter = 100;
    double energy_tol = 1.0e-10;    // |E_n - E_(n-1)| (E_h)
    double gradient_tol = 1.0e-7;   // max |X^T (FPS - SPF) X|
    double lindep_tol = 1.0e-7;     // overlap eigenvalues below this are removed
    int diis_size = 8;
    int threads = 0;                // 0: all hardware threads
    double memory_limit_mb = 1024;  // for storing the two-electron integrals
    // Optional starting orbitals (AO x MO coefficients of a nearby geometry, same basis):
    // the occupied ones are orthonormalised in the current overlap metric and give the
    // starting density instead of the initial guess (geometry optimisations, finite differences).
    const Matrix* guess_orbitals = nullptr;
    // 14.0: the initial guess without starting orbitals: "sad" (default), the superposition of
    // the densities of the neutral atoms, each from an SCF of the atom in its basis functions
    // (and its ECP) with the electrons of its open shell spread evenly over the degenerate
    // orbitals; or "core", the orbitals of the core Hamiltonian (earlier versions).
    std::string guess = "sad";
    // Optional uniform electric field F (atomic units, about the origin): adds F . r to the
    // core Hamiltonian of the electrons and -F . sum_A Z_A R_A to the nuclear energy, so that
    // E(F) = E(0) - mu . F + ... (finite-field dipole moments; energies only, no gradients).
    std::array<double, 3> field{0.0, 0.0, 0.0};
    // 15.0: the reference determinant. "RHF" (closed-shell singlet), "UHF" (different orbitals for
    // the alpha and the beta electrons) or "ROHF" (doubly and singly occupied orbitals of one set);
    // run_scf chooses the routine. The spin multiplicity 2S + 1 fixes n_alpha - n_beta = 2S.
    std::string reference = "RHF";
    int multiplicity = 1;
    // UHF: starting beta orbitals (with guess_orbitals for the alpha ones); nullptr: the alpha ones.
    const Matrix* guess_orbitals_b = nullptr;
    // Open shells: a level shift (E_h) of the virtual orbitals in the first iterations, which
    // stabilises the aufbau of atoms with degenerate open shells; 0: automatic.
    double level_shift = 0.0;
};

struct RHFIteration {
    int iteration = 0;
    double energy = 0.0;      // total energy (E_h)
    double delta_e = 0.0;
    double gradient = 0.0;    // max |orbital gradient|
    bool diis = false;
};

struct RHFResult {
    bool converged = false;
    int iterations = 0;
    double electronic_energy = 0.0;
    double nuclear_repulsion = 0.0;
    double total_energy = 0.0;
    double kinetic_energy = 0.0;           // Tr(P T)
    double nuclear_attraction_energy = 0.0;  // Tr(P V)
    double two_electron_energy = 0.0;      // Tr(P G)/2
    double virial_ratio = 0.0;             // -V/T
    int n_basis = 0;
    int n_mo = 0;          // after removal of linear dependencies
    int n_occupied = 0;
    double smallest_overlap_eigenvalue = 0.0;
    int removed_functions = 0;
    std::vector<double> orbital_energies;  // ascending, n_mo values
    Matrix C;              // AO x MO coefficients
    Matrix P;              // total density
    Matrix H;              // core Hamiltonian
    Matrix S;              // overlap
    std::vector<RHFIteration> history;
    bool eri_stored = false;
    std::size_t unique_integrals = 0;
    int threads = 1;
    double integral_seconds = 0.0;  // one-electron + stored ERIs (or Schwarz setup)
    double scf_seconds = 0.0;
    bool guess_from_orbitals = false;  // the SCF started from RHFOptions::guess_orbitals
    std::string guess;                 // 14.0: "orbitals", "sad" or "core"

    // 15.0, open shells. For RHF: reference "RHF", n_alpha = n_beta = n_occupied, Cb empty.
    // UHF: C and orbital_energies are the alpha orbitals, Cb and orbital_energies_b the beta ones;
    // n_occupied = n_beta. ROHF: one set C, the first n_beta orbitals doubly and the next
    // n_alpha - n_beta singly occupied; orbital_energies are the eigenvalues of the Roothaan
    // effective Fock matrix (PySCF's choice of the blocks). P = Pa + Pb in every case.
    std::string reference = "RHF";
    int multiplicity = 1;
    int n_alpha = 0, n_beta = 0;
    Matrix Cb;
    std::vector<double> orbital_energies_b;
    Matrix Pa, Pb;                      // alpha and beta densities (open shells; RHF: P/2 each)
    Matrix Fa, Fb;                      // alpha and beta Fock matrices (AO) at convergence (open shells)
    double s2 = 0.0;                    // <S^2> of the determinant
    double s2_exact = 0.0;              // S(S+1)
    bool open_shell() const { return reference != "RHF"; }
    // The beta orbitals: Cb for UHF, C otherwise.
    const Matrix& beta_orbitals() const { return Cb.rows() ? Cb : C; }
    const std::vector<double>& beta_orbital_energies() const {
        return orbital_energies_b.empty() ? orbital_energies : orbital_energies_b;
    }
};

// Progress callback: called after each SCF iteration (may be empty).
using RHFProgress = std::function<void(const RHFIteration&)>;

// Throws std::runtime_error for impossible input (e.g. odd electron count).
// Non-convergence is reported through RHFResult::converged.
RHFResult run_rhf(const Molecule& mol, const BasisSet& basis, const RHFOptions& options,
                  const RHFProgress& progress = RHFProgress());

// 15.0: unrestricted Hartree-Fock (Pople-Nesbet equations; DIIS on both Fock matrices with one
// set of coefficients) and restricted open-shell Hartree-Fock (one set of orbitals; the Roothaan
// effective Fock matrix with PySCF's blocks: (F_a + F_b)/2 in the closed-closed, open-open,
// virtual-virtual and closed-virtual blocks, F_b closed-open, F_a open-virtual). Both start from
// the SAD density split evenly into alpha and beta, from given orbitals, or from the core
// Hamiltonian, and report <S^2>. Throws for impossible electron counts and multiplicities.
RHFResult run_uhf(const Molecule& mol, const BasisSet& basis, const RHFOptions& options,
                  const RHFProgress& progress = RHFProgress());
RHFResult run_rohf(const Molecule& mol, const BasisSet& basis, const RHFOptions& options,
                   const RHFProgress& progress = RHFProgress());
// RHF, UHF or ROHF according to options.reference.
RHFResult run_scf(const Molecule& mol, const BasisSet& basis, const RHFOptions& options,
                  const RHFProgress& progress = RHFProgress());

// Numbers of alpha and beta electrons for `electrons` and the multiplicity 2S + 1; throws with an
// explanation if they do not fit (e.g. an even number of electrons in a doublet).
std::pair<int, int> alpha_beta_electrons(int electrons, int multiplicity);
// The superposition of atomic densities of the molecule (the SAD guess of the SCF), total density.
Matrix sad_guess_density(const Molecule& mol, const BasisSet& basis, int threads, double memory_mb);
// Canonical name of a reference ("RHF", "UHF", "ROHF") for any spelling, or "".
std::string canonical_reference_name(const std::string& name);

}  // namespace macqc
