// Gaussian basis sets: Gaussian94 (.gbs) parser, the basis-set library of ΨLab and
// contracted shells (spherical harmonics; 14.0: Cartesian d for the 6-31G family) placed on
// the atoms of a molecule.
#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "macqc/molecule.h"

namespace macqc {

// Highest angular momentum the integral engine supports (h functions, l = 5;
// the tables are sized for one more).
constexpr int kMaxAngularMomentum = 6;

// One contracted shell as written in a basis-set file (before normalisation).
struct ShellData {
    int l = 0;
    std::vector<double> exponents;
    std::vector<double> coefficients;  // for normalised primitives (file convention)
};

struct ElementBasis {
    int Z = 0;
    std::vector<ShellData> shells;
    std::shared_ptr<const AtomEcp> ecp;  // 14.0: the effective core potential, if the file has one
};

struct BasisFile {
    std::string header;                   // the '!' comment lines at the top
    std::map<int, ElementBasis> elements;  // by atomic number
};

// Parses Gaussian94 text (as written by basissetexchange.org). SP shells are split
// into an S and a P shell; a scale factor other than 1 multiplies the exponents
// by its square. 14.0: effective core potentials ("<El> 0", "<El>-ECP L n_core", then
// the local part U_L and U_l - U_L for l = 0 .. L-1, each "<label>", the number of terms
// and lines "n exponent coefficient"; only n = 2 is supported). Throws
// std::runtime_error("line N: ...") on malformed input.
BasisFile parse_gaussian94(const std::string& text);

// A basis set of ΨLab's library (14.0).
struct BasisSetInfo {
    const char* name;    // canonical display name, e.g. "6-311++G(2d,2p)"
    const char* file;    // Gaussian94 file in the basis directory, e.g. "6-311ppg_2d_2p_.gbs"
    const char* family;  // "Pople 6-31G", "Pople 6-311G", "Karlsruhe def2", "Dunning aug-cc-pVXZ",
                         // "Dunning aug-cc-pCVXZ"
    bool cartesian_d;    // d shells as six Cartesian functions (6D): the 6-31G family
};

// The library in display order: the 6-31G family (7 sets), the 6-311G family (9), def2 (9),
// aug-cc-pVXZ (4) and aug-cc-pCVXZ (4).
const std::vector<BasisSetInfo>& basis_set_library();

// The names of the basis sets offered by ΨLab, in display order.
const std::vector<std::string>& supported_basis_sets();

// The library entry for any spelling of a name, or nullptr. Spellings: any case and blanks;
// Gaussian's polarisation notation for the Pople sets, 6-31G(d) = 6-31G*, 6-31+G(d,p) =
// 6-31+G**; def2SVP = def2-SVP.
const BasisSetInfo* basis_set_info(const std::string& name);

// Canonical display name ("aug-cc-pVDZ") for any spelling ("AUG-CC-PVDZ"), or "".
std::string canonical_basis_name(const std::string& name);

// File name in the basis directory, e.g. "aug-cc-pvdz.gbs", "6-31pgss.gbs" (the names of
// Psi4's library: + -> p, * -> s, ( , ) -> _).
std::string basis_file_name(const std::string& canonical_name);

// The kind of functions of a basis set for reports: "Cartesian d functions (6D), spherical f
// functions (7F)" for the 6-31G family, else "spherical harmonics (5D, 7F)".
std::string basis_function_types(const std::string& name);

// Default basis directory: $PSILAB_BASIS_DIR if set, else the compile-time default
// (MACQC_BASIS_DIR, the repository's basis/ folder), else "basis".
std::string default_basis_directory();

// A contracted shell of Gaussians on an atom: real solid harmonics (spherical, "pure",
// 2l+1 functions m = -l .. l) or, 14.0, Cartesian functions x^i y^j z^k (i + j + k = l,
// (l+1)(l+2)/2 functions in the order xx, xy, xz, yy, yz, zz for d; each normalised), as
// the d shells of the 6-31G family are used (Gaussian's 6D).
struct Shell {
    int l = 0;
    int atom = 0;                    // index into Molecule::atoms
    std::array<double, 3> center{};  // bohr
    std::vector<double> exponents;
    // Contraction coefficients that already include the normalisation of each
    // primitive (that of x^l exp(-a r^2)) and of the contracted function.
    std::vector<double> coefficients;
    int first_function = 0;  // index of the first function of the shell
    bool pure = true;        // false: Cartesian functions

    int size() const { return pure ? 2 * l + 1 : (l + 1) * (l + 2) / 2; }
    int cartesian_size() const { return (l + 1) * (l + 2) / 2; }
};

// Basis-set assignment of one atom, for reports.
struct AtomBasisInfo {
    std::string symbol;
    std::string basis;        // basis set actually used for this atom
    std::string contraction;  // e.g. "[4s3p2d]"
    std::string primitives;   // e.g. "(10s5p2d)" (distinct exponents per l)
    int functions = 0;
    std::string note;         // e.g. why another basis set was used for H
    int ecp_electrons = 0;    // 14.0: core electrons replaced by the atom's ECP (0: all electrons)
    std::string ecp;          // 14.0: "ECP28" for a potential with 28 core electrons ("" without)
};

struct BasisSet {
    std::string name;  // requested basis set (canonical)
    bool cartesian_d = false;  // 14.0: its d shells are Cartesian (6D)
    std::vector<Shell> shells;
    std::vector<AtomBasisInfo> atoms;
    int n_functions = 0;
    int max_l = 0;
    int n_primitives = 0;  // sum over shells of exponents x (2l+1)
    std::string source;    // header line(s) of the file(s) read
    // 14.0: the ECP of every atom (nullptr: all electrons), in the order of the molecule's atoms
    // (empty for a basis set built by hand: no ECPs)
    std::vector<std::shared_ptr<const AtomEcp>> ecps;
    int ecp_electrons = 0;  // core electrons replaced by ECPs, all atoms
};

// 14.0: the molecule with the ECPs of `basis` on its atoms (Atom::ecp; atoms without an ECP in
// the basis lose one they had). Every calculation with `basis` uses this molecule: its nuclear
// charges are Z less the ECP core electrons, and its nuclear attraction includes the ECPs.
Molecule attach_ecps(const Molecule& mol, const BasisSet& basis);

// Throws std::logic_error unless the ECPs on the atoms of `mol` are those of `basis`.
void check_ecps(const Molecule& mol, const BasisSet& basis);

// Builds the basis for the molecule from the files in `directory`. For H, the
// core-valence sets aug-cc-pCVXZ are not defined (H has no core); the matching
// aug-cc-pVXZ functions are used and a note is recorded. The d shells of the 6-31G family
// are Cartesian (Shell::pure false), all other shells spherical. Throws std::runtime_error
// with an explanation if a basis set is not available for an element.
BasisSet build_basis(const Molecule& mol, const std::string& basis_name, const std::string& directory);

// Same, from already parsed files: `files` maps canonical basis names to contents.
BasisSet build_basis(const Molecule& mol, const std::string& basis_name,
                     const std::map<std::string, BasisFile>& files);

// Normalised contraction coefficients (primitive and contraction normalisation).
std::vector<double> normalized_coefficients(int l, const std::vector<double>& exponents,
                                            const std::vector<double>& coefficients);

// Letter of an angular momentum (0 -> 's', 5 -> 'h').
char angular_letter(int l);

}  // namespace macqc
