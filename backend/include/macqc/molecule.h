// Molecules of ΨLab: H2, LiH and H2O (3.0) built from bond length (and angle); H2O2 (7.0) built
// from R(O-O), R(O-H), the O-O-H angle and the H-O-O-H dihedral angle; NH3 (7.0) built from
// R(N-H) and the H-N-H angle (C3v; 120 degrees: planar, D3h). ΨLab 9.0: any closed-shell molecule
// of any elements from Cartesian coordinates ("custom", read from .xyz or .gjf files); the basis sets
// cover H-Kr (14.0: the def2 sets with their ECPs up to Rn and Th-Lr; see kMaxAtomicNumber).
#pragma once

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace macqc {

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrInAngstrom = 0.529177210544;

// 14.0: a semi-local effective core potential (ECP) of an atom, as in Gaussian's input:
//   U(r) = U_L(r) + sum_{l < L} sum_m |lm> [U_l(r) - U_L(r)] <lm|,
// each radial part a sum of terms c r^(n-2) exp(-a r^2) (r from the atom; ΨLab supports n = 2,
// the Gaussian terms of the Stuttgart-Cologne potentials of the def2 sets).
struct EcpTerm {
    int power = 2;            // n
    double exponent = 0.0;    // a (bohr^-2)
    double coefficient = 0.0; // c (E_h)
};

struct AtomEcp {
    int core_electrons = 0;                     // electrons replaced by the potential
    int local_l = 0;                            // L of the local part U_L
    std::vector<EcpTerm> local;                 // U_L
    std::vector<std::vector<EcpTerm>> projected;  // [l] = U_l - U_L, l = 0 .. L-1
};

struct Atom {
    int Z = 0;
    std::string symbol;
    std::array<double, 3> position{};  // bohr
    // 14.0: the ECP of the atom in the basis set of the calculation (nullptr: all electrons);
    // set by attach_ecps (basis.h) when the basis is built.
    std::shared_ptr<const AtomEcp> ecp;

    // The nuclear charge the electrons see: Z, less the core electrons of an ECP.
    int nuclear_charge() const { return ecp ? Z - ecp->core_electrons : Z; }
    int core_electrons() const { return ecp ? ecp->core_electrons : 0; }
};

struct Molecule {
    std::string name;  // "H2", "LiH", "H2O", "H2O2", "NH3", "HeH+"; "custom" (9.0)
    std::vector<Atom> atoms;
    int charge = 0;
    // 9.0: a molecule given by Cartesian coordinates (a structure file). Its atoms may come in any
    // order, so the code paths written for the built-in molecules (Z-matrix coordinates of the
    // numerical derivatives, symmetry-adapted displacements) are replaced by general ones.
    bool custom = false;
    std::string title;  // custom: the name of the structure file's comment/title line (may contain LaTeX)

    // Electrons treated explicitly (14.0: without the core electrons of ECPs) and the nuclear
    // repulsion of the nuclear charges the electrons see (Z less the ECP core electrons).
    int electrons() const;
    int core_electrons() const;  // 14.0: replaced by ECPs
    double nuclear_repulsion() const;
    // Molecular formula in Hill order ("CO2", "CH4", "H2O"; with the charge: "NH4+", "SO4 2-").
    std::string formula() const;
    // Name for the output: the formula for custom molecules, else `name`.
    std::string display_name() const;
};

// Element symbol of Z = 1..118 ("H", "He", "Li", ..., "Og"); throws otherwise.
std::string element_symbol(int Z);
// Atomic number of an element symbol (case-insensitive) or of a number given as text ("8"), 0 if
// unknown. 9.0: every element, Z = 1 .. 118.
int atomic_number(const std::string& symbol);
// Largest atomic number (118). A structure may contain any element; which elements a calculation
// can use depends on the basis set: the all-electron aug-cc-pVXZ sets of the Basis Set Exchange
// cover H-Ar and Sc-Kr (Sc-Zn up to aug-cc-pVQZ: their 5Z sets need i functions), aug-cc-pCVXZ
// Li-Ar (aug-cc-pCV5Z: B-Ar); 14.0: 6-31G, 6-31G* and the def2 sets H-Kr, the diffuse 6-31G sets
// H-Ar, the 6-311G family H-Ca and Ga-Kr (diffuse ones up to Ca). Heavier elements need effective
// core potentials (aug-cc-pVXZ-PP, def2-ECP), which ΨLab does not implement; such molecules are
// refused with that explanation.
constexpr int kMaxAtomicNumber = 118;

// Canonical molecule name ("H2", "LiH", "H2O", "H2O2", "NH3", "HeH+"; "custom") for any spelling, or "".
std::string canonical_molecule_name(const std::string& name);

// 15.0: molecules and atoms of the molecule menu that are given by their Cartesian coordinates:
// the atoms H-Ar in their ground states and H2+, H3+, H5+, N2, O2, F2, Cl2, NO, NO2 at their
// equilibrium structures. They run as custom molecules with these atoms, charge and multiplicity.
struct MoleculePreset {
    std::string key;        // "O2", "NO2", "H5+", "N" (atoms: the element symbol)
    std::string title;      // name in LaTeX notation, e.g. "O$_2$", "H$_5^+$"
    int charge = 0;
    int multiplicity = 1;
    std::string state;      // ground state, e.g. "X 3Sigma_g-", "4S"
    std::vector<std::pair<std::string, std::array<double, 3>>> atoms;  // symbol, position in angstrom
    std::string source;     // where the structure comes from
    bool atom = false;      // a single atom
};
const std::vector<MoleculePreset>& molecule_presets();
// The preset of a name in any case ("o2", "H5+", "h5plus", "Cl", "chlorine atom"), or nullptr.
const MoleculePreset* find_molecule_preset(const std::string& name);
// A molecule of the general path: built-in, custom or a preset.
inline bool general_molecule_name(const std::string& name) {
    return !canonical_molecule_name(name).empty() || find_molecule_preset(name) != nullptr;
}

// Default geometries (experimental structures), in bohr / degrees:
// H2 r_e = 0.74144 Å; LiH r_e = 1.5949 Å; H2O r(OH) = 0.9572 Å, angle 104.52°;
// H2O2 (gas phase) r(OO) = 1.4556 Å, r(OH) = 0.9670 Å, angle O-O-H 102.32°, dihedral
// H-O-O-H 113.70°; NH3 r(NH) = 1.012 Å, angle H-N-H 106.7°. default_bond_length() is the X-H
// bond (H2O2: O-H).
double default_bond_length(const std::string& molecule);  // bohr
constexpr double kDefaultWaterAngle = 104.52;              // degrees
constexpr double kDefaultPeroxideOOAngstrom = 1.4556;
constexpr double kDefaultPeroxideAngle = 102.32;           // O-O-H, degrees
constexpr double kDefaultPeroxideDihedral = 113.70;        // H-O-O-H, degrees
constexpr double kDefaultAmmoniaAngle = 106.7;             // H-N-H, degrees
double default_oo_bond_length();                            // bohr

// Default of the `angle` input of a molecule: H-O-H (H2O), O-O-H (H2O2), H-N-H (NH3), else 104.52.
double default_angle(const std::string& molecule);

// Builds the molecule. H2: atoms at z = -R/2, +R/2. LiH: Li at the origin, H at z = R. HeH+ (7.0, for
// the general path with correlation-consistent basis sets): He at the origin, H at z = R, charge +1.
// H2O: O at the origin, the hydrogens in the xz plane at distance R, symmetric about
// the z axis, enclosing `angle_degrees`.
// H2O2: atoms O1, O2, H3 (bonded to O1), H4 (bonded to O2); the O-O bond along x, centred at
// the origin, the C2 axis along z; R = R(O-H), `angle_degrees` = angle O-O-H, `r_oo_bohr` =
// R(O-O) (<= 0: default), `dihedral_degrees` = dihedral angle H3-O1-O2-H4 (0: cis, 180: trans).
// NH3: N at the origin, the C3 axis along z, the hydrogens at z > 0 (planar for 120 degrees),
// H2 in the xz plane; R = R(N-H), `angle_degrees` = angle H-N-H (60-120 degrees).
// Throws std::invalid_argument on bad input.
Molecule make_molecule(const std::string& name, double bond_length_bohr, double angle_degrees,
                       double r_oo_bohr = 0.0, double dihedral_degrees = kDefaultPeroxideDihedral);

// 9.0: a custom molecule from atoms (Z and positions in bohr; the symbols are filled in), its
// charge and title. Throws std::invalid_argument for no atoms, unknown elements, non-finite
// coordinates or two atoms closer than 0.3 bohr (0.16 Å). Open shells are not rejected here
// (prepare_molecular explains that RHF needs an even number of electrons).
Molecule make_custom_molecule(const std::vector<Atom>& atoms, int charge, const std::string& title);

}  // namespace macqc
