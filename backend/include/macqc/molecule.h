// Molecules of ΨLab: H2, LiH and H2O (3.0) built from bond length (and angle); H2O2 (7.0) built
// from R(O-O), R(O-H), the O-O-H angle and the H-O-O-H dihedral angle; NH3 (7.0) built from
// R(N-H) and the H-N-H angle (C3v; 120 degrees: planar, D3h).
#pragma once

#include <array>
#include <string>
#include <vector>

namespace macqc {

// 1 bohr in angstrom (CODATA 2022).
constexpr double kBohrInAngstrom = 0.529177210544;

struct Atom {
    int Z = 0;
    std::string symbol;
    std::array<double, 3> position{};  // bohr
};

struct Molecule {
    std::string name;  // "H2", "LiH", "H2O", "H2O2", "NH3"
    std::vector<Atom> atoms;
    int charge = 0;

    int electrons() const;
    double nuclear_repulsion() const;
};

// Element symbol of Z = 1..10 ("H", "He", "Li", ..., "Ne"); throws otherwise.
std::string element_symbol(int Z);
// Atomic number of an element symbol (case-insensitive), 0 if unknown.
int atomic_number(const std::string& symbol);

// Canonical molecule name ("H2", "LiH", "H2O", "H2O2", "NH3", "HeH+") for any spelling, or "".
std::string canonical_molecule_name(const std::string& name);

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

}  // namespace macqc
