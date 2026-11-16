// Point groups of ΨLab's molecules (7.0; 9.0: any molecule): the symmetry operations that map a
// molecule onto itself (rotations C2 ... C8, reflections, the inversion and S3 ... S8, about axes
// and planes through the centre of mass), the point group they form (C1, Cs, Ci, Cn, Cnv, Cnh, Sn,
// Dn, Dnh, Dnd, T, Td, Th, O, Oh, I, Ih; linear molecules C∞v and D∞h), the rotational symmetry
// number (the number of proper rotations, E included; 1 for C∞v and 2 for D∞h) and the symmetry
// species of normal modes from their characters (reduction formula, so a degenerate pair gives E
// and an accidental near-degeneracy a sum such as "A1+E"). Character tables: C1, Cs, Ci, C2, C2v,
// C2h, C3, C3v, D3h (7.0) and D2, D2h, D3, D3d, Td (9.0); linear molecules Sigma and Pi species.
// Other groups are named, with their symmetry number, but their modes get no species.
//
// Conventions (Mulliken): for a planar C2v molecule the molecular plane is sigma_v'(yz), so the
// in-plane antisymmetric modes are B2; in D3h sigma_h is perpendicular to the C3 axis.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "macqc/molecule.h"

namespace macqc {

struct SymmetryOperation {
    std::string kind;                                  // "E", "C2", "C3", "sigma", "i", "S3"
    std::string symmetry_class;                        // class of the point group, e.g. "sigma_v'"
    std::array<std::array<double, 3>, 3> R{};          // about the centre of mass
    std::vector<int> permutation;                      // atom a is moved onto atom permutation[a]
};

struct PointGroup {
    std::string name = "C1";
    std::vector<SymmetryOperation> operations;
    int symmetry_number = 1;
    bool has_table = true;   // false: a group without a character table here (no species)
    bool consistent = true;  // the operations found form the named group (its order)
    bool linear = false;     // 9.0: C∞v or D∞h (axis below; operations E and, for D∞h, i)
    std::array<double, 3> axis{};
};

// `tolerance`: largest allowed distance (bohr) between an image and an atom. The default (0)
// tries 2e-3 bohr (optimised geometries, also transition states, are symmetric only to about
// 1e-4 bohr) and tightens to 5e-4, 1e-4 and 2e-5 bohr while the operations found do not form
// one of the groups above.
PointGroup point_group(const Molecule& mol, double tolerance = 0.0);

// True if the group has degenerate species (E, T, Pi): near-degenerate modes are then classified
// together.
bool has_degenerate_species(const PointGroup& group);

// Species of a set of (degenerate) modes given as mass-weighted, orthonormal 3N vectors
// (sqrt(m) times the Cartesian displacement); "" if the group has no table here.
std::string mode_species(const PointGroup& group, const std::vector<std::vector<double>>& modes);

}  // namespace macqc
