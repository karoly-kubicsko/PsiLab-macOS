// Redundant internal coordinates (ΨLab 5.0): bond stretches, bends and linear bends, the
// Wilson B matrix, the generalised inverse of G = B B', internal gradients and Hessians, and
// the iterative back-transformation of internal steps to Cartesian coordinates
// (Pulay & Fogarasi, J. Chem. Phys. 96, 2856 (1992); Peng, Ayala, Schlegel & Frisch,
// J. Comput. Chem. 17, 49 (1996)).
//
// Primitives: every pair of atoms closer than 1.3 times the sum of their covalent radii is a
// bond; every pair of bonds at a common atom is a bend, replaced by two linear-bend
// components (in two perpendicular planes through the A...C axis) when the angle exceeds 175
// degrees; every chain A-B-C-D of three bonds is a dihedral (torsion), left out while one of
// its two bends is above 175 degrees (7.0). Units: bohr for stretches, radians for bends and
// dihedrals (IUPAC sign, in (-pi, pi]; differences are taken modulo 2 pi); the linear-bend
// components are the projections (u + v) . e_k of the unit vectors from the central atom
// (dimensionless, close to the deviation from 180 degrees in radians).
//
// Coordinates that must be present (the constrained or scanned coordinate of a relaxed scan)
// can be passed to the constructor: they are added if the connectivity does not contain them,
// and a required bend stays a bend above 175 degrees.
#pragma once

#include <array>
#include <string>
#include <vector>

#include "macqc/linalg.h"
#include "macqc/molecule.h"

namespace macqc {

enum class InternalType { Stretch, Bend, LinearBend, Dihedral };

struct InternalCoordinate {
    InternalType type = InternalType::Stretch;
    std::vector<int> atoms;          // 2 (stretch), 3 (bend: A, centre, C) or 4 (dihedral A-B-C-D)
    int component = 0;               // linear bend: 0 or 1
    std::array<double, 3> axis{};    // linear bend: the fixed direction e_k
    std::string name() const;        // "R", "A", "La", "Lb", "D"
    std::string definition(const Molecule& mol) const;  // e.g. "R(O1-H2)", "A(H2-O1-H3)", "D(H3-O1-O2-H4)"
    // Same type and atoms (also in reverse order); linear bends: also the same component.
    bool same_as(const InternalCoordinate& other) const;
    bool angular() const { return type != InternalType::Stretch; }
};

// Value of one coordinate at a geometry (bohr or rad; dihedrals in (-pi, pi]).
double internal_value(const Molecule& mol, const InternalCoordinate& c);
// Its row of the Wilson B matrix, dq/dx (3 n_atoms).
std::vector<double> wilson_b_row(const Molecule& mol, const InternalCoordinate& c);

class InternalCoordinates {
public:
    InternalCoordinates() = default;
    // Builds the primitives from the connectivity of `mol`, plus `required` coordinates.
    explicit InternalCoordinates(const Molecule& mol, const std::vector<InternalCoordinate>& required = {});

    const std::vector<InternalCoordinate>& coordinates() const { return coords_; }
    std::size_t size() const { return coords_.size(); }
    // Index of a coordinate (see InternalCoordinate::same_as), or -1.
    int find(const InternalCoordinate& c) const;
    // True if any bend became linear (> 175 deg) or a linear bend bent back (< 170 deg), or a
    // dihedral lost or regained its definition: the set must be rebuilt for this geometry.
    bool needs_rebuild(const Molecule& mol) const;

    std::vector<double> values(const Molecule& mol) const;
    Matrix wilson_b(const Molecule& mol) const;  // n_q x 3 n_atoms
    // Second derivatives d2q_i/dx dx' contracted with c_i: sum_i c_i d2q_i (3N x 3N),
    // by central differences of the B matrix.
    Matrix b_derivative_contracted(const Molecule& mol, const std::vector<double>& c) const;
    // q1 - q0, dihedral differences reduced to (-pi, pi].
    std::vector<double> difference(const std::vector<double>& q1, const std::vector<double>& q0) const;

private:
    friend std::vector<double> wilson_b_row(const Molecule& mol, const InternalCoordinate& c);
    std::vector<InternalCoordinate> coords_;
    std::vector<InternalCoordinate> required_;
    std::vector<InternalCoordinate> omitted_dihedrals_;  // undefined (near-linear) when built
};

// An angle reduced to (-pi, pi].
double wrap_angle(double a);

// Generalised inverse of a symmetric positive semidefinite matrix (eigenvalues below
// tol * largest are treated as zero); `rank` receives the number of kept eigenvalues.
Matrix generalized_inverse(const Matrix& G, double tol = 1.0e-8, int* rank = nullptr);

// Internal-coordinate gradient g_q = G^- B g_x.
std::vector<double> internal_gradient(const Matrix& B, const Matrix& Ginv, const std::vector<double>& gx);

struct BackTransform {
    Molecule molecule;
    int iterations = 0;
    double residual = 0.0;  // largest |q(x) - q_target| at the end
    bool converged = false;
};

// Cartesian geometry whose internal coordinates are q0 + dq (iterative: x += B' G^- (q_t - q(x)),
// first iteration exact to first order). Stops when the Cartesian change is below 1e-10 bohr.
BackTransform back_transform(const InternalCoordinates& ic, const Molecule& mol, const std::vector<double>& dq,
                             int max_iter = 50);

// Covalent radius (bohr) used for the connectivity (Pyykkö-Atsumi single-bond radii).
double covalent_radius(int Z);

}  // namespace macqc
