// Molecular orbitals for display (ΨLab 8.0): the atoms, the basis and sets of orbitals as AO
// coefficients (RHF canonical orbitals; Full-CI natural orbitals; the two orbitals of the
// Szabo-Ostlund HeH+ calculation), their values at points, on grids and in planes, and their
// isosurfaces psi = +c and psi = -c (marching tetrahedra).
#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/molecule.h"

namespace macqc {

struct RHFResult;
struct FCIResultGeneral;

struct OrbitalInfo {
    int number = 0;            // 1-based number in its full set (energy order or occupation order)
    double energy = 0.0;       // E_h (canonical orbitals)
    bool has_energy = false;   // false for natural orbitals
    double occupation = 0.0;   // 2 or 0 (RHF), the natural occupation (Full-CI)
    std::string label;         // "HOMO", "HOMO-1", "LUMO", "LUMO+2" or ""
};

struct OrbitalSet {
    std::string kind;    // "canonical" or "natural"
    std::string title;   // e.g. "RHF canonical orbitals"
    std::vector<OrbitalInfo> orbitals;
    std::vector<std::vector<double>> coefficients;  // per orbital: the AO coefficients
};

struct OrbitalData {
    std::vector<Atom> atoms;    // bohr
    std::string basis_name;
    int n_basis = 0;
    std::vector<Shell> shells;  // spherical or (14.0) Cartesian shells; coefficients include all normalisation
    std::vector<OrbitalSet> sets;
};

// RHF canonical orbitals (all occupied ones and up to `virtuals` virtual ones) and, with a Full
// CI result that has natural orbitals, as many Full-CI natural orbitals (largest occupation
// first).
OrbitalData make_orbital_data(const Molecule& mol, const BasisSet& basis, const RHFResult& rhf,
                              const FCIResultGeneral* fci = nullptr, int virtuals = 10);

// The two orbitals of the Szabo-Ostlund HeH+ calculation: STO-nG s functions on atom A (He, at
// the origin) and atom B (H, at z = R), with the coefficients C (2 x 2, columns = orbitals).
OrbitalData make_sto_orbital_data(int n_gaussians, double R, double zeta1, double zeta2, double ZA, double ZB,
                                  const std::vector<double>& C, const std::vector<double>& energies);

// Values of all basis functions (spherical or Cartesian, in basis order) at the point r.
void basis_function_values(const std::vector<Shell>& shells, int n_basis, const std::array<double, 3>& r,
                           double* values);

// Value at r of the orbital with AO coefficients c.
double orbital_value(const OrbitalData& d, const std::vector<double>& c, const std::array<double, 3>& r);

// The points origin + i a + j b + k c for i < n[0] (fastest), j < n[1], k < n[2].
struct GridSpec {
    std::array<double, 3> origin{0.0, 0.0, 0.0};
    std::array<std::array<double, 3>, 3> axes{{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}}};
    std::array<int, 3> n{1, 1, 1};

    std::size_t size() const {
        return static_cast<std::size_t>(n[0]) * static_cast<std::size_t>(n[1]) * static_cast<std::size_t>(n[2]);
    }
    std::size_t index(int i, int j, int k) const {
        return (static_cast<std::size_t>(k) * static_cast<std::size_t>(n[1]) + static_cast<std::size_t>(j)) *
                   static_cast<std::size_t>(n[0]) +
               static_cast<std::size_t>(i);
    }
    std::array<double, 3> point(double i, double j, double k) const {
        std::array<double, 3> p{};
        for (int d = 0; d < 3; ++d) p[d] = origin[d] + i * axes[0][d] + j * axes[1][d] + k * axes[2][d];
        return p;
    }
};

// Orbital values on a grid (order of GridSpec::index); threads = 0: all hardware threads.
std::vector<double> orbital_grid(const OrbitalData& d, const std::vector<double>& c, const GridSpec& g,
                                 int threads = 0);

struct Mesh {
    std::vector<float> positions;  // 3 per vertex (bohr)
    std::vector<float> normals;    // 3 per vertex, unit vectors pointing out of the enclosed region
    std::vector<int> indices;      // 3 per triangle, counter-clockwise seen from outside
    std::size_t vertex_count() const { return positions.size() / 3; }
    std::size_t triangle_count() const { return indices.size() / 3; }
};

// The surface f = level around the region f >= level, for f sampled on a grid with axes
// along x, y and z (marching tetrahedra: six tetrahedra per cube around its main diagonal, so
// that neighbouring cubes split their common faces alike). Vertices on the same grid edge are
// shared; the normals are -grad f (central differences, interpolated along the edge).
Mesh isosurface(const std::vector<double>& f, const GridSpec& g, double level);

struct OrbitalSurfaces {
    Mesh positive;          // psi = +isovalue
    Mesh negative;          // psi = -isovalue
    GridSpec grid;          // the grid of the surfaces
    double max_abs = 0.0;   // largest |psi| on the grid
};

// Isosurfaces of an orbital: a coarse grid (spacing 0.5 bohr, the atoms +- 12 bohr) finds where
// |psi| reaches the isovalue; a grid of at most `max_points` points covers that region with a
// margin of 1 bohr (spacing at least 0.04 bohr).
OrbitalSurfaces orbital_surfaces(const OrbitalData& d, const std::vector<double>& c, double isovalue,
                                 int max_points = 400000, int threads = 0);

// The "orbitals" object of the output JSON (api.cpp) and back.
std::string orbitals_to_json(const OrbitalData& d);
OrbitalData orbitals_from_json(const std::string& text);

}  // namespace macqc
