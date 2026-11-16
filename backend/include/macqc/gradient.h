// Energies and nuclear gradients at arbitrary Cartesian geometries for every method of
// ΨLab (5.0): the building block of geometry optimisations and numerical Hessians.
//
// Analytic gradients: RHF (derivative integrals, energy-weighted density, exact term for
// removed near-linear dependencies); CCSD and CCSD(T) (Lambda equations, relaxed densities,
// Z-vector); Full CI (relaxed densities). CCSDT and CCSDTQ, and the correlated methods when
// near-linear dependencies were removed or the analytic gradient does not fit into memory,
// use central differences of energies along (symmetry-adapted) internal coordinates; the
// result says which was used and why.
#pragma once

#include <array>
#include <functional>
#include <optional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/cc.h"
#include "macqc/cc_gradient.h"
#include "macqc/ci.h"
#include "macqc/correlated_gradient.h"
#include "macqc/molecule.h"
#include "macqc/rhf.h"

namespace macqc {

struct ElectronicOptions {
    int threads = 0;                 // 0: all cores
    double memory_mb = 1024.0;       // AO two-electron integrals kept in memory up to this size
    double fci_memory_mb = 4096.0;   // memory limit of Full CI, coupled cluster and correlated gradients
    int scf_max_iter = 100;
    double scf_energy_tol = 1.0e-10;
    double scf_gradient_tol = 1.0e-9;
    double lindep_tol = 1.0e-7;
    int cc_max_iter = 150;
    double cc_amplitude_tol = 1.0e-9;
    double cc_energy_tol = 1.0e-10;
    double fci_residual_tol = 1.0e-8;
    double numerical_step = 0.005;   // bohr (stretches) and rad (bends) of numerical gradients
    std::array<double, 3> field{0.0, 0.0, 0.0};  // uniform electric field (a.u.), energies only
};

// Callbacks for progress and for the text output (each may be empty).
struct ElectronicCallbacks {
    RHFProgress scf;
    CCProgress cc;
    FCIProgress fci;
    LambdaProgress lambda;
    CPHFCallback cphf;                               // Z-vector / CPHF iterations
    std::function<void(const std::string&)> note;    // one line of text (e.g. a numerical step)
};

// Energy of one method at one geometry.
struct SinglePoint {
    Molecule molecule;
    std::string method;           // canonical: RHF, CCSD, CCSD(T), CCSDT, CCSDTQ, FCI
    double energy = 0.0;          // total energy of the method (E_h)
    RHFResult rhf;
    std::optional<CCResult> cc;
    std::optional<FCIResultGeneral> fci;
    std::vector<std::string> warnings;
    bool converged = false;       // SCF and the correlated method
    double seconds = 0.0;
};

struct GradientResult {
    SinglePoint point;            // the calculation at the geometry itself
    std::vector<double> gradient; // 3 n_atoms (E_h/bohr)
    bool analytic = false;
    std::string description;      // how it was computed (for the output and the UI)
    std::string numerical_reason; // why numerical (empty if analytic)
    int energy_evaluations = 0;   // numerical: energies computed (besides the reference)
    std::vector<std::string> internal_coordinates;  // numerical: the coordinates displaced
    std::array<double, 3> dipole{};  // e bohr: SCF dipole (RHF) or relaxed (correlated, analytic)
    bool has_dipole = false;
    std::optional<RHFGradientResult> rhf_gradient;
    std::optional<CorrelatedGradientResult> correlated;
    bool lambda_converged = false;
    int lambda_iterations = 0;
    double lambda_residual = 0.0;
    std::vector<double> l1, l2;   // 10.0: the CCSD / CCSD(T) multipliers (for the analytic Hessian)
    double seconds = 0.0;
};

// True for the methods with an analytic gradient: RHF, CCSD, CCSD(T), FCI.
bool analytic_gradient_method(const std::string& canonical_method);

// Memory (bytes) of an analytic gradient of `method`, or 0 for RHF.
double analytic_gradient_bytes(const std::string& canonical_method, int n_occupied, int n_mo, int n_basis);

SinglePoint run_single_point(const Molecule& mol, const BasisSet& basis, const std::string& method,
                             const ElectronicOptions& options, const Matrix* guess_orbitals = nullptr,
                             const ElectronicCallbacks& callbacks = ElectronicCallbacks());

// `force_numerical`: central differences of energies even where an analytic gradient exists
// (tests and comparisons).
GradientResult run_gradient(const Molecule& mol, const BasisSet& basis, const std::string& method,
                            const ElectronicOptions& options, const Matrix* guess_orbitals = nullptr,
                            const ElectronicCallbacks& callbacks = ElectronicCallbacks(), bool force_numerical = false);

// Natural internal coordinates of the molecules of ΨLab: a diatomic has R; H2O (atom 0 = O)
// has R(O-H1), R(O-H2) and the angle H1-O-H2; H2O2 (atoms O1, O2, H3, H4) has R(O1-O2),
// R(O1-H3), R(O2-H4), A(O2-O1-H3), A(O1-O2-H4) and the dihedral D(H3-O1-O2-H4); NH3 (atoms N1,
// H2, H3, H4) has R(N1-H2), R(N1-H3), R(N1-H4), A(H2-N1-H3), A(H2-N1-H4) and the improper
// torsion D(H4-N1-H2-H3) (a Z-matrix; 180 degrees at the planar geometry) (7.0).
// values() in bohr and rad.
struct SimpleInternals {
    std::vector<std::string> names;
    std::vector<double> values;
    std::vector<std::vector<double>> B;  // Wilson B matrix rows (dq/dx), 3 n_atoms columns
};
SimpleInternals simple_internals(const Molecule& mol);

// 9.0: numerical derivatives of custom molecules with three or more atoms use Cartesian
// directions instead of the Z-matrix coordinates above (their atoms may come in any order).
bool uses_cartesian_displacements(const Molecule& mol);
// Orthonormal Cartesian directions (3 n_atoms components each) that span the complement of the
// three translations and the rigid rotations about the centroid: 3N - 6 directions, 3N - 5 for a
// linear molecule (*linear). Energies and gradients are invariant along the translations and
// rotations, so central differences along these directions give the whole gradient (2(3N - 6)
// energies) and the whole Hessian at a stationary point (2(3N - 6) gradients).
std::vector<std::vector<double>> vibrational_directions(const Molecule& mol, bool* linear = nullptr);
// The translations (3) and rotations (2 or 3) as orthonormal Cartesian directions u, with the
// derivative of the dipole moment along each (translations change the dipole by charge times
// the shift, a rotation by theta about n through the centroid c by theta n x (mu - charge c))
// and, given the gradient g, the Hessian times u, which invariance fixes exactly: H u = 0 for a
// translation, and for the rotation about n, H (n x (x_a - c))_a = (n x g_a)_a.
struct RigidMotion {
    std::vector<std::vector<double>> directions;
    std::vector<std::array<double, 3>> dipole_derivatives;
    std::vector<std::vector<double>> hessian_times;  // with a gradient: H u for every direction
};
RigidMotion rigid_motions(const Molecule& mol, const std::array<double, 3>& dipole,
                          const std::vector<double>* gradient = nullptr);
// The geometry with the given internal-coordinate values (same orientation of the first bond).
Molecule geometry_from_internals(const Molecule& reference, const std::vector<double>& values);

}  // namespace macqc
