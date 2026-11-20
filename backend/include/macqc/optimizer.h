// Geometry optimisation in redundant internal coordinates (ΨLab 5.0).
//
// Each step: energy and Cartesian gradient (analytic where available, see gradient.h);
// internal gradient g_q = G^- B g_x; projected gradient and Hessian (P = G G^-,
// H~ = P H P + 1000 (1 - P)); rational-function (RFO) step from the lowest eigenvector of the
// augmented Hessian [[H~, g~], [g~', 0]], scaled down to the trust radius; iterative
// back-transformation to Cartesian coordinates; BFGS update of the internal Hessian with the
// actual change of the internal coordinates and gradients (skipped when the curvature
// condition dq.dg > 0 fails). The trust radius adapts to the ratio of the actual and the
// predicted energy change (Fletcher). Initial Hessian: the analytic RHF Hessian transformed
// to internal coordinates (including the gradient term) when the method is RHF and the basis
// is small, otherwise Lindh's model Hessian (Chem. Phys. Lett. 241, 423 (1995)) or
// Schlegel's (Theor. Chim. Acta 66, 333 (1984)).
//
// Convergence, as in Gaussian: maximum and RMS internal force (E_h/bohr, E_h/rad) and maximum
// and RMS internal displacement of the next step (bohr, rad), plus the energy change since
// the previous step; converged when all are met, or when the forces are 100 times below
// their thresholds.
//
// Transition states (7.0, Gaussian Opt=TS): the Hessian is computed at the first point
// (CalcFC: analytic RHF, or central differences of analytic gradients or of energies, see
// hessian.h) or at every point (CalcAll) and transformed to internal coordinates; the step is
// the partitioned RFO step of Baker (J. Comput. Chem. 7, 385 (1986)): in the eigenvectors v_i
// (eigenvalues b_i, gradient components F_i) of the projected Hessian, the followed mode k (the
// lowest at the first step, later the one with the largest overlap with the previous followed
// mode) gets -F_k / (b_k - lambda_p) with lambda_p = b_k/2 + sqrt(b_k^2/4 + F_k^2) (a maximisation)
// and every other mode -F_i / (b_i - lambda_n) with lambda_n the lowest eigenvalue of their
// augmented Hessian (a minimisation); the step is scaled down to the trust radius, which stays
// fixed (0.3 bohr/rad, as in Gaussian for transition states). Between computed Hessians the
// Hessian is updated by Bofill's formula (Int. J. Quantum Chem. 52, 1257 (1994)): the
// Murtagh-Sargent (SR1) and Powell (PSB) updates mixed with phi = (xi.s)^2 / ((xi.xi)(s.s)),
// xi = y - H s. 15.0: also RecalcFC=N (Gaussian: the Hessian computed at the first step and every
// N-th step thereafter, updated in between). The eigenvalue test (on by default, Gaussian's EigenTest) stops the search when
// the Hessian does not have exactly one negative eigenvalue.
//
// 15.0, minimisations: computed Hessians as for transition states (Gaussian Opt=CalcFC, CalcAll,
// RecalcFC=N), made positive definite for the RFO step (negative eigenvalues replaced by their
// magnitude), BFGS updates in between. Optionally (algorithm = "gediis") the hybrid GEDIIS/RFO step
// of Li and Frisch (J. Chem. Theory Comput. 2, 835 (2006)), which Gaussian uses by default for
// minimisations: once the RMS force is below 1e-2, the energy-represented DIIS interpolation of the
// last (up to five) points, the coefficients c_i >= 0, sum c_i = 1, minimising
//   E(c) = sum_i c_i E_i - 1/2 sum_ij c_i c_j (g_i - g_j).(q_i - q_j),
// (used while the RMS force is below 1e-2 and the RMS RFO step above 2.5e-3, the thresholds of
// the hybrid method), followed by the RFO step from the interpolated point with the interpolated gradient
// (q* = sum c_i q_i, g* = sum c_i g_i); a step that is longer than the trust radius, goes uphill
// or does not move away from the current point falls back to the plain RFO step.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "macqc/basis.h"
#include "macqc/gradient.h"
#include "macqc/internal_coordinates.h"
#include "macqc/molecule.h"

namespace macqc {

struct ConvergenceCriteria {
    double max_force = 4.5e-4, rms_force = 3.0e-4, max_displacement = 1.8e-3, rms_displacement = 1.2e-3;
    double energy = 1.0e-6;
};
// "normal" (Gaussian's default) or "tight" (Gaussian's tight: 1.5e-5, 1.0e-5, 6.0e-5, 4.0e-5; 1e-8 E_h).
ConvergenceCriteria convergence_criteria(bool tight);

// A frozen internal coordinate (Gaussian: Opt=ModRedundant ... F; the scanned coordinate of a
// relaxed scan, 7.0). The optimiser keeps it at `target` with the projector
// P~ = P - P C (C' P C)^-1 C' P (Peng, Ayala, Schlegel & Frisch 1996): the RFO step and the
// convergence test use the projected gradient P~ g, and every step adds the change
// P C (C' P C)^-1 (target - q_c) that drives the coordinate to its target. A bend frozen at
// exactly 180 degrees is imposed through its two linear-bend components (both 0).
struct CoordinateConstraint {
    InternalCoordinate coordinate;  // Stretch, Bend or Dihedral
    double target = 0.0;            // bohr or rad
};

struct OptimizationOptions {
    bool tight = false;
    int max_steps = 50;
    double trust_radius = 0.3;   // initial (bohr, rad)
    double max_trust = 0.5, min_trust = 1.0e-3;
    // auto (10.0: the analytic Hessian of RHF up to analytic_hessian_max_functions, of Full CI,
    // CCSD and CCSD(T) where it fits into memory; otherwise Lindh) | analytic | lindh | schlegel
    std::string initial_hessian = "auto";
    int analytic_hessian_max_functions = 100;  // "auto": analytic RHF Hessian up to this size
    // Called before the energy and gradient of each step (for the output file).
    std::function<void(int, const Molecule&)> on_step_start;
    // 7.0: frozen coordinates; SCF guess for the first step; a Hessian in internal coordinates
    // to start from (e.g. that of the previous point of a scan), used when `hessian_coordinates`
    // equals the coordinate set built at the first step.
    std::vector<CoordinateConstraint> constraints;
    const Matrix* guess_orbitals = nullptr;
    const Matrix* initial_hessian_matrix = nullptr;
    std::vector<InternalCoordinate> hessian_coordinates;
    std::string initial_hessian_label = "Hessian of the previous scan point";
    // 7.0 transition-state search (see above)
    bool transition_state = false;
    std::string ts_hessian = "calcfc";   // calcfc (computed at the first point) | calcall (every point) |
                                         // recalcfc (15.0: the first point and every recalc_every-th step)
    // 15.0, minimisations: "" (the initial Hessian of initial_hessian, then BFGS) | calcfc | calcall |
    // recalcfc (computed Hessians as for transition states, see above)
    std::string computed_hessian;
    int recalc_every = 0;                // RecalcFC=N
    std::string algorithm = "rfo";       // minimisations: rfo | gediis (hybrid GEDIIS/RFO)
    bool eigen_test = true;              // stop unless the Hessian has exactly one negative eigenvalue
    bool numerical_hessian = false;      // 10.0: computed Hessians by differences even where analytic
};

struct ConvergenceItem {
    std::string name;
    double value = 0.0, threshold = 0.0;
    bool met = false;
    bool applicable = true;  // the energy change has no value at the first step
};

struct OptimizationStep {
    int step = 0;
    Molecule geometry;                    // where energy and gradient were computed
    double energy = 0.0, delta_e = 0.0;   // delta_e: change since the previous step
    std::string gradient_description;     // analytic / numerical, how
    std::vector<std::string> names, definitions;  // of the internal coordinates
    std::vector<double> q, gq, dq;        // values, (projected) gradient, step taken from here
    std::vector<double> cartesian_gradient;
    std::string hessian_update;           // how the Hessian of this step was obtained
    double trust_radius = 0.0, step_norm = 0.0, rfo_lambda = 0.0, predicted_change = 0.0;
    bool step_scaled = false;
    double ratio = 0.0;                   // actual / predicted change of the previous step
    bool has_ratio = false;
    int backtransform_iterations = 0;
    double backtransform_residual = 0.0;
    std::vector<ConvergenceItem> convergence;
    bool converged = false;
    std::vector<std::string> warnings;
    double seconds = 0.0;
    // 7.0: frozen coordinates (index into names, deviation from the target at this geometry)
    std::vector<int> constrained;
    double constraint_deviation = 0.0;
    // 7.0 transition-state step: Hessian eigenvalues (non-redundant space, ascending), the number
    // of negative ones, the followed mode (its eigenvalue, gradient component, overlap with the
    // previously followed mode, components in the internal coordinates), lambda_p and lambda_n
    bool ts_step = false;
    std::vector<double> hessian_eigenvalues;
    int negative_eigenvalues = 0;
    double followed_eigenvalue = 0.0, followed_gradient = 0.0, followed_overlap = 1.0;
    std::vector<double> followed_mode;
    double lambda_p = 0.0, lambda_n = 0.0;
    bool hessian_computed = false;       // CalcFC first step, CalcAll or RecalcFC
    // 15.0: the kind of step (RFO, P-RFO, GEDIIS + RFO) and the GEDIIS coefficients (of the
    // points of the steps gediis_steps, the last being this one)
    std::string step_type = "RFO";
    std::vector<double> gediis_coefficients;
    std::vector<int> gediis_steps;
    double gediis_energy = 0.0;          // the interpolated energy E(c)
    std::string gediis_note;             // why GEDIIS was not used at this step (when tried)
    std::string hessian_description;     // how the computed Hessian was obtained
};

struct GeometricParameter {
    std::string name;       // e.g. "R(O1–H2)"
    double value = 0.0;     // bohr or degrees
    bool angle = false;
};

struct OptimizationResult {
    bool converged = false;
    int steps = 0;
    std::string message;
    std::string initial_hessian;
    ConvergenceCriteria criteria;
    Molecule initial, final_geometry;
    double initial_energy = 0.0, final_energy = 0.0;
    std::vector<OptimizationStep> history;
    GradientResult final_point;           // energy, gradient, orbitals at the final geometry
    std::vector<GeometricParameter> initial_parameters, final_parameters;
    bool symmetric_start = false;         // H2O: started C2v
    bool symmetric_end = false;           // H2O: C2v kept
    double symmetry_deviation = 0.0;      // |R1 - R2| at the end (bohr)
    bool removed_functions_changed = false;
    double seconds = 0.0;
    // 7.0: the internal coordinates and the (updated) Hessian of the last step, for the next
    // point of a scan
    std::vector<InternalCoordinate> final_coordinates;
    Matrix final_hessian;
    double final_constraint_deviation = 0.0;
    // 7.0 transition states
    bool transition_state = false;
    std::string ts_hessian;
    std::string computed_hessian;        // 15.0: calcfc | calcall | recalcfc (or empty)
    int recalc_every = 0;
    std::string algorithm = "rfo";
    int gediis_steps = 0;                // steps taken with GEDIIS
    int hessian_evaluations = 0;         // computed Hessians
    int final_negative_eigenvalues = -1; // of the last Hessian (computed or updated)
    bool eigen_test_failed = false;
};

using OptimizationCallback = std::function<void(const OptimizationStep&)>;

// `basis` may be built at any geometry of the molecule (shells are moved to each geometry).
OptimizationResult optimize_geometry(const Molecule& start, const BasisSet& basis, const std::string& method,
                                     const ElectronicOptions& electronic, const OptimizationOptions& options,
                                     const ElectronicCallbacks& callbacks = ElectronicCallbacks(),
                                     const OptimizationCallback& on_step = OptimizationCallback());

// Bond lengths (bohr) and angles (degrees) of the internal coordinates of a geometry.
std::vector<GeometricParameter> geometric_parameters(const Molecule& mol);

// The geometry with the constrained coordinates set to their targets, all other internal
// coordinates changed as little as possible (minimal change in the non-redundant space; for a
// non-redundant set the others stay fixed, as in a Z-matrix), in substeps of at most 0.2 bohr /
// 20 degrees. *residual: the largest remaining |q - target|.
Molecule apply_constraints(const Molecule& mol, const std::vector<CoordinateConstraint>& constraints,
                           double* residual = nullptr);

// Model Hessians in internal coordinates (diagonal).
std::vector<double> lindh_model_hessian(const Molecule& mol, const InternalCoordinates& ic);
std::vector<double> schlegel_model_hessian(const Molecule& mol, const InternalCoordinates& ic);

}  // namespace macqc
