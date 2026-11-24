// ΨLab 20.0: nuclear motion. Rovibrational energy levels of the water isotopologues H2(16)O,
// H2(17)O, H2(18)O and D2(16)O on the CVRQD potential energy surfaces.
//
// The method is that of the program cpp_DOPI_expval_array.cpp (reference/nuclear_motion):
//
//   - Jacobi coordinates: R1 = r(O-H1), R2 = the distance of H2 from the centre of mass of O-H1,
//     theta = the angle between them; the body-fixed z axis along R1;
//   - a discrete variable representation (DVR): Hermite DVR points in R1 and in R2 on
//     Re +- Rmax (the largest Hermite root at Rmax), Legendre DVR points in cos(theta) (the
//     associated Legendre functions P_l^K, l = K ... K + n - 1, for each K);
//   - the Hamiltonian
//       H = T(R1) + T(R2) + [1/(2 mu1 R1^2) + 1/(2 mu2 R2^2)] j^2 + V
//           + for J > 0: [J(J+1) - 2K^2] / (2 mu1 R1^2) and the Coriolis coupling of K and K +- 1,
//     in a parity-adapted (Wang) basis of K: K = 0 ... J or K = 1 ... J;
//   - the potential measured from its minimum, so that H is positive definite;
//   - a sparse (CSR) matrix, its lowest eigenvalues by shift-invert Lanczos (every step solves
//     H x = q by preconditioned conjugate gradients or MINRES), the Jacobi-Davidson method, the
//     Davidson method, or a dense diagonalisation for small matrices;
//   - expectation values of the bond lengths r(O-H1), r(O-H2) and of the angle H-O-H.
#pragma once

#include <functional>
#include <string>
#include <vector>

#include "macqc/cvrqd.h"

namespace macqc {

struct NuclearMotionInput {
    CvrqdIsotopologue isotopologue = CvrqdIsotopologue::H2_16O;
    int J = 0;                  // total angular momentum
    int parity_block = 0;       // J > 0: 0 = K from 0 to J, 1 = K from 1 to J ("npar" of the original)
    int radial_points = 20;     // Hermite DVR points in R1 and in R2 ("ns1")
    int angular_points = 24;    // Legendre DVR points in cos(theta) ("np")
    double r_centre = 2.0;      // bohr: centre of the radial grids ("Re")
    double r_half_width = 0.8;  // bohr: the outermost points at r_centre +- r_half_width ("Rmax")
    int levels = 10;            // number of energy levels

    // "lanczos" (shift-invert), "jacobi-davidson", "davidson" or "dense"
    std::string eigensolver = "lanczos";
    // Lanczos: "pcg" (preconditioned conjugate gradients) or "minres"
    std::string linear_solver = "pcg";
    // "ilu" (symmetric Gauss-Seidel factors of H, as in the original), "jacobi" (diagonal) or "none"
    std::string preconditioner = "ilu";

    double eigenvalue_tolerance = 1.0e-10;  // E_h: change of the eigenvalues between checks
    double residual_tolerance = 1.0e-8;     // E_h: ||H x - E x|| of the eigenvectors
    double linear_tolerance = 1.0e-12;      // residual of H x = q (Lanczos)
    int max_iterations = 1000;              // Lanczos steps or (Jacobi-)Davidson iterations
    int max_linear_iterations = 10000;      // per linear solve
    int jd_inner_iterations = 15;           // MINRES steps of the Jacobi-Davidson correction equation
    int subspace_size = 0;                  // (Jacobi-)Davidson subspace before a restart (0: automatic)

    bool expectation_values = true;
    int threads = 0;  // 0: all cores
};

// Expectation values of a distance r (bohr): <r>, <r^2>^(1/2), <r^3>^(1/3), <r^-1>^(-1),
// <r^-2>^(-1/2) and <r^-3>^(-1/3) (as printed by the original program).
struct DistanceMoments {
    double mean = 0, rms = 0, cubic = 0, inverse = 0, inverse2 = 0, inverse3 = 0;
};

struct NuclearLevel {
    double energy = 0;      // E_h, from the minimum of the potential
    double residual = 0;    // ||H x - E x|| (E_h)
    DistanceMoments r1, r2;  // r(O-H1), r(O-H2) (bohr)
    double theta = 0;       // <angle H-O-H> (radians)
    double cos_theta = 0;   // <cos HOH>
    std::vector<double> k_weights;  // J > 0: the weight of each K block (sum 1)
};

struct NuclearSolverStep {
    int iteration = 0;
    double lowest = 0;      // E_h
    double change = 0;      // largest change of the wanted eigenvalues (Lanczos) or residual (Davidson)
    int converged = 0;      // converged levels
    int linear_iterations = 0;  // Lanczos: iterations of this step's linear solve
};

struct NuclearMotionResult {
    NuclearMotionInput input;
    std::string isotopologue_name;   // e.g. "H2(16)O"
    double masses[3] = {0, 0, 0};    // u: H1 (or D), H2 (or D), O
    double mu1 = 0, mu2 = 0;         // reduced masses of R1 and R2 (m_e)
    double potential_minimum = 0;    // E_h (the zero of the energies)
    double minimum_r = 0, minimum_angle = 0;  // bohr, radians (symmetric minimum)
    std::vector<double> r_grid;      // bohr (R1 and R2 share the points)
    std::vector<int> k_values;       // the K of each block
    std::vector<std::vector<double>> cos_grid;  // cos(theta) points of each K block
    long long dimension = 0;
    long long nonzeros = 0;
    double hamiltonian_seconds = 0;
    double solver_seconds = 0;
    int iterations = 0;              // Lanczos steps or (Jacobi-)Davidson iterations
    long long matrix_vector_products = 0;
    long long linear_iterations = 0; // all iterations of the linear solves
    bool converged = false;
    std::string solver_description;
    std::vector<NuclearSolverStep> history;
    std::vector<NuclearLevel> levels;
    std::vector<std::string> notes;
    std::string report;              // the output file
};

using NuclearOutputSink = std::function<void(const std::string&)>;

// Runs the calculation; throws std::invalid_argument for an input it cannot run and Cancelled
// when the current job is cancelled. `sink` receives the output file while it is written.
NuclearMotionResult run_nuclear_motion(const NuclearMotionInput& input, const NuclearOutputSink& sink = {});

// The name of an isotopologue ("H2(16)O", "H2(17)O", "H2(18)O", "D2(16)O") and back.
std::string isotopologue_name(CvrqdIsotopologue iso);
CvrqdIsotopologue isotopologue_from_name(const std::string& name);

// Memory (MB) and dimension of a calculation before it runs (for the app's estimate).
struct NuclearMotionEstimate {
    long long dimension = 0;
    long long nonzeros = 0;
    double memory_mb = 0;
    std::string refusal;  // empty if the calculation can run
};
NuclearMotionEstimate estimate_nuclear_motion(const NuclearMotionInput& input);

}  // namespace macqc
