// ΨLab 20.0/21.0: nuclear motion. Rovibrational energy levels of the water isotopologues H2(16)O,
// H2(17)O, H2(18)O and D2(16)O on the CVRQD potential energy surfaces.
//
// The method is that of the program cpp_DOPI_expval_array.cpp (reference/nuclear_motion):
//
//   - Jacobi coordinates: R1 = r(O-H1), R2 = the distance of H2 from the centre of mass of O-H1,
//     theta = the angle between them; or (21.0) Radau coordinates: R1, R2 = the distances of H1
//     and H2 from the canonical point K = alpha O + (1 - alpha) (m1 H1 + m2 H2)/(m1 + m2),
//     alpha = (m(O)/M)^1/2, theta = the angle H1-K-H2 (the original's nkoord = 2, which it had not
//     written out); both orthogonal, the body-fixed z axis along R1;
//   - a discrete variable representation (DVR): Hermite DVR points in R1 and in R2 on
//     Re +- Rmax (the largest Hermite root at Rmax), Legendre DVR points in cos(theta) (the
//     associated Legendre functions P_l^K, l = K ... K + n - 1, for each K);
//   - the Hamiltonian
//       H = T(R1) + T(R2) + [1/(2 mu1 R1^2) + 1/(2 mu2 R2^2)] j^2 + V
//           + for J > 0: [J(J+1) - 2K^2] / (2 mu1 R1^2) and the Coriolis coupling of K and K +- 1,
//     in a parity-adapted (Wang) basis of K: K = 0 ... J or K = 1 ... J;
//   - the potential measured from its minimum, so that H is positive definite;
//   - a sparse (CSR) matrix, its lowest eigenvalues by shift-invert Lanczos, the Jacobi-Davidson
//     method, the Davidson method, (21.0) inverse (subspace) iteration, Rayleigh-quotient iteration,
//     or a dense diagonalisation for small matrices; the linear systems of every iterative method
//     by preconditioned conjugate gradients, MINRES, (21.0) Chebyshev iteration or BiCG, with the
//     ILU (symmetric Gauss-Seidel), Jacobi or no preconditioner;
//   - (21.0) optionally the eigenvalues only (the original's sv = 0); the start vector of the
//     thesis's Algorithm 5 (the original's psikiind): a product of Gaussians in the bond lengths
//     and the bond angle;
//   - expectation values of the bond lengths r(O-H1), r(O-H2) and of the angle H-O-H; (21.0)
//     cuts and densities of the eigenvectors on the grid (WavefunctionPlot).
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
    double r_centre = 2.0;      // bohr: centre of the radial grids ("Re"; bohr in the original too)
    double r_half_width = 0.8;  // bohr: the outermost points at r_centre +- r_half_width ("Rmax")
    int levels = 10;            // number of energy levels
    std::string coordinates = "jacobi";  // 21.0: "jacobi" or "radau" (the original's nkoord 1, 2)
    bool eigenvectors = true;   // 21.0: false: the eigenvalues only (the original's sv = 0)
    bool gaussian_start = false;  // 21.0: the start vector of Algorithm 5 (psikiind)

    // "lanczos" (shift-invert), "jacobi-davidson", "davidson", "inverse-power" (21.0),
    // "rayleigh-quotient" (21.0) or "dense"
    std::string eigensolver = "lanczos";
    // the linear systems of every iterative eigensolver: "pcg" (preconditioned conjugate
    // gradients), "minres", "chebyshev" (21.0) or "bicg" (21.0)
    std::string linear_solver = "pcg";
    // "ilu" (symmetric Gauss-Seidel factors of H, as in the original), "jacobi" (diagonal) or "none";
    // for a shifted matrix H - theta (Davidson, Jacobi-Davidson, Rayleigh quotient) the factors
    // of H - theta
    std::string preconditioner = "ilu";

    double eigenvalue_tolerance = 1.0e-10;  // E_h: change of the eigenvalues between checks
    double residual_tolerance = 1.0e-8;     // E_h: ||H x - E x|| of the eigenvectors
    double linear_tolerance = 1.0e-12;      // residual of H x = q (Lanczos)
    int max_iterations = 1000;              // Lanczos steps or (Jacobi-)Davidson iterations
    int max_linear_iterations = 10000;      // per linear solve
    int jd_inner_iterations = 15;           // linear-solver steps of the Jacobi-Davidson correction equation
    int davidson_inner_iterations = 0;      // 21.0: Davidson: linear-solver steps of (H - theta) t = -r
                                            // (0: the preconditioner alone, the classical method)
    int subspace_size = 0;                  // (Jacobi-)Davidson subspace before a restart (0: automatic)

    bool expectation_values = true;
    int threads = 0;  // 0: all cores
};

// Expectation values of a distance r (bohr): <r>, <r^2>^(1/2), <r^3>^(1/3), <r^-1>^(-1),
// <r^-2>^(-1/2) and <r^-3>^(-1/3) (as printed by the original program).
struct DistanceMoments {
    double mean = 0, rms = 0, cubic = 0, inverse = 0, inverse2 = 0, inverse3 = 0;
};

// 21.0: the eigenvector of a level on the DVR grid. psi = c / (s sqrt(w)) at the grid points (the
// DVR coefficient c over the sign s and the square root of the quadrature weight w of the DVR
// function), scaled so that its largest value is +1, at the point (block, i, j, l) of the largest
// |psi|: cuts through that point in the planes (R1, R2), (R1, theta) and (R2, theta); and the
// probability densities |psi|^2 integrated over the third coordinate (rho_r1r2 also over the K
// blocks; the planes with theta in the K block `block`): bohr^-2 and bohr^-1 rad^-1. Row-major:
// r1r2: ns x ns (R1 major), r1t: ns x np (R1 major), r2t: ns x np (R2 major).
struct WavefunctionPlot {
    int block = 0;
    int i = 0, j = 0, l = 0;
    std::vector<float> psi_r1r2, psi_r1t, psi_r2t;
    std::vector<float> rho_r1r2, rho_r1t, rho_r2t;
};

struct NuclearLevel {
    double energy = 0;      // E_h, from the minimum of the potential
    double residual = 0;    // ||H x - E x|| (E_h); NaN: not known (Lanczos, eigenvalues only)
    DistanceMoments r1, r2;  // r(O-H1), r(O-H2) (bohr)
    double theta = 0;       // <angle H-O-H> (radians)
    double cos_theta = 0;   // <cos HOH>
    std::vector<double> k_weights;  // J > 0: the weight of each K block (sum 1)
    WavefunctionPlot plot;  // 21.0 (with the eigenvectors, at most the first 60 levels)
};

struct NuclearSolverStep {
    int iteration = 0;
    double lowest = 0;      // E_h: the lowest level, or (Rayleigh quotient) the level being refined
    double change = 0;      // largest change of the wanted eigenvalues (Lanczos) or residual (the others)
    int converged = 0;      // converged levels
    int linear_iterations = 0;  // iterations of this step's linear solves
    int level = 0;          // 21.0: the level being refined (1, 2, ...; 0: all at once)
    double linear_residual = 0;  // 21.0: the residual of this step's (last) linear solve
    double seconds = 0;     // 21.0: since the start of the eigensolver
};

// 21.0: the residual norms of one linear solve, iteration by iteration (Convergence tab).
struct LinearSolveTrace {
    int step = 0;           // the eigensolver iteration it belongs to (0: none recorded)
    std::vector<double> residuals;
};

// 20.0: a picture of the sparsity of the Hamiltonian: rows and columns first ... first + span - 1
// (from 0) in bins x bins pixels (a pixel is one element when span <= bins). density: the fraction
// of the elements of each pixel that are non-zero (0 ... 255, 0 only when none is); magnitude: the
// largest |H_ij| of the pixel on a logarithmic scale, 1 ... 255 for 10^log_min ... 10^log_max E_h
// (0: no non-zero element). Row-major, row 0 at the top.
struct MatrixPicture {
    int bins = 0;
    long long first = 0, span = 0;
    std::vector<unsigned char> density, magnitude;
    double log_min = 0, log_max = 0;
};

struct NuclearMotionResult {
    NuclearMotionInput input;
    std::string isotopologue_name;   // e.g. "H2(16)O"
    double masses[3] = {0, 0, 0};    // u: H1 (or D), H2 (or D), O
    double mu1 = 0, mu2 = 0;         // reduced masses of R1 and R2 (m_e)
    double g1 = 0, g2 = 0;           // 21.0: the Sutcliffe-Tennyson parameters of the coordinates
    bool has_vectors = false;        // 21.0: eigenvectors computed (and kept)
    double start_rayleigh = 0;       // 21.0: Gaussian start vector: its Rayleigh quotient (E_h)
    double start_overlap = -1;       // 21.0: |<q1|x1>| with the lowest eigenvector (-1: unknown)
    double start_mu[3] = {0, 0, 0};     // 21.0: the centres of the Gaussians: r1, r2 (bohr), angle (rad)
    double start_sigma[3] = {0, 0, 0};  // 21.0: their widths (harmonic estimates)
    double chebyshev_min = 0, chebyshev_max = 0;  // 21.0: Chebyshev: the interval of M^-1 A
    LinearSolveTrace first_solve, last_solve;     // 21.0
    double potential_minimum = 0;    // E_h (the zero of the energies)
    double minimum_r = 0, minimum_angle = 0;  // bohr, radians (symmetric minimum)
    std::vector<double> r_grid;      // bohr (R1 and R2 share the points)
    std::vector<double> r_weights;   // 21.0: bohr, the quadrature weights of the radial points
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
    MatrixPicture matrix_whole;      // the whole Hamiltonian (at most 1024 x 1024 pixels)
    MatrixPicture matrix_detail;     // its leading rows and columns, one pixel per element
    std::vector<std::string> notes;
    std::string report;              // the output file
    // 21.0, for the export (not in the JSON result): the normalised eigenvectors (DVR coefficients,
    // empty with the eigenvalues only), the Gaussian start vector q1 (empty without it), and the
    // quadrature weights and signs of the angular DVR functions of each K block
    std::vector<std::vector<double>> vectors;
    std::vector<double> start_vector;
    std::vector<std::vector<double>> cos_weights, cos_signs;
};

// 21.0: the table of the energy levels as in the output file: n, E (E_h), E (cm-1), E - E1 (cm-1),
// ||Hx-Ex|| (E_h), with its two header lines.
std::string nuclear_level_table(const NuclearMotionResult& result);

// 21.0: writes the eigenvalues (the level table), the eigenvectors and/or the Gaussian start vector
// of a result as text to `path` (any combination, at least one). The vectors: one row per basis
// function (its index, K, the grid indices, R1, R2, theta, the sign and the weight of its DVR
// function), then q1 and c_1 ... c_m. Throws std::invalid_argument when a part is not available,
// std::runtime_error when the file cannot be written.
void write_nuclear_export(const NuclearMotionResult& result, const std::string& path, bool eigenvalues, bool eigenvectors,
                          bool start_vector);

using NuclearOutputSink = std::function<void(const std::string&)>;
// 22.0: text written after the results and before the termination banner (the archive of the
// results, which File > Open Output File... of the app reads back).
using NuclearAppendixWriter = std::function<std::string(const NuclearMotionResult&)>;

// Runs the calculation; throws std::invalid_argument for an input it cannot run and Cancelled
// when the current job is cancelled. `sink` receives the output file while it is written; 22.0:
// `appendix` adds its text before the termination banner.
NuclearMotionResult run_nuclear_motion(const NuclearMotionInput& input, const NuclearOutputSink& sink = {},
                                       const NuclearAppendixWriter& appendix = {});

// 22.0: the termination banner of the output file (as that of the electronic-structure jobs):
// status "normal", "unconverged", "cancelled" or "error" (with `message`).
std::string nuclear_termination_text(const std::string& status, double seconds, const std::string& message = "");

// 22.0: nuclear-motion input files of the terminal version (psilab INPUT > OUTPUT): one
// "key = value" per line, '#' starts a comment (keys and values as in nuclear_input_text).
// is_nuclear_input_text: the text asks for a nuclear-motion calculation (calculation = nuclear,
// or molecule/isotopologue = H2(16)O, H2(17)O, H2(18)O, D2(16)O).
bool is_nuclear_input_text(const std::string& text);
// Throws std::runtime_error("line N: ...") on an unknown key, a key given twice or a bad value.
NuclearMotionInput parse_nuclear_input_text(const std::string& text);
// The input as such a file (the "Input" section of the output file; it runs the same calculation).
std::string nuclear_input_text(const NuclearMotionInput& input);
// A commented example input file (psilab --example nuclear).
std::string nuclear_example_input();

// 22.0: the archive section of the output file: the complete results as JSON (as the app receives
// them, without the output text) between "@@PSILAB-ARCHIVE-BEGIN 1" and "@@PSILAB-ARCHIVE-END",
// with "archive": {"kind": "nuclear-motion", ...} (api.cpp).
std::string nuclear_archive_text(const NuclearMotionResult& result);

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
