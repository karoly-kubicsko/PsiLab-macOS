// ΨLab 14.0: integrals of semi-local effective core potentials (ECPs) over Cartesian Gaussian
// shells.
//
// An ECP on the atom C is U(r) = U_L(r) + sum_{l<L} sum_m |lm> [U_l(r) - U_L(r)] <lm| with radial
// parts sum_k d_k exp(-z_k r^2) (r = |r - C|; the Gaussian terms n = 2 of Gaussian's format, those
// of the Stuttgart-Cologne potentials of the def2 sets).
//  - Local part: <a| exp(-z |r - C|^2) |b> is a three-centre overlap (Obara-Saika).
//  - Semi-local parts: the projection of a shell function onto the real spherical harmonic Y_lm
//    on the sphere of radius r about C is, with a = A - C and k = 2 alpha |a| r,
//      int dOmega phi(r) Y_lm = 4 pi exp(-alpha (r - |a|)^2) sum_{pqs,lambda} B_pqs r^(p+q+s)
//                               i~_lambda(k) sum_mu Y_lambda,mu(a^) int dOmega Y_lambda,mu Y_lm x^p y^q z^s
//    (the binomial expansion of (r - a)^c and the expansion of exp(2 alpha r a.r^) in modified
//    spherical Bessel functions, i~_lambda(k) = exp(-k) i_lambda(k); McMurchie and Davidson,
//    J. Comput. Phys. 44, 289 (1981)); the angular integrals are exact, the radial integral of the
//    two projections times exp(-z r^2) r^2 is done for each pair of primitives and term by
//    Gauss-Legendre quadrature over +-9 widths about the peak of the Gaussian product (cf.
//    Flores-Moreno et al., J. Comput. Chem. 27, 1009 (2006)).
// The ECP is added to the nuclear attraction of its atom (nuclear_cartesian), so every derivative
// of the integrals with respect to A, B and (by translational invariance) C follows.
#pragma once

#include <vector>

#include "macqc/basis.h"
#include "macqc/molecule.h"

namespace macqc {

// Highest angular momentum of a shell in ECP integrals (second derivatives of h functions: 7).
constexpr int kEcpMaxShellL = 9;

// Adds <a| U_C |b> of the ECP of `center` (nothing if it has none) for the Cartesian components of
// the contracted shells A and B (coefficients as given) to cart[ncart(A.l)][ncart(B.l)].
// `points`: Gauss-Legendre points of the radial quadrature (0: the default, 64; tests).
void add_ecp_cartesian(const Shell& A, const Shell& B, const Atom& center, std::vector<double>& cart, int points = 0);

// exp(-x) i_n(x), n = 0 .. nmax, of the modified spherical Bessel functions of the first kind
// (x >= 0; nmax <= 2 kEcpMaxShellL). Exposed for tests.
void scaled_bessel_i(int nmax, double x, double* out);

}  // namespace macqc
