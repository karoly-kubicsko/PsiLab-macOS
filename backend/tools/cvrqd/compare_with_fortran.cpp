// Compares the C++ translation of a CVRQD surface (backend/src/cvrqd.cpp) with the gfortran build
// of the original Fortran routine at 300000 random bond geometries. Build one program per
// isotopologue (the H2(16)O and H2(18)O files define the same subroutine names), e.g.
//   gfortran -c reference/nuclear_motion/pot_h216o_cvrqd.f -o f.o
//   g++ -std=c++17 -O2 -ffp-contract=off -Ibackend/include -DFN=potv_h216o_cvrqd_ -DISO=H2_16O \
//       -DNAME='"H2_16O"' backend/tools/cvrqd/compare_with_fortran.cpp backend/src/cvrqd.cpp f.o -lgfortran
// (ΨLab 20.0: the largest difference was 0 for all four surfaces.)
#include <cmath>
#include <cstdio>
#include <random>
#include "macqc/cvrqd.h"
extern "C" void FN(double*, double*, double*, double*);
int main() {
    std::mt19937_64 rng(7);
    std::uniform_real_distribution<double> r(1.0, 4.5), c(-0.999, 0.999);
    double maxabs = 0, maxrel = 0;
    for (int n = 0; n < 300000; ++n) {
        double r1 = r(rng), r2 = r(rng), x = c(rng);
        if (n == 0) { r1 = 1.81; r2 = 1.81; x = std::cos(104.5 * M_PI / 180); }
        double a1 = r1, a2 = r2, ax = x, vf;
        FN(&vf, &a1, &a2, &ax);
        double vc = macqc::cvrqd_potential(macqc::CvrqdIsotopologue::ISO, r1, r2, x);
        double d = std::fabs(vf - vc);
        if (d > maxabs) maxabs = d;
        double rel = d / std::fmax(std::fabs(vf), 1e-3);
        if (rel > maxrel) maxrel = rel;
        if (n == 0) printf("%-7s V(1.81,1.81,104.5 deg): Fortran %.16f  C++ %.16f\n", NAME, vf, vc);
    }
    printf("%-7s 300000 random geometries (R 1-4.5 bohr): max |dV| = %.2e Eh, max rel = %.2e\n", NAME, maxabs, maxrel);
}
