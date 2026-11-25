// CVRQD potential energy surfaces of water isotopologues (ΨLab 20.0, nuclear-motion module).
#pragma once

namespace macqc {

enum class CvrqdIsotopologue { H2_16O = 1, H2_17O = 2, D2_16O = 3, H2_18O = 4 };

// The CVRQD potential (hartree) at bond coordinates: r1 = r(O-H1), r2 = r(O-H2) in bohr and
// the cosine of the H-O-H angle. Identical in arithmetic to the gfortran build of the
// Fortran routines potv_<isotopologue>_cvrqd (see cvrqd.cpp).
double cvrqd_potential(CvrqdIsotopologue iso, double r1, double r2, double cos_theta);

}  // namespace macqc
