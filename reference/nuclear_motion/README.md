# Nuclear-motion reference sources (ΨLab 20.0)

- `cpp_DOPI_expval_array.cpp`: the original rovibrational program of the nuclear-motion module
  (DVR Hamiltonian of a triatomic in Jacobi coordinates, Lanczos, Jacobi–Davidson, conjugate
  gradients, MINRES, expectation values). Line endings converted from CR to LF; otherwise as
  received. ΨLab's implementation is `backend/src/nuclear_motion.cpp`.
- `pot_h216o_cvrqd.f`, `pot_h217o_cvrqd.f`, `pot_h218o_cvrqd.f`, `pot_d216o_cvrqd.f`: the CVRQD
  potential energy surfaces of the water isotopologues (Fortran 77; line endings converted to
  LF). `backend/src/cvrqd.cpp` is generated from them by `backend/tools/cvrqd/gen.py`; the four
  files share one surface and differ only in the adiabatic (BODC) term, which the generated
  file selects by isotopologue. `backend/tools/cvrqd/compare_with_fortran.cpp` compares the
  translation with the gfortran build (identical at 300000 random geometries for all four).
