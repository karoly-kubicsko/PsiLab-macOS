#!/usr/bin/env python3
"""Reference values for the molecules of ΨLab 7.0, H2O2 and NH3, and for HeH+ with
correlation-consistent basis sets (PySCF + geomeTRIC).

At ΨLab's default geometries (H2O2: r(OO) 1.4556 Å, r(OH) 0.9670 Å, O-O-H 102.32°, H-O-O-H 113.70°;
NH3: r(NH) 1.012 Å, H-N-H 106.7°, built as in macqc/molecule.cpp):
  - RHF energies and analytic gradients in aug-cc-pVDZ, aug-cc-pVTZ and aug-cc-pCVDZ (aug-cc-pVDZ
    for H in the core-valence set, as in ΨLab);
  - CCSD(T)/aug-cc-pVDZ energies, all electrons correlated, and the CCSD(T) gradient of NH3 (with
    the CCSD(T) Lambda equations, see pyscf_optimization_reference.py).
Tight optimisations (PySCF analytic gradients, geomeTRIC): RHF/aug-cc-pVDZ of both molecules and
CCSD(T)/aug-cc-pVDZ of NH3; at the RHF minima the harmonic wavenumbers from PySCF's analytic
Hessian with the AME 2020 masses.

HeH+ (charge +1, He at the origin, H at z = R = 1.4632 bohr, as in macqc/molecule.cpp; stage 7.4):
RHF energies and gradients in aug-cc-pVDZ, aug-cc-pVTZ and aug-cc-pCVDZ (aug-cc-pVDZ for He and H,
which have no core-valence sets), the Full CI energy in aug-cc-pVTZ (PySCF's FCI solver), the
CCSD/aug-cc-pVTZ minimum (exact for two electrons; analytic CCSD gradients, geomeTRIC) and its
harmonic wavenumber from a five-point second derivative of CCSD energies (AME 2020 masses).

Usage (needs: pip install pyscf geometric):
  python3 reference/pyscf_molecules_reference.py           # -> backend/tests/data/molecules_reference.json
  python3 reference/pyscf_molecules_reference.py --check   # recompute and compare
"""
import argparse
import json
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
import pyscf_optimization_reference as orf  # noqa: E402
import pyscf_scan_reference as sr  # noqa: E402
import numpy  # noqa: E402
import pyscf  # noqa: E402
from pyscf import cc, gto, scf  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402
from pyscf.grad import ccsd_t as ccsd_t_grad  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "molecules_reference.json"
B = rr.BOHR_IN_ANGSTROM
BASES = ["aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pCVDZ"]


def atoms(name):
    if name == "H2O2":
        g = sr.H2O2
        return sr.h2o2_atoms(g["ROO"], g["R"], g["angle"], g["dihedral"])
    return sr.nh3_atoms(sr.NH3["R"], sr.NH3["angle"])


def build(at, basis):
    bas = {}
    for s, _ in at:
        b, _used = rr.element_basis(basis, s)
        bas[s] = b
    return gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, verbose=0)


def rhf(mol):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.direct_scf_tol = 1e-14
    mf.kernel()
    return mf


def ccsd_t(mol, gradient):
    mf = rhf(mol)
    mc = cc.CCSD(mf)
    mc.conv_tol = 1e-12
    mc.conv_tol_normt = 1e-9
    mc.max_cycle = 300
    mc.kernel()
    e = mc.e_tot + mc.ccsd_t()
    g = ccsd_t_grad.Gradients(mc).kernel() if gradient else None
    return e, g


HEH_R = 1.4632
MASS = {"H": 1.007825031898, "He": 4.00260325413}  # AME 2020
HARTREE_J, AMU_KG, BOHR_M, C_CM = 4.3597447222071e-18, 1.66053906660e-27, 5.29177210903e-11, 2.99792458e10


def heh(r, basis):
    used = "aug-cc-pV" + basis[len("aug-cc-pCV"):] if basis.startswith("aug-cc-pCV") else basis
    at = [("He", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, r))]
    bas = {s: rr.element_basis(used, s)[0] for s, _ in at}
    return gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, charge=1, verbose=0)


def heh_ccsd(r, basis):
    mf = rhf(heh(r, basis))
    mc = cc.CCSD(mf)
    mc.conv_tol = 1e-12
    mc.conv_tol_normt = 1e-9
    mc.kernel()
    return mc


def heh_reference():
    from pyscf import fci
    out = {"points": [], "optimization": None}
    for basis in BASES:
        mol = heh(HEH_R, basis)
        mf = rhf(mol)
        g = mf.nuc_grad_method().kernel()
        out["points"].append(dict(basis=basis, method="RHF", nbf=mol.nao, energy=float(mf.e_tot),
                                  gradient=[[float(c) for c in row] for row in g]))
    mol = heh(HEH_R, "aug-cc-pVTZ")
    mf = rhf(mol)
    e_fci = fci.FCI(mf).kernel()[0]
    out["points"].append(dict(basis="aug-cc-pVTZ", method="FCI", nbf=mol.nao, energy=float(e_fci), gradient=None))
    mc = heh_ccsd(HEH_R, "aug-cc-pVTZ")
    mol_eq = optimize(mc, maxsteps=100, **orf.CONV)
    x = mol_eq.atom_coords(unit="Bohr")
    re = float(abs(x[1][2] - x[0][2]))
    h = 0.005
    e = [heh_ccsd(re + k * h, "aug-cc-pVTZ").e_tot for k in (-2, -1, 0, 1, 2)]
    k2 = (-e[0] + 16 * e[1] - 30 * e[2] + 16 * e[3] - e[4]) / (12 * h * h)  # E_h / bohr^2
    mu = MASS["He"] * MASS["H"] / (MASS["He"] + MASS["H"]) * AMU_KG
    w = (k2 * HARTREE_J / BOHR_M**2 / mu) ** 0.5 / (2 * 3.141592653589793 * C_CM)
    out["optimization"] = dict(basis="aug-cc-pVTZ", method="CCSD", R=re * B, energy=float(e[2]), wavenumber=float(w))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    doc = {"description": "H2O2 and NH3 at ΨLab's default geometries: RHF energies and gradients (E_h, E_h/bohr), "
                          "CCSD(T)/aug-cc-pVDZ energies (NH3 also the gradient); tight optimisations (bond lengths in "
                          "angstrom, angles in degrees) and RHF wavenumbers (cm^-1) at the RHF minima; generated by "
                          "reference/pyscf_molecules_reference.py",
           "program": f"PySCF {pyscf.__version__}, geomeTRIC",
           "settings": "all electrons correlated; geomeTRIC thresholds " + json.dumps(orf.CONV),
           "bohr_in_angstrom": B, "points": [], "optimizations": []}
    worst = {"energy": 0.0, "gradient": 0.0, "geometry": 0.0, "frequency": 0.0}
    for name in ("H2O2", "NH3"):
        for basis in BASES:
            t0 = time.time()
            mol = build(atoms(name), basis)
            mf = rhf(mol)
            g = mf.nuc_grad_method().kernel()
            p = dict(molecule=name, basis=basis, method="RHF", nbf=mol.nao, energy=float(mf.e_tot),
                     gradient=[[float(c) for c in row] for row in g])
            doc["points"].append(p)
            print(f"{name:5s} {basis:13s} RHF      E = {p['energy']:.10f}  nbf {mol.nao}  max|g| "
                  f"{abs(g).max():.6f}  ({time.time() - t0:.1f} s)", flush=True)
        t0 = time.time()
        mol = build(atoms(name), "aug-cc-pVDZ")
        e, g = ccsd_t(mol, name == "NH3")
        p = dict(molecule=name, basis="aug-cc-pVDZ", method="CCSD(T)", nbf=mol.nao, energy=float(e),
                 gradient=None if g is None else [[float(c) for c in row] for row in g])
        doc["points"].append(p)
        print(f"{name:5s} aug-cc-pVDZ   CCSD(T)  E = {e:.10f}  ({time.time() - t0:.1f} s)", flush=True)
    t0 = time.time()
    doc["HeH+"] = heh_reference()
    hr = doc["HeH+"]
    print("HeH+  " + "  ".join(f"{p['basis']} {p['method']} {p['energy']:.10f}" for p in hr["points"]) +
          f"  CCSD/aug-cc-pVTZ minimum R = {hr['optimization']['R']:.6f} Å, E = {hr['optimization']['energy']:.10f}, "
          f"w = {hr['optimization']['wavenumber']:.2f} cm-1  ({time.time() - t0:.1f} s)", flush=True)
    for name, method in (("H2O2", "RHF"), ("NH3", "RHF"), ("NH3", "CCSD(T)")):
        t0 = time.time()
        mol = build(atoms(name), "aug-cc-pVDZ")
        mol_eq = optimize(orf.solver(mol, method), maxsteps=200, **orf.CONV)
        x = mol_eq.atom_coords(unit="Bohr")
        o = dict(molecule=name, basis="aug-cc-pVDZ", method=method, energy=float(orf.energy(mol_eq, method)),
                 geometry=[[float(c) for c in row] for row in x])
        o.update(sr.parameters(name, x))
        if method == "RHF":
            o["frequencies"] = sr.frequencies(mol_eq)
        o["seconds"] = round(time.time() - t0, 1)
        doc["optimizations"].append(o)
        keys = [k for k in ("ROO", "R1", "R2", "R3", "angle1", "angle2", "angle3", "dihedral") if k in o]
        print(f"{name:5s} aug-cc-pVDZ   {method:8s} optimised E = {o['energy']:.10f}  " +
              ", ".join(f"{k} {o[k]:.6f}" for k in keys) +
              ("  freq " + " ".join(f"{f:.2f}" for f in o["frequencies"]) if "frequencies" in o else "") +
              f"  ({o['seconds']} s)", flush=True)
    if args.check and stored:
        for a, b in zip(doc["points"], stored["points"]):
            worst["energy"] = max(worst["energy"], abs(a["energy"] - b["energy"]))
            if a["gradient"] is not None:
                worst["gradient"] = max(worst["gradient"], float(numpy.abs(numpy.array(a["gradient"]) -
                                                                             numpy.array(b["gradient"])).max()))
        for a, b in zip(doc["optimizations"], stored["optimizations"]):
            for k in ("ROO", "R1", "R2", "R3", "angle1", "angle2", "angle3", "dihedral"):
                if k in a:
                    worst["geometry"] = max(worst["geometry"], abs(a[k] - b[k]))
            for f, g in zip(a.get("frequencies", []), b.get("frequencies", [])):
                worst["frequency"] = max(worst["frequency"], abs(f - g))
        if stored.get("HeH+"):
            for a, b in zip(doc["HeH+"]["points"], stored["HeH+"]["points"]):
                worst["energy"] = max(worst["energy"], abs(a["energy"] - b["energy"]))
            worst["geometry"] = max(worst["geometry"], abs(doc["HeH+"]["optimization"]["R"] - stored["HeH+"]["optimization"]["R"]))
            worst["frequency"] = max(worst["frequency"], abs(doc["HeH+"]["optimization"]["wavenumber"] -
                                                             stored["HeH+"]["optimization"]["wavenumber"]))
        print("largest differences from the stored values:", json.dumps(worst))
        ok = worst["energy"] < 1e-9 and worst["gradient"] < 1e-8 and worst["geometry"] < 2e-5 and worst["frequency"] < 0.05
        return 0 if ok else 1
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
