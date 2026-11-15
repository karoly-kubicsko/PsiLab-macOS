#!/usr/bin/env python3
"""Reference values for the molecular orbitals of ΨLab 7.0, computed with PySCF.

Three checks that do not depend on how a program orders or signs its real spherical
harmonics, or on the arbitrary sign of an orbital:

1. Basis functions (l = 0 ... 5): for every contracted shell, the kernel
   K(r, r') = sum_m phi_lm(r) phi_lm(r') at pairs of points. It is the same for any
   orthonormal set of real harmonics, so it checks the radial part, the normalisation and
   the angular functions of the shells. Atoms O and Li (aug-cc-pV5Z), N (aug-cc-pCV5Z), He and
   H (aug-cc-pV5Z), placed off the origin; l goes up to 5 (h).
2. RHF canonical orbitals of H2O, NH3 and H2O2 (aug-cc-pVDZ, the default geometries of ΨLab):
   orbital energies and the values of all occupied orbitals and the lowest four virtual ones
   at points around the molecule (compared up to the sign; a degenerate set, e.g. the e
   orbitals of NH3, through the sum of the squares of its orbitals).
3. Full CI natural orbitals of LiH (aug-cc-pVDZ): the natural occupations and the Full CI
   electron density sum_k n_k chi_k(r)^2 at the points (independent of the phases).

Usage:
  python3 reference/pyscf_orbitals_reference.py           # -> backend/tests/data/orbitals_reference.json
  python3 reference/pyscf_orbitals_reference.py --check   # recompute and compare with the stored file

Needs: pip install pyscf
"""
import argparse
import json
import math
import pathlib
import sys

import numpy
import pyscf
from pyscf import fci, gto, scf

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "orbitals_reference.json"
B = rr.BOHR_IN_ANGSTROM
H2O2 = dict(ROO=1.4556 / B, R=0.9670 / B, angle=102.32, dihedral=113.70)
NH3 = dict(R=1.012 / B, angle=106.7)


def h2o2_atoms(roo, r, ang, dih):
    a, h = math.radians(ang), math.radians(dih) / 2
    c, s = math.cos(a), math.sin(a)
    return [("O", (-roo / 2, 0.0, 0.0)), ("O", (roo / 2, 0.0, 0.0)),
            ("H", (-roo / 2 + r * c, r * s * math.sin(h), r * s * math.cos(h))),
            ("H", (roo / 2 - r * c, -r * s * math.sin(h), r * s * math.cos(h)))]


def nh3_atoms(r, ang):
    c2 = max(0.0, (2.0 * math.cos(math.radians(ang)) + 1.0) / 3.0)
    ct, st = math.sqrt(c2), math.sqrt(1.0 - c2)
    return [("N", (0.0, 0.0, 0.0))] + [("H", (r * st * math.cos(p), r * st * math.sin(p), r * ct))
                                       for p in (0.0, 2.0 * math.pi / 3.0, 4.0 * math.pi / 3.0)]


MOLECULES = {
    "H2O": rr.atoms("H2O", rr.GEOMETRIES["H2O"]["R"], rr.GEOMETRIES["H2O"]["angle"]),
    "NH3": nh3_atoms(NH3["R"], NH3["angle"]),
    "H2O2": h2o2_atoms(H2O2["ROO"], H2O2["R"], H2O2["angle"], H2O2["dihedral"]),
}

SHELL_CASES = [("O", "aug-cc-pV5Z"), ("Li", "aug-cc-pV5Z"), ("N", "aug-cc-pCV5Z"),
               ("He", "aug-cc-pV5Z"), ("H", "aug-cc-pV5Z")]
SHELL_CENTER = (0.31, -0.22, 0.17)


def build(at, basis):
    bas = {s: rr.element_basis(basis, s)[0] for s, _ in at}
    spin = sum(gto.charge(s) for s, _ in at) % 2  # single atoms (shell kernels only)
    return gto.M(atom=[[s, xyz] for s, xyz in at], basis=bas, unit="Bohr", cart=False, spin=spin, verbose=0)


def points_around(at, n, seed):
    """n points within 3 bohr of the atoms (some near the nuclei, some in the tails)."""
    rng = numpy.random.default_rng(seed)
    xyz = numpy.array([p for _, p in at])
    pts = []
    for k in range(n):
        a = xyz[k % len(xyz)]
        d = rng.normal(size=3)
        d /= numpy.linalg.norm(d)
        r = [0.05, 0.3, 0.8, 1.5, 2.5, 3.5][k % 6] * (0.8 + 0.4 * rng.random())
        pts.append(a + r * d)
    return numpy.array(pts)


def shell_kernels():
    out = []
    rng = numpy.random.default_rng(7)
    for sym, basis in SHELL_CASES:
        mol = build([(sym, SHELL_CENTER)], basis)
        pts = numpy.array(SHELL_CENTER) + numpy.array(
            [r * v / numpy.linalg.norm(v) for r, v in zip([0.07, 0.4, 0.9, 1.6, 2.7, 4.0],
                                                          rng.normal(size=(6, 3)))])
        ao = mol.eval_gto("GTOval_sph", pts)  # points x AO
        loc = mol.ao_loc_nr()
        shells = []
        for b in range(mol.nbas):
            l = int(mol.bas_angular(b))
            nctr = int(mol.bas_nctr(b))
            exps = [float(x) for x in mol.bas_exp(b)]
            coefs = mol.bas_ctr_coeff(b)  # nprim x nctr
            width = 2 * l + 1
            for c in range(nctr):
                cols = ao[:, loc[b] + c * width: loc[b] + (c + 1) * width]
                K = cols @ cols.T
                shells.append(dict(l=l, exponents=exps, coefficients=[float(x) for x in coefs[:, c]],
                                   kernel=[[float(x) for x in row] for row in K]))
        out.append(dict(element=sym, basis=basis, center=list(SHELL_CENTER), points=pts.tolist(), shells=shells))
        print(f"{sym:2s} {basis:13s} {len(shells)} shells, l up to {max(s['l'] for s in shells)}")
    return out


def rhf(mol):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.direct_scf_tol = 1e-14
    mf.kernel()
    assert mf.converged
    return mf


def canonical_cases():
    out = []
    for name, at in MOLECULES.items():
        mol = build(at, "aug-cc-pVDZ")
        mf = rhf(mol)
        nocc = mol.nelectron // 2
        pts = points_around(at, 24, {"H2O": 1, "NH3": 2, "H2O2": 3}[name])
        ao = mol.eval_gto("GTOval_sph", pts)
        n = nocc + 4
        values = ao @ mf.mo_coeff[:, :n]
        out.append(dict(molecule=name, basis="aug-cc-pVDZ", energy=float(mf.e_tot), occupied=nocc,
                        orbitalEnergies=[float(e) for e in mf.mo_energy[:n]], points=pts.tolist(),
                        values=[[float(x) for x in row] for row in values.T]))
        print(f"{name:5s} RHF/aug-cc-pVDZ E = {mf.e_tot:.10f}, {n} orbitals at {len(pts)} points")
    return out


def natural_case():
    at = rr.atoms("LiH", rr.GEOMETRIES["LiH"]["R"])
    mol = build(at, "aug-cc-pVDZ")
    mf = rhf(mol)
    solver = fci.FCI(mf)
    solver.conv_tol = 1e-12
    e, civec = solver.kernel()
    dm_mo = solver.make_rdm1(civec, mol.nao, mol.nelectron)
    occ = numpy.sort(numpy.linalg.eigvalsh(dm_mo))[::-1]
    dm_ao = mf.mo_coeff @ dm_mo @ mf.mo_coeff.T
    pts = points_around(at, 18, 4)
    ao = mol.eval_gto("GTOval_sph", pts)
    rho = numpy.einsum("pi,ij,pj->p", ao, dm_ao, ao)
    print(f"LiH   FCI/aug-cc-pVDZ E = {e:.10f}, occupations {occ[:4]}")
    return dict(molecule="LiH", basis="aug-cc-pVDZ", energy=float(e), occupations=[float(x) for x in occ],
                points=pts.tolist(), density=[float(x) for x in rho])


def compute():
    return dict(pyscf=pyscf.__version__, shells=shell_kernels(), canonical=canonical_cases(),
                natural=natural_case())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    data = compute()
    if not args.check:
        OUT.write_text(json.dumps(data, indent=1) + "\n")
        print(f"wrote {OUT}")
        return
    old = json.loads(OUT.read_text())
    worst = 0.0
    for a, b in zip(old["shells"], data["shells"]):
        for sa, sb in zip(a["shells"], b["shells"]):
            worst = max(worst, float(numpy.max(numpy.abs(numpy.array(sa["kernel"]) - numpy.array(sb["kernel"])))
                                     / (1.0 + numpy.max(numpy.abs(sb["kernel"])))))
    for a, b in zip(old["canonical"], data["canonical"]):
        worst = max(worst, abs(a["energy"] - b["energy"]))
        ea, eb = numpy.array(a["orbitalEnergies"]), numpy.array(b["orbitalEnergies"])
        worst = max(worst, float(numpy.max(numpy.abs(ea - eb))))
        va, vb = numpy.array(a["values"]), numpy.array(b["values"])
        # the density of the occupied orbitals (independent of rotations among degenerate ones)
        nocc = a["occupied"]
        rho_a, rho_b = 2 * (va[:nocc] ** 2).sum(0), 2 * (vb[:nocc] ** 2).sum(0)
        worst = max(worst, float(numpy.max(numpy.abs(rho_a - rho_b)) / numpy.max(rho_b)))
        # single orbitals up to the sign, where they are not degenerate
        for k in range(len(eb)):
            if any(abs(eb[k] - eb[q]) < 1e-5 for q in range(len(eb)) if q != k):
                continue
            scale = float(numpy.max(numpy.abs(vb[k])))
            diff = min(float(numpy.max(numpy.abs(va[k] - vb[k]))), float(numpy.max(numpy.abs(va[k] + vb[k]))))
            worst = max(worst, diff / scale * 1e-2)
    na, nb = old["natural"], data["natural"]
    worst = max(worst, abs(na["energy"] - nb["energy"]),
                float(numpy.max(numpy.abs(numpy.array(na["density"]) - numpy.array(nb["density"])))))
    print(f"largest difference from the stored references: {worst:.2e}")
    if worst > 1e-7:
        sys.exit(1)


if __name__ == "__main__":
    main()
