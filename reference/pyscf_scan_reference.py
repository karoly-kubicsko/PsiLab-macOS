#!/usr/bin/env python3
"""Constrained optimisations for the relaxed scans of ΨLab 7.0 (PySCF + geomeTRIC).

Each point of a relaxed scan is a constrained minimum: the scanned bond or angle is fixed and
all other coordinates are optimised. Here every point is computed independently: PySCF gives the
energies and analytic gradients (RHF; CCSD(T) with the CCSD(T) Lambda equations, see
pyscf_optimization_reference.py), geomeTRIC (an independent optimiser in
translation-rotation-internal coordinates) minimises with the coordinate set by a "$set"
constraint, tight thresholds. Same basis files, all electrons correlated. The starting
geometries are ΨLab's defaults with the scanned coordinate at its value.

H2O with the angle frozen at exactly 180 degrees (ΨLab imposes it through the two linear-bend
components) is computed without geomeTRIC: the linear, symmetric molecule (D-infinity-h, which
the constrained minimum keeps by symmetry) is minimised along the symmetric stretch by Newton
steps on the PySCF analytic gradient.

Usage (needs: pip install pyscf geometric):
  python3 reference/pyscf_scan_reference.py           # -> backend/tests/data/scan_reference.json
  python3 reference/pyscf_scan_reference.py --check   # recompute and compare
"""
import argparse
import json
import math
import pathlib
import sys
import tempfile
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
import pyscf_optimization_reference as orf  # noqa: E402
import numpy  # noqa: E402
import pyscf  # noqa: E402
from pyscf import gto  # noqa: E402
from pyscf.geomopt.geometric_solver import optimize  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "scan_reference.json"
B = rr.BOHR_IN_ANGSTROM
# (molecule, basis, method, coordinate as in ΨLab, values in angstrom or degrees)
CASES = [
    ("H2O", "aug-cc-pVDZ", "RHF", "B 1 2", [0.85, 1.2]),
    ("H2O", "aug-cc-pVDZ", "RHF", "A 2 1 3", [90.0, 150.0, 170.0, 180.0]),
    ("H2O", "aug-cc-pVDZ", "CCSD(T)", "B 1 2", [1.1]),
]


def start_atoms(name, coord, value):
    kind = coord.split()[0]
    geo = rr.GEOMETRIES[name]
    if kind == "A":
        return rr.atoms(name, geo["R"], value)
    at = rr.atoms(name, geo["R"], geo["angle"])
    if kind == "B":  # move H2 along O-H2 to the value
        o = numpy.array(at[0][1])
        h = numpy.array(at[1][1])
        u = (h - o) / numpy.linalg.norm(h - o)
        at[1] = ("H", tuple(o + u * value / B))
    return at


def build(atoms, basis):
    bas = {s: rr.element_basis(basis, s)[0] for s, _ in atoms}
    return gto.M(atom=[[s, xyz] for s, xyz in atoms], basis=bas, unit="Bohr", cart=False, verbose=0)


def internal(x, kind, idx):
    p = [x[i - 1] for i in idx]
    if kind == "B":
        return float(numpy.linalg.norm(p[1] - p[0]) * B)
    if kind == "A":
        u, v = p[0] - p[1], p[2] - p[1]
        return math.degrees(math.acos(numpy.dot(u, v) / (numpy.linalg.norm(u) * numpy.linalg.norm(v))))
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    e = b1 / numpy.linalg.norm(b1)
    v, w = b0 - numpy.dot(b0, e) * e, b2 - numpy.dot(b2, e) * e
    return math.degrees(math.atan2(numpy.dot(numpy.cross(e, v), w), numpy.dot(v, w)))


def parameters(name, x):
    return dict(R1=internal(x, "B", [1, 2]), R2=internal(x, "B", [1, 3]), angle=internal(x, "A", [2, 1, 3]))


def linear_gradient(r, basis):
    """Energy and dE/dr of linear symmetric H2O (O at the origin, H at z = +r and -r), RHF."""
    from pyscf import scf
    mol = build([("O", (0.0, 0.0, 0.0)), ("H", (0.0, 0.0, r)), ("H", (0.0, 0.0, -r))], basis)
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-9
    mf.kernel()
    g = mf.nuc_grad_method().kernel()
    return mol, float(mf.e_tot), float(g[1][2] - g[2][2])


def linear_water(basis, method):
    """Linear symmetric H2O (H-O-H = 180 degrees) at the minimum along the symmetric stretch."""
    if method != "RHF":
        raise ValueError("linear_water: RHF only")
    r = 0.93 / B
    for _ in range(30):
        mol, e, d1 = linear_gradient(r, basis)
        if abs(d1) < 1e-10:
            break
        h = 1e-4
        d2 = (linear_gradient(r + h, basis)[2] - d1) / h
        r -= d1 / d2
    return mol, e


def run_point(name, basis, method, coord, value):
    t0 = time.time()
    tok = coord.split()
    kind, idx = tok[0], [int(t) for t in tok[1:]]
    if kind == "A" and value == 180.0:
        mol_eq, e = linear_water(basis, method)
        x = mol_eq.atom_coords(unit="Bohr")
        out = dict(value=value, energy=e, measured=180.0, geometry=[[float(c) for c in row] for row in x])
        out.update(R1=float(x[1][2] * B), R2=float(-x[2][2] * B), angle=180.0)
        out["seconds"] = round(time.time() - t0, 1)
        return out
    mol = build(start_atoms(name, coord, value), basis)
    word = {"B": "distance", "A": "angle", "D": "dihedral"}[kind]
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write("$set\n" + word + " " + " ".join(map(str, idx)) + f" {value:.10f}\n")
        cfile = f.name
    mol_eq = optimize(orf.solver(mol, method), maxsteps=200, constraints=cfile, **orf.CONV)
    x = mol_eq.atom_coords(unit="Bohr")
    out = dict(value=value, energy=float(orf.energy(mol_eq, method)), measured=internal(x, kind, idx),
               geometry=[[float(c) for c in row] for row in x])
    out.update(parameters(name, x))
    out["seconds"] = round(time.time() - t0, 1)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--tolerance", type=float, default=2e-5)
    args = ap.parse_args()
    stored = json.loads(OUT.read_text()) if OUT.exists() else None
    cases, worst = [], 0.0
    for name, basis, method, coord, values in CASES:
        pts = []
        for v in values:
            p = run_point(name, basis, method, coord, v)
            pts.append(p)
            keys = [k for k in ("R1", "R2", "angle") if k in p]
            line = (f"{name:4s} {basis} {method:8s} {coord} = {v:8.3f}: E = {p['energy']:.10f}  " +
                    ", ".join(f"{k} {p[k]:.6f}" for k in keys) + f"  ({p['seconds']} s)")
            if args.check and stored:
                ref = next(c for c in stored["cases"] if (c["molecule"], c["method"], c["coordinate"]) == (name, method, coord))
                rp = next(q for q in ref["points"] if q["value"] == v)
                d = max(abs(rp[k] - p[k]) for k in keys)
                worst = max(worst, d)
                line += f"  |diff| {d:.1e}"
            print(line, flush=True)
        cases.append(dict(molecule=name, basis=basis, method=method, coordinate=coord, points=pts))
    if args.check:
        print(f"largest difference from the stored values: {worst:.1e}")
        return 0 if worst <= args.tolerance else 1
    doc = {
        "description": "Constrained minima (scan points) from PySCF analytic gradients and geomeTRIC with a $set "
                       "constraint (H2O at 180 degrees: Newton minimisation of the linear molecule along the symmetric "
                       "stretch): bond lengths in angstrom, angles in degrees, energies in E_h, geometry in bohr; "
                       "generated by reference/pyscf_scan_reference.py",
        "program": f"PySCF {pyscf.__version__}, geomeTRIC",
        "settings": "all electrons correlated; geomeTRIC thresholds " + json.dumps(orf.CONV),
        "bohr_in_angstrom": B,
        "cases": cases,
    }
    OUT.write_text(json.dumps(doc, indent=1) + "\n")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
