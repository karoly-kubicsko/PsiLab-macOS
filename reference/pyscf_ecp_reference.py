#!/usr/bin/env python3
"""Reference values for the elements beyond Kr of ΨLab 14.0: the def2 sets with their effective core
potentials (ECPs) and all-electron iodine in 6-311G*, computed with PySCF from the same files
(basis/*.gbs: the basis functions and the ECP blocks, read here into PySCF's format).

ECPs covered: def2-ECP with 28 (I, Xe), 46 (Cs) and 60 (Au, Pb) core electrons, the 4f-in-valence
ECP28 of the lanthanides (Lu: projectors up to g), and the small-core ECP60MWB of the actinides
(Th in def2-mTZVPP).

Stored (backend/tests/data/ecp_reference.json, read by backend/tests/test_basis_sets.cpp):
  - RHF energies; CCSD(T) of HI (all valence electrons of the ECP calculation correlated);
  - analytic RHF gradients (HI, CsF, PbO, LuF3 away from their minima) and the CCSD gradient of HI;
  - RHF Hessians of HI and PbO.

Usage (needs: pip install pyscf):
  python3 reference/pyscf_ecp_reference.py           # -> backend/tests/data/ecp_reference.json
  python3 reference/pyscf_ecp_reference.py --check   # recompute and compare
"""
import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
from pyscf import cc, gto, scf  # noqa: E402
from pyscf.gto.basis import parse_gaussian  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "ecp_reference.json"

FILES = {"def2-SVP": "def2-svp.gbs", "def2-TZVP": "def2-tzvp.gbs", "def2-TZVPP": "def2-tzvpp.gbs",
         "def2-TZVPPD": "def2-tzvppd.gbs", "def2-QZVPP": "def2-qzvpp.gbs", "def2-mTZVP": "def2-mtzvp.gbs",
         "def2-mTZVPP": "def2-mtzvpp.gbs", "6-311G*": "6-311gs.gbs"}

# name: (title, charge, XYZ atom lines in angstrom)
STRUCTURES = {
    "HI": ("HI", 0, "H 0 0 0\nI 0 0 1.609"),
    "HI-stretched": ("HI", 0, "H 0 0 0\nI 0 0 1.70"),
    "Xe": ("Xe", 0, "Xe 0 0 0"),
    "CsF": ("CsF", 0, "Cs 0 0 0\nF 0 0 2.40"),
    "AuH": ("AuH", 0, "Au 0 0 0\nH 0 0 1.524"),
    "PbO": ("PbO", 0, "Pb 0 0 0\nO 0 0 1.95"),
    "LuF3": ("LuF$_3$", 0, "Lu 0 0 0\nF 2.00 0 0\nF -0.98 1.70 0\nF -0.99 -1.71 0.05"),
    "ThO2": ("ThO$_2$", 0, "Th 0 0 0\nO 1.66 0 1.06\nO -1.66 0 1.06"),
    "ICl": ("ICl", 0, "I 0 0 0\nCl 0 0 2.321"),
}

ENERGIES = [
    ("HI", "def2-SVP", "RHF"), ("HI", "def2-TZVP", "RHF"), ("HI", "def2-QZVPP", "RHF"),
    ("HI", "def2-TZVPPD", "RHF"), ("HI", "def2-mTZVP", "RHF"), ("HI", "6-311G*", "RHF"),
    ("Xe", "def2-QZVPP", "RHF"), ("CsF", "def2-TZVP", "RHF"), ("AuH", "def2-TZVPP", "RHF"),
    ("PbO", "def2-SVP", "RHF"), ("LuF3", "def2-SVP", "RHF"), ("ThO2", "def2-mTZVPP", "RHF"),
    ("ICl", "def2-SVP", "RHF"),
    ("HI", "def2-SVP", "CCSD(T)"),
]
GRADIENTS = [
    ("HI-stretched", "def2-TZVP", "RHF"), ("CsF", "def2-SVP", "RHF"), ("PbO", "def2-SVP", "RHF"),
    ("LuF3", "def2-SVP", "RHF"), ("HI-stretched", "def2-SVP", "CCSD"),
]
HESSIANS = [("HI", "def2-SVP"), ("PbO", "def2-SVP")]


def read_ecp(path, symbol):
    """The ECP block of `symbol` in a Gaussian94 file as PySCF's [ncore, [[l, [r^0, r^1, r^2 terms]], ...]]
    (l = -1: the local part), or None."""
    lines = [l.strip() for l in open(path) if l.strip() and not l.lstrip().startswith("!")]
    for i, line in enumerate(lines):
        t = line.split()
        if t[0].upper() == symbol.upper() + "-ECP":
            L, ncore = int(t[1]), int(t[2])
            k = i + 1
            parts = []
            for part in range(L + 1):
                k += 1  # label
                nterm = int(lines[k]); k += 1
                terms = [[], [], []]
                for _ in range(nterm):
                    n, a, c = lines[k].split(); k += 1
                    terms[int(n)].append([float(a.replace("D", "E")), float(c.replace("D", "E"))])
                parts.append([-1 if part == 0 else part - 1, terms])
            return [ncore, parts]
    return None


def parse_atoms(text):
    return [[s, [float(x), float(y), float(z)]] for s, x, y, z in (l.split() for l in text.strip().splitlines())]


def build(name, basis):
    _title, charge, text = STRUCTURES[name]
    at = parse_atoms(text)
    path = str(rr.BASIS_DIR / FILES[basis])
    symbols = {s for s, _ in at}
    bas = {s: parse_gaussian.load(path, s) for s in symbols}
    ecp = {s: e for s in symbols if (e := read_ecp(path, s)) is not None}
    atoms_bohr = [[s, [c / rr.BOHR_IN_ANGSTROM for c in xyz]] for s, xyz in at]
    return gto.M(atom=atoms_bohr, basis=bas, ecp=ecp, unit="Bohr", cart=False, charge=charge, spin=0, verbose=0)


def rhf(mol):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-8
    mf.direct_scf_tol = 1e-14
    mf.max_cycle = 300
    mf.kernel()
    if not mf.converged:
        mf = scf.newton(mf)
        mf.kernel()
    if not mf.converged:
        raise SystemExit("RHF did not converge")
    return mf


def compute():
    data = {"structures": {}, "energies": [], "gradients": [], "hessians": []}
    for name, (title, charge, text) in STRUCTURES.items():
        data["structures"][name] = {"title": title, "charge": charge,
                                    "xyz": f"{len(text.strip().splitlines())}\n{title}\n{text}\n"}
    for name, basis, method in ENERGIES:
        mol = build(name, basis)
        mf = rhf(mol)
        e = mf.e_tot
        if method == "CCSD(T)":
            c = cc.CCSD(mf)
            c.conv_tol = 1e-11
            c.conv_tol_normt = 1e-8
            c.kernel()
            e = c.e_tot + c.ccsd_t()
        core = sum(mol.atom_nelec_core(i) for i in range(mol.natm))
        data["energies"].append({"molecule": name, "basis": basis, "method": method, "energy": e,
                                 "basisFunctions": mol.nao, "electrons": mol.nelectron, "coreElectrons": core})
        print(f"{name:13s} {basis:12s} {method:8s} {e:.10f}  ({mol.nao} functions, {mol.nelectron} electrons, "
              f"{core} in ECPs)", flush=True)
    for name, basis, method in GRADIENTS:
        mol = build(name, basis)
        mf = rhf(mol)
        if method == "RHF":
            g = mf.nuc_grad_method().kernel()
            e = mf.e_tot
        else:
            c = cc.CCSD(mf)
            c.conv_tol = 1e-11
            c.conv_tol_normt = 1e-9
            c.kernel()
            g = c.nuc_grad_method().kernel()
            e = c.e_tot
        data["gradients"].append({"molecule": name, "basis": basis, "method": method, "energy": e,
                                  "gradient": [float(v) for v in g.ravel()]})
        print(f"{name:13s} {basis:12s} {method:8s} gradient, |g|max {abs(g).max():.6f}", flush=True)
    for name, basis in HESSIANS:
        mol = build(name, basis)
        mf = rhf(mol)
        h = mf.Hessian().kernel()
        n = mol.natm
        H = np.zeros((3 * n, 3 * n))
        for a in range(n):
            for b in range(n):
                H[3 * a:3 * a + 3, 3 * b:3 * b + 3] = h[a, b]
        data["hessians"].append({"molecule": name, "basis": basis, "method": "RHF", "hessian": H.ravel().tolist()})
        print(f"{name:13s} {basis:12s} RHF      Hessian", flush=True)
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    data = compute()
    if a.check:
        ref = json.loads(OUT.read_text())
        worst = 0.0
        for x, y in zip(data["energies"], ref["energies"]):
            worst = max(worst, abs(x["energy"] - y["energy"]))
        for key, field in (("gradients", "gradient"), ("hessians", "hessian")):
            for x, y in zip(data[key], ref[key]):
                worst = max(worst, max(abs(u - v) for u, v in zip(x[field], y[field])))
        print(f"largest difference from the stored values: {worst:.2e}")
        if worst > 1e-8 or data["structures"] != ref["structures"]:
            raise SystemExit("reference values changed")
        return
    OUT.write_text(json.dumps(data, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
