#!/usr/bin/env python3
"""Reference values for the custom molecules of ΨLab 9.0 (Cartesian coordinates, any closed-shell
molecule of the elements H-Ar and Ga-Kr), computed with PySCF from the same basis-set files.

The structures are given here as XYZ text (angstrom, as a user's .xyz file would give them) and
stored with the results, so the C++ test (backend/tests/test_custom.cpp) reads the same atoms:
  - RHF energies in aug-cc-pVDZ (and two in aug-cc-pCVDZ, where H takes aug-cc-pVDZ, as in
    ΨLab) of molecules with the elements added in 9.0 (Be, B, C, F, Na, Si, P, S, Cl, Cu, Zn, Br, Kr),
    a cation and an anion, a single atom, and molecules whose atoms are listed in an unusual order;
  - CCSD(T)/aug-cc-pVDZ energies (all electrons correlated) of HCl and CH4;
  - the analytic RHF/aug-cc-pVDZ gradient of HOCl near (not at) its minimum.

Usage (needs: pip install pyscf):
  python3 reference/pyscf_custom_reference.py           # -> backend/tests/data/custom_reference.json
  python3 reference/pyscf_custom_reference.py --check   # recompute and compare
"""
import argparse
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
from pyscf import cc, gto, scf  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "custom_reference.json"

# name: (title line, charge, XYZ atom lines in angstrom)
STRUCTURES = {
    "CO2": ("CO$_2$", 0, "C 0 0 0\nO 0 0 1.16\nO 0 0 -1.16"),
    "CH4": ("CH$_4$", 0, "C 0 0 0\nH 0.629118 0.629118 0.629118\nH -0.629118 -0.629118 0.629118\n"
                         "H -0.629118 0.629118 -0.629118\nH 0.629118 -0.629118 -0.629118"),
    "C2H4": ("C$_2$H$_4$", 0, "C 0 0 0.6695\nC 0 0 -0.6695\nH 0 0.9289 1.2321\nH 0 -0.9289 1.2321\n"
                              "H 0 0.9289 -1.2321\nH 0 -0.9289 -1.2321"),
    "HCN": ("HCN", 0, "H 0 0 -1.0655\nC 0 0 0\nN 0 0 1.1532"),
    "BeH2": ("BeH$_2$", 0, "Be 0 0 0\nH 0 0 1.326\nH 0 0 -1.326"),
    "BH3": ("BH$_3$", 0, "B 0 0 0\nH 1.19 0 0\nH -0.595 1.030570 0\nH -0.595 -1.030570 0"),
    "LiF": ("LiF", 0, "Li 0 0 0\nF 0 0 1.5639"),
    "NaH": ("NaH", 0, "Na 0 0 0\nH 0 0 1.8874"),
    "SiH4": ("SiH$_4$", 0, "Si 0 0 0\nH 0.854364 0.854364 0.854364\nH -0.854364 -0.854364 0.854364\n"
                           "H -0.854364 0.854364 -0.854364\nH 0.854364 -0.854364 -0.854364"),
    "PH3": ("PH$_3$", 0, "P 0 0 0.1272\nH 1.1937 0 -0.5937\nH -0.59685 1.033775 -0.5937\nH -0.59685 -1.033775 -0.5937"),
    "H2S": ("H$_2$S", 0, "H 0.9616 0 0.9269\nS 0 0 0\nH -0.9616 0 0.9269"),
    "HCl": ("HCl", 0, "H 0 0 0\nCl 0 0 1.2746"),
    "HBr": ("HBr", 0, "Br 0 0 0\nH 0 0 1.4145"),
    "Kr": ("Kr", 0, "Kr 0 0 0"),
    "Zn": ("Zn", 0, "Zn 0 0 0"),
    "ZnH2": ("ZnH$_2$", 0, "Zn 0 0 0\nH 0 0 1.535\nH 0 0 -1.535"),
    "CuH": ("CuH", 0, "Cu 0 0 0\nH 0 0 1.4626"),
    "NH4+": ("NH$_4^+$", 1, "N 0 0 0\nH 0.590 0.590 0.590\nH -0.590 -0.590 0.590\nH -0.590 0.590 -0.590\n"
                           "H 0.590 -0.590 -0.590"),
    "OH-": ("OH$^-$", -1, "O 0 0 0\nH 0 0 0.9697"),
    "HOCl": ("HOCl", 0, "O 0 0 0\nH 0.9257 0 -0.2137\nCl 0 0 1.67"),
}

CASES = [(n, "aug-cc-pVDZ", "RHF") for n in STRUCTURES] + [
    ("HCl", "aug-cc-pCVDZ", "RHF"),
    ("SiH4", "aug-cc-pCVDZ", "RHF"),
    ("HCl", "aug-cc-pVDZ", "CCSD(T)"),
    ("CH4", "aug-cc-pVDZ", "CCSD(T)"),
]


def parse_atoms(text):
    out = []
    for line in text.strip().splitlines():
        s, x, y, z = line.split()
        out.append([s, [float(x), float(y), float(z)]])
    return out


def build(name, basis):
    _title, charge, text = STRUCTURES[name]
    at = parse_atoms(text)
    bas = {}
    for s, _ in at:
        b, _used = rr.element_basis(basis, s)
        if b is None:
            raise SystemExit(f"{basis} is not defined for {s}")
        bas[s] = b
    # angstrom -> bohr with ΨLab's factor, so the geometries are bit-for-bit the same
    atoms_bohr = [[s, [c / rr.BOHR_IN_ANGSTROM for c in xyz]] for s, xyz in at]
    return gto.M(atom=atoms_bohr, basis=bas, unit="Bohr", cart=False, charge=charge, spin=0, verbose=0)


def rhf(mol):
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-8
    mf.direct_scf_tol = 1e-14
    mf.max_cycle = 200
    mf.kernel()
    if not mf.converged:
        raise SystemExit("RHF did not converge")
    return mf


def compute():
    data = {"structures": {}, "energies": [], "gradients": []}
    for name, (title, charge, text) in STRUCTURES.items():
        data["structures"][name] = {"title": title, "charge": charge,
                                    "xyz": f"{len(text.strip().splitlines())}\n{title}\n{text}\n"}
    for name, basis, method in CASES:
        mol = build(name, basis)
        mf = rhf(mol)
        e = mf.e_tot
        if method == "CCSD(T)":
            c = cc.CCSD(mf)
            c.conv_tol = 1e-11
            c.conv_tol_normt = 1e-8
            c.kernel()
            e = c.e_tot + c.ccsd_t()
        data["energies"].append({"molecule": name, "basis": basis, "method": method, "energy": e,
                                 "basisFunctions": mol.nao})
        print(f"{name:6s} {basis:13s} {method:8s} {e:.10f}  ({mol.nao} functions)", flush=True)
    mol = build("HOCl", "aug-cc-pVDZ")
    mf = rhf(mol)
    g = mf.nuc_grad_method().kernel()
    data["gradients"].append({"molecule": "HOCl", "basis": "aug-cc-pVDZ", "method": "RHF",
                              "gradient": [float(v) for v in g.ravel()]})
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
        for x, y in zip(data["gradients"], ref["gradients"]):
            worst = max(worst, max(abs(u - v) for u, v in zip(x["gradient"], y["gradient"])))
        print(f"largest difference from the stored values: {worst:.2e}")
        if worst > 1e-8 or data["structures"] != ref["structures"]:
            raise SystemExit("reference values changed")
        return
    OUT.write_text(json.dumps(data, indent=1) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
