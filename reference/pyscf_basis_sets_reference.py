#!/usr/bin/env python3
"""Reference values for the basis sets added in ΨLab 14.0 (the Pople 6-31G and 6-311G families and
the Karlsruhe def2 sets), computed with PySCF from the same basis-set files (basis/*.gbs).

Function types, as Gaussian (and Psi4, Q-Chem) use these sets by default and as ΨLab does:
  - the 6-31G family: Cartesian d functions (6D), spherical f functions (7F); PySCF runs with
    cart=True; the one case with both (Zn in 6-31G*: 6D 7F) is computed from PySCF's Cartesian
    integrals transformed to spherical functions in the f shells only;
  - the 6-311G family and def2: spherical harmonics (5D 7F), cart=False.

Stored (backend/tests/data/basis_sets_reference.json, read by backend/tests/test_basis_sets.cpp):
  - RHF energies of H2O in all 25 sets and of molecules with other elements (Li, C, N, Cl, K,
    Ca, Zn, Br; f and g functions; 6D with K, Br and the 6D 7F Zn case);
  - CCSD(T) energies (all electrons correlated) in 6-31G*, 6-311++G** and def2-SVP;
  - analytic RHF gradients (6-31+G**, 6-31G* with Cl, def2-TZVP) and the CCSD gradient of H2O in
    6-31G* away from the minimum; the RHF Hessian and dipole moment of H2O in 6-31G*;
  - values of the RHF HOMO and LUMO of H2O in 6-31G* at points (MO values do not depend on the
    normalisation convention of the Cartesian functions).

Usage (needs: pip install pyscf):
  python3 reference/pyscf_basis_sets_reference.py           # -> backend/tests/data/basis_sets_reference.json
  python3 reference/pyscf_basis_sets_reference.py --check   # recompute and compare
"""
import argparse
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pyscf_rhf_reference as rr  # noqa: E402
from pyscf import ao2mo, cc, gto, scf  # noqa: E402
from pyscf.gto.basis import parse_gaussian  # noqa: E402

OUT = rr.ROOT / "backend" / "tests" / "data" / "basis_sets_reference.json"

# ΨLab's file names (Psi4's convention: + -> p, * -> s, ( , ) -> _) and the 6D sets.
FILES = {
    "6-31G": "6-31g.gbs", "6-31G*": "6-31gs.gbs", "6-31+G": "6-31pg.gbs", "6-31+G*": "6-31pgs.gbs",
    "6-31+G**": "6-31pgss.gbs", "6-31++G*": "6-31ppgs.gbs", "6-31++G**": "6-31ppgss.gbs",
    "6-311G": "6-311g.gbs", "6-311G*": "6-311gs.gbs", "6-311+G": "6-311pg.gbs", "6-311+G*": "6-311pgs.gbs",
    "6-311+G**": "6-311pgss.gbs", "6-311++G*": "6-311ppgs.gbs", "6-311++G**": "6-311ppgss.gbs",
    "6-311++G(2d,2p)": "6-311ppg_2d_2p_.gbs", "6-311++G(3df,3pd)": "6-311ppg_3df_3pd_.gbs",
    "def2-SVP": "def2-svp.gbs", "def2-mTZVP": "def2-mtzvp.gbs", "def2-mTZVPP": "def2-mtzvpp.gbs",
    "def2-TZVP": "def2-tzvp.gbs", "def2-TZVPP": "def2-tzvpp.gbs", "def2-TZVPPD": "def2-tzvppd.gbs",
    "def2-QZVP": "def2-qzvp.gbs", "def2-QZVPP": "def2-qzvpp.gbs", "def2-QZVPPD": "def2-qzvppd.gbs",
}


def cartesian_d(basis):
    return basis.startswith("6-31") and not basis.startswith("6-311")


# name: (title, charge, XYZ atom lines in angstrom)
STRUCTURES = {
    "H2O": ("H$_2$O", 0, "O 0 0 0\nH 0.75695033 0 0.58588228\nH -0.75695033 0 0.58588228"),
    "H2O-distorted": ("H$_2$O", 0, "O 0 0 0.02\nH 0.80 0.03 0.55\nH -0.73 0 0.61"),
    "LiH": ("LiH", 0, "Li 0 0 0\nH 0 0 1.5949"),
    "NH3": ("NH$_3$", 0, "N 0 0 0.1173\nH 0.9377 0 -0.2737\nH -0.46885 0.81207 -0.2737\nH -0.46885 -0.81207 -0.2737"),
    "CH4": ("CH$_4$", 0, "C 0 0 0\nH 0.629118 0.629118 0.629118\nH -0.629118 -0.629118 0.629118\n"
                         "H -0.629118 0.629118 -0.629118\nH 0.629118 -0.629118 -0.629118"),
    "HCl": ("HCl", 0, "H 0 0 0\nCl 0 0 1.2746"),
    "HCl-stretched": ("HCl", 0, "H 0 0 0\nCl 0 0 1.35"),
    "KH": ("KH", 0, "K 0 0 0\nH 0 0 2.2423"),
    "CaH2": ("CaH$_2$", 0, "Ca 0 0 0\nH 0 0 2.04\nH 0 0 -2.04"),
    "ZnH2": ("ZnH$_2$", 0, "Zn 0 0 0\nH 0 0 1.535\nH 0 0 -1.535"),
    "HBr": ("HBr", 0, "Br 0 0 0\nH 0 0 1.4145"),
    "HeH+": ("HeH$^+$", 1, "He 0 0 0\nH 0 0 0.7743"),
}

ENERGIES = [("H2O", b, "RHF") for b in FILES] + [
    ("LiH", "6-31++G**", "RHF"),
    ("NH3", "6-311++G(3df,3pd)", "RHF"),
    ("CH4", "6-31+G*", "RHF"),
    ("CH4", "def2-TZVPPD", "RHF"),
    ("HCl", "6-31G*", "RHF"),          # 6D on Cl
    ("HCl", "6-311G*", "RHF"),         # 5D on Cl (the BSE marks these d functions Cartesian; Gaussian: 5D)
    ("HCl", "def2-QZVP", "RHF"),
    ("KH", "6-31G", "RHF"),            # Cartesian valence d of K
    ("CaH2", "6-311G", "RHF"),         # spherical valence d of Ca
    ("ZnH2", "6-31G*", "RHF"),         # 6D 7F
    ("ZnH2", "def2-TZVPP", "RHF"),     # g functions on Zn
    ("HBr", "6-31G*", "RHF"),
    ("HBr", "def2-SVP", "RHF"),
    ("HeH+", "6-311G", "RHF"),       # He is defined in 6-311G (not in 6-311++G*)
    ("H2O", "6-31G*", "CCSD(T)"),
    ("H2O", "def2-SVP", "CCSD(T)"),
    ("NH3", "6-311++G**", "CCSD(T)"),
]

GRADIENTS = [
    ("H2O-distorted", "6-31+G**", "RHF"),
    ("HCl-stretched", "6-31G*", "RHF"),
    ("H2O-distorted", "def2-TZVP", "RHF"),
    ("H2O-distorted", "6-31G*", "CCSD"),
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
    path = str(rr.BASIS_DIR / FILES[basis])
    bas = {s: parse_gaussian.load(path, s) for s, _ in at}
    atoms_bohr = [[s, [c / rr.BOHR_IN_ANGSTROM for c in xyz]] for s, xyz in at]
    return gto.M(atom=atoms_bohr, basis=bas, unit="Bohr", cart=cartesian_d(basis), charge=charge, spin=0, verbose=0)


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


def rhf_6d7f(mol):
    """RHF with Cartesian d and spherical f functions: PySCF's Cartesian integrals (cart=True)
    transformed to spherical functions in the f shells only (PySCF's own Cartesian -> spherical
    map); returns (energy, number of functions)."""
    blocks = []
    for ib in range(mol.nbas):
        l = mol.bas_angular(ib)
        for _ in range(mol.bas_nctr(ib)):
            blocks.append(gto.cart2sph(l, normalized="sp") if l == 3 else np.eye((l + 1) * (l + 2) // 2))
    nc = sum(b.shape[0] for b in blocks)
    ns = sum(b.shape[1] for b in blocks)
    C = np.zeros((nc, ns))
    r = c = 0
    for b in blocks:
        C[r:r + b.shape[0], c:c + b.shape[1]] = b
        r += b.shape[0]
        c += b.shape[1]
    S = C.T @ mol.intor("int1e_ovlp") @ C
    H = C.T @ (mol.intor("int1e_kin") + mol.intor("int1e_nuc")) @ C
    eri = ao2mo.incore.full(mol.intor("int2e", aosym="s8"), C)
    mf = scf.RHF(mol)
    mf.get_ovlp = lambda *a: S
    mf.get_hcore = lambda *a: H
    mf._eri = eri
    mf.init_guess = "1e"
    mf.conv_tol = 1e-12
    mf.conv_tol_grad = 1e-8
    mf.max_cycle = 300
    mf.kernel(dm0=mf.init_guess_by_1e())
    if not mf.converged:
        raise SystemExit("RHF (6D 7F) did not converge")
    return mf.e_tot, ns


def compute():
    data = {"structures": {}, "energies": [], "gradients": []}
    for name, (title, charge, text) in STRUCTURES.items():
        data["structures"][name] = {"title": title, "charge": charge,
                                    "xyz": f"{len(text.strip().splitlines())}\n{title}\n{text}\n"}
    for name, basis, method in ENERGIES:
        mol = build(name, basis)
        has_f = any(mol.bas_angular(i) == 3 for i in range(mol.nbas))
        if mol.cart and has_f:
            e, nbf = rhf_6d7f(mol)
            assert method == "RHF"
        else:
            mf = rhf(mol)
            e, nbf = mf.e_tot, mol.nao
            if method == "CCSD(T)":
                c = cc.CCSD(mf)
                c.conv_tol = 1e-11
                c.conv_tol_normt = 1e-8
                c.kernel()
                e = c.e_tot + c.ccsd_t()
        data["energies"].append({"molecule": name, "basis": basis, "method": method, "energy": e,
                                 "basisFunctions": nbf, "cartesianD": cartesian_d(basis)})
        print(f"{name:14s} {basis:18s} {method:8s} {e:.10f}  ({nbf} functions)", flush=True)

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
        print(f"{name:14s} {basis:18s} {method:8s} gradient, |g|max {abs(g).max():.6f}", flush=True)

    # RHF Hessian, dipole moment and orbital values of H2O in 6-31G* (6D)
    mol = build("H2O", "6-31G*")
    mf = rhf(mol)
    h = mf.Hessian().kernel()  # [atom][atom][3][3]
    n = mol.natm
    H = np.zeros((3 * n, 3 * n))
    for a in range(n):
        for b in range(n):
            H[3 * a:3 * a + 3, 3 * b:3 * b + 3] = h[a, b]
    dip = mf.dip_moment(unit="AU", verbose=0)
    pts = np.array([[0.3, 0.2, 0.5], [1.0, -0.4, 1.2], [-0.7, 0.6, -0.3], [0.0, 0.0, 1.5], [1.43, 0.1, 1.0]])
    ao = mol.eval_gto("GTOval_cart", pts)
    homo = mol.nelectron // 2 - 1
    values = {"HOMO": (ao @ mf.mo_coeff[:, homo]).tolist(), "LUMO": (ao @ mf.mo_coeff[:, homo + 1]).tolist()}
    data["h2o_631gs"] = {"hessian": H.ravel().tolist(), "dipole": [float(v) for v in dip], "points": pts.tolist(),
                         "orbitals": values, "orbitalEnergies": mf.mo_energy[:homo + 2].tolist()}
    print("H2O 6-31G*: Hessian, dipole", dip, flush=True)
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
