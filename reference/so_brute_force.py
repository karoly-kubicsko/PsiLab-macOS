#!/usr/bin/env python3
"""Brute-force reference energies of the correlated methods of ΨLab 15.0, in a determinant space.

Independent of the equations of ΨLab's spin-orbital engine (backend/tools/so_equations.py): every
method is evaluated from its definition with sparse matrices of the Hamiltonian and of the cluster
(excitation) operators over all determinants with the electrons' numbers of alpha and beta spin, in
the spin orbitals of a PySCF SCF (RHF or UHF) from ΨLab's basis-set files:

  MPn       Rayleigh-Schroedinger perturbation theory with H0 = sum_p eps_p n_p (canonical orbitals),
            intermediate normalisation, to fifth order;
  CI-n      the lowest eigenvalue of H in the determinants up to n-fold excited from the reference
            (the root with the largest weight of the reference);
  QCI-n     the projections of H + [H, T] + selected double commutators (the connected quadratic
            terms, (H T_a T_b)_C = [[H, T_a], T_b], (H T_2^2/2)_C = [[H, T_2], T_2]/2): singles T1 T2,
            doubles T2^2/2, triples T2 T3, quadruples T2 T4 and T2^2/2; E = <0|H (T1 + T2)|0>;
  CC-n      <mu| exp(-T) H exp(T) |0> = 0;
  (T)       t3 = -<T|[H, T2]|0>/D3, E_T = sum t3 (<T|[H, T2]|0> + f <T|T1 H|0>), f = 1 (CCSD(T)) or 2
            (QCISD(T), CISD(T)).

Used for the reference values of macqc_so_tests (backend/tests/data/so_reference.json):
  python3 reference/so_brute_force.py           # -> backend/tests/data/so_reference.json
"""
import itertools
import json
import math
import pathlib
import sys

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
from pyscf import ao2mo, gto, scf
from pyscf.gto.basis import parse_gaussian

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "backend" / "tests" / "data" / "so_reference.json"
BOHR = 0.529177210544


def popcount(x):
    return bin(x).count("1")


class Space:
    """Determinants over 2n spin orbitals (alpha 0..n-1, beta n..2n-1) with na alpha and nb beta electrons."""

    def __init__(self, n, na, nb):
        self.n, self.na, self.nb = n, na, nb
        alpha = [sum(1 << i for i in c) for c in itertools.combinations(range(n), na)]
        beta = [sum(1 << (n + i) for i in c) for c in itertools.combinations(range(n), nb)]
        self.dets = [a | b for a in alpha for b in beta]
        self.index = {d: k for k, d in enumerate(self.dets)}
        self.ref_det = sum(1 << i for i in range(na)) | sum(1 << (n + i) for i in range(nb))
        self.ref = self.index[self.ref_det]
        self.occ = [i for i in range(2 * n) if self.ref_det >> i & 1]
        self.vir = [i for i in range(2 * n) if not self.ref_det >> i & 1]

    @staticmethod
    def apply(ops, det):
        sign = 1
        for orb, creation in reversed(ops):
            bit = 1 << orb
            occupied = bool(det & bit)
            if creation == occupied:
                return None
            if popcount(det & (bit - 1)) % 2:
                sign = -sign
            det ^= bit
        return det, sign

    def level(self, det):
        return popcount(det & ~self.ref_det)

    def excitations(self, rank):
        """unique excitations (A, I) of the given rank that conserve the spin projection"""
        n = self.n
        out = []
        for A in itertools.combinations(self.vir, rank):
            for I in itertools.combinations(self.occ, rank):
                if sum(a >= n for a in A) == sum(i >= n for i in I):
                    out.append((A, I))
        return out

    def ket(self, A, I):
        """index and sign of a+_A a_I |0>"""
        ops = [(a, True) for a in A] + [(i, False) for i in reversed(I)]
        det, sign = self.apply(ops, self.ref_det)
        return self.index[det], sign


def hamiltonian_matrix(space, h, V):
    """sparse <D|H|D'> for spin-orbital h and antisymmetrised V[p,q,r,s] = <pq||rs>"""
    N = 2 * space.n
    rows, cols, vals = [], [], []
    for col, det in enumerate(space.dets):
        occ = [i for i in range(N) if det >> i & 1]
        for q in occ:
            for p in range(N):
                if h[p, q] == 0.0:
                    continue
                r = space.apply([(p, True), (q, False)], det)
                if r and r[0] in space.index:
                    rows.append(space.index[r[0]]); cols.append(col); vals.append(r[1] * h[p, q])
        for r_, s in itertools.combinations(occ, 2):
            for p, q in itertools.combinations(range(N), 2):
                v = V[p, q, r_, s]
                if v == 0.0:
                    continue
                res = space.apply([(p, True), (q, True), (s, False), (r_, False)], det)
                if res and res[0] in space.index:
                    rows.append(space.index[res[0]]); cols.append(col); vals.append(res[1] * v)
    D = len(space.dets)
    return sps.csr_matrix((vals, (rows, cols)), shape=(D, D))


def cluster_matrix(space, amps):
    """sparse matrix of T = sum_(A,I) t a+_A a_I over the determinant space"""
    rows, cols, vals = [], [], []
    for (A, I), t in amps.items():
        if t == 0.0:
            continue
        ops = [(a, True) for a in A] + [(i, False) for i in reversed(I)]
        for col, det in enumerate(space.dets):
            r = space.apply(ops, det)
            if r and r[0] in space.index:
                rows.append(space.index[r[0]]); cols.append(col); vals.append(r[1] * t)
    D = len(space.dets)
    return sps.csr_matrix((vals, (rows, cols)), shape=(D, D))


class Problem:
    def __init__(self, atoms, basis_file, charge, mult, ref, cart=False):
        bas = {s: parse_gaussian.load(str(ROOT / "basis" / basis_file), s) for s, _ in atoms}
        mol = gto.M(atom=atoms, basis=bas, unit="Bohr", charge=charge, spin=mult - 1, cart=cart, verbose=0)
        mf = (scf.UHF if ref == "UHF" else scf.RHF)(mol)
        mf.conv_tol = 1e-13
        mf.conv_tol_grad = 1e-10
        mf.max_cycle = 400
        mf.kernel()
        self.mol, self.mf = mol, mf
        self.e_scf = mf.e_tot
        n = mol.nao
        if ref == "UHF":
            Ca, Cb = mf.mo_coeff
            ea, eb = mf.mo_energy
        else:
            Ca = Cb = mf.mo_coeff
            ea = eb = mf.mo_energy
        na, nb = mol.nelec
        self.space = Space(n, na, nb)
        N = 2 * n
        hcore = mf.get_hcore()
        h = np.zeros((N, N))
        h[:n, :n] = Ca.T @ hcore @ Ca
        h[n:, n:] = Cb.T @ hcore @ Cb
        eri = mol.intor("int2e")
        def mo(C1, C2, C3, C4):
            return np.einsum("pqrs,pi,qj,rk,sl->ijkl", eri, C1, C2, C3, C4, optimize=True)
        g = np.zeros((N, N, N, N))  # chemists' (pq|rs) over spin orbitals
        g[:n, :n, :n, :n] = mo(Ca, Ca, Ca, Ca)
        g[n:, n:, n:, n:] = mo(Cb, Cb, Cb, Cb)
        g[:n, :n, n:, n:] = mo(Ca, Ca, Cb, Cb)
        g[n:, n:, :n, :n] = g[:n, :n, n:, n:].transpose(2, 3, 0, 1)
        phys = g.transpose(0, 2, 1, 3)  # <pq|rs> = (pr|qs)
        self.V = phys - phys.transpose(0, 1, 3, 2)
        self.h = h
        self.eps = np.concatenate([ea, eb])
        self.H = hamiltonian_matrix(self.space, h, self.V)
        e_ref = self.H[self.space.ref, self.space.ref]
        self.e_nuc = mol.energy_nuc()
        self.e_ref = e_ref + self.e_nuc
        D = len(self.space.dets)
        self.Hbar = self.H - e_ref * sps.identity(D, format="csr")
        self.ref_vec = np.zeros(D)
        self.ref_vec[self.space.ref] = 1.0
        # diagonal of H0: sum of the orbital energies of the occupied spin orbitals of each determinant
        self.h0 = np.array([sum(self.eps[i] for i in range(N) if d >> i & 1) for d in self.space.dets])

    # -- perturbation theory ----------------------------------------------------------------------
    def mpn(self, nmax=5):
        sp = self.space
        D = len(sp.dets)
        e0 = self.h0[sp.ref]
        V = self.H - sps.diags(self.h0)
        denom = e0 - self.h0
        denom[sp.ref] = 1.0
        psi = [self.ref_vec.copy()]
        E = [e0, V[sp.ref] @ self.ref_vec]
        E[1] = float(E[1][0]) if hasattr(E[1], "__len__") else float(E[1])
        energies = {}
        for k in range(1, nmax):
            rhs = V @ psi[k - 1]
            for j in range(1, k + 1):
                rhs = rhs - E[j] * psi[k - j]
            new = rhs / denom
            new[sp.ref] = 0.0
            psi.append(new)
            E.append(float((V @ psi[k])[sp.ref]))
        total = self.e_nuc + e0 + E[1]
        for k in range(2, nmax + 1):
            total += E[k]
            energies[f"MP{k}"] = total
        return energies

    # -- truncated CI -------------------------------------------------------------------------------
    def ci(self, rank):
        sp = self.space
        keep = [k for k, d in enumerate(sp.dets) if sp.level(d) <= rank]
        Hs = self.H[keep][:, keep].toarray()
        w, v = np.linalg.eigh(Hs)
        r = keep.index(sp.ref)
        k = int(np.argmax(np.abs(v[r, :])))
        c = v[:, k] / v[r, k]
        full = np.zeros(len(sp.dets))
        full[keep] = c
        return w[k] + self.e_nuc, full

    # -- amplitude equations (QCI, CC) ----------------------------------------------------------------
    def amplitudes_to_matrices(self, x, layout):
        T = {}
        for rank, exc in layout.items():
            amps = {e: x[rank][k] for k, e in enumerate(exc)}
            T[rank] = cluster_matrix(self.space, amps)
        return T

    def project(self, vec, layout):
        out = {}
        for rank, exc in layout.items():
            r = np.zeros(len(exc))
            for k, (A, I) in enumerate(exc):
                idx, sign = self.space.ket(A, I)
                r[k] = sign * vec[idx]
            out[rank] = r
        return out

    def denominators(self, layout):
        out = {}
        for rank, exc in layout.items():
            out[rank] = np.array([sum(self.eps[a] for a in A) - sum(self.eps[i] for i in I) for A, I in exc])
        return out

    def solve(self, residual, layout, tol=1e-11, maxit=400):
        x = {r: np.zeros(len(e)) for r, e in layout.items()}
        Dn = self.denominators(layout)
        hist_x, hist_e = [], []
        for it in range(maxit):
            R, E = residual(x)
            step = {r: -R[r] / Dn[r] for r in x}
            err = max(np.max(np.abs(s)) for s in step.values()) if step else 0.0
            if err < tol:
                return x, E
            xn = {r: x[r] + step[r] for r in x}
            flat = np.concatenate([xn[r] for r in sorted(x)])
            ef = np.concatenate([step[r] for r in sorted(x)])
            hist_x.append(flat); hist_e.append(ef)
            hist_x, hist_e = hist_x[-8:], hist_e[-8:]
            m = len(hist_x)
            if m >= 2:
                B = -np.ones((m + 1, m + 1)); B[m, m] = 0.0
                for i in range(m):
                    for j in range(m):
                        B[i, j] = hist_e[i] @ hist_e[j]
                rhs = np.zeros(m + 1); rhs[m] = -1.0
                try:
                    c = np.linalg.solve(B, rhs)
                    flat = sum(c[i] * hist_x[i] for i in range(m))
                except np.linalg.LinAlgError:
                    pass
            pos = 0
            for r in sorted(x):
                x[r] = flat[pos:pos + len(x[r])]
                pos += len(x[r])
        raise RuntimeError("amplitude equations did not converge")

    def layout(self, nmax):
        return {r: self.space.excitations(r) for r in range(1, nmax + 1)}

    def qci(self, nmax):
        layout = self.layout(nmax)
        Hb, ref = self.Hbar, self.ref_vec

        def comm(A, vec):  # [A, T] acting on vec needs T; helper below
            raise NotImplementedError

        def residual(x):
            T = self.amplitudes_to_matrices(x, layout)
            h0 = Hb @ ref
            lin = {r: Hb @ (T[r] @ ref) - T[r] @ h0 for r in T}  # [H, T_r]|0>
            def dc(a, b):  # [[H, T_a], T_b]|0>
                Ta, Tb = T[a], T[b]
                HTa = lambda v: Hb @ (Ta @ v) - Ta @ (Hb @ v)  # noqa: E731
                return HTa(Tb @ ref) - Tb @ HTa(ref)
            R = {}
            for m in range(1, nmax + 1):
                vec = sum(lin.values())
                if m == 2:
                    vec = vec + h0 + 0.5 * dc(2, 2)
                if m == 1:
                    vec = vec + h0 + dc(1, 2)
                if m >= 3 and m in T:
                    vec = vec + dc(2, m)
                if m == 4:
                    vec = vec + 0.5 * dc(2, 2)
                R[m] = self.project(vec, {m: layout[m]})[m]
            E = float(ref @ (Hb @ (T[1] @ ref + T[2] @ ref)))
            return R, E

        x, E = self.solve(residual, layout)
        return E + self.e_ref, x, layout

    def cc(self, nmax):
        layout = self.layout(nmax)
        Hb, ref = self.Hbar, self.ref_vec

        def expm(T, v, s):
            out, term = v.copy(), v.copy()
            for k in range(1, 2 * nmax + 4):
                term = s * (T @ term) / k
                out = out + term
            return out

        def residual(x):
            T = self.amplitudes_to_matrices(x, layout)
            Tt = sum(T.values())
            vec = expm(Tt, Hb @ expm(Tt, ref, 1.0), -1.0)
            R = {m: self.project(vec, {m: layout[m]})[m] for m in layout}
            return R, float(vec[self.space.ref])

        x, E = self.solve(residual, layout)
        return E + self.e_ref, x, layout

    def triples(self, x, layout, f):
        """(T) correction with singles x[1] and doubles x[2] (any normalisation convention: intermediate)"""
        sp = self.space
        lay = {1: layout[1], 2: layout[2]}
        T = self.amplitudes_to_matrices({1: x[1], 2: x[2]}, lay)
        Hb, ref = self.Hbar, self.ref_vec
        X = Hb @ (T[2] @ ref) - T[2] @ (Hb @ ref)
        Y = T[1] @ (Hb @ ref)
        exc3 = sp.excitations(3)
        P = self.project(X, {3: exc3})[3]
        Q = self.project(Y, {3: exc3})[3]
        D3 = self.denominators({3: exc3})[3]
        t3 = -P / D3
        return float(t3 @ P), float(t3 @ Q) * f


def run_case(name, atoms, basis_file, charge, mult, ref, methods, cart=False):
    p = Problem(atoms, basis_file, charge, mult, ref, cart)
    out = {"name": name, "basis_file": basis_file, "charge": charge, "multiplicity": mult, "reference": ref,
           "atoms": [[s, list(xyz)] for s, xyz in atoms], "scf": p.e_scf, "determinants": len(p.space.dets),
           "energies": {}}
    if any(m.startswith("MP") for m in methods):
        out["energies"].update({k: v for k, v in p.mpn(5).items() if k in methods})
    for rank, mname in ((2, "CISD"), (3, "CISDT"), (4, "CISDTQ")):
        if mname in methods or (mname == "CISD" and "CISD(T)" in methods):
            e, c = p.ci(rank)
            out["energies"][mname] = e
            if "CISD(T)" in methods and rank == 2:
                lay = p.layout(2)
                amps = p.project(c, lay)
                e4, e5 = p.triples(amps, lay, 2.0)
                out["energies"]["CISD(T)"] = e + e4 + e5
    for rank, mname in ((2, "QCISD"), (3, "QCISDT"), (4, "QCISDTQ")):
        if mname in methods or (rank == 2 and "QCISD(T)" in methods):
            e, x, lay = p.qci(rank)
            out["energies"][mname] = e
            if rank == 2 and "QCISD(T)" in methods:
                e4, e5 = p.triples(x, lay, 2.0)
                out["energies"]["QCISD(T)"] = e + e4 + e5
    for rank, mname in ((2, "CCSD"), (3, "CCSDT"), (4, "CCSDTQ")):
        if mname in methods or (rank == 2 and "CCSD(T)" in methods):
            e, x, lay = p.cc(rank)
            out["energies"][mname] = e
            if rank == 2 and "CCSD(T)" in methods:
                e4, e5 = p.triples(x, lay, 1.0)
                out["energies"]["CCSD(T)"] = e + e4 + e5
    for k, v in out["energies"].items():
        print(f"  {name:12s} {ref} {k:9s} {v:.12f}", flush=True)
    return out


ALL = ["MP2", "MP3", "MP4", "MP5", "CISD", "CISD(T)", "CISDT", "CISDTQ", "QCISD", "QCISD(T)", "QCISDT", "QCISDTQ",
       "CCSD", "CCSD(T)", "CCSDT", "CCSDTQ"]


def main():
    R = 1.0 / BOHR
    cases = [
        # closed shell: H2O stretched (strong correlation), STO-3G
        ("H2O", [["O", (0, 0, 0)], ["H", (1.4 * R * 0.7570, 0, 1.4 * R * 0.5859)], ["H", (-1.4 * R * 0.7570, 0, 1.4 * R * 0.5859)]],
         "sto-3g.gbs", 0, 1, "RHF", ALL),
        # open shells: OH radical (doublet) and CH2 (triplet), UHF, STO-3G
        ("OH", [["O", (0, 0, 0)], ["H", (0, 0, 0.98 * R)]], "sto-3g.gbs", 0, 2, "UHF", ALL),
        ("CH2", [["C", (0, 0, 0)], ["H", (0, 0.9906 * R, 0.4537 * R)], ["H", (0, -0.9906 * R, 0.4537 * R)]],
         "sto-3g.gbs", 0, 3, "UHF", ALL),
        # more virtual orbitals: LiH 6-31G (Cartesian d not involved)
        ("LiH", [["Li", (0, 0, 0)], ["H", (0, 0, 1.75 * R)]], "6-31g.gbs", 0, 1, "RHF", ALL),
    ]
    out = []
    for c in cases:
        print(c[0], flush=True)
        out.append(run_case(*c))
    if "--check" not in sys.argv:
        OUT.write_text(json.dumps({"cases": out, "bohr_in_angstrom": BOHR}, indent=1) + "\n")
        print("wrote", OUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
