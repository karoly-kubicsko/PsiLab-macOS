#!/usr/bin/env python3
"""Spin-orbital equations of ΨLab 15.0's dense spin-orbital engine (so_engine.cpp).

Usage:  python3 backend/tools/so_equations.py           (verify and write backend/include/macqc/so_diagrams.inc)
        python3 backend/tools/so_equations.py --check   (verify only)

Every correlated method of the engine -- coupled cluster, quadratic CI, truncated CI, Moller-Plesset
perturbation theory to fifth order and the perturbative triples -- is assembled from two families of
terms of <mu_m| H_N X_1 ... X_k |0> (H_N = F_N + V_N normal ordered, the X_j excitation operators
T_n = sum t(a1..an; i1..in) a1+ .. an+ in .. i1):

  C  connected: one Hamiltonian vertex (F or V) and k <= 4 cluster vertices, each joined to the
     vertex by cp particle and ch hole lines (the diagrams of backend/tools/cc_diagrams.py, up to
     rank 4, enumerated by the same code);
  D  disconnected (truncated CI): a vertex with external lines only, f_ai {a+ i} or
     1/4 <ab||ij> {a+ b+ j i}, times one cluster vertex with no line to it (target rank n + 1 or
     n + 2). They make <mu|H_N C|0> of linear CI complete together with the connected terms with
     k <= 1.

Unlike cc_diagrams.inc (packed antisymmetric tensors, diagonal Fock elements in denominators), the
terms here are written for dense antisymmetric spin-orbital tensors and use the full Fock matrix, so
that the residual of a method is a polynomial in the amplitudes and linear in (f, <pq||rs>): the
engine differentiates it automatically (analytic gradients and Hessians). For each term the output
gives

  * the Hamiltonian block, "F:<upper>|<lower>" or "V:<upper>|<lower>", upper = v^x o^(c-x)
    (external particles, then contracted holes), lower = v^y o^(c-y) (contracted particles, then
    external holes);
  * the evaluation chain: one einsum per cluster vertex (letters: lower case contracted lines,
    upper case external lines; amplitudes t[particles..., holes...] with the contracted lines
    first), giving z with the external particle letters A, B, C, D and the hole letters I, J, K, L
    in group order (vertex, then cluster after cluster);
  * the antisymmetriser as cosets: R[ABCD..IJKL..] += coef * sum_c sign_c * z relabelled by c;
  * coef = sign / (n_v prod cp! ch!) with n_v the number of permutations of identical clusters.

The signs are fitted to a brute-force evaluation in a determinant space (random Fock matrix with its
diagonal, random <pq||rs> and amplitudes): <mu| e^-T H_N e^T |0> for C and <mu| H_N C |0> for D; the
fit must be exact to rounding.
"""
import itertools
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import cc_diagrams as ccd  # noqa: E402

MAX_RANK = 4
ROOT = pathlib.Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "backend" / "include" / "macqc" / "so_diagrams.inc"
PART = "ABCD"
HOLE = "IJKL"


class Term:
    def __init__(self, kind, vertex, x, y, clusters):
        self.kind = kind            # 'C' or 'D'
        self.vertex = vertex        # 'F' or 'V'
        self.x, self.y = x, y
        self.c = 1 if vertex == 'F' else 2
        self.clusters = list(clusters)  # (n, cp, ch)
        self.target = x + sum(n - cp for n, cp, ch in clusters)
        self.sign = 0
        self.chain = []
        self.z_letters = ""
        self.cosets = []
        self.block = ""

    @property
    def symmetry(self):
        counts = {}
        for t in self.clusters:
            counts[t] = counts.get(t, 0) + 1
        s = 1
        for k in counts.values():
            s *= math.factorial(k)
        return s

    @property
    def magnitude(self):
        m = self.symmetry
        for n, cp, ch in self.clusters:
            m *= math.factorial(cp) * math.factorial(ch)
        return 1.0 / m

    def __repr__(self):
        ts = " ".join(f"T{n}({cp},{ch})" for n, cp, ch in self.clusters)
        return f"{self.kind} R{self.target} <- {self.vertex}[x={self.x},y={self.y}] {ts}"


def connected_terms():
    out = []
    for d in ccd.enumerate_diagrams(MAX_RANK):
        out.append(Term('C', d.vertex, d.x, d.y, d.clusters))
    return out


def disconnected_terms():
    out = []
    for vertex, c in (('F', 1), ('V', 2)):
        for n in range(1, MAX_RANK + 1):
            if n + c <= MAX_RANK:
                out.append(Term('D', vertex, c, 0, [(n, 0, 0)]))
    return out


def dense_order_cost(t, order, no=10, nv=80):
    op, oh = t.y, t.c - t.x
    ep, eh = t.x, t.c - t.y
    flops = peak = 0.0
    for (n, cp, ch) in order:
        rest = nv ** (op - cp) * no ** (oh - ch) * nv ** ep * no ** eh
        flops += rest * nv ** cp * no ** ch * nv ** (n - cp) * no ** (n - ch)
        op, oh, ep, eh = op - cp, oh - ch, ep + n - cp, eh + n - ch
        peak = max(peak, nv ** op * no ** oh * nv ** ep * no ** eh)
    return flops, peak


def choose_order(t):
    best = None
    for perm in sorted(set(itertools.permutations(t.clusters))):
        flops, peak = dense_order_cost(t, perm)
        key = (flops + 4.0 * peak, perm)
        if best is None or key < best:
            best = key
    t.clusters = list(best[1])


def build_chain(t):
    """Letters and einsum specs of the evaluation chain, and the cosets."""
    lower = iter("abcdefghmnopqrstuvwxyz")  # contracted lines (no clash with A-L)
    ext_p = iter(PART)
    ext_h = iter(HOLE)
    vx_p = [next(ext_p) for _ in range(t.x)]
    vx_h = [next(ext_h) for _ in range(t.c - t.y)]
    cont_p = [[next(lower) for _ in range(cp)] for (n, cp, ch) in t.clusters]
    cont_h = [[next(lower) for _ in range(ch)] for (n, cp, ch) in t.clusters]
    cl_p = [[next(ext_p) for _ in range(n - cp)] for (n, cp, ch) in t.clusters]
    cl_h = [[next(ext_h) for _ in range(n - ch)] for (n, cp, ch) in t.clusters]
    upper = vx_p + sum(cont_h, [])
    low = sum(cont_p, []) + vx_h
    t.block = f"{t.vertex}:{'v' * t.x + 'o' * (t.c - t.x)}|{'v' * t.y + 'o' * (t.c - t.y)}"
    cur = "".join(upper + low)
    t.h_letters = cur
    open_h = list(cur)  # letters of the running intermediate: open vertex lines, then externals
    ext_so_far = [c for c in cur if c.isupper()]
    vertex_open = [c for c in cur if c.islower()]
    t.chain = []
    for j, (n, cp, ch) in enumerate(t.clusters):
        tl = "".join(cont_p[j] + cl_p[j] + cont_h[j] + cl_h[j])
        vertex_open = [c for c in vertex_open if c not in cont_p[j] and c not in cont_h[j]]
        ext_so_far = ext_so_far + cl_p[j] + cl_h[j]
        out = "".join(vertex_open + ext_so_far)
        t.chain.append((n, f"{cur},{tl}->{out}"))
        cur = out
    # z letters in group order: particles then holes
    groups_p = [vx_p] + cl_p
    groups_h = [vx_h] + cl_h
    zp = sum(groups_p, [])
    zh = sum(groups_h, [])
    t.z_letters = "".join(zp + zh)
    if t.chain:
        last_n, spec = t.chain[-1]
        lhs = spec.split("->")[0]
        t.chain[-1] = (last_n, lhs + "->" + t.z_letters)
    # cosets of the antisymmetriser
    m = len(zp)
    t.cosets = []
    for pp, sp in group_cosets([len(g) for g in groups_p], m, PART):
        for hp, sh in group_cosets([len(g) for g in groups_h], m, HOLE):
            t.cosets.append((pp + hp, sp * sh))


def group_cosets(sizes, m, letters):
    """Assignments of the m target letters to consecutive groups (sorted within each group):
    (letters of the z axes, sign of the permutation)."""
    out = []
    sizes = [s for s in sizes]

    def rec(gi, avail, acc):
        if gi == len(sizes):
            perm = [letters.index(c) for c in acc]
            out.append(("".join(acc), ccd.perm_sign(perm)))
            return
        for comb in itertools.combinations(avail, sizes[gi]):
            rest = [a for a in avail if a not in comb]
            rec(gi + 1, rest, acc + list(comb))

    rec(0, list(letters[:m]), [])
    return out


def evaluate(t, f, V, amps, no, nv):
    """Dense value of the term (without coef and sign) at all index tuples: R[A.., I..]."""
    o, v = slice(0, no), slice(no, no + nv)
    blk = {'o': o, 'v': v}
    up, lo = t.block[2:].split("|")
    H = (f if t.vertex == 'F' else V)[tuple(blk[c] for c in up + lo)]
    z = H
    for n, spec in t.chain:
        z = np.einsum(spec, z, amps[n], optimize=True)
    target = PART[:t.target] + HOLE[:t.target]
    R = np.zeros(z.shape if t.target else ())
    if t.target == 0:
        return np.array(z)
    for zl, sgn in t.cosets:
        R = R + sgn * np.einsum(zl + "->" + target, z)
    return R


def dense_amplitudes(amplitudes, no, nv, max_rank):
    t = {}
    for n in range(1, max_rank + 1):
        dense = np.zeros((nv,) * n + (no,) * n)
        for (r, A, I), val in amplitudes.items():
            if r != n:
                continue
            Ap = [a - no for a in A]
            for pa in itertools.permutations(range(n)):
                for pi in itertools.permutations(range(n)):
                    dense[tuple(Ap[k] for k in pa) + tuple(I[k] for k in pi)] = ccd.perm_sign(pa) * ccd.perm_sign(pi) * val
        t[n] = dense
    return t


def brute_ci(dets, f, V, amplitudes, max_rank):
    """<mu| H_N C |0> (linear), C = sum of the amplitudes as one excitation operator."""
    no, N = dets.no, dets.N
    h = f - np.einsum("piqi->pq", V[:, :no, :, :no])
    ndet = len(dets.dets)
    Hm = np.zeros((ndet, ndet))
    for col, det in enumerate(dets.dets):
        for p in range(N):
            for q in range(N):
                if h[p, q] != 0.0:
                    r = dets.apply([(p, True), (q, False)], det)
                    if r:
                        Hm[dets.index[r[0]], col] += r[1] * h[p, q]
        occs = [i for i in range(N) if det >> i & 1]
        for r_, s in itertools.combinations(occs, 2):
            for p, q in itertools.combinations(range(N), 2):
                res = dets.apply([(p, True), (q, True), (s, False), (r_, False)], det)
                if res:
                    Hm[dets.index[res[0]], col] += res[1] * V[p, q, r_, s]
    Hm -= Hm[dets.ref, dets.ref] * np.eye(ndet)
    vec = np.zeros(ndet)
    for (n, A, I), amp in amplitudes.items():
        det, sign = dets.apply(dets.excitation_ops(A, I), (1 << no) - 1)
        vec[dets.index[det]] += sign * amp
    phi = Hm @ vec
    res = {0: phi[dets.ref]}
    for n in range(1, max_rank + 1):
        for A in itertools.combinations(range(no, N), n):
            for I in itertools.combinations(range(no), n):
                det, sign = dets.apply(dets.excitation_ops(A, I), (1 << no) - 1)
                res[(n, A, I)] = sign * phi[dets.index[det]]
    return res


def fit(terms, kind, no=6, nv=7, samples=3, seed=11, fixed=()):
    rng = np.random.default_rng(seed)
    dets = ccd.Determinants(no, nv)
    rows = {m: [] for m in range(MAX_RANK + 1)}
    rhs = {m: [] for m in range(MAX_RANK + 1)}
    for _ in range(samples):
        f, V, _t, amplitudes = ccd.random_problem(rng, no, nv, MAX_RANK)
        amps = dense_amplitudes(amplitudes, no, nv, MAX_RANK)
        if kind == 'C':
            brute = ccd.brute_force_residuals(dets, f, V, amplitudes, MAX_RANK)
        else:
            brute = brute_ci(dets, f, V, amplitudes, MAX_RANK)
        vals = {id(t): evaluate(t, f, V, amps, no, nv) * t.magnitude for t in terms + list(fixed)}
        for m in range(MAX_RANK + 1):
            ts = [t for t in terms if t.target == m]
            fx = [t for t in fixed if t.target == m]
            idx = [()] if m == 0 else [tuple(A) + tuple(I) for A, I in ccd.sorted_tuples(m, no, nv)]
            for k, ix in enumerate(idx):
                rows[m].append([float(vals[id(t)][ix]) for t in ts])
                key = 0 if m == 0 else (m, tuple(a + no for a in ix[:m]), tuple(ix[m:]))
                b = brute[key] - sum(t.sign * float(vals[id(t)][ix]) for t in fx)
                rhs[m].append(b)
    worst = 0.0
    for m in range(MAX_RANK + 1):
        ts = [t for t in terms if t.target == m]
        if not ts:
            continue
        M = np.array(rows[m])
        b = np.array(rhs[m])
        coef, *_ = np.linalg.lstsq(M, b, rcond=None)
        for t, c in zip(ts, coef):
            if abs(abs(c) - 1.0) > 1e-8:
                raise SystemExit(f"coefficient {c:.6f} is not +-1 for {t}")
            t.sign = 1 if c > 0 else -1
        err = np.max(np.abs(M @ np.array([t.sign for t in ts], dtype=float) - b))
        worst = max(worst, err)
        print(f"{kind} R{m}: {len(ts):3d} terms, {M.shape[0]:5d} equations, max error {err:.2e}")
    if worst > 1e-10:
        raise SystemExit("the terms do not reproduce the brute-force values")


def write_include(terms):
    lines = [
        "// Generated by backend/tools/so_equations.py -- do not edit.",
        "// Dense spin-orbital terms <mu_m| H X_1 .. X_k |0> verified against brute-force evaluations in a",
        "// determinant space. Fields: kind ('C' connected, 'D' disconnected), target rank, Hamiltonian block,",
        "// coefficient, number of clusters, clusters {rank, cp, ch} in evaluation order, letters of the",
        "// Hamiltonian block, the einsum of each step, the letters of z, and the cosets (z letters, sign).",
    ]
    for t in sorted(terms, key=lambda t: (t.kind, t.target, t.block, len(t.clusters), t.clusters)):
        cl = ", ".join(f"{{{n}, {cp}, {ch}}}" for n, cp, ch in t.clusters)
        steps = ", ".join(f'"{spec}"' for _, spec in t.chain)
        cos = " ".join(f"{zl}{'+' if s > 0 else '-'}" for zl, s in t.cosets)
        coef = t.sign * t.magnitude
        lines.append(
            f"{{'{t.kind}', {t.target}, \"{t.block}\", {coef!r}, {len(t.clusters)}, {{{cl}}}, \"{t.h_letters}\", "
            f"{{{steps}}}, \"{t.z_letters}\", \"{cos}\"}},")
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"wrote {OUTPUT.relative_to(ROOT)} ({len(terms)} terms)")


def main():
    conn = connected_terms()
    disc = disconnected_terms()
    for t in conn + disc:
        choose_order(t)
        build_chain(t)
    print(f"{len(conn)} connected and {len(disc)} disconnected terms up to rank {MAX_RANK}")
    fit(conn, 'C')
    # linear CI: the connected terms with at most one cluster are fixed, the disconnected ones fitted
    linear = [t for t in conn if len(t.clusters) == 1]
    fit(disc, 'D', fixed=linear)
    if "--check" not in sys.argv:
        write_include(conn + disc)


if __name__ == "__main__":
    main()
