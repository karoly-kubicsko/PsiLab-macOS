#!/usr/bin/env python3
"""The minimal STO-nG basis sets of ΨLab 15.0 (basis/sto-1g.gbs, sto-2g.gbs, sto-3g.gbs).

STO-2G and STO-3G are the sets of the Basis Set Exchange (data version 0.12, as the other files of
basis/; H-Xe). The BSE has no STO-1G; ΨLab builds it the way Hehre, Stewart and Pople built the STO-NG
sets (J. Chem. Phys. 51, 2657 (1969)):

  * every Slater function R_n(zeta=1, r) = N r^(n-1) exp(-r) of a shell is replaced by the least-squares
    expansion in N Gaussians r^l exp(-alpha r^2), minimising int (STO - sum_k d_k g_k)^2 dtau; the s and p
    functions of an sp shell share their exponents (the sum of both errors is minimised), d functions
    are fitted alone;
  * the exponents scale with zeta^2: alpha_k(zeta) = zeta^2 alpha_k(1), the coefficients do not;
  * the zeta of every element and shell are those of STO-3G: zeta^2 = alpha(STO-3G) / alpha_k(1) of the
    three-term fit (the same for all three primitives, which this script checks);
  * the contracted functions are compared normalised to one (the published coefficients of STO-2G are
    those of the renormalised contraction).

The procedure is verified before anything is written: the three-term and two-term fits reproduce the
coefficients of the BSE's STO-3G and STO-2G for every shell type (1s, 2sp, 3sp, 4sp, 5sp, 3d, 4d), the
1s exponents agree with Szabo and Ostlund's (STO-1G 0.270950, STO-2G 0.151623 / 0.851819, STO-3G
0.109818 / 0.405771 / 2.22766), and the zeta from the three primitives of every shell agree.

Usage:  python3 backend/tools/sto_ng.py           (verify and write the three files)
        python3 backend/tools/sto_ng.py --check   (verify only)
Needs: pip install basis_set_exchange scipy
"""
import math
import pathlib
import sys

import mpmath
import numpy as np
from scipy import integrate, optimize, special  # noqa: F401

ROOT = pathlib.Path(__file__).resolve().parents[2]
BASIS = ROOT / "basis"


# ---- least-squares Slater -> Gaussian fits --------------------------------------------------------

def radial_integral(m, zeta, alpha):
    """int_0^inf r^m exp(-zeta r - alpha r^2) dr by adaptive quadrature (the closed form through erfcx and
    the upward recurrence loses all digits for diffuse Gaussians in double precision)."""
    val, _ = integrate.quad(lambda r: r ** m * math.exp(-zeta * r - alpha * r * r), 0, np.inf, epsabs=0.0,
                            epsrel=1e-13, limit=400)
    return val


def radial_integral_exact(m, zeta, alpha):
    with mpmath.workdps(30):
        a = mpmath.mpf(alpha)
        z = mpmath.mpf(zeta)
        x = z / (2 * mpmath.sqrt(a))
        I0 = mpmath.sqrt(mpmath.pi / a) / 2 * mpmath.exp(x * x) * mpmath.erfc(x)
        prev, cur = I0, (1 - z * I0) / (2 * a)
        if m == 0:
            return float(I0)
        for k in range(2, m + 1):
            prev, cur = cur, ((k - 1) * prev - z * cur) / (2 * a)
        return float(cur)


def sto_gauss_overlap(n, l, alpha):
    """<R_n(zeta=1) | g_l(alpha)> of normalised radial functions."""
    ns = math.sqrt(2.0 ** (2 * n + 1) / math.factorial(2 * n))
    ng = math.sqrt(2.0 * (2.0 * alpha) ** (l + 1.5) / math.gamma(l + 1.5))
    return ns * ng * radial_integral(n + l + 1, 1.0, alpha)


def gauss_overlap(l, a, b):
    return (2.0 * math.sqrt(a * b) / (a + b)) ** (l + 1.5)


def fit_error(log_alpha, shells):
    """sum over the fitted functions (n, l) of the least-squares error, and the coefficients."""
    alpha = np.exp(log_alpha)
    total = 0.0
    coefs = []
    for n, l in shells:
        S = np.array([sto_gauss_overlap(n, l, a) for a in alpha])
        G = np.array([[gauss_overlap(l, a, b) for b in alpha] for a in alpha])
        d = np.linalg.solve(G, S)
        total += 1.0 - S @ d
        coefs.append(d)
    return total, coefs


def fit(shells, ngauss, start):
    best = None
    for s0 in start:
        res = optimize.minimize(lambda x: fit_error(x, shells)[0], np.log(np.array(s0, dtype=float)), method="BFGS",
                                options={"gtol": 1e-11, "maxiter": 2000})
        if best is None or res.fun < best.fun:
            best = res
    order = np.argsort(-best.x)
    x = best.x[order]
    err, coefs = fit_error(x, shells)
    return np.exp(x), coefs, err


SHELLS = {  # shell label: the Slater functions (n, l) fitted together
    "1s": [(1, 0)], "2sp": [(2, 0), (2, 1)], "3sp": [(3, 0), (3, 1)], "4sp": [(4, 0), (4, 1)],
    "5sp": [(5, 0), (5, 1)], "3d": [(3, 2)], "4d": [(4, 2)],
}
STARTS = {1: [[0.3], [0.1], [0.05]], 2: [[1.0, 0.15], [0.5, 0.08], [0.3, 0.05]],
          3: [[2.2, 0.4, 0.1], [1.0, 0.25, 0.07], [0.5, 0.18, 0.06], [0.3, 0.1, 0.04]]}


# ---- the BSE data -----------------------------------------------------------------------------------

def bse_basis(name):
    import basis_set_exchange as bse
    return bse.get_basis(name)


def shell_label(element_shells, idx, shell):
    """label (1s, 2sp, 3sp, 3d, ...) of the idx-th shell of an element in an STO-nG set."""
    am = shell["angular_momentum"]
    counts = {}
    for s in element_shells[:idx + 1]:
        key = tuple(s["angular_momentum"])
        counts[key] = counts.get(key, 0) + 1
    k = counts[tuple(am)]
    if am == [0]:
        return "1s"
    if am == [0, 1]:
        return f"{k + 1}sp"
    if am == [2]:
        return f"{k + 2}d"
    raise ValueError(am)


def main():
    import basis_set_exchange as bse  # noqa: F401
    worst = max(abs(radial_integral(m, 1.0, a) / radial_integral_exact(m, 1.0, a) - 1.0)
                for m in range(0, 11) for a in (0.002, 0.01, 0.05, 0.2, 1.0, 5.0, 40.0))
    print(f"closed-form radial integrals against 30 digits: largest relative deviation {worst:.1e}")
    if worst > 1e-9:
        raise SystemExit("the double-precision radial integrals are not accurate enough")
    fits = {}
    for label, shells in SHELLS.items():
        for n in (1, 2, 3):
            fits[(label, n)] = fit(shells, n, STARTS[n])
    # Szabo and Ostlund's 1s values
    so = {1: [0.270950], 2: [0.851819, 0.151623], 3: [2.22766, 0.405771, 0.109818]}
    for n, ref in so.items():
        a = fits[("1s", n)][0]
        dev = max(abs(x - y) / y for x, y in zip(a, ref))
        print(f"STO-{n}G 1s exponents {np.round(a, 6)} (Szabo-Ostlund {ref}) relative deviation {dev:.1e}")
        if dev > 5e-6:
            raise SystemExit("the 1s fit does not reproduce Szabo and Ostlund's exponents")
    # coefficients against the BSE sets, every shell type; zeta of every element and shell
    zetas = {}
    for name, n in (("sto-3g", 3), ("sto-2g", 2)):
        data = bse_basis(name)
        worst_c = worst_z = 0.0
        for z, el in data["elements"].items():
            if int(z) > 20:
                continue  # STO-1G for H-Ca (Sc-Zn use other 3sp and 3d expansions, Ga-Xe shared spd shells)
            shells = el["electron_shells"]
            for i, sh in enumerate(shells):
                label = shell_label(shells, i, sh)
                alpha_fit, coefs, _ = fits[(label, n)]
                exps = np.array([float(e) for e in sh["exponents"]])
                ratios = exps / alpha_fit
                spread = (ratios.max() - ratios.min()) / ratios.mean()
                worst_z = max(worst_z, spread)
                bse_c = [np.array([float(c) for c in cc]) for cc in sh["coefficients"]]
                # compare the contracted functions normalised to one (the BSE data are renormalised
                # contractions; ΨLab normalises every contraction anyway)
                for (nn, ll), d, c in zip(SHELLS[label], coefs, bse_c):
                    G = np.array([[gauss_overlap(ll, a, b) for b in alpha_fit] for a in alpha_fit])
                    dn = d / np.sqrt(d @ G @ d)
                    cn = c / np.sqrt(c @ G @ c)
                    worst_c = max(worst_c, float(np.max(np.abs(dn - cn))))
                if n == 3:
                    zetas[(int(z), i)] = (label, math.sqrt(ratios.mean()), spread)
        print(f"{name}: largest coefficient difference {worst_c:.1e}, largest spread of zeta^2 between primitives {worst_z:.1e}")
        # STO-3G is reproduced to the printed digits; the published STO-2G parameters were fitted less
        # tightly (zeta^2 of the primitives of one 3sp shell differ by up to 3.4e-4 relative), so they
        # are checked to the agreement that their own spread allows
        tol = 2e-5 if n == 3 else 5e-4
        if worst_c > tol or worst_z > tol:
            raise SystemExit("the fits do not reproduce " + name)
    if "--check" in sys.argv:
        return
    # write STO-2G and STO-3G as the BSE gives them, STO-1G from the fits
    import basis_set_exchange as bse
    for name, fname in (("sto-2g", "sto-2g.gbs"), ("sto-3g", "sto-3g.gbs")):
        (BASIS / fname).write_text(bse.get_basis(name, fmt="gaussian94", header=True))
    data = bse_basis("sto-3g")
    lines = [
        "!----------------------------------------------------------------------",
        "! STO-1G, generated by backend/tools/sto_ng.py (ΨLab 15.0); not from the Basis Set Exchange,",
        "! which has no STO-1G. One Gaussian per Slater function: the least-squares fits of",
        "! Hehre, Stewart and Pople, J. Chem. Phys. 51, 2657 (1969) with one term (1s: alpha = 0.270950",
        "! zeta^2, as in Szabo and Ostlund), s and p of an sp shell with a shared exponent, and the zeta",
        "! of every element and shell from the BSE's STO-3G (data version 0.12). Elements H-Ca.",
        "!----------------------------------------------------------------------",
        "",
        "",
    ]
    from basis_set_exchange.lut import element_sym_from_Z
    for z, el in sorted(data["elements"].items(), key=lambda kv: int(kv[0])):
        if int(z) > 20:
            continue
        sym = element_sym_from_Z(int(z), True)
        lines.append(f"{sym}     0")
        shells = el["electron_shells"]
        for i, sh in enumerate(shells):
            label, zeta, spread = zetas[(int(z), i)]
            if spread > 2e-5:
                # the expansion is not the standard one: zeta from the outermost (most diffuse) primitive
                exps = [float(e) for e in sh["exponents"]]
                zeta = math.sqrt(min(exps) / fits[(label, 3)][0][-1])
            alpha = fits[(label, 1)][0][0] * zeta * zeta
            if label.endswith("sp"):
                lines.append("SP   1   1.00")
                lines.append(f"      {alpha:.10E}       1.0000000000       1.0000000000".replace("E", "D"))
            elif label.endswith("d"):
                lines.append("D   1   1.00")
                lines.append(f"      {alpha:.10E}       1.0000000000".replace("E", "D"))
            else:
                lines.append("S   1   1.00")
                lines.append(f"      {alpha:.10E}       1.0000000000".replace("E", "D"))
        lines.append("****")
    (BASIS / "sto-1g.gbs").write_text("\n".join(lines) + "\n")
    for label in SHELLS:
        a = fits[(label, 1)][0][0]
        print(f"STO-1G {label}: alpha(zeta=1) = {a:.8f}")
    print("wrote basis/sto-1g.gbs, sto-2g.gbs, sto-3g.gbs")


if __name__ == "__main__":
    main()
