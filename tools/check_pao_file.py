# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/pao-basis.md and
# src/modules/two-center-overlap-kinetic.md on a real PAO file.
#
# The PAO file itself is not part of this repository; pass its path.
# Run:
#   uv run --with numpy --with scipy python tools/check_pao_file.py path/to/C6.0.pao
#
# What is checked (all with this book's conventions, independent of OpenMX code):
#   1. mesh: x = ln r uniform, extent relative to the cutoff radius
#   2. radial functions vanish beyond the cutoff radius in the table
#   3. orthonormality  int r^2 R_{l,z} R_{l,z'} dr = delta  for the leading multiplicities
#   4. S_ss(tau), T_ss(tau) for the first s function: k-space evaluation on the
#      default momentum grid (uniform k in [1e-6, 60], 901 points, trapezoid;
#      R~(k) from a uniform-r trapezoid on [r_min, r_c+0.5] with 900 intervals)
#      against a dense reference (direct real-space 2D integration).
import sys
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn


def read_pao(path):
    lines = open(path, errors="replace").read().splitlines()
    kv = {}
    for ln in lines:
        t = ln.split("#")[0].split()
        if len(t) >= 2 and t[0] in ("PAO.Lmax", "PAO.Mul", "grid.num.output",
                                     "radial.cutoff.pao", "AtomSpecies"):
            kv.setdefault(t[0], t[1])
    lmax, mul, n = int(kv["PAO.Lmax"]), int(kv["PAO.Mul"]), int(kv["grid.num.output"])
    rc = float(kv.get("radial.cutoff.pao", 5.0))
    tabs = {}
    for l in range(lmax + 1):
        tag = f"<pseudo.atomic.orbitals.L={l}"
        i0 = next(i for i, s in enumerate(lines) if s.strip().startswith(tag))
        data = np.array([[float(v) for v in lines[i0 + 1 + i].split()] for i in range(n)])
        tabs[l] = data
    x, r = tabs[0][:, 0], tabs[0][:, 1]
    R = {l: tabs[l][:, 2:2 + mul] for l in tabs}
    return dict(lmax=lmax, mul=mul, rc=rc, x=x, r=r, R=R)


def radial_fn(p, l, z):
    """R(r) with the book's out-of-table convention (l=0: constant below r_min; l>0: 0)."""
    r, y = p["r"], p["R"][l][:, z]
    sp = CubicSpline(r, y)

    def f(q):
        q = np.asarray(q, float)
        out = np.where(q > r[-1], 0.0, sp(np.clip(q, r[0], r[-1])))
        below = q < r[0]
        out = np.where(below, y[0] if l == 0 else 0.0, out)
        return out
    return f


def main(path):
    p = read_pao(path)
    x, r, rc = p["x"], p["r"], p["rc"]
    dx = np.diff(x)
    print(f"Lmax={p['lmax']} Mul={p['mul']} rc={rc} points={len(r)}")
    print(f"ln-mesh uniform: max|dx-mean|={np.max(abs(dx - dx.mean())):.2e}, r_min={r[0]:.3e}, r_max={r[-1]:.4f}")
    print(f"exp(x)==r: max rel err {np.max(abs(np.exp(x) - r) / r):.2e}")
    # 2. zero beyond cutoff
    beyond = r > rc
    mx = max(np.max(abs(p["R"][l][beyond, :])) for l in p["R"])
    print(f"max |R| for r>rc in table: {mx:.2e}")
    # 3. orthonormality on the log mesh (int r^3 R R dx), first 4 multiplicities
    w = np.full_like(x, dx.mean()); w[0] *= 0.5; w[-1] *= 0.5
    for l in range(p["lmax"] + 1):
        A = p["R"][l][:, :4]
        O = (A * (r ** 3 * w)[:, None]).T @ A
        print(f"l={l}: max|O-I| (4x4) = {np.max(abs(O - np.eye(4))):.2e}")
    # 4. two-center s-s
    f = radial_fn(p, 0, 0)
    NR, NK, kmax, kmin = 900, 900, 60.0, 1e-6
    rmin, rmax = r[0], rc + 0.5
    rg = np.linspace(rmin, rmax, NR + 1); hr = rg[1] - rg[0]
    wr = np.full(NR + 1, hr); wr[0] *= 0.5; wr[-1] *= 0.5
    ktab = np.arange(NK) * kmax / NK
    Rt_tab = np.array([np.sum(wr * rg ** 2 * spherical_jn(0, k * rg) * f(rg)) for k in ktab])
    Rt = CubicSpline(ktab, Rt_tab)
    kg = np.linspace(kmin, kmax, NR + 1); hk = kg[1] - kg[0]
    wk = np.full(NR + 1, hk); wk[0] *= 0.5; wk[-1] *= 0.5
    Rk = np.where(kg > ktab[-1], 0.0, Rt(kg))
    # Real-space reference: S(tau) = 2pi int r^2 dr int dmu f(r) f(|r - tau|) Y00^2 ...
    # For s-s with real Y00 = 1/sqrt(4pi): S = (1/2) int r^2 dr f(r) int_{-1}^{1} dmu f(sqrt(r^2+t^2-2 r t mu))
    from numpy.polynomial.legendre import leggauss
    mu, wmu = leggauss(200)
    rr = np.linspace(0, rc, 6001); wrr = np.full_like(rr, rr[1]); wrr[0] *= .5; wrr[-1] *= .5
    for tau in (0.0, 1.0, 2.5, 4.0, 8.0, 11.9, 12.0):
        # 8 * G(00;00;00) * Y00 * I = 8 * (1/sqrt(4pi)) * (1/sqrt(4pi)) * I = (2/pi) I
        Sk = (2 / np.pi) * np.sum(wk * kg ** 2 * spherical_jn(0, kg * max(tau, 1e-10)) * Rk * Rk)
        Tk = (1 / np.pi) * np.sum(wk * kg ** 4 * spherical_jn(0, kg * max(tau, 1e-10)) * Rk * Rk)
        d = np.sqrt(np.maximum(rr[:, None] ** 2 + tau ** 2 - 2 * rr[:, None] * tau * mu[None, :], 0))
        inner = (f(d) * wmu[None, :]).sum(1)
        Sref = 0.5 * np.sum(wrr * rr ** 2 * f(rr) * inner)
        print(f"tau={tau:5.2f}  S_k={Sk:.10f}  S_ref={Sref:.10f}  diff={Sk - Sref:.1e}  T_k={Tk:.8f}")


if __name__ == "__main__":
    main(sys.argv[1])
