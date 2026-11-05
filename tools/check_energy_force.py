# SPDX-License-Identifier: MIT
"""Numerical checks for the energy/force/stress chapter (src/modules/energy-force.md).

1. Pulay identity: d/dx sum_n f_n eps_n = Tr[rho dH/dx] - Tr[E dS/dx]
   for the generalized eigenproblem H c = eps S c with fixed occupations.
2. Hartree stress: dE_H/de_{gh} = delta_{gh} E_H + Int (d dn/de) dV_H
   + 4 pi Omega sum_G |dn(G)|^2 G_g G_h / G^4,
   for rigid, neutral Gaussian charge clouds whose centres follow the strain.
3. Screened pair integral Int n_I V_J on a cylindrical grid (uniform z planes,
   Gauss-Legendre in the cylindrical radius) against the analytic erf result,
   and the cutoff rescaling factor.

Run: uv run --with numpy --with scipy python tools/check_energy_force.py
"""
import numpy as np
from scipy.linalg import eigh
from scipy.special import erf

rng = np.random.default_rng(1)


# ---------------------------------------------------------------- 1. Pulay
def check_pulay():
    n = 8
    a = rng.normal(size=(n, n))
    h0 = a + a.T
    h1 = rng.normal(size=(n, n))
    h1 = h1 + h1.T
    b = rng.normal(size=(n, n)) * 0.2
    s0 = np.eye(n) + 0.5 * (b + b.T)
    s1 = rng.normal(size=(n, n)) * 0.1
    s1 = s1 + s1.T
    f = np.array([2, 2, 2, 1.3, 0.4, 0, 0, 0], float)  # fixed occupations

    def band(x):
        e, c = eigh(h0 + x * h1, s0 + x * s1)
        return (f * e).sum(), e, c

    x0, dx = 0.1, 1e-5
    ep = band(x0 + dx)[0]
    em = band(x0 - dx)[0]
    fd = (ep - em) / (2 * dx)
    _, e, c = band(x0)
    rho = (c * f) @ c.T
    edm = (c * (f * e)) @ c.T
    an = np.sum(rho * h1) - np.sum(edm * s1)
    print(f"[Pulay] FD={fd:.10f} analytic={an:.10f} diff={abs(fd-an):.2e}")
    assert abs(fd - an) < 1e-7


# ------------------------------------------------------ 2. Hartree stress
def gauss_ft(q, a):
    """FT of normalised Gaussian (a/pi)^{3/2} exp(-a r^2)."""
    return np.exp(-q**2 / (4 * a))


def check_hartree_stress():
    A0 = np.array([[6.0, 0.3, 0.0], [0.5, 5.5, 0.2], [0.1, -0.4, 7.0]])  # rows = a_i
    tfrac = np.array([[0.1, 0.2, 0.3], [0.6, 0.55, 0.7]])
    # each site: +Z narrow Gaussian - Z wide Gaussian (neutral)
    Zs, a1, a2 = [1.0, 2.0], [1.5, 2.5], [0.6, 0.8]
    nmax = 9
    m = np.array([[i, j, k] for i in range(-nmax, nmax + 1)
                  for j in range(-nmax, nmax + 1)
                  for k in range(-nmax, nmax + 1) if (i, j, k) != (0, 0, 0)], float)

    def setup(eps):
        A = A0 @ (np.eye(3) + eps).T          # a_i' = (1+eps) a_i
        Om = abs(np.linalg.det(A))
        B = 2 * np.pi * np.linalg.inv(A).T    # rows = b_i
        G = m @ B
        g = np.linalg.norm(G, axis=1)
        T = tfrac @ A                        # Cartesian centres (follow strain)
        return A, Om, G, g, T

    def dn_G(Om, G, g, T, dfdq=False):
        s = np.zeros(len(g), complex)
        for Z, p, q, t in zip(Zs, a1, a2, T):
            ph = np.exp(-1j * G @ t)
            if dfdq:
                f = Z * (-g / (2 * p) * gauss_ft(g, p) + g / (2 * q) * gauss_ft(g, q))
            else:
                f = Z * (gauss_ft(g, p) - gauss_ft(g, q))
            s += f * ph
        return s / Om

    def energy(eps):
        A, Om, G, g, T = setup(eps)
        d = dn_G(Om, G, g, T)
        return 0.5 * Om * np.sum(4 * np.pi * np.abs(d) ** 2 / g**2)

    E0 = energy(np.zeros((3, 3)))
    A, Om, G, g, T = setup(np.zeros((3, 3)))
    d = dn_G(Om, G, g, T)
    dV = 4 * np.pi * d / g**2
    dfq = dn_G(Om, G, g, T, dfdq=True)     # derivative of radial FT wrt |G|, /Omega
    h = 1e-6
    worst = 0.0
    for gi in range(3):
        for hi in range(gi, 3):
            e = np.zeros((3, 3)); e[gi, hi] = e[hi, gi] = 0.5 * h if gi != hi else h
            # off-diagonal: symmetric strain e_gh = e_hg = h/2 gives the symmetric component
            fd = (energy(e) - energy(-e)) / (2 * h)
            # analytic pieces (pulled-back density, see chapter)
            # d dn(G)/de_gh = -delta dn(G) + dfq * d|G|/de_gh, d|G|/de = -G_g G_h/|G|
            ddn = -(gi == hi) * d + dfq * (-G[:, gi] * G[:, hi] / g)
            term_vol = (gi == hi) * E0
            term_dn = Om * np.real(np.sum(np.conj(ddn) * dV))
            term_G = 4 * np.pi * Om * np.sum(np.abs(d) ** 2 * G[:, gi] * G[:, hi] / g**4)
            an = term_vol + term_dn + term_G
            worst = max(worst, abs(fd - an))
            print(f"[Hartree stress] ({gi},{hi}) FD={fd:+.8f} analytic={an:+.8f}")
    print(f"[Hartree stress] max diff={worst:.2e}")
    assert worst < 1e-6


# ------------------------------------------- 3. cylindrical pair integral
def check_pair_integral():
    aI, aJ, ZI, ZJ, rc = 1.2, 0.9, 4.0, 6.0, 4.5   # rc: support used for n_I
    def nI(r):
        return ZI * (aI / np.pi) ** 1.5 * np.exp(-aI * r**2)
    def VJ(r):
        r = np.maximum(r, 1e-12)
        return ZJ * erf(np.sqrt(aJ) * r) / r
    def exact(R):
        mu = aI * aJ / (aI + aJ)
        return ZI * ZJ * erf(np.sqrt(mu) * R) / R
    # grid: z uniform in [-rc, rc] with spacing ~ pi/sqrt(9600); cylindrical radius by GL(150)
    dx0 = np.pi / np.sqrt(16 * 600.0)
    nd = 2 * int(rc / dx0) + 1
    z = np.linspace(-rc, rc, nd); dz = z[1] - z[0]
    xg, wg = np.polynomial.legendre.leggauss(150)
    rho = 0.5 * rc * (xg + 1); wr = 0.5 * rc * wg * 2 * np.pi * rho
    Zg, Rg = np.meshgrid(z, rho, indexing="ij")
    W = np.outer(np.full(nd, dz), wr)
    dens = nI(np.sqrt(Zg**2 + Rg**2))
    def quad(R):
        return np.sum(W * dens * VJ(np.sqrt((Zg - R) ** 2 + Rg**2)))
    for R in [0.0, 1.0, 2.5, 4.0]:
        q, ex = quad(R), exact(R) if R > 0 else ZI * ZJ * 2 * np.sqrt(aI * aJ / (aI + aJ) / np.pi)
        print(f"[pair integral] R={R:.2f} quad={q:.10f} exact={ex:.10f} rel={abs(q-ex)/abs(ex):.1e}")
        assert abs(q - ex) / abs(ex) < 1e-6
    s = ZI * ZJ / 9.0 / quad(9.0)
    print(f"[pair integral] cutoff rescale factor at R=9.0: {s:.12f}")


if __name__ == "__main__":
    check_pulay()
    check_hartree_stress()
    check_pair_integral()
    print("all checks passed")
