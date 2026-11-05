# SPDX-License-Identifier: MIT
"""Numerical checks for the dispersion chapter (src/modules/exx-vdw.md).

1. Unit conversion of the DFT-D2 parameters of Grimme, JCC 27, 1787 (2006):
   C6 in J nm^6 mol^-1 -> Hartree bohr^6 (and Rydberg bohr^6), R0 in Angstrom -> bohr.
2. DFT-D2 periodic pair sum: analytic gradient and strain derivative (stress)
   against central finite differences, including per-atom periodicity flags.
3. Image range: ceil(rcut/|a_k|) versus ceil(rcut/h_k) (h_k = interplanar
   spacing) in an oblique cell; the former can miss pairs inside rcut.
4. DFT-D3-type energy with a synthetic (made-up) reference C6 table:
   coordination numbers, Gaussian CN interpolation of C6 (with the max-shift of
   the exponents), zero and Becke-Johnson damping; analytic gradient including
   the dC6/dCN chain rule, and stress, against finite differences.

No parameter table of the original DFT-D3 release is used here; the reference
C6 values below are random numbers and only serve to test the formulas.

5. Two-electron (Coulomb) integral between two product densities centred at
   different points, by the reciprocal-space partial-wave formula of the
   chapter (book's Fourier convention, complex Y_lm, Gaunt coefficients),
   for s-s and p-s Gaussian densities, with the bare and the erfc-screened
   kernel, against closed forms.

Run: uv run --with numpy --with scipy --with sympy python tools/check_dispersion.py
"""
import itertools
import math

import numpy as np

BOHR_ANG = 0.529177210903          # CODATA 2018
HARTREE_J = 4.3597447222071e-18
AVOGADRO = 6.02214076e23

rng = np.random.default_rng(7)


# ------------------------------------------------------------ 1. units
def check_units():
    # 1 J nm^6 mol^-1 in Hartree bohr^6
    fac = (1.0 / (HARTREE_J * AVOGADRO)) * (1e-9 / (BOHR_ANG * 1e-10)) ** 6
    print(f"[units] 1 J nm^6/mol = {fac:.6f} Ha bohr^6 = {2*fac:.6f} Ry bohr^6")
    # (element, C6 [J nm^6/mol], R0 [Angstrom]) from Grimme 2006, Table 1,
    # compared with Ry bohr^6 / bohr values quoted in the chapter.
    cases = [("H", 0.14, 1.001, 4.857, 1.892),
             ("C", 1.75, 1.452, 60.710, 2.744),
             ("Si", 9.23, 1.716, 320.200, 3.243)]
    for el, c6, r0, c6_ry, r0_b in cases:
        got_c6 = 2 * fac * c6
        got_r0 = r0 / BOHR_ANG
        assert abs(got_c6 - c6_ry) / c6_ry < 2e-3, (el, got_c6)
        assert abs(got_r0 - r0_b) < 2e-3, (el, got_r0)
        print(f"[units] {el}: C6 = {got_c6:9.3f} Ry bohr^6 ({c6_ry}), R0 = {got_r0:.3f} bohr ({r0_b})")


# ------------------------------------------------------------ helpers
def image_range(cell, rcut, mode):
    a = cell
    n = []
    for k in range(3):
        if mode == "length":
            n.append(int(math.ceil(rcut / np.linalg.norm(a[k]))))
        else:
            i, j = (k + 1) % 3, (k + 2) % 3
            normal = np.cross(a[i], a[j])
            normal /= np.linalg.norm(normal)
            n.append(int(math.ceil(rcut / abs(normal @ a[k]))))
    return n


def images(cell, rcut, mode="height"):
    n = image_range(cell, rcut, mode)
    for n1, n2, n3 in itertools.product(*(range(-m, m + 1) for m in n)):
        yield (n1, n2, n3), n1 * cell[0] + n2 * cell[1] + n3 * cell[2]


def strained(pos, cell, eps):
    m = np.eye(3) + eps
    return pos @ m.T, cell @ m.T


# ------------------------------------------------------------ 2. D2
def d2_energy(pos, cell, Z, flags, c6, r0, s6=0.75, d=20.0, rcut=30.0, grad=False):
    nat = len(pos)
    E = 0.0
    g = np.zeros_like(pos)
    st = np.zeros((3, 3))
    imgs = list(images(cell, rcut))
    for I in range(nat):
        for J in range(nat):
            C6 = math.sqrt(c6[Z[I]] * c6[Z[J]])
            R0 = r0[Z[I]] + r0[Z[J]]
            for nvec, L in imgs:
                is0 = nvec == (0, 0, 0)
                if not is0 and flags[J] == 0:
                    continue
                w = 1.0 if (not is0 and flags[I] == 0 and flags[J] == 1) else 0.5
                rv = pos[I] - pos[J] - L
                r2 = rv @ rv
                if not (0.1 < r2 <= rcut * rcut):
                    continue
                r = math.sqrt(r2)
                ex = math.exp(-d * (r / R0 - 1.0))
                f = s6 / (1.0 + ex)
                phi = -C6 * f / r**6
                E += w * phi
                # dphi/dr
                dphi = -C6 * f / r**6 * (ex * (d / R0) / (1.0 + ex) - 6.0 / r)
                # term w*phi(|R_I - R_J - L|): derivative wrt R_I and R_J
                g[I] += w * dphi * rv / r
                g[J] -= w * dphi * rv / r
                st += w * dphi * np.outer(rv, rv) / r
    return (E, g, st) if grad else E


def fd_check(efun, pos, cell, label, h=1e-4, hs=1e-7, tol=1e-7):
    # hs: strain step. It is kept small because a strain moves every pair
    # distance (~ r*hs) and a pair crossing the hard cutoff makes the energy jump.
    E, g, st = efun(pos, cell, grad=True)
    gfd = np.zeros_like(g)
    for a in range(len(pos)):
        for k in range(3):
            p = pos.copy(); p[a, k] += h
            m = pos.copy(); m[a, k] -= h
            gfd[a, k] = (efun(p, cell) - efun(m, cell)) / (2 * h)
    sfd = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            eps = np.zeros((3, 3)); eps[i, j] += hs / 2; eps[j, i] += hs / 2
            pp, cp = strained(pos, cell, eps)
            pm, cm = strained(pos, cell, -eps)
            sfd[i, j] = (efun(pp, cp) - efun(pm, cm)) / (2 * hs)
    sym = 0.5 * (st + st.T)
    eg = np.max(np.abs(g - gfd))
    es = np.max(np.abs(sym - sfd))
    print(f"[{label}] E = {E:.10f}  max|grad-FD| = {eg:.2e}  max|stress-FD| = {es:.2e}  "
          f"sum of gradients = {np.max(np.abs(g.sum(0))):.1e}")
    assert eg < tol and es < tol, label
    return E, g, st


def check_d2():
    cell = np.array([[7.0, 0.0, 0.0], [1.5, 6.5, 0.0], [0.7, -0.9, 8.0]])
    pos = rng.uniform(0, 1, size=(4, 3)) @ cell
    Z = [0, 1, 1, 2]
    c6 = {0: 4.857 / 2, 1: 60.710 / 2, 2: 320.200 / 2}   # Hartree bohr^6
    r0 = {0: 1.892, 1: 2.744, 2: 3.243}                 # bohr
    for flags, label in [([1, 1, 1, 1], "D2 all periodic"), ([0, 1, 1, 0], "D2 mixed flags")]:
        fd_check(lambda p, c, grad=False: d2_energy(p, c, Z, flags, c6, r0, grad=grad),
                 pos, cell, label)


# ------------------------------------------------------------ 3. image range
def check_image_range():
    cell = np.array([[6.0, 0.0, 0.0], [5.5, 2.0, 0.0], [0.0, 0.0, 6.0]])
    rcut = 12.0
    nl = image_range(cell, rcut, "length")
    nh = image_range(cell, rcut, "height")
    big = [m + 6 for m in nh]
    def count(nr):
        c = 0
        for n1, n2, n3 in itertools.product(*(range(-m, m + 1) for m in nr)):
            L = n1 * cell[0] + n2 * cell[1] + n3 * cell[2]
            if 0.1 < L @ L <= rcut**2:
                c += 1
        return c
    cl, ch, cb = count(nl), count(nh), count(big)
    print(f"[image range] |a|-based n = {nl} -> {cl} images,  height-based n = {nh} -> {ch},  reference {cb}")
    assert ch == cb and cl < cb


# ------------------------------------------------------------ 4. D3-type
K1, K2, K3 = 16.0, 4.0 / 3.0, 4.0
ALP6, ALP8 = 14.0, 16.0


def make_d3_params(nspec=3):
    rcov = rng.uniform(0.6, 1.4, nspec) / BOHR_ANG * K2      # scaled radii (bohr)
    q = rng.uniform(2.0, 6.0, nspec)                          # "r2r4" like factors
    nref = [1, 3, 2][:nspec]
    cnref = [np.sort(rng.uniform(0, 4, n)) for n in nref]
    c6ref = {}
    for a in range(nspec):
        for b in range(a, nspec):
            tab = rng.uniform(5, 60, size=(nref[a], nref[b]))
            c6ref[(a, b)] = tab
            c6ref[(b, a)] = tab.T
    r0ab = rng.uniform(4.0, 7.0, size=(nspec, nspec))
    r0ab = 0.5 * (r0ab + r0ab.T)
    return dict(rcov=rcov, q=q, cnref=cnref, c6ref=c6ref, r0ab=r0ab)


def c6_of_cn(P, si, sj, cni, cnj):
    tab = P["c6ref"][(si, sj)]
    ai = cni - P["cnref"][si][:, None]
    bj = cnj - P["cnref"][sj][None, :]
    arg = -K3 * (ai**2 + bj**2)
    mask = tab > 1e-12
    shift = arg[mask].max()
    Lw = np.where(mask, np.exp(arg - shift), 0.0)
    Zs = (tab * Lw).sum(); W = Lw.sum()
    c6 = Zs / W
    dZi = (tab * Lw * (-2 * K3 * ai)).sum(); dWi = (Lw * (-2 * K3 * ai)).sum()
    dZj = (tab * Lw * (-2 * K3 * bj)).sum(); dWj = (Lw * (-2 * K3 * bj)).sum()
    return c6, (dZi * W - dWi * Zs) / W**2, (dZj * W - dWj * Zs) / W**2


def d3_energy(pos, cell, spec, P, damp, s6=1.0, s8=0.722, sr6=1.217, a1=0.4289, a2=4.4407,
              rcut=25.0, cncut=15.0, grad=False):
    nat = len(pos)
    imgs_cn = list(images(cell, cncut))
    imgs = list(images(cell, rcut))
    # coordination numbers and their derivative pieces
    CN = np.zeros(nat)
    cn_terms = []
    for I in range(nat):
        for J in range(nat):
            rc = P["rcov"][spec[I]] + P["rcov"][spec[J]]
            for _, L in imgs_cn:
                rv = pos[I] - pos[J] - L
                r2 = rv @ rv
                if 0.1 < r2 < cncut**2:
                    r = math.sqrt(r2)
                    ex = math.exp(-K1 * (rc / r - 1.0))
                    c = 1.0 / (1.0 + ex)
                    CN[I] += c
                    dc = -c * c * ex * K1 * rc / r2
                    cn_terms.append((I, J, rv, r, dc))
    E = 0.0
    g = np.zeros_like(pos)
    st = np.zeros((3, 3))
    dE_dCN = np.zeros(nat)
    for I in range(nat):
        for J in range(nat):
            si, sj = spec[I], spec[J]
            c6, dc6i, dc6j = c6_of_cn(P, si, sj, CN[I], CN[J])
            c8c6 = 3.0 * P["q"][si] * P["q"][sj]
            G = 0.0
            for _, L in imgs:
                rv = pos[I] - pos[J] - L
                r2 = rv @ rv
                if not (0.1 < r2 < rcut**2):
                    continue
                r = math.sqrt(r2)
                if damp == "zero":
                    r0 = P["r0ab"][si, sj]
                    p6 = (r / (sr6 * r0)) ** (-ALP6); f6 = 1 / (1 + 6 * p6)
                    p8 = (r / (1.0 * r0)) ** (-ALP8); f8 = 1 / (1 + 6 * p8)
                    g6 = s6 * f6 / r**6
                    g8 = s8 * c8c6 * f8 / r**8
                    dg6 = g6 * (6.0 / r) * (-1 + ALP6 * p6 * f6)
                    dg8 = g8 * (2.0 / r) * (-4 + 3 * ALP8 * p8 * f8)
                else:
                    R = a1 * math.sqrt(c8c6) + a2
                    t6 = r**6 + R**6; t8 = r**8 + R**8
                    g6 = s6 / t6; g8 = s8 * c8c6 / t8
                    dg6 = -s6 * 6 * r**5 / t6**2
                    dg8 = -s8 * c8c6 * 8 * r**7 / t8**2
                Gp = g6 + g8
                dGp = dg6 + dg8
                G += Gp
                # E = -1/2 sum_IJL C6 * Gp  (all atoms periodic)
                E += -0.5 * c6 * Gp
                g[I] += -0.5 * c6 * dGp * rv / r
                g[J] -= -0.5 * c6 * dGp * rv / r
                st += -0.5 * c6 * dGp * np.outer(rv, rv) / r
            dE_dCN[I] += -0.5 * dc6i * G
            dE_dCN[J] += -0.5 * dc6j * G
    for I, J, rv, r, dc in cn_terms:
        # CN_I depends on R_I - R_J - L
        g[I] += dE_dCN[I] * dc * rv / r
        g[J] -= dE_dCN[I] * dc * rv / r
        st += dE_dCN[I] * dc * np.outer(rv, rv) / r
    return (E, g, st) if grad else E


def check_d3():
    cell = np.array([[6.5, 0.0, 0.0], [1.0, 6.0, 0.0], [0.5, 0.8, 7.0]])
    pos = rng.uniform(0, 1, size=(4, 3)) @ cell
    spec = [0, 1, 1, 2]
    P = make_d3_params()
    for damp in ("zero", "bj"):
        fd_check(lambda p, c, grad=False: d3_energy(p, c, spec, P, damp, grad=grad),
                 pos, cell, f"D3-type {damp}", tol=2e-7)
    # far-from-reference CN: the shifted weights stay finite and the
    # result tends to the reference closest in CN space
    c6, _, _ = c6_of_cn(P, 1, 1, 40.0, 40.0)
    tab = P["c6ref"][(1, 1)]
    ref = tab[-1, -1]
    print(f"[D3 C6 far from refs] C6(CN=40,40) = {c6:.6f}, nearest reference = {ref:.6f}")
    assert abs(c6 - ref) < 1e-8


# ------------------------------------------------------------ 5. ERI formula
def sph(l, m, theta, phi):
    from scipy.special import sph_harm_y
    return sph_harm_y(l, m, theta, phi)


def gaunt_c(l1, m1, l2, m2, l3, m3):
    """G(l1 m1; l2 m2, l3 m3) = int conj(Y_l1m1) Y_l2m2 Y_l3m3 dOmega."""
    from sympy.physics.wigner import gaunt
    return float((-1) ** m1 * gaunt(l1, l2, l3, -m1, m2, m3))


def eri_partial_wave(PF, PG, Rvec, omega=None, rmax=12.0, nr=3000, kmax=40.0, nk=4000):
    """(F|G) = int int F(r) G(r') v(r - r') with F centred at 0, G at Rvec.
    PF, PG: dict {(l, m): radial function P(r)} so that F = sum P_lm(r) Y_lm."""
    from scipy.special import spherical_jn
    r = np.linspace(1e-6, rmax, nr); wr = np.full(nr, r[1] - r[0]); wr[[0, -1]] *= 0.5
    k = np.linspace(1e-8, kmax, nk); wk = np.full(nk, k[1] - k[0]); wk[[0, -1]] *= 0.5
    def bt(P, l):
        return np.array([np.sum(wr * r**2 * spherical_jn(l, kk * r) * P(r)) for kk in k])
    R = np.linalg.norm(Rvec)
    th = math.acos(Rvec[2] / R); ph = math.atan2(Rvec[1], Rvec[0])
    screen = np.ones_like(k) if omega is None else 1.0 - np.exp(-k**2 / (4 * omega**2))
    tot = 0.0 + 0.0j
    for (l, m), P in PF.items():
        tP = bt(P, l)
        for (lp, mp), Q in PG.items():
            tQ = bt(Q, lp)
            for lam in range(abs(l - lp), l + lp + 1):
                mu = m - mp
                if abs(mu) > lam:
                    continue
                G = gaunt_c(l, m, lp, mp, lam, mu)
                if abs(G) < 1e-14:
                    continue
                rad = np.sum(wk * np.conj(tP) * tQ * spherical_jn(lam, k * R) * screen)
                tot += 32 * math.pi * (1j) ** (l - lp - lam) * G * np.conj(sph(lam, mu, th, ph)) * rad
    return tot


def check_eri():
    from scipy.special import erf
    a, b = 1.3, 0.8
    Rvec = np.array([0.7, -1.1, 1.6]); R = np.linalg.norm(Rvec)
    y00 = 1 / math.sqrt(4 * math.pi)
    Fs = {(0, 0): lambda r: np.exp(-a * r * r) / y00}
    Gs = {(0, 0): lambda r: np.exp(-b * r * r) / y00}
    mu = math.sqrt(a * b / (a + b)); qa = (math.pi / a) ** 1.5; qb = (math.pi / b) ** 1.5
    ss = lambda Rv, g=mu: qa * qb * erf(g * np.linalg.norm(Rv)) / np.linalg.norm(Rv)
    got = eri_partial_wave(Fs, Gs, Rvec)
    print(f"[ERI s-s] partial-wave {got.real:.10f} (imag {got.imag:.1e})  analytic {ss(Rvec):.10f}")
    assert abs(got - ss(Rvec)) < 1e-6 * abs(ss(Rvec))
    # p-s: F = x exp(-a r^2) = r exp(-a r^2) sqrt(2 pi/3) (Y_1,-1 - Y_11)
    c = math.sqrt(2 * math.pi / 3)
    Fp = {(1, -1): lambda r: c * r * np.exp(-a * r * r), (1, 1): lambda r: -c * r * np.exp(-a * r * r)}
    got = eri_partial_wave(Fp, Gs, Rvec)
    h = 1e-5; e = np.array([h, 0, 0])
    ref = -(1 / (2 * a)) * (ss(Rvec + e) - ss(Rvec - e)) / (2 * h)
    print(f"[ERI p-s] partial-wave {got.real:.10f} (imag {got.imag:.1e})  analytic {ref:.10f}")
    assert abs(got - ref) < 1e-6 * abs(ref)
    # erfc-screened kernel: 4 pi / k^2 (1 - exp(-k^2 / 4 w^2))
    w = 0.6
    gam = 1 / math.sqrt(1 / a + 1 / b + 1 / w**2)
    got = eri_partial_wave(Fs, Gs, Rvec, omega=w)
    ref = ss(Rvec) - ss(Rvec, gam)
    print(f"[ERI s-s erfc] partial-wave {got.real:.10f}  analytic {ref:.10f}")
    assert abs(got - ref) < 1e-6 * abs(ref)


if __name__ == "__main__":
    check_units()
    check_d2()
    check_image_range()
    check_d3()
    check_eri()
    print("all checks passed")
