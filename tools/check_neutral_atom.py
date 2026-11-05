# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/neutral-atom.md.
#
#  (1) Sphere-sphere Coulomb integral in k-space,
#      int n_I V_H,J = (2/pi) int dk j0(kR) N_I(k) N_J(k),  N = 4 pi int r^2 j0 n,
#      against the closed form for two Gaussian charge clouds.
#  (2) "Product density x spherical potential" two-center formula
#      <phi_a| V(|r - tau|) |phi_b>, both orbitals on the same center,
#      = 8 sum_LM G(lm; l'm', LM) Y*_LM(tau^) int dk k^2 j_L(k tau) P_L(k) V(k)
#      against direct 3D grid integration.
#  (3) Projector expansion of a short-range potential with radial functions
#      orthogonalised in the V-weighted metric: exactness for functions in the
#      span (on-site case), and convergence of a genuine three-center element
#      with L_max, against direct 3D grid integration.
#  (4) The smooth cutoff polynomial (C2 step over a 1 bohr buffer).
#
# Model: n(r) = Z (alpha/pi)^{3/2} exp(-alpha r^2), V_H = Z erf(sqrt(alpha) r)/r,
#        V_loc = -Z erf(sqrt(beta) r)/r, V_na = V_loc + V_H (short range).
# Run (takes several minutes because of the 3D grid in part 3):
#   uv run --with numpy --with scipy --with sympy python tools/check_neutral_atom.py

import numpy as np
from scipy.special import spherical_jn, sph_harm_y, erf
from scipy.integrate import quad
from sympy.physics.wigner import gaunt

np.set_printoptions(precision=8, suppress=True)

# ---------------------------------------------------------------- model
Z, alpha, beta = 4.0, 0.9, 3.0


def vna(r):
    r = np.maximum(r, 1e-12)
    return Z * (erf(np.sqrt(alpha) * r) - erf(np.sqrt(beta) * r)) / r


def Y(l, m, th, ph):
    return sph_harm_y(l, m, th, ph)


def G(l1, m1, l2, m2, L, M):
    # int conj(Y_l1m1) Y_l2m2 Y_LM
    return float((-1) ** m1 * gaunt(l1, l2, L, -m1, m2, M))


def Sreal(l, m, th, ph):
    # real spherical harmonics (any orthonormal real combination is fine here)
    if m == 0:
        return Y(l, 0, th, ph).real
    if m > 0:
        return np.sqrt(2) * (-1) ** m * Y(l, m, th, ph).real
    return np.sqrt(2) * (-1) ** m * Y(l, -m, th, ph).imag


# ---------------------------------------------------------------- (1)
print("(1) sphere-sphere Coulomb integral")
a1, a2, Z1, Z2 = 0.7, 1.6, 3.0, 5.0


def Nk(Zc, a, k):  # 3D Fourier transform of normalised Gaussian charge
    return Zc * np.exp(-k * k / (4 * a))


for R in [0.5, 1.5, 3.0]:
    num = 2 / np.pi * quad(lambda k: spherical_jn(0, k * R) * Nk(Z1, a1, k) * Nk(Z2, a2, k),
                           0, 60, limit=400)[0]
    g = np.sqrt(a1 * a2 / (a1 + a2))
    exact = Z1 * Z2 * erf(g * R) / R
    print(f"  R={R:4.1f}  k-space={num:.12f}  exact={exact:.12f}  diff={num-exact:.1e}")

# ---------------------------------------------------------------- 3D grid
h = 0.11
x = np.arange(-7.5, 7.5 + 1e-9, h)
X, Yg, Zg = np.meshgrid(x, x, x, indexing="ij")
dV = h ** 3


def sph(Xc, Yc, Zc):
    r = np.sqrt(Xc * Xc + Yc * Yc + Zc * Zc) + 1e-300
    return r, np.arccos(np.clip(Zc / r, -1, 1)), np.arctan2(Yc, Xc)


def orb(l, m, a, c, real=True):
    r, th, ph = sph(X - c[0], Yg - c[1], Zg - c[2])
    ang = Sreal(l, m, th, ph) if real else Y(l, m, th, ph)
    return r ** l * np.exp(-a * r * r) * ang


# ---------------------------------------------------------------- (2)
print("(2) product density x spherical potential (two orbitals on one center)")
tau = np.array([0.6, -0.4, 1.1])  # V centered at tau, orbitals at the origin
t = np.linalg.norm(tau)
tth, tph = np.arccos(tau[2] / t), np.arctan2(tau[1], tau[0])
rV, _, _ = sph(X - tau[0], Yg - tau[1], Zg - tau[2])
Vgrid = vna(rV)
ea, eb = 0.8, 1.2


def Vk(k):  # radial transform int r^2 j0(kr) V(r)
    return quad(lambda r: r * r * spherical_jn(0, k * r) * vna(r), 0, 14, limit=300)[0]


kk = np.linspace(1e-6, 25, 1201)
Vtab = np.array([Vk(k) for k in kk])


def Ptab(L, la, lb):
    return np.array([quad(lambda r: r * r * spherical_jn(L, k * r) * r ** (la + lb)
                          * np.exp(-(ea + eb) * r * r), 0, 12, limit=200)[0] for k in kk])


for (la, ma, lb, mb) in [(0, 0, 0, 0), (1, 0, 0, 0), (1, 1, 1, 1), (2, -1, 1, 0)]:
    s = 0
    for L in range(abs(la - lb), la + lb + 1):
        M = mb - ma
        if abs(M) > L:
            continue
        # conj(Y_la ma) Y_lb mb = sum_LM G(la ma; lb mb, L M') conj(Y_LM') with M' = ma - mb
        g = G(la, ma, lb, mb, L, ma - mb)
        if g == 0:
            continue
        P = Ptab(L, la, lb)
        I = np.trapezoid(kk * kk * spherical_jn(L, kk * t) * P * Vtab, kk)
        s += 8 * g * np.conj(Y(L, ma - mb, tth, tph)) * I
    fa = orb(la, ma, ea, [0, 0, 0], real=False)
    fb = orb(lb, mb, eb, [0, 0, 0], real=False)
    direct = np.sum(np.conj(fa) * Vgrid * fb) * dV
    print(f"  {(la,ma,lb,mb)}  formula={complex(np.round(s,8))}  grid={complex(np.round(direct,8))}")

# ---------------------------------------------------------------- (3)
print("(3) projector expansion with V-weighted Gram-Schmidt")
rr = np.linspace(1e-6, 12, 6001)
Vr = vna(rr)
seeds_exp = [0.5, 0.9, 1.6, 2.8, 4.5, 7.0]


def rint(f):
    return np.trapezoid(f * rr * rr, rr)


def projectors(l, nrad):
    seeds = [rr ** l * np.exp(-a * rr * rr) for a in seeds_exp[:nrad]]
    bars, cs = [], []
    for s in seeds:
        b = s.copy()
        for bb, cc in zip(bars, cs):
            b -= bb * rint(bb * Vr * s) / cc
        bars.append(b)
        cs.append(rint(b * Vr * b))
    return bars, cs


def proj_element(fa, fb, Lmax, nrad):
    rK, thK, phK = sph(X, Yg, Zg)  # projector center K at the origin
    VK = vna(rK)
    tot = 0.0
    for l in range(Lmax + 1):
        bars, cs = projectors(l, nrad)
        for m in range(-l, l + 1):
            ang = Sreal(l, m, thK, phK)
            for b, c in zip(bars, cs):
                beta_fun = VK * np.interp(rK, rr, b) * ang
                tot += np.sum(fa * beta_fun) * dV * np.sum(beta_fun * fb) * dV / c
    return tot


rK, _, _ = sph(X, Yg, Zg)
VK = vna(rK)
# (3a) on-site: phi_a on K with an exponent in the seed set -> exact
fa = orb(1, 1, 0.9, [0, 0, 0])
fb = orb(0, 0, 1.0, [1.3, 0.5, -0.8])
exact = np.sum(fa * VK * fb) * dV
approx = proj_element(fa, fb, Lmax=3, nrad=3)
print(f"  i=K case: exact={exact:.8f}  projector(Lmax=3,N=3)={approx:.8f}")
# (3b) genuine three-center: phi_a, phi_b off K
fa = orb(0, 0, 1.0, [1.0, 0.3, 0.2])
fb = orb(1, 0, 1.1, [-0.4, 0.9, -0.7])
exact = np.sum(fa * VK * fb) * dV
print(f"  three-center exact={exact:.8f}")
for Lmax in [0, 2, 4, 6]:
    for nrad in [2, 4, 6]:
        print(f"    Lmax={Lmax} N={nrad}  projector={proj_element(fa, fb, Lmax, nrad):.8f}")

# ---------------------------------------------------------------- (4)
print("(4) smooth cutoff polynomial f(x)=-x^3(10+15x+6x^2), x=R-Rcut, buffer 1")
A = np.array([[-1, 1, -1], [3, -4, 5], [-6, 12, -20]], float)  # f, f', f'' at x=-1
c3, c4, c5 = np.linalg.solve(A, [1, 0, 0])
print("  solved coefficients c3,c4,c5 =", c3, c4, c5, " expected -10,-15,-6")
