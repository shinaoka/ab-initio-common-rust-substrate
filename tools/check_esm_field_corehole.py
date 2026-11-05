# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/esm-field-corehole.md:
#   1. field unit conversion GV/m -> Hartree/(e bohr)
#   2. ESM Green's functions (vacuum|vacuum, metal|metal (+bias), vacuum|metal):
#      spectral "periodic part + homogeneous correction" vs direct quadrature,
#      and the closed-form interval integral for vacuum|vacuum
#   3. sawtooth potential and its energy
#   4. core-hole projectors |R Phi_JM><R Phi_JM| in the real-harmonic basis
#   5. binding-energy arithmetic of the manual examples
# Run: uv run --with numpy --with scipy python tools/check_esm_field_corehole.py

import math
import numpy as np
from scipy.special import factorial

ok_all = True


def check(name, err, tol):
    global ok_all
    flag = err <= tol
    ok_all &= flag
    print(f"[{'OK ' if flag else 'NG '}] {name}: err={err:.3e} tol={tol:.1e}")


HARTREE_EV = 27.211386245988
BOHR_M = 0.529177210903e-10

# ---------------------------------------------------------------- 1. units
au_field = HARTREE_EV / BOHR_M          # V/m per Hartree/(e bohr)
gv = 1.0e9 / au_field
check("1 GV/m in Hartree/(e bohr) = 1.944690e-3", abs(gv - 1.944690e-3), 1e-9)
check("factor used (0.1944688e-2) vs CODATA", abs(gv - 0.1944688e-2) / gv, 2e-6)

# ---------------------------------------------------------------- 2. ESM
L = 24.0
z0 = L / 2
N = 192
h = L / N
zg = (np.arange(N) - N // 2) * h               # grid on [-z0, z0)
m_idx = np.fft.fftfreq(N, d=1.0 / N)           # integers
g = 2 * np.pi * m_idx / L


def density(z):
    # slab-like: positive lump plus smaller negative lump (net charge != 0)
    return (np.exp(-((z - 1.0) / 1.2) ** 2) - 0.4 * np.exp(-((z + 2.0) / 0.9) ** 2))


nz = density(zg)
nm = np.fft.fft(nz) / N * np.cos(np.pi * m_idx)  # n(z) = sum_m n_m exp(i g z), z measured from the centre
Q = np.sum(nz) * h


def green(bc, Gp, z, zp, dV=0.0):
    zl, zr = np.minimum(z, zp), np.maximum(z, zp)
    z1 = z0
    if Gp > 0:
        if bc == "on1":
            return 2 * np.pi / Gp * np.exp(-Gp * np.abs(z - zp))
        if bc == "on2":
            return 4 * np.pi * np.sinh(Gp * (zl + z1)) * np.sinh(Gp * (z1 - zr)) / (Gp * np.sinh(2 * Gp * z1))
        if bc == "on3":
            return 2 * np.pi / Gp * (np.exp(-Gp * np.abs(z - zp)) - np.exp(Gp * (z + zp - 2 * z1)))
    else:
        if bc == "on1":
            return -2 * np.pi * np.abs(z - zp)
        if bc == "on2":
            return 4 * np.pi * (zl + z1) * (z1 - zr) / (2 * z1)
        if bc == "on3":
            return 4 * np.pi * (z1 - zr)
    raise ValueError


def direct(bc, Gp, z, dV=0.0):
    zf = np.linspace(-z0, z0, 40001)
    nf = density(zf)
    w = np.full_like(zf, zf[1] - zf[0]); w[0] *= 0.5; w[-1] *= 0.5
    V = np.array([np.sum(w * green(bc, Gp, zi, zf) * nf) for zi in z])
    if bc == "on2" and Gp == 0:
        V = V + dV * (z0 - z) / (2 * z0)
    return V


def spectral(bc, Gp, z, dV=0.0):
    """periodic solution + homogeneous correction fixed by the boundary conditions."""
    G2 = Gp ** 2 + g ** 2
    mask = G2 > 0
    E = np.exp(1j * np.outer(z, g))
    Ep = np.exp(1j * np.outer(np.array([-z0, z0]), g))
    coef = np.zeros(N, complex); coef[mask] = 4 * np.pi * nm[mask] / G2[mask]
    Vp = (E @ coef).real
    Vb = (Ep @ coef).real                      # values at -z0, +z0 (equal: periodic)
    dVb = (Ep @ (1j * g * coef)).real          # derivatives at -z0, +z0
    n0 = nm[0].real if Gp == 0 else 0.0
    if Gp == 0:
        # particular solution of the uniform part: -2 pi n0 z^2
        Vp = Vp - 2 * np.pi * n0 * z ** 2
        Vb = Vb - 2 * np.pi * n0 * z0 ** 2
        dVb = dVb + np.array([4 * np.pi * n0 * z0, -4 * np.pi * n0 * z0])
        f = [lambda x: np.ones_like(x), lambda x: x]
        df = [lambda x: np.zeros_like(x), lambda x: np.ones_like(x)]
    else:
        f = [lambda x: np.exp(Gp * x), lambda x: np.exp(-Gp * x)]
        df = [lambda x: Gp * np.exp(Gp * x), lambda x: -Gp * np.exp(-Gp * x)]
    lo, hi = np.array([-z0]), np.array([z0])
    A = np.zeros((2, 2)); r = np.zeros(2)
    if bc == "on1":
        if Gp > 0:   # V' = -Gp V at +z0, V' = +Gp V at -z0
            for k in range(2):
                A[0, k] = df[k](hi)[0] + Gp * f[k](hi)[0]
                A[1, k] = df[k](lo)[0] - Gp * f[k](lo)[0]
            r = [-(dVb[1] + Gp * Vb[1]), -(dVb[0] - Gp * Vb[0])]
        else:        # field -2piQ / +2piQ and V(z0)+V(-z0) = -4 pi z0 Q
            for k in range(2):
                A[0, k] = df[k](hi)[0]
                A[1, k] = f[k](hi)[0] + f[k](lo)[0]
            r = [-2 * np.pi * Q - dVb[1], -4 * np.pi * z0 * Q - Vb[0] - Vb[1]]
    elif bc == "on2":   # V(+-z0) = 0
        for k in range(2):
            A[0, k] = f[k](hi)[0]; A[1, k] = f[k](lo)[0]
        r = [-Vb[1], -Vb[0]]
    elif bc == "on3":   # V(z0) = 0, vacuum below
        for k in range(2):
            A[0, k] = f[k](hi)[0]
            A[1, k] = (df[k](lo)[0] - Gp * f[k](lo)[0]) if Gp > 0 else df[k](lo)[0]
        r = [-Vb[1], -((dVb[0] - Gp * Vb[0]) if Gp > 0 else dVb[0])]
    c = np.linalg.solve(A, r)
    V = Vp + c[0] * f[0](z) + c[1] * f[1](z)
    if bc == "on2" and Gp == 0:
        V = V + dV * (z0 - z) / (2 * z0)
    return V


def closed_on1(Gp, z):
    """closed-form interval integral for vacuum|vacuum (book Eqs.)."""
    s = np.cos(np.pi * m_idx)                   # (-1)^m
    out = np.zeros_like(z, dtype=complex)
    for k in range(N):
        gk, nk, sk = g[k], nm[k], s[k]
        if Gp > 0:
            G2 = Gp ** 2 + gk ** 2
            out += nk * (4 * np.pi * np.exp(1j * gk * z) / G2
                         - 2 * np.pi * sk / (Gp * G2) * ((Gp + 1j * gk) * np.exp(Gp * (z - z0))
                                                        + (Gp - 1j * gk) * np.exp(-Gp * (z + z0))))
        elif gk == 0:
            out += -2 * np.pi * nk * (z ** 2 + z0 ** 2)
        else:
            out += nk * (4 * np.pi * np.exp(1j * gk * z) / gk ** 2 - 4j * np.pi * sk * z / gk
                         - 4 * np.pi * sk / gk ** 2)
    return out.real


zt = np.linspace(-z0 + 0.5, z0 - 0.5, 41)
for bc in ("on1", "on2", "on3"):
    for Gp in (0.0, 0.25, 1.3):
        dV = 0.3 if bc == "on2" else 0.0
        err = np.max(np.abs(spectral(bc, Gp, zt, dV) - direct(bc, Gp, zt, dV)))
        check(f"ESM {bc} Gp={Gp}: spectral vs direct Green's function", err, 1e-6)
for Gp in (0.0, 0.25, 1.3):
    err = np.max(np.abs(closed_on1(Gp, zt) - direct("on1", Gp, zt)))
    check(f"ESM on1 Gp={Gp}: closed-form interval integral vs direct", err, 1e-6)

# Green's function symmetry and boundary values
zs = np.linspace(-z0, z0, 7)
for bc in ("on1", "on2", "on3"):
    for Gp in (0.0, 0.7):
        Gm = np.array([[green(bc, Gp, a, b) for b in zs] for a in zs])
        check(f"G symmetric {bc} Gp={Gp}", np.max(np.abs(Gm - Gm.T)), 1e-12)
check("on2: G(+-z1, z') = 0", max(abs(green("on2", 0.7, z0, 1.0)), abs(green("on2", 0.7, -z0, 1.0)),
                                  abs(green("on2", 0.0, z0, 1.0)), abs(green("on2", 0.0, -z0, 1.0))), 1e-12)
check("on3: G(z1, z') = 0", max(abs(green("on3", 0.7, z0, 1.0)), abs(green("on3", 0.0, z0, 1.0))), 1e-12)
# induced charge: on3 Gp=0 -> field at the metal = -4 pi Q (all lines end on the metal)
V = spectral("on3", 0.0, np.array([z0 - 1e-4, z0]))
check("on3 Gp=0: dV/dz at metal = -4 pi Q", abs((V[1] - V[0]) / 1e-4 + 4 * np.pi * Q), 1e-3)
# bias: on2 with zero density -> V = dV at -z1 and 0 at +z1
check("on2 bias term: V(-z1)=dV, V(z1)=0",
      abs(0.3 * (z0 - (-z0)) / (2 * z0) - 0.3) + abs(0.3 * (z0 - z0) / (2 * z0)), 1e-15)

# ---------------------------------------------------------------- 3. sawtooth
a = 20.0; Ng = 80; E = 0.01; c = 7.3
xg = np.arange(Ng) * a / Ng
Vef = E * (xg - c)
jump = Vef[-1] + E * a / Ng - Vef[0]
check("sawtooth jump across the cell face = E a", abs(jump - E * a), 1e-12)
# energy of a narrow normalized lump at x0 ~ E (x0 - c)
x0 = 11.0
lump = np.exp(-((xg - x0) / 0.8) ** 2); lump /= lump.sum() * (a / Ng)
Uef = np.sum(lump * Vef) * (a / Ng)
check("Uef of unit lump = E (x0 - c)", abs(Uef - E * (x0 - c)), 1e-10)

# ---------------------------------------------------------------- 4. core-hole projectors
# real harmonics in the book's order, built from polynomials
def real_harm(l, x, y, z):
    if l == 0:
        return np.array([np.ones_like(x)])
    if l == 1:
        return np.array([x, y, z])
    if l == 2:
        return np.array([3 * z * z - 1, x * x - y * y, x * y, x * z, y * z])
    if l == 3:
        return np.array([z * (5 * z * z - 3), x * (5 * z * z - 1), y * (5 * z * z - 1), z * (x * x - y * y),
                         x * y * z, x * (x * x - 3 * y * y), y * (3 * x * x - y * y)])


def complex_harm(l, th, ph):
    out = []
    for m in range(-l, l + 1):
        am = abs(m)
        # associated Legendre with Condon-Shortley phase
        from scipy.special import lpmv
        P = lpmv(am, l, np.cos(th))
        nrm = math.sqrt((2 * l + 1) / (4 * math.pi) * factorial(l - am) / factorial(l + am))
        Y = nrm * P * np.exp(1j * am * ph)
        if m < 0:
            Y = (-1) ** am * np.conj(Y)
        out.append(Y)
    return np.array(out)


xg_, wg_ = np.polynomial.legendre.leggauss(40)
th = np.arccos(xg_); ph = np.linspace(0, 2 * np.pi, 81)[:-1]
TH, PH = np.meshgrid(th, ph, indexing="ij")
W = np.outer(wg_, np.full(ph.size, 2 * np.pi / ph.size))
X, Y_, Z = np.sin(TH) * np.cos(PH), np.sin(TH) * np.sin(PH), np.cos(TH)


def cg_spinor(l, J2, M2):
    """coefficients of Phi_JM on (Y_lm alpha, Y_lm beta), J2=2J, M2=2M."""
    c = np.zeros((2, 2 * l + 1), complex)  # [spin][m+l]
    if J2 == 2 * l + 1:
        m = (M2 - 1) // 2
        if -l <= m <= l:
            c[0, m + l] = math.sqrt((l + m + 1) / (2 * l + 1))
        if -l <= m + 1 <= l:
            c[1, m + 1 + l] = math.sqrt((l - m) / (2 * l + 1))
    else:
        m = (M2 + 1) // 2
        if -l <= m - 1 <= l:
            c[0, m - 1 + l] = math.sqrt((l - m + 1) / (2 * l + 1))
        if -l <= m <= l:
            c[1, m + l] = -math.sqrt((l + m) / (2 * l + 1))
    return c


diag_expect = {
    (1, 3, 3): {(0, 0): 0.5, (0, 1): 0.5},                       # p J=3/2 M=3/2: x up, y up
    (2, 5, 1): {(0, 0): 0.6, (1, 3): 0.2, (1, 4): 0.2},           # d J=5/2 M=1/2
    (2, 5, 3): {(0, 3): 0.4, (0, 4): 0.4, (1, 1): 0.1, (1, 2): 0.1},
    (1, 1, -1): {(0, 0): 1 / 3, (0, 1): 1 / 3, (1, 2): 1 / 3},  # p J=1/2 M=-1/2
}
for l in (0, 1, 2, 3):
    S = real_harm(l, X, Y_, Z)
    S = S / np.sqrt(np.einsum("kij,kij,ij->k", S, S, W))[:, None, None]
    Yc = complex_harm(l, TH, PH)
    U = np.einsum("aij,bij,ij->ab", S, np.conj(Yc), W)     # S_a = sum_b U_ab Y_b
    check(f"U^(l={l}) unitary", np.max(np.abs(U @ U.conj().T - np.eye(2 * l + 1))), 1e-10)
    total = np.zeros((2 * (2 * l + 1),) * 2, complex)
    states = [(2 * l + 1, M2) for M2 in range(2 * l + 1, -2 * l - 2, -2)]
    if l > 0:
        states += [(2 * l - 1, M2) for M2 in range(2 * l - 1, -2 * l, -2)]
    check(f"number of (J,M) states for l={l} = 2(2l+1)", abs(len(states) - 2 * (2 * l + 1)), 0)
    for (J2, M2) in states:
        c = cg_spinor(l, J2, M2)
        # coefficient on real basis: Y_b = sum_a conj(U_ab) S_a
        creal = np.concatenate([np.conj(U) @ c[0], np.conj(U) @ c[1]])
        P = np.outer(creal, np.conj(creal))
        total += P
        key = (l, J2, M2)
        if key in diag_expect:
            d = np.real(np.diag(P))
            exp = np.zeros_like(d)
            for (s, a_), v in diag_expect[key].items():
                exp[s * (2 * l + 1) + a_] = v
            check(f"diag of projector l={l} J={J2}/2 M={M2}/2", np.max(np.abs(d - exp)), 1e-10)
    check(f"sum over (J,M) projectors = identity (l={l})", np.max(np.abs(total - np.eye(total.shape[0]))), 1e-10)

# ---------------------------------------------------------------- 5. binding energies (manual 60.2-60.3)
eb = (-66.084858926233 - (-76.787732114928)) * HARTREE_EV
check("C2H2 C-1s: Ef(N-1)-Ei(N) = 291.24 eV", abs(eb - 291.24), 0.005)
eb = (-10489.553360141708 - (-10499.900104007471)) * HARTREE_EV
check("TiC C-1s (metal, neutral final state) = 281.55 eV", abs(eb - 281.55), 0.005)
mu0 = -0.201641583
eb32 = (-34816.628201335407 - (-34820.483255130872) + mu0) * HARTREE_EV
eb12 = (-34816.601864921540 - (-34820.483255130872) + mu0) * HARTREE_EV
check("Si 2p3/2 = 99.41 eV", abs(eb32 - 99.41), 0.005)
check("Si 2p1/2 = 100.13 eV", abs(eb12 - 100.13), 0.005)

print("ALL OK" if ok_all else "SOME CHECKS FAILED")
