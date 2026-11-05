# SPDX-License-Identifier: MIT
"""Numerical checks for the nonlocal pseudopotential chapter.

Run:
    uv run --with numpy --with scipy --with sympy python tools/check_nonlocal.py

Checks (all use only the conventions written in this book):
  1. Angular-momentum matrices in the real spherical harmonics of
     src/modules/pao-basis.md (l = 1, 2, 3): L^2 = l(l+1), [Lx, Ly] = i Lz.
  2. j-projectors P_+ = ((l+1) + L.sigma)/(2l+1), P_- = (l - L.sigma)/(2l+1):
     idempotent, complementary, traces 2j+1, and equal to sum_M |Phi_JM><Phi_JM|
     built from Clebsch-Gordan spinors (complex Y_lm, Condon-Shortley phase).
  3. Hermiticity and time-reversal symmetry of an assembled SO block.
  4. j-averaging: spin trace of the j-dependent operator.
  5. MBK diagonalisation: B^{-1} form == sum_k |beta_k> sign(b_k) <beta_k|.
  6. Damping function f(t) = 10 t^3 - 15 t^4 + 6 t^5 (t = rc - r, 1 bohr buffer).
  7. Gaussian three-centre test: V_ab = eps * S(a,K) S(K,b) and its gradient
     with respect to the projector centre vs. finite differences.
  8. SO-factor scaling of the projectors: residual L.sigma when eps_+ != eps_-.
"""
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

try:
    from scipy.special import sph_harm_y as _sph

    def Ylm(l, m, theta, phi):
        return _sph(l, m, theta, phi)
except ImportError:  # older scipy
    from scipy.special import sph_harm as _sph

    def Ylm(l, m, theta, phi):
        return _sph(m, l, phi, theta)

x, y, z = sp.symbols("x y z", real=True)
r2 = x * x + y * y + z * z
POLY = {
    0: [sp.Integer(1)],
    1: [x, y, z],
    2: [3 * z**2 - r2, x**2 - y**2, x * y, x * z, y * z],
    3: [z * (5 * z**2 - 3 * r2), x * (5 * z**2 - r2), y * (5 * z**2 - r2),
        z * (x**2 - y**2), x * y * z, x * (x**2 - 3 * y**2), y * (3 * x**2 - y**2)],
}

# product quadrature on the unit sphere (exact for the degrees used here)
u, wu = leggauss(24)
NPHI = 48
phi = np.linspace(0, 2 * np.pi, NPHI, endpoint=False)
U_, PHI = np.meshgrid(u, phi, indexing="ij")
W = np.outer(wu, np.full(NPHI, 2 * np.pi / NPHI))
TH = np.arccos(U_)
X = np.sin(TH) * np.cos(PHI)
Y = np.sin(TH) * np.sin(PHI)
Z = U_


def ev(e):
    return sp.lambdify((x, y, z), e, "numpy")(X, Y, Z) * np.ones_like(X)


def real_harmonics(l):
    S = []
    exprs = []
    for p in POLY[l]:
        v = ev(p)
        c = np.sqrt(np.sum(W * v * v))
        S.append(v / c)
        exprs.append(p / c)
    return S, exprs


def L_matrices(l):
    S, ex = real_harmonics(l)
    n = 2 * l + 1
    M = {k: np.zeros((n, n), complex) for k in "xyz"}
    for nu, f in enumerate(ex):
        Lf = {
            "x": -sp.I * (y * sp.diff(f, z) - z * sp.diff(f, y)),
            "y": -sp.I * (z * sp.diff(f, x) - x * sp.diff(f, z)),
            "z": -sp.I * (x * sp.diff(f, y) - y * sp.diff(f, x)),
        }
        for k in "xyz":
            v = ev(sp.expand(Lf[k]))
            for mu in range(n):
                M[k][mu, nu] = np.sum(W * S[mu] * v)
    return M, S


SIG = {
    "x": np.array([[0, 1], [1, 0]], complex),
    "y": np.array([[0, -1j], [1j, 0]]),
    "z": np.array([[1, 0], [0, -1]], complex),
}


def LdotS(M):
    # spin index outer, orbital index inner: (sigma, mu)
    return sum(np.kron(SIG[k], M[k]) for k in "xyz")


def projectors(l, M):
    n = 2 * l + 1
    I = np.eye(2 * n)
    Ls = LdotS(M)
    return ((l + 1) * I + Ls) / (2 * l + 1), (l * I - Ls) / (2 * l + 1)


def check(name, ok):
    print(f"[{'OK' if ok else 'NG'}] {name}")
    if not ok:
        raise SystemExit(1)


def main():
    rng = np.random.default_rng(1)
    for l in (1, 2, 3):
        n = 2 * l + 1
        M, S = L_matrices(l)
        L2 = sum(M[k] @ M[k] for k in "xyz")
        check(f"l={l}: L^2 = l(l+1)", np.allclose(L2, l * (l + 1) * np.eye(n), atol=1e-12))
        check(f"l={l}: [Lx,Ly] = i Lz",
              np.allclose(M["x"] @ M["y"] - M["y"] @ M["x"], 1j * M["z"], atol=1e-12))
        check(f"l={l}: Lz, L- purely imaginary / real-antisym structure",
              np.allclose(M["z"].real, 0, atol=1e-12))

        Pp, Pm = projectors(l, M)
        check(f"l={l}: P+^2=P+, P-^2=P-, P+P-=0",
              np.allclose(Pp @ Pp, Pp) and np.allclose(Pm @ Pm, Pm) and np.allclose(Pp @ Pm, 0))
        check(f"l={l}: tr P+ = 2l+2, tr P- = 2l",
              np.isclose(np.trace(Pp).real, 2 * l + 2) and np.isclose(np.trace(Pm).real, 2 * l))

        # U^{(l)}: S_mu = sum_m U[mu,m] Y_lm
        ms = list(range(-l, l + 1))
        Yc = [Ylm(l, m, TH, PHI) for m in ms]
        Umat = np.array([[np.sum(W * np.conj(Yc[j]) * S[mu]) for j in range(n)] for mu in range(n)])
        check(f"l={l}: U unitary", np.allclose(Umat @ Umat.conj().T, np.eye(n), atol=1e-10))

        # Clebsch-Gordan spinors |l 1/2 j M> in the complex basis (spin outer, m inner)
        def cg_proj(j):
            P = np.zeros((2 * n, 2 * n), complex)
            twoM_list = range(-int(2 * j), int(2 * j) + 1, 2)
            for twoM in twoM_list:
                Mj = twoM / 2
                v = np.zeros(2 * n, complex)
                for s, ms_ in ((0, 0.5), (1, -0.5)):
                    m = Mj - ms_
                    if abs(m) > l or abs(m - round(m)) > 1e-9:
                        continue
                    m = int(round(m))
                    if np.isclose(j, l + 0.5):
                        c = np.sqrt((l + 0.5 + 2 * ms_ * Mj) / (2 * l + 1))
                    else:
                        c = -2 * ms_ * np.sqrt((l + 0.5 - 2 * ms_ * Mj) / (2 * l + 1))
                    v[s * n + (m + l)] = c
                P += np.outer(v, v.conj())
            return P

        # Y_m = sum_mu conj(U[mu,m]) S_mu, so coefficients map by conj(U)
        T = np.kron(np.eye(2), Umat.conj())
        for j, P in ((l + 0.5, Pp), (l - 0.5, Pm)):
            Pc = T @ cg_proj(j) @ T.conj().T
            check(f"l={l}, j={j}: P_j equals Clebsch-Gordan projector", np.allclose(Pc, P, atol=1e-10))

        # assembled SO block for one third atom: <phi_a|beta_mu> = A[a,mu] (real)
        na, nb = 4, 3
        A = rng.normal(size=(na, n))
        B = rng.normal(size=(nb, n))
        Aa = rng.normal(size=(na, n))  # j=l-1/2 overlaps
        Bb = rng.normal(size=(nb, n))
        ep, em = 0.7, -1.3

        def block(Ap, Bp, Am, Bm):
            # H[(sigma,a),(sigma',b)]
            Hp = np.kron(np.eye(2), Ap) @ Pp @ np.kron(np.eye(2), Bp.T)
            Hm = np.kron(np.eye(2), Am) @ Pm @ np.kron(np.eye(2), Bm.T)
            return ep * Hp + em * Hm

        Hab = block(A, B, Aa, Bb)
        Hba = block(B, A, Bb, Aa)
        check(f"l={l}: H_ab = (H_ba)^dagger", np.allclose(Hab, Hba.conj().T))
        # time reversal Theta = i sigma_y K acting on spin
        Ty = np.kron(1j * SIG["y"], np.eye(na))
        Tyb = np.kron(1j * SIG["y"], np.eye(nb))
        check(f"l={l}: time reversal (i s_y) H* (i s_y)^dag = H",
              np.allclose(Ty @ Hab.conj() @ Tyb.conj().T, Hab))

        # j-averaging: half spin trace of the j-dependent operator
        bp = rng.normal(size=n)
        bm = rng.normal(size=n)
        half = 0.5 * (np.trace(Pp.reshape(2, n, 2, n), axis1=0, axis2=2))
        check(f"l={l}: (1/2) Tr_spin P+ = (l+1)/(2l+1) * 1",
              np.allclose(half, (l + 1) / (2 * l + 1) * np.eye(n)))
        half_m = 0.5 * (np.trace(Pm.reshape(2, n, 2, n), axis1=0, axis2=2))
        check(f"l={l}: (1/2) Tr_spin P- = l/(2l+1) * 1", np.allclose(half_m, l / (2 * l + 1) * np.eye(n)))

        # SO factor: beta_pm -> bbar +- f dbeta ; L.sigma coefficient
        f_so = 0.0
        bbar, db = 0.5 * (bp + bm), 0.5 * (bp - bm)
        b1, b2 = bbar + f_so * db, bbar - f_so * db
        Ls = LdotS(M)
        H = ep * np.kron(np.eye(2), np.outer(b1, b1)) @ Pp + em * np.kron(np.eye(2), np.outer(b2, b2)) @ Pm
        # L.sigma part: coefficient (ep - em)/(2l+1) on |bbar><bbar|
        G = np.kron(np.eye(2), np.outer(bbar, bbar))
        expect = G @ (((l + 1) * ep + l * em) * np.eye(2 * n) + (ep - em) * Ls) / (2 * l + 1)
        check(f"l={l}: SO factor 0 leaves (eps+ - eps-)/(2l+1) L.sigma", np.allclose(H, expect))

    # 5. MBK diagonalisation
    nb_ = 3
    Bm = rng.normal(size=(nb_, nb_))
    Bm = Bm + Bm.T
    chi = rng.normal(size=(50, nb_))
    V1 = chi @ np.linalg.inv(Bm) @ chi.T
    bvals, Uv = np.linalg.eigh(Bm)
    beta = chi @ Uv / np.sqrt(np.abs(bvals))
    V2 = beta @ np.diag(np.sign(bvals)) @ beta.T
    check("MBK: B^{-1} form equals diagonal +-1 form", np.allclose(V1, V2))

    # 6. damping function
    t = sp.symbols("t")
    fpoly = 10 * t**3 - 15 * t**4 + 6 * t**5
    ok = (fpoly.subs(t, 1) == 1 and sp.diff(fpoly, t).subs(t, 1) == 0
          and sp.diff(fpoly, t, 2).subs(t, 1) == 0 and fpoly.subs(t, 0) == 0
          and sp.diff(fpoly, t).subs(t, 0) == 0 and sp.diff(fpoly, t, 2).subs(t, 0) == 0)
    check("damping: f(0)=f'(0)=f''(0)=0, f(1)=1, f'(1)=f''(1)=0", bool(ok))

    # 7. Gaussian three-centre test (s orbitals, s projector)
    def s_overlap(a, b, tau):
        mu = a * b / (a + b)
        return (2 * np.sqrt(a * b) / (a + b)) ** 1.5 * np.exp(-mu * tau @ tau)

    def grad_s_overlap(a, b, tau):
        mu = a * b / (a + b)
        return -2 * mu * tau * s_overlap(a, b, tau)

    al, be, ga, eps = 0.8, 1.3, 2.0, -1.0
    RA = np.array([0.0, 0.0, 0.0])
    RB = np.array([1.1, 0.3, -0.2])
    RK = np.array([0.4, -0.5, 0.6])

    def V(RK):
        return eps * s_overlap(al, ga, RK - RA) * s_overlap(ga, be, RB - RK)

    gK = eps * (grad_s_overlap(al, ga, RK - RA) * s_overlap(ga, be, RB - RK)
                - s_overlap(al, ga, RK - RA) * grad_s_overlap(ga, be, RB - RK))
    h = 1e-5
    fd = np.array([(V(RK + h * e) - V(RK - h * e)) / (2 * h) for e in np.eye(3)])
    check("Gaussian 3-centre: dV/dR_K analytic vs finite difference", np.allclose(gK, fd, atol=1e-9))

    # direct numerical check of the 3-centre value on a grid
    g = np.linspace(-6, 6, 121)
    dx = g[1] - g[0]
    Xg, Yg, Zg = np.meshgrid(g, g, g, indexing="ij")

    def gauss(a, R):
        return (2 * a / np.pi) ** 0.75 * np.exp(-a * ((Xg - R[0]) ** 2 + (Yg - R[1]) ** 2 + (Zg - R[2]) ** 2))

    pa, pb, pk = gauss(al, RA), gauss(be, RB), gauss(ga, RK)
    Sak = np.sum(pa * pk) * dx**3
    Skb = np.sum(pk * pb) * dx**3
    check("Gaussian 3-centre: S(a,K) S(K,b) grid vs analytic",
          np.isclose(eps * Sak * Skb, V(RK), rtol=1e-6))
    print("all checks passed")


if __name__ == "__main__":
    main()
