# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/analysis.md (post-processing analyses):
#  1. Mulliken sum rule and per-state Mulliken weights (non-orthogonal 1D chain).
#  2. Fuzzy-cell (Becke) weights: partition of unity.
#  3. Linear tetrahedron DOS (6 tetrahedra per cube) for a simple-cubic s band:
#     normalization and integrated DOS vs. direct counting; Gaussian DOS normalization.
#  4. Berry-phase polarization with lattice-only Bloch phases and the e^{-i b.r}
#     link operator: equality with the Resta supercell formula, gauge invariance,
#     and the Thouless pump (integer charge per cycle) for the Rice-Mele chain.
#  5. Chern number by the Fukui-Hatsugai-Suzuki lattice field strength (QWZ model).
#  6. Exchange coupling (Liechtenstein form used by the band mode):
#     relation to the force theorem (mixed second derivative of the grand potential
#     under spin rotations), and pole-expansion (continued fraction) form
#     J(i0,jR) summed over R == periodic-sum form.
# Run: uv run --with numpy --with scipy python tools/check_analysis.py

import itertools
import numpy as np
from scipy.linalg import eigh

rng = np.random.default_rng(7)
ok = True


def check(name, err, tol):
    global ok
    flag = err <= tol
    ok &= bool(flag)
    print(f"{'OK ' if flag else 'NG '} {name}: err={err:.3e} tol={tol:.0e}")


def fermi(e, mu, kT):
    x = np.clip((e - mu) / kT, -200.0, 200.0)
    return 1.0 / (1.0 + np.exp(x))


# ---------------------------------------------------------------------------
# 1D chain with two atoms per cell, one orbital each; L in {-1,0,1}.
# O(k) = sum_L e^{i k L} O(L), k = 2 pi q  (lattice-only phase, as in the book).
# ---------------------------------------------------------------------------
def chain_model(eA, eB, t1, t2, s1, s2):
    """H(L), S(L) for A(0)-B(1/2) chain: intra-cell A-B (t1,s1), inter-cell B-A(L=+1) (t2,s2)."""
    H = {L: np.zeros((2, 2)) for L in (-1, 0, 1)}
    S = {L: np.zeros((2, 2)) for L in (-1, 0, 1)}
    H[0][:] = [[eA, t1], [t1, eB]]
    S[0][:] = [[1.0, s1], [s1, 1.0]]
    H[1][1, 0] = t2
    S[1][1, 0] = s2           # <B(0)|A(L=+1)>
    H[-1] = H[1].T.copy()
    S[-1] = S[1].T.copy()
    return H, S


def bloch(O, q):
    return sum(np.exp(2j * np.pi * q * L) * O[L] for L in O)


# 1. Mulliken -----------------------------------------------------------------
H, S = chain_model(-0.3, 0.4, -1.0, -0.6, 0.15, 0.1)
Nk = 12
qs = np.arange(Nk) / Nk
kT, mu = 0.02, 0.05
rho = {L: np.zeros((2, 2)) for L in (-1, 0, 1)}
Ne = 0.0
gross_from_states = np.zeros(2)
for q in qs:
    e, c = eigh(bloch(H, q), bloch(S, q))
    f = fermi(e, mu, kT)
    Ne += f.sum() / Nk
    for L in rho:
        rho[L] += np.real(np.exp(2j * np.pi * q * L)
                          * np.einsum("n,an,bn->ab", f, c.conj(), c)) / Nk
    # per-state, per-orbital Mulliken weight w_a = Re[conj(c_a) (S c)_a]
    w = np.real(c.conj() * (bloch(S, q) @ c))
    check("per-state Mulliken weights sum to 1", np.max(abs(w.sum(0) - 1)), 1e-12)
    gross_from_states += (w * f).sum(1) / Nk
gross = np.array([sum((rho[L] * S[L])[a].sum() for L in rho) for a in range(2)])
check("Mulliken sum rule sum_I M_I = N_e", abs(gross.sum() - Ne), 1e-12)
check("Mulliken gross pop. == sum_states f*w", np.max(abs(gross - gross_from_states)), 1e-12)

# 2. Fuzzy-cell weights ---------------------------------------------------------
def becke_s(mu_, napply):
    f = mu_
    for _ in range(napply):
        f = 1.5 * f - 0.5 * f**3
    return 0.5 * (1.0 - f)


atoms = rng.uniform(-2, 2, (5, 3))
pts = rng.uniform(-3, 3, (200, 3))
for napply in (3, 6):
    P = np.ones((len(pts), len(atoms)))
    for a, b in itertools.permutations(range(len(atoms)), 2):
        Rab = np.linalg.norm(atoms[a] - atoms[b])
        mu_ = (np.linalg.norm(pts - atoms[a], axis=1) - np.linalg.norm(pts - atoms[b], axis=1)) / Rab
        P[:, a] *= becke_s(mu_, napply)
    W = P / P.sum(1, keepdims=True)
    check(f"fuzzy weights partition of unity (p applied {napply} times)", np.max(abs(W.sum(1) - 1)), 1e-14)


# 3. Tetrahedron DOS ------------------------------------------------------------
def tetra_dos(et, e):
    e1, e2, e3, e4 = et
    out = np.zeros_like(e)
    m = (e > e1) & (e < e2)
    out[m] = 3 * (e[m] - e1) ** 2 / ((e2 - e1) * (e3 - e1) * (e4 - e1))
    m = (e > e2) & (e < e3)
    x = e[m] - e2
    out[m] = (3 * (e2 - e1) + 6 * x - 3 * (e3 - e1 + e4 - e2) * x * x / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1))
    m = (e > e3) & (e < e4)
    out[m] = 3 * (e4 - e[m]) ** 2 / ((e4 - e1) * (e4 - e2) * (e4 - e3))
    return out


TET = [(0, 1, 2, 5), (1, 2, 3, 5), (2, 3, 5, 7), (0, 2, 4, 5), (2, 4, 5, 6), (2, 5, 6, 7)]
Nm = 12
g = np.arange(Nm) / Nm
E = -2 * sum(np.cos(2 * np.pi * x) for x in np.meshgrid(g, g, g, indexing="ij"))
egrid = np.linspace(-6.5, 6.5, 2601)
dos = np.zeros_like(egrid)
for i, j, k in itertools.product(range(Nm), repeat=3):
    corner = [E[(i + a) % Nm, (j + b) % Nm, (k + c) % Nm] for a in (0, 1) for b in (0, 1) for c in (0, 1)]
    for t in TET:
        dos += tetra_dos(np.sort([corner[v] for v in t]), egrid)
dos /= 6 * Nm**3
de = egrid[1] - egrid[0]
check("tetrahedron DOS integrates to 1", abs(dos.sum() * de - 1), 2e-3)
Ef = -1.0
nf = 200
gf = (np.arange(nf) + 0.5) / nf
Ef_count = np.mean(-2 * sum(np.cos(2 * np.pi * x) for x in np.meshgrid(gf, gf, gf, indexing="ij")) < Ef)
check("tetrahedron N(E_F) vs fine counting", abs(dos[egrid < Ef].sum() * de - Ef_count), 5e-3)
a = 0.1
gdos = sum(np.exp(-((egrid - e0) / a) ** 2) / (a * np.sqrt(np.pi)) for e0 in E.ravel()[::7]) / len(E.ravel()[::7])
check("Gaussian DOS integrates to 1", abs(gdos.sum() * de - 1), 1e-6)


# 4. Berry phase ----------------------------------------------------------------
def rice_mele(phase, d0=0.3, t0=1.0, dt=0.4):
    eA, eB = d0 * np.cos(phase), -d0 * np.cos(phase)
    t1, t2 = -(t0 + dt * np.sin(phase)), -(t0 - dt * np.sin(phase))
    return chain_model(eA, eB, t1, t2, 0.0, 0.0)[0]


tau = np.array([0.0, 0.5])   # fractional positions of A, B (lattice constant 1)


def berry_pol(H, Nk, random_gauge=False):
    """Electronic polarization (in units of e*a, per spin) of the lowest band,
    lattice-only Bloch phases, link M = c(k)^+ T c(k+b), T = diag(e^{-i b tau})."""
    b = 2 * np.pi / Nk
    T = np.diag(np.exp(-1j * b * tau))
    cs = []
    for n in range(Nk):
        _, c = eigh(bloch(H, n / Nk))
        c = c[:, :1]
        if random_gauge:
            c = c * np.exp(2j * np.pi * rng.uniform())
        cs.append(c)
    cs.append(cs[0])          # periodic gauge: lattice-only phase => c(k+G) = c(k)
    prod = 1.0 + 0j
    for n in range(Nk):
        prod *= np.linalg.det(cs[n].conj().T @ T @ cs[n + 1])
    phi = -np.angle(prod)
    return -phi / (2 * np.pi)   # P_el = -(e/2pi) phi * a  (electron charge -e)


def resta_pol(H, Ncell):
    """Gamma-point supercell of Ncell cells, P = -(1/2pi) Im ln det <Psi|e^{i 2pi x/(N a)}|Psi>."""
    n = 2 * Ncell
    Hs = np.zeros((n, n))
    x = np.zeros(n)
    for c in range(Ncell):
        x[2 * c:2 * c + 2] = c + tau
        for L in H:
            c2 = (c + L) % Ncell
            Hs[2 * c:2 * c + 2, 2 * c2:2 * c2 + 2] += H[L]
    _, C = eigh(Hs)
    occ = C[:, :Ncell]
    U = occ.conj().T @ np.diag(np.exp(2j * np.pi * x / Ncell)) @ occ
    return -np.angle(np.linalg.det(U)) / (2 * np.pi)


def wrap(x):
    return (x + 0.5) % 1.0 - 0.5


H_rm = rice_mele(0.7)
p1 = berry_pol(H_rm, 15)   # odd N: the cyclic block permutation in the supercell det has sign +1
p2 = resta_pol(H_rm, 15)
check("Berry string == Resta supercell (mod 1)", abs(wrap(p1 - p2)), 1e-10)
check("Berry phase gauge invariant", abs(wrap(berry_pol(H_rm, 15, True) - p1)), 1e-12)
# Thouless pump: accumulate continuous change over one cycle
phases = np.linspace(0, 2 * np.pi, 201)
P = [berry_pol(rice_mele(p), 24) for p in phases]
dP = sum(wrap(P[i + 1] - P[i]) for i in range(len(P) - 1))
check("Thouless pump: |Delta P| over a cycle is an integer (1)", abs(abs(dP) - 1.0), 1e-6)


# 5. Chern number (FHS) -----------------------------------------------------------
def qwz(kx, ky, u):
    return np.array([[u + np.cos(kx) + np.cos(ky), np.sin(kx) - 1j * np.sin(ky)],
                     [np.sin(kx) + 1j * np.sin(ky), -(u + np.cos(kx) + np.cos(ky))]])


def chern_fhs(u, N=24):
    vec = np.empty((N, N, 2), complex)
    for i, j in itertools.product(range(N), repeat=2):
        _, c = np.linalg.eigh(qwz(2 * np.pi * i / N, 2 * np.pi * j / N, u))
        vec[i, j] = c[:, 0]

    def U(a, b):
        z = np.vdot(a, b)
        return z / abs(z)
    F = 0.0
    for i, j in itertools.product(range(N), repeat=2):
        i1, j1 = (i + 1) % N, (j + 1) % N
        F += np.angle(U(vec[i, j], vec[i1, j]) * U(vec[i1, j], vec[i1, j1])
                      * U(vec[i1, j1], vec[i, j1]) * U(vec[i, j1], vec[i, j]))
    return F / (2 * np.pi)


for u, cexp in ((-1.0, 1), (1.0, -1), (2.5, 0)):
    c = chern_fhs(u)
    check(f"FHS Chern number QWZ u={u} -> |C|={abs(cexp)}", abs(abs(c) - abs(cexp)), 1e-9)


# 6. Exchange coupling ------------------------------------------------------------
def omega_gc(e, mu, kT):
    return -kT * np.logaddexp(0.0, -(e - mu) / kT)


# 6a. cluster, force theorem
kT, mu = 0.01, -0.4
H0 = np.array([[0.0, -1.0, -0.3], [-1.0, 0.2, -0.8], [-0.3, -0.8, -0.1]])
S0 = np.array([[1, 0.2, 0.05], [0.2, 1, 0.1], [0.05, 0.1, 1]])
Dl = np.array([-1.5, -1.2, 0.0])          # Delta_I = H_up - H_dn (on-site)
sx = np.array([[0, 1], [1, 0]])
sz = np.array([[1, 0], [0, -1]])


def omega_rot(theta):
    n = 3
    Hs = np.kron(np.eye(2), H0).astype(complex)
    for I in range(n):
        m = np.sin(theta[I]) * sx + np.cos(theta[I]) * sz
        for s1, s2 in itertools.product(range(2), repeat=2):
            Hs[s1 * n + I, s2 * n + I] += Dl[I] / 2 * m[s1, s2]
    return omega_gc(eigh(Hs, np.kron(np.eye(2), S0), eigvals_only=True), mu, kT).sum()


def j_band_formula(eu, Cu, ed, Cd, Pi, Pj, mu, kT):
    fu, fd = fermi(eu, mu, kT), fermi(ed, mu, kT)
    Mj = Cu.conj().T @ Pj @ Cd
    Mi = Cd.conj().T @ Pi @ Cu
    return 0.25 * np.sum((fu[:, None] - fd[None, :]) / (ed[None, :] - eu[:, None]) * Mj * Mi.T)


i, j, h = 0, 1, 1e-3
th = lambda a, b: np.array([a, b, 0.0])
mixed = (omega_rot(th(h, h)) - omega_rot(th(h, -h)) - omega_rot(th(-h, h)) + omega_rot(th(-h, -h))) / (4 * h * h)
eu, Cu = eigh(H0 + np.diag(Dl / 2), S0)
ed, Cd = eigh(H0 - np.diag(Dl / 2), S0)
Pi = np.diag(np.where(np.arange(3) == i, Dl, 0.0))
Pj = np.diag(np.where(np.arange(3) == j, Dl, 0.0))
Jb = j_band_formula(eu, Cu, ed, Cd, Pi, Pj, mu, kT).real
check("force theorem: -d2Omega/dth_i dth_j == 2 J (band-mode prefactor 1/4)", abs(-mixed - 2 * Jb), 1e-6)

# 6b. periodic chain: pole sum of J(i0,jR) over R == periodic-sum formula
Hc, Sc = chain_model(0.0, 0.1, -1.0, -0.7, 0.1, 0.05)
Dc = np.array([-1.2, -0.9])
Hu = {L: Hc[L] + (np.diag(Dc / 2) if L == 0 else 0) for L in Hc}
Hd = {L: Hc[L] - (np.diag(Dc / 2) if L == 0 else 0) for L in Hc}
Nk, kT, mu = 10, 0.05, -0.3
qs = (np.arange(Nk) * 2 + 1) / (2 * Nk) - 0.5        # mesh used by the band mode
Pi = np.diag([Dc[0], 0.0])
Pj = np.diag([0.0, Dc[1]])
data = []
Jpsum = 0.0
for q in qs:
    eu, Cu = eigh(bloch(Hu, q), bloch(Sc, q))
    ed, Cd = eigh(bloch(Hd, q), bloch(Sc, q))
    data.append((q, eu, Cu, ed, Cd))
    Jpsum += j_band_formula(eu, Cu, ed, Cd, Pi, Pj, mu, kT).real / Nk
# continued-fraction poles of the Fermi function
Np = 60
M = 2 * Np
A = np.zeros((M, M))
for a in range(M - 1):
    A[a, a + 1] = A[a + 1, a] = -0.5
lam, V = eigh(A, np.diag(2 * np.arange(1, M + 1) - 1.0))
zp = -1 / lam[:Np]
Rp = -V[0, :Np] ** 2 * zp**2 / 4
xs = np.linspace(-30, 30, 13)
fapprox = 0.5 + sum(Rp[p] * 2 * xs / (xs**2 + zp[p] ** 2) for p in range(Np))
check("continued-fraction poles reproduce the Fermi function", np.max(abs(fapprox - 1 / (1 + np.exp(xs)))), 1e-12)
Jsum = 0.0
for R in range(-Nk // 2, Nk // 2):
    JR = 0.0
    for p in range(Np):
        Ep = mu + 1j * kT * zp[p]
        Gu = np.zeros((2, 2), complex)   # G^up_{ji}(E; -R) : sum_k e^{+i k R} C (E-e)^-1 C^+
        Gd = np.zeros((2, 2), complex)   # G^dn_{ij}(E; +R) : sum_k e^{-i k R} ...
        for q, eu, Cu, ed, Cd in data:
            ph = np.exp(2j * np.pi * q * R)
            Gu += ph * (Cu / (Ep - eu)) @ Cu.conj().T
            Gd += np.conj(ph) * (Cd / (Ep - ed)) @ Cd.conj().T
        JR += Rp[p] * np.real(np.trace(Pi @ Gd @ Pj @ Gu))
    Jsum += kT / 2 * JR / Nk**2
check("sum_R J_pole(i0,jR) == J_periodic-sum", abs(Jsum - Jpsum), 1e-10)

print("ALL OK" if ok else "SOME CHECKS FAILED")
