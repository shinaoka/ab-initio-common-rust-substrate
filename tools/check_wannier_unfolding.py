# SPDX-License-Identifier: MIT
"""Numerical checks for the chapter "Wannier functions, band unfolding and
localized-orbital analysis".

Everything here is a toy model written from the formulas in the chapter
(book notation).  No OpenMX code is used.

  1. Band unfolding (LCAO, non-orthogonal basis) on a 1D chain:
     - perfect supercell: weights are 0 or 1 and reproduce the primitive bands
     - sum rule: sum over the n unfolded k of W/n = 1 (also with defects)
     - the variant with the conjugate on the other coefficient breaks the sum rule
  2. Closest Wannier functions (polar decomposition of the projection matrix):
     - F[U] = mean (sigma-1)^2, F[U] <= F[B] for random partial unitaries
     - gauge invariance of the tight-binding parameters
     - full window + all orbitals: Wannier-interpolated bands = original bands
     - polynomial window: C^1 at the inner edge, C^2 at the outer edge
  3. Maximally localized Wannier functions (Marzari-Vanderbilt, 1D):
     - Omega_I is gauge invariant
     - steepest descent lowers the spread; complete set -> Omega -> 0
  4. Orbital optimization: analytic gradient of the band energy with respect
     to contraction coefficients equals a finite difference.

Run:  uv run --with numpy --with scipy python tools/check_wannier_unfolding.py
"""
import numpy as np
from scipy.linalg import eigh, svd, expm

rng = np.random.default_rng(1234)
OK = True


def check(name, cond, info=""):
    global OK
    OK &= bool(cond)
    print(f"[{'ok' if cond else 'NG'}] {name} {info}")


# ---------------------------------------------------------------------------
# 1D Gaussian-orbital model: s and p orbitals, extended-Hueckel Hamiltonian
# ---------------------------------------------------------------------------
ALPHA = 0.7
XG = np.linspace(-40, 40, 16001)
DX = XG[1] - XG[0]


def orb(t, x):
    g = np.exp(-ALPHA * x * x)
    return g / np.sqrt(np.sum(g * g) * DX) if t == 0 else \
        x * g / np.sqrt(np.sum((x * g) ** 2) * DX)


def overlap(t1, t2, d):
    """<phi_t1(x)|phi_t2(x-d)>"""
    return np.sum(orb(t1, XG) * orb(t2, XG - d)) * DX


S_CACHE = {}


def S2(t1, t2, d):
    key = (t1, t2, round(d, 10))
    if key not in S_CACHE:
        S_CACHE[key] = overlap(t1, t2, d)
    return S_CACHE[key]


ONSITE = {"A": (-0.30, 0.10), "B": (-0.45, 0.05)}  # (s, p) on-site energies
KWH = 1.75
RCUT = 9.0


def matrices(pos, spec, cell, k):
    """H(k), S(k) with the book convention: phase exp(+i k L) on the ket image."""
    nat = len(pos)
    n = 2 * nat
    H = np.zeros((n, n), complex)
    S = np.zeros((n, n), complex)
    nimg = int(RCUT / cell) + 2
    for i in range(nat):
        for j in range(nat):
            for m in range(-nimg, nimg + 1):
                L = m * cell
                d = pos[j] + L - pos[i]
                if abs(d) > RCUT:
                    continue
                ph = np.exp(1j * k * L)
                for a in range(2):
                    for b in range(2):
                        s = S2(a, b, d)
                        if i == j and m == 0:
                            h = ONSITE[spec[i]][a] if a == b else 0.0
                        else:
                            h = 0.5 * KWH * (ONSITE[spec[i]][a] + ONSITE[spec[j]][b]) * s
                        S[2 * i + a, 2 * j + b] += s * ph
                        H[2 * i + a, 2 * j + b] += h * ph
    return H, S


def solve(H, S):
    e, C = eigh(H, S)  # C^dagger S C = 1
    return e, C


# ---------------------------------------------------------------------------
# 1. band unfolding
# ---------------------------------------------------------------------------
print("== 1. band unfolding ==")
a = 2.2      # reference-cell length
n = 3        # number of reference cells in the supercell
A = n * a
x0 = 0.3


def unfold_weights(pos, spec, group, kref, origin, conj_on_ket=False):
    """Return eigenvalues and W/n (per state) at the reference-cell k = kref.

    W = sum_{N home, ket image (J,L)} sum_{M: group(M)=group(J)}
        conj(C_N) C_M exp(i k (r(J,L) - r'(M))) S_{N,(J,L)}   (orbital-matched)
    """
    H, S = matrices(pos, spec, A, kref)  # Bloch sum at k itself (= K + G)
    e, C = solve(H, S)
    nat = len(pos)
    rcell = lambda x: np.floor((x - origin) / a) * a
    W = np.zeros(len(e), complex)
    nimg = int(RCUT / A) + 2
    for i in range(nat):
        for j in range(nat):
            for m in range(-nimg, nimg + 1):
                L = m * A
                d = pos[j] + L - pos[i]
                if abs(d) > RCUT:
                    continue
                rJ = rcell(pos[j] + L)
                for M in range(nat):
                    if group[M] != group[j]:
                        continue
                    ph = np.exp(1j * kref * (rJ - rcell(pos[M])))
                    for aa in range(2):
                        for bb in range(2):
                            s = S2(aa, bb, d)
                            if conj_on_ket:
                                W += np.conj(C[2 * M + bb]) * C[2 * i + aa] * ph * s
                            else:
                                W += np.conj(C[2 * i + aa]) * C[2 * M + bb] * ph * s
    return e, W / n


# perfect supercell
pos = [x0 + j * a for j in range(n)]
spec = ["A"] * n
group = [1] * n
Kpt = 0.37 * 2 * np.pi / A
ks = [Kpt + 2 * np.pi * g / A for g in range(n)]
ws = []
for kk in ks:
    e, w = unfold_weights(pos, spec, group, kk, origin=-0.1)
    ws.append(w)
ws = np.array(ws)
check("perfect: imaginary part of W/n ~ 0", np.max(np.abs(ws.imag)) < 1e-10,
      f"max={np.max(np.abs(ws.imag)):.1e}")
wr = ws.real
check("perfect: W/n in {0,1}", np.max(np.minimum(np.abs(wr), np.abs(wr - 1))) < 1e-8,
      f"dev={np.max(np.minimum(np.abs(wr), np.abs(wr - 1))):.1e}")
check("perfect: sum_j W/n = 1", np.max(np.abs(wr.sum(0) - 1)) < 1e-10)
# energies with weight 1 at k_j equal primitive-cell bands at k_j
maxdev = 0.0
e_sc, _ = solve(*matrices(pos, spec, A, Kpt))
for jj, kk in enumerate(ks):
    ep, _ = solve(*matrices([x0], ["A"], a, kk))
    sel = e_sc[wr[jj] > 0.5]
    maxdev = max(maxdev, np.max(np.abs(np.sort(sel) - np.sort(ep))))
check("perfect: unfolded energies = primitive bands", maxdev < 1e-10, f"dev={maxdev:.1e}")
# origin independence in the perfect case
ws2 = np.array([unfold_weights(pos, spec, group, kk, origin=-0.73)[1] for kk in ks])
check("perfect: weights independent of origin", np.max(np.abs(ws2 - ws)) < 1e-10)

# imperfect supercell: one atom substituted and displaced
pos_d = [x0, x0 + a + 0.25, x0 + 2 * a]
spec_d = ["A", "B", "A"]
group_d = [1, 1, 1]
wd = np.array([unfold_weights(pos_d, spec_d, group_d, kk, origin=-0.1)[1] for kk in ks])
check("defect: sum_j W/n = 1", np.max(np.abs(wd.sum(0) - 1)) < 1e-10,
      f"dev={np.max(np.abs(wd.sum(0) - 1)):.1e}")
check("defect: some weights fractional", np.any((wd.real > 0.05) & (wd.real < 0.95)))
wv = np.array([unfold_weights(pos_d, spec_d, group_d, kk, origin=-0.1, conj_on_ket=True)[1]
               for kk in ks])
check("variant with conj on relabelled coefficient breaks sum rule",
      np.max(np.abs(wv.sum(0) - 1)) > 1e-3, f"dev={np.max(np.abs(wv.sum(0) - 1)):.2e}")

# ---------------------------------------------------------------------------
# 2. closest Wannier functions
# ---------------------------------------------------------------------------
print("== 2. closest Wannier functions ==")
ALPHA = 1.2                   # better-conditioned overlap for this part
S_CACHE.clear()
cell = 3.0
pos_c = [0.0, 1.5]
spec_c = ["A", "B"]          # 2 atoms x (s,p) = 4 basis functions
nb = 4
Nk = 9                        # odd grid
kgrid = 2 * np.pi * (np.arange(Nk) / Nk) / cell


def poly_window(x, o0, i0, i1, o1):
    c = [1.0, 0.0, -6.0, 8.0, -3.0]
    P = lambda y: sum(ci * y ** p for p, ci in enumerate(c))
    if x <= o0 or x > o1:
        return 0.0
    if x <= i0:
        return P((i0 - x) / (i0 - o0))
    if x <= i1:
        return 1.0
    return P((x - i1) / (o1 - i1))


def fermi_window(x, e0, e1, kT0, kT1, delta):
    y0 = np.clip((e0 - x) / kT0, -700, 700)
    y1 = np.clip((x - e1) / kT1, -700, 700)
    return 1.0 / (np.exp(y0) + 1) / (np.exp(y1) + 1) + delta


def cwf(Qcols, wfun, phases=None):
    """Return list of (A, U, e) per k and hopping t(R) for R=-(Nk-1)/2..(Nk-1)/2."""
    data = []
    for ik, k in enumerate(kgrid):
        H, S = matrices(pos_c, spec_c, cell, k)
        e, C = solve(H, S)
        if phases is not None:
            C = C * phases[ik]
        w = np.array([wfun(x) for x in e])
        Amat = (w[:, None]) * (C.conj().T @ S @ Qcols)  # a_{k mu,p} = w <psi|Q>
        Wl, sig, Vh = svd(Amat, full_matrices=False)
        U = Wl @ Vh
        data.append((Amat, U, e, sig, C))
    Rs = np.arange(-(Nk - 1) // 2, (Nk - 1) // 2 + 1)
    t = {}
    for R in Rs:
        acc = 0
        for ik, k in enumerate(kgrid):
            _, U, e, _, _ = data[ik]
            acc = acc + np.exp(-1j * k * R * cell) * (U.conj().T @ np.diag(e) @ U)
        t[R] = acc / Nk
    return data, t


# guiding functions: s orbitals of both atoms (AO guiding), window on lower bands
Q = np.zeros((nb, 2))
Q[0, 0] = 1.0
Q[2, 1] = 1.0
win = lambda x: fermi_window(x, -2.0, 0.0, 0.05, 0.3, 1e-12)
data, t = cwf(Q, win)
F_svd = np.mean([np.sum((d[3] - 1) ** 2) for d in data])
F_dir = np.mean([np.linalg.norm(d[0] - d[1]) ** 2 for d in data])
check("F[U] = mean sum (sigma-1)^2", abs(F_svd - F_dir) < 1e-12, f"F={F_svd:.6f}")
worst = np.inf
for _ in range(200):
    Fb = 0.0
    for d in data:
        X = rng.normal(size=(nb, 2)) + 1j * rng.normal(size=(nb, 2))
        Bq, _ = np.linalg.qr(X)
        Fb += np.linalg.norm(d[0] - Bq) ** 2
    worst = min(worst, Fb / Nk - F_dir)
check("F[U] <= F[B] for random partial unitaries", worst >= -1e-12, f"min gap={worst:.3e}")
phases = [np.exp(2j * np.pi * rng.random(nb)) for _ in kgrid]
_, t2 = cwf(Q, win, phases)
dev = max(np.max(np.abs(t[R] - t2[R])) for R in t)
check("hopping t(R) invariant under Bloch phases", dev < 1e-12, f"dev={dev:.1e}")
dev = max(np.max(np.abs(t[R].imag)) for R in t)
check("collinear: t(R) real", dev < 1e-10, f"max|Im|={dev:.1e}")
# full window, all orbitals -> interpolation exact on the grid
data_f, tf = cwf(np.eye(nb), lambda x: 1.0)
dev = 0.0
for ik, k in enumerate(kgrid):
    Hk = sum(np.exp(1j * k * R * cell) * tf[R] for R in tf)
    dev = max(dev, np.max(np.abs(np.linalg.eigvalsh(Hk) - data_f[ik][2])))
check("all orbitals, w=1: interpolated = original on grid", dev < 1e-10, f"dev={dev:.1e}")
kmid = 0.5 * (kgrid[1] + kgrid[2])
Hk = sum(np.exp(1j * kmid * R * cell) * tf[R] for R in tf)
emid = solve(*matrices(pos_c, spec_c, cell, kmid))[0]
dev = np.max(np.abs(np.linalg.eigvalsh(Hk) - emid))
check("all orbitals, w=1: interpolation off grid (Nk=9, approximate)", dev < 2e-2,
      f"dev={dev:.1e}")
# Loewdin equivalence: U = A (A^dag A)^{-1/2} when A full rank
d0 = data_f[0]                # full window: A has full rank
Aq = d0[0]
ev, V = np.linalg.eigh(Aq.conj().T @ Aq)
U_l = Aq @ V @ np.diag(ev ** -0.5) @ V.conj().T
check("polar factor = Loewdin A S^-1/2", np.max(np.abs(U_l - d0[1])) < 1e-8)
# polynomial window smoothness
o0, i0, i1, o1 = -3.0, -2.0, 0.5, 2.0
h = 1e-6
for x, order, name in [(i0, 1, "inner-low"), (i1, 1, "inner-high"),
                       (o0, 2, "outer-low"), (o1, 2, "outer-high")]:
    f = lambda y: poly_window(y, o0, i0, i1, o1)
    d1 = abs((f(x + h) - f(x)) / h - (f(x) - f(x - h)) / h)
    ok = abs(f(x + h) - f(x - h)) < 1e-5 and d1 < 1e-4
    check(f"poly window continuous (C^{order}) at {name}", ok)

# ---------------------------------------------------------------------------
# 3. MLWF spread (Marzari-Vanderbilt), 1D orthonormal point-orbital chain
# ---------------------------------------------------------------------------
print("== 3. maximally localized Wannier functions ==")
Nk = 16
cellm = 1.0
tau = np.array([0.0, 0.5])
kk = 2 * np.pi * np.arange(Nk) / Nk / cellm
bvec = 2 * np.pi / (Nk * cellm)
bs = [bvec, -bvec]
wb = 1.0 / (2 * bvec ** 2)        # 1D: sum_b w_b b b = 1


def ssh(k, t1=1.0, t2=0.6):
    hk = t1 + t2 * np.exp(1j * k * cellm)  # book convention (lattice phase only)
    return np.array([[0, hk], [np.conj(hk), 0]])


Ck = [np.linalg.eigh(ssh(k))[1] for k in kk]


def Mmat(U, ib):
    """M_mn(k,b) = <u_mk|u_n,k+b> for point orbitals: sum_i c*_i c_i exp(-i b tau_i)."""
    out = []
    for ik in range(Nk):
        sh = 1 if bs[ib] > 0 else -1
        jk = (ik + sh) % Nk
        G = (kk[ik] + bs[ib]) - kk[jk]  # 0 or a reciprocal vector
        cj = Ck[jk] @ U[jk]
        # c(k+b) = c(k+b-G) for the lattice-only phase convention
        ci = Ck[ik] @ U[ik]
        ph = np.exp(-1j * bs[ib] * tau)
        out.append(ci.conj().T @ (ph[:, None] * cj))
    return out


def spread(U):
    Ms = [Mmat(U, ib) for ib in range(2)]
    nw = U[0].shape[1]
    rbar = np.zeros(nw)
    for ib in range(2):
        for M in Ms[ib]:
            rbar += -wb * bs[ib] * np.angle(np.diag(M)) / Nk
    OmI = OmD = OmOD = 0.0
    for ib in range(2):
        for M in Ms[ib]:
            OmI += wb * (nw - np.sum(np.abs(M) ** 2)) / Nk
            OmOD += wb * (np.sum(np.abs(M) ** 2) - np.sum(np.abs(np.diag(M)) ** 2)) / Nk
            q = np.angle(np.diag(M)) + bs[ib] * rbar
            OmD += wb * np.sum(q ** 2) / Nk
    return OmI, OmD, OmOD, Ms, rbar


def rand_unitary(n):
    X = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return np.linalg.qr(X)[0]


# one occupied band: Omega_I gauge invariant
U1 = [np.array([[1.0 + 0j]]) for _ in kk]
U1r = [np.array([[np.exp(2j * np.pi * rng.random())]]) for _ in kk]
Ck_full = Ck
Ck = [c[:, :1] for c in Ck_full]
a1 = spread(U1)
a2 = spread(U1r)
check("1 band: Omega_I gauge invariant", abs(a1[0] - a2[0]) < 1e-12, f"Omega_I={a1[0]:.6f}")
# two bands (complete set): steepest descent drives Omega -> 0, centers -> tau
Ck = Ck_full
U = [rand_unitary(2) for _ in kk]
O0 = sum(spread(U)[:3])
alpha = 0.5
wtot = 2 * wb
for it in range(3000):
    OmI, OmD, OmOD, Ms, rbar = spread(U)
    G = [np.zeros((2, 2), complex) for _ in kk]
    for ib in range(2):
        for ik, M in enumerate(Ms[ib]):
            Rm = M * np.conj(np.diag(M))[None, :]
            T = (M / np.diag(M)[None, :]) * (np.angle(np.diag(M)) + bs[ib] * rbar)[None, :]
            AR = 0.5 * (Rm - Rm.conj().T)
            ST = (T + T.conj().T) / 2j
            G[ik] += 4 * wb * (AR - ST)
    U = [U[ik] @ expm(alpha / (4 * wtot) * G[ik]) for ik in range(Nk)]
O1 = sum(spread(U)[:3])
OmI, OmD, OmOD, Ms, rbar = spread(U)
check("SD lowers the spread", O1 < O0, f"{O0:.4f} -> {O1:.2e}")
check("complete set: Omega_I = 0", abs(OmI) < 1e-12)
check("complete set: Omega -> 0", O1 < 1e-5, f"Omega={O1:.1e}")
cent = np.sort(np.mod(rbar + 0.25, cellm) - 0.25)   # fold to [-0.25, 0.75)
check("centers -> orbital positions", np.max(np.abs(cent - np.sort(tau))) < 1e-3,
      f"{cent}")

# ---------------------------------------------------------------------------
# 4. orbital optimization gradient (finite cluster, fixed Hamiltonian)
# ---------------------------------------------------------------------------
print("== 4. orbital optimization gradient ==")
nat, npr = 3, 3           # atoms, primitive functions per atom (one l, m-free toy)
Np = nat * npr
X = rng.normal(size=(Np, Np))
Sp = X @ X.T / Np + np.eye(Np)
Sp = Sp / np.sqrt(np.outer(np.diag(Sp), np.diag(Sp)))
Hp = rng.normal(size=(Np, Np))
Hp = (Hp + Hp.T) / 2
ncont = 2                 # contracted functions per atom
nocc = 2


def contraction(a):
    """a[I, alpha, p] -> block matrix (Np x nat*ncont)."""
    Am = np.zeros((Np, nat * ncont))
    for I in range(nat):
        for al in range(ncont):
            Am[I * npr:(I + 1) * npr, I * ncont + al] = a[I, al]
    return Am


def band_energy(a):
    Am = contraction(a)
    e, C = eigh(Am.T @ Hp @ Am, Am.T @ Sp @ Am)
    return 2 * e[:nocc].sum(), e, C


a0 = rng.normal(size=(nat, ncont, npr))
E0, e, C = band_energy(a0)
rho = C[:, :nocc] @ C[:, :nocc].T                 # per spin
edm = C[:, :nocc] @ np.diag(e[:nocc]) @ C[:, :nocc].T
Am = contraction(a0)
HA, SA = Hp @ Am, Sp @ Am                        # <chi|H|phi>, <chi|phi>
grad = np.zeros_like(a0)
for I in range(nat):
    for al in range(ncont):
        i = I * ncont + al
        rows = slice(I * npr, (I + 1) * npr)
        # 2 (spin) * 2 * sum_j (rho_ij <chi|H|phi_j> - E_ij <chi|phi_j>)
        grad[I, al] = 4 * (HA[rows] @ rho[i] - SA[rows] @ edm[i])
fd = np.zeros_like(a0)
h = 1e-6
for idx in np.ndindex(a0.shape):
    ap, am = a0.copy(), a0.copy()
    ap[idx] += h
    am[idx] -= h
    fd[idx] = (band_energy(ap)[0] - band_energy(am)[0]) / (2 * h)
check("dE/da analytic = finite difference", np.max(np.abs(grad - fd)) < 1e-6,
      f"dev={np.max(np.abs(grad - fd)):.1e}")

print("ALL OK" if OK else "SOME CHECKS FAILED")
raise SystemExit(0 if OK else 1)
