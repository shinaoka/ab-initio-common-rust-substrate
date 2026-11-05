# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/kpoints-bloch.md and
# src/modules/generalized-eigen.md, using a random model of localized
# orbitals (real H(L), S(L) with O(-L) = O(L)^T) on a small 3D lattice.
# Run: uv run --with numpy --with scipy python tools/check_kpoints_eigen.py

import itertools
import numpy as np
import scipy.linalg as sla

rng = np.random.default_rng(1)
N = 6                      # orbitals per cell
LV = [L for L in itertools.product((-1, 0, 1), repeat=3)]


def random_real_model(scale, diag_shift):
    """O(L) real, O(-L) = O(L)^T, O(0) symmetric."""
    O = {}
    for L in LV:
        if L in O:
            continue
        mL = tuple(-x for x in L)
        A = scale * rng.standard_normal((N, N))
        if L == mL:
            A = A + A.T + diag_shift * np.eye(N)
            O[L] = A
        else:
            O[L] = A
            O[mL] = A.T.copy()
    return O


Hr = random_real_model(0.3, 0.0)
Sr = random_real_model(0.02, 0.0)            # small overlaps between different orbitals
Sr[(0, 0, 0)] += np.eye(N)                    # S(0) ~ I  => S(k) well conditioned
tau = rng.uniform(0, 1, (N, 3))               # fractional "atomic positions" of the orbitals


def bloch(O, q, with_tau=False):
    """O(k) = sum_L e^{i 2pi q.n} O(L)  (phase on lattice vector only)."""
    M = np.zeros((N, N), complex)
    for L, A in O.items():
        M += np.exp(2j * np.pi * np.dot(q, L)) * A
    if with_tau:   # alternative gauge: phase includes orbital positions
        D = np.diag(np.exp(2j * np.pi * tau @ q))
        M = D.conj() @ M @ D
    return M


def check(name, val, tol):
    ok = val < tol
    print(f"{'OK ' if ok else 'NG '} {name}: {val:.3e} (tol {tol:.0e})")
    assert ok


# ---------------------------------------------------------------- Bloch sum
q = np.array([0.13, -0.27, 0.41])
Hk, Sk = bloch(Hr, q), bloch(Sr, q)
check("H(k) Hermitian", np.abs(Hk - Hk.conj().T).max(), 1e-13)
check("H(-k) = conj H(k)", np.abs(bloch(Hr, -q) - Hk.conj()).max(), 1e-13)
check("H(k+G) = H(k) (lattice-only phase)", np.abs(bloch(Hr, q + [1, 0, 0]) - Hk).max(), 1e-12)
e1 = sla.eigh(Hk, Sk, eigvals_only=True)
e2 = sla.eigh(bloch(Hr, q, True), bloch(Sr, q, True), eigvals_only=True)
check("gauge change keeps eigenvalues", np.abs(e1 - e2).max(), 1e-11)

# ---------------------------------------------------------------- k meshes


def mesh_regular(Ns):
    ax = [[0.0] if n == 1 else [-0.5 + (2 * l + 1) / (2 * n) for l in range(n)] for n in Ns]
    return ax


def mesh_gcenter(Ns):
    return [[0.0] if n == 1 else [l / n for l in range(n)] for n in Ns]


def reduce_inversion(Ns, kind):
    """Keep one representative of each {k, -k} pair; weight 2 (1 if k = -k mod G)."""
    def partner(l, n):
        if kind == "regular":
            return n - 1 - l
        return 0 if l == 0 else n - l
    seen, out = set(), []
    for idx in itertools.product(*[range(n) for n in Ns]):
        if idx in seen:
            continue
        p = tuple(partner(l, n) for l, n in zip(idx, Ns))
        seen.add(idx)
        seen.add(p)
        out.append((idx, 1 if p == idx else 2))
    return out


for kind, gen in (("regular", mesh_regular), ("gcenter", mesh_gcenter)):
    for Ns in [(2, 3, 4), (3, 3, 3), (4, 4, 1), (1, 1, 1)]:
        ax = gen(Ns)
        red = reduce_inversion(Ns, kind)
        assert sum(w for _, w in red) == np.prod(Ns)
        # band-energy-like sum (any function with f(k) = f(-k)): full vs reduced mesh
        f = lambda k: np.sum(np.minimum(sla.eigh(bloch(Hr, k), bloch(Sr, k), eigvals_only=True), 0.0))
        full = sum(f(np.array(k)) for k in itertools.product(*ax))
        redsum = sum(w * f(np.array([ax[d][idx[d]] for d in range(3)])) for idx, w in red)
        check(f"{kind} {Ns}: reduced mesh sum = full mesh sum (n_red={len(red)})",
              abs(full - redsum), 1e-10)
        n_trim = sum(1 for _, w in red if w == 1)
        print(f"    self-inverse points (weight 1): {n_trim}")

# GCenter mesh == Gamma point of the N1xN2xN3 supercell (folding check)
Ns = (2, 2, 1)
ax = mesh_gcenter(Ns)
ek = np.sort(np.concatenate([sla.eigh(bloch(Hr, np.array(k)), bloch(Sr, np.array(k)), eigvals_only=True)
                             for k in itertools.product(*ax)]))
cells = list(itertools.product(range(Ns[0]), range(Ns[1]), range(Ns[2])))
M = len(cells) * N
Hs, Ss = np.zeros((M, M)), np.zeros((M, M))
for a, ca in enumerate(cells):
    for b, cb in enumerate(cells):
        for L, A in Hr.items():
            d = tuple((cb[i] + L[i] - ca[i]) % Ns[i] for i in range(3))
            if d == (0, 0, 0):
                Hs[a*N:(a+1)*N, b*N:(b+1)*N] += A
                Ss[a*N:(a+1)*N, b*N:(b+1)*N] += Sr[L]
check("GCenter k-mesh = supercell Gamma (folding)", np.abs(np.sort(sla.eigh(Hs, Ss, eigvals_only=True)) - ek).max(), 1e-10)

# density matrix from k-mesh: rho(L) = sum_k w_k sum_n f_n Re[e^{ik.L} c*_i c_j]; check E_band = sum_L tr(rho(L) H(L)^T)
Ns = (3, 2, 2)
ax = mesh_gcenter(Ns)
red = reduce_inversion(Ns, "gcenter")
nocc = 2
rho = {L: np.zeros((N, N)) for L in LV}
Eb = 0.0
for idx, w in red:
    k = np.array([ax[d][idx[d]] for d in range(3)])
    e, C = sla.eigh(bloch(Hr, k), bloch(Sr, k))
    Eb += w * e[:nocc].sum()
    for L in LV:
        ph = np.exp(2j * np.pi * np.dot(k, L))
        rho[L] += w * np.real(ph * (C[:, :nocc].conj() @ C[:, :nocc].T))
tot = np.prod(Ns)
Eb /= tot
EbL = sum(np.sum(rho[L] * Hr[L]) for L in LV) / tot
NeL = sum(np.sum(rho[L] * Sr[L]) for L in LV) / tot
check("E_band from rho(L) (inversion-reduced mesh)", abs(Eb - EbL), 1e-11)
check("N_el from rho(L)", abs(NeL - nocc), 1e-11)

# ------------------------------------------------------ generalized eigen


def solve(H, S, rule, nev=None, cut=1e-10, drop=1e-9):
    lam, U = np.linalg.eigh(S)
    if rule == "A":
        lt = np.where(lam < 0, cut, lam)
    elif rule == "B":
        lt = np.maximum(lam, cut)
    elif rule == "C":
        keep = lam >= drop
        lam, U = lam[keep], U[:, keep]
        lt = lam
    X = U / np.sqrt(lt)
    Hp = X.conj().T @ H @ X
    Hp = (Hp + Hp.conj().T) / 2
    e, D = np.linalg.eigh(Hp)
    if nev is not None:
        e, D = e[:nev], D[:, :nev]
    return e, X @ D, U, lt


def eta(H, S, e, C):
    nH, nS = np.linalg.norm(H, 2), np.linalg.norm(S, 2)
    r = H @ C - S @ C * e
    return (np.linalg.norm(r, axis=0) / ((nH + np.abs(e) * nS) * np.linalg.norm(C, axis=0))).max()


Hk, Sk = bloch(Hr, q), bloch(Sr, q)
ref = sla.eigh(Hk, Sk, eigvals_only=True)          # LAPACK zhegvd-type reference
for rule in "ABC":
    e, C, U, lt = solve(Hk, Sk, rule, nev=4)
    check(f"well-conditioned S, rule {rule}: |eps - zhegv|", np.abs(e - ref[:4]).max(), 1e-12)
    check(f"well-conditioned S, rule {rule}: eta", eta(Hk, Sk, e, C), 1e-13)
    check(f"well-conditioned S, rule {rule}: C^H S C - I", np.abs(C.conj().T @ Sk @ C - np.eye(4)).max(), 1e-12)

# near-singular S: two nearly linearly dependent basis functions
n = 8
B = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
B[:, 1] = B[:, 0] + 1e-7 * rng.standard_normal(n)   # gives lambda_min ~ 1e-14
B[:, 3] = B[:, 2] + 3e-6 * rng.standard_normal(n)   # gives lambda ~ 1e-11
S = B.conj().T @ B
S = S / np.linalg.norm(S, 2)
A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
H = B.conj().T @ (A + A.conj().T) @ B
lam = np.linalg.eigvalsh(S)
print("    near-singular S eigenvalues:", np.array2string(lam[:3], precision=2))
for rule in "ABC":
    e, C, U, lt = solve(H, S, rule)
    if rule in "AB":
        Seff = (U * lt) @ U.conj().T
        check(f"near-singular S, rule {rule}: eta wrt S_eff", eta(H, Seff, e, C), 1e-12)
        kappa = lt.max() / lt.min()
        orth = np.abs(C.conj().T @ Seff @ C - np.eye(n)).max()
        # orthonormality error grows like eps * kappa(S_eff); only meaningful when kappa is moderate
        check(f"near-singular S, rule {rule}: C^H S_eff C - I <= 1e3 eps kappa", orth / (1e3 * 2.2e-16 * kappa), 1.0)
        print(f"    rule {rule}: kappa(S_eff) = {kappa:.1e}, C^H S_eff C - I = {orth:.1e}, "
              f"||S_eff - S|| = {np.linalg.norm(Seff - S, 2):.2e}, max|eps| = {np.abs(e).max():.1e}")
    else:
        P = U / np.sqrt(lt)
        check("near-singular S, rule C: projected residual X^H(HC - SC eps)",
              np.abs(P.conj().T @ (H @ C - S @ C * e)).max() / np.linalg.norm(H, 2), 1e-12)
        check("near-singular S, rule C: C^H S C - I", np.abs(C.conj().T @ S @ C - np.eye(len(e))).max(), 1e-8)
        print(f"    rule C: dim = {len(e)} (dropped {n - len(e)})")
    print(f"    rule {rule}: lowest 3 eps = {np.array2string(e[:3], precision=6)}")

# ------------------------------------------------------ noncollinear layout
# H2 = [[H_uu, H_ud], [H_ud^H, H_dd]] (spin-block-major), S2 = I_2 (x) S
Huu, Hdd = Hk, bloch(Hr, q) + 0.1 * np.eye(N)
Hud = 0.05 * (rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N)))
H2 = np.block([[Huu, Hud], [Hud.conj().T, Hdd]])
S2 = np.kron(np.eye(2), Sk)
lam, U = np.linalg.eigh(Sk)
X = U / np.sqrt(lam)
X2 = np.kron(np.eye(2), X)
Hp = np.block([[X.conj().T @ Huu @ X, X.conj().T @ Hud @ X], [X.conj().T @ Hud.conj().T @ X, X.conj().T @ Hdd @ X]])
check("noncollinear: blockwise X^H H_ss' X = X2^H H2 X2", np.abs(Hp - X2.conj().T @ H2 @ X2).max(), 1e-12)
e, D = np.linalg.eigh(Hp)
check("noncollinear: eigenvalues vs zhegv(H2, S2)", np.abs(e - sla.eigh(H2, S2, eigvals_only=True)).max(), 1e-11)
C = X2 @ D
check("noncollinear: C^H S2 C - I", np.abs(C.conj().T @ S2 @ C - np.eye(2 * N)).max(), 1e-12)
print("all checks passed")
