# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/on-methods.md (DC, Krylov, DC-LNO).
# Model: a 1D ring of atoms with two orbitals per atom, non-orthogonal
# overlap, hopping and overlap up to second neighbours ("overlap range" = 2
# sites), alternating on-site energies so that the chain has a gap.
# Run: uv run --with numpy --with scipy python tools/check_on_methods.py

import numpy as np
import scipy.linalg as sla

rng = np.random.default_rng(7)
NAT, NO = 48, 2            # atoms, orbitals per atom
RNG_OV = 2                 # |i-j| <= 2 : first neighbours (non-zero H, S)
BETA = 1.0 / (300.0 * 8.617251324e-5 / 27.2113845) / 4.0   # ~1200 K in 1/Hartree
DROP = 1.0e-9              # S eigenvalue threshold (rule C)


def check(name, val, tol, less=True):
    ok = (val < tol) if less else (val > tol)
    rel = "<" if less else ">"
    print(f"{'OK ' if ok else 'NG '} {name}: {val:.3e} ({rel} {tol:.1e})")
    assert ok


# ------------------------------------------------------------------ model
def blocks():
    h = {0: None}
    s = {0: np.eye(NO)}
    h[1] = np.array([[-0.10, 0.06], [-0.06, 0.08]])
    h[2] = np.array([[-0.02, 0.01], [-0.01, 0.01]])
    s[1] = np.array([[0.12, -0.05], [0.05, -0.03]])
    s[2] = np.array([[0.01, -0.005], [0.005, -0.002]])
    return h, s


hb, sb = blocks()
N = NAT * NO
H = np.zeros((N, N))
S = np.zeros((N, N))
for i in range(NAT):
    eps = np.array([-0.40, 0.35]) + (0.05 if i % 2 else -0.05)
    H[i*NO:(i+1)*NO, i*NO:(i+1)*NO] = np.diag(eps)
    S[i*NO:(i+1)*NO, i*NO:(i+1)*NO] = sb[0]
    for d in (1, 2):
        j = (i + d) % NAT              # ring: no edge states
        if True:
            H[i*NO:(i+1)*NO, j*NO:(j+1)*NO] = hb[d]
            H[j*NO:(j+1)*NO, i*NO:(i+1)*NO] = hb[d].T
            S[i*NO:(i+1)*NO, j*NO:(j+1)*NO] = sb[d]
            S[j*NO:(j+1)*NO, i*NO:(i+1)*NO] = sb[d].T
NEL = float(NAT * 2)       # 2 electrons per atom (spin-degenerate): half of the orbitals filled
assert np.linalg.eigvalsh(S).min() > 0.5


def orb(i):
    return list(range(i*NO, (i+1)*NO))


def fermi(e, mu):
    x = np.clip((e - mu) * BETA, -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(x))


def bisect(count, lo=-10.0, hi=10.0, tol=1e-12):
    for _ in range(1000):
        mu = 0.5 * (lo + hi)
        d = NEL - count(mu)
        if abs(d) < tol:
            break
        if d >= 0:
            lo = mu
        else:
            hi = mu
    return mu


# ------------------------------------------------------------- exact
ev, C = sla.eigh(H, S)
mu_ex = bisect(lambda m: 2 * fermi(ev, m).sum())
rho_ex = (C * fermi(ev, mu_ex)) @ C.T
edm_ex = (C * (fermi(ev, mu_ex) * ev)) @ C.T
check("exact: 2 Tr(rho S) = Nel", abs(2 * np.trace(rho_ex @ S) - NEL), 1e-9)


def dist(i, j):
    d = abs(i - j) % NAT
    return min(d, NAT - d)


def fnan(i):
    return [j for j in range(NAT) if dist(i, j) <= RNG_OV]


def cluster(i, R):
    """central atom first, then the other atoms sorted by distance."""
    js = sorted([j for j in range(NAT) if dist(i, j) <= R], key=lambda j: (dist(i, j), j))
    return js


def solve_rule_c(Hc, Sc):
    lam, U = np.linalg.eigh(Sc)
    keep = lam >= DROP
    X = U[:, keep] / np.sqrt(lam[keep])
    e, d = np.linalg.eigh(X.T @ Hc @ X)
    return e, X @ d


# ------------------------------------------------------------- DC
def local_problems_dc(R):
    """Per atom: eigenvalues, residues r[mu, a, (j,b)] for the central rows and
    first-neighbour columns, and Mulliken weights kappa[mu]."""
    out = []
    for i in range(NAT):
        at = cluster(i, R)
        idx = sum((orb(j) for j in at), [])
        e, c = solve_rule_c(H[np.ix_(idx, idx)], S[np.ix_(idx, idx)])
        cen = c[:NO, :]                                   # central atom rows
        cols = [k for k, j in enumerate(at) if j in fnan(i)]
        colidx = sum(([k*NO + b for b in range(NO)] for k in cols), [])
        glob = sum((orb(at[k]) for k in cols), [])
        res = np.einsum("am,bm->mab", cen, c[colidx, :])  # (mu, a, jb)
        Sblk = S[np.ix_(orb(i), glob)]                    # S_{i a, j b}
        kappa = np.einsum("mab,ab->m", res, Sblk)
        out.append((i, e, res, kappa, glob))
    return out


def dm_from_local(loc):
    mu = bisect(lambda m: 2 * sum((fermi(e, m) * k).sum() for _, e, _, k, _ in loc))
    rho = np.zeros((N, N))
    edm = np.zeros((N, N))
    for i, e, res, k, glob in loc:
        f = fermi(e, mu)
        rho[np.ix_(orb(i), glob)] = np.einsum("m,mab->ab", f, res)
        edm[np.ix_(orb(i), glob)] = np.einsum("m,mab->ab", f * e, res)
    return mu, rho, edm


mask = np.zeros((N, N), bool)
for i in range(NAT):
    for j in fnan(i):
        mask[np.ix_(orb(i), orb(j))] = True

print("-- DC")
mu_full, rho_full, edm_full = dm_from_local(local_problems_dc(NAT))
check("DC with whole system as cluster == exact rho", np.abs(rho_full - rho_ex)[mask].max(), 1e-10)
check("DC with whole system as cluster == exact EDM", np.abs(edm_full - edm_ex)[mask].max(), 1e-10)
errs = []
for R in (3, 5, 7, 9, 11):
    mu, rho, edm = dm_from_local(local_problems_dc(R))
    ne = 2 * np.sum(rho[mask] * S.T[mask])              # sum_i sum_j rho_ij S_ji
    check(f"DC R={R}: electron count with common mu", abs(ne - NEL), 1e-9)
    errs.append(np.abs(rho - rho_ex)[mask].max())
    print(f"    R={R:2d}  max|rho_DC - rho_exact| = {errs[-1]:.2e}  mu-mu_ex = {mu - mu_ex:+.2e}")
check("DC error decays with R (gapped chain): err(11)/err(3)", errs[-1] / errs[0], 1e-6)
rho7 = dm_from_local(local_problems_dc(3))[1]
asym = max(np.abs(rho7[np.ix_(orb(i), orb(j))] - rho7[np.ix_(orb(j), orb(i))].T).max()
           for i in range(NAT) for j in fnan(i))
check("DC density matrix is not symmetric (rho_ij != rho_ji^T)", asym, 1e-8, less=False)


# ------------------------------------------------------------- Krylov
def krylov_basis(Hc, Sc, W0, nblk, invS):
    """S-orthonormal block Krylov basis of span{W0, A W0, ..., A^(nblk-1) W0}, A = invS Hc."""
    def s_orth(V):
        G = V.T @ Sc @ V
        lam, X = np.linalg.eigh(G)
        il = np.where(lam > DROP, 1.0 / np.sqrt(np.maximum(np.abs(lam), DROP)), 0.0)
        return V @ (X * il)
    blocks_ = [s_orth(W0)]
    for _ in range(nblk - 1):
        V = invS @ (Hc @ blocks_[-1])
        for Ub in blocks_:                    # classical block Gram-Schmidt in the S metric
            V = V - Ub @ (Ub.T @ Sc @ V)
        blocks_.append(s_orth(V))
    U = np.hstack(blocks_)
    lam, X = np.linalg.eigh(U.T @ Sc @ U)       # final re-orthonormalisation
    il = np.where(lam > DROP, 1.0 / np.sqrt(np.maximum(np.abs(lam), DROP)), 0.0)
    return U @ (X * il), int(np.sum(lam <= DROP))


print("-- Krylov")
i0, R = NAT // 2, 9
at = cluster(i0, R)
idx = sum((orb(j) for j in at), [])
Hc, Sc = H[np.ix_(idx, idx)], S[np.ix_(idx, idx)]
ncore = 3                                       # central + two nearest neighbours
W0 = np.zeros((len(idx), ncore * NO))
W0[: ncore * NO, :] = np.eye(ncore * NO)
invS = np.linalg.inv(Sc)
q = 3
U, nz = krylov_basis(Hc, Sc, W0, q + 1, invS)
live = np.abs(U).sum(axis=0) > 0
check("Krylov basis is S-orthonormal", np.abs(U[:, live].T @ Sc @ U[:, live] - np.eye(live.sum())).max(), 1e-11)
# core/buffer decomposition of the projected Hamiltonian (core = first-neighbour atoms)
nc = sum(NO for j in at if j in fnan(i0))
uc, ub = U[:nc], U[nc:]
HK = uc.T @ Hc[:nc, :nc] @ uc + uc.T @ Hc[:nc, nc:] @ ub + ub.T @ Hc[nc:, :nc] @ uc + ub.T @ Hc[nc:, nc:] @ ub
check("core/buffer sum == U^T H U", np.abs(HK - U.T @ Hc @ U).max(), 1e-12)
# moments: W0^T S A^p W0 reproduced up to p = 2q+1
A = invS @ Hc
t = U.T @ Sc @ W0
HKl = (U.T @ Hc @ U)[np.ix_(live, live)]
tl = t[live]
mom_err = []
for p in range(0, 2 * q + 3):
    ex = W0.T @ Sc @ np.linalg.matrix_power(A, p) @ W0
    ap = tl.T @ np.linalg.matrix_power(HKl, p) @ tl
    mom_err.append(np.abs(ex - ap).max() / max(1.0, np.abs(ex).max()))
check(f"moments p<=2q+1={2*q+1} reproduced", max(mom_err[: 2 * q + 2]), 1e-10)
check(f"moment p=2q+2={2*q+2} not reproduced (error ratio to p<=2q+1)",
      mom_err[2 * q + 2] / max(mom_err[: 2 * q + 2]), 1e3, less=False)
# full Krylov dimension -> identical to DC on the same cluster
Ufull, _ = krylov_basis(Hc, Sc, W0, 40, invS)
live = np.abs(Ufull).sum(axis=0) > 0
eK = np.linalg.eigvalsh(Ufull[:, live].T @ Hc @ Ufull[:, live])
eD, _ = solve_rule_c(Hc, Sc)
check("Krylov with full dimension: eigenvalues == DC cluster", np.abs(np.sort(eK)[: len(eD)] - eD).max(), 1e-9)
# approximate S^{-1} from a Krylov space of S (Euclidean orthonormal basis)
V = [np.linalg.qr(W0)[0]]
for _ in range(60):
    Y = Sc @ V[-1]
    for Vb in V:
        Y = Y - Vb @ (Vb.T @ Y)
    Qb, Rb = np.linalg.qr(Y)
    if np.abs(Rb).max() < 1e-12:
        break
    V.append(Qb)
Vm = np.hstack(V)
Q = Vm @ np.linalg.pinv(Vm.T @ Sc @ Vm, rcond=1e-14) @ Vm.T
check("Krylov S^{-1} with full dimension == S^{-1} on the start block", np.abs(Q @ Sc @ W0 - W0).max(), 1e-9)

# ------------------------------------------------------------- LNO
print("-- DC-LNO")
Lam = [sum(rho_ex[np.ix_(orb(i), orb(j))] @ S[np.ix_(orb(j), orb(i))] for j in fnan(i)) for i in range(NAT)]
check("2 sum_i tr Lambda_i = Nel", abs(2 * sum(np.trace(L) for L in Lam) - NEL), 1e-9)
lno = []
for L in Lam:
    w, v = np.linalg.eig(L)
    o = np.argsort(-w.real)
    v = v[:, o].real
    v /= np.linalg.norm(v, axis=0)            # Euclidean normalisation
    lno.append((w.real[o], v))
pops = np.array([p for p, _ in lno])
print(f"    LNO populations (atom 0,1): {pops[0]}, {pops[1]}")
check("Lambda eigenvalues are real here", max(np.abs(np.linalg.eigvals(L).imag).max() for L in Lam), 1e-12)


def dclno_local(R, Rpao, nlno):
    out = []
    for i in range(NAT):
        at = cluster(i, R)
        T_blocks, idx = [], []
        for j in at:
            idx += orb(j)
            T_blocks.append(np.eye(NO) if dist(i, j) <= Rpao else lno[j][1][:, :nlno])
        T = sla.block_diag(*T_blocks)
        Hc = T.T @ H[np.ix_(idx, idx)] @ T
        Sc = T.T @ S[np.ix_(idx, idx)] @ T
        lam, Us = np.linalg.eigh(Sc)
        X = Us / np.sqrt(np.abs(lam))           # rule D: |lambda|^{-1/2}, no threshold
        e, d = np.linalg.eigh(X.T @ Hc @ X)
        c = X @ d                               # coefficients in the (PAO, LNO) basis
        cols = [k for k, j in enumerate(at) if j in fnan(i)]
        colidx = sum(([k*NO + b for b in range(NO)] for k in cols), [])  # PAO atoms come first
        glob = sum((orb(at[k]) for k in cols), [])
        res = np.einsum("am,bm->mab", c[:NO], c[colidx])
        kappa = np.einsum("mab,ab->m", res, S[np.ix_(orb(i), glob)])
        out.append((i, e, res, kappa, glob))
    return out


R, Rpao = 7, 3
rho_dc = dm_from_local(local_problems_dc(R))[1]
rho_l2 = dm_from_local(dclno_local(R, Rpao, 2))[1]
check("DC-LNO with all LNOs == DC", np.abs(rho_l2 - rho_dc)[mask].max(), 1e-10)
rho_l1 = dm_from_local(dclno_local(R, Rpao, 1))[1]
e_l1 = np.abs(rho_l1 - rho_ex)[mask].max()
e_dc = np.abs(rho_dc - rho_ex)[mask].max()
e_small = np.abs(dm_from_local(local_problems_dc(Rpao))[1] - rho_ex)[mask].max()
print(f"    err DC(R={Rpao}) = {e_small:.2e}, DC-LNO(1 LNO) = {e_l1:.2e}, DC(R={R}) = {e_dc:.2e}")
# informational only: in this simple model one LNO per outer atom changes the error little

# ------------------------------------------------------------- buffer rule of DC-LNO
def n_pao_snan(dists, r0, rhop, b):
    """number of second-neighbour atoms that keep PAOs (sorted distances)."""
    sf = 0.1
    for _ in range(400):
        rc = r0 + sf * (rhop - r0)
        cnt = int(np.sum(np.asarray(dists) <= rc))
        if int(len(dists) * b) < cnt:
            return cnt
        sf *= 1.05
    return cnt


d = np.sort(rng.uniform(3.0, 8.0, 40))
check("buffer=1.0 keeps PAOs on all second neighbours", abs(n_pao_snan(d, 2.9, 8.0, 1.0) - 40), 0.5)
check("buffer=0.0 still keeps PAOs on >=1 second neighbour", n_pao_snan(d, 2.9, 8.0, 0.0), 0.5, less=False)
n02 = n_pao_snan(d, 2.9, 8.0, 0.2)
check("buffer=0.2 keeps just over 20% of second neighbours", abs(n02 - 9) , 1.5)
print("all checks passed")
