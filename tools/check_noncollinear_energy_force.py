# SPDX-License-Identifier: MIT
"""Numerical checks for src/modules/noncollinear-energy-force.md.

All conventions are the ones written in the book:
  * real-space blocks O^{ss'}_{ab}(L) = <a s, 0 | O | b s', L>,
    Bloch sum O(k) = sum_L exp(+i k L) O(L), spinor index p = s N + a;
  * density matrix rho^{ss'}_{ab}(L) = sum_k w_k sum_n f_n exp(+i k L)
    conj(c_{s a, n}) c_{s' b, n};
  * grid spin density n_{ss'}(r) = sum_n f_n psi_{s}(r) conj(psi_{s'}(r)),
    m = (2 Re n12, -2 Im n12, n11 - n22).

Checks
  1. Band energy == sum_{ss'} sum_{L,ab} rho^{ss'} H^{ss'} (complex, SO-like H),
     and its real/imag decomposition (Re Re - Im Im, factor 2 on up-down).
  2. Pulay identity for the 2N x 2N spinor problem with S^(2) = I_2 (x) S:
     d(sum f e) = Tr[rho dH] - Tr[E dS], and only Re(E_uu + E_dd) is needed.
  3. Grid identity: sum_r Tr[V(r) n(r)] == sum rho^{ss'} H_grid^{ss'}, with
     n_{ss'}(r) = sum rho^{s's}_{ab}(L) phi_a phi_b  (index swap!).
  4. Non-collinear XC: V^{ss'} = dE/dn_{s's} for E = sum f(n'_up, n'_dn)
     (local-frame model functional), by finite differences on Hermitian n.
  5. Global SU(2) rotation without SO: eigenvalues, Re(rho_uu + rho_dd) and the
     local-frame densities are invariant; atomic moments rotate as vectors.
  6. Time reversal with SO: for a time-reversal-symmetric H with complete Kramers
     pairs, rho^{ud}_{ab}(L) = -rho^{ud}_{ba}(-L); the Hermitian-part Mulliken
     moments vanish, while the row-only definition (existing implementation)
     can give spurious per-atom moments in low symmetry.
  7. Constraint target N0 = 1/2 (N sigma0 + dN n0.sigma) equals the rotated
     diag(N_up, N_dn) and Ucs = v ||N - N0||_F^2 vanishes when aligned.
  8. Collinear reduction: if H^{ud} = 0 and H^{uu}, H^{dd} real, the spinor
     band energy and density equal two collinear problems.

Run: uv run --with numpy --with scipy python tools/check_noncollinear_energy_force.py
"""
import numpy as np
from scipy.linalg import eigh, expm

rng = np.random.default_rng(7)
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]])
SZ = np.array([[1, 0], [0, -1]], complex)
LS = [-1, 0, 1]          # lattice vectors of a 1D chain
NK = 12                  # exact for |L| <= 1 products up to |L| <= 2
KS = 2 * np.pi * np.arange(NK) / NK
W = np.full(NK, 1.0 / NK)


def herm_blocks(n, spin_offdiag=True, complex_diag=True, scale=1.0):
    """Random spinor blocks H[L] (2n x 2n) with H[-L] = H[L]^dagger."""
    H = {}
    a = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n)) * complex_diag
    h0 = 0.5 * (a + a.conj().T)
    h1 = (rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n)) * complex_diag) * 0.4
    for h in (h0, h1):
        if not spin_offdiag:
            h[:n, n:] = 0
            h[n:, :n] = 0
    H[0] = scale * h0
    H[1] = scale * h1
    H[-1] = scale * h1.conj().T
    return H


def real_overlap(n):
    b = rng.normal(size=(n, n)) * 0.1
    s0 = np.eye(n) + 0.5 * (b + b.T)
    s1 = rng.normal(size=(n, n)) * 0.05
    return {0: s0, 1: s1, -1: s1.T}


def bloch(O, k):
    return sum(np.exp(1j * k * L) * O[L] for L in O)


def solve(H, S, n, nocc):
    """Fill the lowest nocc spinor bands at every k. Returns rho, E (2n x 2n per L), Eband."""
    rho = {L: np.zeros((2 * n, 2 * n), complex) for L in range(-2, 3)}
    edm = {L: np.zeros((2 * n, 2 * n), complex) for L in range(-2, 3)}
    eband = 0.0
    for k, w in zip(KS, W):
        hk = bloch(H, k)
        sk = np.kron(np.eye(2), bloch(S, k))
        e, c = eigh(hk, sk)
        c = c[:, :nocc]
        e = e[:nocc]
        eband += w * e.sum()
        P = c.conj() @ c.T                      # P[p,q] = sum conj(c_p) c_q
        Pe = (c.conj() * e) @ c.T
        for L in rho:
            rho[L] += w * np.exp(1j * k * L) * P
            edm[L] += w * np.exp(1j * k * L) * Pe
    return rho, edm, eband


def spin_block(M, n, s, t):
    return M[s * n:(s + 1) * n, t * n:(t + 1) * n]


def contract(rho, H):
    return sum(np.sum(rho[L] * H[L]) for L in H)


# ------------------------------------------------------------------ 1
def check_energy_contraction():
    n, nocc = 4, 3
    H = herm_blocks(n)
    S = real_overlap(n)
    rho, _, eb = solve(H, S, n, nocc)
    e = contract(rho, H)
    # real/imag decomposition with the up-down factor 2
    dec = 0.0
    for L in H:
        r, h = rho[L], H[L]
        for s in (0, 1):
            rr, hh = spin_block(r, n, s, s), spin_block(h, n, s, s)
            dec += np.sum(rr.real * hh.real - rr.imag * hh.imag)
        rr, hh = spin_block(r, n, 0, 1), spin_block(h, n, 0, 1)
        dec += 2 * np.sum(rr.real * hh.real - rr.imag * hh.imag)
    print(f"[1] Eband={eb:.12f} sum rho H={e.real:.12f} (Im {e.imag:.1e}) decomposed={dec:.12f}")
    assert abs(eb - e.real) < 1e-10 and abs(e.imag) < 1e-10 and abs(dec - eb) < 1e-10
    # Hermiticity of rho: rho^{s's}_{ba}(-L) = conj rho^{ss'}_{ab}(L)
    err = max(np.abs(rho[-L] - rho[L].conj().T).max() for L in (0, 1, 2))
    print(f"    Hermiticity of rho(L): {err:.1e}")
    assert err < 1e-12


# ------------------------------------------------------------------ 2
def check_pulay():
    n, nocc = 4, 3
    H0, H1 = herm_blocks(n), herm_blocks(n, scale=0.3)
    S0, S1 = real_overlap(n), real_overlap(n)
    S1 = {L: S1[L] - (np.eye(n) if L == 0 else 0) for L in S1}  # perturbation only

    def at(x):
        H = {L: H0[L] + x * H1[L] for L in H0}
        S = {L: S0[L] + x * S1[L] for L in S0}
        return solve(H, S, n, nocc)

    x0, dx = 0.05, 1e-5
    fd = (at(x0 + dx)[2] - at(x0 - dx)[2]) / (2 * dx)
    rho, edm, _ = at(x0)
    an_full = contract(rho, H1).real - contract(edm, {L: np.kron(np.eye(2), S1[L]) for L in S1}).real
    an_diag = contract(rho, H1).real - sum(
        np.sum((spin_block(edm[L], n, 0, 0) + spin_block(edm[L], n, 1, 1)).real * S1[L]) for L in S1)
    print(f"[2] FD={fd:.10f} Tr[rho dH]-Tr[E dS]={an_full:.10f} with Re(E_uu+E_dd) only={an_diag:.10f}")
    assert abs(fd - an_full) < 1e-7 and abs(an_full - an_diag) < 1e-12


# ------------------------------------------------------------------ 3
def check_grid_identity():
    n, nocc, npt = 3, 2, 50
    H = herm_blocks(n)
    S = real_overlap(n)
    rho, _, _ = solve(H, S, n, nocc)
    Lset = range(-2, 3)
    cen = rng.uniform(0, 1, size=n)          # orbital centres in the unit cell (length 1)
    wid = rng.uniform(4, 9, size=n)
    xs = (np.arange(npt) + 0.5) / npt        # grid points of the cell
    Rs = range(-4, 5)

    def phi(x):  # phi_a(x) for all a, x array -> (len(x), n)
        return np.exp(-wid * (x[:, None] - cen) ** 2) * (1 + 0.3 * (x[:, None] - cen))

    v = rng.normal(size=(npt, 2, 2)) + 1j * rng.normal(size=(npt, 2, 2))
    v = 0.5 * (v + np.conj(np.transpose(v, (0, 2, 1))))
    # n_{ss'}(r) = sum_R sum_L sum_ab rho^{s's}_{ab}(L) phi_a(r-R) phi_b(r-R-L)
    nr = np.zeros((npt, 2, 2), complex)
    hg = {L: np.zeros((2, 2, n, n), complex) for L in Lset}
    for R in Rs:
        pa = phi(xs - R)
        for L in Lset:
            pb = phi(xs - R - L)
            for s in (0, 1):
                for t in (0, 1):
                    nr[:, s, t] += np.einsum('ra,ab,rb->r', pa, spin_block(rho[L], n, t, s), pb)
                    hg[L][s, t] += np.einsum('ra,r,rb->ab', pa, v[:, s, t], pb)
    herm = np.abs(nr - np.conj(np.transpose(nr, (0, 2, 1)))).max()
    e_grid = np.einsum('rst,rts->', v, nr)
    e_mat = sum(np.sum(spin_block(rho[L], n, s, t) * hg[L][s, t]) for L in Lset for s in (0, 1) for t in (0, 1))
    print(f"[3] n(r) Hermitian {herm:.1e}; sum Tr[V n]={e_grid.real:.12f} sum rho Hgrid={e_mat.real:.12f}")
    assert herm < 1e-12 and abs(e_grid - e_mat) < 1e-10


# ------------------------------------------------------------------ 4
def local_frame(nm):
    ntot = np.trace(nm).real
    m = np.array([2 * nm[0, 1].real, -2 * nm[0, 1].imag, (nm[0, 0] - nm[1, 1]).real])
    am = np.linalg.norm(m)
    return 0.5 * (ntot + am), 0.5 * (ntot - am), m / am


def check_xc_potential():
    def f(nu, nd):  # model local-frame energy density (LSDA-like)
        return -0.7 * (nu ** (4 / 3) + nd ** (4 / 3)) + 0.05 * nu * nd

    def dfu(nu, nd):
        return -0.7 * 4 / 3 * nu ** (1 / 3) + 0.05 * nd

    def dfd(nu, nd):
        return -0.7 * 4 / 3 * nd ** (1 / 3) + 0.05 * nu

    for _ in range(5):
        n0 = rng.uniform(0.5, 1.0)
        m = rng.normal(size=3) * 0.3
        nm = 0.5 * (n0 * np.eye(2) + m[0] * SX + m[1] * SY + m[2] * SZ)
        # sigma_y[0,1] = -i, so n12 = (m_x - i m_y)/2: the density convention
        nu, nd, mh = local_frame(nm)
        assert np.allclose(mh, m / np.linalg.norm(m))
        vb, dv = 0.5 * (dfu(nu, nd) + dfd(nu, nd)), 0.5 * (dfu(nu, nd) - dfd(nu, nd))
        th, ph = np.arccos(mh[2]), np.arctan2(mh[1], mh[0])
        V = np.array([[vb + dv * np.cos(th), dv * np.sin(th) * np.exp(-1j * ph)],
                      [dv * np.sin(th) * np.exp(1j * ph), vb - dv * np.cos(th)]])
        # directional derivative along a random Hermitian dn: dE = Tr[V dn]
        d = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        d = 0.5 * (d + d.conj().T)
        h = 1e-6
        ep = f(*local_frame(nm + h * d)[:2])
        em = f(*local_frame(nm - h * d)[:2])
        fd = (ep - em) / (2 * h)
        an = np.trace(V @ d).real
        assert abs(fd - an) < 1e-7, (fd, an)
    print(f"[4] Tr[V_xc dn] == dE_xc (local frame model), last FD={fd:.10f} analytic={an:.10f}")


# ------------------------------------------------------------------ 5
def su2(axis, angle):
    axis = axis / np.linalg.norm(axis)
    return expm(-0.5j * angle * (axis[0] * SX + axis[1] * SY + axis[2] * SZ))


def atomic_moments(rho, S, n, atoms, hermitian_part=False):
    """Per-atom 2x2 Mulliken matrix in the density convention, returns list of (N, m).

    Default: the existing implementation's definition: A^{st} = sum_{i in I} sum rho^{st}_{i,Jj} S_{i,Jj},
    N^{ss} = Re A^{ss}, N^{ud} = conj(A^{ud}), N^{du} = A^{ud}.
    hermitian_part=True: N = conj((A + A^dagger)/2) (rotation covariant even with SO).
    """
    out = []
    for idx in atoms:
        A = np.zeros((2, 2), complex)
        for L in S:
            for s in (0, 1):
                for t in (0, 1):
                    blk = spin_block(rho[L], n, s, t)
                    A[s, t] += np.sum(blk[idx, :] * S[L][idx, :])
        if hermitian_part:
            Nm = (0.5 * (A + A.conj().T)).conj()
        else:
            Nm = np.array([[A[0, 0].real, np.conj(A[0, 1])], [A[0, 1], A[1, 1].real]])
        m = np.array([2 * Nm[0, 1].real, -2 * Nm[0, 1].imag, (Nm[0, 0] - Nm[1, 1]).real])
        out.append((np.trace(Nm).real, m))
    return out


def check_spin_rotation():
    n, nocc = 4, 3
    # spin-dependent but SO-free H: H = h0 (x) 1 + sum_i b_i (x) sigma_i with real symmetric blocks
    S = real_overlap(n)

    def realsym(scale):
        a = rng.normal(size=(n, n)) * scale
        return {0: a + a.T, 1: 0.3 * rng.normal(size=(n, n)) * scale}

    h0 = realsym(1.0)
    bs = [realsym(0.3) for _ in range(3)]
    for d in [h0] + bs:
        d[-1] = d[1].T
    # with only real-symmetric scalar parts and sigma couplings, this is SO-free only if
    # the b_i share one spatial structure; use b_i = c_i * b0 (collinear-type exchange field)
    b0 = bs[0]
    cvec = rng.normal(size=3)

    def build(cv):
        return {L: np.kron(np.eye(2), h0[L]) + sum(np.kron(P, c * b0[L]) for P, c in zip((SX, SY, SZ), cv))
                for L in h0}

    H = build(cvec)
    U = su2(rng.normal(size=3), 1.234)
    R = np.array([[0.5 * np.trace(Pj @ U @ Pi @ U.conj().T).real for Pi in (SX, SY, SZ)] for Pj in (SX, SY, SZ)])
    H2 = {L: np.kron(U, np.eye(n)) @ H[L] @ np.kron(U, np.eye(n)).conj().T for L in H}
    rho1, _, eb1 = solve(H, S, n, nocc)
    rho2, _, eb2 = solve(H2, S, n, nocc)
    sc1 = sum(np.sum((spin_block(rho1[L], n, 0, 0) + spin_block(rho1[L], n, 1, 1)).real * h0[L]) for L in h0)
    sc2 = sum(np.sum((spin_block(rho2[L], n, 0, 0) + spin_block(rho2[L], n, 1, 1)).real * h0[L]) for L in h0)
    m1 = atomic_moments(rho1, S, n, [[0, 1], [2, 3]])
    m2 = atomic_moments(rho2, S, n, [[0, 1], [2, 3]])
    # psi -> U psi, so n = (N + m.sigma)/2 -> U n U^dagger and m_j -> sum_i R_ji m_i
    # with R_ji = Tr(sigma_j U sigma_i U^dagger)/2
    errs = [np.abs(R @ a[1] - b[1]).max() for a, b in zip(m1, m2)]
    print(f"[5] SU(2) rotation: dEband={abs(eb1-eb2):.1e} dTr[Re(rho_uu+rho_dd) h0]={abs(sc1-sc2):.1e}"
          f" moment rotation err={max(errs):.1e} |m|={[round(float(np.linalg.norm(a[1])), 6) for a in m1]}")
    assert abs(eb1 - eb2) < 1e-10 and abs(sc1 - sc2) < 1e-10 and max(errs) < 1e-10
    # with an SO-like (generic complex) H the row-only definition is not covariant,
    # while the Hermitian part is; record the size of the difference
    Hg = herm_blocks(n)
    Hg2 = {L: np.kron(U, np.eye(n)) @ Hg[L] @ np.kron(U, np.eye(n)).conj().T for L in Hg}
    g1, g2 = solve(Hg, S, n, nocc)[0], solve(Hg2, S, n, nocc)[0]
    eh = max(np.abs(R @ a[1] - b[1]).max() for a, b in zip(atomic_moments(g1, S, n, [[0, 1], [2, 3]], True),
                                                         atomic_moments(g2, S, n, [[0, 1], [2, 3]], True)))
    ec = max(np.abs(R @ a[1] - b[1]).max() for a, b in zip(atomic_moments(g1, S, n, [[0, 1], [2, 3]]),
                                                         atomic_moments(g2, S, n, [[0, 1], [2, 3]])))
    print(f"    generic complex H: covariance error Hermitian-part def={eh:.1e}, row-only def={ec:.1e}")
    assert eh < 1e-10


# ------------------------------------------------------------------ 6
def check_time_reversal():
    n = 4
    S = real_overlap(n)
    H = herm_blocks(n)
    # project onto time-reversal-symmetric part: T H(L) T^-1 with T = (i sigma_y (x) 1) K
    Ty = np.kron(1j * SY, np.eye(n))
    Hs = {L: 0.5 * (H[L] + Ty @ H[L].conj() @ Ty.conj().T) for L in H}
    rho, _, _ = solve(Hs, S, n, 2 * 2)  # even number of bands: Kramers pairs complete at TRIM; check moments
    ms = atomic_moments(rho, S, n, [[0, 1], [2, 3]], hermitian_part=True)
    mx = max(np.abs(m).max() for _, m in ms)
    mr = max(np.abs(m).max() for _, m in atomic_moments(rho, S, n, [[0, 1], [2, 3]]))
    # time reversal makes rho^{ud}_{ab}(L) antisymmetric in the orbital pair
    anti = max(np.abs(spin_block(rho[L], n, 0, 1) + spin_block(rho[-L], n, 0, 1).T).max() for L in (0, 1))
    imtr = sum(np.sum((spin_block(rho[L], n, 0, 0) + spin_block(rho[L], n, 1, 1)).imag * S[L]) for L in S)
    print(f"[6] time reversal (SO-like H): rho_ud antisymmetric err={anti:.1e}; max |m_I| Hermitian-part"
          f" def={mx:.1e}; row-only def={mr:.1e} (spurious); sum Im(rho_uu+rho_dd) S={abs(imtr):.1e}")
    assert anti < 1e-10 and mx < 1e-8 and abs(imtr) < 1e-10


# ------------------------------------------------------------------ 7
def check_constraint_target():
    for _ in range(5):
        th, ph = rng.uniform(0, np.pi), rng.uniform(-np.pi / 2, 3 * np.pi / 2)
        nu, nd = rng.uniform(3, 6), rng.uniform(0, 3)
        n0 = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
        # rotation form: N0 = sum_s N_s conj(U_s i) U_s j
        U = np.array([[np.exp(0.5j * ph) * np.cos(th / 2), np.exp(-0.5j * ph) * np.sin(th / 2)],
                      [-np.exp(0.5j * ph) * np.sin(th / 2), np.exp(-0.5j * ph) * np.cos(th / 2)]])
        N0 = nu * np.outer(U[0].conj(), U[0]) + nd * np.outer(U[1].conj(), U[1])
        ref = (0.5 * ((nu + nd) * np.eye(2) + (nu - nd) * (n0[0] * SX + n0[1] * SY + n0[2] * SZ)))
        assert np.allclose(N0, ref)
        _, _, mh = local_frame(N0)
        assert np.allclose(mh, n0)
    print("[7] constraint target N0 = (N s0 + dN n0.s)/2 in the density convention; axis recovered")


# ------------------------------------------------------------------ 8
def check_collinear_reduction():
    n, nocc_up, nocc_dn = 4, 2, 1
    S = real_overlap(n)

    def realsym():
        a = rng.normal(size=(n, n))
        b = 0.3 * rng.normal(size=(n, n))
        return {0: a + a.T, 1: b, -1: b.T}

    hu, hd = realsym(), realsym()
    # make the two spins well separated so lowest 3 spinor bands = 2 up + 1 down at every k
    hd = {L: hd[L] + (8.0 * S[L] if L in S else 0) for L in hd}
    hu = {L: hu[L] - (8.0 * S[L] if L in S else 0) for L in hu}
    H = {L: np.block([[hu[L], np.zeros((n, n))], [np.zeros((n, n)), hd[L]]]) for L in hu}
    # choose fillings so that the spinor filling is the union (check the eigenvalue ordering)
    ebu = ebd = 0.0
    for k, w in zip(KS, W):
        eu = eigh(bloch(hu, k), bloch(S, k), eigvals_only=True)
        ed = eigh(bloch(hd, k), bloch(S, k), eigvals_only=True)
        ebu += w * eu[:nocc_up].sum()
        allv = np.sort(np.concatenate([eu, ed]))
        assert np.isclose(allv[:nocc_up].sum(), eu[:nocc_up].sum())
    rho, _, eb = solve(H, S, n, nocc_up)
    ud = max(np.abs(spin_block(rho[L], n, 0, 1)).max() for L in rho)
    print(f"[8] collinear reduction: Eband spinor={eb:.12f} collinear up={ebu:.12f} max|rho_ud|={ud:.1e}")
    assert abs(eb - ebu) < 1e-10 and ud < 1e-12


if __name__ == "__main__":
    check_energy_contraction()
    check_pulay()
    check_grid_identity()
    check_xc_potential()
    check_spin_rotation()
    check_time_reversal()
    check_constraint_target()
    check_collinear_reduction()
    print("all checks passed")
