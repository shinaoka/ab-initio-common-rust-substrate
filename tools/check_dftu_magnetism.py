# SPDX-License-Identifier: MIT
"""Numerical checks for the DFT+U / constraint / Zeeman chapter
(src/modules/dftu-magnetism.md).

Everything here is derived from the formulas written in the chapter; no
OpenMX code is used.

 1. Coulomb tensor U_abcd = <ab|V|cd> in the real harmonics of the PAO chapter
    (built from Slater integrals and Gaunt coefficients): reality, the
    averages U and J, the F^k <- (U, J) relations and rotational invariance.
 2. Energy / potential consistency (v = dE/dn, finite differences) for
    the simplified (Dudarev) form and the general form with sFLL, sAMF,
    cFLL and cAMF, collinear and non-collinear.  sFLL with J=0 equals the
    simplified form with U = F0.
 3. Dual occupation: sum rule, identity Tr(rho H_U) = sum v n, and the
    fixed-rho overlap-derivative force formula (dual and full), with a
    complex non-collinear toy model.
 4. Spin-direction constraint (on / on2) and spin Zeeman term: closed forms
    of the penalty, the effective field h = dE/dm and v = h . sigma.
 5. Orbital-polarization occupation rule conserves the electron number.
 6. Yukawa kernel: Legendre expansion of exp(-lam R)/R and its lam -> 0 limit.
 7. Tesla -> atomic unit.

Run: uv run --with numpy --with scipy --with sympy python tools/check_dftu_magnetism.py
"""
import numpy as np
from scipy.special import sph_harm_y, spherical_in, spherical_kn, eval_legendre
from scipy.spatial.transform import Rotation
from sympy.physics.wigner import wigner_3j

rng = np.random.default_rng(7)
SIG = [np.array([[0, 1], [1, 0]], complex),
       np.array([[0, -1j], [1j, 0]], complex),
       np.array([[1, 0], [0, -1]], complex)]
results = []


def report(name, err, tol):
    ok = err <= tol
    results.append(ok)
    print(f"[{'OK ' if ok else 'NG '}] {name}: {err:.3e} (tol {tol:.0e})")


# ------------------------------------------------------------ real harmonics
def real_from_complex(l):
    """T[a, m+l] with S_a = sum_m T[a, m] Y_lm, in the order of the PAO chapter."""
    T = np.zeros((2 * l + 1, 2 * l + 1), complex)
    s = 1 / np.sqrt(2)
    rows = [(0, 'c')]
    if l == 1:
        rows = [(1, 'c'), (1, 's'), (0, 'c')]   # x, y, z
    elif l == 2:
        rows += [(2, 'c'), (2, 's'), (1, 'c'), (1, 's')]
    else:
        for m in range(1, l + 1):
            rows += [(m, 'c'), (m, 's')]
    for a, (m, t) in enumerate(rows):
        if m == 0:
            T[a, l] = 1
        elif t == 'c':
            T[a, l - m] = s
            T[a, l + m] = (-1) ** m * s
        else:
            T[a, l - m] = 1j * s
            T[a, l + m] = -1j * (-1) ** m * s
    return T


def real_harm(l, theta, phi):
    T = real_from_complex(l)
    Y = np.array([sph_harm_y(l, m, theta, phi) for m in range(-l, l + 1)])
    return (T @ Y).real


def check_real_harmonics_signs():
    # S_{1,x} must be +sqrt(3/4pi) x, d_xz positive at x=z>0, etc.
    th, ph = np.pi / 3, 0.4
    x, y, z = np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)
    p = real_harm(1, th, ph)
    c = np.sqrt(3 / (4 * np.pi))
    err = np.abs(p - c * np.array([x, y, z])).max()
    d = real_harm(2, th, ph)
    ref = np.array([np.sqrt(5 / (16 * np.pi)) * (3 * z * z - 1),
                    np.sqrt(15 / (16 * np.pi)) * (x * x - y * y),
                    np.sqrt(15 / (4 * np.pi)) * x * y,
                    np.sqrt(15 / (4 * np.pi)) * x * z,
                    np.sqrt(15 / (4 * np.pi)) * y * z])
    err = max(err, np.abs(d - ref).max())
    report("real harmonics p,d match the PAO-chapter table", err, 1e-14)


# ------------------------------------------------------------ Coulomb tensor
def gaunt(l, k, m1, m2):
    """c^k(l m1, l m2) = sqrt(4pi/(2k+1)) Int Y*_{l m1} Y_{k,m1-m2} Y_{l m2}."""
    return float((-1) ** m1 * (2 * l + 1) * wigner_3j(l, k, l, 0, 0, 0)
                 * wigner_3j(l, k, l, -m1, m1 - m2, m2))


def coulomb_complex(l, F):
    n = 2 * l + 1
    V = np.zeros((n,) * 4)
    ms = range(-l, l + 1)
    for i1, m1 in enumerate(ms):
        for i2, m2 in enumerate(ms):
            for i3, m3 in enumerate(ms):
                for i4, m4 in enumerate(ms):
                    if m1 + m2 != m3 + m4:
                        continue
                    V[i1, i2, i3, i4] = sum(gaunt(l, k, m1, m3) * gaunt(l, k, m4, m2) * F[k]
                                            for k in range(0, 2 * l + 1, 2))
    return V


def coulomb_real(l, F):
    T = real_from_complex(l)
    V = coulomb_complex(l, F)
    U = np.einsum('am,bn,cp,dq,mnpq->abcd', T.conj(), T.conj(), T, T, V)
    return U


def slater_from_UJ(l, U, J, ratio=0.625):
    F = {0: U}
    if l == 1:
        F[2] = 5 * J
    elif l == 2:
        F[2] = 14 * J / (1 + ratio)
        F[4] = ratio * F[2]
    elif l == 3:
        F[2] = 6435 * J / (286 + 195 * 0.668 + 250 * 0.494)
        F[4] = 0.668 * F[2]
        F[6] = 0.494 * F[2]
    return F


def check_coulomb():
    for l in (1, 2, 3):
        Uin, Jin = 0.25, 0.035
        F = slater_from_UJ(l, Uin, Jin)
        U = coulomb_real(l, F)
        report(f"l={l}: U_abcd real", np.abs(U.imag).max(), 1e-14)
        U = U.real
        n = 2 * l + 1
        Uav = np.einsum('abab->', U) / n**2
        UmJ = (np.einsum('abab->', U) - np.einsum('abba->', U)) / (2 * l * n)
        report(f"l={l}: mean U_abab = F0", abs(Uav - Uin), 1e-14)
        report(f"l={l}: J from averages = input J", abs((Uav - UmJ) - Jin), 1e-14)
        # symmetries <ab|cd> = <ba|dc> = <cd|ab> (real orbitals: also <cb|ad>)
        err = max(np.abs(U - U.transpose(1, 0, 3, 2)).max(), np.abs(U - U.transpose(2, 3, 0, 1)).max(),
                  np.abs(U - U.transpose(2, 1, 0, 3)).max())
        report(f"l={l}: tensor symmetries", err, 1e-14)
        # rotational invariance: D from least squares on sample points
        R = Rotation.from_rotvec([0.3, -0.7, 0.5])
        pts = rng.normal(size=(200, 3))
        pts /= np.linalg.norm(pts, axis=1)[:, None]
        rp = R.apply(pts)

        def ang(p):
            return np.arccos(p[:, 2]), np.arctan2(p[:, 1], p[:, 0])
        A = real_harm(l, *ang(pts)).T          # (npts, n)
        B = real_harm(l, *ang(rp)).T
        D = np.linalg.lstsq(A, B, rcond=None)[0].T   # S(R r) = D S(r)
        err_d = np.abs(D @ D.T - np.eye(n)).max()
        Ur = np.einsum('ia,jb,kc,ld,abcd->ijkl', D, D, D, D, U)
        report(f"l={l}: D orthogonal / U rotation invariant", max(err_d, np.abs(Ur - U).max()), 1e-11)
    # J for d from F2, F4
    F = {0: 0.2, 2: 0.30, 4: 0.19}
    U = coulomb_real(2, F).real
    J = np.einsum('abab->', U) / 25 - (np.einsum('abab->', U) - np.einsum('abba->', U)) / 20
    report("d: J = (F2+F4)/14", abs(J - (0.30 + 0.19) / 14), 1e-14)
    F = {0: 0.2, 2: 0.30, 4: 0.19, 6: 0.14}
    U = coulomb_real(3, F).real
    J = np.einsum('abab->', U) / 49 - (np.einsum('abab->', U) - np.einsum('abba->', U)) / 42
    report("f: J = (286F2+195F4+250F6)/6435", abs(J - (286 * .3 + 195 * .19 + 250 * .14) / 6435), 1e-14)


# ------------------------------------------------- energy functionals (NC form)
# n has shape (2, 2, m, m): n[s, t] = n^{st}, Hermitian as a 2m x 2m matrix.
def flat(n):
    m = n.shape[-1]
    return n.transpose(0, 2, 1, 3).reshape(2 * m, 2 * m)


def unflat(x, m):
    return x.reshape(2, m, 2, m).transpose(0, 2, 1, 3)


def e_dudarev(n, U):
    x = flat(n)
    return (0.5 * U * (np.trace(x) - np.trace(x @ x))).real


def v_dudarev(n, U):
    m = n.shape[-1]
    return unflat(U * (0.5 * np.eye(2 * m) - flat(n).T), m)


def e_hf(n, Ut):
    N = n[0, 0] + n[1, 1]
    e = 0.5 * np.einsum('ijkl,ik,jl->', Ut, N, N)
    e -= 0.5 * sum(np.einsum('ijlk,ik,jl->', Ut, n[s, t], n[t, s]) for s in range(2) for t in range(2))
    return e.real


def v_hf(n, Ut):
    N = n[0, 0] + n[1, 1]
    v = np.zeros_like(n)
    for s in range(2):
        v[s, s] += np.einsum('ijkl,jl->ik', Ut, N)
        for t in range(2):
            v[s, t] -= np.einsum('ijlk,jl->ik', Ut, n[t, s])
    return v


def tr2(n):
    return np.array([[np.trace(n[s, t]) for t in range(2)] for s in range(2)])


def e_general(n, Ut, U, J, dc):
    m = n.shape[-1]
    Ns = tr2(n)
    N = (Ns[0, 0] + Ns[1, 1]).real
    if dc in ('sFLL', 'cFLL'):
        e = e_hf(n, Ut)
        if dc == 'sFLL':
            e -= 0.5 * U * N * (N - 1) - 0.5 * J * ((Ns @ Ns).trace().real - N)
        else:
            e -= 0.5 * U * N * (N - 1) - 0.25 * J * N * (N - 2)
        return e
    return e_hf(amf_shift(n, dc), Ut)


def amf_shift(n, dc):
    m = n.shape[-1]
    Ns = tr2(n)
    nt = n.copy()
    if dc == 'sAMF':
        for s in range(2):
            for t in range(2):
                nt[s, t] = n[s, t] - Ns[s, t] / m * np.eye(m)
    else:  # cAMF
        N = Ns[0, 0] + Ns[1, 1]
        for s in range(2):
            nt[s, s] = n[s, s] - N / (2 * m) * np.eye(m)
    return nt


def v_general(n, Ut, U, J, dc):
    m = n.shape[-1]
    Ns = tr2(n)
    N = (Ns[0, 0] + Ns[1, 1])
    if dc in ('sAMF', 'cAMF'):
        return v_hf(amf_shift(n, dc), Ut)
    v = v_hf(n, Ut)
    I = np.eye(m)
    for s in range(2):
        for t in range(2):
            if dc == 'sFLL':
                d = U * (N - 0.5) * (s == t) - J * (Ns[t, s] - 0.5 * (s == t))
            else:
                d = (U * (N - 0.5) - J * (N / 2 - 0.5)) * (s == t)
            v[s, t] -= d * I
    return v


def rand_occ(m, collinear):
    a = rng.normal(size=(2 * m, 2 * m)) + 1j * rng.normal(size=(2 * m, 2 * m))
    if collinear:
        a = a.real
        a[:m, m:] = 0
        a[m:, :m] = 0
    h = a @ a.conj().T
    w, q = np.linalg.eigh(h)
    occ = rng.uniform(0, 1, size=2 * m)
    x = (q * occ) @ q.conj().T
    return unflat(x, m)


def fd_check(efun, vfun, n, collinear, label):
    m = n.shape[-1]
    worst = 0.0
    for _ in range(4):
        d = rand_occ(m, collinear) - rand_occ(m, collinear)
        eps = 1e-5
        num = (efun(n + eps * d) - efun(n - eps * d)) / (2 * eps)
        ana = np.sum(vfun(n) * d).real
        worst = max(worst, abs(num - ana) / max(1.0, abs(ana)))
    report(label, worst, 1e-8)


def check_functionals():
    l = 2
    m = 2 * l + 1
    U, J = 0.22, 0.03
    Ut = coulomb_real(l, slater_from_UJ(l, U, J)).real
    for coll in (True, False):
        n = rand_occ(m, coll)
        tag = "collinear" if coll else "non-collinear"
        fd_check(lambda x: e_dudarev(x, U), lambda x: v_dudarev(x, U), n, coll, f"Dudarev v=dE/dn ({tag})")
        for dc in ('sFLL', 'sAMF', 'cFLL', 'cAMF'):
            fd_check(lambda x, dc=dc: e_general(x, Ut, U, J, dc), lambda x, dc=dc: v_general(x, Ut, U, J, dc),
                     n, coll, f"general {dc} v=dE/dn ({tag})")
    # sFLL with J=0 equals Dudarev with U=F0
    Ut0 = coulomb_real(l, slater_from_UJ(l, U, 0.0)).real
    n = rand_occ(m, False)
    err = abs(e_general(n, Ut0, U, 0.0, 'sFLL') - e_dudarev(n, U))
    err = max(err, np.abs(v_general(n, Ut0, U, 0.0, 'sFLL') - v_dudarev(n, U)).max())
    report("sFLL(J=0) == Dudarev(U=F0): energy and potential", err, 1e-12)
    # Dudarev energy vanishes for idempotent (integer) occupation
    q = np.linalg.qr(rng.normal(size=(2 * m, 2 * m)) + 1j * rng.normal(size=(2 * m, 2 * m)))[0]
    x = q[:, :4] @ q[:, :4].conj().T
    report("Dudarev E=0 for idempotent n", abs(e_dudarev(unflat(x, m), U)), 1e-14)
    # Dudarev energy and potential invariant under a global spin rotation
    n = rand_occ(m, False)
    Rs = np.array([[np.cos(.3), -np.sin(.3) * np.exp(-0.2j)], [np.sin(.3) * np.exp(0.2j), np.cos(.3)]])
    W = np.kron(Rs, np.eye(m))
    nr = unflat(W @ flat(n) @ W.conj().T, m)
    err = max(abs(e_dudarev(nr, U) - e_dudarev(n, U)),
              abs(e_general(nr, Ut, U, J, 'sFLL') - e_general(n, Ut, U, J, 'sFLL')),
              abs(e_general(nr, Ut, U, J, 'sAMF') - e_general(n, Ut, U, J, 'sAMF')))
    report("spin-rotation invariance (Dudarev, sFLL, sAMF)", err, 1e-12)


# ------------------------------------------------------- dual / full occupation
def toy_system(nc):
    """Atom A: s + d (6 orbitals, U on the d block); atom B: 3 orbitals."""
    na, nb = 6, 3
    N = na + nb
    b = rng.normal(size=(N, N)) * 0.25
    S0 = np.eye(N) + (b + b.T) / 2
    S0[:na, :na] = np.eye(na)        # on-site blocks: orthonormal on one atom
    S0[na:, na:] = np.eye(nb)
    dS = np.zeros((N, N))
    g = rng.normal(size=(na, nb)) * 0.1
    dS[:na, na:] = g
    dS[na:, :na] = g.T              # only the inter-atomic block depends on R
    ns = 2 if nc else 1
    Nt = ns * N
    h = rng.normal(size=(Nt, Nt)) + (1j * rng.normal(size=(Nt, Nt)) if nc else 0)
    h = (h + h.conj().T) / 2
    S2 = np.kron(np.eye(ns), S0)
    w, c = np.linalg.eigh(np.linalg.solve(np.linalg.cholesky(S2), h) @ np.linalg.inv(np.linalg.cholesky(S2)).conj().T)
    c = np.linalg.solve(np.linalg.cholesky(S2).conj().T, c)
    f = np.zeros(Nt)
    f[: Nt // 2] = 1.0
    f[Nt // 2] = 0.37
    rho = np.conj((c * f) @ c.conj().T)    # rho_ab = sum f conj(c_a) c_b
    return na, nb, S0, dS, rho, f.sum()


def dual_occ(rho, S2):
    return 0.5 * (rho @ S2 + S2 @ rho)


def full_occ(rho, S2):
    return S2 @ rho @ S2


def check_dual():
    for nc in (False, True):
        na, nb, S0, dS, rho, Ne = toy_system(nc)
        ns = 2 if nc else 1
        N = na + nb
        tag = "NC" if nc else "collinear"
        S2 = np.kron(np.eye(ns), S0)
        nfull = dual_occ(rho, S2)
        # sum rule: sum over all on-site blocks of Tr n = N_e
        report(f"dual sum rule Tr n = N_e ({tag})", abs(np.trace(nfull).real - Ne), 1e-12)
        idx = [s * N + i for s in range(ns) for i in range(1, na)]   # d block of atom A
        U = 0.2

        def energy(S, occ):
            x = occ(rho, np.kron(np.eye(ns), S))[np.ix_(idx, idx)]
            return (0.5 * U * (np.trace(x) - np.trace(x @ x))).real * (2 if not nc else 1)

        for occ, name in ((dual_occ, "dual"), (full_occ, "full")):
            x = occ(rho, S2)[np.ix_(idx, idx)]
            v = U * (0.5 * np.eye(len(idx)) - x.T)          # v = dE/dn (per spin set)
            V = np.zeros((ns * N, ns * N), complex)
            V[np.ix_(idx, idx)] = v
            HU = 0.5 * (V @ S2 + S2 @ V) if name == "dual" else S2 @ V @ S2
            g = 2 if not nc else 1
            # identity: sum_ab H_ab rho_ab = sum v_ab n_ab
            err = abs(np.sum(HU * rho) - np.sum(v * x))
            report(f"{name}: sum H_U rho = sum v n ({tag})", err, 1e-12)
            dS2 = np.kron(np.eye(ns), dS)
            dHU = 0.5 * (V @ dS2 + dS2 @ V) if name == "dual" else dS2 @ V @ S2 + S2 @ V @ dS2
            ana = g * np.sum(dHU * rho).real
            eps = 1e-6
            num = (energy(S0 + eps * dS, occ) - energy(S0 - eps * dS, occ)) / (2 * eps)
            report(f"{name}: dE_U/dR at fixed rho = sum rho dH_U/dR ({tag})", abs(num - ana), 1e-8)


# --------------------------------------------------------- constraint / Zeeman
def spin_T(n):
    """T_{st} = sum_i n^{ts}_{ii} (charge-density convention, Hermitian 2x2)."""
    return tr2(n).T


def moment(T):
    return np.array([np.trace(s @ T).real for s in SIG])


def check_constraint():
    m = 5
    n = rand_occ(m, False)
    V = 0.3
    th0, ph0 = 0.7, 2.1
    e0 = np.array([np.sin(th0) * np.cos(ph0), np.sin(th0) * np.sin(ph0), np.cos(th0)])
    # m-vector convention: m = (2Re T01, -2Im T01, T00-T11)
    T = spin_T(n)
    mv = moment(T)
    err = abs(mv[0] - 2 * T[0, 1].real) + abs(mv[1] + 2 * T[0, 1].imag) + abs(mv[2] - (T[0, 0] - T[1, 1]).real)
    report("m = (2Re T01, -2Im T01, T00-T11)", err, 1e-14)
    # TN2 Eq.16 rotation matrix: U^dag diag(a,b) U = (a+b)/2 + (a-b)/2 e.sigma
    Um = np.array([[np.exp(.5j * ph0) * np.cos(th0 / 2), np.exp(-.5j * ph0) * np.sin(th0 / 2)],
                   [-np.exp(.5j * ph0) * np.sin(th0 / 2), np.exp(-.5j * ph0) * np.cos(th0 / 2)]])
    a, b = 3.1, 1.4
    T0 = Um.conj().T @ np.diag([a, b]) @ Um
    ref = 0.5 * (a + b) * np.eye(2) + 0.5 * (a - b) * sum(e0[c] * SIG[c] for c in range(3))
    report("U^dag diag U is aligned with e(theta,phi) in the T convention", np.abs(T0 - ref).max(), 1e-14)

    def target(nn, mode, M0=1.7):
        TT = spin_T(nn)
        Ntot = np.trace(TT).real
        mm = np.linalg.norm(moment(TT)) if mode == 'on' else M0
        return 0.5 * (Ntot * np.eye(2) + mm * sum(e0[c] * SIG[c] for c in range(3)))

    def e_cs(nn, mode):
        D = spin_T(nn) - target(nn, mode)
        return V * np.sum(np.abs(D) ** 2)

    for mode in ('on', 'on2'):
        mv = moment(spin_T(n))
        if mode == 'on':
            closed = V * mv @ mv * (1 - mv @ e0 / np.linalg.norm(mv))
            h = V * (2 * mv - mv / np.linalg.norm(mv) * (mv @ e0) - np.linalg.norm(mv) * e0)
        else:
            closed = 0.5 * V * np.sum((mv - 1.7 * e0) ** 2)
            h = V * (mv - 1.7 * e0)
        report(f"constraint {mode}: V||T-T0||^2 closed form", abs(e_cs(n, mode) - closed), 1e-13)
        hs = sum(h[c] * SIG[c] for c in range(3))
        vfun = (lambda x, hs=hs: np.einsum('st,ij->stij', hs, np.eye(m)))
        # v must be evaluated at the current n; finite difference around n
        fd_check(lambda x, mode=mode: e_cs(x, mode), lambda x: vfun(x), n, False,
                 f"constraint {mode}: v = (h.sigma) on the diagonal")
    # Zeeman spin: E = b.s, s = m/2, b = -B e0 ; v = (b/2).sigma
    Bm = 2.0e-4
    bvec = -Bm * e0
    ez = lambda x: bvec @ moment(spin_T(x)) / 2
    vz = lambda x: np.einsum('st,ij->stij', sum(bvec[c] * SIG[c] for c in range(3)) / 2, np.eye(m))
    fd_check(ez, vz, n, False, "Zeeman spin: v = (b/2).sigma")


# ---------------------------------------------------- orbital polarization rule
def polarize(eigs):
    """Occupations for the enhancement: floor(N) highest get 1, next gets frac."""
    N = eigs.sum()
    k = int(np.floor(N))
    order = np.argsort(eigs)[::-1]
    occ = np.zeros_like(eigs)
    occ[order[:k]] = 1.0
    if k < len(eigs):
        occ[order[k]] = N - k
    return occ


def check_polarize():
    worst = 0.0
    for _ in range(20):
        e = rng.uniform(0, 1, size=10)
        occ = polarize(e)
        worst = max(worst, abs(occ.sum() - e.sum()))
    report("orbital-polarization occupations conserve N", worst, 1e-12)


# ------------------------------------------------------------- Yukawa kernel
def check_yukawa():
    lam = 0.9
    r1, r2, cg = 0.8, 2.3, 0.37
    R = np.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * cg)
    exact = np.exp(-lam * R) / R
    ser = sum((2 * k + 1) * lam * spherical_in(k, lam * r1) * (2 / np.pi) * spherical_kn(k, lam * r2)
              * eval_legendre(k, cg) for k in range(60))
    report("Yukawa: exp(-lam R)/R = lam sum (2k+1) i_k k~_k P_k", abs(ser - exact), 1e-12)
    lam = 1e-4
    err = max(abs((2 * k + 1) * lam * spherical_in(k, lam * r1) * (2 / np.pi) * spherical_kn(k, lam * r2)
                  - r1**k / r2**(k + 1)) / (r1**k / r2**(k + 1)) for k in range(4))
    report("Yukawa kernel -> r<^k/r>^(k+1) as lam -> 0", err, 1e-3)


def check_units():
    au_tesla = 2.35051757e5      # CODATA 2018 atomic unit of magnetic flux density
    report("1 T = 1/2.35051742e5 a.u. vs CODATA (relative)", abs(2.35051742e5 / au_tesla - 1), 1e-6)


if __name__ == "__main__":
    check_real_harmonics_signs()
    check_coulomb()
    check_functionals()
    check_dual()
    check_constraint()
    check_polarize()
    check_yukawa()
    check_units()
    print("ALL OK" if all(results) else "SOME CHECKS FAILED")
