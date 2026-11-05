# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/density-mixing.md:
#  1. Fermi-Dirac occupations and chemical-potential bisection (two-level system).
#  2. Real-space density matrix of a 1D two-orbital chain with overlap, using the
#     phase convention rho_ab(L) = sum_k w_k sum_n f_n Re[e^{i k.L} c_an^* c_bn],
#     and the invariants: electron count, band energy, Hermiticity,
#     idempotency rho S rho = rho at T -> 0 for an insulator.
#  3. Mixing on linear model problems: simple, adaptive weight, Pulay (DIIS),
#     Kerker preconditioning of a Thomas-Fermi-like response.
#  4. Coordinate-based extrapolation coefficients (constrained / unconstrained LS).
# Run: uv run --with numpy --with scipy python tools/check_density_mixing.py

import numpy as np
from scipy.linalg import eigh

rng = np.random.default_rng(1)
ok = True


def check(name, err, tol):
    global ok
    flag = err <= tol
    ok &= flag
    print(f"{'OK ' if flag else 'NG '} {name}: err={err:.3e} tol={tol:.0e}")


def fermi(e, mu, beta):
    x = np.clip((e - mu) * beta, -60.0, 60.0)
    return 1.0 / (1.0 + np.exp(x))


def bisect_mu(eps, wk, nelec, beta, spin_deg, lo=-20.0, hi=20.0, tol=1e-12):
    # eps[k, n], wk[k] normalized to 1
    for _ in range(1000):
        mu = 0.5 * (lo + hi)
        n = spin_deg * np.sum(wk[:, None] * fermi(eps, mu, beta))
        d = nelec - n
        if d >= 0:
            lo = mu
        else:
            hi = mu
        if abs(d) < tol:
            break
    return mu


kB = 8.617251324e-5 / 27.2113845  # Hartree/K

# ---- 1. two-level system --------------------------------------------------
e1, e2 = -0.30, 0.10
eps = np.array([[e1, e2]])
wk = np.array([1.0])
for T in (300.0, 3000.0):
    beta = 1.0 / (kB * T)
    mu = bisect_mu(eps, wk, 2.0, beta, 2.0)
    # exact solution is the midpoint, but N(mu) is flat inside the gap: the
    # bisection stops on |N - N0| < 1e-12, so mu itself is only fixed loosely.
    n2 = 2 * (fermi(e1, mu, beta) + fermi(e2, mu, beta))
    check(f"two-level N=2 electron count at T={T}", abs(n2 - 2.0), 1e-11)
    print(f"    mu - (e1+e2)/2 = {mu - 0.5*(e1+e2):.3e} Ha (inside the gap: {e1 < mu < e2})")
    ok &= e1 < mu < e2
    # N=3: closed form f2 = 1/2 -> mu = e2 exactly? no: f1+f2 = 1.5
    mu3 = bisect_mu(eps, wk, 3.0, beta, 2.0)
    n3 = 2 * (fermi(e1, mu3, beta) + fermi(e2, mu3, beta))
    check(f"two-level N=3 electron count at T={T}", abs(n3 - 3.0), 1e-11)
# N=3 with nondegenerate levels, analytic: solve f1 + f2 = 1.5
beta = 1.0 / (kB * 3000.0)
from scipy.optimize import brentq
mu_exact = brentq(lambda m: 2 * (fermi(e1, m, beta) + fermi(e2, m, beta)) - 3.0, -1, 1, xtol=1e-15)
mu3 = bisect_mu(eps, wk, 3.0, beta, 2.0)
check("two-level N=3 analytic mu (3000 K)", abs(mu3 - mu_exact), 1e-9)

# ---- 2. 1D two-orbital chain with overlap --------------------------------
# orbitals a=0 (on-site e_a), b=1 (e_b); hoppings inside cell and to next cell
ea, eb = -0.4, 0.3
t_in, t_out = -0.20, -0.12
s_in, s_out = 0.15, 0.08


def HS_L():
    H = {0: np.array([[ea, t_in], [t_in, eb]]),
         1: np.array([[0.0, 0.0], [t_out, 0.0]]),   # <phi_a(0)|H|phi_b(L=1)>: b-orbital of cell 0 to a of cell 1
         }
    S = {0: np.array([[1.0, s_in], [s_in, 1.0]]),
         1: np.array([[0.0, 0.0], [s_out, 0.0]])}
    H[-1] = H[1].T.copy(); S[-1] = S[1].T.copy()
    return H, S


Hl, Sl = HS_L()
Nk = 12
ks = 2 * np.pi * np.arange(Nk) / Nk  # Gamma-centred grid (closed under k -> -k), lattice constant 1
wk = np.full(Nk, 1.0 / Nk)


def bloch(O, k):  # O(k) = sum_L e^{ikL} O(L)
    return sum(np.exp(1j * k * L) * O[L] for L in O)


eigs, vecs = [], []
for k in ks:
    e, c = eigh(bloch(Hl, k), bloch(Sl, k))
    eigs.append(e); vecs.append(c)
eigs = np.array(eigs)
nelec = 2.0  # one band per spin: insulator
for T, tolidem in ((1.0, 1e-10), (3000.0, None)):
    beta = 1.0 / (kB * T)
    mu = bisect_mu(eigs, wk, nelec, beta, 2.0)
    Lrange = range(-(Nk // 2) + 1, Nk // 2 + 1)
    rho = {L: np.zeros((2, 2)) for L in Lrange}
    edm = {L: np.zeros((2, 2)) for L in Lrange}
    for ik, k in enumerate(ks):
        f = fermi(eigs[ik], mu, beta)
        c = vecs[ik]
        for L in Lrange:
            ph = np.exp(1j * k * L)
            M = np.einsum('n,an,bn->ab', f, c.conj(), c)
            Me = np.einsum('n,an,bn->ab', f * eigs[ik], c.conj(), c)
            rho[L] += wk[ik] * np.real(ph * M)
            edm[L] += wk[ik] * np.real(ph * Me)
    N = 2 * sum(np.sum(rho[L] * Sl[L]) for L in Sl)
    check(f"chain electron count (T={T})", abs(N - nelec), 1e-12)
    Eb = 2 * sum(np.sum(rho[L] * Hl[L]) for L in Hl)
    Eb_k = 2 * np.sum(wk[:, None] * fermi(eigs, mu, beta) * eigs)
    check(f"chain band energy sum rho H = sum f eps (T={T})", abs(Eb - Eb_k), 1e-12)
    # EDM: sum E_ab(L) S_ab(L) = band energy as well
    Eb2 = 2 * sum(np.sum(edm[L] * Sl[L]) for L in Sl)
    check(f"chain sum EDM*S = band energy (T={T})", abs(Eb2 - Eb_k), 1e-12)
    herm = max(np.max(np.abs(rho[L] - rho[-L].T)) for L in Lrange if -L in rho)
    check(f"chain Hermiticity rho(L)=rho(-L)^T (T={T})", herm, 1e-14)
    if tolidem:
        # cyclic convolution over the Nk-cell supercell
        def g(O, L):
            Lm = ((L + Nk // 2 - 1) % Nk) - Nk // 2 + 1
            return O.get(Lm, np.zeros((2, 2)))
        err = 0.0
        for L in Lrange:
            acc = np.zeros((2, 2))
            for L1 in Lrange:
                for L2 in (-1, 0, 1):
                    acc += g(rho, L1) @ Sl[L2] @ g(rho, L - L1 - L2)
            err = max(err, np.max(np.abs(acc - rho[L])))
        check("chain idempotency rho S rho = rho (T=1 K)", err, tolidem)


# time-reversal folding: k and -k merged into one point of weight 2 (real part only)
beta = 1.0 / (kB * 3000.0)
mu = bisect_mu(eigs, wk, nelec, beta, 2.0)
def build_rho(klist, wlist, L):
    acc = np.zeros((2, 2))
    for k, w in zip(klist, wlist):
        e, c = eigh(bloch(Hl, k), bloch(Sl, k))
        f = fermi(e, mu, beta)
        acc += w * np.real(np.exp(1j * k * L) * np.einsum('n,an,bn->ab', f, c.conj(), c))
    return acc
half = [k for k in ks if 0 <= k <= np.pi + 1e-12]
whalf = [(1.0 if (abs(k) < 1e-12 or abs(k - np.pi) < 1e-12) else 2.0) / Nk for k in half]
err = max(np.max(np.abs(build_rho(ks, wk, L) - build_rho(half, whalf, L))) for L in (-2, -1, 0, 1, 2))
check("time-reversal folding of k points", err, 1e-12)

# noncollinear spinor with uncoupled spins reproduces the collinear rho (complex form, no Re)
def rho_nc(L):
    acc = np.zeros((4, 4), complex)
    for k, w in zip(ks, wk):
        H2 = np.kron(np.eye(2), bloch(Hl, k)); S2 = np.kron(np.eye(2), bloch(Sl, k))
        e, c = eigh(H2, S2)
        f = fermi(e, mu, beta)
        acc += w * np.exp(1j * k * L) * np.einsum('n,pn,qn->pq', f, c.conj(), c)
    return acc
err = 0.0
for L in (-1, 0, 1):
    r = rho_nc(L); rc = build_rho(ks, wk, L)
    err = max(err, np.max(np.abs(r[:2, :2] - rc)), np.max(np.abs(r[2:, 2:] - rc)), np.max(np.abs(r[:2, 2:])))
    # Hermiticity rho^{s s'}_{ab}(L) = conj(rho^{s' s}_{ba}(-L))
    err = max(err, np.max(np.abs(r - rho_nc(-L).conj().T)))
check("noncollinear (uncoupled) rho equals collinear rho; NC Hermiticity", err, 1e-12)

# ---- 3. mixing on linear models ------------------------------------------
n = 40
lam = np.concatenate([np.linspace(-0.9, 0.95, n - 4), [-3.0, -2.5, 0.99, 0.995]])
Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
Jac = Q @ np.diag(lam) @ Q.T          # x_out - x* = J (x_in - x*)
xs = rng.normal(size=n)


def F(x):
    return xs + Jac @ (x - xs)


def run_simple(w, iters=20000):
    x = np.zeros(n)
    for i in range(iters):
        r = F(x) - x
        if np.linalg.norm(r) < 1e-10:
            return i
        x = x + w * r
    return None


# stability bound: |1 - w(1-lam)| < 1 -> w < 2/(1-lam_min) = 0.5
print("simple w=0.45 iters:", run_simple(0.45), " w=0.55 iters:", run_simple(0.55))
ok &= run_simple(0.55) is None


def adapt(w, n0, n1, wmin, wmax):
    if n0 < n1:
        w = w * n1 / max(n1 - n0, 1e-9)
    elif n1 < n0:
        w = w * n1 / (n1 + n0)
    return min(max(w, wmin), wmax)


def run_pulay(hist=8, start=3, w0=0.1, iters=200, coef_table=True):
    x = np.zeros(n)
    X, R = [], []
    w = w0; nprev = None
    for it in range(1, iters + 1):
        r = F(x) - x
        nr = r @ r
        if np.sqrt(nr) < 1e-10:
            return it
        X.insert(0, x.copy()); R.insert(0, r.copy())
        X, R = X[:hist], R[:hist]
        if it < start or len(R) < 2:
            if nprev is not None:
                w = adapt(w, nr, nprev, 1e-3, 0.4)
            x = x + w * r
        else:
            m = len(R)
            A = np.array([[Ri @ Rj for Rj in R] for Ri in R])
            B = np.zeros((m + 1, m + 1)); B[:m, :m] = A; B[:m, m] = -1; B[m, :m] = -1
            rhs = np.zeros(m + 1); rhs[m] = -1
            sol = np.linalg.solve(B, rhs)
            c = sol[:m]
            # closed form check
            Ai1 = np.linalg.solve(A, np.ones(m))
            assert np.allclose(c, Ai1 / Ai1.sum(), atol=1e-8 * max(1, np.abs(c).max()))
            xbar = sum(ci * Xi for ci, Xi in zip(c, X))
            rbar = sum(ci * Ri for ci, Ri in zip(c, R))
            a = 1.0
            if coef_table:
                a = 0.5 if nr >= 1e-1 else 0.6 if nr >= 1e-2 else 0.7 if nr >= 1e-3 else 0.8 if nr >= 1e-4 else 1.0
            x = xbar + a * rbar
        nprev = nr
    return None


it_p = run_pulay(hist=n + 2, start=3, coef_table=False)
print("Pulay (full history, linear problem) iterations:", it_p)
ok &= it_p is not None and it_p <= n + 5
it_p8 = run_pulay(hist=8, iters=2000)
print("Pulay (history 8, OpenMX-like coefficient table) iterations:", it_p8)
ok &= it_p8 is not None

# Kerker on a 1D periodic model: response n_out - n* = -(qTF^2/q^2) (n_in - n*), q != 0
M = 64; Lbox = 200.0
q = 2 * np.pi * np.fft.fftfreq(M, d=Lbox / M)
qTF = 1.0
nstar = rng.normal(size=M); nstar -= nstar.mean()


def Fq(nin_q):
    out = nin_q.copy()
    mask = q != 0
    out[mask] = np.fft.fft(nstar)[mask] - (qTF ** 2 / q[mask] ** 2) * (nin_q[mask] - np.fft.fft(nstar)[mask])
    return out


def run_kerker(w, q0, modified=False, iters=500, plain=False):
    nq = np.zeros(M, complex)
    for it in range(iters):
        r = Fq(nq) - nq
        if np.linalg.norm(r) / np.sqrt(M) < 1e-10:
            return it
        if plain:
            g = np.ones(M)
        elif modified:   # weight w*sqrt((q^2+(0.2 q0)^2)/(q^2+q0^2)), 1 at q=0
            g = np.where(q == 0, 1.0, np.sqrt((q ** 2 + (0.2 * q0) ** 2) / (q ** 2 + q0 ** 2)))
        else:          # textbook Kerker
            g = np.where(q == 0, 1.0, q ** 2 / (q ** 2 + q0 ** 2))
        nq = nq + w * g * r
    return None


itk = run_kerker(1.0, qTF)
print("Kerker w=1, q0=qTF iterations (expect 1-2):", itk)
ok &= itk is not None and itk <= 2
print("simple mixing w=0.3 on Kerker model (expect divergence/None):", run_kerker(0.3, 0.0, plain=True))
ok &= run_kerker(0.3, 0.0, plain=True) is None
# bounded weight w*sqrt((q^2+(0.2q0)^2)/(q^2+q0^2)) -> 0.2 w at q -> 0: the slowest mode
# contracts by 1 - 0.2 w (1 + qTF^2/qmin^2); stable only if that stays above -1.
qmin = 2 * np.pi / Lbox
print("bounded-Kerker factor at qmin:", 1 - 0.2 * 0.3 * (1 + qTF**2 / qmin**2), "(< -1 means divergence)")
# Kerker with q0 = qTF: modes contract by |1 - w (q^2+qTF^2)/(q^2+q0^2)| = 1-w
print("Kerker w=0.5 iterations:", run_kerker(0.5, qTF))

# ---- 4. extrapolation coefficients ---------------------------------------
natom = 5
a0, v0, acc = rng.normal(size=3 * natom), rng.normal(size=3 * natom), rng.normal(size=3 * natom)
traj = [a0 + v0 * t + 0.5 * acc * t * t for t in (0.0, 0.1, 0.2, 0.3)]
hist = [traj[2], traj[1], traj[0]]   # most recent first
xnew = traj[3]
A = np.array([[hi @ hj for hj in hist] for hi in hist])
B = np.array([hi @ xnew for hi in hist])
c_unc = np.linalg.solve(A, B)
check("extrapolation unconstrained LS -> (3,-3,1)", np.max(np.abs(c_unc - [3, -3, 1])), 1e-8)
Bc = np.zeros((4, 4)); Bc[:3, :3] = A; Bc[:3, 3] = 1; Bc[3, :3] = 1
c_con = np.linalg.solve(Bc, np.concatenate([B, [1.0]]))[:3]
check("extrapolation constrained LS -> (3,-3,1)", np.max(np.abs(c_con - [3, -3, 1])), 1e-8)


# ---- 5. automatic Kerker factor -------------------------------------------
def auto_q0(a1, a2, a3):
    A = np.array([a1, a2, a3]); B = 2 * np.pi * np.linalg.inv(A).T
    b = np.linalg.norm(B, axis=1)
    bbar = b.mean(); db = (abs(b[1]-b[0]) + abs(b[2]-b[1]) + abs(b[0]-b[2])) / 3
    gamma = 0.5 / b.min() * (2 * db / bbar + 1)
    return gamma * b.min()
check("auto Kerker q0 = 0.5 for a cubic cell", abs(auto_q0([7, 0, 0], [0, 7, 0], [0, 0, 7]) - 0.5), 1e-14)
print("auto Kerker q0 for a 5x5x40 bohr cell:", auto_q0([5, 0, 0], [0, 5, 0], [0, 0, 40]))

print("ALL OK" if ok else "SOME CHECKS FAILED")
