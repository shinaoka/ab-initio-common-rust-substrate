# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/geometry-md.md:
#  1. Unit conversions (fs, amu, m/s, GPa, K) used by the MD chapter.
#  2. Velocity-Verlet NVE on harmonic oscillators: bounded energy error O(dt^2).
#  3. The "full kick + Taylor drift" scheme used by the velocity-scaling family
#     (no scaling): characteristic roots have modulus sqrt(1 - h^2/2), i.e.
#     the scheme is dissipative and only first-order accurate.
#  4. Nose-Hoover (velocity-Verlet type: explicit predictor for zeta, implicit
#     corrector solved by Newton) conserves E + K + Q zeta^2/2 + g kT eta with
#     dzeta/dt = (2K - g kT)/Q.
#  5. Nose-Hoover with dzeta/dt = (2K - g kT)/(3Q) (NPT_NH form), integrated
#     accurately: the extended energy is conserved with 3Q zeta^2/2, not Q.
#  6. Velocity-scaling factor and random initial velocities (zero momentum).
#  7. Optimizers on a quadratic model: GDIIS coefficients, BFGS secant
#     conditions, RF step = (H - lambda)^-1 (-g), convergence of SD / DIIS /
#     BFGS / RF / EF variants with the 0.5 bohr step cap.
#  8. Schlegel-type bond-stretch model Hessian: translational invariance.
#  9. Parrinello-Rahman geometry: Sigma = Omega h^{-T}, d Omega / d h = Sigma,
#     and the derivative of Sigma used in the Wentzcovitch-type equation.
# 10. NEB tangent and projected force on a 2D model potential; DIIS+BFGS
#     NEB optimisation converges, spacing equalised, top image near saddle.
# 11. Periodic wrap of fractional coordinates.
# Run: uv run --with numpy --with scipy python tools/check_geometry_md.py

import numpy as np

rng = np.random.default_rng(7)
ok = True


def check(name, err, tol):
    global ok
    flag = bool(err <= tol)
    ok &= flag
    print(f"{'OK ' if flag else 'NG '} {name}: err={err:.3e} tol={tol:.0e}")


# ---- 1. units ---------------------------------------------------------------
fs = 41.3411                       # atomic time units per fs (as used)
amu = 1.660538921 / 0.000910938291  # electron masses per amu
mps = 0.4571028e-6                 # atomic velocity per m/s
gpa = 0.3398827e-4                 # Hartree/bohr^3 per GPa
kB = 8.617251324e-5 / 27.2113845   # Hartree/K
# independent CODATA-like values
au_time_s = 2.418884326e-17
bohr_m = 0.52917721e-10
check("1 fs in a.u. of time", abs(fs - 1e-15 / au_time_s) / fs, 1e-5)
check("1 m/s in a.u. of velocity", abs(mps - au_time_s / bohr_m) / mps, 1e-6)
hartree_J = 4.35974e-18
check("1 GPa in Hartree/bohr^3 (constant as used, 3e-5 rel. off)", abs(gpa - 1e9 * bohr_m**3 / hartree_J) / gpa, 1e-4)
check("amu / m_e", abs(amu - 1822.888) / amu, 1e-6)


# ---- 2. velocity Verlet -----------------------------------------------------
def vv_energy_drift(omega, m, dt, nsteps):
    x, v = 1.0, 0.0
    a = -omega**2 * x
    e0 = 0.5 * m * v * v + 0.5 * m * omega**2 * x * x
    emax = 0.0
    for _ in range(nsteps):
        vh = v + 0.5 * dt * a
        x = x + dt * vh
        a = -omega**2 * x
        v = vh + 0.5 * dt * a
        e = 0.5 * m * v * v + 0.5 * m * omega**2 * x * x
        emax = max(emax, abs(e - e0) / e0)
    return emax


m = 12.0 * amu
omega = 2 * np.pi / (33.0 * fs)       # ~1000 cm^-1 mode, period 33 fs
d1 = vv_energy_drift(omega, m, 0.5 * fs, 4000)
d2 = vv_energy_drift(omega, m, 0.25 * fs, 8000)
check("velocity Verlet: relative energy error bounded (dt=0.5 fs)", d1, 1e-2)
check("velocity Verlet: error ratio for dt/2 is ~4", abs(d1 / d2 - 4.0), 0.2)


# ---- 3. scheme of the velocity-scaling family -------------------------------
def vs_scheme(omega, dt, nsteps):
    # step 1: x <- x + v dt + a dt^2/2, v <- v + a dt
    # step n>1: v <- v + a_n dt ; x <- x + v dt + a_n dt^2/2
    x, v = 1.0, 0.0
    amp = []
    for n in range(1, nsteps + 1):
        a = -omega**2 * x
        if n == 1:
            x = x + v * dt + 0.5 * a * dt * dt
            v = v + a * dt
        else:
            v = v + a * dt
            x = x + v * dt + 0.5 * a * dt * dt
        amp.append(x * x + (v / omega) ** 2)
    return np.array(amp)


h = omega * 0.5 * fs
roots = np.roots([1.0, 1.5 * h * h - 2.0, 1.0 - 0.5 * h * h])
check("VS scheme: |z| = sqrt(1 - h^2/2)", abs(abs(roots[0]) - np.sqrt(1 - 0.5 * h * h)), 1e-12)
amp = vs_scheme(omega, 0.5 * fs, 2000)
ratio = (amp[-1] / amp[999]) ** (1.0 / 1000)
check("VS scheme: simulated decay per step = h^2/2 (relative)", abs((1 - ratio) - 0.5 * h * h) / (0.5 * h * h), 5e-2)


# ---- 4. Nose-Hoover (velocity Verlet, implicit corrector) -------------------
def nh_run(Q, T0, dt, nsteps, ndof=12):
    m = np.full(ndof, 2.0 * amu)
    w = np.linspace(0.6, 1.4, ndof) * omega
    x = rng.normal(size=ndof) * 0.05
    v = rng.normal(size=ndof) * np.sqrt(kB * 600.0 / m)
    g = ndof
    zeta = 0.0
    eta = 0.0

    def acc(x):
        return -w**2 * x

    def kin(v):
        return 0.5 * np.sum(m * v * v)

    def pot(x):
        return 0.5 * np.sum(m * w * w * x * x)

    def hext():
        return pot(x) + kin(v) + 0.5 * Q * zeta**2 + g * kB * T0 * eta

    a = acc(x)
    H = []
    for n in range(nsteps):
        # first half: v(t) -> v(t+dt/2), x -> x(t+dt), zeta predictor, eta
        vh = v + 0.5 * dt * (a - zeta * v)
        x = x + dt * vh
        G = (kin(v) - 0.5 * g * kB * T0) / Q          # half of dzeta/dt
        zeta_pred = zeta + G * dt
        zeta = zeta_pred
        eta = eta + zeta_pred * dt
        a = acc(x)
        # second half: solve v = vh + dt/2 (a - zeta v), zeta = zeta_pred + G(v) dt
        vv = vh.copy()
        for _ in range(50):
            r1 = vv - vh - 0.5 * dt * (a - zeta * vv)
            r2 = zeta - zeta_pred - (kin(vv) - 0.5 * g * kB * T0) * dt / Q
            # Newton on (vv, zeta)
            J11 = 1.0 + 0.5 * dt * zeta
            J12 = 0.5 * dt * vv
            J21 = -m * vv * dt / Q
            # eliminate vv: dv = (-r1 - J12 dz)/J11
            dz = (-r2 + np.sum(J21 * r1 / J11)) / (1.0 - np.sum(J21 * J12 / J11))
            dv = (-r1 - J12 * dz) / J11
            vv += dv
            zeta += dz
            if np.sum(dv * dv) + dz * dz < 1e-30:
                break
        v = vv
        H.append(hext())
    return np.array(H), kin(v)


Q = 20.0 * amu  # NH.Mass.HeatBath default 20 amu bohr^2
H1, _ = nh_run(Q, 300.0, 0.5 * fs, 4000)
H2, _ = nh_run(Q, 300.0, 0.25 * fs, 8000)
e1 = np.max(np.abs(H1 - H1[0]))
e2 = np.max(np.abs(H2 - H2[0]))
check("NH: extended energy conserved (relative to kT)", e1 / (kB * 300.0), 5e-2)
check("NH: extended energy error O(dt^2)", abs(e1 / e2 - 4.0), 1.0)


# ---- 5. NPT_NH thermostat form: dzeta/dt = (2K - g kT)/(3Q) -----------------
from scipy.integrate import solve_ivp


def nh_exact(Q, T0, tmax, ndof=12):
    m = np.full(ndof, 2.0 * amu)
    w = np.linspace(0.6, 1.4, ndof) * omega
    g = ndof

    def rhs(t, y):
        x, v, zeta, eta = y[:ndof], y[ndof:2*ndof], y[2*ndof], y[2*ndof+1]
        K = 0.5 * np.sum(m * v * v)
        T = 2 * K / (g * kB)
        dz = (g / 3.0) * kB * (T - T0) / Q
        return np.concatenate([v, -w**2 * x - zeta * v, [dz, zeta]])

    y0 = np.concatenate([np.zeros(ndof), rng.normal(size=ndof) * np.sqrt(kB * 600.0 / m), [0.0, 0.0]])
    sol = solve_ivp(rhs, (0, tmax), y0, rtol=1e-11, atol=1e-14, dense_output=False, max_step=tmax / 2000)
    x, v, zeta, eta = sol.y[:ndof], sol.y[ndof:2*ndof], sol.y[2*ndof], sol.y[2*ndof+1]
    U = 0.5 * np.sum((m * w * w)[:, None] * x * x, axis=0)
    K = 0.5 * np.sum(m[:, None] * v * v, axis=0)
    hq = U + K + 0.5 * Q * zeta**2 + g * kB * T0 * eta
    h3q = U + K + 1.5 * Q * zeta**2 + g * kB * T0 * eta
    return np.ptp(hq), np.ptp(h3q)


dq, d3q = nh_exact(Q, 300.0, 2000 * fs)
print(f"    NPT_NH form: spread of extended energy with Q: {dq:.3e}, with 3Q: {d3q:.3e} Ha")
check("NPT_NH form: conserved quantity uses 3Q (spread relative to kT)", d3q / (kB * 300.0), 1e-6)
check("NPT_NH form: the Q version is not conserved", float(dq < 1e3 * d3q), 0.0)


# ---- 6. velocity scaling, initial velocities --------------------------------
nat = 8
mass = rng.uniform(1, 60, nat) * amu
vel = rng.normal(size=(nat, 3)) * 1e-3
K = 0.5 * np.sum(mass[:, None] * vel**2)
T = K / (1.5 * kB * nat)
T0, alpha = 500.0, 0.3
s = np.sqrt((T0 + (T - T0) * alpha) / T)
K2 = 0.5 * np.sum(mass[:, None] * (s * vel) ** 2)
check("VS: scaled temperature = T0 + (T - T0) alpha", abs(K2 / (1.5 * kB * nat) - (T0 + (T - T0) * alpha)), 1e-9)
d = rng.uniform(-0.5, 0.5, (nat, 3))
d /= np.linalg.norm(d, axis=1)[:, None]
v0 = d * np.sqrt(3 * kB * 300.0 / mass)[:, None]
Pbar = np.sum(mass[:, None] * v0, axis=0) / nat
v0 -= Pbar[None, :] / mass[:, None]
check("initial velocities: total momentum zero", np.max(np.abs(np.sum(mass[:, None] * v0, axis=0))), 1e-12)


# ---- 7. optimizers on a quadratic -------------------------------------------
n = 9
A = rng.normal(size=(n, n))
Hq = A @ A.T / n + 0.05 * np.eye(n)   # positive definite, Hartree/bohr^2
xmin = rng.normal(size=n) * 0.3
def grad(x):
    return Hq @ (x - xmin)
def energy(x):
    return 0.5 * (x - xmin) @ Hq @ (x - xmin)

# GDIIS coefficients: bordered system == closed form
G = rng.normal(size=(4, n))
B = G @ G.T
M = np.zeros((5, 5)); M[:4, :4] = B / np.max(np.abs(B)); M[:4, 4] = 1; M[4, :4] = 1
rhs = np.zeros(5); rhs[4] = 1
c = np.linalg.solve(M, rhs)[:4]
cref = np.linalg.solve(B, np.ones(4)); cref /= cref.sum()
check("GDIIS: bordered system = A^-1 1 / (1^T A^-1 1)", np.max(np.abs(c - cref)), 1e-10)

# BFGS secant conditions
s_ = rng.normal(size=n); y_ = Hq @ s_
H0 = np.eye(n)
Hn = H0 + np.outer(y_, y_) / (y_ @ s_) - np.outer(H0 @ s_, H0 @ s_) / (s_ @ H0 @ s_)
check("BFGS (Hessian): H_new s = y", np.max(np.abs(Hn @ s_ - y_)), 1e-12)
Bi = np.eye(n)
sy = s_ @ y_
Bn = Bi + (sy + y_ @ Bi @ y_) / sy**2 * np.outer(s_, s_) - (np.outer(Bi @ y_, s_) + np.outer(s_, y_ @ Bi)) / sy
check("BFGS (inverse): Hinv_new y = s", np.max(np.abs(Bn @ y_ - s_)), 1e-12)
check("BFGS: inverse update = inverse of Hessian update", np.max(np.abs(Bn - np.linalg.inv(Hn))), 1e-9)

# RF step
g0 = grad(np.zeros(n))
aug = np.zeros((n + 1, n + 1)); aug[:n, :n] = Hq; aug[:n, n] = g0; aug[n, :n] = g0
w_, V = np.linalg.eigh(aug)
step = V[:n, 0] / V[n, 0]
lam = w_[0]
check("RF: step = -(H - lambda)^-1 g", np.max(np.abs(step + np.linalg.solve(Hq - lam * np.eye(n), g0))), 1e-10)
check("RF: lambda = g.s", abs(lam - g0 @ step), 1e-10)
check("RF: lambda below lowest Hessian eigenvalue", max(0.0, lam - np.linalg.eigvalsh(Hq)[0]), 0.0)

CAP = 0.5


def cap(xnew, xold):
    d = xnew - xold
    mx = np.max(np.abs(d))
    return xold + d * (CAP / mx) if mx > CAP else xnew


def run_opt(kind, hist=3, crit=3e-4, maxit=300):
    x = np.zeros(n)
    X, Gs = [], []
    Hm = np.eye(n); Hinv = np.eye(n)
    for it in range(1, maxit + 1):
        g = grad(x)
        if np.max(np.abs(g)) < crit:
            return it
        X.insert(0, x.copy()); Gs.insert(0, g.copy())
        X, Gs = X[:hist], Gs[:hist]
        k = len(X)
        Bm = np.array([[a @ b for b in Gs] for a in Gs])
        Mx = np.zeros((k + 1, k + 1)); Mx[:k, :k] = Bm / np.max(np.abs(Bm))
        Mx[:k, k] = 1; Mx[k, :k] = 1
        r = np.zeros(k + 1); r[k] = 1
        cc = np.linalg.solve(Mx, r)[:k]
        xt = sum(ci * xi for ci, xi in zip(cc, X))
        gt = sum(ci * gi for ci, gi in zip(cc, Gs))
        if k > 1:
            s_ = X[0] - X[1]; y_ = Gs[0] - Gs[1]
            if kind in ("RF", "EF"):
                sy = max(s_ @ y_, 1e-10) if kind == "RF" else s_ @ y_
                sHs = max(s_ @ Hm @ s_, 1e-10)
                if kind == "RF" or sy > 0:
                    Hm = Hm + np.outer(y_, y_) / sy - np.outer(Hm @ s_, Hm @ s_) / sHs
            if kind == "BFGS":
                sy = s_ @ y_
                Hinv = Hinv + (sy + y_ @ Hinv @ y_) / sy**2 * np.outer(s_, s_) \
                    - (np.outer(Hinv @ y_, s_) + np.outer(s_, y_ @ Hinv)) / sy
        if kind == "DIIS":
            rr = np.linalg.norm(gt) / 3.0
            mu = -0.2 if rr > 1e-2 else -0.3 if rr > 1e-3 else -0.4 if rr > 1e-4 else -0.5 if rr > 1e-5 else -0.6 if rr > 1e-6 else -0.7
            xn = xt + mu * gt
        elif kind == "BFGS":
            xn = xt - Hinv @ gt
        elif kind == "RF":
            aug = np.zeros((n + 1, n + 1)); aug[:n, :n] = Hm; aug[:n, n] = gt; aug[n, :n] = gt
            w_, V = np.linalg.eigh(aug)
            xn = xt + V[:n, 0] / V[n, 0]
        elif kind == "EF":
            ev, U = np.linalg.eigh(Hm)
            ev = np.maximum(ev, 0.005)
            xn = xt - U @ ((U.T @ gt) / ev)
        x = cap(xn, x)
    return maxit


for kind in ("DIIS", "BFGS", "RF", "EF"):
    it = run_opt(kind)
    check(f"optimizer {kind}: converges on a 9-dim quadratic (steps={it})", float(it >= 300), 0.0)

# steepest descent with the adaptive prefactor
x = np.zeros(n)
BohrR = 0.529177249
s_init = (2 * fs) ** 2 / amu
smin, smax = 0.04 * s_init, 15 * s_init
sc = min(max(1.0 / (1.5 * BohrR), smin), smax)
E_hist = []
for it in range(1, 3001):
    g = grad(x)
    gmax = np.max(np.abs(g))
    if gmax < 3e-4:
        break
    if it > 1:
        E = energy(x)
        if len(E_hist) >= 1 and E > E_hist[-1]:
            sc /= 2
        elif len(E_hist) >= 2 and E_hist[-1] < E_hist[-2] and E < E_hist[-1] and it % 4 == 1:
            sc *= 2.5
        sc = min(max(sc, smin), smax)
        E_hist.append(E)
    stepv = sc * g
    mx = np.max(np.abs(stepv))
    x = x - stepv * (CAP / mx if mx > CAP else 1.0)
check(f"optimizer SD: converges (steps={it})", float(it >= 3000), 0.0)


# ---- 8. model Hessian: translational invariance -----------------------------
pos = rng.normal(size=(5, 3)) * 2.0
Bpar = 1.0
Hm = np.zeros((15, 15))
for I in range(5):
    for J in range(5):
        if I == J:
            continue
        r = pos[I] - pos[J]; d = np.linalg.norm(r); u = r / d
        k = 1.734 / (d - Bpar) ** 3 if d > Bpar + 0.5 else 1.734 / 0.5**3
        blk = k * np.outer(u, u)
        Hm[3*I:3*I+3, 3*I:3*I+3] += blk
        Hm[3*I:3*I+3, 3*J:3*J+3] -= blk
t = np.tile(rng.normal(size=3), 5)
check("model Hessian (pair part): H t = 0 for a uniform translation", np.max(np.abs(Hm @ t)), 1e-12)
check("model Hessian: symmetric", np.max(np.abs(Hm - Hm.T)), 1e-12)


# ---- 9. Parrinello-Rahman geometry ------------------------------------------
h = rng.normal(size=(3, 3)) + 3 * np.eye(3)    # rows = lattice vectors
Om = np.linalg.det(h)
Sig = np.array([np.cross(h[1], h[2]), np.cross(h[2], h[0]), np.cross(h[0], h[1])])
check("Sigma = Omega h^{-T}", np.max(np.abs(Sig - Om * np.linalg.inv(h).T)), 1e-10)
eps = 1e-6
dOm = np.zeros((3, 3))
for i in range(3):
    for j in range(3):
        hp = h.copy(); hp[i, j] += eps
        hm = h.copy(); hm[i, j] -= eps
        dOm[i, j] = (np.linalg.det(hp) - np.linalg.det(hm)) / (2 * eps)
check("d Omega / d h_ij = Sigma_ij", np.max(np.abs(dOm - Sig)), 1e-7)


def sig_of(h):
    return np.array([np.cross(h[1], h[2]), np.cross(h[2], h[0]), np.cross(h[0], h[1])])


err = 0.0
for i in range(3):
    for a in range(3):
        hp = h.copy(); hp[i, a] += eps
        hm = h.copy(); hm[i, a] -= eps
        num = (sig_of(hp) - sig_of(hm)) / (2 * eps)
        ana = (Sig[i, a] * Sig - np.outer(Sig[:, a], Sig[i])) / Om
        err = max(err, np.max(np.abs(num - ana)))
check("d sigma_j / d a_{i alpha} = (sigma_{i alpha} sigma_j - sigma_{j alpha} sigma_i)/Omega", err, 1e-7)


# ---- 10. NEB on a 2D model ---------------------------------------------------
def V2(p):
    x, y = p
    return (x * x - 1) ** 2 + 2.0 * (y - 0.3 * np.sin(np.pi * x)) ** 2
def dV2(p):
    x, y = p
    dy = 4.0 * (y - 0.3 * np.sin(np.pi * x))
    return np.array([4 * x * (x * x - 1) - dy * 0.3 * np.pi * np.cos(np.pi * x), dy])


def tangent(Rm, R, Rp, Vm, V, Vp):
    tp, tm = Rp - R, R - Rm
    if Vm < V < Vp:
        t = tp
    elif Vp < V < Vm:
        t = tm
    else:
        dmax = max(abs(Vp - V), abs(Vm - V)); dmin = min(abs(Vp - V), abs(Vm - V))
        t = dmax * tp + dmin * tm if Vm <= Vp else dmin * tp + dmax * tm
    return t / np.linalg.norm(t)


def neb_grad(Rs, k):
    Vs = [V2(r) for r in Rs]
    out = np.zeros_like(Rs)
    for p in range(1, len(Rs) - 1):
        t = tangent(Rs[p-1], Rs[p], Rs[p+1], Vs[p-1], Vs[p], Vs[p+1])
        g = dV2(Rs[p])
        gp = g - (g @ t) * t
        d1 = np.linalg.norm(Rs[p+1] - Rs[p]); d2 = np.linalg.norm(Rs[p] - Rs[p-1])
        out[p] = gp - k * (d1 - d2) * t
    return out


P = 9
Rs = np.array([[-1.0 + 2.0 * i / (P + 1), 0.0] for i in range(P + 2)])
Rs[0] = [-1.0, 0.3 * np.sin(-np.pi)]; Rs[-1] = [1.0, 0.0]
k = 0.4
hist_x, hist_g = [], []
Hinv = np.eye(2 * P)
converged = False
for it in range(1, 400):
    G = neb_grad(Rs, k)
    maxF = max(np.linalg.norm(G[p]) for p in range(1, P + 1))
    if maxF < 1e-4:
        converged = True
        break
    xv = Rs[1:-1].ravel().copy(); gv = G[1:-1].ravel().copy()
    hist_x.insert(0, xv); hist_g.insert(0, gv)
    hist_x, hist_g = hist_x[:3], hist_g[:3]
    kk = len(hist_x)
    Bm = np.array([[a @ b for b in hist_g] for a in hist_g])
    Mx = np.zeros((kk + 1, kk + 1)); Mx[:kk, :kk] = Bm / np.max(np.abs(Bm)); Mx[:kk, kk] = 1; Mx[kk, :kk] = 1
    rhs = np.zeros(kk + 1); rhs[kk] = 1
    cc = np.linalg.solve(Mx, rhs)[:kk]
    xt = sum(c_ * a for c_, a in zip(cc, hist_x)); gt = sum(c_ * a for c_, a in zip(cc, hist_g))
    if kk > 1:
        s_ = hist_x[0] - hist_x[1]; y_ = hist_g[0] - hist_g[1]; sy = s_ @ y_
        if sy > 1e-14:
            Hinv = Hinv + (sy + y_ @ Hinv @ y_) / sy**2 * np.outer(s_, s_) - (np.outer(Hinv @ y_, s_) + np.outer(s_, y_ @ Hinv)) / sy
    f = 0.1 if maxF > 0.7 else 0.2 if maxF > 0.3 else 0.3 if maxF > 0.1 else 0.4 if maxF > 0.05 else 0.5 if maxF > 0.02 else 0.6 if maxF > 0.01 else 0.7 if maxF > 0.005 else 0.8 if maxF > 0.002 else 0.9 if maxF > 0.001 else 1.0
    xn = xt - f * Hinv @ gt
    d = xn - xv; mx = np.max(np.abs(d))
    if mx > 0.2:
        xn = xv + d * 0.2 / mx
    Rs[1:-1] = xn.reshape(P, 2)
check(f"NEB DIIS+BFGS: converged (iterations={it})", float(not converged), 0.0)
dist = np.array([np.linalg.norm(Rs[p] - Rs[p-1]) for p in range(1, P + 2)])
check("NEB: springs equalise the image spacing (std/mean)", np.std(dist[1:-1]) / np.mean(dist[1:-1]), 2e-2)
check("NEB: highest image energy close to the saddle (1.0)", abs(max(V2(r) for r in Rs) - 1.0), 0.05)


# ---- 11. wrap ----------------------------------------------------------------
def wrap(f):
    out = f.copy()
    for i in range(3):
        if out[i] > 1.0:
            out[i] -= np.floor(out[i])
        elif out[i] < 0.0:
            out[i] += np.floor(-out[i]) + 1.0
    return out


f = np.array([1.25, -0.25, 1.0])
check("wrap: 1.25 -> 0.25, -0.25 -> 0.75, 1.0 kept", np.max(np.abs(wrap(f) - np.array([0.25, 0.75, 1.0]))), 1e-15)

print("ALL OK" if ok else "SOME CHECKS FAILED")
