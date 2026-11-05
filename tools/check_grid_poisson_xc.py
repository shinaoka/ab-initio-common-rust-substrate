# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/grid-poisson-xc.md:
#   1. grid-size rule (nearest 7-smooth number in log scale) and cluster box rule
#   2. FFT Poisson solver: plane-wave density (periodic) and Gaussian charge
#      with the spherical Coulomb cutoff (isolated)
#   3. 27-point quadratic-fit gradient on a skewed grid
#   4. LDA (PZ81, PW92 + Slater) and PBE written from the papers vs libxc
#   5. non-collinear local frame: V = V0 + dV m_hat.sigma
# Run: uv run --with numpy --with scipy --with pyscf python tools/check_grid_poisson_xc.py

import math
import itertools
import numpy as np
from scipy.special import erf

ok_all = True
def check(name, err, tol):
    global ok_all
    flag = err <= tol
    ok_all &= flag
    print(f"[{'OK ' if flag else 'NG '}] {name}: err={err:.3e} tol={tol:.1e}")

# ---------------------------------------------------------------- 1. grid size
def smooth_numbers(nmax):
    out = []
    for n in range(1, nmax + 1):
        m = n
        for p in (2, 3, 5, 7):
            while m % p == 0:
                m //= p
        if m == 1:
            out.append(n)
    return out
SMOOTH = smooth_numbers(5000)

def grid_count(E_ry, a):
    """a: 3x3 rows = lattice vectors (bohr). returns N_i and used cutoffs (Ry)."""
    V = np.dot(a[0], np.cross(a[1], a[2]))
    b = 2 * np.pi * np.array([np.cross(a[1], a[2]), np.cross(a[2], a[0]), np.cross(a[0], a[1])]) / V
    N = []
    for i in range(3):
        target = 2.0 * math.sqrt(E_ry) / np.linalg.norm(b[i])
        N.append(min(SMOOTH, key=lambda n: abs(math.log(n) - math.log(target))))
    used = [(N[i] * np.linalg.norm(b[i]) / 2) ** 2 for i in range(3)]
    return N, used

print("grid examples (Ry -> N, used cutoff Ry):")
for E, a in [(150, np.eye(3) * 10.0), (220, np.eye(3) * 10.0), (200, np.diag([5.0, 7.3, 21.0])),
             (220, 5.4 * 0.5 * np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))]:
    N, used = grid_count(E, a)
    print("  E=%g  N=%s  used=%s" % (E, N, ["%.2f" % u for u in used]))

def cluster_box(extent, E_ry):
    """extent: (max - min) incl. +-rc, per Cartesian axis (bohr)."""
    h = math.pi / math.sqrt(E_ry + 0.001)
    L, N = [], []
    for e in extent:
        Lr = 1.10 * e
        j = int(Lr / h)
        while not (j in SMOOTH and h * j >= Lr):
            j += 1
        L.append(h * j); N.append(j)
    return L, N
L, Nc = cluster_box([12.0, 12.0, 14.3], 150.0)
Ng, _ = grid_count(150.0, np.diag(L))
print("  cluster box L=%s N(box)=%s N(rule)=%s" % (["%.4f" % x for x in L], Nc, Ng))
check("cluster box grid count consistent with general rule", float(np.max(np.abs(np.array(Nc) - Ng))), 0)

# ------------------------------------------------------------- 2. Poisson FFT
def poisson_fft(dn, a, kernel="periodic", Rc=None):
    Ns = dn.shape
    V = abs(np.dot(a[0], np.cross(a[1], a[2])))
    b = 2 * np.pi * np.array([np.cross(a[1], a[2]), np.cross(a[2], a[0]), np.cross(a[0], a[1])]) / np.dot(a[0], np.cross(a[1], a[2]))
    m = [np.fft.fftfreq(n, 1.0 / n) for n in Ns]          # 0..N/2-1, -N/2..-1
    # (for even N, fftfreq puts -N/2; same as the convention k<N/2 -> k else k-N)
    M1, M2, M3 = np.meshgrid(*m, indexing="ij")
    G = M1[..., None] * b[0] + M2[..., None] * b[1] + M3[..., None] * b[2]
    G2 = np.sum(G * G, axis=-1)
    nG = np.fft.fftn(dn) / dn.size
    with np.errstate(divide="ignore", invalid="ignore"):
        if kernel == "periodic":
            K = 4 * np.pi / G2
            K[0, 0, 0] = 0.0
        else:
            K = 4 * np.pi / G2 * (1 - np.cos(np.sqrt(G2) * Rc))
            K[0, 0, 0] = 2 * np.pi * Rc ** 2
    return np.real(np.fft.ifftn(K * nG) * dn.size), V

# periodic: dn = cos(G.r) in a skewed cell
a = np.array([[8.0, 0, 0], [2.0, 7.0, 0], [1.0, 1.5, 9.0]])
Ns = (24, 25, 27)
fr = [np.arange(n) / n for n in Ns]
F1, F2, F3 = np.meshgrid(*fr, indexing="ij")
R = F1[..., None] * a[0] + F2[..., None] * a[1] + F3[..., None] * a[2]
V = np.dot(a[0], np.cross(a[1], a[2]))
b = 2 * np.pi * np.array([np.cross(a[1], a[2]), np.cross(a[2], a[0]), np.cross(a[0], a[1])]) / V
Gv = 2 * b[0] - 1 * b[1] + 3 * b[2]
dn = np.cos(R @ Gv) + 0.3
VH, _ = poisson_fft(dn, a)
check("periodic Poisson, cos(G.r)+const (q=0 dropped)", np.max(np.abs(VH - 4 * np.pi / (Gv @ Gv) * np.cos(R @ Gv))), 1e-12)

# isolated: Gaussian charge with the spherical cutoff Rc = L/2 (cubic cell)
Lc, n = 24.0, 96
a = np.eye(3) * Lc
x = np.arange(n) * Lc / n
X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
c = Lc / 2
r = np.sqrt((X - c) ** 2 + (Y - c) ** 2 + (Z - c) ** 2)
alpha = 1.5
rho = (alpha / np.pi) ** 1.5 * np.exp(-alpha * r ** 2)       # unit charge
VH, Vol = poisson_fft(rho, a, kernel="cutoff", Rc=Lc / 2)
dV = Vol / rho.size
with np.errstate(invalid="ignore", divide="ignore"):
    Vex = np.where(r > 1e-12, erf(np.sqrt(alpha) * r) / r, 2 * np.sqrt(alpha / np.pi))
mask = r < Lc / 2 - 4.0
check("Gaussian + spherical Coulomb cutoff: potential (r < L/2 - 4)", np.max(np.abs(VH - Vex)[mask]), 1e-8)
EH = 0.5 * dV * np.sum(rho * VH)
check("Gaussian self energy sqrt(alpha/(2 pi))", abs(EH - math.sqrt(alpha / (2 * np.pi))), 1e-8)
print("  charge on grid dV*sum(rho) - 1 = %.2e" % (dV * rho.sum() - 1))
VHp, _ = poisson_fft(rho - rho.mean(), a)
print("  periodic (no cutoff) Gaussian energy = %.6f, isolated = %.6f" % (0.5 * dV * np.sum(rho * VHp), math.sqrt(alpha / (2 * np.pi))))

# odd N: index kappa -> k = kappa if kappa < floor(N/2) else kappa - N, real part taken
def poisson_mapped(dn, a):
    Ns = dn.shape
    Vc = np.dot(a[0], np.cross(a[1], a[2]))
    bb = 2 * np.pi * np.array([np.cross(a[1], a[2]), np.cross(a[2], a[0]), np.cross(a[0], a[1])]) / Vc
    m = [np.array([k if k < n // 2 else k - n for k in range(n)], float) for n in Ns]
    M1, M2, M3 = np.meshgrid(*m, indexing="ij")
    G = M1[..., None] * bb[0] + M2[..., None] * bb[1] + M3[..., None] * bb[2]
    G2 = np.sum(G * G, axis=-1)
    with np.errstate(divide="ignore"):
        K = 4 * np.pi / G2
    K[0, 0, 0] = 0.0
    return np.real(np.fft.ifftn(K * np.fft.fftn(dn)))
for n_odd, expect_equal in ((16, True), (15, False)):
    a = np.eye(3) * 10.0
    x = np.arange(n_odd) * 10.0 / n_odd
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    g = np.exp(-3.0 * ((X - 5) ** 2 + (Y - 5) ** 2 + (Z - 5) ** 2)); g -= g.mean()
    v_sym, _ = poisson_fft(g, a)
    v_map = poisson_mapped(g, a)
    dE = 0.5 * (1000.0 / g.size) * np.sum((v_sym - v_map) * g)
    if expect_equal:
        check("N=16: mapped frequencies = symmetric choice", np.max(np.abs(v_sym - v_map)), 1e-12)
    else:
        print("  N=15: E_dH(symmetric) - E_dH(mapped) = %.2e Ha" % dE)

# ------------------------------------------------- 3. 27-point gradient fit
def grad_fit_weights(g):
    rows, disp = [], []
    for i, j, k in itertools.product((-1, 0, 1), repeat=3):
        d = i * g[0] + j * g[1] + k * g[2]
        disp.append((i, j, k))
        rows.append([d[0] ** p * d[1] ** q * d[2] ** s for p in range(3) for q in range(3) for s in range(3)])
    A = np.array(rows)
    inv = np.linalg.inv(A)              # coefficients = inv @ f
    # monomial index p*9+q*3+s : x -> (1,0,0)=9, y -> 3, z -> 1
    return disp, inv[[9, 3, 1], :]
gsk = np.array([[0.20, 0.0, 0.0], [0.07, 0.19, 0.0], [0.03, 0.05, 0.21]])
disp, W = grad_fit_weights(gsk)
rng = np.random.default_rng(0)
cq = rng.normal(size=10)
def poly2(r):
    x, y, z = r
    return (cq[0] + cq[1] * x + cq[2] * y + cq[3] * z + cq[4] * x * x + cq[5] * y * y + cq[6] * z * z
            + cq[7] * x * y + cq[8] * y * z + cq[9] * z * x)
def grad2(r):
    x, y, z = r
    return np.array([cq[1] + 2 * cq[4] * x + cq[7] * y + cq[9] * z,
                     cq[2] + 2 * cq[5] * y + cq[7] * x + cq[8] * z,
                     cq[3] + 2 * cq[6] * z + cq[8] * y + cq[9] * x])
r0 = np.array([0.3, -0.2, 0.5])
f = np.array([poly2(r0 + i * gsk[0] + j * gsk[1] + k * gsk[2]) for (i, j, k) in disp])
check("27-point fit exact for total degree <= 2 (skewed grid)", np.max(np.abs(W @ f - grad2(r0))), 1e-10)
_, Wo = grad_fit_weights(np.diag([0.2, 0.25, 0.3]))
cd = np.zeros(27); cd[disp.index((1, 0, 0))] = 1 / 0.4; cd[disp.index((-1, 0, 0))] = -1 / 0.4
check("orthogonal grid: x-derivative = central difference", np.max(np.abs(Wo[0] - cd)), 1e-10)
fg = np.array([math.exp(-0.7 * np.sum((r0 + i * gsk[0] + j * gsk[1] + k * gsk[2]) ** 2)) for (i, j, k) in disp])
ge = -1.4 * r0 * math.exp(-0.7 * r0 @ r0)
print("  Gaussian gradient rel. error (h~0.2 bohr): %.2e" % (np.linalg.norm(W @ fg - ge) / np.linalg.norm(ge)))

# ---------------------------------------------------------------- 4. XC
CRS = (3 / (4 * np.pi)) ** (1 / 3)
def rs_of(n): return CRS * n ** (-1 / 3)
def ex_unif(n):           # Slater exchange per particle, unpolarized
    return -0.75 * (3 / np.pi) ** (1 / 3) * n ** (1 / 3)
def fz(z): return ((1 + z) ** (4 / 3) + (1 - z) ** (4 / 3) - 2) / (2 ** (4 / 3) - 2)
def ex_pol(nu, nd):       # spin scaling
    n = nu + nd
    return (nu * ex_unif(2 * nu) + nd * ex_unif(2 * nd)) / n
def ec_pz_single(rs, pol):
    g, b1, b2, A, B, C, D = ((-0.1423, 1.0529, 0.3334, 0.0311, -0.048, 0.0020, -0.0116) if pol == 0
                             else (-0.0843, 1.3981, 0.2611, 0.01555, -0.0269, 0.0007, -0.0048))
    return np.where(rs >= 1, g / (1 + b1 * np.sqrt(rs) + b2 * rs), A * np.log(rs) + B + C * rs * np.log(rs) + D * rs)
def ec_pz(nu, nd):
    n = nu + nd; z = (nu - nd) / n; rs = rs_of(n)
    return ec_pz_single(rs, 0) + fz(z) * (ec_pz_single(rs, 1) - ec_pz_single(rs, 0))
PW = {0: (0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294),
      1: (0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517),
      2: (0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671)}
def G_pw(rs, i):
    A, a1, b1, b2, b3, b4 = PW[i]
    s = np.sqrt(rs)
    return -2 * A * (1 + a1 * rs) * np.log(1 + 1 / (2 * A * (b1 * s + b2 * rs + b3 * rs * s + b4 * rs * rs)))
def ec_pw(nu, nd):
    n = nu + nd; z = (nu - nd) / n; rs = rs_of(n)
    fpp0 = 1.709921
    ec0, ec1, mac = G_pw(rs, 0), G_pw(rs, 1), -G_pw(rs, 2)
    return ec0 + mac * fz(z) / fpp0 * (1 - z ** 4) + (ec1 - ec0) * fz(z) * z ** 4
def pbe(nu, nd, gu, gd):   # gu, gd: gradient vectors (3, ...) ; returns exc per particle
    kappa, mu = 0.804, 0.2195149727645171
    beta = 0.06672455060314922; gam = (1 - math.log(2)) / np.pi ** 2
    def ex_s(n2, g2):   # exchange per particle for unpolarized density n2 with |grad| g2
        kf = (3 * np.pi ** 2 * n2) ** (1 / 3); s = g2 / (2 * kf * n2)
        return ex_unif(n2) * (1 + kappa - kappa / (1 + mu * s * s / kappa))
    n = nu + nd
    ex = (nu * ex_s(2 * nu, 2 * np.linalg.norm(gu, axis=0)) + nd * ex_s(2 * nd, 2 * np.linalg.norm(gd, axis=0))) / n
    z = (nu - nd) / n; rs = rs_of(n)
    phi = 0.5 * ((1 + z) ** (2 / 3) + (1 - z) ** (2 / 3))
    kf = (3 * np.pi ** 2 * n) ** (1 / 3); ks = np.sqrt(4 * kf / np.pi)
    t = np.linalg.norm(gu + gd, axis=0) / (2 * phi * ks * n)
    ecu = ec_pw(nu, nd)
    A = beta / gam / (np.exp(-ecu / (gam * phi ** 3)) - 1)
    H = gam * phi ** 3 * np.log(1 + beta / gam * t * t * (1 + A * t * t) / (1 + A * t * t + A * A * t ** 4))
    return ex + ecu + H

try:
    from pyscf.dft import libxc
    have_libxc = True
except Exception:
    have_libxc = False
    print("pyscf/libxc not available; skipping XC comparison")

if have_libxc:
    nlist = np.array([1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0])
    # unpolarized
    for name, mine in [("LDA,PZ", lambda n: ex_unif(n) + ec_pz(n / 2, n / 2)),
                       ("LDA,PW", lambda n: ex_unif(n) + ec_pw(n / 2, n / 2))]:
        e = libxc.eval_xc(name, nlist, spin=0, deriv=1)[0]
        check(f"{name} unpolarized eps_xc vs libxc", np.max(np.abs(mine(nlist) - e) / np.abs(e)), 1e-6)
    # polarized
    nu = np.array([0.08, 0.5, 0.003, 2.0]); nd = np.array([0.02, 0.5, 0.0005, 0.7])
    for name, mine in [("LDA,PZ", lambda u, d: ex_pol(u, d) + ec_pz(u, d)),
                       ("LDA,PW", lambda u, d: ex_pol(u, d) + ec_pw(u, d))]:
        e = libxc.eval_xc(name, np.array([nu, nd]), spin=1, deriv=1)[0]
        check(f"{name} polarized eps_xc vs libxc", np.max(np.abs(mine(nu, nd) - e) / np.abs(e)), 1e-6)
    # PBE (spin polarized, gradients along x only for simplicity + a general one)
    gu = np.array([[0.05, 0.3, 0.001, 1.0], [0.0, 0.1, 0.002, 0.2], [0.0, 0.0, 0.0, -0.4]])
    gd = np.array([[0.01, 0.2, 0.0001, 0.3], [0.0, -0.1, 0.0, 0.0], [0.0, 0.05, 0.0002, 0.1]])
    rho = np.zeros((2, 4, 4)); rho[0, 0] = nu; rho[1, 0] = nd; rho[0, 1:] = gu; rho[1, 1:] = gd
    e = libxc.eval_xc("PBE,PBE", rho, spin=1, deriv=1)[0]
    check("PBE polarized eps_xc vs libxc", np.max(np.abs(pbe(nu, nd, gu, gd) - e) / np.abs(e)), 1e-6)
    # reference table for the chapter (unpolarized, Hartree)
    print("reference values (libxc %s), unpolarized:" % libxc.__version__)
    for name in ("LDA,PZ", "LDA,PW"):
        e, v = libxc.eval_xc(name, nlist, spin=0, deriv=1)[:2]
        for ni, ei, vi in zip(nlist, e, v[0]):
            print(f"  {name:7s} n={ni:8.1e}  eps={ei:.10f}  v={vi:.10f}")
    rho = np.zeros((4, 3)); rho[0] = [0.01, 0.1, 1.0]; rho[1] = [0.01, 0.1, 0.5]
    e, v = libxc.eval_xc("PBE,PBE", rho, spin=0, deriv=1)[:2]
    for k in range(3):
        print(f"  PBE     n={rho[0,k]:.2f} |grad n|={rho[1,k]:.2f}  eps={e[k]:.10f}  vrho={v[0][k]:.10f}  vsigma={v[1][k]:.10f}")
    rho = np.array([[0.08], [0.02]])
    e, v = libxc.eval_xc("LDA,PZ", rho, spin=1, deriv=1)[:2]
    print(f"  LDA,PZ  n_up=0.08 n_dn=0.02  eps={e[0]:.10f}  v_up={v[0][0,0]:.10f}  v_dn={v[0][0,1]:.10f}")
    e, v = libxc.eval_xc("LDA,PW", rho, spin=1, deriv=1)[:2]
    print(f"  LDA,PW  n_up=0.08 n_dn=0.02  eps={e[0]:.10f}  v_up={v[0][0,0]:.10f}  v_dn={v[0][0,1]:.10f}")

# ------------------------------------------------- 5. non-collinear frame
s0 = np.eye(2); sx = np.array([[0, 1], [1, 0]]); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0, -1.0])
n_tot, m = 0.7, np.array([0.1, -0.2, 0.15])
rho2 = 0.5 * (n_tot * s0 + m[0] * sx + m[1] * sy + m[2] * sz)
mx, my, mz = 2 * rho2[0, 1].real, -2 * rho2[0, 1].imag, (rho2[0, 0] - rho2[1, 1]).real
check("m recovered from rho_12 = (m_x - i m_y)/2", np.max(np.abs(np.array([mx, my, mz]) - m)), 1e-14)
ev = np.linalg.eigvalsh(rho2)
am = np.linalg.norm(m)
check("local densities = (n +- |m|)/2", np.max(np.abs(np.sort(ev) - np.array([(n_tot - am) / 2, (n_tot + am) / 2]))), 1e-14)
vu, vd = -0.9, -0.7
mh = m / am
Vmat = 0.5 * (vu + vd) * s0 + 0.5 * (vu - vd) * (mh[0] * sx + mh[1] * sy + mh[2] * sz)
th = math.acos(mh[2]); ph = math.atan2(mh[1], mh[0])
dv = 0.5 * (vu - vd)
V12 = dv * math.sin(th) * (math.cos(ph) - 1j * math.sin(ph))
check("V_12 = dV sin(theta) e^{-i phi}", abs(Vmat[0, 1] - V12), 1e-14)
# E_xc derivative consistency: Tr(V d rho) = vu d n_up + vd d n_dn for a small change
d = 1e-6 * np.array([[0.3, 0.1 - 0.2j], [0.1 + 0.2j, -0.1]])
def loc(r):
    e = np.linalg.eigvalsh(r); return e[1], e[0]
u1, d1 = loc(rho2 + d); u0, d0 = loc(rho2)
check("first-order: Tr(V drho) = v_up dn_up + v_dn dn_dn", abs(np.trace(Vmat @ d).real - (vu * (u1 - u0) + vd * (d1 - d0))), 1e-12)

print("ALL OK" if ok_all else "SOME CHECKS FAILED")
