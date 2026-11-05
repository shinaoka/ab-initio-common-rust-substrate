# SPDX-License-Identifier: MIT
# Numerical checks for src/modules/negf.md (NEGF electron transport).
# A small non-orthogonal tight-binding model (lead cell with 2 orbitals,
# central region = a few lead cells with a modified on-site block) is used to
# check: the pole expansion of the Fermi function, the surface Green's function
# by decimation, the self-energy, the equilibrium density matrix by poles,
# the non-equilibrium correction, transmission (Caroli formula), eigenchannels,
# and unit conversions.
# Run: uv run --with numpy --with scipy python tools/check_negf.py

import numpy as np
import scipy.linalg as sla

rng = np.random.default_rng(7)
ok_all = True


def check(name, val, tol):
    global ok_all
    flag = val <= tol
    ok_all &= flag
    print(f"{'OK ' if flag else 'NG '} {name}: {val:.3e} (tol {tol:.0e})")


# ---------------------------------------------------------------------------
# 1. Pole expansion of the Fermi function (continued fraction, tech note 4)
# ---------------------------------------------------------------------------
def fermi_poles(n):
    """Return z_p > 0 and R_p (p = 1..n) such that
    1/(1+e^x) ~ 1/2 + sum_p R_p [1/(x - i z_p) + 1/(x + i z_p)].
    Derivation: truncated continued fraction of depth M = 2n gives
    f = 1/2 - (x/4) [(D + i x B)^{-1}]_{11}, D = diag(1,3,...,2M-1),
    B tridiagonal with 1/2 on the off-diagonals.  With B v = lam D v,
    v^T D v = 1, the poles are x = +- i/lam, and R = -(v_1^2 / 4) z^2."""
    M = 2 * n
    D = np.diag(2.0 * np.arange(1, M + 1) - 1.0)
    B = np.zeros((M, M))
    for i in range(M - 1):
        B[i, i + 1] = B[i + 1, i] = 0.5
    lam, V = sla.eigh(B, D)
    pos = lam > 0
    lam, V = lam[pos], V[:, pos]
    z = 1.0 / lam
    R = -0.25 * V[0, :] ** 2 * z ** 2
    order = np.argsort(z)
    return z[order], R[order]


def fermi_approx(x, z, R):
    x = np.asarray(x, dtype=complex)
    s = 0.5 + 0 * x
    for zp, rp in zip(z, R):
        s += rp * (1.0 / (x - 1j * zp) + 1.0 / (x + 1j * zp))
    return s


x = np.linspace(-60, 60, 2001)
fex = 1.0 / (1.0 + np.exp(x))
for n in (20, 40, 80):
    z, R = fermi_poles(n)
    err = np.max(np.abs(fermi_approx(x, z, R) - fex))
    print(f"    poles N={n}: max |f_approx - f| on |x|<=60 = {err:.2e}, "
          f"z_1={z[0]:.6f} (pi={np.pi:.6f}), z_N={z[-1]:.3e}")
z, R = fermi_poles(150)
check("Fermi pole expansion, N=150, |x|<=60", np.max(np.abs(fermi_approx(x, z, R) - fex)), 1e-12)
# the first poles coincide with Matsubara poles (2p-1) pi
check("first 20 poles = Matsubara (2p-1)pi (rel)", np.max(np.abs(z[:20] / ((2 * np.arange(1, 21) - 1) * np.pi) - 1)), 1e-10)
check("first 20 residues = -1 (rel)", np.max(np.abs(R[:20] + 1.0)), 1e-10)

# ---------------------------------------------------------------------------
# 2. Model: lead cell (2 orbitals), nonorthogonal
# ---------------------------------------------------------------------------
n = 2
H00 = np.array([[-0.20, -0.30], [-0.30, -0.10]])
H01 = np.array([[0.02, 0.00], [-0.25, 0.03]])      # coupling cell j -> j+1 (to the right)
S00 = np.array([[1.0, 0.10], [0.10, 1.0]])
S01 = np.array([[0.00, 0.00], [0.08, 0.01]])
H10, S10 = H01.T, S01.T                             # coupling cell j -> j-1


def surface_green(zz, h00, s00, h_out, s_out, h_in, s_in, eps=1e-12, itmax=600):
    """Decimation (renormalisation) for the surface cell of a semi-infinite lead.
    h_out: coupling surface cell -> next cell away from C, h_in: reverse."""
    e_s = zz * s00 - h00
    e_b = e_s.copy()
    a = -(zz * s_out - h_out)
    b = -(zz * s_in - h_in)
    g_old = np.linalg.inv(e_s)
    for it in range(1, itmax):
        ge = np.linalg.inv(e_b)
        agb = a @ ge @ b
        bga = b @ ge @ a
        e_s = e_s - agb
        e_b = e_b - agb - bga
        a = a @ ge @ a
        b = b @ ge @ b
        g = np.linalg.inv(e_s)
        if np.max(np.abs(g - g_old)) < eps:
            return g, it
        g_old = g
    return g, itmax


# Dyson self-consistency of the surface Green's function
for E in (-0.6, -0.15, 0.3):
    zz = E + 1e-6j
    gR, it = surface_green(zz, H00, S00, H01, S01, H10, S10)
    res = np.linalg.inv(zz * S00 - H00 - (zz * S01 - H01) @ gR @ (zz * S10 - H10)) - gR
    check(f"right-lead surface G Dyson residual E={E:+.2f} ({it} it)", np.max(np.abs(res)), 1e-9)

# scalar chain: analytic surface Green function
eps0, t = 0.0, -0.5
for E in (-0.7, 0.2, 1.5):
    zz = E + 1e-9j
    g, it = surface_green(zz, np.array([[eps0]]), np.eye(1), np.array([[t]]), np.zeros((1, 1)),
                          np.array([[t]]), np.zeros((1, 1)))
    w = zz - eps0
    ga = (w - np.sqrt(w - 2 * t) * np.sqrt(w + 2 * t)) / (2 * t * t)
    check(f"scalar chain g_s analytic E={E:+.2f}", abs(g[0, 0] - ga), 1e-7)

# ---------------------------------------------------------------------------
# 3. Extended central region C = L0 | C0 | R0 (5 lead cells; middle one perturbed)
# ---------------------------------------------------------------------------
ncell = 5
N = n * ncell
HC = np.zeros((N, N)); SC = np.zeros((N, N))
for c in range(ncell):
    sl = slice(n * c, n * c + n)
    HC[sl, sl] = H00; SC[sl, sl] = S00
    if c + 1 < ncell:
        sr = slice(n * c + n, n * c + 2 * n)
        HC[sl, sr] = H01; SC[sl, sr] = S01
        HC[sr, sl] = H10; SC[sr, sl] = S10
mid = slice(n * 2, n * 3)
HC_pert = HC.copy(); HC_pert[mid, mid] += np.array([[0.15, 0.05], [0.05, -0.08]])
# couplings C -> L1 (only first cell couples) and C -> R1 (only last cell)
HCL = np.zeros((N, n)); SCL = np.zeros((N, n))
HCL[:n, :] = H10; SCL[:n, :] = S10
HCR = np.zeros((N, n)); SCR = np.zeros((N, n))
HCR[-n:, :] = H01; SCR[-n:, :] = S01


def self_energies(zz, shiftL=0.0, shiftR=0.0):
    # bias: lead H -> H + V S (applied to all lead blocks incl. couplings)
    gL, _ = surface_green(zz, H00 + shiftL * S00, S00, H10 + shiftL * S10, S10, H01 + shiftL * S01, S01)
    gR, _ = surface_green(zz, H00 + shiftR * S00, S00, H01 + shiftR * S01, S01, H10 + shiftR * S10, S10)
    tL = zz * SCL - (HCL + shiftL * SCL)
    tR = zz * SCR - (HCR + shiftR * SCR)
    SigL = tL @ gL @ (np.conj(zz) * SCL - (HCL + shiftL * SCL)).conj().T
    SigR = tR @ gR @ (np.conj(zz) * SCR - (HCR + shiftR * SCR)).conj().T
    return SigL, SigR


def greenC(zz, H, shiftL=0.0, shiftR=0.0):
    SigL, SigR = self_energies(zz, shiftL, shiftR)
    return np.linalg.inv(zz * SC - H - SigL - SigR), SigL, SigR


def transmission(E, H, eta=1e-6, shiftL=0.0, shiftR=0.0):
    Gr, SLr, SRr = greenC(E + 1j * eta, H, shiftL, shiftR)
    Ga, SLa, SRa = greenC(E - 1j * eta, H, shiftL, shiftR)
    GamL = 1j * (SLr - SLa); GamR = 1j * (SRr - SRa)
    return np.trace(GamL @ Gr @ GamR @ Ga).real, Gr, GamL, GamR


# perfect chain: T(E) = number of right-moving channels (from lead band structure)
def nchannels(E):
    ks = np.linspace(-np.pi, np.pi, 4001)
    bands = []
    for k in ks:
        Hk = H00 + H01 * np.exp(1j * k) + H10 * np.exp(-1j * k)
        Sk = S00 + S01 * np.exp(1j * k) + S10 * np.exp(-1j * k)
        bands.append(sla.eigh(Hk, Sk, eigvals_only=True))
    bands = np.array(bands)
    cnt = 0
    for b in range(n):
        s = np.sign(bands[:, b] - E)
        cnt += np.sum(s[1:] != s[:-1])
    return cnt // 2


errT = 0.0
for E in (-0.75, -0.45, -0.1, 0.25, 0.6):
    T, *_ = transmission(E, HC, eta=1e-8)
    errT = max(errT, abs(T - nchannels(E)))
    print(f"    perfect chain E={E:+.2f}: T={T:.8f}, channels={nchannels(E)}")
check("perfect chain T = channel count (eta=1e-8)", errT, 1e-5)
# a finite eta acts as an absorber inside C: 1 - T grows linearly with eta
for eta in (1e-4, 1e-5, 1e-6):
    print(f"    perfect chain E=-0.10, eta={eta:.0e}: 1-T = {1 - transmission(-0.1, HC, eta=eta)[0]:.3e}")

# transmission symmetry T_LR = T_RL and 0 <= eigenchannel <= 1; channel sum = T
E = -0.1
T, Gr, GamL, GamR = transmission(E, HC_pert)
T_RL = np.trace(GamR @ Gr @ GamL @ Gr.conj().T).real
check("T_LR = T_RL", abs(T - T_RL), 1e-10)
# eigenchannels (Loewdin basis): A_L = G Gam_L G^+, T = A_L Gam_R
lamS, US = np.linalg.eigh(SC)
Sh = US @ np.diag(np.sqrt(lamS)) @ US.T
Shi = US @ np.diag(1 / np.sqrt(lamS)) @ US.T
AL = Sh @ (Gr @ GamL @ Gr.conj().T) @ Sh
GR_o = Shi @ GamR @ Shi
a, Ua = np.linalg.eigh(AL)
AL_half = Ua @ np.diag(np.sqrt(np.clip(a, 0, None)))
Tt = AL_half.conj().T @ GR_o @ AL_half
tch = np.sort(np.linalg.eigvalsh(Tt))[::-1]
print(f"    perturbed chain E={E:+.2f}: T={T:.8f}, channels={tch[:3]}")
check("sum of channel eigenvalues = T", abs(tch.sum() - T), 1e-10)
check("channel eigenvalues within [0,1]", max(0.0, -tch.min(), tch.max() - 1), 1e-10)
# uniform shift: leads and C all shifted by V*S  =>  T(E + V) = T(E)
Vs = 0.07
T_shift, *_ = transmission(E + Vs, HC_pert + Vs * SC, shiftL=Vs, shiftR=Vs)
check("uniform shift H->H+VS: T(E+V) = T(E)", abs(T_shift - T), 1e-10)

# ---------------------------------------------------------------------------
# 4. Equilibrium density matrix by poles; zeroth moment
# ---------------------------------------------------------------------------
kT = 0.01  # Hartree (about 3158 K: keeps the real-axis comparison cheap)
beta = 1.0 / kT
mu = -0.1


def rho_poles(mu, H, npole=150, Rbig=1e12, shiftL=0.0, shiftR=0.0):
    z, Rp = fermi_poles(npole)
    G_iR, *_ = greenC(1j * Rbig, H, shiftL, shiftR)
    rho = 0.5 * (1j * Rbig) * G_iR
    for zp, rp in zip(z, Rp):
        Gp, *_ = greenC(mu + 1j * zp / beta, H, shiftL, shiftR)
        rho += (-rp / beta) * (Gp + Gp.conj().T)
    return rho


rho150 = rho_poles(mu, HC_pert)
rho300 = rho_poles(mu, HC_pert, npole=300)
check("rho(150 poles) vs rho(300 poles)", np.max(np.abs(rho150 - rho300)), 1e-11)
check("rho Hermitian", np.max(np.abs(rho150 - rho150.conj().T)), 1e-12)
# zeroth moment differs from S_C^{-1} when S_CL != 0
G_iR, *_ = greenC(1j * 1e12, HC_pert)
m0 = 1j * 1e12 * G_iR
print(f"    |iR G(iR) - S_C^-1|_max = {np.max(np.abs(m0 - np.linalg.inv(SC))):.3e} (nonzero because S_CL != 0)")

# isolated system (no leads): poles reproduce sum_n c_n f(e_n) c_n^+
Hiso = HC_pert[:4, :4]; Siso = SC[:4, :4]
ev, C = sla.eigh(Hiso, Siso)
rho_ex = C @ np.diag(1 / (1 + np.exp(beta * (ev - mu)))) @ C.T
z, Rp = fermi_poles(150)
rho_iso = 0.5 * np.linalg.inv(Siso).astype(complex)
for zp, rp in zip(z, Rp):
    Gp = np.linalg.inv((mu + 1j * zp / beta) * Siso - Hiso)
    rho_iso += (-rp / beta) * (Gp + Gp.conj().T)
check("isolated: pole rho = Fermi sum", np.max(np.abs(rho_iso - rho_ex)), 1e-12)

# real-axis reference: rho = -(1/pi) Im_m int G(E+i eta) f dE  (Im_m X = (X-X^+)/(2i))
eta = 2e-3
Es = np.arange(-1.6, mu + 40 * kT, 2e-4)
acc = np.zeros((N, N), dtype=complex)
for Ei in Es:
    Gi, *_ = greenC(Ei + 1j * eta, HC_pert)
    acc += (Gi - Gi.conj().T) / (2j) * (1 / (1 + np.exp(beta * (Ei - mu))))
rho_real = -acc * 2e-4 / np.pi
check("pole rho vs real-axis rho (eta=2e-3, O(eta) error)", np.max(np.abs(rho150 - rho_real)), 5e-3)

# ---------------------------------------------------------------------------
# 5. Non-equilibrium: reference-lead independence
#    rho = rho_eq(mu_low) + int G Gam_high G^+ (f_high - f_low) dE / 2pi
#        = rho_eq(mu_high) - int G Gam_low G^+ (f_high - f_low) dE / 2pi
# ---------------------------------------------------------------------------
V = 0.04
muL, muR = mu - V / 2, mu + V / 2
dE = 5e-4
eta_neq = 1e-4
Es = np.arange(mu - 0.4, mu + 0.4, dE)
accR = np.zeros((N, N), dtype=complex); accL = np.zeros((N, N), dtype=complex)
for Ei in Es:
    w = 1 / (1 + np.exp(beta * (Ei - muR))) - 1 / (1 + np.exp(beta * (Ei - muL)))
    if abs(w) < 1e-8:
        continue
    Gr_, SL_, SR_ = greenC(Ei + 1j * eta_neq, HC_pert, -V / 2, V / 2)
    GamL_ = 1j * (SL_ - SL_.conj().T); GamR_ = 1j * (SR_ - SR_.conj().T)
    accR += Gr_ @ GamR_ @ Gr_.conj().T * w
    accL += Gr_ @ GamL_ @ Gr_.conj().T * w
rho_a = rho_poles(muL, HC_pert, shiftL=-V / 2, shiftR=V / 2) + accR * dE / (2 * np.pi)
rho_b = rho_poles(muR, HC_pert, shiftL=-V / 2, shiftR=V / 2) - accL * dE / (2 * np.pi)
check("non-eq rho independent of reference lead", np.max(np.abs(rho_a - rho_b)), 2e-3)

# ---------------------------------------------------------------------------
# 6. Units
# ---------------------------------------------------------------------------
e = 1.602176634e-19; Eh = 4.3597447222071e-18; hbar = 1.054571817e-34; h = 2 * np.pi * hbar
I_au = e * Eh / hbar
print(f"    1 a.u. of current e*Eh/hbar = {I_au:.7e} A (implementation constant 6.6236178e-3)")
check("current unit", abs(I_au - 6.6236178e-3) / I_au, 1e-6)
G0 = e * e / h
print(f"    G0 = e^2/h = {G0:.6e} S ; 2e^2/h = {2*G0:.6e} S")

print("ALL OK" if ok_all else "SOME CHECKS FAILED")
