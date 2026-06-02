"""
SSEE Local Void — Hubble Tension via phi-DM structure suppression

Hypothesis: SSEE two-sector power spectrum predicts larger voids than LCDM.
If we inhabit a local underdensity delta ~ -0.46 within R ~ 200 Mpc/h,
the outflow velocity Delta_H0 = H0 * |delta| * f(Omega_m) / 3 ~ KAL ~ 5.52 km/s/Mpc
could explain the Hubble tension gap without new physics.

Method:
  1. Compute P(k) for LCDM and SSEE two-sector (CAMB + WDM mixing)
  2. Compute sigma(R) = top-hat rms fluctuation for R = 8..300 Mpc/h
  3. Compute P(delta < delta_thresh) under Gaussian approx
  4. Compute implied Delta_H0 from void depth
  5. Check whether KAL = beta + pi = 5.521 emerges naturally

Verification mandate: all intermediates printed before use.
"""

import numpy as np
import camb
from scipy.integrate import quad
from scipy.special import erfc

# ── algebraic constants ────────────────────────────────────────────────────────
phi   = (1 + 5**0.5) / 2
pi_   = np.pi
beta  = (phi + pi_) / 2                    # BIAL = 2.3798
KAL   = beta + pi_                         # 5.5214  h/Mpc → km/s/Mpc (same numeric)
H0_alg = 3 * (phi + pi_)**2               # 67.96

print("=" * 65)
print("SSEE Local Void — Hubble Tension test")
print("=" * 65)
print(f"\n  phi  = {phi:.10f}")
print(f"  pi   = {pi_:.10f}")
print(f"  beta = {beta:.10f}   (BIAL)")
print(f"  KAL  = {KAL:.10f}   (target Delta_H0 km/s/Mpc)")
print(f"  H0   = {H0_alg:.6f}  (algebraic)")
print(f"  KAL/H0 = {KAL/H0_alg:.6f}  (= fractional local excess)")

# ── SSEE parameters ────────────────────────────────────────────────────────────
w0_SSEE   = -0.840
wa_SSEE   = -0.670
OMM_MIRA  = 0.3199
OMB_H2    = 0.02237
ns_SSEE   = 1 - phi**(-7)
As_SSEE   = 2.1e-9
alpha_WDM = 1.6561    # calibrated h/Mpc
f_phiDM   = 0.5
k_fs      = 0.493

# ── LCDM Planck 2018 ───────────────────────────────────────────────────────────
H0_LCDM  = 67.36
OMM_LCDM = 0.3153
OMB_H2_L = 0.02237
ns_LCDM  = 0.9649
As_LCDM  = 2.1e-9


def T_WDM(k, alpha, mu=1.12):
    """Viel+2005 WDM transfer function."""
    return (1 + (alpha * k)**(2*mu))**(-5/mu)


def sigma_R(k_arr, pk_arr, R):
    """Top-hat rms sigma(R) in Mpc/h units."""
    def integrand(lnk):
        kk = np.exp(lnk)
        x  = kk * R
        if x < 1e-6:
            W = 1.0
        else:
            W = 3 * (np.sin(x) - x * np.cos(x)) / x**3
        Pk = np.interp(kk, k_arr, pk_arr)
        return kk**3 * Pk * W**2
    lnk_min = np.log(max(k_arr[0], 1e-4))
    lnk_max = np.log(min(k_arr[-1], 5.0))
    res, _ = quad(integrand, lnk_min, lnk_max, limit=2000,
                  epsabs=0, epsrel=1e-5)
    return np.sqrt(res / (2 * pi_**2))


def run_camb(H0, ombh2, omch2, w, wa, ns, As):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2,
                       mnu=0.06, omk=0, tau=0.0544)
    pars.set_dark_energy(w=w, wa=wa, dark_energy_model='ppf')
    pars.InitPower.set_params(ns=ns, As=As)
    pars.set_matter_power(redshifts=[0.0], kmax=10.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    kh, _, pk = results.get_matter_power_spectrum(minkh=5e-4, maxkh=5.0, npoints=600)
    return kh, pk[0]


print("\n  Running CAMB (LCDM + SSEE) ...")

omch2_L = OMM_LCDM * (H0_LCDM/100)**2 - OMB_H2_L
kL, pkL = run_camb(H0_LCDM, OMB_H2_L, omch2_L, -1.0, 0.0, ns_LCDM, As_LCDM)

omch2_S = OMM_MIRA * (H0_alg/100)**2 - OMB_H2
kS, pkS = run_camb(H0_alg, OMB_H2, omch2_S, w0_SSEE, wa_SSEE, ns_SSEE, As_SSEE)

# SSEE two-sector: apply WDM mixing (Paper 6 formula)
T2 = T_WDM(kS, alpha_WDM)
pk_2s = pkS * ((1 - f_phiDM) + f_phiDM * T2)**2

print("  CAMB done.")

# ── sigma(R) over a range of void radii ───────────────────────────────────────
R_vals = np.array([8, 20, 50, 80, 100, 130, 150, 200, 250, 300])

print("\n── sigma(R) — LCDM vs SSEE two-sector ────────────────────────")
print(f"  {'R (Mpc/h)':<12}  {'sigma_LCDM':<13}  {'sigma_SSEE':<13}  {'ratio'}")

sig_L   = []
sig_S   = []
for R in R_vals:
    sL = sigma_R(kL, pkL,  R)
    sS = sigma_R(kS, pk_2s, R)
    sig_L.append(sL)
    sig_S.append(sS)
    print(f"  {R:<12}  {sL:<13.5f}  {sS:<13.5f}  {sS/sL:.4f}")

sig_L = np.array(sig_L)
sig_S = np.array(sig_S)

# ── Gaussian void probability ─────────────────────────────────────────────────
# P(delta < delta_thresh) = 0.5 * erfc(delta_thresh / (sqrt(2) * sigma))
# For underdensity: delta_thresh < 0 -> P is > 0.5 for sigma > 0

print("\n── P(delta < -0.46) — required void depth for KAL resolution ──")
delta_thresh = -0.46   # ~46% underdensity
print(f"  Threshold: delta = {delta_thresh}")
print(f"\n  {'R (Mpc/h)':<12}  {'P_LCDM':<14}  {'P_SSEE':<14}  {'ratio P'}")

for i, R in enumerate(R_vals):
    x_L = -delta_thresh / (np.sqrt(2) * sig_L[i])
    x_S = -delta_thresh / (np.sqrt(2) * sig_S[i])
    pL  = 0.5 * erfc(x_L)
    pS  = 0.5 * erfc(x_S)
    ratio = pS / pL if pL > 0 else float('inf')
    print(f"  {R:<12}  {pL:<14.4e}  {pS:<14.4e}  {ratio:.3f}x")

# ── Implied Delta_H0 for different void depths at R=200 Mpc/h ─────────────────
print("\n── Delta_H0 from local void at R = 200 Mpc/h ─────────────────")
print(f"  Formula: Delta_H0 = H0 * |delta| * f(Omega_m) / 3")
print(f"  f(Omega_m,CMB=0.320) = (0.320)^0.55 = {0.320**0.55:.4f}")
f_growth = 0.320**0.55
print(f"  f(Omega_m,dyn=0.160) = (0.160)^0.55 = {0.160**0.55:.4f}")
f_growth_dyn = 0.160**0.55

print(f"\n  {'|delta|':<10}  {'DH0 (CMB Omega_m)':<22}  {'DH0 (dyn Omega_m)'}")
for delta_abs in [0.20, 0.30, 0.40, 0.46, 0.50, 0.55, 0.60, KAL*3/(H0_alg*f_growth)]:
    dH_cmb = H0_alg * delta_abs * f_growth / 3
    dH_dyn = H0_alg * delta_abs * f_growth_dyn / 3
    marker = "  ← KAL" if abs(dH_cmb - KAL) < 0.15 else ""
    print(f"  {delta_abs:<10.4f}  {dH_cmb:<22.4f}  {dH_dyn:.4f}{marker}")

# ── KAL derivation check ──────────────────────────────────────────────────────
print("\n── Algebraic check: KAL = H0 * delta_KAL * f / 3 ────────────")
delta_KAL_cmb = 3 * KAL / (H0_alg * f_growth)
delta_KAL_dyn = 3 * KAL / (H0_alg * f_growth_dyn)
print(f"  KAL = {KAL:.5f}")
print(f"  Required |delta| (f from Omega_m,CMB) = {delta_KAL_cmb:.4f}")
print(f"  Required |delta| (f from Omega_m,dyn) = {delta_KAL_dyn:.4f}")

# sigma at R=200 Mpc/h
i200 = list(R_vals).index(200)
print(f"\n  sigma(R=200) LCDM  = {sig_L[i200]:.5f}")
print(f"  sigma(R=200) SSEE  = {sig_S[i200]:.5f}")
print(f"  delta_KAL / sigma_LCDM = {delta_KAL_cmb/sig_L[i200]:.3f} sigma")
print(f"  delta_KAL / sigma_SSEE = {delta_KAL_cmb/sig_S[i200]:.3f} sigma")

# ── Comparison with SH0ES ─────────────────────────────────────────────────────
H0_shoes  = 73.04
H0_planck = 67.36
gap_obs   = H0_shoes - H0_alg
print("\n── Observed gap vs KAL ────────────────────────────────────────")
print(f"  H0(SH0ES)   = {H0_shoes}")
print(f"  H0(SSEE alg) = {H0_alg:.4f}")
print(f"  Gap observed = {gap_obs:.4f}  km/s/Mpc")
print(f"  KAL          = {KAL:.4f}  km/s/Mpc")
print(f"  KAL - gap    = {KAL - gap_obs:+.4f}  km/s/Mpc  ({(KAL/gap_obs - 1)*100:+.1f}%)")

# ── Physical summary ──────────────────────────────────────────────────────────
print("\n── Physical summary ───────────────────────────────────────────")
pS200 = 0.5 * erfc(-delta_KAL_cmb / (np.sqrt(2) * sig_S[i200]))
pL200 = 0.5 * erfc(-delta_KAL_cmb / (np.sqrt(2) * sig_L[i200]))
print(f"  P(inhabit void with |delta|>={delta_KAL_cmb:.3f}, R=200 Mpc/h):")
print(f"    LCDM:  {pL200:.4e}")
print(f"    SSEE:  {pS200:.4e}")
if pL200 > 0:
    print(f"    SSEE/LCDM ratio: {pS200/pL200:.3f}x")
print(f"\n  KAL = beta + pi = {beta:.4f} + {pi_:.4f} = {KAL:.4f}")
print(f"  If delta_void = {delta_KAL_cmb:.4f} in R=200 Mpc/h sphere:")
print(f"    Delta_H0 = {H0_alg * delta_KAL_cmb * f_growth / 3:.4f} km/s/Mpc")
print(f"    H0_local = H0_alg + Delta_H0 = {H0_alg + H0_alg * delta_KAL_cmb * f_growth / 3:.4f}")
print(f"    H0_SH0ES = {H0_shoes}")
print("=" * 65)
