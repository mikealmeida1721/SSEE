"""
SSEE non-linear S8 correction via CAMB Halofit
Compares linear vs non-linear sigma8, S8 for:
  - LCDM (Planck 2018)
  - SSEE MIRA-sector (Omega_m=0.3199, w0=-0.840, wa=-0.670)
  - SSEE two-sector + WDM (Paper 6 mixing formula: P_mix = P*((1-f)+f*T)^2)

Verification mandate: all intermediate values printed before use.
"""

import numpy as np
import camb
from scipy.integrate import quad

# ── algebraic constants ────────────────────────────────────────────────────────
phi  = (1 + 5**0.5) / 2
pi_  = np.pi

H0_SSEE  = 3*(phi + pi_)**2          # 67.96 km/s/Mpc
w0_SSEE  = -0.840
wa_SSEE  = -0.670
OMM_MIRA = 0.3199                    # MIRA-mapped Omega_m
OMB_H2   = 0.02237
ns_SSEE  = 1 - phi**(-7)            # 0.9656
As_SSEE  = 2.1e-9

# WDM: Paper 6 calibration (alpha_WDM=1.6561, f_phiDM=0.5)
alpha_WDM = 1.6561   # h/Mpc  (calibrated: sigma8_eff=0.737 vs CLASS)
f_phiDM   = 0.5      # Omega_phiDM / Omega_m,total  (≈ 0.160/0.320)
k_fs      = 0.493    # h/Mpc  free-streaming scale

# LCDM Planck 2018
H0_LCDM  = 67.36
OMM_LCDM = 0.3153
OMB_H2_L = 0.02237
ns_LCDM  = 0.9649
As_LCDM  = 2.1e-9


def T_WDM(k, alpha, mu=1.12):
    """Viel+2005 WDM transfer function."""
    return (1 + (alpha * k)**(2*mu))**(-5/mu)


def sigma8_from_pk(k_arr, pk_arr, R=8.0):
    """Top-hat sigma(R) from P(k) [h/Mpc, (Mpc/h)^3] via Peebles 1980."""
    def integrand(lnk):
        kk = np.exp(lnk)
        x  = kk * R
        W  = 3*(np.sin(x) - x*np.cos(x)) / x**3
        Pk = np.interp(kk, k_arr, pk_arr)
        return kk**3 * Pk * W**2
    lnk_min = np.log(max(k_arr[0],  1e-4))
    lnk_max = np.log(min(k_arr[-1], 5.0))
    res, err = quad(integrand, lnk_min, lnk_max, limit=2000,
                    epsabs=0, epsrel=1e-5)
    return np.sqrt(res / (2 * pi_**2))


def run_camb(H0, ombh2, omch2, w, wa, ns, As, nonlinear=False):
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ombh2, omch2=omch2,
                       mnu=0.06, omk=0, tau=0.0544)
    pars.set_dark_energy(w=w, wa=wa, dark_energy_model='ppf')
    pars.InitPower.set_params(ns=ns, As=As)
    pars.set_matter_power(redshifts=[0.0], kmax=10.0)
    pars.NonLinear = (camb.model.NonLinear_both if nonlinear
                      else camb.model.NonLinear_none)
    results = camb.get_results(pars)
    kh, _, pk = results.get_matter_power_spectrum(minkh=5e-4, maxkh=5.0, npoints=600)
    return kh, pk[0]


# ── 1. LCDM ────────────────────────────────────────────────────────────────────
omch2_L = OMM_LCDM * (H0_LCDM/100)**2 - OMB_H2_L
kL, pkL_lin = run_camb(H0_LCDM, OMB_H2_L, omch2_L, -1.0, 0.0,
                       ns_LCDM, As_LCDM, nonlinear=False)
kL, pkL_nl  = run_camb(H0_LCDM, OMB_H2_L, omch2_L, -1.0, 0.0,
                       ns_LCDM, As_LCDM, nonlinear=True)

s8_lcdm_lin = sigma8_from_pk(kL, pkL_lin)
s8_lcdm_nl  = sigma8_from_pk(kL, pkL_nl)

# ── 2. SSEE MIRA-sector ────────────────────────────────────────────────────────
omch2_S = OMM_MIRA * (H0_SSEE/100)**2 - OMB_H2
kS, pkS_lin = run_camb(H0_SSEE, OMB_H2, omch2_S, w0_SSEE, wa_SSEE,
                       ns_SSEE, As_SSEE, nonlinear=False)
kS, pkS_nl  = run_camb(H0_SSEE, OMB_H2, omch2_S, w0_SSEE, wa_SSEE,
                       ns_SSEE, As_SSEE, nonlinear=True)

s8_ssee_lin = sigma8_from_pk(kS, pkS_lin)
s8_ssee_nl  = sigma8_from_pk(kS, pkS_nl)

# ── 3. SSEE two-sector: Paper 6 mixing formula ────────────────────────────────
# P_mix(k) = P_lin(k) * ((1 - f_phiDM) + f_phiDM * T_WDM(k))^2
# Same formula used in calibrate_wdm_alpha.py to calibrate alpha=1.6561
T2_S = T_WDM(kS, alpha_WDM)
pk_2sec_lin = pkS_lin * ((1 - f_phiDM) + f_phiDM * T2_S)**2

# Non-linear two-sector: apply WDM mixing to non-linear MIRA P(k)
# (conservative: Halofit applied to background, then WDM suppresses small scales)
pk_2sec_nl  = pkS_nl  * ((1 - f_phiDM) + f_phiDM * T2_S)**2

s8_2sec_lin = sigma8_from_pk(kS, pk_2sec_lin)
s8_2sec_nl  = sigma8_from_pk(kS, pk_2sec_nl)

# S8 = sigma8 * sqrt(Omega_m / 0.3)
def S8(s8, omm):
    return s8 * np.sqrt(omm / 0.3)

# ── 4. Observational targets ───────────────────────────────────────────────────
kids_s8, kids_err = 0.766, 0.020   # Asgari+2021 KiDS-1000
des_s8,  des_err  = 0.776, 0.017   # Abbott+2022 DES Y3

def t(val, ref, err):
    return (val - ref) / err

# ── Print verification ─────────────────────────────────────────────────────────
print("=" * 65)
print("SSEE non-linear S8 — CAMB Halofit verification")
print("=" * 65)

print("\n── Input parameters (verify before use) ──────────────────────")
print(f"  SSEE:  H0={H0_SSEE:.4f}  Omega_m={OMM_MIRA}  w0={w0_SSEE}  wa={wa_SSEE}")
print(f"  LCDM:  H0={H0_LCDM}   Omega_m={OMM_LCDM}   w0=-1.0   wa=0.0")
print(f"  WDM:   alpha={alpha_WDM} h/Mpc   f_phiDM={f_phiDM}   k_fs={k_fs} h/Mpc")
print(f"  WDM mixing: P_mix = P * ((1-f) + f*T_WDM)^2  [Paper 6 formula]")

print("\n── T_WDM spot check ───────────────────────────────────────────")
for kk in [0.1, 0.3, 0.493, 1.0, 3.0]:
    t_val = T_WDM(kk, alpha_WDM)
    mix   = ((1-f_phiDM) + f_phiDM*t_val)**2
    print(f"  k={kk:.3f}: T={t_val:.4f}  mix_factor={mix:.4f}  "
          f"suppression={100*(1-mix):.1f}%")

print("\n── sigma8 results ─────────────────────────────────────────────")
print(f"  {'Model':<20}  {'σ₈ linear':<12}  {'σ₈ non-lin':<12}  {'Halofit boost'}")
print(f"  {'LCDM':<20}  {s8_lcdm_lin:<12.4f}  {s8_lcdm_nl:<12.4f}  "
      f"+{(s8_lcdm_nl/s8_lcdm_lin-1)*100:.1f}%")
print(f"  {'SSEE MIRA':<20}  {s8_ssee_lin:<12.4f}  {s8_ssee_nl:<12.4f}  "
      f"+{(s8_ssee_nl/s8_ssee_lin-1)*100:.1f}%")
print(f"  {'SSEE two-sec':<20}  {s8_2sec_lin:<12.4f}  {s8_2sec_nl:<12.4f}  "
      f"+{(s8_2sec_nl/s8_2sec_lin-1)*100:.1f}%")

print("\n── S8 = sigma8 × sqrt(Omega_m/0.3) ───────────────────────────")
rows = [
    ("LCDM lin",      S8(s8_lcdm_lin, OMM_LCDM)),
    ("LCDM non-lin",  S8(s8_lcdm_nl,  OMM_LCDM)),
    ("SSEE lin",      S8(s8_ssee_lin, OMM_MIRA)),
    ("SSEE non-lin",  S8(s8_ssee_nl,  OMM_MIRA)),
    ("SSEE 2-sec lin",S8(s8_2sec_lin, OMM_MIRA)),
    ("SSEE 2-sec nl", S8(s8_2sec_nl,  OMM_MIRA)),
]
print(f"  {'Model':<20}  {'S8':<8}  {'KiDS σ':>8}  {'DES σ':>8}")
for label, val in rows:
    tk = t(val, kids_s8, kids_err)
    td = t(val, des_s8,  des_err)
    print(f"  {label:<20}  {val:.4f}    {tk:+.2f}σ    {td:+.2f}σ")

print("\n── Key comparison: CLASS two-sector result vs CAMB ───────────")
print(f"  CLASS two-sector sigma8_eff (Paper 6) = 0.7370  (calibrated)")
print(f"  CAMB  two-sector sigma8 lin            = {s8_2sec_lin:.4f}")
print(f"  CAMB  two-sector sigma8 non-lin        = {s8_2sec_nl:.4f}")
print(f"  Consistency: {'✅ within 2%' if abs(s8_2sec_lin - 0.737) < 0.02 else '⚠️ check'}")

print("\n── Non-linear boost suppression by WDM ───────────────────────")
boost_mira   = s8_ssee_nl / s8_ssee_lin
boost_2sec   = s8_2sec_nl / s8_2sec_lin
print(f"  Halofit boost, MIRA single-sector:  +{(boost_mira-1)*100:.1f}%")
print(f"  Halofit boost, two-sector + WDM:    +{(boost_2sec-1)*100:.1f}%")
print(f"  WDM reduces non-linear boost by:     {(boost_mira - boost_2sec)*100:.2f} pp")

print("=" * 65)
