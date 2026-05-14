"""
SSEE-V3.6 — Fase 2b/3: Calibración numérica del alpha WDM
Encuentra el alpha que reproduce sigma_8_eff = 0.737 (Paper 6) con CLASS P(k).

Usa el filtro top-hat estándar (R=8 Mpc/h) para sigma_8 correcta.
Fase 3 encontró que la proxy anterior (sin filtro) subestimaba la supresión
necesaria. Este script usa la fórmula sigma_8 estándar de cosmología.
"""
import numpy as np
from scipy.optimize import brentq

ssee_pk = np.loadtxt("output/ssee_v36__pk.dat",        comments="#")
lcdm_pk = np.loadtxt("output/lcdm_planck2018__pk.dat", comments="#")
k_s, P_s = ssee_pk[:,0], ssee_pk[:,1]
k_l, P_l = lcdm_pk[:,0], lcdm_pk[:,1]
P_l_i = np.interp(k_s, k_l, P_l)

# Planck 2018: sigma_8 = 0.811 (best fit)
# Paper 6 SSEE 2-sector: sigma_8_eff = 0.737
sigma8_planck  = 0.811
sigma8_target  = 0.737   # Paper 6 two-sector

mu       = 1.12
k_fs     = 0.493
f_phiDM  = 0.5   # Omega_phiDM / Omega_m = 0.160/0.320

def T_wdm(k, alpha):
    return (1 + (alpha * k)**(2*mu))**(-5/mu)

def W_tophat(k, R=8.0):
    x = k * R
    return 3.0 * (np.sin(x) - x * np.cos(x)) / x**3

def sigma8_real(k, P):
    """Sigma_8 estándar con filtro top-hat R=8 Mpc/h (Peebles 1980)."""
    W = W_tophat(k)
    return np.sqrt(np.trapezoid(k**2 * P * W**2 / (2 * np.pi**2), k))

sig_lcdm = sigma8_real(k_s, P_l_i)
sig_ssee  = sigma8_real(k_s, P_s)
# Ratio target relativo a la LCDM que CLASS computa (no sigma8_planck=0.811)
ratio_target = sigma8_target / sigma8_planck   # 0.9087 vs Planck
ratio_target_class = ratio_target * sigma8_planck / sig_lcdm  # ajuste al ref CLASS

print(f"ΛCDM  sigma_8 (CLASS top-hat): {sig_lcdm:.4f}")
print(f"SSEE+MIRA sigma_8:             {sig_ssee:.4f}  (ratio={sig_ssee/sig_lcdm:.4f})")
print(f"Target sigma_8 (Paper 6):      {sigma8_target}  ratio_vs_Planck={ratio_target:.4f}")
print(f"Target ratio vs CLASS LCDM:    {ratio_target_class:.4f}")
print()

def ratio_diff(alpha):
    T2 = T_wdm(k_s, alpha)
    P_mix = P_s * ((1-f_phiDM) + f_phiDM * T2)**2
    return sigma8_real(k_s, P_mix) / sig_lcdm - ratio_target_class

# Scan para ver el rango
print("── Alpha scan (sigma_8 real con top-hat) ────────────────────────")
alphas = np.logspace(-1, 1, 25)
r_prev = ratio_diff(alphas[0])
bracket = None
for a in alphas:
    T2 = T_wdm(k_s, a)
    P_mix = P_s * ((1-f_phiDM) + f_phiDM * T2)**2
    r = sigma8_real(k_s, P_mix) / sig_lcdm
    print(f"  alpha={a:.4f}  →  sigma_8={sigma8_real(k_s,P_mix):.4f}  ratio={r:.4f}")
    if r < ratio_target_class and bracket is None:
        bracket = a

print()

# Encontrar raíz
if bracket is None:
    print("⚠️  Target ratio no alcanzado — σ₈ SSEE MIRA > target incluso sin WDM")
    print(f"   σ₈ base SSEE MIRA = {sig_ssee:.4f} vs target {sigma8_target}")
    print(f"   Diferencia física: SSEE con Ω_m=0.320 y EoS SSEE da más poder que ΛCDM")
    print(f"   El WDM con f_phiDM=0.5 puede suprimir hasta ~30% en k>k_fs")
    # Encontrar alpha que maximiza la supresión vs. target
    try:
        alpha_cal = brentq(ratio_diff, 0.1, 20.0, xtol=1e-5)
        T2_cal = T_wdm(k_s, alpha_cal)
        P_cal  = P_s * ((1-f_phiDM) + f_phiDM * T2_cal)**2
        sig_cal = sigma8_real(k_s, P_cal)
        print(f"\n✅  Calibrated alpha = {alpha_cal:.6f} h/Mpc")
        print(f"   sigma_8 = {sig_cal:.4f}  ratio={sig_cal/sig_lcdm:.4f}  (target={ratio_target_class:.4f})")
        print(f"   T_WDM(k_fs)  = {T_wdm(k_fs, alpha_cal):.4f}")
        print(f"   Supresión en k_fs: {(1-T_wdm(k_fs, alpha_cal)**2)*100:.1f}%")
        np.save("output/alpha_calibrated.npy", alpha_cal)
    except Exception as e:
        print(f"\nNo se encontró raíz: {e}")
        # Guardar alpha que da la supresión más cercana al target
        alphas_fine = np.logspace(-1, 2, 200)
        ratios = []
        for a in alphas_fine:
            T2 = T_wdm(k_s, a)
            P_mix = P_s * ((1-f_phiDM) + f_phiDM * T2)**2
            ratios.append(sigma8_real(k_s, P_mix) / sig_lcdm)
        ratios = np.array(ratios)
        idx_min = np.argmin(np.abs(ratios - ratio_target_class))
        alpha_best = alphas_fine[idx_min]
        T2_best = T_wdm(k_s, alpha_best)
        P_best  = P_s * ((1-f_phiDM) + f_phiDM * T2_best)**2
        sig_best = sigma8_real(k_s, P_best)
        print(f"\n⚠️  Mejor alpha disponible = {alpha_best:.4f} h/Mpc")
        print(f"   sigma_8 = {sig_best:.4f}  ratio={sig_best/sig_lcdm:.4f}  (target={ratio_target_class:.4f})")
        print(f"   Δσ₈ residual: {sig_best - sigma8_target:.4f}")
        np.save("output/alpha_calibrated.npy", alpha_best)
        np.save("output/alpha_bestfit.npy", np.array([alpha_best, sig_best, ratio_target_class]))
else:
    try:
        alpha_cal = brentq(ratio_diff, 0.1, bracket * 2, xtol=1e-5)
        T2_cal = T_wdm(k_s, alpha_cal)
        P_cal  = P_s * ((1-f_phiDM) + f_phiDM * T2_cal)**2
        sig_cal = sigma8_real(k_s, P_cal)
        print(f"\n✅  Calibrated alpha = {alpha_cal:.6f} h/Mpc")
        print(f"   sigma_8 = {sig_cal:.4f}  ratio={sig_cal/sig_lcdm:.4f}  (target={ratio_target_class:.4f})")
        print(f"   T_WDM(k_fs)  = {T_wdm(k_fs, alpha_cal):.4f}")
        print(f"   Supresión en k_fs: {(1-T_wdm(k_fs, alpha_cal)**2)*100:.1f}%")
        np.save("output/alpha_calibrated.npy", alpha_cal)
    except Exception as e:
        print(f"Root not found: {e}")
