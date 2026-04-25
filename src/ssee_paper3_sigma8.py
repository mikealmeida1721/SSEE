"""
SSEE-V3.6 — Paper 3: σ8 and S8 from matter power spectrum
Computes:
  1. σ8 directly from CAMB with SSEE parameters (w0/wa, H0_CMB=67.075)
  2. Growth index correction: γ_SSEE=φ⁻¹=0.618 vs γ_ΛCDM=0.55
  3. σ8_IS, S8_IS after IS-viscosity growth modification
  4. f(z)σ8(z) for RSD comparison

Key tension: S8 (KiDS: 0.759, DES: 0.773) vs Planck (0.832)
SSEE prediction with γ=φ⁻¹ should lower S8 toward weak-lensing values.
"""

import numpy as np
from scipy.integrate import quad
import os

phi = (1 + 5**0.5) / 2
pi  = np.pi

# SSEE algebraic parameters
Omega_  = pi + phi
beta_   = (pi + phi) / 2
P_sc    = Omega_ + phi
Kv      = phi + pi + Omega_
Tr      = 3 * (phi + beta_)
Mv      = phi + pi + Kv
MIRA    = (phi + (phi + pi) / 2) / 2

w0_s    = -Tr / Mv          # -0.8399
wa_s    = -P_sc / Kv        # -0.6700
Ofld    = Tr / Mv           # Omega_DE

# H0 from CMB chi2 minimum (sesión 5)
H0_s    = 67.075
ob_s    = 0.02237
Omm_cmb = (1 - Ofld) * MIRA  # 0.31993
oc_s    = Omm_cmb * (H0_s / 100)**2 - ob_s

ns_s    = 1.0 - (1.0 / phi)**7   # 0.96556
As_s    = np.exp(3.044) * 1e-10
tau_s   = 0.054

# Growth indices
gamma_lcdm = 0.55
gamma_ssee = 1.0 / phi  # φ⁻¹ = 0.6180

print(f"\nSSEE — σ8 / S8 analysis")
print("=" * 60)
print(f"w0={w0_s:.5f}  wa={wa_s:.5f}  H0={H0_s:.3f}")
print(f"Ωm,CMB={Omm_cmb:.5f}  Ωc h²={oc_s:.5f}")
print(f"γ_SSEE=φ⁻¹={gamma_ssee:.5f}  γ_ΛCDM={gamma_lcdm:.4f}")


# ── 1. CAMB: σ8 from matter power spectrum ──────────────────────────────────
def compute_sigma8_camb(H0, oc, w0, wa, label):
    import camb
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=ob_s, omch2=oc,
                       mnu=0.06, omk=0, tau=tau_s, YHe=0.2454)
    pars.set_dark_energy(w=w0, wa=wa, dark_energy_model='ppf')
    pars.InitPower.set_params(ns=ns_s, As=float(As_s))
    pars.set_matter_power(redshifts=[0.0, 0.5, 1.0], kmax=10.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    sig8_z = results.get_sigma8()    # sorted high-z first by CAMB
    sig8 = sig8_z[-1]               # z=0 is last after CAMB re-sort
    # f σ8: growth rate × σ8(z)
    derived = results.get_derived_params()
    return sig8, sig8_z, results


print("\n--- σ8 from CAMB (linear, z=0) ---")
s8_ssee, s8z_ssee, res_ssee = compute_sigma8_camb(H0_s, oc_s, w0_s, wa_s, "SSEE")
s8_lcdm, s8z_lcdm, res_lcdm = compute_sigma8_camb(67.36, 0.1200, -1.0, 0.0, "ΛCDM")

Omm_lcdm = (0.1200 + 0.02237) / (67.36 / 100)**2  # ≈ 0.315

print(f"  ΛCDM (H0=67.36, Ωm=0.315): σ8 = {s8_lcdm:.5f}")
print(f"  SSEE (H0={H0_s}, Ωm={Omm_cmb:.4f}): σ8 = {s8_ssee:.5f}")

S8_lcdm_camb = s8_lcdm * (Omm_lcdm / 0.3)**0.5
S8_ssee_camb = s8_ssee * (Omm_cmb / 0.3)**0.5
print(f"\n  S8 = σ8×(Ωm/0.3)^0.5:")
print(f"    ΛCDM: {S8_lcdm_camb:.5f}")
print(f"    SSEE (CAMB, sin γ corrección): {S8_ssee_camb:.5f}")


# ── 2. Growth index correction: γ=φ⁻¹ vs γ=0.55 ────────────────────────────
def growth_integral(gamma, Omm):
    """Compute D(z=0)/D_uncorrected by integrating f=Omega_m(z)^gamma.
    Returns ratio D_gamma / D_0.55 at z=0, normalized to D=1 at z=0 for ΛCDM.
    We integrate from z=2 (matter dominated, both ~same) to z=0.
    """
    def Omm_z(z):
        # DE density: Omega_fld × (1+z)^(3*(1+w0+wa)) * exp(−3*wa*z/(1+z))
        # For simplicity use SSEE w0/wa values
        a = 1.0 / (1.0 + z)
        w_eff = w0_s + wa_s * (1 - a)
        rho_de = Ofld * (1 + z)**(3 * (1 + w0_s + wa_s)) * np.exp(-3 * wa_s * z / (1 + z))
        rho_m  = Omm * (1 + z)**3
        # ignore radiation at low z
        return rho_m / (rho_m + rho_de)

    def integrand(z):
        f = Omm_z(z)**gamma
        return f / (1 + z)

    # Integrate ln D from z_hi to 0: D(0) = exp(integral)
    z_hi = 3.0
    val, _ = quad(integrand, 0, z_hi)
    return val


ln_D_lcdm = growth_integral(gamma_lcdm, Omm_lcdm)
ln_D_ssee = growth_integral(gamma_ssee, Omm_cmb)

# σ8 ∝ D(z=0): ratio of growth factors
# D_ssee / D_lcdm evaluated from same starting point
ratio_D = np.exp(ln_D_ssee - ln_D_lcdm)

# Correct SSEE σ8: CAMB already uses proper EOM for w0/wa, so the γ correction
# is an additional IS modification on top of the w0/wa background
# σ8_IS = σ8_CAMB × (D_γ=0.618 / D_γ=0.55_effective)
# We compute ratio relative to SSEE's own γ_eff from CAMB

# CAMB's effective growth index for w0/wa dark energy ≈ 0.55 + 0.05(1+w0)
gamma_camb_eff = 0.55 + 0.05 * (1 + w0_s)  # ≈ 0.542
ln_D_camb_eff = growth_integral(gamma_camb_eff, Omm_cmb)
ratio_IS = np.exp(ln_D_ssee - ln_D_camb_eff)

sigma8_IS = s8_ssee * ratio_IS
S8_IS = sigma8_IS * (Omm_cmb / 0.3)**0.5

print(f"\n--- Growth index correction ---")
print(f"  γ_ΛCDM = {gamma_lcdm:.3f}  →  ln D = {ln_D_lcdm:.5f}")
print(f"  γ_CAMB_eff(SSEE w0/wa) = {gamma_camb_eff:.4f}  →  ln D = {ln_D_camb_eff:.5f}")
print(f"  γ_SSEE = φ⁻¹ = {gamma_ssee:.5f}  →  ln D = {ln_D_ssee:.5f}")
print(f"  IS correction ratio D(γ=φ⁻¹)/D(γ_eff) = {ratio_IS:.6f}")
print(f"  σ8_IS = {s8_ssee:.5f} × {ratio_IS:.5f} = {sigma8_IS:.5f}")
print(f"  S8_IS = σ8_IS × (Ωm/0.3)^0.5 = {S8_IS:.5f}")


# ── 3. f(z)σ8(z) for RSD ────────────────────────────────────────────────────
print(f"\n--- fσ8(z) for RSD surveys ---")
z_rsd = [0.295, 0.510, 0.706, 0.934, 1.317, 2.330]  # DESI DR2 redshifts
print(f"  {'z':>6}  {'Ωm(z)':>8}  {'f_ΛCDM':>8}  {'f_SSEE':>8}  {'fσ8_SSEE':>10}")
print(f"  {'-'*55}")

for z in z_rsd:
    a = 1.0 / (1.0 + z)
    rho_de = Ofld * (1 + z)**(3 * (1 + w0_s + wa_s)) * np.exp(-3 * wa_s * z / (1 + z))
    rho_m  = Omm_cmb * (1 + z)**3
    Omm_z_val = rho_m / (rho_m + rho_de)
    f_lcdm = Omm_z_val**gamma_lcdm
    f_ssee = Omm_z_val**gamma_ssee
    # rough σ8(z) ~ σ8_IS × D(z)/D(0), D(z)/D(0) ≈ (1+z)^(-(f-1)) crudely
    # Use growth integral properly
    ln_Dz, _ = quad(lambda zp: Omm_cmb**gamma_ssee / (1 + zp), 0, z)
    Dz_ratio  = np.exp(-ln_Dz)
    fsig8_ssee = f_ssee * sigma8_IS * Dz_ratio
    print(f"  {z:>6.3f}  {Omm_z_val:>8.5f}  {f_lcdm:>8.4f}  {f_ssee:>8.4f}  {fsig8_ssee:>10.5f}")


# ── 4. Summary ───────────────────────────────────────────────────────────────
print(f"\n{'=' * 60}")
print(f"  RESUMEN S8")
print(f"  {'Modelo':<30}  {'σ8':>7}  {'S8':>7}")
print(f"  {'-'*48}")
print(f"  {'ΛCDM (CAMB)':<30}  {s8_lcdm:>7.5f}  {S8_lcdm_camb:>7.5f}")
print(f"  {'SSEE (CAMB, γ_eff)':<30}  {s8_ssee:>7.5f}  {S8_ssee_camb:>7.5f}")
print(f"  {'SSEE+IS (γ=φ⁻¹=0.618)':<30}  {sigma8_IS:>7.5f}  {S8_IS:>7.5f}")
print(f"  {'-'*48}")
print(f"  {'Planck 2018 (obs)':<30}  {'0.8120':>7}  {'0.8320':>7}  ±0.013")
print(f"  {'KiDS-1000 (obs)':<30}  {'—':>7}  {'0.7590':>7}  ±0.024")
print(f"  {'DES Y3 (obs)':<30}  {'—':>7}  {'0.7730':>7}  ±0.025")
print(f"  {'DES+KiDS joint':<30}  {'—':>7}  {'0.7660':>7}  ±0.020")

tension_planck = (S8_IS - 0.832) / 0.013
tension_kids   = (S8_IS - 0.759) / 0.024
tension_des    = (S8_IS - 0.773) / 0.025
print(f"\n  Tensión SSEE+IS vs Planck S8: {tension_planck:+.2f}σ")
print(f"  Tensión SSEE+IS vs KiDS-1000: {tension_kids:+.2f}σ")
print(f"  Tensión SSEE+IS vs DES Y3:    {tension_des:+.2f}σ")
print(f"{'=' * 60}")

# Save
out = os.path.join(os.path.dirname(__file__), "..", "results", "sigma8_s8_results.txt")
with open(out, "w") as f:
    f.write("SSEE-V3.6 — σ8 / S8 analysis\n")
    f.write("=" * 60 + "\n")
    f.write(f"w0={w0_s:.5f}  wa={wa_s:.5f}  H0={H0_s:.3f}  Ωm={Omm_cmb:.5f}\n")
    f.write(f"γ_SSEE=φ⁻¹={gamma_ssee:.5f}\n\n")
    f.write(f"ΛCDM (CAMB):             σ8={s8_lcdm:.5f}  S8={S8_lcdm_camb:.5f}\n")
    f.write(f"SSEE (CAMB, no IS):      σ8={s8_ssee:.5f}  S8={S8_ssee_camb:.5f}\n")
    f.write(f"SSEE+IS (γ=φ⁻¹):         σ8={sigma8_IS:.5f}  S8={S8_IS:.5f}\n\n")
    f.write(f"Observado:\n")
    f.write(f"  Planck 2018: S8=0.832±0.013\n")
    f.write(f"  KiDS-1000:   S8=0.759±0.024\n")
    f.write(f"  DES Y3:      S8=0.773±0.025\n\n")
    f.write(f"Tensión SSEE+IS vs Planck: {tension_planck:+.2f}σ\n")
    f.write(f"Tensión SSEE+IS vs KiDS:   {tension_kids:+.2f}σ\n")
    f.write(f"Tensión SSEE+IS vs DES:    {tension_des:+.2f}σ\n")
print(f"\n  Guardado en {out}")
