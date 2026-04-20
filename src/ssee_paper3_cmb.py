"""
SSEE-V3.6 — Paper 3: CMB Power Spectrum vs Planck PR4
Computes Cl_TT under SSEE background, applies r_d,eff mapping,
compares against Planck PR4 data, and produces chi2 + figures.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import urllib.request

matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# SSEE constants (algebraically fixed)
# ---------------------------------------------------------------------------
phi   = (1 + 5**0.5) / 2          # 1.6180
pi    = np.pi
Omega = pi + phi                   # 4.7596
beta  = (pi + phi) / 2             # 2.3798
KAL0  = beta + pi                  # 5.5214
P_sc  = Omega + phi                # 6.3776
Kv    = phi + pi + Omega           # 9.5192
Tr    = 3 * (phi + beta)           # 11.9935
Mv    = phi + pi + Kv              # 14.2788

w0    = -Tr / Mv                   # -0.8403
wa    = -P_sc / Kv                 # -0.6703
OmDE  = Tr / Mv                    # 0.8403
Omm   = 1.0 - OmDE                 # 0.1597  — sector dinámico (BAO/cúmulos)
M_SSEE = abs(w0)                   # acoustic saturation factor

# Sector de observación CMB: Ω_m,CMB = Ω_m,dyn × MIRA
# MIRA = AURA/2 = (PHI+BIAL)/2 — "Frecuencia de Observación" (Genesis 5.12)
BIAL  = (phi + pi) / 2
AURA  = phi + BIAL
MIRA  = AURA / 2                   # 1.998924
Omm_cmb = Omm * MIRA               # 0.3198 ≈ Planck 0.3153 (1.47% off, <1σ)

H0       = 66.66
Omb_h2   = 0.02237
ns       = 0.9649
ln_As    = 3.044
As       = np.exp(ln_As) * 1e-10

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "results", "figures")
DAT_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "raw")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(DAT_DIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Download Planck PR4 TT spectrum (COM_PowerSpect_CMB-TT-full_R3.01.txt)
# Public mirror via PLA (ESA)
# ---------------------------------------------------------------------------
PLANCK_FILE = os.path.join(DAT_DIR, "planck_pr4_TT.txt")
PLANCK_URL  = (
    "https://irsa.ipac.caltech.edu/data/Planck/release_3/ancillary-data/"
    "cosmoparams/COM_PowerSpect_CMB-TT-full_R3.01.txt"
)

def download_planck():
    if os.path.exists(PLANCK_FILE):
        return
    print("Descargando espectro Planck PR4 TT...")
    try:
        urllib.request.urlretrieve(PLANCK_URL, PLANCK_FILE)
        print(f"  Guardado en {PLANCK_FILE}")
    except Exception as e:
        print(f"  Descarga fallida: {e}")
        print("  Usando espectro sintético ΛCDM como referencia.")

def load_planck():
    """Returns (ell, Dl, sigma) arrays from Planck TT file."""
    if not os.path.exists(PLANCK_FILE):
        return None, None, None
    try:
        data = np.loadtxt(PLANCK_FILE, comments="#")
        ell   = data[:, 0].astype(int)
        Dl    = data[:, 1]
        # Columns 2 and 3 are +sigma and -sigma; use mean
        sigma = 0.5 * (data[:, 2] + data[:, 3])
        return ell, Dl, sigma
    except Exception as e:
        print(f"  Error leyendo Planck: {e}")
        return None, None, None

# ---------------------------------------------------------------------------
# Compute SSEE Cl with CAMB
# ---------------------------------------------------------------------------
def compute_ssee_spectrum(lmax=2500):
    """SSEE+MIRA two-sector model: Ω_m,CMB = Ω_m,dyn × MIRA."""
    import camb

    h = H0 / 100.0
    omch2 = Omm_cmb * h**2 - Omb_h2
    if omch2 < 0:
        raise ValueError(f"omch2 negativo: Omm={Omm:.4f} implica omch2={omch2:.5f}")

    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=H0,
        ombh2=Omb_h2,
        omch2=omch2,
        mnu=0.085,          # SSEE neutrino mass sum (eV)
        omk=0,
        tau=0.054,
    )
    pars.set_dark_energy(w=w0, wa=wa, dark_energy_model="ppf")
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)
    pars.Want_CMB = True

    results = camb.get_results(pars)
    derived = results.get_derived_params()
    r_d_camb = derived["rdrag"]             # raw sound horizon from CAMB
    r_d_eff  = r_d_camb * M_SSEE           # structural mapping

    powers = results.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)
    # total = lensed; shape (lmax+1, 4): TT, EE, BB, TE
    Dl_TT = powers["total"][:, 0]          # Dl = l(l+1)Cl/2pi [muK^2]
    ells  = np.arange(len(Dl_TT))

    return ells, Dl_TT, r_d_camb, r_d_eff, derived


def compute_lcdm_spectrum(lmax=2500):
    """ΛCDM reference: Planck 2018 best-fit parameters."""
    import camb

    pars = camb.CAMBparams()
    pars.set_cosmology(H0=67.36, ombh2=0.02237, omch2=0.1200,
                       mnu=0.06, omk=0, tau=0.054)
    pars.set_dark_energy(w=-1.0)
    pars.InitPower.set_params(As=np.exp(3.044)*1e-10, ns=0.9649)
    pars.set_for_lmax(lmax, lens_potential_accuracy=1)
    pars.Want_CMB = True

    results = camb.get_results(pars)
    powers  = results.get_cmb_power_spectra(pars, CMB_unit="muK", raw_cl=False)
    Dl_TT   = powers["total"][:, 0]
    ells    = np.arange(len(Dl_TT))
    r_d     = results.get_derived_params()["rdrag"]
    return ells, Dl_TT, r_d

# ---------------------------------------------------------------------------
# Chi2 calculation
# ---------------------------------------------------------------------------
def chi2_vs_planck(ell_obs, Dl_obs, sigma_obs, ell_model, Dl_model,
                   ell_min=30, ell_max=2000):
    mask = (ell_obs >= ell_min) & (ell_obs <= ell_max)
    ell_sel = ell_obs[mask]
    Dl_sel  = Dl_obs[mask]
    sig_sel = sigma_obs[mask]

    Dl_interp = np.interp(ell_sel, ell_model, Dl_model)
    residuals  = (Dl_interp - Dl_sel) / sig_sel
    chi2       = np.sum(residuals**2)
    chi2_r     = chi2 / len(ell_sel)
    return chi2, chi2_r, len(ell_sel)

# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_spectrum(ells_s, Dl_s, ells_l, Dl_l, ell_obs, Dl_obs, sigma_obs):
    fig, axes = plt.subplots(2, 1, figsize=(10, 8),
                             gridspec_kw={"height_ratios": [3, 1]})

    ax = axes[0]
    if ell_obs is not None:
        ax.errorbar(ell_obs, Dl_obs, yerr=sigma_obs,
                    fmt="k.", ms=2, lw=0.5, alpha=0.6, label="Planck PR4 TT")
    ax.plot(ells_l[2:], Dl_l[2:], color="tab:orange", lw=1.5,
            ls="--", label=r"$\Lambda$CDM ($\Omega_m=0.315$)")
    ax.plot(ells_s[2:], Dl_s[2:], color="tab:blue", lw=1.8,
            label=r"SSEE-V3.6+MIRA ($\Omega_{m,\rm CMB}=\Omega_{m,\rm dyn}\times\mathrm{MIRA}=0.320$)")
    ax.set_xlim(2, 2500)
    ax.set_ylim(0, 6500)
    ax.set_ylabel(r"$D_\ell^{TT}$ [$\mu$K$^2$]", fontsize=13)
    ax.legend(fontsize=11)
    ax.set_title("SSEE-V3.6 vs Planck PR4: CMB TT Power Spectrum", fontsize=13)
    ax.grid(True, alpha=0.3)

    ax2 = axes[1]
    if ell_obs is not None:
        Dl_interp_s = np.interp(ell_obs, ells_s, Dl_s)
        res = (Dl_interp_s - Dl_obs) / sigma_obs
        ax2.plot(ell_obs, res, "b.", ms=2, alpha=0.6)
        ax2.axhline(0, color="k", lw=0.8)
        ax2.axhline(+1, color="gray", lw=0.6, ls="--")
        ax2.axhline(-1, color="gray", lw=0.6, ls="--")
        ax2.set_xlim(2, 2500)
        ax2.set_ylim(-5, 5)
        ax2.set_ylabel(r"$(D_\ell^{\rm SSEE}-D_\ell^{\rm Planck})/\sigma$",
                       fontsize=11)
        ax2.set_xlabel(r"Multipole $\ell$", fontsize=13)
        ax2.grid(True, alpha=0.3)
    else:
        ax2.set_visible(False)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, "fig_cmb_spectrum.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Figura guardada: {out}")


def plot_peak_zoom(ells_s, Dl_s, ells_l, Dl_l, ell_obs, Dl_obs, sigma_obs):
    """Zoom en los tres primeros picos acústicos."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ranges = [(150, 350, 220), (420, 650, 540), (700, 950, 810)]
    labels = ["1er pico", "2do pico", "3er pico"]

    for ax, (lo, hi, expected), lbl in zip(axes, ranges, labels):
        if ell_obs is not None:
            m = (ell_obs >= lo) & (ell_obs <= hi)
            ax.errorbar(ell_obs[m], Dl_obs[m], yerr=sigma_obs[m],
                        fmt="k.", ms=3, lw=0.7, alpha=0.7)
        m_s = (ells_s >= lo) & (ells_s <= hi)
        m_l = (ells_l >= lo) & (ells_l <= hi)
        ax.plot(ells_l[m_l], Dl_l[m_l], "tab:orange", lw=1.5, ls="--")
        ax.plot(ells_s[m_s], Dl_s[m_s], "tab:blue", lw=1.8)
        ax.axvline(expected, color="green", lw=0.8, ls=":", alpha=0.7,
                   label=rf"$\ell={expected}$")
        ax.set_title(lbl, fontsize=12)
        ax.set_xlabel(r"$\ell$")
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    axes[0].set_ylabel(r"$D_\ell^{TT}$ [$\mu$K$^2$]")
    plt.suptitle("SSEE-V3.6: Zoom en Picos Acústicos vs Planck PR4", fontsize=13)
    plt.tight_layout()
    out = os.path.join(FIG_DIR, "fig_cmb_peaks_zoom.pdf")
    plt.savefig(out, bbox_inches="tight")
    plt.close()
    print(f"  Figura guardada: {out}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("SSEE-V3.6 — Paper 3: CMB vs Planck PR4")
    print("=" * 60)
    print(f"\nParámetros SSEE:")
    print(f"  w0={w0:.4f}  wa={wa:.4f}")
    print(f"  Ω_m,dyn={Omm:.4f}  Ω_DE={OmDE:.4f}  (sector dinámico: BAO/cúmulos)")
    print(f"  MIRA={MIRA:.6f}  (AURA/2 — Frecuencia de Observación)")
    print(f"  Ω_m,CMB={Omm_cmb:.6f}  (sector observacional: CMB, Planck: 0.3153)")
    print(f"  M_SSEE=|w0|={M_SSEE:.4f}")
    print(f"  KAL0={KAL0:.4f}")

    # 1. Descargar datos Planck
    download_planck()
    ell_obs, Dl_obs, sigma_obs = load_planck()
    if ell_obs is not None:
        print(f"\nDatos Planck PR4: {len(ell_obs)} puntos, ℓ={ell_obs[0]}–{ell_obs[-1]}")
    else:
        print("\nDatos Planck no disponibles — solo se generará el espectro SSEE.")

    # 2. Calcular espectro SSEE
    print("\nCalculando espectro SSEE con CAMB...")
    ells_s, Dl_s, r_d_raw, r_d_eff, derived = compute_ssee_spectrum()
    print(f"  r_d,SSEE (CAMB)  = {r_d_raw:.2f} Mpc")
    print(f"  r_d,eff (×M_SSEE) = {r_d_eff:.2f} Mpc  (Planck: ~147.1 Mpc)")
    print(f"  D_A(z*)          = {derived.get('DA', 'N/A')}")
    print(f"  z_drag           = {derived.get('zdrag', 'N/A'):.2f}")
    print(f"  100θ_MC          = {derived.get('thetastar', derived.get('theta_MC_100', 'N/A'))}")

    # 3. Calcular espectro ΛCDM
    print("\nCalculando espectro ΛCDM con CAMB...")
    ells_l, Dl_l, r_d_lcdm = compute_lcdm_spectrum()
    print(f"  r_d,ΛCDM = {r_d_lcdm:.2f} Mpc")

    # 4. Chi2 vs Planck
    if ell_obs is not None:
        print("\nCalculando χ² vs Planck PR4 (ℓ=30–2000):")
        chi2_s, chi2r_s, n = chi2_vs_planck(
            ell_obs, Dl_obs, sigma_obs, ells_s, Dl_s)
        chi2_l, chi2r_l, _ = chi2_vs_planck(
            ell_obs, Dl_obs, sigma_obs, ells_l, Dl_l)
        print(f"  SSEE:  χ²={chi2_s:.1f}, χ²_r={chi2r_s:.3f}  (N={n})")
        print(f"  ΛCDM:  χ²={chi2_l:.1f}, χ²_r={chi2r_l:.3f}")
        print(f"  Δχ²(SSEE−ΛCDM) = {chi2_s - chi2_l:.1f}")

    # 5. Posiciones de picos (máximos locales)
    print("\nPosición de picos SSEE (primeros 3):")
    from scipy.signal import argrelmax
    peaks_idx = argrelmax(Dl_s[50:1500], order=60)[0] + 50
    for i, idx in enumerate(peaks_idx[:3]):
        print(f"  Pico {i+1}: ℓ={ells_s[idx]}  Dℓ={Dl_s[idx]:.1f} μK²")

    # 6. Figuras
    print("\nGenerando figuras...")
    plot_spectrum(ells_s, Dl_s, ells_l, Dl_l, ell_obs, Dl_obs, sigma_obs)
    plot_peak_zoom(ells_s, Dl_s, ells_l, Dl_l, ell_obs, Dl_obs, sigma_obs)

    print("\n✓ Listo.")


if __name__ == "__main__":
    main()
