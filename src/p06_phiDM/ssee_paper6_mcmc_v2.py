#!/usr/bin/env python3
"""
Paper 6 — MCMC two-sector φ-DM v2 (RIGUROSO, emulador CLASS validado).
======================================================================
Reemplaza ssee_paper6_mcmc.py (RETIRADO: σ₈ single-k roto → S₈=0.820 espurio).

Física correcta:
  - σ₈(z), S₈ vienen de CLASS REAL (grilla p6_grid.npz) + integral top-hat.
  - Amplitud: floto A_s (primordial, prior Planck), NO σ₈.  σ₈ ∝ √A_s exacto.
    Historia física: el CMB fija A_s; el free-streaming φ-DM DERIVA el σ₈ tardío.
  - Ω_φDM, m_φ: priors PLANOS (test honesto — el dato decide, no el prior).

Parámetros:  θ = (Ω_φDM, m_φ [eV], A_s9 ≡ 10⁹A_s)
Datos:       fσ₈×6 (RSD) + S₈ KiDS-1000 + S₈ DES Y3
Salida:      posterior, tensión vs punto algebraico forward (0.149, 40.70),
             σ₈_eff/S₈_eff, χ² forward vs posterior, ΔBIC honesto vs ΛCDM.
Validación:  re-corre CLASS directo en ~10 puntos del posterior, reporta error emulador.
"""
import numpy as np, emcee, os, sys, subprocess, glob, time
from scipy.interpolate import RegularGridInterpolator
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner

GRIDNPZ = "/mnt/datos/SSEE_data/mcmc/p6_grid.npz"
OUTDIR  = "/mnt/datos/SSEE_data/mcmc/"
os.makedirs(OUTDIR, exist_ok=True)

# ── Cargar grilla CLASS ───────────────────────────────────────────────────────
G = np.load(GRIDNPZ)
g_phidm, g_mphi, z_dense = G["grid_phidm"], G["grid_mphi"], G["z_dense"]
s8z_grid = G["s8z_grid"]              # (nP, nM, nZ) σ₈(z) @ A_s=2.1e-9
As_ref   = float(G["As_ref"])*1e9    # → en unidades de 10⁹
Om_m     = float(G["Om_m"])          # 0.30889 fijo
M_PHI_CAN, OM_PHIDM_CAN = float(G["m_phi_can"]), float(G["om_phidm_can"])

# Interpolador σ₈(Ω_φDM, m_φ, z)  (lineal, sobre grilla CLASS real)
INTERP = RegularGridInterpolator((g_phidm, g_mphi, z_dense), s8z_grid,
                                 bounds_error=False, fill_value=None)

# ── Datos ─────────────────────────────────────────────────────────────────────
fsig8_data = np.array([
    [0.067, 0.423, 0.055], [0.150, 0.490, 0.145], [0.380, 0.497, 0.045],
    [0.510, 0.458, 0.038], [0.610, 0.436, 0.034], [1.480, 0.462, 0.045]])
S8_KIDS, S8_KIDS_E = 0.758, 0.0225   # KiDS-1000 Asgari+2021 (S8 directo)
S8_DES,  S8_DES_E  = 0.776, 0.017    # DES Y3 Amon/Secco+2022
GAMMA = 0.55                          # índice de crecimiento (RSD)
S8_FAC = np.sqrt(Om_m/0.3)           # S8 = σ8(0)·1.0147

# Prior Planck sobre A_s9 (primordial, externo): 10⁹A_s = 2.101 ± 0.030
AS9_MU, AS9_SIG = 2.101, 0.030

# ── Observables del modelo (vía emulador) ─────────────────────────────────────
def sigma8_curve(Om_phiDM, m_phi, As9):
    """σ₈(z_dense) escalado por amplitud:  σ₈ ∝ √(A_s/A_s_ref)."""
    pts = np.column_stack([np.full_like(z_dense, Om_phiDM),
                           np.full_like(z_dense, m_phi), z_dense])
    s8 = INTERP(pts)
    return s8 * np.sqrt(As9/As_ref)

def observables(Om_phiDM, m_phi, As9):
    s8z = sigma8_curve(Om_phiDM, m_phi, As9)
    s8_0 = s8z[0]                      # z_dense[0] = 0
    S8   = s8_0 * S8_FAC
    # f(z) = -(1+z) dlnσ8/dz  (crecimiento escala-dependiente capturado por CLASS)
    lns8 = np.log(s8z)
    dlnz = np.gradient(lns8, z_dense)
    fz   = -(1.0 + z_dense) * dlnz
    fsig8_z = fz * s8z
    return s8_0, S8, (z_dense, fsig8_z)

def chi2_model(Om_phiDM, m_phi, As9):
    s8_0, S8, (zc, fs) = observables(Om_phiDM, m_phi, As9)
    chi2 = 0.0
    for z, val, err in fsig8_data:
        pred = np.interp(z, zc, fs)
        chi2 += ((pred - val)/err)**2
    chi2 += ((S8 - S8_KIDS)/S8_KIDS_E)**2
    chi2 += ((S8 - S8_DES )/S8_DES_E )**2
    return chi2

# ── ΛCDM de referencia (CMB-fijo) sobre los MISMOS datos ──────────────────────
def chi2_lcdm():
    Om_l, s8_l = 0.3153, 0.811
    # crecimiento lineal ΛCDM: D(a) ∝ H(a)∫₀ᵃ da'/(a'H')³  (Heath, integral DESDE a≈0)
    from scipy.integrate import cumulative_trapezoid
    a = np.linspace(1e-4, 1.0, 5000)
    Ha = np.sqrt(Om_l/a**3 + (1-Om_l))
    integ = 1.0/(a*Ha)**3
    Dun = Ha * cumulative_trapezoid(integ, a, initial=0)
    D_of_a = Dun / Dun[-1]            # D(a=1)=1
    chi2 = 0.0
    for z, val, err in fsig8_data:
        az = 1.0/(1+z)
        D  = np.interp(az, a, D_of_a)
        Omz = Om_l*(1+z)**3 / (Om_l*(1+z)**3 + (1-Om_l))
        fs  = Omz**GAMMA * s8_l * D
        chi2 += ((fs - val)/err)**2
    S8_l = s8_l*np.sqrt(Om_l/0.3)     # 0.831
    chi2 += ((S8_l - S8_KIDS)/S8_KIDS_E)**2
    chi2 += ((S8_l - S8_DES )/S8_DES_E )**2
    return chi2, S8_l

# ── Posterior ─────────────────────────────────────────────────────────────────
B_PHIDM = (g_phidm[0], g_phidm[-1])   # plano, rango de grilla
B_MPHI  = (g_mphi[0],  g_mphi[-1])
def log_prob(theta):
    Om_phiDM, m_phi, As9 = theta
    if not (B_PHIDM[0] < Om_phiDM < B_PHIDM[1]): return -np.inf
    if not (B_MPHI[0]  < m_phi    < B_MPHI[1]):  return -np.inf
    if not (1.9 < As9 < 2.3):                    return -np.inf
    lp = -0.5*((As9-AS9_MU)/AS9_SIG)**2          # prior Planck A_s; resto plano
    return lp - 0.5*chi2_model(Om_phiDM, m_phi, As9)

# ── Validación: CLASS directo vs emulador en puntos del posterior ─────────────
def validate_emulator(samples, n=10):
    sys.path.insert(0, os.path.dirname(__file__))
    from ssee_paper6_mcmc_grid import run_class_multiz
    idx = np.random.default_rng(0).choice(len(samples), n, replace=False)
    errs = []
    print("\n── Validación emulador vs CLASS directo ──────────────────────────")
    for i in idx:
        opd, mp, As9 = samples[i]
        s8_emu = sigma8_curve(opd, mp, As9)[0]
        res = run_class_multiz(f"val_{i}", opd, mp, np.array([0.0]))
        if res is None: continue
        s8_cls = res[0][0.0]*np.sqrt(As9/As_ref)
        err = abs(s8_emu - s8_cls)/s8_cls*100
        errs.append(err)
        print(f"  Ω_φDM={opd:.3f} m_φ={mp:5.1f} A_s9={As9:.3f}  "
              f"σ8_emu={s8_emu:.4f} σ8_CLASS={s8_cls:.4f}  err={err:.3f}%")
    if errs:
        print(f"  → error emulador: medio {np.mean(errs):.3f}%  máx {np.max(errs):.3f}%")
    return errs

# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    ndim, nwalkers, nsteps, burnin = 3, 60, 12000, 3000
    print("="*70)
    print("  SSEE Paper 6 — MCMC two-sector φ-DM v2 (emulador CLASS)")
    print(f"  θ=(Ω_φDM, m_φ, A_s9)  priors PLANOS en φ-DM; Planck en A_s")
    print(f"  Ω_φDM∈{B_PHIDM}  m_φ∈({B_MPHI[0]:.1f},{B_MPHI[1]:.1f})  walkers={nwalkers} steps={nsteps}")
    print("="*70)

    rng = np.random.default_rng(42)
    p0 = np.column_stack([
        rng.uniform(0.08, 0.22, nwalkers),
        rng.uniform(20, 70, nwalkers),
        rng.normal(AS9_MU, AS9_SIG, nwalkers)])
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob)
    t0 = time.time()
    sampler.run_mcmc(p0, nsteps, progress=True)
    print(f"  muestreo: {time.time()-t0:.0f}s")
    np.save(OUTDIR+"p6_mcmc_v2_chain.npy", sampler.get_chain())

    try:
        tau = sampler.get_autocorr_time(tol=0)
        print(f"  τ={tau}  N_eff={((nsteps-burnin)*nwalkers/tau).astype(int)}")
    except Exception as e: print(f"  (autocorr {e})")
    print(f"  accept={np.mean(sampler.acceptance_fraction):.3f}")

    flat = sampler.get_chain(discard=burnin, thin=10, flat=True)
    def q(a): p=np.percentile(a,[16,50,84]); return p[1],p[2]-p[1],p[1]-p[0]
    Om_q, mp_q, As_q = q(flat[:,0]), q(flat[:,1]), q(flat[:,2])

    # derivados
    s8_chain = np.array([observables(*s)[0] for s in flat])
    S8_chain = s8_chain*S8_FAC
    s8_q, S8_q = q(s8_chain), q(S8_chain)

    # tensión vs punto algebraico forward
    ten_Om = abs(Om_q[0]-OM_PHIDM_CAN)/np.mean(Om_q[1:])
    ten_mp = abs(mp_q[0]-M_PHI_CAN)/np.mean(mp_q[1:])

    # χ²: forward algebraico vs posterior median
    chi2_fwd = chi2_model(OM_PHIDM_CAN, M_PHI_CAN, AS9_MU)
    chi2_med = chi2_model(Om_q[0], mp_q[0], As_q[0])
    chi2_l, S8_l = chi2_lcdm()
    N = len(fsig8_data)+2
    # ΔBIC honesto: SSEE k=3 (floteados) vs ΛCDM k=1 (solo amplitud); conservador
    dBIC = (chi2_med + 3*np.log(N)) - (chi2_l + 1*np.log(N))
    dchi2_fwd = chi2_fwd - chi2_l   # forward (k=0) vs ΛCDM (k=0): comparación limpia

    print("\n── Posterior ─────────────────────────────────────────────────────")
    print(f"  Ω_φDM = {Om_q[0]:.4f} +{Om_q[1]:.4f} -{Om_q[2]:.4f}  (algebraico {OM_PHIDM_CAN:.4f}, Δ={ten_Om:.2f}σ)")
    print(f"  m_φ   = {mp_q[0]:.2f} +{mp_q[1]:.2f} -{mp_q[2]:.2f} eV  (algebraico {M_PHI_CAN:.2f}, Δ={ten_mp:.2f}σ)")
    print(f"  A_s9  = {As_q[0]:.3f} +{As_q[1]:.3f} -{As_q[2]:.3f}  (Planck {AS9_MU})")
    print(f"  σ₈_eff = {s8_q[0]:.4f} +{s8_q[1]:.4f} -{s8_q[2]:.4f}  (forward 0.747)")
    print(f"  S₈_eff = {S8_q[0]:.4f} +{S8_q[1]:.4f} -{S8_q[2]:.4f}  (forward 0.758; KiDS {S8_KIDS})")
    print(f"\n── Bondad de ajuste ──────────────────────────────────────────────")
    print(f"  χ² forward(alg) = {chi2_fwd:.2f} / {N}")
    print(f"  χ² posterior    = {chi2_med:.2f} / {N}  (χ²_r={chi2_med/(N-3):.3f})")
    print(f"  χ² ΛCDM         = {chi2_l:.2f} / {N}   (S₈_ΛCDM={S8_l:.3f})")
    print(f"  Δχ²(forward−ΛCDM) = {dchi2_fwd:+.2f}   ΔBIC(post k=3 − ΛCDM k=1) = {dBIC:+.2f}")

    # figuras
    fig = corner.corner(flat, labels=[r"$\Omega_{\phi\rm DM}$", r"$m_\phi$ [eV]", r"$10^9A_s$"],
        truths=[OM_PHIDM_CAN, M_PHI_CAN, AS9_MU], truth_color="darkorange",
        quantiles=[0.16,0.5,0.84], show_titles=True, title_fmt=".3f", color="steelblue", bins=40, smooth=1.0)
    fig.suptitle(r"SSEE P6 $\phi$-DM v2 — emulador CLASS, priors planos", y=1.01)
    fig.savefig(OUTDIR+"fig_paper6_mcmc_v2_corner.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR+"fig_paper6_mcmc_v2_corner.png", dpi=150, bbox_inches="tight")

    # validación emulador
    validate_emulator(flat, n=10)
    print("\n✓ MCMC v2 completado.")
