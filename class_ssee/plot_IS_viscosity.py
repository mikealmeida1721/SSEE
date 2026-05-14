"""
SSEE-V3.6 — Fase 3: IS Viscosity CLASS Validation
Compara cs2=1 (nomira) vs cs2=0 (IS) y valida sigma_8=0.702 de Paper 5.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# ─── Cargar P(k) de CLASS ────────────────────────────────────────────────────
nomira  = np.loadtxt("output/ssee_v36_nomira__pk.dat",  comments="#")
IS_pk   = np.loadtxt("output/ssee_v36_IS__pk.dat",      comments="#")
lcdm_pk = np.loadtxt("output/lcdm_planck2018__pk.dat",  comments="#")

k_n, P_n = nomira[:,0],  nomira[:,1]
k_i, P_i = IS_pk[:,0],   IS_pk[:,1]
k_l, P_l = lcdm_pk[:,0], lcdm_pk[:,1]

# Interpolar al grid del IS run (rango común)
k_ref  = k_i
P_n_i  = np.interp(k_ref, k_n, P_n)
P_l_i  = np.interp(k_ref, k_l, P_l)

# ─── Función sigma_8 con filtro top-hat ──────────────────────────────────────
def W_tophat(k, R=8.0):
    x = k * R
    return 3.0 * (np.sin(x) - x * np.cos(x)) / x**3

def sigma8_CLASS(k, P, R=8.0):
    W  = W_tophat(k, R)
    integrand = k**2 * P * W**2 / (2 * np.pi**2)
    return np.sqrt(np.trapezoid(integrand, k))

sig_nomira = sigma8_CLASS(k_ref, P_n_i)
sig_IS     = sigma8_CLASS(k_ref, P_i)
sig_lcdm   = sigma8_CLASS(k_ref, P_l_i)

print("=" * 62)
print("  SSEE Fase 3 — IS Viscosity: sigma_8 CLASS")
print("=" * 62)
print(f"  σ₈  ΛCDM   (CLASS):  {sig_lcdm:.4f}")
print(f"  σ₈  nomira (cs2=1):  {sig_nomira:.4f}  ratio={sig_nomira/sig_lcdm:.4f}")
print(f"  σ₈  IS     (cs2≈0):  {sig_IS:.4f}  ratio={sig_IS/sig_lcdm:.4f}")
print()
print(f"  Paper 5 predicción:  G=0.866 → σ₈=0.702")
print(f"  CLASS IS ratio:      {sig_IS/sig_lcdm:.4f}")
diff = abs(sig_IS/sig_lcdm - 0.866)
print(f"  Δ vs Paper 5:        {diff:.4f}  ({diff/0.866*100:.1f}%)")
print()

# ─── ODE de crecimiento (idéntico a Paper 5 / ssee_paper6_verification.py) ───
Om_CDM = 0.111568
Om_b   = 0.048432
Om_m   = Om_CDM + Om_b   # 0.160
OmDE   = 1.0 - Om_m
w0, wa = -0.8399, -0.6699

def f_DE(a):
    return a**(-3*(1+w0+wa)) * np.exp(-3*wa*(1-a))

def E2_ssee(a):
    return max(Om_m * a**(-3) + OmDE * f_DE(a), 1e-30)

def E2_lcdm(a, Om_lcdm=0.3153):
    return max(Om_lcdm * a**(-3) + (1 - Om_lcdm), 1e-30)

def growth_rhs(lna, Y, E2_func, Om_func):
    a  = np.exp(lna)
    e2 = E2_func(a)
    Om = Om_func(a) / e2
    # dlnE/dlna ≈ (3/2) for matter dom; general:
    da  = 1e-5
    dlnE2_dlna = (E2_func(a + da) - E2_func(a - da)) / (2*da) * a / e2
    h_par = 0.5 * dlnE2_dlna
    D, Dp = Y
    return [Dp, -(2 + h_par)*Dp + 1.5*Om*D]

lna_span = (np.log(0.01), 0.0)
IC = [0.01, 0.01]
opts = dict(method="DOP853", rtol=1e-9, atol=1e-11, dense_output=True)

sol_ssee  = solve_ivp(growth_rhs, lna_span, IC,
                      args=(E2_ssee, lambda a: Om_m * a**(-3)), **opts)
sol_lcdm  = solve_ivp(growth_rhs, lna_span, IC,
                      args=(E2_lcdm, lambda a: 0.3153 * a**(-3)), **opts)

G_ODE = sol_ssee.y[0,-1] / sol_lcdm.y[0,-1]
print(f"  G (ODE growth factor):  {G_ODE:.4f}  [Paper 5: 0.866]")
print(f"  σ₈_ODE:                 {G_ODE * 0.811:.4f}  [Paper 5: 0.702]")
print()

# ─── Tabla resumen ────────────────────────────────────────────────────────────
header = f"{'Modelo':<20} {'σ₈':>8} {'G=σ₈/σ₈_ΛCDM':>14} {'Δ Paper5 (%)':>14}"
print(header)
print("-"*60)
rows = [
    ("ΛCDM (CLASS)",    sig_lcdm,   sig_lcdm/sig_lcdm),
    ("SSEE nomira cs2=1", sig_nomira, sig_nomira/sig_lcdm),
    ("SSEE IS cs2≈0",  sig_IS,     sig_IS/sig_lcdm),
    ("Paper 5 (ODE)",  G_ODE*0.811, G_ODE),
]
for name, s8, G in rows:
    d = (G - 0.866)/0.866*100 if name != "ΛCDM (CLASS)" else 0.0
    print(f"  {name:<18} {s8:>8.4f} {G:>14.4f} {d:>13.1f}%")
print()

# ─── Figura 1: P(k) comparación ───────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                gridspec_kw={"height_ratios": [2, 1], "hspace": 0.05})

ax1.loglog(k_ref, P_l_i,  "k-",  lw=1.8, label=r"$\Lambda$CDM Planck 2018")
ax1.loglog(k_ref, P_n_i,  "b--", lw=1.8, label=r"SSEE nomira ($c_s^2=1$, $\Omega_m=0.160$)")
ax1.loglog(k_ref, P_i,    "r-",  lw=2.0, label=r"SSEE IS ($c_s^2\approx0$, Paper 5)")

ax1.axvline(0.493, ls=":", color="gray", lw=1.2, label=r"$k_{\rm fs}=0.493$ h/Mpc")
ax1.set_ylabel(r"$P(k)$ [$h^{-3}$ Mpc$^3$]", fontsize=12)
ax1.set_title("SSEE-V3.6 Fase 3 — IS Viscosity: CLASS P(k)", fontsize=13)
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3, which="both")

ratio_IS    = P_i    / P_l_i
ratio_nomira = P_n_i / P_l_i
ax2.semilogx(k_ref, ratio_IS,    "r-",  lw=2.0, label=r"IS ($c_s^2\approx0$) / $\Lambda$CDM")
ax2.semilogx(k_ref, ratio_nomira,"b--", lw=1.8, label=r"nomira ($c_s^2=1$) / $\Lambda$CDM")
ax2.axhline(1.0,   ls="-",  color="k",    lw=0.8)
ax2.axhline(0.866**2, ls=":", color="red", lw=1.0,
            label=fr"Paper 5: $G^2={0.866**2:.3f}$")
ax2.axvline(0.493, ls=":", color="gray", lw=1.2)
ax2.set_xlabel(r"$k$ [$h$ Mpc$^{-1}$]", fontsize=12)
ax2.set_ylabel(r"$P_{\rm SSEE}/P_{\Lambda\rm CDM}$", fontsize=11)
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3, which="both")
ax2.set_ylim(0.3, 1.3)

fig.savefig("output/ssee_IS_Pk_comparison.pdf", bbox_inches="tight", dpi=150)
fig.savefig("output/ssee_IS_Pk_comparison.png", bbox_inches="tight", dpi=150)
print("  → output/ssee_IS_Pk_comparison.pdf")

# ─── Figura 2: cs2 scan — σ₈ vs cs2 ─────────────────────────────────────────
cs2_vals = np.array([0.001, 0.01, 0.1, 0.5, 1.0])
sig_vals  = np.array([sig_IS, np.nan, np.nan, np.nan, sig_nomira])

fig2, ax = plt.subplots(figsize=(7, 4))
ax.scatter([0.001, 1.0], [sig_IS, sig_nomira], color=["red", "blue"],
           zorder=5, s=80, label="CLASS runs")
ax.axhline(0.702, ls="--", color="red",   lw=1.5, label=r"Paper 5: $\sigma_8=0.702$")
ax.axhline(0.811, ls="--", color="black", lw=1.0, label=r"$\Lambda$CDM: $\sigma_8=0.811$")
ax.set_xscale("log")
ax.set_xlabel(r"$c_s^2$", fontsize=12)
ax.set_ylabel(r"$\sigma_8$", fontsize=12)
ax.set_title(r"SSEE $\sigma_8$ vs $c_s^2$ (IS: $c_s^2\to0$)", fontsize=12)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

fig2.savefig("output/ssee_IS_sigma8_cs2.pdf", bbox_inches="tight", dpi=150)
fig2.savefig("output/ssee_IS_sigma8_cs2.png", bbox_inches="tight", dpi=150)
print("  → output/ssee_IS_sigma8_cs2.pdf")
print()
print("Fase 3 completada. ✅")
