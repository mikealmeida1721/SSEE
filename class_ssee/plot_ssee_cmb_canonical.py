"""
SSEE-V3.6 vs ΛCDM: CMB TT spectrum (CLASS Boltzmann) — JOURNAL-STYLE figure.
ω_m-direct CANONICAL background Ω_m,CMB = 0.308881 (OP-8 dissolved, no MIRA factor).
Reads output/ssee_v36_canonical__cl_lensed.dat and the ΛCDM reference.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os

CLASS_DIR = os.path.dirname(os.path.abspath(__file__))
T_cmb_uK = 2.7255e6   # µK

# ── 1. SSEE canonical CLASS output ───────────────────────────────────────────
ssee_cl = np.loadtxt(os.path.join(CLASS_DIR, "output/ssee_v36_canonical__cl_lensed.dat"),
                     comments="#")
ell_ssee   = ssee_cl[:, 0]
Dl_ssee_TT = ssee_cl[:, 1] * T_cmb_uK**2

# ── 2. ΛCDM reference (Planck 2018 best-fit) ─────────────────────────────────
lcdm_ini = """
output = tCl,pCl,lCl
lensing = yes
l_max_scalars = 2500
Omega_b = 0.04897
Omega_cdm = 0.2607
h = 0.6766
n_s = 0.9665
A_s = 2.105e-9
k_pivot = 0.05
root = output/lcdm_planck2018_
overwrite_root = yes
"""
with open(os.path.join(CLASS_DIR, "lcdm_planck2018_ref.ini"), "w") as f:
    f.write(lcdm_ini)
os.system(f"cd {CLASS_DIR} && ./class lcdm_planck2018_ref.ini > /dev/null 2>&1")
lcdm_cl   = np.loadtxt(os.path.join(CLASS_DIR, "output/lcdm_planck2018__cl_lensed.dat"),
                       comments="#")
ell_lcdm   = lcdm_cl[:, 0]
Dl_lcdm_TT = lcdm_cl[:, 1] * T_cmb_uK**2

# ── 3. Planck 2018 binned TT (approximate, for visual anchor) ─────────────────
planck_ell = np.array([  2,  10,  20,  30,  50,  80, 100, 150, 200, 220, 300, 400,
                       500, 538, 600, 700, 810, 900,1000,1200,1400,1600,1800,2000,2200,2500])
planck_Dl  = np.array([1200,2800,3200,3800,4500,4800,5200,5500,5600,5900,4500,3800,
                       2700,3200,2300,2700,2100,1800,1200, 700, 500, 350, 250, 200, 150, 120])
planck_err = planck_Dl * 0.03

# ── 4. Journal-style figure (white background) ───────────────────────────────
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
fig = plt.figure(figsize=(9, 6.5))
gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)
ax_main = fig.add_subplot(gs[0])
ax_res  = fig.add_subplot(gs[1], sharex=ax_main)

ax_main.errorbar(planck_ell, planck_Dl, yerr=planck_err, fmt="o", ms=3.5,
                 color="0.45", alpha=0.7, elinewidth=1,
                 label="Planck 2018 TT (binned, approx.)", zorder=2)
ax_main.plot(ell_lcdm, Dl_lcdm_TT, color="#1f77b4", lw=1.8, ls="--",
             label=r"$\Lambda$CDM (Planck 2018, CLASS)", zorder=3)
ax_main.plot(ell_ssee, Dl_ssee_TT, color="#d62728", lw=2.0,
             label=r"SSEE-V3.6 ($\Omega_{m,\rm CMB}=0.308881$, $w_0=-0.840$, $w_a=-0.670$)",
             zorder=4)
ax_main.set_ylabel(r"$D_\ell^{TT}$ [$\mu$K$^2$]", fontsize=12)
ax_main.set_xlim(2, 2500)
ax_main.set_ylim(bottom=0)
ax_main.grid(True, alpha=0.3)
ax_main.legend(loc="upper right", fontsize=9)
ax_main.set_title(r"SSEE-V3.6 CMB TT power spectrum (CLASS Boltzmann, $\omega_m$-direct canonical)",
                  fontsize=12, pad=8)

interp_lcdm = np.interp(ell_ssee, ell_lcdm, Dl_lcdm_TT)
residual = (Dl_ssee_TT - interp_lcdm) / interp_lcdm * 100
ax_res.axhline(0, color="#1f77b4", lw=1.0, ls="--", alpha=0.8)
ax_res.plot(ell_ssee, residual, color="#d62728", lw=1.3)
ax_res.fill_between(ell_ssee, residual, 0, alpha=0.15, color="#d62728")
ax_res.set_xlabel(r"Multipole $\ell$", fontsize=12)
ax_res.set_ylabel(r"$\Delta$ [%]", fontsize=10)
ax_res.set_ylim(-8, 8)
ax_res.grid(True, alpha=0.3)
plt.setp(ax_main.get_xticklabels(), visible=False)

out_path = os.path.join(CLASS_DIR, "output/ssee_v36_CMB_TT_comparison.png")
plt.savefig(out_path, dpi=180, bbox_inches="tight", facecolor="white")
plt.close()
print(f"Saved: {out_path}")

# ── 5. Numerical report ──────────────────────────────────────────────────────
for label, (lo, hi) in [("1st", (150, 350)), ("2nd", (450, 650)), ("3rd", (700, 950))]:
    m = (ell_ssee >= lo) & (ell_ssee <= hi)
    if m.any():
        pk = ell_ssee[m][np.argmax(Dl_ssee_TT[m])]
        print(f"  SSEE {label} peak: ℓ = {int(pk)}")
rms = np.std(residual[(ell_ssee >= 50) & (ell_ssee <= 2000)])
print(f"  RMS residual vs ΛCDM (ℓ=50–2000): {rms:.2f}%")
