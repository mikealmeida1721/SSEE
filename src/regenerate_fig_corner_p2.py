"""
Regenera fig_corner_ssee_professional.pdf con la cadena MCMC-MIRA actualizada
(reemplaza el corner del prior Planck por el del prior MIRA).
"""
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import corner
import os

# Cargar cadena P2-MIRA
d = np.load("results/logs/mcmc_paper2_mira_ckpt.npz")
chain = d["chain"]
H0 = chain[:,0]
ob = chain[:,1]

# Estadísticas
H0_med = np.median(H0); H0_std = np.std(H0)
H0_p16, H0_p84 = np.percentile(H0, [16, 84])
ob_med = np.median(ob); ob_std = np.std(ob)

print(f"H0: {H0_med:.4f}  +{H0_p84-H0_med:.4f}/-{H0_med-H0_p16:.4f}")
print(f"Ob_h2: {ob_med:.5f}  ±{ob_std:.5f}")
print(f"N samples: {len(chain)}")

# Corner figure — mismo estilo que original
fig = corner.corner(
    chain,
    labels=[r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]", r"$\Omega_b h^2$"],
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 11},
    label_kwargs={"fontsize": 12},
    truths=[67.068, 0.02218],         # prior MIRA + BBN
    truth_color="#1f77b4",
    color="#333333",
    hist_kwargs={"density": True},
    levels=(0.393, 0.865),             # 1σ, 2σ contours (2D)
    plot_datapoints=False,
    fill_contours=True,
)
fig.suptitle(
    f"SSEE-V3.6 posterior (MIRA prior, $N_{{\\rm walkers}}\\times N_{{\\rm steps}}=100\\times25000$). "
    f"Blue lines: MIRA $H_0$=67.068 km/s/Mpc, BBN $\\Omega_b h^2$=0.02218.",
    y=1.02, fontsize=11,
)
os.makedirs("results/figures", exist_ok=True)
out = "results/figures/fig_corner_ssee_professional.pdf"
fig.savefig(out, bbox_inches="tight")
plt.close(fig)
print(f"Saved: {out}")
