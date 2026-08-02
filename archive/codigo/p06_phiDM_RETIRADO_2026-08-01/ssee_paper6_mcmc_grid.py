import os
#!/usr/bin/env python3
"""
Paper 6 — Generador de grilla CLASS para el MCMC two-sector φ-DM (rediseño riguroso).
=====================================================================================
Reemplaza el σ₈ "single-k" de juguete (ssee_paper6_mcmc.py, RETIRADO) por σ₈ y
fσ₈(z) calculados con CLASS REAL + integral top-hat (fórmula de Peebles).

Estrategia (gold-standard, validable):
  - Grilla 2D sobre (Ω_φDM, m_φ).  El background SSEE (H₀, w₀, wₐ, Ω_m,CMB) es FIJO.
  - Ω_m,CMB = 0.308881 fijo (anclado por CMB, independiente).  Al flotar Ω_φDM,
    el sector frío Ω_cdm = Ω_m − Ω_φDM − Ω_b se ajusta (la partición es lo que se prueba).
  - m_φ fija T_φ vía relic constraint (idéntico al forward ssee_paper6_canonical_particle.py).
  - A_s NO va en la grilla: σ₈ ∝ √A_s exacto en teoría lineal → se escala analíticamente
    en el MCMC (cero error de interpolación en la dirección de amplitud).
  - Para cada punto: 1 corrida CLASS con P(k) a múltiples z → σ₈(z_i) top-hat.
    fσ₈(z) = f(z)·σ₈(z), con f(z) = −(1+z) dlnσ₈/dz (diferencia finita en grilla z densa).

Anti-ambigüedad de archivos multi-z de CLASS: σ₈(z) es monótona decreciente, así que
ordeno σ₈ de todos los archivos y mapeo σ₈↓ ↔ z↑.  No dependo del orden de CLASS.

Uso:
  python ssee_paper6_mcmc_grid.py --test     # solo punto canónico, valida σ₈≈0.747
  python ssee_paper6_mcmc_grid.py --run       # grilla completa → .npz en HDD
"""
import numpy as np, subprocess, os, sys, glob, time

# ── Rutas portables (HDD para outputs pesados si existe; NUNCA root SSD) ──────
# CLASS relativo al repo; outputs vía SSEE_DATA_DIR (default: HDD si existe, si no results/)
_REPO   = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D       = os.path.join(_REPO, "class_ssee") + os.sep
CLASS   = os.path.join(D, "class")
_SSEE_DATA = os.environ.get("SSEE_DATA_DIR") or ("/mnt/datos/SSEE_data" if os.path.isdir("/mnt/datos") else os.path.join(_REPO, "results", "data"))
OUTDIR  = os.path.join(_SSEE_DATA, "mcmc/grid_class") + os.sep   # P(k) temporales de la grilla
GRIDNPZ = os.path.join(_SSEE_DATA, "mcmc/p6_grid.npz")           # producto final
os.makedirs(OUTDIR, exist_ok=True)

# ── Constantes (φ,π) y fondo ω_m-directo — idéntico al forward canónico ───────
phi   = (1 + 5**0.5) / 2
pi    = np.pi
OMEGA = pi + phi
KAL   = (pi + phi)/2 + pi
TRIAL = 3*(phi + (pi + phi)/2)
N_S   = 1 - phi**-7                       # 0.96556
H     = 0.67962
h2    = H**2
T_NU  = 0.71611

om_b_h2 = (pi - phi)/(3*OMEGA**2)         # 0.02242
Omega_b = om_b_h2 / h2                     # 0.04854
om_c_h2 = KAL*om_b_h2*N_S                  # 0.11951
R2      = OMEGA/(KAL*TRIAL)
Smnu    = R2*0.960318
om_nu_h2= Smnu/93.14
om_m_h2 = om_b_h2 + om_c_h2 + om_nu_h2     # 0.14267
Om_m    = om_m_h2 / h2                     # 0.308881  Ω_m,CMB (fijo, anclado CMB)

# Punto canónico (forward): Ω_φDM=0.149, m_φ=40.70 eV
SOLAR   = (pi + phi)/2 + KAL               # φ+2π = 7.9012
KRYSTOS_V = phi + pi + OMEGA               # 9.5192 — padres {φ,π,Ω}, NO 2Ω (colapso)
MULT    = SOLAR**2 * KRYSTOS_V             # 594.28
M_PHI_CAN = Smnu*MULT                      # 40.70 eV
OM_PHIDM_CAN = Om_m - 0.160                # 0.149

# ── Redshifts: 6 RSD + densos para derivar f(z) ──────────────────────────────
Z_RSD = np.array([0.067, 0.150, 0.380, 0.510, 0.610, 1.480])
# grilla z densa (incluye 0 y cubre hasta 1.55 para diferencia finita en z=1.48)
Z_DENSE = np.unique(np.concatenate([[0.0], np.linspace(0.02, 1.55, 30)]))

def sigma8_tophat(k, P, R=8.0):
    """σ₈ = sqrt[(1/2π²)∫k²P(k)W²(kR)dk], W top-hat (Peebles 1980)."""
    x = k*R
    W = 3*(np.sin(x) - x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(k**2 * P * W**2 / (2*np.pi**2), k))

def run_class_multiz(tag, Om_phiDM, m_phi, z_list):
    """Corre CLASS two-sector; devuelve dict {z: σ₈(z)} usando mapeo σ₈↓↔z↑."""
    om_ncdm_h2 = Om_phiDM * h2
    T_phi = T_NU*(om_ncdm_h2*93.14/m_phi)**(1.0/3.0)
    Nur   = max(0.01, 3.046 - (T_phi/T_NU)**4)
    Om_cdm = Om_m - Om_phiDM - Omega_b     # sector frío restante
    if Om_cdm <= 0:
        return None
    z_sorted = np.sort(z_list)             # ascendente
    z_str = ",".join(f"{z:.4f}" for z in z_sorted)
    root = f"{tag}_"
    lines = [
        "output = mPk", "P_k_max_h/Mpc = 5", f"z_pk = {z_str}",
        f"h = {H}", f"Omega_b = {Omega_b:.6f}", f"Omega_cdm = {Om_cdm:.6f}",
        "Omega_fld = 0.8399", "w0_fld = -0.8399", "wa_fld = -0.6699", "cs2_fld = 1.0",
        f"n_s = {N_S:.5f}", "A_s = 2.1e-9", "k_pivot = 0.05",
        "N_ncdm = 1", f"m_ncdm = {m_phi:.5f}", f"T_ncdm = {T_phi:.6f}", f"N_ur = {Nur:.4f}",
        "write_warnings = no", "overwrite_root = yes", f"root = {OUTDIR}{root}",
    ]
    ini = f"/tmp/{tag}.ini"
    open(ini, "w").write("\n".join(lines) + "\n")
    # limpiar outputs viejos de este tag (total y _cb)
    for f in glob.glob(OUTDIR + f"{root}*.dat"):
        os.remove(f)
    r = subprocess.run([CLASS, ini], cwd=D, capture_output=True, text=True, timeout=600)
    # solo P(k) TOTAL: '<root>z{i}_pk.dat' (excluye '_pk_cb.dat' cold+baryons)
    files = sorted(f for f in glob.glob(OUTDIR + f"{root}*_pk.dat")
                   if not f.endswith("_pk_cb.dat"))
    if len(files) < len(z_sorted):
        print(f"  [FALLO {tag}] archivos={len(files)} esperados={len(z_sorted)}  {r.stderr[-200:]}")
        return None
    s8_vals = []
    for f in files:
        d = np.loadtxt(f)
        if d.ndim != 2 or not np.all(np.isfinite(d)):
            return None
        s8_vals.append(sigma8_tophat(d[:, 0], d[:, 1]))
    s8_vals = np.array(s8_vals)
    if len(s8_vals) != len(z_sorted):
        # tomar los len(z) mayores si CLASS dejó archivos extra
        s8_vals = np.sort(s8_vals)[::-1][:len(z_sorted)]
    # mapeo: σ₈ descendente ↔ z ascendente
    s8_desc = np.sort(s8_vals)[::-1]
    return {float(z): float(s8) for z, s8 in zip(z_sorted, s8_desc)}, T_phi, Nur, Om_cdm

def test_canonical():
    print("="*70)
    print("  TEST punto canónico — debe reproducir σ₈(0)≈0.747 del forward")
    print("="*70)
    print(f"  Ω_m,CMB={Om_m:.5f}  Ω_φDM={OM_PHIDM_CAN:.5f}  m_φ={M_PHI_CAN:.4f} eV")
    t0 = time.time()
    res = run_class_multiz("test_can", OM_PHIDM_CAN, M_PHI_CAN, Z_DENSE)
    if res is None:
        print("  ✗ CLASS falló"); return
    s8z, T_phi, Nur, Om_cdm = res
    print(f"  T_φ={T_phi:.5f}  N_ur={Nur:.4f}  Ω_cdm={Om_cdm:.5f}  ({time.time()-t0:.1f}s)")
    print(f"  σ₈(z=0)   = {s8z[0.0]:.4f}   (forward CLASS: 0.7470)")
    s8_0 = s8z[0.0]
    print(f"  S₈(z=0)   = {s8_0*np.sqrt(Om_m/0.3):.4f}   (forward: 0.758)")
    ok = abs(s8_0 - 0.7470) < 0.01
    print(f"  {'✓ COINCIDE' if ok else '✗ NO coincide'} con el forward (tol 0.01)")
    return ok

def run_grid():
    """Barre (Ω_φDM, m_φ); guarda σ₈(z_dense) por punto en .npz (A_s=2.1e-9 ref)."""
    # rangos: cubren el punto algebraico (0.149, 41 eV) y dejan que el dato deambule.
    # Ω_φDM acotado por Ω_cdm = Ω_m − Ω_φDM − Ω_b > 0  →  Ω_φDM < 0.260
    grid_phidm = np.linspace(0.04, 0.25, 22)
    grid_mphi  = np.logspace(np.log10(8.0), np.log10(100.0), 20)  # log: k_fs ∝ m_φ
    nP, nM, nZ = len(grid_phidm), len(grid_mphi), len(Z_DENSE)
    s8z_grid = np.full((nP, nM, nZ), np.nan)
    Tphi_grid = np.full((nP, nM), np.nan)
    print("="*70)
    print(f"  GRILLA CLASS  {nP}×{nM} = {nP*nM} puntos × {nZ} z   (A_s ref=2.1e-9)")
    print(f"  Ω_φDM ∈ [{grid_phidm[0]:.3f}, {grid_phidm[-1]:.3f}]   "
          f"m_φ ∈ [{grid_mphi[0]:.1f}, {grid_mphi[-1]:.1f}] eV")
    print("="*70)
    t0 = time.time(); done = 0; failed = 0
    for i, opd in enumerate(grid_phidm):
        for j, mp in enumerate(grid_mphi):
            res = run_class_multiz(f"g_{i}_{j}", opd, mp, Z_DENSE)
            if res is None:
                failed += 1; continue
            s8z, T_phi, Nur, Om_cdm = res
            zs = np.sort(Z_DENSE)
            s8z_grid[i, j, :] = np.array([s8z[float(z)] for z in zs])
            Tphi_grid[i, j] = T_phi
            done += 1
        el = time.time()-t0
        print(f"  fila {i+1}/{nP}  (Ω_φDM={opd:.3f})  ok={done} fail={failed}  {el:.0f}s")
    np.savez_compressed(
        GRIDNPZ, grid_phidm=grid_phidm, grid_mphi=grid_mphi,
        z_dense=np.sort(Z_DENSE), s8z_grid=s8z_grid, Tphi_grid=Tphi_grid,
        Om_m=Om_m, h=H, Omega_b=Omega_b, As_ref=2.1e-9,
        m_phi_can=M_PHI_CAN, om_phidm_can=OM_PHIDM_CAN)
    print(f"\n  ✓ Grilla guardada: {GRIDNPZ}  ({done} ok, {failed} fail, {time.time()-t0:.0f}s)")
    # sanity: σ₈(0) en el vértice más cercano al canónico
    ic = np.argmin(abs(grid_phidm - OM_PHIDM_CAN))
    jc = np.argmin(abs(grid_mphi - M_PHI_CAN))
    print(f"  σ₈(0) vértice≈canónico (Ω_φDM={grid_phidm[ic]:.3f}, m_φ={grid_mphi[jc]:.1f}) "
          f"= {s8z_grid[ic, jc, 0]:.4f}  (forward 0.747)")

if __name__ == "__main__":
    if "--test" in sys.argv:
        test_canonical()
    elif "--run" in sys.argv:
        run_grid()
    else:
        print("Usa --test primero. Luego --run para la grilla completa.")
