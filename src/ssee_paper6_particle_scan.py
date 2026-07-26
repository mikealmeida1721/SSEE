"""
PASO 1 (FIT explícito) — partícula SUELTA: ¿qué m_φ prefieren los datos?
========================================================================
Método de Mike (2026-06-18): el fondo ω_m-directo está FIJO (Ω_m=0.308881,
Ω_φDM=0.14889). Lo único que mueve S8 es la masa m_φ. Soltamos m_φ como
parámetro libre, corremos CLASS two-sector para una grilla, y encontramos el
m_φ que los datos prefieren (S8 = KiDS-1000 = 0.758 → 0.00σ).
Esto es un FIT a S8 (declarado), NO forward. El paso 2 (separado) pregunta si
el diccionario SSEE tiene un multiplicador simple en ese m_φ preferido.

RESULTADO CANÓNICO de esta exploración (2026-07-10, ν-closure C=93.14):
m_φ = 40.70 eV = Σm_ν·(SOLAR²·KRYSTOS), S8 = 0.758 (0.04σ KiDS-1000).
Las referencias a 42.47/41.02 abajo son candidatos SUPERSEDED del linaje previo.
Fuente autoritativa: VERIFICATION_LEDGER.md + ssee_paper6_canonical_particle.py.
"""
import numpy as np, subprocess, os
from scipy.interpolate import interp1d

_REPO=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D=os.path.join(_REPO,"class_ssee")+os.sep; OUT=D+"output"+os.sep; CLASS=os.path.join(D,"class")
T_NU=0.71611; H=0.67962
phi=(1+5**0.5)/2; pi=np.pi; OMEGA=pi+phi
KAL=(pi+phi)/2+pi; N_S=1-phi**-7
# fondo ω_m-directo FIJO
Omega_b=(pi-phi)/(3*OMEGA**2)/H**2
om_b_h2=(pi-phi)/(3*OMEGA**2); om_c_h2=KAL*om_b_h2*N_S
Smnu=0.06849; om_nu_h2=Smnu/93.14   # C_nu=93.14 PDG (era 0.06902)
Om_m=(om_b_h2+om_c_h2+om_nu_h2)/H**2     # 0.308881
Om_cold=0.160; Om_hot=Om_m-Om_cold       # 0.14889
Om_cdm=Om_cold-Omega_b
om_ncdm_h2=Om_hot*H**2                    # FIJO (abundancia fija; T_φ se ajusta por m_φ)

KIDS_S8=0.758; KIDS_ERR=0.024            # KiDS-1000 (Asgari+2021) — EXTERNO
S8_target_sig8 = KIDS_S8/np.sqrt(Om_m/0.3)   # σ8 que daría S8=KiDS

def sigma8(k,P,R=8.0):
    x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(k**2*P*W**2/(2*np.pi**2),k))

def run_class(m_phi, T_phi, Nur, tag):
    lines=["output = mPk","P_k_max_h/Mpc = 5","z_pk = 0.0",
           f"h = {H}", f"Omega_b = {Omega_b}",
           "Omega_fld = 0.8399","w0_fld = -0.8399","wa_fld = -0.6699","cs2_fld = 1.0",
           "n_s = 0.96556","A_s = 2.1e-9","k_pivot = 0.05",
           "write_warnings = no","overwrite_root = yes", f"root = output/{tag}_",
           f"Omega_cdm = {Om_cdm:.5f}", f"N_ur = {Nur:.4f}",
           "N_ncdm = 1", f"m_ncdm = {m_phi:.4f}", f"T_ncdm = {T_phi:.5f}"]
    open(f"/tmp/{tag}.ini","w").write("\n".join(lines)+"\n")
    subprocess.run([CLASS,f"/tmp/{tag}.ini"],cwd=D,capture_output=True,text=True,timeout=600)
    f=OUT+f"{tag}__pk.dat"
    if not os.path.exists(f): return None
    d=np.loadtxt(f); k,P=d[:,0],d[:,1]
    if not(np.all(np.isfinite(k)) and np.all(np.isfinite(P))): return None
    return sigma8(k,P)

def kfs_analytic(m_phi, T_phi):
    return 0.018*np.sqrt(Om_m)*m_phi*(T_NU/T_phi)

print("="*78)
print("  PASO 1 (FIT) — partícula SUELTA: ¿qué m_φ prefieren los datos (S8)?")
print("="*78)
print(f"  Fondo FIJO: Ω_m=0.308881, Ω_φDM=0.14889, Σm_ν=0.069 eV (todo SSEE)")
print(f"  Objetivo (EXTERNO): KiDS-1000 S8={KIDS_S8}±{KIDS_ERR} → σ8_target={S8_target_sig8:.4f}")
print(f"  Referencia: m_φ=42.47 (PYROS·VITA·MIKA) dio S8=0.765 (0.24σ)")
print("-"*78)
print(f"  {'m_φ(eV)':>8} {'mult M':>8} {'k_fs':>7} {'σ8':>8} {'S8':>8} {'σ(KiDS)':>9}")

grid = [28, 32, 35, 38, 40, 42.47, 45, 50]
rows=[]
for m in grid:
    T_phi=T_NU*(om_ncdm_h2*93.14/m)**(1.0/3.0)
    Nur=max(0.01,3.046-(T_phi/T_NU)**4)
    s8=run_class(m, T_phi, Nur, f"scan_{int(m*100)}")
    if s8 is None:
        print(f"  {m:8.2f}  [FALLO CLASS]"); continue
    S8=s8*np.sqrt(Om_m/0.3); M=m/Smnu; kfs=kfs_analytic(m,T_phi)
    tens=abs(S8-KIDS_S8)/KIDS_ERR
    mark=" <- PYROS·VITA·MIKA" if abs(m-42.47)<0.1 else ""
    print(f"  {m:8.2f} {M:8.1f} {kfs:7.3f} {s8:8.4f} {S8:8.4f} {tens:8.2f}σ{mark}")
    rows.append((m,M,kfs,s8,S8,tens))

# interpolar m_φ donde S8 = KiDS (data-preferred)
rows=np.array(rows)
m_arr,S8_arr=rows[:,0],rows[:,4]
if S8_arr.min()<=KIDS_S8<=S8_arr.max():
    f_inv=interp1d(S8_arr,m_arr)
    m_best=float(f_inv(KIDS_S8))
    M_best=m_best/Smnu
    print("-"*78)
    print(f"  DATA-PREFERRED (S8=KiDS, 0.00σ):  m_φ = {m_best:.3f} eV  →  multiplicador M = {M_best:.2f}")
    print(f"  (vs 615.33 actual → m_φ=42.47, 0.24σ). El paso 2 busca M={M_best:.1f} en el diccionario.")
else:
    print(f"  KiDS S8={KIDS_S8} fuera del rango escaneado [{S8_arr.min():.3f},{S8_arr.max():.3f}]")
print("="*78)
