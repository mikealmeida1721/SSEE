"""
PARTÍCULA CANÓNICA φ-DM (Vía 2 + multiplicador SOLAR²·KRYSTOS).
Cadena forward, CERO fiteo (reframe ω_m-directo 2026-06-18, OP-8 cerrado;
mecanismo g²·v SOLAR²·KRYSTOS adoptado 2026-06-19, OP-17):
    R2 = Ω/(KAL·TRIAL)
    Σm_ν = R2 × ω_b × C_ν/τ_Π  (ν-closure C=93.14; ≈ R2 × 0.9530 eV)
    m_φ = Σm_ν × (SOLAR²·KRYSTOS)         [= 594.28; antes PYROS·VITA·MIKA=615.33]
    T_φ SALE del relic constraint.
Escenario de fondo: NO hay factor materia. Ω_m,CMB = ω_m/h² = 0.30889; el sector
φ-DM es la diferencia Ω_φDM = Ω_m,CMB − Ω_m,dyn = 0.149 (antes 0.160 vía MIRA).
Calcula las DOS predicciones falseables (se ACTUALIZAN con el escenario nuevo):
    (1) alpha  — OUTPUT de CLASS: fit Viel al cociente P_partícula/P_frío.
                 (NO se impone; emerge del Boltzmann real).
    (2) k_fs   — Dodelson-Widrow free-streaming a partir de m_φ.
"""
import numpy as np, subprocess, os
from scipy.optimize import curve_fit

_REPO=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))   # portable: raíz del repo
D=os.path.join(_REPO,"class_ssee")+os.sep; OUT=D+"output"+os.sep; CLASS=os.path.join(D,"class")
T_NU=0.71611; H=0.67962
# --- constantes (φ,π) ---
phi=(1+5**0.5)/2; pi=np.pi
OMEGA=pi+phi                  # Ω_DNAV
KAL=(pi+phi)/2+pi             # KAL0
TRIAL=3*(phi+(pi+phi)/2)      # TRIAL
AURA=(3*phi+pi)/2             # AURA = 2*MIRA
N_S=1-phi**-7                 # n_s = 0.96556
# multiplicador canónico ω_m-directo (SOLAR²·KRYSTOS — mecanismo g²·v, OP-9/17)
# SOLAR = BIAL+KAL = φ+2π (linaje radiativo);  KRYSTOS_V = φ+π+Ω (padres, anclado por wₐ)
BIAL=(pi+phi)/2               # 2.3798
SOLAR=BIAL+KAL                # 7.9012 (= φ+2π)
KRYSTOS_V=phi+pi+OMEGA        # 9.5192 — padres {φ,π,Ω}, NO 2Ω (colapso); gemelo IGNIS=π+PYROS

# --- fondo ω_m-directo (sin factor materia) ---
Omega_b=(pi-phi)/(3*OMEGA**2)/H**2    # ω_b/h² algebraico = 0.04854
om_b_h2=(pi-phi)/(3*OMEGA**2)         # 0.02242
om_c_h2=KAL*om_b_h2*N_S               # 0.11951 (forward, ya en Paper 1)
R2=OMEGA/(KAL*TRIAL)
# Σm_ν = R₂·ω_b·C/(τ_Π H₀), transparente. C=93.14 eV = constante de cierre ν
# PRECISA (N_eff=3.046, la de Planck); demostrada desde 1ros principios en
# derive_nu_closure.py (deriva 94.06 instantáneo; 93.14 con corrección e⁺e⁻).
# Univaluada en toda la suite (antes convivía con 94.07; auditoría 2026-07-10).
C_nu=93.14
tau_Pi=KAL*OMEGA/TRIAL            # = KAL/(3·Ω_DE) = 2.191, número puro
Smnu=R2*om_b_h2*C_nu/tau_Pi       # = 0.0685 eV (antes 0.06902 con C≈93.85 baked)
om_nu_h2=Smnu/93.14
om_m_h2=om_b_h2+om_c_h2+om_nu_h2      # 0.14267
Om_m=om_m_h2/H**2                     # 0.30889  Ω_m,CMB (derivado)
Om_cold=0.160                         # sector dinámico (CDM frío + b), DESI
Om_hot=Om_m-Om_cold                   # 0.149  Ω_φDM (diferencia; antes 0.160 vía MIRA)
Om_cdm=Om_cold-Omega_b
om_ncdm_h2=Om_hot*H**2
mu=1.12

MULT=SOLAR**2*KRYSTOS_V
m_phi=Smnu*MULT

T_phi=T_NU*(om_ncdm_h2*93.14/m_phi)**(1.0/3.0)
Nur=max(0.01,3.046-(T_phi/T_NU)**4)

def load(f):
    d=np.loadtxt(f); return d[:,0],d[:,1]
def sigma8(k,P,R=8.0):
    x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(k**2*P*W**2/(2*np.pi**2),k))
def viel_mix(k,a,f):
    Tv=(1+(a*k)**(2*mu))**(-5.0/mu); return ((1-f)+f*Tv)**2

def run(tag, ncdm, z=0.0):
    lines=["output = mPk","P_k_max_h/Mpc = 5", f"z_pk = {z}",
           f"h = {H}", f"Omega_b = {Omega_b}",
           "Omega_fld = 0.8399","w0_fld = -0.8399","wa_fld = -0.6699","cs2_fld = 1.0",
           "n_s = 0.96556","A_s = 2.1e-9","k_pivot = 0.05",
           "write_warnings = no","overwrite_root = yes", f"root = output/{tag}_"]
    if ncdm:
        lines += [f"Omega_cdm = {Om_cdm:.5f}", f"N_ur = {Nur:.4f}",
                  "N_ncdm = 1", f"m_ncdm = {m_phi:.4f}", f"T_ncdm = {T_phi:.5f}"]
    else:
        lines += [f"Omega_cdm = {Om_m-Omega_b:.5f}", "N_ur = 3.046", "N_ncdm = 0"]
    open(f"/tmp/{tag}.ini","w").write("\n".join(lines)+"\n")
    r=subprocess.run([CLASS,f"/tmp/{tag}.ini"],cwd=D,capture_output=True,text=True,timeout=600)
    f=OUT+f"{tag}__pk.dat"
    if not os.path.exists(f):
        print(f"  [FALLO {tag}] {r.stderr[-300:]}"); return None
    k,P=load(f)
    if not(np.all(np.isfinite(k))and np.all(np.isfinite(P))):
        print(f"  [ROTO {tag}]"); return None
    return k,P

print("="*74)
print("  PARTÍCULA CANÓNICA φ-DM — alpha (OUTPUT) y k_fs")
print("="*74)
print(f"  R2 = Ω/(KAL·TRIAL)        = {R2:.6f}")
print(f"  Σm_ν                      = {Smnu:.5f} eV")
print(f"  multiplicador SOLAR²·KRYSTOS_V = {MULT:.4f}")
print(f"  m_φ                       = {m_phi:.4f} eV")
print(f"  Ω_m,CMB (ω_m/h², sin factor) = {Om_m:.5f}")
print(f"  Ω_φDM = Ω_m,CMB − 0.160      = {Om_hot:.5f}")
print(f"  T_φ (relic)               = {T_phi:.5f}  (= {T_phi/T_NU:.4f} T_ν)")
print(f"  N_ur                      = {Nur:.4f}")
print("-"*74)

kp,Pp = run("can_part", True)
kc,Pc = run("can_cold", False)

# --- alpha OUTPUT: fit Viel al cociente P_partícula/P_frío ---
kg=np.logspace(np.log10(0.05),np.log10(2.0),400)
rg=np.clip(np.interp(kg,kp,Pp)/np.interp(kg,kc,Pc),1e-6,1.2)
popt,_=curve_fit(viel_mix,kg,rg,p0=[1.6,0.5],bounds=([0.05,0.0],[60,1.0]),maxfev=40000)
alpha_out, frac_out = popt

# --- k_fs free-streaming ---
# (a) MEDIDA por CLASS: k donde el cociente P_part/P_frío cae a 0.5 (half-mode).
#     Robusto, independiente de fórmula. Es la escala de borrado real del Boltzmann.
half = np.interp(0.5, rg[::-1], kg[::-1]) if rg.min()<0.5 else np.nan
# (b) ANALÍTICO no-relativista (forma neutrino, corregida por T más frío):
#     k_fs ≈ 0.018 · Ω_m^(1/2) · (m_φ/eV) · (T_ν/T_φ)   h/Mpc
#     La partícula es más fría que el ν (T_φ/T_ν<1) → menos velocidad → k_fs mayor.
kfs_an = 0.018*np.sqrt(Om_m)*m_phi*(T_NU/T_phi)

print(f"  alpha (OUTPUT fit Viel)   = {alpha_out:.4f} Mpc/h   (frac={frac_out:.3f})")
print(f"      [φ = {phi:.4f} — comparar]")
print(f"  k_half (CLASS, ratio=0.5) = {half:.4f} h/Mpc   <-- borrado medido por el Boltzmann")
print(f"  k_fs (analítico no-rel)   = {kfs_an:.4f} h/Mpc   <-- predicción falseable DESI/Euclid")
s8=sigma8(kp,Pp); s8c=sigma8(kc,Pc)
print(f"  σ8 (partícula)            = {s8:.4f}")
print(f"  σ8 (frío, techo)          = {s8c:.4f}")
print("="*74)
