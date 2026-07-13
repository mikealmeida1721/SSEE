"""
FILTRO HONESTO sobre la particula phi-DM de ~40 eV (las 3 vias de Mike).
Dos sectores: 0.160 CDM frio + 0.160 ncdm caliente (la particula).
m_phi = Sum_mnu * [(pi+phi)^4 + (phi+2pi)^2]  con Sum_mnu = R * 0.960318 eV.
T_phi SALE del relic constraint. CERO fiteo.

Dos pruebas que faltaban:
  (1) CMB TT: posiciones de los 3 primeros picos vs LCDM (un ~40 eV caliente
      podria desplazarlos). Tambien chi2 RMS relativo vs LCDM.
  (2) fsigma8(z): crecimiento via derivada de sigma8(z) de CLASS,
      comparado contra BOSS/DESI (refs ~0.45-0.47).
"""
import numpy as np, subprocess, os

_REPO=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D=os.path.join(_REPO,"class_ssee")+os.sep; OUT=D+"output"+os.sep; CLASS=os.path.join(D,"class")
T_NU=0.71611; H=0.67962; Om_m=0.3200; Omega_b=0.04897
Om_cold=0.160; Om_hot=0.160
Om_cdm=Om_cold-Omega_b
om_ncdm_h2=Om_hot*H**2

phi=(1+5**0.5)/2; pi=np.pi
OMEGA=pi+phi; MAR=phi+2*pi; KAL=(pi+phi)/2+pi
TRIAL=3*(phi+(pi+phi)/2); PYROS=OMEGA+phi
MULT=OMEGA**4+MAR**2          # 575.635

vias={"V1_estructura": pi/(KAL*MAR),
      "V2_estabilidad": OMEGA/(KAL*TRIAL),
      "V3_fuego":       KAL/(PYROS*TRIAL)}

def load_pk(f):
    d=np.loadtxt(f); return d[:,0],d[:,1]
def sigma8(k,P,R=8.0):
    x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3
    return np.sqrt(np.trapezoid(k**2*P*W**2/(2*np.pi**2),k))

def base_lines(tag, ncdm_block, want_cmb, z):
    L=[f"h = {H}", f"Omega_b = {Omega_b}",
       "Omega_fld = 0.8399","w0_fld = -0.8399","wa_fld = -0.6699","cs2_fld = 1.0",
       "n_s = 0.96556","A_s = 2.1e-9","k_pivot = 0.05",
       "write_warnings = no","overwrite_root = yes", f"root = output/{tag}_"]
    if want_cmb:
        L=["output = tCl,pCl,lCl,mPk","lensing = yes","l_max_scalars = 2700",
           "P_k_max_h/Mpc = 5", f"z_pk = {z}"]+L
    else:
        L=["output = mPk","P_k_max_h/Mpc = 5", f"z_pk = {z}"]+L
    return L+ncdm_block

def run(tag, m_phi, T_phi, Nur, want_cmb=False, z=0.0):
    if m_phi is None:   # LCDM-like single cold sector
        nb=[f"Omega_cdm = {0.320-Omega_b:.5f}", "N_ur = 3.046", "N_ncdm = 0"]
    else:
        nb=[f"Omega_cdm = {Om_cdm:.5f}", f"N_ur = {Nur:.4f}",
            "N_ncdm = 1", f"m_ncdm = {m_phi:.4f}", f"T_ncdm = {T_phi:.5f}"]
    open(f"/tmp/{tag}.ini","w").write("\n".join(base_lines(tag,nb,want_cmb,z))+"\n")
    r=subprocess.run([CLASS,f"/tmp/{tag}.ini"],cwd=D,capture_output=True,text=True,timeout=900)
    out={}
    pk=OUT+f"{tag}__pk.dat"
    if os.path.exists(pk):
        k,P=load_pk(pk)
        if np.all(np.isfinite(k)) and np.all(np.isfinite(P)): out["pk"]=(k,P)
    if want_cmb:
        cl=OUT+f"{tag}__cl_lensed.dat"
        if os.path.exists(cl):
            d=np.loadtxt(cl)
            if np.all(np.isfinite(d)): out["cl"]=(d[:,0],d[:,1])  # l, TT (l(l+1)Cl/2pi)
    if not out: print(f"  [FALLO {tag}] {r.stderr[-200:]}")
    return out

def peaks(l,TT):
    """3 primeros maximos locales en TT por encima de l=100."""
    m=l>=100; l=l[m]; TT=TT[m]
    pk=[]
    for i in range(2,len(TT)-2):
        if TT[i]>TT[i-1] and TT[i]>TT[i+1] and TT[i]>TT[i-2] and TT[i]>TT[i+2]:
            pk.append(l[i])
            if len(pk)==3: break
    return pk

def f_growth(tag, m_phi, T_phi, Nur, z_lo, z_hi):
    a=run(f"{tag}_lo", m_phi,T_phi,Nur, z=z_lo)["pk"]
    b=run(f"{tag}_hi", m_phi,T_phi,Nur, z=z_hi)["pk"]
    s_lo=sigma8(*a); s_hi=sigma8(*b)
    a_lo=1/(1+z_lo); a_hi=1/(1+z_hi)
    f=(np.log(s_lo)-np.log(s_hi))/(np.log(a_lo)-np.log(a_hi))
    return f, s_lo

print("="*78)
print("  FILTRO CMB + fsigma8 — particula phi-DM ~40 eV (3 vias, dos sectores)")
print(f"  multiplicador (pi+phi)^4+(phi+2pi)^2 = {MULT:.4f}")
print("="*78)

# --- LCDM de referencia (mismo background fld, un solo sector frio) ---
ref=run("ref_lcdm", None,None,None, want_cmb=True, z=0.0)
l_r,TT_r=ref["cl"]; pk_r=peaks(l_r,TT_r)
print(f"  REF (1 sector frio, Om_cdm=0.271): picos TT = {pk_r}")
print("-"*78)

for nom,R in vias.items():
    Smnu=R*0.960318
    m_phi=Smnu*MULT
    T_phi=T_NU*(om_ncdm_h2*93.14/m_phi)**(1/3)
    Nur=max(0.01,3.046-(T_phi/T_NU)**4)
    o=run(nom, m_phi,T_phi,Nur, want_cmb=True, z=0.0)
    if "cl" not in o: continue
    l,TT=o["cl"]; pk=peaks(l,TT)
    # chi2 RMS relativo vs ref en el rango con datos (l 30-2500)
    msk=(l>=30)&(l<=2500)
    rms=np.sqrt(np.mean(((np.interp(l[msk],l,TT)-np.interp(l[msk],l_r,TT_r))/
                          np.interp(l[msk],l_r,TT_r))**2))*100
    s8=sigma8(*o["pk"]); S8=s8*np.sqrt(Om_m/0.3)
    f0,_=f_growth(nom,m_phi,T_phi,Nur,0.0,0.1)
    f05,s05=f_growth(nom,m_phi,T_phi,Nur,0.5,0.6)
    fs8_05=f05*s05
    print(f"  {nom}:  m_phi={m_phi:.2f} eV  T_phi/T_nu={T_phi/T_NU:.4f}")
    print(f"     CMB picos TT      = {pk}   (ref {pk_r})")
    print(f"     RMS(TT) vs ref    = {rms:.2f}%")
    print(f"     sigma8 / S8       = {s8:.4f} / {S8:.4f}")
    print(f"     fsigma8(z=0.5)    = {fs8_05:.4f}   [BOSS ref ~0.46]")
    print("-"*78)
print("="*78)
