#!/usr/bin/env python3
"""
Robustness check (Mike, 2026-06-14): rerun REAL de CPL anclando el prior de H0
en varios valores, para responder con MUESTREO DIRECTO (sin el problema de ESS
de la reponderación) a qué H0 lleva el ajuste libre w0,wa hacia la predicción
algebraica de SSEE (-0.840, -0.670), y si ese H0 es SH0ES (73) o no.

Es la versión rigurosa del notebook lab_w0wa_degeneracy.ipynb (que solo reponderaba).
Corre CPL (5 params) con el prior de H0 centrado en cada ancla; el resto del
likelihood (DESI BAO + Planck Om/Obh2 + priors w0,wa) es idéntico a
ssee_paper2_mcmc.py.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import emcee
from ssee_core import KAL0, OMEGA_M_DYN as OM_DYN

C_KM = 2.998e5

# ── DESI DR2 (copiado de ssee_paper2_mcmc.py) ──
DESI_Z    = [0.295,0.510,0.510,0.706,0.706,0.930,0.930,1.317,1.317,1.491,1.491,2.330,2.330]
DESI_TYPE = ["DV_rd","DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd",
             "DM_rd","DH_rd","DM_rd","DH_rd","DM_rd","DH_rd"]
DESI_OBS  = np.array([7.93,13.62,20.08,16.85,19.50,21.71,17.88,27.79,13.82,30.21,13.23,39.71,8.52])
DESI_SIGMA= np.array([0.15,0.25,0.60,0.32,0.55,0.28,0.35,0.69,0.42,0.79,0.55,0.94,0.17])
RHO_PAIRS = {(1,2):-0.44,(3,4):-0.45,(5,6):-0.44,(7,8):-0.43,(9,10):-0.42,(11,12):-0.45}
DESI_COV  = np.diag(DESI_SIGMA**2)
for (i,j),rho in RHO_PAIRS.items():
    DESI_COV[i,j]=rho*DESI_SIGMA[i]*DESI_SIGMA[j]; DESI_COV[j,i]=DESI_COV[i,j]
DESI_COV_INV = np.linalg.inv(DESI_COV)

# ── Prior Planck Om/Obh2 (la parte de H0 se ANCLA por escenario) ──
PLANCK_OM, PLANCK_OBH2 = (0.3153,0.0073), (0.02237,0.00015)
SIGMA_H0 = 0.54

def f_de_cpl(z,w0,wa):
    a=1.0/(1.0+z); return (1+z)**(3*(1+w0+wa))*np.exp(-3*wa*(1-a))
def E_cpl(z,Om,w0,wa):
    return np.sqrt(Om*(1+z)**3+(1-Om)*f_de_cpl(z,w0,wa))
def DC(zmax,Om,w0,wa,n=300):
    zz=np.linspace(0,zmax,n); return np.trapezoid(1.0/E_cpl(zz,Om,w0,wa),zz)
def rd(ob_h2,om_h2):
    return 147.27*(om_h2/0.1432)**(-0.255)*(ob_h2/0.02237)**(-0.134)
def predict(H0,Om,w0,wa,rdv):
    out=[]
    for z,q in zip(DESI_Z,DESI_TYPE):
        dm=(C_KM/H0)*DC(z,Om,w0,wa); dh=C_KM/(H0*E_cpl(z,Om,w0,wa))
        out.append(dm/rdv if q=="DM_rd" else dh/rdv if q=="DH_rd" else (z*dm**2*dh)**(1/3)/rdv)
    return np.array(out)

def make_lpost(H0_anchor):
    """log-posterior CPL con el prior de H0 centrado en H0_anchor."""
    def lpost(th):
        H0,Om,w0,wa,ob=th
        if not(40<H0<100 and 0.15<Om<0.55 and -2.5<w0<0.5 and -3.0<wa<2.0 and 0.015<ob<0.030):
            return -np.inf
        om_h2=Om*(H0/100)**2
        lp =-0.5*((H0-H0_anchor)/SIGMA_H0)**2                       # ancla H0
        lp+=-0.5*((Om-PLANCK_OM[0])/PLANCK_OM[1])**2                # Planck Om
        lp+=-0.5*((ob-PLANCK_OBH2[0])/PLANCK_OBH2[1])**2            # Planck Obh2
        lp+=-0.5*((w0+1.0)/0.5)**2-0.5*(wa/1.0)**2                  # priors débiles w0,wa
        r=predict(H0,Om,w0,wa,rd(ob,om_h2))-DESI_OBS
        return lp-0.5*(r@DESI_COV_INV@r)
    return lpost

def run(anchor,nw=80,ns=4000,nb=1000):
    rng=np.random.default_rng(42)
    th0=np.array([anchor,0.315,-0.8,-0.6,0.02237])
    sc =np.array([0.4,0.01,0.1,0.25,0.0003])
    pos=th0+rng.standard_normal((nw,5))*sc
    s=emcee.EnsembleSampler(nw,5,make_lpost(anchor))
    pos,_,_=s.run_mcmc(pos,nb,progress=False); s.reset()
    s.run_mcmc(pos,ns,progress=False)
    ch=s.get_chain(flat=True)
    med=np.median(ch,axis=0); acc=np.mean(s.acceptance_fraction)
    return med,acc

if __name__=="__main__":
    SSEE=np.array([-0.840,-0.670])
    print(f"{'ancla H0':>10} | {'w0':>7} {'wa':>7} {'H0post':>7} {'Om':>6} | dist a SSEE | acc")
    print("-"*72)
    for name,a in [("MIRA",67.037),("Planck",67.36),("H_alg",67.962),("H0*~68.8",68.8),("SH0ES",73.0)]:
        med,acc=run(a)
        d=np.hypot(med[2]-SSEE[0],med[3]-SSEE[1])
        print(f"{name:>10} | {med[2]:+.3f} {med[3]:+.3f} {med[0]:7.2f} {med[1]:6.3f} | "
              f"{d:.3f} (eucl) | {acc:.2f}")
    print(f"{'SSEE alg':>10} | {SSEE[0]:+.3f} {SSEE[1]:+.3f}     —      —   |  0 (objetivo)")
