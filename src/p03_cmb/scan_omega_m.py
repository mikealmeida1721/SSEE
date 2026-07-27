"""
Scan chi2_CMB vs Omega_m a H=67.962 FIJO (todos los demas ingredientes SSEE).
Compara factor pi/phi (0.31076) vs factor (phi+0.1*pi) (0.30915) y busca el minimo.
"""
import numpy as np, math
import ssee_paper3_cobaya_unified as P3

# 2026-07-26: mnu estaba HARDCODEADO en 0.069 —el residuo de C_ν=94.07, retirado
# el 2026-07-25— mientras el resto de la suite ya usaba 0.06849. No era un log
# rancio: era el FUENTE rancio, y produce el barrido que justifica el ancla H₀.
# Ahora se toma del núcleo, para que no pueda volver a divergir.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), ".."))
from ssee_core import SUM_MNU_EV as _MNU                      # noqa: E402
H0=67.962; ombh2=0.02242; mnu=_MNU
phi=(1+5**0.5)/2; pi=math.pi
f_piphi = pi/phi                 # 1.94161 -> 0.31076
f_alma  = phi + 0.1*pi           # 1.93216 -> 0.30915

print(f"factor pi/phi      = {f_piphi:.5f} -> Om = {0.160*f_piphi:.5f}")
print(f"factor phi+0.1*pi  = {f_alma:.5f} -> Om = {0.160*f_alma:.5f}\n")

grid = [0.3000, 0.3050, 0.3070, 0.30915, 0.31076, 0.3130, 0.3160, 0.3200]
print(f"{'Om_cmb':>8} {'omch2':>9} {'omega_m':>9} {'chi2':>10}")
print("-"*40)
best=(None,1e9)
for Om in grid:
    omch2 = Om*(H0/100)**2 - ombh2
    om_m  = Om*(H0/100)**2
    c = P3.evaluate_model(H0, ombh2, omch2, P3.w0_ssee, P3.wa_ssee,
                          P3.As_ssee, P3.ns_ssee, P3.tau_ssee, mnu=mnu, quiet=True)
    tag=""
    if abs(Om-0.30915)<1e-4: tag=" <- phi+0.1pi"
    if abs(Om-0.31076)<1e-4: tag=" <- pi/phi (actual)"
    print(f"{Om:8.5f} {omch2:9.5f} {om_m:9.5f} {c:10.3f}{tag}")
    if c<best[1]: best=(Om,c)
print(f"\nMIN CMB en Om = {best[0]:.5f}  (chi2 = {best[1]:.3f})")
print(f"LCDM ref = 1003.76")
