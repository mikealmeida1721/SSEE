#!/usr/bin/env python3
"""
ssee_mira_disformal_test.py — ¿MIRA es geometría o contenido?
==============================================================
El factor MIRA lleva Ω_m de 0.160 a 0.320 en el sector CMB. Pregunta:
¿esa diferencia es un efecto puramente geométrico (que un acoplamiento
disformal del fotón PODRÍA reproducir) o física de fuente pre-recombinación
(que NO podría, porque actúa sobre el fotón ya emitido)?

DISCRIMINANTE: las razones de altura de picos P2/P1 y P3/P1 son INVARIANTES
bajo un reescalado puro del eje ℓ — y todo efecto de geometría/distancia
(incluido el disformal) actúa como un reescalado del eje ℓ.

  · Si P2/P1 y P3/P1 coinciden entre Ω_m=0.160 y 0.320  → diferencia
    geométrica → un mecanismo disformal PODRÍA reproducir MIRA. PUERTA ABIERTA.
  · Si difieren → la física de fuente (driving acústico, ISW temprano,
    fijada en la recombinación) cambió → ningún efecto de propagación lo
    reproduce. PUERTA CERRADA: MIRA exige contenido, no solo geometría.

CONTROL: se varía SOLO Ω_m (vía Ω_cdm). Ω_b h², H0, w0, wa, ns, As fijos,
así la razón par/impar (baryon-driven) y la amplitud no contaminan el test.
Las razones de pico son además invariantes bajo el factor e^-2τ de reionización.
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ssee_core import W0, WA, N_S, OMEGA_M_DYN, OMEGA_M_CMB_MIRA

import numpy as np

try:
    import camb
except ImportError:
    print("ERROR: camb no instalado — 'pip install camb' antes de correr.")
    sys.exit(1)

H0    = 67.0
OMBH2 = 0.02237          # fijo — baryones idénticos en ambos casos
NS    = N_S
AS    = 2.1e-9
h     = H0 / 100.0


def get_TT(Om):
    """Espectro TT Dℓ = ℓ(ℓ+1)Cℓ/2π [µK²] para una Ω_m dada."""
    omch2 = Om * h**2 - OMBH2
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=OMBH2, omch2=omch2, mnu=0.06, omk=0)
    pars.set_dark_energy(w=W0, wa=WA, dark_energy_model='ppf')
    pars.InitPower.set_params(ns=NS, As=AS)
    pars.set_for_lmax(2500, lens_potential_accuracy=0)
    results = camb.get_results(pars)
    powers = results.get_cmb_power_spectra(pars, CMB_unit='muK')
    return powers['total'][:, 0]


from scipy.signal import find_peaks


def acoustic_peaks(tt):
    """Los picos acústicos como máximos locales reales (no argmax en ventana
    fija: las ventanas fijas fallan cuando los picos se corren con Ω_m)."""
    seg = tt[100:1400]
    idx, _ = find_peaks(seg, distance=120, prominence=100)
    return sorted((100 + int(i), float(seg[i])) for i in idx)


print("=" * 68)
print("TEST MIRA — ¿geometría o contenido? (razones de altura de picos)")
print("=" * 68)
print(f"  Control: H0={H0}, Ωb h²={OMBH2}, ns={NS:.5f}, w0={W0:.4f}, wa={WA:.4f}")
print(f"  Se varía solo Ω_m.  Ω_m,dyn={OMEGA_M_DYN:.4f}  vs  Ω_m,CMB={OMEGA_M_CMB_MIRA:.4f}")
print("-" * 68)

results = {}
for Om, name in [(OMEGA_M_DYN, "0.160"), (OMEGA_M_CMB_MIRA, "0.320")]:
    tt = get_TT(Om)
    pk = acoustic_peaks(tt)
    if len(pk) < 3:
        print(f"\n  ERROR Ω_m={name}: solo {len(pk)} picos detectados: {pk}")
        sys.exit(1)
    p1, p2, p3 = pk[0], pk[1], pk[2]
    r21 = p2[1] / p1[1]
    r31 = p3[1] / p1[1]
    results[name] = dict(p1=p1, p2=p2, p3=p3, r21=r21, r31=r31)
    print(f"\n  Ω_m = {name}")
    print(f"    Pico 1:  ℓ={p1[0]:4d}   Dℓ={p1[1]:8.1f} µK²")
    print(f"    Pico 2:  ℓ={p2[0]:4d}   Dℓ={p2[1]:8.1f} µK²")
    print(f"    Pico 3:  ℓ={p3[0]:4d}   Dℓ={p3[1]:8.1f} µK²")
    print(f"    Razón P2/P1 = {r21:.4f}")
    print(f"    Razón P3/P1 = {r31:.4f}")

a, b = results["0.160"], results["0.320"]

# Posiciones: ¿un reescalado único alinea los 3 picos? (consistencia geométrica)
s1 = b["p1"][0] / a["p1"][0]
s2 = b["p2"][0] / a["p2"][0]
s3 = b["p3"][0] / a["p3"][0]

# Discriminante: cambio relativo de las razones de altura
d21 = abs(b["r21"] - a["r21"]) / a["r21"] * 100
d31 = abs(b["r31"] - a["r31"]) / a["r31"] * 100

print("\n" + "-" * 68)
print("  Reescalado de posiciones (geométrico si s1≈s2≈s3):")
print(f"    s1={s1:.4f}  s2={s2:.4f}  s3={s3:.4f}   dispersión={np.std([s1,s2,s3]):.4f}")
print("\n  DISCRIMINANTE — cambio de las razones de altura (invariantes geométricos):")
print(f"    Δ(P2/P1) = {d21:5.1f} %")
print(f"    Δ(P3/P1) = {d31:5.1f} %")

print("\n" + "=" * 68)
worst = max(d21, d31)
if worst < 5.0:
    print("  VEREDICTO: PUERTA ABIERTA")
    print("  Las razones de pico coinciden → la diferencia 0.160↔0.320 es")
    print("  geométrica. Un mecanismo disformal PODRÍA reproducir MIRA.")
elif worst > 15.0:
    print("  VEREDICTO: PUERTA CERRADA")
    print("  Las razones de pico difieren fuerte → la diferencia está en la")
    print("  física de fuente (pre-recombinación). Ningún efecto de propagación")
    print("  /disformal la reproduce. MIRA exige contenido real, no solo geometría.")
else:
    print("  VEREDICTO: AMBIGUO")
    print(f"  Cambio de razones {worst:.1f}% (zona 5–15%). Requiere análisis fino:")
    print("  un disformal podría cubrir parte, pero no es limpio.")
print("=" * 68)
