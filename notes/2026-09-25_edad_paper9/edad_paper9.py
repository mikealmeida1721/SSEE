#!/usr/bin/env python3
"""Edad del universo — Paper 9 (nota 2026-09-25_edad_paper9).

Reproduce el calculo del Paper 9 y muestra que cambia con Omega_m=0.308881.
Uso: python3 edad_paper9.py   (requiere numpy, scipy)
"""
import numpy as np
from scipy.integrate import quad

W0, WA = -0.840, -0.670  # (w0, wa) SSEE


def E(z, Om, lcdm=False):
    Ode = 1.0 - Om
    if lcdm:
        fde = 1.0
    else:
        fde = (1 + z) ** (3 * (1 + W0 + WA)) * np.exp(-3 * WA * z / (1 + z))
    return np.sqrt(Om * (1 + z) ** 3 + Ode * fde)


def age_gyr(Om, H0, lcdm=False):
    I, _ = quad(lambda z: 1.0 / ((1 + z) * E(z, Om, lcdm=lcdm)),
                0, np.inf, limit=300)
    return 977.8 / H0 * I  # 1/H0 en Gyr


casos = [
    ("Paper9 actual (s_m en E(z))", 0.160050, 73.04, False),
    ("Corregida",                   0.308881, 73.04, False),
    ("Corregida, H_global SSEE",    0.308881, 67.962, False),
    ("LCDM, mismo H0",              0.308881, 73.04, True),
    ("LCDM Planck",                 0.308881, 67.4,  True),
]
for nombre, Om, H0, lcdm in casos:
    print(f"{nombre:28s} Om={Om:.6f}  H0={H0:6.3f}  ->  t0 = {age_gyr(Om, H0, lcdm):.2f} Gyr")
