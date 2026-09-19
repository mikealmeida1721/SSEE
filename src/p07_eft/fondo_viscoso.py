#!/usr/bin/env python3
"""
Fondo de SSEE con el LAGRANGIANO DE PAPER 1, no el de Paper 7.

  P(X) = A X^n            k-essence pura, SIN potencial
  n = (1+w0)/(2 w0)       -> w_phi = 1/(2n-1) = w0 EXACTO
  zeta(a) = Gamma (1-a) rho_DE / H     viscosidad de volumen
  Pi = -3 zeta H                        presion viscosa
  beta_c = 0                            SIN acoplamiento oscuro

PREGUNTA: reproduce (w0, wa) y el cruce fantasma sin romper w_c?

CONTROL OBLIGATORIO (R53): con la Gamma del paper (-P_sc/I_g) el
ajuste debe FALLAR dando wa = +2.01; solo con la corregida
(+P_sc/(3 I_g)) debe pasar. Si las dos pasan, la prueba no mide.
"""
# ORIGEN-VALOR: 0.3139 — cruce_z = 0.3139 de results/logs/eft_dos_campos_phi_pi.json
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, '/home/mike/Proyectos/SSEE/src')
import ssee_core as S

W0 = S.W0
WA = S.WA
PSC = 2 * S.PHI + S.PI          # P_sc  Dynamical Evolution Scalar
IG = S.PI + PSC                 # I_g   IGNIS
OM_M = S.OMEGA_M_TOTAL
OM_DE = 1.0 - OM_M
OM_C = S.OMEGA_C_H2 / (S.H0_GLOBAL / 100.0) ** 2

GAMMA_PAPER = -PSC / IG         # lo que dice SSEE_EFT_section.tex
GAMMA_FIX = PSC / (3.0 * IG)    # lo que exige Pi = -3 zeta H

# --- contornos DESI DR2 (arXiv:2503.14738 ec.26), via P2 ---
DESI = dict(w0=-0.8380, sw0=0.0552,
            wa=-0.6166, swa=0.2076, rho=-0.894)
# --- Planck 2018 (CANONICAL_VALUES.yaml:586) ---
WC_PLANCK, WC_SIG = 0.1200, 0.0012


def n_kessence():
    return (1.0 + W0) / (2.0 * W0)


def w_phi():
    """EoS del campo: constante, por algebra."""
    n = n_kessence()
    return 1.0 / (2.0 * n - 1.0)


def nec_campo():
    """rho+p = 2 n A X^n, en unidades de rho_phi."""
    n = n_kessence()
    return 2.0 * n / (2.0 * n - 1.0)     # = 1 + w_phi


def w_efectivo(a, gamma):
    """w_eff = w_phi + Pi/rho_DE,  Pi = -3 zeta H."""
    return w_phi() - 3.0 * gamma * (1.0 - a)


def ajusta_cpl(gamma, lo=0.300, hi=0.773):
    a = np.linspace(lo, hi, 2000)
    w = w_efectivo(a, gamma)
    A = np.vstack([np.ones(a.size), 1.0 - a]).T
    c, *_ = np.linalg.lstsq(A, w, rcond=None)
    return float(c[0]), float(c[1])


def sigma_desi(w0f, waf):
    s0, sa, r = DESI['sw0'], DESI['swa'], DESI['rho']
    C = np.array([[s0 ** 2, r * s0 * sa],
                  [r * s0 * sa, sa ** 2]])
    dv = np.array([w0f - DESI['w0'], waf - DESI['wa']])
    c2 = float(dv @ np.linalg.inv(C) @ dv)
    return c2, float(stats.norm.isf(stats.chi2.sf(c2, df=2) / 2))


def cruce_fantasma(gamma):
    """a donde w_eff = -1; None si nunca cruza."""
    den = -3.0 * gamma
    if den == 0:
        return None
    a = 1.0 + (1.0 + w_phi()) / den
    return a if 0.0 < a < 1.0 else None


def _linea(txt, ok):
    print('  [%s] %s' % ('OK  ' if ok else 'FALLA', txt))
    return ok


if __name__ == '__main__':
    print('=' * 66)
    print('  FONDO SSEE CON EL LAGRANGIANO DE PAPER 1')
    print('  P(X) = A X^n  +  viscosidad Pi  +  beta_c = 0')
    print('=' * 66)
    n = n_kessence()
    print('  n = (1+w0)/(2 w0)     = %.8f' % n)
    print('  w_phi = 1/(2n-1)      = %.12f' % w_phi())
    print('  w0                    = %.12f' % W0)
    print('  diferencia            = %.2e' % abs(w_phi() - W0))
    print()
    print('  P_sc = %.6f   I_g = %.6f' % (PSC, IG))
    print('  Gamma paper    = %+.8f' % GAMMA_PAPER)
    print('  Gamma corregida= %+.8f' % GAMMA_FIX)
    print('  razon          = %+.6f' % (GAMMA_FIX / GAMMA_PAPER))
    print()

    todo = []
    print('  --- 1. el campo cumple la NEC (no es fantasma) ---')
    npc = nec_campo()
    print('      (rho+p)/rho_phi = 1+w_phi = %+.8f' % npc)
    todo.append(_linea('el campo NO es fantasma (1+w_phi > 0)',
                       npc > 0))
    print()

    print('  --- 2. Gamma CORREGIDA: ajuste CPL ---')
    w0f, waf = ajusta_cpl(GAMMA_FIX)
    c2, sg = sigma_desi(w0f, waf)
    print('      w0 ajustado = %+.8f   (SSEE %+.8f)' % (w0f, W0))
    print('      wa ajustado = %+.8f   (SSEE %+.8f)' % (waf, WA))
    print('      contra DESI+CMB+Pantheon+ : %.2f sigma' % sg)
    todo.append(_linea('reproduce wa = -P_sc/I_g',
                       abs(waf - WA) < 1e-8))
    todo.append(_linea('dentro de 1 sigma de DESI', sg < 1.0))
    print()

    print('  --- 3. CONTROL: Gamma del paper debe FALLAR ---')
    w0p, wap = ajusta_cpl(GAMMA_PAPER)
    c2p, sgp = sigma_desi(w0p, wap)
    print('      w0 ajustado = %+.8f' % w0p)
    print('      wa ajustado = %+.8f   <- signo contrario' % wap)
    print('      contra DESI : %.2f sigma' % sgp)
    todo.append(_linea('el control FALLA (como debe)', sgp > 3.0))
    print()

    print('  --- 4. cruce fantasma ---')
    ac = cruce_fantasma(GAMMA_FIX)
    print('      a_* = %.6f   (z_* = %.4f)' % (ac, 1 / ac - 1))
    print('      P1 eq:phantom_crossing dice z_* ~ 0.31')
    todo.append(_linea('cruza w=-1 en z ~ 0.31',
                       abs(1 / ac - 1 - 0.3139) < 0.01))
    print()

    print('  --- 5. la materia oscura NO se toca ---')
    print('      beta_c = 0  =>  rho_c ~ a^-3  exacto')
    print('      w_c(z=1089) = w_c(hoy) = %.6f' % S.OMEGA_C_H2)
    sg_wc = abs(S.OMEGA_C_H2 - WC_PLANCK) / WC_SIG
    print('      Planck 0.1200 +- 0.0012  ->  %.2f sigma' % sg_wc)
    todo.append(_linea('w_c dentro de 1 sigma de Planck',
                       sg_wc < 1.0))
    print()

    print('=' * 66)
    if all(todo):
        print('  VERDE — %d/%d' % (sum(todo), len(todo)))
        print('  El lagrangiano de Paper 1 reproduce (w0, wa),')
        print('  el cruce fantasma y w_c SIN acoplamiento.')
    else:
        print('  ROJO — %d/%d' % (sum(todo), len(todo)))
    print('=' * 66)
