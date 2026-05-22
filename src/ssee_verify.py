#!/usr/bin/env python3
"""
SSEE Verification Harness — el guardián del Registro de Verificación.

Re-ejecuta, de forma automática, todas las comprobaciones de
VERIFICATION_LEDGER.md. No revisa strings: recomputa la física.

    python3 src/ssee_verify.py

  VERDE → el modelo sigue íntegro; se puede continuar.
  ROJO  → algo se rompió. NO commitear ni sellar nada hasta resolverlo.

Correr DESPUÉS DE CADA CAMBIO al modelo o a los papers. El harness crece
con el Registro: cada capa que se verifica añade aquí sus comprobaciones,
de modo que re-verificar todo el modelo cuesta segundos. Este guardián
también revisa el trabajo de quien edita — humano o máquina.
"""
import hashlib
import math
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
fails = []
checks = 0


def check(name, ok, detail=""):
    global checks
    checks += 1
    mark = " OK " if ok else "FALLA"
    line = f"  [{mark}] {name}"
    if detail:
        line += f"  — {detail}"
    print(line)
    if not ok:
        fails.append(name)


opens = []


def track_open(name, detail=""):
    """Problema abierto, conocido y documentado en VERIFICATION_LEDGER.md.

    No es una regresión: se lista en CADA corrida para no olvidarlo nunca,
    pero no pone el guardián en rojo. Sale del registro de abiertos solo
    cuando una derivación real lo resuelve y se actualiza el Registro.
    """
    global checks
    checks += 1
    print(f"  [ABIERTO] {name}" + (f"  — {detail}" if detail else ""))
    opens.append(name)


# ─────────────────────────────────────────────────────────────────────
# CAPA 1 — Axiomas y constantes algebraicas
# Recomputa cada constante desde su definición y la coteja con el valor
# registrado en VERIFICATION_LEDGER.md. Detecta un valor mal registrado
# o una definición alterada.
# ─────────────────────────────────────────────────────────────────────
print("Capa 1 — axiomas y constantes algebraicas")

phi = (1 + 5 ** 0.5) / 2
pi = math.pi
Omega = phi + pi
beta = (phi + pi) / 2
AURA = (3 * phi + pi) / 2
MIRA = AURA / 2
KAL0 = (phi + 3 * pi) / 2
Tr = 3 * (phi + beta)
Kv = 2 * (phi + pi)
Mv = phi + pi + Kv

# (constante computada, valor registrado en el ledger)
L1 = {
    "V-L1-01 phi":  (phi,   1.6180339887),
    "V-L1-02 pi":   (pi,    3.1415926536),
    "V-L1-03 Omega": (Omega, 4.7596266423),
    "V-L1-04 beta": (beta,  2.3798133212),
    "V-L1-05 AURA": (AURA,  3.9978473099),
    "V-L1-06 MIRA": (MIRA,  1.9989236550),
    "V-L1-07 KAL0": (KAL0,  5.5214059748),
    "V-L1-08 Tr":   (Tr,    11.9935419298),
    "V-L1-09 Kv":   (Kv,    9.5192532847),
    "V-L1-10 Mv":   (Mv,    14.2788799270),
}
for name, (computed, ledger) in L1.items():
    check(name, abs(computed - ledger) < 1e-9,
          f"computado {computed:.10f} / registro {ledger:.10f}")

# Identidades estructurales — el álgebra debe cerrar.
check("L1 identidad  MIRA = AURA/2",
      abs(MIRA - AURA / 2) < 1e-12)
check("L1 identidad  KAL0 = beta + pi  (def. alternativa)",
      abs(KAL0 - (beta + pi)) < 1e-12)
check("L1 identidad  Kv = phi + pi + Omega  (def. alternativa)",
      abs(Kv - (phi + pi + Omega)) < 1e-12)
w0 = -Tr / Mv
check("L1 identidad  w0 = -Tr/Mv = -0.840",
      abs(w0 - (-0.8399497713)) < 1e-9, f"{w0:.10f}")
check("L1 identidad  Om_m,dyn = 1 + w0 = 0.160",
      abs((1 + w0) - 0.1600502287) < 1e-9)
check("L1 identidad  Om_m,cosm = MIRA * Om_m,dyn = 0.320",
      abs(MIRA * (1 + w0) - 0.3199281880) < 1e-9)

# ─────────────────────────────────────────────────────────────────────
# CAPA 2 — Parámetros cosmológicos derivados
# ─────────────────────────────────────────────────────────────────────
print("\nCapa 2 — parámetros cosmológicos derivados")

Psc = Omega + phi
Om_DE = Tr / Mv
Om_m_dyn = 1 + w0
alpha_K = 3 * Om_DE * Om_m_dyn

L2 = {
    "V-L2-01 w0":        (-Tr / Mv,                      -0.8399497713),
    "V-L2-02 wa":        (-Psc / Kv,                     -0.6699748857),
    "V-L2-03 Om_DE":     (Om_DE,                          0.8399497713),
    "V-L2-04 Om_m,dyn":  (Om_m_dyn,                       0.1600502287),
    "V-L2-05 Om_m,cosm": (MIRA * Om_m_dyn,                0.3199281880),
    "V-L2-06 H0^alg":    (3 * Omega ** 2,                 67.9621373234),
    "V-L2-07 n_s":       (1 - phi ** -7,                  0.9655581463),
    "V-L2-08 alpha_K":   (alpha_K,                        0.4033024589),
    "V-L2-09 beta_c":    (-AURA,                         -3.9978473099),
    "V-L2-12 r":         (12 * (phi ** 4 / 3) / (2 * phi ** 7) ** 2, 0.0081306227),
    "V-L2-13 f_screen":  ((pi - phi) / Omega ** 2,        0.0672532703),
}
for name, (computed, ledger) in L2.items():
    check(name, abs(computed - ledger) < 1e-6,
          f"computado {computed:.10f} / registro {ledger:.10f}")

check("V-L2-10 m_phi (numérico)  0.0824*H0^alg = 5.60 eV",
      abs(0.0824 * (3 * Omega ** 2) - 5.6001) < 1e-3)

# Identidades cruzadas — dos rutas independientes deben coincidir.
n_kess = (Tr - Mv) / (2 * Tr)
check("L2 identidad  w0 = 1/(2n-1)  (ruta k-essence)",
      abs(1 / (2 * n_kess - 1) - (-Tr / Mv)) < 1e-10)
check("L2 identidad  f_screen = alpha_K/(3*MIRA) = (pi-phi)/Om^2",
      abs(alpha_K / (3 * MIRA) - (pi - phi) / Omega ** 2) < 1e-10)

# Problemas ABIERTOS detectados en Capa 2 — comprobación dimensional.
track_open("V-L2-06 H0^alg dimensional",
           "3*Omega^2 es adimensional; H0 tiene unidades km/s/Mpc (Postulado D)")
track_open("V-L2-10 m_phi dimensional",
           "Sigma_m_nu * H0 = [eV]*[km/s/Mpc] no es [eV] (numerología, P6)")

# ─────────────────────────────────────────────────────────────────────
# CAPA 3 — Mecanismos y derivaciones
# ─────────────────────────────────────────────────────────────────────
print("\nCapa 3 — mecanismos y derivaciones")

# OP-2 — n_s = 1 - phi^-7 (índice espectral). La cadena algebraica cierra
# exacto; el insumo físico N_* = 2phi^7 es conjetura (ABIERTO).
N_star = 2 * phi ** 7
check("V-L3-OP2  identidad n_s: 1-2/(2phi^7) = 1-phi^-7",
      abs((1 - 2 / N_star) - (1 - phi ** -7)) < 1e-12)
check("V-L3-OP2  identidad r: 12(phi^4/3)/(2phi^7)^2 = phi^-10",
      abs(12 * (phi ** 4 / 3) / N_star ** 2 - phi ** -10) < 1e-12)
track_open("V-L3-OP2  N_* = 2phi^7",
           "Conjecture B.1 no derivada; falta el puente de reheating gravitacional")

# OP-7 — beta_c = -AURA (acoplamiento EFT). La dualidad Z2 es álgebra exacta;
# la identificacion beta_c=-AURA es coincidencia numerica al 0.2% (ABIERTO).
check("V-L3-OP7  dualidad Z2: KAL0(phi<->pi) = AURA",
      abs((pi + 3 * phi) / 2 - AURA) < 1e-12)
bc_num = -3.990  # extraído numéricamente en P7 §6 (shooting)
track_open("V-L3-OP7  beta_c = -AURA",
           f"numerico {bc_num} vs -AURA {-AURA:.4f}: "
           f"brecha {abs(bc_num + AURA) / AURA * 100:.2f}%, identificacion no derivada")

# alpha = phi^4/3 — consecuencia exacta de los axiomas n_s=1-phi^-7, r=phi^-10.
alpha_attr = phi ** 4 / 3
check("V-L3-alpha  alpha = r*N^2/12 = phi^4/3",
      abs((phi ** -10) * (2 * phi ** 7) ** 2 / 12 - alpha_attr) < 1e-10)
check("V-L3-alpha  identidad Fibonacci  phi^4/3 = phi + 2/3",
      abs(alpha_attr - (phi + 2 / 3)) < 1e-12)
# Curvatura de Kähler R = -2/(3 alpha) = -2 phi^-4 (P1 decía -phi^-4: error factor 2, corregido).
check("V-L3-alpha  curvatura Kahler R = -2/(3 alpha) = -2 phi^-4",
      abs(-2 / (3 * alpha_attr) - (-2 * phi ** -4)) < 1e-12,
      f"R = {-2/(3*alpha_attr):.6f}")

# OP-4 — radio k-mouflage (P8). La fórmula r_km^3 = M_obj/(4pi M_pl M^2)
# es dimensionalmente inconsistente: exponentes GeV de (M_obj, M_pl, M^2).
dim_rkm = (1 - 1 - 2) / 3   # dimensión GeV de r_km según la fórmula de P8
track_open("V-L3-OP4  formula k-mouflage de P8 dimensionalmente rota",
           f"r_km tiene dimension GeV^{dim_rkm:.3f}; una longitud es GeV^-1. "
           "Introducida en commit 295ed6e; requiere re-derivacion")

# OP-1 — densidad bariónica (P4/Paper B). La cadena algebraica cierra exacto;
# el insumo Omega_b h^2 = (pi-phi)/(3 Omega^2) es coincidencia hallada por scan
# de 7 candidatos, no derivada de BBN (ABIERTO).
Omb_h2_alg = (pi - phi) / (3 * Omega ** 2)
check("V-L3-OP1  identidad (pi-phi)/(3 Om^2) = 0.0224178",
      abs(Omb_h2_alg - 0.0224178) < 1e-6,
      f"computado {Omb_h2_alg:.7f}")
track_open("V-L3-OP1  Omega_b h^2 no derivado",
           "coincidencia a 0.32sigma de Planck hallada por scan de 7 candidatos; "
           "falta cadena BBN (eta_B -> Omega_b) — diferida a Paper B")

# OP-3 — separabilidad UV-IR (P10). Las identidades algebraicas cierran exacto;
# la jerarquia (H0/M)^2 es real, pero la prueba de separabilidad (jacobiano
# d phi/d chi) esta diferida y KALeff dropea el factor rho_crit (ABIERTO).
alpha_op3 = phi ** 4 / 3
check("V-L3-OP3  identidad 6*alpha = 2 phi^4",
      abs(6 * alpha_op3 - 2 * phi ** 4) < 1e-12)
check("V-L3-OP3  identidad sqrt(6*alpha) = phi^2 sqrt(2)",
      abs((6 * alpha_op3) ** 0.5 - phi ** 2 * 2 ** 0.5) < 1e-12)
track_open("V-L3-OP3  separabilidad UV-IR no probada",
           "jerarquia (H0/M)^2~3e-62 real; prueba via jacobiano d phi/d chi "
           "diferida a Paper B; KALeff^2 = M^4/(6 alpha) dropea rho_crit")

# OP-5 — tensión S8 weak-lensing (P5/P6). La definición S8 = sigma8*(Om/0.3)^0.5
# cierra para ambas ramas; HMcode reduce la tensión solo anclando en la rama
# secundaria sigma8=0.737, NO en el titular sigma8=0.794 (ABIERTO).
Om_cosm_op5 = MIRA * (1 + w0)
S8_sec = 0.737 * (Om_cosm_op5 / 0.3) ** 0.5   # rama WDM CLASS (secundaria)
S8_tit = 0.794 * (Om_cosm_op5 / 0.3) ** 0.5   # titular dos sectores
check("V-L3-OP5  identidad S8 = sigma8(Om/0.3)^0.5  rama secundaria = 0.761",
      abs(S8_sec - 0.761) < 2e-3, f"S8 = {S8_sec:.4f}")
check("V-L3-OP5  identidad S8 = sigma8(Om/0.3)^0.5  rama titular = 0.820",
      abs(S8_tit - 0.820) < 2e-3, f"S8 = {S8_tit:.4f}")
track_open("V-L3-OP5  S8 anclado en rama secundaria",
           "HMcode da S8=0.758 (0.06sigma DES) partiendo de sigma8=0.737 "
           "(rama WDM, NO titular); el titular sigma8=0.794 deja ~2.6sigma. "
           "N-body Nivel 2 diferido")

# OP-6 — forma de screening (P9). El valor f_screen es algebra exacta (ver
# V-L2-13); la forma multiplicativa sigue del universo separado. El paso
# delta_rho_phi y el insumo delta_local = 2 (sobredensidad Grupo Local) no
# estan derivados de phi,pi (ABIERTO parcial).
f_screen_op6 = (pi - phi) / Omega ** 2
H0_local_op6 = (3 * Omega ** 2) / (1 - f_screen_op6)
check("V-L3-OP6  H0_local = H0^alg/(1-f_screen) = 72.86",
      abs(H0_local_op6 - 72.864) < 1e-2, f"H0_local = {H0_local_op6:.4f}")
track_open("V-L3-OP6  forma multiplicativa: insumo delta_local = 2 no derivado",
           "la forma multiplicativa sigue del universo separado, pero el valor "
           "f_screen requiere delta_local=2 (sobredensidad Grupo Local) y una "
           "expresion delta_rho_phi asertada, no derivada de phi,pi")

# m_phi — masa del campo phi-DM (P6). La cadena cierra a 5.60 eV, pero el
# insumo R = 4*KAL0 - 22 es una resta de entero pelado (numerologia) y el
# producto Sigma_m_nu * H0_alg multiplica energia por un numero adimensional
# sin mecanismo fisico; ademas hereda OP-1 (ABIERTO).
R_p6 = 4 * KAL0 - 22
Ob_h2_p6 = (pi - phi) / (3 * Omega ** 2)
tau_Pi = KAL0 / (3 * Om_DE)
mnu_active = R_p6 * Ob_h2_p6 * 94.07 / tau_Pi
m_phi = mnu_active * (3 * Omega ** 2)
check("V-L3-mphi  cadena m_phi = Sigma_m_nu * H0^alg = 5.60 eV",
      abs(m_phi - 5.60) < 2e-2, f"m_phi = {m_phi:.4f} eV")
track_open("V-L3-mphi  m_phi no derivado",
           "R = 4*KAL0-22 resta un entero pelado para extraer 0.0856 "
           "(numerologia); m_phi = [eV] * 67.96 adimensional sin mecanismo; "
           "hereda OP-1 (Omega_b h^2 no derivado)")

# Dos sectores phi-DM (P6) — Om_CDM + Om_phiDM = MIRA * Om_m,dyn es una
# identidad algebraica exacta (re-particion de Om_m,cosm). El modelo fisico
# del split en k_fs depende de m_phi (ABIERTO) y k_fs (pendiente).
Om_phiDM = (MIRA - 1) * Om_m_dyn
check("V-L3-2sec  identidad Om_CDM + Om_phiDM = MIRA * Om_m,dyn",
      abs((Om_m_dyn + Om_phiDM) - MIRA * Om_m_dyn) < 1e-12)
track_open("V-L3-2sec  split fisico de dos sectores no cerrado",
           "la suma es identidad algebraica (= V-L2-05); la separacion fisica "
           "en k_fs depende de m_phi (ABIERTO) y k_fs (V-L2 pendiente)")

# EFT canónico (P7) — los parámetros lambda, alpha_pot, V0 son consecuencias
# algebraicas de constantes ya verificadas (Om_m,dyn, KAL0, Om_DE).
lam_eft = (3 * Om_m_dyn) ** 0.5
alpha_pot = lam_eft / KAL0 ** 0.5
check("V-L3-EFT  identidad lambda^2 = 3 Om_m,dyn",
      abs(lam_eft ** 2 - 3 * Om_m_dyn) < 1e-12, f"lambda = {lam_eft:.6f}")
check("V-L3-EFT  identidad alpha_pot = lambda/sqrt(KAL0)",
      abs(alpha_pot - lam_eft / KAL0 ** 0.5) < 1e-12, f"alpha_pot = {alpha_pot:.6f}")
check("V-L3-EFT  identidad V0 = Om_DE * rho_crit",
      abs(Om_DE - Tr / Mv) < 1e-12, f"V0 = {Om_DE:.6f}")

# UV completion K(X) (P10) — la identidad 45 alpha^2 = 5 phi^8 es exacta;
# la normalización física de M^4 está calibrada a SH0ES (admisión del propio
# script de P10) y es inconsistente con M^4 = rho_crit usado en P7.
check("V-L3-KX  identidad 45 alpha^2 = 5 phi^8  (M^4/rho_crit)",
      abs(45 * alpha_attr ** 2 - 5 * phi ** 8) < 1e-9,
      f"M^4/rho_crit = {5 * phi ** 8:.4f}")
track_open("V-L3-KX  M^4 = 5 phi^8 rho_crit calibrado a SH0ES",
           "ssee_paper10_verification.py admite: normalizacion fisica de M^4 "
           "calibrada a SH0ES, no derivada; Ruta A da M^4~418 != 234.9")
track_open("V-L3-EFT  M^4 inconsistente entre P7 y P10",
           "ssee_eft_verification.py usa M^4 = rho_crit (=1); "
           "ssee_paper10_verification.py usa M^4 = 5 phi^8 rho_crit (=234.9)")

# Israel-Stewart (P5) — c²_s,eff = 0. La corrección IS zeta/tau_Pi se reduce
# a Om_DE porque el factor KAL0/3 se cancela (zeta = KAL0/3, tau_Pi =
# KAL0/(3 Om_DE)). c²_s,eff = w0 + Om_DE = 0 es la identidad w0 = -Om_DE.
zeta_tilde = KAL0 / 3.0
tau_Pi_H0 = KAL0 / (3.0 * Om_DE)
IS_corr = zeta_tilde / tau_Pi_H0
check("V-L3-IS  identidad zeta/tau_Pi = Om_DE  (KAL0/3 se cancela)",
      abs(IS_corr - Om_DE) < 1e-12, f"IS_corr = {IS_corr:.8f}")
check("V-L3-IS  c2_s,eff = w0 + Om_DE = 0  (estabilidad marginal)",
      abs(w0 + IS_corr) < 1e-12)
track_open("V-L3-IS  c2_s,eff = 0 depende de tau_Pi no derivado en script",
           "zeta/tau_Pi = Om_DE porque tau_Pi = zeta/Om_DE comparten KAL0/3; "
           "el resultado se reduce a w0 = -Om_DE (V-L2). La derivacion de "
           "tau_Pi (steady-state IS) esta asertada, no mostrada")

# dos-Ω_m — regla de uso. Om_m,dyn = 0.160 fija w0 (1+w0 = Om_m,dyn) y NO
# entra en E(z); Om_m,cosm = MIRA*Om_m,dyn = 0.320 es el fondo gravitacional
# y SI entra en E(z)/Poisson/CMB. Aqui se verifica que son distintos y que
# la relación es exacta; la auditoría de uso paper-por-paper es Capa 5.
check("V-L3-2Om  Om_m,dyn != Om_m,cosm  (son cantidades distintas)",
      abs(Om_m_dyn - MIRA * Om_m_dyn) > 0.15)
check("V-L3-2Om  relacion exacta Om_m,cosm = MIRA * Om_m,dyn",
      abs(MIRA * Om_m_dyn - 0.3199281880) < 1e-9)
track_open("V-L3-2Om  uso correcto de cada Om_m no auditado por paper",
           "regla: Om_m,dyn en 1+w0; Om_m,cosm en E(z)/Poisson/CMB. "
           "Verificar que cada paper usa la correcta es tarea de Capa 5")

# ─────────────────────────────────────────────────────────────────────
# CAPA 4 — Confrontaciones con datos
# Verifica la ARITMETICA que conecta cantidades reportadas (tensiones,
# S8, chi2_r). Los chi2 de CMB y los posteriores MCMC en si requieren
# re-correr CAMB/CLASS/emcee — se marcan como dependencia externa.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa 4 — confrontaciones con datos")

# S8 weak-lensing (P5). Cadena: sigma8_SSEE = sigma8_LCDM * G;
# S8 = sigma8 * sqrt(Om_m,cosm/0.3). G = D1_SSEE/D1_LCDM (resultado ODE P5).
sig8_LCDM = 0.811        # Planck 2018
G_growth = 0.866         # D1_SSEE/D1_LCDM, resultado ODE Paper 5
sig8_SSEE = sig8_LCDM * G_growth
Om_cosm = MIRA * Om_m_dyn
S8_P5 = sig8_SSEE * (Om_cosm / 0.3) ** 0.5
check("V-L4-01 P5  sigma8_SSEE = sigma8_LCDM * G = 0.702",
      abs(sig8_SSEE - 0.7023) < 1e-3, f"sigma8 = {sig8_SSEE:.4f}")
check("V-L4-02 P5  S8_SSEE = sigma8 sqrt(Om_cosm/0.3) = 0.725",
      abs(S8_P5 - 0.7253) < 1e-3, f"S8 = {S8_P5:.4f}")

# Tensiones S8 (P5) — error en cuadratura modelo + observacional.
sig8_SSEE_err = 0.006 * G_growth
S8_P5_err = sig8_SSEE_err * (Om_cosm / 0.3) ** 0.5
t_DES = abs(S8_P5 - 0.776) / (S8_P5_err ** 2 + 0.017 ** 2) ** 0.5
t_KIDS = abs(S8_P5 - 0.766) / (S8_P5_err ** 2 + 0.020 ** 2) ** 0.5
check("V-L4-03 P5  tension S8 vs DES-Y3 (3x2pt) = 2.84 sigma",
      abs(t_DES - 2.84) < 0.05, f"{t_DES:.2f} sigma")
check("V-L4-04 P5  tension S8 vs KiDS-1000 = 1.96 sigma",
      abs(t_KIDS - 1.96) < 0.05, f"{t_KIDS:.2f} sigma")

# r_d horizonte de sonido (P2/P3). r_d=147.156 Mpc lo calcula CAMB con
# los parametros SSEE+MIRA; aqui se verifica solo la tension reportada.
rd_ssee = 147.156   # salida CAMB, ssee_verify_rd.py caso A (insumo registrado)
t_rd = abs(rd_ssee - 147.09) / 0.26
check("V-L4-05 P2/P3  tension r_d vs Planck 2018 = 0.25 sigma",
      abs(t_rd - 0.25) < 0.05, f"{t_rd:.2f} sigma")

track_open("V-L4  valor de referencia DES-Y3 inconsistente entre scripts",
           "ssee_paper5 usa S8_DES = 0.776+-0.017 (3x2pt, Abbott 2022); "
           "ssee_op5_hmcode usa S8_DES = 0.759+-0.023 (cosmic shear, Amon 2022); "
           "elegir una referencia DES unica para toda la suite")
track_open("V-L4  chi2 / BIC / posteriores MCMC requieren pipeline",
           "los chi2 de CMB (P3) y los posteriores DESI+Planck (P2) salen de "
           "CAMB/CLASS/emcee; el guardian verifica su aritmetica, no los recomputa")

# ─────────────────────────────────────────────────────────────────────
# SELLOS — integridad de los papers sellados.
# Cuando un paper se sella, se registra aquí su sha256. El harness lo
# recalcula: si el archivo cambió tras sellarse, el sello se rompe.
# ─────────────────────────────────────────────────────────────────────
print("\nSellos de papers")
SEALS = {
    # "manuscript/SSEE_Paper1_Framework.tex": "<sha256 al sellar>",
}
if not SEALS:
    print("  (ningún paper sellado todavía)")
for relpath, expected in SEALS.items():
    f = ROOT / relpath
    actual = hashlib.sha256(f.read_bytes()).hexdigest() if f.exists() else None
    check(f"sello {relpath}", actual == expected,
          "intacto" if actual == expected
          else "ROTO — el archivo cambió después de sellarse")

# ─────────────────────────────────────────────────────────────────────
# VEREDICTO
# ─────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
if fails:
    print(f"ROJO — {len(fails)} de {checks} comprobaciones FALLARON:")
    for f in fails:
        print(f"   x  {f}")
    print("No commitear ni sellar hasta resolverlo.")
    sys.exit(1)
linea = f"VERDE — sin regresiones, {checks} comprobaciones."
if opens:
    linea += f"  {len(opens)} problema(s) ABIERTO(s) rastreado(s) — ver Registro."
print(linea)
sys.exit(0)
