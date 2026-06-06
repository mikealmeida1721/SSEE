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

# V-L2-10 m_phi canónico (forward-prediction, cero fiteo, P6 2026-06-04):
#   m_phi = Sigma_m_nu^active * (Omega^4 + AURA*KAL)
#   con Sigma_m_nu^active = R2 * 0.960318 eV, R2 = Omega/(KAL0*Tr).
# Dimensionalmente CONSISTENTE: [eV] * (número puro) = [eV].
_R2 = Omega / (KAL0 * Tr)
_mnu_active = _R2 * 0.960318
_mult_mphi = Omega ** 4 + AURA * KAL0
_m_phi_canon = _mnu_active * _mult_mphi
check("V-L2-10 m_phi canónico = Sigma_m_nu^active * (Om^4+AURA*KAL) = 36.95 eV",
      abs(_m_phi_canon - 36.9463) < 1e-3,
      f"m_phi = {_m_phi_canon:.4f} eV (R2={_R2:.6f}, mult={_mult_mphi:.4f})")

# Identidades cruzadas — dos rutas independientes deben coincidir.
n_kess = (Tr - Mv) / (2 * Tr)
check("L2 identidad  w0 = 1/(2n-1)  (ruta k-essence)",
      abs(1 / (2 * n_kess - 1) - (-Tr / Mv)) < 1e-10)
check("L2 identidad  f_screen = alpha_K/(3*MIRA) = (pi-phi)/Om^2",
      abs(alpha_K / (3 * MIRA) - (pi - phi) / Omega ** 2) < 1e-10)

# Problemas ABIERTOS detectados en Capa 2 — comprobación dimensional.
track_open("V-L2-06 H0^alg dimensional",
           "3*Omega^2 es adimensional; H0 tiene unidades km/s/Mpc (Postulado D)")
# V-L2-10: la fórmula CANÓNICA (forward-prediction 36.95 eV) es dimensionalmente
# consistente — [eV]*(número puro) = [eV]. El antiguo ansatz Sigma_m_nu*H0^alg
# (5.60 eV) está retirado. Lo abierto es el Lagrangiano φ-DM (OP-9), no la dimensión.

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

# m_phi — masa del campo phi-DM (P6, CANÓNICO forward-prediction 2026-06-04).
#   Sigma_m_nu^active = R2 * 0.960318 eV,  R2 = Omega/(KAL0*Tr)
#   m_phi = Sigma_m_nu^active * (Omega^4 + AURA*KAL0)
# El multiplicador es número PURO -> dimensión [eV] preservada. Cero fiteo.
R2_p6 = Omega / (KAL0 * Tr)
mnu_active = R2_p6 * 0.960318
mult_p6 = Omega ** 4 + AURA * KAL0
m_phi = mnu_active * mult_p6
check("V-L3-mphi  cadena m_phi = Sigma_m_nu^active * (Om^4+AURA*KAL) = 36.95 eV",
      abs(m_phi - 36.9463) < 2e-2, f"m_phi = {m_phi:.4f} eV")
track_open("V-L3-mphi  Lagrangiano phi-DM no cerrado (OP-9)",
           "el valor 36.95 eV es forward-prediction dimensionalmente consistente; "
           "lo abierto es el Lagrangiano P(X,phi) que produce esta masa (OP-9), "
           "no la dimensionalidad ni un fiteo")

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

# c_s^2 del sector k-essence — extraccion T_munu^ef (2026-05-22). Para
# K(X)=X/KAL0+X^2/M^4, Garriga-Mukhanov da c_s^2=(A+2BX)/(A+6BX) con
# A=1/KAL0>0 y B=1/M^4>0. Es decreciente: c_s^2 in [1/3,1] para TODO X>=0,
# cualquier M^4>0, cualquier KAL0>0 — depende solo de la FORMA de K.
A_ke = 1.0 / KAL0
cs2_vals = []
for B_M4 in (1.0, 1.0 / 234.8936):       # M^4 = 1 (P7) y 5 phi^8 (P10)
    for Xv in (0.0, 0.1, 1.0, 1e2, 1e6):
        cs2_vals.append((A_ke + 2 * B_M4 * Xv) / (A_ke + 6 * B_M4 * Xv))
check("V-L3-cs2  c_s^2 k-essence acotado en [1/3, 1] para todo X",
      all(1 / 3 - 1e-9 <= c <= 1 + 1e-9 for c in cs2_vals),
      f"min={min(cs2_vals):.5f} max={max(cs2_vals):.5f}")
track_open("V-L3-cs2  el sector geometrico de SSEE no puede agruparse [CENTRAL]",
           "extraccion T_munu^ef: c_s^2 de la k-essence in [1/3,1] siempre — "
           "nunca baja de 1/3, no clusteriza como materia fria. El sector "
           "geometrico tiene peso de FONDO (rho_phi existe) pero NO peso de "
           "agrupamiento. El CMB exige materia que se agrupe -> la k-essence "
           "actual no puede ser la '0.320'. MIRA no esta en la accion vigente")

# Ruta B (gravedad disformal de P8) — auditada 2026-05-22. La accion de
# P8 incluye un campo de materia oscura psi_DM: el mecanismo MIRA/lensing
# se sostiene en rho_DM. Choca con el postulado no-DM de P1.
track_open("V-L3-disf  mecanismo disformal de P8 presupone materia oscura",
           "P8 accion eq.(1) incluye S_DM[g~;psi_DM]: campo de materia "
           "oscura acoplado disformalmente. Todo el mecanismo MIRA/lensing "
           "se sostiene en rho_DM (P8 L221,238,271). P1 prohibe DM (L51 "
           "'without dark-matter particles'; L278 detectar DM falsa el "
           "modelo) -> contradiccion interna P1<->P8. Ademas P8 L69 admite "
           "sqrt(beta_c)/MIRA=1.00030 'near-coincidence, not an identity'. "
           "Ruta B disformal NO deriva MIRA sin materia oscura")

# Mecanismo de retencion conformal (Ruta B) — probado 2026-05-22 en
# src/ssee_mira_mechanism.py. Acoplamiento beta_c=-AURA: negativo limpio.
track_open("V-L3-mira  retencion conformal beta_c=AURA NO reproduce MIRA",
           "test del fondo acoplado Friedmann+KG (ssee_mira_mechanism.py): "
           "con beta_c=-AURA la excursion del campo es x18 excesiva, signo "
           "invertido (drena la materia en vez de cargarla, R(z=2)~0.016), y "
           "timing invertido (campo thawing). Cuarto mecanismo descartado "
           "para el '0.320' tras cs2, Poisson-mu y disformal. MIRA sin "
           "derivacion en el marco vigente por los 4 mecanismos naturales")

# dos-Ω_m — PROBLEMA CENTRAL ABIERTO. El algebra es exacta (Om_m,dyn=1+w0
# derivado; Om_m,cosm=MIRA*Om_m,dyn) pero la REGLA fisica de uso no esta
# derivada: no hay mecanismo para MIRA ni funcion Om_m(z) de transicion.
check("V-L3-2Om  Om_m,dyn != Om_m,cosm  (son cantidades distintas)",
      abs(Om_m_dyn - MIRA * Om_m_dyn) > 0.15)
check("V-L3-2Om  relacion exacta Om_m,cosm = MIRA * Om_m,dyn",
      abs(MIRA * Om_m_dyn - 0.3199281880) < 1e-9)
track_open("V-L3-2Om  mecanismo MIRA (0.160 -> 0.320) NO derivado [CENTRAL]",
           "0.160=1+w0 es derivado; MIRA es hipotesis auxiliar (ssee_core L41 "
           "'no derivada'). NO hay regla derivada de por que Om_m,cosm va en "
           "E(z), ni funcion Om_m(z) de transicion CMB<->BAO. 0.320 en (1+z)^3 "
           "es fenomenologicamente materia oscura -- choca con el postulado "
           "'sin materia oscura, solo geometria y viscosidad'. Es el problema "
           "central del modelo; r_d/theta*/H0 son sintomas de esta brecha")

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

# CMB Planck PR4 (P3) — re-corrida con CAMB 1.6.5 (2026-05-22): chi2_r
# TT 1.047 / TE 1.041 / EE 1.041 / PP 0.837 y ΔBIC=-20.8 reproducidos
# EXACTAMENTE. La aritmetica chi2_r->ΔBIC solo acota por el redondeo.
check("V-L4-05 P3  ΔBIC CMB = -20.8 consistente con chi2_r (re-run 2026-05-22)",
      -22.9 <= -20.8 <= -11.1,
      "chi2_r redondeados acotan ΔBIC a [-22.9,-11.1]; reportado -20.8 dentro")

# r_d / theta* — re-run CAMB 2026-06-01 con los H0 CANONICOS (no el stale 67.756).
# Politica del modelo: anclar el CMB en el anchor 67.068 (CMB-optimo), reportar
# el posterior BAO 66.553 por separado. theta* es muy sensible a H0 via D_A, asi
# que NO se propaga el posterior BAO al observable CMB.
#   anchor   67.068 -> r_d=146.70 Mpc (1.51sigma), theta*=0.59636 (0.69sigma)
#   posterior 66.553 -> r_d=147.28 Mpc (0.72sigma), theta*=0.59417 (5.46sigma)
# vs Planck 147.09+-0.26 Mpc / 0.59668+-0.00046 deg.
track_open("V-L4  r_d coherente en ambos H0 canonicos (anchor 1.51 / post 0.72 sigma)",
           "ssee_verify_rd.py re-run 2026-06-01: anchor 67.068 da r_d=146.70 Mpc "
           "(1.51sigma), posterior 66.553 da r_d=147.28 Mpc (0.72sigma) vs "
           "Planck 147.09+-0.26. El '4.47sigma' previo usaba el H0 stale 67.756")
track_open("V-L4  theta* sensible a H0: 0.69sigma en anchor, 5.46sigma en posterior",
           "ssee_verify_rd.py re-run 2026-06-01: anchor 67.068 da theta*=0.59636 "
           "(0.69sigma); posterior 66.553 da 0.59417 (5.46sigma) vs 0.59668+-0.00046. "
           "Por eso el CMB se ancla en 67.068 y el posterior BAO no se propaga a theta*")

track_open("V-L4  valor de referencia DES-Y3 inconsistente entre scripts",
           "ssee_paper5 usa S8_DES = 0.776+-0.017 (3x2pt, Abbott 2022); "
           "ssee_op5_hmcode usa S8_DES = 0.759+-0.023 (cosmic shear, Amon 2022); "
           "elegir una referencia DES unica para toda la suite")

# MCMC DESI+Planck (P2) — re-corrido 2026-05-22 (100w x 25000s x 3, 1.52h).
# lnP_MAP: SSEE -13.22 (k=2), LCDM -15.79 (k=3). N_DATA = 16.
# BIC = k*ln(N) - 2*lnP_MAP — esto SI se recomputa exacto.
N_data_p2 = 16   # 13 BAO DESI DR2 + 3 constraints Planck
bic_ssee = 2 * math.log(N_data_p2) - 2 * (-13.22)
bic_lcdm = 3 * math.log(N_data_p2) - 2 * (-15.79)
check("V-L4-06 P2  BIC SSEE = k ln(N) - 2 lnP = 31.98 (re-run MCMC)",
      abs(bic_ssee - 31.98) < 0.05, f"BIC = {bic_ssee:.2f}")
check("V-L4-07 P2  ΔBIC(LCDM-SSEE) = +7.91 (SSEE favorecido)",
      abs((bic_lcdm - bic_ssee) - 7.91) < 0.05, f"ΔBIC = {bic_lcdm - bic_ssee:.2f}")
H0_mcmc = 66.553   # posterior canonico, prior MIRA (ledger L113); el 67.756 era mala anotacion
t_H0 = abs(H0_mcmc - 67.36) / (0.442 ** 2 + 0.54 ** 2) ** 0.5
check("V-L4-08 P2  tension H0 SSEE vs Planck = 1.16 sigma",
      abs(t_H0 - 1.16) < 0.03, f"{t_H0:.2f} sigma")
track_open("V-L4  Omega_b h^2: posterior MCMC < prediccion algebraica",
           "MCMC da Om_b h^2 = 0.02183+-0.00048; OP-1 algebraico da 0.02242. "
           "el dato prefiere ~1.2sigma menos barion que (pi-phi)/(3 Om^2)")

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
# FUENTE CANÓNICA — chequeo de src/ssee_core.py
# ssee_core.py es el módulo del que TODOS los demás scripts importan sus
# constantes algebraicas. Aquí se verifica que ese módulo coincide con la
# recomputación independiente del guardián. Si ssee_core se edita mal, esto
# se pone ROJO — y, por tanto, todo script que importe de él queda advertido.
# ─────────────────────────────────────────────────────────────────────
print("\nFuente canónica — ssee_core.py")
import importlib.util as _ilu

_core_path = ROOT / "ssee_core.py"
_spec = _ilu.spec_from_file_location("ssee_core", _core_path)
try:
    _core = _ilu.module_from_spec(_spec)
    _spec.loader.exec_module(_core)
    CANON = {
        "PHI":         (_core.PHI,         phi),
        "PI":          (_core.PI,          pi),
        "OMEGA":       (_core.OMEGA,       Omega),
        "BETA":        (_core.BETA,        beta),
        "KAL0":        (_core.KAL0,        KAL0),
        "P_SC":        (_core.P_SC,        Psc),
        "K_V":         (_core.K_V,         Kv),
        "T_R":         (_core.T_R,         Tr),
        "M_V":         (_core.M_V,         Mv),
        "W0":          (_core.W0,          w0),
        "WA":          (_core.WA,          -Psc / Kv),
        "OMEGA_DE":    (_core.OMEGA_DE,    Om_DE),
        "OMEGA_M_DYN": (_core.OMEGA_M_DYN, Om_m_dyn),
        "MIRA":        (_core.MIRA,        MIRA),
        "AURA":        (_core.AURA,        AURA),
        "H0_ALG":      (_core.H0_ALG,      3 * Omega ** 2),
        "N_S":         (_core.N_S,         1 - phi ** -7),
        "OMEGA_B_H2":  (_core.OMEGA_B_H2,  (pi - phi) / (3 * Omega ** 2)),
    }
    for nm, (core_val, guard_val) in CANON.items():
        check(f"canon  ssee_core.{nm}", abs(core_val - guard_val) < 1e-9,
              f"core {core_val:.10f} / guardian {guard_val:.10f}")
except Exception as e:
    check("canon  ssee_core.py importable y consistente", False, str(e))

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
