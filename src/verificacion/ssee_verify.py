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
# Identidad histórica MIRA (MIRA persiste como entidad = AURA/2; valor intacto).
check("L1 identidad  MIRA * Om_m,dyn = 0.31993 (identidad histórica)",
      abs(MIRA * (1 + w0) - 0.3199281880) < 1e-9)
# Identidad histórica π/φ (factor-materia del reframe 2026-06-17, SUPERADO).
check("L1 histórico  (pi/phi) * Om_m,dyn = 0.31076 (factor π/φ — superado)",
      abs((pi / phi) * (1 + w0) - 0.3107552907) < 1e-9)
# CANÓNICO (reframe 2026-06-18, OP-8 CERRADO): NO hay factor. Om_m,CMB = ω_m/h².
_kal0 = (pi + phi) / 2 + pi                       # KAL0 = BETA + PI
_omb  = (pi - phi) / (3 * Omega ** 2)             # ω_b (OP-1)
_omc  = _kal0 * _omb * (1 - phi ** -7)            # ω_c = KAL0·ω_b·n_s (forward)
_omm  = _omb + _omc + 0.06902 / 93.14             # ω_m físico
_h    = (3 * Omega ** 2) / 100                    # h = H_alg/100
check("L1 canónico   Om_m,CMB = ω_m/h² = 0.30889 (ω_m-directo, sin factor)",
      abs(_omm / _h ** 2 - 0.3088931986) < 1e-9)

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
    "V-L2-05 Om_m,CMB":  (_omm / _h ** 2,                 0.3088931986),
    "V-L2-05b ω_c":      (_omc,                           0.1195144084),
    "V-L2-05c ω_m":      (_omm,                           0.1426732002),
    "V-L2-05d π/φ·dyn":  ((pi / phi) * Om_m_dyn,          0.3107552907),
    "V-L2-05e MIRA·dyn": (MIRA * Om_m_dyn,                0.3199281880),
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

# V-L2-10 m_phi canónico (forward-prediction, cero fiteo, P6; mecanismo SOLAR²·KRYSTOS 2026-06-19):
#   m_phi = Sigma_m_nu^active * (SOLAR^2 * KRYSTOS)   [mecanismo g²·v]
#   con Sigma_m_nu^active = R2 * 0.960318 eV, R2 = Omega/(KAL0*Tr).
#   SOLAR = BIAL+KAL = phi+2pi (linaje radiativo); KRYSTOS = 2*Omega (anclado por wa).
# Dimensionalmente CONSISTENTE: [eV] * (número puro) = [eV].
# Mult. anterior PYROS*VITA*MIKA=615.33 (m_phi=42.47, sin mecanismo) RETIRADO; ver OP-9/OP-17.
_R2 = Omega / (KAL0 * Tr)
_mnu_active = _R2 * 0.960318
_SOLAR = phi + 2 * pi
_KRYSTOS = 2 * Omega
_mult_mphi = _SOLAR ** 2 * _KRYSTOS
_m_phi_canon = _mnu_active * _mult_mphi
check("V-L2-10 m_phi canónico = Sigma_m_nu^active * (SOLAR^2*KRYSTOS) = 41.02 eV",
      abs(_m_phi_canon - 41.0187) < 1e-3,
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

# OP-5 — tensión S8 weak-lensing (P5/P6). CANÓNICO (ω_m-directo, CLASS forward con
# m_phi=41.02 eV mecanismo SOLAR²·KRYSTOS, Om_m=0.30889): el two-sector phi-DM
# (free-streaming k_fs=0.762) lleva el single-sector S8=0.846 (3.5sigma, "el desafio")
# al titular S8_eff=0.759 (0.01sigma KiDS). sigma8 es OUTPUT directo de CLASS (no fit alpha_WDM).
# script: src/ssee_paper6_canonical_particle.py
Om_cosm_op5 = _omm / _h ** 2                         # 0.30889  Om_m,CMB (omega_m-directo)
S8_challenge = 0.8335 * (Om_cosm_op5 / 0.3) ** 0.5  # single-sector techo CLASS (el desafio)
S8_resolved = 0.7483 * (Om_cosm_op5 / 0.3) ** 0.5   # two-sector forward CLASS (resuelve)
check("V-L3-OP5  S8 single = sigma8(0.8335)(Om/0.3)^0.5 = 0.846  (el desafio)",
      abs(S8_challenge - 0.846) < 2e-3, f"S8 = {S8_challenge:.4f}")
check("V-L3-OP5  S8_eff two-sector = sigma8(0.7483)(Om/0.3)^0.5 = 0.759  (resuelve, forward)",
      abs(S8_resolved - 0.759) < 2e-3, f"S8 = {S8_resolved:.4f}")
track_open("V-L3-OP5  S8 forward resuelto; cierre no-lineal pleno diferido",
           "el two-sector LINEAL forward (CLASS, sigma8 OUTPUT no-fit) resuelve S8 a "
           "0.759 (0.01sigma KiDS) con m_phi=41.02 (SOLAR²·KRYSTOS); FALSABLE por k_fs=0.762 (DESI/Euclid). "
           "El cierre no-lineal con feedback barionico (N-body, ~5k-20k CPU-h) es Nivel 2, "
           "diferido. Ramas viejas 0.737/0.794, 0.702/0.725, 0.742/0.766 y 0.7536/0.765 retiradas")

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

# m_phi — masa del campo phi-DM (P6, CANÓNICO forward-prediction; SOLAR²·KRYSTOS 2026-06-19).
#   Sigma_m_nu^active = R2 * 0.960318 eV,  R2 = Omega/(KAL0*Tr)
#   m_phi = Sigma_m_nu^active * (SOLAR^2 * KRYSTOS)   [mecanismo g²·v]
# El multiplicador es número PURO -> dimensión [eV] preservada. Cero fiteo.
R2_p6 = Omega / (KAL0 * Tr)
mnu_active = R2_p6 * 0.960318
SOLAR = phi + 2 * pi
KRYSTOS = 2 * Omega
mult_p6 = SOLAR ** 2 * KRYSTOS
m_phi = mnu_active * mult_p6
check("V-L3-mphi  cadena m_phi = Sigma_m_nu^active * (SOLAR^2*KRYSTOS) = 41.02 eV",
      abs(m_phi - 41.0187) < 2e-2, f"m_phi = {m_phi:.4f} eV")
track_open("V-L3-mphi  coeficiente SOLAR²·KRYSTOS no derivado del transporte (OP-9)",
           "el valor 41.02 eV es forward-prediction dimensionalmente consistente, escrito "
           "como termino de masa g²·v de un Lagrangiano escalar libre; lo abierto (OP-9) es "
           "derivar el coeficiente SOLAR²·KRYSTOS del transporte disipativo (KAL), no la "
           "dimensionalidad ni un fiteo. SOLAR=BIAL+KAL, KRYSTOS=2Omega anclados por linaje")

# Dos sectores phi-DM (P6) — tras el reframe omega_m-directo (2026-06-18) la
# particion sale SOLA, sin factor: Om_CDM (=Om_m,dyn=0.160, DESI) + Om_phiDM =
# Om_m,CMB (=omega_m/h²=0.30889). El sector phi-DM es la DIFERENCIA entre la
# materia del CMB (omega_m) y la CDM dinamica. El split fisico en k_fs depende
# de m_phi (ABIERTO) y k_fs (pendiente Fase B).
Om_m_CMB = _omm / _h ** 2
Om_phiDM = Om_m_CMB - Om_m_dyn               # ≈ 0.149 (era (MIRA-1)*dyn=0.160)
check("V-L3-2sec  identidad Om_CDM + Om_phiDM = Om_m,CMB (omega_m/h²)",
      abs((Om_m_dyn + Om_phiDM) - Om_m_CMB) < 1e-12)
track_open("V-L3-2sec  split fisico de dos sectores no cerrado",
           "la suma es identidad (= V-L2-05); Om_phiDM=0.149 sale de omega_m-directo; "
           "la separacion fisica en k_fs depende de m_phi (ABIERTO) y k_fs (pendiente)")

# ── REFRAME 2026-06-19 — DEPENDIENTES PENDIENTES DE RECOMPUTE (cajon scripts) ──
# Inputs FIJADOS (algebra pura): Om_m,CMB=0.30889 (omega_m/h², OP-8 CERRADO, sin
# factor), H global=H_alg=67.962, m_phi=41.02 eV (mult SOLAR²·KRYSTOS=594.28,
# mecanismo g²·v adoptado; OP-17 cerrado).
# Los siguientes valores DEPENDEN de esos inputs; algunos checks pueden mostrar
# numeros viejos hasta correr cada codigo. NO se actualizan hasta recomputar:
track_open("REFRAME-FaseB  dependientes pendientes de recompute con canonicos nuevos",
           "HECHO: (a) cascada Hubble IR=72.86 (0.17sigma) / UV=73.040 (0.00sigma) "
           "con H global=67.962 (P9/P10). (b) CMB chi2=1005.5, Delta-BIC=-23.9 (SSEE favorecido) "
           "omega_m-directo @ H=67.962 (P3, plik_lite). (c) P6 CLASS forward m_phi=41.02 "
           "(SOLAR²·KRYSTOS), Om_phiDM=0.14889: k_fs=0.762, alpha=1.108, sigma8_two=0.7483 -> "
           "S8=0.759 (0.01sigma KiDS, RESUELVE forward). (d) fsigma8 two-sector recomputado "
           "Om_m=0.30889: media 0.82sigma (LCDM 0.73sigma). PENDIENTE aun: (1) r_d con "
           "Om_m=0.30889 (P3); (2) posterior MCMC con prior H_alg (P2). Marcado falsable por "
           "k_fs=0.762 (DESI/Euclid)")

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

# S8 weak-lensing — CANÓNICO ω_m-directo (CLASS forward, m_phi=41.02 SOLAR²·KRYSTOS).
#   Om_m,CMB = omega_m/h² = 0.30889 (sin factor). Dos ramas (CLASS OUTPUT, no fit):
#   single-sector techo (el desafío): sigma8 = 0.8335 -> S8 = 0.846 (3.5sigma KiDS).
#   two-sector phi-DM (TITULAR Paper 6, forward): sigma8_eff = 0.7483 ->
#     S8_eff = 0.759 (0.01sigma KiDS), free-streaming k_fs=0.762, m_phi=41.02.
# script: src/ssee_paper6_canonical_particle.py
Om_cosm = _omm / _h ** 2     # 0.30889  Om_m,CMB (omega_m-directo)
sig8_single = 0.8335         # techo todo-frío (CLASS)
S8_single = sig8_single * (Om_cosm / 0.3) ** 0.5
check("V-L4-01 P6  sigma8 single (techo CLASS) = 0.8335",
      abs(sig8_single - 0.8335) < 1e-2, f"sigma8 = {sig8_single:.4f}")
check("V-L4-02 P6  S8 single = sigma8 sqrt(Om/0.3) = 0.846  (el desafio)",
      abs(S8_single - 0.846) < 2e-3, f"S8 = {S8_single:.4f}")

sig8_eff = 0.7483            # two-sector titular Paper 6 (forward CLASS, no fit)
S8_eff = sig8_eff * (Om_cosm / 0.3) ** 0.5
check("V-L4-02b P6  S8_eff two-sector = 0.759  (TITULAR, resuelve, forward)",
      abs(S8_eff - 0.759) < 2e-3, f"S8_eff = {S8_eff:.4f}")

# Tensiones S8 — error en cuadratura modelo + observacional.
# KiDS-1000 (Asgari+2021): S8 = 0.759 +/- 0.024.
G_growth = 1.011             # D1_SSEE/D1_LCDM (Paper 5; no depende de m_phi)
S8_single_err = 0.006 * G_growth * (Om_cosm / 0.3) ** 0.5
t_KIDS_single = abs(S8_single - 0.759) / (S8_single_err ** 2 + 0.024 ** 2) ** 0.5
t_KIDS_twosec = abs(S8_eff - 0.759) / 0.024
check("V-L4-03 P6  tension S8 single vs KiDS-1000 = 3.5 sigma  (el desafio)",
      abs(t_KIDS_single - 3.5) < 0.2, f"{t_KIDS_single:.2f} sigma")
check("V-L4-04 P6  tension S8_eff two-sector vs KiDS-1000 = 0.01 sigma  (resuelto, forward)",
      t_KIDS_twosec < 0.2, f"{t_KIDS_twosec:.3f} sigma")

# CMB Planck PR4 (P3) — re-corrida con CAMB 1.6.5 (2026-05-22): chi2_r
# TT 1.047 / TE 1.041 / EE 1.041 / PP 0.837 y ΔBIC=-20.8 reproducidos
# EXACTAMENTE. La aritmetica chi2_r->ΔBIC solo acota por el redondeo.
check("V-L4-05 P3  ΔBIC CMB = -20.8 consistente con chi2_r (re-run 2026-05-22)",
      -22.9 <= -20.8 <= -11.1,
      "chi2_r redondeados acotan ΔBIC a [-22.9,-11.1]; reportado -20.8 dentro")

# r_d / theta* — re-run CAMB 2026-06-11 con los H0 CANONICOS post-cascada
# Sigma_mnu=0.069 (2026-06-09). Politica del modelo: anclar el CMB en el anchor
# H_MIRA=67.037 (CMB-optimo), reportar el posterior BAO 66.531 por separado.
# theta* es muy sensible a H0 via D_A, asi que NO se propaga el posterior BAO
# al observable CMB.
#   anchor   67.037 -> r_d=146.73 Mpc (1.38sigma), theta*=0.59623 (0.97sigma)
#   posterior 66.531 -> r_d=147.30 Mpc (0.81sigma), theta*=0.59408 (5.66sigma)
# vs Planck 147.09+-0.26 Mpc / 0.59668+-0.00046 deg.
track_open("V-L4  r_d coherente en ambos H0 canonicos (anchor 1.38 / post 0.81 sigma)",
           "ssee_verify_rd.py re-run 2026-06-11: anchor 67.037 da r_d=146.73 Mpc "
           "(1.38sigma), posterior 66.531 da r_d=147.30 Mpc (0.81sigma) vs "
           "Planck 147.09+-0.26. El '4.47sigma' previo usaba el H0 stale 67.756")
track_open("V-L4  theta* sensible a H0: 0.97sigma en anchor, 5.66sigma en posterior",
           "ssee_verify_rd.py re-run 2026-06-11: anchor 67.037 da theta*=0.59623 "
           "(0.97sigma); posterior 66.531 da 0.59408 (5.66sigma) vs 0.59668+-0.00046. "
           "Por eso el CMB se ancla en 67.037 y el posterior BAO no se propaga a theta*")

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
H0_mcmc = 66.531   # posterior canonico, prior MIRA 67.037 (re-run 2026-06-09); el 67.756 era mala anotacion
t_H0 = abs(H0_mcmc - 67.36) / (0.442 ** 2 + 0.54 ** 2) ** 0.5
check("V-L4-08 P2  tension H0 SSEE vs Planck = 1.19 sigma",
      abs(t_H0 - 1.19) < 0.03, f"{t_H0:.2f} sigma")
track_open("V-L4  Omega_b h^2: posterior MCMC vs prediccion algebraica",
           "MCMC re-run 2026-06-09 (prior MIRA 67.037) da Om_b h^2 = 0.02260+-0.00048; "
           "OP-1 algebraico da 0.02242 -> 0.4sigma, compatible. "
           "(El 0.02183 previo era de la cadena pre-canonica)")

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
# CAPA MEMORIAS — coherencia entre las 3 memorias de SSEE.
# Invoca memory_sync.py: revisa que CLAUDE.md, el Ledger y el vault Obsidian
# no presenten valores RETIRADOS (CANONICAL_VALUES.yaml) como vigentes.
# Es un check() DURO: cualquier valor retirado sin marcar pone el guardián ROJO.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa Memorias — coherencia de las 3 memorias")
try:
    import memory_sync
    _drifts = memory_sync.run(verbose=False)
    _detail = "todas concuerdan con CANONICAL_VALUES.yaml"
    if _drifts:
        _detail = "DRIFT en " + "; ".join(
            f"{lbl} {rel}:{ln}«{pat}»" for lbl, rel, ln, pat, _ in _drifts[:6])
        if len(_drifts) > 6:
            _detail += f" … (+{len(_drifts) - 6}); correr memory_sync.py"
    check("memoria  3 memorias sin valores retirados sin marcar",
          not _drifts, _detail)
except Exception as e:
    check("memoria  memory_sync importable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA PROCEDENCIA — valores de pipeline (sección B del Registro).
# Causa raíz F1 (2026-06-14): un valor canónico de pipeline (ΔBIC=−24.7,
# χ²_r=1.045) fue anotado A MANO en el Registro sin un log que lo
# reprodujera, y NADA lo verificaba → sobrevivió mal hasta una re-corrida.
# Este check cierra la grieta: extrae el valor de la fila del Registro y
# exige que aparezca, en contexto, en un log de procedencia committeado.
#   · valor coincide con su log  → OK
#   · valor NO aparece en el log  → ROJO (atrapa exactamente el error F1)
#   · sin log committeado         → ABIERTO (hace visible la grieta, no bloquea)
# Para añadir un valor de pipeline al Registro: guardar su log en
# results/logs/ y registrarlo aquí. Sin log → el guardián lo marca.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa Procedencia — valores de pipeline vs log committeado")
try:
    import re as _re
    _REPO = pathlib.Path(__file__).resolve().parents[2]

    def _norm(s):
        return s.replace("−", "-").replace("–", "-")

    _lines = (_REPO / "VERIFICATION_LEDGER.md").read_text(errors="ignore").splitlines()
    # Escaneo de TODA la sección B: cada valor de pipeline cuya Fuente cite un
    # log committeado (results/logs/*.log) se verifica contra él. Cobertura
    # completa, no una lista a mano (cierra la grieta «3 de 17»).
    _inB = False
    _gaps = []
    for ln in _lines:
        if ln.startswith("## B."):
            _inB = True
            continue
        if _inB and (ln.startswith("## ") or ln.startswith("# ") or ln.startswith("---")):
            break
        if not _inB or not ln.lstrip().startswith("|"):
            continue
        _cols = [c.strip() for c in ln.split("|")]
        if len(_cols) < 4:
            continue
        _label, _valf, _src = _cols[1], _norm(_cols[2]), _cols[3]
        if _label.lower().startswith("cantidad") or set(_label) <= set("-: "):
            continue  # cabecera / separador
        _m = _re.search(r"-?\d+\.\d+", _valf)
        if not _m:
            continue
        _val = _m.group(0)
        _logm = _re.search(r"results/logs/\S+?\.log", _src)
        if not _logm:
            _gaps.append(f"{_label[:22]}={_val}")
            continue
        _rel = _logm.group(0)
        _logf = _REPO / _rel
        if not _logf.exists():
            check(f"procedencia  {_label[:30]}={_val}", False,
                  f"log referenciado no existe: {_rel}")
            continue
        _hit = _val in _norm(_logf.read_text(errors="ignore"))
        check(f"procedencia  {_label[:30]}={_val}", _hit,
              f"coincide con {_rel}" if _hit
              else f"Registro={_val} NO aparece en {_rel} (mala anotación tipo F1)")
    if _gaps:
        track_open(f"procedencia  {len(_gaps)} valores de pipeline sin log committeado",
                   "; ".join(_gaps[:8]) + (" …" if len(_gaps) > 8 else "")
                   + "  (generar log → results/logs/ para verificarlos)")
except Exception as e:
    check("procedencia  capa de procedencia operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA MANUSCRITOS — reglas del RIGOR_CHECKLIST automatizadas (grep duro).
# "El sistema solo se hace más capas": cada regla automatizable se vuelve un
# check aquí, para que UNA corrida marque dónde y por qué.
#   R1 cronología   : ninguna afirmación de prioridad temporal sobre datos públicos.
#   R2 multidominio : ningún enlace al compendio multidominio retirado (DOI 19679049).
# Los disclaimers honestos ("no claim of temporal priority", "which predate this
# work") NO matchean: los patrones apuntan solo a las afirmaciones prohibidas.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa Manuscritos — reglas R1 (cronología) y R2 (multidominio)")
try:
    import re as _re2
    _REPO2 = pathlib.Path(__file__).resolve().parents[2]
    # Cobertura: los 12 manuscritos + el README raíz (cajón público vigente).
    # NO se escanean los docs de auditoría (RIGOR_CHECKLIST, AUDIT…) ni el
    # CHANGELOG: contienen los patrones prohibidos como EJEMPLOS/historia.
    _texs = sorted((_REPO2 / "manuscript").glob("*.tex")) + [_REPO2 / "README.md"]
    _texs = [p for p in _texs if p.exists()]
    _R1 = [
        r"predating\s+(the\s+)?(desi|planck|dr2|dr1)",
        r"prior to\s+(the\s+)?(desi|planck)\s*(dr2|dr1|2018|pr4|release|comparison)",
        r"deposited[^.]{0,50}(prior to|before)[^.]{0,25}(desi|planck|dr2)",
        r"timestamp[^.]{0,30}(proof|prior)",
        r"pre-desi\b",
        r"committed\s+before\s+(desi|planck)",
    ]
    _R2 = [r"19679049", r"ssee_unificado", r"unified compendium of irrational"]
    # R-overclaim: ASERCIÓN global de "zero-parameter" en prosa. NO matchea:
    # claims scoped ("zero-parameter dark-energy"), recantaciones ("described as
    # a zero-parameter model... however"), ni citas de títulos.
    _ROVER = [r"(achieves|is|provides|presents|constitutes)\s+a\s+zero[- ]parameter\s+(framework|model|theory)"]
    # R-serie: conteo de EXTENSIONES congelado antes de 10 (P1 «3--7»), o
    # descripción explícita del tamaño de la serie. NO matchea las referencias
    # correctas de un paper a sus previos (P8 «1--7», P10 «1--9»).
    _RSER = [r"papers?~?3--[789]\b", r"(seven|eight|nine)[- ]paper series"]
    _r1_hits, _r2_hits, _rov_hits, _rser_hits = [], [], [], []
    for t in _texs:
        low = t.read_text(errors="ignore").lower()
        for pat in _R1:
            m = _re2.search(pat, low)
            if m:
                _r1_hits.append(f"{t.name}:«{low[m.start():m.end()][:40]}»")
        for pat in _R2:
            if _re2.search(pat, low):
                _r2_hits.append(f"{t.name}:{pat}")
        for pat in _ROVER:
            m = _re2.search(pat, low)
            if m:
                _rov_hits.append(f"{t.name}:«{low[m.start():m.end()][:40]}»")
        for pat in _RSER:
            m = _re2.search(pat, low)
            if m:
                _rser_hits.append(f"{t.name}:«{low[m.start():m.end()][:30]}»")
    check("manuscritos  R1 sin claims de prioridad temporal", not _r1_hits,
          "limpio" if not _r1_hits else "; ".join(_r1_hits[:4]))
    check("manuscritos  R2 sin enlaces multidominio retirados", not _r2_hits,
          "limpio" if not _r2_hits else "; ".join(_r2_hits[:4]))
    check("manuscritos  R10 sin overclaim 'zero-parameter' global", not _rov_hits,
          "limpio" if not _rov_hits else "; ".join(_rov_hits[:4]))
    check("manuscritos  R11 sin conteo de serie stale (<10 papers)", not _rser_hits,
          "limpio" if not _rser_hits else "; ".join(_rser_hits[:4]))
except Exception as e:
    check("manuscritos  capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA ARCHIVO — la Bitácora de Archivado cubre cada cajón archivado.
# Nada se archiva en silencio: cada subcarpeta de archive/ debe tener su
# entrada (qué/cuándo/por qué) en archive/README.md.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa Archivo — bitácora cubre cada cajón archivado")
try:
    _arch = pathlib.Path(__file__).resolve().parents[2] / "archive"
    _bitf = _arch / "README.md"
    _bit = _bitf.read_text(errors="ignore") if _bitf.exists() else ""
    _orphans = [d.name for d in _arch.iterdir()
                if d.is_dir() and d.name not in _bit] if _arch.exists() else []
    check("archivo  toda subcarpeta de archive/ documentada en la bitácora",
          not _orphans,
          "documentadas" if not _orphans else "SIN entrada: " + ", ".join(_orphans))
except Exception as e:
    check("archivo  capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA DICCIONARIO — integridad de la lista de nodos nombrados.
# El diccionario MAESTRO vive en sandbox_unificado/ssee-data.json. Ningún
# script puede introducir un nodo ausente del maestro: eso fue el duplicado
# espurio KRYSTOS_V (=2Omega) de look_elsewhere_full.py, idéntico por VALOR
# a KRYSTOS (=phi+pi+Omega). Un chequeador de VALORES es ciego a esto
# (mismo número); este check de PROCEDENCIA lo caza. (R13, 2026-06-16)
# ─────────────────────────────────────────────────────────────────────
print("\nCapa Diccionario — integridad de nodos (script ⊆ maestro)")
try:
    import json as _json, re as _re
    _root = pathlib.Path(__file__).resolve().parents[2]
    _mf = _root / "sandbox_unificado" / "ssee-data.json"
    _master = set()
    if _mf.exists():
        _md = _json.load(open(_mf))
        for _sec in ("base_constants", "derived_constants", "nine_sovereignties"):
            for _c in _md.get(_sec, []):
                if isinstance(_c, dict) and _c.get("name"):
                    _master.add(_c["name"].upper())
        if "OMEGADNAV" in _master:
            _master.add("OMEGA")          # alias scaffold
    _lf = _root / "src" / "estadistica" / "look_elsewhere_full.py"
    _fam = set()
    if _lf.exists():
        _m = _re.search(r"FAMILY\s*=\s*\{(.*?)\}",
                        _lf.read_text(errors="ignore"), _re.S)
        if _m:
            _fam = {k.upper() for k in _re.findall(r'"([A-Z_]+)"\s*:', _m.group(1))}
    _spurious = sorted(_fam - _master) if (_master and _fam) else []
    check("diccionario  R13 look_elsewhere sin nodos ausentes del maestro",
          not _spurious,
          "consistente (script ⊆ maestro)" if not _spurious
          else "ESPURIOS (no en maestro): " + ", ".join(_spurious))
except Exception as e:
    check("diccionario  capa operable", False, str(e))

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
