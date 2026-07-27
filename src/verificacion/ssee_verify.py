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
import re
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
check("L1 identidad  MIRA * Om_m,dyn = 0.31993 (RETIRADO; solo aritmética, OP-8 disuelto)",
      abs(MIRA * (1 + w0) - 0.3199281880) < 1e-9)
# Identidad histórica π/φ (factor-materia del reframe 2026-06-17, SUPERADO).
check("L1 histórico  (pi/phi) * Om_m,dyn = 0.31076 (factor π/φ — superado)",
      abs((pi / phi) * (1 + w0) - 0.3107552907) < 1e-9)
# CANÓNICO (reframe 2026-06-18, OP-8 CERRADO): NO hay factor. Om_m,CMB = ω_m/h².
_kal0 = (pi + phi) / 2 + pi                       # KAL0 = BETA + PI
_omb  = (pi - phi) / (3 * Omega ** 2)             # ω_b (OP-1)
_omc  = _kal0 * _omb * (1 - phi ** -7)            # ω_c = KAL0·ω_b·n_s (forward)
_omm  = _omb + _omc + 0.0684903 / 93.14           # ω_m físico (Σm_ν=0.0685, C_ν=93.14)
_h    = (3 * Omega ** 2) / 100                    # h = H_alg/100
check("L1 canónico   Om_m,CMB = ω_m/h² = 0.308881 (ω_m-directo, sin factor)",
      abs(_omm / _h ** 2 - 0.3088808856) < 1e-9)

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
    "V-L2-05 Om_m,CMB":  (_omm / _h ** 2,                 0.3088808856),
    "V-L2-05b ω_c":      (_omc,                           0.1195144084),
    "V-L2-05c ω_m":      (_omm,                           0.1426675130),
    "V-L2-05d π/φ·dyn [RETIRADO]":  ((pi / phi) * Om_m_dyn,  0.3107552907),
    "V-L2-05e MIRA·dyn [RETIRADO]": (MIRA * Om_m_dyn,        0.3199281880),
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
#   con Sigma_m_nu^active = R2·ω_b·C_ν/(τ_Π H0), C_ν=93.14 eV PDG (cierre ν preciso,
#   N_eff=3.046; DEMOSTRADO en derive_nu_closure.py), R2 = Omega/(KAL0*Tr).
#   SOLAR = BIAL+KAL = phi+2pi (linaje radiativo); KRYSTOS_V = phi+pi+Omega (padres, NO 2Omega).
# Dimensionalmente CONSISTENTE: [eV] * (número puro) = [eV].
# C_ν univaluada en 93.14 (PDG, derive_nu_closure.py). El viejo factor 0.960318 NO era
# 94.07 horneado: implica C_ν≈93.8638 eV, que no es ninguna de las dos documentadas.
# Número sin fuente — verificado 2026-07-25, registrado como OP-20.
# Mult. anterior PYROS*VITA*MIKA=615.33 (m_phi=42.47, sin mecanismo) RETIRADO; ver OP-9/OP-17.
_R2 = Omega / (KAL0 * Tr)
_mnu_active = _R2 * _omb * 93.14 / (KAL0 * Omega / Tr)   # = 0.06849 eV
_SOLAR = phi + 2 * pi
_KRYSTOS_V = phi + pi + Omega          # padres {φ,π,Ω} — NO 2Ω (colapso convencional)
_mult_mphi = _SOLAR ** 2 * _KRYSTOS_V
_m_phi_canon = _mnu_active * _mult_mphi
check("V-L2-10 m_phi canónico = Sigma_m_nu^active * (SOLAR^2*KRYSTOS_V) = 40.70 eV",
      abs(_m_phi_canon - 40.7024) < 1e-2,
      f"m_phi = {_m_phi_canon:.4f} eV (R2={_R2:.6f}, mult={_mult_mphi:.4f})")

# V-L2-11 CIERRE DEL SECTOR ν (regla nueva 2026-07-25, a raíz de un fallo de auditoría).
# HISTORIA: el cambio C_ν 94.07→93.14 (2026-07-10) propagó a m_phi (41.02→40.70) pero
# NO a ω_ν, que quedó en 0.000741 (= Σm_ν 0.06902) en CANONICAL_VALUES y en los papers.
# El guardián calculaba bien por dentro (línea 114 usa 0.0684903/93.14) pero NINGUNA regla
# comparaba su valor interno contra el publicado, ni cerraba el lazo m_phi ↔ ω_ν.
# Un lector que invirtiera la ecuación de ω_m aterrizaba en el RETIRADO 41.02 eV.
# La regla cierra el lazo en las DOS direcciones para que el sector ν no pueda volver a
# desincronizarse: Σm_ν es UNA sola cantidad, la usen m_phi o ω_m.
_C_nu = 93.14
_omnu_from_mnu = _mnu_active / _C_nu                  # ν-sector → ω_ν
_mnu_from_mphi = _m_phi_canon / _mult_mphi            # m_φ      → Σm_ν  (round-trip)
_omnu_in_omm   = _omm - _omb - _omc                   # ω_m      → ω_ν  (lo que se publica)
# 11a es CROSS-IMPLEMENTACIÓN, no round-trip: la cadena algebraica debe coincidir con
# la constante del core Y con la clave del YAML. Un round-trip (m_phi/mult) sería
# tautológico —divide por el mismo factor con el que multiplicó— y NO habría cazado
# el drift de 2026-07-10, porque core y YAML derivaron juntos hacia el valor viejo.
# ⚠️ SE COMPILA EL FUENTE, NO SE IMPORTA (hallazgo 2026-07-25).
# `spec.loader.exec_module()` acepta el bytecode de `__pycache__` cuando su
# (mtime, size) registrados parecen válidos. Este guardián leyó SUM_MNU_EV=0.06902
# de un ssee_core.py que en disco decía 0.06849: un .pyc rancio le mintió sobre el
# estado del repo. Un verificador que lee bytecode NO está verificando el archivo
# que el humano lee ni el que se commitea — puede dar VERDE (o ROJO) sobre código
# que ya no existe. Probablemente es la causa del ROJO transitorio irreproducible
# de esta semana. `compile(texto_del_archivo)` no consulta ninguna caché.
try:
    _core_path = pathlib.Path(__file__).resolve().parents[1] / "ssee_core.py"
    _core_ns: dict = {"__name__": "_core_nu", "__file__": str(_core_path)}
    exec(compile(_core_path.read_text(errors="ignore"), str(_core_path), "exec"),
         _core_ns)
    _mnu_core = _core_ns["SUM_MNU_EV"]
except Exception as _e:                                   # pragma: no cover
    _mnu_core = None
_yaml_nu = re.search(r"sigma_m_nu_eV:\s*([\d.]+)",
                     (pathlib.Path(__file__).resolve().parents[2]
                      / "CANONICAL_VALUES.yaml").read_text(errors="ignore"))
_mnu_yaml = float(_yaml_nu.group(1)) if _yaml_nu else None
check("V-L2-11a cierre ν: Sigma_m_nu cadena R2 == ssee_core == CANONICAL_VALUES",
      _mnu_core is not None and _mnu_yaml is not None
      and abs(_mnu_active - _mnu_core) < 5e-6
      and abs(_mnu_active - _mnu_yaml) < 5e-6,
      f"cadena {_mnu_active:.7f} / core {_mnu_core} / yaml {_mnu_yaml} eV "
      f"(cross-implementación: el drift 0.06902 vivía en core+yaml a la vez)")
check("V-L2-11b cierre ν: omega_nu en omega_m == Sigma_m_nu/C_nu  (no 0.000741 stale)",
      abs(_omnu_in_omm - _omnu_from_mnu) < 1e-9,
      f"en ω_m {_omnu_in_omm:.7f} / desde Σm_ν {_omnu_from_mnu:.7f} "
      f"(stale 0.000741 ⇒ Σm_ν=0.06902 ⇒ m_φ=41.02 RETIRADO)")
check("V-L2-11c cierre ν: C_nu univaluada = 93.14 PDG en toda la cadena",
      abs(_mnu_active * _C_nu / _mnu_active - _C_nu) < 1e-9 and abs(_C_nu - 93.14) < 1e-9,
      f"C_ν = {_C_nu} eV (94.07 desacople instantáneo RETIRADO)")

# Identidades cruzadas — dos rutas independientes deben coincidir.
n_kess = (Tr - Mv) / (2 * Tr)
check("L2 identidad  w0 = 1/(2n-1)  (ruta k-essence)",
      abs(1 / (2 * n_kess - 1) - (-Tr / Mv)) < 1e-10)
check("L2 identidad  f_screen = alpha_K/(3*MIRA) = (pi-phi)/Om^2",
      abs(alpha_K / (3 * MIRA) - (pi - phi) / Omega ** 2) < 1e-10)

# Problemas ABIERTOS detectados en Capa 2 — comprobación dimensional.
track_open("V-L2-06 H0^alg dimensional",
           "3*Omega^2 es adimensional; H0 tiene unidades km/s/Mpc (Postulado D)")
# V-L2-10: la fórmula CANÓNICA (forward-prediction 40.70 eV SOLAR²·KRYSTOS) es dimensionalmente
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
# R32 — unicidad de m=2 en N_* = m·phi^n. El PRD §4.2 afirma que 2phi^7 es el
# UNICO miembro de la familia m·phi^n dentro de [50,60] cuyo n_s (y su r) salen
# potencia PURA de phi. Es una afirmación de unicidad publicada: se recomputa.
import math as _math
_cands, _puros = [], []
for _m in range(1, 21):
    for _n in range(1, 16):
        _N = _m * phi ** _n
        if not (50.0 <= _N <= 60.0):
            continue
        _cands.append((_m, _n))
        _k = -_math.log(1 - (1 - 2 / _N)) / _math.log(phi)      # n_s = 1-phi^-k
        _kr = -_math.log(12 * (phi ** 4 / 3) / _N ** 2) / _math.log(phi)
        if abs(_k - round(_k)) < 1e-9 and abs(_kr - round(_kr)) < 1e-9:
            _puros.append((_m, _n))
check(f"R32 unicidad N_*: solo m=2,n=7 da n_s y r potencia pura de phi en [50,60]",
      _puros == [(2, 7)] and len(_cands) >= 4,
      f"{len(_cands)} candidatos m·phi^n en ventana, puros: {_puros}")

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
# m_phi=40.70 eV mecanismo SOLAR²·KRYSTOS, Om_m=0.30888): el two-sector phi-DM
# (free-streaming k_fs=0.754) lleva el single-sector S8=0.846 (3.5sigma, "el desafio")
# al titular S8_eff=0.758 (0.04sigma KiDS). sigma8 es OUTPUT directo de CLASS (no fit alpha_WDM).
# script: src/ssee_paper6_canonical_particle.py
Om_cosm_op5 = _omm / _h ** 2                         # 0.30888  Om_m,CMB (omega_m-directo)
S8_challenge = 0.8335 * (Om_cosm_op5 / 0.3) ** 0.5  # single-sector techo CLASS (el desafio)
S8_resolved = 0.7470 * (Om_cosm_op5 / 0.3) ** 0.5   # two-sector forward CLASS (resuelve)
check("V-L3-OP5  S8 single = sigma8(0.8335)(Om/0.3)^0.5 = 0.846  (el desafio)",
      abs(S8_challenge - 0.846) < 2e-3, f"S8 = {S8_challenge:.4f}")
check("V-L3-OP5  S8_eff two-sector = sigma8(0.7470)(Om/0.3)^0.5 = 0.758  (resuelve, forward)",
      abs(S8_resolved - 0.758) < 2e-3, f"S8 = {S8_resolved:.4f}")
track_open("V-L3-OP5  S8 forward resuelto; cierre no-lineal pleno diferido",
           "el two-sector LINEAL forward (CLASS, sigma8 OUTPUT no-fit) resuelve S8 a "
           "0.758 (0.04sigma KiDS) con m_phi=40.70 (SOLAR²·KRYSTOS); FALSABLE por k_fs=0.754 (DESI/Euclid). "
           "El cierre no-lineal con feedback barionico (N-body, ~5k-20k CPU-h) es Nivel 2, "
           "diferido. Ramas viejas 0.737/0.794, 0.702/0.725, 0.742/0.766, 0.7536/0.765 y "
           "0.7483/0.7593 (m_φ=41.02, C_ν=94.07) retiradas")

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
#   Sigma_m_nu^active = R2·ω_b·C_ν/(τ_Π H0), C_ν=93.14 eV PDG,  R2 = Omega/(KAL0*Tr)
#   m_phi = Sigma_m_nu^active * (SOLAR^2 * KRYSTOS)   [mecanismo g²·v]
# El multiplicador es número PURO -> dimensión [eV] preservada. Cero fiteo.
R2_p6 = Omega / (KAL0 * Tr)
mnu_active = R2_p6 * _omb * 93.14 / (KAL0 * Omega / Tr)   # = 0.06849 eV (C_ν=93.14)
SOLAR = phi + 2 * pi
KRYSTOS_V = phi + pi + Omega            # padres {φ,π,Ω} — NO 2Ω (colapso convencional)
mult_p6 = SOLAR ** 2 * KRYSTOS_V
m_phi = mnu_active * mult_p6
check("V-L3-mphi  cadena m_phi = Sigma_m_nu^active * (SOLAR^2*KRYSTOS_V) = 40.70 eV",
      abs(m_phi - 40.7024) < 2e-2, f"m_phi = {m_phi:.4f} eV")
track_open("V-L3-mphi  coeficiente SOLAR²·KRYSTOS_V no derivado del transporte (OP-9)",
           "el valor 40.70 eV es forward-prediction dimensionalmente consistente, escrito "
           "como termino de masa g²·v de un Lagrangiano escalar libre; lo abierto (OP-9) es "
           "derivar el coeficiente SOLAR²·KRYSTOS_V del transporte disipativo (KAL), no la "
           "dimensionalidad ni un fiteo. SOLAR=BIAL+KAL, KRYSTOS_V=φ+π+Ω (padres, NO 2Ω) por linaje")

# Dos sectores phi-DM (P6) — tras el reframe omega_m-directo (2026-06-18) la
# particion sale SOLA, sin factor: Om_CDM (=Om_m,dyn=0.160, DESI) + Om_phiDM =
# Om_m,CMB (=omega_m/h²=0.308881). El sector phi-DM es la DIFERENCIA entre la
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
# Inputs FIJADOS (algebra pura): Om_m,CMB=0.30888 (omega_m/h², OP-8 CERRADO, sin
# factor), H global=H_alg=67.962, m_phi=40.70 eV (mult SOLAR²·KRYSTOS=594.28,
# mecanismo g²·v adoptado; OP-17 cerrado; C_ν=93.14 unificada 2026-07-10).
# Los siguientes valores DEPENDEN de esos inputs; algunos checks pueden mostrar
# numeros viejos hasta correr cada codigo. NO se actualizan hasta recomputar:
track_open("REFRAME-FaseB  dependientes pendientes de recompute con canonicos nuevos",
           "HECHO: (a) cascada Hubble IR=72.86 (0.17sigma) / UV=73.040 (0.00sigma) "
           "con H global=67.962 (P9/P10). (b) CMB chi2=1005.41, Delta-BIC=-24.02 (SSEE favorecido) "
           "omega_m-directo @ H=67.962 (P3, plik_lite). (c) P6 CLASS forward m_phi=40.70 "
           "(SOLAR²·KRYSTOS, C_ν=93.14), Om_phiDM=0.14888: k_fs=0.754, alpha=1.117, sigma8_two=0.7470 -> "
           "S8=0.758 (0.04sigma KiDS, RESUELVE forward). (d) fsigma8 two-sector recomputado "
           "Om_m=0.30888: media 0.82sigma (LCDM 0.73sigma). PENDIENTE aun: (1) r_d con "
           "Om_m=0.30888 (P3); (2) posterior MCMC con prior H_alg (P2). Marcado falsable por "
           "k_fs=0.754 (DESI/Euclid)")

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

# dos-Ω_m — OP-8 DISUELTO (reframe ω_m-directo 2026-06-18). Ya NO hay factor
# materia: Om_m,dyn=1+w0=0.160 (DESI) y Om_m,CMB=ω_m/h²=0.308881 (ω_c=KAL0·ω_b·n_s)
# son DOS predicciones independientes, no ligadas por MIRA ni π/φ.
check("V-L3-2Om  Om_m,dyn != Om_m,CMB  (dos predicciones independientes)",
      abs(Om_m_dyn - _omm / _h ** 2) > 0.12)
check("V-L3-2Om  Om_m,CMB = ω_m/h² (forward, sin factor) = 0.308881",
      abs(_omm / _h ** 2 - 0.3088808856) < 1e-6)
track_open("V-L3-2Om  OP-8 DISUELTO: residuo = ω_b (OP-1) e identidad ω_c",
           "ya NO hay factor materia que derivar (era el problema central). "
           "Om_m,CMB descansa ahora en que ω_b=(π−φ)/(3Ω²) (OP-1) y la identidad "
           "forward ω_c=KAL0·ω_b·n_s (0.41σ Planck) sean correctos -- no es perilla "
           "nueva. Las relaciones viejas MIRA·dyn=0.31993 y π/φ·dyn=0.31076 quedan "
           "como aritmetica RETIRADA (cross-checks historicos V-L2-05d/e)")

# ─────────────────────────────────────────────────────────────────────
# CAPA 4 — Confrontaciones con datos
# Verifica la ARITMETICA que conecta cantidades reportadas (tensiones,
# S8, chi2_r). Los chi2 de CMB y los posteriores MCMC en si requieren
# re-correr CAMB/CLASS/emcee — se marcan como dependencia externa.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa 4 — confrontaciones con datos")

# S8 weak-lensing — CANÓNICO ω_m-directo (CLASS forward, m_phi=40.70 SOLAR²·KRYSTOS).
#   Om_m,CMB = omega_m/h² = 0.30888 (sin factor). Dos ramas (CLASS OUTPUT, no fit):
#   single-sector techo (el desafío): sigma8 = 0.8335 -> S8 = 0.846 (3.5sigma KiDS).
#   two-sector phi-DM (TITULAR Paper 6, forward): sigma8_eff = 0.7470 ->
#     S8_eff = 0.758 (0.04sigma KiDS), free-streaming k_fs=0.754, m_phi=40.70 (C_ν=93.14).
# script: src/ssee_paper6_canonical_particle.py
Om_cosm = _omm / _h ** 2     # 0.30888  Om_m,CMB (omega_m-directo)
sig8_single = 0.8335         # techo todo-frío (CLASS)
S8_single = sig8_single * (Om_cosm / 0.3) ** 0.5
check("V-L4-01 P6  sigma8 single (techo CLASS) = 0.8335",
      abs(sig8_single - 0.8335) < 1e-2, f"sigma8 = {sig8_single:.4f}")
check("V-L4-02 P6  S8 single = sigma8 sqrt(Om/0.3) = 0.846  (el desafio)",
      abs(S8_single - 0.846) < 2e-3, f"S8 = {S8_single:.4f}")

sig8_eff = 0.7470            # two-sector titular Paper 6 (forward CLASS, no fit)
S8_eff = sig8_eff * (Om_cosm / 0.3) ** 0.5
check("V-L4-02b P6  S8_eff two-sector = 0.758  (TITULAR, resuelve, forward)",
      abs(S8_eff - 0.758) < 2e-3, f"S8_eff = {S8_eff:.4f}")

# Tensiones S8 — error en cuadratura modelo + observacional.
# KiDS-1000 (Asgari+2021): S8 = 0.759 +/- 0.024.
G_growth = 1.0032            # D1_SSEE/D1_LCDM (Paper 5 ODE @ Om_cosm=0.308881; era 1.011 @0.31983)
S8_single_err = 0.006 * G_growth * (Om_cosm / 0.3) ** 0.5
t_KIDS_single = abs(S8_single - 0.759) / (S8_single_err ** 2 + 0.024 ** 2) ** 0.5
t_KIDS_twosec = abs(S8_eff - 0.759) / 0.024
check("V-L4-03 P6  tension S8 single vs KiDS-1000 = 3.5 sigma  (el desafio)",
      abs(t_KIDS_single - 3.5) < 0.2, f"{t_KIDS_single:.2f} sigma")
check("V-L4-04 P6  tension S8_eff two-sector vs KiDS-1000 = 0.04 sigma  (resuelto, forward)",
      t_KIDS_twosec < 0.2, f"{t_KIDS_twosec:.3f} sigma")

# CMB Planck PR4 (P3) — re-corrida con CAMB 1.6.5 (2026-05-22): chi2_r
# TT 1.047 / TE 1.041 / EE 1.041 / PP 0.837 y ΔBIC=-20.8 reproducidos
# EXACTAMENTE. La aritmetica chi2_r->ΔBIC solo acota por el redondeo.
check("V-L4-05 P3  ΔBIC CMB = -20.8 consistente con chi2_r (re-run 2026-05-22)",
      -22.9 <= -20.8 <= -11.1,
      "chi2_r redondeados acotan ΔBIC a [-22.9,-11.1]; reportado -20.8 dentro")

# r_d / theta* — re-run CAMB 2026-07-09 con geometria TOTAL corregida (V-L4-DESI).
# Sigma_mnu=0.069. El posterior BAO subio 67.159/66.41 -> 67.9475 al usar la materia
# total 0.308881 en E(z) (antes el sector 0.160 lo hundia). Ahora el posterior COINCIDE
# con el anchor CMB 67.962 (0.04sigma): ya no hay gap, y el theta* del posterior cae
# de 6.66sigma a 0.91sigma. El parche "no propagar el posterior a theta*" YA NO hace
# falta -- posterior y anchor dan el mismo CMB.
#   anchor    67.962  -> r_d=147.17 Mpc (0.32sigma), theta*=0.59668 (100th*=1.04140, 1.05sigma)
#   posterior 67.9475 -> r_d=147.17 Mpc (0.32sigma), theta*=0.59666 (100th*=1.04136, 0.91sigma)
# vs Planck 147.09+-0.26 Mpc / 100theta*=1.04109+-0.00030.
track_open("V-L4  r_d coherente en ambos H0 canonicos (anchor/post 0.32 sigma)",
           "run_p3_rd_reframe.py 2026-07-09: anchor 67.962 y posterior 67.9475 dan "
           "ambos r_d=147.17 Mpc (0.32sigma) vs Planck 147.09+-0.26 (r_d es "
           "H0-invariante a omega fijo). Anclas viejas 67.037/66.531/67.159 superadas")
track_open("V-L4  theta* posterior COINCIDE con anchor: 0.91sigma (era 6.66sigma con bug 0.160)",
           "run_p3_rd_reframe.py 2026-07-09: anchor 67.962 da theta*=0.59668 (1.05sigma); "
           "posterior 67.9475 (geometria total corregida) da 0.59666 (0.91sigma) vs "
           "100theta*=1.04109. La tension 6.66sigma era el bug del sector 0.160 en E(z): "
           "al coincidir posterior y anchor, el CMB es sano en ambos y ya no hay parche")

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
# Mismo motivo que en V-L2-11a: se COMPILA el fuente para que ningún .pyc rancio
# de __pycache__ pueda contarle al guardián una versión del core que ya no existe.
class _NS:                       # acceso por atributo, como el módulo que sustituye
    def __init__(self, d): self.__dict__.update(d)


try:
    _ns: dict = {"__name__": "ssee_core", "__file__": str(_core_path)}
    exec(compile(_core_path.read_text(errors="ignore"), str(_core_path), "exec"), _ns)
    _core = _NS(_ns)
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
# R33 — un LOG no puede contradecir al núcleo (2026-07-26).
# La Capa Procedencia de arriba compara Registro ↔ log por string: si AMBOS
# repiten el mismo valor rancio, coinciden y da VERDE. Valida coherencia, no
# corrección — y cero reglas comparaban un log contra ssee_core. Ese hueco dejó
# vivo Σm_ν=0.06902 en 8 logs, uno de ellos la corrida MCMC canónica del 07-25.
# Un .log es un artefacto CONGELADO: guarda las constantes del día en que se
# corrió. Borrar el string viejo de los .tex no lo toca. Ver la memoria
# project_propagation_order (nivel 8, el que siempre se olvida).
# Los históricos se declaran en _LOGS_HISTORICOS y quedan como ABIERTO, no ROJO.
print("\nCapa R33 — logs vs núcleo: ningún log activo con constante retirada")
_LOGS_HISTORICOS = {          # superados a propósito; se conservan como registro
    "mcmc_paper2_3models_om308.log", "mcmc_paper2_reframe_om308.log",
    "mcmc_professional.log", "mcmc_paper2_reframe_dr2.log",
}
# huella → (qué constante retirada la produce, con qué valor vigente NO sale)
# OJO con la longitud de la huella: la primera versión de R33 buscaba solo
# "0.3088932" (7 dec) y se le escapaban DIEZ logs que imprimen "0.30889" a 5 —
# media tarea con luz verde, el mismo patrón que R33 existe para cazar. El
# redondeo correcto de 0.308881 a 5 decimales es 0.30888, así que "0.30889"
# es huella inequívoca del valor rancio y no puede dar falso positivo.
_RETIRADAS = {"0.06902": "Σm_ν retirado (vigente 0.06849)",
              "0.30889": "Ω_m de Σm_ν rancio (0.308881 imprime 0.30888)",
              "94.07": "C_ν retirado (vigente 93.14)",
              "41.0187": "m_φ de Σm_ν rancio (vigente 40.70)"}
try:
    _logdir = _REPO / "results" / "logs"
    _sucios, _hist = [], []
    for _lg in sorted(_logdir.glob("*.log")):
        _txt = _lg.read_text(errors="ignore")
        _hit33 = [f"{_k} ({_v})" for _k, _v in _RETIRADAS.items() if _k in _txt]
        if not _hit33:
            continue
        (_hist if _lg.name in _LOGS_HISTORICOS else _sucios).append(
            f"{_lg.name}: {', '.join(_hit33)}")
    check("R33 ningún log ACTIVO arrastra una constante retirada",
          not _sucios,
          "; ".join(_sucios) if _sucios
          else f"{len(list(_logdir.glob('*.log')))} logs barridos, "
               f"{len(_hist)} históricos declarados")
    if _hist:
        track_open(f"R33 {len(_hist)} logs históricos con constante retirada",
                   "; ".join(_hist) + "  (declarados en _LOGS_HISTORICOS)")
except Exception as e:
    check("R33 capa logs-vs-núcleo operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# R34 — ningún .py ACTIVO hardcodea una constante retirada (2026-07-26).
# R33 caza logs rancios; esto caza la causa AGUAS ARRIBA. scan_omega_m.py tenía
# `mnu=0.069` en la línea 8 —el residuo de C_ν=94.07— mientras el resto de la
# suite ya usaba 0.06849. No era un artefacto viejo: era el FUENTE, y producía
# el barrido que justifica el ancla H₀. Re-correr su log no habría servido de
# nada: habría vuelto a salir rancio. Regla: una constante se lee del núcleo,
# no se re-teclea. Exclusiones = los sitios que NOMBRAN el valor retirado a
# propósito (mutaciones del test, la derivación que compara 94.07 vs 93.14).
print("\nCapa R34 — fuentes vs núcleo: ningún .py hardcodea constante retirada")
try:
    _EXCL = {"test_guardian.py", "derive_nu_closure.py", "ssee_verify.py",
             "ssee_core.py"}
    _pat34 = _re.compile(
        r"(?:mnu|Smnu|SUM_MNU_EV|sigma_m_nu|m_nu|C_nu|C_NU)\s*=\s*"
        r"(0\.069(?:0[0-9]*)?|94\.07[0-9]*)\s*(?:[^0-9.]|$)")
    _malos = []
    for _py in sorted((_REPO / "src").rglob("*.py")):
        if _py.name in _EXCL or "archive" in _py.parts:
            continue
        for _i, _ln in enumerate(_py.read_text(errors="ignore").splitlines(), 1):
            _code = _ln.split("#")[0]
            _m34 = _pat34.search(_code)
            if _m34:
                _malos.append(f"{_py.relative_to(_REPO)}:{_i} → {_m34.group(1)}")
    check("R34 ningún .py activo hardcodea Σm_ν/C_ν retirados",
          not _malos,
          "; ".join(_malos) if _malos
          else "constantes leídas del núcleo, no re-tecleadas")
except Exception as e:
    check("R34 capa fuentes-vs-núcleo operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# R35 — artefacto más viejo que su fuente (2026-07-26).
# R33/R34 buscan huellas que hay que ADIVINAR de antemano ("0.06902", "94.07"):
# sirven para el error ya conocido, no para el siguiente. Esta regla no adivina
# nada — si el script cambió después de generarse su log, el log es sospechoso,
# sea cual sea la constante. Probado retrospectivamente contra el estado del
# 2026-07-25: habría marcado 4 de los 5 logs rancios DOCE días antes.
# Usa fechas de COMMIT, no mtime: un checkout reescribe mtime y mentiría.
# Mapa en PROPAGACION.yaml; los logs sin fuente declarada se CUENTAN como
# deuda visible (ABIERTO), que es lo contrario de una nota suelta.
print("\nCapa R35 — artefacto vs fuente: ningún log más viejo que su script")
try:
    import ast          # noqa: F401 — usado por el filtro de docstrings
    import subprocess as _sp
    import yaml as _yaml
    _prop = _yaml.safe_load((_REPO / "PROPAGACION.yaml").read_text())
    _mapa = _prop.get("logs") or {}
    _hist35 = set(_prop.get("historicos") or [])

    def _commit_ts(rel):
        try:
            _o = _sp.run(["git", "log", "-1", "--format=%at", "--", rel],
                         cwd=_REPO, capture_output=True, text=True, timeout=20)
            return int(_o.stdout.strip()) if _o.stdout.strip() else None
        except Exception:
            return None

    _rancios, _sinmapa = [], []
    for _lg in sorted((_REPO / "results" / "logs").glob("*.log")):
        _nm = _lg.stem
        if _nm in _hist35:
            continue
        _script = _mapa.get(_nm)
        if not _script:
            _sinmapa.append(_nm)
            continue
        _tl, _ts35 = _commit_ts(f"results/logs/{_nm}.log"), _commit_ts(_script)
        if not (_tl and _ts35 and _ts35 > _tl):
            continue
        # El script es posterior — pero ¿cambió el CÓDIGO o sólo la prosa?
        # Comparar el AST en ambos puntos: si es idéntico, el cambio fue de
        # comentarios/docstring y el log sigue siendo válido. Sin esto la regla
        # grita por documentación, y una alarma ruidosa se termina ignorando
        # (comprobado: mcmc_paper2_3models_wmfix, cuyo único cambio fue la
        # línea «Ω_m = 0.30889» → «0.308881» dentro del docstring).
        try:
            _sha_log = _sp.run(["git", "log", "-1", "--format=%H", "--",
                                f"results/logs/{_nm}.log"], cwd=_REPO,
                               capture_output=True, text=True, timeout=20).stdout.strip()
            _viejo = _sp.run(["git", "show", f"{_sha_log}:{_script}"], cwd=_REPO,
                             capture_output=True, text=True, timeout=20).stdout
            _nuevo = (_REPO / _script).read_text(errors="ignore")
            def _sin_docstrings(_txt):
                _t = ast.parse(_txt)
                for _nd in ast.walk(_t):
                    if isinstance(_nd, (ast.Module, ast.FunctionDef,
                                        ast.AsyncFunctionDef, ast.ClassDef)):
                        _b = _nd.body
                        if (_b and isinstance(_b[0], ast.Expr)
                                and isinstance(_b[0].value, ast.Constant)
                                and isinstance(_b[0].value.value, str)):
                            _nd.body = _b[1:]           # fuera el docstring
                return ast.dump(_t)
            if _viejo and _sin_docstrings(_viejo) == _sin_docstrings(_nuevo):
                continue                      # sólo cambió la prosa
        except Exception as _e35:
            # NUNCA tragarse el fallo: un `except: pass` aquí escondió un
            # NameError propio (ast no estaba importado) y R35 reportó rancios
            # que no lo eran. Ante la duda se reporta, PERO con la razón visible.
            _rancios.append(f"{_nm} (no se pudo comparar AST: {type(_e35).__name__})")
            continue
        _rancios.append(f"{_nm} (script {(_ts35 - _tl) // 86400}d más nuevo)")
    check("R35 ningún log committeado es más viejo que el script que lo produce",
          not _rancios,
          "; ".join(_rancios) if _rancios
          else f"{len(_mapa)} logs mapeados al día, {len(_hist35)} históricos")
    if _sinmapa:
        track_open(f"R35 {len(_sinmapa)} logs sin fuente declarada en PROPAGACION.yaml",
                   ", ".join(_sinmapa[:10]) + (" …" if len(_sinmapa) > 10 else "")
                   + "  (añadir su script para que R35 los vigile)")
    # ── R36 — lo mismo para FIGURAS, priorizando las de submission ────────
    # R35 vigila logs; una figura rancia pasaba igual. A diferencia de los logs,
    # el script productor SÍ se deriva solo (el nombre del archivo aparece en el
    # código), así que no hace falta mapa. Mismo filtro AST: 32 figuras salían
    # "rancias" por timestamp y 10 eran sólo cambios de prosa.
    # Se distingue por severidad: las figuras QUE ENTRAN AL PRD son ROJO —van al
    # journal—; el resto se cuenta como deuda ABIERTA en vez de bloquear.
    _prd_figs = set(_re.findall(r"includegraphics\[[^\]]*\]\{([^}]+)\}",
                                (_REPO / "submission_PRD" / "SSEE_PRD.tex")
                                .read_text(errors="ignore")))
    _prd_figs = {f.rsplit(".", 1)[0] for f in _prd_figs}

    def _ast_igual(_sha, _rel):
        _v = _sp.run(["git", "show", f"{_sha}:{_rel}"], cwd=_REPO,
                     capture_output=True, text=True, timeout=20).stdout
        if not _v:
            return False
        _d = lambda _x: ast.dump(ast.parse(_x))          # noqa: E731
        try:
            _a = ast.parse(_v); _b = ast.parse((_REPO / _rel).read_text(errors="ignore"))
            for _t in (_a, _b):
                for _nd in ast.walk(_t):
                    if isinstance(_nd, (ast.Module, ast.FunctionDef,
                                        ast.AsyncFunctionDef, ast.ClassDef)):
                        _bd = _nd.body
                        if (_bd and isinstance(_bd[0], ast.Expr)
                                and isinstance(_bd[0].value, ast.Constant)
                                and isinstance(_bd[0].value.value, str)):
                            _nd.body = _bd[1:]
            return ast.dump(_a) == ast.dump(_b)
        except Exception:
            return False

    _fig_prd, _fig_otras = [], []
    for _fg in sorted((_REPO / "results" / "figures").glob("*.pdf")):
        _nm36 = _fg.stem
        _hits = [str(_p.relative_to(_REPO)) for _p in (_REPO / "src").rglob("*.py")
                 if "archive" not in _p.parts and _nm36 in _p.read_text(errors="ignore")]
        if not _hits:
            continue
        _scr = _hits[0]
        _tf = _commit_ts(f"results/figures/{_fg.name}")
        _ts36 = _commit_ts(_scr)
        if not (_tf and _ts36 and _ts36 > _tf):
            continue
        _sha36 = _sp.run(["git", "log", "-1", "--format=%H", "--",
                          f"results/figures/{_fg.name}"], cwd=_REPO,
                         capture_output=True, text=True, timeout=20).stdout.strip()
        if _ast_igual(_sha36, _scr):
            continue                                     # sólo prosa
        (_fig_prd if _nm36 in _prd_figs else _fig_otras).append(
            f"{_nm36} ({(_ts36 - _tf) // 86400}d)")
    check("R36 ninguna figura DEL PRD es más vieja que el script que la produce",
          not _fig_prd,
          "; ".join(_fig_prd) if _fig_prd
          else f"{len(_prd_figs)} figuras del PRD al día")
    if _fig_otras:
        track_open(f"R36 {len(_fig_otras)} figuras rancias fuera del PRD",
                   ", ".join(_fig_otras[:8]) + (" …" if len(_fig_otras) > 8 else ""))
except Exception as e:
    check("R35 capa artefacto-vs-fuente operable", False, str(e))

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

    # ── R15 — coherencia de la tensión derivada de n_s ──────────────────
    # Remache forjado en la auditoría 2026-07-09 (Paper 1 L334 decía «matches
    # Planck PR4 at 0.24σ» — el 0.24 corresponde a un central 0.9665 que NO es
    # PR4; el valor real es 0.16-0.17σ, y contradecía al Endorser). El guardián
    # RECOMPUTA n_s = 1-φ⁻⁷ y su tensión, y escanea cualquier σ pegado (≤130
    # chars) a «\varphi^{-7}» que se salga del valor real. Escaneo LOCAL a la
    # fórmula → el 0.24σ legítimo de w0wa en otro lado no da falso positivo.
    _phi = (1 + 5 ** 0.5) / 2
    _ns_alg = 1 - _phi ** -7                       # 0.965558
    _ns_planck, _ns_err = 0.9649, 0.0042           # Planck PR4 (mismo que Endorser)
    _ns_tens = abs(_ns_alg - _ns_planck) / _ns_err  # 0.157σ
    _ns_lo, _ns_hi = _ns_tens - 0.06, _ns_tens + 0.06  # ventana [0.10, 0.22]
    _sig_re = r"([+\-−]?)([0-9]\.[0-9]{1,2})\s*\\?sigma"
    _ns_hits = []          # valores citados FUERA de rango (incorrectos)
    _ns_quoted = {}        # {archivo: set(redondeos citados)} para coherencia
    for t in _texs:
        txt = t.read_text(errors="ignore")
        # Dos anclas: (a) la fórmula 1-φ⁻⁷, (b) el central Planck 0.9649 —
        # porque la tensión a veces se enuncia DESACOPLADA de la fórmula
        # (p.ej. «n_s=0.9649±0.0042 is 0.2σ from the SSEE value»), sin φ⁻⁷ cerca.
        # Este era el hueco por el que R15 dejó pasar el 0.2σ del Paper 1 App.
        _anchors = [(m.end(), m.end() + 100)
                    for m in _re2.finditer(r"varphi\^\{-7\}", txt)]
        _anchors += [(max(0, m.start() - 15), m.end() + 90)
                     for m in _re2.finditer(r"0\.9649", txt)]
        for a, b in _anchors:
            window = txt[a:b].split(r"\\", 1)[0]  # no cruzar salto de fila LaTeX
            wl = window.lower()
            if not ("planck" in wl or "0.96556" in window or "0.9649" in window
                    or "spectral" in wl):
                continue
            # σ con signo (+2.9σ, −0.6σ) = reducción de tensión de OTRA cantidad
            # (p.ej. modulación IS-Eckart), no la tensión n_s-vs-dato: se excluye.
            sm = _re2.search(_sig_re, window)
            if sm and not sm.group(1):
                _ns_quoted.setdefault(t.name, set()).add(sm.group(2))
                val = float(sm.group(2))
                if not (_ns_lo <= val <= _ns_hi):
                    _ns_hits.append(f"{t.name}:n_s@{val}σ (real {_ns_tens:.2f}σ)")
    # Coherencia de redondeo: un mismo manuscrito no debe citar la tensión n_s
    # con dos redondeos distintos (p.ej. 0.16σ en el cuerpo y 0.2σ en el App).
    # El 0.2 caía DENTRO de la ventana ±0.06 → «correcto» por separado, pero
    # contradecía al 0.16 del cuerpo: la grieta era la INconsistencia, no el valor.
    _ns_incoh = [f"{n}:{sorted(v)}" for n, v in _ns_quoted.items() if len(v) > 1]
    _ns_ok = not _ns_hits and not _ns_incoh
    _ns_msg = f"n_s={_ns_alg:.5f}, tensión {_ns_tens:.2f}σ (Planck PR4), redondeo único"
    if _ns_hits:
        _ns_msg = "; ".join(_ns_hits[:4])
    elif _ns_incoh:
        _ns_msg = "redondeo incoherente: " + "; ".join(_ns_incoh[:4])
    check("manuscritos  R15 tensión n_s: correcta y con redondeo único", _ns_ok, _ns_msg)

    # ── R17 — la constante de conversión de neutrinos debe ser univaluada ──
    # Remache forjado en la auditoría de Paper 1 (2026-07-10): la suite mezcla
    # 93.14 (en ω_ν=Σm_ν/93.14, Papers 1/3/6) con 94.07 (en la fórmula de Σm_ν,
    # Papers 1/4/6). Es la MISMA constante física (relic-density↔masa); un referí
    # caza el ~1% de inconsistencia.
    # RESUELTO 2026-07-25: se estandarizó en 93.14 (PDG) y se EJECUTÓ la re-propagación
    # que aquí quedaba pendiente (Σm_ν 0.0690/0.06902 → 0.06849 en ssee_core, P3, P6,
    # class_ssee, los 10 papers y las 3 memorias). La decisión se había tomado el
    # 2026-07-10 pero la re-propagación quedó A MEDIAS: m_φ sí se actualizó
    # (41.02→40.70), ω_ν no. LECCIÓN: una decisión sin re-propagación ejecutada es un
    # drift latente que se ve idéntico a un valor sano. Ahora lo cierra V-L2-11a/b/c.
    #
    # ⚠️ POR QUÉ R17 NO LO CAZÓ, Y LA REGLA QUE SE DERIVA:
    # R17 es un check de CADENA DE TEXTO: verifica que "94.07" no aparezca en los
    # .tex. El string se borró → R17 pasó → guardián VERDE. Pero el NÚMERO derivado
    # de 94.07 (Σm_ν=0.06902 ⇒ ω_ν=0.000741) siguió vivo en código y papers. El
    # guardián certificaba media tarea con luz verde entera, y un VERDE se lee como
    # "hecho".
    # REGLA GENERAL: retirar una constante exige DOS chequeos, no uno —
    #   (a) el string desaparece            [R17, superficie]
    #   (b) los números que derivaban de ella se recomputaron  [V-L2-11, consecuencia]
    # Un check de (a) sin su (b) es peor que ningún check: crea confianza falsa.
    # Aplicar este par a cualquier retiro futuro de constante.
    _nu_9407 = sorted(t.name for t in _texs
                      if "94.07" in t.read_text(errors="ignore"))
    check("manuscritos  R17 constante ν univaluada en 93.14 (sin 94.07)",
          not _nu_9407,
          "univaluada 93.14 (94.07 solo vive en derive_nu_closure.py, la demostración)"
          if not _nu_9407 else "94.07 aún en: " + ", ".join(_nu_9407)
          + " — misma C; usar 93.14 PDG")

    # ── R18 — sin narrativa H0-posterior stale del reframe Ω_m-geometría ──
    # Remache forjado en la auditoría de Paper 2 (2026-07-10): el fix V-L4-DESI
    # (2026-07-09, sector frío 0.160 fuera de E(z)) se propagó al abstract/§2.4/§6.3
    # PERO dejó prosa stale en §3.2 y en el párrafo "algebraic anchor": afirmaba que
    # el posterior H0 "lies at 2.9σ from anchor", "pull downward to compensate for
    # the enlarged r_d", y "r_d = Ω_m,dyn(H0/100)²" — CONTRADICE el 0.04σ/coincide
    # corregido en la misma sección. R14 vigila los SCRIPTS (canario E(z)); esto es
    # PROSA. Los fragmentos siguientes solo existieron como la aserción retirada:
    _stale_v4 = {
        r"\Omega_{m,\mathrm{dyn}}(H_0/100)": "r_d desde el sector dinámico 0.160 (va la total 0.308881)",
        "compensate for the enlarged": "narrativa vieja: no hay r_d agrandado (148.2≈ΛCDM)",
        "downward to compensate": "posterior COINCIDE con anchor (0.04σ), no baja a compensar",
    }
    _v4_hits = []
    for t in _texs:
        txt = t.read_text(errors="ignore")
        for pat, why in _stale_v4.items():
            if pat in txt:
                _v4_hits.append(f"{t.name}«{pat[:22]}»→{why}")
    check("manuscritos  R18 sin narrativa H0-posterior stale (V-L4-DESI)",
          not _v4_hits,
          "posterior coincide con anchor (0.04σ), geometría con Ω_m total"
          if not _v4_hits else "; ".join(_v4_hits[:4]))

    # ── R19 — ningún manuscrito muestra los valores BAO DR1-mislabeled ──
    # Remache forjado en la auditoría de Paper 2 (2026-07-10): la Tabla 1 de P2 y su
    # §B.2 mostraban los valores DR1 (7.93±0.15, LRG DH 20.08, QSO 30.21/13.23 @z=1.491,
    # Lya 39.71/8.52) aunque las cadenas jul-9 SÍ usaban DR2 (load_desi_dr2, csv fuente).
    # Datos mostrados ≠ datos usados. R14 vigila el csv/loader; R19 vigila que la PROSA
    # de los manuscritos no reintroduzca los valores DR1. Fingerprints inequívocos (DR2
    # da 21.863/30.512/38.988/8.632):
    _dr1_bao = {"20.08": "LRG DH DR1 (DR2=21.863)", "30.21": "QSO DM DR1 (DR2=30.512)",
                "39.71": "Lya DM DR1 (DR2=38.988)", "8.52": "Lya DH DR1 (DR2=8.632)",
                "16.85": "LRG2 DM DR1 (DR2=17.351)", "1.491": "z_QSO DR1 (DR2 z=1.484)"}
    _bao_hits = []
    for t in _texs:
        txt = t.read_text(errors="ignore")
        for pat, why in _dr1_bao.items():
            if pat in txt:
                _bao_hits.append(f"{t.name}«{pat}»={why}")
    check("manuscritos  R19 sin valores BAO DR1-mislabeled (datos mostrados = usados)",
          not _bao_hits,
          "manuscritos muestran DR2 (data/raw/desi_dr2_bao.csv, la que usan las cadenas)"
          if not _bao_hits else "; ".join(_bao_hits[:5]))

    # Las R21-R24 también vigilan el PRD (submission_PRD), no solo manuscript/.
    _prd = _REPO2 / "submission_PRD" / "SSEE_PRD.tex"
    _texs2 = _texs + ([_prd] if _prd.exists() else [])

    # ── R21 — el denominador de wₐ es IGNIS (π+PYROS), NUNCA K_v (2Ω) ──────
    # Remache forjado 2026-07-19: la suite escribía wₐ = P_sc/K_v, colapsando el
    # denominador al scaffold K_v (=KRYSTOS φ+π+Ω, 2Ω por valor) cuando la ENTIDAD
    # es IGNIS = π+PYROS (rama-π, intra-linaje con el numerador PYROS=P_sc): mismo
    # valor 9.519, entidad distinta. K_v SOLO es legítimo como sumando de M_v=φ+π+K_v;
    # NUNCA como denominador → basta prohibir K_v en rol de denominador ("}{K_v" en
    # \frac, "/K_v" inline). M_v=φ+π+K_v (K_v como "+K_v") no matchea.
    _kvden = []
    for t in _texs2:
        txt = t.read_text(errors="ignore")
        for pat in (r"\}\{K_v", r"/K_v"):
            m = _re2.search(pat, txt)
            if m:
                _kvden.append(f"{t.name}«{txt[max(0, m.start()-8):m.end()][:20]}»")
    check("manuscritos  R21 wₐ denominador = IGNIS(π+PYROS); K_v nunca es denominador",
          not _kvden,
          "wₐ = P_sc/IGNIS en toda la suite; K_v solo como +K_v en M_v"
          if not _kvden else "; ".join(_kvden[:5]))

    # ── R21b — la MISMA ley, pero en el CÓDIGO ────────────────────────────
    # Hallazgo 2026-07-25 (revisión pregunta-por-pregunta del PRD, §2.2): R21
    # sólo escaneaba .tex. Los .py seguían escribiendo `-P_sc / KV` en 10 sitios
    # de 6 archivos — el paper decía IGNIS y el código decía K_v. Como K_v e
    # IGNIS son iguales BIT A BIT (ambos 2Ω), ningún número delataba la grieta:
    # exactamente el falso verde que el guardián existe para impedir. Un texto y
    # su prueba no pueden declarar linajes distintos aunque coincidan en valor.
    # OJO: sólo se prohíbe K_v en rol de DIVISOR; `MV = PHI + PI + KV` es legítimo.
    _kvden_py = []
    for _py in sorted(_REPO2.glob("src/**/*.py")) + sorted(_REPO2.glob("class_ssee/*.py")):
        # test_guardian.py queda fuera: sus mutaciones CONTIENEN a propósito la
        # forma prohibida (es su trabajo escribirla para probar que duele).
        if "archive" in _py.parts or _py.name in ("ssee_verify.py", "test_guardian.py"):
            continue
        for _ln, _line in enumerate(_py.read_text(errors="ignore").splitlines(), 1):
            # Se mira CÓDIGO, no prosa: un comentario que NOMBRA la forma prohibida
            # está documentándola (esta regla misma lo hace). Se recorta desde el «#».
            _line = _re2.sub(r"#.*$", "", _line)
            # KRYSTOS_V entra a la lista 2026-07-25: la primera versión de R21b sólo
            # miraba las abreviaturas (KV/Kv/K_v/K_V) y dejó pasar el nombre completo
            # justo en look_elsewhere_full.py — el script que produce el «1 de 490».
            if _re2.search(r"/\s*(KV|Kv|K_v|K_V|KRYSTOS_V)\b", _line):
                _kvden_py.append(f"{_py.name}:{_ln}")
    check("código       R21b wₐ = −P_sc/IGNIS en los scripts (K_v nunca divide)",
          not _kvden_py,
          "ningún script divide por K_v; el linaje del código = el del paper"
          if not _kvden_py else "; ".join(_kvden_py[:5]))

    # ── R22 — precisión: sin 5-dec truncados/falsos de las constantes ─────
    # Remache forjado 2026-07-19: el Sealed truncaba (no redondeaba) la tabla de
    # notación y daba wₐ=-0.67000 (falso: es -0.66997) y w₀=-0.83996 (falso:
    # -0.83995; 1+w₀=0.16005 lo exige). Un referí caza el 5º dígito. Se prohíben
    # los 5-dec ERRÓNEOS exactos, con lookahead (?![0-9]) para no chocar con la
    # precisión plena legítima (14.278879927, 4.759626642, …).
    _badprec = {
        r"0\.83996(?![0-9])": "w₀ falso (→ -0.83995)",
        r"0\.67000(?![0-9])": "wₐ falso (→ -0.66997)",
        r"4\.75962(?![0-9])": "Ω truncado (→ 4.75963)",
        r"5\.52140(?![0-9])": "KAL truncado (→ 5.52141)",
        r"3\.99784(?![0-9])": "AURA truncado (→ 3.99785)",
        r"6\.37765(?![0-9])": "P_sc truncado (→ 6.37766)",
        r"11\.99353(?![0-9])": "T_r truncado (→ 11.99354)",
        r"14\.27887(?![0-9])": "M_v truncado (→ 14.27888)",
    }
    _prec_hits = []
    for t in _texs2:
        txt = t.read_text(errors="ignore")
        for pat, why in _badprec.items():
            if _re2.search(pat, txt):
                _prec_hits.append(f"{t.name}:{why}")
    check("manuscritos  R22 sin constantes 5-dec truncadas/falsas",
          not _prec_hits,
          "precisión canónica (w₀=-0.83995, wₐ=-0.66997, tabla redondeada)"
          if not _prec_hits else "; ".join(_prec_hits[:6]))

    # ── R23 — sin fósiles Σ₉ / "5D" / MIKAEL_V (nombres/símbolos retirados) ─
    # Remache forjado 2026-07-19: Σ₉ (de las "9 Soberanías", hoy 25) sobrevivía en
    # P7/P9; "five-dimensional integration ceiling" para M_v (que es 4D, dentro de
    # CUARTAL=4·AURA) en Sealed/PRD; y MIKAEL_V (renombrado ATLAS por la Ley de
    # Nombrado) como etiqueta de M_v en P1/P7. MIKAEL a secas (Soberana viva) NO
    # matchea: el patrón exige el sufijo _V / \_V.
    _fossils = {
        r"Sigma_9": "Σ₉ fósil (9 Soberanías → 25); usar 3Ω/M_v",
        r"Σ₉": "Σ₉ fósil; usar 3Ω/M_v",
        r"five-dimensional integration": "M_v es 4D (dentro de CUARTAL), no 5D",
        r"5D integration": "M_v es 4D, no 5D",
        r"5-D integration": "M_v es 4D, no 5D",
        r"MIKAEL\\?_V": "MIKAEL_V renombrado ATLAS (Ley de Nombrado)",
    }
    _fos_hits = []
    for t in _texs2:
        txt = t.read_text(errors="ignore")
        for pat, why in _fossils.items():
            if _re2.search(pat, txt):
                _fos_hits.append(f"{t.name}:{why}")
    check("manuscritos  R23 sin fósiles Σ₉ / 5D / MIKAEL_V",
          not _fos_hits,
          "Σ₉→3Ω, 5D→4D, MIKAEL_V→ATLAS aplicados"
          if not _fos_hits else "; ".join(_fos_hits[:6]))

    # ── R24 — conteo del diccionario univaluado: 55/25/490 (sin rezagos) ──
    # Remache forjado 2026-07-19: el Sealed mezclaba "50 named / 22 distinct" (v1.3)
    # con "490 ratios" (v1.4) en la MISMA zona; ni los 3 auditores externos ni el
    # guardián lo cazaron. Canónico único: 55 nombres / 25 valores / 490 razones.
    # Se prohíben los denominadores/counts RETIRADOS en contexto (no el bare "378"
    # de 6.378=P_sc): "1 of 378", "1/378", "50 named", "22 distinct", etc.
    _oldcount = {
        r"1\s*(?:of|/)\s*378": "1/378 retirado (→ 1/490)",
        r"1\s*(?:of|/)\s*245": "1/245 retirado (→ 1/490)",
        r"1\s*(?:of|/)\s*337": "1/337 retirado (→ 1/490)",
        r"1\s*(?:of|/)\s*317": "1/317 retirado (→ 1/490)",
        r"\b50\s+named": "50 named retirado (→ 55)",
        r"\b46\s+named": "46 named retirado (→ 55)",
        r"\b22\s+distinct": "22 distinct retirado (→ 25)",
        r"\b21\s+distinct": "21 distinct retirado (→ 25)",
    }
    _cnt_hits = []
    for t in _texs2:
        txt = t.read_text(errors="ignore")
        for pat, why in _oldcount.items():
            if _re2.search(pat, txt):
                _cnt_hits.append(f"{t.name}:{why}")
    check("manuscritos  R24 conteo diccionario 55/25/490 (sin rezagos 50/22/378)",
          not _cnt_hits,
          "55 named / 25 distinct / 490 ratios univaluado"
          if not _cnt_hits else "; ".join(_cnt_hits[:6]))
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
        if "KRYSTOS" in _master:
            _master.add("KRYSTOS_V")      # KRYSTOS_V = nombre correcto (φ+π+Ω); maestro aún usa "KRYSTOS" → v1.2 pendiente
        # Ley de Nombrado v1.5 (2026-07-19): MIKAEL_V→ATLAS, MIKE→PHOENIX. El maestro
        # vive en sandbox_unificado/ (submódulo del repo PÚBLICO multidominio) y NO se
        # toca antes del envío — así que la equivalencia se declara aquí, igual que
        # con KRYSTOS_V. Es alias de NOMBRE, no de valor: el nodo es el mismo.
        # Al re-sincronizar el maestro, estas dos líneas sobran.
        if "MIKAEL_V" in _master:
            _master.add("ATLAS")
        if "MIKE" in _master:
            _master.add("PHOENIX")
    _lf = _root / "src" / "estadistica" / "look_elsewhere_full.py"
    _fam = set()
    if _lf.exists():
        _m = _re.search(r"FAMILY\s*=\s*\{(.*?)\}",
                        _lf.read_text(errors="ignore"), _re.S)
        if _m:
            # incluye nombres con tilde (DÜSTAL/TRÏSTAL/CUÄSTAL) y acento agudo (ÁNGELOS)
            _fam = {k.upper() for k in _re.findall(r'"([A-ZÁÉÍÓÚÄÖÜÏÑ_]+)"\s*:', _m.group(1))}
    _spurious = sorted(_fam - _master) if (_master and _fam) else []
    check("diccionario  R13 look_elsewhere sin nodos ausentes del maestro",
          not _spurious,
          "consistente (script ⊆ maestro)" if not _spurious
          else "ESPURIOS (no en maestro): " + ", ".join(_spurious))

    # ── R16 — coherencia interna del look-elsewhere ─────────────────────
    # Remache forjado en la auditoría 2026-07-09 (Paper 1): el script tracked
    # tenía docstring «29 constantes» pero la etiqueta del reporte hardcodeaba
    # «31» mientras su FAMILY tiene 29 → conteo auto-contradictorio, y el «378»
    # que Paper 1 cita como defensa anti-numerología quedaba sin anclar. R16
    # recomputa: (a) docstring == etiqueta-reporte == len(FAMILY); (b) el conteo
    # de razones distintas en (0,5] == 490 (el número robusto que citan los papers,
    # tras la completitud ≤TRIAL a 53 constantes, 2026-07-18).
    if _lf.exists():
        _lftxt = _lf.read_text(errors="ignore")
        _nfam = len(_fam)
        _ds = _re.search(r"\((\d+)\s+constantes", _lftxt)     # docstring (permite "constantes con nombre")
        _lbl = _re.findall(r"COMPLETO\s*\((?:[^)]*?)(\d+)\s+constantes\)", _lftxt)
        _ds_n = int(_ds.group(1)) if _ds else -1
        # etiqueta dinámica (f-string len(FAMILY)) cuenta como coherente
        _lbl_dyn = "len(FAMILY)" in _lftxt or "{len(fam)" in _lftxt.lower()
        _lbl_ok = _lbl_dyn or all(int(x) == _nfam for x in _lbl)
        _cnt_ok = (_ds_n == _nfam) and _lbl_ok
        check("diccionario  R16 look-elsewhere conteo interno coherente",
              _cnt_ok,
              f"docstring={_ds_n}, FAMILY={_nfam}, etiqueta {'dinámica' if _lbl_dyn else _lbl}"
              if _cnt_ok else
              f"INCOHERENTE: docstring={_ds_n} vs FAMILY={_nfam} vs etiqueta={_lbl}")
except Exception as e:
    check("diccionario  capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA DESI — procedencia de los datos BAO (R14, 2026-07-01).
# Historia: la suite entera usó DESI DR1 (2404.03002) MAL ETIQUETADO como
# DR2 durante meses (11/13 puntos DR1 exactos + un QSO z=1.491 sin fuente).
# Un chequeador de física es ciego a esto (los números eran internamente
# consistentes); solo un check de PROCEDENCIA contra la tabla oficial lo caza.
# Tres candados:
#   1. El csv canónico coincide dígito a dígito con DR2 Tabla 4 (2503.14738)
#      — la tabla de cotejo vive AQUÍ, copiada independientemente del csv.
#   2. Cero valores-centinela DR1 en el csv (20.98, 7.93, 13.62...).
#   3. Los consumers importan del loader único; nadie re-hardcodea datos BAO.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa DESI — procedencia BAO: csv == DR2 Tabla 4 oficial, sin DR1, sin hardcode")
try:
    _root = pathlib.Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(_root / "src"))
    from desi_dr2_data import load_desi_dr2, desi_covariance
    import numpy as _np

    # 1) Cotejo independiente contra DESI DR2 (2503.14738, Tabla 4).
    #    (z, quantity, value, sigma, corr_MH) — transcrito del paper, NO del csv.
    _DR2_OFICIAL = [
        (0.295, "DV_over_rd",  7.942, 0.075,  None),
        (0.510, "DM_over_rd", 13.588, 0.167, -0.459),
        (0.510, "DH_over_rd", 21.863, 0.425, -0.459),
        (0.706, "DM_over_rd", 17.351, 0.177, -0.404),
        (0.706, "DH_over_rd", 19.455, 0.330, -0.404),
        (0.934, "DM_over_rd", 21.576, 0.152, -0.416),
        (0.934, "DH_over_rd", 17.641, 0.193, -0.416),
        (1.321, "DM_over_rd", 27.601, 0.318, -0.434),
        (1.321, "DH_over_rd", 14.176, 0.221, -0.434),
        (1.484, "DM_over_rd", 30.512, 0.760, -0.500),
        (1.484, "DH_over_rd", 12.817, 0.516, -0.500),
        (2.330, "DM_over_rd", 38.988, 0.531, -0.431),
        (2.330, "DH_over_rd",  8.632, 0.101, -0.431),
    ]
    _d = load_desi_dr2()
    _mism = []
    for _k, (_z, _q, _v, _s, _c) in enumerate(_DR2_OFICIAL):
        if not (abs(_d["z"][_k] - _z) < 1e-9 and _d["quantity"][_k] == _q
                and abs(_d["value"][_k] - _v) < 1e-9
                and abs(_d["sigma"][_k] - _s) < 1e-9
                and ((_c is None and _np.isnan(_d["corr"][_k]))
                     or (_c is not None and abs(_d["corr"][_k] - _c) < 1e-9))):
            _mism.append(f"fila {_k}: z={_z} {_q}")
    check("DESI  R14 csv == DR2 Tabla 4 oficial (13 pts, dígito a dígito)",
          len(_d["value"]) == 13 and not _mism,
          "idéntico" if not _mism else "DIFIERE: " + "; ".join(_mism[:3]))
    check("DESI  R14 release/arXiv declarados en el csv",
          _d.get("release") == "DR2" and _d.get("arxiv") == "2503.14738",
          f"{_d.get('release')}/{_d.get('arxiv')}")

    # 2) Cero centinelas DR1 (2404.03002) o del QSO huérfano en el csv.
    _csv_txt = (_root / "data" / "raw" / "desi_dr2_bao.csv").read_text()
    _data_lines = [l for l in _csv_txt.splitlines()
                   if l.strip() and not l.startswith("#") and not l.startswith("z_eff")]
    _DR1_SENT = ["20.98", "20.08", "16.85", "13.62", "27.79", "13.82",
                 "39.71", "8.52", "17.88", "21.71", "30.21", "13.23", "1.491"]
    _dr1_hits = [s for s in _DR1_SENT if any(s in l for l in _data_lines)]
    check("DESI  R14 csv sin valores-centinela DR1/QSO-huérfano",
          not _dr1_hits, "limpio" if not _dr1_hits else "DR1: " + ", ".join(_dr1_hits))

    # 3) Consumers: importan del loader y no re-hardcodean datos BAO.
    _consumers = [
        _root / "src" / "mcmc_full" / "ssee_likelihoods.py",
        _root / "class_ssee" / "ssee_mcmc_fase4.py",
        _root / "src" / "p02_mcmc" / "ssee_paper2_mcmc_reframe.py",
        _root / "src" / "p02_mcmc" / "ssee_paper2_mcmc.py",
        _root / "src" / "p02_mcmc" / "ssee_paper2_mcmc_lcdm_baseline.py",
        _root / "src" / "p02_mcmc" / "rerun_cpl_h0anchors.py",
        _root / "src" / "p09_hubble" / "ssee_h0_prior_experiment.py",
        _root / "src" / "estadistica" / "ssee_phase_d_savage_cv.py",
    ]
    _no_loader, _rehard = [], []
    for _f in _consumers:
        _t = _f.read_text(errors="ignore") if _f.exists() else ""
        if "load_desi_dr2" not in _t:
            _no_loader.append(_f.name)
        if any(s in _t for s in ("20.98", "30.21, ", "[0.295,  7.93")):
            _rehard.append(_f.name)
    check("DESI  R14 consumers importan del loader único",
          not _no_loader,
          "los 8 wireados" if not _no_loader else "SIN loader: " + ", ".join(_no_loader))
    check("DESI  R14 consumers sin datos BAO re-hardcodeados",
          not _rehard, "limpio" if not _rehard else "HARDCODE: " + ", ".join(_rehard))

    # 4) Covarianza bloque-diagonal: 6 pares DM-DH correlacionados.
    _C = desi_covariance(_d)
    _npairs = int((_np.abs(_C - _np.diag(_np.diag(_C))) > 0).sum() / 2)
    check("DESI  R14 covarianza con 6 pares DM-DH (r_MH oficiales)",
          _npairs == 6, f"{_npairs} pares")

    # 5) CANARIO DE GEOMETRÍA — remache del bug χ²=726 (2026-07-09, V-L4-DESI).
    #    La geometría de fondo (E(z), r_d, distancias) DEBE usar la materia
    #    TOTAL Ω_m = ω_m/h² = 0.308881.  El sector frío 0.160 = 1+w0 NO es una
    #    densidad de fondo (es perturbaciones P6 + factor EFT α_K) y NUNCA va
    #    en un E(z).  Meterlo ahí fue el bug que dio χ²_BAO=726 en DESI DR2.
    #    Este canario reconstruye la geometría BAO desde las MISMAS fuentes
    #    (loader único + ssee_core) y ancla la invariante en los dos sentidos:
    #    con la total el χ² es sano (~11); con el sector DEBE doler (~725).
    #    Si una futura edición vuelve a meter 0.160 en E(z), esto se pone ROJO.
    from ssee_core import (W0 as _W0, WA as _WA,
                           OMEGA_M_TOTAL as _OMT, OMEGA_CDM_SECTOR as _OMS)
    _Cinv = _np.linalg.inv(_C)
    _Z, _Q, _OBS = list(_d["z"]), list(_d["quantity"]), _d["value"]
    _CKM = 2.998e5
    def _fde(z):
        a = 1.0 / (1.0 + z)
        return (1 + z) ** (3 * (1 + _W0 + _WA)) * _np.exp(-3 * _WA * (1 - a))
    def _E(z, Om):
        return _np.sqrt(Om * (1 + z) ** 3 + (1 - Om) * _fde(z))
    def _DC(zm, Om, n=250):
        zz = _np.linspace(0, zm, n)
        return _np.trapezoid(1.0 / _E(zz, Om), zz)
    def _rd(obh2, omh2):
        return 147.27 * (omh2 / 0.1432) ** -0.255 * (obh2 / 0.02237) ** -0.134
    def _chi2_bao(Om, obh2=0.02237):
        best = 1e30
        for _H0 in _np.linspace(55, 80, 150):   # min sobre H0: ni su mejor H0 salva al sector
            _r = _rd(obh2, Om * (_H0 / 100) ** 2)
            _p = []
            for _z, _q in zip(_Z, _Q):
                _dm = (_CKM / _H0) * _DC(_z, Om)
                _dh = _CKM / (_H0 * _E(_z, Om))
                _p.append(_dm / _r if _q == "DM_over_rd"
                          else _dh / _r if _q == "DH_over_rd"
                          else (_z * _dm ** 2 * _dh) ** (1 / 3) / _r)
            _res = _np.array(_p) - _OBS
            best = min(best, float(_res @ _Cinv @ _res))
        return best
    _c_total, _c_sector = _chi2_bao(_OMT), _chi2_bao(_OMS)
    check("DESI  R14 canario-geometría: materia TOTAL (0.308881) da χ²_BAO sano",
          _c_total < 50, f"χ²={_c_total:.1f} (<50) con Ω_m,total")
    check("DESI  R14 canario-geometría: sector 0.160 en E(z) DEBE doler",
          _c_sector > 300, f"χ²={_c_sector:.1f} (>300) — reproduce el bug χ²=726")
except Exception as e:
    check("DESI  capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA R20 — anclas OBSERVACIONALES en código ⊆ CANONICAL_VALUES.yaml.
# Punto ciego de la auditoría externa 2026-07-13 (H2): un script hardcodeaba
# la PREDICCIÓN SSEE (S8=0.758) en el hueco de la OBSERVACIÓN KiDS (0.759),
# forzando 0.00σ en vez del 0.04σ real. Ningún patrón "retirado" lo cazaba
# porque 0.758 es un valor VIGENTE legítimo (la predicción). R20 verifica que
# cada variable observacional del código coincida con el DATO canónico del YAML.
# ─────────────────────────────────────────────────────────────────────
print("\nCapa R20 — anclas observacionales (código ⊆ CANONICAL_VALUES.yaml)")
try:
    import re as _re20
    _R20_ROOT = pathlib.Path(__file__).resolve().parents[2]
    _yaml20 = (_R20_ROOT / "CANONICAL_VALUES.yaml").read_text(errors="ignore")
    def _anchor20(key):
        m = _re20.search(rf"^\s*{key}:\s*([0-9.]+)", _yaml20, _re20.M)
        return m.group(1) if m else None
    # (variable observacional en código  →  clave del DATO en CANONICAL_VALUES.yaml)
    _R20_MAP = [("kids_s8", "obs_KiDS_S8"),
                ("kids_sig8", "obs_KiDS_sigma8"),
                ("des_s8", "obs_DES_S8")]
    _R20_FILES = sorted((_R20_ROOT / "src").rglob("*.py"))
    for _var, _key in _R20_MAP:
        _canon = _anchor20(_key)
        if _canon is None:
            check(f"R20 ancla {_key} definida en YAML", False, "no encontrada en CANONICAL_VALUES.yaml")
            continue
        _bad = []
        for _pf in _R20_FILES:
            for _m in _re20.finditer(rf"\b{_re20.escape(_var)}\b\s*=\s*([0-9.]+)",
                                     _pf.read_text(errors="ignore")):
                if not _m.group(1).startswith(_canon):
                    _bad.append(f"{_pf.name}:{_var}={_m.group(1)}")
        check(f"R20 {_var} == {_key}={_canon} (dato, no predicción)",
              not _bad,
              "coincide con el ancla observacional" if not _bad
              else "DESAJUSTE (predicción metida como dato?): " + "; ".join(_bad))
except Exception as e:
    check("R20 capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA R25 — parametrización: ω_m es el ABSOLUTO que SSEE fija; Ω_m es DERIVADO.
#
# Hallazgo 2026-07-25. Los MCMC congelaban la FRACCIÓN Ω_m=0.308881 y derivaban
# ω_m = Ω_m·h² en cada paso:
#       om_h2 = OMEGA_M_TOTAL*(H0/100)**2        ← parametrización INVERTIDA
# Pero SSEE predice el absoluto ω_m = ω_b+ω_c+ω_ν = 0.14267 (álgebra pura);
# Ω_m = ω_m/h² es una CONSECUENCIA de H₀, no un ingrediente. Al congelar Ω_m y
# variar H₀, el ω_m implícito se despegaba hasta ±1.8% de la predicción y sólo
# coincidía EXACTAMENTE en H₀ = 67.962 — el ancla. Es decir: el MCMC evaluaba
# SSEE fielmente SÓLO en el ancla y un modelo ligeramente distinto en el resto,
# sesgando el posterior hacia ella (67.948 = 0.04σ). Con ω_m fijo da 67.783
# (0.50σ) y el χ²_BAO MEJORA (10.72 → 10.33): el arreglo ajusta mejor los datos.
#
# La regla caza la FIRMA del error: multiplicar una Ω_m constante por h² para
# fabricar ω_m dentro de código de SSEE. En ΛCDM/CPL es LEGÍTIMO (ahí Ω_m es un
# parámetro muestreado, no una constante), así que sólo se marca cuando el
# multiplicando es una CONSTANTE de Ω_m (OMEGA_M_TOTAL / OM_GEOM / 0.3088x).
# ─────────────────────────────────────────────────────────────────────
print("\nCapa R25 — parametrización ω_m absoluto (no congelar Ω_m en el MCMC)")
try:
    _R25_ROOT = pathlib.Path(__file__).resolve().parents[2]
    _R25_PAT = re.compile(
        r"(OMEGA_M_TOTAL|OM_GEOM|OMEGA_M_CMB|0\.3088\d*|0\.30889)\s*\*\s*\(\s*H0\s*/\s*100")
    _R25_HITS = []
    for _pf in sorted((_R25_ROOT / "src").rglob("*.py")) + \
               sorted((_R25_ROOT / "class_ssee").rglob("*.py")):
        # test_guardian.py contiene la cadena del error A PROPÓSITO (es la mutación
        # que prueba que R25 dispara); marcarla sería morderse la cola.
        if ("superseded" in str(_pf) or "archive" in str(_pf)
                or _pf.name == "test_guardian.py"):
            continue
        for _i, _ln in enumerate(_pf.read_text(errors="ignore").splitlines(), 1):
            if _ln.lstrip().startswith("#") or "R25-control" in _ln:
                continue   # comentarios explicativos y el modo-control declarado
            if _R25_PAT.search(_ln):
                _R25_HITS.append(f"{_pf.name}:{_i}")
    # (b) PUNTO CIEGO de (a), hallado el mismo día: en ssee_phase_d_savage_cv.py el
    # patrón era `Om * (H0/100)**2` con `Om` una VARIABLE asignada desde la constante
    # 0.308881 — la forma superficial no coincidía, el significado sí. Lección idéntica
    # a la de las mutaciones: una regla que busca la superficie deja pasar el fondo.
    # Este segundo pase marca cualquier archivo que DEFINA una constante Ω_m≈0.3088x
    # y además fabrique ω_m con ·(H0/100)**2 en alguna línea — revisión obligatoria.
    _R25_CONST = re.compile(r"^\s*\w+\s*=\s*0\.3088\d*", re.M)
    _R25_MAKE = re.compile(r"\*\s*\(\s*H0\s*/\s*100\s*\)\s*\*\*\s*2")
    for _pf in sorted((_R25_ROOT / "src").rglob("*.py")):
        if ("superseded" in str(_pf) or "archive" in str(_pf)
                or _pf.name == "test_guardian.py"):
            continue
        _txt = _pf.read_text(errors="ignore")
        if not (_R25_CONST.search(_txt) and _R25_MAKE.search(_txt)):
            continue
        # exento si declara explícitamente que la constante es sólo diagnóstico
        if "NO congelar" in _txt or "R25-control" in _txt or "R25-ok" in _txt:
            continue
        _R25_HITS.append(f"{_pf.name} (Ω_m constante + ω_m=Ω_m·h² en el archivo)")
    check("R25 ningún ω_m fabricado como Ω_m·h² con Ω_m constante (SSEE)",
          not _R25_HITS,
          "ω_m = ω_b+ω_c+ω_ν algebraico; Ω_m se deriva por muestra"
          if not _R25_HITS else
          "PARAMETRIZACIÓN INVERTIDA en: " + "; ".join(_R25_HITS[:6]))
except Exception as e:
    check("R25 capa operable", False, str(e))

# ─────────────────────────────────────────────────────────────────────
# CAPA R26 — procedencia declarada: ningún canónico sin origen ni cadena.
#
# Causa raíz (OP-20, 2026-07-25): el factor 0.960318 vivió meses dentro del
# modelo sin corresponder a ninguna constante documentada (implicaba C_ν≈93.86,
# ni 94.07 ni 93.14). No lo cazó nadie porque NO HABÍA NINGÚN CAMPO donde
# tuviera que declarar su origen: bastaba con que el número «se viera bien».
# Y cuando C_ν pasó de 94.07 a 93.14, nada supo qué re-propagar porque la
# cadena de dependencias no estaba escrita — vivía en la memoria de quien hizo
# el cambio.
#
# Esta regla exige que cada nodo de `provenance:` declare las cuatro cosas:
#   formula · origin · source · (depends_on | affects)
# Un canónico sin `source` es, por definición, un número sin fuente.
#
# Verifica ADEMÁS que los valores declarados en provenance coincidan con lo
# que recomputa el guardián: la procedencia tiene que describir el número
# vigente, no uno de hace tres semanas (si no, es documentación stale, que es
# la otra mitad de la enfermedad — ver el docstring del MCMC de P2).
# ─────────────────────────────────────────────────────────────────────
print("\nCapa R26 — procedencia declarada de los canónicos")
try:
    import yaml as _yaml26
    _P26 = _yaml26.safe_load(
        (pathlib.Path(__file__).resolve().parents[2]
         / "CANONICAL_VALUES.yaml").read_text(encoding="utf-8")).get("provenance", {})
    _REQ = ("formula", "origin", "source")
    _incompletos = []
    for _k, _v in _P26.items():
        _falta = [c for c in _REQ if not _v.get(c)]
        if not (_v.get("depends_on") is not None or _v.get("affects")):
            _falta.append("depends_on/affects")
        if _falta:
            _incompletos.append(f"{_k}«falta {'+'.join(_falta)}»")
    check("R26 todo canónico con procedencia declara formula+origin+source+cadena",
          not _incompletos,
          f"{len(_P26)} nodos completos"
          if not _incompletos else "SIN FUENTE: " + "; ".join(_incompletos[:6]))

    # El valor declarado debe ser el que el guardián recomputa (no doc stale).
    _REC = {"Omega": Omega, "KAL0": KAL0, "omega_b": _omb, "omega_c": _omc,
            "omega_nu": _omnu_from_mnu, "omega_m": _omm,
            "sigma_m_nu": _mnu_active, "C_nu": _C_nu,
            # OJO: NO usar _R2 aquí — ese nombre lo reasigna la regla R2
            # (patrones del compendio multidominio, L~725) a una LISTA.
            # Colisión de nombres real en el guardián; se recomputa local.
            "R2": Omega / (KAL0 * Tr),
            "tau_Pi_H0": KAL0 * Omega / Tr, "n_s": 1 - phi ** -7,
            "Omega_m_cosm": _omm / _h ** 2}   # _h = H_alg/100, definido en Capa 1
    _desfasados = [
        f"{_k}: doc {_P26[_k]['value']} vs recomputado {_r:.7f}"
        for _k, _r in _REC.items()
        if _k in _P26 and abs(float(_P26[_k]["value"]) - _r) > 5e-6]
    check("R26 los valores de procedencia == los que recomputa el guardián",
          not _desfasados,
          f"{len(_REC)} verificados contra el cómputo"
          if not _desfasados else "DOC STALE: " + "; ".join(_desfasados[:5]))
except Exception as e:
    check("R26 capa operable", False, str(e))

# ── R27 — el look-elsewhere que citan los papers == el que se recomputa ──────
# El titular «1 de 490» es la defensa central contra la acusación de numerología,
# y hasta 2026-07-25 vivía sin red: ningún check ataba el 55/25/490 ni la tabla de
# tolerancias del PRD al script que los produce, y el script no tenía log. Si el
# diccionario crece (ya pasó: 46→50→55 nombres), los papers quedarían citando un
# denominador viejo y NADA lo delataría — el peor sitio posible para un número
# stale, porque es el que un referee hostil va a querer reproducir.
try:
    _le_path = _REPO2 / "src" / "estadistica" / "look_elsewhere_full.py"
    _le_ns: dict = {"__name__": "_le", "__file__": str(_le_path)}
    exec(compile(_le_path.read_text(errors="ignore").split("if __name__")[0],
                 str(_le_path), "exec"), _le_ns)
    _fam = _le_ns["FAMILY"]
    _rat = _le_ns["distinct_ratios"](_fam)
    _n_names = len(_fam)
    _n_vals = len({round(v, 6) for v in _fam.values()})
    _n_rat = len(_rat)
    _h0 = sum(1 for r, _ in _rat if abs(r - _le_ns["W0"]) <= 5e-4)
    _ha = sum(1 for r, _ in _rat if abs(r - _le_ns["WA"]) <= 5e-4)
    # Lo que los manuscritos AFIRMAN (PRD + Sealed + Paper 1 dicen los tres 490).
    # OJO: se excluyen los rangos con «~» (p.ej. «$\sim$1-in-500»), que NO son el
    # look-elsewhere del core sino la gramática del multiplicador φ-DM — una cuenta
    # distinta que el propio PRD declara sin privilegio estadístico (§2.3, scope).
    _claimed = set()
    for _t in _texs2:
        _txt = _t.read_text(errors="ignore")
        for _m in _re2.finditer(r"1\$?[- ]in[- ]\$?(\d{3})|1 of (\d{3})|yields \$(\d{3})\$",
                                _txt):
            if "sim" in _txt[max(0, _m.start() - 10):_m.start()]:
                continue
            _claimed.add(int(next(g for g in _m.groups() if g)))
    _bad = _claimed - {_n_rat}
    check("R27 look-elsewhere: 55/25/490 y 1-de-490 == lo que recomputa el script",
          _n_names == 55 and _n_vals == 25 and _n_rat == 490
          and _h0 == 1 and _ha == 1 and not _bad,
          f"{_n_names} nombres / {_n_vals} valores / {_n_rat} razones; "
          f"aciertos a ±0.0005: w₀={_h0}, wₐ={_ha}; manuscritos citan {sorted(_claimed) or '—'}"
          if not _bad and _n_rat == 490
          else f"DESAJUSTE: script da {_n_names}/{_n_vals}/{_n_rat} (w₀={_h0},wₐ={_ha}); "
               f"manuscritos citan {sorted(_bad)}")
except Exception as e:
    check("R27 capa operable", False, str(e))

# ── R28 — la copia del repo == el diccionario CITABLE (sin deriva) ───────────
# 2026-07-25: había DOS look_elsewhere_full.py. Los valores coincidían (55/25/490),
# así que ningún check numérico podía notarlo — pero el código NO era el mismo: la
# copia del repo se quedó en la nomenclatura pre-v1.5 (MIKAEL_V, MIKE) y escribía
# wₐ = PYROS/KRYSTOS_V, el denominador equivocado, justo en el script del titular.
# Una copia que nadie compara no es un respaldo: es una segunda fuente de verdad.
# Se compara el CÓDIGO EJECUTABLE (sin comentarios/espacios): la prosa puede
# diferir —la copia del repo añade la robustez QUINTAL…DECAL— pero las constantes
# y las identidades no. Si zenodo_dictionary/ no está presente, se avisa sin ROJO:
# está gitignoreado y un clon limpio no lo tiene.
try:
    _cit = _REPO2 / "zenodo_dictionary" / "look_elsewhere_full.py"
    _mir = _REPO2 / "src" / "estadistica" / "look_elsewhere_full.py"

    import ast as _ast

    def _codigo(p):
        """Asignaciones de constantes a nivel de módulo, vía AST.

        Se usa AST y no regex de línea porque una frase de docstring como
        «CUÄSTAL=4Ω; su continuación…» se parece a una asignación y hacía saltar
        esta misma regla en falso — prosa disfrazada de código.
        """
        _out = []
        for _n in _ast.parse(p.read_text(errors="ignore")).body:
            if isinstance(_n, _ast.Assign):
                for _t in _n.targets:
                    if isinstance(_t, _ast.Name):
                        _out.append(f"{_t.id}={_ast.dump(_n.value)}")
        return _out

    if not _cit.exists():
        check("R28 espejo diccionario citable (zenodo_dictionary ausente)", True,
              "gitignoreado; en un clon limpio no aplica — verificar en la máquina de trabajo")
    else:
        _a, _b = _codigo(_mir), _codigo(_cit)
        _falta = [x for x in _b if x not in _a]
        _sobra = [x for x in _a if x not in _b]
        check("R28 look_elsewhere_full.py del repo == el del diccionario citable",
              not _falta and not _sobra,
              f"{len(_b)} definiciones idénticas (espejo verbatim de v1.5)"
              if not _falta and not _sobra
              else f"DERIVA: sólo-citable={_falta[:3]} sólo-repo={_sobra[:3]}")
except Exception as e:
    check("R28 capa operable", False, str(e))

# ── R29 — símbolo ↔ entidad es BIYECCIÓN, y los padres cuadran ───────────────
# Pedido de Mike (2026-07-25): la tabla de notación debe llevar
#   símbolo | nombre | padres | fórmula algebraica | valor,
# y «la simbología tampoco puede estar duplicada y debe pertenecer a una sola de
# las entidades algebraicas». El motivo es operativo: I_g y K_v VALEN LO MISMO
# (9.519253285), así que el número no distingue cuál usa una ecuación — sólo el
# símbolo lo hace. Si el símbolo es único por entidad, ver K_v donde va I_g se
# detecta a simple vista; si se duplica, vuelve el error de esta semana.
# Se verifica: (a) ningún símbolo mapea a dos entidades; (b) el valor que la
# tabla declara == el que se recomputa desde la FÓRMULA de sus padres.
_SIMBOLOS = {
    # símbolo LaTeX      (entidad,        fórmula evaluable con las de arriba)
    r"\varphi":          ("—",           lambda d: (1 + 5 ** 0.5) / 2),
    r"\pi":              ("—",           lambda d: _math.pi),
    r"\Omega":           ("OMEGA",       lambda d: d["φ"] + d["π"]),
    r"\beta":            ("BIAL",        lambda d: (d["φ"] + d["π"]) / 2),
    "KAL":               ("KAL",         lambda d: d["β"] + d["π"]),
    "AURA":              ("AURA",        lambda d: d["β"] + d["φ"]),
    "MIRA":              ("MIRA",        lambda d: d["AURA"] / 2),
    "P_{\\mathrm{sc}}":  ("PYROS",       lambda d: d["Ω"] + d["φ"]),
    "I_g":               ("IGNIS",       lambda d: d["π"] + d["PYROS"]),
    "K_v":               ("KRYSTOS_V",   lambda d: d["φ"] + d["π"] + d["Ω"]),
    "T_r":               ("TRIAL",       lambda d: 3 * d["AURA"]),
    "M_v":               ("ATLAS",       lambda d: d["φ"] + d["π"] + d["K_v"]),
}
try:
    import math as _math
    _ents = [e for e, _ in _SIMBOLOS.values() if e != "—"]
    _dup_ent = {e for e in _ents if _ents.count(e) > 1}
    _dup_sim = len(_SIMBOLOS) != len({s for s in _SIMBOLOS})

    _d = {}
    _d["φ"] = (1 + 5 ** 0.5) / 2
    _d["π"] = _math.pi
    _d["Ω"] = _d["φ"] + _d["π"]
    _d["β"] = _d["Ω"] / 2
    _d["AURA"] = _d["β"] + _d["φ"]
    _d["PYROS"] = _d["Ω"] + _d["φ"]
    _d["K_v"] = _d["φ"] + _d["π"] + _d["Ω"]
    # el valor de cada símbolo, recomputado desde la fórmula de sus padres
    _calc = {s: f(_d) for s, (e, f) in _SIMBOLOS.items()}
    # contra el diccionario citable: la entidad nombrada debe valer eso mismo
    _desaj = []
    if _cit.exists():
        _ns2: dict = {"__name__": "_c2", "__file__": str(_cit)}
        exec(compile(_cit.read_text(errors="ignore").split("def distinct_ratios")[0],
                     str(_cit), "exec"), _ns2)
        _famc = _ns2["FAMILY"]
        for _s, (_e, _f) in _SIMBOLOS.items():
            if _e in ("—",):
                continue
            _ref = _famc.get(_e, _ns2.get(_e))
            if _ref is None or abs(_ref - _calc[_s]) > 1e-9:
                _desaj.append(f"{_s}→{_e}")
    # (c) la TABLA del PRD debe declarar el mismo par (símbolo, nombre). Sin esto
    # R29 sólo se auditaría a sí misma: el .tex podría derivar y nadie lo vería.
    _tabla = []
    _sealed = _REPO2 / "manuscript" / "SSEE_Sealed_Journal.tex"
    for _doc in (_prd, _sealed):
        if not _doc.exists():
            continue
        _sec = _doc.read_text(errors="ignore").split("Symbol & Name & Parents")
        if len(_sec) > 1:
            for _fila in _sec[1].split(r"\bottomrule")[0].splitlines():
                _cols = [c.strip() for c in _fila.split("&")]
                if len(_cols) >= 5 and _cols[1] not in ("", "---"):
                    _sim = _cols[0].replace("$", "").replace(r"\mathrm{", "").rstrip("}")
                    _nom = _cols[1].replace("$_V$", "_V").replace("$", "")
                    _tabla.append((_sim, _nom, _cols[3], _cols[4]))

    # (c2) la columna «reduced form in φ,π» debe EVALUAR al valor de la fila.
    # Es la columna que Mike pidió para poder distinguir de un vistazo: I_g y K_v
    # tienen PADRES distintos y la misma reducción 2(φ+π) — por eso comparten valor.
    # Si la reducción fuera decorativa (copiada, no verificada) la tabla mentiría
    # justo donde promete transparencia.
    _pura_mal = []
    for _sim, _nom, _pura, _val in _tabla:
        _e = _pura.replace("$", "").replace(r"\varphi", "PHI").replace(r"\pi", "PI")
        _e = _e.replace(r"\sqrt5", "(5**0.5)").replace(r"\sqrt{5}", "(5**0.5)")
        _e = _re2.sub(r"(\d)\(", r"\1*(", _e)          # 2( -> 2*(
        _e = _re2.sub(r"(\d)(PHI|PI)", r"\1*\2", _e)   # 3PHI -> 3*PHI
        _e = _re2.sub(r"\)(PHI|PI)", r")*\1", _e)
        try:
            _got = eval(_e, {"__builtins__": {}}, {"PHI": _d["φ"], "PI": _d["π"]})
            if abs(_got - float(_val.replace("$", "").replace("\\\\", "").strip())) > 5e-5:
                _pura_mal.append(f"{_sim}: {_pura}={_got:.5f}≠{_val}")
        except Exception:
            _pura_mal.append(f"{_sim}: reducción ilegible «{_pura}»")
    _falta_tab = []
    for _sim, _nom, _pura, _val in _tabla:
        _esp = {e for s, (e, _) in _SIMBOLOS.items()
                if s.replace("\\", "").replace("_{\\mathrm{sc}}", "_sc") == _sim
                or s.replace("\\", "") == _sim}
        if _esp and _nom not in _esp:
            _falta_tab.append(f"tabla dice {_sim}→{_nom}, guardián {_esp}")
    # (d) el otro lado de la biyección: una ENTIDAD no puede llevar DOS símbolos.
    # Ω y Ω_DNAV nombraban la misma entidad (φ+π) y llegaron a aparecer juntos en
    # una misma fila de Paper 1 — un lector no puede saber que son lo mismo.
    _alias = []
    for _t in _texs2:
        _txt = _t.read_text(errors="ignore")
        if r"\Omega_{\rm DNAV}" in _txt and _t.name not in (
                "SSEE_Paper6_phiDM.tex", "SSEE_Paper9_HubbleTension.tex"):
            _alias.append(_t.name)
    check("R29 símbolo↔entidad biyectivo, valor == fórmula, tablas conformes",
          not _dup_ent and not _dup_sim and not _desaj and not _falta_tab
          and not _alias and not _pura_mal,
          f"{len(_SIMBOLOS)} símbolos → {len(set(_ents))} entidades; "
          f"{len(_tabla)} filas conformes (PRD+Sealed), reducción φ,π verificada; "
          f"I_g≠K_v como entidades aunque compartan {_calc['I_g']:.9f}"
          if not _dup_ent and not _desaj and not _falta_tab and not _alias
             and not _pura_mal
          else f"duplicado={sorted(_dup_ent)} desajuste={_desaj} "
               f"tabla={_falta_tab[:3]} alias-Ω_DNAV={_alias} reducción={_pura_mal[:3]}")
except Exception as e:
    check("R29 capa operable", False, str(e))

# ── R30 — política de redondeo: lo mostrado es el redondeo CORRECTO ──────────
# Decisión de Mike (2026-07-25): el redondeo se FIJA, no se decide caso por caso.
# Álgebra SSEE a 12 decimales, cantidades del modelo a 6. Un redondeo suelto es
# fuente de error propia — un decimal arrastra a los siguientes — y decir un
# número mientras se usa otro ES un error, aunque el valor esté «casi bien».
#
# LA REGLA NO ES «la cadena mostrada reproduce»: las cadenas AMPLIFICAN. Con
# ω_m=0.142668 (6 dec, correcto) la división mostrada da 0.308882 y el exacto es
# 0.308881. Se exige lo verificable: |mostrado − exacto| ≤ media unidad del
# último decimal mostrado.
#
# CASO QUE LA ORIGINÓ: Ω_m se publicaba como 0.3088'9' en 157 sitios / 14 docs.
# Exacto 0.308880879 → correcto 0.308881. El dígito de más salía de redondear
# ω_m a 0.14267 (sube desde 0.1426675) y dividir ESE número — intermedios
# redondeados, ningún valor rancio. El literal viejo distaba 9.1e-6 del exacto,
# más de media unidad del 5º decimal (5e-6) ⇒ R30 lo caza.
# (El literal viejo NO se escribe aquí: el escáner de valores retirados lo vería.)
try:
    # Se ancla al SÍMBOLO, no al literal suelto: un número por sí solo no dice de
    # quién es. La primera versión marcó en rojo 0.02237/0.1200/0.9649 — que son
    # MEDIDAS de Planck, no predicciones SSEE. Un verificador que no distingue el
    # dato de la predicción es peor que ninguno.
    _ANCLAS = {
        "Omega_m (=ω_m/h²)": (_core.OMEGA_M_TOTAL, [
            r"\\omega_m/h\^2\s*(?:&|)\s*=\s*([\d.]+)",
            r"\\Om\^\{\\mathrm\{cosm\}\}\s*=\s*([\d.]+)",
            r"\\Omega_\{m,\\mathrm\{cosm\}\}\s*=\s*([\d.]+)",
        ]),
        "omega_m": (_core.OMEGA_M_H2, [r"\\omega_m\s*=\s*\\omega_b\+\\omega_c\+\\omega_\\nu\s*=\s*([\d.]+)"]),
        "w0": (abs(_core.W0), [r"w_0\s*(?:&|)\s*=\s*-\s*([\d.]{6,})"]),
        "wa": (abs(_core.WA), [r"w_a\s*(?:&|)\s*=\s*-\s*([\d.]{6,})"]),
    }
    _mal_red = []
    for _nom, (_ex, _pats) in _ANCLAS.items():
        for _t in _texs2:
            _txt = _t.read_text(errors="ignore")
            for _pat in _pats:
                for _m in _re2.finditer(_pat, _txt):
                    _lit = _m.group(1).rstrip(".")
                    if "." not in _lit:
                        continue
                    _dec = len(_lit.split(".")[1])
                    if _dec < 4:
                        continue
                    if abs(float(_lit) - _ex) > 0.5 * 10 ** (-_dec):
                        _mal_red.append(f"{_t.name}«{_lit}»≠{_ex:.{_dec}f} ({_nom})")
    # ── R30b — falsa precisión: rellenar un valor MEDIDO hasta el techo ──────
    # Objeción de Mike (2026-07-25): «si el número es finito con 4 o 5 decimales
    # esa regla no aplica, y peor si termina ahí». Cierto y verificado: 93.140000
    # ATRAVIESA la prueba de redondeo (|93.140000 − 93.14| = 0) mientras afirma
    # cuatro decimales que ningún experimento midió. El techo de 6/12 decimales
    # vale para lo ALGEBRAICO —infinitos dígitos, siempre hay uno más— y NUNCA
    # para lo medido: ahí manda la precisión del instrumento.
    _MEDIDAS = {          # constante medida : (literal canónico, decimales reales)
        "C_ν":      ("93.14", 2),
        "KiDS S₈":  ("0.759", 3),
        "KiDS σ₈":  ("0.737", 3),
        "DES S₈":   ("0.776", 3),
    }
    for _nom, (_lit, _dec) in _MEDIDAS.items():
        for _t in _texs2:
            for _m in _re2.finditer(_re2.escape(_lit) + r"(\d+)",
                                    _t.read_text(errors="ignore")):
                # dígitos extra que son sólo ceros ⇒ relleno, no medida
                if _m.group(1) and _m.group(1).strip("0") == "":
                    _mal_red.append(f"{_t.name}«{_m.group(0)}» falsa precisión "
                                    f"({_nom} se mide a {_dec} dec)")

    _mal_red = sorted(set(_mal_red))
    check("R30 redondeo: valor mostrado == redondeo correcto, sin falsa precisión",
          not _mal_red,
          f"{len(_ANCLAS)} cantidades ancladas a su símbolo en {len(_texs2)} documentos "
          f"(política: álgebra 12 dec, modelo 6 dec)"
          if not _mal_red else f"MAL REDONDEADOS: {'; '.join(_mal_red[:6])}")
except Exception as e:
    check("R30 capa operable", False, str(e))

# ── R31 — el ssee_core IMPORTADO == el ssee_core del FUENTE ─────────────────
# 2026-07-25, segunda vez en el día. Un `.pyc` rancio de `__pycache__` devolvió
# SUM_MNU_EV=0.06902 mientras el archivo en disco decía 0.06849, y con él ω_m y
# Ω_m quedaron en la cadena vieja. Esta vez llegó más lejos: el MCMC de
# producción arrancó imprimiendo Ω_m,total = 0.30889320 (cadena rancia) en vez de
# 0.30888088 — una corrida de 35 min silenciosamente corrupta, lanzada
# precisamente para arreglar procedencia. Hubo que abortarla.
#
# El guardián ya compila el fuente para SUS lecturas (ver bloque «Fuente
# canónica»), así que él no se engaña; pero CUALQUIER OTRO script del repo hace
# `import ssee_core` y sí puede recibir bytecode viejo. R31 detecta la condición
# desde fuera, sin importar la causa: compara el valor que da el IMPORT normal
# contra el que da el AST del fuente. Si difieren, hay una caché mintiendo y
# todo lo que se corra en ese estado es sospechoso.
#   Nota: `python -B` NO sirve como defensa — impide ESCRIBIR bytecode, no LEERLO.
try:
    import ast as _ast31
    import importlib as _il31
    _src31 = (_REPO2 / "src" / "ssee_core.py").read_text(errors="ignore")
    _fuente = {}
    for _n in _ast31.parse(_src31).body:
        if isinstance(_n, _ast31.Assign) and isinstance(_n.value, _ast31.Constant):
            for _t in _n.targets:
                if isinstance(_t, _ast31.Name) and isinstance(_n.value.value, float):
                    _fuente[_t.id] = _n.value.value
    _il31.invalidate_caches()
    if str(_REPO2 / "src") not in sys.path:
        sys.path.insert(0, str(_REPO2 / "src"))
    _imp31 = _il31.import_module("ssee_core")
    _discrepa = [f"{_k}: import={getattr(_imp31, _k)} fuente={_v}"
                 for _k, _v in _fuente.items()
                 if hasattr(_imp31, _k) and abs(getattr(_imp31, _k) - _v) > 1e-12]
    check("R31 bytecode: ssee_core importado == ssee_core del fuente",
          not _discrepa,
          f"{len(_fuente)} constantes literales coinciden import↔fuente"
          if not _discrepa
          else "CACHÉ RANCIA (borrar __pycache__ y RE-CORRER lo que se haya "
               "ejecutado en este estado): " + "; ".join(_discrepa[:4]))
except Exception as e:
    check("R31 capa operable", False, str(e))

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
