#!/usr/bin/env python3
"""
test_guardian.py — el vigilante vigilado ("who watches the watcher").

El guardián (`ssee_verify.py`) protege contra el FALSO VERDE (decir VERDE
cuando debería ser ROJO). Pero un check puede dejar de detectar en silencio
tras un refactor → el guardián daría falso verde sin que nadie lo note.

Este test vuelve PERMANENTES las pruebas negativas: para cada check crítico,
perturba un cajón a un estado que DEBE disparar ROJO, corre el guardián, y
exige que (a) salga con código ≠ 0 y (b) aparezca el check correcto. Luego
restaura el archivo (siempre, vía finally). Si un check deja de cazar lo suyo,
ESTE test falla → la grieta se delata en vez de pasar callada.

    .venv/bin/python3 src/verificacion/test_guardian.py
    → TODOS VERDES: el guardián detecta lo que debe.
    → ROJO: un check se debilitó; arreglarlo antes de confiar en un VERDE.
"""
import contextlib
import pathlib
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]
GUARD = "src/verificacion/ssee_verify.py"
PY = sys.executable
_passed, _failed = 0, 0


def run_guardian():
    r = subprocess.run([PY, GUARD], cwd=REPO, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr


@contextlib.contextmanager
def perturbed(relpath, transform):
    """Perturba un archivo y SIEMPRE lo restaura (bytes originales)."""
    p = REPO / relpath
    original = p.read_bytes()
    try:
        p.write_text(transform(original.decode("utf-8", errors="ignore")),
                     encoding="utf-8")
        yield
    finally:
        p.write_bytes(original)


def expect_red(name, relpath, transform, must_contain):
    global _passed, _failed
    with perturbed(relpath, transform):
        rc, out = run_guardian()
    low = out
    ok = rc != 0 and all(k in low for k in must_contain)
    print(f"  [{'PASS' if ok else 'FALLA'}] {name}")
    if not ok:
        print(f"        esperaba ROJO con {must_contain}; rc={rc}")
    _passed += ok
    _failed += (not ok)


# ── Test positivo: estado limpio → VERDE (exit 0) ───────────────────────
print("Vigilante vigilado — el guardián debe detectar cada falla simulada\n")
rc0, out0 = run_guardian()
_ok0 = rc0 == 0 and "VERDE" in out0
print(f"  [{'PASS' if _ok0 else 'FALLA'}] estado limpio → VERDE (exit 0)")
_passed += _ok0
_failed += (not _ok0)

# ── Tests negativos: cada perturbación DEBE dar ROJO ────────────────────
# F1 — valor de pipeline que NO coincide con su log de procedencia (R9)
expect_red("R9 procedencia: ΔBIC −24.7 no está en el log → ROJO",
           "VERIFICATION_LEDGER.md",
           lambda t: t.replace("−35.0 (SSEE favorecido)", "−24.7 (SSEE favorecido)"),
           ["procedencia", "ROJO"])

# R1 — claim de prioridad temporal
expect_red("R1 cronología: 'predating DESI DR2' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n% test predating DESI DR2\n",
           ["R1", "ROJO"])

# R2 — enlace al compendio multidominio retirado
expect_red("R2 multidominio: DOI 19679049 → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n% test zenodo.19679049\n",
           ["R2", "ROJO"])

# R10 — overclaim 'zero-parameter' global
expect_red("R10 overclaim: 'is a zero-parameter framework' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\nSSEE is a zero-parameter framework.\n",
           ["R10", "ROJO"])

# R11 — conteo de serie congelado
expect_red("R11 conteo serie: 'Papers~3--7' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\nextensions in Papers~3--7.\n",
           ["R11", "ROJO"])

# R21 — wₐ con denominador K_v en vez de IGNIS
expect_red("R21 wₐ/K_v: 'P_{sc}/K_v' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n$-P_{sc}/K_v$\n",
           ["R21", "ROJO"])

# R22 — precisión falsa (wₐ = -0.67000)
expect_red("R22 precisión: '-0.67000' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n$w_a=-0.67000$\n",
           ["R22", "ROJO"])

# R23 — fósil Σ₉ / MIKAEL_V / 5D
expect_red("R23 fósil: '\\Sigma_9' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n$\\Sigma_9$\n",
           ["R23", "ROJO"])

# R24 — conteo del diccionario retirado (1 of 378)
expect_red("R24 conteo: '1 of 378' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n1 of 378 ratios.\n",
           ["R24", "ROJO"])

# ── Cobertura añadida 2026-07-25 ────────────────────────────────────────
# Motivo: el drift de Σm_ν sobrevivió 15 días en VERDE porque R17 verificaba
# el STRING "94.07" y no el NÚMERO derivado. Al medir la cobertura del propio
# auto-test salió que de 16 reglas R con nombre, sólo 9 tenían mutación: las
# otras 7 se asumían buenas — el mismo falso verde, un nivel más arriba.
# Estas mutaciones cierran esa brecha. R14 (procedencia DESI) es la más
# importante del repo: existe porque los BAO venían mal etiquetados DR1/DR2.

# R15 — dos redondeos distintos de la MISMA tensión n_s en un manuscrito
expect_red("R15 n_s redondeo incoherente: 0.16σ y 0.2σ juntos → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           # DEBE caer en la ventana de R15: ancla varphi^{-7} + mención Planck.
           lambda t: t + ("\n$n_s=1-\\varphi^{-7}$ is $0.16\\sigma$ from Planck.\n"
                          "\n$n_s=1-\\varphi^{-7}$ is $0.2\\sigma$ from Planck.\n"),
           ["R15", "ROJO"])

# R17 — la constante ν vuelve a bifurcarse (94.07 reaparece en un manuscrito)
expect_red("R17 constante ν bifurcada: '94.07' reaparece → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n$\\Sigma m_\\nu = 94.07\\,\\mathrm{eV}\\,\\Omega_\\nu h^2$\n",
           ["R17", "ROJO"])

# R18 — narrativa H0-posterior stale (el r_d agrandado que V-L4-DESI retiró)
expect_red("R18 narrativa stale: 'downward to compensate' → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\nThe posterior shifts downward to compensate.\n",
           ["R18", "ROJO"])

# R19 — un valor BAO de DR1 se cuela en un manuscrito (el bug del mislabel)
expect_red("R19 BAO DR1-mislabeled: D_H=20.08 (DR2=21.863) → ROJO",
           "manuscript/SSEE_Paper7_EFT.tex",
           lambda t: t + "\n$D_H/r_d = 20.08$ for the LRG sample.\n",
           ["R19", "ROJO"])

# R14 — procedencia DESI: alterar un dígito del csv oficial DEBE doler
expect_red("R14 procedencia DESI: csv alterado ≠ Tabla 4 oficial → ROJO",
           "data/raw/desi_dr2_bao.csv",
           # OJO: "21.863" aparece antes en un COMENTARIO de cabecera; hay que
           # mutar la FILA DE DATOS o la mutación no prueba nada (lección propia).
           lambda t: t.replace("0.510,LRG1,DH_over_rd,21.863",
                               "0.510,LRG1,DH_over_rd,21.860", 1),
           ["R14", "ROJO"])

# R20 — meter la PREDICCIÓN SSEE (0.758) en el hueco del DATO KiDS (0.759)
expect_red("R20 ancla observacional: KiDS 0.759 → 0.758 (predicción) → ROJO",
           "CANONICAL_VALUES.yaml",
           lambda t: t.replace("obs_KiDS_S8:         0.759", "obs_KiDS_S8:         0.758", 1),
           ["R20", "ROJO"])

# V-L2-11 — el cierre del sector ν: desincronizar core vs YAML DEBE doler.
# Ésta es la mutación que, de haber existido el 2026-07-10, habría delatado
# el drift en el momento de escribirse la regla.
expect_red("V-L2-11 cierre ν: Σm_ν del core desincronizado del YAML → ROJO",
           "src/ssee_core.py",
           lambda t: t.replace("SUM_MNU_EV = 0.06849", "SUM_MNU_EV = 0.06902", 1),
           ["L2-11", "ROJO"])

# R25 — reintroducir la parametrización invertida (Ω_m congelado · h² para
# fabricar ω_m dentro de código SSEE) DEBE doler. Es el error que sesgó el
# posterior de H₀ hacia el ancla (0.04σ) hasta 2026-07-25.
expect_red("R25 parametrización: ω_m = Ω_m·h² con Ω_m constante → ROJO",
           "src/p02_mcmc/ssee_paper2_mcmc_reframe.py",
           lambda t: t + "\n_bad_om = OMEGA_M_TOTAL*(H0/100)**2\n",
           ["R25", "ROJO"])

# Escáner de valores retirados: un valor RETIRADO presentado como vigente
expect_red("memorias: valor retirado (ω_ν=0.000741) sin marcar → ROJO",
           "VERIFICATION_LEDGER.md",
           lambda t: t + "\n\nEl valor vigente de omega_nu es 0.000741 en toda la suite.\n",
           ["memoria", "ROJO"])


# ── Archivo: cajón sin entrada en la Bitácora → ROJO (caso especial dir) ─
def test_archive_orphan():
    global _passed, _failed
    d = REPO / "archive" / "_test_orphan_dir"
    try:
        d.mkdir(exist_ok=True)
        rc, out = run_guardian()
    finally:
        if d.exists():
            d.rmdir()
    ok = rc != 0 and "archivo" in out
    print(f"  [{'PASS' if ok else 'FALLA'}] R12 archivo: cajón sin bitácora → ROJO")
    _passed += ok
    _failed += (not ok)


test_archive_orphan()

# ── Veredicto ───────────────────────────────────────────────────────────
print(f"\n{'='*55}")
total = _passed + _failed
if _failed:
    print(f"ROJO — {_failed}/{total} pruebas del vigilante FALLARON.")
    print("Un check del guardián se debilitó: arreglarlo antes de confiar en un VERDE.")
    sys.exit(1)
print(f"VERDE — el guardián detecta las {total} fallas simuladas. Vigilante íntegro.")
sys.exit(0)
