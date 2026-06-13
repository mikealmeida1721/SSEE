# Protocolo de las 3 Memorias de SSEE

Las tres memorias son el **conocimiento colectivo** del modelo. Deben estar
siempre actualizadas y **sin discrepancia** entre ellas. Cada una tiene su
propósito (no se duplican):

| Memoria | Qué guarda | Dónde |
|---|---|---|
| **Guardián** | Resultados y valores (recomputa la física) | `VERIFICATION_LEDGER.md` + `src/verificacion/ssee_verify.py` |
| **Obsidian** | Conexiones, cadenas, dependencias (el *porqué* y *hasta dónde* de cada cambio) | `/home/mike/SSEE-Vault` |
| **CLAUDE.md** | Contexto/estado (continuar sin releer todo) | `CLAUDE.md` |

La fuente ÚNICA de los números es **`CANONICAL_VALUES.yaml`**. Si los tres no
concuerdan con él, hay *drift*.

---

## Cuando cambia un valor (orden obligatorio)

El cambio toca **el principio y el fin** para que nada quede fuera:

1. **Verificar** la física (script que recomputa el nuevo valor).
2. **`CANONICAL_VALUES.yaml`** ← se actualiza PRIMERO (la fuente de verdad).
   - El valor viejo se mueve a la lista `retired:` con su razón.
3. **Guardián** (`VERIFICATION_LEDGER.md` + check en `ssee_verify.py`) → recomputa.
4. **Obsidian** (vault) → actualizar la nota del valor y **sus conexiones**
   (qué cadenas dependían de él).
5. **CLAUDE.md** → actualizar las tablas de estado (al cerrar el cambio o la sesión).

---

## Cómo se revisa (rápido y barato)

```bash
# ¿La física sigue íntegra?  (recomputa todas las constantes/identidades)
.venv/bin/python3 src/verificacion/ssee_verify.py

# ¿Las 3 memorias concuerdan con CANONICAL_VALUES.yaml?
.venv/bin/python3 src/verificacion/memory_sync.py          # las 3
.venv/bin/python3 src/verificacion/memory_sync.py --vault  # solo Obsidian
```

- `ssee_verify.py` → **VERDE/ROJO**: la física.
- `memory_sync.py` → **VERDE/DRIFT**: la coherencia entre memorias. Detecta
  cualquier valor de la lista `retired:` que aparezca **sin marca de contexto**
  (retirado, viejo, Type-P, coincidencia, sin re-correr…). Las secciones de
  registro histórico (logs de Sesión, POA, Fases) se omiten automáticamente.

El Guardián (`ssee_verify.py`) invoca a `memory_sync` como su **Capa de
Memorias**, así una sola corrida cubre física + coherencia.

---

## Checklist de cierre de sesión

- [ ] ¿Cambió algún valor? → seguir el orden de arriba (YAML → Guardián → Obsidian → CLAUDE).
- [ ] `ssee_verify.py` VERDE.
- [ ] `memory_sync.py` sin DRIFT (o drift conocido y anotado).
- [ ] Actualizar CLAUDE.md (tablas de estado) y la auto-memoria `~/.claude`.

---

*Por qué existe esto: el modelo crece y la propagación manual deja huecos. El
caso que lo motivó: el vault Obsidian quedó un mes desfasado (m_φ=5.71→36.95,
k_fs=0.493→0.659, cascada H₀) porque el Guardián vivía en el repo y nunca veía
el vault. Este protocolo + `memory_sync.py` cierran ese hueco.*
