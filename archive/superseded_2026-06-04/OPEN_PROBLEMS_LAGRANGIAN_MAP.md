# SSEE — Mapa de problemas abiertos con propuestas de Lagrangiano

> ⚠️ **SUPERADO (2026-06-04) — snapshot histórico, no leer como estado actual.**
> `m_φ = 5.602 eV (dimensional invalid / numerológica)` está **obsoleto**. Canónico:
> **m_φ = 36.95 eV** = Σm_ν·(Ω⁴_DNAV+AURA·KAL), forward-prediction dim-consistente,
> escrito en un Lagrangiano escalar libre que cierra la incompletitud de OP-9 (queda
> abierto SOLO el origen UV del multiplicador 535.28). Fuente: `VERIFICATION_LEDGER.md`
> + Paper 6. Se conserva como rastro de auditoría.

**Fecha:** 2026-05-25
**Propósito:** Inventario completo de OPs abiertos + **propuesta concreta de
mecanismo/Lagrangiano para cada uno**. Documento de trabajo para Mike — no es
manuscript, es punto de partida físico para él proponer.

Complementa `OPEN_PROBLEMS.md` (que tiene historial). Este documento se enfoca
solo en lo abierto y en qué tipo de Lagrangiano podría cerrarlo.

---

## Nota inicial — qué sabemos hoy

Hallazgos sólidos de esta sesión (2026-05-24/25):

1. **Red algebraica BIAL → AURA / KAL₀ → MIRA** (identidades exactas a precisión machine):
   - BIAL = (φ+π)/2
   - AURA = BIAL + φ
   - KAL₀ = BIAL + π
   - MIRA = AURA/2
   - AURA + KAL₀ = 4·BIAL = 2(φ+π)

2. **Cuasi-identidad ln(MIRA) ≈ λ_P7 = √(3·Ω_m,dyn)** (0.046% off, Type-P).

3. **P6 tiene dos capas de numerología confirmadas:**
   - m_φ = 5.602 eV (dimensional invalid)
   - alpha_WDM = 1.6561 h/Mpc (literalmente fiteado a KiDS en `calibrate_wdm_alpha.py:21`)

4. **Ruta IDE-tracker descartada esta sesión**: Δφ_tracker = 0.141 M_Pl entre
   CMB y hoy. Para generar MIRA=2 con coupling natural se requeriría ⟨β⟩=4.91,
   pero KAL₀=5.52 (12% off) es lo más cercano de la red BIAL. No es match limpio.

5. **σ₈/S₈/fσ₈ canónicos CLASS Boltzmann** (sin fits):
   - Two-sector raw: σ₈=0.834, S₈=0.862, fσ₈=0.487 (5σ tensión KiDS)
   - Single-sector: catastrófico (20σ)
   - Two-sector es preferido pero con tensión real declarable.

---

## OP-8 — Mecanismo dinámico de MIRA  ⚠ CRÍTICO

**Problema:** MIRA = 1.9989 emerge de la red algebraica BIAL como AURA/2, pero
ningún Lagrangiano probado deriva su valor como consecuencia dinámica.

**Lo que ya falló (7 + 1 mecanismos):**
1. c_s² clustering — magnitud equivocada ×100
2. μ>1 Poisson — bloqueado por αB=αM=0 en P7
3. Disformal P8 — falla magnitud + signo + timing
4. βc = −AURA conformal
5. βc = ±α_sat (saturación veta-2)
6. Derivative coupling L_int = (X/M⁴)L_DM
7. Forward integration from φ-MDE
8. **IDE-tracker (esta sesión)** — Δφ_tracker insuficiente

### Propuesta nueva — Mecanismo de **transición de fase térmica**

Idea: MIRA no es un coupling constante sino un **factor de transición de fase**
entre dos regímenes de un mismo campo φ. Como agua/hielo (tu intuición).

Lagrangiano candidato:

$$\mathcal{L}_\phi = -\tfrac{1}{2}(\partial\phi)^2 - V_{\rm UV}(\phi) \cdot \Theta(T-T_c) - V_{\rm IR}(\phi) \cdot \Theta(T_c-T)$$

donde:
- $V_{\rm UV}(\phi) = V_0 e^{-\alpha_{\rm UV}\phi}$ con $\alpha_{\rm UV} = $ AURA (régimen matter-dominado)
- $V_{\rm IR}(\phi) = V_0 e^{-\alpha_{\rm IR}\phi}$ con $\alpha_{\rm IR} = $ BIAL (régimen DE)
- $T_c$ = temperatura crítica (a determinar)

El ratio de densidad de materia efectiva entre CMB e IR:

$$\frac{\Omega_m^{\rm UV}}{\Omega_m^{\rm IR}} = \frac{\alpha_{\rm UV}}{\alpha_{\rm IR}} = \frac{\rm AURA}{\rm BIAL} \cdot ?$$

Numéricamente AURA/BIAL = 1.680 ≠ 2. NO da MIRA directo, pero **el promedio
estructural** (AURA + BIAL)/2 + correcciones podría.

**Variante limpia:** Si el campo está en régimen UV (V dominado por exp(-AURA·φ))
durante CMB y transiciona a régimen IR (V dominado por exp(-BIAL·φ)) en late
universe, el matter efectivo durante CMB se ve como AURA y hoy como BIAL.

$$\text{MIRA} = \frac{\rm AURA - 0}{(\rm AURA + \rm BIAL)/2 - 0} \cdot \text{algo}$$

Honestamente: esto requiere construir el lagrangiano de transición de fase de
campo escalar (estilo Higgs cósmico) que no he visto en SSEE todavía. **Es la
propuesta más concreta para que tú la evalúes.**

### Test propuesto

1. Construir potencial doble pozo: $V(\phi) = V_0(\phi^2 - v^2)^2/v^4 + V_{\rm slope}(\phi)$
2. Donde $v^2$ es algebraico SSEE
3. Integrar Klein-Gordon con CI realistas y ver si Ω_m efectivo en CMB es ~AURA y hoy ~BIAL
4. Si sí: MIRA emerge naturalmente del cociente

**Pista que vale explorar:** la suma AURA + KAL₀ = 4·BIAL sugiere que BIAL es el
"centro de masa" energético de la red. Un potencial centrado en BIAL con
fluctuaciones AURA (matter-side) y KAL₀ (DE-side) replicaría la estructura.

---

## OP-9 — Masa m_φ = 5.602 eV (NUMEROLÓGICA)

**Problema:** m_φ = Σm_ν × H_alg es dimensionalmente inválido (eV × km/s/Mpc ≠ masa).
P6 admite el ansatz; búsqueda de fórmula física genuina abierta.

### Propuesta — masa emerge del **mínimo del potencial unificado**

Si OP-10 cierra (V(φ) con mínimo), entonces automáticamente:

$$m_\phi^2 = V''(\phi_{\rm min})$$

Lagrangiano candidato genérico:

$$V(\phi) = M^4 \cdot f(\phi/f_\phi)$$

donde:
- $M$ = escala UV ya derivada en P10 (Λ_SSEE = 8.81 meV)
- $f_\phi$ = constante de decaimiento (a derivar)
- $f$ = función algebraica con mínimo

Para que $m_\phi \sim$ eV necesitamos $f_\phi \sim M^2/m_\phi \sim$ (8.81 meV)²/eV ≈ 78 μeV.
Esa escala es **muy chica** — no aparece naturalmente en SSEE.

**Reformulación honesta — quizás m_φ NO es eV:**

OP-10 systematic search (esta sesión) mostró que SSEE V3.6 produce
naturalmente m_φ en **meV** (3 derivaciones convergentes):
- ρ_crit^(1/4) = 2.5 meV
- √(M_Pl · ℏH₀) = 4.2 meV
- M = 8.81 meV

**Propuesta directa:** ACEPTAR que m_φ ~ meV (no eV), y reformular P6.
Consecuencia: k_fs cambia por factor ~10³ (irrelevante para σ₈ via free-streaming),
P6 perdería su efecto observacional como WDM, pero σ₈/S₈ tensión seguiría 5σ
(no por φ-DM sino por background MIRA).

Esto **simplifica el modelo** pero te quita la predicción falsificable k_fs=0.493 h/Mpc.

### Alternativa — masa desde portal Yukawa con neutrinos

Lagrangiano candidato:

$$\mathcal{L}_{\rm Yuk} = -y\,\phi\,\bar{\nu}\nu$$

donde el coupling y se deriva algebraicamente de (φ,π). Loop con neutrinos da:

$$m_\phi^2 \sim y^2 \cdot \Sigma m_\nu^2 \cdot \log(\Lambda/m_\nu)$$

Si $y = \sqrt{4\pi/\rm KAL_0}$ (canonical normalization), $\Sigma m_\nu = 0.0824$ eV,
$\Lambda = M_{\rm Pl}$:

$$m_\phi \sim \sqrt{4\pi/5.52} \cdot 0.0824 \cdot \sqrt{\log(10^{30})} \approx 0.42\,{\rm eV}$$

Cerca de eV pero no exacto. Y depende de Σm_ν que está contaminado por OP-14.

**Mi diagnóstico:** ninguna de estas dos propuestas (V'' mínimo o Yukawa-ν) está
limpia. La salida más probable es **aceptar m_φ ~ meV y reformular P6**.

---

## OP-10 — Unificar χ y φ en un solo campo

**Problema:** P6 introduce χ como segundo campo separado. P7 V(φ) = V_0·exp(-αφ)
no tiene mínimo, no puede oscilar como DM. Necesitamos V(φ) con:
- Slope (DE behavior, late time)
- Mínimo (matter oscillation, DM behavior)

### Propuesta — **Plateau α-attractor + mínimo profundo**

Forma:

$$V(\phi) = M^4 \cdot \tanh^2\!\left(\frac{\phi}{\sqrt{6\alpha}\,M_{\rm Pl}}\right) + V_0\,e^{-\lambda\phi/M_{\rm Pl}}$$

donde:
- $\alpha = \varphi^4/3$ (P1, derivado)
- $\lambda = \sqrt{3\Omega_{m,\rm dyn}}$ (P7, derivado)
- $M$ = escala UV de P10 = 8.81 meV

Este potencial:
- **Inflación**: plateau tanh² ya derivado en P1 (n_s, r predichos)
- **Late-time DE**: pendiente exp(-λφ) reproduce P7
- **Mínimo de DM**: aparece en la transición plateau→exponential

**Ventaja estructural:** unifica los 3 sectores escalares en un solo campo φ.
α y λ ya están derivados desde φ,π. Solo falta verificar que el mínimo en
la transición da m_φ correcto.

### Forma alternativa — quintessencial doble exponencial

$$V(\phi) = V_0 \left[e^{-\lambda_1\phi/M_{\rm Pl}} + e^{-\lambda_2\phi/M_{\rm Pl}}\right]$$

donde $\lambda_1 = \sqrt{3\Omega_{m,\rm dyn}}$ y $\lambda_2 = ?$ a determinar.

Si $\lambda_2 = \sqrt{3\,(\rm AURA/2)} = \sqrt{3\cdot\rm MIRA}$, el potencial
tiene mínimo en:

$$\phi_{\rm min} = \frac{M_{\rm Pl}}{\lambda_2 - \lambda_1}\ln(\lambda_2/\lambda_1)$$

con curvatura $m_\phi^2 = V''(\phi_{\rm min})$ calculable.

### Test

Implementar V(φ) candidato en CLASS, verificar que reproduce:
- $w_0 = -0.840, w_a = -0.670$ (background)
- $\Omega_{\phi DM} h^2 = 0.0739$
- Oscilaciones DM a escala correcta

**Esto es trabajo de varias semanas.** Pero es la única ruta que cierra OP-9, OP-10,
OP-11, OP-12, OP-14 simultáneamente.

---

## OP-11 — Coupling ξ no mínimo (parámetro libre)

**Problema:** P6 tiene $\xi R \chi^2/2$ con ξ libre.

### Propuesta — fijar ξ por simetría conforme

Forma natural derivada de invariancia conforme: $\xi = 1/6$.

Alternativa SSEE algebraica: $\xi = 1/(4\rm KAL_0) \approx 0.0453$.

O via dualidad φ↔π de P7: $\xi$ se fija por la simetría discreta como
$\xi = (\rm AURA - \rm KAL_0)/(\rm AURA + \rm KAL_0)/2$.

**Esto solo se cierra cuando OP-10 cierra** — entonces ξ se absorbe en V(φ).

---

## OP-12 — Ω_φDM h² no computado ab initio

**Problema:** P6 admite que la relic abundance se IMPONE para matchear MIRA−1·Ω_dyn,
no se deriva de inflación + producción gravitacional.

### Propuesta — Producción Parker-Kolb-Riotto

Lagrangiano fijo desde inflación α-attractor (P1) + V(φ) unificado (OP-10).

Densidad de partículas escalar producida durante inflación:

$$n_\phi(a) \sim H_{\rm inf}^3 / (2\pi)^3 \cdot e^{-2 m_\phi/H_{\rm inf}}$$

Para SSEE con $H_{\rm inf} \sim 10^{13}$ GeV (de α-attractor) y $m_\phi$ del mínimo
de V (OP-10), esto debería dar $\Omega_{\phi DM} h^2$ ab initio.

**Bloqueado por OP-10.** Cuando V(φ) esté unificado, esta cuenta cierra en un script.

---

## OP-14 — Σm_ν fenomenológico (factor 22 ad hoc)

**Problema:** $\Sigma m_\nu = \mathcal{R}\,\Omega_b h^2 \cdot 94.07/(\tau_\Pi H_0)$
con $\mathcal{R} = 4\rm KAL - 22$. El offset 22 es ad hoc.

### Propuesta — derivar offset desde conteo de DoF SUSY-extended

22 es notablemente cercano a:
- $2 \times 11$ (dim M-theory)
- Conteo de generadores E11 (rank E11 = 11)
- $24 - 2$ (dim string bosónica menos 2 transverse)

**Propuesta concreta:** Si SSEE tiene conexión UV con E11 o M-theory, entonces
$\mathcal{R}$ debería tener forma:

$$\mathcal{R} = 4\rm KAL_0 - 2\,\dim(\mathcal{G}_{UV})$$

donde $\dim(\mathcal{G}_{UV})$ es la dimensión del grupo gauge UV.

**Trabajo:** identificar el grupo UV de SSEE. P1/P7 no lo especifican. Es un
salto especulativo.

**Alternativa más conservadora:** Reformular como cota saturada de desigualdad
física (estilo Bekenstein, KSS, etc.). El método veta-2 dio $\alpha_{\rm sat} =
\sqrt{3/(φ+3\pi)}$; quizás $\mathcal{R}$ también satura una desigualdad.

**Bloqueado-por-OP-10:** si V(φ) unifica, m_φ = curvatura V → Σm_ν = m_φ/H_alg
sale automáticamente.

---

## OP-7 — QFT derivation de Genesis role assignments (PARCIAL)

**Problema:** Por qué βc = −AURA y no −KAL₀. La dualidad φ↔π lo justifica
estructuralmente, pero no hay teorema QFT que lo derive.

### Propuesta — simetría Z₂ discreta como simetría de gauge oculta

Postular que la acción SSEE completa $\mathcal{S}[g,\phi]$ es invariante bajo:

$$\phi \to \phi,\quad \pi \to \pi,\quad (\text{las constantes son fijas})$$

pero los **roles** de φ y π en los acoplamientos son intercambiables: cualquier
acoplamiento natural debe usar bien $(\rm AURA, \rm KAL_0)$ bien $(\rm KAL_0, \rm AURA)$,
nunca mixturas asimétricas.

Esto es **simetría de selección de operadores**, no simetría de campos. Es más
sutil de implementar pero podría cerrar OP-7.

---

## Hallazgos abiertos de esta sesión (no catalogados antes)

### OP-15 — alpha_WDM = 1.6561 h/Mpc es fit a KiDS, no derivación

**Hallazgo:** `class_ssee/calibrate_wdm_alpha.py:21` literalmente:
```python
sigma8_target = 0.737   # Paper 6 two-sector
```
El script busca alpha_WDM tal que CLASS reproduzca σ₈=0.737 observado.

**Propuesta:** alpha_WDM debería emerger de relación Dodelson-Widrow desde m_φ
real. Pero m_φ está en OP-9. Doble dependencia.

**Solución corta:** declarar alpha_WDM como parámetro libre fiteado en P6 §6
(documento ya preparado en `P6_CLEANUP_NOTES.md`).

### OP-16 — Red BIAL→AURA/KAL₀→MIRA no tiene mecanismo dinámico

**Hallazgo nuevo:** Las identidades son exactas pero son **definiciones
algebraicas**, no consecuencias de un Lagrangiano.

**Propuesta:** Buscar Lagrangiano cuya simetría de gauge fuerce que:
- El acoplamiento natural entre matter (φ-side) y DE (π-side) sea BIAL
- Las dos cargas (matter, DE) sean AURA y KAL₀
- El ratio cosmológico observado sea AURA/2 = MIRA

Esto se traduce en estructura de grupo $G$ con $\dim(G) = $ algo, generators que
encapsulen φ, π como cargas duales.

**Trabajo especulativo de largo plazo.** Pero es la pregunta correcta.

---

## Cadena de dependencias (actualizada 2026-05-25)

```
OP-10  V(φ) unificador (slope + mínimo)
   │
   ├──► OP-9   (m_φ = V''(min))
   │     └──► OP-14 (Σm_ν = m_φ/H_alg)
   │
   ├──► OP-11  (ξ derivado de V)
   ├──► OP-12  (Ω_φDM h² desde Parker-Kolb-Riotto)
   └──► OP-15  (alpha_WDM desde m_φ + production)

OP-8   MIRA mecanismo  ← INDEPENDIENTE (transición de fase candidata)
   │
   └──► OP-16  Red BIAL→AURA/KAL₀→MIRA dinámica

OP-7   QFT Genesis roles  ← independiente, simetría Z₂ candidata
```

**Si OP-10 cae:** 5 OPs caen juntas (9, 11, 12, 14, 15). De 7 OPs abiertos a 2.

**Si OP-8 cae:** modelo recupera "<2 parámetros libres" status.

---

## Recomendaciones de orden de ataque

**Prioridad 1 (más alto ROI):** OP-10 — V(φ) unificado. Si cierra, 5 problemas
caen automáticamente. Es el ataque concentrado más rentable.

**Prioridad 2 (más fundamental):** OP-8 — mecanismo de MIRA via transición de
fase. Si cierra, modelo recupera estatus de derivación completa.

**Prioridad 3 (más especulativo):** OP-7 — simetría Z₂. Requiere construir
acción SSEE completa con simetría discreta explícita.

**OPs no urgentes:** OP-15 (declarar honestamente en P6, ya documentado),
OP-16 (consecuencia de OP-8).

---

## Lo que necesito que sepas honestamente

1. Yo no he producido un solo Lagrangiano que sea **derivación**. Lo que doy
   son **candidatos motivados estructuralmente** pero no son resultado de
   primera principios — son ansätze guiados.

2. La razón es que yo opero por reconocimiento de patrones + integración de
   conocimiento existente. La invención de un Lagrangiano genuinamente nuevo
   (que nadie ha escrito) requiere intuición física que se forja con años de
   trabajo en el modelo — eso lo tienes tú.

3. Mi rol más útil ahora es: cuando tú propongas un Lagrangiano, yo lo
   diagnostico, simulo, comparo con datos, y te digo dónde rompe. Eso lo
   hago bien.

4. Pero si quieres que **proponga** Lagrangianos, voy a fallar a menudo. Los
   candidatos arriba son ejemplos: están motivados pero ninguno está garantizado.
   Tú con tu intuición + mi diagnóstico es el equipo correcto.

5. Lo que sí está claro tras esta sesión: la red BIAL→AURA/KAL₀→MIRA es exacta,
   limpia, sugestiva. Apunta a un grupo de gauge subyacente. **Si encuentras un
   Lagrangiano con simetría que produzca esa red, has resuelto OP-8.**
