# SN Encore: tercera H₀ por supernova con lente — nota, NO confirmación

**Fuente:** Pierel et al. / Suyu et al. 2025, arXiv:2509.12301 (retrasos y H₀) y
arXiv:2509.12319 (comparación de modelos de lente). Traída por Max el 2026-09-19.

## El dato

| | valor |
|---|---|
| Objeto | SN Ia *Encore*, z = 1.949, lente por el cúmulo MACS J0138.0−2155 |
| Retraso | Δt₁b,₁ₐ = −39.8 (+3.9/−3.3) días |
| **H₀** | **66.9 (+11.2 / −8.1) km/s/Mpc**, análisis a ciegas, 7 modelos de lente |
| Lo que viene | *Requiem* (misma galaxia): su imagen tardía daría H₀ al **~2–3 %** |

## Contra qué se lee (distancia en σ, usando el brazo del lado que toca)

| referencia | valor | distancia |
|---|---|---|
| H_global UV (P10) | 67.962142 | **0.09σ** 🟢 |
| H_global IR (P9) | 68.13 | **0.11σ** 🟢 |
| Planck | 67.36 | 0.04σ 🟢 |
| SH0ES | 73.04 | 0.55σ 🟢 |

Con ±10 **no discrimina nada**: todo cae dentro de 1σ.

## Corrección a la lectura de Max

La nota decía: «si Requiem confirma la convergencia [hacia ~68], la cascada de
f_screen pierde su razón de ser». **Es al revés.**

La cascada de P9/P10 dice: la escalera LOCAL (cefeidas → SH0ES) está
apantallada y lee de más; lo GLOBAL se expande a H_global = 73.04·(1−f_screen)
≈ 67.96–68.13. Un retraso de lente a z ≈ 2 mide la distancia a lo largo de
toda la expansión, sin cefeidas: es una sonda **global**. Entonces:

- **Requiem da ~68 al 2–3 %** ⟹ la sonda global coincide con H_global y no con
  SH0ES: es lo que la cascada **predice**. Es un acierto, no un problema.
- **Requiem da ~73 al 2–3 %** ⟹ la expansión global es la de SH0ES y no hay
  nada que apantallar: **eso sí falsa f_screen.**

Con σ ≈ 1.4–2.0 km/s/Mpc, la separación 73.04 − 67.96 = 5.08 queda a
**2.5–3.7σ**: Requiem es una prueba falsable real de la cascada. Hay que
registrarla como predicción **antes** de que salga el dato (R24).

## Lo que va en contra y hay que tener a la vista

Los cuásares con lente (H0LiCOW, Wong et al. 2020) dieron 73.3 (+1.7/−1.8),
del lado de SHOES. Si la cascada dice que las sondas globales leen ~68, ese
resultado es el que la tensiona, y cualquier texto que cite a Encore a favor
tiene que citar también a H0LiCOW. (Los análisis posteriores de TDCOSMO cambian
con el perfil de masa asumido; **no verificado aquí**, hay que leerlos antes de
usarlos.)

## Estado

Nota de interés. No se toca ningún paper ni número. Pendiente, si Mike lo
decide: registrar en P9 la predicción «H₀ de Requiem = 67.96 ± su σ; 73 lo
falsa» con fecha anterior al dato.

## Actualización (mismo día, tras el intercambio con Max)

Max corrigió su lectura («invertí la conclusión») y aportó TDCOSMO IV
(Birrer+2020): 74.5 (+5.6/−6.1) solo, 67.4 (+4.1/−3.2) con SLACS — verificados.
Hay uno más reciente, **TDCOSMO-2025** (arXiv:2506.03023): **71.6 (+3.9/−3.3)**
en ΛCDM plano, con SLACS y SL2S en acuerdo.

**Lo que faltaba en la conversación:** todas esas H₀ se infieren con la FORMA de
expansión de ΛCDM. Con el fondo de SSEE (w₀, wₐ del núcleo) la misma medición
da una H₀ un **1.6–1.9 % más baja**. Medido en
`src/p09_hubble/h0_lente_fondo_ssee.py` → `results/logs/h0_lente_fondo_ssee.log`
(aproximación declarada: Ω_m igual en los dos, efecto de forma a primer orden).
Distancias a H_global, leída en ΛCDM → leída con el fondo SSEE: TDCOSMO-2025
1.10σ → 0.70σ; H0LiCOW 2.97σ → 2.24σ (sigue siendo el más incómodo).

**Consecuencia para la predicción de Requiem:** hay que fecharla con las dos
lecturas, porque el corrimiento (~1.6 %) es del tamaño de su error (2–3 %):
- analizada con el fondo de SSEE: H₀ = H_global;
- analizada en ΛCDM plano (como se publicará): H_global / 0.9838 — el valor
  exacto está en el log (`H_global_leido_en_LCDM` del sistema MACS J0138).
Falsación: un análisis ΛCDM de Requiem cerca de 73 (al 2–3 %).
