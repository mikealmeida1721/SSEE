# Veredicto de Claude — push de Nova (`3a9115c`, `78fa756`, `25c0b89`, `312eb91`, `2ed2a61`, `32f4fa9`, `3d57269`)

Fecha: 2026-09-25 noche. Rama: `fix/etiquetas-h0-ir`.
Protocolo §0 (decisión de Mike): verifico el push de Nova. Lo corrí yo; no leí
su reporte y lo di por bueno.

## Qué verifiqué, y cómo

- **`3a9115c` — convención de símbolo `H_0^{glob,IR}`.** ✓ **Y es mejor que lo
  que yo había hecho**: Nova renombra el SÍMBOLO, yo sólo había cambiado la
  etiqueta de texto. Con el régimen dentro del símbolo es imposible
  confundirlos, que es la propia regla de un-símbolo-por-entidad. Descarté mi
  versión y adopté la suya en los 17 sitios restantes.
- **`32f4fa9` — nota de la edad.** ✓ Reproducido **por mi cuenta**, integral
  independiente: SSEE 13.733 Gyr (él dice 13.73) y ΛCDM 13.796 (él dice 13.80).
  Control del otro lado: el ΛCDM cae a 0.01 % del publicado (13.797). Su
  corrección de categoría —edad = cantidad global ⟹ H_global, no H_local— es
  correcta.
- **`312eb91` — σ₈ 0.7446 → 0.8153.** ✓ `sigma8_ssee_unif: 0.8153` es lo que
  dice `CANONICAL_VALUES.yaml:55`. El razonamiento escrito en el comentario
  (A_s se paga una vez en el CMB y no se vuelve a cobrar en el crecimiento) es
  el correcto y está documentado en el propio script.
- **`78fa756` — 20 correcciones editoriales.** ✓ en lo muestreado: P3 l.1045
  lleva `H_0^{\mathrm{glob,IR}}`; los tres β_c de P7 (l.308, 354, 788) están en
  contexto de retiro; el de P8 es el disformal `+AURA`, deliberado y
  desambiguado.
- **`25c0b89` — relabel s_m/s_DE + convención Δ.** ✓ en el relabel.
- **`3d57269` — protocolo.** ✓ Recoge fielmente la decisión de Mike, incluido
  el borrado de la rama tras el merge.

## Observaciones — ninguna bloquea, todas ya corregidas por mí

Las cuatro son **la misma patología**, y es la que Mike lleva tiempo señalando:
un barrido por documento no ve lo que vive en otro documento ni en los cajones.

1. **β_c = −AURA seguía vivo en Paper 1.** `78fa756` lo retiró en Paper 4, pero
   el Registro Predictivo de Paper 1 lo seguía listando como Retrodicción con
   «Δ=0.2 %», 18 días después de que Paper 7 lo retirara. Corregido en
   `866c1be`. El guardián estaba verde: R58 mira si el valor aparece, no si un
   documento lo presenta como vigente mientras otro lo retira.
2. **`s_m=0.160` redondeado en Paper 2**, de la propia relabel. Lo cazó R44.
   Corregido en `866c1be`.
3. **La convención Δ que declaró contradecía a una clave canónica.** No es
   error suyo: el YAML tenía LAS DOS convenciones a la vez, y dos claves
   (`delta_AIC_ssee_lcdm`, `delta_BIC_ssee_lcdm`) llevaban el nombre de una
   dirección y el signo de la otra. Comprobado por aritmética: con Δχ²=+2.65 y
   Δk=4 sobre N=225, AIC = −5.35 y BIC = −19.01. Raíz corregida en el YAML y
   propagada a 9 sitios.
4. **El reporte A1–A19 decía «0 de los 19 cubiertos por la rama»**, y A11/A12
   ya los había cerrado él mismo en `25c0b89`. Y el radio del δc que da como
   3 sitios **es de 7**: le faltan `README.md:314` (la portada, que lo daba como
   resultado), `AUDIT.md`, `CHANGELOG.md` y `ssee_press_schechter.py`, que
   además calculaba con él. Corregidos los 7 en el commit de A1–A4.

## Para el expediente de firmas (va para Nova)

El conteo de halos masivos a z ≳ 10 estaba anotado **a favor** del modelo
(enhancement 1.05×–1.89×). Medido con el δc que el modelo realmente produce,
apunta **ligeramente en contra** en el extremo de masa alta (0.892 a
3×10¹² M☉, z=10; 0.778 a z=15) y es **neutro** a las masas que JWST mide (0.99).
Hay que asentarlo con ese signo antes de cruzarlo con la literatura.

## Veredicto: VERDE
