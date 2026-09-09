# Con el fondo clavado, KiDS pide el MISMO A_s en los dos modelos: 0.45σ

> 🔧 **REESCRITO 2026-09-09 tras la corrección de Mike.** La primera versión
> titulaba *«ΛCDM está más tenso que SSEE»* y sacaba al frente 3.63σ contra
> 4.65σ. **Ese encuadre está mal, por tres razones que él señaló:**
>
> 1. **Que ΛCDM rígido salga tenso ya se sabe** — es la tensión S₈ publicada.
>    No hacían falta 15 horas para eso.
> 2. **Compara como si el fondo de SSEE pudiera ser libre.** No puede: está
>    fijado por álgebra. La fila de ΛCDM-libre no es una alternativa que SSEE
>    tenga; es otro modelo.
> 3. **No son la misma tensión.** La de SSEE está **localizada en A_s**, porque
>    con todo lo demás fijo A_s es el único sitio donde puede aparecer. La del
>    modelo estándar se reporta en **S₈**, un estadístico comprimido que **no
>    puede decir dónde** está el problema. Llamarlas «la misma» borra justo lo
>    que hace informativa a la nuestra.
>
> **Para lo que se lanzó esta corrida** era para ver **dónde cae el A_s de ΛCDM**
> cuando se le clava su fondo y se le suelta solo la amplitud — y si cae cerca
> del que saca SSEE con el suyo. Eso es lo que mide, y es lo que va arriba.

**Corrida:** `src/p06_growth/cobaya_kids.py lcdmfijo` · 4 cadenas MPI ·
**14.91 h** · convergió (R−1 medias = **0.025014**, criterio 0.03; bordes
0.1176) · 24 805 pasos aceptados ·
**Log:** `results/logs/kids_lcdm_fondofijo_reparto.log` ·
**No toca ningún paper.**

## 1. El criterio estaba escrito ANTES de correr

En el propio código, `cobaya_kids.py` L160-165:

> *«Si ΛCDM-fijo también sale en tensión, la tensión está en el DATO y ninguno
> de los dos modelos la fabrica. Si solo SSEE sale en tensión, entonces su fondo
> algebraico es el sospechoso.»*

Comiteado el **2026-09-08 03:41** (`01d4369`). La cadena arrancó ~18 h después.
**Pre-registro verificado en git, no afirmado.**

## 2. La tabla, con su control

**Control primero (R24):** el medidor de degeneración recupera **D = 7.083** en
una degeneración plantada de r=0.99, donde el álgebra da 7.089. **PASA.**

ΛCDM-fijo va clavado en Planck 2018 (ω_b=0.02237, ω_c=0.1200, h=0.6736,
n_s=0.9649) — verificado en el código, L167.

| | logA que pide KiDS | **D** | χ²_min | libres | **vs el CMB** |
|---|---|---|---|---|---|
| **SSEE**, fondo fijo por álgebra | 2.8633 ± 0.0500 | **1.16** | 265.440 | 9 | **3.63σ** |
| **ΛCDM**, fondo fijo en Planck | 2.8327 ± 0.0456 | **1.13** | 265.591 | 9 | **4.65σ** |
| ΛCDM, fondo **libre** | 3.4774 ± 0.4905 | 4.73 | 262.746 | 13 | −0.88σ |

(El fondo cósmico pide 3.0448 en SSEE y 3.0451 en ΛCDM — el mismo número.)

## 2bis. EL RESULTADO — que en la primera versión quedó enterrado

**Con el fondo clavado y solo la amplitud suelta, KiDS pide lo mismo en los dos
modelos:**

| | logA que pide KiDS |
|---|---|
| SSEE, fondo fijo por álgebra | 2.8633 ± 0.0500 |
| ΛCDM, fondo fijo en Planck | 2.8327 ± 0.0456 |
| **diferencia** | **0.0306 = 0.45σ — coinciden** |

Y el déficit contra el A_s que pide **su propio** fondo cósmico es el mismo
agujero:

| | pide el CMB | pide KiDS | déficit | en σ₈ |
|---|---|---|---|---|
| SSEE | 3.0448 | 2.8633 | 0.1815 | **−9.1%** |
| ΛCDM | 3.0451 | 2.8327 | 0.2124 | **−10.6%** |

Los dos déficits difieren en **0.46σ**: es **el mismo agujero**.

### El veredicto que esto da

**El fondo está bien — el de los dos.** Si le clavas a KiDS el fondo que su
propio CMB prefiere y le sueltas solo la amplitud, **siempre pide esta cantidad
de A_s**, y da igual de qué modelo venga el fondo. El número que sale es una
propiedad de **KiDS**, no del modelo que se le ponga delante.

Por tanto el agujero **no lo fabrica ningún fondo**. Lo que hay es que **el CMB
y KiDS no perciben la misma cantidad de materia**. Y eso ya no es cuestión de
amplitud: es cuestión de **cómo mira cada uno**. Uno pesa el total en un
instante temprano; el otro ve cómo se cae la materia a lo largo del camino. Si
hay materia que no se agrupa, el primero la cuenta y el segundo no la ve.

**Eso es una hipótesis, no un resultado, y tiene su prueba: la cola #19** —
distinguir si la supresión es **plana** (amplitud, todas las escalas por igual)
o **con escala** (una fuga por debajo de cierto corte). KiDS mide ξ± en un rango
de escalas, así que el dato puede decidirlo. **Verificado: la #19 NO está
lanzada, ni escrita.** La conclusión estaba, la prueba no.

## 3. Tres lecturas, y la primera desinfla algo que me habría gustado decir

**(a) La rigidez es la causa, no SSEE.** ΛCDM con el fondo clavado da **D = 1.13**,
prácticamente igual que el 1.16 de SSEE, y con el mismo socio (`halo_A`). O sea
que *«KiDS mide A_s»* **no es una virtud de SSEE: es lo que le pasa a cualquier
modelo cuyo fondo esté fijo.** Lo específico de SSEE es que el suyo está fijo
**por álgebra** y no por decreto — pero la mejora en la medición es genérica.
Esto había que decirlo antes que lo siguiente.

**(b) El criterio pre-registrado se cumple: la tensión está en el DATO.** El
fondo algebraico de SSEE **no la fabrica**, porque el fondo de Planck produce el
mismo agujero (0.46σ de diferencia entre los dos déficits).

⚠️ **Lo que NO hay que leer aquí es «SSEE gana 3.63 contra 4.65».** Los dos
números salen del mismo agujero visto con barras algo distintas, y esa
diferencia de 0.45σ es ruido, no mérito. Además compararlos así insinúa que
SSEE podría haber salido con el fondo libre, y no puede: el suyo está fijado por
álgebra. **La lectura correcta es la coincidencia, no la carrera.**

**(c) Soltar el fondo no resuelve nada: apaga el instrumento.** Con ΛCDM libre no
hay tensión (−0.88σ)… porque la barra es **diez veces más ancha** (0.49 contra
0.046) y `logA` se enreda con `ω_c` a r=−0.95, con D=4.73. Ese acuerdo no lo
compra el dato: lo compra la indeterminación. Es el mismo espejismo que ya
cazamos en BOSS, y por eso la columna D va siempre al lado del número.

## 4. Lo que NO afirmo

- **Que SSEE resuelva la tensión.** 3.63σ sigue siendo tensión. La reduce
  respecto a ΛCDM-rígido; no la cierra.
- **Que sea la comparación definitiva ΛCDM vs SSEE.** ΛCDM-rígido no es el ΛCDM
  que nadie defiende: su versión honesta tiene el fondo libre, y ésa gasta 13
  libres en vez de 9. La fila (b) iguala la **rigidez**, no la libertad — que es
  justo lo que hacía falta para preguntar *de quién es la tensión*, y nada más.
- **Que el 3.63σ sea el número final.** Sigue siendo una **marginal** (cadena)
  contra un **perfil** (minimización del CMB). Eso es la cola #22, y hasta
  entonces la cifra no es citable con esa precisión.

## 5. Qué queda pendiente de esta corrida

Escribir la fila en Paper 6 **no está autorizado** hasta que Mike lea esto. Y
antes hay que emparejar las reglas (#22), porque el 3.63σ y el 4.65σ arrastran
los dos el mismo defecto de emparejamiento — aunque su **diferencia** sí es
robusta, porque ambos se midieron igual.
