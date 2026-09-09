# La celda que faltaba: con la MISMA rigidez, ΛCDM está MÁS tenso que SSEE

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

## 3. Tres lecturas, y la primera desinfla algo que me habría gustado decir

**(a) La rigidez es la causa, no SSEE.** ΛCDM con el fondo clavado da **D = 1.13**,
prácticamente igual que el 1.16 de SSEE, y con el mismo socio (`halo_A`). O sea
que *«KiDS mide A_s»* **no es una virtud de SSEE: es lo que le pasa a cualquier
modelo cuyo fondo esté fijo.** Lo específico de SSEE es que el suyo está fijo
**por álgebra** y no por decreto — pero la mejora en la medición es genérica.
Esto había que decirlo antes que lo siguiente.

**(b) Con el mismo protocolo, la brecha de SSEE es MENOR.** Mismos datos, mismos
9 libres, mismo canal, mismos priors: **3.63σ contra 4.65σ.** Y el ajuste
empata (Δχ²_min = 0.151 a favor de SSEE, que es cero sobre 225 puntos). SSEE no
ajusta mejor; **necesita estirar menos la amplitud.**

Por el criterio pre-registrado, esto responde la pregunta que lo motivó:
**la tensión está en el DATO. El fondo algebraico de SSEE no la fabrica** — y de
hecho la deja más pequeña que el fondo de Planck.

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
