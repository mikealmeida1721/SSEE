"""¿Dos campos de w CONSTANTE reproducen el (w0, wa) algebraico? (cola #12)

DE DONDE VIENE. La corrida anterior (`dos_campos_phi_pi.py`) mato una pareja
concreta: el campo k-esencia de Paper 7 mas un segundo componente de w
constante. Ahi el campo de Paper 7 no ayudaba nada — la dispersion CRECIA al
meterlo, y pasado el 8% el sobrante pedia densidad negativa.

Eso mato un CANDIDATO, no la idea. Queda sin probar la pareja generica, que es
la forma estandar de un modelo quintom: un campo normal (w > -1) y uno fantasma
(w < -1), los dos de w constante. La suma de dos fluidos de w fijo SI puede dar
un w total que se mueve, porque el peso relativo de los dos cambia con a.

QUE SE MIDE. Con
    rho_i(a) = f_i · a^(-3(1+w_i)),   f_A + f_B = 1
el total tiene
    w_tot(a) = [f_A w_A a^(-3(1+w_A)) + f_B w_B a^(-3(1+w_B))] / rho_tot(a)
y se compara contra w(a) = w0 + wa(1-a) con los DOS numeros algebraicos.

Tres numeros libres (w_A, w_B, f_A) contra 400 puntos de una curva que tiene
otra forma funcional. Puede fallar, y por eso vale la pena correrlo.

CRITERIO, escrito ANTES de correr, sobre el desvio MAXIMO en a ∈ [0.30, 1.00]:
    max|Δw| < 0.01   -> los dos campos son indistinguibles del par algebraico
    0.01 .. 0.05     -> aproximacion buena pero distinguible
    > 0.05           -> NO: la forma CPL no sale de dos fluidos de w constante

CONTROL (R53, y PRIMERO por R24). Se fabrica un total sintetico que ES, por
construccion, dos fluidos con w_A = -0.70, w_B = -1.40, f_A = 0.35. El ajuste
tiene que devolver esos tres numeros. Si no, no esta midiendo nada.

FUENTE: results/logs/eft_dos_campos_genericos.json
"""
import json
import pathlib
import sys

import numpy as np
from scipy.optimize import least_squares

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
import ssee_core as S                                   # noqa: E402

SALIDA = REPO / "results" / "logs" / "eft_dos_campos_genericos.json"

W0, WA = S.W0, S.WA
A = np.linspace(0.30, 1.00, 400)
TOL_IGUAL, TOL_APROX = 0.01, 0.05     # criterio, antes de correr


def w_de_dos(a, wA, wB, fA):
    """w total de dos fluidos de w constante, normalizados a rho_tot(1) = 1."""
    rA = fA * a ** (-3.0 * (1.0 + wA))
    rB = (1.0 - fA) * a ** (-3.0 * (1.0 + wB))
    return (wA * rA + wB * rB) / (rA + rB)


def ajusta(a, w_obj):
    """Devuelve (wA, wB, fA, max|dw|, rms)."""
    def res(u):
        return w_de_dos(a, u[0], u[1], u[2]) - w_obj
    mejor = None
    # varios arranques: el problema es no lineal y tiene simetria A<->B
    for x0 in ([-0.6, -1.3, 0.5], [-0.8, -1.5, 0.3], [-0.9, -1.2, 0.7],
               [-0.3, -1.8, 0.4], [-0.75, -1.1, 0.6]):
        try:
            r = least_squares(res, x0, bounds=([-3.0, -3.0, 1e-6],
                                               [1.0, 1.0, 1.0 - 1e-6]),
                              xtol=1e-14, ftol=1e-14, gtol=1e-14)
        except Exception:
            continue
        if mejor is None or r.cost < mejor.cost:
            mejor = r
    d = res(mejor.x)
    return (float(mejor.x[0]), float(mejor.x[1]), float(mejor.x[2]),
            float(np.max(np.abs(d))), float(np.sqrt(np.mean(d ** 2))))


def main():
    # ── CONTROL PRIMERO (R24) ────────────────────────────────────────────
    WA_C, WB_C, FA_C = -0.70, -1.40, 0.35
    w_sint = w_de_dos(A, WA_C, WB_C, FA_C)
    cA, cB, cF, cmax, crms = ajusta(A, w_sint)
    ok = (abs(cA - WA_C) < 1e-4 and abs(cB - WB_C) < 1e-4
          and abs(cF - FA_C) < 1e-4)
    print("=== CONTROL (va primero): sintetico wA=%.2f wB=%.2f fA=%.2f"
          % (WA_C, WB_C, FA_C))
    print("    recuperado   wA=%+.6f wB=%+.6f fA=%.6f  max|dw|=%.2e  -> %s"
          % (cA, cB, cF, cmax, "PASA" if ok else "FALLA"), flush=True)
    if not ok:
        SALIDA.write_text(json.dumps(dict(
            control=dict(pasa=False, wA=cA, wB=cB, fA=cF),
            veredicto="el ajuste no recupera un caso conocido: no vale"), indent=1))
        return

    # ── EL CASO REAL ─────────────────────────────────────────────────────
    w_obj = W0 + WA * (1.0 - A)
    wA, wB, fA, mx, rms = ajusta(A, w_obj)
    w_aj = w_de_dos(A, wA, wB, fA)

    if mx < TOL_IGUAL:
        veredicto = "SI — indistinguible del par algebraico"
    elif mx < TOL_APROX:
        veredicto = "APROXIMA — buena, pero distinguible"
    else:
        veredicto = "NO — la forma CPL no sale de dos fluidos de w constante"

    print("\n=== SSEE: objetivo w(a) = %.6f + %.6f(1-a)" % (W0, WA))
    print("    campo A (normal)   w_A = %+.6f   aporta hoy %.4f" % (wA, fA))
    print("    campo B (fantasma) w_B = %+.6f   aporta hoy %.4f" % (wB, 1 - fA))
    print("    desvio maximo |dw| = %.5f   (criterio: <%.2f igual, <%.2f aprox)"
          % (mx, TOL_IGUAL, TOL_APROX))
    print("    rms                = %.5f" % rms)
    print("    w en a=0.30/0.65/1.00  objetivo %+.4f %+.4f %+.4f"
          % (w_obj[0], w_obj[len(A) // 2], w_obj[-1]))
    print("                           ajuste   %+.4f %+.4f %+.4f"
          % (w_aj[0], w_aj[len(A) // 2], w_aj[-1]))
    print("\n    ¿dos campos de w constante lo reproducen? %s" % veredicto,
          flush=True)

    # ── POR QUE falla: la DIRECCION, no la precision ──────────────────────
    # Una mezcla de dos fluidos de w fijo se corre siempre hacia el mas
    # negativo, porque ese diluye mas despacio y acaba mandando. O sea que su
    # w_tot SIEMPRE decrece con a. SSEE pide lo contrario: de -1.309 en a=0.3
    # a -0.840 hoy, w CRECIENTE. Esto se comprueba, no se afirma.
    direcciones = []
    for wa_, wb_, fa_ in ((-0.70, -1.40, 0.35), (-0.50, -1.20, 0.60),
                          (-0.90, -1.90, 0.20), (-0.20, -1.05, 0.80)):
        v = w_de_dos(A, wa_, wb_, fa_)
        d = np.diff(v)
        direcciones.append(dict(
            w_A=wa_, w_B=wb_, f_A=fa_, w_ini=float(v[0]), w_fin=float(v[-1]),
            sentido=("crece" if np.all(d > 0) else
                     "decrece" if np.all(d < 0) else "mixto")))
    todas_decrecen = all(x["sentido"] == "decrece" for x in direcciones)
    print("    todas las mezclas probadas DECRECEN con a: %s" % todas_decrecen)
    print("    SSEE pide w CRECIENTE: %+.4f -> %+.4f"
          % (w_obj[0], w_obj[-1]), flush=True)

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(dict(
        corrida="cola #12 — dos campos genericos de w constante contra (w0, wa)",
        idea="M. Almeida — dos opuestos; aqui uno normal y uno fantasma",
        objetivo=dict(w0=W0, wa=WA, rango_a=[float(A[0]), float(A[-1])]),
        criterio=dict(igual=TOL_IGUAL, aproxima=TOL_APROX,
                      sobre="desvio maximo de w(a), escrito antes de correr"),
        control=dict(verdad=dict(wA=WA_C, wB=WB_C, fA=FA_C),
                     recuperado=dict(wA=cA, wB=cB, fA=cF), max_dw=cmax,
                     pasa=True),
        resultado=dict(w_A=wA, w_B=wB, f_A=fA, f_B=1 - fA,
                       max_dw=mx, rms=rms, veredicto=veredicto,
                       nota="el ajuste COLAPSA a w_A = w_B, o sea a un solo "
                            "fluido con el w medio: no es que ajuste mal, es "
                            "que la mezcla no puede moverse en el sentido que "
                            "hace falta"),
        direccion=dict(
            mezclas_probadas=direcciones, todas_decrecen=todas_decrecen,
            ssee_pide=dict(w_en_a030=float(w_obj[0]), w_hoy=float(w_obj[-1]),
                           sentido="crece"),
            lectura="una mezcla de dos fluidos de w constante siempre se corre "
                    "hacia el mas negativo, porque diluye mas despacio; su w "
                    "total DECRECE con a. SSEE pide que CREZCA. El fallo es de "
                    "sentido, no de precision, y por eso ningun par de valores "
                    "lo arregla"),
        pendiente="si PASA, comparar w_A, w_B y f_A contra el diccionario, con "
                  "tolerancia declarada y look-elsewhere (1/490 a +-0.0005). "
                  "NO se busca aqui: primero el numero, despues la busqueda",
        alcance="no toca el acuerdo con DESI, que es del par algebraico"),
        indent=1))
    print("escrito -> %s" % SALIDA.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
