"""Si los dos campos SE PASAN ENERGIA, ¿sale el (w0, wa) algebraico? (cola #14)

DE DONDE VIENE. Dos corridas anteriores cerraron la via sin acoplar:
  · `dos_campos_phi_pi.py` — el campo de Paper 7 mas un compañero: el de Paper 7
    ESTORBA, y pasado f1 = 0.08 el sobrante pediria densidad negativa.
  · `dos_campos_genericos.py` — un campo normal mas uno fantasma, los dos de w
    constante: el ajuste COLAPSA a w_A = w_B. Y la razon es estructural: en una
    mezcla de w fijos manda al final el que diluye mas despacio, el mas
    negativo, asi que w_tot SIEMPRE DECRECE. SSEE lo pide CRECIENTE
    (-1.3089 -> -0.8399). Es fallo de SENTIDO, no de precision.

Queda una sola puerta: que los dos campos se pasen energia. Si el flujo va
hacia el componente MENOS negativo, su peso crece con el tiempo y w_tot puede
SUBIR. Eso es lo que aqui se mide.

QUE SE INTEGRA. Con N = ln a y un intercambio Q = lam · H · rho_tot:
    drho_1/dN = -3(1+w_1) rho_1 + lam · rho_tot
    drho_2/dN = -3(1+w_2) rho_2 - lam · rho_tot
El total sigue conservandose, asi que sigue siendo UNA energia oscura; lo unico
que cambia es como se reparte por dentro. Se compara
    w_tot(a) = (w_1 rho_1 + w_2 rho_2)/rho_tot
contra w(a) = w0 + wa(1-a) con los DOS numeros algebraicos.

LA PREGUNTA QUE DE VERDAD IMPORTA (la hizo Mike): ¿hace falta un campo FANTASMA
de verdad, o basta con que se toquen? Un fantasma real tiene energia cinetica
negativa y arrastra una inestabilidad conocida (ghost), asi que no tenerlo es
mejor fisica. Por eso se corren DOS escenarios:
    LIBRE      w_1, w_2 ∈ [-3, +1]   — puede salir fantasma si le conviene
    SIN_GHOST  w_1, w_2 ∈ [-1, +1]   — NINGUNO puede ser fantasma
Si SIN_GHOST tambien cumple, la respuesta es que el acoplamiento hace el trabajo
y el fantasma NO hace falta.

CRITERIO, escrito ANTES de correr, sobre el desvio MAXIMO en a ∈ [0.30, 1.00]:
    max|Δw| < 0.01   -> CUMPLE
    0.01 .. 0.05     -> aproxima, distinguible
    > 0.05           -> NO

CONTROL (R53, y PRIMERO por R24). Dos controles, no uno:
  (a) con lam = 0 el integrador tiene que reproducir la formula cerrada de dos
      fluidos sin acoplar, la misma que uso la corrida #12;
  (b) se fabrica un total con un acoplamiento CONOCIDO (w1=-0.60, w2=-1.30,
      lam=0.25, f1=0.40) y el ajuste tiene que devolver esos cuatro numeros.
Si cualquiera falla, no se mide nada y se aborta.

NO SE BUSCA NADA EN EL DICCIONARIO AQUI. Primero el numero, despues la busqueda,
y con el riesgo look-elsewhere declarado (490 razones: 1/490 a +-0.0005).

FUENTE: results/logs/eft_dos_campos_acoplados.json
"""
# ORIGEN-VALOR: 0.0005 — tolerancia DECLARADA de la busqueda look-elsewhere (1 de 490 razones), no es medida
import json
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
import ssee_core as S                                   # noqa: E402

SALIDA = REPO / "results" / "logs" / "eft_dos_campos_acoplados.json"

W0, WA = S.W0, S.WA
A = np.linspace(0.30, 1.00, 400)
N = np.log(A)
TOL_CUMPLE, TOL_APROX = 0.01, 0.05     # criterio, antes de correr


def evoluciona(w1, w2, lam, f1, n_grid=N):
    """rho_1(a), rho_2(a) hacia atras desde hoy. rho_tot(1) = 1."""
    def rhs(n, y):
        rt = y[0] + y[1]
        return [-3.0 * (1.0 + w1) * y[0] + lam * rt,
                -3.0 * (1.0 + w2) * y[1] - lam * rt]
    s = solve_ivp(rhs, (0.0, float(n_grid.min())), [f1, 1.0 - f1],
                  rtol=1e-11, atol=1e-14, dense_output=True)
    if not s.success:
        return None, None
    y = s.sol(n_grid)
    return y[0], y[1]


def w_total(w1, w2, lam, f1, n_grid=N):
    r1, r2 = evoluciona(w1, w2, lam, f1, n_grid)
    if r1 is None:
        return None
    rt = r1 + r2
    # exigencia fisica: las dos densidades positivas en todo el rango
    if np.any(r1 <= 0) or np.any(r2 <= 0) or np.any(~np.isfinite(rt)):
        return None
    return (w1 * r1 + w2 * r2) / rt


def ajusta(w_obj, w_min):
    """Ajusta (w1, w2, lam, f1). w_min = -3 libre, -1 sin fantasmas."""
    def res(u):
        v = w_total(u[0], u[1], u[2], u[3])
        return np.full(len(N), 1e3) if v is None else v - w_obj

    mejor = None
    for x0 in ([-0.6, -1.3, 0.25, 0.40], [-0.9, -1.1, 0.10, 0.60],
               [-0.4, -1.5, 0.50, 0.30], [-1.0, -0.8, -0.20, 0.50],
               [-0.7, -0.95, 0.05, 0.70], [-0.99, -0.99, 0.30, 0.50]):
        u0 = [max(x0[0], w_min), max(x0[1], w_min), x0[2], x0[3]]
        try:
            r = least_squares(res, u0,
                              bounds=([w_min, w_min, -5.0, 1e-4],
                                      [1.0, 1.0, 5.0, 1.0 - 1e-4]),
                              xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=8000)
        except Exception:
            continue
        if mejor is None or r.cost < mejor.cost:
            mejor = r
    d = res(mejor.x)
    return dict(w_1=float(mejor.x[0]), w_2=float(mejor.x[1]),
                lam=float(mejor.x[2]), f_1=float(mejor.x[3]),
                max_dw=float(np.max(np.abs(d))),
                rms=float(np.sqrt(np.mean(d ** 2))))


def veredicto_de(mx):
    if mx < TOL_CUMPLE:
        return "CUMPLE"
    if mx < TOL_APROX:
        return "aproxima, distinguible"
    return "NO"


def main():
    # ── CONTROL (a): lam = 0 tiene que dar la formula cerrada ────────────
    w1, w2, f1 = -0.70, -1.40, 0.35
    v = w_total(w1, w2, 0.0, f1)
    rA = f1 * A ** (-3.0 * (1.0 + w1))
    rB = (1.0 - f1) * A ** (-3.0 * (1.0 + w2))
    cerrada = (w1 * rA + w2 * rB) / (rA + rB)
    d_a = float(np.max(np.abs(v - cerrada)))
    ok_a = d_a < 1e-8
    print("=== CONTROL (a): con lam=0 el integrador vs la formula cerrada")
    print("    desvio maximo = %.2e   -> %s"
          % (d_a, "PASA" if ok_a else "FALLA"), flush=True)

    # ── CONTROL (b): recuperar un acoplamiento CONOCIDO ──────────────────
    V = dict(w_1=-0.60, w_2=-1.30, lam=0.25, f_1=0.40)
    w_sint = w_total(V["w_1"], V["w_2"], V["lam"], V["f_1"])
    rec = ajusta(w_sint, -3.0)
    ok_b = all(abs(rec[k] - V[k]) < 1e-3 for k in V)
    print("=== CONTROL (b): recuperar w1=%.2f w2=%.2f lam=%.2f f1=%.2f"
          % (V["w_1"], V["w_2"], V["lam"], V["f_1"]))
    print("    recuperado    w1=%+.5f w2=%+.5f lam=%+.5f f1=%.5f  -> %s"
          % (rec["w_1"], rec["w_2"], rec["lam"], rec["f_1"],
             "PASA" if ok_b else "FALLA"), flush=True)

    if not (ok_a and ok_b):
        SALIDA.write_text(json.dumps(dict(
            control=dict(lam_cero=dict(desvio=d_a, pasa=ok_a),
                         recupera_conocido=dict(recuperado=rec, pasa=ok_b)),
            veredicto="el metodo no pasa sus controles: no mide nada"), indent=1))
        return

    # ── EL CASO REAL, dos escenarios ─────────────────────────────────────
    w_obj = W0 + WA * (1.0 - A)
    esc = {}
    for nombre, wmin in (("LIBRE", -3.0), ("SIN_GHOST", -1.0)):
        r = ajusta(w_obj, wmin)
        r["veredicto"] = veredicto_de(r["max_dw"])
        r["hay_fantasma"] = bool(min(r["w_1"], r["w_2"]) < -1.0 + 1e-6)
        esc[nombre] = r
        print("\n=== %s  (w minimo permitido: %+.1f)" % (nombre, wmin))
        print("    w_1 = %+.6f     w_2 = %+.6f" % (r["w_1"], r["w_2"]))
        print("    lam = %+.6f     f_1 = %.6f  (aporte de 1 hoy)"
              % (r["lam"], r["f_1"]))
        print("    desvio maximo = %.5f   rms = %.5f   -> %s"
              % (r["max_dw"], r["rms"], r["veredicto"]))
        print("    ¿alguno es fantasma (w < -1)? %s"
              % ("SI" if r["hay_fantasma"] else "NO"), flush=True)

    hace_falta = (esc["LIBRE"]["veredicto"] == "CUMPLE"
                  and esc["SIN_GHOST"]["veredicto"] != "CUMPLE")
    print("\n=== ¿HACE FALTA UN CAMPO FANTASMA DE VERDAD?")
    if esc["SIN_GHOST"]["veredicto"] == "CUMPLE":
        print("    NO — con los dos campos por encima de -1 ya cumple: el")
        print("    trabajo lo hace el ACOPLAMIENTO, no un fantasma.", flush=True)
    elif esc["LIBRE"]["veredicto"] == "CUMPLE":
        print("    SI — solo cumple dejando que uno baje de -1.", flush=True)
    else:
        print("    la pregunta no se contesta: ningun escenario cumple.",
              flush=True)

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(dict(
        corrida="cola #14 — dos campos ACOPLADOS contra el (w0, wa) algebraico",
        idea="M. Almeida — dos opuestos que se TOCAN, no que solo conviven",
        intercambio="Q = lam · H · rho_tot; el total se conserva",
        objetivo=dict(w0=W0, wa=WA, rango_a=[float(A[0]), float(A[-1])]),
        criterio=dict(cumple=TOL_CUMPLE, aproxima=TOL_APROX,
                      sobre="desvio maximo de w(a), escrito antes de correr"),
        control=dict(lam_cero=dict(desvio=d_a, pasa=ok_a),
                     recupera_conocido=dict(verdad=V, recuperado=rec, pasa=ok_b)),
        escenarios=esc,
        fantasma_necesario=hace_falta,
        pendiente="si algun escenario CUMPLE, comparar lam, w_1, w_2 y f_1 "
                  "contra el diccionario, con tolerancia declarada y "
                  "look-elsewhere (1/490 a +-0.0005). NO se busca aqui",
        alcance="no toca el acuerdo con DESI, que es del par algebraico"),
        indent=1))
    print("\nescrito -> %s" % SALIDA.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
