"""¿Hacen falta DOS campos, y el segundo esta en el algebra? (idea de M. Almeida)

DE DONDE SALE. Un solo campo escalar minimamente acoplado NO PUEDE cruzar
w = -1 (barrera fantasma). SSEE la cruza: con (w0, wa) = (-0.839950, -0.669975)
el cruce cae en z = 0.3139. Y el lagrangiano de Paper 7, k-esencia
K = c1 X + c2 X^2, da w0 EXACTO pero wa = +0.4135, signo contrario, porque su
punto de reposo es w = -1 y el campo esta CAYENDO hacia el.

La salida conocida es quintom: dos campos, uno normal y uno fantasma. Encaja
con el principio de Mike (dos opuestos que crean lo que ninguno hace solo), y
encaja con que su gramatica tenga dos leyes, phi (la copia) y pi (la no
auto-suma).

EL PELIGRO, y por eso el diseno es este. Con dos campos LIBRES se ajusta
cualquier par (w0, wa), asi que la prueba obvia no puede fallar y no demuestra
nada. Aqui NO se ajusta el resultado: se DESPEJA lo que el segundo componente
tendria que ser, y despues se mira si eso esta en el algebra.

COMO. Cero libertad en el objetivo:
  · el TOTAL tiene que dar w(a) = w0 + wa(1-a) con los dos numeros algebraicos;
  · el componente 1 es el campo de Paper 7, con u fijado por el algebra: su
    rho_1(a) sale de integrar du/dN = -6u(1+2u)/(1+6u), sin ningun ajuste;
  · el componente 2 es lo que sobra: rho_2 = rho_tot - rho_1, y su presion
    p_2 = p_tot - p_1. Su ecuacion de estado queda DESPEJADA, no elegida.

El unico numero libre es f1, la fraccion que el campo 1 aporta hoy. Se fija
exigiendo lo mas simple que puede ser un segundo campo: **que su w sea
CONSTANTE**. Eso es una condicion sobre todo el rango, no un punto, asi que
puede no tener solucion — y ahi es donde la prueba puede fallar.

QUE SE MIRA DESPUES. El f1 que resulte y el w2 que resulte, contra el
diccionario. Se reporta la tolerancia y el riesgo de acierto por azar (490
razones: 1/490 a +-0.0005).

CONTROL (R53, y va PRIMERO por R24). Se fabrica un total SINTETICO que es, por
construccion, el campo 1 con f1 = 0.6 mas un segundo componente de w constante
= -1.3. El despeje tiene que devolver esos dos numeros. Si no los devuelve, el
metodo no despeja nada y se aborta.

FUENTE: results/logs/eft_dos_campos_phi_pi.json
"""
# ORIGEN: results/logs/eft_dos_campos_phi_pi.json   (cruce_z 0.3139, wa_que_da 0.4135)
# ORIGEN-VALOR: 0.0005 — tolerancia DECLARADA para comparar contra el diccionario de 490 razones (la eleccion fija el riesgo 1/490)
import json
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
import ssee_core as S                                   # noqa: E402

SALIDA = REPO / "results" / "logs" / "eft_dos_campos_phi_pi.json"

W0, WA = S.W0, S.WA
U0 = -0.522735                     # c2 X/c1 fijado por el algebra (Paper 7)
A = np.linspace(0.30, 1.00, 400)   # rango donde se exige la constancia


def campo1(a_grid, u_hoy=U0):
    """rho_1(a) y p_1(a) del k-esencia de Paper 7, normalizados a rho_1(1)=1."""
    f = lambda N, y: [-6 * y[0] * (1 + 2 * y[0]) / (1 + 6 * y[0])]
    N = np.log(a_grid)
    s = solve_ivp(f, (0.0, N.min()), [u_hoy], rtol=1e-11, atol=1e-14,
                  dense_output=True)
    u = s.sol(N)[0]
    rho = u * (1 + 3 * u)          # rho = c1^2 u(1+3u)/c2, factor comun fuera
    p = u * (1 + u)                # p   = c1^2 u(1+u)/c2
    n = u_hoy * (1 + 3 * u_hoy)
    return rho / n, p / n


def total_cpl(a_grid, w0, wa):
    """rho_DE(a)/rho_DE(1) y p/rho_DE(1) para w(a) = w0 + wa(1-a)."""
    rho = a_grid ** (-3.0 * (1.0 + w0 + wa)) * np.exp(-3.0 * wa * (1.0 - a_grid))
    w = w0 + wa * (1.0 - a_grid)
    return rho, w * rho


def despeja(a_grid, r_tot, p_tot, r1, p1, f1):
    """Componente 2 = lo que sobra. Devuelve (w2(a), rho2(a))."""
    r2 = r_tot - f1 * r1
    p2 = p_tot - f1 * p1
    with np.errstate(divide="ignore", invalid="ignore"):
        w2 = np.where(np.abs(r2) > 1e-12, p2 / r2, np.nan)
    return w2, r2


def dispersion(f1, a_grid, r_tot, p_tot, r1, p1):
    """Cuanto se aparta de constante el w2 despejado (desviacion tipica)."""
    w2, r2 = despeja(a_grid, r_tot, p_tot, r1, p1, f1)
    if not np.all(np.isfinite(w2)) or np.any(r2 <= 0):
        return 1e6
    return float(np.std(w2))


def resuelve(a_grid, r_tot, p_tot, r1, p1):
    """Barrido grueso + refinado local.

    FIX 2026-09-08: aqui habia un `minimize_scalar` acotado, y FALLABA su
    propio control: el minimo es muy estrecho (dispersion 3e-1 alrededor,
    6e-16 dentro) y ademas hay un acantilado donde rho_2 se vuelve negativa y
    la funcion salta a 1e6. Brent acotado se iba al borde. Un barrido no se
    pierde en un valle estrecho, y el refinado local le da la precision.
    """
    g = np.linspace(1e-4, 0.9999, 2000)
    v = np.array([dispersion(x, a_grid, r_tot, p_tot, r1, p1) for x in g])
    i = int(np.argmin(v))
    lo = g[max(i - 1, 0)]
    hi = g[min(i + 1, len(g) - 1)]
    o = minimize_scalar(dispersion, bounds=(lo, hi), method="bounded",
                        args=(a_grid, r_tot, p_tot, r1, p1),
                        options=dict(xatol=1e-12))
    f1 = float(o.x) if o.fun <= v[i] else float(g[i])
    w2, r2 = despeja(a_grid, r_tot, p_tot, r1, p1, f1)
    return f1, float(np.mean(w2)), float(np.std(w2))


def main():
    r1, p1 = campo1(A)

    # ── CONTROL PRIMERO (R24): total sintetico con f1 y w2 conocidos ─────
    F1_C, W2_C = 0.6, -1.3
    r2c = A ** (-3.0 * (1.0 + W2_C))
    r2c = r2c / r2c[-1] * (1.0 - F1_C)
    rt_c = F1_C * r1 + r2c
    pt_c = F1_C * p1 + W2_C * r2c
    f1c, w2c, sdc = resuelve(A, rt_c, pt_c, r1, p1)
    ok = abs(f1c - F1_C) < 1e-3 and abs(w2c - W2_C) < 1e-3
    print("=== CONTROL (va primero): total sintetico f1=%.3f  w2=%.3f" % (F1_C, W2_C))
    print("    despejado    f1=%.6f  w2=%.6f  dispersion=%.2e  -> %s"
          % (f1c, w2c, sdc, "PASA" if ok else "FALLA"), flush=True)
    if not ok:
        SALIDA.write_text(json.dumps(dict(
            control=dict(pasa=False, f1=f1c, w2=w2c),
            veredicto="el despeje no recupera un caso conocido: no vale"), indent=1))
        return

    # ── EL CASO REAL ─────────────────────────────────────────────────────
    rt, pt = total_cpl(A, W0, WA)
    f1, w2, sd = resuelve(A, rt, pt, r1, p1)

    # Diagnostico: ¿el f1 esta en el borde porque el barrido falla, o porque
    # la funcion es monotona? Se guarda la forma para poder distinguirlo.
    forma = [dict(f1=float(x), dispersion=dispersion(float(x), A, rt, pt, r1, p1))
             for x in (1e-4, 0.01, 0.05, 0.10, 0.30, 0.60)]
    tope = None
    for x in np.linspace(1e-4, 0.999, 400):
        if np.any(despeja(A, rt, pt, r1, p1, float(x))[1] <= 0):
            tope = float(x)
            break
    w2v, r2v = despeja(A, rt, pt, r1, p1, f1)
    print("\n=== SSEE: total = CPL(w0=%.6f, wa=%.6f)" % (W0, WA))
    print("    f1 (fraccion del campo de Paper 7 hoy) = %.6f" % f1)
    print("    f2 = 1 - f1                            = %.6f" % (1 - f1))
    print("    w2 medio del componente que sobra      = %+.6f" % w2)
    print("    dispersion de w2 (0 = constante)       = %.4e" % sd)
    print("    w2 en a=0.30 / a=0.65 / a=1.00         = %+.4f %+.4f %+.4f"
          % (w2v[0], w2v[len(A) // 2], w2v[-1]), flush=True)

    const = sd < 0.02
    print("\n    ¿el sobrante puede ser un campo de w CONSTANTE? %s"
          % ("SI (dispersion < 0.02)" if const else "NO — se mueve demasiado"),
          flush=True)

    out = dict(
        corrida="dos campos: se DESPEJA el segundo, no se ajusta",
        idea="M. Almeida — un solo campo no cruza w=-1; ¿el segundo esta en el algebra?",
        objetivo=dict(w0=W0, wa=WA, cruce_z=0.3139),
        campo1=dict(lagrangiano="k-esencia Paper 7, K = c1 X + c2 X^2",
                    u_algebraico=U0, w0_que_da=S.W0, wa_que_da=0.4135),
        control=dict(criterio="recuperar un total sintetico con f1=0.6, w2=-1.3",
                     f1=f1c, w2=w2c, dispersion=sdc, pasa=True),
        resultado=dict(f1=f1, f2=1 - f1, w2_medio=w2, w2_dispersion=sd,
                       w2_constante=bool(const),
                       w2_en=dict(a030=float(w2v[0]),
                                  a065=float(w2v[len(A) // 2]),
                                  a100=float(w2v[-1]))),
        diagnostico=dict(
            forma_dispersion=forma,
            f1_donde_rho2_se_vuelve_negativa=tope,
            lectura="la dispersion CRECE con f1 desde el borde: el campo de "
                    "Paper 7 no ayuda nada, y pasado ese tope el sobrante "
                    "tendria densidad negativa. El f1 en el borde no es un "
                    "fallo del barrido, es la respuesta"),
        pendiente="comparar f1, f2 y w2 contra el diccionario, con la "
                  "tolerancia declarada y el riesgo look-elsewhere (1/490 a "
                  "+-0.0005). NO se ha hecho aqui a proposito: primero el "
                  "numero, despues la busqueda",
        alcance="no toca el acuerdo con DESI, que es del par algebraico")
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(out, indent=1))
    print("\nescrito -> %s" % SALIDA.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
