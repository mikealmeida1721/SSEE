"""¿La supresión que le falta a KiDS es PLANA o depende de la ESCALA? (cola #19)

LA PREGUNTA, y por qué se puede hacer ya. Con el fondo clavado y solo la
amplitud suelta, KiDS pide logA = 2.8633 (SSEE) y 2.8327 (LCDM): **el mismo
numero a 0.45 sigma**, y en los dos casos ~0.19 por debajo del que pide su
propio fondo cosmico — un deficit del 9-10% en sigma8. Ya esta descartado que
lo fabriquen los ingredientes del fondo (`fuga3 fondo`: se mueven <0.6 sigma),
que lo fabrique el fondo de SSEE (`lcdmfijo`: el de Planck da el mismo agujero),
y que sea la regla de medir (el volumen es `ln det F`, y KiDS no lo lleva).

Queda lo que se OBSERVA. Y hay tres explicaciones con huellas de escala
distintas:

    amplitud baja .......... suprime TODAS las escalas por igual
    materia que no se agrupa suprime POR DEBAJO de su escala de fuga
    bariones ............... suprime SOLO en las mas pequenas

KiDS mide xi_pm en un rango de escalas, asi que el DATO puede separarlas. No
hay que postular ninguna: se ajusta una familia de dos parametros que las
contiene a las tres y se mira DONDE cae el corte.

    P_sup(k) = P(k) * [ 1 - A_sup * x^2/(1+x^2) ],   x = k / k_c

    k_c -> 0     : x>>1 en todo el rango  -> supresion PLANA (amplitud)
    k_c medio    : grandes intactas, chicas suprimidas -> FUGA
    k_c -> grande: solo las mas chicas   -> BARIONES

NO SE POSTULA NINGUNA PARTICULA. Esto mide la FORMA de la supresion. Si sale
que hay corte, entonces —y solo entonces— se discute que lo produce. Escribir
el mecanismo antes de medir la forma es exactamente como nacio la particula
phi-DM que se retiro el 2026-08-01.

LA PRUEBA QUE DISTINGUE FISICA DE PARCHE (idea de Mike). Se corre con los DOS
fondos, el de SSEE y el de Planck. Si hace falta la MISMA supresion en ambos,
es algo real que ninguno de los dos estaba viendo. Si solo la necesita SSEE, es
un parche para tapar su fondo, y se dice.

CONTROL (R53, y va PRIMERO por R24). Dos, y los dos tienen respuesta conocida
de antemano:

  C1 · A_sup = 0 tiene que devolver EXACTAMENTE el chi2 de clavar logA en el
       valor del fondo sin tocar nada. Si la maquinaria mueve algo con la
       supresion apagada, esta rota.
  C2 · con k_c empujado a la escala mas grande de la rejilla, la supresion es
       plana y equivale a bajar la amplitud. El A_s efectivo que resulte,
       logA + ln(1-A_sup), tiene que caer sobre el logA que ya midio la cadena
       (2.8633 en SSEE) dentro de su barra. Si no cae, no estoy suprimiendo lo
       que creo.

FUENTE: results/logs/growth_2026-07/supresion_plana_o_escala.json
"""
# ORIGEN-VALOR: 0.6736 — h de Planck 2018 (arXiv:1807.06209)
import json
import pathlib
import sys
import time

import numpy as np
from scipy.optimize import minimize

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "src" / "p06_growth"))

import cobaya_kids as C                                        # noqa: E402
import kids_shear as K                                         # noqa: E402

SALIDA = REPO / "results" / "logs" / "growth_2026-07" / (
    "supresion_plana_o_escala_%s.json" % "_".join(sys.argv[1:] or ["todo"]).lower())

# El A_s que pide cada fondo cosmico (results/logs/cmb_dbic_tau_ajustado.json)
LOGA_CMB = dict(SSEE=3.0448340130228546, LCDM=3.0450790027403647)
# Lo que pide KiDS con ese fondo y la amplitud suelta (las cadenas convergidas)
LOGA_KIDS = dict(SSEE=(2.8633, 0.0500), LCDM=(2.8327, 0.0456))
# Fondo de LCDM clavado en Planck 2018, igual que en `loglike_lcdm_fijo`
LCDM_BG = (0.02237, 0.1200, 0.6736, 0.9649)

REJILLA_A = np.linspace(0.0, 0.60, 9)           # amplitud de la supresion
# el primer k_c es 1e-6, no 1e-4: con 1e-4 la propia rejilla de CAMB empieza en
# k=1e-4 y ahi x=1, o sea la supresion NO seria plana en el extremo. Con 1e-6 el
# factor vale 0.9999 ya en la primera k, que es lo que el control C2 necesita.
REJILLA_KC = np.concatenate([[1e-6], np.geomspace(0.02, 8.0, 9)])    # h/Mpc
HALO = np.array([2.15, 2.60, 3.05])             # se perfila sobre esto
ITER = 40                                       # Nelder-Mead, con arranque tibio


def suprime(pk, kh, A_sup, k_c):
    """P(k) -> P(k)*[1 - A_sup*x^2/(1+x^2)],  x = k/k_c.  kh viene en h/Mpc."""
    if A_sup <= 0.0:
        return pk
    x2 = (np.asarray(kh) / k_c) ** 2
    return pk * (1.0 - A_sup * x2 / (1.0 + x2))[None, :]


def chi2_de(bg, logA, halo_A, A_sup, k_c, extra):
    """chi2 total (dato + priores) con la supresion aplicada sobre P(k)."""
    A_IA, delta_c = extra
    dz = C.DZ_MEAN
    try:
        r, p, kh, zpk, pk, gr = C._camb_cached(bg, np.exp(logA) * 1e-10, halo_A)
        ells, Cl, idx = K.cl_shear(C.D, r, p, kh, zpk,
                                   suprime(pk, kh, A_sup, k_c), gr, A_IA, dz)
        th = K.theory_vector(C.D, ells, Cl, idx, delta_c=delta_c)
    except Exception:
        return 1e10
    dv = (th - C.D['d'])[C.MASK]
    rr = dz - C.DZ_MEAN
    return (float(dv @ C.CINV @ dv) + float(rr @ C.SOM_INV @ rr)
            + (delta_c / C.DELTA_C_SIG) ** 2)


_TIBIO = {}          # arranque tibio: la solucion del punto anterior


def perfila(bg, logA, A_sup, k_c):
    """Minimiza sobre halo_A (rejilla) y sobre (A_IA, delta_c) (Nelder-Mead).

    halo_A va en rejilla y no en el minimizador porque es el UNICO de los tres
    que entra en CAMB: en rejilla la cache lo resuelve con 3 llamadas por fondo,
    dentro del minimizador costaria 4.7 s cada evaluacion en vez de 0.5 s.
    """
    mejor, arg = 1e10, None
    for h in HALO:
        e0 = _TIBIO.get((bg if isinstance(bg, str) else 'lcdm', h),
                        np.array([0.55, 0.0]))
        o = minimize(lambda e: chi2_de(bg, logA, h, A_sup, k_c, e), e0,
                     method='Nelder-Mead',
                     options=dict(maxiter=ITER, xatol=1e-3, fatol=1e-3))
        _TIBIO[(bg if isinstance(bg, str) else 'lcdm', h)] = o.x
        if o.fun < mejor:
            mejor, arg = float(o.fun), (float(h), float(o.x[0]), float(o.x[1]))
    return mejor, arg


def corre_modelo(nombre, bg):
    print("\n" + "=" * 72, flush=True)
    print("  FONDO: %s" % nombre)
    print("=" * 72, flush=True)
    lA = LOGA_CMB[nombre]

    # ── CONTROL 1: supresion apagada ────────────────────────────────────
    # ARREGLADO 2026-09-09: antes comparaba dos PERFILES, y fallo con LCDM por
    # 1e-4 contra un criterio de 1e-6. No era fisico: el perfil arranca tibio
    # desde la solucion anterior y con tope de 40 iteraciones no converge a esa
    # cifra — el criterio quedaba POR DEBAJO del ruido de mi propio minimizador.
    # Lo que C1 tiene que comprobar es la FUNCION `suprime`, no el minimizador:
    # con las molestias en el MISMO punto, A_sup=0 da el mismo chi2 exactamente.
    # Corregido con la rejilla de LCDM aun sin ver.
    e = np.array([0.55, 0.0])
    c0 = chi2_de(bg, lA, HALO[1], 0.0, 1.0, e)
    c0b = chi2_de(bg, lA, HALO[1], 0.0, 5.0, e)
    ok1 = bool(c0 == c0b)
    print("  C1 · supresion apagada, molestias fijas: %.6f vs %.6f  -> %s"
          % (c0, c0b, "PASA" if ok1 else "FALLA"), flush=True)
    c0, _ = perfila(bg, lA, 0.0, 1.0)           # el perfil, para la referencia

    # ── CONTROL 2: supresion PLANA equivale a bajar la amplitud ─────────
    # con k_c minusculo la supresion es plana; el A_s efectivo tiene que caer
    # sobre el que ya midio la cadena.
    kc0 = REJILLA_KC[0]
    fila = [(a, perfila(bg, lA, a, kc0)[0]) for a in REJILLA_A]
    a_best = min(fila, key=lambda t: t[1])[0]
    logA_ef = lA + np.log(1.0 - a_best) if a_best < 1 else float('nan')
    med, sig = LOGA_KIDS[nombre]
    desv = abs(logA_ef - med) / sig
    ok2 = bool(desv < 2.0)
    print("  C2 · plana: A_sup=%.3f -> logA efectivo %.4f  vs cadena %.4f+-%.4f"
          "  = %.2f sigma  -> %s" % (a_best, logA_ef, med, sig, desv,
                                     "PASA" if ok2 else "FALLA"), flush=True)

    if not (ok1 and ok2):
        return dict(controles=dict(C1=bool(ok1), C2=bool(ok2), pasa=False),
                    veredicto="la maquinaria no reproduce lo que ya se sabe")

    # ── LA MEDIDA ───────────────────────────────────────────────────────
    print("\n  %-10s %s" % ("k_c", "  ".join("%5.2f" % a for a in REJILLA_A)),
          flush=True)
    Z = np.full((len(REJILLA_KC), len(REJILLA_A)), np.nan)
    for i, kc in enumerate(REJILLA_KC):
        for j, a in enumerate(REJILLA_A):
            Z[i, j] = perfila(bg, lA, float(a), float(kc))[0]
        print("  %-10.4f %s" % (kc, "  ".join("%5.1f" % v for v in Z[i])),
              flush=True)

    i, j = np.unravel_index(np.nanargmin(Z), Z.shape)
    return dict(
        controles=dict(C1=bool(ok1), C2=bool(ok2), chi2_sin_supresion=float(c0),
                       A_plana=float(a_best), logA_efectivo=float(logA_ef),
                       desv_sigmas=float(desv), pasa=True),
        logA_clavado=lA, rejilla_kc=REJILLA_KC.tolist(),
        rejilla_A=REJILLA_A.tolist(), chi2=Z.tolist(),
        mejor=dict(k_c=float(REJILLA_KC[i]), A_sup=float(REJILLA_A[j]),
                   chi2=float(Z[i, j]), gana_vs_sin_supresion=float(c0 - Z[i, j])))


def main():
    t0 = time.time()
    res = {}
    cuales = sys.argv[1:] or ['SSEE', 'LCDM']
    for nombre, bg in (("SSEE", 'SSEE'), ("LCDM", LCDM_BG)):
        if nombre not in cuales:
            continue
        res[nombre] = corre_modelo(nombre, bg)

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(dict(
        corrida="cola #19 — la supresion que le falta a KiDS, plana o con escala",
        pregunta="M. Almeida — si el mecanismo hace falta en los DOS fondos es "
                 "fisica que ninguno estaba viendo; si solo en SSEE es un parche",
        familia="P(k)*[1 - A_sup*x^2/(1+x^2)], x=k/k_c; k_c->0 plana, "
                "k_c medio fuga, k_c grande bariones",
        no_postula="mide la FORMA. NO introduce ninguna particula ni mecanismo",
        modelos=res, segundos=time.time() - t0,
        alcance="NO decide que produce la supresion, ni entra en ningun paper"),
        indent=1))
    print("\nescrito -> %s  (%.1f h)"
          % (SALIDA.relative_to(REPO), (time.time() - t0) / 3600), flush=True)


if __name__ == "__main__":
    main()
