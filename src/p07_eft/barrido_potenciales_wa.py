"""¿Existe un potencial V(phi) cuyo atractor de w_a = -0.669975? (OP-24)

QUE CONTESTA, y la distincion importa. El par (w0, wa) = (-0.840, -0.670) sale
del ALGEBRA de SSEE y ajusta el dato de DESI a 0.24 sigma. Eso esta medido y no
se discute aqui. Lo que se pregunta es otra cosa: **el propio modelo dice que la
energia oscura es un campo escalar** (Paper 7). Si se toma ese campo, se le pone
un potencial y se le deja evolucionar, ¿sale w_a = -0.670?

Hoy NO. Las tres formas probadas dan -0.093 (atractor exponencial), -0.211
(lambda = 1.0205) y +0.406 (acoplado, signo CONTRARIO). O sea que el modelo
concuerda con el dato pero **no consigo mismo**: el algebra dice una cosa y la
dinamica de su propio campo dice otra. Eso es lo que hay que cerrar, y no lo
cierra ningun dato nuevo.

QUE SE BARRE. Familias estandar de quintaesencia, cada una con su parametro:
  exp        V = V0 exp(-lambda phi)                  [Wetterich, Ferreira-Joyce]
  ley_pot    V = V0 phi^(-n)                          [Ratra-Peebles]
  doble_exp  V = V0 [exp(-a phi) + exp(-b phi)]       [Barreiro-Copeland-Nunes]
  pNGB       V = V0 [1 + cos(phi/f)]                  [Frieman et al., axion-like]
  colina     V = V0 [1 - (phi/mu)^2]                  [hilltop]

En cada punto se dispara V0 hasta que Omega_phi(hoy) = 0.691119 (= 1 - Om_m,CMB
algebraico). Despues se leen w0 = w(a=1) y wa = -dw/da|_{a=1}, que es la
definicion CPL con la que DESI reporta.

CRITERIO, escrito ANTES de correr: se busca |w0 + 0.839950| < 0.02 **y**
|wa + 0.669975| < 0.05 a la vez. Si ninguna familia entra, OP-24 se refuerza y
la conclusion es que la energia oscura de SSEE **no es un campo con potencial**
en estas familias — que es un resultado, no un fracaso.

CONTROL (R53). El exponencial con lambda -> 0 tiene que dar w = -1 exacto (es
el limite de constante cosmologica). Si no lo da, el integrador esta mal y el
barrido entero no vale. Se corre PRIMERO (R24) y aborta si falla.

FUENTE: results/logs/eft_barrido_potenciales_wa.json
"""
import json
import pathlib
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

REPO = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
import ssee_core as S                                   # noqa: E402

SALIDA = REPO / "results" / "logs" / "eft_barrido_potenciales_wa.json"

OM_M = S.OMEGA_M_CMB              # 0.308881, materia hoy
OM_R = 9.1e-5                     # radiacion hoy (fotones + neutrinos sin masa)
OM_DE = 1.0 - OM_M - OM_R         # lo que tiene que salir del disparo
W0_OBJ, WA_OBJ = S.W0, S.WA       # -0.839950, -0.669975
TOL_W0, TOL_WA = 0.02, 0.05       # criterio, escrito antes de correr

N_INI = -14.0                     # ln a inicial (z ~ 1.2e6)


def _pot(fam, par):
    """Devuelve (V(phi), dV/dphi) en unidades de V0 = 1."""
    if fam == "exp":
        lam = par
        return (lambda f: np.exp(-lam * f), lambda f: -lam * np.exp(-lam * f))
    if fam == "ley_pot":
        n = par
        return (lambda f: (f + 1.0) ** (-n),
                lambda f: -n * (f + 1.0) ** (-n - 1.0))
    if fam == "doble_exp":
        a, b = par, 0.1
        return (lambda f: np.exp(-a * f) + np.exp(-b * f),
                lambda f: -a * np.exp(-a * f) - b * np.exp(-b * f))
    if fam == "pNGB":
        fdec = par
        return (lambda f: 1.0 + np.cos(f / fdec),
                lambda f: -np.sin(f / fdec) / fdec)
    if fam == "colina":
        mu = par
        return (lambda f: max(1.0 - (f / mu) ** 2, 1e-12),
                lambda f: -2.0 * f / mu ** 2)
    raise ValueError(fam)


def evoluciona(fam, par, V0, N_fin=0.0, n_out=4000):
    """Integra el campo en e-folds. Devuelve (N, w_phi, Omega_phi)."""
    V, dV = _pot(fam, par)

    def rhs(N, y):
        f, fp = y                       # phi y dphi/dN
        a = np.exp(N)
        rm = 3.0 * OM_M * a ** -3       # 3 H0^2 Om (unidades M_pl=1, H0=1)
        rr = 3.0 * OM_R * a ** -4
        Vv = V0 * V(f)
        den = 3.0 - 0.5 * fp ** 2
        if den <= 1e-8:
            return [0.0, 0.0]
        H2 = (rm + rr + Vv) / den
        rf = 0.5 * H2 * fp ** 2 + Vv
        pf = 0.5 * H2 * fp ** 2 - Vv
        rt = rm + rr + rf
        pt = rr / 3.0 + pf
        HpH = -1.5 * (1.0 + pt / rt)
        return [fp, -(3.0 + HpH) * fp - V0 * dV(f) / H2]

    s = solve_ivp(rhs, (N_INI, N_fin), [1e-4, 0.0], rtol=1e-9, atol=1e-12,
                  dense_output=True, method="LSODA")
    if not s.success:
        return None
    N = np.linspace(N_INI, N_fin, n_out)
    f, fp = s.sol(N)
    a = np.exp(N)
    rm, rr = 3.0 * OM_M * a ** -3, 3.0 * OM_R * a ** -4
    Vv = V0 * np.array([V(x) for x in f])
    H2 = (rm + rr + Vv) / (3.0 - 0.5 * fp ** 2)
    rf = 0.5 * H2 * fp ** 2 + Vv
    pf = 0.5 * H2 * fp ** 2 - Vv
    return N, pf / rf, rf / (rm + rr + rf)


def dispara(fam, par):
    """Ajusta V0 para que Omega_phi(hoy) = OM_DE. Devuelve (w0, wa) o None."""
    def falta(lv):
        r = evoluciona(fam, par, 10.0 ** lv)
        return 1e3 if r is None else r[2][-1] - OM_DE
    try:
        lv = brentq(falta, -6.0, 3.0, xtol=1e-10)
    except Exception:
        return None
    r = evoluciona(fam, par, 10.0 ** lv)
    if r is None:
        return None
    N, w, om = r
    a = np.exp(N)
    # CPL alrededor de a=1: w0 = w(1), wa = -dw/da en a=1
    m = a > 0.9
    c = np.polyfit(a[m], w[m], 2)
    w0 = float(np.polyval(c, 1.0))
    wa = float(-np.polyval(np.polyder(c), 1.0))
    return w0, wa, float(10.0 ** lv), float(om[-1])


def main():
    print("=== CONTROL (R24, va primero): exponencial con lambda -> 0", flush=True)
    c = dispara("exp", 1e-4)
    if c is None:
        raise SystemExit("el control no integra: el barrido no vale")
    ok = abs(c[0] + 1.0) < 1e-3 and abs(c[1]) < 1e-2
    print("    w0 = %+.6f  wa = %+.6f   (debe ser -1 y 0)  -> %s"
          % (c[0], c[1], "PASA" if ok else "FALLA"), flush=True)
    if not ok:
        SALIDA.write_text(json.dumps(dict(
            control=dict(pasa=False, w0=c[0], wa=c[1]),
            veredicto="el integrador no reproduce el limite de constante "
                      "cosmologica: el barrido no vale"), indent=1))
        return

    FAM = {"exp": np.linspace(0.05, 1.6, 24),
           "ley_pot": np.linspace(0.1, 6.0, 20),
           "doble_exp": np.linspace(0.2, 4.0, 16),
           "pNGB": np.linspace(0.3, 3.0, 20),
           "colina": np.linspace(1.0, 12.0, 20)}

    res, mejor = {}, None
    print("\n%-11s %8s %10s %10s %9s" % ("familia", "param", "w0", "wa", "dist"),
          flush=True)
    for fam, pars in FAM.items():
        res[fam] = []
        for par in pars:
            r = dispara(fam, float(par))
            if r is None:
                continue
            w0, wa, V0, omf = r
            d = np.hypot((w0 - W0_OBJ) / TOL_W0, (wa - WA_OBJ) / TOL_WA)
            fila = dict(param=float(par), w0=w0, wa=wa, V0=V0, Omega_phi=omf,
                        cumple=bool(abs(w0 - W0_OBJ) < TOL_W0
                                    and abs(wa - WA_OBJ) < TOL_WA),
                        distancia=float(d))
            res[fam].append(fila)
            if mejor is None or d < mejor[0]:
                mejor = (d, fam, fila)
            print("%-11s %8.3f %+10.5f %+10.5f %9.2f%s"
                  % (fam, par, w0, wa, d, "  <== CUMPLE" if fila["cumple"] else ""),
                  flush=True)

    cumplen = [(f, x) for f, v in res.items() for x in v if x["cumple"]]
    print("\nobjetivo: w0 = %+.6f  wa = %+.6f" % (W0_OBJ, WA_OBJ), flush=True)
    print("familias que cumplen las DOS a la vez: %d" % len(cumplen), flush=True)
    if mejor:
        print("mas cercano: %s param=%.3f  w0=%+.5f  wa=%+.5f"
              % (mejor[1], mejor[2]["param"], mejor[2]["w0"], mejor[2]["wa"]),
              flush=True)

    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(dict(
        corrida="OP-24 — barrido de potenciales buscando wa = -0.669975",
        objetivo=dict(w0=W0_OBJ, wa=WA_OBJ, tol_w0=TOL_W0, tol_wa=TOL_WA,
                      criterio="escrito antes de correr; hay que cumplir LAS DOS"),
        control=dict(criterio="exp con lambda->0 da w=-1, wa=0",
                     w0=c[0], wa=c[1], pasa=True),
        Omega_DE_objetivo=OM_DE,
        familias=res,
        n_cumplen=len(cumplen),
        mas_cercano=None if mejor is None else dict(familia=mejor[1], **mejor[2]),
        alcance="mide si la energia oscura de SSEE puede describirse como un "
                "campo con potencial en estas familias. NO toca el acuerdo con "
                "DESI, que es del par algebraico y esta medido aparte"),
        indent=1))
    print("escrito -> %s" % SALIDA.relative_to(REPO), flush=True)


if __name__ == "__main__":
    main()
