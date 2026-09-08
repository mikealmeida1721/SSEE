"""El techo sigma8 = 0.8335: que es, y su log, que nunca existio.

POR QUE EXISTE (2026-09-08). Lo pidio Mike al preguntar cuales afirmaciones
tienen respaldo de verdad. Este era el ULTIMO valor de la seccion B del
Registro sin log committeado, y ademas estaba MAL ATRIBUIDO: la columna
Fuente citaba `ssee_paper5_IS_perturbations.py`, que no contiene el numero ni
lo calcula. Sale de una corrida de CLASS, `output/can_cold__pk.dat`.

QUE CONTESTA EL NUMERO. Cuanta estructura predice el fondo de SSEE si la
amplitud primordial A_s se HEREDA de Planck en vez de ajustarse. No es una
prediccion del modelo: A_s es uno de los dos libres del sector CMB (k=2), asi
que clavarlo al valor de Planck es importar una inferencia hecha dentro de
LCDM, y con ella la discrepancia Planck-cizalla. Por eso el "3.5 sigma" que
salia de aqui era un ARTEFACTO del condicionamiento, no una tension del
modelo. Con A_s libre contra el dato crudo de KiDS: S8 = 0.7555 +/- 0.0192,
0.11 sigma (Paper 6, R3). El techo se conserva como DIAGNOSTICO bajo una
condicion declarada, no como prediccion.

COMO. sigma8^2 = (1/2pi^2) Int k^2 P(k) W^2(kR) dk, R = 8 Mpc/h,
W(x) = 3(sin x - x cos x)/x^3 (Peebles 1980). CLASS entrega 119 puntos, que
no bastan para el filtro: se interpola en log-log a 4000.

CONTROL (R53) — NO CONCLUYENTE, y se dice por que. La idea era pasar el mismo
integrador por un P(k) de referencia LCDM: si cae donde Planck pone su sigma8
(0.8111 +/- 0.006), el metodo queda validado y no solo el numero de SSEE. Sale
0.8221, un 1.4% alto. La causa esta identificada y NO es el integrador: las
corridas LCDM guardadas no llevan neutrinos masivos (sin `Omega_ncdm`), y sin
ellos sobra poder a escala de 8 Mpc/h justo en ese orden. Ademas el .ini que
lleva ese nombre pide `output = tCl,pCl,lCl`, o sea ni siquiera genera mPk: el
fichero vino de otra corrida no identificada.

Asi que este log respalda que el 0.8335 SALE de `can_cold__pk.dat` con el
metodo declarado, y NO respalda todavia que el metodo este calibrado. Para
cerrarlo hace falta correr CLASS con la baseline de Planck (Sum m_nu = 0.06 eV)
y comprobar que da 0.811. Queda anotado como pendiente, no como verde.
"""
import json
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
PK_SSEE = REPO / "class_ssee" / "output" / "can_cold__pk.dat"
PK_LCDM = REPO / "class_ssee" / "output" / "ref_lcdm__pk.dat"
SALIDA = REPO / "results" / "logs" / "p5_techo_sigma8_As_fijo.json"
R_TOPHAT = 8.0            # Mpc/h
N_INTERP = 4000


def sigma8(archivo):
    k, P = np.loadtxt(archivo, unpack=True)
    lk = np.log(k)
    kk = np.exp(np.linspace(lk[0], lk[-1], N_INTERP))
    PP = np.exp(np.interp(np.log(kk), lk, np.log(P)))
    x = kk * R_TOPHAT
    W = 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3
    return float(np.sqrt(np.trapezoid(kk ** 2 * PP * W ** 2, kk) / (2 * np.pi ** 2)))


def main():
    s_ssee = sigma8(PK_SSEE)
    s_lcdm = sigma8(PK_LCDM)
    OM = 0.3088808787787524            # Omega_m,CMB = omega_m/h^2 (algebraico)
    S8 = s_ssee * (OM / 0.3) ** 0.5
    out = {
        "corrida": "techo sigma8/S8 con A_s FIJADO a Planck (diagnostico, NO prediccion)",
        "que_contesta": "cuanta estructura da el fondo de SSEE si A_s se hereda "
                        "de Planck en vez de ajustarse; A_s es libre en el modelo",
        "fuente_Pk": str(PK_SSEE.relative_to(REPO)),
        "metodo": "top-hat R=8 Mpc/h, Peebles 1980; interp log-log a %d puntos"
                  % N_INTERP,
        "sigma8_techo": s_ssee,
        "Omega_m_CMB": OM,
        "S8_techo": S8,
        "control_LCDM": {
            "estado": "NO CONCLUYENTE",
            "fuente_Pk": str(PK_LCDM.relative_to(REPO)),
            "sigma8": s_lcdm,
            "sigma8_Planck2018": 0.8111,
            "nota": "sale 1.4% alto. Causa identificada y ajena al integrador: "
                    "las corridas LCDM guardadas no llevan neutrinos masivos, y "
                    "sin ellos sobra poder a 8 Mpc/h en ese orden. Para cerrarlo "
                    "hay que correr CLASS con la baseline de Planck "
                    "(Sum m_nu = 0.06 eV) y comprobar que da 0.811.",
        },
        "superado_por": {
            "valor": "S8 = 0.7555 +/- 0.0192 (0.11 sigma vs KiDS-1000)",
            "log": "results/logs/growth_2026-07/R3_ssee_kids_S8.json",
            "razon": "A_s libre contra el dato crudo; el 3.5 sigma que salia "
                     "del techo era artefacto de fijar A_s",
        },
    }
    SALIDA.parent.mkdir(parents=True, exist_ok=True)
    SALIDA.write_text(json.dumps(out, indent=1))
    print("sigma8 techo (A_s fijo) = %.6f   [Registro 0.8335]" % s_ssee)
    print("S8 techo                = %.6f   [Registro 0.846]" % S8)
    print("CONTROL sigma8 LCDM ref = %.6f" % s_lcdm)
    print("escrito ->", SALIDA.relative_to(REPO))


if __name__ == "__main__":
    main()
