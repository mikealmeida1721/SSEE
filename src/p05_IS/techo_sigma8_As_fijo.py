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

CONTROL (R53) — QUE MIDE Y QUE NO. Lo senalo Mike: comparar el 0.8221 de la
referencia LCDM contra el 0.8111 que publica Planck NO es un control de nada
de SSEE, porque son dos cosas con ingredientes distintos. El control que si
dice algo es SSEE contra LCDM con EL MISMO integrador y el mismo tratamiento:

    SSEE  0.833368  /  LCDM  0.822068  =  1.0137   (SSEE da 1.4% mas)

y ese 1.4% es el contenido fisico del techo. La comparacion contra el numero
publicado de Planck queda como diagnostico del METODO, y ahi sale alta porque
ninguna de las dos corridas guardadas lleva neutrinos masivos; sin ellos sobra
poder a 8 Mpc/h justo en ese orden. Al faltar en las DOS, la razon de arriba
no se ve afectada, pero el valor absoluto de cada una si.

LO QUE ESTE LOG **NO** ESTABLECE, y hay que decirlo entero:

  · `can_cold__pk.dat` no tenia .ini. No hay fichero de configuracion que
    diga con que se corrio, asi que no puedo certificar sus ingredientes.
  · Estaba FUERA del repositorio: `class_ssee/output/` esta en .gitignore, y
    el fichero nunca se commiteo. Vivia solo en el disco de Mike. Se copia
    aqui, a `results/logs/p5_techo_evidencia/`, para que exista de verdad.
  · Su gemelo `can_part__pk.dat`, del mismo minuto, es la corrida de DOS
    SECTORES con la particula, retirada el 2026-08-01. El techo nacio como
    su termino de comparacion.
  · El fondo canonico de SSEE SI lleva neutrinos masivos (Sum m_nu = 0.06849
    eV, w_nu = 0.000735 dentro de w_m = 0.14267). Si esta corrida no los
    lleva, entonces 0.8335 no es el techo del modelo canonico sino el de una
    variante sin ellos.

PENDIENTE REAL para cerrarlo: un .ini VERSIONADO con el fondo canonico y sus
neutrinos, corrido con CLASS, mas la misma referencia LCDM con la baseline de
Planck. Hasta entonces el numero se cita como diagnostico de procedencia
incompleta, no como cantidad certificada.
"""
import json
import pathlib

import numpy as np

REPO = pathlib.Path(__file__).resolve().parents[2]
EVID = REPO / "results" / "logs" / "p5_techo_evidencia"
PK_SSEE = EVID / "can_cold__pk.dat"          # copiado de class_ssee/output/,
PK_LCDM = EVID / "lcdm_planck2018__pk.dat"   # que esta en .gitignore
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
        "control_SSEE_vs_LCDM_mismo_metodo": {
            "fuente_Pk": str(PK_LCDM.relative_to(REPO)),
            "sigma8_LCDM": s_lcdm,
            "razon_SSEE_sobre_LCDM": s_ssee / s_lcdm,
            "nota": "ESTE es el control que dice algo: el mismo integrador y el "
                    "mismo tratamiento en los dos. SSEE da 1.4% mas estructura.",
        },
        "diagnostico_del_metodo": {
            "estado": "NO CONCLUYENTE",
            "sigma8_Planck2018_publicado": 0.8111,
            "nota": "la referencia sale 1.4% alta contra el valor publicado, "
                    "porque ninguna de las dos corridas guardadas lleva "
                    "neutrinos masivos. Al faltar en las DOS, la razon de "
                    "arriba no se ve afectada; el valor absoluto de cada una si.",
        },
        "lo_que_NO_esta_establecido": [
            "can_cold__pk.dat no tiene .ini: no se puede certificar con que se corrio",
            "estaba fuera del repo (class_ssee/output/ esta en .gitignore); "
            "copiado a results/logs/p5_techo_evidencia/ para que exista",
            "su gemelo can_part__pk.dat, del mismo minuto, es la corrida de dos "
            "sectores con la particula, retirada el 2026-08-01",
            "el fondo canonico de SSEE lleva Sum m_nu = 0.06849 eV; si esta "
            "corrida no los lleva, 0.8335 no es el techo del modelo canonico",
        ],
        "pendiente_para_cerrarlo": ".ini versionado con el fondo canonico y sus "
                                   "neutrinos + la misma referencia LCDM con la "
                                   "baseline de Planck",
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
