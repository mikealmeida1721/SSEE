"""
DESI DR2 BAO — cargador de la FUENTE ÚNICA DE VERDAD (data/raw/desi_dr2_bao.csv).

TODOS los scripts que usen los puntos BAO de DESI DR2 deben importar de aquí,
NO hardcodear los valores. Así se evita el drift csv↔código que introdujo datos
DR1 mal etiquetados como DR2 (ver cabecera del csv; corregido 2026-07-01).

Uso:
    import sys; sys.path.insert(0, "<repo>/src")
    from desi_dr2_data import load_desi_dr2, desi_covariance
    d = load_desi_dr2()          # dict con arrays z, value, sigma, quantity, tracer, corr
    C = desi_covariance(d)       # covarianza 13x13 bloque-diagonal (usa corr_MH)

Formato 'type' para likelihoods legacy: 0=DM/rd, 1=DH/rd, 2=DV/rd.
"""
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
# repo root = padre de src/
_CSV = os.path.normpath(os.path.join(_HERE, "..", "data", "raw", "desi_dr2_bao.csv"))

_QMAP = {"DM_over_rd": 0, "DH_over_rd": 1, "DV_over_rd": 2}


def load_desi_dr2(csv_path=None):
    """Lee el csv canónico y devuelve arrays alineados. Falla ruidoso si algo no cuadra."""
    path = csv_path or _CSV
    z, val, sig, corr, quant, trac, typ = [], [], [], [], [], [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("z_eff"):
                continue
            parts = line.split(",")
            if len(parts) < 6:
                raise ValueError(f"Fila DESI mal formada (esperaba >=6 columnas): {line!r}")
            z_eff, tracer, quantity, value, sigma, corr_MH = parts[:6]
            if quantity not in _QMAP:
                raise ValueError(f"Cantidad DESI desconocida: {quantity!r} en {line!r}")
            z.append(float(z_eff)); trac.append(tracer); quant.append(quantity)
            val.append(float(value)); sig.append(float(sigma))
            corr.append(float(corr_MH) if corr_MH != "" else np.nan)
            typ.append(_QMAP[quantity])
    d = {
        "z": np.array(z), "value": np.array(val), "sigma": np.array(sig),
        "corr": np.array(corr), "quantity": quant, "tracer": trac,
        "type": np.array(typ), "release": "DR2", "arxiv": "2503.14738",
    }
    if len(d["value"]) != 13:
        raise ValueError(f"DESI DR2 debe tener 13 puntos, csv tiene {len(d['value'])}")
    return d


def desi_covariance(d):
    """Covarianza 13x13 bloque-diagonal: bloques 2x2 (DM,DH) por tracer usando corr_MH.

    BGS (DV, isotrópico) queda diagonal. Cada par (DM,DH) del mismo tracer usa
    su correlación oficial r_MH de la Tabla 4 de DESI DR2 (2503.14738).
    """
    n = len(d["value"])
    C = np.diag(d["sigma"] ** 2).astype(float)
    for i in range(n - 1):
        same_tracer = d["tracer"][i] == d["tracer"][i + 1]
        is_dm_dh = d["type"][i] == 0 and d["type"][i + 1] == 1
        if same_tracer and is_dm_dh and not np.isnan(d["corr"][i]):
            cov = d["corr"][i] * d["sigma"][i] * d["sigma"][i + 1]
            C[i, i + 1] = C[i + 1, i] = cov
    return C


if __name__ == "__main__":
    d = load_desi_dr2()
    print(f"DESI {d['release']} ({d['arxiv']}): {len(d['value'])} puntos")
    for i in range(len(d["value"])):
        c = d["corr"][i]
        print(f"  z={d['z'][i]:.3f}  {d['tracer'][i]:10s} {d['quantity'][i]:11s} "
              f"= {d['value'][i]:7.3f} ± {d['sigma'][i]:.3f}  "
              f"corr={'—' if np.isnan(c) else f'{c:+.3f}'}")
    C = desi_covariance(d)
    npairs = int((np.abs(C - np.diag(np.diag(C))) > 0).sum() / 2)
    print(f"Covarianza 13x13 construida; pares DM-DH correlacionados: {npairs} (esperado 6)")
