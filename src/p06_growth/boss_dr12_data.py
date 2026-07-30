#!/usr/bin/env python3
"""
Lector de los multipolos CRUDOS de BOSS DR12 (Beutler et al. 2017, 1607.03150).

Por que crudos: el `fsigma8_rsd.csv` que usabamos son valores YA derivados con una
cosmologia fiducial LCDM dentro (efecto Alcock-Paczynski cocido en el numero). Con
P0/P2/P4 + covarianza + ventana, el AP se aplica al MODELO, no al dato, y el test
del fondo de SSEE es limpio.

Fuente: https://data.sdss.org/sas/dr12/boss/papers/clustering/
        Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz  (3.2 MB, 2016-07-08)
Local:  /mnt/datos/SSEE_data/boss_dr12/public_material_RSD/

Fiducial DECLARADO en los headers: Omega_m=0.31, Omega_nu=0, w=-1, H0=100 (unid. h).
z_eff = 0.38 (z1), 0.51 (z2), 0.61 (z3).

Formatos:
  multipolos : k_min  k_max  P(k)  sigma          (~40 bins, dk=0.01)
  covarianza : i  j  k_i  k_j  C_ij               (P0+P2+P4 conjuntos, mocks PATCHY)
  ventana    : s  dummy  RR0 RR2 RR4 RR6 RR8      (espacio de configuracion)
"""
from pathlib import Path
import numpy as np

BOSS = Path('/mnt/datos/SSEE_data/boss_dr12/public_material_RSD')
ZBINS = {'z1': 0.38, 'z2': 0.51, 'z3': 0.61}
CAPS = ('NGC', 'SGC')
MULTIPOLES = {0: 'monopole', 2: 'quadrupole', 4: 'hexadecapole'}
# mocks usados en cada covarianza (va en el nombre del fichero)
NMOCKS = {'NGC': 2045, 'SGC': 2048}


def _numeric_rows(path):
    """filas puramente numericas — los headers son texto libre sin '#' fiable."""
    out = []
    for line in path.read_text(errors='ignore').splitlines():
        parts = line.split()
        if not parts:
            continue
        try:
            out.append([float(p) for p in parts])
        except ValueError:
            continue
    return out


def load_multipole(zbin, cap, ell):
    """devuelve (k_centro, P_ell, sigma) del fichero crudo."""
    f = BOSS / (f'Beutleretal_pk_{MULTIPOLES[ell]}_DR12_{cap}_{zbin}'
                f'_prerecon_120.dat')
    rows = [r for r in _numeric_rows(f) if len(r) == 4]
    a = np.array(rows)
    return 0.5 * (a[:, 0] + a[:, 1]), a[:, 2], a[:, 3]


def load_cov(zbin, cap):
    """
    Covarianza densa de (P0,P2,P4) concatenados, mas su k y su ell por elemento.

    Los indices del fichero viven en una rejilla de 60 bins POR MULTIPOLO
    (P0: 1-60, P2: 61-120, P4: 121-180) y solo aparecen los bins usados
    (aqui 1-14, 61-74, 121-129). Se compactan a posiciones consecutivas.
    """
    n = NMOCKS[cap]
    hits = list(BOSS.glob(f'Beutleretal_cov_patchy_{zbin}_{cap}_*_{n}_60.dat'))
    if not hits:
        raise FileNotFoundError(f'sin covarianza para {zbin} {cap}')
    rows = np.array([r for r in _numeric_rows(hits[0]) if len(r) == 5])
    raw_i = rows[:, 0].astype(int)
    raw_j = rows[:, 1].astype(int)
    used = np.unique(np.concatenate([raw_i, raw_j]))          # ya ordenado
    pos = {v: p for p, v in enumerate(used)}
    N = len(used)
    C = np.zeros((N, N))
    kk = np.zeros(N)
    ell = np.zeros(N, dtype=int)
    for a, b, ka, _kb, c in rows:
        p, q = pos[int(a)], pos[int(b)]
        C[p, q] = C[q, p] = c                                  # simetrizar
        kk[p] = ka
        ell[p] = 2 * ((int(a) - 1) // 60)                      # 0, 2 o 4
    return C, kk, ell


def load_window(zbin, cap):
    """ventana en espacio de configuracion: s, RR0, RR2, RR4, RR6, RR8."""
    f = BOSS / f'Beutleretal_window_{zbin}_{cap}.dat'
    a = np.array([r for r in _numeric_rows(f) if len(r) == 7])
    return a[:, 0], a[:, 2:]


if __name__ == '__main__':
    print('BOSS DR12 — inventario y verificacion de lectura')
    print('=' * 74)
    for zbin, zeff in ZBINS.items():
        for cap in CAPS:
            k0, P0, s0 = load_multipole(zbin, cap, 0)
            k2, P2, s2 = load_multipole(zbin, cap, 2)
            k4, P4, s4 = load_multipole(zbin, cap, 4)
            C, kc, ellc = load_cov(zbin, cap)
            s, RR = load_window(zbin, cap)
            print(f'\n  {zbin} {cap}  (z_eff = {zeff})')
            print(f'    P0: {len(k0):3d} bins  k = {k0[0]:.4f} .. {k0[-1]:.4f} '
                  f'h/Mpc   P0(k_min) = {P0[0]:.1f}')
            print(f'    P2: {len(k2):3d} bins   P4: {len(k4):3d} bins')
            print(f'    cov: {C.shape[0]}x{C.shape[1]}  '
                  f'(esperado {len(k0)+len(k2)+len(k4)} si cubre P0+P2+P4)')
            print(f'    ventana: {len(s)} filas, s = {s[0]:.3f} .. {s[-1]:.1f} '
                  f'Mpc/h, {RR.shape[1]} multipolos RR')
            # CONTROL DE LECTURA: sigma del fichero de datos vs sqrt(diag C),
            # emparejando por (ell, k) — no por posicion.
            d = np.sqrt(np.diag(C))
            for ell, kd, sig in ((0, k0, s0), (2, k2, s2), (4, k4, s4)):
                m = ellc == ell
                idx = [int(np.argmin(abs(kd - kv))) for kv in kc[m]]
                r = d[m] / sig[idx]
                print(f'    control P{ell}: {m.sum():2d} bins  '
                      f'sqrt(C_ii)/sigma_fichero = {r.mean():.4f} '
                      f'+- {r.std():.4f}   (1.0000 = lectura consistente)')
            ev = np.linalg.eigvalsh(C)
            print(f'    cov definida positiva: {ev.min() > 0}  '
                  f'(autovalor minimo {ev.min():.3e})')
