"""
SSEE-V3.6 — Likelihoods custom para el MCMC FULL con Cobaya.

Dos likelihoods que complementan a Planck (plik_lite_native + lensing.native):

  - DESIDR2BAO : 13 puntos BAO de DESI DR2 (Abdul-Karim et al. 2025,
                 arXiv:2503.14738), idénticos a los del MCMC Fase 4.
  - FSigma8    : 6 encuestas fσ8 (z ≤ 0.57), idénticas a Paper 6 / Fase 4,
                 pero usando el fsigma8 que CAMB calcula nativamente en cada
                 paso (no un σ8_eff fijo como hacía Fase 4).

Ambas usan covarianza diagonal (errores independientes), igual que Fase 4.
"""

import os
import sys
import numpy as np
from cobaya.likelihood import Likelihood

# Fuente única de verdad para DESI DR2 (evita hardcodear → drift DR1/DR2)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from desi_dr2_data import load_desi_dr2, desi_covariance  # noqa: E402

C_KMS = 299792.458  # velocidad de la luz [km/s]


class DESIDR2BAO(Likelihood):
    """DESI DR2 BAO — 13 puntos (2503.14738 Tabla 4), cargados de la FUENTE ÚNICA
    data/raw/desi_dr2_bao.csv vía desi_dr2_data. Covarianza bloque-diagonal 2x2
    (DM,DH) por tracer con las correlaciones oficiales r_MH. type: 0=DM/rd,
    1=DH/rd, 2=DV/rd. NO hardcodear valores aquí (ver historia del error en el csv)."""

    def initialize(self):
        self.d = load_desi_dr2()
        self.zuniq = np.unique(self.d["z"])
        self.zidx = {z: i for i, z in enumerate(self.zuniq)}
        self.Cinv = np.linalg.inv(desi_covariance(self.d))
        self.log.info(f"DESI {self.d['release']} BAO ({self.d['arxiv']}): "
                      f"{len(self.d['value'])} puntos, {len(self.zuniq)} z únicos, "
                      f"covarianza bloque-diagonal (r_MH oficiales).")

    def get_requirements(self):
        return {
            "comoving_radial_distance": {"z": self.zuniq},
            "Hubble": {"z": self.zuniq, "units": "km/s/Mpc"},
            "rdrag": None,
        }

    def logp(self, **params_values):
        DM_arr = self.provider.get_comoving_radial_distance(self.zuniq)  # plano → D_M
        H_arr = self.provider.get_Hubble(self.zuniq, units="km/s/Mpc")
        rd = self.provider.get_param("rdrag")
        pred = np.empty(len(self.d["value"]))
        for k in range(len(pred)):
            z = self.d["z"][k]
            i = self.zidx[z]
            DM = DM_arr[i]
            DH = C_KMS / H_arr[i]
            tipo = self.d["type"][k]
            if tipo == 0:                       # DM/rd
                pred[k] = DM / rd
            elif tipo == 1:                     # DH/rd
                pred[k] = DH / rd
            else:                               # DV/rd = [DM^2 · z · DH]^(1/3) / rd
                pred[k] = (DM * DM * z * DH) ** (1.0 / 3.0) / rd
        r = pred - self.d["value"]
        return -0.5 * float(r @ self.Cinv @ r)


class FSigma8(Likelihood):
    """fσ8 — set canónico Paper 5/6 (6dFGRS, SDSS MGS, BOSS DR12, eBOSS DR16).

    Usa el fσ8(z) nativo de CAMB (crecimiento lineal consistente con los
    parámetros del paso), no un σ8_eff fijo. Covarianza diagonal.
    Idéntico a FSIG8_* de ssee_mcmc_fase4.py / verification.py (consistencia
    con toda la suite; el set anterior z=0.320/0.570 era no canónico).
    """

    def initialize(self):
        self.z = np.array([0.067, 0.150, 0.380, 0.510, 0.610, 1.480])
        self.obs = np.array([0.423, 0.490, 0.497, 0.458, 0.436, 0.462])
        self.err = np.array([0.055, 0.145, 0.045, 0.038, 0.034, 0.045])
        self.log.info(f"fσ8 cargado: {len(self.z)} encuestas canónicas (Paper 5/6, z ≤ 1.48).")

    def get_requirements(self):
        return {"fsigma8": {"z": self.z}}

    def logp(self, **params_values):
        pred = self.provider.get_fsigma8(self.z)
        chi2 = np.sum(((pred - self.obs) / self.err) ** 2)
        return -0.5 * chi2
