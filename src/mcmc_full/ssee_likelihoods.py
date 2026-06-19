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

import numpy as np
from cobaya.likelihood import Likelihood

C_KMS = 299792.458  # velocidad de la luz [km/s]


class DESIDR2BAO(Likelihood):
    """DESI DR2 BAO — 13 puntos. tipo: 0=DM/rd, 1=DH/rd, 2=DV/rd."""

    def initialize(self):
        # z, valor observado, sigma, tipo
        self.data = np.array([
            [0.295,  7.93, 0.15, 2],   # BGS   DV/rd
            [0.510, 13.62, 0.25, 0],   # LRG1  DM/rd
            [0.510, 20.98, 0.61, 1],   # LRG1  DH/rd
            [0.706, 16.85, 0.32, 0],   # LRG2  DM/rd
            [0.706, 20.08, 0.60, 1],   # LRG2  DH/rd
            [0.930, 21.71, 0.28, 0],   # LRG3  DM/rd
            [0.930, 17.88, 0.35, 1],   # LRG3  DH/rd
            [1.317, 27.79, 0.69, 0],   # ELG   DM/rd
            [1.317, 13.82, 0.42, 1],   # ELG   DH/rd
            [1.491, 30.21, 0.79, 0],   # QSO   DM/rd
            [1.491, 13.23, 0.55, 1],   # QSO   DH/rd
            [2.330, 39.71, 0.94, 0],   # Lya   DM/rd
            [2.330,  8.52, 0.17, 1],   # Lya   DH/rd
        ])
        self.zuniq = np.unique(self.data[:, 0])
        self.zidx = {z: i for i, z in enumerate(self.zuniq)}
        self.log.info(f"DESI DR2 BAO cargado: {len(self.data)} puntos, "
                      f"{len(self.zuniq)} redshifts únicos.")

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
        chi2 = 0.0
        for z, obs, sig, tipo in self.data:
            i = self.zidx[z]
            DM = DM_arr[i]
            DH = C_KMS / H_arr[i]
            if tipo == 0:
                pred = DM / rd
            elif tipo == 1:
                pred = DH / rd
            else:  # DV/rd = [DM^2 · z · DH]^(1/3) / rd
                pred = (DM * DM * z * DH) ** (1.0 / 3.0) / rd
            chi2 += ((pred - obs) / sig) ** 2
        return -0.5 * chi2


class FSigma8(Likelihood):
    """fσ8 — 6 encuestas (6dFGRS, SDSS MGS, BOSS DR12 ×2, eBOSS, WiggleZ).

    Usa el fσ8(z) nativo de CAMB (crecimiento lineal consistente con los
    parámetros del paso), no un σ8_eff fijo. Covarianza diagonal.
    """

    def initialize(self):
        self.z = np.array([0.067, 0.150, 0.320, 0.380, 0.510, 0.570])
        self.obs = np.array([0.423, 0.490, 0.427, 0.477, 0.458, 0.426])
        self.err = np.array([0.055, 0.070, 0.056, 0.051, 0.038, 0.048])
        self.log.info(f"fσ8 cargado: {len(self.z)} encuestas (z ≤ 0.57).")

    def get_requirements(self):
        return {"fsigma8": {"z": self.z}}

    def logp(self, **params_values):
        pred = self.provider.get_fsigma8(self.z)
        chi2 = np.sum(((pred - self.obs) / self.err) ** 2)
        return -0.5 * chi2
