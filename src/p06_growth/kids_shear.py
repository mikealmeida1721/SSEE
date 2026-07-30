#!/usr/bin/env python3
"""
Pipeline de cizalla cosmica KiDS-1000 xi_pm — evaluador de likelihood independiente.

OBJETIVO: evaluar chi^2 de un modelo cosmologico arbitrario (incluido SSEE)
contra los 225 puntos de xi_pm de KiDS-1000, SIN pasar por un posterior LCDM.

CONTROL NEGATIVO OBLIGATORIO: reproducir like = -130.157350 en el punto de
maxima posterior de la cadena oficial (chain/maxpost_multinest_start_C.txt).
Hasta que eso pase, ningun resultado de este script es utilizable.

Referencias:
  Asgari et al. 2021, A&A 645, A104 (arXiv:2007.15633) — medidas xi_pm KiDS-1000
  Heymans et al. 2021, A&A 646, A140 (arXiv:2007.15632) — 3x2pt / cosmologia
  Hildebrandt et al. 2021, A&A 647, A124 (arXiv:2007.15635) — n(z) SOM gold
  Mead et al. 2015, MNRAS 454, 1958 — HMcode
  Bridle & King 2007, NJP 9, 444 — alineamiento intrinseco NLA
"""
import numpy as np
from astropy.io import fits
from scipy.integrate import simpson
from scipy.interpolate import InterpolatedUnivariateSpline as IUS
from scipy.special import jv
import camb

DATA = ('/mnt/datos/SSEE_data/kids1000/KiDS1000_cosmis_shear_data_release/'
        'data_fits/xipm_KIDS1000_BlindC_with_m_bias_V1.0.0A_ugriZYJHKs_photoz_'
        'SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid.fits')

NZBINS = 5
# cortes de escala del analisis fiducial (K1K_CorrelationFunctions.data)
KEEP_XIP = (0.5, 300.0)
KEEP_XIM = (4.0, 300.0)
# NLA: C1 * rho_crit en unidades de h^2 Msol/Mpc^3 -> adimensional (Bridle & King)
C1_RHOCRIT = 0.0134


def load_data():
    """Vector de datos, covarianza y n(z). Devuelve todo crudo, sin cortes."""
    f = fits.open(DATA)
    C = np.array(f['COVMAT'].data, dtype=float)
    xp, xm, nz = f['xiP'].data, f['xiM'].data, f['NZ_SOURCE'].data

    def vec(t):
        # orden del fits: tal cual viene; se respeta para casar con COVMAT
        return (np.array(t['BIN1'], int), np.array(t['BIN2'], int),
                np.array(t['ANG'], float), np.array(t['VALUE'], float))

    b1p, b2p, angp, valp = vec(xp)
    b1m, b2m, angm, valm = vec(xm)

    z_low = np.array(nz['Z_LOW'], float)
    z_high = np.array(nz['Z_HIGH'], float)
    z_mid = np.array(nz['Z_MID'], float)
    dz = z_high - z_low
    # OJO: las columnas BINi suman 1 (probabilidad por bin), NO son densidad.
    # Se convierten a densidad dividiendo por dz.
    nofz = np.array([np.array(nz[f'BIN{i+1}'], float) / dz for i in range(NZBINS)])

    d = np.concatenate([valp, valm])
    return dict(C=C, d=d,
                b1=np.concatenate([b1p, b1m]), b2=np.concatenate([b2p, b2m]),
                ang=np.concatenate([angp, angm]),
                kind=np.array(['p'] * len(valp) + ['m'] * len(valm)),
                z=z_mid, nofz=nofz, dz=dz)


def scale_mask(D):
    """True donde el punto entra en el ajuste fiducial (225 de 270)."""
    m = np.zeros(len(D['d']), bool)
    isp = D['kind'] == 'p'
    m[isp] = (D['ang'][isp] >= KEEP_XIP[0]) & (D['ang'][isp] <= KEEP_XIP[1])
    m[~isp] = (D['ang'][~isp] >= KEEP_XIM[0]) & (D['ang'][~isp] <= KEEP_XIM[1])
    return m


def shift_nz(z, nofz, deltaz):
    """n_biased(z) = n(z - delta_z), renormalizada. Convencion CosmoSIS."""
    out = np.zeros_like(nofz)
    for i in range(NZBINS):
        s = IUS(z, nofz[i], k=1, ext=1)
        n = s(z - deltaz[i])
        n[n < 0] = 0.0
        out[i] = n / simpson(n, x=z)
    return out


def run_camb(omch2, ombh2, h0, ns, As, mnu=0.06, w=-1.0, wa=0.0,
             halo_A=2.6, zmax=6.0, kmax=20.0, nz_pk=100):
    """P(k,z) no lineal con HMcode-2015 en la variante de 1 parametro de KiDS:
       c_min = halo_A ; eta_0 = 0.98 - 0.12*c_min  (Mead+2015 ec.30, valores KiDS)."""
    p = camb.CAMBparams()
    p.set_cosmology(H0=h0 * 100.0, ombh2=ombh2, omch2=omch2, mnu=mnu, omk=0.0)
    p.set_dark_energy(w=w, wa=wa, dark_energy_model='ppf')
    p.InitPower.set_params(As=As, ns=ns)
    zs = np.linspace(0.0, zmax, nz_pk)[::-1]
    p.set_matter_power(redshifts=zs, kmax=kmax, nonlinear=True)
    p.NonLinearModel = camb.nonlinear.Halofit()
    p.NonLinearModel.set_params(halofit_version='mead2015',
                                HMCode_A_baryon=halo_A,
                                HMCode_eta_baryon=0.98 - 0.12 * halo_A)
    r = camb.get_results(p)
    kh, z_pk, pk = r.get_matter_power_spectrum(minkh=1e-4, maxkh=kmax, npoints=400)
    # crecimiento D(z)/D(0) desde el P(k) LINEAL a k grande (escala segura)
    lin = camb.get_matter_power_interpolator(p, nonlinear=False, hubble_units=True,
                                             k_hunit=True, kmax=1.0, zmax=zmax)
    zg = np.linspace(0.0, zmax, 200)
    growth = np.sqrt(lin.P(zg, 0.01) / lin.P(0.0, 0.01))
    return r, p, kh, z_pk, pk, (zg, np.asarray(growth).ravel())


def cl_shear(D, res, pars, kh, z_pk, pk, growth, A_IA, deltaz, nell=60,
             ell_min=1.0, ell_max=3e4, nchi=400):
    """C_ell^ij de cizalla (GG + GI + II) por aproximacion de Limber extendida."""
    h = pars.h
    nofz = shift_nz(D['z'], D['nofz'], deltaz)
    z = D['z']

    # geometria comovil (Mpc, no Mpc/h)
    chi_of_z = IUS(z, res.comoving_radial_distance(z), k=3)
    chi = chi_of_z(z)
    chi_max = chi[-1]

    Om = (pars.omch2 + pars.ombh2 + pars.omnuh2) / h**2
    H0_c = 100.0 * h / 299792.458          # 1/Mpc

    # eficiencia de lente q_i(chi) = 1.5 Om H0^2 (chi/a) int_chi^inf n(chi') (chi'-chi)/chi'
    q = np.zeros((NZBINS, len(chi)))
    for i in range(NZBINS):
        n_chi = nofz[i] / np.gradient(chi, z)     # n(chi) = n(z) dz/dchi
        for j, c in enumerate(chi):
            sel = chi >= c
            if sel.sum() < 2:
                continue
            q[i, j] = simpson(n_chi[sel] * (chi[sel] - c) / chi[sel], x=chi[sel])
        q[i] *= 1.5 * Om * H0_c**2 * chi * (1.0 + z)

    # P(k,z) interpolado en (k [1/Mpc], z)
    from scipy.interpolate import RectBivariateSpline
    o = np.argsort(z_pk)                 # CAMB puede devolver z en cualquier orden
    zs_asc = np.asarray(z_pk)[o]
    pk_asc = np.asarray(pk)[o, :]
    spl = RectBivariateSpline(zs_asc, np.log(kh), np.log(pk_asc))

    def P(k_phys, zz):
        """P(k) en Mpc^3 con k en 1/Mpc."""
        kk = np.clip(k_phys / h, kh[0] * 1.001, kh[-1] * 0.999)
        zz = np.clip(zz, zs_asc[0], zs_asc[-1])
        return np.exp(spl.ev(zz, np.log(kk))) / h**3

    # factor de alineamiento intrinseco NLA (Bridle & King 2007):
    #   F(z) = -A_IA * C1*rho_crit * Om / D(z),  D normalizado a D(0)=1
    zg, Dg = growth
    Dz_spl = IUS(zg, Dg, k=3)
    F_IA = -A_IA * C1_RHOCRIT * Om / np.maximum(Dz_spl(z), 1e-6)

    ells = np.logspace(np.log10(ell_min), np.log10(ell_max), nell)
    npair = NZBINS * (NZBINS + 1) // 2
    Cl = np.zeros((npair, nell))
    idx = {}
    p_ = 0
    for i in range(NZBINS):
        for j in range(i, NZBINS):
            idx[(i + 1, j + 1)] = p_
            p_ += 1

    n_chi_all = np.array([nofz[i] / np.gradient(chi, z) for i in range(NZBINS)])

    good = chi > 1e-3
    for a, ell in enumerate(ells):
        k_phys = (ell + 0.5) / np.maximum(chi, 1e-3)
        Pk = P(k_phys, z)          # vectorizado sobre chi
        for i in range(NZBINS):
            for j in range(i, NZBINS):
                # GG + GI + IG + II con el kernel NLA
                W_i = q[i] + F_IA * n_chi_all[i]
                W_j = q[j] + F_IA * n_chi_all[j]
                integ = W_i * W_j / np.maximum(chi, 1e-3)**2 * Pk
                Cl[idx[(i + 1, j + 1)], a] = simpson(integ[good], x=chi[good])
    return ells, Cl, idx


def xi_pm(ells, Cl, theta_arcmin, kind):
    """xi_pm(theta) = 1/(2pi) int dl l C_l J_{0,4}(l theta)."""
    th = theta_arcmin * (np.pi / 180.0 / 60.0)
    order = 0 if kind == 'p' else 4
    lgrid = np.logspace(np.log10(ells[0]), np.log10(ells[-1]), 4000)
    spl = IUS(np.log(ells), np.log(np.maximum(Cl, 1e-40)), k=3)
    Cl_f = np.exp(spl(np.log(lgrid)))
    return simpson(lgrid * Cl_f * jv(order, lgrid * th), x=lgrid) / (2 * np.pi)


def theory_vector(D, ells, Cl, idx, delta_c=0.0):
    t = np.zeros(len(D['d']))
    for n in range(len(D['d'])):
        p_ = idx[(min(D['b1'][n], D['b2'][n]), max(D['b1'][n], D['b2'][n]))]
        t[n] = xi_pm(ells, Cl[p_], D['ang'][n], D['kind'][n])
        if D['kind'][n] == 'p':
            t[n] += delta_c**2          # termino c aditivo (xip_add_c_term = 1)
    return t


def chi2(D, theory, mask):
    dv = (theory - D['d'])[mask]
    Ci = np.linalg.inv(D['C'][np.ix_(mask, mask)])
    return float(dv @ Ci @ dv)


if __name__ == '__main__':
    D = load_data()
    m = scale_mask(D)
    print(f'puntos totales {len(D["d"])}  usados {m.sum()} '
          f'(xi+ {(m & (D["kind"]=="p")).sum()}, xi- {(m & (D["kind"]=="m")).sum()})')
