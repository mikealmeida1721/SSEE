#!/usr/bin/env python3
"""
Modelo de multipolos RSD para los datos CRUDOS de BOSS DR12.

Por que se necesita: el fsigma8 publicado ya lleva dentro la cosmologia fiducial
LCDM (Omega_m=0.31) a traves del efecto Alcock-Paczynski. Aqui el AP se aplica al
MODELO -- se deforma la prediccion de SSEE hacia la rejilla fiducial en la que
Beutler midio -- de modo que el dato queda intacto.

Piezas, en el orden en que actuan:

  1. P_lin(k) de CAMB con el fondo que se le pase (SSEE algebraico o LCDM).
  2. Kaiser (b1 + f mu^2)^2 con amortiguamiento de dispersion de velocidades
     (Lorentziano al cuadrado, el mismo funcional que usa Beutler para el FoG).
  3. Alcock-Paczynski via alpha_par, alpha_perp contra el fiducial DECLARADO en
     los headers del fichero (Omega_m=0.31, w=-1, H0=100).
  4. Proyeccion a multipolos por cuadratura Gauss-Legendre en mu.
  5. Convolucion con la ventana del survey en espacio de configuracion
     (Wilson et al. 2017 eq. 20-22): P_ell -> xi_ell -> mezcla RR -> P_ell.

Unidades: TODO en h/Mpc y (Mpc/h)^3, como los ficheros de BOSS.
"""
# ORIGEN-VALOR: 2997.92458 — c/(100 km/s/Mpc) = 299792.458/100 Mpc/h
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import InterpolatedUnivariateSpline as Spline
from mcfit import P2xi, xi2P

# ---------- fiducial DECLARADO en los headers de Beutler ----------
FID_OM = 0.31
FID_W = -1.0
C_KMS = 299792.458

# rejilla interna: ancha, porque la transformada de Hankel necesita colas
K_MODEL = np.logspace(-5, 2, 2048)          # h/Mpc

_GL_X, _GL_W = leggauss(120)                # cuadratura en mu, [-1,1]


# ==================================================================
#  Fondo
# ==================================================================
def _E(z, Om, w0=-1.0, wa=0.0):
    """H(z)/H0 con materia + DE con w(a)=w0+wa(1-a). Sin radiacion (z<1)."""
    a = 1.0 / (1.0 + z)
    ode = (1.0 - Om) * a ** (-3 * (1 + w0 + wa)) * np.exp(-3 * wa * (1 - a))
    return np.sqrt(Om * (1 + z) ** 3 + ode)


def _DM(z, Om, w0=-1.0, wa=0.0, n=2000):
    """distancia comovil transversa en Mpc/h (c/H0 = c/(100 h km/s/Mpc) = 299792.458/100 = 2997.92458 Mpc/h)."""
    zz = np.linspace(0.0, z, n)
    return 2997.92458 * np.trapezoid(1.0 / _E(zz, Om, w0, wa), zz)


def alphas(z, Om, w0=-1.0, wa=0.0):
    """(alpha_par, alpha_perp) del fondo dado CONTRA el fiducial de Beutler.

    alpha_par  = [H_fid / H] ,  alpha_perp = [D_M / D_M_fid]
    (r_d NO entra: estos datos son full-shape pre-reconstruccion; el bariónico
    se absorbe en la forma, no se usa la escala BAO como regla.)
    """
    a_par = _E(z, FID_OM, FID_W, 0.0) / _E(z, Om, w0, wa)
    a_per = _DM(z, Om, w0, wa) / _DM(z, FID_OM, FID_W, 0.0)
    return a_par, a_per


# ==================================================================
#  Multipolos con AP
# ==================================================================
def multipoles(k_out, plin, f, b1, sig_v, a_par=1.0, a_per=1.0):
    """
    P_0, P_2, P_4 evaluados en k_out (rejilla del OBSERVADOR, o sea fiducial).

    plin : callable P_lin(k) en (Mpc/h)^3, k en h/Mpc, al z_eff del bin.
    sig_v: dispersion de velocidades en Mpc/h.
    """
    mu_o = _GL_X[None, :]
    k_o = np.asarray(k_out)[:, None]

    # rejilla observada -> rejilla verdadera (AP)
    F = a_par / a_per
    fac = np.sqrt(1.0 + mu_o ** 2 * (1.0 / F ** 2 - 1.0))
    k_t = k_o / a_per * fac
    mu_t = mu_o / F / fac

    kaiser = (b1 + f * mu_t ** 2) ** 2
    fog = 1.0 / (1.0 + 0.5 * (k_t * mu_t * sig_v) ** 2) ** 2
    P = kaiser * fog * plin(k_t) / (a_per ** 2 * a_par)

    out = []
    for ell in (0, 2, 4):
        L = np.polynomial.legendre.legval(_GL_X, [0] * ell + [1])
        out.append((2 * ell + 1) / 2.0 * (P * (_GL_W * L)[None, :]).sum(axis=1))
    return out


# ==================================================================
#  Ventana del survey
# ==================================================================
def window_normalised(s, RR, s_fit=(2.0, 20.0)):
    """
    RR crudos -> W_ell(s) con la normalizacion W_0(s->0) = 1.

    Dos cosas que hay que hacer bien, las dos aprendidas a base de que fallara:

    1. Los bines de s son LOGARITMICOS (ds/s = 2.30e-3 constante), asi que el
       conteo de una capa esferica de campo uniforme va como s^2 ds, NO como
       s^2. Dividir por s^2 a secas deja una pendiente espuria de un factor s
       (W_0 subiria a ~13 en s=200 en vez de decaer hacia 0).

    2. La normalizacion es el valor en s -> 0, y ese limite no se alcanza
       promediando a lo bruto: por debajo de s ~ 5 Mpc/h hay menos de 1500 pares
       y el ruido de conteo pasa del 20%, mientras que por encima la ventana YA
       esta cayendo (0.997 en s=10, 0.932 en s=50). El rango tiene que quedar en
       medio, y la eleccion NO es inocua: moverlo cambia la escala un 1.3-3.7%,
       que se propaga entera a la amplitud y por tanto a f*sigma8.

       CRITERIO (medido, no elegido a ojo): el rango bueno es el que deja W_0
       valiendo 1 cerca del origen. Medido sobre los 6 conjuntos:

           [2,20] -> W_0(8) = 1.0005 +- 0.0057   <-- elegido
           [3,30] -> W_0(8) = 1.0069
           [5,50] -> W_0(8) = 1.0192
           [2,60] -> W_0(8) = 1.0123

       Se probo tambien extrapolar con un ajuste pesado g = A + B s^2: sale PEOR
       (dispersion 2-3.8% y W_0 > 1 donde ya deberia bajar), porque el peso de
       Poisson carga en los pares lejanos y la caida no es cuadratica. Retirado.

    LIMITACION CONOCIDA: queda una sistematica de ~2% en la escala de la ventana,
    heredada de esta eleccion. Es menor que el sesgo de Kaiser (9-15%) pero NO es
    cero, y hay que tratarla en la fase de rigor (marginalizar la normalizacion).
    """
    ds = np.gradient(s)
    vol = s ** 2 * ds
    m = (s > s_fit[0]) & (s < s_fit[1]) & (vol > 0)
    A = np.median(RR[m, 0] / vol[m])
    return RR / (A * vol[:, None])


# mezcla de multipolos de Wilson et al. 2017, eq. 20-22 (truncada en ell=4)
def _mix(x0, x2, x4, W):
    """W = (W0, W2, W4, W6, W8) evaluados en la misma s que x_ell."""
    W0, W2, W4, W6, W8 = W
    y0 = x0 * W0 + x2 * (W2 / 5.0) + x4 * (W4 / 9.0)
    y2 = (x0 * W2
          + x2 * (W0 + 2.0 / 7.0 * W2 + 2.0 / 7.0 * W4)
          + x4 * (2.0 / 7.0 * W2 + 100.0 / 693.0 * W4 + 25.0 / 143.0 * W6))
    y4 = (x0 * W4
          + x2 * (18.0 / 35.0 * W2 + 20.0 / 77.0 * W4 + 45.0 / 143.0 * W6)
          + x4 * (W0 + 20.0 / 77.0 * W2 + 162.0 / 1001.0 * W4
                  + 20.0 / 143.0 * W6 + 490.0 / 2431.0 * W8))
    return y0, y2, y4


# Indices de la ley de potencias con que se rellena el modelo fuera de su
# rejilla, ANTES de la transformada. Ver la NOTA CRITICA de Window.
N_LO = 1.0        # k -> 0 : P ~ k^n_s, medido +0.56/+1.00/+0.96 en P0/P2/P4
N_HI = -2.0       # k -> oo: cola amortiguada; medido +0.33..-1.09, que es el
                  # regimen de los contraterminos y no una ley fisica


class Window:
    """Convolucion con la ventana, precalculando los splines de W_ell.

    NOTA CRITICA (2026-09-07) --- POR QUE EL RELLENO LLEVA INDICE FIJO.

    Hasta hoy las tres transformadas se llamaban con `extrap=True`, que en
    mcfit rellena el vector hasta la longitud de la convolucion con una ley de
    potencias cuyo INDICE SE MIDE DEL PROPIO VECTOR DE ENTRADA. Eso tiene dos
    consecuencias, ambas medidas:

      1. La ventana deja de ser un OPERADOR LINEAL. Con el modelo evaluado en
         el mejor punto de la cadena, la identidad
             m(a0,a2,sn) = m(0,0,0) + a0*T0 + a2*T2 + sn*Tn
         fallaba con error relativo 0.2 a 1.2, pese a que los multipolos
         CRUDOS que salen de velocileptors la cumplen a 1.8e-14. Toda la no
         linealidad la metia el relleno.
      2. Con contraterminos grandes la cola cambia de indice y el FFTLog
         devuelve valores de hasta 1e80 (medido con th=[1.2,0,0,0,0,-3000]).
         Eso son precipicios reales del modelo, no del optimizador.

    Con indice FIJO el relleno es `F[0]*(x/x0)^N_LO` y `F[-1]*(x/x1)^N_HI`,
    o sea LINEAL en F: la ventana entera vuelve a ser un operador lineal
    (verificado a 4e-12) y desaparecen los precipicios. Eso habilita
    marginalizar analiticamente los 18 parametros que entran lineales.

    Alternativas descartadas, con su medida en el mejor punto de la cadena
    (chi2 de referencia 268.13, y el modelo NO se re-ajusto):
        extrap=False (relleno con ceros)      chi2 = 690.86   lineal
        extrap='const' (valor del extremo)    chi2 = 170058   lineal
        extrap=True   (el que habia)          chi2 = 268.13   NO lineal

    LIMITACION CONOCIDA: el chi2 al MISMO punto depende de N_HI
    (-1 -> 247.5, -1.5 -> 298.2, -2 -> 366.6, -3 -> 483.1). El indice es por
    tanto una eleccion de modelado que hay que declarar, y que se controla
    re-ajustando con varios valores: lo que tiene que ser estable no es el
    chi2 a parametros congelados sino logA ajustado. Se aplica IDENTICA a los
    dos fondos, asi que en la diferencia entre modelos se cancela en gran parte.
    """

    def __init__(self, s, RR, integral_constraint=True,
                 n_lo=N_LO, n_hi=N_HI):
        self.W = window_normalised(s, RR)
        self.s_grid = s
        self.ic = integral_constraint
        self.n_lo = float(n_lo)
        self.n_hi = float(n_hi)
        self._spl = [Spline(s, self.W[:, i], k=1, ext=1) for i in range(5)]
        self._pad = {}

    def _relleno(self, tr):
        """Factores de relleno de esta transformada. Se cachean por (rejilla,
        multipolo): la clave lleva la rejilla, NO el id() del objeto, porque
        mcfit crea uno nuevo en cada llamada y los id() se reciclan."""
        x = tr.x
        cl = (len(x), float(x[0]), float(x[-1]), tr.N)
        if cl not in self._pad:
            X = tr._x_
            n = (tr.N - len(x)) // 2
            self._pad[cl] = ((X[:n] / x[0]) ** self.n_lo,
                             (X[len(x) + n:] / x[-1]) ** self.n_hi)
        return self._pad[cl]

    def _a_xi(self, k_model, P, ell):
        tr = P2xi(k_model, l=ell, lowring=False)
        fl, fh = self._relleno(tr)
        F = np.concatenate([P[0] * fl, P, P[-1] * fh])
        return tr(F, extrap=False)

    def __call__(self, k_model, P0, P2, P4):
        xs, x0 = self._a_xi(k_model, P0, 0)
        _, x2 = self._a_xi(k_model, P2, 2)
        _, x4 = self._a_xi(k_model, P4, 4)
        Wg = [sp(xs) for sp in self._spl]

        y0, y2, y4 = _mix(x0, x2, x4, Wg)

        if self.ic:
            # Restriccion integral: la densidad media se estima del PROPIO
            # survey, asi que la potencia del modo de la ventana se anula.
            # En espacio de configuracion (Wilson+2017 eq. 19) equivale a
            # restar a cada xi_ell la ventana escalada por la media pesada
            # del monopolo convolucionado.
            m = (xs >= self.s_grid[0]) & (xs <= self.s_grid[-1])
            num = np.trapezoid((y0 * Wg[0] * xs ** 2)[m], xs[m])
            den = np.trapezoid((Wg[0] * xs ** 2)[m], xs[m])
            ic = num / den
            y0 = y0 - ic * Wg[0]
            y2 = y2 - ic * Wg[1]
            y4 = y4 - ic * Wg[2]

        ks, q0 = xi2P(xs, l=0, lowring=False)(y0, extrap=False)
        _, q2 = xi2P(xs, l=2, lowring=False)(y2, extrap=False)
        _, q4 = xi2P(xs, l=4, lowring=False)(y4, extrap=False)
        return ks, q0.real, q2.real, q4.real
