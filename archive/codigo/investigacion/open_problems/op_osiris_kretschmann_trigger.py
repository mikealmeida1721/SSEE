"""
Filtro de Osiris — primer test físico del disparador de Kretschmann.

Pregunta: el acoplamiento K/(2 M_*^2) phi^2 cambia la masa efectiva de phi
solo donde el escalar de Kretschmann K supera el umbral K_th = m_phi^2 M_*^2.
Cerca de un PBH de masa-Ceres, ¿qué fracción del volumen del universo queda
"disparada"? Si esa fracción es ~0, el mecanismo NO puede convertir el grueso
de DM en DE (ni viceversa) y la unificacion via PBH falla.

Unidades naturales (hbar = c = 1), todo en eV salvo donde se indique.
Autor de referencia: Mike Almeida (SSEE). Cuenta honesta, sin fiteo.
"""
import numpy as np

# ---- constantes fundamentales ----
phi = (1 + np.sqrt(5)) / 2
M_Pl_red = 2.435e27      # masa de Planck reducida [eV]
M_Pl_np  = 1.2209e28     # masa de Planck NO reducida [eV]  (G = 1/M_Pl_np^2)
Lambda   = 8.81e-3       # Lambda_SSEE = M_* = escala UV de DE [eV]
H0       = 1.44e-33      # H0 ~ 67 km/s/Mpc en eV

# ---- conversiones ----
eV_per_m_inv = 1.973e-7  # 1/m  = 1.973e-7 eV   (hbar c = 1.973e-7 eV*m)
m_per_eV_inv = 1 / eV_per_m_inv * 1e0   # no usado directamente
kg_to_eV = 5.6096e35     # 1 kg = 5.6096e35 eV
M_sun_kg = 1.989e30
Ceres_kg = 9.1e20

def m_to_eVinv(x_m):      # longitud en metros -> eV^-1
    return x_m / 1.973e-7
def eVinv_to_m(x):        # eV^-1 -> metros
    return x * 1.973e-7

print("="*64)
print("FILTRO DE OSIRIS — disparador de Kretschmann (cuenta honesta)")
print("="*64)

# ---- 1. masa critica de PBH (re-derivacion del numero de Mike) ----
# Condicion: K(r_s) = Lambda^4  con  m_phi = M_* = Lambda
# => M_BH,crit = (3/4)^(1/4) * M_Pl^2 / Lambda
M_BH_crit_eV = (3/4)**0.25 * M_Pl_red**2 / Lambda
M_BH_crit_kg = M_BH_crit_eV / kg_to_eV
print(f"\n[1] Masa critica de PBH (3/4)^1/4 M_Pl^2/Lambda:")
print(f"    M_BH,crit = {M_BH_crit_kg:.3e} kg = {M_BH_crit_kg/M_sun_kg:.3e} M_sun")
print(f"    (masa de Ceres = {Ceres_kg:.2e} kg = {Ceres_kg/M_sun_kg:.2e} M_sun)")

# ---- 2. radio de Schwarzschild del PBH de masa Ceres ----
M_BH = Ceres_kg * kg_to_eV          # usamos Ceres como PBH representativo
r_s = 2 * M_BH / M_Pl_np**2          # eV^-1
print(f"\n[2] PBH de masa-Ceres:  r_s = {eVinv_to_m(r_s)*1e6:.2f} micrometros")

# ---- 3. umbral de disparo y radio de influencia ----
# m_phi = M_* = Lambda  => K_th = m_phi^2 M_*^2 = Lambda^4
K_th = Lambda**4
# K(r) = 48 G^2 M^2 / r^6 = 12 r_s^2 / r^6   (Schwarzschild)
# resolver K(r_inf) = K_th
r_inf = (12 * r_s**2 / K_th)**(1/6)   # eV^-1
print(f"\n[3] Umbral K_th = Lambda^4 = {K_th:.3e} eV^4")
print(f"    Radio de influencia r_inf = {eVinv_to_m(r_inf)*1e6:.2f} micrometros")
print(f"    (mas alla de r_inf, K < umbral -> phi NO disparado)")

# ---- 4. separacion media entre PBH si los PBH SON la materia oscura ----
rho_crit_gcm3 = 8.5e-30                       # g/cm^3
Omega_DM = 0.26
rho_DM_gcm3 = Omega_DM * rho_crit_gcm3
# convertir a eV^4:  1 g/cm^3 ...   mejor: numero por cm^3
rho_DM_eV_cm3 = rho_DM_gcm3 * 1e-3 * kg_to_eV  # eV por cm^3
n_PBH_cm3 = rho_DM_eV_cm3 / M_BH               # PBH por cm^3 (cosmic mean)
d_sep_cm = n_PBH_cm3**(-1/3)
print(f"\n[4] Si los PBH masa-Ceres = TODA la DM (densidad cosmica media):")
print(f"    n_PBH = {n_PBH_cm3:.3e} / cm^3")
print(f"    separacion media = {d_sep_cm:.3e} cm = {d_sep_cm/3.086e18:.3e} pc")

# ---- 5. FRACCION DE VOLUMEN disparada ----
r_inf_cm = eVinv_to_m(r_inf) * 100
V_inf_cm3 = (4/3) * np.pi * r_inf_cm**3
f_fill = V_inf_cm3 * n_PBH_cm3
print(f"\n[5] >>> FRACCION DE VOLUMEN con K > umbral <<<")
print(f"    f_fill = {f_fill:.3e}")
print(f"    (fraccion del espacio donde el disparador de Osiris actua)")

# ---- 6. Kretschmann cosmologico de fondo hoy ----
K_cosmo = 24 * H0**4   # de Sitter: K = 24 H^4 (orden de magnitud)
print(f"\n[6] Kretschmann cosmologico de fondo HOY:")
print(f"    K_cosmo ~ 24 H0^4 = {K_cosmo:.3e} eV^4")
print(f"    K_th                = {K_th:.3e} eV^4")
print(f"    K_cosmo / K_th      = {K_cosmo/K_th:.3e}")
print(f"    -> el fondo cosmico {'SUPERA' if K_cosmo>K_th else 'NO alcanza'} el umbral por "
      f"factor {K_th/K_cosmo:.1e}")

print("\n" + "="*64)
print("VEREDICTO")
print("="*64)
print(f" Disparador PBH cubre {f_fill:.1e} del volumen -> "
      f"{'suficiente' if f_fill>0.01 else 'DESPRECIABLE'} para conversion masiva.")
print(f" Fondo cosmico esta {K_th/K_cosmo:.0e}x por debajo del umbral -> "
      f"phi {'puede' if K_cosmo>K_th else 'NUNCA'} es DE en gran escala con este K.")
