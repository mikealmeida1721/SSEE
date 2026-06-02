"""
OP-9 ruta 1 — Misalignment: ¿un campo de meV producido FRIO da Omega_DM ~ 0.26?

Mecanismo (estandar, sin inventar nada):
  - Campo phi con V = 1/2 m^2 phi^2. Al principio H >> m: la friccion de Hubble
    lo tiene CONGELADO en un valor inicial phi_i (como un pendulo sostenido).
  - Cuando H baja hasta H ~ m, phi se suelta y empieza a OSCILAR.
  - Desde ahi se comporta como materia oscura FRIA: rho ~ a^-3, sin velocidad
    termica => NO hace free-streaming => SI forma galaxias.

Pregunta honesta: ¿que valor inicial phi_i se necesita para Omega_DM=0.26?
¿Es un numero natural de SSEE o hay que ajustarlo a mano (= parametro libre)?

Unidades naturales (hbar=c=1), eV. Cuenta sin fiteo.
"""
import numpy as np

print("="*64)
print("OP-9 RUTA 1 — Misalignment de un campo de meV (cuenta honesta)")
print("="*64)

# ---- constantes ----
M_Pl_red = 2.435e27      # masa de Planck reducida [eV]
M_Pl_np  = 1.2209e28     # masa de Planck NO reducida [eV]
H0       = 1.44e-33      # eV
T0       = 2.35e-4       # temperatura CMB hoy [eV] (2.725 K)
Omega_DM = 0.26
m_phi    = 8.81e-3       # CANDIDATO: m_phi = Lambda_SSEE = 8.81 meV [eV]

# ---- 1. densidad de DM hoy que hay que reproducir ----
rho_crit = 3 * H0**2 * M_Pl_red**2          # eV^4
rho_DM   = Omega_DM * rho_crit
print(f"\n[1] Objetivo:  rho_crit = {rho_crit:.3e} eV^4  (rho_crit^1/4 = {rho_crit**0.25*1e3:.2f} meV)")
print(f"    rho_DM (Omega=0.26) = {rho_DM:.3e} eV^4")

# ---- 2. cuando empieza a oscilar:  H(T_osc) = m_phi ----
# radiacion: H = 1.66 sqrt(g*) T^2 / M_Pl_np   =>   T_osc = sqrt(m M_Pl /(1.66 sqrt g*))
g_star = 106.75          # todos los grados de libertad del SM (T_osc resultara ~TeV)
T_osc = np.sqrt(m_phi * M_Pl_np / (1.66 * np.sqrt(g_star)))
print(f"\n[2] Oscila cuando H = m_phi:")
print(f"    T_osc = {T_osc:.3e} eV = {T_osc/1e12:.2f} TeV   (g* = {g_star})")
print(f"    (mucho antes de BBN; el campo lleva FRIO desde entonces)")

# ---- 3. dilucion por expansion (conservacion de entropia) ----
g_s_osc = 106.75
g_s_0   = 3.91
dilucion = (g_s_0 / g_s_osc) * (T0 / T_osc)**3      # (a_osc/a_0)^3
print(f"\n[3] Dilucion (a_osc/a_0)^3 = {dilucion:.3e}")

# ---- 4. phi_i NECESARIO para dar Omega_DM=0.26 ----
# rho_DM = (1/2 m^2 phi_i^2) * dilucion   =>   phi_i = sqrt(2 rho_DM /(m^2 dilucion))
phi_i = np.sqrt(2 * rho_DM / (m_phi**2 * dilucion))
print(f"\n[4] >>> phi_i NECESARIO <<<")
print(f"    phi_i = {phi_i:.3e} eV = {phi_i/1e9:.3e} GeV")
print(f"    phi_i / M_Pl_red = {phi_i/M_Pl_red:.3e}")

# ---- 5. ¿es un numero natural de SSEE? comparacion con escalas conocidas ----
print(f"\n[5] ¿Es phi_i un valor natural? Comparacion con escalas:")
print(f"    M_Pl_red                 = {M_Pl_red:.3e} eV   (phi_i/M_Pl = {phi_i/M_Pl_red:.2e})")
print(f"    sqrt(M_Pl * H0)          = {np.sqrt(M_Pl_red*H0):.3e} eV")
print(f"    geometrica sqrt(M_Pl*meV)= {np.sqrt(M_Pl_red*m_phi):.3e} eV")
print(f"    (M_Pl^2 * m_phi)^(1/3)   = {(M_Pl_red**2 * m_phi)**(1/3):.3e} eV")
print(f"    sqrt(M_Pl * 1 GeV)       = {np.sqrt(M_Pl_red*1e9):.3e} eV")

print("\n" + "="*64)
print("LECTURA")
print("="*64)
print(f" - El campo de meV producido por misalignment es FRIO: oscila desde T~TeV,")
print(f"   rho ~ a^-3, sin velocidad termica. NO free-streaming -> SI forma galaxias.")
print(f"   => El miedo de OP-9 (meV se desparrama) se DISUELVE para produccion fria.")
print(f" - Costo: hace falta phi_i ~ {phi_i/M_Pl_red:.1e} M_Pl = {phi_i/1e9:.1e} GeV.")
print(f"   Si ese numero sale de (phi,pi) sin ajuste -> 0 parametros libres, GANAMOS.")
print(f"   Si hay que ponerlo a mano -> es 1 parametro (mejor que la formula rota, pero no 0).")
