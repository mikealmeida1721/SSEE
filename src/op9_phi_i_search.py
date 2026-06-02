"""
OP-9 ruta 1 — ¿el phi_i del misalignment sale de (phi,pi) sin ajuste?

Objetivo (de op9_misalignment_relic.py):
    phi_i = 2.857e21 eV = 1.17e-6 M_Pl_red    para dar Omega_DM=0.26 con m_phi=8.81 meV.

Pregunta honesta: ¿hay una combinacion NATURAL de las constantes SSEE
(phi, pi, Omega, beta, KAL, ...) y/o escalas fisicas (M_Pl, H0, Lambda, rho_c^1/4)
que reproduzca phi_i SIN exponentes feos? Si si -> 0 parametros. Si no -> 1 parametro.

Referee hostil: solo cuenta como "hallazgo" si el error es < 1% Y los exponentes
son enteros chicos / razones simples. Forzar con lambda~60000 = numerologia, se rechaza.
"""
import numpy as np
import itertools

print("="*64)
print("OP-9 RUTA 1 — ¿phi_i es natural en SSEE? (busqueda honesta)")
print("="*64)

# ---- el objetivo, recomputado limpio ----
M_Pl_red = 2.435e27
M_Pl_np  = 1.2209e28
H0       = 1.44e-33
T0       = 2.35e-4
Omega_DM = 0.26
m_phi    = 8.81e-3
g_star   = 106.75
g_s_0    = 3.91

rho_crit = 3 * H0**2 * M_Pl_red**2
rho_DM   = Omega_DM * rho_crit
T_osc    = np.sqrt(m_phi * M_Pl_np / (1.66 * np.sqrt(g_star)))
dilucion = (g_s_0 / g_star) * (T0 / T_osc)**3
phi_i    = np.sqrt(2 * rho_DM / (m_phi**2 * dilucion))

print(f"\n[OBJETIVO]  phi_i = {phi_i:.4e} eV")
print(f"            phi_i / M_Pl_red = {phi_i/M_Pl_red:.4e}")
print(f"            log10(phi_i)     = {np.log10(phi_i):.4f}")

# ---- constantes SSEE adimensionales ----
phi  = (1 + np.sqrt(5)) / 2
pi   = np.pi
Omega = pi + phi
beta  = (pi + phi)/2
KAL   = beta + pi
P_sc  = Omega + phi
Kv    = phi + pi + Omega
Tr    = 3*(phi + beta)
Mv    = phi + pi + Kv

ssee_consts = {
    "phi": phi, "pi": pi, "Omega": Omega, "beta": beta, "KAL": KAL,
    "P_sc": P_sc, "Kv": Kv, "Tr": Tr, "Mv": Mv, "sqrt5": np.sqrt(5),
}

# ---- escalas fisicas dimensionales (eV) ----
scales = {
    "M_Pl_red": M_Pl_red, "M_Pl_np": M_Pl_np, "H0": H0,
    "Lambda": m_phi, "rho_c_4th": rho_crit**0.25, "T0": T0,
}

# ====================================================================
# Estrategia 1: phi_i = (const_adim) * (escala_dim)^a * (escala_dim2)^b
# con a,b en {-2,-1,-1/2,1/2,1,2} y a+b... no exigimos dim, solo eV total.
# Pero queremos que el RESULTADO tenga unidades de eV (phi_i es eV).
# Las escalas ya son eV, asi que a+b+... = 1 para que sea eV^1.
# ====================================================================
print("\n" + "="*64)
print("[1] phi_i = C_adim * escala1^a * escala2^b   (suma exponentes = 1)")
print("="*64)

target = phi_i
exps = [-2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2]
hits = []

scale_items = list(scales.items())
for (n1, v1), (n2, v2) in itertools.product(scale_items, repeat=2):
    for a in exps:
        for b in exps:
            if abs(a + b - 1.0) > 1e-9:   # exigir eV^1 total
                continue
            dim_part = v1**a * v2**b
            C_needed = target / dim_part
            # ¿C_needed coincide con alguna constante SSEE o razon simple?
            for cn, cv in ssee_consts.items():
                for p in [-2,-1,1,2]:
                    pred = cv**p * dim_part
                    err = abs(pred/target - 1)
                    if err < 0.02:
                        hits.append((err, f"{cn}^{p} * {n1}^{a} * {n2}^{b}", pred))
            # tambien C_needed crudo, para inspeccion
            if 0.5 < C_needed < 50:   # constante "natural" de orden 1-50
                hits.append((abs(0), f"[C={C_needed:.4f}] * {n1}^{a} * {n2}^{b}  (C crudo)", target))

hits = sorted(set((round(h[0],6), h[1]) for h in hits), key=lambda x: x[0])
print(f"\nCandidatos con error < 2% (o C crudo de orden 1-50):")
for err, desc in hits[:25]:
    print(f"   err={err*100:5.2f}%   phi_i = {desc}")

# ====================================================================
# Estrategia 2: razon phi_i/M_Pl = 1.17e-6 como potencia de phi o combo
# ====================================================================
print("\n" + "="*64)
print("[2] phi_i/M_Pl = 1.17e-6  como potencia de constantes")
print("="*64)
ratio = phi_i / M_Pl_red
print(f"    ratio = {ratio:.4e},  log10 = {np.log10(ratio):.4f}")
for cn, cv in ssee_consts.items():
    n = np.log(ratio) / np.log(cv)
    print(f"    ratio = {cn}^({n:.3f})    [{'ENTERO?' if abs(n-round(n))<0.05 else 'no'}]")

# geometrico: ratio ~ (H0/M_Pl)^x ?
for cn, cv in [("H0/M_Pl", H0/M_Pl_red), ("Lambda/M_Pl", m_phi/M_Pl_red),
               ("T0/M_Pl", T0/M_Pl_red), ("rho_c4/M_Pl", rho_crit**0.25/M_Pl_red)]:
    x = np.log(ratio)/np.log(cv)
    print(f"    ratio = ({cn})^({x:.4f})    [{'ENTERO?' if abs(x-round(x))<0.05 else 'no'}]")

print("\n" + "="*64)
print("LECTURA (referee)")
print("="*64)
print(" - Si arriba hay un match con error <1% y exponentes simples -> GOLAZO (0 params).")
print(" - Si solo hay 'C crudo' raros o exponentes feos -> phi_i es 1 parametro libre.")
print(" - phi_i/M_Pl ~ 1e-6 es el numero a explicar. Ojo: NO inventar exponente.")
