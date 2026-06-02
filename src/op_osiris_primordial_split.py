"""
Filtro de Osiris — ULTIMA prueba: ¿el disparo temporal deja una separacion DM/DE?

Idea de la rendija: en el universo bebe K era gigante en todo el espacio, pero
las regiones un pelin mas densas (delta ~ 1e-5, las semillas de las galaxias)
tenian K un pelin mayor, asi que cruzaron el umbral K_th un instante DESPUES que
el vacio. Si phi rodaba mientras estaba disparado y se congelo al apagarse el
disparador, una region que estuvo disparada mas tiempo rodo mas lejos.

Pregunta honesta: ¿esa diferencia de "rodada" produce un contraste de orden 1
(galaxias necesitan delta ~ 1 ... 1e6) o solo un contraste minusculo ~1e-5?

Unidades naturales. Cuenta sin fiteo.
"""
import numpy as np

print("="*64)
print("OSIRIS — ULTIMA PRUEBA: ¿el disparo temporal separa DM de DE?")
print("="*64)

# ---- 1. contraste de densidad primordial ----
# Fluctuaciones de inflacion (potencial gravitatorio) ~ escala-invariante
delta_primordial = 1e-5   # delta_rho/rho de las semillas de estructura
print(f"\n[1] Contraste de densidad primordial (semillas de galaxias):")
print(f"    delta = delta_rho/rho ~ {delta_primordial:.0e}")

# ---- 2. como se traduce a contraste en K ----
# El escalar de Kretschmann cosmologico va como la densidad al cuadrado: K ~ rho^2
# (los terminos de FRW son ~ (G rho)^2). Entonces:
#   delta_K / K = 2 * delta_rho/rho
delta_K = 2 * delta_primordial
print(f"\n[2] K ~ rho^2  =>  delta_K/K = 2*delta = {delta_K:.0e}")
print(f"    Una zona densa tiene K un {delta_K*100:.4f}% mas alto que el promedio.")

# ---- 3. retraso temporal en cruzar el umbral ----
# K(a) ~ a^-n  (n=6 materia, n=8 radiacion). El disparo se apaga cuando K=K_th.
# Una zona con K un factor (1+delta_K) mayor cruza el umbral a un factor de
# escala (1+delta_a) mayor:  (1+delta_K)(a(1+delta_a))^-n = a^-n
#   => (1+delta_a)^-n = 1/(1+delta_K)  => delta_a ≈ delta_K / n
n_rad = 8.0
delta_a = delta_K / n_rad
print(f"\n[3] K ~ a^-{n_rad:.0f} (radiacion). Retraso en apagarse el disparador:")
print(f"    delta_a/a = delta_K/n = {delta_a:.2e}")
print(f"    (la zona densa sigue disparada un {delta_a*100:.5f}% mas de expansion)")

# ---- 4. cuanto rueda phi de mas en ese tiempo extra ----
# Mientras esta disparado, phi rueda. La distancia rodada es proporcional al
# tiempo (o factor de escala) que paso disparado. Una diferencia relativa
# delta_a en el tiempo disparado => diferencia relativa ~delta_a en lo rodado.
# Por tanto el contraste IMPRESO en phi entre zona densa y vacio es ~ delta_a.
contraste_phi = delta_a
print(f"\n[4] Contraste impreso en phi (denso vs vacio):  ~{contraste_phi:.2e}")

# ---- 5. lo que se NECESITA para formar galaxias ----
contraste_necesario = 1.0   # delta ~ 1 (entrada al regimen no-lineal); galaxias 1e5-1e6
print(f"\n[5] Contraste NECESARIO para que phi se agrupe como DM:  ~{contraste_necesario:.0f}")
print(f"    (las galaxias hoy tienen delta ~ 1e5 ... 1e6)")

deficit = contraste_necesario / contraste_phi
print(f"\n    >>> Osiris imprime {contraste_phi:.0e}, se necesita {contraste_necesario:.0f}")
print(f"    >>> FALTA un factor {deficit:.0e}  (~{np.log10(deficit):.0f} ordenes de magnitud)")

# ---- 6. el clavo final: phi sigue siendo meV-ligero hoy ----
print(f"\n[6] Ademas: tras apagarse el disparador, phi vuelve a ser meV-ligero")
print(f"    en TODO el espacio. Un campo de meV hace free-streaming (se desparrama)")
print(f"    y borra cualquier grumo pequeno. El imprint de {contraste_phi:.0e} se")
print(f"    suaviza aun mas con el tiempo, en vez de crecer.")

print("\n" + "="*64)
print("VEREDICTO FINAL")
print("="*64)
print(f" Disparo temporal imprime contraste ~{contraste_phi:.0e} en phi.")
print(f" Se necesita ~1 para formar estructura. Falta {np.log10(deficit):.0f} ordenes.")
print(f" Y phi meV se desparrama, no se agrupa. -> Osiris NO separa DM de DE.")
