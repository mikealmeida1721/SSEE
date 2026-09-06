"""Gamma del EFT_section: que w_a sale de zeta, y con que factor."""
import sys, os
sys.path.insert(0, "/home/mike/Proyectos/SSEE/src")
from ssee_core import PI, P_SC, W0, WA, IGNIS

I_g = IGNIS
wa_alg = -P_SC / I_g
print(f"  w_a algebraico = -P_sc/I_g = {wa_alg:.12f}")
print(f"  WA en el core             = {WA:.12f}")
print(f"  dif                       = {abs(wa_alg-WA):.2e}\n")

# Continuidad viscosa:
#   rho' + 3H(rho + p + Pi) = 0,  Pi = -3 zeta H
# CPL: p = w rho,  w = w0 + wa(1-a)
# El k-essence puro ya da w0 => la parte viscosa aporta wa(1-a)rho:
#   Pi = wa (1-a) rho
# El paper escribe zeta = Gamma (1-a) rho / H  =>  Pi = -3 Gamma (1-a) rho
print("  Pi que hace falta  =  w_a (1-a) rho")
print("  Pi del paper       = -3 Gamma (1-a) rho")
print("  => Gamma = -w_a/3\n")

G_paper = -P_SC / I_g
G_fix   =  P_SC / (3.0 * I_g)
print(f"  Gamma del paper   = -P_sc/I_g      = {G_paper:+.12f}")
print(f"  Gamma correcto    = +P_sc/(3 I_g)  = {G_fix:+.12f}")
print(f"  cociente                          = {G_paper/G_fix:+.6f}\n")

print("  w_a que devuelve cada uno (w_a = -3 Gamma):")
print(f"    con el del paper  {-3*G_paper:+.12f}   (objetivo {WA:+.12f})")
print(f"    con el corregido  {-3*G_fix:+.12f}   dif {abs(-3*G_fix-WA):.2e}")
