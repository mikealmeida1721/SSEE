"""Prueba de MUTACIÓN del guardián.

Inyecta un defecto CONOCIDO en un manuscrito, corre el guardián y comprueba que
se pone ROJO; después restaura el archivo y verifica que quedó idéntico.

POR QUÉ EXISTE (2026-07-29). Mike preguntó lo que había que preguntar: si las
reglas se acumulan, ¿cómo sabemos que cada una funciona, que no se estorban entre
sí, y que no estamos poniendo parches donde tocaba refinar? Un VERDE prueba que
nada saltó; NO prueba que algo habría saltado. Esto último sólo se demuestra
rompiendo el documento a propósito.

Y la primera corrida encontró el agujero: un valor sencillamente MAL
(«n_s = 1−φ⁻⁷ & 0.965123», a 4.4e-4 del exacto) pasaba VERDE. R38 lo descartaba
por su ventana de identificación —lo leía como «esta celda no habla de n_s»— y
R30 no lo veía porque sus patrones esperan «fórmula = valor» y el valor vivía en
otra celda. Dos reglas, cada una con su motivo, y el error grueso en medio.
Se refinó R38 (identificación DECLARADA por el par símbolo|fórmula) en vez de
añadir una regla más.

REGLA DE USO: toda regla nueva entra con su caso aquí. Si añadir el caso no
enrojece al guardián, la regla no sirve todavía."""
import subprocess, pathlib, sys
REPO = pathlib.Path("/home/mike/Proyectos/SSEE")
TEX  = REPO / "manuscript/SSEE_Paper1_Framework.tex"
orig = TEX.read_text()

CASOS = [
 ("R42 igualdad sin unidad",
  "anchor $H_0=3(\\varphi+\\pi)^2\\,\\kmsu$ (derived, Paper~9)",
  "anchor $H_0=3(\\varphi+\\pi)^2$ (derived, Paper~9)"),
 ("R43 potencia de φ truncada con «=»",
  "$2\\varphi^7=58.068884$", "$2\\varphi^7=58.07$"),
 ("R44 Ω_m,dyn a 3 decimales",
  "$\\Omega_{m,\\rm dyn}=0.160050$ (DESI)", "$\\Omega_{m,\\rm dyn}=0.160$ (DESI)"),
 ("R44 SOLAR²·K_v a 2 decimales",
  "\\mathrm{KRYSTOS}_V=594.279999$", "\\mathrm{KRYSTOS}_V=594.28$"),
 ("R45 OP resuelto citado como abierto",
  "with the remaining open problems OP-9 and OP-11",
  "with the remaining open problems OP-9/11/14"),
 ("R37 constante SSEE a 5 decimales",
  "$w_0 = -T_r/M_v$ & $-0.839950$", "$w_0 = -T_r/M_v$ & $-0.83995$"),
 ("valor MAL, error 4e-4 (el que se escapó)",
  "$n_s = 1-\\varphi^{-7}$ & $0.965558$", "$n_s = 1-\\varphi^{-7}$ & $0.965123$"),
 ("valor MAL, error 1e-5 (último dígito)",
  "$n_s = 1-\\varphi^{-7}$ & $0.965558$", "$n_s = 1-\\varphi^{-7}$ & $0.965548$"),
 ("valor MAL, error 3e-2 (grosero)",
  "$n_s = 1-\\varphi^{-7}$ & $0.965558$", "$n_s = 1-\\varphi^{-7}$ & $0.935558$"),
 ("valor MAL en w_0 (otra fila, otra fórmula)",
  "$w_0 = -T_r/M_v$ & $-0.839950$", "$w_0 = -T_r/M_v$ & $-0.812340$"),
 ("CONTROL: valor CORRECTO no debe disparar",
  "$n_s = 1-\\varphi^{-7}$ & $0.965558$", "$n_s = 1-\\varphi^{-7}$ & $0.965558$"),
 ("CONFLICTO: R42 enmascarado por la exención de mención",
  "anchor $H_0=3(\\varphi+\\pi)^2\\,\\kmsu$ (derived, Paper~9)",
  "This is not a fitted quantity: the anchor $H_0=3(\\varphi+\\pi)^2$ (derived, Paper~9)"),
]

def verde():
    r = subprocess.run([sys.executable, str(REPO/"src/verificacion/ssee_verify.py")],
                       capture_output=True, text=True, cwd=REPO)
    return "VERDE" in r.stdout.split("=" * 40)[-1]

print("línea base:", "VERDE" if verde() else "ROJO")
print()
for nombre, viejo, nuevo in CASOS:
    if viejo not in orig:
        print(f"  [?] {nombre:52s} — patrón ancla no encontrado"); continue
    TEX.write_text(orig.replace(viejo, nuevo, 1))
    ok = not verde()          # queremos que se ponga ROJO
    TEX.write_text(orig)
    print(f"  [{'DETECTA' if ok else ' PASA  '}] {nombre}")
assert TEX.read_text() == orig, "¡el archivo no volvió a su estado original!"
print("\narchivo restaurado idéntico ✓")
