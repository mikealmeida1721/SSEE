"""REGISTRO DE REGLAS — la única declaración de qué vigila el guardián.

POR QUÉ EXISTE (2026-07-29). Mike: «que el guardián me dé verde no significa que
esté bien, y eso es preocupante, porque lo creamos para tener la confianza de que
está bien». Tenía razón: un VERDE demostraba que nada saltó, no que algo habría
saltado. Y la primera prueba de mutación destapó un valor sencillamente MAL que
pasaba verde porque la exención de una regla caía justo en el hueco de otra.

CUATRO MANERAS DE ESTAR VERDE SIN ESTAR BIEN, y qué las detecta:

  vacía        la regla no examina ningún candidato    → caso de MUTACIÓN obligatorio
  enmascarada  la exención de otra tapa el hallazgo     → cada exención dice QUIÉN cubre
  duplicada    dos reglas, misma intención y ámbito     → colisión (intencion, ambito)
  huérfana     regla añadida sin poder verificarla      → sin entrada aquí, no entra

CONTRATO PARA AÑADIR UNA REGLA (esto es el método que faltaba):

  1. ¿Existe ya una regla con esta misma `intencion` y `ambito`? Si sí, se REFINA
     esa. No se crea una nueva «igual pero con cambios».
  2. Se declara aquí: intencion, ambito, exenciones (cada una con `cubierta_por`)
     y al menos un caso de `mutacion`.
  3. Se corre `mutacion_guardian.py`: el caso debe poner el guardián en ROJO.
     Si no enrojece, la regla no sirve todavía — no se commitea.
  4. `meta_guardian.py` comprueba 1–3 y que ninguna capa del código quede sin
     entrada aquí (ni al revés).

`ambito` responde «¿sobre qué mira?»; `intencion`, «¿qué exige?». Dos reglas
pueden compartir una de las dos, nunca las dos.
"""

# Anclas de mutación: texto REAL del Paper 1 sobre el que se inyecta el defecto.
# Si un ancla deja de existir, la prueba avisa en vez de pasar en silencio.
TEX_MUTACION = "manuscript/SSEE_Paper1_Framework.tex"

REGLAS = {
    # ── capas de FÍSICA ──────────────────────────────────────────────────────
    # No viven en un manuscrito: afirman identidades algebraicas. Se prueban
    # tocando el archivo donde vive lo que afirman.
    "canon": dict(
        capa="Fuente canónica — ssee_core.py",
        intencion="nucleo-vs-recomputacion",
        ambito="src/ssee_core.py contra el álgebra rehecha por el guardián",
        # Es el eslabón que conecta las capas de física con el MODELO: la Capa 2
        # compara su propia recomputación contra literales; sola no notaría que
        # ssee_core derivó. Esta sí, porque rehace el álgebra por su cuenta.
        exenciones=[],
        archivo="src/ssee_core.py",
        mutacion=[("una constante del núcleo alterada",
                   "W0          = -T_R / M_V", "W0          = -T_R / M_V * 1.0001")],
    ),
    "V-L2": dict(
        capa="2 — parámetros cosmológicos derivados",
        intencion="identidad-algebraica",
        ambito="álgebra de fondo (w0, wa, Ω_DE, Ω_m, H_alg, n_s, α_K)",
        # La mutación prueba VIVEZA —que la comprobación compara de verdad y está
        # clavada a un número— no que el modelo sea correcto. El puente al modelo
        # lo pone «canon»; el puente al registro publicado, R20/R26.
        exenciones=[],
        archivo="src/verificacion/ssee_verify.py",
        # Historia de este caso, que enseña más que el caso: la primera versión
        # movía el literal 7.7e-7 y NO enrojecía, porque la capa toleraba 1e-6
        # absoluto. Al medir el residuo REAL (peor: 3.9e-9) resultó que la holgura
        # era 256× — y dentro cabía un literal equivocado: r estaba escrito
        # 0.0081306227 en vez de 0.0081306188. Se corrigió el literal y se apretó
        # la tolerancia a 1e-9. Ahora basta con mover el último dígito.
        mutacion=[("el valor esperado de w0 movido en su última cifra",
                   '"V-L2-01 w0":        (-Tr / Mv,                      -0.8399497713)',
                   '"V-L2-01 w0":        (-Tr / Mv,                      -0.8399497813)')],
    ),
    "V-L3": dict(
        capa="3 — mecanismos y derivaciones",
        intencion="identidad-algebraica",
        ambito="mecanismos (n_s/r, f_screen, m_φ, EFT, IS, dos sectores)",
        exenciones=[],
        archivo="src/verificacion/ssee_verify.py",
        prefijos=["V-L3"],
        # Se muta la AFIRMACIÓN, no la etiqueta: un primer caso añadía «[MUT]» al
        # nombre del check y por supuesto no cambiaba nada — el rótulo no es lo
        # que se comprueba. Y un segundo rompía f_screen, que resultó ser un
        # check de la Capa 2 («L2 identidad»), no de la 3. Lo dijo la atribución,
        # que existe justamente para eso.
        mutacion=[("la identidad EFT lambda^2 = 3 Om_m,dyn rota",
                   "abs(lam_eft ** 2 - 3 * Om_m_dyn) < 1e-12",
                   "abs(lam_eft ** 2 - 4 * Om_m_dyn) < 1e-12")],
    ),
    "V-L4": dict(
        capa="4 — confrontaciones con datos",
        intencion="dato-vs-prediccion",
        ambito="σ8/S8, tensiones KiDS, ΔBIC de P2/P3/P6",
        exenciones=[],
        archivo="CANONICAL_VALUES.yaml",
        prefijos=["V-L4"],
        # Se muta el YAML, no el guardián: desde 2026-07-29 estas comprobaciones
        # LEEN el valor de la fuente en vez de llevarlo escrito. Antes σ8 se
        # asignaba «= 0.8335» y se comprobaba «abs(σ8 − 0.8335) < 1e-2» — una
        # tautología imposible de fallar. La mutación es la única forma en que
        # eso salía a la luz: el check estaba VERDE y no miraba nada.
        mutacion=[("el techo σ8 del YAML alterado",
                   "sigma8_single_ceiling: 0.8335",
                   "sigma8_single_ceiling: 0.7000")],
    ),
    "DESI": dict(
        capa="DESI — procedencia BAO: csv == DR2 Tabla 4 oficial",
        intencion="procedencia-del-dato",
        ambito="data/raw/desi_dr2_bao.csv contra la tabla oficial 2503.14738",
        exenciones=[],
        archivo="data/raw/desi_dr2_bao.csv",
        prefijos=["DESI"],
        mutacion=[("un punto BAO alterado respecto a la tabla oficial",
                   "0.295", "0.296")],
    ),
    "R30": dict(
        capa="R30 — política de redondeo",
        intencion="redondeo-correcto",
        ambito="prosa .tex, anclas símbolo↔patrón",
        exenciones=[("valores con menos de 4 decimales (no discriminan)", None),
                    ("medidas de Planck/SH0ES, que no son predicciones", None)],
        mutacion=[("valor anclado que no es el redondeo correcto",
                   "$\\Omega_{m,\\rm CMB}=\\omega_m/h^2=0.308881$ at the global anchor",
                   "$\\Omega_{m,\\rm CMB}=\\omega_m/h^2=0.308999$ at the global anchor")],
    ),
    "R41": dict(
        capa="R41 — coherencia de precisión",
        intencion="una-precision-por-simbolo",
        ambito="prosa .tex, cantidades CON unidades",
        exenciones=[("precisión declarada en el texto («to three decimals»)", None)],
        mutacion=[("misma cantidad dimensional con dos precisiones",
                   "$H_0 = 3(\\varphi+\\pi)^2\\,\\kmsu \\approx 67.962\\,\\kmsu$",
                   "$H_0 = 3(\\varphi+\\pi)^2\\,\\kmsu \\approx 67.9621\\,\\kmsu$")],
    ),
    "R37": dict(
        capa="R37 — igualdades de constantes SSEE a 6 decimales",
        intencion="precision-igualdad",
        ambito="prosa .tex, constantes de _C37",
        exenciones=[("cantidades con unidades", "R41")],
        # El caso va en PROSA, no en una fila de tabla: la primera versión usaba
        # «$w_0 = -T_r/M_v$ & $-0.83995$» y la prueba de atribución avisó de que
        # quien lo cazaba era R38 —cuyo ámbito son las tablas— y no R37. Un
        # defecto detectado por la regla equivocada deja a la propia sin probar.
        mutacion=[("constante SSEE en prosa a 5 decimales",
                   "\\times0.96556 = 0.119514$",
                   "\\times0.96556 = 0.11951$")],
    ),
    "R38": dict(
        capa="R38 — igualdades que cruzan celdas de tabla",
        intencion="precision-igualdad",
        ambito="filas de tabla .tex",
        exenciones=[("≈/≃/∼ en la columna 1", "R40"),
                    ("unidad pegada al valor", "R41"),
                    ("\\ldots explícito", None)],
        mutacion=[("valor MAL, error 4e-4 (el que se escapó)",
                   "$n_s = 1-\\varphi^{-7}$ & $0.965558$",
                   "$n_s = 1-\\varphi^{-7}$ & $0.965123$"),
                  ("valor MAL, error 1e-5 (último dígito)",
                   "$n_s = 1-\\varphi^{-7}$ & $0.965558$",
                   "$n_s = 1-\\varphi^{-7}$ & $0.965548$"),
                  ("valor MAL, error 3e-2 (grosero)",
                   "$n_s = 1-\\varphi^{-7}$ & $0.965558$",
                   "$n_s = 1-\\varphi^{-7}$ & $0.935558$"),
                  ("valor MAL en otra fila y otra fórmula",
                   "$w_0 = -T_r/M_v$ & $-0.839950$",
                   "$w_0 = -T_r/M_v$ & $-0.812340$")],
    ),
    "R40": dict(
        capa="R40 — tipo de relación: «=» exacto vs «≈» truncado",
        intencion="tipo-de-relacion",
        ambito="prosa .tex, símbolo↔fórmula↔decimal",
        exenciones=[],
        mutacion=[("igualdad exacta escondida bajo un «≈»",
                   "\\item $\\boldsymbol{\\Omega} = \\pi+\\varphi \\approx 4.7596$",
                   "\\item $\\boldsymbol{\\Omega}$ ($\\pi+\\varphi \\approx 4.7596$)")],
    ),
    "R42": dict(
        capa="R42 — tipo dimensional: número puro vs cantidad física",
        intencion="tipo-dimensional",
        ambito="prosa .tex, H_0 y cocientes adimensionales",
        exenciones=[("mención negada («we do not write it as…»)", None)],
        mutacion=[("igualdad sin unidad",
                   "anchor $H_0=3(\\varphi+\\pi)^2\\,\\kmsu$ (derived, Paper~9)",
                   "anchor $H_0=3(\\varphi+\\pi)^2$ (derived, Paper~9)"),
                  ("CONFLICTO: ¿la exención de mención enmascara?",
                   "anchor $H_0=3(\\varphi+\\pi)^2\\,\\kmsu$ (derived, Paper~9)",
                   "This is not a fitted quantity: the anchor "
                   "$H_0=3(\\varphi+\\pi)^2$ (derived, Paper~9)")],
    ),
    "R43": dict(
        capa="R43 — potencias de φ: el decimal tras «=» es el valor",
        intencion="precision-igualdad",
        ambito="prosa .tex, expresiones φ^n",
        exenciones=[("≈/≃ (truncamiento declarado)", "R30")],
        mutacion=[("potencia de φ truncada con «=»",
                   "$2\\varphi^7=58.068884$", "$2\\varphi^7=58.07$")],
    ),
    "R44": dict(
        capa="R44 — constantes de la lectura con «=» a 6 decimales",
        intencion="precision-igualdad",
        ambito="prosa .tex, constantes incorporadas por la lectura",
        exenciones=[("documentos aún no leídos (deuda contada, sólo baja)", None)],
        mutacion=[("Ω_m,dyn a 3 decimales",
                   "$\\Omega_{m,\\rm dyn}=0.160050$ (DESI)",
                   "$\\Omega_{m,\\rm dyn}=0.160$ (DESI)"),
                  ("SOLAR²·K_v a 2 decimales",
                   "\\mathrm{KRYSTOS}_V=594.279999$",
                   "\\mathrm{KRYSTOS}_V=594.28$")],
    ),
    "R45": dict(
        capa="R45 — estado de los OP: la prosa concuerda con OPEN_PROBLEMS.md",
        intencion="coherencia-con-registro",
        ambito="prosa .tex, menciones OP-N",
        exenciones=[("documentos aún no leídos (deuda contada)", None)],
        mutacion=[("OP resuelto citado como abierto",
                   "with the remaining open problems OP-9 and OP-11",
                   "with the remaining open problems OP-9/11/14")],
    ),
}

# Capas que EXISTEN en el código y todavía NO tienen entrada arriba. No es una
# lista de perdón: es la deuda visible, y sólo puede bajar. Mientras una capa
# esté aquí, su VERDE no está demostrado — puede ser verde por vacío.
SIN_COBERTURA = [
    "Capa Memorias — coherencia de las 3 memorias",
    "Capa Procedencia — valores de pipeline vs log committeado",
    "Capa R33 — logs vs núcleo: ningún log activo con constante retirada",
    "Capa R39 — compuerta de publicación: conoce todas las portadas",
    "Capa R34 — fuentes vs núcleo: ningún .py hardcodea constante retirada",
    "Capa R35 — artefacto vs fuente: ningún log más viejo que su script",
    "Capa Manuscritos — reglas R1 (cronología) y R2 (multidominio)",
    "Capa Archivo — bitácora cubre cada cajón archivado",
    "Capa Diccionario — integridad de nodos (script ⊆ maestro)",
    "Capa R20 — anclas observacionales (código ⊆ CANONICAL_VALUES.yaml)",
    "Capa R25 — parametrización ω_m absoluto (no congelar Ω_m en el MCMC)",
    "Capa R26 — procedencia declarada de los canónicos",
]
