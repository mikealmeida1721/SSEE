# Sondeos de la partícula φ-DM — retirados el 2026-08-01

Estos nueve scripts investigaban la partícula `m_φ = 40.70 eV` y el segundo
sector φ-DM: su cota de masa, su free-streaming, su huella en KiDS, su papel en
el cierre de Ω. **Todos ellos quedaron sin objeto** cuando la partícula se
retiró: la resta que definía su densidad mezclaba una densidad medida con
`1+w₀`, que es un número de la ecuación de estado, y la tensión S₈ que la
motivaba no existe contra el dato crudo (MCMC R3: 0.11σ).

**Por qué se archivan y no se borran.** Son el registro de cómo se puso a
prueba la hipótesis, incluida la medición que la excluyó (`cota_masa.py` /
`fig_cota.py`: la cizalla cruda exige `m_φ > 70.3 eV`, muy por encima de los
40.70 predichos). Ese trabajo es la razón por la que la retirada está
justificada, y borrarlo dejaría la conclusión sin su prueba.

**Por qué se archivan y no se marcan.** Seguían dentro de `src/`, que es la
superficie viva que barre el guardián (R60), con líneas como `M_PHI = 40.70` y
`OM_PHI = 0.14889` escritas como valores en uso — 24 sitios. Marcarlas una a
una habría dejado nueve scripts muertos en el cajón de los vivos. `archive/`
es el cajón de lo retirado, y el barrido lo excluye por eso.

Movidos el 2026-09-19. Ninguno era importado por código vivo (comprobado).
