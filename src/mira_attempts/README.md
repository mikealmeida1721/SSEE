# MIRA Derivation Attempts — Archive

This directory contains the 12 scripts produced during the systematic
attempts to derive the **MIRA dynamical mechanism** (OP-8 in
`OPEN_PROBLEMS.md`). All 7 mechanisms tested here have been ruled out;
the scripts are preserved as scientific evidence of *what was tried and
why it failed*, not as production code.

See `VERIFICATION_LEDGER.md` §V-L3-mira and §V-L3-saturacion for the
formal record. See memory `project_veta2_saturation.md` for the
session notes.

## Mechanisms tested and ruled out

| # | Script(s) | Mechanism | Failure mode |
|---|---|---|---|
| 1 | (sound speed clustering — Paper 5) | $c^2_s$ clustering | magnitude (effect 100× too small) |
| 2 | (μ>1 Poisson — P7) | Modified Poisson | structural (forces $\alpha_B=\alpha_M=0$) |
| 3 | (disformal — P8) | Disformal coupling | axiomatic (requires $\psi_{\rm DM}$) |
| 4 | `ssee_mira_mechanism.py`, `ssee_alpha_saturation*.py` | Conformal $\beta_c=-\mathrm{AURA}$ | magnitude+sign+timing |
| 5 | `ssee_mira_saturated_*.py` | Saturated $\alpha=\sqrt{3/(\varphi+3\pi)}$ backward | dynamical connectivity |
| 6 | `ssee_mira_derivative_test.py` | Derivative coupling $(X/M^2)\mathcal{L}_{\rm DM}$ | UV-incompatible scale |
| 7 | `ssee_mira_phimde_forward_test.py` | Forward from $\varphi$-MDE | matter drainage in 7 e-folds |

Additional diagnostic scripts:
- `ssee_mira_decay_test.py` — algebraic-$\tau$ scan
- `ssee_mira_projection_test.py` — 3D→2D projection hypothesis
- `ssee_mira_saturated_diagnostic.py` — R(z=1100) signature matrix

## What this means

Seven candidate mechanisms have been formally ruled out. The MIRA
algebraic value $(3\varphi+\pi)/4 \simeq 1.9989$ is now Postulate~M
(see Paper~1 §1.3). The *value* is structural; the *dynamical
mechanism* that realises this ratio in a Friedmann background remains
open (OP-8).
