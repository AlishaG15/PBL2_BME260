# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

BME 260L PBL2 course project: an ODE model of *V. cholerae* infection. Its states are bacterial load in the small-intestine lumen, the activated epithelial fraction, whole-body fluid volume, and whole-body potassium, with ORS treatment. The current equations are in `SPEC.md` and supersede all three prototypes. Right now it is three standalone prototypes. `SPEC.md` describes a planned restructure into a `src/cholera/` package. **Read `SPEC.md` before any structural work.** It defines the target layout, design requirements, tests, and unresolved modeling questions. Its workflow says to confirm the layout and the open-question defaults with Jake in plan mode before writing code.

## Running

No packaging, dependency file, or tests exist yet. Dependencies are numpy, scipy, and matplotlib.

```bash
python3 "Prelim Modeling PBL2.py"      # writes cholera_si_timecourse.png and cholera_si_sweep.png to cwd
python3 cholera_equations_model.py     # prints a single dN_SI/dt evaluation
jupyter notebook jakes_barebones_cholera_model.ipynb
```

Once the SPEC package exists, the tests will be run with `pytest` (for a single test, `pytest tests/test_model.py::test_name`).

## The three prototypes and how they relate

They are now legacy references only: the equations in `SPEC.md` (no toxin state T, a new K⁺ state, bolus at t = 0, linear immune term) supersede them.

- **`jakes_barebones_cholera_model.ipynb`**: the older 4-state model `y = [N_SI, T, E, V]`. Outflow `Q(E) = Q0 + (Q_max − Q0)·E` couples toxin effect to washout and volume loss. ORS is `eta·Q·u(t − t_R)`. With `dynamic_VL`, washout uses the state `V` in place of the constant `V_L`. Its parameter values are unitless placeholders.
- **`Prelim Modeling PBL2.py`**: 1-state `N_SI` model. Its physiological values (hours, cells, mL, mL/h) are candidates for later presets. Bacteria arrive as a bolus `N_intake·(1 − f_acid)` at `t_delay`, and `N_SI = 0` before that. It has an optional sinusoidal `Q(t)` and a one-parameter sweep. Its docstring's mention of mol is stale and should be ignored.
- **`cholera_equations_model.py`**: 1-state model as a frozen dataclass plus helper functions, with a continuous intake term (superseded by the t = 0 bolus).

Naming differs across files for the same quantities: `f_planktonic` / `f_free`, `N_max` / `N_SI_max`, `V_lumen` / `V_L`. The new code uses SPEC.md's names: `f_planktonic`, `N_SI_max`, `V_lumen`, and `k_E_max` / `K_half_N` for the doc's K_max / K_N, so they aren't confused with potassium.

## Numerical gotchas (from SPEC.md)

- `N_SI` spans ~1 to 1e13 cells, while E is in [0, 1] and V, K⁺ are ~10³–10⁴. A single scalar `atol` is wrong for at least one state, so use per-state atol.
- Two discontinuities: the ORS step at `t_R` (split the integration there) and the corner in `Q_stool = max(0, Q_SI − A_max)` (terminal event when Q_SI crosses A_max, then restart). Don't fall back to `max_step`.
- Physiological parameter values are deliberately unassigned during the skeleton phase. Don't invent them; leave sources as `TODO`.

## Repo notes

- PNGs in the root are committed outputs of the prototypes. The SPEC plans to move generated figures to a gitignored `figures/` and move the prototypes unchanged into `legacy/`.
- The cross-check for the refactor is the 1-state reduction (k_imm = k_E_max = 0) run with the prelim values. It should match the prelim script with `Q_amp=0`, shifted by its 20-min delay.
