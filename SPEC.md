# SPEC — PBL2 (BME 260L): V. cholerae small-intestine model

## Goal
Restructure the three prototype scripts into one small Python package. It should simulate the full 4-state cholera model, sweep parameters, and make the figures for the PBL2 deliverable. It replaces:
- `cholera_equations_model.py`: 1-state N_SI equation with continuous intake, implemented as dataclass + helper functions
- `Prelim Modeling PBL2.py`: 1-state model with a delayed bolus, oscillating Q(t), and sweeps
- `jakes_barebones_cholera_model.ipynb`: earlier 4-state model (N_SI, T, E, V) with ORS

The **source of truth for the equations** is the "Updated Eqns" document (Nov 2026), transcribed below. It supersedes the notebook. The toxin state T is gone, and potassium is a new state.

**Current phase: code skeleton only.** Physiological parameter values are not assigned yet (see *Parameters*).

## Model
States: $y = [N_{SI}, E, V, K^+_{body}]$

$$\frac{dN_{SI}}{dt} = r_{growth} N_{SI}\left(1-\frac{N_{SI}}{N_{SI,max}}\right) - k_{imm} N_{SI} - f_{planktonic}\, Q_{SI}(t)\,\frac{N_{SI}}{V_{lumen}}$$
$$\frac{dE}{dt} = k_{E,max}\frac{N_{SI}}{K_{half,N} + N_{SI}}(1-E) - r_{turnover} E$$
$$\frac{dV}{dt} = \mathrm{ORS}(t) + r_{liquid,intake} - Q_{stool}(t) - U_{liquid}$$
$$\frac{dK^+_{body}}{dt} = r_{K,intake} + c_{K,ORS}\,\mathrm{ORS}(t) - c_{K,stool}\,Q_{stool}(t) - U_K$$

where
- $Q_{SI}(t) = Q_{SI,0} + (Q_{SI,max}-Q_{SI,0})\,E(t)$: secretory outflow from the small intestine
- $Q_{stool}(t) = \max(0,\ Q_{SI}(t) - A_{max})$: what colonic absorption can't take back
- $\mathrm{ORS}(t) = \eta\, Q_{stool}(t)\, u(t-t_R)$

Initial conditions: $N_{SI}(0) = N_{intake}(1-f_{cleared})$, $E(0)=0$, $V(0)=14\,000$ mL, $K^+_{body}(0)\approx 3500$ mmol.

Notes:
- **Intake** is a bolus at $t=0$. There is no $t_{delay}$, so $t=0$ is arrival in the small intestine.
- **Starting assumptions** (defaults, not hard-coded): $r_{liquid,intake}=U_{liquid}$, $r_{K,intake}=U_K$, and $c_{K,stool}$ is constant (independent of $E$). With these, $V$ and $K^+$ stay constant at baseline. After $t_R$, $dV/dt = -(1-\eta)Q_{stool}$ and $dK^+/dt = (\eta c_{K,ORS} - c_{K,stool})Q_{stool}$.
- **Naming:** the doc's $K_{max}$ and $K_N$ are `k_E_max` and `K_half_N` in code, so they aren't confused with potassium.
- **Removed** relative to the prototypes: toxin $T$, `dynamic_VL` (washout always uses constant $V_{lumen}$; $V$ is whole-body fluid), the oscillating $Q(t)$, and $t_{delay}$.

## Package layout
```
PBL2_BME260/
├── SPEC.md
├── CLAUDE.md
├── pyproject.toml           # numpy, scipy, matplotlib; pytest as dev dep
├── src/cholera/
│   ├── params.py            # Params dataclass, units in field comments, validation
│   ├── model.py             # Q_SI_of, Q_stool_of, ors, rhs(t, y, p) — pure functions
│   ├── simulate.py          # simulate(p, t_end) -> Result; handles t_R split and A_max event
│   ├── sweep.py             # sweep(p, name, values) -> list[Result]; 2-D sweep optional
│   └── plots.py             # timecourse(result), sweep_plot(results), no plt.show() inside
├── scripts/
│   ├── run_baseline.py      # writes figures/ timecourse
│   └── run_sweeps.py        # f_planktonic, r_growth, eta, t_R sweeps
├── tests/
│   ├── test_model.py
│   └── test_simulate.py
├── figures/                 # generated, gitignored
└── legacy/                  # move the three original files here unchanged
```

## Parameters
One frozen dataclass, `Params`. Every field has a `# unit, source` comment, and source is `TODO` for now.

| Field | Unit | Meaning |
|---|---|---|
| `r_growth` | 1/h | logistic growth rate |
| `N_SI_max` | cells | carrying capacity |
| `k_imm` | 1/h | linear immune clearance |
| `f_planktonic` | – | free-floating fraction, [0, 1] |
| `V_lumen` | mL | SI lumen volume (constant) |
| `k_E_max` | 1/h | max epithelial activation rate |
| `K_half_N` | cells | half-saturation of activation |
| `r_turnover` | 1/h | epithelial turnover |
| `Q_SI_0`, `Q_SI_max` | mL/h | baseline / max SI outflow |
| `A_max` | mL/h | max colonic absorption |
| `r_liquid_intake`, `U_liquid` | mL/h | fluid intake / urine output |
| `r_K_intake`, `U_K` | mmol/h | K⁺ intake / urinary loss |
| `c_K_ORS`, `c_K_stool` | mmol/mL | K⁺ concentration in ORS / stool |
| `eta` | – | ORS replacement fraction, [0, 1] |
| `t_R` | h | ORS start time |
| `N_intake` | cells | ingested dose |
| `f_cleared` | – | fraction killed by stomach acid, [0, 1] |
| `V0` | mL | initial whole-body fluid (default 14 000, from doc) |
| `K0` | mmol | initial body K⁺ (default 3500, from doc) |

- **No physiological values yet.** Fields other than `V0` and `K0` have no defaults. Tests build their own `Params` from a fixture with simple numbers. Presets (e.g. `PHYSIOLOGICAL`) come later, once values have sources. The prelim script's values ($r_g=1.55$ h⁻¹, $N_{max}=3\times10^{13}$, $V_L=43$ mL, $Q_0=240$ mL/h, $f_{acid}=0.99$, $N_{intake}=10^8$) are candidates.
- **Validation:** fractions in [0, 1]; volumes, capacities and half-saturations > 0; rates and flows ≥ 0; $Q_{SI,max} \ge Q_{SI,0}$.
- `with_(**overrides)` for sweeps.

## Design requirements
- **Units:** time in h, $N$ in cells, volume in mL, flow in mL/h, potassium in mmol. Do not use mol.
- **Pure `rhs`:** no globals, no plotting, no I/O. Clamp $N_{SI} \ge 0$ and $E \in [0,1]$ only in the RHS evaluation, not by mutating solver state.
- **Solver:** `solve_ivp`, LSODA, rtol 1e-8. Use per-state atol, because $N$ spans ~1–1e13 while $E \in [0,1]$ and $V$, $K^+$ are ~10³–10⁴.
- **Discontinuities:** integrate piecewise rather than relying on `max_step`:
  - split at $t_R$ (ORS step);
  - add a `solve_ivp` event for $Q_{SI}(E) = A_{max}$ (the corner in $\max(0,\cdot)$), i.e. $E = (A_{max}-Q_{SI,0})/(Q_{SI,max}-Q_{SI,0})$. Make it terminal and restart from the event state. Skip it when $A_{max} \le Q_{SI,0}$ or $A_{max} \ge Q_{SI,max}$, since no crossing is possible.
- **Results** are a dataclass with `t, N_SI, E, V, K` plus derived `Q_SI, Q_stool, ORS` arrays, so plots never recompute model logic.
- Plot functions return `Figure`. Scripts save to `figures/` at 150 dpi. Scripts read a params object that stays `TODO` until values are assigned.

## Tests (pytest)
1. **Washout only** ($r_{growth}=k_{imm}=k_{E,max}=0$): $N_{SI}$ decays at rate $f_{planktonic}Q_{SI,0}/V_{lumen}$, matching the analytic solution.
2. **Logistic only:** $N_{SI}\to N_{SI,max}$.
3. **Immune only:** $N_{SI}$ decays at rate $k_{imm}$.
4. **No bacteria** ($N_{intake}=0$, $Q_{SI,0} \le A_{max}$): $E$ stays 0, and $V$, $K^+$ stay at $V_0$, $K_0$.
5. `Q_stool_of` is 0 below $A_{max}$, equals $Q_{SI}-A_{max}$ above, and is continuous at the threshold.
6. **ORS** is off before $t_R$. After it (with balanced intake/urine), $dV/dt = -(1-\eta)Q_{stool}$ and $dK^+/dt = (\eta c_{K,ORS}-c_{K,stool})Q_{stool}$. With $\eta=1$, $V$ is constant after $t_R$.
7. The $A_{max}$ event fires when $E$ crosses the threshold, and results match a fine-`max_step` reference solve.
8. Parameter validation raises on bad fractions, zero volume, and $Q_{SI,max} < Q_{SI,0}$.
9. $E$ stays in [0, 1] and $N_{SI} \ge 0$ across the sweep grid.

## Open questions
- **$\eta > 1$:** validated as ≤ 1 for now. Over-replacement would make $V$ rise after $t_R$. Relax the check if that's wanted.
- **$V \le 0$:** with no ORS and large $Q_{stool}$, $V$ can be driven negative. Decide whether to stop the run with a terminal event (e.g. at a fatal-dehydration threshold) or just flag it.
- **Parameter values and sources:** everything in the table except `V0` and `K0`.
- **Oscillating $Q_{SI,0}$:** dropped from the new equations. Re-add as an option only if needed.

## Workflow for Claude Code
1. Read this file.
2. Plan mode: confirm the plan with Jake before writing code.
3. Move the originals to `legacy/`, then scaffold `params` → `model` → `simulate` → tests → plots/scripts.
4. Validate with the tests above. As a cross-check, the 1-state reduction ($k_{imm}=k_{E,max}=0$, so $Q_{SI}=Q_{SI,0}$) run with the prelim values should match `Prelim Modeling PBL2.py` with `Q_amp=0`, shifted by its $t_{delay}$.
