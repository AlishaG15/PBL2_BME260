

"""Units
-----
Time in hours, N in cells (the notes use mol; convert with N / AVOGADRO).
Q in mL/hour, V in mL.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

avogadrosnum = 6.022e23


# Parameters, still need to check values to be sure everything matches

params = dict(
    r_growth=1.55,       # 1/h   net in-gut growth rate (doubling ~1 h)
    N_max=3e13,         # cells carrying capacity of the SI lumen
    f_planktonic=0.2,   # fraction of cells free-floating (vs. mucus-attached)
    V_lumen=43,      # mL    SI lumen volume
    Q0=240,           # mL/h  mean luminal flow rate, had to reduce to get feasible graphs
    Q_amp=0.7,          # relative amplitude of flow oscillation (0 = constant)
    Q_period=1.5,       # h     period of flow oscillation
    I_loss=0.0,         # cells/h  immmune clearance rate, set to 0 as immune response is not defined 
    N_intake=1e8,       # cells ingested
    f_acid=0.99,        # fraction killed by stomach acid
    t_delay=20 / 60,    # h     stomach -> SI transit delay (20 min)
)

t_end = 72.0


def Q(t, p):
    Qt = p["Q0"] * (1.0 + p["Q_amp"] * np.sin(2 * np.pi * t / p["Q_period"]))
    return Qt


def immune_loss(t, N, p):
    return min(p["I_loss"], N)


def rhs(t, y, p):
    N = max(y[0], 0.0)
    growth = p["r_growth"] * N * (1.0 - N / p["N_max"])
    washout = p["f_planktonic"] * Q(t, p) * N / p["V_lumen"]
    return [growth - immune_loss(t, N, p) - washout]


def simulate(p=params, t_end=t_end, n_points=1000):
    t_eval = np.linspace(0, t_end, n_points)
    N = np.zeros_like(t_eval)  # no bacteria in the SI before t_delay

    N0 = p["N_intake"] * (1.0 - p["f_acid"])
    mask = t_eval >= p["t_delay"]
    sol = solve_ivp(
        rhs,
        (p["t_delay"], t_end),
        [N0],
        args=(p,),
        t_eval=t_eval[mask],
        method="LSODA",
        rtol=1e-8,
        atol=1e-6,
    )
    N[mask] = np.maximum(sol.y[0], 0.0)
    return t_eval, N


def plot(t, N, p=params):
    fig, axes = plt.subplots(1, 1, figsize=(12, 4.5))

    plt.plot(t, N, lw=2)
    plt.xlabel("Time (h)")
    plt.ylabel("N_SI (cells)")
    plt.title("Cells in Small Intestine Lumen")
    plt.axvline(p["t_delay"], color="gray", ls="--", label="arrival in SI")
    plt.legend()

    fig.suptitle("V. cholerae load in small intestine lumen")
    fig.tight_layout()
    return fig


def parameter_sweep(p=params, param="r_growth", values=(0.3, 0.7, 1.4, 2.0)):
    """Plot N(t) for several values of one parameter."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for v in values:
        q = {**p, param: v}
        t, N = simulate(q)
        ax.semilogy(t, np.maximum(N, 1e-3), lw=2, label=f"{param} = {v}")
    ax.set_xlabel("Time (h)")
    ax.set_ylabel("N_SI (cells)")
    ax.set_ylim(1, p["N_max"] * 10)
    ax.set_title(f"Sensitivity to {param}")
    ax.legend()
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    t, N = simulate()
    
    plot(t, N).savefig("cholera_si_timecourse.png", dpi=150)
    parameter_sweep(param="f_planktonic", values=(0.1, 0.2, 0.5, 0.7, 0.9, 1.0)).savefig(
        "cholera_si_sweep.png", dpi=150)
    plt.show()