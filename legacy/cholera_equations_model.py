"""General equations and helper functions for a cholera biological model.

The equations in this module use the names supplied by the modeler:

- N_SI: non-specific immune population or concentration
- r_growth: growth rate
- N_SI_max: carrying capacity
- N_intake: material entering the lumen
- f_cleared_by_stomach_acid: fraction cleared before entering the lumen
- f_planktonic: fraction used in the planktonic removal term
- Q(t): planktonic material flow rate
- V_L: lumen volume

Time-dependent terms may be supplied as a scalar or a function of time.
"""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Union


ScalarOrTimeFunction = Union[float, Callable[[float], float]]


@dataclass(frozen=True)
class CholeraModelParameters:
    """All constant parameters used by the preliminary NSI equation."""

    r_growth: float
    N_SI_max: float
    f_cleared_by_stomach_acid: float
    N_intake: ScalarOrTimeFunction
    f_planktonic: float
    Q: ScalarOrTimeFunction
    V_L: float

    def __post_init__(self) -> None:
        if self.N_SI_max == 0:
            raise ValueError("N_SI_max must be non-zero.")
        if self.V_L == 0:
            raise ValueError("V_L must be non-zero.")
        if not 0.0 <= self.f_cleared_by_stomach_acid <= 1.0:
            raise ValueError("f_cleared_by_stomach_acid must be between 0 and 1.")
        if not 0.0 <= self.f_planktonic <= 1.0:
            raise ValueError("f_planktonic must be between 0 and 1.")


def value_at_time(value: ScalarOrTimeFunction, time: float) -> float:
    """Return a scalar or evaluate a time-dependent model term."""
    return value(time) if callable(value) else value


def logistic_growth(
    population: float,
    growth_rate: float,
    carrying_capacity: float,
) -> float:
    """Calculate logistic growth: r * N * (1 - N / N_max)."""
    if carrying_capacity == 0:
        raise ValueError("carrying_capacity must be non-zero.")
    return growth_rate * population * (1.0 - population / carrying_capacity)


def planktonic_removal(
    flow_rate: float,
    fraction_planktonic: float,
    population: float,
    lumen_volume: float,
) -> float:
    """Calculate f_planktonic * Q(t) * N_SI / V_L."""
    if lumen_volume == 0:
        raise ValueError("lumen_volume must be non-zero.")
    return fraction_planktonic * flow_rate * population / lumen_volume


def dNSI_dt(
    time: float,
    N_SI: float,
    parameters: CholeraModelParameters,
) -> float:
    """Return dNSI/dt for the preliminary equation.

    dNSI/dt = r_growth * N_SI * (1 - N_SI / N_SI_max)
              + (1 - f_cleared_by_stomach_acid) * N_intake
              - f_planktonic * Q(t) * N_SI / V_L

    N_intake and Q may be constants or functions of time.
    """
    growth = logistic_growth(N_SI, parameters.r_growth, parameters.N_SI_max)
    intake_value = value_at_time(parameters.N_intake, time)
    flow_value = value_at_time(parameters.Q, time)
    intake = (
        1.0 - parameters.f_cleared_by_stomach_acid
    ) * intake_value
    removal = planktonic_removal(
        flow_value,
        parameters.f_planktonic,
        N_SI,
        parameters.V_L,
    )
    return growth + intake - removal


if __name__ == "__main__":
    # Define the model parameters in one place, then pass them to the equation.
    model_parameters = CholeraModelParameters(
        r_growth=0.2,
        N_SI_max=100.0,
        f_cleared_by_stomach_acid=0.0,
        N_intake=10.0,
        f_planktonic=0.1,
        Q=lambda time: 10.0,
        V_L=5.0,
    )

    initial_N_SI = 50.0
    time = 1.0
    print(dNSI_dt(time, initial_N_SI, model_parameters))

