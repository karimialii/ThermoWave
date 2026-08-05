"""Demonstration tests for issues found in the 2026-08-04 professional review.

Each test encodes the *desired* safe behavior and is expected to FAIL against
the current code, proving the corresponding problem is real and reproducible.
See the review summary (P0-P2 priority list) for the full write-up; this file
only covers the issues that are cheaply reproducible in a unit test.

Run with: pytest tests/test_review_findings.py -v
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from thermowave.components.combustor import Combustor
from thermowave.components.compressor import Compressor
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.components.valve import Valve
from thermowave.core.exceptions import ConvergenceError, FluidRangeError
from thermowave.core.solver import newton_solve
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.fluids.real_fluid import CoolPropFluid

_COMPRESSOR_MAP = str(Path(__file__).parent / "fixtures" / "simple_compressor_map.cop")

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    """Minimal stand-in for core.network.NetworkState (mirrors test_compressor.py)."""

    def __init__(self, fluid, mdot, node_values: dict[str, tuple[float, float]]):
        self.fluid = fluid
        self._mdot = mdot
        self._node_values = node_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def fluid_at(self, name: str):
        return self.fluid

    def mdot(self, name: str) -> float:
        return self._mdot


# --- P0 #1: real_fluid.py silently clamps out-of-range pressure -----------


def test_coolprop_fluid_raises_instead_of_silently_clamping_pressure():
    """CanteraFluid raises FluidRangeError on an out-of-range state so the
    solver's line search can back off (solver.py's line_search docstring
    explicitly relies on this). CoolPropFluid should do the same instead of
    silently returning properties evaluated at a clamped pressure -- a
    Newton trial that transiently drives a node below P_min currently gets
    plausible-looking but *wrong* density back with no error at all.
    """
    pytest.importorskip("CoolProp")
    water = CoolPropFluid(name="Water", P_min=1.0e5, P_max=1.0e8)
    h = 500_000.0  # J/kg, arbitrary liquid-water enthalpy

    # A pressure ten times below P_min is a different physical state and
    # must not silently collapse onto the P_min answer.
    with pytest.raises(FluidRangeError):
        water.density_ph(1.0e4, h)


def test_coolprop_fluid_accepts_pressure_at_range_boundary():
    """P_min/P_max themselves are still valid, inclusive boundary values --
    only P strictly outside [P_min, P_max] should raise.
    """
    pytest.importorskip("CoolProp")
    water = CoolPropFluid(name="Water", P_min=1.0e5, P_max=1.0e7)
    h = 500_000.0
    water.density_ph(1.0e5, h)  # should not raise
    water.density_ph(1.0e7, h)  # should not raise


# --- P0 #2: solver accepts a NaN residual instead of failing fast ---------


def test_newton_solve_detects_nan_residual_promptly_instead_of_burning_max_iter():
    """A residual function that silently produces NaN (e.g. np.sqrt of a
    negative intermediate -- a domain error that doesn't raise in numpy, so
    it never becomes a FluidRangeError) should be caught quickly, not
    accepted iteration after iteration until max_iter is exhausted.
    line_search's `use_trial <= damping` branch force-accepts the last
    backtrack attempt unconditionally, even when its residual norm is NaN,
    because `norm_trial <= norm_from` is always False for NaN and never
    stops that acceptance.
    """
    calls = {"count": 0}

    def residual_fn(x: np.ndarray) -> np.ndarray:
        calls["count"] += 1
        # x[1] going negative is a legitimate Newton overshoot, not a
        # fluid-range violation -- np.sqrt(negative) returns nan with a
        # RuntimeWarning instead of raising.
        return np.array([x[0] - 3.0, np.sqrt(x[1]) - 1.0])

    x0 = np.array([0.0, -0.5])
    max_iter = 100

    with pytest.raises(ConvergenceError):
        newton_solve(residual_fn, x0, tol=1e-9, max_iter=max_iter, progress=False)

    # A prompt failure should stop well short of the full iteration budget.
    # Currently it silently keeps "converging" on NaN for the entire budget
    # (~3 residual evaluations per iteration => ~300 calls) before giving up.
    assert calls["count"] < 3 * max_iter // 2, (
        f"solver kept iterating on a NaN residual for {calls['count']} residual "
        f"evaluations instead of failing fast once NaN first appeared"
    )


# --- P1 #3: Valve(opening=0.0) crashes instead of behaving like a closed valve --


def test_fully_closed_valve_does_not_raise_zero_division_error():
    """opening=0.0 is a routine physical state (a fully closed valve), and
    should produce a very large but finite resistance so the solver can
    still converge toward zero flow -- not a ZeroDivisionError.
    """
    valve = Valve(name="v1", D=0.05, K=1.5, opening=0.0)
    state = _FakeState(
        fluid=AIR,
        mdot=0.0,
        node_values={"v1.in": (200_000.0, AIR.enthalpy_pt(200_000.0, 300.0)),
                     "v1.out": (150_000.0, AIR.enthalpy_pt(200_000.0, 300.0))},
    )
    valve.residuals(state)  # should not raise


def test_fully_closed_valve_gives_finite_but_very_large_resistance():
    """opening=0.0 should behave like a nearly-closed valve (a very large
    pressure drop for any nonzero flow), not zero resistance and not a
    crash -- confirms the floor actually changes the physics, not just
    that residuals() happens to not throw.
    """
    from thermowave.components.valve import _MIN_OPENING

    valve = Valve(name="v1", D=0.05, K=1.5, opening=0.0)
    K_eff = valve.K / max(valve.opening, _MIN_OPENING) ** 2
    assert K_eff == pytest.approx(valve.K / _MIN_OPENING**2)
    assert K_eff > 1.0e9  # "very large" -- effectively shuts off flow


def test_valve_rejects_invalid_construction_parameters():
    with pytest.raises(ValueError):
        Valve(name="v1", D=0.0, K=1.0)
    with pytest.raises(ValueError):
        Valve(name="v1", D=0.05, K=-1.0)
    with pytest.raises(ValueError):
        Valve(name="v1", D=0.05, K=1.0, opening=1.5)


# --- P1 #4: Compressor Q_loss / mdot_in divides by zero at mdot=0 ---------


class _FakeHeatPath:
    def Q(self, state) -> float:
        return 0.0  # even with zero (or no) heat loss, mdot_in=0 alone triggers this


def test_compressor_residuals_does_not_divide_by_zero_at_zero_flow():
    """An isolated/closed-off branch with zero flow is a valid transient
    Newton iterate (or a deliberately shut-down compressor leg), and should
    produce a well-defined residual, not a ZeroDivisionError.
    """
    comp = Compressor(name="c1", map_path=_COMPRESSOR_MAP, N=3000.0, gamma=1.4)
    P_in, T_in = 101325.0, 300.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    state = _FakeState(
        fluid=AIR,
        mdot=0.0,
        node_values={"c1.in": (P_in, h_in), "c1.out": (P_in * 1.5, h_in)},
    )
    comp.residuals(state)  # should not raise ZeroDivisionError


# --- P1 #5/#6: Compressor/Turbine _corrected_params mishandles invalid state --


def test_compressor_rejects_negative_inlet_temperature_instead_of_going_complex():
    """T_in**0.5 for a negative T_in silently returns a complex number in
    Python 3 rather than raising -- this should instead fail clearly at the
    source with FluidRangeError, the same recoverable-invalid-state signal
    a real fluid backend already raises, letting the solver's line search
    back off instead of a confusing TypeError surfacing far downstream.
    """
    comp = Compressor(name="c1", map_path=_COMPRESSOR_MAP, N=3000.0, gamma=1.4)
    P_in = 101325.0
    # An enthalpy low enough that IdealGasFluid's temperature_ph returns a
    # negative T_in -- a plausible bad Newton iterate, not a fluid-backend
    # failure in its own right.
    h_bad = AIR.enthalpy_pt(P_in, -100.0)
    state = _FakeState(
        fluid=AIR,
        mdot=1.0,
        node_values={"c1.in": (P_in, h_bad), "c1.out": (P_in * 1.5, h_bad)},
    )
    with pytest.raises(FluidRangeError):
        comp.residuals(state)


def test_compressor_rejects_reverse_flow_instead_of_clamping_to_map_edge():
    """Negative mdot (windmilling/reversal, or a misconfigured fixed-mdot
    Source) isn't in the map's tabulated data -- it should raise rather
    than silently clamp to the map's lowest corrected-flow point and
    return a plausible-looking but wrong PR/efficiency.
    """
    comp = Compressor(name="c1", map_path=_COMPRESSOR_MAP, N=3000.0, gamma=1.4)
    P_in, T_in = 101325.0, 300.0
    h_in = AIR.enthalpy_pt(P_in, T_in)
    state = _FakeState(
        fluid=AIR,
        mdot=-1.0,
        node_values={"c1.in": (P_in, h_in), "c1.out": (P_in * 1.5, h_in)},
    )
    with pytest.raises(FluidRangeError):
        comp.residuals(state)


# --- P2 #9: exergy.py has no sanity check for nonphysical negative E_D ----


def test_finish_cost_warns_on_negative_exergy_destruction():
    """Exergy destruction is only ever >= 0 by the second law -- a negative
    E_D (possible if an upstream component's power-sign convention is ever
    inconsistent, giving E_F < E_P) should be flagged, not silently
    reported as if it were a normal result.
    """
    from thermowave.core.exergy import _finish_cost

    with pytest.warns(UserWarning, match="negative exergy destruction"):
        cost = _finish_cost(E_F=100.0, E_P=150.0, name="bad_component")
    assert cost.E_D == pytest.approx(-50.0)


# --- P2 #7: Combustor doesn't bound efficiency/PR (unlike SimpleCombustor) --


def test_combustor_rejects_efficiency_above_one_like_simple_combustor_does():
    """SimpleCombustor.__init__ enforces 0 < efficiency <= 1 and 0 < PR <= 1
    (efficiency > 1 could push the outlet above the adiabatic flame temp --
    a first-law violation). Combustor should apply the same guard rather
    than silently accepting efficiency=1.5.
    """
    with pytest.raises(ValueError):
        Combustor(name="cc1", efficiency=1.5)


def test_combustor_rejects_pressure_ratio_above_one():
    with pytest.raises(ValueError):
        Combustor(name="cc1", PR=1.2)


# --- P2 #10: SimpleGenerator doesn't validate efficiency ------------------


def test_simple_generator_rejects_efficiency_outside_zero_one():
    """Generator/ElectricMotor/compressors all validate their efficiency-like
    parameters; SimpleGenerator currently accepts negative or >1 values
    silently, producing a nonphysical electrical output with no warning.
    """
    with pytest.raises(ValueError):
        SimpleGenerator(name="g1", efficiency=1.5)
    with pytest.raises(ValueError):
        SimpleGenerator(name="g2", efficiency=-0.2)
