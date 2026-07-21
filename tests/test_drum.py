import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.drum import Drum  # noqa: E402
from thermowave.fluids.ideal_gas import IdealGasFluid  # noqa: E402
from thermowave.fluids.real_fluid import CoolPropFluid  # noqa: E402

WATER = CoolPropFluid(name="Water")


class _FakeState:
    def __init__(self, fluid, node_values, mdots, params):
        self.fluid = fluid
        self._node_values = node_values
        self._mdots = mdots
        self._params = params

    def node(self, name):
        return self._node_values[name]

    def fluid_at(self, name):
        return self.fluid

    def mdot(self, name):
        return self._mdots[name]

    def param(self, name):
        return self._params[name]


def test_rejects_non_positive_volume():
    with pytest.raises(ValueError, match="V must be > 0"):
        Drum(name="d1", V=0.0, P0=1.0e6, fluid=WATER)


def test_rejects_level0_out_of_range():
    with pytest.raises(ValueError, match="level0"):
        Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=1.0)


def test_rejects_non_two_phase_fluid():
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    with pytest.raises(ValueError, match="two-phase-capable"):
        Drum(name="d1", V=2.0, P0=1.0e6, fluid=air)


def test_ports_with_and_without_riser():
    d_riser = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=True)
    assert set(d_riser.ports()) == {"feed_in", "riser_in", "steam_out", "water_out"}
    d_plain = Drum(name="d2", V=2.0, P0=1.0e6, fluid=WATER, has_riser=False)
    assert set(d_plain.ports()) == {"feed_in", "steam_out", "water_out"}


def test_differential_parameters_seeded_from_p0_and_level0():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=0.5)
    dp = d.differential_parameters()
    assert set(dp) == {"P", "h"}
    assert dp["P"] == 1.0e6
    assert math.isclose(dp["h"], d.h0, rel_tol=1e-12)


def test_residuals_zero_at_saturated_outlets_hand_calc():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=False)
    P = 1.0e6
    h_drum = d.h0
    h_g = WATER.saturated_vapor_enthalpy(P)
    h_f = WATER.saturated_liquid_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.feed_in": (P, WATER.enthalpy_pt(P, 430.0)),
            "d1.steam_out": (P, h_g),
            "d1.water_out": (P, h_f),
        },
        mdots={"d1.feed_in": 1.0, "d1.steam_out": 0.2, "d1.water_out": 0.8},
        params={"d1.P": P, "d1.h": h_drum},
    )
    residuals = d.residuals(state)
    assert len(residuals) == 4
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-4)


def test_state_derivative_matches_hand_solved_2x2_system():
    d = Drum(name="d1", V=3.0, P0=1.0e6, fluid=WATER, has_riser=True)
    P, h = 1.0e6, d.h0
    h_feed = WATER.enthalpy_pt(P, 430.0)
    h_riser = WATER.enthalpy_pq(P, 0.15)
    h_g = WATER.saturated_vapor_enthalpy(P)
    h_f = WATER.saturated_liquid_enthalpy(P)
    mdot_feed, mdot_riser, mdot_steam, mdot_water = 1.0, 8.0, 0.9, 8.1

    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.feed_in": (P, h_feed),
            "d1.riser_in": (P, h_riser),
            "d1.steam_out": (P, h_g),
            "d1.water_out": (P, h_f),
        },
        mdots={
            "d1.feed_in": mdot_feed, "d1.riser_in": mdot_riser,
            "d1.steam_out": mdot_steam, "d1.water_out": mdot_water,
        },
        params={"d1.P": P, "d1.h": h},
    )

    # Rebuild the exact 2x2 (mirrors Drum.state_derivative / Tank).
    V = 3.0
    rho = WATER.density_ph(P, h)
    eps_P = max(abs(P) * 1.0e-6, 1.0)
    eps_h = max(abs(h) * 1.0e-6, 1.0)
    drho_dP = (WATER.density_ph(P + eps_P, h) - rho) / eps_P
    drho_dh = (WATER.density_ph(P, h + eps_h) - rho) / eps_h
    mdot_in = mdot_feed + mdot_riser
    energy_in = mdot_feed * (h_feed - h) + mdot_riser * (h_riser - h)
    mdot_out = mdot_steam + mdot_water
    energy_out = mdot_steam * (h_g - h) + mdot_water * (h_f - h)
    a11, a12, b1 = V * drho_dP, V * drho_dh, mdot_in - mdot_out
    a21, a22, b2 = -V, rho * V, energy_in - energy_out
    det = a11 * a22 - a12 * a21
    dPdt = (b1 * a22 - a12 * b2) / det
    dhdt = (a11 * b2 - a21 * b1) / det

    result = d.state_derivative(state)
    assert math.isclose(result["P"], dPdt, rel_tol=1e-9)
    assert math.isclose(result["h"], dhdt, rel_tol=1e-9)


def test_state_derivative_zero_at_balanced_inflow_outflow():
    # Feed at the drum's own enthalpy, steam+water drawn at the mass quality
    # that balances both mass and energy -> both derivatives vanish.
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=False)
    P, h = 1.0e6, d.h0
    h_g = WATER.saturated_vapor_enthalpy(P)
    h_f = WATER.saturated_liquid_enthalpy(P)
    x0 = (h - h_f) / (h_g - h_f)
    mdot_feed = 1.0
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.feed_in": (P, h),  # feed at drum enthalpy
            "d1.steam_out": (P, h_g),
            "d1.water_out": (P, h_f),
        },
        mdots={"d1.feed_in": mdot_feed, "d1.steam_out": x0, "d1.water_out": 1.0 - x0},
        params={"d1.P": P, "d1.h": h},
    )
    deriv = d.state_derivative(state)
    assert math.isclose(deriv["P"], 0.0, abs_tol=1e-6)
    assert math.isclose(deriv["h"], 0.0, abs_tol=1e-6)


def test_report_metrics_level_reflects_liquid_volume_fraction():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=0.6, has_riser=False)
    P, h = 1.0e6, d.h0  # h0 was seeded from level0=0.6
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.feed_in": (P, h), "d1.steam_out": (P, h), "d1.water_out": (P, h),
        },
        mdots={"d1.feed_in": 1.0, "d1.steam_out": 0.0, "d1.water_out": 1.0},
        params={"d1.P": P, "d1.h": h},
    )
    metrics = d.report_metrics(state)
    assert math.isclose(metrics["level [-]"], 0.6, abs_tol=1e-3)
    assert math.isclose(metrics["P [Pa]"], P, rel_tol=1e-12)


def _level_from_ph(fluid, P, h):
    h_f = fluid.saturated_liquid_enthalpy(P)
    h_g = fluid.saturated_vapor_enthalpy(P)
    x = min(max((h - h_f) / (h_g - h_f), 0.0), 1.0)
    v_avg = 1.0 / fluid.density_ph(P, h)
    v_f = 1.0 / fluid.density_ph(P, h_f)
    return (1.0 - x) * v_f / v_avg


def test_net_liquid_inflow_raises_level_and_condenses_pressure():
    # An adiabatic drum's level is a pure integrator (no steady-state
    # restoring force -- real drums need level control, which is why a plain
    # steady solve is singular in h_drum). Its correct *transient* signs:
    # injecting subcooled liquid faster than draw condenses vapor (dP/dt < 0,
    # the classic "shrink" effect) and raises the liquid level (dh/dt < 0).
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=0.5, has_riser=False)
    P, h = 1.0e6, d.h0
    h_g = WATER.saturated_vapor_enthalpy(P)
    h_f = WATER.saturated_liquid_enthalpy(P)
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.feed_in": (P, WATER.enthalpy_pt(P, 430.0)),  # subcooled liquid feed
            "d1.steam_out": (P, h_g),
            "d1.water_out": (P, h_f),
        },
        # feed 2.0 in, draw 0.5 steam + 0.5 water -> net mass accumulation.
        mdots={"d1.feed_in": 2.0, "d1.steam_out": 0.5, "d1.water_out": 0.5},
        params={"d1.P": P, "d1.h": h},
    )
    deriv = d.state_derivative(state)
    assert deriv["h"] < 0.0  # average enthalpy falling => quality down => level rising
    assert deriv["P"] < 0.0  # subcooled feed condenses vapor => pressure shrink


def test_manual_backward_euler_step_raises_level():
    # One backward-Euler step of the drum's own (P, h) under net liquid
    # inflow must increase the reported level -- the transient use case,
    # exercised directly at the component level (a full network needs a level
    # controller to be well-posed; see examples/28).
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=0.5, has_riser=False)
    P, h = 1.0e6, d.h0
    level_before = _level_from_ph(WATER, P, h)

    h_g = WATER.saturated_vapor_enthalpy(P)
    h_f = WATER.saturated_liquid_enthalpy(P)

    def deriv_at(P_cur, h_cur):
        state = _FakeState(
            fluid=WATER,
            node_values={
                "d1.feed_in": (P_cur, WATER.enthalpy_pt(P_cur, 430.0)),
                "d1.steam_out": (P_cur, WATER.saturated_vapor_enthalpy(P_cur)),
                "d1.water_out": (P_cur, WATER.saturated_liquid_enthalpy(P_cur)),
            },
            mdots={"d1.feed_in": 2.0, "d1.steam_out": 0.5, "d1.water_out": 0.5},
            params={"d1.P": P_cur, "d1.h": h_cur},
        )
        return d.state_derivative(state)

    # Explicit Euler a few small steps (enough to move level measurably).
    dt = 0.02
    P_cur, h_cur = P, h
    for _ in range(10):
        r = deriv_at(P_cur, h_cur)
        P_cur += dt * r["P"]
        h_cur += dt * r["h"]
    level_after = _level_from_ph(WATER, P_cur, h_cur)
    assert level_after > level_before
