import math

import pytest

pytest.importorskip("CoolProp")

from thermowave.components.drum import Drum  # noqa: E402
from thermowave.components.sink import Sink  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.exceptions import NetworkTopologyError  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
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


# ---------------------------------------------------------------------------
# P_target / water_out_mdot: extra closing residuals for the free (P_drum,
# riser recirculation mdot) degrees of freedom a typical surrounding
# topology otherwise leaves unclosed (see Drum's own docstring).
# ---------------------------------------------------------------------------


def test_p_target_and_water_out_mdot_add_no_extra_residuals_by_default():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=True)
    P = 1.0e6
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (P, WATER.saturated_vapor_enthalpy(P)),
            "d1.water_out": (P, WATER.saturated_liquid_enthalpy(P)),
        },
        mdots={"d1.water_out": 5.0},
        params={"d1.P": P, "d1.h": d.h0},
    )
    assert len(d.residuals(state)) == 4


def test_p_target_adds_fifth_residual_pinning_drum_pressure():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=False, P_target=1.2e6)
    P_drum = 1.0e6  # deliberately off-target
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (P_drum, WATER.saturated_vapor_enthalpy(P_drum)),
            "d1.water_out": (P_drum, WATER.saturated_liquid_enthalpy(P_drum)),
        },
        mdots={},
        params={"d1.P": P_drum, "d1.h": d.h0},
    )
    residuals = d.residuals(state)
    assert len(residuals) == 5
    assert math.isclose(residuals[-1], P_drum - 1.2e6, rel_tol=1e-12)

    # Zero exactly when P_drum sits at the target.
    state_at_target = _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (1.2e6, WATER.saturated_vapor_enthalpy(1.2e6)),
            "d1.water_out": (1.2e6, WATER.saturated_liquid_enthalpy(1.2e6)),
        },
        mdots={},
        params={"d1.P": 1.2e6, "d1.h": d.h0},
    )
    assert math.isclose(d.residuals(state_at_target)[-1], 0.0, abs_tol=1e-9)


def test_water_out_mdot_adds_residual_pinning_downcomer_flow():
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=True, water_out_mdot=8.0)
    P = 1.0e6
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (P, WATER.saturated_vapor_enthalpy(P)),
            "d1.water_out": (P, WATER.saturated_liquid_enthalpy(P)),
        },
        mdots={"d1.water_out": 6.0},  # deliberately off-target
        params={"d1.P": P, "d1.h": d.h0},
    )
    residuals = d.residuals(state)
    assert len(residuals) == 5
    assert math.isclose(residuals[-1], 6.0 - 8.0, rel_tol=1e-12)


def test_p_target_and_water_out_mdot_together_add_both_residuals():
    d = Drum(
        name="d1", V=2.0, P0=1.0e6, fluid=WATER, has_riser=True,
        P_target=1.1e6, water_out_mdot=8.0,
    )
    P = 1.0e6
    state = _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (P, WATER.saturated_vapor_enthalpy(P)),
            "d1.water_out": (P, WATER.saturated_liquid_enthalpy(P)),
        },
        mdots={"d1.water_out": 8.0},
        params={"d1.P": P, "d1.h": d.h0},
    )
    residuals = d.residuals(state)
    assert len(residuals) == 6
    assert math.isclose(residuals[-2], P - 1.1e6, rel_tol=1e-12)
    assert math.isclose(residuals[-1], 0.0, abs_tol=1e-12)


def test_manual_backward_euler_step_raises_level():
    # One backward-Euler step of the drum's own (P, h) under net liquid
    # inflow must increase the reported level -- the transient use case,
    # exercised directly at the component level (a full network needs a level
    # controller to be well-posed; see examples/28).
    d = Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level0=0.5, has_riser=False)
    P, h = 1.0e6, d.h0
    level_before = _level_from_ph(WATER, P, h)

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


# ---------------------------------------------------------------------------
# level_target: the steady-state closure for the drum's own `h`, and the
# one-sided density partials that keep state_derivative() single-regime near
# the dome boundary. Both were found together while diagnosing why the SEGS
# boiler-internals subsystem stalled -- see Drum's own docstring.
# ---------------------------------------------------------------------------


def test_rejects_level_target_out_of_range():
    with pytest.raises(ValueError, match="level_target must be in"):
        Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level_target=0.0)
    with pytest.raises(ValueError, match="level_target must be in"):
        Drum(name="d1", V=2.0, P0=1.0e6, fluid=WATER, level_target=1.0)


def _drum_state(d, P, h, mdot_water=5.0):
    return _FakeState(
        fluid=WATER,
        node_values={
            "d1.steam_out": (P, WATER.saturated_vapor_enthalpy(P)),
            "d1.water_out": (P, WATER.saturated_liquid_enthalpy(P)),
        },
        mdots={"d1.water_out": mdot_water, "d1.steam_out": 1.0},
        params={"d1.P": P, "d1.h": h},
    )


def test_level_target_adds_residual_zero_at_the_targeted_level():
    P = 1.0e6
    d = Drum(
        name="d1", V=2.0, P0=P, fluid=WATER, level0=0.5, has_riser=False,
        level_target=0.5,
    )
    residuals = d.residuals(_drum_state(d, P, d.h0))
    assert len(residuals) == 5
    # d.h0 was constructed from level0=0.5, so the level residual is zero there.
    assert math.isclose(residuals[-1], 0.0, abs_tol=1e-9)
    # ... and it agrees with what report_metrics reports as the level.
    assert math.isclose(
        d.report_metrics(_drum_state(d, P, d.h0))["level [-]"], 0.5, abs_tol=1e-9
    )


def test_level_target_residual_is_strictly_monotone_in_h_across_the_dome():
    # This is the whole point of the closure: it must give the solver a
    # NON-ZERO d(residual)/d(h_drum) everywhere, including outside the dome.
    # A version that clamped quality to [0, 1] (the way report_metrics does)
    # would be exactly flat for h < h_f / h > h_g -- i.e. the Jacobian column
    # for `h` would go back to zero in precisely the region an overshooting
    # Newton step lands in, which is what made the SEGS boiler subsystem
    # stall in the first place.
    P = 1.0e6
    d = Drum(name="d1", V=2.0, P0=P, fluid=WATER, has_riser=False, level_target=0.5)
    h_f = WATER.saturated_liquid_enthalpy(P)
    h_g = WATER.saturated_vapor_enthalpy(P)
    span = h_g - h_f
    hs = [h_f - 0.05 * span, h_f, h_f + 0.25 * span, h_g, h_g + 0.05 * span]
    values = [d.residuals(_drum_state(d, P, h))[-1] for h in hs]
    assert all(b < a for a, b in zip(values, values[1:])), values


def test_state_derivative_does_not_flip_sign_across_the_saturation_boundary():
    # A finite-difference density partial whose step STRADDLES h_f blends two
    # regimes whose drho/dh differ by ~14x, which flips the SIGN of the 2x2
    # solve's determinant (measured at this pressure, one J/kg below h_f:
    # det = -6.4e-2 straddling vs. +9.6e-2 one-sided) and so of dP/dt and
    # dh/dt themselves -- a discontinuity sitting exactly where a boiler drum
    # lives. Isolated here by giving the drum a pure MASS imbalance with zero
    # energy imbalance (every inlet and outlet carrying the drum's own h), so
    # dP/dt reduces to b1 * rho * V / det and its sign is exactly det's. The
    # magnitude still changes a lot across the boundary -- liquid and
    # two-phase compressibility genuinely differ ~170x -- but never the sign.
    P = 9.997764e6  # the SEGS boiler drum's own pressure, where this was found
    d = Drum(name="d1", V=10.0, P0=P, fluid=WATER, has_riser=False)
    h_f = WATER.saturated_liquid_enthalpy(P)
    signs = set()
    for offset in (-50.0, -5.0, -1.0, -0.1, 0.0, 0.1, 1.0, 5.0, 50.0):
        h = h_f + offset
        state = _FakeState(
            fluid=WATER,
            node_values={
                "d1.feed_in": (P, h),
                "d1.steam_out": (P, h),
                "d1.water_out": (P, h),
            },
            # Net mass accumulation, held identical at every h -> dP/dt > 0.
            mdots={"d1.feed_in": 2.0, "d1.steam_out": 0.5, "d1.water_out": 0.5},
            params={"d1.P": P, "d1.h": h},
        )
        signs.add(math.copysign(1.0, d.state_derivative(state)["P"]))
    assert signs == {1.0}, f"dP/dt changed sign across h_f: {signs}"


# ---------------------------------------------------------------------------
# End-to-end regression: a Drum wired into a real Network.solve() (dt=None,
# the exact mode the singularity this whole file is about only exists in --
# see Drum's own docstring on level_target). Every test above exercises
# residuals()/state_derivative() directly against a _FakeState double, which
# checks the closure's math but never actually runs it through the Newton
# solver -- this is the test that would have caught the original stall.
# ---------------------------------------------------------------------------


def _steady_drum_network(level_target: float | None):
    """A minimal steady network: Source -> Drum -> two Sinks (has_riser=False).

    Feed enters well subcooled (430 K at 1 MPa, ~23 K below saturation), so
    sustaining any steam draw at all needs an external heat input -- there's
    no riser return here to supply it. heat_loss is set negative (i.e. heat
    ADDED, per Drum's own doc: "positive = lost") to exactly the value that
    balances the drum's own mass/energy equations at level_target=0.5 with a
    0.1/1.9 kg/s steam/water split -- chosen by hand-solving the same 2x2
    (b1, b2) system Drum.state_derivative() itself solves, not tuned by
    trial and error against the solver.

    Only P_target and level_target close this drum's own free (P, h) --
    water_out_mdot is deliberately left None: the steam/water split falls
    out of the mass/energy balance on its own once h is closed, and adding
    a third closing residual on top of an already-square 8-unknown system
    would over-determine it (P_drum, h_drum, and the two outlet mdots are
    exactly 4 of those unknowns; P_target + level_target + the two
    automatic steady-state balance equations (state_derivative()==0 for P
    and h) already close all 4).
    """
    P0 = 1.0e6
    net = Network(fluid=WATER)
    src = Source(name="src", P=P0, T=430.0, mdot=2.0)
    drum = Drum(
        name="drum", V=2.0, P0=P0, fluid=WATER, level0=0.5, has_riser=False,
        P_target=P0, level_target=level_target, heat_loss=-402371.42020333465,
    )
    steam_sink = Sink(name="steam_sink")
    water_sink = Sink(name="water_sink")
    for component in (src, drum, steam_sink, water_sink):
        net.add_component(component)
    net.connect(src, "out", drum, "feed_in")
    net.connect(drum, "steam_out", steam_sink, "in")
    net.connect(drum, "water_out", water_sink, "in")
    return net, drum


def test_level_target_closes_a_real_steady_network_end_to_end():
    net, drum = _steady_drum_network(level_target=0.5)
    result = net.solve(progress=False)
    assert result.converged is True

    from thermowave.core.network import NetworkState

    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    metrics = drum.report_metrics(state)
    assert math.isclose(metrics["level [-]"], 0.5, abs_tol=1e-6)


def test_without_level_target_the_same_steady_network_is_not_solvable():
    # The whole point of level_target: omitting it leaves this drum's own
    # differential `h` genuinely unclosed (see Drum's docstring) -- here
    # that shows up as a non-square system (8 unknowns, 7 equations) that
    # Network.solve() rejects outright, rather than a silent numerical
    # stall. Either way, it must NOT converge -- this is what pins the
    # regression test above to level_target specifically, not just "some
    # working Drum network."
    net, _drum = _steady_drum_network(level_target=None)
    with pytest.raises(NetworkTopologyError):
        net.solve(progress=False)
