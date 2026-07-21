import math

from thermowave.components.junction import Junction
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


class _FakeState:
    def __init__(self, fluid, node_values, mdot_values):
        self.fluid = fluid
        self._node_values = node_values
        self._mdot_values = mdot_values

    def node(self, name: str) -> tuple[float, float]:
        return self._node_values[name]

    def mdot(self, name: str) -> float:
        return self._mdot_values[name]


def test_junction_ports_splitter():
    j = Junction(name="j1", n_inlets=1, n_outlets=2)
    assert j.ports() == {"in0": "j1.in0", "out0": "j1.out0", "out1": "j1.out1"}


def test_junction_ports_merger():
    j = Junction(name="j1", n_inlets=2, n_outlets=1)
    assert j.ports() == {"in0": "j1.in0", "in1": "j1.in1", "out0": "j1.out0"}


def test_junction_rejects_mismatched_split_fractions_length():
    try:
        Junction(name="j1", n_inlets=1, n_outlets=2, split_fractions=[1.0])
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_junction_splitter_residuals_matches_hand_calc():
    P_in, h_in = 300000.0, 400000.0
    mdot_in = 2.0
    fractions = [0.3, 0.7]
    j = Junction(name="j1", n_inlets=1, n_outlets=2, split_fractions=fractions)

    # Correct guesses: uniform pressure, mixed (== inlet, single stream) enthalpy,
    # mdot split by the given fractions.
    state = _FakeState(
        fluid=AIR,
        node_values={
            "j1.in0": (P_in, h_in),
            "j1.out0": (P_in, h_in),
            "j1.out1": (P_in, h_in),
        },
        mdot_values={
            "j1.in0": mdot_in,
            "j1.out0": fractions[0] * mdot_in,
            "j1.out1": fractions[1] * mdot_in,
        },
    )
    residuals = j.residuals(state)
    assert len(residuals) == 6  # 3 per outlet (P, h, mdot) * 2 outlets
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-9)


def test_junction_merger_mixes_enthalpy_by_mass_weighted_average():
    P = 300000.0
    h1, h2 = 400000.0, 600000.0
    mdot1, mdot2 = 1.0, 3.0
    j = Junction(name="j1", n_inlets=2, n_outlets=1)

    expected_h_mix = (mdot1 * h1 + mdot2 * h2) / (mdot1 + mdot2)
    state = _FakeState(
        fluid=AIR,
        node_values={
            "j1.in0": (P, h1),
            "j1.in1": (P, h2),
            "j1.out0": (P, expected_h_mix),
        },
        mdot_values={
            "j1.in0": mdot1,
            "j1.in1": mdot2,
            "j1.out0": mdot1 + mdot2,
        },
    )
    residuals = j.residuals(state)
    assert len(residuals) == 3
    for r in residuals:
        assert math.isclose(r, 0.0, abs_tol=1e-6)


def test_junction_merger_mass_residual_reflects_imbalance():
    P, h = 300000.0, 400000.0
    j = Junction(name="j1", n_inlets=2, n_outlets=1)

    state = _FakeState(
        fluid=AIR,
        node_values={
            "j1.in0": (P, h),
            "j1.in1": (P, h),
            "j1.out0": (P, h),
        },
        mdot_values={"j1.in0": 1.0, "j1.in1": 1.0, "j1.out0": 2.5},
    )
    _p_res, _h_res, mass_residual = j.residuals(state)
    assert math.isclose(mass_residual, 0.5, abs_tol=1e-9)


def test_junction_residuals_does_not_raise_on_all_zero_inflow():
    # A plausible mid-Newton-iteration (or transient t=0) guess: every
    # inlet's mdot is exactly 0, which would otherwise divide-by-zero
    # computing h_mix.
    P, h = 300000.0, 400000.0
    j = Junction(name="j1", n_inlets=2, n_outlets=1)

    state = _FakeState(
        fluid=AIR,
        node_values={
            "j1.in0": (P, h),
            "j1.in1": (P, h),
            "j1.out0": (P, h),
        },
        mdot_values={"j1.in0": 0.0, "j1.in1": 0.0, "j1.out0": 0.0},
    )
    residuals = j.residuals(state)
    assert all(math.isfinite(r) for r in residuals)


def test_junction_warm_start_pairs_pairs_first_inlet_with_every_outlet():
    j = Junction(name="j1", n_inlets=2, n_outlets=3)
    assert j.warm_start_pairs() == [("in0", "out0"), ("in0", "out1"), ("in0", "out2")]
