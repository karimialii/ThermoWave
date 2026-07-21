import math

import pytest

from thermowave.components.multi_pass_heat_exchanger import MultiPassHeatExchanger
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network, NetworkState
from thermowave.fluids.ideal_gas import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def test_multi_pass_hx_rejects_non_positive_ua():
    with pytest.raises(ValueError, match="UA must be > 0"):
        MultiPassHeatExchanger(name="hx1", UA=0.0, PR_hot=0.98, PR_cold=0.97)


def test_multi_pass_hx_rejects_bad_n_passes():
    with pytest.raises(ValueError, match="n_passes"):
        MultiPassHeatExchanger(name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97, n_passes=0)


def test_multi_pass_hx_rejects_unknown_arrangement():
    with pytest.raises(ValueError, match="arrangement"):
        MultiPassHeatExchanger(
            name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97, arrangement="radial"
        )


def test_multi_pass_hx_ports_returns_four_named_ports():
    hx = MultiPassHeatExchanger(name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97)
    assert hx.ports() == {
        "hot_in": "hx1.hot_in", "hot_out": "hx1.hot_out",
        "cold_in": "hx1.cold_in", "cold_out": "hx1.cold_out",
    }


def test_multi_pass_hx_internal_nodes_scale_with_n_passes():
    hx1 = MultiPassHeatExchanger(name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97, n_passes=1)
    hx3 = MultiPassHeatExchanger(name="hx3", UA=500.0, PR_hot=0.98, PR_cold=0.97, n_passes=3)
    assert hx1.internal_nodes() == []
    assert len(hx3.internal_nodes()) == 4  # 2 mid-nodes per stream, 2 streams


def _run(n_passes, arrangement="counterflow", UA=800.0):
    hot_src = Source(name="hsrc", P=110000.0, T=800.0, mdot=0.4)
    cold_src = Source(name="csrc", P=300000.0, T=400.0, mdot=0.4)
    hx = MultiPassHeatExchanger(
        name="hx", UA=UA, PR_hot=0.98, PR_cold=0.97, n_passes=n_passes, arrangement=arrangement,
    )
    hot_snk = Sink(name="hsnk")
    cold_snk = Sink(name="csnk")

    network = Network(fluid=AIR)
    for component in (hot_src, cold_src, hx, hot_snk, cold_snk):
        network.add_component(component)
    network.connect(hot_src, "out", hx, "hot_in")
    network.connect(hx, "hot_out", hot_snk, "in")
    network.connect(cold_src, "out", hx, "cold_in")
    network.connect(hx, "cold_out", cold_snk, "in")

    result = network.solve(tol=1e-10, max_iter=300, verbose=False)
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    return result, hx.report_metrics(state)


def test_multi_pass_hx_single_pass_converges_and_transfers_heat_hot_to_cold():
    result, metrics = _run(n_passes=1)
    assert result.converged
    assert metrics["power [W]"] > 0.0
    assert 0.0 < metrics["effectiveness [-]"] < 1.0


def test_multi_pass_hx_counterflow_effectiveness_is_invariant_to_n_passes():
    # Uniform mdot/cp/UA along the exchanger: subdividing an exact
    # counterflow closed-form solution along a fixed flow direction doesn't
    # change it (see the component's own docstring) -- this is the intended,
    # verified behavior, not a bug.
    _, metrics_1 = _run(n_passes=1, arrangement="counterflow")
    _, metrics_4 = _run(n_passes=4, arrangement="counterflow")
    assert math.isclose(
        metrics_1["effectiveness [-]"], metrics_4["effectiveness [-]"], rel_tol=1e-6
    )


def test_multi_pass_hx_parallel_effectiveness_is_invariant_to_n_passes():
    _, metrics_1 = _run(n_passes=1, arrangement="parallel")
    _, metrics_4 = _run(n_passes=4, arrangement="parallel")
    assert math.isclose(
        metrics_1["effectiveness [-]"], metrics_4["effectiveness [-]"], rel_tol=1e-6
    )


def test_multi_pass_hx_counterflow_beats_parallel_at_same_ua():
    _, metrics_counter = _run(n_passes=1, arrangement="counterflow")
    _, metrics_parallel = _run(n_passes=1, arrangement="parallel")
    assert metrics_counter["effectiveness [-]"] > metrics_parallel["effectiveness [-]"]


def test_multi_pass_hx_counterflow_effectiveness_matches_balanced_closed_form():
    # Equal mdot on both sides + same fluid -> C_hot == C_cold (Cr == 1),
    # where counterflow's effectiveness has the simple closed form
    # eff = NTU/(1+NTU) -- an independent hand-calc cross-check of the
    # general (1-exp(...))/(1-Cr*exp(...)) formula's Cr==1 branch.
    UA = 800.0
    _, metrics = _run(n_passes=1, arrangement="counterflow", UA=UA)
    C_min = 0.4 * AIR.cp(110000.0, 800.0)  # both streams: mdot=0.4, same air
    NTU = UA / C_min
    expected_eff = NTU / (1.0 + NTU)
    # rel_tol loosened slightly to absorb _smooth_min's small (eps/2) bias
    # right at C_hot == C_cold -- see its own docstring for why that
    # smoothing exists at all.
    assert math.isclose(metrics["effectiveness [-]"], expected_eff, rel_tol=1e-4)


# ---------------------------------------------------------------------------
# "shell_and_tube": genuine reversing-header multi-pass F-correction factor
# ---------------------------------------------------------------------------


def test_effectiveness_one_shell_pass_matches_incropera_worked_example():
    # Classic textbook check (Incropera): NTU=1.5, Cr=0.5 -> eff ~ 0.64 for
    # a single 1-2 shell-and-tube pass (Bowman/Mueller/Nagle relation).
    eff = MultiPassHeatExchanger._effectiveness_one_shell_pass(1.5, 0.5)
    assert math.isclose(eff, 0.6385, abs_tol=1e-3)


def test_effectiveness_shell_and_tube_reduces_to_one_shell_pass_at_n_equals_one():
    eff1 = MultiPassHeatExchanger._effectiveness_one_shell_pass(1.5, 0.5)
    eff_n1 = MultiPassHeatExchanger._effectiveness_shell_and_tube(1.5, 0.5, 1)
    assert math.isclose(eff1, eff_n1, rel_tol=1e-12)


def test_effectiveness_shell_and_tube_increases_toward_counterflow_limit_with_n():
    def counterflow(NTU, Cr):
        exp_term = math.exp(-NTU * (1.0 - Cr))
        return (1.0 - exp_term) / (1.0 - Cr * exp_term)

    NTU, Cr = 1.5, 0.5
    cf_limit = counterflow(NTU, Cr)
    effs = [
        MultiPassHeatExchanger._effectiveness_shell_and_tube(NTU, Cr, N)
        for N in (1, 2, 4, 8, 16, 32)
    ]
    # Strictly increasing, and every value stays below the true counterflow
    # limit it's approaching -- unlike counterflow/parallel (mathematically
    # invariant to n_passes), more shell passes genuinely helps here.
    assert all(a < b for a, b in zip(effs, effs[1:]))
    assert all(e < cf_limit for e in effs)
    assert math.isclose(effs[-1], cf_limit, rel_tol=1e-3)


def test_effectiveness_shell_and_tube_matches_incropera_balanced_closed_form():
    # Cr == 1 special-cased branch: eff = N*eff_1 / (1 + (N-1)*eff_1), with
    # eff_1 evaluated at NTU_1 = NTU/N (one shell's own share of the total
    # NTU), same as the function under test.
    NTU, N = 1.5, 3
    eff1 = MultiPassHeatExchanger._effectiveness_one_shell_pass(NTU / N, 1.0)
    expected = N * eff1 / (1.0 + (N - 1.0) * eff1)
    actual = MultiPassHeatExchanger._effectiveness_shell_and_tube(NTU, 1.0, N)
    assert math.isclose(actual, expected, rel_tol=1e-12)


def test_multi_pass_hx_shell_and_tube_no_internal_nodes():
    hx1 = MultiPassHeatExchanger(
        name="hx1", UA=500.0, PR_hot=0.98, PR_cold=0.97, n_passes=1, arrangement="shell_and_tube"
    )
    hx4 = MultiPassHeatExchanger(
        name="hx4", UA=500.0, PR_hot=0.98, PR_cold=0.97, n_passes=4, arrangement="shell_and_tube"
    )
    assert hx1.internal_nodes() == []
    assert hx4.internal_nodes() == []


def test_multi_pass_hx_shell_and_tube_converges_and_transfers_heat():
    result, metrics = _run(n_passes=2, arrangement="shell_and_tube")
    assert result.converged
    assert metrics["power [W]"] > 0.0
    assert 0.0 < metrics["effectiveness [-]"] < 1.0


def test_multi_pass_hx_shell_and_tube_effectiveness_genuinely_improves_with_n_passes():
    # Unlike counterflow/parallel (n_passes is mathematically a no-op) and
    # crossflow (not guaranteed to improve), shell_and_tube's n_passes is a
    # real physical multi-pass lever -- more shell passes should measurably
    # raise effectiveness at the same total UA.
    _, metrics_1 = _run(n_passes=1, arrangement="shell_and_tube")
    _, metrics_4 = _run(n_passes=4, arrangement="shell_and_tube")
    assert metrics_4["effectiveness [-]"] > metrics_1["effectiveness [-]"]


def test_multi_pass_hx_shell_and_tube_at_n1_lies_between_parallel_and_counterflow():
    _, metrics_counter = _run(n_passes=1, arrangement="counterflow")
    _, metrics_parallel = _run(n_passes=1, arrangement="parallel")
    _, metrics_shell = _run(n_passes=1, arrangement="shell_and_tube")
    assert metrics_parallel["effectiveness [-]"] < metrics_shell["effectiveness [-]"]
    assert metrics_shell["effectiveness [-]"] < metrics_counter["effectiveness [-]"]
