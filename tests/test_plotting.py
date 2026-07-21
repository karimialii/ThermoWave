import math

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from thermowave.components.compressor import Compressor
from thermowave.components.pipe import Pipe
from thermowave.components.shaft import Shaft
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network, NetworkState
from thermowave.core.plotting import ThermoPlot, plot_series
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.maps.characteristic_map import CharacteristicMap

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)


def _build_turboshaft(N0: float, inertia: float = 0.05):
    src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
    comp = Compressor(name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None)
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-300000.0)
    turb = Turbine(name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None)
    shaft = Shaft(
        name="shaft", components=[comp, turb], signs=[-1.0, 1.0],
        efficiency=0.98, inertia=inertia, dynamic=True, N0=N0,
    )
    snk = Sink(name="snk")

    network = Network(fluid=AIR)
    for component in (src, comp, heater, turb, shaft, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", snk, "in")
    return network, comp, turb, shaft


def test_transient_result_plot_returns_axes_with_expected_data():
    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    history = network.solve_transient(
        duration=0.1, dt=0.05, tol=1e-8, max_iter=400, damping=0.3,
    )

    ax = history.plot((comp, "power [W]"), (turb, "power [W]"), show=False)
    lines = ax.get_lines()
    assert len(lines) == 2
    assert list(lines[0].get_xdata()) == history.times

    expected_turb_power = [
        turb.report_metrics(
            NetworkState(
                fluid=step.fluid, node_P=step.node_P, node_h=step.node_h,
                node_mdot=step.node_mdot, params=step.params,
            )
        )["power [W]"]
        for step in history.steps
    ]
    assert list(lines[1].get_ydata()) == pytest.approx(expected_turb_power)


def test_plot_series_rejects_mismatched_lengths():
    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    with pytest.raises(ValueError, match="must match"):
        plot_series(x=[0.0, 1.0], results=[result], series=[(comp, "power [W]")], show=False)


def test_plot_series_raises_on_component_with_no_report_metrics():
    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    src = next(c for c in network.components if c.name == "src")
    with pytest.raises(ValueError, match="no report_metrics"):
        plot_series(x=[0.0], results=[result], series=[(src, "P [Pa]")], show=False)


def test_plot_series_custom_label_used_on_line():
    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    ax = plot_series(
        x=[0.0], results=[result], series=[(comp, "power [W]", "compressor power")],
        show=False,
    )
    assert ax.get_legend_handles_labels()[1] == ["compressor power"]


def test_plot_series_onto_existing_axes_does_not_create_new_figure():
    import matplotlib.pyplot as plt

    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    result = network.solve(tol=1e-8, max_iter=400, damping=0.3)
    fig, ax = plt.subplots()
    returned = plot_series(
        x=[0.0], results=[result], series=[(comp, "power [W]")], ax=ax, show=False,
    )
    assert returned is ax


def test_thermoplot_line_plots_expected_data():
    ax = ThermoPlot(title="t", xlabel="x", ylabel="y").line(
        [0, 1, 2], [0, 1, 4], label="squares"
    ).finish(show=False)
    lines = ax.get_lines()
    assert len(lines) == 1
    assert list(lines[0].get_xdata()) == [0, 1, 2]
    assert list(lines[0].get_ydata()) == [0, 1, 4]
    assert ax.get_legend_handles_labels()[1] == ["squares"]


def test_thermoplot_scatter_produces_path_collection():
    ax = ThermoPlot().scatter([0, 1, 2], [3, 2, 1], label="pts").finish(show=False)
    assert len(ax.collections) == 1


def test_thermoplot_no_legend_when_no_labels():
    ax = ThermoPlot().line([0, 1], [0, 1]).finish(show=False)
    assert ax.get_legend() is None


def test_thermoplot_twin_axis_shares_x_has_own_y():
    plot = ThermoPlot(xlabel="x", ylabel="power [W]")
    plot.line([0, 1, 2], [10, 20, 30], label="power")
    twin = plot.twin_axis(ylabel="efficiency [-]")
    twin.line([0, 1, 2], [0.5, 0.6, 0.7], label="efficiency")
    ax = plot.finish(show=False)

    assert twin.ax is not ax
    assert twin.ax.get_shared_x_axes().joined(twin.ax, ax)
    assert list(twin.ax.get_lines()[0].get_ydata()) == [0.5, 0.6, 0.7]
    assert twin.ax.get_ylabel() == "efficiency [-]"
    legend_labels = {text.get_text() for text in ax.get_legend().get_texts()}
    assert legend_labels == {"power", "efficiency"}


def test_thermoplot_bar_single_group():
    ax = ThermoPlot().bar(["a", "b", "c"], [1.0, 2.0, 3.0]).finish(show=False)
    heights = sorted(patch.get_height() for patch in ax.patches)
    assert heights == [1.0, 2.0, 3.0]
    assert [t.get_text() for t in ax.get_xticklabels()] == ["a", "b", "c"]


def test_thermoplot_bar_grouped_side_by_side():
    plot = ThermoPlot()
    plot.bar(["a", "b"], [1.0, 2.0], label="before")
    plot.bar(["a", "b"], [1.5, 2.5], label="after")
    ax = plot.finish(show=False)
    assert len(ax.patches) == 4
    assert sorted(round(p.get_height(), 2) for p in ax.patches) == [1.0, 1.5, 2.0, 2.5]
    assert ax.get_legend_handles_labels()[1] == ["before", "after"]


def _simple_map() -> CharacteristicMap:
    return CharacteristicMap.from_file("tests/fixtures/simple_compressor_map.cop")


def test_thermoplot_map_plots_one_line_per_speed():
    char_map = _simple_map()
    n_speeds = len(char_map.speed_lines("pressure_ratio"))
    ax = ThermoPlot().map(char_map, kind="pressure_ratio").finish(show=False)
    assert len(ax.get_lines()) == n_speeds


def test_thermoplot_map_with_operating_points_adds_overlay_line():
    char_map = _simple_map()
    n_speeds = len(char_map.speed_lines("pressure_ratio"))
    ax = ThermoPlot().map(
        char_map, kind="pressure_ratio", operating_points=[(1.0, 2.0), (2.0, 2.5)],
    ).finish(show=False)
    assert len(ax.get_lines()) == n_speeds + 1
    overlay = ax.get_lines()[-1]
    assert list(overlay.get_xdata()) == [1.0, 2.0]
    assert list(overlay.get_ydata()) == [2.0, 2.5]


def test_thermoplot_map_unknown_kind_raises():
    char_map = _simple_map()
    with pytest.raises(ValueError, match="Unknown kind"):
        ThermoPlot().map(char_map, kind="bogus")


def test_thermoplot_transient_matches_plot_series_shape():
    network, comp, turb, shaft = _build_turboshaft(N0=55000.0)
    history = network.solve_transient(
        duration=0.1, dt=0.05, tol=1e-8, max_iter=400, damping=0.3,
    )
    ax = ThermoPlot(ylabel="power [W]").transient(
        history, [(comp, "power [W]"), (turb, "power [W]")]
    ).finish(show=False)
    assert len(ax.get_lines()) == 2
    assert ax.get_xlabel() == "t [s]"
