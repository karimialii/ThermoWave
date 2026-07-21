"""ThermoPlot: the themed, chainable plotting class, showing three chart
types beyond the plain line-vs-time plots seen in earlier examples.

Sweeps a Compressor across a few inlet mass flows (fixed N), collecting one
SolveResult per case, then uses ThermoPlot to build:

  1. A turbomachinery map: the compressor's own T100 Comp.cop pressure-ratio
     map (one line per iso-speed curve) with the swept cases' actual
     (corrected mass flow, PR) operating points overlaid as a traversed
     path — ThermoPlot.map().
  2. A twin-axis plot: power [W] and isentropic efficiency [-] against inlet
     mass flow on independent y-axes, since they're on very different
     scales — ThermoPlot.twin_axis().
  3. A bar chart comparing outlet temperature across the swept cases —
     ThermoPlot.bar().

Every plot_series()/TransientResult.plot() call elsewhere in examples/ and
tutorials/ already gets this same theme for free (they're now thin wrappers
around ThermoPlot) — this example is only about the extra chart types.

Run: .venv/bin/python examples/23_advanced_plotting.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.constants import PA_PER_BAR
from thermowave.core.network import Network, NetworkState
from thermowave.core.plotting import ThermoPlot
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)

N_SHAFT = 65000.0
mdots = [0.50, 0.55, 0.60, 0.63, 0.66]


def metrics_for(comp, result):
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params,
    )
    return comp.report_metrics(state)


results = []
comps = []
for mdot in mdots:
    src = Source(name="src", P=101325.0, T=288.15, mdot=mdot)
    comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=N_SHAFT)
    snk = Sink(name="snk")

    network = Network(fluid=air)
    for component in (src, comp, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    results.append(network.solve(tol=1e-8, max_iter=200, damping=0.5))
    comps.append(comp)

# Corrected mass flow B for each case, from the same non-dimensional group
# Compressor.residuals() itself uses (A = N/sqrt(T_in), B = mdot*sqrt(T_in)/P_in).
operating_points = []
outlet_temps = []
for comp, result, mdot in zip(comps, results, mdots):
    P_in, h_in = result.node_P[comp._inlet_node], result.node_h[comp._inlet_node]
    T_in = air.temperature_ph(P_in, h_in)
    B = mdot * T_in**0.5 / (P_in / PA_PER_BAR)
    operating_points.append((B, metrics_for(comp, result)["PR [-]"]))

    P_out, h_out = result.node_P[comp._outlet_node], result.node_h[comp._outlet_node]
    outlet_temps.append(air.temperature_ph(P_out, h_out))

# --- 1. Turbomachinery map with the traversed operating path ---------------
ThermoPlot(title="Compressor map with swept operating points").map(
    comp.map, kind="pressure_ratio", operating_points=operating_points,
).finish(show=False, save_path="advanced_plotting_map.png")
print("Saved advanced_plotting_map.png")

# --- 2. Twin-axis: power and efficiency against mass flow -------------------
power_plot = ThermoPlot(
    title="Compressor power and efficiency vs. mass flow",
    xlabel="mdot [kg/s]", ylabel="power [W]",
)
power_plot.line(
    mdots, [metrics_for(comp, result)["power [W]"] for comp, result in zip(comps, results)],
    label="power [W]",
)
efficiency_plot = power_plot.twin_axis(ylabel="eta_s [-]")
efficiency_plot.line(
    mdots, [metrics_for(comp, result)["eta_s [-]"] for comp, result in zip(comps, results)],
    label="eta_s [-]", color="#D55E00",
)
power_plot.finish(show=False, save_path="advanced_plotting_twin_axis.png")
print("Saved advanced_plotting_twin_axis.png")

# --- 3. Bar chart: outlet temperature per case ------------------------------
ThermoPlot(
    title="Compressor outlet temperature by case", ylabel="T_out [K]",
).bar([f"{m:.2f} kg/s" for m in mdots], outlet_temps).finish(
    show=False, save_path="advanced_plotting_bar.png"
)
print("Saved advanced_plotting_bar.png")
