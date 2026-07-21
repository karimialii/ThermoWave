"""Combined-cycle showcase: a T100-class microturbine topped onto a steam
bottoming cycle — steady state, load sweep, emissions, spool-lag transient,
and a results-annotated cycle diagram, all from one script.

This is the package's "everything at once" example. One run produces:

  1. Full-load steady state of the gas topping cycle (real T100 compressor/
     turbine maps, Cantera equilibrium combustor, mass flow solved from the
     maps rather than dictated), plus a steam bottoming cycle (CoolProp
     water, Pump -> SimpleEvaporator -> SteamTurbine -> SimpleCondenser)
     sized from the gas exhaust through an HRSG energy balance with an
     explicit pinch-point check.
  2. A part-load sweep (fuel flow swept at grid-synchronous shaft speed,
     warm-starting each solve from the previous one) -> compressor map with
     the traversed operating path, power/efficiency curves, and equilibrium
     NO/CO2 emissions vs. load.
  3. A fuel-step transient on a dynamic shaft (rotor inertia) -> the
     classic compressor/spool lag: fuel jumps instantly, shaft speed,
     pressure ratio, and power follow with a visible time constant.
  4. A schematic drawing of the whole plant, annotated with the actual
     solved full-load numbers (station temperatures/pressures/flows and
     the power split).

Sizing is based on the published Turbec T100 microturbine, whose actual
compressor/turbine characteristic maps ship with this package ("T100
Comp.cop"/"T100 Turb.tur"): ~100 kWe class, pressure ratio ~4.3-4.5,
nominal speed 70,000 rev/min, turbine inlet ~950 degC, turbine outlet
limited to ~645 degC (see e.g. Zornek et al., Applied Energy 2015,
https://elib.dlr.de/99014/ and the Ensola T100 technical description).
The production T100 reaches ~30% electrical efficiency with an internal
recuperator; here the recuperator is deliberately omitted so the exhaust
leaves hot (~590 degC) and the heat goes to a heat-recovery steam
generator (HRSG) instead — the combined-cycle arrangement — recovering a
comparable overall efficiency with extra power instead of hot water.

The gas and steam sides are two Network objects (a network has one
boundary fluid — air here, water there) coupled the way sequential
combined-cycle solvers do it: the solved gas-side exhaust state feeds an
HRSG energy balance that sizes the steam mass flow, then the steam network
solves with that flow. The pinch check afterwards confirms the coupling is
thermodynamically feasible, not just energy-consistent.

Requires the cantera, coolprop, and plot extras:
    pip install thermowave[cantera,coolprop,plot]

Runtime is a few minutes: every Newton iteration of every solve runs a
Cantera equilibrium calculation inside the combustor.

Run: .venv/bin/python examples/25_combined_cycle_showcase.py
"""

from thermowave.components.combustor import Combustor
from thermowave.components.compressor import Compressor
from thermowave.components.pump import Pump
from thermowave.components.schedule import Schedule
from thermowave.components.shaft import Shaft
from thermowave.components.simple_condenser import SimpleCondenser
from thermowave.components.simple_evaporator import SimpleEvaporator
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.steam_turbine import SteamTurbine
from thermowave.components.turbine import Turbine
from thermowave.core.constants import PA_PER_BAR
from thermowave.core.network import Network, NetworkState
from thermowave.core.plotting import ThermoPlot
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.fluids.real_fluid import CoolPropFluid

# --- Shared plant definition -------------------------------------------------

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
water = CoolPropFluid(name="Water")
GAMMA = 1005.0 / (1005.0 - 287.05)

P_AMBIENT = 101325.0   # Pa
T_AMBIENT = 288.15     # K, ISO conditions
N_GRID = 65000.0       # rev/min, grid-synchronous shaft speed for the sweep
LHV_CH4 = 50.0e6       # J/kg, methane lower heating value (for efficiency only —
#                        the combustor itself uses Cantera equilibrium, not LHV)
ETA_GEN_GT = 0.95      # gas-turbine generator (incl. power electronics)
ETA_GEN_ST = 0.96      # steam-turbine generator

FUEL_FULL_LOAD = 0.012          # kg/s CH4 — full load at N_GRID
FUEL_SWEEP = [0.007, 0.008, 0.009, 0.010, 0.011, 0.012]

# Steam bottoming cycle sizing
P_BOILER = 15.0e5      # Pa
P_COND = 1.0e4         # Pa (0.1 bar, ~46 degC condensing)
SUPERHEAT = 250.0      # K above saturation -> ~448 degC live steam
T_STACK = 420.0        # K HRSG stack temperature (well above water dewpoint)
PINCH_MIN = 10.0       # K minimum allowed evaporator pinch


def build_gas_turbine(mdot_fuel, N):
    """Gas topping cycle. Source(mdot=None) + Sink(P=ambient) leave total
    mass flow a Newton unknown: the flow that actually passes is whatever
    the compressor/turbine maps imply at shaft speed N — the physical
    closure a real machine has, instead of dictating the flow up front."""
    src = Source(name="src", P=P_AMBIENT, T=T_AMBIENT, mdot=None, mdot_guess=0.63)
    comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=GAMMA, N=N)
    comb = Combustor(name="comb", PR=0.96, mdot_fuel=mdot_fuel, fuel="CH4")
    turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=GAMMA, N=N)
    snk = Sink(name="snk", P=P_AMBIENT)

    network = Network(fluid=air)
    for component in (src, comp, comb, turb, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", comb, "in")
    network.connect(comb, "out", turb, "in")
    network.connect(turb, "out", snk, "in")
    return network, comp, comb, turb


def metrics_for(component, result):
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    return component.report_metrics(state)


def gt_operating_point(result, comp, comb, turb):
    """Everything downstream needs from one solved gas-turbine state."""
    cm, bm, tm = (metrics_for(c, result) for c in (comp, comb, turb))
    mdot_air = result.node_mdot["comp.in"]
    mdot_exhaust = result.node_mdot["turb.out"]
    T_exhaust = air.temperature_ph(result.node_P["turb.out"], result.node_h["turb.out"])
    W_gt_shaft = tm["power [W]"] - cm["power [W]"]
    return {
        "mdot_air": mdot_air,
        "mdot_fuel": bm["mdot_fuel [kg/s]"],
        "mdot_exhaust": mdot_exhaust,
        "PR": cm["PR [-]"],
        "TIT": bm["T_out [K]"],
        "T_exhaust": T_exhaust,
        "W_comp": cm["power [W]"],
        "W_turb": tm["power [W]"],
        "W_gt_shaft": W_gt_shaft,
        "W_gt_elec": W_gt_shaft * ETA_GEN_GT,
        "fuel_power": bm["mdot_fuel [kg/s]"] * LHV_CH4,
        "X_NO": bm.get("X_NO [-]", 0.0),
        "X_CO": bm.get("X_CO [-]", 0.0),
        "X_CO2": bm.get("X_CO2 [-]", 0.0),
    }


def hrsg_steam_mdot(gt):
    """HRSG energy balance: heat recovered cooling the exhaust from turbine
    outlet to the stack temperature sizes the steam flow the boiler can
    sustain. This is the script-level coupling between the two networks —
    the same sequential (gas solve -> heat duty -> steam solve) scheme
    dedicated combined-cycle programs use."""
    Q_hrsg = gt["mdot_exhaust"] * (
        air.enthalpy_pt(P_AMBIENT, gt["T_exhaust"]) - air.enthalpy_pt(P_AMBIENT, T_STACK)
    )
    # Feedwater -> live steam specific enthalpy rise, evaluated at the same
    # states the steam network below will solve to.
    T_cond = water.temperature_ph(P_COND, water.enthalpy_pq(P_COND, 0.0))
    h_feed = water.enthalpy_pq(P_COND, 0.0)
    h_feed_pumped = h_feed + (P_BOILER - P_COND) / water.density_ph(P_COND, h_feed) / 0.75
    T_sat = water.saturation_temperature(P_BOILER)
    h_steam = water.enthalpy_pt(P_BOILER, T_sat + SUPERHEAT)
    mdot_steam = Q_hrsg / (h_steam - h_feed_pumped)

    # Evaporator pinch check: the gas temperature where the water side is at
    # saturated liquid (economizer done, boiling about to start) must still
    # sit PINCH_MIN above the saturation temperature — an energy balance
    # alone can't tell you the heat can actually flow.
    h_sat_liquid = water.enthalpy_pq(P_BOILER, 0.0)
    Q_economizer = mdot_steam * (h_sat_liquid - h_feed_pumped)
    h_gas_pinch = air.enthalpy_pt(P_AMBIENT, T_STACK) + Q_economizer / gt["mdot_exhaust"]
    T_gas_pinch = air.temperature_ph(P_AMBIENT, h_gas_pinch)
    pinch = T_gas_pinch - T_sat
    return Q_hrsg, mdot_steam, pinch, T_cond


def build_steam_cycle(mdot_steam):
    """Steam bottoming cycle, the unrolled closed Rankine loop from
    tutorials/05_rankine_steam_cycle (Source pins feedwater, Sink terminates
    the condensate — the condenser outlet lands back on the feed state)."""
    h_feed = water.enthalpy_pq(P_COND, 0.0)
    T_feed = water.temperature_ph(P_COND, h_feed) - 1.0  # slightly subcooled
    feed = Source(name="feed", P=P_COND, T=T_feed, mdot=mdot_steam)
    pump = Pump(name="pump", P_out=P_BOILER, eta=0.75)
    boiler = SimpleEvaporator(name="hrsg", superheat=SUPERHEAT)
    steam_turb = SteamTurbine(name="steam_turb", P_out=P_COND, eta_s=0.80)
    condenser = SimpleCondenser(name="condenser", outlet_quality=0.0)
    drain = Sink(name="drain")

    network = Network(fluid=water)
    for component in (feed, pump, boiler, steam_turb, condenser, drain):
        network.add_component(component)
    network.connect(feed, "out", pump, "in")
    network.connect(pump, "out", boiler, "in")
    network.connect(boiler, "out", steam_turb, "in")
    network.connect(steam_turb, "out", condenser, "in")
    network.connect(condenser, "out", drain, "in")
    return network, pump, boiler, steam_turb, condenser


def solve_combined_cycle(mdot_fuel, N, warm_start=None):
    """One combined-cycle operating point: gas solve, HRSG balance, steam
    solve. Returns everything the charts below need."""
    gas_network, comp, comb, turb = build_gas_turbine(mdot_fuel, N)
    gas_result = gas_network.solve(tol=1e-7, max_iter=1000, damping=0.3, warm_start=warm_start)
    gt = gt_operating_point(gas_result, comp, comb, turb)

    Q_hrsg, mdot_steam, pinch, T_cond = hrsg_steam_mdot(gt)
    steam_network, pump, boiler, steam_turb, condenser = build_steam_cycle(mdot_steam)
    steam_result = steam_network.solve(tol=1e-6, max_iter=200)

    W_st = metrics_for(steam_turb, steam_result)["power [W]"]
    W_pump = metrics_for(pump, steam_result)["power [W]"]
    Q_cond = -metrics_for(condenser, steam_result)["power [W]"]
    gt.update({
        "Q_hrsg": Q_hrsg, "mdot_steam": mdot_steam, "pinch": pinch,
        "W_st_shaft": W_st - W_pump,
        "W_st_elec": (W_st - W_pump) * ETA_GEN_ST,
        "Q_cond": Q_cond,
        "gas_result": gas_result, "steam_result": steam_result,
        "components": (comp, comb, turb, pump, boiler, steam_turb, condenser),
    })
    gt["W_cc_elec"] = gt["W_gt_elec"] + gt["W_st_elec"]
    gt["eta_gt"] = gt["W_gt_elec"] / gt["fuel_power"]
    gt["eta_cc"] = gt["W_cc_elec"] / gt["fuel_power"]
    return gt


def banner(text):
    print(f"\n{'=' * 72}\n{text}\n{'=' * 72}")


# --- 1. Full-load steady state ----------------------------------------------

banner("1. Full-load steady state (design point)")
full = solve_combined_cycle(FUEL_FULL_LOAD, N_GRID)

print("\nGas topping cycle report:")
full["gas_result"].print_report()
print("\nSteam bottoming cycle report:")
full["steam_result"].print_report()

if full["pinch"] < PINCH_MIN:
    raise RuntimeError(
        f"HRSG evaporator pinch {full['pinch']:.1f} K is below the {PINCH_MIN:.0f} K "
        f"minimum — the steam conditions ask for heat the exhaust can't deliver."
    )

print(f"""
Combined-cycle design point summary
  shaft speed                : {N_GRID:9.0f} rev/min (grid-synchronous)
  air mass flow (from maps)  : {full['mdot_air']:9.3f} kg/s
  fuel flow (CH4)            : {full['mdot_fuel'] * 1000:9.2f} g/s   ({full['fuel_power'] / 1e3:6.1f} kW LHV)
  compressor pressure ratio  : {full['PR']:9.2f}
  turbine inlet temperature  : {full['TIT']:9.1f} K ({full['TIT'] - 273.15:6.1f} degC)
  exhaust (HRSG inlet) temp  : {full['T_exhaust']:9.1f} K ({full['T_exhaust'] - 273.15:6.1f} degC)
  HRSG duty                  : {full['Q_hrsg'] / 1e3:9.1f} kW  (stack at {T_STACK - 273.15:.0f} degC)
  HRSG evaporator pinch      : {full['pinch']:9.1f} K
  steam flow                 : {full['mdot_steam']:9.3f} kg/s at {P_BOILER / PA_PER_BAR:.0f} bar
  GT electrical output       : {full['W_gt_elec'] / 1e3:9.1f} kW
  ST electrical output       : {full['W_st_elec'] / 1e3:9.1f} kW
  combined electrical output : {full['W_cc_elec'] / 1e3:9.1f} kW
  GT-only efficiency         : {full['eta_gt'] * 100:9.1f} %
  combined-cycle efficiency  : {full['eta_cc'] * 100:9.1f} %
""")

# --- 2. Part-load sweep: map path, power/efficiency, emissions ----------------

banner("2. Part-load sweep (fuel swept at fixed grid speed, warm-started)")
# Continuation: walk DOWN from the already-converged full-load solution,
# warm-starting each solve from its neighbor — far off the design point a
# cold Newton start can fall outside the maps, but each small step from a
# converged neighbor stays inside. Then flip back to ascending order.
sweep = []
warm = full["gas_result"]
for fuel in reversed(FUEL_SWEEP):
    point = solve_combined_cycle(fuel, N_GRID, warm_start=warm)
    warm = point["gas_result"]
    sweep.append(point)
sweep.reverse()
for point in sweep:
    print(
        f"  fuel {point['mdot_fuel'] * 1000:5.1f} g/s -> GT {point['W_gt_elec'] / 1e3:6.1f} kWe + "
        f"ST {point['W_st_elec'] / 1e3:5.1f} kWe = {point['W_cc_elec'] / 1e3:6.1f} kWe, "
        f"eta_cc {point['eta_cc'] * 100:4.1f} %, TIT {point['TIT']:6.1f} K, "
        f"NO {point['X_NO'] * 1e6:6.1f} ppm"
    )

load_pct = [p["W_cc_elec"] / sweep[-1]["W_cc_elec"] * 100 for p in sweep]

# Compressor map with the traversed operating path (corrected flow B, PR) —
# at fixed N the path climbs the same speed line as backpressure rises.
operating_points = []
for point in sweep:
    B = point["mdot_air"] * T_AMBIENT**0.5 / (P_AMBIENT / PA_PER_BAR)
    operating_points.append((B, point["PR"]))
comp_full = full["components"][0]
ThermoPlot(
    title="Compressor map: part-load operating path at grid speed",
).map(
    comp_full.map, kind="pressure_ratio", operating_points=operating_points,
    operating_label="Load sweep (60 -> 100 %)",
).finish(show=False, save_path="combined_cycle_map.png")
print("Saved combined_cycle_map.png")

power_plot = ThermoPlot(
    title="Combined-cycle power and efficiency vs. load",
    xlabel="load [% of full-load output]", ylabel="electrical power [kW]",
)
power_plot.line(load_pct, [p["W_gt_elec"] / 1e3 for p in sweep], label="GT output")
power_plot.line(load_pct, [p["W_st_elec"] / 1e3 for p in sweep], label="ST output")
power_plot.line(load_pct, [p["W_cc_elec"] / 1e3 for p in sweep], label="combined output")
eff_plot = power_plot.twin_axis(ylabel="efficiency [%]")
eff_plot.line(
    load_pct, [p["eta_cc"] * 100 for p in sweep],
    label="combined efficiency", color="#D55E00", linestyle="--",
)
eff_plot.line(
    load_pct, [p["eta_gt"] * 100 for p in sweep],
    label="GT-only efficiency", color="#CC79A7", linestyle="--",
)
power_plot.finish(show=False, save_path="combined_cycle_load_sweep.png")
print("Saved combined_cycle_load_sweep.png")

# Equilibrium emissions vs. load: NO climbs steeply with load because flame
# temperature rises with fuel flow at fixed speed — the classic thermal-NO
# trend; CO2 tracks fuel-air ratio almost linearly.
emissions_plot = ThermoPlot(
    title="Equilibrium emissions vs. load",
    xlabel="load [% of full-load output]", ylabel="NO [ppmv]",
)
emissions_plot.line(load_pct, [p["X_NO"] * 1e6 for p in sweep], label="NO [ppmv]")
co2_plot = emissions_plot.twin_axis(ylabel="CO2 [mol %]")
co2_plot.line(
    load_pct, [p["X_CO2"] * 100 for p in sweep],
    label="CO2 [mol %]", color="#D55E00",
)
emissions_plot.finish(show=False, save_path="combined_cycle_emissions.png")
print("Saved combined_cycle_emissions.png")

# Full-load energy breakdown: where every kW of the fuel's LHV input ends up.
other = (
    full["fuel_power"] - full["W_cc_elec"] - full["Q_cond"]
    - full["mdot_exhaust"]
    * (air.enthalpy_pt(P_AMBIENT, T_STACK) - air.enthalpy_pt(P_AMBIENT, T_AMBIENT))
)
ThermoPlot(
    title="Full-load energy breakdown", ylabel="power [kW]",
).bar(
    ["fuel (LHV)", "GT elec", "ST elec", "condenser", "stack", "other losses"],
    [
        full["fuel_power"] / 1e3, full["W_gt_elec"] / 1e3, full["W_st_elec"] / 1e3,
        full["Q_cond"] / 1e3,
        full["mdot_exhaust"]
        * (air.enthalpy_pt(P_AMBIENT, T_STACK) - air.enthalpy_pt(P_AMBIENT, T_AMBIENT)) / 1e3,
        other / 1e3,
    ],
).finish(show=False, save_path="combined_cycle_energy_balance.png")
print("Saved combined_cycle_energy_balance.png")

# --- 3. Transient: fuel step and the compressor/spool lag ---------------------

banner("3. Fuel-step transient: compressor/spool lag")

# The transient plant swaps the boundary closure: inlet flow pinned, shaft
# speed freed as a true differential state (rotor inertia). A Schedule steps
# the fuel command instantly at t = 2 s; the spool cannot follow instantly —
# the gap between the two curves IS the compressor lag.
src_t = Source(name="src", P=P_AMBIENT, T=T_AMBIENT, mdot=0.63)
comp_t = Compressor(name="comp", map_path="T100 Comp.cop", gamma=GAMMA, N=None)
comb_t = Combustor(name="comb", PR=0.96, mdot_fuel=0.0085, fuel="CH4")
turb_t = Turbine(name="turb", map_path="T100 Turb.tur", gamma=GAMMA, N=None)
shaft_t = Shaft(
    name="shaft", components=[comp_t, turb_t], signs=[-1.0, 1.0],
    efficiency=0.98, inertia=0.05, dynamic=True, N0=64000.0,
)
snk_t = Sink(name="snk")
fuel_step = Schedule(
    name="fuel_step", target=comb_t, attr="mdot_fuel",
    breakpoints=[(0.0, 0.0085), (2.0, 0.0085), (2.001, 0.0105)],
    interpolation="step",
)

transient_network = Network(fluid=air)
for component in (src_t, comp_t, comb_t, turb_t, shaft_t, snk_t, fuel_step):
    transient_network.add_component(component)
transient_network.connect(src_t, "out", comp_t, "in")
transient_network.connect(comp_t, "out", comb_t, "in")
transient_network.connect(comb_t, "out", turb_t, "in")
transient_network.connect(turb_t, "out", snk_t, "in")

history = transient_network.solve_transient(
    duration=10.0, dt=0.1, tol=1e-8, max_iter=600, damping=0.3,
)
N_series = history.diff_history["shaft.N"]
print(f"  N before step: {N_series[0]:8.0f} rev/min")
print(f"  N at t=10 s  : {N_series[-1]:8.0f} rev/min")

# Fuel command (instant) against shaft speed (lagging) on one figure. The
# command comes from the Schedule's own profile — the combustor's live
# mdot_fuel attribute is whatever the Schedule last set it to, so reading
# its report_metrics() after the run would show only the final value.
fuel_command = [fuel_step.value_at(t) * 1000 for t in history.times]
lag_plot = ThermoPlot(
    title="Compressor lag: instant fuel step, inertia-limited spool response",
    xlabel="t [s]", ylabel="fuel flow [g/s]",
)
lag_plot.line(history.times, fuel_command, label="fuel command [g/s]")
speed_plot = lag_plot.twin_axis(ylabel="shaft speed [rev/min]")
speed_plot.line(
    history.times, N_series, label="shaft speed", color="#D55E00",
)
lag_plot.finish(show=False, save_path="combined_cycle_spool_lag.png")
print("Saved combined_cycle_spool_lag.png")

# The same lag seen by the machines: compressor PR and both shaft powers.
response_plot = ThermoPlot(
    title="Fuel step: compressor and turbine response",
    xlabel="t [s]", ylabel="power [kW]",
)
for component, label in ((comp_t, "compressor power"), (turb_t, "turbine power")):
    y = [metrics_for(component, step)["power [W]"] / 1e3 for step in history.steps]
    response_plot.line(history.times, y, label=label)
pr_plot = response_plot.twin_axis(ylabel="compressor PR [-]")
pr_plot.line(
    history.times,
    [metrics_for(comp_t, step)["PR [-]"] for step in history.steps],
    label="compressor PR", color="#009E73", linestyle="--",
)
response_plot.finish(show=False, save_path="combined_cycle_transient_response.png")
print("Saved combined_cycle_transient_response.png")

# Emissions during the step: equilibrium NO jumps with the fuel step (flame
# temperature jumps first) and keeps evolving as the spool catches up. The
# combustor reads its live mdot_fuel attribute, so replay each step with the
# fuel the Schedule actually commanded at that time before reading metrics.
NO_series = []
for t, step in zip(history.times, history.steps):
    comb_t.mdot_fuel = fuel_step.value_at(t)
    NO_series.append(metrics_for(comb_t, step).get("X_NO [-]", 0.0) * 1e6)
transient_emissions = ThermoPlot(
    title="NO during the fuel step",
    xlabel="t [s]", ylabel="NO [ppmv]",
)
transient_emissions.line(history.times, NO_series, label="NO [ppmv]")
transient_emissions.finish(show=False, save_path="combined_cycle_transient_emissions.png")
print("Saved combined_cycle_transient_emissions.png")

# --- 4. Cycle diagram, annotated with the solved full-load numbers ------------

banner("4. Cycle diagram")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle  # noqa: E402

AIR_C, GAS_C, EXH_C, STEAM_C, WATER_C = "#0072B2", "#D55E00", "#E69F00", "#009E73", "#56B4E9"

fig, ax = plt.subplots(figsize=(13.0, 7.5))
ax.set_xlim(0, 13)
ax.set_ylim(0, 7.5)
ax.axis("off")
ax.set_title(
    "T100-class combined cycle — solved full-load state", fontsize=14, fontweight="bold"
)


def arrow(x0, y0, x1, y1, color, lw=2.2):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
        color=color, linewidth=lw, zorder=3,
    ))


def label(x, y, text, size=8.5, color="#333333", ha="center", weight="normal"):
    ax.text(x, y, text, fontsize=size, color=color, ha=ha, va="center",
            fontweight=weight, zorder=6)


# Gas turbine: compressor and turbine trapezoids joined by a shaft.
ax.add_patch(Polygon([(1.6, 4.0), (1.6, 6.0), (2.8, 5.6), (2.8, 4.4)],
                     closed=True, facecolor="#B8CCE0", edgecolor="#1F2933", zorder=4))
label(2.2, 5.0, "C", 13, weight="bold")
ax.add_patch(Polygon([(4.6, 4.4), (4.6, 5.6), (5.8, 6.0), (5.8, 4.0)],
                     closed=True, facecolor="#F0C4A8", edgecolor="#1F2933", zorder=4))
label(5.2, 5.0, "T", 13, weight="bold")
ax.plot([2.8, 4.6], [5.0, 5.0], color="#1F2933", linewidth=3, zorder=2)  # shaft
ax.plot([5.8, 6.9], [5.0, 5.0], color="#1F2933", linewidth=3, zorder=2)
gen1 = Circle((7.2, 5.0), 0.32, facecolor="#FFD966", edgecolor="#1F2933", zorder=4)
ax.add_patch(gen1)
label(7.2, 5.0, "G", 11, weight="bold")
label(7.2, 4.4, f"{full['W_gt_elec'] / 1e3:.0f} kWe", 9, weight="bold")

# Combustor between them.
comb_circle = Circle((3.7, 6.3), 0.45, facecolor="#F4A0A0", edgecolor="#1F2933", zorder=4)
ax.add_patch(comb_circle)
label(3.7, 6.3, "CC", 10, weight="bold")
arrow(3.7, 7.3, 3.7, 6.8, GAS_C)
label(4.5, 7.15, f"CH4 {full['mdot_fuel'] * 1000:.1f} g/s", 8.5, GAS_C)

# Gas path.
arrow(0.4, 5.0, 1.5, 5.0, AIR_C)
label(0.9, 5.35, f"air {full['mdot_air']:.2f} kg/s\n{T_AMBIENT - 273.15:.0f} degC", 8, AIR_C)
arrow(2.8, 5.15, 3.35, 5.95, AIR_C)
label(2.5, 6.35, f"{full['PR']:.1f} bar", 8, AIR_C)
arrow(4.05, 5.95, 4.6, 5.15, GAS_C)
label(4.9, 6.35, f"TIT {full['TIT'] - 273.15:.0f} degC", 8, GAS_C)
arrow(5.8, 4.6, 5.8, 3.4, EXH_C)  # turbine exhaust down to HRSG
label(6.75, 4.0, f"{full['T_exhaust'] - 273.15:.0f} degC\n{full['mdot_exhaust']:.2f} kg/s",
      8, EXH_C)

# HRSG box with coil.
ax.add_patch(Rectangle((4.9, 1.9), 3.4, 1.5, facecolor="#F5F0E6",
                       edgecolor="#1F2933", zorder=4))
label(6.6, 3.15, "HRSG", 10, weight="bold")
coil_x = [5.2 + 0.28 * i for i in range(11)]
coil_y = [2.35 if i % 2 == 0 else 2.85 for i in range(11)]
ax.plot(coil_x, coil_y, color=STEAM_C, linewidth=2, zorder=5)
label(6.6, 2.1, f"{full['Q_hrsg'] / 1e3:.0f} kW", 8.5, EXH_C)
arrow(4.9, 2.65, 3.9, 2.65, EXH_C)  # stack
arrow(3.9, 2.65, 3.9, 1.5, EXH_C)
label(3.35, 2.0, f"stack\n{T_STACK - 273.15:.0f} degC", 8, EXH_C)

# Steam loop: HRSG -> steam turbine -> condenser -> pump -> HRSG.
ax.add_patch(Polygon([(9.5, 2.2), (9.5, 3.1), (10.7, 3.5), (10.7, 1.8)],
                     closed=True, facecolor="#B7E0CE", edgecolor="#1F2933", zorder=4))
label(10.1, 2.65, "ST", 12, weight="bold")
ax.plot([10.7, 11.5], [2.65, 2.65], color="#1F2933", linewidth=3, zorder=2)
gen2 = Circle((11.8, 2.65), 0.32, facecolor="#FFD966", edgecolor="#1F2933", zorder=4)
ax.add_patch(gen2)
label(11.8, 2.65, "G", 11, weight="bold")
label(11.8, 2.05, f"{full['W_st_elec'] / 1e3:.0f} kWe", 9, weight="bold")

arrow(8.3, 2.95, 9.5, 2.95, STEAM_C)  # live steam
label(8.9, 3.4, f"steam {full['mdot_steam']:.2f} kg/s\n"
      f"{P_BOILER / PA_PER_BAR:.0f} bar / "
      f"{water.saturation_temperature(P_BOILER) + SUPERHEAT - 273.15:.0f} degC",
      8, STEAM_C)
ax.add_patch(Rectangle((9.6, 0.5), 1.6, 0.75, facecolor="#DCEBF5",
                       edgecolor="#1F2933", zorder=4))
label(10.4, 0.87, "condenser", 8.5, weight="bold")
label(10.4, 0.25, f"{full['Q_cond'] / 1e3:.0f} kW rejected", 8, WATER_C)
arrow(10.2, 1.8, 10.4, 1.3, STEAM_C)
pump_circle = Circle((8.6, 0.87), 0.3, facecolor="#DCEBF5", edgecolor="#1F2933", zorder=4)
ax.add_patch(pump_circle)
label(8.6, 0.87, "P", 10, weight="bold")
arrow(9.55, 0.87, 8.95, 0.87, WATER_C)
arrow(8.6, 1.2, 8.6, 2.6, WATER_C)
arrow(8.6, 2.6, 8.3, 2.6, WATER_C)
label(7.95, 1.5, f"{P_COND / PA_PER_BAR:.1f} bar", 8, WATER_C)

# Summary box.
ax.add_patch(Rectangle((0.4, 0.4), 4.4, 1.15, facecolor="#FFFFFF",
                       edgecolor="#1F2933", linewidth=1.2, zorder=4))
label(2.6, 1.22, f"Combined output: {full['W_cc_elec'] / 1e3:.0f} kWe", 10,
      weight="bold")
label(2.6, 0.72,
      f"eta GT-only {full['eta_gt'] * 100:.1f} %  ->  combined {full['eta_cc'] * 100:.1f} %",
      9.5)

fig.savefig("combined_cycle_diagram.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("Saved combined_cycle_diagram.png")

banner("Done — 7 figures saved (combined_cycle_*.png)")
