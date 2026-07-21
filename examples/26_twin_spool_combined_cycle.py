"""Utility-scale twin-spool combined cycle: LP + HP compressors and turbines
on separate spools, topped onto a steam bottoming cycle.

The utility-scale counterpart to examples/25_combined_cycle_showcase.py's
microturbine. The gas generator is an LM6000-class twin-spool aeroderivative
(GE's published architecture: two-shaft, overall pressure ratio 29:1, LP
spool at 3600 rev/min synchronous speed so the generator couples directly
with no gearbox, ~130 kg/s airflow, ~40 MWe — see
https://www.geaerospace.com/sites/default/files/datasheet-lm6000.pdf):

    Source -> LPC -> HPC -> Combustor -> HPT -> LPT -> Sink
                                                 |
    LP spool (3600 rev/min, grid): LPC + LPT + generator
    HP spool (free, ~10000 rev/min): HPC + HPT, dynamic Shaft

Real OEM component maps are proprietary, so the four maps here ("LM6000-class
*.cop/.tur") are the packaged T100 map *shapes* scaled to this machine's
published design points — the standard map-scaling practice when OEM data
isn't available; see generate_lm6000_class_maps.py for exactly how they were
made (they are NOT GE data). At full load the solved cycle lands on the
published architecture: 127 kg/s, OPR ~29, HP spool at ~10000 rev/min.

The interesting twin-spool physics this shows:
  - Only the LP spool speed is imposed (grid). The HP spool speed is a free
    differential state: its steady value is wherever HPT power balances HPC
    power, and in the transient it *floats down* at part load exactly like
    the real machine's gas generator.
  - Mass flow is a Newton unknown everywhere (Source(mdot=None) + ambient
    Sink) — airflow follows the maps as load changes.
  - The fuel-step transient shows the HP spool lag: fuel steps instantly,
    the HP spool (and with it airflow and both compressors' operating
    points) takes seconds to follow — with a utility-scale rotor inertia
    this is much slower than the microturbine's.

The steam bottoming cycle is single-pressure (Pump -> SimpleEvaporator ->
SteamTurbine -> SimpleCondenser at 25 bar), coupled through the same
sequential HRSG energy balance + pinch check as example 25. (Real LM6000
combined cycles use dual-pressure HRSGs and reach ~52%; single-pressure
gives up a few points — the gap is HRSG design, not gas-side physics.)

Requires the cantera, coolprop, and plot extras:
    pip install thermowave[cantera,coolprop,plot]

Runtime is a few minutes (Cantera equilibrium inside every Newton
iteration).

Run: .venv/bin/python examples/26_twin_spool_combined_cycle.py
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

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
water = CoolPropFluid(name="Water")
GAMMA = 1005.0 / (1005.0 - 287.05)

P_AMBIENT = 101325.0
T_AMBIENT = 288.15
N_LP = 3600.0            # rev/min, grid-synchronous LP spool (direct drive)
N_HP_GUESS = 10000.0     # rev/min, HP spool design speed (solved, not imposed)
HP_INERTIA = 30.0        # kg*m^2, utility-scale HP rotor
LHV_CH4 = 50.0e6
ETA_GEN = 0.98           # utility synchronous generator

FUEL_FULL_LOAD = 2.0                    # kg/s CH4
FUEL_SWEEP = [1.4, 1.55, 1.7, 1.85, 2.0]

# Steam bottoming cycle
P_BOILER = 25.0e5
P_COND = 8.0e3           # 0.08 bar, ~41 degC condensing
SUPERHEAT = 120.0        # K above saturation
PINCH_DESIGN = 15.0      # K evaporator pinch the HRSG is sized to
T_STACK_MIN = 373.0      # K, stay above water dewpoint in the stack


def build_gas_turbine(mdot_fuel):
    """Twin-spool gas generator. LPC/LPT speeds are imposed (grid); the HP
    spool is a dynamic Shaft whose speed is a differential state — steady
    solves close it at the HPT-drives-HPC torque balance, transients
    integrate it through rotor inertia."""
    src = Source(name="src", P=P_AMBIENT, T=T_AMBIENT, mdot=None, mdot_guess=127.0)
    lpc = Compressor(name="lpc", map_path="LM6000-class LPC.cop", gamma=GAMMA, N=N_LP)
    hpc = Compressor(name="hpc", map_path="LM6000-class HPC.cop", gamma=GAMMA, N=None)
    comb = Combustor(name="comb", PR=0.96, mdot_fuel=mdot_fuel, fuel="CH4")
    hpt = Turbine(name="hpt", map_path="LM6000-class HPT.tur", gamma=GAMMA, N=None)
    lpt = Turbine(name="lpt", map_path="LM6000-class LPT.tur", gamma=GAMMA, N=N_LP)
    hp_shaft = Shaft(
        name="hp_shaft", components=[hpc, hpt], signs=[-1.0, 1.0],
        efficiency=0.99, inertia=HP_INERTIA, dynamic=True, N0=N_HP_GUESS,
    )
    snk = Sink(name="snk", P=P_AMBIENT)

    network = Network(fluid=air)
    for component in (src, lpc, hpc, comb, hpt, lpt, hp_shaft, snk):
        network.add_component(component)
    network.connect(src, "out", lpc, "in")
    network.connect(lpc, "out", hpc, "in")
    network.connect(hpc, "out", comb, "in")
    network.connect(comb, "out", hpt, "in")
    network.connect(hpt, "out", lpt, "in")
    network.connect(lpt, "out", snk, "in")
    return network, lpc, hpc, comb, hpt, lpt


def metrics_for(component, result):
    state = NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )
    return component.report_metrics(state)


def gt_operating_point(result, lpc, hpc, comb, hpt, lpt):
    lm, hm, bm, htm, ltm = (metrics_for(c, result) for c in (lpc, hpc, comb, hpt, lpt))
    T_exhaust = air.temperature_ph(result.node_P["lpt.out"], result.node_h["lpt.out"])
    W_gt_shaft = ltm["power [W]"] - lm["power [W]"]  # LP spool net = generator input
    return {
        "mdot_air": result.node_mdot["lpc.in"],
        "mdot_fuel": bm["mdot_fuel [kg/s]"],
        "mdot_exhaust": result.node_mdot["lpt.out"],
        "N_HP": result.params["hp_shaft.N"],
        "PR_lpc": lm["PR [-]"], "PR_hpc": hm["PR [-]"],
        "OPR": lm["PR [-]"] * hm["PR [-]"],
        "TIT": bm["T_out [K]"],
        "T_exhaust": T_exhaust,
        "W_lpc": lm["power [W]"], "W_hpc": hm["power [W]"],
        "W_hpt": htm["power [W]"], "W_lpt": ltm["power [W]"],
        "W_gt_shaft": W_gt_shaft,
        "W_gt_elec": W_gt_shaft * ETA_GEN,
        "fuel_power": bm["mdot_fuel [kg/s]"] * LHV_CH4,
        "X_NO": bm.get("X_NO [-]", 0.0),
        "X_CO2": bm.get("X_CO2 [-]", 0.0),
    }


def hrsg_steam_mdot(gt):
    """Pinch-based HRSG sizing — the way single-pressure HRSGs are actually
    rated. The evaporator pinch (gas-to-saturation temperature difference
    where boiling starts) is the design constraint: steam flow is whatever
    the gas can boil+superheat *above* the pinch point, and the stack
    temperature then falls out of the economizer balance below it (rather
    than picking a stack temperature and hoping the pinch works out —
    example 25 does the latter and checks; at utility scale the pinch is
    the binding constraint, so size to it directly)."""
    T_sat = water.saturation_temperature(P_BOILER)
    T_gas_pinch = T_sat + PINCH_DESIGN
    h_feed = water.enthalpy_pq(P_COND, 0.0)
    h_feed_pumped = h_feed + (P_BOILER - P_COND) / water.density_ph(P_COND, h_feed) / 0.75
    h_sat_liquid = water.enthalpy_pq(P_BOILER, 0.0)
    h_steam = water.enthalpy_pt(P_BOILER, T_sat + SUPERHEAT)

    # Evaporator + superheater duty comes from the gas above the pinch.
    Q_above = gt["mdot_exhaust"] * (
        air.enthalpy_pt(P_AMBIENT, gt["T_exhaust"]) - air.enthalpy_pt(P_AMBIENT, T_gas_pinch)
    )
    mdot_steam = Q_above / (h_steam - h_sat_liquid)

    # Economizer below the pinch sets the stack temperature.
    Q_economizer = mdot_steam * (h_sat_liquid - h_feed_pumped)
    h_gas_stack = air.enthalpy_pt(P_AMBIENT, T_gas_pinch) - Q_economizer / gt["mdot_exhaust"]
    T_stack = air.temperature_ph(P_AMBIENT, h_gas_stack)
    if T_stack < T_STACK_MIN:
        raise RuntimeError(
            f"HRSG stack temperature {T_stack:.0f} K fell below {T_STACK_MIN:.0f} K — "
            f"steam conditions recover more heat than the stack can safely give up."
        )
    Q_hrsg = Q_above + Q_economizer
    return Q_hrsg, mdot_steam, T_stack


def build_steam_cycle(mdot_steam):
    h_feed = water.enthalpy_pq(P_COND, 0.0)
    T_feed = water.temperature_ph(P_COND, h_feed) - 1.0
    feed = Source(name="feed", P=P_COND, T=T_feed, mdot=mdot_steam)
    pump = Pump(name="pump", P_out=P_BOILER, eta=0.75)
    boiler = SimpleEvaporator(name="hrsg", superheat=SUPERHEAT)
    steam_turb = SteamTurbine(name="steam_turb", P_out=P_COND, eta_s=0.82)
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


def _solve_gas(mdot_fuel, warm_start, warm_fuel, depth=0):
    """Gas-side solve with adaptive continuation: if a fuel step from the
    warm-start point diverges (twin-spool Newton systems are sensitive to
    step size — a too-large step can throw the HP spool off its speed
    line), bisect the step and walk there in halves. Same technique any
    industrial cycle deck uses to traverse an operating line."""
    from thermowave.core.exceptions import ConvergenceError

    components = build_gas_turbine(mdot_fuel)
    try:
        result = components[0].solve(
            tol=1e-6, max_iter=1500, damping=0.2, warm_start=warm_start,
        )
        return components, result
    except ConvergenceError:
        if warm_start is None or warm_fuel is None or depth >= 3:
            raise
        mid_fuel = 0.5 * (warm_fuel + mdot_fuel)
        _, mid_result = _solve_gas(mid_fuel, warm_start, warm_fuel, depth + 1)
        return _solve_gas(mdot_fuel, mid_result, mid_fuel, depth + 1)


def solve_combined_cycle(mdot_fuel, warm_start=None, warm_fuel=None):
    components, gas_result = _solve_gas(mdot_fuel, warm_start, warm_fuel)
    gas_network, lpc, hpc, comb, hpt, lpt = components
    gt = gt_operating_point(gas_result, lpc, hpc, comb, hpt, lpt)

    Q_hrsg, mdot_steam, T_stack = hrsg_steam_mdot(gt)
    steam_network, pump, boiler, steam_turb, condenser = build_steam_cycle(mdot_steam)
    steam_result = steam_network.solve(tol=1e-6, max_iter=200)

    W_st = metrics_for(steam_turb, steam_result)["power [W]"]
    W_pump = metrics_for(pump, steam_result)["power [W]"]
    gt.update({
        "Q_hrsg": Q_hrsg, "mdot_steam": mdot_steam, "T_stack": T_stack,
        "W_st_elec": (W_st - W_pump) * ETA_GEN,
        "Q_cond": -metrics_for(condenser, steam_result)["power [W]"],
        "gas_result": gas_result, "steam_result": steam_result,
        "components": (lpc, hpc, comb, hpt, lpt),
    })
    gt["W_cc_elec"] = gt["W_gt_elec"] + gt["W_st_elec"]
    gt["eta_gt"] = gt["W_gt_elec"] / gt["fuel_power"]
    gt["eta_cc"] = gt["W_cc_elec"] / gt["fuel_power"]
    return gt


def banner(text):
    print(f"\n{'=' * 72}\n{text}\n{'=' * 72}")


# --- 1. Full-load steady state ----------------------------------------------

banner("1. Full-load steady state (design point)")
full = solve_combined_cycle(FUEL_FULL_LOAD)
full["gas_result"].print_report()

print(f"""
Twin-spool combined-cycle design point
  LP spool speed             : {N_LP:9.0f} rev/min (grid, direct drive)
  HP spool speed (solved)    : {full['N_HP']:9.0f} rev/min
  air mass flow (from maps)  : {full['mdot_air']:9.1f} kg/s
  fuel flow (CH4)            : {full['mdot_fuel']:9.2f} kg/s ({full['fuel_power'] / 1e6:6.1f} MW LHV)
  LPC / HPC pressure ratio   : {full['PR_lpc']:6.2f} / {full['PR_hpc']:5.2f}  (OPR {full['OPR']:.1f})
  turbine inlet temperature  : {full['TIT']:9.0f} K ({full['TIT'] - 273.15:5.0f} degC)
  HPT / LPT power            : {full['W_hpt'] / 1e6:6.1f} / {full['W_lpt'] / 1e6:5.1f} MW
  exhaust (HRSG inlet) temp  : {full['T_exhaust']:9.0f} K ({full['T_exhaust'] - 273.15:5.0f} degC)
  HRSG duty / stack temp     : {full['Q_hrsg'] / 1e6:6.1f} MW / {full['T_stack'] - 273.15:.0f} degC (pinch {PINCH_DESIGN:.0f} K)
  steam flow                 : {full['mdot_steam']:9.2f} kg/s at {P_BOILER / PA_PER_BAR:.0f} bar
  GT electrical output       : {full['W_gt_elec'] / 1e6:9.2f} MW
  ST electrical output       : {full['W_st_elec'] / 1e6:9.2f} MW
  combined electrical output : {full['W_cc_elec'] / 1e6:9.2f} MW
  GT-only efficiency         : {full['eta_gt'] * 100:9.1f} %
  combined-cycle efficiency  : {full['eta_cc'] * 100:9.1f} %
""")

# --- 2. Part-load sweep --------------------------------------------------------

banner("2. Part-load sweep (warm-started continuation, lowest load first)")
# Continuation direction matters on these scaled maps: walking UP in fuel
# converges at every step (each solve warm-started from the point below),
# while walking down from full load loses the HP-spool speed line. The
# lowest-load point starts cold — it's within cold-start range; well below
# it the maps run out of speed lines, exactly like a real machine's
# minimum-load limit.
sweep = []
warm, warm_fuel = None, None
for fuel in FUEL_SWEEP:
    point = solve_combined_cycle(fuel, warm_start=warm, warm_fuel=warm_fuel)
    warm, warm_fuel = point["gas_result"], fuel
    sweep.append(point)
for point in sweep:
    print(
        f"  fuel {point['mdot_fuel']:.2f} kg/s -> {point['W_cc_elec'] / 1e6:5.1f} MWe "
        f"(GT {point['W_gt_elec'] / 1e6:5.1f} + ST {point['W_st_elec'] / 1e6:4.1f}), "
        f"eta_cc {point['eta_cc'] * 100:4.1f} %, N_HP {point['N_HP']:5.0f} rev/min, "
        f"OPR {point['OPR']:4.1f}, NO {point['X_NO'] * 1e6:5.0f} ppm"
    )

load_pct = [p["W_cc_elec"] / sweep[-1]["W_cc_elec"] * 100 for p in sweep]

# Both compressor maps with the traversed part-load path. The LPC runs at
# fixed grid speed (path climbs one speed line); the HPC's own spool slows
# down at part load, so its path crosses speed lines — the twin-spool
# signature.
lpc_full, hpc_full = full["components"][0], full["components"][1]
lpc_points = [
    (p["mdot_air"] * T_AMBIENT**0.5 / (P_AMBIENT / PA_PER_BAR), p["PR_lpc"]) for p in sweep
]
hpc_points = []
for p in sweep:
    r = p["gas_result"]
    P_in, h_in = r.node_P["hpc.in"], r.node_h["hpc.in"]
    T_in = air.temperature_ph(P_in, h_in)
    hpc_points.append((r.node_mdot["hpc.in"] * T_in**0.5 / (P_in / PA_PER_BAR), p["PR_hpc"]))

ThermoPlot(title="LPC map: part-load path at grid speed").map(
    lpc_full.map, kind="pressure_ratio", operating_points=lpc_points,
    operating_label="Load sweep",
).finish(show=False, save_path="twin_spool_lpc_map.png")
ThermoPlot(title="HPC map: part-load path (HP spool floats down)").map(
    hpc_full.map, kind="pressure_ratio", operating_points=hpc_points,
    operating_label="Load sweep",
).finish(show=False, save_path="twin_spool_hpc_map.png")
print("Saved twin_spool_lpc_map.png, twin_spool_hpc_map.png")

power_plot = ThermoPlot(
    title="Twin-spool combined cycle: power and efficiency vs. load",
    xlabel="load [% of full-load output]", ylabel="electrical power [MW]",
)
power_plot.line(load_pct, [p["W_gt_elec"] / 1e6 for p in sweep], label="GT output")
power_plot.line(load_pct, [p["W_st_elec"] / 1e6 for p in sweep], label="ST output")
power_plot.line(load_pct, [p["W_cc_elec"] / 1e6 for p in sweep], label="combined output")
eff_plot = power_plot.twin_axis(ylabel="efficiency [%]")
eff_plot.line(
    load_pct, [p["eta_cc"] * 100 for p in sweep],
    label="combined efficiency", color="#D55E00", linestyle="--",
)
eff_plot.line(
    load_pct, [p["eta_gt"] * 100 for p in sweep],
    label="GT-only efficiency", color="#CC79A7", linestyle="--",
)
power_plot.finish(show=False, save_path="twin_spool_load_sweep.png")
print("Saved twin_spool_load_sweep.png")

# Equilibrium NO/CO2 track fuel-air ratio (flame temperature), which on a
# twin-spool PEAKS just below full load: pushing toward full load the fuel
# leads and the HP spool lags, so FAR is richest around ~90% — at full load
# the spool has caught up, swallows full airflow, and the mixture leans
# back out. A single-shaft machine (example 25) shows no such peak.
emissions_plot = ThermoPlot(
    title="Equilibrium emissions vs. load",
    xlabel="load [% of full-load output]", ylabel="NO [ppmv]",
)
emissions_plot.line(load_pct, [p["X_NO"] * 1e6 for p in sweep], label="NO [ppmv]")
co2_plot = emissions_plot.twin_axis(ylabel="CO2 [mol %]")
co2_plot.line(
    load_pct, [p["X_CO2"] * 100 for p in sweep], label="CO2 [mol %]", color="#D55E00",
)
emissions_plot.finish(show=False, save_path="twin_spool_emissions.png")
print("Saved twin_spool_emissions.png")

# --- 3. Transient: fuel step and the HP spool lag ------------------------------

banner("3. Fuel-step transient: HP spool lag")

# Same twin-spool network, but now solve_transient() integrates the HP
# spool through its rotor inertia while a Schedule steps the fuel command.
# Start from the converged 92%-fuel steady state (continuation) so t=0 is a
# genuine part-load equilibrium.
FUEL_LOW, FUEL_HIGH = 1.85, 2.0
transient_network, lpc_t, hpc_t, comb_t, hpt_t, lpt_t = build_gas_turbine(FUEL_LOW)
fuel_step = Schedule(
    name="fuel_step", target=comb_t, attr="mdot_fuel",
    breakpoints=[(0.0, FUEL_LOW), (2.0, FUEL_LOW), (2.001, FUEL_HIGH)],
    interpolation="step",
)
transient_network.add_component(fuel_step)

# Warm-start from the sweep's own converged solution at this same fuel flow
# (stepping down from full load is the direction that doesn't converge).
part_load = next(p for p in sweep if abs(p["mdot_fuel"] - FUEL_LOW) < 1e-9)
initial = transient_network.solve(
    tol=1e-6, max_iter=1500, damping=0.2, warm_start=part_load["gas_result"],
)
history = transient_network.solve_transient(
    duration=12.0, dt=0.1, initial=initial, tol=1e-6, max_iter=1500, damping=0.2,
)
N_series = history.diff_history["hp_shaft.N"]
print(f"  HP spool before step: {N_series[0]:7.0f} rev/min")
print(f"  HP spool at t=12 s  : {N_series[-1]:7.0f} rev/min")

fuel_command = [fuel_step.value_at(t) for t in history.times]
lag_plot = ThermoPlot(
    title="HP spool lag: instant fuel step, inertia-limited response",
    xlabel="t [s]", ylabel="fuel flow [kg/s]",
)
lag_plot.line(history.times, fuel_command, label="fuel command [kg/s]")
speed_plot = lag_plot.twin_axis(ylabel="HP spool speed [rev/min]")
speed_plot.line(history.times, N_series, label="HP spool speed", color="#D55E00")
lag_plot.finish(show=False, save_path="twin_spool_lag.png")
print("Saved twin_spool_lag.png")

# Airflow and both compressors chase the HP spool — the whole gas path lags
# the fuel command, not just the shaft speed.
response_plot = ThermoPlot(
    title="Fuel step: airflow and compressor response",
    xlabel="t [s]", ylabel="air mass flow [kg/s]",
)
response_plot.line(
    history.times, [step.node_mdot["lpc.in"] for step in history.steps],
    label="air mass flow",
)
pr_plot = response_plot.twin_axis(ylabel="pressure ratio [-]")
pr_plot.line(
    history.times, [metrics_for(hpc_t, step)["PR [-]"] for step in history.steps],
    label="HPC pressure ratio", color="#D55E00", linestyle="--",
)
pr_plot.line(
    history.times, [metrics_for(lpc_t, step)["PR [-]"] for step in history.steps],
    label="LPC pressure ratio", color="#009E73", linestyle="--",
)
response_plot.finish(show=False, save_path="twin_spool_transient_response.png")
print("Saved twin_spool_transient_response.png")

# --- 4. Cycle diagram, annotated with the solved full-load numbers ------------

banner("4. Cycle diagram")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Rectangle  # noqa: E402

AIR_C, GAS_C, EXH_C, STEAM_C, WATER_C = "#0072B2", "#D55E00", "#E69F00", "#009E73", "#56B4E9"

fig, ax = plt.subplots(figsize=(14.0, 8.0))
ax.set_xlim(0, 14)
ax.set_ylim(0, 8)
ax.axis("off")
ax.set_title(
    "LM6000-class twin-spool combined cycle — solved full-load state",
    fontsize=14, fontweight="bold",
)


def arrow(x0, y0, x1, y1, color, lw=2.2):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=16,
        color=color, linewidth=lw, zorder=3,
    ))


def label(x, y, text, size=8.5, color="#333333", ha="center", weight="normal"):
    ax.text(x, y, text, fontsize=size, color=color, ha=ha, va="center",
            fontweight=weight, zorder=6)


def machine(x, y, w, h, kind, face, name):
    """kind: 'comp' narrows left->right outlet, 'turb' widens."""
    if kind == "comp":
        pts = [(x, y - h / 2), (x, y + h / 2), (x + w, y + h * 0.32), (x + w, y - h * 0.32)]
    else:
        pts = [(x, y - h * 0.32), (x, y + h * 0.32), (x + w, y + h / 2), (x + w, y - h / 2)]
    ax.add_patch(Polygon(pts, closed=True, facecolor=face, edgecolor="#1F2933", zorder=4))
    label(x + w / 2, y, name, 11, weight="bold")


Y = 5.6  # gas-path centerline
machine(1.3, Y, 1.1, 1.9, "comp", "#B8CCE0", "LPC")
machine(3.0, Y, 1.1, 1.4, "comp", "#8FB3D9", "HPC")
machine(6.1, Y, 1.1, 1.4, "turb", "#F0C4A8", "HPT")
machine(7.8, Y, 1.1, 1.9, "turb", "#E8A87C", "LPT")

# Combustor between HPC and HPT.
comb_circle = Circle((5.15, Y + 1.35), 0.45, facecolor="#F4A0A0", edgecolor="#1F2933", zorder=4)
ax.add_patch(comb_circle)
label(5.15, Y + 1.35, "CC", 10, weight="bold")
arrow(5.15, Y + 2.35, 5.15, Y + 1.85, GAS_C)
label(6.15, Y + 2.2, f"CH4 {full['mdot_fuel']:.2f} kg/s", 8.5, GAS_C)

# Spools: HP shaft (inner, short) and LP shaft (outer, runs through to the
# generator on the right).
ax.plot([4.1, 6.1], [Y, Y], color="#1F2933", linewidth=5, zorder=2, solid_capstyle="butt")
label(5.1, Y - 0.35, f"HP spool {full['N_HP']:.0f} rev/min", 8, ha="center")
ax.plot([2.4, 3.0], [Y, Y], color="#6B7280", linewidth=2.5, zorder=1)
ax.plot([6.95, 7.8], [Y, Y], color="#6B7280", linewidth=2.5, zorder=1)
ax.plot([8.9, 10.0], [Y, Y], color="#6B7280", linewidth=2.5, zorder=1)
gen1 = Circle((10.3, Y), 0.34, facecolor="#FFD966", edgecolor="#1F2933", zorder=4)
ax.add_patch(gen1)
label(10.3, Y, "G", 11, weight="bold")
label(10.3, Y - 0.62, f"{full['W_gt_elec'] / 1e6:.1f} MWe", 9, weight="bold")
label(9.4, Y + 0.3, f"LP spool {N_LP:.0f} rev/min", 8)

# Gas path.
arrow(0.2, Y, 1.2, Y, AIR_C)
label(0.65, Y + 0.42, f"air {full['mdot_air']:.0f} kg/s", 8, AIR_C)
arrow(2.45, Y + 0.75, 2.95, Y + 0.75, AIR_C, lw=1.8)
label(2.7, Y + 1.05, f"{full['PR_lpc']:.1f} bar", 8, AIR_C)
arrow(4.1, Y + 0.6, 4.7, Y + 1.2, AIR_C)
label(4.0, Y + 1.55, f"{full['OPR'] * P_AMBIENT / PA_PER_BAR / 1e0:.0f} bar", 8, AIR_C)
arrow(5.6, Y + 1.2, 6.2, Y + 0.6, GAS_C)
label(6.6, Y + 1.5, f"TIT {full['TIT'] - 273.15:.0f} degC", 8, GAS_C)
arrow(8.9, Y - 0.4, 8.9, 3.3, EXH_C)
label(9.75, 4.35, f"{full['T_exhaust'] - 273.15:.0f} degC\n{full['mdot_exhaust']:.0f} kg/s", 8, EXH_C)

# HRSG.
ax.add_patch(Rectangle((6.6, 1.8), 3.6, 1.5, facecolor="#F5F0E6", edgecolor="#1F2933", zorder=4))
label(8.4, 3.05, "HRSG", 10, weight="bold")
coil_x = [6.9 + 0.3 * i for i in range(11)]
coil_y = [2.25 if i % 2 == 0 else 2.75 for i in range(11)]
ax.plot(coil_x, coil_y, color=STEAM_C, linewidth=2, zorder=5)
label(8.4, 2.0, f"{full['Q_hrsg'] / 1e6:.1f} MW", 8.5, EXH_C)
arrow(6.6, 2.55, 5.6, 2.55, EXH_C)
arrow(5.6, 2.55, 5.6, 1.3, EXH_C)
label(5.0, 1.9, f"stack\n{full['T_stack'] - 273.15:.0f} degC", 8, EXH_C)

# Steam loop.
machine(11.0, 2.55, 1.2, 1.5, "turb", "#B7E0CE", "ST")
ax.plot([12.2, 12.9], [2.55, 2.55], color="#1F2933", linewidth=3, zorder=2)
gen2 = Circle((13.2, 2.55), 0.32, facecolor="#FFD966", edgecolor="#1F2933", zorder=4)
ax.add_patch(gen2)
label(13.2, 2.55, "G", 11, weight="bold")
label(13.2, 1.95, f"{full['W_st_elec'] / 1e6:.1f} MWe", 9, weight="bold")

arrow(10.2, 2.85, 11.0, 2.85, STEAM_C)
label(10.6, 3.35,
      f"{full['mdot_steam']:.1f} kg/s\n{P_BOILER / PA_PER_BAR:.0f} bar / "
      f"{water.saturation_temperature(P_BOILER) + SUPERHEAT - 273.15:.0f} degC",
      8, STEAM_C)
ax.add_patch(Rectangle((10.9, 0.35), 1.6, 0.7, facecolor="#DCEBF5", edgecolor="#1F2933", zorder=4))
label(11.7, 0.7, "condenser", 8.5, weight="bold")
label(11.7, 0.1, f"{full['Q_cond'] / 1e6:.1f} MW rejected", 8, WATER_C)
arrow(11.6, 1.9, 11.7, 1.1, STEAM_C)
pump_circle = Circle((9.8, 0.7), 0.28, facecolor="#DCEBF5", edgecolor="#1F2933", zorder=4)
ax.add_patch(pump_circle)
label(9.8, 0.7, "P", 10, weight="bold")
arrow(10.85, 0.7, 10.1, 0.7, WATER_C)
arrow(9.8, 1.0, 9.8, 2.4, WATER_C)
arrow(9.8, 2.4, 10.2, 2.4, WATER_C)

# Summary box.
ax.add_patch(Rectangle((0.4, 0.4), 4.0, 1.3, facecolor="#FFFFFF",
                       edgecolor="#1F2933", linewidth=1.2, zorder=4))
label(2.4, 1.35, f"Combined output: {full['W_cc_elec'] / 1e6:.1f} MWe", 10.5, weight="bold")
label(2.4, 0.85, f"OPR {full['OPR']:.1f}   eta GT {full['eta_gt'] * 100:.1f} % "
      f"->  combined {full['eta_cc'] * 100:.1f} %", 9)

fig.savefig("twin_spool_diagram.png", dpi=160, bbox_inches="tight")
plt.close(fig)
print("Saved twin_spool_diagram.png")

banner("Done — 7 figures saved (twin_spool_*.png)")
