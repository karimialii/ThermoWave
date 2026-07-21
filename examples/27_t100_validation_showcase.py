"""T100 microturbine: validation against test-rig data + two-PID closed loop.

The package's validation showcase: one clean model of the Turbec T100
recuperated microturbine — its real compressor/turbine maps ("T100
Comp.cop"/"T100 Turb.tur"), a Cantera equilibrium combustor, casing/shaft
thermal masses with real two-way heat-transfer coupling, and the actual
machine's two-loop control architecture — validated point-by-point against
measured steady-state data from the University of Stavanger's T100 rig,
then driven through a multi-step power-setpoint transient.

The control architecture (both sections, same two loops):
  - TOT loop:   fuel flow holds the turbine outlet temperature at the
    machine's 645 degC limit (the T100's own control strategy).
  - Power loop: shaft speed is commanded to meet the electrical power
    setpoint (the power electronics' role in the real machine).
Section 1 uses the ideal steady-state `Controller` for both (infinite-gain:
pins the target exactly); section 2 swaps in two `PIDController`s — the
same loops with finite gains, showing real settling, lag, and coupling.

Section 1 — steady-state validation. 16 measured operating points
(35-70 kW power setpoints at each point's own inlet pressure/temperature)
are re-solved and compared against the rig's measured shaft speed,
compressor outlet temperature (COT), and compressor outlet pressure (COP).
Calibration policy: the compressor map is used exactly as shipped; three
machine-level constants are set once against the data and then never
touched per-point — the turbine map's efficiency factor scaled by
K_TURB = 0.92 (this rig delivers full power only at 100% speed, i.e. it
performs below the nominal map — `factor_overrides` is the package's
built-in mechanism for exactly this calibration), generator+power-
electronics efficiency ETA_GEN = 0.90, and recuperator effectiveness 0.87
(published T100 value). Everything else follows from the maps and physics.
Humidity (the rig's RH/WAR columns) is not modeled — dry air; the largest
recorded water-air ratio is 0.0077, a percent-level effect.

Section 2 — transient. A Schedule steps the power setpoint through
multiple levels (35 -> 49 -> 63 -> 70 -> 56 kW); the two PIDs track it.
Every response plot carries the power setpoint on the right axis. The
compressor lag is visible at each step: speed, airflow, and pressure
ratio trail the instantaneous setpoint by several seconds of closed-loop
settling, while the casing/shaft thermal masses soak over minutes behind
that — three separated timescales from one model.

Validation data source: "T100 Validation Data for Ali.xlsx" (UiS T100
rig, 2025) — embedded below so this example runs standalone. Power_Set
in that sheet is percent of 70 kW; Speed is percent of 70,000 rev/min.

Requires the cantera and plot extras: pip install thermowave[cantera,plot]
Runtime is a few minutes (Cantera equilibrium inside every Newton
iteration; 16 steady solves + a 200-step transient).

Run: .venv/bin/python examples/27_t100_validation_showcase.py
"""

from thermowave.components.combustor import Combustor
from thermowave.components.compressor import Compressor
from thermowave.components.controller import Controller
from thermowave.components.heat_transfer import Conduction, Convection, ThermalMass
from thermowave.components.pid_controller import PIDController
from thermowave.components.schedule import Schedule
from thermowave.components.sensor import Sensor
from thermowave.components.shaft import Shaft
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.components.simple_heat_exchanger import SimpleHeatExchanger
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.exceptions import ConvergenceError
from thermowave.core.network import Network, NetworkState
from thermowave.core.plotting import ThermoPlot
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)

# --- Machine constants (calibrated once against the rig data, see docstring) --
TOT_SETPOINT = 645.0 + 273.15   # K, the T100's turbine-outlet-temperature limit
K_TURB = 0.92                   # turbine map E_fact calibration factor
ETA_GEN = 0.90                  # generator + power electronics
RECUP_EFFECTIVENESS = 0.87      # published T100 recuperator rating
N_RATED = 70000.0               # rev/min = 100% speed
P_RATED = 70000.0               # W = 100% power on this rig

# --- Measured steady-state data: UiS T100 rig ---------------------------------
# (P_set [W], inlet P [bar], inlet T [degC], measured speed [% of 70000],
#  measured COT [degC], measured COP [Pa])
VALIDATION_DATA = [
    (35000, 0.99750375, 8.0455, 81.39, 147.0167, 308813.63),
    (35000, 0.99860475, 4.6335, 81.86, 144.7560, 307723.27),
    (35000, 0.99812325, 10.0325, 81.52, 150.0670, 306827.40),
    (42000, 0.99847450, 4.8705, 84.96, 156.1147, 333136.83),
    (42000, 0.99751275, 9.3380, 85.04, 160.1463, 332673.80),
    (49000, 0.99846050, 6.4455, 88.01, 167.7557, 356290.30),
    (49000, 0.99713500, 11.0087, 89.28, 176.4583, 352741.13),
    (49000, 0.99750275, 9.6875, 87.93, 170.7160, 356397.80),
    (56000, 0.99821425, 10.0275, 90.63, 179.8607, 377560.83),
    (56000, 0.99693667, 11.2450, 91.22, 183.9700, 377158.00),
    (56000, 0.99772775, 13.2635, 91.29, 185.9000, 379682.47),
    (63000, 0.99696800, 11.6100, 94.11, 194.2867, 402188.50),
    (63000, 0.99751975, 11.2805, 93.72, 192.4923, 405360.23),
    (70000, 0.99667533, 11.2270, 96.99, 204.5887, 428851.20),
    (70000, 0.99756625, 9.9190, 96.53, 201.0590, 431206.87),
    (70000, 0.99138000, 18.0000, 100.00, 211.9400, 439000.00),
]


def build_t100(P_amb, T_amb, k_turb=K_TURB):
    """The T100 model both sections share. Flow path:
    Source -> Compressor -> recuperator (cold) -> Combustor -> Turbine ->
    recuperator (hot) -> Sink, with mass flow a Newton unknown (Source
    mdot=None + Sink pinning ambient back-pressure) and each machine's
    outlet fluid convecting into its own casing ThermalMass, the casings
    conducting through a shaft mass, and both casings convecting to
    ambient — the example-24 heat-transfer pattern, coupled two-way into
    the energy balances via heat_path."""
    src = Source(name="src", P=P_amb, T=T_amb, mdot=None, mdot_guess=0.7)
    comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=GAMMA, N=None)
    recup = SimpleHeatExchanger(
        name="recup", effectiveness=RECUP_EFFECTIVENESS, PR_hot=0.98, PR_cold=0.97,
    )
    comb = Combustor(name="comb", PR=0.96, mdot_fuel=None, fuel="CH4")
    turb = Turbine(
        name="turb", map_path="T100 Turb.tur", gamma=GAMMA, N=None,
        factor_overrides={"E_fact": 0.961 * k_turb},
    )
    shaft = Shaft(name="shaft", components=[comp, turb], signs=[-1.0, 1.0], efficiency=0.98)
    gen = SimpleGenerator(name="gen", component=shaft, efficiency=ETA_GEN)
    tot_sensor = Sensor(name="tot_sensor")
    snk = Sink(name="snk", P=P_amb)

    comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=2000.0, T0=T_amb)
    turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=2000.0, T0=T_amb)
    shaft_mass = ThermalMass(name="shaft_mass", thermal_capacitance=500.0, T0=T_amb)
    conv_comp = Convection(name="conv_comp", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3)
    conv_turb = Convection(name="conv_turb", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3)
    cond_ts = Conduction(name="cond_ts", a=turb_casing, b=shaft_mass, k=15.0, A=0.01, L=0.4)
    cond_sc = Conduction(name="cond_sc", a=shaft_mass, b=comp_casing, k=15.0, A=0.01, L=0.4)
    conv_ta = Convection(name="conv_ta", a=turb_casing, b=T_amb, h=10.0, A=1.0)
    conv_ca = Convection(name="conv_ca", a=comp_casing, b=T_amb, h=10.0, A=1.0)
    comp.heat_path = conv_comp
    turb.heat_path = conv_turb
    comp_casing.heat_sources = [(conv_comp, 1.0), (cond_sc, 1.0), (conv_ca, -1.0)]
    turb_casing.heat_sources = [(conv_turb, 1.0), (cond_ts, -1.0), (conv_ta, -1.0)]
    shaft_mass.heat_sources = [(cond_ts, 1.0), (cond_sc, -1.0)]

    network = Network(fluid=air)
    for component in (
        src, comp, recup, comb, turb, shaft, gen, tot_sensor, snk,
        comp_casing, turb_casing, shaft_mass,
        conv_comp, conv_turb, cond_ts, cond_sc, conv_ta, conv_ca,
    ):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", recup, "cold_in")
    network.connect(recup, "cold_out", comb, "in")
    network.connect(comb, "out", turb, "in")
    network.connect(turb, "out", tot_sensor, "tap")
    network.connect(turb, "out", recup, "hot_in")
    network.connect(recup, "hot_out", snk, "in")

    parts = {
        "network": network, "src": src, "comp": comp, "recup": recup, "comb": comb,
        "turb": turb, "shaft": shaft, "gen": gen, "tot_sensor": tot_sensor,
        "comp_casing": comp_casing, "turb_casing": turb_casing, "shaft_mass": shaft_mass,
    }
    return parts


def state_of(result):
    return NetworkState(
        fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
        node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
    )


def banner(text):
    print(f"\n{'=' * 78}\n{text}\n{'=' * 78}")


# =============================================================================
# Section 1 — steady-state validation against the rig data
# =============================================================================

banner("1. Steady-state validation against UiS T100 rig data")


def solve_steady(P_set, P_amb, T_amb, k_turb=K_TURB, warm=None, warm_power=None, depth=0):
    """One steady point: the two ideal Controllers pin generator power and
    TOT exactly. If a power step from the warm-start point diverges, walk
    there in halves (same adaptive continuation as example 26)."""
    parts = build_t100(P_amb, T_amb, k_turb=k_turb)
    ctrl_power = Controller(
        name="ctrl_power", sensor=parts["gen"], quantity="power [W]",
        component=parts["comp"], free_param="N", value=P_set,
    )
    ctrl_tot = Controller(
        name="ctrl_tot", sensor=parts["tot_sensor"], quantity="T [K]",
        component=parts["comb"], free_param="mdot_fuel", value=TOT_SETPOINT,
    )
    parts["network"].add_component(ctrl_power)
    parts["network"].add_component(ctrl_tot)
    try:
        result = parts["network"].solve(tol=1e-7, max_iter=800, damping=0.3, warm_start=warm)
        return parts, result
    except ConvergenceError:
        if warm is None or warm_power is None or depth >= 3:
            raise
        mid = 0.5 * (warm_power + P_set)
        _, mid_result = solve_steady(mid, P_amb, T_amb, k_turb, warm, warm_power, depth + 1)
        return solve_steady(P_set, P_amb, T_amb, k_turb, mid_result, mid, depth + 1)


# Branch seeding: this cycle has two mathematically valid steady states per
# power setpoint (the power-vs-speed curve folds); the physical one — the
# branch the rig actually runs on, reaching full power at ~100% speed — is
# found by a cold solve at full load with the nominal (k_turb=1) map, then
# walking the calibration factor in. Every later point warm-starts from a
# neighbor, staying on that branch.
seed_row = VALIDATION_DATA[13]
warm = None
for k_seed in (1.0, K_TURB):
    _, warm = solve_steady(
        float(seed_row[0]), seed_row[1] * 1e5, seed_row[2] + 273.15, k_turb=k_seed, warm=warm,
    )

results = []
warm_power = float(seed_row[0])
for row in sorted(VALIDATION_DATA, key=lambda r: -r[0]):
    P_set, cip_bar, cit_c, meas_speed, meas_cot, meas_cop = row
    parts, r = solve_steady(
        float(P_set), cip_bar * 1e5, cit_c + 273.15, warm=warm, warm_power=warm_power,
    )
    warm, warm_power = r, float(P_set)
    model_speed = r.params["comp.N"] / N_RATED * 100.0
    model_cot = air.temperature_ph(r.node_P["comp.out"], r.node_h["comp.out"]) - 273.15
    model_cop = r.node_P["comp.out"]
    results.append({
        "P_set": P_set, "cit": cit_c,
        "meas": (meas_speed, meas_cot, meas_cop),
        "model": (model_speed, model_cot, model_cop),
        "err": (
            (model_speed - meas_speed) / meas_speed * 100.0,
            (model_cot - meas_cot) / meas_cot * 100.0,
            (model_cop - meas_cop) / meas_cop * 100.0,
        ),
        "fuel": r.params["comb.mdot_fuel"],
    })
results.sort(key=lambda p: (p["P_set"], p["cit"]))

print(f"\n{'P_set':>6} {'T_in':>6} │ {'speed [%]':^21} │ {'COT [degC]':^21} │ {'COP [bar]':^21}")
print(f"{'[kW]':>6} {'[degC]':>6} │ {'meas':>6} {'model':>6} {'err':>7} │ "
      f"{'meas':>6} {'model':>6} {'err':>7} │ {'meas':>6} {'model':>6} {'err':>7}")
print("─" * 100)
for p in results:
    (ms, mc, mp), (Ms, Mc, Mp), (es, ec, ep) = p["meas"], p["model"], p["err"]
    print(f"{p['P_set'] / 1e3:6.0f} {p['cit']:6.1f} │ {ms:6.2f} {Ms:6.2f} {es:+6.2f}% │ "
          f"{mc:6.1f} {Mc:6.1f} {ec:+6.2f}% │ {mp / 1e5:6.3f} {Mp / 1e5:6.3f} {ep:+6.2f}%")

for i, quantity in enumerate(("speed", "COT", "COP")):
    errors = [abs(p["err"][i]) for p in results]
    print(f"  {quantity:>5}: mean abs error {sum(errors) / len(errors):4.2f} %, "
          f"max {max(errors):4.2f} %")

# --- Validation charts --------------------------------------------------------
P_kw = [p["P_set"] / 1e3 for p in results]

speed_plot = ThermoPlot(
    title="T100 validation: shaft speed vs. power setpoint",
    xlabel="power setpoint [kW]", ylabel="shaft speed [% of 70,000 rev/min]",
)
speed_plot.scatter(P_kw, [p["meas"][0] for p in results], label="measured (UiS rig)",
                   marker="s", s=45, zorder=5)
speed_plot.line(P_kw, [p["model"][0] for p in results], label="ThermoWave model",
                linestyle="--")
speed_plot.finish(show=False, save_path="t100_validation_speed.png")
print("\nSaved t100_validation_speed.png")

cot_plot = ThermoPlot(
    title="T100 validation: compressor outlet state",
    xlabel="power setpoint [kW]", ylabel="COT [degC]",
)
cot_plot.scatter(P_kw, [p["meas"][1] for p in results], label="COT measured",
                 marker="s", s=45, zorder=5)
cot_plot.line(P_kw, [p["model"][1] for p in results], label="COT model", linestyle="--")
cop_plot = cot_plot.twin_axis(ylabel="COP [bar]")
cop_plot.scatter(P_kw, [p["meas"][2] / 1e5 for p in results], label="COP measured",
                 marker="D", s=40, color="#D55E00", zorder=5)
cop_plot.line(P_kw, [p["model"][2] / 1e5 for p in results], label="COP model",
              linestyle="--", color="#D55E00")
cot_plot.finish(show=False, save_path="t100_validation_cot_cop.png")
print("Saved t100_validation_cot_cop.png")

# Per-power-level mean error bars (grouping repeated points at one setpoint).
levels = sorted({p["P_set"] for p in results})
level_labels = [f"{lvl / 1e3:.0f} kW" for lvl in levels]
error_plot = ThermoPlot(
    title="T100 validation: mean error by power level", ylabel="error [%]",
)
for i, (quantity, color) in enumerate(
    (("speed", "#0072B2"), ("COT", "#D55E00"), ("COP", "#009E73"))
):
    means = [
        sum(p["err"][i] for p in results if p["P_set"] == lvl)
        / sum(1 for p in results if p["P_set"] == lvl)
        for lvl in levels
    ]
    error_plot.bar(level_labels, means, label=quantity, color=color)
error_plot.finish(show=False, save_path="t100_validation_error.png")
print("Saved t100_validation_error.png")

# =============================================================================
# Section 2 — transient: two PIDs tracking a multi-step power schedule
# =============================================================================

banner("2. Transient: two PID loops, multiple power setpoints, compressor lag")

T_AMB, P_AMB = 288.15, 101325.0

# Find the true 35 kW operating point at ISO conditions first (same ideal
# Controllers as section 1, warm-started from its 35 kW result), so both
# PIDs' output0 start the transient exactly ON the initial setpoints —
# a real dispatch scenario starts from an operating plant, not a cold guess.
_, iso_35 = solve_steady(35000.0, P_AMB, T_AMB, warm=warm, warm_power=warm_power)
N_35 = iso_35.params["comp.N"]
fuel_35 = iso_35.params["comb.mdot_fuel"]
print(f"  35 kW ISO operating point: N = {N_35:.0f} rev/min, "
      f"fuel = {fuel_35 * 1e3:.2f} g/s")

parts = build_t100(P_AMB, T_AMB)
comp_t, comb_t, gen_t = parts["comp"], parts["comb"], parts["gen"]
tot_sensor_t = parts["tot_sensor"]

pid_fuel = PIDController(
    name="pid_fuel", sensor=tot_sensor_t, quantity="T [K]",
    component=comb_t, free_param="mdot_fuel", setpoint=TOT_SETPOINT,
    Kp=6.0e-6, Ki=3.0e-6, Kd=0.0,
    output0=fuel_35, output_min=0.0012, output_max=0.0080,
)
# Positive gains: on this recuperated, free-airflow cycle more speed means
# more power (the validation curve: 81% speed at 35 kW up to ~97% at 70 kW).
# Note tutorial 04's power loop uses NEGATIVE gains — its fixed-airflow
# cycle has the opposite power-vs-speed slope. Same loop, opposite plant
# sign: get this wrong and the loop runs away instead of settling.
pid_power = PIDController(
    name="pid_power", sensor=gen_t, quantity="power [W]",
    component=comp_t, free_param="N", setpoint=35000.0,
    Kp=0.05, Ki=0.03, Kd=0.0,
    output0=N_35, output_min=50000.0, output_max=70000.0,
)
# (t [s], electrical power setpoint [W]) — step interpolation: each new
# setpoint arrives as an instantaneous dispatch command.
POWER_PROFILE = [
    (0.0, 35000.0), (30.0, 35000.0),
    (30.5, 49000.0), (70.0, 49000.0),
    (70.5, 63000.0), (110.0, 63000.0),
    (110.5, 70000.0), (150.0, 70000.0),
    (150.5, 56000.0), (200.0, 56000.0),
]
sched_power = Schedule(
    name="sched_power", target=pid_power, attr="setpoint",
    breakpoints=POWER_PROFILE, interpolation="step",
)
for component in (pid_fuel, pid_power, sched_power):
    parts["network"].add_component(component)

# t=0: the PID residuals pin N/fuel to output0, which IS the 35 kW point,
# so this initial solve just re-finds the same state on this network.
initial = parts["network"].solve(tol=1e-7, max_iter=800, damping=0.3, warm_start=iso_35)

# Heavier damping than section 1's steady solves: mid-transient the plant
# crosses the map region where PR-vs-flow goes flat (the same region behind
# the steady N(P) plateau), and Newton needs the smaller steps there.
history = parts["network"].solve_transient(
    duration=200.0, dt=1.0, initial=initial, tol=1e-7, max_iter=3000, damping=0.15,
)

power_set_kw = [sched_power.value_at(t) / 1e3 for t in history.times]
power_kw, speed_pct, mdot_series, pr_series, tot_series = [], [], [], [], []
for step in history.steps:
    s = state_of(step)
    power_kw.append(gen_t.report_metrics(s)["power [W]"] / 1e3)
    speed_pct.append(step.params["comp.N"] / N_RATED * 100.0)
    mdot_series.append(step.node_mdot["comp.in"])
    pr_series.append(comp_t.report_metrics(s)["PR [-]"])
    tot_series.append(tot_sensor_t.report_metrics(s)["T [K]"] - 273.15)

print(f"  power at t=0    : {power_kw[0]:6.1f} kW (setpoint {power_set_kw[0]:.0f} kW)")
print(f"  power at t=200 s: {power_kw[-1]:6.1f} kW (setpoint {power_set_kw[-1]:.0f} kW)")
print(f"  TOT excursion   : {min(tot_series):.1f} .. {max(tot_series):.1f} degC "
      f"(setpoint {TOT_SETPOINT - 273.15:.0f} degC)")


def with_power_set_axis(plot):
    """Overlay the commanded power setpoint on the right axis — the user-
    facing convention for every transient chart in this example."""
    right = plot.twin_axis(ylabel="power setpoint [kW]")
    right.line(history.times, power_set_kw, label="power setpoint",
               color="#999999", linestyle=":", marker=None, linewidth=2.0)
    return right


# The visible ringing after the first (35 -> 49 kW) step is genuine
# finite-gain controller behavior, not numerical noise: that step drives
# the plant across the map region where power-vs-speed goes locally flat
# (the same region behind section 1's N(P) plateau), so the loop gain
# effectively spikes and the PID rings before settling. The later steps
# traverse better-behaved map regions and settle smoothly.
power_plot = ThermoPlot(
    title="Power tracking: two-PID response to dispatch steps",
    xlabel="t [s]", ylabel="electrical power [kW]",
)
power_plot.line(history.times, power_kw, label="delivered power")
with_power_set_axis(power_plot)
power_plot.finish(show=False, save_path="t100_transient_power.png")
print("Saved t100_transient_power.png")

lag_plot = ThermoPlot(
    title="Compressor lag: speed and airflow trail each power step",
    xlabel="t [s]", ylabel="shaft speed [%]  /  PR [-] x 20",
)
lag_plot.line(history.times, speed_pct, label="shaft speed [%]")
lag_plot.line(history.times, [pr * 20 for pr in pr_series], label="compressor PR x 20")
lag_plot.line(history.times, [m * 100 for m in mdot_series], label="air mdot [kg/s] x 100")
with_power_set_axis(lag_plot)
lag_plot.finish(show=False, save_path="t100_transient_lag.png")
print("Saved t100_transient_lag.png")

tot_plot = ThermoPlot(
    title="TOT held by the fuel PID through every power step",
    xlabel="t [s]", ylabel="TOT [degC]",
)
tot_plot.line(history.times, tot_series, label="TOT")
tot_plot.line(history.times, [TOT_SETPOINT - 273.15] * len(history.times),
              label="TOT setpoint", linestyle="--", marker=None, color="#000000")
with_power_set_axis(tot_plot)
tot_plot.finish(show=False, save_path="t100_transient_tot.png")
print("Saved t100_transient_tot.png")

casing_plot = ThermoPlot(
    title="Casing and shaft thermal masses soaking behind the load steps",
    xlabel="t [s]", ylabel="temperature [degC]",
)
casing_plot.line(history.times,
                 [T - 273.15 for T in history.diff_history["turb_casing.T"]],
                 label="turbine casing")
casing_plot.line(history.times,
                 [T - 273.15 for T in history.diff_history["shaft_mass.T"]],
                 label="shaft")
casing_plot.line(history.times,
                 [T - 273.15 for T in history.diff_history["comp_casing.T"]],
                 label="compressor casing")
with_power_set_axis(casing_plot)
casing_plot.finish(show=False, save_path="t100_transient_casings.png")
print("Saved t100_transient_casings.png")

# =============================================================================
# Wrap-up
# =============================================================================

banner("Summary")
print(f"""Model: T100 recuperated microturbine — shipped maps, Cantera combustion,
casing/shaft thermal masses (two-way coupled), two-loop control.
Calibration constants (set once): turbine map E_fact x {K_TURB},
generator/PE efficiency {ETA_GEN}, recuperator effectiveness {RECUP_EFFECTIVENESS}.

Steady-state validation vs. {len(results)} measured rig points (35-70 kW):""")
for i, quantity in enumerate(("speed", "COT", "COP")):
    errors = [abs(p["err"][i]) for p in results]
    print(f"  {quantity:>5}: mean abs error {sum(errors) / len(errors):4.2f} %")
print("""
Transient: TOT held by the fuel PID while the power PID chased five
dispatch steps; compressor speed/airflow lag each step by seconds, the
thermal masses by minutes.

Figures: t100_validation_{speed,cot_cop,error}.png,
         t100_transient_{power,lag,tot,casings}.png""")
