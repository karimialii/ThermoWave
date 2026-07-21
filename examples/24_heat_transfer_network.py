"""Heat transfer: ThermalMass, Convection, and Conduction, on a real cycle.

Same fixed-N turboshaft as 19_generator_component.py (Source -> Compressor
-> Pipe(heater) -> Turbine -> Sink), but this time neither machine is
adiabatic: each one's outlet fluid convects heat into its own casing
(a ThermalMass), the two casings conduct through a third ThermalMass
standing in for the shaft between them, and both casings convect to
ambient. This is full two-way coupling, not just bookkeeping: each
machine's own energy residual actually includes the heat it's losing
through heat_path (see Turbine/Compressor's own docstrings), so the
casings' temperatures change what T_out the cycle itself predicts.

thermowave.components.heat_transfer has four pieces:
  - ThermalMass: a solid's own temperature as a differential state
    (dT/dt = Q_net / thermal_capacitance), the node any Convection/
    Conduction/Radiation path attaches to on the solid side.
  - Convection: Q = h*A*(T_a - T_b) -- covers both free and forced
    convection (they differ only in how h is obtained physically).
  - Conduction: Q = (k*A/L)*(T_a - T_b) -- steady 1D conduction.
  - Radiation: Q = emissivity*view_factor*sigma*A*(T_a^4 - T_b^4) --
    general-purpose surface-to-surface physics (see its own docstring for
    why this is the intended primitive a future 1D combustion-chamber
    liner model builds on) -- not used in this example.

Every Q(state) is positive when `a` is hotter than `b` (heat flowing a->b),
the same sign convention Pipe's own heat_loss already uses. Wiring a path
in takes two steps because both endpoints need to know about it: build the
endpoints first, then the path referencing both, then set
<flow_component>.heat_path = path and <thermal_mass>.heat_sources =
[(path, sign)] by plain attribute assignment -- see ThermalMass's own
docstring for why sign is +1.0/-1.0 per source.

Run: .venv/bin/python examples/24_heat_transfer_network.py
"""

from thermowave.components.compressor import Compressor
from thermowave.components.heat_transfer import Conduction, Convection, ThermalMass
from thermowave.components.pipe import Pipe
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.turbine import Turbine
from thermowave.core.network import Network
from thermowave.core.plotting import ThermoPlot
from thermowave.fluids.ideal_gas import IdealGasFluid

air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
gamma = 1005.0 / (1005.0 - 287.05)
T_ambient = 288.15

src = Source(name="src", P=101325.0, T=288.15, mdot=0.63)
comp = Compressor(name="comp", map_path="T100 Comp.cop", gamma=gamma, N=65000.0)
heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=-431000.0)
turb = Turbine(name="turb", map_path="T100 Turb.tur", gamma=gamma, N=65000.0)
snk = Sink(name="snk")

comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=200.0, T0=300.0)
turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=200.0, T0=300.0)
shaft_mass = ThermalMass(name="shaft_mass", thermal_capacitance=50.0, T0=300.0)

conv_comp_fluid = Convection(name="conv_comp_fluid", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3)
conv_turb_fluid = Convection(name="conv_turb_fluid", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3)
cond_turb_shaft = Conduction(name="cond_turb_shaft", a=turb_casing, b=shaft_mass, k=15.0, A=0.01, L=0.5)
cond_shaft_comp = Conduction(name="cond_shaft_comp", a=shaft_mass, b=comp_casing, k=15.0, A=0.01, L=0.5)
conv_turb_ambient = Convection(name="conv_turb_ambient", a=turb_casing, b=T_ambient, h=10.0, A=1.0)
conv_comp_ambient = Convection(name="conv_comp_ambient", a=comp_casing, b=T_ambient, h=10.0, A=1.0)

# Two-step wiring: both components above already exist, so now set the
# two-way back-references plain attribute assignment (see this file's own
# docstring for why this can't all happen in one constructor call).
comp.heat_path = conv_comp_fluid
turb.heat_path = conv_turb_fluid
comp_casing.heat_sources = [(conv_comp_fluid, 1.0), (cond_shaft_comp, 1.0), (conv_comp_ambient, -1.0)]
turb_casing.heat_sources = [(conv_turb_fluid, 1.0), (cond_turb_shaft, -1.0), (conv_turb_ambient, -1.0)]
shaft_mass.heat_sources = [(cond_turb_shaft, 1.0), (cond_shaft_comp, -1.0)]

network = Network(fluid=air)
for component in (
    src, comp, heater, turb, snk,
    comp_casing, turb_casing, shaft_mass,
    conv_comp_fluid, conv_turb_fluid, cond_turb_shaft, cond_shaft_comp,
    conv_turb_ambient, conv_comp_ambient,
):
    network.add_component(component)
network.connect(src, "out", comp, "in")
network.connect(comp, "out", heater, "in")
network.connect(heater, "out", turb, "in")
network.connect(turb, "out", snk, "in")

# Cold start: solve_transient()'s default initial condition is already a
# steady state for every differential state, including these casings
# (Q_net == 0), so start from an explicit cold override instead to see
# them actually warm up.
cold_initial = network.solve(tol=1e-8, max_iter=400, damping=0.3)
for mass_name in ("comp_casing", "turb_casing", "shaft_mass"):
    cold_initial.params[f"{mass_name}.T"] = 300.0

history = network.solve_transient(
    duration=200.0, dt=10.0, tol=1e-8, max_iter=400, damping=0.3, initial=cold_initial,
)

print(f"{'t [s]':>8}  {'turb_casing [K]':>16}  {'shaft_mass [K]':>15}  {'comp_casing [K]':>16}")
for t, step in zip(history.times, history.steps):
    print(
        f"{t:8.1f}  {step.params['turb_casing.T']:16.2f}  "
        f"{step.params['shaft_mass.T']:15.2f}  {step.params['comp_casing.T']:16.2f}"
    )

T_out_final = air.temperature_ph(
    history.steps[-1].node_P["turb.out"], history.steps[-1].node_h["turb.out"]
)
print(f"\nFinal turbine T_out (with heat loss) = {T_out_final:.2f} K")

# Plotting requires the 'plot' extra: pip install thermowave[plot]
ThermoPlot(
    title="Casing/shaft temperatures warming from a cold start",
    xlabel="t [s]", ylabel="T [K]",
).transient(
    history, [
        (turb_casing, "T [K]", "turbine casing"),
        (shaft_mass, "T [K]", "shaft"),
        (comp_casing, "T [K]", "compressor casing"),
    ],
).finish(show=False, save_path="heat_transfer_network.png")
print("Saved heat_transfer_network.png")
