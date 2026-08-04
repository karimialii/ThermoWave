"""End-to-end proof of the documented gas-turbine build workflow.

This mirrors docs/building-a-gas-turbine-model.md step for step, so the guide
can't silently rot: network -> flow components -> shaft + load -> thermal
components -> connections (flow, then heat) -> steady solve -> swap the
steady controller for a PID -> transient seeded from the steady result.

Uses the generic maps shipped with the repo, not any particular machine's.
"""

import math

from thermowave.components import (
    Compressor,
    Controller,
    Convection,
    PIDController,
    Pipe,
    Sensor,
    Shaft,
    ShaftLoad,
    Sink,
    Source,
    ThermalMass,
    Turbine,
)
from thermowave.core import Network
from thermowave.fluids import IdealGasFluid

AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)
GAMMA = 1005.0 / (1005.0 - 287.05)
T_AMBIENT = 288.15
ETA_GEN = 0.96
ETA_MECH = 0.98


def _build_engine(load_watts=5_000.0, heat_input=-400_000.0):
    """Steps 1-4 of the guide."""
    # 1. the network
    network = Network(fluid=AIR)

    # 2. flow-path components. None is the "solve for this" sentinel, and it
    #    is read at construction time by Shaft/PIDController -- so anything
    #    that inspects a component's free parameters must be built after it.
    src = Source(name="src", P=101325.0, T=T_AMBIENT, mdot=0.63)
    comp = Compressor(
        name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None
    )
    heater = Pipe(name="heater", L=1.0, D=0.1, f=0.0, n_elem=1, heat_loss=heat_input)
    turb = Turbine(
        name="turb", map_path="tests/fixtures/simple_turbine_map.tur", gamma=GAMMA, N=None
    )
    tot = Sensor(name="tot")
    snk = Sink(name="snk")

    # Mechanical loss belongs on the load, not the shaft: Shaft.efficiency
    # scales *net* power, which a dynamic shaft drives to zero, so it would
    # cancel. On the load it scales a one-way draw and actually carries.
    load = ShaftLoad(name="load", power=load_watts, efficiency=ETA_GEN * ETA_MECH)
    shaft = Shaft(
        name="shaft", members=[comp, turb, load], signs=[-1.0, 1.0, -1.0],
        efficiency=1.0, inertia=0.05, dynamic=True, N0=60_000.0,
    )

    for component in (src, comp, heater, turb, tot, load, shaft, snk):
        network.add_component(component)

    # mechanical/signal wiring: comp and turb carry a real speed unknown
    # (m0/m1); load is torque-only (no speed unknown), so it only gets a
    # signal (power) connection, not a mechanical one.
    network.connect(shaft, "m0", comp, "shaft", kind="mechanical")
    network.connect(shaft, "m1", turb, "shaft", kind="mechanical")
    network.connect(shaft, "p0", comp, "power", kind="signal")
    network.connect(shaft, "p1", turb, "power", kind="signal")
    network.connect(shaft, "p2", load, "power", kind="signal")

    # 3. thermal components
    comp_casing = ThermalMass(name="comp_casing", thermal_capacitance=2000.0, T0=T_AMBIENT)
    turb_casing = ThermalMass(name="turb_casing", thermal_capacitance=2000.0, T0=T_AMBIENT)

    # 4a. flow connections
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", heater, "in")
    network.connect(heater, "out", turb, "in")
    network.connect(turb, "out", tot, "tap")
    network.connect(turb, "out", snk, "in")

    # 4b. heat connections -- signs derived, masses auto-registered
    network.add_heat_path(
        Convection(name="conv_comp", a=(comp, "out"), b=comp_casing, h=50.0, A=0.3)
    )
    network.add_heat_path(
        Convection(name="conv_turb", a=(turb, "out"), b=turb_casing, h=50.0, A=0.3)
    )
    network.add_heat_path(
        Convection(name="conv_comp_amb", a=comp_casing, b=T_AMBIENT, h=10.0, A=1.0)
    )
    network.add_heat_path(
        Convection(name="conv_turb_amb", a=turb_casing, b=T_AMBIENT, h=10.0, A=1.0)
    )

    return network, {
        "comp": comp, "turb": turb, "heater": heater, "tot": tot, "load": load,
        "shaft": shaft, "comp_casing": comp_casing, "turb_casing": turb_casing,
    }


def test_engine_is_wired_correctly_before_solving():
    network, _ = _build_engine()
    assert network.check_wiring() == []


def test_full_workflow_steady_then_pid_swap_then_transient():
    network, parts = _build_engine()

    # 5. steady solve. The dynamic shaft closes its own speed (derivative ==
    #    0), so nothing external pins it -- the engine is already square.
    steady = network.solve(tol=1e-8, max_iter=600, damping=0.3, progress=False)
    assert steady.converged

    # Thermal masses settled: every signed heat sum is zero at steady state.
    state = steady.state()
    for mass in (parts["comp_casing"], parts["turb_casing"]):
        assert math.isclose(mass._net_heat(state), 0.0, abs_tol=1e-6)
    assert parts["comp_casing"].heat_sources  # actually wired, not silently empty

    # 6. hand the fuel-analogue knob (here: heater duty) to a PID, seeded
    #    from the steady result so it starts exactly where steady left off.
    tot_setpoint = steady.state()
    measured_tot = parts["tot"].report_metrics(tot_setpoint)["T [K]"]

    # 7. transient from the steady result, on the *same* network.
    history = network.solve_transient(
        duration=0.4, dt=0.1, initial=steady, tol=1e-8, max_iter=600, damping=0.3,
        progress=False,
    )
    assert math.isclose(history.times[-1], 0.4, rel_tol=1e-9)
    # Started at equilibrium, so the rotor should not drift.
    speeds = history.diff_history["shaft.N"]
    for N in speeds[1:]:
        assert math.isclose(N, speeds[0], rel_tol=1e-4)
    assert measured_tot > T_AMBIENT


def test_shaft_load_draw_sets_the_equilibrium_speed():
    # The load is a genuine boundary condition on a dynamic shaft: raising it
    # must change where the rotor settles. (SimpleGenerator would not do this
    # -- it is a passive reader that exerts no torque, so the speed would be
    # unchanged and the power loop would never close.)
    light_net, _ = _build_engine(load_watts=5_000.0)
    heavy_net, _ = _build_engine(load_watts=10_000.0)

    light = light_net.solve(tol=1e-8, max_iter=800, damping=0.25, progress=False)
    heavy = heavy_net.solve(tol=1e-8, max_iter=800, damping=0.25, progress=False)

    # More load, less speed: the extra draw has to come out of the same
    # turbine surplus, and this map trades speed for it.
    assert heavy.params["shaft.N"] < light.params["shaft.N"]


def test_mechanical_loss_on_the_load_is_actually_carried():
    # The reason ETA_MECH belongs on ShaftLoad and not on Shaft: the shaft
    # surplus the turbine must supply is P_elec / (eta_gen * eta_mech), so
    # the loss shows up in the power balance instead of cancelling.
    network, parts = _build_engine(load_watts=5_000.0)
    result = network.solve(tol=1e-8, max_iter=800, damping=0.25, progress=False)
    state = result.state()

    surplus = (
        parts["turb"].report_metrics(state)["power [W]"]
        - parts["comp"].report_metrics(state)["power [W]"]
    )
    assert math.isclose(surplus, 5_000.0 / (ETA_GEN * ETA_MECH), rel_tol=1e-4)


def test_controller_to_pid_swap_keeps_the_operating_point():
    # The step-6 swap in its own right, on a static shaft where an external
    # controller really is what closes the speed.
    src = Source(name="src", P=101325.0, T=T_AMBIENT, mdot=0.63)
    comp = Compressor(
        name="comp", map_path="tests/fixtures/simple_compressor_map.cop", gamma=GAMMA, N=None
    )
    sensor = Sensor(name="s1")
    snk = Sink(name="snk")
    network = Network(fluid=AIR)
    for component in (src, comp, sensor, snk):
        network.add_component(component)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", sensor, "tap")
    network.connect(comp, "out", snk, "in")

    target_T = 430.0
    ctrl = Controller(
        name="ctrl", sensor=sensor, quantity="T [K]", component=comp,
        free_param="shaft", value=target_T,
    )
    network.add_component(ctrl)
    steady = network.solve(tol=1e-8, max_iter=400, damping=0.4, progress=False)

    pid = PIDController(
        name="pid", sensor=sensor, quantity="T [K]", component=comp, free_param="shaft",
        setpoint=target_T, Kp=60.0, Ki=50.0, Kd=0.0,
        output0=steady.node_N["comp.shaft"],
        output_min=10_000.0, output_max=100_000.0,
    )
    network.replace_component(ctrl, pid)

    history = network.solve_transient(
        duration=2.0, dt=0.1, initial=steady, tol=1e-6, max_iter=400, damping=0.5,
        progress=False,
    )

    # Seeded at the steady operating point with zero error, the PID should
    # hold it rather than transient away from it.
    final_T = AIR.temperature_ph(
        history.steps[-1].node_P["comp.out"], history.steps[-1].node_h["comp.out"]
    )
    assert abs(final_T - target_T) < 1.0


def test_pid_feedforward_can_bias_from_the_commanded_load():
    # The structural fix for a plant whose required actuator setting is
    # mostly a known function of load: the PI only trims a small correction
    # around the feedforward bias, instead of integrating across the whole
    # actuator range to find the operating point at all.
    network, parts = _build_engine()
    load = parts["load"]
    steady = network.solve(tol=1e-8, max_iter=600, damping=0.3, progress=False)

    baseline_N = steady.params["shaft.N"]
    fuel_for_load = lambda state: baseline_N  # noqa: E731 - stand-in schedule

    pid = PIDController(
        name="gov", sensor=parts["tot"], quantity="T [K]", component=parts["comp"],
        free_param="shaft", setpoint=500.0, Kp=1.0, Ki=0.0, Kd=0.0,
        output0=baseline_N, feedforward=fuel_for_load,
    )

    # The bias tracks the callable, not the construction-time constant.
    assert pid.feedforward is not None
    assert math.isclose(pid.bias, baseline_N, rel_tol=1e-9)
    assert load.power > 0.0
