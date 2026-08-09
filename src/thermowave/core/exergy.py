"""Exergy analysis: a SolveResult-consuming derived view, the same role
core/reporting.py plays for the plain P/T/h/mdot report -- not a solver
change, not a new component, just a way to look at an already-converged
network's own state.

WHAT THIS COVERS
-----------------
- Per-node physical exergy (e_PH, split into thermal e_T and mechanical e_M
  components, the standard SPECO decomposition), for any node
  whose fluid exposes entropy_ph (CoolPropFluid, IdealGasFluid,
  IdealGasMixtureFluid -- see fluids/base_fluid.py's ConstantCpFluid and
  fluids/two_phase.py's supports_entropy()). CanteraFluid nodes raise a
  clear ValueError instead of a silently wrong number -- exergy needs
  entropy, which CanteraFluid doesn't provide here.
- Per-component fuel/product/destruction (E_F, E_P, E_D) costing for the
  component types where "what's fuel and what's product" is unambiguous
  from the component's own physics: turbines/compressors/pumps (work
  producers/consumers) and the two-stream heat exchangers, Condenser
  included (its wf_in/wf_out side is the hot/fuel side, cool_in/cool_out
  the cold/product side). Dispatched by
  isinstance, not report_category() -- see _cost_component()'s own
  docstring for why (Turbine and Compressor share report_category() ==
  "turbomachinery" but have opposite fuel/product roles report_category()
  can't distinguish on its own).

WHAT THIS DELIBERATELY DOES NOT COVER (see exergy_report()'s docstring for
the caller-facing summary of each):
- Combustor/SimpleCombustor: no fuel chemical-exergy/LHV database exists in
  this package, so these get connection-level node exergy at their in/out
  nodes like everything else, but never a component-level E_F/E_P/E_D row.
  This is a deliberate scope boundary, not an oversight.
- Pipe, Valve, Junction, Source, Sink, Sensor, Setpoint, Controller: these
  fix or read boundaries rather than doing anything a fuel/product split
  would describe -- connection-level exergy only.
- Shaft, ShaftLoad, ElectricMotor, Generator, SimpleGenerator: no flow
  ports at all, so no connection exergy either; not separately costed
  since their story is already fully captured by whatever turbomachinery
  component's "power [W]" they read from (costing them too would double
  count the same power).
- SimpleEvaporator/SimpleCondenser (single-stream: no secondary stream is
  modeled at all) can only be costed if the caller supplies the external
  heat source/sink's own temperature via `source_temperatures`; omitted by
  default (not costed, connection-level exergy only) rather than guessed.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import TYPE_CHECKING

from thermowave.components.compressor import Compressor
from thermowave.components.condenser import Condenser
from thermowave.components.electric_motor import ElectricMotor
from thermowave.components.generator import Generator
from thermowave.components.heat_exchanger import HeatExchanger, MultiPassHeatExchanger
from thermowave.components.pump import Pump
from thermowave.components.shaft import Shaft
from thermowave.components.shaft_load import ShaftLoad
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.simple_generator import SimpleGenerator
from thermowave.components.simple_turbine import SimpleTurbine
from thermowave.components.steam_turbine import SteamTurbine
from thermowave.components.turbine import Turbine
from thermowave.core.network import NetworkState, TargetValue
from thermowave.core.reporting import render_metrics_table, table_row
from thermowave.fluids.two_phase import supports_entropy

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.solver import SolveResult
    from thermowave.fluids.base_fluid import BaseFluid

_EXPANDERS = (Turbine, SimpleTurbine, SteamTurbine)  # work OUT is the product
_COMPRESSION = (Compressor, SimpleCompressor, Pump)  # work IN is the fuel
_TWO_STREAM_HX = (HeatExchanger, MultiPassHeatExchanger, Condenser)  # MultiPassHeatExchanger
# is a HeatExchanger subclass, so isinstance already covers it -- listed explicitly for
# clarity. The single-stream SimpleHeatExchanger is NOT here: it has no two streams to
# split into a hot/fuel and a cold/product side at all.
#
# Condenser IS a two-stream exchanger and costs identically (fuel = the hot side's exergy
# drop, product = the cold side's exergy rise); it just names its ports after its two
# streams' roles rather than their temperatures. Costing it is not optional cosmetics: it
# was previously falling through every branch here silently -- no E_F/E_P/E_D row, no
# contribution to the network totals, no warning -- which in the SEGS benchmark quietly
# dropped the plant's single largest heat-rejection component (58.3 MW) out of the balance.
_DEFAULT_TWO_STREAM_PORTS = ("hot_in", "hot_out", "cold_in", "cold_out")
_TWO_STREAM_PORTS: dict[type, tuple[str, str, str, str]] = {
    Condenser: ("wf_in", "wf_out", "cool_in", "cool_out"),
}
_POWER_ONLY = (Shaft, ShaftLoad, ElectricMotor, Generator, SimpleGenerator)
_SINGLE_STREAM_PHASE_CHANGE_PORTS = ("in", "out")  # SimpleEvaporator/SimpleCondenser


@dataclass
class NodeExergy:
    """Specific (e_*, J/kg) and flow-rate (E_*, W) physical exergy at one
    node, split into thermal (e_T) and mechanical (e_M) components:
    e_PH = e_T + e_M, the standard SPECO thermal/mechanical decomposition.
    E_* is None where the node's mdot isn't part of the solve
    (mirrors Sensor.report_metrics()'s own "omit mdot [kg/s] if unavailable"
    behavior, not a separate special case invented here).
    """

    e_PH: float
    e_T: float
    e_M: float
    E_PH: float | None
    E_T: float | None
    E_M: float | None


@dataclass
class ComponentExergyCost:
    """E_F (fuel), E_P (product), E_D (destroyed) for one component. E_P is
    None for a component marked `dissipative` (its whole E_F is instead
    tallied into the network-level E_L, not E_D -- see exergy_report()'s
    docstring). epsilon (exergetic efficiency, E_P/E_F) is None wherever
    E_P is None or E_F is ~0.
    """

    E_F: float
    E_P: float | None
    E_D: float | None
    epsilon: float | None


def _require_entropy_for_exergy(fluid: "BaseFluid", node_name: str) -> None:
    if not supports_entropy(fluid):
        raise ValueError(
            f"Cannot compute exergy at node {node_name!r}: fluid "
            f"{type(fluid).__name__} {getattr(fluid, 'name', '?')!r} has no "
            f"entropy_ph/enthalpy_ps. Exergy analysis needs an entropy-capable "
            f"fluid (CoolPropFluid, IdealGasFluid, or IdealGasMixtureFluid) -- "
            f"CanteraFluid isn't supported here yet."
        )


def node_exergy(state: "NetworkState", node: str, T0: float, P0: float) -> NodeExergy:
    """Physical exergy at `node`, split into thermal/mechanical components:

        e_T = (h(T,P0) - h(T0,P0)) - T0*(s(T,P0) - s(T0,P0))   # temperature deviation, at P0
        e_M = (h(T,P)  - h(T,P0))  - T0*(s(T,P)  - s(T,P0))    # pressure deviation, at T
        e_PH = e_T + e_M

    Uses only temperature_ph/enthalpy_pt/entropy_ph -- no fluid-model API
    beyond what ConstantCpFluid/CoolPropFluid already provide. Raises
    ValueError (via _require_entropy_for_exergy) for a fluid without
    entropy_ph (CanteraFluid today).
    """
    P, h = state.node(node)
    fluid = state.fluid_at(node)
    _require_entropy_for_exergy(fluid, node)

    T = fluid.temperature_ph(P, h)
    s = fluid.entropy_ph(P, h)

    h_T0P0 = fluid.enthalpy_pt(P0, T0)
    s_T0P0 = fluid.entropy_ph(P0, h_T0P0)

    h_TP0 = fluid.enthalpy_pt(P0, T)
    s_TP0 = fluid.entropy_ph(P0, h_TP0)

    e_T = (h_TP0 - h_T0P0) - T0 * (s_TP0 - s_T0P0)
    e_M = (h - h_TP0) - T0 * (s - s_TP0)
    e_PH = e_T + e_M

    try:
        mdot = state.mdot(node)
    except KeyError:
        mdot = None

    return NodeExergy(
        e_PH=e_PH, e_T=e_T, e_M=e_M,
        E_PH=None if mdot is None else mdot * e_PH,
        E_T=None if mdot is None else mdot * e_T,
        E_M=None if mdot is None else mdot * e_M,
    )


def _E_PH(state: "NetworkState", node: str, T0: float, P0: float) -> float:
    """A component's own local port name (e.g. "turb.in") always resolves
    correctly via NetworkState.node()/mdot() -- that's the whole point of
    the existing architecture (components address their own ports by their
    own local ids, not some externally-canonicalized name) -- so exergy at
    a given port is always computed fresh here, not looked up in a
    pre-built table keyed by Network's deduplicated node_order (which is
    only a display list, not every valid node name)."""
    E = node_exergy(state, node, T0, P0).E_PH
    if E is None:
        raise ValueError(
            f"Node {node!r} has no mass flow rate available -- cannot cost a "
            f"component whose exergy balance needs E_PH (mdot * e_PH) there."
        )
    return E


def _two_stream_port_names(component: "BaseComponent") -> tuple[str, str, str, str]:
    """(hot_in, hot_out, cold_in, cold_out) port ids for a two-stream
    exchanger, as ITS OWN class names them. The costing itself is identical
    across these types -- fuel is the hot side's exergy drop, product the
    cold side's rise -- so only the four port ids vary, and they're looked
    up here rather than hardcoded at the one call site."""
    for cls, names in _TWO_STREAM_PORTS.items():
        if isinstance(component, cls):
            return names
    return _DEFAULT_TWO_STREAM_PORTS


def _cost_component(
    component: "BaseComponent",
    state: "NetworkState",
    dissipative: frozenset[str],
    source_temperatures: dict[str, float],
    T0: float,
    P0: float,
) -> ComponentExergyCost | None:
    """Dispatches by isinstance, not report_category(): reporting.py's
    string-based categories only ever decide which TABLE/COLUMNS a
    component's metrics render under, where every component sharing a
    category already agrees on report_metrics()'s key *shape*. Exergy
    costing needs component *roles* instead -- Turbine (fuel = exergy
    drop, product = power out) and Compressor (fuel = power in, product =
    exergy rise) are both report_category() == "turbomachinery" but have
    swapped fuel/product roles no string tag distinguishes. A
    report_category()-style opt-in hook every component declares for
    itself is the natural longer-term answer; isinstance dispatch
    contained entirely in this module is enough for the component types
    this pass actually costs, with zero changes needed to any component
    file to add it.
    """
    ports = component.ports()

    if isinstance(component, _EXPANDERS):
        E_in = _E_PH(state, ports["in"], T0, P0)
        E_out = _E_PH(state, ports["out"], T0, P0)
        E_F = E_in - E_out
        E_P = component.report_metrics(state)["power [W]"]
        return _finish_cost(E_F, E_P, component.name)

    if isinstance(component, _COMPRESSION):
        E_F = component.report_metrics(state)["power [W]"]
        E_in = _E_PH(state, ports["in"], T0, P0)
        E_out = _E_PH(state, ports["out"], T0, P0)
        E_P = E_out - E_in
        return _finish_cost(E_F, E_P, component.name)

    if isinstance(component, _TWO_STREAM_HX):
        hot_in, hot_out, cold_in, cold_out = _two_stream_port_names(component)
        E_hot_in = _E_PH(state, ports[hot_in], T0, P0)
        E_hot_out = _E_PH(state, ports[hot_out], T0, P0)
        E_F = E_hot_in - E_hot_out
        if component.name in dissipative:
            return ComponentExergyCost(E_F=E_F, E_P=None, E_D=None, epsilon=None)
        E_cold_in = _E_PH(state, ports[cold_in], T0, P0)
        E_cold_out = _E_PH(state, ports[cold_out], T0, P0)
        E_P = E_cold_out - E_cold_in
        return _finish_cost(E_F, E_P, component.name)

    if isinstance(component, _POWER_ONLY):
        # No flow ports, hence no exergy of their own to speak of -- their
        # power already shows up as the E_F/E_P of whatever turbomachinery
        # component reads/feeds it (see this module's own docstring).
        return None

    if component.name in source_temperatures:
        # Single-stream phase-change (SimpleEvaporator/SimpleCondenser):
        # no secondary stream is modeled, so the external heat source's own
        # exergy needs its temperature explicitly, via source_temperatures.
        # T_source < T0 (a sub-ambient sink) makes (1 - T0/T_source)
        # negative, which is physically correct for e.g. a refrigeration
        # evaporator -- not a bug, so it is deliberately not special-cased.
        if set(ports) != set(_SINGLE_STREAM_PHASE_CHANGE_PORTS):
            raise ValueError(
                f"{component.name!r} is in source_temperatures but isn't a "
                f"single-stream phase-change component (expected ports "
                f"{_SINGLE_STREAM_PHASE_CHANGE_PORTS}, got {sorted(ports)})."
            )
        T_source = source_temperatures[component.name]
        Q = component.report_metrics(state)["power [W]"]
        E_F = abs(Q) * (1.0 - T0 / T_source)
        E_in = _E_PH(state, ports["in"], T0, P0)
        E_out = _E_PH(state, ports["out"], T0, P0)
        E_P = E_out - E_in
        return _finish_cost(E_F, E_P, component.name)

    return None


def _finish_cost(E_F: float, E_P: float, name: str = "") -> ComponentExergyCost:
    E_D = E_F - E_P
    if E_D < 0.0:
        # Exergy destruction is only ever >= 0 by the second law -- a
        # negative value here means the fuel/product split fed in is
        # nonphysical (most likely an inconsistent power-sign convention
        # upstream), not a real result. Surfacing it as a warning rather
        # than raising: this function only *reports* on an
        # already-converged solve, and the caller may still want the rest
        # of the exergy report even with one bad row in it.
        warnings.warn(
            f"{name or 'component'!r}: negative exergy destruction "
            f"(E_D={E_D:.3g} W) is physically impossible -- check upstream "
            "power-sign conventions",
            stacklevel=3,
        )
    epsilon = E_P / E_F if abs(E_F) > 1e-9 else None
    return ComponentExergyCost(E_F=E_F, E_P=E_P, E_D=E_D, epsilon=epsilon)


def _resolve(items: "list[TargetValue]", state: "NetworkState") -> float:
    return sum(v(state) if callable(v) else v for v in items)


@dataclass
class ExergyReport:
    """Mirrors SolveResult's own spirit: plain attributes, not a
    heavyweight class. node_exergy/component_costs are name -> dataclass
    dicts (component_costs values are None for anything not costed --
    see this module's own docstring for the full list of what that
    covers). E_D is the sum of every costed, non-dissipative component's
    own E_D; E_L is the sum of dissipative components' whole E_F,
    reclassified as a loss to the environment rather than an internal
    irreversibility -- the standard exergy destruction/loss
    distinction. epsilon is E_P/E_F, None if E_F is ~0.
    """

    T0: float
    P0: float
    node_exergy: dict[str, NodeExergy]
    component_costs: dict[str, ComponentExergyCost | None]
    E_F: float
    E_P: float
    E_D: float
    E_L: float
    epsilon: float | None

    def print_report(self) -> None:
        print(format_exergy_report(self))


def exergy_report(
    result: "SolveResult",
    T0: float,
    P0: float,
    fuel: "list[TargetValue]",
    product: "list[TargetValue]",
    dissipative: "list[str]" = (),
    source_temperatures: dict[str, float] | None = None,
) -> ExergyReport:
    """Exergy analysis of a converged SolveResult against ambient (T0, P0).

    fuel/product: each entry is a plain float or a callable
    NetworkState -> float (the same TargetValue type Setpoint/Controller's
    own `value` parameter already accepts) -- e.g.
    `fuel=[lambda s: heater.report_metrics(s)["power [W]"]]`. Whether these
    cover the network's whole boundary (so E_F - E_P lines up with
    E_D + E_L) is the caller's responsibility, the same posture Setpoint's
    own docstring takes about closing a free parameter.

    dissipative: component names (must be HeatExchanger or
    MultiPassHeatExchanger instances in this network) whose cold-side gain
    isn't a useful product -- e.g. a cooler dumping heat to ambient cooling
    water. Their E_P is reported as None and their whole E_F is tallied
    into E_L (a loss to the environment) instead of E_D (an internal
    irreversibility).

    source_temperatures: component name -> external heat source/sink
    temperature [K], for costing a single-stream SimpleEvaporator/
    SimpleCondenser (which models no secondary stream, so its own external
    heat source's exergy is otherwise uncomputable). Omitted by default:
    such components then get connection-level node exergy only, no
    component-level cost.

    See this module's own docstring for the full list of what's NOT
    costed at the component level (Combustor/SimpleCombustor, Pipe/Valve/
    Junction/Source/Sink/Sensor/Setpoint/Controller, and the power-only
    bookkeeping components) and why.
    """
    dissipative_set = frozenset(dissipative)
    source_temperatures = source_temperatures or {}

    component_by_name = {c.name: c for c in result.components}
    for name in dissipative_set:
        target = component_by_name.get(name)
        if target is None or not isinstance(target, _TWO_STREAM_HX):
            raise ValueError(
                f"exergy_report(dissipative=...) names {name!r}, but that isn't a "
                f"HeatExchanger/MultiPassHeatExchanger/Condenser component in this network."
            )

    state = result.state()
    node_exergies = {node: node_exergy(state, node, T0, P0) for node in result.node_order}

    component_costs: dict[str, ComponentExergyCost | None] = {}
    E_D_total = 0.0
    E_L_total = 0.0
    for component in result.components:
        cost = _cost_component(
            component, state, dissipative_set, source_temperatures, T0, P0
        )
        component_costs[component.name] = cost
        if cost is None:
            continue
        if component.name in dissipative_set:
            E_L_total += cost.E_F
        else:
            E_D_total += cost.E_D

    E_F_total = _resolve(fuel, state)
    E_P_total = _resolve(product, state)
    epsilon = E_P_total / E_F_total if abs(E_F_total) > 1e-9 else None

    return ExergyReport(
        T0=T0, P0=P0,
        node_exergy=node_exergies,
        component_costs=component_costs,
        E_F=E_F_total, E_P=E_P_total, E_D=E_D_total, E_L=E_L_total,
        epsilon=epsilon,
    )


_NODE_EXERGY_COLUMNS = ["e_PH [J/kg]", "e_T [J/kg]", "e_M [J/kg]", "E_PH [W]"]
_COMPONENT_COST_COLUMNS = ["E_F [W]", "E_P [W]", "E_D [W]", "epsilon [-]"]


def format_exergy_report(report: ExergyReport) -> str:
    """Render an ExergyReport as a per-node exergy table + per-component
    cost table + network totals, reusing core.reporting's own table
    renderer rather than a second, parallel implementation."""
    node_rows = [
        (
            name,
            {
                "e_PH [J/kg]": ex.e_PH, "e_T [J/kg]": ex.e_T, "e_M [J/kg]": ex.e_M,
                "E_PH [W]": ex.E_PH,
            },
        )
        for name, ex in report.node_exergy.items()
    ]
    component_rows = [
        (
            name,
            {
                "E_F [W]": cost.E_F, "E_P [W]": cost.E_P, "E_D [W]": cost.E_D,
                "epsilon [-]": None if cost.epsilon is None else cost.epsilon * 100.0,
            },
        )
        for name, cost in report.component_costs.items()
        if cost is not None
    ]

    sections = [
        "ThermoWave Exergy Analysis",
        f"  T0 = {report.T0:.2f} K, P0 = {report.P0:.2f} Pa",
        "",
        "Connections",
        *render_metrics_table(node_rows, _NODE_EXERGY_COLUMNS),
    ]
    if component_rows:
        sections += ["", "Components", *render_metrics_table(component_rows, _COMPONENT_COST_COLUMNS)]

    epsilon_str = "-" if report.epsilon is None else f"{report.epsilon * 100.0:.1f}"
    sections += [
        "",
        "Network",
        table_row([
            f"{'E_F [W]':>14}", f"{'E_P [W]':>14}", f"{'E_D [W]':>14}",
            f"{'E_L [W]':>14}", f"{'epsilon [%]':>14}",
        ]),
        table_row([
            f"{report.E_F:>14.4g}", f"{report.E_P:>14.4g}", f"{report.E_D:>14.4g}",
            f"{report.E_L:>14.4g}", f"{epsilon_str:>14}",
        ]),
    ]

    rule_w = max(len(line) for line in sections)
    rule = "═" * rule_w
    return "\n".join([rule, *sections, rule])
