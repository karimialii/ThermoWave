from __future__ import annotations

from typing import TYPE_CHECKING, Callable, Union

import networkx as nx

from thermowave.core.exceptions import NetworkTopologyError

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.solver import SolveResult
    from thermowave.core.transient import TransientResult
    from thermowave.fluids.base_fluid import BaseFluid

# Connection kinds this Network can wire. "flow" merges two ports into one
# shared (P, h, mdot) node; "mechanical" merges two ports into one shared
# shaft-speed (N) node; "signal" merges two ports into one shared scalar
# (e.g. shaft power) resolved from whichever connected component computes
# it. Heat-transfer connections use their own mechanism (see
# Network.add_heat_path()) rather than this list.
_SUPPORTED_CONNECTION_KINDS = {"flow", "mechanical", "signal"}

# A plain float (a fixed design target), or a callable NetworkState -> float
# re-evaluated fresh every residual call, for tying a target to something
# that's itself live/computed rather than a constant -- see Setpoint's own
# docstring (components/setpoint.py) for the original rationale and an
# example. Lives here rather than in setpoint.py so exergy.py's network-level
# fuel/product specification (core/exergy.py) can reuse the exact same type
# without components/* importing from core/exergy.py or vice versa -- both
# setpoint.py and exergy.py depend downward on this module instead.
TargetValue = Union[float, Callable[["NetworkState"], float]]


class NetworkState:
    """Current node/branch state passed to each component's residuals() call."""

    def __init__(
        self,
        fluid: "BaseFluid",
        node_P: dict[str, float],
        node_h: dict[str, float],
        node_mdot: dict[str, float],
        params: dict[str, float] | None = None,
        node_fluid: dict[str, "BaseFluid"] | None = None,
        node_N: dict[str, float] | None = None,
        node_signal: dict[str, float] | None = None,
    ):
        self.fluid = fluid
        self.node_P = node_P
        self.node_h = node_h
        self.node_mdot = node_mdot
        self.params = params if params is not None else {}
        # Shaft speed [rev/min] at every mechanical port this network has
        # (raw component-local port id -> value, already expanded from the
        # canonical/shared value the same way node_P/node_h are) -- see
        # BaseComponent.fixed_mechanical_values()/free_mechanical_ports() and
        # Network._all_mechanical_nodes(). Empty for a network with no
        # mechanical ports at all (the common case today).
        self.node_N = node_N if node_N is not None else {}
        # Scalar signal values (e.g. shaft power [W]) at every signal port
        # some component provides -- raw port id -> value, resolved once per
        # residual evaluation by Network._resolve_node_signal(). Empty for a
        # network with no signal ports.
        self.node_signal = node_signal if node_signal is not None else {}
        # node -> BaseFluid, populated by Network._resolve_node_fluid() for
        # any node downstream of a component whose outlet_fluid() changes
        # composition (e.g. Combustor) -- nodes not present here (the
        # overwhelming majority) fall back to `fluid` via fluid_at() below.
        # Left as a plain public dict (like node_P/node_h/params) rather than
        # always-populated, since resolving it is only worth doing for
        # networks that actually have a composition-changing component.
        self.node_fluid = node_fluid if node_fluid is not None else {}
        # Component-owned scratch space for memoizing an expensive
        # intermediate that more than one hook on the same component needs
        # within one residual evaluation (e.g. Combustor's outlet_fluid()
        # and residuals() both need the same Cantera equilibrium result) --
        # keyed however the owning component likes (e.g. (component.name,
        # "mixture")). A fresh NetworkState is constructed once per Newton
        # iteration/finite-difference perturbation (see Solver.solve()), so
        # this cache is exactly as fresh as it needs to be: never stale
        # within one residual assembly pass, never reused across one.
        self._cache: dict = {}

    def node(self, name: str) -> tuple[float, float]:
        return self.node_P[name], self.node_h[name]

    def mdot(self, name: str) -> float:
        return self.node_mdot[name]

    def param(self, name: str) -> float:
        return self.params[name]

    def N(self, name: str) -> float:
        """Shaft speed [rev/min] at this component's own raw mechanical port
        id `name` -- shares a value with every other port connected to it
        via kind="mechanical", the same way node()/mdot() do for flow ports.
        """
        return self.node_N[name]

    def signal(self, name: str) -> float:
        """Scalar value (e.g. power [W]) at this component's own raw signal
        port id `name`, as computed by whichever connected component
        provides it (see BaseComponent.provided_signal_values()). Raises
        KeyError if nothing in the network provides a value for this port's
        canonical node.
        """
        return self.node_signal[name]

    def fluid_at(self, name: str) -> "BaseFluid":
        """BaseFluid that's actually flowing through node `name` -- the
        network's own default `fluid` unless a component upstream of this
        node changed composition (see BaseComponent.outlet_fluid()) and that
        change was resolved into node_fluid. Prefer this over reading
        `.fluid` directly in any component whose physics reads fluid
        properties at a specific node, so it automatically picks up whatever
        composition actually reaches that node instead of always assuming
        the network's single default fluid.
        """
        return self.node_fluid.get(name, self.fluid)


class Connection:
    """A typed link between one component's port and another's.

    All three kinds merge the two ports into one shared node via the same
    union-find (Network._union()); what differs is which NetworkState
    channel that shared node lives in:
      - "flow": (P, h, mdot) -- the two components are physically joined by
        a stream.
      - "mechanical": N (shaft speed) -- the two components turn together.
      - "signal": a scalar (e.g. shaft power) one connected component
        computes and the other reads (see BaseComponent.
        provided_signal_values()).
    Heat-transfer connections use their own mechanism instead
    (Network.add_heat_path()), since a heat path couples its endpoints
    through their own residuals rather than by merging node state.
    """

    def __init__(
        self,
        from_component: "BaseComponent",
        from_port: str,
        to_component: "BaseComponent",
        to_port: str,
        kind: str = "flow",
    ):
        self.from_component = from_component
        self.from_port = from_port
        self.to_component = to_component
        self.to_port = to_port
        self.kind = kind


def _can_change_composition(component: "BaseComponent") -> bool:
    """Whether this component can alter fluid composition between its inlet
    and outlet — i.e. whether it overrides either composition hook.

    Only Combustor (outlet_fluid) and Junction (merge_fluids) do so among the
    package's components. Comparing the *class* attribute against
    BaseComponent's is the correct test (bound methods would compare unequal
    for every instance); the vars() check additionally catches an override
    monkeypatched onto a single instance.
    """
    from thermowave.components.base_component import BaseComponent

    cls = type(component)
    if cls.outlet_fluid is not BaseComponent.outlet_fluid:
        return True
    if cls.merge_fluids is not BaseComponent.merge_fluids:
        return True
    instance_attrs = vars(component)
    return "outlet_fluid" in instance_attrs or "merge_fluids" in instance_attrs


class _FluidPropagationPlan:
    """Precomputed, topology-only structure for Network._resolve_node_fluid().

    Everything the fluid fixed-point loop needs that does *not* depend on the
    Newton state — each component's ports, its canonical inlet/outlet ids, and
    whether it participates as a merger — resolved once per topology version
    instead of on every one of the ~hundreds-of-thousands of residual
    evaluations a transient run performs. Only the actual outlet_fluid()/
    merge_fluids() calls stay in the per-state loop.
    """

    __slots__ = (
        "has_composition_change",
        "seed_canon",
        "entries",
        "export",
    )

    def __init__(
        self,
        has_composition_change: bool,
        seed_canon: tuple[str, ...],
        entries: list,
        export: tuple[tuple[str, str], ...],
    ):
        self.has_composition_change = has_composition_change
        self.seed_canon = seed_canon
        self.entries = entries
        self.export = export


class _PlanEntry:
    """One component's frozen role in the fluid propagation plan."""

    __slots__ = ("component", "is_merger", "merge_read", "merge_write", "pairs")

    def __init__(self, component, is_merger, merge_read, merge_write, pairs):
        self.component = component
        self.is_merger = is_merger
        # (port_name, canon) for every port, used to assemble merge_fluids'
        # inlet dict; merge_write maps a merged output's port_name -> canon.
        self.merge_read = merge_read
        self.merge_write = merge_write
        # (inlet_port, outlet_port, inlet_canon, outlet_canon) with the
        # "port not in ports" cases already filtered out at build time.
        self.pairs = pairs


class Network:
    """A graph of components sharing a single working fluid and mass flow rate.

    Components are added with add_component(), then wired together with
    connect(from_component, from_port, to_component, to_port), which links one
    component's outlet-like port to another's inlet-like port.
    """

    def __init__(self, fluid: "BaseFluid"):
        self.fluid = fluid
        self.components: list["BaseComponent"] = []
        self.connections: list[Connection] = []
        self.graph = nx.DiGraph()
        self._parent: dict[str, str] = {}  # union-find: port id -> port id
        # Caches keyed on _topology_version (bumped by _bump_topology whenever
        # the port set or any union-find root can change). Both are recomputed
        # lazily and cost <0.01s to rebuild across a whole transient run, but
        # save the O(passes x components x ports) union-find traversal that
        # _resolve_node_fluid otherwise repeats on *every* residual evaluation.
        self._topology_version = 0
        self._canon_cache: dict[str, str] = {}
        self._fluid_plan: "_FluidPropagationPlan | None" = None

    def _bump_topology(self) -> None:
        """Invalidate topology-derived caches after a structural change."""
        self._topology_version += 1
        self._canon_cache.clear()
        self._fluid_plan = None

    def invalidate_caches(self) -> None:
        """Public escape hatch: call after mutating a component *in place*
        (e.g. changing a Junction's inlet count, or monkeypatching
        outlet_fluid onto an instance) once it has already been added. The
        ordinary add_component/connect/remove_component paths invalidate
        automatically; direct attribute mutation on an already-registered
        component does not, since Network can't observe it."""
        self._bump_topology()

    def add_component(self, component: "BaseComponent") -> None:
        self.components.append(component)
        component._network = self
        self._register_nodes(component)

    def _register_nodes(self, component: "BaseComponent") -> None:
        for port_id in component.ports().values():
            self._parent.setdefault(port_id, port_id)
            self.graph.add_node(port_id)
        for port_id in component.mechanical_ports().values():
            self._parent.setdefault(port_id, port_id)
            self.graph.add_node(port_id)
        for port_id in component.signal_ports().values():
            self._parent.setdefault(port_id, port_id)
            self.graph.add_node(port_id)
        for node_name in component.internal_nodes():
            self.graph.add_node(node_name)
        self._bump_topology()

    def add_heat_path(self, path: "BaseComponent") -> "BaseComponent":
        """Wire a Convection/Conduction/Radiation into this network, both
        endpoints included — the heat-transfer counterpart to connect().

        A heat path couples its two endpoints (`path.a` and `path.b`) through
        their own residuals rather than by merging node state, so it isn't a
        connect() kind: there's no shared (P, h) node to union. What it does
        need is a back-reference on each endpoint, which is what this does:

          - a ThermalMass endpoint gets (path, sign) appended to its
            heat_sources, so its dT/dt picks the path up;
          - a (component, port) endpoint gets the path appended to that
            component's heat_path, so its energy residual picks it up;
          - a plain float endpoint (ambient) needs nothing.

        The sign is derived from which endpoint each side actually is, never
        asked for. Q is positive when `a` is hotter than `b` (heat flowing
        a -> b), so:

          - a ThermalMass sums heat *in*: +1.0 as `b` (gains), -1.0 as `a`.
          - a flow component reports heat *lost*: +1.0 as `a` (loses),
            -1.0 as `b` (gains).

        Those two are mirror images of each other, which is exactly the kind
        of bookkeeping that's easy to get backwards by hand — and a wrong
        sign doesn't fail, it silently creates or destroys energy.

        Also registers the path itself (so it appears in the Heat Transfer
        Paths report table) and any ThermalMass endpoint not yet added — a
        mass that never reaches add_component() contributes no differential
        state, and the first thing to notice is a bare KeyError on its own
        temperature deep inside a residual call.

        Idempotent: adding the same path twice is a no-op rather than
        double-counting its Q.
        """
        for attr in ("a", "b"):
            if not hasattr(path, attr):
                raise NetworkTopologyError(
                    f"add_heat_path() expects a heat path with .a/.b endpoints "
                    f"(Convection, Conduction, Radiation), but "
                    f"{getattr(path, 'name', path)!r} has no {attr!r}."
                )
        if not callable(getattr(path, "Q", None)):
            raise NetworkTopologyError(
                f"add_heat_path() expects a heat path exposing Q(state), but "
                f"{getattr(path, 'name', path)!r} doesn't."
            )

        if path in self.components:
            return path  # already wired; re-adding must not double-count Q

        # Endpoint `a` loses Q, endpoint `b` gains it. mass_sign sums heat
        # *into* the mass; loss_sign reports heat *out of* the fluid.
        for endpoint, mass_sign, loss_sign in ((path.a, -1.0, 1.0), (path.b, 1.0, -1.0)):
            self._attach_heat_endpoint(path, endpoint, mass_sign, loss_sign)

        self.add_component(path)
        return path

    def _attach_heat_endpoint(
        self, path: "BaseComponent", endpoint: object, mass_sign: float, loss_sign: float
    ) -> None:
        from thermowave.components.heat_transfer import ThermalMass, normalized_heat_paths

        if isinstance(endpoint, ThermalMass):
            if endpoint not in self.components:
                self.add_component(endpoint)
            if not any(existing is path for existing, _ in endpoint.heat_sources):
                endpoint.heat_sources.append((path, mass_sign))
            return

        if isinstance(endpoint, (int, float)):
            return  # a fixed ambient temperature has nothing to wire back

        try:
            component, port_name = endpoint
        except (TypeError, ValueError) as exc:
            raise NetworkTopologyError(
                f"Heat path {getattr(path, 'name', path)!r} has an endpoint that isn't a "
                f"ThermalMass, a number (ambient), or a (component, port) pair: "
                f"{endpoint!r}"
            ) from exc

        # Raises a named NetworkTopologyError now rather than a bare KeyError
        # later, once this path is already deep inside a residual evaluation.
        self._resolve_port(component, port_name)

        if not hasattr(component, "heat_path"):
            raise NetworkTopologyError(
                f"Heat path {getattr(path, 'name', path)!r} attaches to "
                f"{component.name!r}.{port_name}, but {component.name!r} has no "
                f"heat_path attribute to hang it on — only components whose energy "
                f"residual reads one (Compressor, Turbine, Combustor and their Simple* "
                f"counterparts) can sit on a heat path. Pipe takes a scalar heat_loss "
                f"instead."
            )

        paths = normalized_heat_paths(component.heat_path)
        if not any(existing is path for existing, _ in paths):
            paths.append((path, loss_sign))
        component.heat_path = paths

    def remove_component(self, component: "BaseComponent") -> None:
        """Remove a component, its connections, and any heat-path references
        to it — the inverse of add_component().

        The main use is swapping what closes a free parameter between solves:
        a steady solve pins a compressor's speed with a Controller, then the
        transient run wants a PIDController on the same unknown. Both pin it
        with one residual each, so the Controller has to go first or the
        system is over-determined by exactly one equation.

        Rebuilds the port graph and union-find from the survivors rather than
        unpicking them in place: _parent is a flat child -> parent map with
        path compression and no reverse index, so removing a port id that
        happens to be a union root can't be undone incrementally. Everything
        else the solver reads (_all_nodes, _fixed_node_values, ...) is already
        recomputed from self.components on every call, so there's nothing
        else cached to unwind.
        """
        if component not in self.components:
            raise NetworkTopologyError(
                f"Component {getattr(component, 'name', component)!r} is not in this "
                f"network, so it can't be removed."
            )

        self.components = [c for c in self.components if c is not component]
        self.connections = [
            c
            for c in self.connections
            if c.from_component is not component and c.to_component is not component
        ]
        self._drop_heat_references(component)
        self._rebuild_topology()

    def replace_component(
        self, old: "BaseComponent", new: "BaseComponent"
    ) -> "BaseComponent":
        """remove_component(old) then add_component(new). See
        remove_component() for why the Controller -> PIDController swap needs
        this rather than just adding the new one."""
        self.remove_component(old)
        self.add_component(new)
        return new

    def _drop_heat_references(self, component: "BaseComponent") -> None:
        """Forget a removed component wherever a heat path or thermal mass
        still points at it, so nothing calls Q() on a path that's no longer
        part of the network (or on a path whose endpoint just left)."""
        from thermowave.components.heat_transfer import ThermalMass, normalized_heat_paths

        def endpoint_is_gone(endpoint: object) -> bool:
            if endpoint is component:
                return True
            if isinstance(endpoint, (int, float)) or isinstance(endpoint, ThermalMass):
                return False
            if isinstance(endpoint, tuple) and len(endpoint) == 2:
                return endpoint[0] is component
            return False

        def path_is_gone(path: "BaseComponent") -> bool:
            if path is component:
                return True
            return endpoint_is_gone(getattr(path, "a", None)) or endpoint_is_gone(
                getattr(path, "b", None)
            )

        for remaining in self.components:
            if isinstance(remaining, ThermalMass):
                remaining.heat_sources = [
                    (path, sign) for path, sign in remaining.heat_sources if not path_is_gone(path)
                ]
            if hasattr(remaining, "heat_path") and remaining.heat_path is not None:
                kept = [
                    (path, sign)
                    for path, sign in normalized_heat_paths(remaining.heat_path)
                    if not path_is_gone(path)
                ]
                remaining.heat_path = kept or None

    def _rebuild_topology(self) -> None:
        """Rebuild graph + union-find from self.components/self.connections."""
        self.graph = nx.DiGraph()
        self._parent = {}
        self._bump_topology()
        for component in self.components:
            self._register_nodes(component)
        for connection in self.connections:
            from_id = self._resolve_port(
                connection.from_component, connection.from_port, connection.kind
            )
            to_id = self._resolve_port(
                connection.to_component, connection.to_port, connection.kind
            )
            self._union(from_id, to_id)
            self.graph.add_edge(from_id, to_id, kind=connection.kind)

    def connect(
        self,
        from_component: "BaseComponent",
        from_port: str,
        to_component: "BaseComponent",
        to_port: str,
        kind: str = "flow",
    ) -> Connection:
        if kind not in _SUPPORTED_CONNECTION_KINDS:
            raise NotImplementedError(
                f"Connection kind {kind!r} is not yet supported "
                f"(supported: {sorted(_SUPPORTED_CONNECTION_KINDS)})"
            )
        from_id = self._resolve_port(from_component, from_port, kind)
        to_id = self._resolve_port(to_component, to_port, kind)

        self._union(from_id, to_id)
        self.graph.add_edge(from_id, to_id, kind=kind)

        connection = Connection(from_component, from_port, to_component, to_port, kind)
        self.connections.append(connection)
        return connection

    def _ports_for_kind(self, component: "BaseComponent", kind: str) -> dict[str, str]:
        if kind == "mechanical":
            return component.mechanical_ports()
        if kind == "signal":
            return component.signal_ports()
        return component.ports()

    def _resolve_port(
        self, component: "BaseComponent", port_name: str, kind: str = "flow"
    ) -> str:
        ports = self._ports_for_kind(component, kind)
        kind_label = "" if kind == "flow" else f"{kind} "
        if port_name not in ports:
            raise NetworkTopologyError(
                f"Component {component.name!r} has no {kind_label}port {port_name!r}; "
                f"available {kind_label}ports: {sorted(ports)}"
            )
        return ports[port_name]

    def _find(self, port_id: str) -> str:
        self._parent.setdefault(port_id, port_id)
        root = port_id
        while self._parent[root] != root:
            root = self._parent[root]
        while self._parent[port_id] != root:
            self._parent[port_id], port_id = root, self._parent[port_id]
        return root

    def _union(self, a: str, b: str) -> None:
        root_a, root_b = self._find(a), self._find(b)
        if root_a != root_b:
            self._parent[root_b] = root_a
        self._bump_topology()

    def _canonical(self, port_id: str) -> str:
        """This component port's shared network-node id after connect() merges.

        Memoized on _topology_version: within one topology the canonical root
        of a port never changes, and this is called O(millions) of times per
        solve. _find is still the sole path-compressing mutator; the cache
        only remembers its answer.
        """
        cached = self._canon_cache.get(port_id)
        if cached is not None:
            return cached
        root = self._find(port_id)
        self._canon_cache[port_id] = root
        return root

    def _all_nodes(self) -> list[str]:
        nodes: list[str] = []
        for component in self.components:
            for node_name in component.internal_nodes():
                if node_name not in nodes:
                    nodes.append(node_name)
            for port_id in component.ports().values():
                canonical = self._canonical(port_id)
                if canonical not in nodes:
                    nodes.append(canonical)
        return nodes

    def _port_nodes(self) -> list[str]:
        """_all_nodes() minus every component's solver-internal nodes.

        Used for reporting: internal discretization nodes (e.g. a multi-element
        Pipe's mid-nodes) are still solved and readable via SolveResult.node_P/
        node_h, just not shown in the per-node results table.
        """
        internal = {
            node_name
            for component in self.components
            for node_name in component.internal_nodes()
        }
        return [n for n in self._all_nodes() if n not in internal]

    def _fixed_node_values(self) -> dict[str, tuple[float, float]]:
        fixed: dict[str, tuple[float, float]] = {}
        for component in self.components:
            for port_id, value in component.fixed_node_values(self.fluid).items():
                fixed[self._canonical(port_id)] = value
        return fixed

    def _fixed_node_mdot(self) -> dict[str, float]:
        fixed: dict[str, float] = {}
        for component in self.components:
            for port_id, value in component.fixed_node_mdot().items():
                fixed[self._canonical(port_id)] = value
        return fixed

    def _guess_node_mdot(self) -> dict[str, float]:
        guesses: dict[str, float] = {}
        for component in self.components:
            for port_id, value in component.guess_node_mdot().items():
                guesses[self._canonical(port_id)] = value
        return guesses

    def _fixed_mechanical_values(self) -> dict[str, float]:
        fixed: dict[str, float] = {}
        for component in self.components:
            for port_id, value in component.fixed_mechanical_values().items():
                fixed[self._canonical(port_id)] = value
        return fixed

    def _guess_mechanical_ports(self) -> dict[str, float]:
        guesses: dict[str, float] = {}
        for component in self.components:
            for port_id, value in component.free_mechanical_ports().items():
                canon = self._canonical(port_id)
                guesses.setdefault(canon, value)
        return guesses

    def _all_mechanical_nodes(self) -> list[str]:
        """Every canonical mechanical (shaft-speed) node -- the union of
        every component's fixed_mechanical_values() and
        free_mechanical_ports() keys, canonicalized. A port only becomes a
        mechanical node by appearing in one of those two hooks; there is no
        separate "this port is mechanical" declaration, mirroring how flow
        ports never need one either.
        """
        nodes: list[str] = []
        for canon in {**self._guess_mechanical_ports(), **self._fixed_mechanical_values()}:
            if canon not in nodes:
                nodes.append(canon)
        return nodes

    def _free_mechanical_nodes(self) -> list[str]:
        fixed = self._fixed_mechanical_values()
        return [n for n in self._all_mechanical_nodes() if n not in fixed]

    def _resolve_node_signal(self, state: "NetworkState") -> dict[str, float]:
        """raw_port_id -> value for every signal port some component
        provides, given the rest of `state` (P/h/mdot/N) already resolved --
        a single pass over every component's provided_signal_values(), not a
        fixed-point loop (unlike _resolve_node_fluid()), since a signal
        value is expected to depend only on its own component's already-
        known port state, never on another component's provided signal.
        Every raw port id sharing a canonical node with a provided value
        gets that same value, the same expansion node_P/node_h/node_mdot get
        in Solver.solve()'s unpack().
        """
        canon_signal: dict[str, float] = {}
        for component in self.components:
            for port_id, value in component.provided_signal_values(state).items():
                canon_signal[self._canonical(port_id)] = value
        if not canon_signal:
            return {}
        node_signal: dict[str, float] = {}
        for component in self.components:
            for port_id in component.signal_ports().values():
                canon = self._canonical(port_id)
                if canon in canon_signal:
                    node_signal[port_id] = canon_signal[canon]
        return node_signal

    def _resolve_node_fluid(self, state: "NetworkState") -> dict[str, "BaseFluid"]:
        """node -> BaseFluid, forward-propagated from every fixed boundary
        node (this network's own default `fluid`, e.g. a Source's inlet air)
        through each component's outlet_fluid() hook -- the composition-
        aware analogue of Solver.solve()'s (P, h) warm-start propagation
        loop, but run against the network's *actual* current state (not
        just a guess), so it stays correct through every Newton iteration
        even when a component's outlet composition itself depends on a free
        unknown (e.g. Combustor's mdot_fuel).

        Same fixed-point technique as that warm-start loop: repeating
        len(components)+1 passes over every component handles them being
        registered in any order (not necessarily upstream-to-downstream) --
        a component whose inlet isn't resolved yet on one pass is simply
        revisited on a later pass, until every reachable node has settled.
        Nodes never reached by any fixed-boundary node (there shouldn't be
        any in a valid topology) are simply absent from the result, and
        NetworkState.fluid_at() falls back to `fluid` for those.

        Two propagation mechanisms per component per pass: merge_fluids()
        first (a genuine multi-inlet merge point, e.g. Junction mixing
        several differently-composed streams — needs every inlet at once,
        see its own docstring), then fluid_flow_pairs()/outlet_fluid() for
        any of that component's outlets merge_fluids() didn't resolve (the
        single-inlet-per-pair pass-through/transform path everything else
        uses). A component that opts fully into merge_fluids() returns []
        from fluid_flow_pairs() so the fallback pass-through path never
        races ahead and locks in a wrong single-inlet answer before all of
        that component's actual inlets are known.

        Fast path: when no component in the network can change composition
        (the common case — only Combustor/Junction can), every reachable node
        resolves to `self.fluid`, and every *unreachable* node's fluid_at()
        also falls back to `self.fluid` — so mapping every node to `self.fluid`
        directly is observationally identical to running the loop, without the
        O(passes x components x ports) union-find traversal. See
        _can_change_composition().
        """
        plan = self._fluid_propagation_plan()

        if not plan.has_composition_change:
            fluid = self.fluid
            return {raw: fluid for raw, _canon in plan.export}

        canon_fluid: dict[str, "BaseFluid"] = dict.fromkeys(plan.seed_canon, self.fluid)
        # The fixed point is reached as soon as a full pass adds nothing new:
        # both write paths are first-wins (setdefault / the "already known"
        # guard), and `state` is constant across the call, so an unproductive
        # pass proves every later pass is unproductive too. This caps a chain
        # at 2 passes instead of the worst-case len(components)+1.
        for _ in range(len(plan.entries) + 1):
            added = False
            for entry in plan.entries:
                if entry.is_merger:
                    inlet_fluids = {
                        port_name: canon_fluid[canon]
                        for port_name, canon in entry.merge_read
                        if canon in canon_fluid
                    }
                    merged = entry.component.merge_fluids(state, inlet_fluids)
                    if merged:
                        for port_name, fluid in merged.items():
                            canon = entry.merge_write.get(port_name)
                            if canon is not None and canon not in canon_fluid:
                                canon_fluid[canon] = fluid
                                added = True

                for inlet_port, outlet_port, inlet_canon, outlet_canon in entry.pairs:
                    if outlet_canon in canon_fluid or inlet_canon not in canon_fluid:
                        continue
                    inlet_fluid = canon_fluid[inlet_canon]
                    outlet_fluid = entry.component.outlet_fluid(
                        state, (inlet_port, outlet_port), inlet_fluid
                    )
                    canon_fluid[outlet_canon] = (
                        outlet_fluid if outlet_fluid is not None else inlet_fluid
                    )
                    added = True
            if not added:
                break

        # Every component addresses its own ports by a raw, component-local
        # id (e.g. "pipe.in"), which may have been merged by connect() into
        # a canonical node shared with other components' raw ids (e.g.
        # "cc1.out") -- the same raw_to_canonical expansion Solver.solve()'s
        # own unpack() does for node_P/node_h/node_mdot, so NetworkState.
        # fluid_at() (always called with a component's own raw id) and
        # SolveResult.node_fluid (which callers also index by raw id, e.g.
        # result.node_fluid["pipe.in"]) both resolve correctly regardless of
        # which raw id happened to become the union-find root.
        node_fluid: dict[str, "BaseFluid"] = {}
        for raw, canon in plan.export:
            fluid = canon_fluid.get(canon)
            if fluid is not None:
                node_fluid[raw] = fluid
        return node_fluid

    def _fluid_propagation_plan(self) -> "_FluidPropagationPlan":
        """Build (and cache per topology version) the state-independent
        structure _resolve_node_fluid() walks. See _FluidPropagationPlan."""
        if self._fluid_plan is not None:
            return self._fluid_plan

        from thermowave.components.base_component import BaseComponent

        has_change = any(_can_change_composition(c) for c in self.components)

        seed_canon = tuple(self._fixed_node_values())

        entries: list = []
        export: list[tuple[str, str]] = []
        for component in self.components:
            ports = component.ports()
            is_merger = (
                type(component).merge_fluids is not BaseComponent.merge_fluids
                or "merge_fluids" in vars(component)
            )
            merge_read = tuple(
                (port_name, self._canonical(port_id))
                for port_name, port_id in ports.items()
            )
            merge_write = {
                port_name: self._canonical(port_id)
                for port_name, port_id in ports.items()
            }
            pairs = tuple(
                (
                    inlet_port,
                    outlet_port,
                    self._canonical(ports[inlet_port]),
                    self._canonical(ports[outlet_port]),
                )
                for inlet_port, outlet_port in component.fluid_flow_pairs()
                if inlet_port in ports and outlet_port in ports
            )
            entries.append(
                _PlanEntry(component, is_merger, merge_read, merge_write, pairs)
            )

            for port_id in ports.values():
                export.append((port_id, self._canonical(port_id)))
            for node_name in component.internal_nodes():
                export.append((node_name, node_name))

        self._fluid_plan = _FluidPropagationPlan(
            has_change, seed_canon, entries, tuple(export)
        )
        return self._fluid_plan

    def check_wiring(self) -> list[str]:
        """Problems with how free parameters are targeted, as plain sentences.

        The solver's own squareness check knows only that unknowns and
        equations disagree by some count — it can't say which parameter is
        loose or which two components are fighting over one. This walks the
        same wiring it does and names them, and is appended to that error
        when it fires (see Solver.solve()).

        Who closes what comes from BaseComponent.closes_parameters(), which
        a Setpoint/Controller/PIDController answers with its own target and a
        Shaft answers with its speed-tied members' speeds. Returns [] when
        nothing looks wrong; a network can still be non-square for reasons
        this doesn't model (a component contributing the wrong number of
        residuals for its own ports), so an empty list is not a guarantee.
        """
        problems: list[str] = []

        free_params: set[str] = set()
        for component in self.components:
            declared = set(component.free_parameters())
            # The solver actually discovers unknowns through
            # guess_free_parameters(), while every controller validates
            # against free_parameters() — override one without the other and
            # the two silently disagree about what's even solvable.
            guessed = set(component.guess_free_parameters(self.fluid, 101325.0, 3e5, 1.0))
            if declared != guessed:
                problems.append(
                    f"{component.name!r} declares free parameters {sorted(declared)} from "
                    f"free_parameters() but {sorted(guessed)} from guess_free_parameters() "
                    f"— the solver uses the latter to create unknowns while Setpoint/"
                    f"Controller/PIDController validate against the former, so these must "
                    f"return the same keys."
                )
            free_params |= {f"{component.name}.{key}" for key in declared | guessed}

        closed_by: dict[str, list[str]] = {name: [] for name in free_params}
        for component in self.components:
            target = getattr(component, "component", None)
            if target is not None and target not in self.components:
                problems.append(
                    f"{component.name!r} targets {getattr(target, 'name', target)!r}, but "
                    f"that component is not in this network — it was probably removed (or "
                    f"never added) while its controller stayed."
                )
                continue
            for name in component.closes_parameters():
                closed_by.setdefault(name, []).append(component.name)

        for name, holders in sorted(closed_by.items()):
            if name not in free_params:
                continue
            if not holders:
                problems.append(
                    f"{name} is free (left as None) but nothing closes it — add a "
                    f"Setpoint/Controller/PIDController to drive it, tie it to a Shaft, "
                    f"or give the component a fixed value instead of None."
                )
            elif len(holders) > 1:
                problems.append(
                    f"{name} is closed by {len(holders)} components ({', '.join(sorted(holders))}) "
                    f"but is only one unknown — each contributes its own residual, so all "
                    f"but one must be removed (see Network.remove_component/"
                    f"replace_component, e.g. swapping a steady-state Controller for a "
                    f"transient PIDController)."
                )

        free_mech_nodes = set(self._free_mechanical_nodes())
        mech_closed_by: dict[str, list[str]] = {n: [] for n in free_mech_nodes}
        for component in self.components:
            for canon in component.closes_mechanical_nodes():
                mech_closed_by.setdefault(canon, []).append(component.name)

        for canon, holders in sorted(mech_closed_by.items()):
            if canon not in free_mech_nodes:
                continue
            if not holders:
                problems.append(
                    f"mechanical node {canon} is free (no component fixes its N) but "
                    f"nothing closes it — add a Setpoint/Controller/PIDController "
                    f"targeting it, tie it to a Shaft, or fix N on one of its members "
                    f"instead of leaving it None."
                )
            elif len(holders) > 1:
                problems.append(
                    f"mechanical node {canon} is closed by {len(holders)} components "
                    f"({', '.join(sorted(holders))}) but is only one unknown."
                )

        return problems

    def validate_topology(self) -> None:
        if not self._fixed_node_values():
            raise NetworkTopologyError(
                "Network has no boundary component fixing node state (P, h). "
                "Add a Source."
            )
        # No analogous check for mdot: a Source may legitimately leave it
        # unfixed (mdot=None) so total mass flow is solved for instead,
        # closed by some other residual elsewhere (e.g. a Sink pinning its
        # inlet pressure) — Solver.solve()'s own square-system check (equal
        # unknown/equation counts) is what actually catches a genuinely
        # under- or over-constrained network, with a more specific message.

    def solve(
        self,
        tol: float = 1e-6,
        max_iter: int = 100,
        damping: float = 1.0,
        verbose: bool = False,
        progress: bool = True,
        dt: float | None = None,
        prev_diff_values: dict[str, float] | None = None,
        warm_start: "SolveResult | None" = None,
        jacobian_reuse: int | None = None,
        step_growth: float | None = None,
    ) -> "SolveResult":
        """dt/prev_diff_values are advanced/internal — see Solver.solve()
        and BaseComponent.differential_parameters(). Ordinary steady-state
        callers should leave both at their defaults. warm_start: see
        Solver.solve()'s docstring — an earlier SolveResult (this network's
        or a related one's) to seed x0 from instead of the flat/propagated
        default guess, for staged/continuation solving.

        progress: a fixed, in-place terminal progress bar for the Newton
        iteration budget — on by default (thermowave.core.progress
        .ProgressBar only actually redraws in place on a real terminal, so
        this is safe to leave on in scripts/logs/tests: it just prints one
        final summary line there instead). Set False for total silence.
        verbose adds iteration/residual/step detail to that bar's text; it
        has no effect when progress=False.

        step_growth: off by default (None -> newton_solve()'s own default
        of 1.0) -- see newton_solve()'s docstring for the full contract.
        Lets damped Newton iterations grow toward a full Newton step once a
        trajectory has demonstrated several steps landing safely, which can
        cut iteration count substantially and let Jacobian reuse actually
        trigger -- but ‖F‖ improving doesn't guarantee the new point stays
        physically valid, and this was found to reliably break a stiff
        real-gas network (chemical equilibrium + coupled PID control) that
        plain fixed damping handles fine. Only pass step_growth > 1.0 for a
        network you've specifically verified still converges correctly with
        it enabled — it is not a safe-by-default speed knob."""
        self.validate_topology()
        from thermowave.core.solver import Solver

        return Solver(self).solve(
            tol=tol, max_iter=max_iter, damping=damping, verbose=verbose, progress=progress,
            dt=dt, prev_diff_values=prev_diff_values, warm_start=warm_start,
            jacobian_reuse=jacobian_reuse, step_growth=step_growth,
        )

    def solve_transient(
        self,
        duration: float,
        dt: float,
        initial: "SolveResult | None" = None,
        tol: float = 1e-6,
        max_iter: int = 100,
        damping: float = 1.0,
        verbose: bool = False,
        progress: bool = True,
        adaptive: bool = False,
        rtol: float = 1e-3,
        atol: float = 1e-6,
        dt_min: float | None = None,
        dt_max: float | None = None,
        safety: float = 0.9,
        growth_limit: float = 5.0,
        shrink_limit: float = 0.2,
        max_step_shrinks: int = 10,
        jacobian_reuse: int | None = None,
        step_growth: float | None = None,
    ) -> "TransientResult":
        """Quasi-steady transient over every differential state any
        component in this network declares (e.g. a dynamic Shaft's rotor
        speed) — see thermowave.core.transient.solve_transient() for the
        full contract, including what adaptive/rtol/atol/dt_min/dt_max and
        the rest of the step-size-control knobs do.

        step_growth: forwarded to every per-timestep Network.solve() call —
        see that method's own docstring. Off by default; only enable for a
        network you've verified converges correctly with it on."""
        self.validate_topology()
        from thermowave.core.transient import solve_transient

        return solve_transient(
            self, duration, dt, initial=initial,
            tol=tol, max_iter=max_iter, damping=damping, verbose=verbose, progress=progress,
            adaptive=adaptive, rtol=rtol, atol=atol, dt_min=dt_min, dt_max=dt_max,
            safety=safety, growth_limit=growth_limit, shrink_limit=shrink_limit,
            max_step_shrinks=max_step_shrinks, jacobian_reuse=jacobian_reuse,
            step_growth=step_growth,
        )
