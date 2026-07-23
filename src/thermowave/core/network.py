from __future__ import annotations

from typing import TYPE_CHECKING

import networkx as nx

from thermowave.core.exceptions import NetworkTopologyError

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.solver import SolveResult
    from thermowave.core.transient import TransientResult
    from thermowave.fluids.base_fluid import BaseFluid

# Connection kinds this Network can wire today. Mechanical (shared shaft/speed),
# signal (controller setpoints), and heat-transfer connections are planned but
# not yet implemented — connect() raises NotImplementedError for those kinds.
_SUPPORTED_CONNECTION_KINDS = {"flow"}


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
    ):
        self.fluid = fluid
        self.node_P = node_P
        self.node_h = node_h
        self.node_mdot = node_mdot
        self.params = params if params is not None else {}
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

    Only kind="flow" is implemented: it merges the two ports into a single
    shared (P, h) node, i.e. the two components are physically joined by a
    stream. Future kinds (mechanical, signal, heat) will couple components
    through their own residuals instead of merging node state.
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

    def add_component(self, component: "BaseComponent") -> None:
        self.components.append(component)
        for port_id in component.ports().values():
            self._parent.setdefault(port_id, port_id)
            self.graph.add_node(port_id)
        for node_name in component.internal_nodes():
            self.graph.add_node(node_name)

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
        from_id = self._resolve_port(from_component, from_port)
        to_id = self._resolve_port(to_component, to_port)

        self._union(from_id, to_id)
        self.graph.add_edge(from_id, to_id, kind=kind)

        connection = Connection(from_component, from_port, to_component, to_port, kind)
        self.connections.append(connection)
        return connection

    def _resolve_port(self, component: "BaseComponent", port_name: str) -> str:
        ports = component.ports()
        if port_name not in ports:
            raise NetworkTopologyError(
                f"Component {component.name!r} has no port {port_name!r}; "
                f"available ports: {sorted(ports)}"
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

    def _canonical(self, port_id: str) -> str:
        """This component port's shared network-node id after connect() merges."""
        return self._find(port_id)

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
        """
        canon_fluid: dict[str, "BaseFluid"] = dict.fromkeys(self._fixed_node_values(), self.fluid)
        for _ in range(len(self.components) + 1):
            for component in self.components:
                ports = component.ports()

                inlet_fluids = {
                    port_name: canon_fluid[self._canonical(port_id)]
                    for port_name, port_id in ports.items()
                    if self._canonical(port_id) in canon_fluid
                }
                merged = component.merge_fluids(state, inlet_fluids)
                if merged:
                    for port_name, fluid in merged.items():
                        if port_name not in ports:
                            continue
                        canon = self._canonical(ports[port_name])
                        canon_fluid.setdefault(canon, fluid)

                for pair in component.fluid_flow_pairs():
                    inlet_port, outlet_port = pair
                    if inlet_port not in ports or outlet_port not in ports:
                        continue
                    inlet_canon = self._canonical(ports[inlet_port])
                    outlet_canon = self._canonical(ports[outlet_port])
                    if outlet_canon in canon_fluid or inlet_canon not in canon_fluid:
                        continue
                    inlet_fluid = canon_fluid[inlet_canon]
                    outlet_fluid = component.outlet_fluid(state, pair, inlet_fluid)
                    canon_fluid[outlet_canon] = outlet_fluid if outlet_fluid is not None else inlet_fluid

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
        for component in self.components:
            for port_id in component.ports().values():
                canon = self._canonical(port_id)
                if canon in canon_fluid:
                    node_fluid[port_id] = canon_fluid[canon]
            for node_name in component.internal_nodes():
                if node_name in canon_fluid:
                    node_fluid[node_name] = canon_fluid[node_name]
        return node_fluid

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
        has no effect when progress=False."""
        self.validate_topology()
        from thermowave.core.solver import Solver

        return Solver(self).solve(
            tol=tol, max_iter=max_iter, damping=damping, verbose=verbose, progress=progress,
            dt=dt, prev_diff_values=prev_diff_values, warm_start=warm_start,
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
    ) -> "TransientResult":
        """Quasi-steady transient over every differential state any
        component in this network declares (e.g. a dynamic Shaft's rotor
        speed) — see thermowave.core.transient.solve_transient() for the
        full contract, including what adaptive/rtol/atol/dt_min/dt_max and
        the rest of the step-size-control knobs do."""
        self.validate_topology()
        from thermowave.core.transient import solve_transient

        return solve_transient(
            self, duration, dt, initial=initial,
            tol=tol, max_iter=max_iter, damping=damping, verbose=verbose, progress=progress,
            adaptive=adaptive, rtol=rtol, atol=atol, dt_min=dt_min, dt_max=dt_max,
            safety=safety, growth_limit=growth_limit, shrink_limit=shrink_limit,
            max_step_shrinks=max_step_shrinks,
        )
