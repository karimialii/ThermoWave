from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.core.network import Network
    from thermowave.core.solver import SolveResult


class TransientResult:
    """Time history from Network.solve_transient(): one steady-state SolveResult
    per timestep, plus the time history of every differential state any
    component in the network declared (e.g. a dynamic Shaft's rotor speed).

    times/steps are aligned 1:1 (len(times) == len(steps)); diff_history maps
    each differential state's full name ("component.key") to a value list
    the same length as times/steps. diff_history is a convenience only —
    everything in it is also readable from steps[i].params[full_name], or
    (for anything a component also exposes via report_metrics(), like a
    Shaft's "N [rev/min]") via plot()/print_report() below.
    """

    def __init__(
        self,
        times: list[float],
        steps: list["SolveResult"],
        diff_history: dict[str, list[float]],
    ):
        self.times = times
        self.steps = steps
        self.diff_history = diff_history

    def print_report(self, every: int = 1) -> None:
        """Print one steady-state report per recorded timestep (every Nth by
        default all of them) preceded by a time header line."""
        for i in range(0, len(self.times), every):
            print(f"--- t = {self.times[i]:.4f} s ---")
            self.steps[i].print_report()

    def plot(self, *series, ylabel=None, title=None, ax=None, show=True, save_path=None):
        """Plot one or more component metrics vs. time. Each series entry is
        (component, "metric [unit]") or (component, "metric [unit]", label) —
        see thermowave.core.plotting.ThermoPlot.transient() for the full
        contract. Requires the 'plot' extra (matplotlib):
        pip install thermowave[plot].
        """
        from thermowave.core.plotting import ThermoPlot

        plot = ThermoPlot(title=title, xlabel="t [s]", ylabel=ylabel, ax=ax)
        plot.transient(self, series)
        return plot.finish(show=show, save_path=save_path)


def _differential_full_names(network: "Network") -> list[str]:
    return [
        f"{component.name}.{key}"
        for component in network.components
        for key in component.differential_parameters()
    ]


def solve_transient(
    network: "Network",
    duration: float,
    dt: float,
    initial: "SolveResult | None" = None,
    tol: float = 1e-6,
    max_iter: int = 100,
    damping: float = 1.0,
    verbose: bool = False,
) -> TransientResult:
    """Quasi-steady transient over every differential state declared by any
    component in the network (BaseComponent.differential_parameters() /
    state_derivative() — e.g. a dynamic Shaft's rotor speed, driven by
    inertia and net torque). Nothing about *which* components are dynamic is
    passed in here: it's discovered purely from what the network's own
    components declare, so this function has no notion of "shaft" or
    "setpoint" at all — it's the same machinery for any future component
    with its own time constant.

    At each timestep the network's flow/thermal state (and every algebraic
    free parameter, e.g. anything pinned by a Setpoint/Controller) is solved
    exactly as in an ordinary Network.solve(), except each differential
    state's closing equation is backward-Euler instead of "derivative == 0":
    (value - value_at_previous_step) / dt == state_derivative(). Only
    differential states carry a true time history — there is no fluid
    inertia/accumulation term, so every step is otherwise a genuine steady
    state. This holds well when whatever's driving the network's dynamics
    (rotor inertia, ...) has a time constant much slower than the
    flow/thermal ones, which is typical for gas-turbine spool-up/load-step
    transients.

    initial: a SolveResult to start from (its .params supplies each
    differential state's t=0 value) — e.g. a prior Network.solve() at a
    known operating point, or the last step of a previous solve_transient()
    run, to continue it further in time. Defaults to None, which runs an
    ordinary Network.solve() first to establish a genuine steady-state
    initial condition (every differential state's derivative is exactly
    zero at t=0) before stepping forward.

    Also steps any component exposing a step(state, dt) method once per
    timestep (a PIDController's finite-response control law — see its own
    docstring) — discovered the same way, by what the network's components
    actually expose, not by anything passed in here.

    Raises ValueError if the network has neither differential state nor any
    step()-able component — there would be nothing for this to evolve over
    time (a plain Network.solve() already covers that case).
    """
    full_names = _differential_full_names(network)
    pid_like = [c for c in network.components if hasattr(c, "step") and callable(c.step)]
    if not full_names and not pid_like:
        raise ValueError(
            "solve_transient() found nothing to evolve over time — no component "
            "declares differential_parameters() (e.g. a Shaft(..., dynamic=True)) "
            "and none exposes a step() method (e.g. a PIDController). A plain "
            "Network.solve() already covers a network with nothing time-varying."
        )

    if initial is None:
        initial = network.solve(tol=tol, max_iter=max_iter, damping=damping, verbose=verbose)

    n_steps = max(1, round(duration / dt))
    prev_diff_values = {name: initial.params[name] for name in full_names}

    times: list[float] = [0.0]
    steps: list["SolveResult"] = [initial]
    diff_history: dict[str, list[float]] = {name: [prev_diff_values[name]] for name in full_names}

    from thermowave.core.network import NetworkState

    t = 0.0
    warm_start = initial
    for _ in range(n_steps):
        result = network.solve(
            tol=tol, max_iter=max_iter, damping=damping, verbose=verbose,
            dt=dt, prev_diff_values=prev_diff_values, warm_start=warm_start,
        )
        steps.append(result)
        warm_start = result
        prev_diff_values = {name: result.params[name] for name in full_names}
        for name in full_names:
            diff_history[name].append(prev_diff_values[name])

        if pid_like:
            state = NetworkState(
                fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
                node_mdot=result.node_mdot, params=result.params, node_fluid=result.node_fluid,
            )
            for pid in pid_like:
                pid.step(state, dt)

        t += dt
        times.append(t)

    return TransientResult(times, steps, diff_history)
