from __future__ import annotations

import math
from typing import TYPE_CHECKING

from thermowave.core import reporting
from thermowave.core.exceptions import ConvergenceError

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


def _weighted_rms_error(
    coarse: "SolveResult", fine: "SolveResult", full_names: list[str], rtol: float, atol: float,
) -> float:
    """Richardson step-doubling error norm, backward-Euler's local
    truncation error being O(dt^2) per step: for a step of size h split
    into two half-steps, |coarse - fine| for each differential state
    estimates the fine (more accurate) result's own local error, up to an
    O(1) constant. Combined across states as a weighted RMS, the classic
    embedded-Runge-Kutta convention (Hairer/Wanner) — <=1 means every
    state's error is within its own tolerance band."""
    sum_sq = 0.0
    for name in full_names:
        c, f = coarse.params[name], fine.params[name]
        scale = atol + rtol * max(abs(c), abs(f))
        sum_sq += ((c - f) / scale) ** 2
    return math.sqrt(sum_sq / len(full_names))


def solve_transient(
    network: "Network",
    duration: float,
    dt: float,
    initial: "SolveResult | None" = None,
    tol: float = 1e-6,
    max_iter: int = 100,
    damping: float = 1.0,
    verbose: bool = False,
    adaptive: bool = False,
    rtol: float = 1e-3,
    atol: float = 1e-6,
    dt_min: float | None = None,
    dt_max: float | None = None,
    safety: float = 0.9,
    growth_limit: float = 5.0,
    shrink_limit: float = 0.2,
    max_step_shrinks: int = 10,
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
    actually expose, not by anything passed in here. When adaptive=True,
    step() only ever runs once per *accepted* step, with that step's actual
    dt — never for a trial/rejected step — since a PID's internal integral
    term accumulates and can't be un-applied.

    adaptive: False (the default) reproduces the original fixed-step
    behavior exactly — duration/dt rounded to an integer step count, every
    step the same size. True switches to step-doubling error control: each
    step is tried once at the current size h and once as two h/2 half-steps;
    since backward-Euler's local error is O(dt^2), the difference between
    the two estimates the half-step (more accurate) result's own error, and
    that result is kept (local extrapolation) if the weighted-RMS error
    across every differential state is <= 1 (rtol/atol combine the same way
    as scipy.integrate.solve_ivp: scale = atol + rtol*|value|). h for the
    *next* step is then rescaled by a classic PI-free (integral-only)
    controller, safety*(1/err)**(1/2) (backward-Euler is a 1st-order
    method), clamped to [shrink_limit, growth_limit]*h and to
    [dt_min, dt_max] if given. A rejected step (err > 1, or the Newton solve
    itself fails to converge at that h) retries at the shrunk h without
    advancing time; max_step_shrinks caps retries per step before raising
    ConvergenceError, so a genuinely unsolvable network still fails loudly
    rather than looping forever. dt is used as the *initial* step size.
    Costs roughly 3x the nonlinear solves of a fixed-step run of the same
    duration (one full-step trial + two half-step trials per accepted step)
    plus whatever rejected retries occur — pay that for a network whose time
    constants vary enough (e.g. a fast initial transient settling into a
    slow one) that no single fixed dt is efficient across the whole run.

    Raises ValueError if the network has neither differential state nor any
    step()-able component — there would be nothing for this to evolve over
    time (a plain Network.solve() already covers that case). adaptive=True
    additionally requires at least one differential state (step-doubling's
    error estimate has nothing to compare otherwise) — a network with only
    step()-able components (e.g. a PIDController with no dynamic Shaft/Tank)
    must use adaptive=False.

    verbose: shows one fixed, in-place progress bar over t/duration for the
    whole run (thermowave.core.progress.ProgressBar — never scrolls, turns
    green on completion), not a per-timestep iteration table — every inner
    Network.solve() call the time-marching loop itself makes is forced to
    verbose=False regardless of this flag, since printing a full Newton
    iteration table once per timestep is exactly the scrolling behavior this
    is meant to replace. The one exception is establishing the t=0
    equilibrium when initial=None: that single Network.solve() does honor
    verbose, since it happens once, before the transient bar starts.
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
    if adaptive and not full_names:
        raise ValueError(
            "solve_transient(adaptive=True) needs at least one differential state "
            "to estimate step-doubling error from (e.g. a Shaft(..., dynamic=True) "
            "or a Tank) — a network with only step()-able components (e.g. a "
            "PIDController alone) has nothing for adaptive error control to compare "
            "against; use adaptive=False (a fixed dt) for it instead."
        )

    if initial is None:
        initial = network.solve(tol=tol, max_iter=max_iter, damping=damping, verbose=verbose)

    prev_diff_values = {name: initial.params[name] for name in full_names}

    times: list[float] = [0.0]
    steps: list["SolveResult"] = [initial]
    diff_history: dict[str, list[float]] = {name: [prev_diff_values[name]] for name in full_names}

    from thermowave.core.network import NetworkState

    bar = reporting.new_progress_bar() if verbose else None
    n_steps_taken = 0

    def _apply_step(result: "SolveResult", h: float) -> None:
        nonlocal prev_diff_values, warm_start, t, n_steps_taken
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
                pid.step(state, h)
        t += h
        times.append(t)
        n_steps_taken += 1
        if bar is not None:
            reporting.render_transient_progress(bar, t, duration, n_steps_taken, h)

    def _solve_step(h: float, prev: dict[str, float], seed: "SolveResult") -> "SolveResult":
        # Always quiet: a per-timestep Newton iteration table is exactly the
        # scrolling behavior the transient progress bar above replaces.
        return network.solve(
            tol=tol, max_iter=max_iter, damping=damping, verbose=False,
            dt=h, prev_diff_values=prev, warm_start=seed,
        )

    t = 0.0
    warm_start = initial

    if not adaptive:
        n_steps = max(1, round(duration / dt))
        for _ in range(n_steps):
            try:
                result = _solve_step(dt, prev_diff_values, warm_start)
            except ConvergenceError:
                if bar is not None:
                    reporting.finish_transient_progress(bar, n_steps_taken, t, success=False)
                raise
            _apply_step(result, dt)
        if bar is not None:
            reporting.finish_transient_progress(bar, n_steps_taken, t)
        return TransientResult(times, steps, diff_history)

    h = dt
    while t < duration - 1e-12:
        h = min(h, duration - t)
        shrinks = 0
        while True:
            try:
                coarse = _solve_step(h, prev_diff_values, warm_start)
                half1 = _solve_step(h / 2.0, prev_diff_values, warm_start)
                half_prev = {name: half1.params[name] for name in full_names}
                fine = _solve_step(h / 2.0, half_prev, half1)
            except ConvergenceError:
                err = float("inf")
            else:
                err = _weighted_rms_error(coarse, fine, full_names, rtol, atol)

            factor = safety * err ** (-0.5) if err > 0 else growth_limit
            factor = min(growth_limit, max(shrink_limit, factor))
            h_next = h * factor
            if dt_min is not None:
                h_next = max(h_next, dt_min)
            if dt_max is not None:
                h_next = min(h_next, dt_max)

            if err <= 1.0:
                _apply_step(fine, h)
                h = h_next
                break

            shrinks += 1
            if h_next >= h or (dt_min is not None and h <= dt_min):
                shrinks = max_step_shrinks  # h can't shrink further; stop retrying
            if shrinks >= max_step_shrinks:
                message = (
                    f"solve_transient(adaptive=True) rejected {shrinks} consecutive "
                    f"steps at t={t:.6g} without meeting rtol={rtol:g}/atol={atol:g} "
                    f"(last dt tried: {h:.6g}, last error: {err:.3g}) — the network's "
                    f"dynamics may be stiffer than dt_min allows, or rtol/atol may be "
                    f"tighter than the Newton solve's own tol can support."
                )
                if bar is not None:
                    reporting.finish_transient_progress(bar, n_steps_taken, t, success=False)
                raise ConvergenceError(message)
            h = h_next

    if bar is not None:
        reporting.finish_transient_progress(bar, n_steps_taken, t)
    return TransientResult(times, steps, diff_history)
