from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np

from thermowave.core import reporting
from thermowave.core.exceptions import ConvergenceError, FluidRangeError, NetworkTopologyError

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.network import Network, NetworkState
    from thermowave.fluids.base_fluid import BaseFluid


def _finite_difference_jacobian(
    residual_fn: Callable[[np.ndarray], np.ndarray],
    x: np.ndarray,
    F: np.ndarray,
    eps: float = 1e-6,
) -> np.ndarray:
    n = len(x)
    m = len(F)
    J = np.zeros((m, n))
    for j in range(n):
        x_pert = x.copy()
        step = eps * max(abs(x[j]), 1.0)
        x_pert[j] += step
        F_pert = residual_fn(x_pert)
        J[:, j] = (F_pert - F) / step
    return J


def _ruiz_equilibrate(
    J: np.ndarray, n_passes: int = 3
) -> tuple[np.ndarray, np.ndarray]:
    """Row/column scale vectors (r, c) such that diag(r) @ J @ diag(c) has
    every entry closer to O(1) — Ruiz equilibration, alternating row- and
    column-normalization by max-abs-value a few times until it settles.

    This network's unknowns and residuals mix wildly different physical
    scales in one linear system (pressure ~1e5 Pa, enthalpy ~1e6 J/kg, mass
    flow ~1 kg/s, a compressor's N ~1e4-1e5 rev/min, ...), which — once mass
    flow becomes a free unknown alongside enthalpy (e.g. a Source with
    mdot=None closed by a Sink pinning exit pressure) — can put raw
    (unscaled) finite-difference Jacobian entries many orders of magnitude
    apart, driving its condition number toward the double-precision noise
    floor even when the system isn't actually analytically singular.
    Equilibrating before each linear solve (the same technique general NLP
    solvers like IPOPT use, and what explicit per-stream "scale factors" in
    commercial process simulators approximate by hand) fixes that without
    needing any physical-unit knowledge of which column is which.
    """
    m, n = J.shape
    r = np.ones(m)
    c = np.ones(n)
    A = J.copy()
    # A structurally all-zero row/column (e.g. a residual_fn that's locally
    # constant, or a genuinely singular Jacobian) has nothing for
    # equilibration to normalize against — dividing by the `tiny` floor
    # would otherwise blow the scale factor up toward float64 overflow over
    # a few passes for no benefit (that row/column stays singular either
    # way; np.linalg.solve downstream already raises a clean
    # ConvergenceError for it). Cap row_scale/col_scale instead of letting
    # them run away.
    tiny = 1e-300
    max_scale = 1e100
    for _ in range(n_passes):
        row_scale = np.minimum(1.0 / np.sqrt(np.maximum(np.max(np.abs(A), axis=1), tiny)), max_scale)
        A = A * row_scale[:, None]
        r = np.minimum(r * row_scale, max_scale)
        col_scale = np.minimum(1.0 / np.sqrt(np.maximum(np.max(np.abs(A), axis=0), tiny)), max_scale)
        A = A * col_scale[None, :]
        c = np.minimum(c * col_scale, max_scale)
    return r, c


def newton_solve(
    residual_fn: Callable[[np.ndarray], np.ndarray],
    x0: np.ndarray,
    tol: float = 1e-6,
    max_iter: int = 100,
    damping: float = 1.0,
    clamp_fn: Callable[[np.ndarray], np.ndarray] | None = None,
    step_limit_fn: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None,
    verbose: bool = False,
    progress: bool = True,
    max_jacobian_reuse: int = 4,
    jacobian_reuse_contraction: float = 0.5,
    max_backtracks: int = 8,
    step_growth: float = 1.0,
) -> tuple[np.ndarray, int, float]:
    """Damped Newton-Raphson with a finite-difference Jacobian and a
    backtracking line search.

    progress: show a fixed, in-place terminal progress bar (see
    thermowave.core.progress.ProgressBar) — on by default; live redraws only
    actually appear on a real terminal (a non-tty stream like a pipe/log
    file/pytest's capsys just gets the one-line final summary), so this is
    safe to leave on for library/script use. progress=False suppresses all
    of it, matching the original fully-silent default. verbose only changes
    how much detail that bar's text carries (iteration/residual/step
    numbers) — it has no effect when progress=False.

    max_jacobian_reuse: how many consecutive iterations may reuse the last
    finite-difference Jacobian (and its Ruiz scaling) before a mandatory
    refresh — a modified-Newton / Shamanskii step. Rebuilding the Jacobian
    costs n+1 residual evaluations (the dominant cost of a solve), so while
    the residual is still contracting healthily the *same* Jacobian usually
    still points in a good-enough descent direction, and reusing it turns
    those iterations into a single residual evaluation each. It never changes
    which point is accepted: convergence is still declared only by
    ‖F(x)‖ < tol below, so a stale Jacobian can only lengthen the path, not
    move the solution. A refresh is forced whenever the residual fails to drop
    to <= jacobian_reuse_contraction of the previous iteration's norm (so a
    stalling or diverging step immediately gets a fresh Jacobian), and always
    on the first iteration (which keeps a genuinely singular system raising
    ConvergenceError right away). max_jacobian_reuse=0 disables reuse entirely,
    restoring a plain Newton iteration that rebuilds the Jacobian every step.

    step_limit_fn(x, step) -> step, applied before damping and clamp_fn, caps
    how far a single iteration may move — unlike damping (one scalar factor
    on the whole step vector) or clamp_fn (a hard floor applied after the
    step lands), it can rein in individual components before they overshoot
    at all. Without it, a raw Newton step computed far from the solution can
    send a positive-definite quantity like a mass flow or pressure past zero
    in one iteration; clamp_fn then pins it at its floor, but if it does
    that to several such unknowns simultaneously, they land on the exact
    same floor value — which, for a component whose own physics has a
    derivative that's singular or ill-conditioned exactly when two
    of its inputs coincide (e.g. an effectiveness-NTU heat exchanger's
    Cmin = min(C_hot, C_cold)), reproduces the very problem clamp_fn was
    meant to guard against, one iteration later. This is the same role
    trust-region step bounding plays in commercial equation-oriented process
    solvers (Aspen Plus, HYSYS, ...), which by default cap how much any
    physical variable may change in a single iteration for exactly this
    reason. Default: no limiting.

    damping / max_backtracks / step_growth: damping is no longer applied as
    a flat scalar on every step by itself. Iteration 1 always tries exactly
    `damping * dx`, same as plain fixed-damping Newton; every later
    iteration first tries a *ceiling* that grows by a factor of
    `step_growth` once a step has actually landed safely at the current
    ceiling with no backoff, up to a full Newton step -- and any backoff
    (the ceiling's trial not reducing ‖F‖) resets the ceiling straight back
    down to `damping`. If a trial doesn't reduce ‖F‖, it's retried at half
    that size, then a quarter, and so on, down to `damping` itself, which is
    always the last resort and is accepted unconditionally -- so regardless
    of step_growth, the worst case is exactly plain fixed-damping behavior,
    never worse.

    step_growth defaults to 1.0, which means the ceiling never actually
    grows past `damping` -- i.e. by default this reduces to plain
    fixed-damping Newton, identical to the solver before this line search
    existed. That default is deliberate, not a placeholder: a fixed damping
    (e.g. 0.3) forces linear (geometric, ratio (1-damping)) convergence for
    as long as it's applied, even once the iterate is in Newton's actual
    quadratically-convergent regime, which both wastes iterations directly
    and (since each iteration then only contracts ‖F‖ by ~(1-damping),
    worse than the 0.5 jacobian_reuse_contraction needs for damping=0.3)
    starves the Jacobian-reuse check above of the contraction it needs to
    ever trigger. Growing the ceiling (step_growth > 1.0, e.g. 1.2) can
    remove both costs once a trajectory proves well-behaved -- but ‖F‖
    dropping does not guarantee the new point is well-conditioned or even
    physically valid: on a stiff real-gas network (a combustor's chemical
    equilibrium plus two coupled PID loops was the case this was found on),
    growth eventually drove a trial into a Cantera HP-solve failure
    somewhere along a 65-90-iteration trajectory, reproduced with immediate
    full-step, gradual 1.2x/iteration growth, and even a ceiling hard-capped
    at 2x damping -- pacing the growth more gently only delayed the failure,
    it didn't prevent it. So step_growth > 1.0 is an explicit, per-call
    opt-in for networks/callers who have verified it's safe for their case,
    not a default; enable it and re-verify against your own network before
    relying on it, the same way you'd validate any other solver-tuning
    knob (tol, damping, max_jacobian_reuse) against a new model.
    max_backtracks bounds how many halvings are tried before falling back to
    `damping`; each trial costs one extra residual evaluation (far cheaper
    than a Jacobian rebuild), so this is a small, bounded price for the
    larger steps a proven-safe, step_growth>1.0 trajectory gets to take.
    """
    x = np.array(x0, dtype=float)
    if clamp_fn is not None:
        x = clamp_fn(x)

    bar = None
    if progress:
        if verbose:
            reporting.print_solve_header(len(x), tol, max_iter)
        bar = reporting.new_progress_bar()

    # Jacobian (already Ruiz-scaled) reused across iterations, with its row/
    # column scale vectors — always assigned together so they can't drift out
    # of sync. reuse_count tracks how many iterations the current Jacobian has
    # been reused; prev_norm starts at inf so iteration 1 always refreshes.
    J_scaled: np.ndarray | None = None
    r: np.ndarray | None = None
    c: np.ndarray | None = None
    reuse_count = 0
    prev_norm = float("inf")

    # F/residual_norm are carried across iterations rather than recomputed at
    # the top of each one: the line search below already evaluates
    # residual_fn at whatever point it ends up accepting, so that's reused
    # directly as next iteration's starting F instead of paying for it twice.
    F = residual_fn(x)
    residual_norm = float(np.linalg.norm(F))

    # Line-search ceiling: starts at `damping` -- exactly the old fixed-step
    # behavior -- and only grows (toward a full Newton step) once an
    # iteration succeeds *at* that ceiling without needing to back off.
    # Jumping straight to a full step from iteration 1 was tried and
    # rejected: on a network whose components have real edges (a
    # characteristic map's bounds, a clamp floor two unknowns can collide
    # on), one big early step can land the *next* iteration's Jacobian on a
    # singular point even though the step itself looked like an improvement
    # -- ‖F‖ dropping doesn't guarantee the new point is well-conditioned.
    # Growing gradually only after consecutive steps have already landed
    # safely keeps iteration 1 identical to the original behavior and only
    # takes bigger steps once the trajectory has *demonstrated* it's in a
    # well-behaved region -- any backoff immediately resets the ceiling back
    # down to `damping`, the same instinct as the Jacobian-refresh-on-stall
    # rule above.
    ceiling = damping
    ceiling_growth = step_growth  # see step_growth's own docstring above -- default 1.0 (off)

    def line_search(x_from, F_from, norm_from, dx, start_ceiling):
        """One backtracking line search from (x_from, F_from, norm_from)
        along dx, starting at start_ceiling and halving down to `damping`
        (always accepted unconditionally as the last resort, bounded by
        max_backtracks so a pathological damping<=0 still terminates).
        Returns (x_trial, F_trial, norm_trial, next_ceiling).

        A FluidRangeError from residual_fn (a trial's (P, h, ...) landed
        somewhere a real fluid backend -- Cantera's HP solve, an
        equilibrium calc -- can't evaluate) is treated exactly like a
        residual that didn't improve: back off and try a smaller trial.
        Only the final, most-conservative attempt (use_trial <= damping)
        re-raises it, since there's nothing smaller left to fall back to
        *within this call* -- the caller's rollback-to-last-known-good-point
        handling (below) takes over from there.
        """
        trial = start_ceiling
        backed_off = False
        for attempt in range(max_backtracks):
            use_trial = damping if trial <= damping or attempt == max_backtracks - 1 else trial
            x_trial = x_from + use_trial * dx
            if clamp_fn is not None:
                x_trial = clamp_fn(x_trial)
            try:
                F_trial = residual_fn(x_trial)
            except FluidRangeError:
                if use_trial <= damping:
                    raise
                trial *= 0.5
                backed_off = True
                continue
            norm_trial = float(np.linalg.norm(F_trial))
            if norm_trial <= norm_from or use_trial <= damping:
                break
            trial *= 0.5
            backed_off = True
        next_ceiling = damping if backed_off else min(1.0, start_ceiling * ceiling_growth)
        return x_trial, F_trial, norm_trial, next_ceiling

    # Rollback state: the point/direction that produced the *current* x, so
    # that if the current x turns out to have a singular Jacobian (possible
    # even though it looked like an improvement in ‖F‖ when accepted -- see
    # the docstring above), that step can be undone and retried at a much
    # smaller ceiling instead of failing outright. Only sound to retry
    # something we can actually undo, hence x_prev starting at None.
    x_prev: np.ndarray | None = None
    F_prev: np.ndarray | None = None
    norm_prev_point: float | None = None
    dx_prev: np.ndarray | None = None
    rollback_budget = 6

    for iteration in range(1, max_iter + 1):
        if residual_norm < tol:
            if bar is not None:
                reporting.finish_solve_progress(
                    bar, True, iteration - 1, residual_norm, tol, verbose=verbose
                )
            return x, iteration - 1, residual_norm

        refresh = (
            J_scaled is None
            or max_jacobian_reuse <= 0
            or reuse_count >= max_jacobian_reuse
            or not (residual_norm <= jacobian_reuse_contraction * prev_norm)
        )

        # Everything that can discover "the point we're standing on (x) is
        # actually untenable" lives in this one try: a finite-difference
        # perturbation landing somewhere a fluid backend can't evaluate
        # (FluidRangeError), a Jacobian that doesn't invert (LinAlgError),
        # or the line search itself exhausting every trial down to
        # `damping` without a valid evaluation (FluidRangeError, re-raised
        # by line_search's own last-resort attempt). All three get the same
        # recovery: x wasn't the problem when it was *accepted* -- ‖F‖
        # improving doesn't guarantee the next differentiation or the next
        # step from it stays valid -- so roll back to the last point that
        # demonstrably worked (x_prev) and retry its step at a smaller,
        # forced-conservative ceiling instead of giving up outright.
        try:
            if refresh:
                J = _finite_difference_jacobian(residual_fn, x, F)
                r, c = _ruiz_equilibrate(J)
                J_scaled = (J * r[:, None]) * c[None, :]
                reuse_count = 0
            else:
                reuse_count += 1

            y = np.linalg.solve(J_scaled, -(r * F))
            dx = c * y

            if step_limit_fn is not None:
                dx = step_limit_fn(x, dx)

            x_trial, F_trial, norm_trial, next_ceiling = line_search(x, F, residual_norm, dx, ceiling)
        except (np.linalg.LinAlgError, FluidRangeError) as exc:
            recovered = False
            while x_prev is not None and rollback_budget > 0:
                # Halving the *rollback* ceiling (not `ceiling` itself) on
                # each successive failure means repeated bad luck backs off
                # further each retry, and dx_prev -- computed from a
                # Jacobian at x_prev that *did* invert cleanly -- is reused
                # rather than rederived, since x_prev was never the problem.
                rollback_budget -= 1
                retry_ceiling = damping / (2 ** (6 - rollback_budget))
                try:
                    x, F, residual_norm, _ = line_search(
                        x_prev, F_prev, norm_prev_point, dx_prev, retry_ceiling
                    )
                except FluidRangeError:
                    continue
                recovered = True
                break
            if recovered:
                ceiling = damping
                J_scaled = None  # force a fresh Jacobian at the retried point
                continue
            if bar is not None:
                reporting.finish_solve_progress(
                    bar, False, iteration - 1, residual_norm, tol, verbose=verbose
                )
            raise ConvergenceError(
                f"Solve failed at iteration {iteration}: {exc}"
            ) from exc

        ceiling = next_ceiling
        if bar is not None:
            reporting.render_solve_progress(
                bar, iteration, max_iter, residual_norm,
                float(np.linalg.norm(x_trial - x)), verbose=verbose,
            )

        prev_norm = residual_norm
        x_prev, F_prev, norm_prev_point, dx_prev = x, F, residual_norm, dx
        rollback_budget = 6
        x, F, residual_norm = x_trial, F_trial, norm_trial

    if residual_norm < tol:
        if bar is not None:
            reporting.finish_solve_progress(
                bar, True, max_iter, residual_norm, tol, verbose=verbose
            )
        return x, max_iter, residual_norm

    if bar is not None:
        reporting.finish_solve_progress(bar, False, max_iter, residual_norm, tol, verbose=verbose)
    raise ConvergenceError(
        f"Solver failed to converge after {max_iter} iterations "
        f"(residual norm={residual_norm:.3e}, tol={tol:.3e})"
    )


def _failure_context(
    message: str,
    network: "Network",
    n_unknowns: int,
    n_equations: int,
    dt: float | None,
) -> str:
    """A failed solve's message, plus what the solver knows about the
    network that produced it.

    A ConvergenceError on its own ("Singular Jacobian at iteration 12")
    gives a caller nothing to act on: it says the linear algebra broke, not
    which part of their model broke it. The system size and any
    Network.check_wiring() findings are both cheap and only meaningful
    here, so they get attached at the point of failure rather than left for
    the caller to go find with verbose=True and a second run.
    """
    lines = [
        message,
        "",
        f"System: {n_unknowns} unknown(s), {n_equations} equation(s)"
        + (f", backward-Euler step dt={dt:g}" if dt is not None else ", steady state"),
    ]

    problems = network.check_wiring()
    if problems:
        lines += ["", "Wiring problems found:"]
        lines += [f"  - {p}" for p in problems]

    lines += [
        "",
        "A square system that still won't converge is usually a starting-guess "
        "problem rather than a wiring one. Things that help, roughly in order:",
        "  - warm_start=<an earlier SolveResult> from an easier version of this "
        "network (a lower load, a fixed value where this one solves for it), then "
        "step toward the point you actually want -- continuation.",
        "  - a lower damping= (0.1-0.3) so early iterations can't overshoot off "
        "the edge of a component's map.",
        "  - fixing a free parameter at a plausible value to check the rest of "
        "the network converges at all, then freeing it again with that result as "
        "warm_start.",
    ]
    return "\n".join(lines)


class SolveResult:
    """Outcome of a Network.solve() call."""

    def __init__(
        self,
        converged: bool,
        iterations: int,
        residual_norm: float,
        node_P: dict[str, float],
        node_h: dict[str, float],
        node_mdot: dict[str, float],
        params: dict[str, float],
        fluid: "BaseFluid",
        node_order: list[str],
        components: list["BaseComponent"],
        node_fluid: dict[str, "BaseFluid"] | None = None,
        node_N: dict[str, float] | None = None,
        network: "Network | None" = None,
    ):
        self.converged = converged
        self.iterations = iterations
        self.residual_norm = residual_norm
        self.node_P = node_P
        self.node_h = node_h
        self.node_mdot = node_mdot
        self.params = params
        self.fluid = fluid
        self.node_order = node_order
        self.components = components
        # node -> BaseFluid, resolved once at the converged solution (see
        # Network._resolve_node_fluid()) -- empty for any network with no
        # composition-changing component, in which case every node just
        # falls back to `fluid` via NetworkState.fluid_at().
        self.node_fluid = node_fluid if node_fluid is not None else {}
        # Shaft speed [rev/min] at every mechanical port -- empty for any
        # network with no mechanical ports at all.
        self.node_N = node_N if node_N is not None else {}
        # Kept only so state() can re-resolve node_signal (a signal port's
        # value is computed from the rest of state, not stored as its own
        # array like node_P/node_h are, so it isn't captured anywhere else).
        self._network = network

    def state(self) -> "NetworkState":
        """This result as a NetworkState, for calling a component's own
        report_metrics()/Q() against the converged solution.

        Reading anything a component computes rather than stores (a
        compressor's PR, a heat path's Q, a shaft's net power) means handing
        it a NetworkState, and the only place one otherwise exists is inside
        the solve that just finished. Without this, every caller rebuilds the
        same six-field object by hand.
        """
        from thermowave.core.network import NetworkState

        state = NetworkState(
            fluid=self.fluid,
            node_P=self.node_P,
            node_h=self.node_h,
            node_mdot=self.node_mdot,
            params=self.params,
            node_fluid=self.node_fluid,
            node_N=self.node_N,
        )
        if self._network is not None:
            state.node_signal = self._network._resolve_node_signal(state)
        return state

    def print_report(self) -> None:
        """Print a nicely formatted summary + per-node (P, T, h) table."""
        print(reporting.format_solve_report(self))


class Solver:
    """Assembles a Network's free variables into a Newton-Raphson problem."""

    P_MIN = 1.0e3  # Pa, physical clamp floor to protect fluid property calls
    MDOT_MIN = 1.0e-6  # kg/s, physical clamp floor — mass flow cannot reverse
    N_MIN = 1.0  # rev/min, physical clamp floor — shaft speed cannot reverse

    # step_limit()'s per-iteration bounds — named here (rather than as bare
    # literals in the two nearly-identical blocks that use them) so the two
    # blocks can't silently drift apart, and so tuning any of these has one
    # place to change.
    STEP_LOWER_FACTOR = 0.5  # positive-definite unknowns: at most an 80% drop...
    STEP_UPPER_FACTOR = 2.0  # ...or a 5x jump, per Newton iteration
    H_CAP_FACTOR = 10.0  # enthalpy step cap, as a multiple of H_CAP_FLOOR/|default_h|
    H_CAP_FLOOR = 1.0e5  # J/kg, floor under H_CAP_FACTOR's scale reference
    PARAM_CAP_FACTOR = 10.0  # non-positive free-parameter step cap, as a
    # multiple of PARAM_SCALE_FLOOR/|initial guess|
    PARAM_SCALE_FLOOR = 1.0  # floor under PARAM_CAP_FACTOR's scale reference

    # Fallback (P, h) guess used only when the network has no fixed nodes at
    # all to seed a reference state from (a degenerate edge case) — roughly
    # atmospheric pressure and a modest reference enthalpy.
    DEFAULT_P_FALLBACK = 1.0e5  # Pa
    DEFAULT_H_FALLBACK = 3.0e5  # J/kg
    DEFAULT_N_FALLBACK = 1.0e4  # rev/min

    def __init__(self, network: "Network"):
        self.network = network

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
    ) -> SolveResult:
        """dt/prev_diff_values are advanced/internal knobs used by
        Network.solve_transient() — see its docstring and
        BaseComponent.differential_parameters()/state_derivative() for the
        full contract. Leave both at their defaults (None) for an ordinary
        steady-state solve, which is what Network.solve() does.

        step_growth: forwarded to newton_solve() (see its own docstring for
        the full contract and why the default there is 1.0, i.e. off).
        Leave at the default None here (also off) unless you've verified
        step_growth > 1.0 is safe for *this* network -- it was found to
        drive a stiff real-gas network into solver failures that plain
        fixed damping doesn't, so it's an explicit per-call opt-in, not a
        setting to enable broadly.

        warm_start: an earlier SolveResult (from this network or a related
        one — e.g. the same physical network with fewer free unknowns, one
        component fixed instead of solved for) whose node_P/node_h/node_mdot
        and params override the usual flat/propagated x0 guess wherever a
        node or param name matches. This is what makes staged/continuation
        solving practical: converging a genuinely hard, tightly-coupled
        network (several turbomachinery/heat-exchanger free unknowns at
        once) from a cold start can walk straight into a real fold in the
        solution manifold, or just be needlessly fragile — converging an
        easier version of the same network first (e.g. one target fixed)
        and handing its result in here as warm_start for the full version
        is the standard fix (the same role step-by-step continuation plays
        in commercial process/turbomachinery solvers). Node/param names that
        don't appear in warm_start fall back to the normal guess.
        """
        from thermowave.core.network import NetworkState

        network = self.network
        # A component may have been mutated in place between solves (e.g. a
        # Controller swapped for a PIDController, or a component's port set
        # changed) without going through add_component/connect. Rebuild the
        # topology-derived caches once here so this solve sees current
        # topology; the rebuild is <0.01s and only happens per solve, not per
        # residual evaluation.
        network.invalidate_caches()
        fluid = network.fluid
        fixed_nodes = network._fixed_node_values()
        fixed_mdot = network._fixed_node_mdot()
        guess_mdot = network._guess_node_mdot()
        all_nodes = network._all_nodes()
        port_nodes = network._port_nodes()
        free_nodes = [n for n in all_nodes if n not in fixed_nodes]
        mdot_free_nodes = [n for n in port_nodes if n not in fixed_mdot]

        fixed_mech = network._fixed_mechanical_values()
        mech_free_nodes = network._free_mechanical_nodes()
        # free_params/param_names are computed further down, after the
        # warm-start (P, h) propagation below — components use that
        # propagated inlet guess (via guess_free_parameters()) rather than a
        # hardcoded reference. unpack()/residual_vector() below read these
        # as closures, resolved at call time, so defining them later is fine.
        free_params: dict[str, float] = {}
        param_names: list[str] = []

        # Components address their own ports by a raw, component-local id (e.g.
        # "p1.in"), which connect() may have merged into a shared canonical
        # node id (e.g. "src.out") along with other components' ports. The
        # unknown vector and fixed_nodes are keyed by canonical id; this map
        # lets unpack() expand back out so every component's residuals() can
        # call state.node(<its own raw port id>) and get the right value.
        raw_to_canonical: dict[str, str] = {}
        for component in network.components:
            for port_id in component.ports().values():
                raw_to_canonical[port_id] = network._canonical(port_id)
            for port_id in component.mechanical_ports().values():
                raw_to_canonical[port_id] = network._canonical(port_id)
            for port_id in component.signal_ports().values():
                raw_to_canonical[port_id] = network._canonical(port_id)
            for node_name in component.internal_nodes():
                raw_to_canonical[node_name] = node_name
        # unpack() runs once per residual evaluation (hundreds of thousands of
        # times per transient run); freeze the mapping as a tuple of pairs so
        # each call iterates a flat sequence instead of re-walking the dict.
        raw_to_canonical_items = tuple(raw_to_canonical.items())

        # Differential unknowns: one per (component, key) declared via
        # differential_parameters() (e.g. a dynamic Shaft's rotor speed).
        # Unlike free_parameters(), these are closed automatically by the
        # solver itself (see residual_vector() below) rather than needing a
        # matching residual wired up elsewhere: state_derivative() == 0 for
        # an ordinary steady-state solve (dt=None), or the backward-Euler
        # equation (x - x_prev)/dt - state_derivative() == 0 for a
        # Network.solve_transient() step (dt given). prev_diff_values seeds
        # both the initial guess and (in transient mode) x_prev; components
        # with no prior value (steady state, or a differential parameter new
        # since the last step) fall back to their own declared initial guess.
        prev_diff_values = prev_diff_values or {}
        diff_owner: dict[str, tuple["BaseComponent", str]] = {}
        diff_params: dict[str, float] = {}
        for component in network.components:
            for key, guess in component.differential_parameters().items():
                full_name = f"{component.name}.{key}"
                diff_owner[full_name] = (component, key)
                diff_params[full_name] = prev_diff_values.get(full_name, guess)

        n_pressure = len(free_nodes)
        n_mdot = len(mdot_free_nodes)
        n_mech = len(mech_free_nodes)

        def unpack(
            x: np.ndarray,
        ) -> tuple[
            dict[str, float], dict[str, float], dict[str, float], dict[str, float], dict[str, float]
        ]:
            canonical_P = {n: v[0] for n, v in fixed_nodes.items()}
            canonical_h = {n: v[1] for n, v in fixed_nodes.items()}
            for i, n in enumerate(free_nodes):
                canonical_P[n] = x[2 * i]
                canonical_h[n] = x[2 * i + 1]
            canonical_mdot = dict(fixed_mdot)
            offset = 2 * n_pressure
            for i, n in enumerate(mdot_free_nodes):
                canonical_mdot[n] = x[offset + i]
            canonical_N = dict(fixed_mech)
            offset = 2 * n_pressure + n_mdot
            for i, n in enumerate(mech_free_nodes):
                canonical_N[n] = x[offset + i]
            params = {}
            offset = 2 * n_pressure + n_mdot + n_mech
            for i, name in enumerate(param_names):
                params[name] = x[offset + i]

            # Mechanical/signal port raw ids share this same raw_to_canonical
            # map (see its construction above) but have no entry in
            # canonical_P/canonical_h/canonical_mdot (those dicts only ever
            # cover flow nodes) -- filtered out here exactly like
            # canonical_mdot already is, rather than assuming every raw id
            # is a flow node.
            node_P = {
                raw: canonical_P[canon]
                for raw, canon in raw_to_canonical_items
                if canon in canonical_P
            }
            node_h = {
                raw: canonical_h[canon]
                for raw, canon in raw_to_canonical_items
                if canon in canonical_h
            }
            node_mdot = {
                raw: canonical_mdot[canon]
                for raw, canon in raw_to_canonical_items
                if canon in canonical_mdot
            }
            node_N = {
                raw: canonical_N[canon]
                for raw, canon in raw_to_canonical_items
                if canon in canonical_N
            }
            return node_P, node_h, node_mdot, node_N, params

        def residual_vector(x: np.ndarray) -> np.ndarray:
            node_P, node_h, node_mdot, node_N, params = unpack(x)
            state = NetworkState(
                fluid=fluid,
                node_P=node_P,
                node_h=node_h,
                node_mdot=node_mdot,
                params=params,
                node_N=node_N,
            )
            state.node_fluid = network._resolve_node_fluid(state)
            state.node_signal = network._resolve_node_signal(state)
            residuals: list[float] = []
            for component in network.components:
                residuals.extend(component.residuals(state))
            for full_name, (component, key) in diff_owner.items():
                rate = component.state_derivative(state)[key]
                if dt is None:
                    residuals.append(rate)
                else:
                    prev_value = prev_diff_values.get(full_name, params[full_name])
                    residuals.append((params[full_name] - prev_value) / dt - rate)
            return np.array(residuals, dtype=float)

        def clamp(x: np.ndarray) -> np.ndarray:
            x = x.copy()
            for i in range(n_pressure):
                if x[2 * i] < self.P_MIN:
                    x[2 * i] = self.P_MIN
            offset = 2 * n_pressure
            for i in range(n_mdot):
                if x[offset + i] < self.MDOT_MIN:
                    x[offset + i] = self.MDOT_MIN
            offset = 2 * n_pressure + n_mdot
            for i in range(n_mech):
                if x[offset + i] < self.N_MIN:
                    x[offset + i] = self.N_MIN
            return x

        def step_limit(x: np.ndarray, dx: np.ndarray) -> np.ndarray:
            # P and mdot are both positive-definite; a raw Newton step
            # computed far from the solution can otherwise overshoot one (or
            # several at once) past zero in a single iteration — see
            # newton_solve's step_limit_fn docstring for why that's worse
            # than it sounds once clamp_fn is involved. Cap each to at most
            # an 80% drop or a 5x jump per iteration instead.
            dx = dx.copy()
            positive_definite_indices = [2 * i for i in range(n_pressure)]
            offset = 2 * n_pressure
            positive_definite_indices += [offset + i for i in range(n_mdot)]
            offset = 2 * n_pressure + n_mdot
            positive_definite_indices += [offset + i for i in range(n_mech)]
            for i in positive_definite_indices:
                current = x[i]
                if current <= 0:
                    continue
                lower, upper = self.STEP_LOWER_FACTOR * current, self.STEP_UPPER_FACTOR * current
                proposed = current + dx[i]
                if proposed < lower:
                    dx[i] = lower - current
                elif proposed > upper:
                    dx[i] = upper - current
            # h has no sign constraint (fluid-model-dependent reference
            # state), so an absolute cap scaled to the network's own
            # reference enthalpy stands in for the relative cap above —
            # without it, the same runaway-step failure mode above can hit
            # h just as easily (a huge finite-difference coefficient
            # anywhere in the row is enough), just without a floor to catch
            # it, so it diverges instead of merely landing on a bad tie.
            h_cap = self.H_CAP_FACTOR * max(abs(default_h), self.H_CAP_FLOOR)
            for i in range(n_pressure):
                h_idx = 2 * i + 1
                if abs(dx[h_idx]) > h_cap:
                    dx[h_idx] = h_cap if dx[h_idx] > 0 else -h_cap
            # Free/differential parameters (a compressor's N, a combustor's
            # mdot_fuel, a PID's output, a dynamic shaft's speed, ...) are
            # just as exposed to the runaway-step failure mode above as P,
            # h, and mdot are, but this closure is built before param_names/
            # all_params_guess are populated further down — both are read by
            # name here (not passed in), so this sees their final contents
            # by the time newton_solve actually calls step_limit.
            param_offset = 2 * n_pressure + n_mdot + n_mech
            for i, name in enumerate(param_names):
                idx = param_offset + i
                current = x[idx]
                if current > 0:
                    # Most declared free/differential parameters are
                    # positive-definite physical quantities (shaft speed,
                    # fuel mass flow, ...) — reuse the same relative bound
                    # as P/mdot above.
                    lower, upper = self.STEP_LOWER_FACTOR * current, self.STEP_UPPER_FACTOR * current
                    proposed = current + dx[idx]
                    if proposed < lower:
                        dx[idx] = lower - current
                    elif proposed > upper:
                        dx[idx] = upper - current
                else:
                    # Non-positive (e.g. a PID/controller output that can
                    # legitimately swing negative, or a param still sitting
                    # on a zero seed) has no sign to bound a relative step
                    # against — cap the step to a magnitude scaled off its
                    # own initial guess instead, since params span wildly
                    # different physical scales (rev/min, kg/s, dimensionless
                    # PID output, ...) with no shared reference like h has.
                    scale = max(abs(all_params_guess.get(name, 0.0)), self.PARAM_SCALE_FLOOR)
                    param_cap = self.PARAM_CAP_FACTOR * scale
                    if abs(dx[idx]) > param_cap:
                        dx[idx] = param_cap if dx[idx] > 0 else -param_cap
            return dx

        if fixed_nodes:
            default_P, default_h = next(iter(fixed_nodes.values()))
        else:
            default_P, default_h = self.DEFAULT_P_FALLBACK, self.DEFAULT_H_FALLBACK
        default_mdot = next(iter(fixed_mdot.values())) if fixed_mdot else 1.0

        # A flat default (every free node starts at the network's first fixed
        # (P, h), typically an ambient Source) is a poor initial guess several
        # components downstream — e.g. a turbine after a compressor + heater
        # would start Newton thinking it sees ~1 bar instead of ~4 bar, which
        # can put map-based components so far outside their valid range that
        # the very first Jacobian is singular. Forward-propagate a rough
        # guess instead, via each component's warm_start_pairs() (default:
        # the single "in"/"out" pair) and guess_outlet()/
        # guess_outlet_for_pair() hooks. Repeated passes handle components
        # added out of topological order; extra passes beyond convergence
        # are harmless since already-guessed nodes are skipped.
        canonical_P_guess: dict[str, float] = {n: v[0] for n, v in fixed_nodes.items()}
        canonical_h_guess: dict[str, float] = {n: v[1] for n, v in fixed_nodes.items()}
        # Seeded the same way as P/h above, but from fixed_mdot and (for a
        # Source left with mdot=None) guess_mdot instead of fixed_nodes —
        # then propagated forward by the same pass, via guess_outlet_mdot()
        # (default: pass through unchanged; SimpleCombustor/Combustor add
        # their own fuel guess). Without this, every free mdot node
        # downstream of the seed defaults to the exact same flat guess,
        # which is harmless on its own but can make two streams that are
        # supposed to differ (e.g. a recuperator's hot vs. cold side once
        # one has picked up fuel mass) land on an artificial exact tie —
        # and if a component's physics has a derivative that's only
        # ill-conditioned exactly at such a tie (an effectiveness-NTU heat
        # exchanger's Cmin = min(C_hot, C_cold), even smoothed), that alone
        # can make the very first Jacobian singular.
        canonical_mdot_guess: dict[str, float] = {**fixed_mdot, **guess_mdot}
        default_N = next(iter(fixed_mech.values())) if fixed_mech else self.DEFAULT_N_FALLBACK
        # Free nodes' guesses are seeded below, per component, via
        # guess_free_mechanical_ports() (context-aware; defaults to plain
        # free_mechanical_ports() when a component doesn't override it) --
        # fixed_mech is the only thing seeded here.
        canonical_N_guess: dict[str, float] = dict(fixed_mech)
        for _ in range(len(network.components) + 1):
            for component in network.components:
                ports = component.ports()
                pairs = component.warm_start_pairs()
                for pair in pairs:
                    inlet_port, outlet_port = pair
                    if inlet_port not in ports or outlet_port not in ports:
                        continue
                    inlet_canon = network._canonical(ports[inlet_port])
                    outlet_canon = network._canonical(ports[outlet_port])
                    if (
                        inlet_canon not in canonical_P_guess
                        or outlet_canon in canonical_P_guess
                    ):
                        continue
                    mdot_in_guess = canonical_mdot_guess.get(inlet_canon, default_mdot)
                    if pairs == [("in", "out")]:
                        P_guess, h_guess = component.guess_outlet(
                            canonical_P_guess[inlet_canon],
                            canonical_h_guess[inlet_canon],
                            mdot_in_guess,
                        )
                    else:
                        P_guess, h_guess = component.guess_outlet_for_pair(
                            pair,
                            canonical_P_guess[inlet_canon],
                            canonical_h_guess[inlet_canon],
                            mdot_in_guess,
                        )
                    canonical_P_guess[outlet_canon] = P_guess
                    canonical_h_guess[outlet_canon] = h_guess
                    canonical_mdot_guess.setdefault(
                        outlet_canon, component.guess_outlet_mdot(pair, mdot_in_guess)
                    )
                    for internal_node in component.internal_nodes():
                        canonical_P_guess.setdefault(internal_node, P_guess)
                        canonical_h_guess.setdefault(internal_node, h_guess)

        # Now that every node has a (possibly propagated) warm-start guess,
        # let each component guess its own free parameters using its own
        # inlet's guess — e.g. a turbine's shaft-speed guess computed from
        # the actual warm-started inlet temperature it'll see, not a
        # hardcoded reference that's only right for some positions in the
        # network.
        for component in network.components:
            ports = component.ports()
            if "in" in ports:
                inlet_canon = network._canonical(ports["in"])
                P_in_guess = canonical_P_guess.get(inlet_canon, default_P)
                h_in_guess = canonical_h_guess.get(inlet_canon, default_h)
            else:
                P_in_guess, h_in_guess = default_P, default_h
            for key, guess in component.guess_free_parameters(
                fluid, P_in_guess, h_in_guess, default_mdot
            ).items():
                free_params[f"{component.name}.{key}"] = guess
            for port_id, guess in component.guess_free_mechanical_ports(
                fluid, P_in_guess, h_in_guess, default_mdot
            ).items():
                canonical_N_guess.setdefault(network._canonical(port_id), guess)
        param_names[:] = list(free_params.keys()) + list(diff_owner.keys())
        all_params_guess = {**free_params, **diff_params}

        if warm_start is not None:
            for n in free_nodes:
                if n in warm_start.node_P:
                    canonical_P_guess[n] = warm_start.node_P[n]
                if n in warm_start.node_h:
                    canonical_h_guess[n] = warm_start.node_h[n]
            for n in mdot_free_nodes:
                if n in warm_start.node_mdot:
                    canonical_mdot_guess[n] = warm_start.node_mdot[n]
            for n in mech_free_nodes:
                if n in warm_start.node_N:
                    canonical_N_guess[n] = warm_start.node_N[n]
            for name in param_names:
                if name in warm_start.params:
                    all_params_guess[name] = warm_start.params[name]

        n_unknowns = 2 * n_pressure + n_mdot + n_mech + len(param_names)

        x0 = np.zeros(n_unknowns)
        for i, n in enumerate(free_nodes):
            x0[2 * i] = canonical_P_guess.get(n, default_P)
            x0[2 * i + 1] = canonical_h_guess.get(n, default_h)
        offset = 2 * n_pressure
        for i, n in enumerate(mdot_free_nodes):
            x0[offset + i] = canonical_mdot_guess.get(n, default_mdot)
        offset = 2 * n_pressure + n_mdot
        for i, n in enumerate(mech_free_nodes):
            x0[offset + i] = canonical_N_guess.get(n, default_N)
        offset = 2 * n_pressure + n_mdot + n_mech
        for i, name in enumerate(param_names):
            x0[offset + i] = all_params_guess[name]

        F0 = residual_vector(x0)
        n_equations = len(F0)

        if verbose and progress:
            reporting.print_system_summary(n_unknowns, n_equations)

        if n_unknowns != n_equations:
            message = (
                f"Network is not solvable as wired: {n_unknowns} unknown(s) but "
                f"{n_equations} equation(s) (a Newton solve needs a square system, "
                f"equal counts of both). Check that every free parameter (e.g. a "
                f"component left with N=None) has exactly one matching target "
                f"residual (e.g. a Setpoint), and that nothing is targeted twice."
            )
            problems = network.check_wiring()
            if problems:
                message += "\n\nLikely cause:\n" + "\n".join(f"  - {p}" for p in problems)
            raise NetworkTopologyError(message)

        if n_unknowns == 0:
            node_P, node_h, node_mdot, node_N, params = unpack(x0)
            final_state = NetworkState(
                fluid=fluid, node_P=node_P, node_h=node_h, node_mdot=node_mdot, params=params,
                node_N=node_N,
            )
            node_fluid = network._resolve_node_fluid(final_state)
            return SolveResult(
                True,
                0,
                float(np.linalg.norm(F0)),
                node_P,
                node_h,
                node_mdot,
                params,
                fluid,
                port_nodes,
                network.components,
                node_fluid=node_fluid,
                node_N=node_N,
                network=network,
            )

        try:
            newton_kwargs = {}
            if jacobian_reuse is not None:
                newton_kwargs["max_jacobian_reuse"] = jacobian_reuse
            if step_growth is not None:
                newton_kwargs["step_growth"] = step_growth
            x_sol, iterations, residual_norm = newton_solve(
                residual_vector,
                x0,
                tol=tol,
                max_iter=max_iter,
                damping=damping,
                clamp_fn=clamp,
                step_limit_fn=step_limit,
                verbose=verbose,
                progress=progress,
                **newton_kwargs,
            )
        except ConvergenceError as exc:
            # A bare "Singular Jacobian at iteration 12" says nothing about
            # the network that produced it. Attach the system size and any
            # wiring findings here, where they're actually known -- failure
            # is exactly when that context is worth having, and it costs
            # nothing on the path where the solve works.
            raise ConvergenceError(
                _failure_context(str(exc), network, n_unknowns, n_equations, dt)
            ) from exc
        node_P, node_h, node_mdot, node_N, params = unpack(x_sol)
        final_state = NetworkState(
            fluid=fluid, node_P=node_P, node_h=node_h, node_mdot=node_mdot, params=params,
            node_N=node_N,
        )
        node_fluid = network._resolve_node_fluid(final_state)
        return SolveResult(
            True,
            iterations,
            residual_norm,
            node_P,
            node_h,
            node_mdot,
            params,
            fluid,
            port_nodes,
            network.components,
            node_fluid=node_fluid,
            node_N=node_N,
            network=network,
        )
