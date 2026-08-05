"""Recycle-torn composition cycles -- the EGR-style case Recycle's flow
anchoring alone doesn't cover (see recycle.py's own docstring and
Network.solve()'s recycle_fluid_tol/recycle_fluid_max_iter docstring).

Topology: a fresh Source stream merges (Junction) with a recirculated
stream closed by Recycle -- Source.out -> Junction.in0, Recycle.out ->
Junction.in1, Junction.out0 -> Recycle.in. This is a genuine composition
cycle (Junction's in1 depends, through the loop, on Junction's own out0)
with no reaction involved, so the analytically correct converged
composition at the recirculating node is exactly the fresh stream's own
composition (nothing in the loop consumes or creates any species -- a pure
mass-weighted mixing recycle, Y_recirc = (mdot_fresh*Y_fresh +
mdot_recirc*Y_recirc) / (mdot_fresh + mdot_recirc) has Y_recirc = Y_fresh as
its only fixed point). That makes this a clean, exactly-checkable test of
the tear-stream mechanism itself, distinct from Combustor's own separate
outlet_fluid() restriction (gated on isinstance(CanteraFluid), so a second
composition-changing hop -- recirculating a Combustor's own product back
into itself, or into a Junction blend -- doesn't retrigger real equilibrium
chemistry; a real fully-reactive EGR loop is therefore blocked by that
separate, narrower restriction today, not by the acyclic-propagation issue
fixed here).

Separate file from test_recycle.py because these need the 'cantera' extra
specifically (a real composition-changing component), not 'CoolProp'.
"""
import pytest

pytest.importorskip("cantera")

from thermowave.components.base_component import BaseComponent  # noqa: E402
from thermowave.components.junction import Junction  # noqa: E402
from thermowave.components.recycle import Recycle  # noqa: E402
from thermowave.components.source import Source  # noqa: E402
from thermowave.core.exceptions import ConvergenceError  # noqa: E402
from thermowave.core.network import Network  # noqa: E402
from thermowave.fluids.cantera_fluid import CanteraFluid  # noqa: E402

AIR = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
USED = CanteraFluid(name="used", composition="CO2:0.10, H2O:0.10, O2:0.10, N2:0.70")


def _build_recycle_junction_loop(fluid_guess) -> tuple[Network, Source, Junction, Recycle]:
    net = Network(fluid=AIR)
    src = Source(name="src", P=3.0e5, T=500.0, mdot=1.0)
    junc = Junction(name="junc", n_inlets=2, n_outlets=1)
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=fluid_guess)
    for c in (src, junc, rc):
        net.add_component(c)
    net.connect(src, "out", junc, "in0")
    net.connect(rc, "out", junc, "in1")
    net.connect(junc, "out0", rc, "in")
    return net, src, junc, rc


def test_recycle_junction_loop_without_fluid_guess_never_resolves_recirc_composition():
    """Baseline: today's documented silent-fallback behavior. Junction's
    in1 inlet has no fixed (P, h) boundary feeding it (only Recycle, which
    fixes no (P, h) at all) -- Network._resolve_node_fluid()'s seed set only
    ever reaches the Source side, so merge_fluids() can never fire (it needs
    *both* inlets known at once, see Junction's own docstring) and the whole
    recirculating branch's fluid is never resolved. The flow (P/h/mdot) side
    still solves fine on its own -- this is a composition-only bug."""
    net, src, junc, rc = _build_recycle_junction_loop(fluid_guess=None)

    result = net.solve(tol=1e-8, max_iter=100, progress=False)

    assert result.converged
    # Only the Source-anchored side ever resolves; the recirculating branch
    # (junc.in1, junc.out0, rc.in, rc.out) is silently absent, so fluid_at()
    # on those nodes falls back to the network's plain default -- not
    # necessarily wrong here (see the next test's actual answer), but
    # arrived at by omission, not by resolving the real recirculated blend.
    assert set(result.node_fluid) == {"src.out", "junc.in0"}


def test_recycle_junction_loop_with_fluid_guess_converges_to_analytic_fixed_point():
    """The fix: Recycle(fluid_guess=...) tears the composition cycle open.
    Network.solve()'s outer loop iterates Solver.solve() + Recycle's
    recycle_fluid_delta()/update_fluid_guess() (direct substitution) until
    the guess seeded at the Recycle's outlet matches what Junction's
    merge_fluids() actually sends back around to the Recycle's inlet --
    started from a deliberately wrong guess (USED, not AIR) to prove this
    is genuine iteration, not a lucky first pass."""
    net, src, junc, rc = _build_recycle_junction_loop(fluid_guess=USED)

    result = net.solve(
        tol=1e-8, max_iter=100, progress=False,
        recycle_fluid_tol=1e-6, recycle_fluid_max_iter=30,
    )

    assert result.converged
    # Every node in the loop now resolves, not just the Source-anchored side.
    assert set(result.node_fluid) == {
        "src.out", "junc.in0", "junc.in1", "junc.out0", "rc.in", "rc.out",
    }

    # Analytic fixed point of a non-reactive mixing recycle: the
    # recirculating stream settles at exactly the fresh stream's own
    # composition (nothing in the loop consumes or creates any species).
    Y_recirc = result.node_fluid["rc.out"].mass_fractions()
    Y_fresh = AIR.mass_fractions()
    species = set(Y_recirc) | set(Y_fresh)
    max_gap = max(abs(Y_recirc.get(s, 0.0) - Y_fresh.get(s, 0.0)) for s in species)
    assert max_gap < 1.0e-4

    # The outer loop actually reached its own convergence criterion, not
    # just the ordinary Newton tolerance -- and it started from USED, a
    # materially different composition, so this wasn't a single-pass fluke.
    delta = rc.recycle_fluid_delta(result.state())
    assert delta is not None
    assert delta < 1.0e-6


class _FlipComposition(BaseComponent):
    """Synthetic composition-changing component that always hands back
    whichever of two fixed fluids its inlet *isn't* -- guaranteed to never
    settle, purely to exercise Network.solve()'s ConvergenceError path when
    a Recycle's composition guess genuinely never converges. Not a
    physically meaningful component."""

    def __init__(self, name, fluid_a, fluid_b):
        self.name = name
        self._in = f"{name}.in"
        self._out = f"{name}.out"
        self._fluid_a = fluid_a
        self._fluid_b = fluid_b

    def ports(self):
        return {"in": self._in, "out": self._out}

    def residuals(self, state):
        P_in, h_in = state.node(self._in)
        P_out, h_out = state.node(self._out)
        return [
            P_out - P_in,
            h_out - h_in,
            state.mdot(self._out) - state.mdot(self._in),
        ]

    def outlet_fluid(self, state, pair, inlet_fluid):
        if pair != ("in", "out"):
            return None
        return self._fluid_b if inlet_fluid is self._fluid_a else self._fluid_a


class _TwoSpeciesFluid:
    def __init__(self, Y):
        self._Y = Y

    def mass_fractions(self):
        return dict(self._Y)

    def temperature_ph(self, P, h):
        return 300.0


def test_recycle_composition_loop_raises_when_it_never_settles():
    fluid_a = _TwoSpeciesFluid({"X": 1.0})
    fluid_b = _TwoSpeciesFluid({"Y": 1.0})

    net = Network(fluid=fluid_a)
    flip = _FlipComposition("flip", fluid_a, fluid_b)
    rc = Recycle(name="rc", mdot=1.0, fluid_guess=fluid_a)
    net.add_component(flip)
    net.add_component(rc)
    net.connect(flip, "out", rc, "in")
    net.connect(rc, "out", flip, "in")

    with pytest.raises(ConvergenceError, match="did not converge"):
        net.solve(
            tol=1e-6,
            max_iter=200,
            progress=False,
            recycle_fluid_tol=1e-6,
            recycle_fluid_max_iter=3,
        )
