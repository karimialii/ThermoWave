# Closing a recycle/EGR composition loop

`Recycle` closes a flow loop that has no `Source`/`Sink` boundary — see
[`Recycle`](../components/flow-elements/recycle.md) and the [gas-turbine
tutorial](building-a-gas-turbine-model.md) for the flow side. This guide
covers the separate question: what happens when the fluid recirculating
through that loop has a **composition** of its own, one that depends on
itself through the loop (EGR, a recompression split, any recirculated
stream that's been through a `Combustor` or a `Junction` blend) — and why
that needs one extra argument, `Recycle(fluid_guess=...)`, that plain flow
recirculation doesn't.

A runnable version of everything below is `tests/test_recycle_composition.py`.

---

## Why the flow side and the composition side are different problems

`Network.solve()`'s Newton system is one flat, simultaneous set of
equations — it doesn't care whether the underlying flow graph has a cycle
in it. A `Recycle`-closed Rankine loop with no `Source` anywhere converges
from one ordinary `net.solve()` call (see `tests/test_recycle.py`).

Fluid **composition** doesn't go through that Newton system at all. It's
resolved by a separate pass, `Network._resolve_node_fluid()`, that
forward-propagates composition from the network's own fixed `(P, h)`
boundaries (a `Source`, typically) through each composition-changing
component's `outlet_fluid()`/`merge_fluids()` hook — a `Combustor`'s
reaction products, a `Junction`'s mass-weighted blend. It only ever writes
a node once its inlet is already known.

A loop closed by `Recycle` fixes no `(P, h)` anywhere (that's the whole
point of it — see its own docstring). So if a node's composition depends,
through that loop, on itself, the propagation pass has no seed to start
from: nothing is ever written, and every node in the loop silently falls
back to the network's plain default fluid instead of the real recirculated
composition. No error — the flow side still converges fine, since it
doesn't need composition at all.

`Recycle(fluid_guess=<a BaseFluid>)` gives the propagation pass the seed it
was missing: a starting composition at the `Recycle`'s outlet.
`Network.solve()` then re-solves and compares what actually arrives back at
the `Recycle`'s inlet against that guess, updates the guess, and re-solves
again — direct-substitution ("tear stream") iteration, the same technique
real process simulators use for recycle convergence — until the two agree.

## A minimal example: fresh air merging with a recirculated stream

The smallest network that needs this: a fresh stream (`Source`) merges with
a recirculated stream (`Recycle`) at a `Junction`, and the merged result
feeds straight back around to the `Recycle`'s own inlet.

```python
from thermowave.components import Junction, Recycle, Source
from thermowave.core.network import Network
from thermowave.fluids import CanteraFluid

air = CanteraFluid(name="air", composition="O2:0.21, N2:0.79")
used = CanteraFluid(name="used", composition="CO2:0.10, H2O:0.10, O2:0.10, N2:0.70")

network = Network(fluid=air)
src = Source(name="src", P=3.0e5, T=500.0, mdot=1.0)
junc = Junction(name="junc", n_inlets=2, n_outlets=1)
recycle = Recycle(name="recycle", mdot=1.0, fluid_guess=used)

for c in (src, junc, recycle):
    network.add_component(c)
network.connect(src, "out", junc, "in0")
network.connect(recycle, "out", junc, "in1")
network.connect(junc, "out0", recycle, "in")

result = network.solve(recycle_fluid_tol=1e-6, recycle_fluid_max_iter=30)
```

`fluid_guess=used` is deliberately *wrong* — a different composition than
the fresh stream — to prove this is genuine iteration, not a lucky first
pass. Since nothing in this network reacts (no `Combustor`), the
analytically correct answer is that the recirculating stream settles at
exactly the fresh stream's own composition — a pure mass-weighted mixing
recycle has no other fixed point:

```python
Y_recirc = result.node_fluid["recycle.out"].mass_fractions()
Y_fresh = air.mass_fractions()
# Y_recirc ≈ Y_fresh once converged
```

Two new `solve()` arguments control the outer loop:

- `recycle_fluid_tol` (default `1e-4`) — convergence threshold, max
  absolute per-species mass-fraction difference between successive
  guesses.
- `recycle_fluid_max_iter` (default `25`) — outer iterations before giving
  up. Exceeding it raises `ConvergenceError` naming the last delta reached.

Leave `fluid_guess` unset (the default) for a loop with no
composition-changing component in it — an ordinary Rankine/Brayton cycle
needs none of this, and `Network.solve()` takes the exact same single-call
path as before this feature existed.

## Checking it actually converged

`Recycle.recycle_fluid_delta(state)` is the same check `Network.solve()`
runs internally — call it on the converged result to confirm (or to sanity
check when debugging a `ConvergenceError`):

```python
delta = recycle.recycle_fluid_delta(result.state())
assert delta < 1e-6
```

It returns `None` when there's nothing to check — `fluid_guess` was never
set, or the resolved fluids don't expose `mass_fractions()` — the same
"nothing comparable" fallback `Junction.merge_fluids()` uses rather than
guessing.

## Two traps worth knowing before you hit them

### A `Recycle` directly across a component with a forced state change is over-determined

`Recycle`'s own residuals are an *ideal, lossless pass-through*:
`P_out == P_in`, `h_out == h_in`. Wire it directly across a component whose
own physics force `P`/`h` to differ substantially — a `Combustor` with
`PR < 1` and real heat release, say — and the two requirements contradict
each other outright: `Recycle` insists nothing changes, `Combustor` insists
it does, and there is no solution except the degenerate one (no reaction at
all). Newton doesn't fail cleanly here; it just never converges, sitting at
a fixed nonzero residual norm no matter how many iterations or how much
damping you throw at it.

Put `Recycle` somewhere the loop's *net* pressure/enthalpy change really is
zero — after whatever undoes the forced change (a cooler, a compressor
restoring pressure) — not directly across the component that causes it.
The example above sidesteps this by construction: nothing in that loop
changes `(P, h)` at all, so `Recycle`'s pass-through is trivially
consistent with everything else.

### Composition doesn't propagate through a *second* consecutive `Combustor` pass

`Combustor.outlet_fluid()` only fires when its inlet fluid is *literally* a
`CanteraFluid` instance (see the README's "Composition-aware fluid
propagation" section for why — it's a real-datum consistency requirement,
not an oversight). A `Combustor`'s own reaction product, or a `Junction`
blend, is a different class (`_CanteraCompositionFluid`) — so feeding
either back into a `Combustor`'s inlet across a recycle loop won't retrigger
genuine equilibrium chemistry; `outlet_fluid()` just returns `None`
(pass-through) at that point instead.

This means a genuinely *reactive* EGR loop — recirculated combustion
products actually participating in a second combustion event — isn't
supported yet, independent of the recycle-composition mechanism this guide
covers. Use a `Junction`-only loop (as above) or a single combustion pass
per loop traversal until that gate is loosened.
