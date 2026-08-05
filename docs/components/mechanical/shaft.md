# Shaft

<img src="../../_static/diagrams/shaft.svg" alt="Shaft diagram" class="component-diagram">

A mechanical connection: two or more components sharing one physical shaft.
No flow ports of its own — one mechanical port per speed-tied member
(`"m0"`, `"m1"`, ...) and one signal port per member (`"p0"`, `"p1"`, ...,
for the power balance). Connecting a member's mechanical port pulls its
signal port in automatically — see below.

**Ports:** none (flow) — mechanical: `m0`, `m1`, ... &nbsp;·&nbsp; signal:
`p0`, `p1`, ... &nbsp;·&nbsp; **Parameters:** `members`, `gear_ratios`,
`signs`, `efficiency`, `inertia`, `dynamic` (bool), `N0`

`members` lists every component riding on the shaft, in order. `Shaft` gets
one mechanical port per member that has its own `"shaft"` mechanical port (a
real speed unknown — e.g. `Turbine`/`Compressor`, whether `N` is given or
left free) and one signal port per member (every member, including a
torque-only one like `ShaftLoad`, since power still needs to reach the
balance). Connecting a member's mechanical `"shaft"` port also wires its
`"power"` signal port automatically — a real shaft doesn't have a separate
wire for power either, so one `connect()` call per speed-tied member is
enough:

```python
shaft = Shaft("sh1", members=[turb, comp, load])  # signs inferred: [1.0, -1.0, -1.0]
network.add_component(shaft)
network.connect(shaft, "m0", turb, "shaft", kind="mechanical")
network.connect(shaft, "m1", comp, "shaft", kind="mechanical")
network.connect(shaft, "p2", load, "power", kind="signal")
```

`ShaftLoad` has no `"shaft"` mechanical port (no speed unknown of its own —
a pure torque contribution), so there's no mechanical connect() for it to
piggyback on — it still needs its own explicit signal connection, as above.
`signs` defaults to `None`, inferred from each member's own `shaft_sign()`
(pass it explicitly to override, e.g. for a motor running in regenerative
mode); it has one entry per `members` in the order given, but `gear_ratios`
only has one per *speed-tied* member.

Or skip writing any of these `connect()` calls by hand: `shaft.autowire(network)`
(after `add_component(shaft)`) makes every one of them for you, including
`load`'s.

Two modes:

- **`dynamic=False`** (default, steady-state coupling): every speed-tied
  member's shaft speed is tied to the first speed-tied member's (`"m0"`):

  $$N_i - \text{gear\_ratio}_i \cdot N_0 = 0 \quad \text{for each follower } i$$

  `"m0"`'s own speed is still whatever the connected member made it (fixed
  or a free Newton unknown) — a free one still needs an external
  [`Setpoint`](../control/setpoint.md)/[`Controller`](../control/controller.md)
  targeting that member's own `"shaft"` port (e.g. `free_param="shaft"`).

- **`dynamic=True`**: the shaft owns its own rotor speed as **differential
  state**, integrated from net torque and inertia:

  $$\frac{dN}{dt} = \frac{P_\text{net}/\omega}{\text{inertia}}
  \qquad \omega = \frac{2\pi N}{60}$$

  and every speed-tied member is tied to that shared speed:
  $N_i - \text{gear\_ratio}_i \cdot N_\text{shaft} = 0$.
  `Network.solve()` closes `N_shaft` via `d N/dt = 0` (the torque-balance
  speed); `Network.solve_transient()` integrates it forward.

`efficiency` (mechanical transmission efficiency) scales the *net* reported
shaft power — put it on a `dynamic=True` shaft's load (`ShaftLoad`) instead,
since a dynamic shaft's steady equilibrium is exactly net power `== 0`, and
the factor would silently cancel.

---
Part of [Mechanical & electrical](index.md).
