# Source

<img src="../../_static/diagrams/source.svg" alt="Source diagram" class="component-diagram">

A boundary condition fixing outlet pressure and temperature, and (usually)
mass flow rate.

**Ports:** `out`

**Fixes:** `(P, h)` at `out`, from the given `P`/`T` (converted to SI via
`fluid.enthalpy_pt`). If `mdot` is given, it also fixes `mdot` at `out`.

**Free parameter:** none of its own — but leaving `mdot=None` removes the
mass-flow fixing entirely, turning total network flow into an ordinary
Newton unknown (closed by something else downstream, e.g. a [`Sink`](sink.md)
pinning exit pressure). `mdot_guess` seeds that unknown's initial guess.

**Residuals:** none — a `Source` only fixes values, it doesn't compute
anything.

**See also:** for a genuinely closed loop with no open boundary at all (a
recirculating ORC/Rankine/Brayton cycle), see [`Recycle`](recycle.md)
instead — it anchors mass flow without also fixing `(P, h)`.

---
Part of [Flow elements](index.md).
