# CheckValve

<img src="../../_static/diagrams/check_valve.svg" alt="CheckValve diagram" style="max-width:100%">

A one-way variant of [`Valve`](valve.md): forward flow behaves identically,
reverse flow sees a much stiffer resistance — a real boundary condition
against backflow rather than hoping the network's topology never asks for
it.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `D`, `K`, `opening`,
`reverse_factor` (default 1000.0)

$$
K = \begin{cases} K_0 & v \ge 0 \text{ (forward)} \\ K_0 \cdot \text{reverse\_factor} & v < 0 \text{ (reverse)} \end{cases}
\qquad
\Delta p = \frac{K}{\text{opening}^2}\cdot\frac{\rho\, v\,|v|}{2}
$$

Using `v|v|` instead of `Valve`'s `v²` makes the pressure drop track the
*direction* of flow, not just its magnitude — required for a directional
device. `reverse_factor=1000` means about 32× less reverse flow than the
same `|Δp|` would drive forward (flow scales with `1/√K_eff`).

---
Part of [Flow elements](index.md).
