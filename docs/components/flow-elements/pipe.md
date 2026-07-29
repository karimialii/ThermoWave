# Pipe

<img src="../../_static/diagrams/pipe.svg" alt="Pipe diagram" style="max-width:100%">

Darcy-Weisbach friction pressure drop, with optional heat loss.
Discretized into `n_elem` sub-elements (default 1), each contributing a
momentum and an energy residual from the pipe's own constant inlet `mdot`.

**Ports:** `in`, `out` &nbsp;·&nbsp; **Parameters:** `L` (length), `D`
(diameter), `f` (Darcy friction factor), `n_elem`, `heat_loss` `[W]`
(optional, total across the whole pipe)

**Per element** (length `L/n_elem`, heat `heat_loss/n_elem`):

$$
\rho = \rho(P, h) \qquad v = \frac{\dot m}{\rho A} \qquad A = \frac{\pi D^2}{4}
$$

$$
\Delta p_\text{friction} = f\,\frac{L_\text{elem}}{D}\,\frac{\rho v^2}{2}
$$

$$
P_\text{in} - P_\text{out} - \Delta p_\text{friction} = 0
\qquad\qquad
h_\text{in} - h_\text{out} - \frac{q_\text{elem}}{\dot m} = 0
$$

plus one mass residual tying the outlet `mdot` back to the inlet's.

---
Part of [Flow elements](index.md).
