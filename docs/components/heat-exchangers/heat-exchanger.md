# HeatExchanger

<img src="../../_static/diagrams/heat_exchanger.svg" alt="HeatExchanger diagram" style="max-width:100%">

Two-stream heat exchanger: fixed effectiveness, or effectiveness derived
from geometry (`UA`, flow arrangement) — auto-selected by which one you
give. Give exactly one of `effectiveness` or `UA`; passing both, or
neither, raises `ValueError`.

**Ports:** `hot_in`, `hot_out`, `cold_in`, `cold_out` &nbsp;·&nbsp;
**Parameters:** `PR_hot`, `PR_cold`, plus one of:
- `effectiveness` — fixed-effectiveness mode
- `UA`, `n_passes` (default 1), `arrangement` (`"counterflow"` /
  `"parallel"` / `"crossflow"` / `"shell_and_tube"` / `"custom"`),
  `correlation` (required iff `arrangement="custom"`) — UA/NTU mode

## Fixed-effectiveness mode (`effectiveness=...`)

A datasheet or an existing exchanger's known performance rating:

$$
C_\text{hot} = \dot m_\text{hot}\,c_{p,\text{hot}}
\qquad
C_\text{cold} = \dot m_\text{cold}\,c_{p,\text{cold}}
\qquad
C_\text{min} = \min(C_\text{hot}, C_\text{cold})
$$

$$
Q = \varepsilon \cdot C_\text{min} \cdot (T_\text{hot,in} - T_\text{cold,in})
$$

$$
P_\text{hot,out} - PR_\text{hot}\,P_\text{hot,in} = 0
\qquad
h_\text{hot,out} - \left(h_\text{hot,in} - \frac{Q}{\dot m_\text{hot}}\right) = 0
$$

(and the mirror pair on the cold side, with `+Q/mdot_cold`). `Q` is not
clamped to `≥ 0` — a network wired with the "hot" side actually colder than
the "cold" side just yields negative `Q` (heat flowing the other way).

## UA/NTU mode (`UA=...`)

Effectiveness comes from the standard effectiveness-NTU relations instead
of being a direct input:

$$
C_r = \frac{C_\text{min}}{C_\text{max}} \qquad NTU = \frac{UA}{C_\text{min}}
$$

$$
\varepsilon_\text{counterflow} = \frac{1-e^{-NTU(1-C_r)}}{1-C_r\,e^{-NTU(1-C_r)}}
\qquad
\varepsilon_\text{parallel} = \frac{1-e^{-NTU(1+C_r)}}{1+C_r}
$$

$$
\varepsilon_\text{crossflow} = 1 - \exp\left[\frac{NTU^{0.22}}{C_r}\left(e^{-C_r\,NTU^{0.78}}-1\right)\right]
$$

`n_passes > 1` chains that many stages in series for these three
arrangements — mathematically a no-op for counterflow/parallel (subdividing
an exact closed-form solution doesn't change it) and not guaranteed to help
for crossflow (an approximation on top of an approximation).

`"shell_and_tube"` is the genuine reversing-header case: `n_passes` is the
number of *shell* passes `N` (each internally 2 tube passes, standard TEMA
"1-2N"), combined via the Bowman/Mueller/Nagle F-correction relation:

$$
\varepsilon_1 = \frac{2}{1+C_r+\sqrt{1+C_r^2}\cdot\dfrac{1+e^{-NTU_1\sqrt{1+C_r^2}}}{1-e^{-NTU_1\sqrt{1+C_r^2}}}}
\qquad NTU_1 = \frac{NTU}{N}
$$

combined in series across `N` shells (see the source docstring for the
full combination formula) — unlike the other three arrangements, `n_passes`
here is a real physical lever: effectiveness strictly increases with `N`,
converging toward the true counterflow limit.

`"custom"` plugs in your own `correlation(NTU, Cr) -> effectiveness`
callable for a plate-fin/finned-tube/vendor correlation.

## `MultiPassHeatExchanger`

`MultiPassHeatExchanger(name, UA, PR_hot, PR_cold, n_passes=1,
arrangement="counterflow", correlation=None)` is a thin subclass of
`HeatExchanger` that always uses UA/NTU mode — no math is duplicated, it
just forwards to `HeatExchanger.__init__(..., UA=UA, ...)`. Use it when you
want that guarantee at the type level (e.g. an `isinstance()` check, or
just to make "this exchanger's rating comes from geometry, not a given
effectiveness" explicit at the call site); otherwise `HeatExchanger(UA=...)`
is equivalent.

Both single-fluid and dual-fluid streams share the network's single fluid
model (no distinct hot/cold fluid types yet — same limitation as every
other component here). Both side's pressure drop is a simple fixed ratio
(not a K-factor loss model); in UA/staged mode it's split evenly in log
space across passes.

For a single-stream, `Q`-driven heater/cooler (no second stream at all),
see [`SimpleHeatExchanger`](simple-heat-exchanger.md).

---
Part of [Heat exchangers & phase change](index.md).
