# SimpleHeatExchanger

<img src="../../_static/diagrams/simple_heat_exchanger.svg" alt="SimpleHeatExchanger diagram" style="max-width:100%">

A two-stream, 0D heat exchanger with a fixed effectiveness rating — the way
a datasheet or an existing exchanger's known performance would be given.

**Ports:** `hot_in`, `hot_out`, `cold_in`, `cold_out` &nbsp;·&nbsp;
**Parameters:** `effectiveness`, `PR_hot`, `PR_cold`

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

For effectiveness derived from geometry (`UA` + flow arrangement) instead of
given directly, see [`MultiPassHeatExchanger`](multi-pass-heat-exchanger.md).

---
Part of [Heat exchangers & phase change](index.md).
