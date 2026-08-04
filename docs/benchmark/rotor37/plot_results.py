"""Sweep pressure ratio through ThermoWave's solver and plot the result.

This is a companion to ``rotor37_benchmark.py`` (repo root, gitignored). That
script is a single-point machine-precision check: it hand-calculates the
exit state implied by NASA Rotor 37's published (PR, eta_s) and confirms
ThermoWave's SimpleCompressor reproduces it. A single point doesn't make an
interesting plot on its own.

What *is* useful is seeing where that design point sits on the compressor's
broader operating range. This script builds and solves a fresh ThermoWave
network -- not the textbook formula -- at each of several pressure ratios
around the design point (holding eta_s, inlet state, and mass flow fixed at
their published design values), and plots the solved shaft power and exit
temperature vs. PR, marking NASA's published design point (PR=2.106) on the
curve.

Run directly from anywhere (no CLI args):

    python docs/benchmark/rotor37/plot_results.py

Requires matplotlib (already a dev dependency of this repo).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from thermowave.components import SimpleCompressor, Sink, Source
from thermowave.core import Network
from thermowave.fluids import IdealGasFluid

# --- Published NASA Rotor 37 design point (same as rotor37_benchmark.py) ---
MDOT_DESIGN = 20.19      # kg/s
PR_DESIGN = 2.106        # total-to-total pressure ratio
ETA_DESIGN = 0.877       # adiabatic (isentropic) efficiency
P_IN = 101_325.0         # Pa
T_IN = 288.15            # K
GAMMA_AIR = 1.4


def solve_at_pr(pr: float) -> tuple[float, float]:
    """Build and solve a fresh ThermoWave network at the given PR.

    Returns (shaft power [W], exit temperature [K]), both taken from the
    solver's converged network state -- not the hand-calc formula.
    """
    air = IdealGasFluid(name="air", R=287.05, cp=1005.0)
    net = Network(fluid=air)

    src = Source(name="src", P=P_IN, T=T_IN, mdot=MDOT_DESIGN)
    comp = SimpleCompressor(name="comp", PR=pr, eta_s=ETA_DESIGN, gamma=GAMMA_AIR)
    sink = Sink(name="sink")

    net.add_component(src)
    net.add_component(comp)
    net.add_component(sink)
    net.connect(src, "out", comp, "in")
    net.connect(comp, "out", sink, "in")

    result = net.solve(progress=False)
    if not result.converged:
        raise RuntimeError(f"network did not converge at PR={pr}")

    metrics = comp.report_metrics(result.state())
    return metrics["power [W]"], metrics["T_out [K]"]


def main() -> None:
    pr_values = np.linspace(1.5, 2.5, 41)
    powers = np.empty_like(pr_values)
    t_outs = np.empty_like(pr_values)

    for i, pr in enumerate(pr_values):
        power, t_out = solve_at_pr(float(pr))
        powers[i] = power / 1e6  # MW
        t_outs[i] = t_out

    design_power, design_t_out = solve_at_pr(PR_DESIGN)
    design_power /= 1e6

    fig, (ax_power, ax_temp) = plt.subplots(1, 2, figsize=(11, 4.5))

    ax_power.plot(pr_values, powers, color="#1f6f8b", lw=2, label="ThermoWave (solved)")
    ax_power.scatter(
        [PR_DESIGN], [design_power], color="#d1495b", zorder=5, s=60,
        label=f"Rotor 37 design point\n(PR={PR_DESIGN}, {design_power:.3f} MW)",
    )
    ax_power.set_xlabel("Pressure ratio, PR [-]")
    ax_power.set_ylabel("Shaft power [MW]")
    ax_power.set_title("Shaft power vs. pressure ratio")
    ax_power.grid(True, alpha=0.3)
    ax_power.legend(loc="upper left", fontsize=8)

    ax_temp.plot(pr_values, t_outs, color="#1f6f8b", lw=2, label="ThermoWave (solved)")
    ax_temp.scatter(
        [PR_DESIGN], [design_t_out], color="#d1495b", zorder=5, s=60,
        label=f"Rotor 37 design point\n(PR={PR_DESIGN}, {design_t_out:.1f} K)",
    )
    ax_temp.set_xlabel("Pressure ratio, PR [-]")
    ax_temp.set_ylabel("Exit temperature [K]")
    ax_temp.set_title("Exit temperature vs. pressure ratio")
    ax_temp.grid(True, alpha=0.3)
    ax_temp.legend(loc="upper left", fontsize=8)

    fig.suptitle(
        "NASA Rotor 37: ThermoWave solver sweep at fixed eta_s=0.877, mdot=20.19 kg/s",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    out_path = Path(__file__).parent / "pr_sweep.png"
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
