"""Generate the two result plots embedded in segs_exergy_benchmark.md.

Re-runs the benchmark's own build_network()/solve() (importing straight from
segs_exergy_benchmark.py at the repo root -- see that script's own module
docstring for the full model) and saves two PNGs into this directory:

    turbine_stage_power.png  -- ThermoWave vs. published (Table 3) turbine
                                 stage power, grouped bar chart.
    exergy_efficiency.png    -- per-component exergy destruction (bar) and
                                 exergetic efficiency (line, right axis)
                                 from ThermoWave's own exergy_report().

Rerunnable any time the benchmark script or its numbers change:

    python docs/benchmark/segs_exergy/plot_results.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT))

from segs_exergy_benchmark import TURBINES, _warm_start, build_network  # noqa: E402
from thermowave.core.exergy import exergy_report  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent


def _solve():
    net, c = build_network()
    net.validate_topology()
    warm_start = _warm_start(net)
    result = net.solve(
        tol=1e-7, max_iter=300, damping=0.4, step_growth=1.03, progress=False, warm_start=warm_start,
    )
    assert result.converged, "benchmark network failed to converge"
    return net, c, result


def plot_turbine_stage_power(c, result, out_path: Path) -> None:
    state = result.state()
    stages = list(TURBINES.keys())
    tw_kw = [c[name].report_metrics(state)["power [W]"] / 1000.0 for name in stages]
    pub_kw = [TURBINES[name][4] for name in stages]

    x = np.arange(len(stages))
    width = 0.38

    fig, ax = plt.subplots(figsize=(9, 5))
    bars_tw = ax.bar(x - width / 2, tw_kw, width, label="ThermoWave", color="#2b6cb0")
    bars_pub = ax.bar(x + width / 2, pub_kw, width, label="Paper (Table 3)", color="#c05621")

    ax.set_ylabel("Stage power [kW]")
    ax.set_title("Turbine stage power: ThermoWave vs. SAND95-1293 Table 3")
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)

    for bars in (bars_tw, bars_pub):
        for b in bars:
            h = b.get_height()
            ax.annotate(f"{h:,.0f}", (b.get_x() + b.get_width() / 2, h),
                        textcoords="offset points", xytext=(0, 3),
                        ha="center", fontsize=7, color="#333333")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    net, c, result = _solve()
    state = result.state()

    plot_turbine_stage_power(c, result, OUT_DIR / "turbine_stage_power.png")

    # --- exergy report, same fuel/product definition as the main script ---
    q_boiler_w = c["boiler"].report_metrics(state)["power [W]"]
    q_reheat_w = c["reheater"].report_metrics(state)["power [W]"]

    gross_mech_w = sum(c[name].report_metrics(state)["power [W]"] for name in TURBINES)

    report = exergy_report(
        result, T0=298.15, P0=101325.0,
        fuel=[lambda s: q_boiler_w + q_reheat_w],
        product=[lambda s: gross_mech_w
                  - c["condensatePump"].report_metrics(s)["power [W]"]
                  - c["feedPump"].report_metrics(s)["power [W]"]],
    )

    names = [n for n, cost in report.component_costs.items() if cost is not None and cost.epsilon is not None]
    e_d_kw = [report.component_costs[n].E_D / 1000.0 for n in names]
    epsilon = [report.component_costs[n].epsilon * 100.0 for n in names]

    x = np.arange(len(names))
    fig, ax1 = plt.subplots(figsize=(10, 5.5))
    bars = ax1.bar(x, e_d_kw, color="#9b2c2c", alpha=0.85, label="Exergy destruction $E_D$")
    ax1.set_ylabel("Exergy destruction $E_D$ [kW]")
    ax1.set_xticks(x)
    ax1.set_xticklabels(names, rotation=30, ha="right")
    ax1.grid(axis="y", linestyle="--", alpha=0.4)

    ax2 = ax1.twinx()
    ax2.plot(x, epsilon, color="#2c5282", marker="o", linewidth=2, label="Exergetic efficiency $\\epsilon$")
    ax2.set_ylabel("Exergetic efficiency $\\epsilon$ [%]")

    for xi, h in zip(x, e_d_kw):
        ax1.annotate(f"{h:,.0f}", (xi, h), textcoords="offset points", xytext=(0, 3),
                     ha="center", fontsize=7, color="#333333")

    ax2.set_ylim(0, 115)  # headroom so the epsilon line clears the legend box
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left", framealpha=0.95)

    ax1.set_title("Per-component exergy destruction and exergetic efficiency")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "exergy_efficiency.png", dpi=200)
    plt.close(fig)

    print(f"wrote {OUT_DIR / 'turbine_stage_power.png'}")
    print(f"wrote {OUT_DIR / 'exergy_efficiency.png'}")


if __name__ == "__main__":
    main()
