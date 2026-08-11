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

from segs_exergy_benchmark import (  # noqa: E402
    C0, FWH_NAMES, FWH_UA_STAGE0, TTD_U, TURBINES, WATER, _frac, _warm_start, build_network,
)
from thermowave.core.exergy import exergy_report  # noqa: E402
from thermowave.core.solver import SolveResult  # noqa: E402

OUT_DIR = Path(__file__).resolve().parent

SPLIT_TAPS = [
    ("sp1", "HPT1", "HPT2"), ("sp2", "HPT2", "LPT1"), ("sp4", "LPT2", "LPT3"),
    ("sp5", "LPT3", "LPT4"), ("sp6", "LPT4", "LPT5"),
]


def _solve():
    """Same two-stage continuation as segs_exergy_benchmark.py's own
    __main__ -- see that script's build_network() docstring for why."""
    net0, _ = build_network(free=False)
    result0 = net0.solve(
        tol=1e-5, max_iter=300, damping=0.4, step_growth=1.03, progress=False,
        warm_start=_warm_start(net0),
    )
    net, c = build_network(free=True)
    net.validate_topology()
    warm_start = SolveResult(
        converged=True, iterations=0, residual_norm=0.0,
        node_P=result0.node_P, node_h=result0.node_h, node_mdot=result0.node_mdot,
        params={
            **result0.params,
            **{f"{sp}.frac0": _frac(a, b) for sp, a, b in SPLIT_TAPS},
            **{f"{n}.UA": FWH_UA_STAGE0 for n in ("hps1", "hps2", "lps1", "lps2", "lps3")},
        },
        fluid=WATER, node_order=[], components=[],
    )
    result = net.solve(
        tol=1e-5, max_iter=300, damping=0.3, step_growth=1.0, progress=False, warm_start=warm_start,
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


def plot_fwh_ttd(c, state, out_path: Path) -> None:
    """Each heater's solved ttd_u/ttd_l vs. its own 5 K/10 K design target --
    verification that the free-fraction/free-UA design (see the module
    docstring's FEEDWATER HEATERS section) actually lands there."""
    labels = [label for label, _cond, _sub, _valve in FWH_NAMES]
    ttd_u = [c[cond].report_metrics(state)["pinch [K]"] for _l, cond, _s, _v in FWH_NAMES]
    ttd_l = [
        c[sub].report_metrics(state)["T_hot_out [K]"] - c[sub].report_metrics(state)["T_cold_in [K]"]
        for _l, _cond, sub, _v in FWH_NAMES
    ]

    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(x - width / 2, ttd_u, width, label="ttd_u (condenser)", color="#2b6cb0")
    ax.bar(x + width / 2, ttd_l, width, label="ttd_l (subcooler)", color="#c05621")
    ax.axhline(TTD_U, color="#2b6cb0", linestyle="--", linewidth=1, alpha=0.6)
    ax.axhline(10.0, color="#c05621", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_ylabel("Terminal temperature difference [K]")
    ax.set_title("Feedwater heater ttd: solved vs. design target (5 K / 10 K)")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_extraction_split(state, out_path: Path) -> None:
    """Solved extraction fraction at each feedwater-heater tap vs. Table 3's
    own fixed mFeed0 sequence -- the two only coincide if Table 3's design
    point already sits exactly on the ttd_u/ttd_l targets."""
    labels = [sp for sp, _a, _b in SPLIT_TAPS]
    solved = [state.param(f"{sp}.frac0") for sp, _a, _b in SPLIT_TAPS]
    table3 = [_frac(a, b) for _sp, a, b in SPLIT_TAPS]

    x = np.arange(len(labels))
    width = 0.38
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.bar(x - width / 2, solved, width, label="Solved (this benchmark)", color="#2b6cb0")
    ax.bar(x + width / 2, table3, width, label="Table 3 mFeed0", color="#c05621")
    ax.set_ylabel("Extraction fraction [-]")
    ax.set_title("Extraction split: solved (ttd-driven) vs. Table 3's fixed sequence")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.legend()
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    net, c, result = _solve()
    state = result.state()

    plot_turbine_stage_power(c, result, OUT_DIR / "turbine_stage_power.png")
    plot_fwh_ttd(c, state, OUT_DIR / "fwh_ttd.png")
    plot_extraction_split(state, OUT_DIR / "extraction_split.png")

    # Same fuel/product definition as segs_exergy_benchmark.py's own __main__.
    T0, T_HTF_K = 298.15, 390.0 + C0
    gross_mech_w = sum(c[name].report_metrics(state)["power [W]"] for name in TURBINES)

    report = exergy_report(
        result, T0=T0, P0=101325.0,
        fuel=[
            lambda s: sum(
                abs(c[n].report_metrics(s)["power [W]"]) for n in ("eco", "eva", "sup")
            ) * (1.0 - T0 / T_HTF_K),
            lambda s: abs(c["reheater"].report_metrics(s)["power [W]"]) * (1.0 - T0 / T_HTF_K),
        ],
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
