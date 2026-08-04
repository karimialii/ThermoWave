"""Generate the comparison plots embedded in sco2_recompression_benchmark.md.

Reuses the already-solved network and exergy report from
`sco2_recompression_benchmark.py` (kept at the repo root -- see that
script's own docstring for why) rather than re-deriving anything. Run
from anywhere; it locates the repo root relative to this file's own
path and adds it to `sys.path` so the benchmark script can be imported
directly:

    python docs/benchmark/sco2_recompression/plot_results.py

Writes two PNGs alongside this file:
  - exergy_comparison.png   component exergetic efficiency, TW vs. published
  - validation_summary.png  predicted vs. published c3 temperature, and the
                             recompression-split discrepancy discussed in the
                             .md's "remaining, precisely-identified limitation"
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(REPO_ROOT))

import sco2_recompression_benchmark as bench  # noqa: E402  (path set up above)

COMPONENTS = ("turb", "rec1", "rec2", "cp1", "cp2")


def plot_exergy_comparison(out_path: Path) -> None:
    tw_eps = [bench.report.component_costs[c].epsilon * 100.0 for c in COMPONENTS]
    pub_eps = [bench.PUBLISHED[c][3] for c in COMPONENTS]

    x = range(len(COMPONENTS))
    width = 0.35

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.bar([i - width / 2 for i in x], tw_eps, width, label="ThermoWave", color="#3b6fa0")
    ax.bar([i + width / 2 for i in x], pub_eps, width, label="Published (Table 6c)",
           color="#c97a3d")

    for i, (tw, pub) in enumerate(zip(tw_eps, pub_eps)):
        ax.text(i - width / 2, tw + 0.6, f"{tw:.1f}", ha="center", fontsize=8)
        ax.text(i + width / 2, pub + 0.6, f"{pub:.1f}", ha="center", fontsize=8)

    ax.set_xticks(list(x))
    ax.set_xticklabels(COMPONENTS)
    ax.set_ylabel("Exergetic efficiency, epsilon [%]")
    ax.set_ylim(75, 102)
    ax.set_title("Component exergetic efficiency: ThermoWave vs.\nPenkuhn & Tsatsaronis (2018), Table 6(c)")
    ax.legend(loc="lower left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def plot_validation_summary(out_path: Path) -> None:
    T3_solved = bench.CO2.temperature_ph(*bench.state.node("rec2.cold_out")) - 273.15
    T3_published = 433.63

    frac_solved = bench.frac0
    frac_calibration = bench.FRAC_CALIBRATION

    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))

    ax = axes[0]
    bars = ax.bar(["ThermoWave\n(predicted)", "Published\n(Table 5c)"],
                   [T3_solved, T3_published], color=["#3b6fa0", "#c97a3d"], width=0.5)
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5,
                f"{b.get_height():.1f} C", ha="center", fontsize=9)
    ax.set_ylabel("c3 temperature [C]")
    ax.set_ylim(0, max(T3_solved, T3_published) * 1.2)
    ax.set_title("Cycle-closure check: predicted vs.\npublished c3 (heater inlet)")
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    bars = ax.bar(["ThermoWave's own\nfree-split solve", "Derivation\ncalibration value"],
                   [frac_solved, frac_calibration], color=["#3b6fa0", "#c97a3d"], width=0.5)
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.01,
                f"{b.get_height():.4f}", ha="center", fontsize=9)
    ax.set_ylabel("Recompression fraction (mdot_cp1 / mdot_total)")
    ax.set_ylim(0, 1.0)
    ax.set_title("Recompression-split discrepancy\n(see .md for the UA-convention explanation)")
    ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle("Where ThermoWave's solve agrees with -- and diverges from -- the paper", y=1.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    exergy_path = HERE / "exergy_comparison.png"
    summary_path = HERE / "validation_summary.png"
    plot_exergy_comparison(exergy_path)
    plot_validation_summary(summary_path)
    print(f"wrote {exergy_path}")
    print(f"wrote {summary_path}")
