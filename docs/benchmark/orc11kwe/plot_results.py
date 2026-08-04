"""Grouped bar chart: working-fluid-side vs. oil-side evaporator duty for
the 11 kWe ORC benchmark's two operating points (low / high), annotated
with each case's relative deviation and the paper's own +/-5% heat-balance
closure band as reference lines.

This re-derives Q_wf/Q_oil by running the same network + oil-duty
calculation as orc11kwe_benchmark.py (rather than hardcoding numbers), so
the plot can't silently drift from the benchmark script's actual output.

Run from anywhere; it imports orc11kwe_benchmark.py by path (the script
lives at the repo root and is gitignored -- see that script's docstring):
    python docs/benchmark/orc11kwe/plot_results.py

Writes docs/benchmark/orc11kwe/heat_balance_comparison.png (needs
matplotlib; the coolprop extra is needed transitively via the benchmark
script's Network.solve()).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[3]
BENCHMARK_SCRIPT = REPO_ROOT / "orc11kwe_benchmark.py"
OUT_PNG = Path(__file__).resolve().parent / "heat_balance_comparison.png"

if not BENCHMARK_SCRIPT.exists():
    raise FileNotFoundError(
        f"expected the benchmark script at {BENCHMARK_SCRIPT} (repo root); "
        "it's gitignored there and moves with the repo, so this plot script "
        "must be run against a checkout that still has it."
    )

# Import orc11kwe_benchmark.py as a module *without* re-running its top-level
# side effects twice if it's already been imported -- it prints its own table
# and asserts on import, which is fine (it's the same computation this plot
# reuses), just noisy. We accept the printout so the two stay in lockstep.
spec = importlib.util.spec_from_file_location("orc11kwe_benchmark", BENCHMARK_SCRIPT)
bench = importlib.util.module_from_spec(spec)
sys.modules["orc11kwe_benchmark"] = bench
spec.loader.exec_module(bench)  # runs the benchmark, populates bench.CASES / recomputation below

# Recompute per-case Q_wf/Q_oil/rel_dev exactly as the benchmark script does,
# so we have them as data rather than scraping stdout.
from thermowave.components import Pump, Sink, SimpleEvaporator, Source  # noqa: E402
from thermowave.core import Network  # noqa: E402

labels = []
Q_wf_list = []
Q_oil_list = []
rel_dev_list = []

for label, c in bench.CASES.items():
    net = Network(fluid=bench.R245FA)
    T_sat_cond = bench.R245FA.saturation_temperature(c["p_cond"])
    src = Source(name="src", P=c["p_cond"], T=T_sat_cond - 0.5, mdot=c["mdot_wf"])
    pump = Pump(name="pump", P_out=c["p_evap"], eta=bench.PUMP_ETA_ASSUMED)
    evap = SimpleEvaporator(name="evap", superheat=c["T_superheat"])
    sink = Sink(name="sink")
    for comp in (src, pump, evap, sink):
        net.add_component(comp)
    net.connect(src, "out", pump, "in")
    net.connect(pump, "out", evap, "in")
    net.connect(evap, "out", sink, "in")
    result = net.solve(progress=False)
    assert result.converged

    h_src = result.node_h["src.out"]
    h_evap_out = result.node_h["evap.out"]
    Q_wf = c["mdot_wf"] * (h_evap_out - h_src)

    T_mean_oil_C = (c["T_hf_in_C"] + c["T_hf_out_C"]) / 2.0
    cp_oil = bench.therminol66_cp(T_mean_oil_C)
    Q_oil = c["mdot_hf"] * cp_oil * (c["T_hf_in_C"] - c["T_hf_out_C"])

    rel_dev = (Q_wf - Q_oil) / Q_oil

    labels.append(label.split(" (")[0])  # "low" / "high"
    Q_wf_list.append(Q_wf / 1e3)
    Q_oil_list.append(Q_oil / 1e3)
    rel_dev_list.append(rel_dev)

# --- Plot --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5.5))

x = np.arange(len(labels))
width = 0.32

bars_wf = ax.bar(x - width / 2, Q_wf_list, width, label="Q_wf (working fluid, ThermoWave)",
                  color="#2f6f9f")
bars_oil = ax.bar(x + width / 2, Q_oil_list, width, label="Q_oil (Therminol 66 oil, published)",
                   color="#d68a2c")

ax.bar_label(bars_wf, fmt="%.1f", padding=3, fontsize=9)
ax.bar_label(bars_oil, fmt="%.1f", padding=3, fontsize=9)

# Annotate relative deviation between the two bars in each pair.
for xi, (qwf, qoil, dev) in enumerate(zip(Q_wf_list, Q_oil_list, rel_dev_list)):
    y_top = max(qwf, qoil) + 6.5
    ax.annotate(
        f"{dev:+.1%}",
        xy=(xi, y_top),
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
        color="#1a1a1a",
    )

# +/-5% closure band, drawn relative to the oil-side duty (the reference
# quantity) for each case as short reference bracket lines.
tol = 0.05
for xi, qoil in enumerate(Q_oil_list):
    lo, hi = qoil * (1 - tol), qoil * (1 + tol)
    ax.plot([xi - width, xi + width], [hi, hi], color="gray", linestyle="--", linewidth=1)
    ax.plot([xi - width, xi + width], [lo, lo], color="gray", linestyle="--", linewidth=1)

# Single legend entry for the band (avoid duplicate labels).
ax.plot([], [], color="gray", linestyle="--", linewidth=1, label="paper's own ±5% closure band")

ax.set_xticks(x)
ax.set_xticklabels([f"{lbl}\n({'mdot_hf~1.5 kg/s, T_hf,in=110 C' if lbl == 'low' else 'mdot_hf~3.0 kg/s, T_hf,in=120 C'})"
                     for lbl in labels], fontsize=9)
ax.set_ylabel("Evaporator duty [kW]")
ax.set_title("11 kWe ORC benchmark: evaporator heat-balance closure\n"
              "(ThermoWave working-fluid duty vs. published oil-side duty)")
ax.set_ylim(0, max(Q_wf_list + Q_oil_list) * 1.35)
ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6)

fig.tight_layout()
fig.savefig(OUT_PNG, dpi=200)
print(f"wrote {OUT_PNG}")
