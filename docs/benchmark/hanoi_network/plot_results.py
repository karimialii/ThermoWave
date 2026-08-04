"""Plot nodal head vs. node number for the Hanoi trunk benchmark.

Re-runs `hanoi_trunk_benchmark.py` (imported as a module, so the numbers here
can never drift from what that script actually prints) and plots ThermoWave's
solved heads against the Savic & Walters (1997) reference heads for nodes
2-13, saving the result as nodal_heads.png in this directory.

Run from anywhere:
    python docs/benchmark/hanoi_network/plot_results.py
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_PATH = REPO_ROOT / "hanoi_trunk_benchmark.py"
OUTPUT_PATH = Path(__file__).resolve().parent / "nodal_heads.png"

# Import hanoi_trunk_benchmark.py (repo root, gitignored) as a module rather
# than re-deriving the network here, so this plot can never drift from what
# the benchmark script itself solves and prints.
spec = importlib.util.spec_from_file_location("hanoi_trunk_benchmark", SCRIPT_PATH)
bench = importlib.util.module_from_spec(spec)
sys.modules["hanoi_trunk_benchmark"] = bench
spec.loader.exec_module(bench)

node_ids = bench.NODE_IDS
h_ref = bench.NODE_HEAD_REF_M
h_solved = [
    bench.result.node_P[p_name] / (bench.RHO_NOMINAL * bench.G)
    for p_name in bench.node_pressure_names
]

fig, ax = plt.subplots(figsize=(7, 4.5))

ax.plot(node_ids, h_ref, marker="o", color="#4c72b0", linewidth=2,
        label="Savic & Walters (1997) reference")
ax.plot(node_ids, h_solved, marker="s", color="#dd8452", linewidth=2,
        linestyle="--", label="ThermoWave (Darcy-Weisbach/Colebrook)")

# Highlight nodes 3-5, where deviation is largest (right after the branch
# point at node 3, where flow drops 3x on the same-diameter pipes).
ax.axvspan(3 - 0.4, 5 + 0.4, color="#dd8452", alpha=0.08, zorder=0)
ax.annotate(
    "largest deviation\n(nodes 3-5, after branch point)",
    xy=(4, h_ref[node_ids.index(4)]),
    xytext=(6.3, 66),
    fontsize=9, color="#555555",
    arrowprops=dict(arrowstyle="->", color="#888888", lw=1),
)

ax.set_xlabel("Node")
ax.set_ylabel("Head [m]")
ax.set_title("Hanoi trunk benchmark: nodal head, ThermoWave vs. published reference")
ax.set_xticks(node_ids)
ax.grid(True, alpha=0.3)
ax.legend(loc="upper right", frameon=True)
fig.tight_layout()

fig.savefig(OUTPUT_PATH, dpi=150)
print(f"Saved {OUTPUT_PATH}")
