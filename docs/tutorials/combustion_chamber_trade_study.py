"""Engineering trade study: choosing the primary/dilution air split for the
1D combustion-chamber model in building-a-combustion-chamber-model.md.

The real design question this answers: pick F_PRIMARY (the fraction of the
0.8 kg/s total air routed to the primary/root combustion zone, the rest
being dilution air) to simultaneously

  1. hit the target chamber exit temperature -- 918 C +/- 15 C, the
     measured/design turbine-inlet condition -- and
  2. keep CO emissions under an illustrative regulatory-style limit of
     50 ppmvd @ 15% O2 (dry, oxygen-corrected -- the standard convention
     for reporting gas-turbine emissions, not a raw mole fraction).

Run directly: `python docs/tutorials/combustion_chamber_trade_study.py`
(needs the `cantera` extra: pip install thermowave[cantera]).

IMPORTANT CAVEAT on the NOx numbers this script also computes: they come
from the SAME equilibrium chemistry Combustor already uses for T_out, which
has no notion of finite-rate kinetics or flame residence time/quenching.
Equilibrium NOx overpredicts real engine NOx by roughly 2-3 orders of
magnitude (thousands of ppm here vs. tens of ppm on an actual small gas
turbine) because equilibrium keeps forming NO for as long as the mixture
sits hot, with no "flame residence time" to freeze it early the way a real
flame does. CO happens to survive close to its equilibrium value reasonably
well in practice (which is why it's used here as an order-of-magnitude
design proxy); NOx does not. Treat every NOx number below as illustrative
of the qualitative trend only, never as a compliance-grade prediction -- a
real NOx prediction needs a finite-rate reactor-network model (e.g. a
perfectly-stirred-reactor + plug-flow chain in Cantera) or an empirical
correlation fitted to test data, neither of which this package provides.
"""

from __future__ import annotations

import math

import cantera as ct
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from thermowave.components import (
    Combustor,
    Conduction,
    Convection,
    Junction,
    Pipe,
    Radiation,
    Sink,
    Source,
    ThermalMass,
)
from thermowave.core import Network
from thermowave.fluids import CanteraFluid

# --- Boundary conditions (measured) ----------------------------------------
P_IN = 417_000.0             # 4.17 bar
T_IN = 587.0 + 273.15        # 587 C
MDOT_AIR = 0.8                # kg/s
MDOT_FUEL = 0.0068             # kg/s, 6.8 g/s CH4
T_TARGET_OUT = 918.0 + 273.15  # 918 C, measured chamber exit

# --- Design/engineering constraints -----------------------------------------
T_TARGET_BAND_C = 15.0        # +/- band around the target exit temperature
CO_LIMIT_PPMVD15 = 50.0        # illustrative regulatory-style CO limit

# --- Estimated liner geometry (see the guide's own disclaimer) -------------
D_LINER = 0.10
L_LINER = 0.22
N_STATIONS = 3
L_SEG = L_LINER / N_STATIONS
WALL_THICKNESS = 0.0015
D_ANNULUS = 0.13  # equivalent circular flow-area diameter, NOT the physical
                   # liner-to-casing gap -- see the guide's "traps" section
K_WALL, RHO_WALL, CP_WALL = 15.0, 8220.0, 500.0  # Hastelloy X, approximate
A_RING = math.pi * D_LINER * L_SEG
VOL_RING = math.pi * D_LINER * WALL_THICKNESS * L_SEG
C_RING = RHO_WALL * CP_WALL * VOL_RING
A_COND = math.pi * D_LINER * WALL_THICKNESS

MECHANISM = "gri30.yaml"
_MW = dict(zip(
    ct.Solution(MECHANISM).species_names,
    ct.Solution(MECHANISM).molecular_weights,
))  # g/mol, keyed by species name -- for the mass -> mole conversion below


def build_network(f_primary: float) -> tuple[Network, Combustor, Junction, list[ThermalMass]]:
    """Steps 1-6 of the guide, as a function of the one design variable
    this trade study sweeps."""
    air = CanteraFluid(name="air", mechanism=MECHANISM, composition="O2:0.21,N2:0.79")
    network = Network(fluid=air)

    air_src = Source(name="air_src", P=P_IN, T=T_IN, mdot=MDOT_AIR)
    splitter = Junction(
        name="splitter", n_inlets=1, n_outlets=2,
        split_fractions=[f_primary, 1.0 - f_primary],
    )
    comb = Combustor(
        name="comb", PR=0.98, efficiency=1.0,
        mdot_fuel=MDOT_FUEL, fuel="CH4", mechanism=MECHANISM,
    )
    hotgas = [
        Pipe(name=f"hotgas{i}", L=L_SEG, D=D_LINER, f=0.02, n_elem=1)
        for i in range(N_STATIONS)
    ]
    annulus = [
        Pipe(name=f"annulus{i}", L=L_SEG, D=D_ANNULUS, f=0.02, n_elem=1)
        for i in range(N_STATIONS)
    ]
    liners = [
        ThermalMass(name=f"liner{i}", thermal_capacitance=C_RING, T0=T_IN)
        for i in range(N_STATIONS)
    ]
    mixer = Junction(name="mixer", n_inlets=2, n_outlets=1)
    snk = Sink(name="snk")

    for component in (air_src, splitter, comb, *hotgas, *annulus, *liners, mixer, snk):
        network.add_component(component)

    network.connect(air_src, "out", splitter, "in0")
    network.connect(splitter, "out0", comb, "in")
    network.connect(splitter, "out1", annulus[0], "in")
    network.connect(comb, "out", hotgas[0], "in")
    for i in range(N_STATIONS - 1):
        network.connect(hotgas[i], "out", hotgas[i + 1], "in")
        network.connect(annulus[i], "out", annulus[i + 1], "in")
    network.connect(hotgas[-1], "out", mixer, "in0")
    network.connect(annulus[-1], "out", mixer, "in1")
    network.connect(mixer, "out0", snk, "in")

    for i in range(N_STATIONS):
        network.add_heat_path(
            Radiation(name=f"rad{i}", a=(hotgas[i], "in"), b=liners[i], emissivity=0.7, A=A_RING)
        )
        network.add_heat_path(
            Convection(name=f"conv_hot{i}", a=(hotgas[i], "in"), b=liners[i], h=80.0, A=A_RING)
        )
        network.add_heat_path(
            Convection(name=f"conv_cold{i}", a=liners[i], b=(annulus[i], "in"), h=40.0, A=A_RING)
        )
    for i in range(N_STATIONS - 1):
        network.add_heat_path(
            Conduction(name=f"cond{i}", a=liners[i], b=liners[i + 1], k=K_WALL, A=A_COND, L=L_SEG)
        )

    return network, comb, mixer, liners


def mass_to_mole(mass_fractions: dict[str, float]) -> dict[str, float]:
    """Convert a Cantera-mechanism species mass-fraction dict to mole
    fractions, using the mechanism's own molecular weights."""
    moles = {
        species: y / _MW[species]
        for species, y in mass_fractions.items()
        if species in _MW and y > 0.0
    }
    total = sum(moles.values())
    return {species: n / total for species, n in moles.items()}


def dry_ppmvd_at_15pct_o2(mole_fractions_wet: dict[str, float], species: str) -> float:
    """species's concentration in ppmvd (parts per million, volume, dry) at
    a 15% O2 reference -- the standard gas-turbine emissions-reporting
    convention, so results are comparable to a real regulatory limit
    instead of a raw (and dilution-air-diluted) mole fraction.

    Drops H2O and renormalizes (dry basis), then applies the standard O2
    correction: ppm_corrected = ppm_measured * (20.9 - 15) / (20.9 - O2_measured%).
    """
    dry = {sp: x for sp, x in mole_fractions_wet.items() if sp != "H2O"}
    total_dry = sum(dry.values())
    dry = {sp: x / total_dry for sp, x in dry.items()}
    X_O2_dry_pct = dry.get("O2", 0.0) * 100.0
    X_species_dry = dry.get(species, 0.0)
    correction = (20.9 - 15.0) / max(20.9 - X_O2_dry_pct, 1e-6)
    return X_species_dry * 1.0e6 * correction


def run_trade_study(f_primary_values, warm_start=None):
    """Solve the network at each f_primary (warm-started from the previous
    point -- these are close together, so continuation converges reliably),
    and return one result dict per point."""
    rows = []
    warm = warm_start
    for f_primary in f_primary_values:
        network, comb, mixer, liners = build_network(f_primary)
        result = network.solve(
            tol=1e-8, max_iter=1000, damping=0.2,
            warm_start=warm, verbose=False, progress=False,
        )
        warm = result

        state = result.state()
        P_out, h_out = state.node("mixer.out0")
        T_exit = state.fluid_at("mixer.out0").temperature_ph(P_out, h_out)
        X_mix = mass_to_mole(state.fluid_at("mixer.out0").mass_fractions())
        comb_metrics = comb.report_metrics(state)

        rows.append(dict(
            f_primary=f_primary,
            converged=result.converged,
            T_flame_C=comb_metrics["T_out [K]"] - 273.15,
            T_exit_C=T_exit - 273.15,
            CO_ppmvd15=dry_ppmvd_at_15pct_o2(X_mix, "CO"),
            NOx_ppmvd15=dry_ppmvd_at_15pct_o2(X_mix, "NO"),
        ))
    return rows


def evaluate_constraints(row):
    ok_temperature = abs(row["T_exit_C"] - (T_TARGET_OUT - 273.15)) <= T_TARGET_BAND_C
    ok_co = row["CO_ppmvd15"] <= CO_LIMIT_PPMVD15
    return ok_temperature, ok_co


def plot_trade_study(rows, path):
    f = [r["f_primary"] for r in rows]
    T_flame = [r["T_flame_C"] for r in rows]
    T_exit = [r["T_exit_C"] for r in rows]
    CO = [r["CO_ppmvd15"] for r in rows]
    NOx = [r["NOx_ppmvd15"] for r in rows]
    T_target_C = T_TARGET_OUT - 273.15

    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)

    ax = axes[0]
    ax.plot(f, T_flame, "o-", color="#d62728", label="Primary-zone flame T")
    ax.plot(f, T_exit, "s-", color="#1f77b4", label="Chamber exit T")
    ax.axhspan(
        T_target_C - T_TARGET_BAND_C, T_target_C + T_TARGET_BAND_C,
        color="#1f77b4", alpha=0.12,
        label=f"Target exit T ({T_target_C:.0f} +/- {T_TARGET_BAND_C:.0f} C)",
    )
    ax.set_ylabel("Temperature [C]")
    ax.set_title("Temperature vs primary/dilution air split")
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.6)
    ax.legend(loc="best", fontsize=8)

    ax2 = axes[1]
    ax2.plot(f, CO, "o-", color="#8c564b")
    ax2.axhline(
        CO_LIMIT_PPMVD15, color="#8c564b", linestyle="--", linewidth=1.2,
        label=f"Illustrative CO limit ({CO_LIMIT_PPMVD15:.0f} ppmvd @15%O2)",
    )
    ax2.set_yscale("log")
    ax2.set_ylabel("CO [ppmvd @ 15% O2]\n(log scale)")
    ax2.set_title("CO -- equilibrium estimate (reasonable order-of-magnitude proxy)")
    ax2.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, which="both")
    ax2.legend(loc="best", fontsize=8)

    ax3 = axes[2]
    ax3.plot(f, NOx, "^-", color="#7f7f7f")
    ax3.set_yscale("log")
    ax3.set_xlabel("Primary-zone air fraction  F_PRIMARY  [-]")
    ax3.set_ylabel("\"NOx\" [ppmvd @ 15% O2]\n(log scale)")
    ax3.set_title("Equilibrium NO -- NOT a quantitative emissions prediction (see module docstring)")
    ax3.grid(True, linestyle="--", linewidth=0.5, alpha=0.6, which="both")

    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    f_values = [0.10, 0.12, 0.146, 0.16, 0.18, 0.20, 0.22, 0.26, 0.30]
    rows = run_trade_study(f_values)

    header = f"{'F_PRIMARY':>10} {'T_exit[C]':>10} {'T ok?':>6} {'CO[ppmvd15]':>13} {'CO ok?':>7}"
    print(header)
    print("-" * len(header))
    any_feasible = False
    for row in rows:
        ok_T, ok_CO = evaluate_constraints(row)
        any_feasible = any_feasible or (ok_T and ok_CO)
        print(
            f"{row['f_primary']:>10.3f} {row['T_exit_C']:>10.1f} "
            f"{'yes' if ok_T else 'no':>6} {row['CO_ppmvd15']:>13.1f} "
            f"{'yes' if ok_CO else 'no':>7}"
        )
    print()
    if any_feasible:
        print("At least one F_PRIMARY meets both constraints.")
    else:
        print(
            "No single F_PRIMARY in this sweep meets both the temperature "
            "band and the CO limit -- see the guide's trade-study section "
            "for why, and what a real fix looks like."
        )

    plot_trade_study(rows, "combustion_chamber_trade_study.png")
    print("\nSaved combustion_chamber_trade_study.png")
