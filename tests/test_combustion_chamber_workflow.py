"""End-to-end proof of the documented 1D discretized combustion-chamber
build workflow (docs/tutorials/building-a-combustion-chamber-model.md),
mirrored step for step so the guide can't silently rot.

T100-scale operating point: 0.8 kg/s air at 4.17 bar / 587 C into the
chamber, 6.8 g/s CH4, chamber exit measured at 918 C -- the fuel/air numbers
and the two temperatures are measured; the liner geometry is an engineering
estimate (no OEM geometry is public); the primary/dilution air split comes
from an independent design rule (a stoichiometric primary zone), NOT from
fitting to the measured exit temperature -- see the guide's own step 7 for
why that distinction matters. This model's predicted exit temperature is
then compared against, not tuned to match, the measured 918 C.
"""

import math

import pytest

pytest.importorskip("cantera")

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

P_IN = 417_000.0            # 4.17 bar
T_IN = 587.0 + 273.15       # 587 C
T_TARGET_OUT = 918.0 + 273.15  # 918 C, measured
MDOT_AIR = 0.8
MDOT_FUEL = 0.0068

# Estimated liner geometry -- no OEM T100 combustor geometry is public, see
# the guide's own disclaimer. Representative small-can-combustor sizing for
# this thermal duty (~340 kW fuel input), not a sourced spec.
D_LINER = 0.10
L_LINER = 0.22
N_STATIONS = 3
L_SEG = L_LINER / N_STATIONS
WALL_THICKNESS = 0.0015
D_ANNULUS = 0.13  # equivalent circular flow-area diameter for the dilution
# annulus -- sized for a plausible ~30 m/s duct velocity carrying most of
# MDOT_AIR, NOT the physical liner-to-casing gap width. A too-small D here
# is the single easiest way to make this network fail to converge: at
# D=0.02 the annulus velocity comes out near 1200 m/s (hypersonic, and the
# Newton solve stalls on it) -- see the guide's own "traps" section.

K_WALL, RHO_WALL, CP_WALL = 15.0, 8220.0, 500.0  # Hastelloy X, approximate
A_RING = math.pi * D_LINER * L_SEG
VOL_RING = math.pi * D_LINER * WALL_THICKNESS * L_SEG
C_RING = RHO_WALL * CP_WALL * VOL_RING
A_COND = math.pi * D_LINER * WALL_THICKNESS

# F_PRIMARY as an INDEPENDENT design input -- a stoichiometric primary zone
# (the conventional flame-stability design target), derived only from the
# boundary inputs (MDOT_FUEL, MDOT_AIR) and a textbook AFR, never from the
# measured exit temperature. Swap in the real hardware's air-split fraction
# if you have it; that's a better independent input, not a different kind
# of fit.
AFR_STOICH_CH4 = 17.2  # kg air / kg fuel, stoichiometric methane combustion
F_PRIMARY = (MDOT_FUEL * AFR_STOICH_CH4) / MDOT_AIR  # ~0.146


def _build(f_primary):
    """Steps 1-4 of the guide: network, splitter, combustor, discretized
    liner (hot-gas/annulus Pipe pairs + ThermalMass rings), final mixer."""
    air = CanteraFluid(name="air", mechanism="gri30.yaml", composition="O2:0.21,N2:0.79")
    network = Network(fluid=air)

    air_src = Source(name="air_src", P=P_IN, T=T_IN, mdot=MDOT_AIR)
    splitter = Junction(
        name="splitter", n_inlets=1, n_outlets=2,
        split_fractions=[f_primary, 1.0 - f_primary],
    )
    comb = Combustor(
        name="comb", PR=0.98, efficiency=1.0,
        mdot_fuel=MDOT_FUEL, fuel="CH4", mechanism="gri30.yaml",
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

    return network, mixer, liners


def _exit_temperature(result, network):
    state = result.state()
    P_out, h_out = state.node("mixer.out0")
    return state.fluid_at("mixer.out0").temperature_ph(P_out, h_out)


def test_check_wiring_is_clean():
    network, _, _ = _build(F_PRIMARY)
    assert network.check_wiring() == []


def test_prediction_from_independent_inputs_is_within_a_few_percent_of_measured():
    # F_PRIMARY here is the stoichiometric-primary-zone design value, derived
    # only from MDOT_FUEL/MDOT_AIR -- never fit to T_TARGET_OUT. This checks
    # the forward PREDICTION lands in the right ballpark of the measurement,
    # not that it matches exactly (it shouldn't, and forcing it to would
    # mean the split was reverse-fit instead of independently chosen -- see
    # the guide's own step 7 for why that distinction matters).
    network, mixer, _ = _build(F_PRIMARY)
    result = network.solve(tol=1e-8, max_iter=800, damping=0.2, verbose=False, progress=False)
    assert result.converged

    T_exit = _exit_temperature(result, network)
    # Measured 901.5 C predicted vs 918 C measured (~-16.5 C, ~2%) when this
    # was last run -- a real discrepancy from this model's own
    # simplifications (no film cooling, estimated geometry, an assumed
    # rather than measured air split), not something to tune away.
    assert math.isclose(T_exit, T_TARGET_OUT, rel_tol=0.03)


def test_exit_temperature_increases_with_more_primary_air():
    # Less dilution air (higher primary fraction) should mix to a hotter
    # exit temperature -- monotone in the physically expected direction.
    network_lo, _, _ = _build(0.16)
    network_hi, _, _ = _build(0.20)
    result_lo = network_lo.solve(tol=1e-7, max_iter=800, damping=0.2, verbose=False, progress=False)
    result_hi = network_hi.solve(tol=1e-7, max_iter=800, damping=0.2, verbose=False, progress=False)
    assert result_lo.converged and result_hi.converged
    assert _exit_temperature(result_lo, network_lo) < _exit_temperature(result_hi, network_hi)


def test_liner_ring_temperatures_decrease_downstream():
    # Each ring loses heat to its own dilution-side annulus and conducts to
    # its downstream neighbor -- the hottest ring should be the first one,
    # closest to the combustor's adiabatic-flame outlet.
    network, _, liners = _build(F_PRIMARY)
    result = network.solve(tol=1e-8, max_iter=800, damping=0.2, verbose=False, progress=False)
    assert result.converged
    state = result.state()
    ring_temps = [state.param(f"{liner.name}.T") for liner in liners]
    assert ring_temps == sorted(ring_temps, reverse=True)
