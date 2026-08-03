"""Benchmark: the water/steam Rankine cycle of a 30 MWe SEGS parabolic-
trough solar plant, at its 100%-solar design point, vs. the ORIGINAL
PUBLISHED source: F. Lippke, "Simulation of the Part-Load Behavior of a
30 MWe SEGS Plant", SAND95-1293, Sandia National Laboratories, 1995
(https://www.osti.gov/biblio/95571) -- Table 3 (every component's design
UA/mass-flow/pressure/efficiency) and Fig. 4 (the design heat balance) give
every number this script's topology and component parameters come from.

This plant is also the subject of a TESPy exergy-analysis example repo
(docs/benchmark/Tespy/SEGS_exergy-main), run this session for cross-
checking purposes -- see segs_exergy_benchmark.md for what that run did
and didn't show. Per this benchmark's own scope, only the ORIGINAL PAPER's
numbers are used as the validation target below, not TESPy's.

SCOPE: the water/steam power cycle only, not the solar field/HTF loop
------------------------------------------------------------------------
The oil-heated boiler (economizer+evaporator+superheater) and reheater
cross the system boundary between the solar/HTF side and the water/steam
side. Both endpoints of each are fully specified by the paper's own Fig. 4
design heat balance (T *and* P), so -- exactly like the sco2_recompression
benchmark's own treatment of its external heater/cooler -- there is
nothing to solve there: they're boundary conditions. This network starts
at a Source (boiler outlet steam, 100 bar/371 C) and ends at a Sink
(feedwater returning to the boiler); the reheater is a second, independent
Source/Sink boundary pair (hpt2 exhaust -> Sink; 17.1 bar/371 C -> Source
-> lpt1).

WHY THE FEEDWATER HEATERS AREN'T MultiPassHeatExchanger
------------------------------------------------------------------------
Every extraction feeding a regenerative feedwater heater here is
condensing steam -- exactly the case SimpleCondenser's own docstring
warns about: "cp is effectively infinite during a constant-pressure phase
change, so an effectiveness-NTU (C = mdot*cp) framework can't represent
condensation." Confirmed the hard way this session: wiring these heaters
as MultiPassHeatExchanger reliably crashed CoolProp's cp() lookup, because
several of the real extraction states here (this plant's own LP turbine
exhausts, in particular) land ON or inside the two-phase dome, where cp
isn't well-defined. So each heater here is instead a boundary-crossing
PAIR, one component per side, matching how a real feedwater heater's own
two jobs are actually different physics:

  - hot side: SimpleCondenser(outlet_quality=0) -- the extraction steam
    condenses fully to saturated liquid at its own (solved) mdot/h_in;
    its duty is a genuine, solved result, not assumed.
  - cold side: SimpleEvaporator(duty=<precomputed constant>, despite the
    "evaporator" name -- duty mode just adds a fixed Q to a single stream,
    the feedwater here staying subcooled liquid throughout) -- its duty is
    precomputed ONCE outside the solve (see Q_FIXED_W below), from Table 3's
    own design mass flow and a standard feedwater-heater terminal-
    temperature-difference (TTD = 5 K below the extraction's saturation
    temperature) assumption, using only safe, single-phase-liquid CoolProp
    calls (enthalpy_pt at subcooled states) -- never a cp() call anywhere
    near the dome.

The two sides' duties are therefore NOT forced to match exactly (no shared
UA/effectiveness ties them) -- see the printed "energy balance" section,
which reports both and their gap explicitly, rather than hiding it.

TOPOLOGY (7 turbine stages, 5 regenerative feedwater heaters with a
condensate-cascade drain train, one open feedwater heater/deaerator) --
matches Table 3's own component list and TESPy's SEGS.py wiring:

    Source(100 bar, 371 C, 38.6415 kg/s)
      -> hpt1 -(33.61 bar)-> sp1 --extraction--> hpPreheater2_hot.in
                              \\--main--> hpt2 -(18.58 bar)-> sp2
                                                    --extraction--> merge_hppre1.in0
                                                     \\--main--> Sink (reheater in)

    Source(17.1 bar, 371 C, 32.8068 kg/s)             [reheater out]
      -> lpt1 -(7.977 bar)-> sp3 --extraction(direct steam)--> FWT.in0
                              \\--main--> lpt2 -(2.728 bar)-> sp4
                                          --extraction--> lpPreheater3_hot.in
                                           \\--main--> lpt3 -(0.9625 bar)-> sp5
                                                       --extraction--> merge_lppre2.in0
                                                        \\--main--> lpt4 -(0.2868 bar)-> sp6
                                                                    --extraction--> merge_lppre1.in0
                                                                     \\--main--> lpt5 -(0.08 bar)-> merge_cond.in0

    Drain cascade (each heater's condensed extraction steam falls to the
    next lower-pressure heater, exactly like a real regenerative train):
      hpPreheater2_hot.out -> merge_hppre1.in1  (-> hpPreheater1_hot.in)
      hpPreheater1_hot.out -> FWT.in2                  (deaerator, open FWH)
      lpPreheater3_hot.out -> merge_lppre2.in1  (-> lpPreheater2_hot.in)
      lpPreheater2_hot.out -> merge_lppre1.in1  (-> lpPreheater1_hot.in)
      lpPreheater1_hot.out -> merge_cond.in1           (-> condenser.in)

    Feedwater (cold) train:
      condenser.out -> condensatePump -> lpPreheater1_cold -> ... ->
      lpPreheater3_cold.out -> FWT.in1 (passing feedwater; FWT also
      receives sp3's direct steam extraction FWT.in0, and hpPreheater1_hot's
      drain FWT.in2 -- a genuine 3-stream open feedwater heater/deaerator,
      same as TESPy's own 3-inlet "Feedwater tank" Merge)
      -> feedPump -> hpPreheater1_cold -> hpPreheater2_cold ->
      Sink (boiler feed)

Every Junction split fraction is FIXED (not solved for) -- Table 3's own
turbine mFeed0 sequence gives the exact design-point extraction at each
tap directly (e.g. hpt1's 38.6415 -> 35.7326 kg/s IS the 2.9089 kg/s design
extraction to hpPreheater2), so there's no free-parameter/Controller
machinery needed here, unlike the sco2_recompression benchmark's
recompression split (whose fraction genuinely wasn't known in advance).
Mass conservation is a clean check of the topology, not an assumption: with
these fixed fractions, the drain cascade re-derives the full 38.6415 kg/s
back at the FWT outlet, with zero unaccounted loss -- see the .md.

Run directly (needs the coolprop extra: pip install thermowave[coolprop]):
    python docs/benchmark/segs_exergy/segs_exergy_benchmark.py
"""

from __future__ import annotations

from thermowave.components import (
    Junction,
    Pump,
    Sink,
    SimpleCondenser,
    SimpleEvaporator,
    Source,
    SteamTurbine,
)
from thermowave.core import Network
from thermowave.core.solver import SolveResult
from thermowave.fluids import CoolPropFluid

WATER = CoolPropFluid(name="Water")
BAR = 1.0e5
C0 = 273.15
TTD_K = 5.0  # feedwater-heater terminal temperature difference (standard design assumption)

# --- Table 3 (SAND95-1293, p.7): turbine stages -----------------------
# (mFeed0 [kg/s], pFeed0 [bar], pDrain0 [bar], etas0 [-], genPower [kW])
TURBINES = {
    "hpt1": (38.6415, 100.0000, 33.6100, 0.8376, 7643.538),
    "hpt2": (35.7326, 33.6100, 18.5800, 0.8463, 3480.355),
    "lpt1": (32.8068, 17.1000, 7.9770, 0.8623, 5730.183),
    "lpt2": (30.7936, 7.9770, 2.7280, 0.9170, 6697.607),
    "lpt3": (29.0324, 2.7280, 0.9625, 0.9352, 5045.830),
    "lpt4": (27.4158, 0.9625, 0.2868, 0.8800, 4505.279),
    "lpt5": (26.6117, 0.2868, 0.0800, 0.6445, 2979.683),
}

PR_HOT, PR_COLD = 0.97, 0.98


def _frac(upstream: str, downstream: str) -> float:
    """Fixed extraction fraction, derived from TURBINES' own mFeed0
    sequence: fraction of the upstream stage's flow that peels off to a
    heater's hot side, everything else continues to the next stage."""
    m_up = TURBINES[upstream][0]
    m_down = TURBINES[downstream][0]
    return (m_up - m_down) / m_up


F_HPPRE2 = _frac("hpt1", "hpt2")
F_HPPRE1 = _frac("hpt2", "lpt1")
F_FWT_DIRECT = _frac("lpt1", "lpt2")
F_LPPRE3 = _frac("lpt2", "lpt3")
F_LPPRE2 = _frac("lpt3", "lpt4")
F_LPPRE1 = _frac("lpt4", "lpt5")

# --- Precomputed cold-side duties (see module docstring) --------------
# (extraction pressure [bar], cold-side design mass flow [kg/s] -- Table 3's
# own mCold0). Chained: each stage's assumed T_out becomes the next stage's
# assumed T_in, following the real feedwater path, ending near Table 3's own
# hpPreheater2 exit condition (close to the 371 C/100 bar boiler inlet after
# the (unmodeled) economizer/superheater picks up the rest).
_FWH_STAGES = [
    ("lpPreheater1", 0.2868, 31.0949),
    ("lpPreheater2", 0.9625, 31.0949),
    ("lpPreheater3", 2.7280, 31.0949),
    ("hpPreheater1", 18.5800, 38.8435),
    ("hpPreheater2", 33.6100, 38.8435),
]


def _fwh_cold_duties() -> dict[str, float]:
    duties = {}
    T_in_c = 45.0  # condensatePump outlet, design condenser hotwell temperature
    for name, P_bar, mdot_cold in _FWH_STAGES:
        if name == "hpPreheater1":
            T_in_c = 170.0  # FWT/deaerator outlet, ~Tsat(lpt1 exhaust, 7.977 bar)
        P_liquid = (9.0 if P_bar < 10.0 else 104.1) * BAR
        T_sat_c = WATER.saturation_temperature(P_bar * BAR) - C0
        T_out_c = T_sat_c - TTD_K
        h_in = WATER.enthalpy_pt(P_liquid, T_in_c + C0)
        h_out = WATER.enthalpy_pt(P_liquid, T_out_c + C0)
        duties[name] = mdot_cold * (h_out - h_in)
        T_in_c = T_out_c
    return duties


Q_FIXED_W = _fwh_cold_duties()


def build_network() -> tuple[Network, dict]:
    net = Network(fluid=WATER)

    src_boiler = Source(
        name="boiler_out", P=TURBINES["hpt1"][1] * BAR, T=371.0 + C0,
        mdot=TURBINES["hpt1"][0],
    )
    hpt1 = SteamTurbine(name="hpt1", P_out=TURBINES["hpt1"][2] * BAR, eta_s=TURBINES["hpt1"][3])
    sp1 = Junction(name="sp1", n_inlets=1, n_outlets=2, split_fractions=[F_HPPRE2, 1.0 - F_HPPRE2])

    hpt2 = SteamTurbine(name="hpt2", P_out=TURBINES["hpt2"][2] * BAR, eta_s=TURBINES["hpt2"][3])
    sp2 = Junction(name="sp2", n_inlets=1, n_outlets=2, split_fractions=[F_HPPRE1, 1.0 - F_HPPRE1])
    sink_reheat_in = Sink(name="sink_reheat_in")

    src_reheat = Source(
        name="reheat_out", P=TURBINES["lpt1"][1] * BAR, T=371.0 + C0, mdot=TURBINES["lpt1"][0],
    )
    lpt1 = SteamTurbine(name="lpt1", P_out=TURBINES["lpt1"][2] * BAR, eta_s=TURBINES["lpt1"][3])
    sp3 = Junction(name="sp3", n_inlets=1, n_outlets=2, split_fractions=[F_FWT_DIRECT, 1.0 - F_FWT_DIRECT])

    lpt2 = SteamTurbine(name="lpt2", P_out=TURBINES["lpt2"][2] * BAR, eta_s=TURBINES["lpt2"][3])
    sp4 = Junction(name="sp4", n_inlets=1, n_outlets=2, split_fractions=[F_LPPRE3, 1.0 - F_LPPRE3])

    lpt3 = SteamTurbine(name="lpt3", P_out=TURBINES["lpt3"][2] * BAR, eta_s=TURBINES["lpt3"][3])
    sp5 = Junction(name="sp5", n_inlets=1, n_outlets=2, split_fractions=[F_LPPRE2, 1.0 - F_LPPRE2])

    lpt4 = SteamTurbine(name="lpt4", P_out=TURBINES["lpt4"][2] * BAR, eta_s=TURBINES["lpt4"][3])
    sp6 = Junction(name="sp6", n_inlets=1, n_outlets=2, split_fractions=[F_LPPRE1, 1.0 - F_LPPRE1])

    lpt5 = SteamTurbine(name="lpt5", P_out=TURBINES["lpt5"][2] * BAR, eta_s=TURBINES["lpt5"][3])

    merge_cond = Junction(name="merge_cond", n_inlets=2, n_outlets=1)
    condenser = SimpleCondenser(name="condenser", PR=1.0, outlet_quality=0.0)
    condensate_pump = Pump(name="condensatePump", P_out=9.0 * BAR, eta=0.7125)

    lp1_hot = SimpleCondenser(name="lpPreheater1_hot", PR=PR_HOT, outlet_quality=0.0)
    lp1_cold = SimpleEvaporator(name="lpPreheater1_cold", PR=PR_COLD, duty=Q_FIXED_W["lpPreheater1"])
    merge_lppre1 = Junction(name="merge_lppre1", n_inlets=2, n_outlets=1)

    lp2_hot = SimpleCondenser(name="lpPreheater2_hot", PR=PR_HOT, outlet_quality=0.0)
    lp2_cold = SimpleEvaporator(name="lpPreheater2_cold", PR=PR_COLD, duty=Q_FIXED_W["lpPreheater2"])
    merge_lppre2 = Junction(name="merge_lppre2", n_inlets=2, n_outlets=1)

    lp3_hot = SimpleCondenser(name="lpPreheater3_hot", PR=PR_HOT, outlet_quality=0.0)
    lp3_cold = SimpleEvaporator(name="lpPreheater3_cold", PR=PR_COLD, duty=Q_FIXED_W["lpPreheater3"])

    fwt = Junction(name="FWT", n_inlets=3, n_outlets=1)
    feed_pump = Pump(name="feedPump", P_out=104.1 * BAR, eta=0.7125)

    hp1_hot = SimpleCondenser(name="hpPreheater1_hot", PR=PR_HOT, outlet_quality=0.0)
    hp1_cold = SimpleEvaporator(name="hpPreheater1_cold", PR=PR_COLD, duty=Q_FIXED_W["hpPreheater1"])
    merge_hppre1 = Junction(name="merge_hppre1", n_inlets=2, n_outlets=1)

    hp2_hot = SimpleCondenser(name="hpPreheater2_hot", PR=PR_HOT, outlet_quality=0.0)
    hp2_cold = SimpleEvaporator(name="hpPreheater2_cold", PR=PR_COLD, duty=Q_FIXED_W["hpPreheater2"])
    sink_boiler_feed = Sink(name="sink_boiler_feed")

    components = dict(
        boiler_out=src_boiler, hpt1=hpt1, sp1=sp1, hpt2=hpt2, sp2=sp2,
        sink_reheat_in=sink_reheat_in, reheat_out=src_reheat, lpt1=lpt1, sp3=sp3,
        lpt2=lpt2, sp4=sp4, lpt3=lpt3, sp5=sp5, lpt4=lpt4, sp6=sp6, lpt5=lpt5,
        merge_cond=merge_cond, condenser=condenser, condensatePump=condensate_pump,
        lpPreheater1_hot=lp1_hot, lpPreheater1_cold=lp1_cold, merge_lppre1=merge_lppre1,
        lpPreheater2_hot=lp2_hot, lpPreheater2_cold=lp2_cold, merge_lppre2=merge_lppre2,
        lpPreheater3_hot=lp3_hot, lpPreheater3_cold=lp3_cold,
        FWT=fwt, feedPump=feed_pump,
        hpPreheater1_hot=hp1_hot, hpPreheater1_cold=hp1_cold, merge_hppre1=merge_hppre1,
        hpPreheater2_hot=hp2_hot, hpPreheater2_cold=hp2_cold,
        sink_boiler_feed=sink_boiler_feed,
    )
    for c in components.values():
        net.add_component(c)

    net.connect(src_boiler, "out", hpt1, "in")
    net.connect(hpt1, "out", sp1, "in0")
    net.connect(sp1, "out0", hp2_hot, "in")
    net.connect(sp1, "out1", hpt2, "in")

    net.connect(hpt2, "out", sp2, "in0")
    net.connect(sp2, "out0", merge_hppre1, "in0")
    net.connect(sp2, "out1", sink_reheat_in, "in")

    net.connect(src_reheat, "out", lpt1, "in")
    net.connect(lpt1, "out", sp3, "in0")
    net.connect(sp3, "out0", fwt, "in0")
    net.connect(sp3, "out1", lpt2, "in")

    net.connect(lpt2, "out", sp4, "in0")
    net.connect(sp4, "out0", lp3_hot, "in")
    net.connect(sp4, "out1", lpt3, "in")

    net.connect(lpt3, "out", sp5, "in0")
    net.connect(sp5, "out0", merge_lppre2, "in0")
    net.connect(sp5, "out1", lpt4, "in")

    net.connect(lpt4, "out", sp6, "in0")
    net.connect(sp6, "out0", merge_lppre1, "in0")
    net.connect(sp6, "out1", lpt5, "in")

    net.connect(lpt5, "out", merge_cond, "in0")
    net.connect(lp1_hot, "out", merge_cond, "in1")
    net.connect(merge_cond, "out0", condenser, "in")
    net.connect(condenser, "out", condensate_pump, "in")

    net.connect(condensate_pump, "out", lp1_cold, "in")
    net.connect(merge_lppre1, "out0", lp1_hot, "in")
    net.connect(lp2_hot, "out", merge_lppre1, "in1")

    net.connect(lp1_cold, "out", lp2_cold, "in")
    net.connect(merge_lppre2, "out0", lp2_hot, "in")
    net.connect(lp3_hot, "out", merge_lppre2, "in1")

    net.connect(lp2_cold, "out", lp3_cold, "in")

    net.connect(lp3_cold, "out", fwt, "in1")
    net.connect(hp1_hot, "out", fwt, "in2")
    net.connect(fwt, "out0", feed_pump, "in")

    net.connect(feed_pump, "out", hp1_cold, "in")
    net.connect(merge_hppre1, "out0", hp1_hot, "in")
    net.connect(hp2_hot, "out", merge_hppre1, "in1")

    net.connect(hp1_cold, "out", hp2_cold, "in")
    net.connect(hp2_cold, "out", sink_boiler_feed, "in")

    return net, components


def _warm_start(net: Network) -> SolveResult:
    """A hand-built warm start for the cold (feedwater) chain: without it,
    Newton's own flat/propagated default guess for these nodes can land
    inside or right at the edge of the two-phase dome for at least one
    intermediate node (a purely numerical hazard of guess propagation
    compounding a fixed +4e5 J/kg swing per hop across several chained
    components -- see this session's own debugging, summarized in the
    .md), which crashes CoolProp's cp()/temperature_ph() before Newton
    even gets a first step. These are physically reasonable subcooled-
    liquid temperatures (see _fwh_cold_duties()'s own T chain above), not
    a fitted answer -- only the ENERGY level matters for a warm start, not
    precision.
    """
    node_T_c = {
        "condensatePump.out": 45.0, "lpPreheater1_cold.in": 45.0, "lpPreheater1_cold.out": 63.0,
        "lpPreheater2_cold.in": 63.0, "lpPreheater2_cold.out": 93.0,
        "lpPreheater3_cold.in": 93.0, "lpPreheater3_cold.out": 125.0,
        "FWT.out0": 125.0, "feedPump.out": 170.0,
        "hpPreheater1_cold.in": 170.0, "hpPreheater1_cold.out": 204.0,
        "hpPreheater2_cold.in": 204.0, "hpPreheater2_cold.out": 235.0,
    }
    node_h = {net._canonical(raw): WATER.enthalpy_pt(9.0 * BAR, T_c + C0) for raw, T_c in node_T_c.items()}
    return SolveResult(
        converged=True, iterations=0, residual_norm=0.0,
        node_P={}, node_h=node_h, node_mdot={}, params={},
        fluid=WATER, node_order=[], components=[],
    )


if __name__ == "__main__":
    net, c = build_network()
    net.validate_topology()
    warm_start = _warm_start(net)
    result = net.solve(
        tol=1e-7, max_iter=300, damping=0.4, step_growth=1.03, progress=False, warm_start=warm_start,
    )
    print(f"converged: {result.converged}, iterations: {result.iterations}")
    result.print_report()

    state = result.state()

    print("\n--- Turbine stage power vs. Table 3 genPower ---")
    turb_objs = {n: c[n] for n in TURBINES}
    tw_power_w_total = 0.0
    pub_power_kw_total = 0.0
    for name, (mFeed0, pFeed0, pDrain0, etas0, genPower_kw) in TURBINES.items():
        tw_w = turb_objs[name].report_metrics(state)["power [W]"]
        tw_power_w_total += tw_w
        pub_power_kw_total += genPower_kw
        d = (tw_w / 1000.0 - genPower_kw) / genPower_kw * 100.0
        print(f"{name:6s}  TW={tw_w/1000.0:9.2f} kW   paper={genPower_kw:9.2f} kW   d={d:+6.2f}%")
    tw_power_kw_total = tw_power_w_total / 1000.0
    print(f"{'TOTAL':6s}  TW={tw_power_kw_total:9.2f} kW   paper={pub_power_kw_total:9.2f} kW   "
          f"d={(tw_power_kw_total - pub_power_kw_total) / pub_power_kw_total * 100.0:+6.2f}%")

    print("\n--- Mass balance closure ---")
    mdot_boiler = TURBINES["hpt1"][0]
    mdot_feed_return = state.mdot("hpPreheater2_cold.out")
    print(f"boiler outlet mdot:         {mdot_boiler:.4f} kg/s")
    print(f"feedwater return to boiler: {mdot_feed_return:.4f} kg/s "
          f"(should match, closing the drain cascade)")

    print("\n--- Feedwater heater energy balance (hot side solved, cold side fixed-duty) ---")
    for name in ("hpPreheater1", "hpPreheater2", "lpPreheater1", "lpPreheater2", "lpPreheater3"):
        q_hot = c[f"{name}_hot"].report_metrics(state)["power [W]"]  # negative: heat rejected
        q_cold = c[f"{name}_cold"].report_metrics(state)["power [W]"]
        gap_pct = (abs(q_hot) - q_cold) / q_cold * 100.0
        print(f"{name:14s} Q_hot={-q_hot/1000:9.1f} kW   Q_cold={q_cold/1000:9.1f} kW   gap={gap_pct:+6.1f}%")

    condensate_pump_kw = c["condensatePump"].report_metrics(state)["power [W]"] / 1000.0
    feed_pump_kw = c["feedPump"].report_metrics(state)["power [W]"] / 1000.0
    print(f"\ncondensatePump power: {condensate_pump_kw:.2f} kW  (Table 2 published: 190 kW*)")
    print(f"feedPump power:       {feed_pump_kw:.2f} kW  (Table 2 published: 880 kW*, "
          f"combined condensate+feed = 1070 kW)")

    P_boiler_in, h_boiler_in = TURBINES["hpt1"][1] * BAR, WATER.enthalpy_pt(TURBINES["hpt1"][1] * BAR, 371.0 + C0)
    P_feed_ret, h_feed_ret = state.node("hpPreheater2_cold.out")
    q_boiler_w = mdot_boiler * (h_boiler_in - h_feed_ret)

    P_rh_in, h_rh_in = state.node("sink_reheat_in.in")
    q_reheat_w = TURBINES["lpt1"][0] * (
        WATER.enthalpy_pt(TURBINES["lpt1"][1] * BAR, 371.0 + C0) - h_rh_in
    )

    gross_mech_w = tw_power_w_total
    gross_eff_steam_cycle = gross_mech_w / (q_boiler_w + q_reheat_w) * 100.0
    print(f"\nboiler duty:   {q_boiler_w/1e6:.3f} MW")
    print(f"reheater duty: {q_reheat_w/1e6:.3f} MW")
    print(f"steam-cycle gross mechanical efficiency (turbine power / (boiler+reheat duty)): "
          f"{gross_eff_steam_cycle:.1f}%")
    print("(paper's headline 38.2% gross ELECTRIC efficiency uses a wider boundary -- "
          "solar-absorbed power, including the HTF loop this benchmark doesn't model; "
          "see the .md for the correct comparison.)")

    generator_eff = 0.97  # Fig. 4's own kW -> kWe factor
    net_kw = (gross_mech_w * generator_eff - (condensate_pump_kw + feed_pump_kw) * 1000.0) / 1000.0
    print(f"\ngross mechanical power: {gross_mech_w/1000:.1f} kW")
    print(f"gross electrical power (x{generator_eff} generator): {gross_mech_w*generator_eff/1000:.1f} kW "
          f"(paper Fig. 4: 34985 kWe)")
    print(f"net electrical power (- condensate/feed pump only, no HTF/cooling-tower pumps modeled): "
          f"{net_kw:.1f} kW (paper's full-plant net: 31400 kWe)")

    print("\n--- ThermoWave exergy analysis (own value-add; the paper has no exergy data) ---")
    from thermowave.core.exergy import exergy_report

    report = exergy_report(
        result, T0=298.15, P0=101325.0,
        fuel=[lambda s: q_boiler_w + q_reheat_w],
        product=[lambda s: gross_mech_w
                  - c["condensatePump"].report_metrics(s)["power [W]"]
                  - c["feedPump"].report_metrics(s)["power [W]"]],
    )
    report.print_report()
