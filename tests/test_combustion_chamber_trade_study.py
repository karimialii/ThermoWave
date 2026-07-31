"""Proof of the combustion-chamber trade study
(docs/tutorials/combustion_chamber_trade_study.py, backing the guide's own
"Engineering trade study" section) -- imports the actual script as a module
instead of duplicating its network-building/ppm-conversion logic a second
time, so the guide's script can't drift from what's tested.
"""

import importlib.util
import math
import sys
from pathlib import Path

import pytest

pytest.importorskip("cantera")

_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "docs" / "tutorials" / "combustion_chamber_trade_study.py"


def _load_trade_study_module():
    spec = importlib.util.spec_from_file_location("combustion_chamber_trade_study", _SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


cts = _load_trade_study_module()


def test_mass_to_mole_normalizes_and_weights_by_molar_mass():
    # A 50/50 mass mix of CO2 (M=44.01) and N2 (M=28.01) is NOT a 50/50
    # mole mix -- the heavier species (CO2) contributes fewer moles per kg.
    Y = {"CO2": 0.5, "N2": 0.5}
    X = cts.mass_to_mole(Y)
    assert math.isclose(sum(X.values()), 1.0, rel_tol=1e-9)
    assert X["N2"] > X["CO2"]  # lighter species: more moles per unit mass
    expected_ratio = cts._MW["CO2"] / cts._MW["N2"]  # actual mechanism MW, not a hand rounded one
    assert math.isclose(X["N2"] / X["CO2"], expected_ratio, rel_tol=1e-9)


def test_dry_ppmvd_drops_water_before_normalizing():
    # Wet vs dry basis must give different answers when H2O is present --
    # otherwise the "dry" conversion isn't doing anything.
    wet_with_water = {"CO": 0.001, "O2": 0.15, "N2": 0.70, "H2O": 0.149}
    wet_no_water = {"CO": 0.001, "O2": 0.15, "N2": 0.849}
    ppm_with_water = cts.dry_ppmvd_at_15pct_o2(wet_with_water, "CO")
    ppm_no_water = cts.dry_ppmvd_at_15pct_o2(wet_no_water, "CO")
    assert not math.isclose(ppm_with_water, ppm_no_water, rel_tol=1e-6)


def test_dry_ppmvd_at_reference_o2_is_a_no_op_correction():
    # When the (dry) measured O2 already equals the 15% reference, the O2
    # correction factor is exactly 1 -- ppm should equal the raw dry mole
    # fraction * 1e6, unchanged.
    X = {"CO": 0.0001, "O2": 0.15, "N2": 0.8499}
    ppm = cts.dry_ppmvd_at_15pct_o2(X, "CO")
    assert math.isclose(ppm, 0.0001 * 1.0e6, rel_tol=1e-6)


def test_trade_study_reproduces_the_documented_temperature_co_conflict():
    # The guide's whole point: a rich-enough primary zone to hit the exit-
    # temperature band leaves CO far over the limit, and a lean-enough
    # primary zone to clear the CO limit overshoots the temperature band --
    # no single point in this small sub-sweep satisfies both. Three points
    # (not the full 9-point sweep) to keep this test's runtime reasonable
    # while still locking in the documented conflict as a real regression
    # check, not just a one-off chat result.
    rows = cts.run_trade_study([0.16, 0.20, 0.30])
    assert all(row["converged"] for row in rows)

    by_f = {round(row["f_primary"], 3): row for row in rows}

    # F=0.16: within the temperature band, but CO is nowhere near the limit.
    ok_T_016, ok_CO_016 = cts.evaluate_constraints(by_f[0.16])
    assert ok_T_016 and not ok_CO_016

    # F=0.30: clears the CO limit, but overshoots the temperature band.
    ok_T_030, ok_CO_030 = cts.evaluate_constraints(by_f[0.3])
    assert ok_CO_030 and not ok_T_030

    # CO drops monotonically as the primary zone leans out across this range.
    assert by_f[0.16]["CO_ppmvd15"] > by_f[0.20]["CO_ppmvd15"] > by_f[0.30]["CO_ppmvd15"]

    # No point in this sub-sweep satisfies both -- the documented conclusion.
    assert not any(all(cts.evaluate_constraints(row)) for row in rows)
