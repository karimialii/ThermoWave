import math

import pytest

from thermowave.maps.characteristic_map import CharacteristicMap

_MAP_PATH = "tests/fixtures/simple_compressor_map.cop"


def test_default_factors_match_from_file_with_no_overrides():
    with_defaults = CharacteristicMap.from_file(_MAP_PATH)
    explicit_none = CharacteristicMap.from_file(_MAP_PATH, factor_overrides=None)
    A, B = with_defaults.mid_speed(), 1.0
    assert math.isclose(
        with_defaults.pressure_ratio(A, B), explicit_none.pressure_ratio(A, B)
    )


def test_b_fact_override_scales_corrected_mass_flow_axis():
    base = CharacteristicMap.from_file(_MAP_PATH)
    scaled = CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"B_fact": 2.0})
    A = base.mid_speed()
    B = 1.0
    # scaled map's B-axis values are all 2x the base map's, so querying at
    # 2*B on the scaled map should match querying at B on the base map.
    assert math.isclose(
        scaled.pressure_ratio(A, 2.0 * B), base.pressure_ratio(A, B), rel_tol=1e-9
    )


def test_c_fact_override_scales_pressure_ratio_values():
    base = CharacteristicMap.from_file(_MAP_PATH)
    scaled = CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"C_fact": 1.1})
    A, B = base.mid_speed(), 1.0
    assert math.isclose(
        scaled.pressure_ratio(A, B), 1.1 * base.pressure_ratio(A, B), rel_tol=1e-9
    )


def test_e_fact_override_scales_efficiency_values():
    base = CharacteristicMap.from_file(_MAP_PATH)
    scaled = CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"E_fact": 0.9})
    A, B = base.mid_speed(), 1.0
    assert math.isclose(
        scaled.efficiency(A, B), 0.9 * base.efficiency(A, B), rel_tol=1e-9
    )


def test_a_fact_override_scales_speed_axis():
    base = CharacteristicMap.from_file(_MAP_PATH)
    scaled = CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"A_fact": 2.0})
    assert math.isclose(scaled.mid_speed(), 2.0 * base.mid_speed(), rel_tol=1e-9)


def test_unknown_factor_key_raises():
    with pytest.raises(ValueError, match="Unknown map conversion factor"):
        CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"D_fact": 1.5})


def test_partial_override_leaves_other_factors_at_file_default():
    base = CharacteristicMap.from_file(_MAP_PATH)
    scaled = CharacteristicMap.from_file(_MAP_PATH, factor_overrides={"C_fact": 1.2})
    A, B = base.mid_speed(), 1.0
    # E_fact untouched -> efficiency should be identical to the base map.
    assert math.isclose(scaled.efficiency(A, B), base.efficiency(A, B), rel_tol=1e-9)


def _map_text(pr_block: str) -> str:
    """A minimal single-angle map around a caller-supplied PR data block."""
    return f"""\
       23
************************ Compressor characteristics *********************
Compressor
************************* Conversion factors *************************
       1       (A_fact) To convert to  [rev/s]/sqrt[K]
       1     (B_fact) To convert to  [kg/s]*sqrt[K]/Bar
       1    (C_fact) To convert to  pressure ratio
       1   (E_fact) To convert to  fraction
********************** Pressure Ratio vs Non-Dimensional Mass Flow
           0     Angle
{pr_block}-999999999 *** End of angle 0
-999999999999 *** Terminator code
********************** Efficiency Or Corrected Work vs Non-Dimensional Mass Flow
           0     Angle
          50
           1            2            3   -999999999
        0.70         0.80         0.75   -999999999
          70
           1            2            3   -999999999
        0.75         0.85         0.80   -999999999
-999999999 *** End of angle 0
-999999999999 *** Terminator code
"""


# Two speed lines whose choke flows differ 2x (3.0 vs 6.0), each ending in the
# near-vertical choke leg these maps really use: the last two tabulated mass
# flows are 0.001 apart while PR collapses. This is the shape that made the
# T100 compressor map produce a spurious near-discontinuity partway between
# speed lines -- see _interpolate's docstring.
_MISMATCHED_CHOKE = _map_text(
    """          50
           1            2            3        3.001   -999999999
         2.0          2.2          2.4          1.0   -999999999
          70
           2            4            6        6.002   -999999999
         3.0          3.2          3.4          1.0   -999999999
"""
)

# Same curves, but both lines choking at the same mass flow, so aligning on
# fraction-of-choke is a no-op and the blend is the plain raw-B one.
_MATCHED_CHOKE = _map_text(
    """          50
           1            2            3        3.001   -999999999
         2.0          2.2          2.4          1.0   -999999999
          70
           1            2            3        3.001   -999999999
         3.0          3.2          3.4          1.0   -999999999
"""
)


def test_speed_lines_are_blended_at_matched_fraction_of_choke():
    """At an intermediate speed, each line must be read at the same fraction of
    its OWN choke flow -- not at the same raw corrected mass flow."""
    cmap = CharacteristicMap.from_text(_MISMATCHED_CHOKE)
    # Midway in speed, choke flow blends to (3.0 + 6.0)/2 = 4.5. Querying at
    # B = 3.0 is then 2/3 of the way along both curves: the lower line's 2.0
    # (PR 2.2) and the upper line's 4.0 (PR 3.2), blended 50/50 -> 2.7.
    assert math.isclose(cmap.pressure_ratio(60.0, 3.0), 2.7, rel_tol=1e-9)
    # And the blended choke flow itself lands on both lines' collapsed ends.
    assert math.isclose(cmap.pressure_ratio(60.0, 4.5), 2.4 / 2 + 3.4 / 2, rel_tol=1e-9)


def test_lower_line_choke_does_not_cut_into_intermediate_speeds():
    """Regression: the lower line's choke collapse must not appear as a cliff at
    intermediate speeds, at a mass flow where the higher line is nowhere near
    choke.

    Before choke alignment this was the T100 model's central defect -- the
    interpolated characteristic dropped near-vertically at the LOWER line's
    choke flow for every speed in between, which froze mass flow, put two
    spurious folds in the steady speed-vs-power curve, and made a
    load-following transient across that band diverge.
    """
    cmap = CharacteristicMap.from_text(_MISMATCHED_CHOKE)
    step = 1.0e-3
    # Straddle the lower line's choke flow (3.0)...
    across_choke = abs(
        cmap.pressure_ratio(60.0, 3.0 + step) - cmap.pressure_ratio(60.0, 3.0)
    )
    # ...and compare with an equal-sized step on the gentle part of the curve.
    reference = abs(
        cmap.pressure_ratio(60.0, 2.0 + step) - cmap.pressure_ratio(60.0, 2.0)
    )
    assert across_choke < 5.0 * reference, (
        f"pressure ratio moved {across_choke:.3e} across the lower line's choke "
        f"flow versus {reference:.3e} over the same step elsewhere -- the lower "
        f"line's choke leg is cutting into an intermediate speed"
    )


def test_matched_choke_flows_reduce_to_plain_blend():
    """When both lines choke at the same flow, normalising changes nothing."""
    cmap = CharacteristicMap.from_text(_MATCHED_CHOKE)
    for B in (1.0, 1.5, 2.0, 2.5, 3.0):
        expected = 0.5 * (
            cmap.pressure_ratio(50.0, B) + cmap.pressure_ratio(70.0, B)
        )
        assert math.isclose(cmap.pressure_ratio(60.0, B), expected, rel_tol=1e-9)


def test_speeds_outside_the_map_range_use_the_nearest_line_directly():
    cmap = CharacteristicMap.from_text(_MISMATCHED_CHOKE)
    # Below the lowest / above the highest speed there is no second line to
    # align against, so the nearest line is read at the raw B.
    assert math.isclose(cmap.pressure_ratio(10.0, 2.0), 2.2, rel_tol=1e-9)
    assert math.isclose(cmap.pressure_ratio(99.0, 4.0), 3.2, rel_tol=1e-9)
