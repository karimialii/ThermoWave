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
