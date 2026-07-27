"""Turbomachinery and torque-speed maps, re-exported for a one-line import.

    from thermowave.maps import CharacteristicMap
"""

from thermowave.maps.characteristic_map import CharacteristicMap
from thermowave.maps.torque_speed_map import TorqueSpeedMap

__all__ = ["CharacteristicMap", "TorqueSpeedMap"]
