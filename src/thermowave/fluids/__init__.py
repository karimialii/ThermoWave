"""Fluid property models, re-exported for a one-line import.

    from thermowave.fluids import IdealGasFluid, CoolPropFluid

CanteraFluid and CoolPropFluid import their optional extra inside __init__,
not at module level, so importing this works without either installed --
the ImportError still arrives when you construct one.
"""

from thermowave.fluids.base_fluid import BaseFluid, ConstantCpFluid
from thermowave.fluids.cantera_fluid import CanteraFluid
from thermowave.fluids.humid_air import HumidAirFluid
from thermowave.fluids.ideal_gas import IdealGasFluid
from thermowave.fluids.ideal_gas_mixture import IdealGasMixtureFluid
from thermowave.fluids.real_fluid import CoolPropFluid

__all__ = [
    "BaseFluid",
    "CanteraFluid",
    "ConstantCpFluid",
    "CoolPropFluid",
    "HumidAirFluid",
    "IdealGasFluid",
    "IdealGasMixtureFluid",
]
