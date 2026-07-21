from thermowave.fluids.base_fluid import ConstantCpFluid


class IdealGasFluid(ConstantCpFluid):
    """Analytic ideal-gas fluid model with constant specific heat.

    Enthalpy is referenced to h=0 at T=0: h = cp * T.
    """

    def __init__(self, name: str, R: float, cp: float):
        self.name = name
        self.R = R
        self._cp = cp
