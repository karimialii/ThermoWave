from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class ShaftLoad(BaseComponent):
    """Constant-power mechanical load on a Shaft — the electrical demand a
    generator (plus its power electronics) actually places on the rotor.

    Passive readers like SimpleGenerator/Generator only *report* electrical
    power derived from what the turbine happens to produce; nothing pushes
    back on the shaft, so a dynamic Shaft's equilibrium is always the
    free-spool point (turbine power == compressor power) regardless of any
    generator hanging off it. ShaftLoad is the other way around: it *is* a
    torque on the shaft. List it in a Shaft's components (sign -1.0 for a
    draw) and its demanded power enters the shaft power balance directly —
    the steady equilibrium speed becomes wherever turbine power covers
    compressor power *plus* this demand, and a transient integrates any
    imbalance through the rotor inertia. This is the grid-dispatch model of
    a generator: power electronics command a power draw, the shaft speed is
    whatever the machine settles at (exactly how a variable-speed
    microturbine like the T100 operates).

    power is the ELECTRICAL demand [W]; the shaft-side draw reported into
    the balance is power / efficiency (generator + power-electronics
    losses land on the shaft as extra torque). power is a plain mutable
    attribute, so a Schedule can step it through a dispatch profile during
    solve_transient(), and a plain assignment re-dispatches between steady
    solves.

    Contributes no ports and no residuals of its own — it has no speed
    unknown (the Shaft skips the speed-tie residual for torque-only members
    like this; see Shaft's own docstring) and no flow connection. Sign
    convention note: report_metrics()["power [W]"] is positive for a draw —
    the Shaft's own signs list (-1.0) is what makes it subtract, the same
    way a Compressor's positive power is subtracted there.
    """

    def __init__(self, name: str, power: float, efficiency: float = 1.0):
        if power < 0.0:
            raise ValueError(f"ShaftLoad {name!r}: power must be >= 0, got {power}")
        if not (0.0 < efficiency <= 1.0):
            raise ValueError(
                f"ShaftLoad {name!r}: efficiency must be in (0, 1], got {efficiency}"
            )
        self.name = name
        self.power = power
        self.efficiency = efficiency

    def ports(self) -> dict[str, str]:
        return {}

    def report_category(self) -> str:
        return "shaft_load"

    def residuals(self, state: "NetworkState") -> list[float]:
        return []

    def report_metrics(self, state: "NetworkState") -> dict[str, float]:
        return {
            "power [W]": self.power / self.efficiency,
            "power_elec [W]": self.power,
            "eta [-]": self.efficiency,
        }
