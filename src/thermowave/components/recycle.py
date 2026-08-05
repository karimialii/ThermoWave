from __future__ import annotations

from typing import TYPE_CHECKING

from thermowave.components.base_component import BaseComponent

if TYPE_CHECKING:
    from thermowave.core.network import NetworkState


class Recycle(BaseComponent):
    """Mass-flow anchor for a genuinely closed loop -- no Source/Sink needed.

    Every flow component's own mass residual (mdot_out - mdot_in = 0) already
    chains mdot all the way around a loop; what a loop with no Source in it
    is missing isn't a *relationship* between its mdot unknowns (those are
    already tied together), it's an anchor for their shared *magnitude* --
    conservation alone only pins mdot up to an arbitrary constant. Likewise,
    P and h at every node in a real closed cycle (ORC, Rankine, Brayton) are
    already fully determined by the loop's own component physics -- a
    SimpleCondenser computes its own outlet (P, h) from subcool/quality, a
    Pump raises P by a given ratio, etc. -- so nothing needs a Source's other
    job (fixing (P, h) as a hard boundary condition) either.

    Recycle supplies exactly the piece that's missing and nothing else: drop
    it anywhere in a loop (typically condenser outlet -> pump inlet) to pin
    the loop's overall mass-flow scale, with every other node's state left
    entirely to the rest of the network's own physics:

        net.connect(condenser, "out", recycle, "in")
        net.connect(recycle, "out", pump, "in")

    Ports: in, out.

    Fixes: mdot at out only -- never (P, h). This is deliberately the
    opposite half of Source's job: Source fixes (P, h) always and mdot
    optionally, Recycle fixes mdot always and never touches (P, h).

    Residuals: exactly 2 -- P_out - P_in = 0, h_out - h_in = 0, an ideal
    lossless pass-through. NOT the 3-residual pattern most flow components
    use: there's deliberately no mass-conservation residual here, since
    fixed_node_mdot() already fixes this component's own outlet mdot
    directly -- adding a mass residual on top of that would try to pin the
    same already-fixed unknown with a second equation.
    """

    def __init__(self, name: str, mdot: float):
        if mdot <= 0.0:
            raise ValueError(f"Recycle {name!r}: mdot must be > 0, got {mdot}")
        self.name = name
        self.mdot = mdot
        self._inlet_node = f"{name}.in"
        self._outlet_node = f"{name}.out"

    def ports(self) -> dict[str, str]:
        return {"in": self._inlet_node, "out": self._outlet_node}

    def residuals(self, state: "NetworkState") -> list[float]:
        P_in, h_in = state.node(self._inlet_node)
        P_out, h_out = state.node(self._outlet_node)
        return [P_out - P_in, h_out - h_in]

    def fixed_node_mdot(self) -> dict[str, float]:
        return {self._outlet_node: self.mdot}
