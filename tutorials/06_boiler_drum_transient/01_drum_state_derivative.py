"""Tutorial 6, step 1: reading a Drum's own differential state directly.

A Drum's (P, h) are differential state (see BaseComponent.
differential_parameters()/state_derivative(), the same mechanism tutorial
3's dynamic Shaft uses for rotor speed): dP/dt and dh/dt come from mass and
energy conservation on the drum's saturated liquid/vapor inventory, given
whatever's flowing in and out of its four ports (feed_in, riser_in,
steam_out, water_out).

This step calls state_derivative() directly — no Network, no solve, no
time-stepping — just to see what it actually computes at two hand-picked
states: one where inflow and outflow are in exact mass/energy balance
(both derivatives should come out at essentially zero), and one where
outflow exceeds inflow (the drum should be draining: pressure and stored
enthalpy both falling). Understanding this one function is really all you
need before letting solve_transient() call it automatically hundreds of
times, which step 2 does.

What you'll learn:
  - state_derivative(state) takes anything shaped like a NetworkState (here,
    a tiny hand-written stand-in) and returns {"P": dP/dt, "h": dh/dt}.
  - A drum has no steady-state level — feed exactly matching draw gives
    dP/dt == dh/dt == 0 (equilibrium), not "level holds at some computed
    value"; there's no such value to compute.
  - Draining (draw > feed) shows up as both derivatives going negative:
    less stored mass (falling P) and a shift toward more liquid in the
    remaining inventory (falling average h).

Requires the 'coolprop' extra: pip install thermowave[coolprop]

Run: .venv/bin/python tutorials/06_boiler_drum_transient/01_drum_state_derivative.py
"""

from thermowave.components.drum import Drum
from thermowave.fluids.real_fluid import CoolPropFluid

water = CoolPropFluid(name="Water")

P0 = 1.0e6  # 10 bar
drum = Drum(name="drum", V=1.0, P0=P0, fluid=water, level0=0.5, has_riser=True)

h_f = water.saturated_liquid_enthalpy(P0)
h_g = water.saturated_vapor_enthalpy(P0)
h_fg = h_g - h_f
H_FEED = water.enthalpy_pt(P0, 440.0)  # subcooled feedwater, fixed boundary state
CIRC = 8.0  # kg/s downcomer/riser circulation


class _DrumState:
    """The minimum a caller needs to hand state_derivative(): node (P, h)
    and mdot for each of the drum's own ports, plus its current (P, h)."""

    def __init__(self, mdot_feed, mdot_riser, mdot_steam, mdot_water, h_riser):
        self.P = P0
        self.h = drum.h0
        self._mdot = {
            "drum.feed_in": mdot_feed, "drum.riser_in": mdot_riser,
            "drum.steam_out": mdot_steam, "drum.water_out": mdot_water,
        }
        self._h_riser = h_riser

    def node(self, name):
        return {
            "drum.feed_in": (self.P, H_FEED),
            "drum.riser_in": (self.P, self._h_riser),
            "drum.steam_out": (self.P, water.saturated_vapor_enthalpy(self.P)),
            "drum.water_out": (self.P, water.saturated_liquid_enthalpy(self.P)),
        }[name]

    def fluid_at(self, name):
        return water

    def mdot(self, name):
        return self._mdot[name]

    def param(self, name):
        return {"drum.P": self.P, "drum.h": self.h}[name]


# Riser quality chosen so mass AND energy exactly balance at 1.0 kg/s
# feed/steam and 8.0 kg/s circulation.
mdot_feed = 1.0
x_riser_balanced = (h_g - H_FEED) / (CIRC * h_fg)
h_riser_balanced = h_f + x_riser_balanced * h_fg

balanced = _DrumState(
    mdot_feed=mdot_feed, mdot_riser=CIRC, mdot_steam=mdot_feed, mdot_water=CIRC,
    h_riser=h_riser_balanced,
)
deriv_balanced = drum.state_derivative(balanced)
print("Balanced (feed == steam draw):")
print(f"  dP/dt = {deriv_balanced['P']:.4f} Pa/s, dh/dt = {deriv_balanced['h']:.4f} J/(kg*s)")

draining = _DrumState(
    mdot_feed=mdot_feed, mdot_riser=CIRC, mdot_steam=1.3 * mdot_feed, mdot_water=CIRC,
    h_riser=h_riser_balanced,
)
deriv_draining = drum.state_derivative(draining)
print("\nDraining (steam draw 30% above feed):")
print(f"  dP/dt = {deriv_draining['P']:.2f} Pa/s, dh/dt = {deriv_draining['h']:.2f} J/(kg*s)")
print("  (both negative: pressure falling, stored enthalpy falling -> level falling)")

print(
    "\nStep 2 integrates dP/dt and dh/dt forward through exactly this kind "
    "of imbalance over time, to see the drum's actual level and pressure "
    "trajectory during a steam-demand step."
)
