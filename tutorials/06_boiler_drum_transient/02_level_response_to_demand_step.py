"""Tutorial 6, step 2: the drum's actual level trajectory over time.

Step 1 looked at state_derivative() at two frozen snapshots. This
integrates it forward through real time: starting from the balanced state,
run for 5 seconds (nothing moves — feed matches draw), then step steam
demand up 30% and watch pressure and level respond continuously.

This uses a small explicit Euler step directly (stable and step-size-
independent for this system) rather than calling Network.solve_transient()
on a full network — the same (P, h) ODE either way, just without the
overhead of a surrounding network and level controller a real closed-loop
drum would need to be well-posed as a full Network.solve_transient() case.

What you'll learn:
  - Reading a differential-state component's own state_derivative()
    directly, integrated by hand, is a legitimate way to explore its
    dynamics in isolation — you don't need a full network to understand
    one component's time response.
  - A steam-demand step shows up as a smooth pressure and level decline,
    not an instantaneous jump — the drum's stored mass provides a real
    buffer, which is exactly why real drums exist (to give a control system
    time to react to load changes).

Requires the 'coolprop' extra: pip install thermowave[coolprop]

Run: .venv/bin/python tutorials/06_boiler_drum_transient/02_level_response_to_demand_step.py
"""

from thermowave.components.drum import Drum
from thermowave.fluids.real_fluid import CoolPropFluid

water = CoolPropFluid(name="Water")

P0 = 1.0e6  # 10 bar drum pressure
drum = Drum(name="drum", V=1.0, P0=P0, fluid=water, level0=0.5, has_riser=True)

# Fixed boundary-stream enthalpies (a real feedwater/riser stream's state
# doesn't follow the drum pressure).
h_f = water.saturated_liquid_enthalpy(P0)
h_g = water.saturated_vapor_enthalpy(P0)
h_fg = h_g - h_f
H_FEED = water.enthalpy_pt(P0, 440.0)  # subcooled feedwater
CIRC = 8.0  # kg/s downcomer/riser circulation
MDOT_FEED = 1.0  # kg/s feedwater (= steam production at balance)
# Riser quality chosen so the drum starts in exact energy balance (flat
# baseline until the demand step) — the same balanced state step 1 checked.
x_riser = (h_g - H_FEED) / (CIRC * h_fg)
H_RISER = h_f + x_riser * h_fg


class _DrumState:
    """Minimal NetworkState stand-in: prescribes the drum's four port flows
    and reads its differential (P, h)."""

    def __init__(self, P, h, mdot_steam):
        self.P = P
        self.h = h
        self._mdot = {
            "drum.feed_in": MDOT_FEED,
            "drum.riser_in": CIRC,
            "drum.steam_out": mdot_steam,
            "drum.water_out": CIRC,
        }

    def node(self, name):
        return {
            "drum.feed_in": (self.P, H_FEED),
            "drum.riser_in": (self.P, H_RISER),
            "drum.steam_out": (self.P, water.saturated_vapor_enthalpy(self.P)),
            "drum.water_out": (self.P, water.saturated_liquid_enthalpy(self.P)),
        }[name]

    def fluid_at(self, name):
        return water

    def mdot(self, name):
        return self._mdot[name]

    def param(self, name):
        return {"drum.P": self.P, "drum.h": self.h}[name]


def level(P, h):
    hf = water.saturated_liquid_enthalpy(P)
    hg = water.saturated_vapor_enthalpy(P)
    x = min(max((h - hf) / (hg - hf), 0.0), 1.0)
    v_avg = 1.0 / water.density_ph(P, h)
    v_f = 1.0 / water.density_ph(P, hf)
    return (1.0 - x) * v_f / v_avg


duration, dt, step_time = 40.0, 0.01, 5.0
n_steps = int(duration / dt)

times, pressures, levels = [], [], []
P, h = P0, drum.h0
for i in range(n_steps + 1):
    t = i * dt
    mdot_steam = MDOT_FEED if t < step_time else 1.3 * MDOT_FEED  # +30% demand step
    times.append(t)
    pressures.append(P)
    levels.append(level(P, h))
    deriv = drum.state_derivative(_DrumState(P, h, mdot_steam))
    P += dt * deriv["P"]
    h += dt * deriv["h"]

print(f"{'t [s]':>7}  {'P [bar]':>9}  {'level [-]':>10}")
for i in range(0, len(times), n_steps // 10):
    print(f"{times[i]:7.1f}  {pressures[i] / 1e5:9.3f}  {levels[i]:10.4f}")

print(
    f"\nSteam-demand step at t={step_time:.0f}s: "
    f"level {levels[0]:.4f} -> {levels[-1]:.4f}, "
    f"P {pressures[0] / 1e5:.2f} -> {pressures[-1] / 1e5:.2f} bar (drum draining)."
)

# Plotting requires the 'plot' extra: pip install thermowave[plot]. Saved
# headlessly (not shown) so the example stays non-interactive.
try:
    import matplotlib.pyplot as plt
except ImportError:
    print("\n(install thermowave[plot] to render drum_transient.png)")
else:
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 6))
    ax1.plot(times, levels, color="tab:blue")
    ax1.axvline(step_time, color="grey", ls="--", alpha=0.6)
    ax1.set_ylabel("drum level [-]")
    ax1.set_title("Steam drum response to a +30% steam-demand step")
    ax1.grid(True, alpha=0.3)
    ax2.plot(times, [p / 1e5 for p in pressures], color="tab:red")
    ax2.axvline(step_time, color="grey", ls="--", alpha=0.6)
    ax2.set_xlabel("t [s]")
    ax2.set_ylabel("drum pressure [bar]")
    ax2.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig("drum_transient.png", bbox_inches="tight")
    print("\nSaved drum_transient.png")

print(
    "\nThis completes the tutorial series. From here: combine this drum "
    "with tutorial 5's Rankine cycle (the drum would sit upstream of that "
    "cycle's boiler), or revisit tutorial 4's PID loops to design a level "
    "controller that holds this drum's level steady through the same "
    "demand step instead of letting it drain."
)
