"""More fluid models: IdealGasMixtureFluid, CanteraFluid, and seeing what a
Combustor's equilibrium chemistry actually produces.

Three things this shows:

1. IdealGasMixtureFluid — a calorically-perfect ideal gas blended from named
   species (mass fractions), instead of a single hand-picked (R, cp) pair.
   No extra dependency. Good for e.g. humid air or a fixed flue-gas blend
   when constant-cp accuracy is good enough.

2. CanteraFluid — the same idea, but backed by Cantera's real (temperature-
   dependent) thermo instead of constant cp, for when that accuracy matters.
   Requires the 'cantera' extra.

3. Combustor.product_composition(state) / report_metrics()'s "X_<species>
   [-]" entries — the equilibrium combustion product mole fractions Combustor
   already computes internally to find T_out, exposed directly instead of
   discarded. This is real chemistry (CO2, H2O, excess O2, ppm-level NO, ...),
   just not propagated onward through the network's own single shared fluid
   model (see Combustor's docstring for why).

Run: .venv/bin/python examples/22_fluid_mixtures_and_combustion_products.py
"""

from thermowave.components.combustor import Combustor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.core.network import Network
from thermowave.fluids.cantera_fluid import CanteraFluid
from thermowave.fluids.ideal_gas_mixture import IdealGasMixtureFluid

# --- 1. IdealGasMixtureFluid: humid air (mass fractions) -------------------
humid_air = IdealGasMixtureFluid(
    name="humid_air", composition={"N2": 0.750, "O2": 0.228, "H2O": 0.022},
)
print(f"Humid air: R = {humid_air.R:.2f} J/(kg*K), cp = {humid_air._cp:.2f} J/(kg*K)")
print(f"  h(300 K) = {humid_air.enthalpy_pt(101325.0, 300.0):.1f} J/kg")

dry_air = IdealGasMixtureFluid(name="dry_air", composition={"N2": 0.767, "O2": 0.233})
print(f"Dry air:   R = {dry_air.R:.2f} J/(kg*K), cp = {dry_air._cp:.2f} J/(kg*K)")

# --- 2. CanteraFluid: real-gas thermo for the same nominal air mixture ------
cantera_air = CanteraFluid(name="cantera_air", composition="O2:0.21, N2:0.79")
h_300 = cantera_air.enthalpy_pt(101325.0, 300.0)
h_1000 = cantera_air.enthalpy_pt(101325.0, 1000.0)
cp_300 = cantera_air.cp(101325.0, 300.0)
cp_1000 = cantera_air.cp(101325.0, 1000.0)
print(f"\nCanteraFluid air: cp(300 K) = {cp_300:.1f} J/(kg*K), cp(1000 K) = {cp_1000:.1f} J/(kg*K)")
print(
    f"  (constant-cp models can't show this rise — real air's cp increases "
    f"noticeably with T; ideal-gas mixture cp above was fixed at ~{dry_air._cp:.0f})"
)

# --- 3. Combustor product composition --------------------------------------
src = Source(name="src", P=300000.0, T=500.0, mdot=1.0)
comb = Combustor(name="cc1", PR=0.96, mdot_fuel=0.02, fuel="CH4")
snk = Sink(name="snk")

network = Network(fluid=cantera_air)
for component in (src, comb, snk):
    network.add_component(component)
network.connect(src, "out", comb, "in")
network.connect(comb, "out", snk, "in")

result = network.solve(tol=1e-6, max_iter=200)
result.print_report()

from thermowave.core.network import NetworkState

state = NetworkState(
    fluid=result.fluid, node_P=result.node_P, node_h=result.node_h,
    node_mdot=result.node_mdot, params=result.params,
)
products = comb.product_composition(state)
print("\nFull equilibrium product composition (mole fractions):")
for species, fraction in sorted(products.items(), key=lambda kv: -kv[1]):
    print(f"  {species:>6}: {fraction:.6f}")
