"""Physical/unit constants and cross-component tuning defaults shared by more
than one module — kept here so each value has exactly one definition instead
of being copy-pasted (and able to silently drift) across the files that use
it.
"""

STANDARD_ATMOSPHERE_PA = 101325.0  # Pa, standard atmosphere

PA_PER_BAR = 1.0e5  # Pa per bar, for characteristic-map corrected-flow units

N_GUESS_T_FALLBACK = 300.0  # K, fallback inlet temperature for a map-based
# Compressor/Turbine's free shaft-speed initial guess when no warm-started
# inlet context is available yet

MDOT_FUEL_GUESS_FRACTION = 0.02  # ~2% fuel-air ratio, a plausible gas-turbine
# fallback initial guess for a combustor's free mdot_fuel

STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m^2*K^4), for Radiation's Q ~ T^4 law
