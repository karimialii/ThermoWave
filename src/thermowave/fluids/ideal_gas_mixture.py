import math

from thermowave.fluids.base_fluid import ConstantCpFluid

_R_UNIVERSAL = 8.314462618  # J/(mol*K)

# name -> (molar mass [kg/mol], molar cp [J/(mol*K)]). cp values are rough
# constants (representative of ~300-1000 K), the same "calorically perfect
# gas" simplification IdealGasFluid already makes for a single species, just
# per-species here — good enough to blend a mixture's R/cp, not a substitute
# for temperature-dependent real thermo (see CanteraFluid for that).
_DEFAULT_SPECIES: dict[str, tuple[float, float]] = {
    "N2": (0.028014, 29.1),
    "O2": (0.031998, 29.4),
    "CO2": (0.044010, 37.1),
    "H2O": (0.018015, 33.6),
    "Ar": (0.039948, 20.8),
    "CH4": (0.016043, 35.7),
    "CO": (0.028010, 29.1),
    "H2": (0.002016, 28.8),
}


class IdealGasMixtureFluid(ConstantCpFluid):
    """Calorically-perfect ideal-gas mixture: mass-fraction-weighted R and cp
    from a small built-in species table, generalizing IdealGasFluid's single
    fixed (R, cp) pair to an arbitrary named mixture (humid air, flue gas, a
    fuel-air blend, ...) with no extra dependency (no CoolProp/Cantera). For
    temperature-dependent real-gas thermo of an arbitrary mixture instead,
    see CanteraFluid.

    composition: {species_name: mass_fraction}, must sum to ~1.0. Species
    must be one of IdealGasMixtureFluid.SPECIES, or pass extra_species (same
    (molar_mass [kg/mol], molar_cp [J/(mol*K)]) form) to add/override
    entries — e.g. a custom fuel not in the built-in table.
    """

    SPECIES = dict(_DEFAULT_SPECIES)

    def __init__(
        self,
        name: str,
        composition: dict[str, float],
        extra_species: dict[str, tuple[float, float]] | None = None,
    ):
        species = {**self.SPECIES, **(extra_species or {})}
        total = sum(composition.values())
        if not math.isclose(total, 1.0, abs_tol=1e-3):
            raise ValueError(
                f"IdealGasMixtureFluid {name!r}: composition mass fractions sum to "
                f"{total:.4f}, expected 1.0"
            )
        unknown = set(composition) - set(species)
        if unknown:
            raise ValueError(
                f"IdealGasMixtureFluid {name!r}: unknown species {sorted(unknown)}; "
                f"pass extra_species=... or use one of {sorted(species)}"
            )

        R = 0.0
        cp = 0.0
        for species_name, mass_fraction in composition.items():
            molar_mass, molar_cp = species[species_name]
            R += mass_fraction * (_R_UNIVERSAL / molar_mass)
            cp += mass_fraction * (molar_cp / molar_mass)

        self.name = name
        self.composition = dict(composition)
        self.R = R
        self._cp = cp
