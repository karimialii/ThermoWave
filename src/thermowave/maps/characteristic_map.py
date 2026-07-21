from __future__ import annotations

from pathlib import Path

import numpy as np

_SENTINEL = -999999999.0
_PR_SECTION_HEADER = "Pressure Ratio vs Non-Dimensional Mass Flow"
_EFF_SECTION_HEADER = "Efficiency Or Corrected Work vs Non-Dimensional Mass Flow"
_FACTOR_KEYS = ("A_fact", "B_fact", "C_fact", "E_fact")


class _SpeedLine:
    """One iso-speed curve: dependent variable (PR or efficiency) vs
    non-dimensional mass flow, at a single corrected speed."""

    def __init__(self, speed: float, mass_flow: list[float], value: list[float]):
        self.speed = speed
        self.mass_flow = mass_flow
        self.value = value


class CharacteristicMap:
    """A Flownex-style turbomachine characteristic map (.cop/.tur format):
    iso-speed curves of pressure ratio and efficiency vs non-dimensional mass
    flow. Same file format for compressors and turbines.

    Non-dimensional groups (matching the file's own header comments):
        A = N / sqrt(T_in)             corrected speed        [rev/s / sqrt(K)]
        B = mdot * sqrt(T_in) / P_in   corrected mass flow     [kg/s * sqrt(K) / Bar]
        C = P_out / P_in               pressure ratio          [-]

    (The file's inline comment literally reads "C = P03/P04", i.e. inlet
    over outlet, but every actual data value in these maps is > 1 for a
    compressor and its own C_fact line says "to convert to pressure ratio" —
    so C is read here as the standard P_out/P_in ratio, trusting the data
    and the conversion-factor line over what looks like a header typo.)

    Only supports single-angle maps (both known example files have exactly
    one "Angle 0" block per section); a variable-geometry map with multiple
    angle blocks isn't handled.

    The file itself specifies four conversion factors (A_fact for speed,
    B_fact for mass flow, C_fact for pressure ratio, E_fact for efficiency)
    that scale the raw tabulated values into the units above — pass
    factor_overrides to from_file()/from_text() to replace any of them
    (e.g. to calibrate a map against test data) without editing the map
    file itself; factors not named in the override dict keep the file's own
    value, and passing nothing at all reproduces the file exactly as before.
    """

    def __init__(self, pr_lines: list[_SpeedLine], eff_lines: list[_SpeedLine]):
        self._pr_lines = sorted(pr_lines, key=lambda line: line.speed)
        self._eff_lines = sorted(eff_lines, key=lambda line: line.speed)

    @classmethod
    def from_file(
        cls, path: str, factor_overrides: dict[str, float] | None = None
    ) -> "CharacteristicMap":
        return cls.from_text(Path(path).read_text(), factor_overrides=factor_overrides)

    @classmethod
    def from_text(
        cls, text: str, factor_overrides: dict[str, float] | None = None
    ) -> "CharacteristicMap":
        factors = _parse_conversion_factors(text)
        if factor_overrides:
            unknown = set(factor_overrides) - set(_FACTOR_KEYS)
            if unknown:
                raise ValueError(
                    f"Unknown map conversion factor(s) {sorted(unknown)}; "
                    f"expected a subset of {list(_FACTOR_KEYS)}"
                )
            factors.update(factor_overrides)

        pr_lines = _parse_section(
            text,
            _PR_SECTION_HEADER,
            a_fact=factors["A_fact"],
            b_fact=factors["B_fact"],
            y_fact=factors["C_fact"],
        )
        eff_lines = _parse_section(
            text,
            _EFF_SECTION_HEADER,
            a_fact=factors["A_fact"],
            b_fact=factors["B_fact"],
            y_fact=factors["E_fact"],
        )
        return cls(pr_lines, eff_lines)

    def pressure_ratio(self, A: float, B: float) -> float:
        """Interpolated pressure ratio at corrected speed A, corrected mass flow B."""
        return _interpolate(self._pr_lines, A, B)

    def efficiency(self, A: float, B: float) -> float:
        """Interpolated efficiency (fraction) at corrected speed A, mass flow B."""
        return _interpolate(self._eff_lines, A, B)

    def mid_speed(self) -> float:
        """Midpoint of the map's corrected-speed range A — a reasonable Newton
        initial guess when speed is solved for instead of given directly."""
        speeds = [line.speed for line in self._pr_lines]
        return (min(speeds) + max(speeds)) / 2.0

    def speed_lines(
        self, kind: str
    ) -> list[tuple[float, list[float], list[float]]]:
        """Raw iso-speed curves for plotting: one (speed, mass_flow, value)
        tuple per corrected speed, sorted by speed. kind is "pressure_ratio"
        or "efficiency", selecting which section's curves to return."""
        if kind == "pressure_ratio":
            lines = self._pr_lines
        elif kind == "efficiency":
            lines = self._eff_lines
        else:
            raise ValueError(
                f"Unknown kind {kind!r}; expected 'pressure_ratio' or 'efficiency'"
            )
        return [(line.speed, line.mass_flow, line.value) for line in lines]


def _parse_conversion_factors(text: str) -> dict[str, float]:
    factors = {"A_fact": 1.0, "B_fact": 1.0, "C_fact": 1.0, "E_fact": 1.0}
    for line in text.splitlines():
        stripped = line.strip()
        for key in factors:
            if f"({key})" in stripped:
                factors[key] = float(stripped.split()[0])
    return factors


def _parse_section(
    text: str, section_header: str, a_fact: float, b_fact: float, y_fact: float
) -> list[_SpeedLine]:
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if section_header in line)
    except StopIteration:
        raise ValueError(f"Map is missing section {section_header!r}") from None

    i = start + 1
    while "Angle" not in lines[i]:
        i += 1
    i += 1  # past the "0  Angle" line

    speed_lines: list[_SpeedLine] = []
    while True:
        stripped = lines[i].strip()
        if stripped.startswith(str(int(_SENTINEL))):
            break
        speed = float(stripped) * a_fact
        mass_flow = _parse_row(lines[i + 1], b_fact)
        value = _parse_row(lines[i + 2], y_fact)
        speed_lines.append(_SpeedLine(speed=speed, mass_flow=mass_flow, value=value))
        i += 3

    return speed_lines


def _parse_row(line: str, fact: float) -> list[float]:
    values: list[float] = []
    for token in line.split():
        value = float(token)
        if value <= _SENTINEL:
            break
        values.append(value * fact)
    return values


def _interpolate(speed_lines: list[_SpeedLine], A: float, B: float) -> float:
    if not speed_lines:
        raise ValueError("Map has no speed lines to interpolate")

    if len(speed_lines) == 1 or A <= speed_lines[0].speed:
        return _interpolate_1d(speed_lines[0], B)
    if A >= speed_lines[-1].speed:
        return _interpolate_1d(speed_lines[-1], B)

    for lower, upper in zip(speed_lines, speed_lines[1:]):
        if lower.speed <= A <= upper.speed:
            y_lower = _interpolate_1d(lower, B)
            y_upper = _interpolate_1d(upper, B)
            weight = (A - lower.speed) / (upper.speed - lower.speed)
            return y_lower + weight * (y_upper - y_lower)

    raise AssertionError("unreachable: A within the map's overall speed range")


def _interpolate_1d(line: _SpeedLine, B: float) -> float:
    b, y = line.mass_flow, line.value
    b_clamped = min(max(B, b[0]), b[-1])
    return float(np.interp(b_clamped, b, y))
