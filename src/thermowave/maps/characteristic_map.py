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

    Variable-geometry maps are supported: a section may have several "Angle"
    blocks (one per vane/nozzle angle), each a complete set of iso-speed
    curves like the single-angle case. `pressure_ratio()`/`efficiency()`
    take an optional `angle=` (defaulting to 0.0, transparent for a
    single-angle map regardless of what angle its one block is actually
    tagged with) and interpolate across angle the same way they already
    interpolate across speed: bracket the two nearest angle blocks, run the
    existing speed/choke-fraction interpolation independently against each,
    then blend the two results linearly by angle fraction. Outside the map's
    angle range (or for a single-angle map) the nearest/only angle block is
    used directly, unchanged.

    The file itself specifies four conversion factors (A_fact for speed,
    B_fact for mass flow, C_fact for pressure ratio, E_fact for efficiency)
    that scale the raw tabulated values into the units above — pass
    factor_overrides to from_file()/from_text() to replace any of them
    (e.g. to calibrate a map against test data) without editing the map
    file itself; factors not named in the override dict keep the file's own
    value, and passing nothing at all reproduces the file exactly as before.
    """

    def __init__(
        self,
        pr_blocks: dict[float, list[_SpeedLine]],
        eff_blocks: dict[float, list[_SpeedLine]],
    ):
        self._pr_blocks = {
            angle: sorted(lines, key=lambda line: line.speed)
            for angle, lines in pr_blocks.items()
        }
        self._eff_blocks = {
            angle: sorted(lines, key=lambda line: line.speed)
            for angle, lines in eff_blocks.items()
        }

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

        pr_blocks = _parse_section(
            text,
            _PR_SECTION_HEADER,
            a_fact=factors["A_fact"],
            b_fact=factors["B_fact"],
            y_fact=factors["C_fact"],
        )
        eff_blocks = _parse_section(
            text,
            _EFF_SECTION_HEADER,
            a_fact=factors["A_fact"],
            b_fact=factors["B_fact"],
            y_fact=factors["E_fact"],
        )
        return cls(pr_blocks, eff_blocks)

    def pressure_ratio(self, A: float, B: float, angle: float = 0.0) -> float:
        """Interpolated pressure ratio at corrected speed A, corrected mass
        flow B, and (for a variable-geometry map) vane/nozzle angle. angle is
        ignored when the map has only one angle block (the default 0.0 is
        transparent regardless of what angle that one block is tagged with)."""
        return _interpolate_angle(self._pr_blocks, A, B, angle)

    def efficiency(self, A: float, B: float, angle: float = 0.0) -> float:
        """Interpolated efficiency (fraction) at corrected speed A, mass flow
        B, and (for a variable-geometry map) vane/nozzle angle. Same angle
        semantics as pressure_ratio()."""
        return _interpolate_angle(self._eff_blocks, A, B, angle)

    def mid_speed(self, angle: float = 0.0) -> float:
        """Midpoint of the corrected-speed range A of the block nearest
        `angle` — a reasonable Newton initial guess when speed is solved for
        instead of given directly."""
        lines = _nearest_block(self._pr_blocks, angle)
        speeds = [line.speed for line in lines]
        return (min(speeds) + max(speeds)) / 2.0

    def speed_lines(
        self, kind: str, angle: float = 0.0
    ) -> list[tuple[float, list[float], list[float]]]:
        """Raw iso-speed curves for plotting: one (speed, mass_flow, value)
        tuple per corrected speed, sorted by speed, from the angle block
        nearest `angle`. kind is "pressure_ratio" or "efficiency", selecting
        which section's curves to return."""
        if kind == "pressure_ratio":
            blocks = self._pr_blocks
        elif kind == "efficiency":
            blocks = self._eff_blocks
        else:
            raise ValueError(
                f"Unknown kind {kind!r}; expected 'pressure_ratio' or 'efficiency'"
            )
        lines = _nearest_block(blocks, angle)
        return [(line.speed, line.mass_flow, line.value) for line in lines]

    def angles(self) -> list[float]:
        """Sorted list of vane/nozzle angles this map has blocks for — a
        single-angle map returns a single-element list."""
        return sorted(self._pr_blocks)


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
) -> dict[float, list[_SpeedLine]]:
    lines = text.splitlines()
    try:
        start = next(i for i, line in enumerate(lines) if section_header in line)
    except StopIteration:
        raise ValueError(f"Map is missing section {section_header!r}") from None

    i = start + 1
    while "Angle" not in lines[i]:
        i += 1

    blocks: dict[float, list[_SpeedLine]] = {}
    while True:
        angle = float(lines[i].strip().split()[0])
        i += 1  # past the "<angle>  Angle" line

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
        blocks[angle] = speed_lines

        i += 1  # past "-999999999 *** End of angle <angle>"
        if "Terminator" in lines[i]:
            break
        while "Angle" not in lines[i]:  # blank/comment lines between blocks, if any
            i += 1

    return blocks


def _parse_row(line: str, fact: float) -> list[float]:
    values: list[float] = []
    for token in line.split():
        value = float(token)
        if value <= _SENTINEL:
            break
        values.append(value * fact)
    return values


def _nearest_block(
    blocks: dict[float, list[_SpeedLine]], angle: float
) -> list[_SpeedLine]:
    """The block whose own angle is closest to `angle` (clamped at the ends)
    -- used where a single representative curve set is wanted (mid_speed(),
    speed_lines()), as opposed to _interpolate_angle()'s blend across two
    bracketing blocks."""
    angles = sorted(blocks)
    if angle <= angles[0]:
        return blocks[angles[0]]
    if angle >= angles[-1]:
        return blocks[angles[-1]]
    closest = min(angles, key=lambda a: abs(a - angle))
    return blocks[closest]


def _interpolate_angle(
    blocks: dict[float, list[_SpeedLine]], A: float, B: float, angle: float
) -> float:
    """Interpolate across vane/nozzle angle blocks at (A, B), same bracket-
    and-blend pattern as _interpolate()'s speed bracketing: run the existing
    speed/choke-fraction interpolation independently against each of the two
    bracketing angle blocks, then blend the two results linearly by angle
    fraction. A single-angle map (or angle outside the map's angle range)
    uses the nearest/only block directly, unchanged."""
    angles = sorted(blocks)
    if len(angles) == 1 or angle <= angles[0]:
        return _interpolate(blocks[angles[0]], A, B)
    if angle >= angles[-1]:
        return _interpolate(blocks[angles[-1]], A, B)

    for lower, upper in zip(angles, angles[1:]):
        if lower <= angle <= upper:
            weight = (angle - lower) / (upper - lower)
            y_lower = _interpolate(blocks[lower], A, B)
            y_upper = _interpolate(blocks[upper], A, B)
            return y_lower + weight * (y_upper - y_lower)

    raise AssertionError("unreachable: angle within the map's overall angle range")


def _interpolate(speed_lines: list[_SpeedLine], A: float, B: float) -> float:
    """Interpolate between iso-speed curves at corrected speed A, mass flow B.

    Between two bracketing speed lines, each line is evaluated at the same
    FRACTION of its own maximum (choke) mass flow rather than at the same raw
    corrected mass flow B, then the two results are blended by speed.

    Why this matters, concretely. Every speed line in these maps ends in a
    near-vertical choke leg: its last few tabulated mass flows are spaced ~1e-3
    apart while the pressure ratio collapses. Each line chokes at its own mass
    flow, and those differ a lot between lines — in the T100 compressor map the
    pair bracketing 60,000 rev/min choke at B = 9.978 and B = 12.308, 23% apart.
    Evaluating both at the same raw B therefore puts the lower line hundreds of
    times further along its curve than the upper one: at B = 9.98 the lower line
    is exactly on its choke knee (PR falling 2.744 -> 2.152 over dB = 0.003,
    roughly 680x steeper than its normal leg) while the upper line is still on
    its gentle leg. The blend inherits that collapse, so the interpolated
    characteristic has a spurious near-discontinuity at the LOWER line's choke
    flow, at every speed in between — a mass flow where the real machine at the
    higher speed is nowhere near choke.

    In the T100 model that artifact froze air mass flow at 0.6044 kg/s across
    59,500-61,000 rev/min while speed rose 1,500 rev/min, produced two spurious
    folds in the steady speed-vs-power curve (two valid steady states at 65 kW,
    3,400 rev/min apart), and made a load-following transient across that band
    diverge. Aligning on fraction-of-choke removes all of it: the folds go away,
    mass flow becomes smooth and monotone, and the modelled speed-vs-power slope
    goes from 96 rev/min per kW on that spurious upper branch to 249, against
    ~212 measured on the rig.

    This is the same normalisation as the beta-line (curve-parameter) map
    representation used by GasTurb/NPSS/PROOSIS, reduced to the minimum needed
    here: a shared 0..1 coordinate along each curve so corresponding points on
    neighbouring speed lines are blended with each other.

    Outside the map's speed range (and for a single-line map) the nearest line is
    used directly at the raw B, unchanged — there is no second line to align to.
    """
    if not speed_lines:
        raise ValueError("Map has no speed lines to interpolate")

    if len(speed_lines) == 1 or A <= speed_lines[0].speed:
        return _interpolate_1d(speed_lines[0], B)
    if A >= speed_lines[-1].speed:
        return _interpolate_1d(speed_lines[-1], B)

    for lower, upper in zip(speed_lines, speed_lines[1:]):
        if lower.speed <= A <= upper.speed:
            weight = (A - lower.speed) / (upper.speed - lower.speed)
            b_lower, b_upper = lower.mass_flow[-1], upper.mass_flow[-1]
            b_choke = b_lower + weight * (b_upper - b_lower)
            # A degenerate line (single point, or a zero choke flow) has no
            # scale to normalise against; fall back to the raw-B blend.
            if b_choke > 0.0 and b_lower > 0.0 and b_upper > 0.0:
                fraction = B / b_choke
                y_lower = _interpolate_1d(lower, fraction * b_lower)
                y_upper = _interpolate_1d(upper, fraction * b_upper)
            else:
                y_lower = _interpolate_1d(lower, B)
                y_upper = _interpolate_1d(upper, B)
            return y_lower + weight * (y_upper - y_lower)

    raise AssertionError("unreachable: A within the map's overall speed range")


def _interpolate_1d(line: _SpeedLine, B: float) -> float:
    b, y = line.mass_flow, line.value
    b_clamped = min(max(B, b[0]), b[-1])
    return float(np.interp(b_clamped, b, y))
