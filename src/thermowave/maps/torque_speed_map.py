from __future__ import annotations

from pathlib import Path

import numpy as np


class TorqueSpeedMap:
    """A generator's torque-vs-speed characteristic curve: one table of
    (shaft speed [rev/min], torque [N*m]) points, linearly interpolated
    (and clamped at the ends, same as CharacteristicMap) between them.

    Much simpler than the Flownex .cop/.tur format (no iso-speed families,
    no non-dimensional groups) since a generator's rating is normally given
    directly as a single speed-torque curve, not corrected for inlet
    conditions the way a compressor/turbine map is.

    File format: one "speed torque" pair per line (whitespace-separated);
    blank lines and lines starting with "#" are ignored.
    """

    def __init__(self, speed: list[float], torque: list[float]):
        if len(speed) != len(torque) or len(speed) == 0:
            raise ValueError("speed and torque must be the same non-empty length")
        order = np.argsort(speed)
        self._speed = np.asarray(speed, dtype=float)[order]
        self._torque = np.asarray(torque, dtype=float)[order]

    @classmethod
    def from_file(cls, path: str) -> "TorqueSpeedMap":
        return cls.from_text(Path(path).read_text())

    @classmethod
    def from_text(cls, text: str) -> "TorqueSpeedMap":
        speed: list[float] = []
        torque: list[float] = []
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) != 2:
                continue
            speed.append(float(parts[0]))
            torque.append(float(parts[1]))
        return cls(speed, torque)

    def torque(self, N: float) -> float:
        """Interpolated (clamped) torque [N*m] at shaft speed N [rev/min]."""
        N_clamped = min(max(N, self._speed[0]), self._speed[-1])
        return float(np.interp(N_clamped, self._speed, self._torque))

    def mid_speed(self) -> float:
        """Midpoint of the map's speed range — a reasonable Newton initial
        guess when speed is solved for instead of given directly."""
        return (self._speed[0] + self._speed[-1]) / 2.0
