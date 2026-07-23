from __future__ import annotations

import sys
from typing import TextIO

_GREEN = "\033[32m"
_RED = "\033[31m"
_BOLD = "\033[1m"
_RESET = "\033[0m"


class ProgressBar:
    """A single-line, in-place terminal progress bar for anything with a
    knowable [0, 1] completion fraction (a Newton iteration budget, a
    transient run's elapsed/duration) — render() overwrites the same
    terminal line via '\\r' + ANSI erase-to-end-of-line rather than printing
    a new one each call, so a long iteration count or a fine transient dt
    never scrolls the terminal. finish() prints the final state once, in
    green on success or red otherwise, and moves to a new line so
    whatever's printed next starts clean.

    Only redraws in place on a real terminal (stream.isatty()); anything
    else (piped output, a log file, pytest's capsys) gets a single plain
    line from finish() only — intermediate render() calls are a no-op, since
    replaying every redraw into a non-terminal sink would just dump a wall
    of '\\r'-separated text with no benefit over the final summary.
    """

    def __init__(self, width: int = 30, stream: TextIO | None = None):
        self.width = width
        self.stream = stream if stream is not None else sys.stdout
        self.interactive = bool(getattr(self.stream, "isatty", lambda: False)())

    def _bar(self, fraction: float) -> str:
        fraction = min(1.0, max(0.0, fraction))
        filled = int(round(fraction * self.width))
        return "█" * filled + "─" * (self.width - filled)

    def render(self, fraction: float, text: str) -> None:
        if not self.interactive:
            return
        pct = f"{min(1.0, max(0.0, fraction)) * 100:5.1f}%"
        line = f"[{self._bar(fraction)}] {pct}  {text}"
        self.stream.write("\r" + line + "\033[K")
        self.stream.flush()

    def finish(self, text: str, success: bool = True) -> None:
        line = f"[{self._bar(1.0)}] 100.0%  {text}"
        if self.interactive:
            color = _GREEN if success else _RED
            self.stream.write("\r" + f"{color}{_BOLD}{line}{_RESET}" + "\033[K" + "\n")
        else:
            self.stream.write(line + "\n")
        self.stream.flush()
