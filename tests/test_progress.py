import io

from thermowave.core.progress import ProgressBar


class _FakeStream(io.StringIO):
    """io.StringIO reports isatty() == False; ProgressBar's in-place
    redraws are gated on isatty(), so a fake terminal needs to override it
    to exercise that path."""

    def __init__(self, interactive: bool):
        super().__init__()
        self._interactive = interactive

    def isatty(self) -> bool:
        return self._interactive


def test_progress_bar_render_is_noop_on_non_interactive_stream():
    stream = _FakeStream(interactive=False)
    bar = ProgressBar(stream=stream)

    bar.render(0.5, "halfway")

    assert stream.getvalue() == ""


def test_progress_bar_finish_prints_one_plain_line_on_non_interactive_stream():
    stream = _FakeStream(interactive=False)
    bar = ProgressBar(stream=stream)

    bar.render(0.3, "in progress")  # no-op, shouldn't leak into finish's line
    bar.finish("Converged in 5 iterations", success=True)

    out = stream.getvalue()
    assert out.count("\n") == 1
    assert out.endswith("\n")
    assert "Converged in 5 iterations" in out
    assert "100.0%" in out
    assert "\033[" not in out  # no ANSI codes on a non-terminal sink


def test_progress_bar_render_redraws_in_place_on_interactive_stream():
    stream = _FakeStream(interactive=True)
    bar = ProgressBar(width=10, stream=stream)

    bar.render(0.0, "start")
    bar.render(0.5, "halfway")
    bar.render(1.0, "almost done")

    out = stream.getvalue()
    # Every render is a '\r' + erase-to-end-of-line redraw, never a newline —
    # that's what keeps it from scrolling the terminal.
    assert "\n" not in out
    assert out.count("\r") == 3
    assert "\033[K" in out
    assert "halfway" in out
    assert "almost done" in out


def test_progress_bar_render_reflects_fraction_in_bar_and_percentage():
    stream = _FakeStream(interactive=True)
    bar = ProgressBar(width=10, stream=stream)

    bar.render(0.5, "")
    out = stream.getvalue()

    assert "50.0%" in out
    assert "█████─────" in out  # half-filled, 10-wide


def test_progress_bar_render_clamps_out_of_range_fractions():
    stream = _FakeStream(interactive=True)
    bar = ProgressBar(width=4, stream=stream)

    bar.render(-0.5, "below zero")
    bar.render(2.0, "above one")

    out = stream.getvalue()
    assert "0.0%" in out
    assert "100.0%" in out
    assert "────" in out  # the -0.5 render: fully empty
    assert "████" in out  # the 2.0 render: fully filled


def test_progress_bar_finish_colors_green_on_success_on_interactive_stream():
    stream = _FakeStream(interactive=True)
    bar = ProgressBar(stream=stream)

    bar.finish("Converged", success=True)

    out = stream.getvalue()
    assert "\033[32m" in out  # green
    assert "\033[0m" in out  # reset
    assert out.endswith("\n")


def test_progress_bar_finish_colors_red_on_failure_on_interactive_stream():
    stream = _FakeStream(interactive=True)
    bar = ProgressBar(stream=stream)

    bar.finish("Failed to converge", success=False)

    out = stream.getvalue()
    assert "\033[31m" in out  # red
    assert "Failed to converge" in out
