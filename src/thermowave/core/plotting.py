from __future__ import annotations

from typing import TYPE_CHECKING, Sequence

if TYPE_CHECKING:
    from thermowave.components.base_component import BaseComponent
    from thermowave.core.solver import SolveResult
    from thermowave.core.transient import TransientResult
    from thermowave.maps.characteristic_map import CharacteristicMap

# A plot series: (component, "metric [unit]" key from that component's own
# report_metrics(), optional display label — defaults to "component.metric").
Series = tuple["BaseComponent", str] | tuple["BaseComponent", str, str]

# Okabe-Ito qualitative palette: distinguishable under the common forms of
# color-blindness, unlike matplotlib's default tab10 cycle.
_PALETTE = [
    "#0072B2",  # blue
    "#E69F00",  # orange
    "#009E73",  # green
    "#D55E00",  # vermillion
    "#CC79A7",  # reddish purple
    "#56B4E9",  # sky blue
    "#F0E442",  # yellow
    "#000000",  # black
]

_RC = {
    "axes.prop_cycle": None,  # set below, needs cycler
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.edgecolor": "#333333",
    "axes.labelcolor": "#222222",
    "axes.titlesize": 12,
    "axes.titleweight": "bold",
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
    "font.family": "sans-serif",
    "font.size": 9,
    "legend.frameon": False,
    "xtick.color": "#333333",
    "ytick.color": "#333333",
}


def _import_matplotlib():
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError(
            "Plotting requires the 'plot' extra: pip install thermowave[plot]"
        ) from exc
    return plt


def _theme_rc_params():
    from cycler import cycler

    rc = dict(_RC)
    rc["axes.prop_cycle"] = cycler(color=_PALETTE)
    return rc


def _metrics_for(result: "SolveResult", component: "BaseComponent") -> dict[str, float]:
    """component.report_metrics() evaluated against one SolveResult's state.

    Each component/result pair is looked up freshly (no cross-call cache) —
    NetworkState construction is just dict aliasing, and report_metrics() is
    typically a handful of arithmetic ops, so this is cheap even across many
    steps/components; a cache would only pay off well past the sizes this
    solver targets.
    """
    from thermowave.core.network import NetworkState

    state = NetworkState(
        fluid=result.fluid,
        node_P=result.node_P,
        node_h=result.node_h,
        node_mdot=result.node_mdot,
        params=result.params,
        node_fluid=result.node_fluid,
    )
    metrics = component.report_metrics(state)
    if metrics is None:
        raise ValueError(
            f"Component {component.name!r} has no report_metrics() to plot "
            f"(report_category() is None or it returns None)."
        )
    return metrics


class ThermoPlot:
    """A themed, chainable wrapper around one matplotlib Axes.

    Mutating calls (line/scatter/series/transient/map/bar) return self, so
    calls chain: ThermoPlot(title="...").line(x, y1).line(x, y2).finish().
    The theme (colorblind-safe palette, light background, trimmed spines,
    subtle gridlines) is applied only to the figure this instance owns, via
    a scoped rc_context — it never mutates matplotlib's global rcParams, so
    it can't bleed into unrelated plots a caller makes elsewhere in the same
    process.

    .fig/.ax are public attributes for anything this class doesn't wrap
    directly — call matplotlib on them like any other Axes.
    """

    def __init__(
        self,
        title: str | None = None,
        xlabel: str = "",
        ylabel: str | None = None,
        figsize: tuple[float, float] = (7.0, 4.5),
        ax=None,
    ):
        plt = _import_matplotlib()
        self._plt = plt
        self._owns_figure = ax is None
        if ax is None:
            with plt.rc_context(_theme_rc_params()):
                self.fig, self.ax = plt.subplots(figsize=figsize)
        else:
            self.fig, self.ax = ax.figure, ax

        self.ax.set_xlabel(xlabel)
        if ylabel is not None:
            self.ax.set_ylabel(ylabel)
        if title is not None:
            self.ax.set_title(title)
        self._has_labels = False
        self._twins: list["ThermoPlot"] = []

    def _note_label(self, label) -> None:
        if label is not None:
            self._has_labels = True

    def line(self, x, y, label: str | None = None, **kwargs) -> "ThermoPlot":
        """Plot an XY line series. Accepts any matplotlib Line2D kwargs."""
        kwargs.setdefault("marker", "o")
        kwargs.setdefault("markersize", 3)
        self.ax.plot(x, y, label=label, **kwargs)
        self._note_label(label)
        return self

    def scatter(self, x, y, label: str | None = None, **kwargs) -> "ThermoPlot":
        """Plot an XY scatter series (no connecting line)."""
        self.ax.scatter(x, y, label=label, **kwargs)
        self._note_label(label)
        return self

    def series(
        self,
        x: Sequence[float],
        results: Sequence["SolveResult"],
        series: Sequence[Series],
        **kwargs,
    ) -> "ThermoPlot":
        """Plot one or more component metrics against a shared x-axis.

        x/results must be the same length (one x-value per already-solved
        SolveResult — e.g. TransientResult.times/steps, or your own manual
        sweep over some parameter). Each series entry is (component,
        metric_key) or (component, metric_key, label); metric_key must be a
        key returned by that component's own report_metrics(), i.e. anything
        already visible in a printed report table is plottable with no
        extra wiring.
        """
        if len(x) != len(results):
            raise ValueError(
                f"x has {len(x)} values but results has {len(results)}; must match"
            )
        for entry in series:
            component, metric_key = entry[0], entry[1]
            label = entry[2] if len(entry) > 2 else f"{component.name}.{metric_key}"
            y = [_metrics_for(result, component)[metric_key] for result in results]
            self.line(x, y, label=label, **kwargs)
        return self

    def twin_axis(self, ylabel: str | None = None) -> "ThermoPlot":
        """A second ThermoPlot sharing this one's x-axis but with its own
        independent y-axis (matplotlib's twinx()) — for plotting two metrics
        with different units/scales on one figure, e.g. power [W] against
        efficiency [-]. Any series/line/scatter calls on the returned
        ThermoPlot land on the twin axis; finish() still needs to be called
        only on the original (owning the figure)."""
        twin_ax = self.ax.twinx()
        twin = ThermoPlot.__new__(ThermoPlot)
        twin._plt = self._plt
        twin._owns_figure = False
        twin.fig = self.fig
        twin.ax = twin_ax
        twin._has_labels = False
        twin._twins = []
        if ylabel is not None:
            twin_ax.set_ylabel(ylabel)
        twin_ax.grid(False)
        self._twins.append(twin)
        return twin

    def transient(
        self,
        transient_result: "TransientResult",
        series: Sequence[Series],
        **kwargs,
    ) -> "ThermoPlot":
        """Plot one or more component metrics vs. time from a TransientResult.
        Same series convention as series()/plot_series()."""
        self.ax.set_xlabel(self.ax.get_xlabel() or "t [s]")
        return self.series(transient_result.times, transient_result.steps, series, **kwargs)

    def map(
        self,
        char_map: "CharacteristicMap",
        kind: str = "pressure_ratio",
        operating_points: Sequence[tuple[float, float]] | None = None,
        operating_label: str = "Operating points",
    ) -> "ThermoPlot":
        """Plot a CharacteristicMap's iso-speed curves (pressure ratio or
        efficiency vs. corrected mass flow), one line per corrected speed.

        If operating_points is given (a list of (B, value) pairs — corrected
        mass flow and pressure ratio/efficiency, matching kind), they're
        overlaid as a distinct scatter, connected by a thin line in the
        order given so a traversed operating path (e.g. across a sweep or a
        transient) is visible against the map's own curves.
        """
        for speed, mass_flow, value in char_map.speed_lines(kind):
            self.ax.plot(
                mass_flow, value, linestyle="-", marker=None,
                alpha=0.7, label=f"N={speed:.3g}",
            )
        self._has_labels = True

        if operating_points:
            xs = [p[0] for p in operating_points]
            ys = [p[1] for p in operating_points]
            self.ax.plot(
                xs, ys, linestyle="-", linewidth=1.5, color="#000000",
                marker="D", markersize=5, markerfacecolor="#D55E00",
                markeredgecolor="#000000", label=operating_label, zorder=5,
            )
            self._has_labels = True

        if not self.ax.get_xlabel():
            self.ax.set_xlabel("Corrected mass flow B [-]")
        if self.ax.get_ylabel() is None or self.ax.get_ylabel() == "":
            self.ax.set_ylabel(
                "Pressure ratio [-]" if kind == "pressure_ratio" else "Efficiency [-]"
            )
        return self

    def bar(
        self,
        categories: Sequence[str],
        values: Sequence[float],
        label: str | None = None,
        **kwargs,
    ) -> "ThermoPlot":
        """A bar chart comparing values across discrete categories (e.g. a
        metric across several swept cases or components). Calling bar()
        more than once on the same ThermoPlot draws grouped bars side by
        side (one group of bars per call, redrawn together so groups stay
        centered on each category tick)."""
        if not hasattr(self, "_bar_groups"):
            self._bar_groups: list[tuple[Sequence[float], str | None, dict]] = []
            self._bar_categories = list(categories)
        elif list(categories) != self._bar_categories:
            raise ValueError("bar() categories must match across calls on one ThermoPlot")
        self._bar_groups.append((values, label, kwargs))
        self._note_label(label)

        for container in list(self.ax.containers):
            container.remove()

        import numpy as np

        positions = np.arange(len(self._bar_categories))
        n_groups = len(self._bar_groups)
        width = 0.8 / n_groups
        for i, (group_values, group_label, group_kwargs) in enumerate(self._bar_groups):
            group_offset = (i - (n_groups - 1) / 2.0) * width
            self.ax.bar(
                positions + group_offset, group_values, width=width,
                label=group_label, **group_kwargs,
            )
        self.ax.set_xticks(positions)
        self.ax.set_xticklabels(self._bar_categories)
        return self

    def finish(self, show: bool = True, save_path: str | None = None):
        """Finalize: add a legend (combining entries from any twin_axis()
        axes into one legend on the primary axes) if any series had a
        label, then save and/or show (only meaningful if this ThermoPlot
        owns its figure — i.e. wasn't built via twin_axis() or ax=), and
        return the Axes."""
        any_labels = self._has_labels or any(twin._has_labels for twin in self._twins)
        if any_labels:
            handles, labels = self.ax.get_legend_handles_labels()
            for twin in self._twins:
                twin_handles, twin_labels = twin.ax.get_legend_handles_labels()
                handles += twin_handles
                labels += twin_labels
            self.ax.legend(handles, labels)
        self.ax.grid(True, alpha=0.3)

        if self._owns_figure:
            if save_path is not None:
                self._plt.savefig(save_path, bbox_inches="tight")
            if show:
                self._plt.show()
        return self.ax


def plot_series(
    x: Sequence[float],
    results: Sequence["SolveResult"],
    series: Sequence[Series],
    xlabel: str = "",
    ylabel: str | None = None,
    title: str | None = None,
    ax=None,
    show: bool = True,
    save_path: str | None = None,
):
    """Plot one or more component metrics against a shared x-axis.

    Thin backward-compatible wrapper around ThermoPlot.series(); see that
    method for the series/x/results contract. Prefer
    ThermoPlot(...).series(...).finish(...) directly in new code for access
    to the rest of ThermoPlot (twin_axis, scatter, map, bar, ...).
    """
    plot = ThermoPlot(title=title, xlabel=xlabel, ylabel=ylabel, ax=ax)
    plot.series(x, results, series)
    return plot.finish(show=show, save_path=save_path)
