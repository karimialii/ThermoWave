from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from thermowave.core.solver import SolveResult

_ITER_W = 8
_NORM_W = 18


def _iteration_rule(joint: str) -> str:
    return "─" * _ITER_W + joint + "─" * _NORM_W + joint + "─" * _NORM_W


def print_system_summary(n_unknowns: int, n_equations: int) -> None:
    if n_unknowns == n_equations:
        print(
            f"Network system: {n_unknowns} unknowns, {n_equations} equations "
            f"-> square, solvable"
        )
    else:
        print(
            f"Network system: {n_unknowns} unknowns, {n_equations} equations "
            f"-> NOT square, not solvable"
        )


def print_iteration_header(n_unknowns: int, tol: float, max_iter: int) -> None:
    print(f"Newton-Raphson solve: {n_unknowns} unknowns, tol={tol:.2e}, max_iter={max_iter}")
    print(_iteration_rule("┬"))
    print(f"{'iter':^{_ITER_W}}│{'residual norm':^{_NORM_W}}│{'step norm':^{_NORM_W}}")
    print(_iteration_rule("┼"))


def print_iteration_row(iteration: int, residual_norm: float, step_norm: float) -> None:
    print(f"{iteration:^{_ITER_W}d}│{residual_norm:^{_NORM_W}.4e}│{step_norm:^{_NORM_W}.4e}")


def print_iteration_footer(
    converged: bool, iterations: int, residual_norm: float, tol: float
) -> None:
    print(_iteration_rule("┴"))
    if converged:
        print(
            f"Converged in {iterations} iterations "
            f"(residual norm = {residual_norm:.4e} < tol = {tol:.2e})"
        )
    else:
        print(
            f"Failed to converge after {iterations} iterations "
            f"(residual norm = {residual_norm:.4e}, tol = {tol:.2e})"
        )


# One column layout + title per report_category() — components sharing a
# category are shown in one table, in this order; a category with no
# components present (e.g. "generator" until a Generator component exists)
# is simply omitted rather than printed empty.
_CATEGORY_TABLES: list[tuple[str, str, list[str]]] = [
    ("turbomachinery", "TurboMachinery", ["power [W]", "eta_s [-]", "PR [-]", "N [rev/min]"]),
    (
        "heat_exchanger",
        "Heat Exchangers",
        ["power [W]", "PR_hot [-]", "PR_cold [-]", "T_hot_in [K]", "T_cold_in [K]"],
    ),
    (
        "combustor",
        "Combustors",
        [
            "power [W]", "mdot_fuel [kg/s]", "PR [-]", "T_out [K]",
            "X_CO2 [-]", "X_H2O [-]",
        ],
    ),
    ("controller", "Controllers", ["target [-]", "measured [-]", "error [-]"]),
    (
        "phase_change",
        "Phase-Change HX",
        ["power [W]", "PR [-]", "x_out [-]", "T_sat [K]", "T_out [K]", "dT_sat [K]", "pinch [K]"],
    ),
    (
        "nozzle",
        "Nozzles",
        [
            "V [m/s]", "Mach [-]", "mdot [kg/s]", "PR [-]", "choked [-]",
            "Mach_exit [-]", "P_exit_ideal [Pa]",
        ],
    ),
    (
        "shaft",
        "Shafts",
        ["N [rev/min]", "power [W]", "eta [-]", "inertia [kg*m^2]"],
    ),
    (
        "tank",
        "Tanks",
        ["P [Pa]", "T [K]", "mdot_in [kg/s]", "mdot_out [kg/s]", "V [m^3]"],
    ),
    (
        "drum",
        "Steam Drums",
        ["P [Pa]", "T_sat [K]", "level [-]", "mdot_steam [kg/s]", "mdot_water [kg/s]", "V [m^3]"],
    ),
    ("sensor", "Sensors", ["P [Pa]", "T [K]", "mdot [kg/s]"]),
    ("generator", "Generators", ["power [W]", "eta [-]", "N [rev/min]"]),
    ("motor", "Motors", ["power [W]", "eta [-]", "N [rev/min]"]),
    ("thermal_mass", "Thermal Masses", ["T [K]", "Q_net [W]"]),
    ("heat_transfer", "Heat Transfer Paths", ["Q [W]", "T_a [K]", "T_b [K]"]),
]


def _table_row(cells: list[str]) -> str:
    return "│".join(cells)


def _format_node_table(result: "SolveResult") -> list[str]:
    node_col_w = max([len("node")] + [len(n) for n in result.node_order]) + 2
    num_w = 14
    columns = ["P [Pa]", "T [K]", "h [J/kg]", "mdot [kg/s]"]

    header_row = _table_row(
        [f"{'node':<{node_col_w}}"] + [f"{c:>{num_w}}" for c in columns]
    )
    sep = "─" * node_col_w + "┼" + "┼".join("─" * num_w for _ in columns)

    lines = [header_row, sep]
    for name in result.node_order:
        P = result.node_P[name]
        h = result.node_h[name]
        T = result.node_fluid.get(name, result.fluid).temperature_ph(P, h)
        mdot = result.node_mdot.get(name)
        mdot_cell = "-" if mdot is None else f"{mdot:.4f}"
        lines.append(
            _table_row(
                [f"{name:<{node_col_w}}"]
                + [
                    f"{P:>{num_w}.2f}",
                    f"{T:>{num_w}.2f}",
                    f"{h:>{num_w}.2f}",
                    f"{mdot_cell:>{num_w}}",
                ]
            )
        )
    return lines


def _render_metrics_table(
    rows: list[tuple[str, dict[str, float]]], columns: list[str]
) -> list[str]:
    name_col_w = max([len("component")] + [len(name) for name, _ in rows]) + 2
    num_w = 14

    def cell(value: float | None) -> str:
        if value is None:
            return f"{'-':>{num_w}}"
        return f"{value:>{num_w}.4g}"

    header_row = _table_row(
        [f"{'component':<{name_col_w}}"] + [f"{col:>{num_w}}" for col in columns]
    )
    sep = "─" * name_col_w + "┼" + "┼".join("─" * num_w for _ in columns)

    lines = [header_row, sep]
    for name, metrics in rows:
        lines.append(
            _table_row([f"{name:<{name_col_w}}"] + [cell(metrics.get(col)) for col in columns])
        )
    return lines


def _format_component_tables(result: "SolveResult") -> list[tuple[str, list[str]]]:
    """One (title, lines) table per report_category() present among the
    network's components — e.g. TurboMachinery, Heat Exchangers, Controllers
    (Source, Sink, and Pipe report no category and are left out entirely).
    Returns [] if no component has anything to report.
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

    rows_by_category: dict[str, list[tuple[str, dict[str, float]]]] = {}
    for component in result.components:
        category = component.report_category()
        if category is None:
            continue
        metrics = component.report_metrics(state)
        if metrics is not None:
            rows_by_category.setdefault(category, []).append((component.name, metrics))

    tables: list[tuple[str, list[str]]] = []
    for category, title, columns in _CATEGORY_TABLES:
        rows = rows_by_category.get(category)
        if rows:
            tables.append((title, _render_metrics_table(rows, columns)))
    return tables


def format_solve_report(result: "SolveResult") -> str:
    """Render a SolveResult as a summary + components table + per-node table."""
    status = "CONVERGED" if result.converged else "NOT CONVERGED"
    summary_lines = [
        "ThermoWave Network Solution",
        f"  status      : {status}",
        f"  iterations  : {result.iterations}",
        f"  residual    : {result.residual_norm:.4e}",
    ]

    component_tables = _format_component_tables(result)
    node_lines = _format_node_table(result)

    sections = list(summary_lines)
    for title, lines in component_tables:
        sections += ["", title, *lines]
    sections += ["", "Nodes", *node_lines]

    rule_w = max(len(line) for line in sections)
    rule = "═" * rule_w
    return "\n".join([rule, *sections, rule])
