import numpy as np
import pytest

from thermowave.core.exceptions import ConvergenceError
from thermowave.core.solver import _finite_difference_jacobian, newton_solve


def test_finite_difference_jacobian_matches_analytic_for_linear_system():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([2.0 * x[0] + 3.0 * x[1] - 1.0, -x[0] + 5.0 * x[1] + 2.0])

    x = np.array([1.0, 1.0])
    F = residual_fn(x)
    J = _finite_difference_jacobian(residual_fn, x, F)
    expected_J = np.array([[2.0, 3.0], [-1.0, 5.0]])
    assert np.allclose(J, expected_J, atol=1e-4)


def test_newton_solve_converges_to_known_linear_root():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 3.0, 2.0 * x[1] + 4.0])

    x0 = np.array([0.0, 0.0])
    x_sol, iterations, residual_norm = newton_solve(residual_fn, x0, tol=1e-9)
    assert np.allclose(x_sol, [3.0, -2.0], atol=1e-6)
    assert iterations >= 1
    assert residual_norm < 1e-9


def test_newton_solve_converges_to_known_nonlinear_root():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] ** 2 - 4.0])

    x0 = np.array([3.0])
    x_sol, _iterations, residual_norm = newton_solve(residual_fn, x0, tol=1e-9)
    assert np.allclose(x_sol, [2.0], atol=1e-6)
    assert residual_norm < 1e-9


def test_newton_solve_raises_convergence_error_when_residual_never_vanishes():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([1.0, 1.0])

    x0 = np.array([0.0, 0.0])
    with pytest.raises(ConvergenceError):
        newton_solve(residual_fn, x0, tol=1e-9, max_iter=3)


def test_newton_solve_raises_convergence_error_on_singular_jacobian():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - x[1], x[0] - x[1]])

    x0 = np.array([0.0, 1.0])
    with pytest.raises(ConvergenceError):
        newton_solve(residual_fn, x0, tol=1e-9, max_iter=5)


def test_newton_solve_applies_clamp_fn_each_iteration():
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 100.0])

    def clamp_fn(x: np.ndarray) -> np.ndarray:
        return np.maximum(x, np.array([0.0]))

    x0 = np.array([50.0])
    x_sol, _iterations, _residual_norm = newton_solve(
        residual_fn, x0, tol=1e-9, clamp_fn=clamp_fn
    )
    assert np.allclose(x_sol, [100.0], atol=1e-6)


from thermowave.components.base_component import BaseComponent
from thermowave.components.pipe import Pipe
from thermowave.components.simple_compressor import SimpleCompressor
from thermowave.components.sink import Sink
from thermowave.components.source import Source
from thermowave.components.valve import Valve
from thermowave.core.network import Network, NetworkState
from thermowave.core.solver import Solver
from thermowave.fluids.ideal_gas import IdealGasFluid

_AIR = IdealGasFluid(name="air", R=287.05, cp=1005.0)


def _simple_network() -> Network:
    network = Network(fluid=_AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    pipe = Pipe(name="p1", L=5.0, D=0.2, f=0.02)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")
    return network


def test_solver_solve_returns_converged_result_for_single_pipe_network():
    network = _simple_network()

    result = Solver(network).solve(tol=1e-6, max_iter=50)

    assert result.converged is True
    assert result.node_P["src.out"] == 101325.0  # fixed by Source, untouched
    assert result.node_P["p1.out"] < 101325.0  # pressure drops across the pipe
    assert result.node_mdot["src.out"] == 1.0
    assert result.node_mdot["p1.out"] == pytest.approx(1.0)


def test_newton_solve_verbose_prints_progress_summary(capsys):
    # capsys's captured stream isn't a real terminal, so ProgressBar's
    # in-place '\r' redraws no-op (see test_progress.py for that path) —
    # this only checks the header + final one-line summary, which print
    # unconditionally regardless of interactivity.
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 3.0, 2.0 * x[1] + 4.0])

    x0 = np.array([0.0, 0.0])
    newton_solve(residual_fn, x0, tol=1e-9, verbose=True)

    out = capsys.readouterr().out
    assert "Newton-Raphson solve: 2 unknowns" in out
    assert "residual norm" in out
    assert "Converged in" in out


def test_newton_solve_quiet_by_default_prints_nothing(capsys):
    def residual_fn(x: np.ndarray) -> np.ndarray:
        return np.array([x[0] - 3.0])

    newton_solve(residual_fn, np.array([0.0]), tol=1e-9)

    assert capsys.readouterr().out == ""


def test_network_solve_verbose_prints_progress_summary(capsys):
    network = _simple_network()

    network.solve(tol=1e-6, max_iter=50, verbose=True)

    out = capsys.readouterr().out
    assert "Newton-Raphson solve" in out
    assert "Converged in" in out


def test_solve_result_print_report_shows_status_and_node_table(capsys):
    network = _simple_network()

    result = Solver(network).solve(tol=1e-6, max_iter=50)
    result.print_report()

    out = capsys.readouterr().out
    assert "CONVERGED" in out
    assert "src.out" in out
    assert "p1.out" in out
    assert "P [Pa]" in out
    assert "T [K]" in out
    assert "h [J/kg]" in out
    assert "mdot [kg/s]" in out
    assert "1.0000" in out  # mdot value, repeated per row


def test_solve_result_print_report_omits_pipe_internal_nodes(capsys):
    network = Network(fluid=_AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    pipe = Pipe(name="p1", L=6.0, D=0.2, f=0.02, n_elem=3)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(pipe)
    network.add_component(snk)
    network.connect(src, "out", pipe, "in")
    network.connect(pipe, "out", snk, "in")

    result = Solver(network).solve(tol=1e-6, max_iter=50)

    # Internal discretization nodes are solved and readable...
    assert "p1__mid1" in result.node_P
    assert "p1__mid2" in result.node_P

    # ...but must not appear in the printed table, only the component ports.
    result.print_report()
    out = capsys.readouterr().out
    assert "src.out" in out
    assert "p1.out" in out
    assert "p1__mid1" not in out
    assert "p1__mid2" not in out


def test_solve_result_print_report_omits_components_table_when_nothing_reportable(capsys):
    # Source/Pipe/Sink is the whole network: none of them implement
    # report_metrics(), so the "Components" table should not appear at all.
    network = _simple_network()

    result = Solver(network).solve(tol=1e-6, max_iter=50)
    result.print_report()

    out = capsys.readouterr().out
    assert "Components" not in out
    assert "Nodes" in out


def test_solve_result_print_report_shows_components_table_for_valve_and_compressor():
    network = Network(fluid=_AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    gamma = 1005.0 / (1005.0 - 287.05)
    comp = SimpleCompressor(name="c1", PR=2.0, eta_s=0.8, gamma=gamma)
    valve = Valve(name="v1", D=0.2, K=1.0)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(comp)
    network.add_component(valve)
    network.add_component(snk)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", valve, "in")
    network.connect(valve, "out", snk, "in")

    result = Solver(network).solve(tol=1e-6, max_iter=50)

    state = NetworkState(
        fluid=result.fluid,
        node_P=result.node_P,
        node_h=result.node_h,
        node_mdot=result.node_mdot,
    )
    comp_metrics = comp.report_metrics(state)
    valve_metrics = valve.report_metrics(state)

    assert comp_metrics["eta_s [-]"] == 0.8
    assert comp_metrics["PR [-]"] == 2.0
    assert comp_metrics["power [W]"] > 0  # compressor consumes work

    assert "eta_s [-]" not in valve_metrics  # not applicable to a valve
    assert valve_metrics["power [W]"] == pytest.approx(0.0, abs=1e-6)  # isenthalpic


def test_solve_result_print_report_component_table_lists_names_not_source_sink(capsys):
    network = Network(fluid=_AIR)
    src = Source(name="src", P=101325.0, T=300.0, mdot=1.0)
    gamma = 1005.0 / (1005.0 - 287.05)
    comp = SimpleCompressor(name="c1", PR=2.0, eta_s=0.8, gamma=gamma)
    snk = Sink(name="snk")
    network.add_component(src)
    network.add_component(comp)
    network.add_component(snk)
    network.connect(src, "out", comp, "in")
    network.connect(comp, "out", snk, "in")

    result = Solver(network).solve(tol=1e-6, max_iter=50)
    result.print_report()

    out = capsys.readouterr().out
    assert "TurboMachinery" in out
    assert "power [W]" in out
    assert "eta_s [-]" in out
    assert "PR [-]" in out
    components_section = out.split("TurboMachinery")[1].split("Nodes")[0]
    assert "c1" in components_section
    assert "src" not in components_section
    assert "snk" not in components_section


class _PHBoundaryComponent(BaseComponent):
    """Test double: pins P and h via residuals() instead of
    fixed_node_values() (a legitimate alternative — fixed_node_values() is
    just a convenience for the common "P/h known outright" case). Used to
    build a network where NOTHING fixes (P, h) via fixed_node_values(),
    exercising Solver.solve()'s DEFAULT_P_FALLBACK/DEFAULT_H_FALLBACK path
    — the initial (P, h) guess used only when fixed_nodes is empty.

    Note this path is unreachable through Network.solve() (its own
    validate_topology() call rejects a network with no fixed_node_values()
    boundary before Solver ever runs), so this test goes through
    Solver(network).solve() directly, same as this file's other tests.
    """

    def __init__(self, name: str, node_name: str, P_target: float, h_target: float, mdot: float):
        self.name = name
        self._node = node_name
        self.P_target = P_target
        self.h_target = h_target
        self.mdot = mdot

    def ports(self) -> dict[str, str]:
        return {"out": self._node}

    def fixed_node_mdot(self) -> dict[str, float]:
        return {self._node: self.mdot}

    def residuals(self, state: "NetworkState") -> list[float]:
        P, h = state.node(self._node)
        return [P - self.P_target, h - self.h_target]


def test_solver_uses_default_fallback_guess_when_network_has_no_fixed_ph_nodes():
    network = Network(fluid=_AIR)
    boundary = _PHBoundaryComponent(
        name="b", node_name="b.out", P_target=180000.0, h_target=250000.0, mdot=1.0,
    )
    network.add_component(boundary)
    assert network._fixed_node_values() == {}  # confirms the fallback path is exercised

    result = Solver(network).solve(tol=1e-9, max_iter=50)

    assert result.converged is True
    assert result.node_P["b.out"] == pytest.approx(180000.0, rel=1e-9)
    assert result.node_h["b.out"] == pytest.approx(250000.0, rel=1e-9)
