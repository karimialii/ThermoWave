"""Network, solver, and result types, re-exported for a one-line import.

    from thermowave.core import Network, ConvergenceError
"""

from thermowave.core.exceptions import (
    ConvergenceError,
    FluidRangeError,
    NetworkTopologyError,
)
from thermowave.core.network import Connection, Network, NetworkState
from thermowave.core.progress import ProgressBar
from thermowave.core.settings import Settings
from thermowave.core.solver import SolveResult, Solver
from thermowave.core.transient import TransientResult

__all__ = [
    "Connection",
    "ConvergenceError",
    "FluidRangeError",
    "Network",
    "NetworkState",
    "NetworkTopologyError",
    "ProgressBar",
    "Settings",
    "SolveResult",
    "Solver",
    "TransientResult",
]
