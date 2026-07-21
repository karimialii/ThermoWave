class ConvergenceError(Exception):
    """Raised when the Newton-Raphson solver fails to converge."""


class FluidRangeError(Exception):
    """Raised when a fluid property call falls outside a physically valid range."""


class NetworkTopologyError(Exception):
    """Raised when a network's component graph is invalid or under/over-determined."""
