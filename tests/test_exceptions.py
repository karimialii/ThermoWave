import pytest

from thermowave.core.exceptions import (
    ConvergenceError,
    FluidRangeError,
    NetworkTopologyError,
)


@pytest.mark.parametrize(
    "exc_type",
    [ConvergenceError, FluidRangeError, NetworkTopologyError],
)
def test_exception_is_exception_subclass_with_message(exc_type):
    exc = exc_type("something went wrong")
    assert isinstance(exc, Exception)
    assert str(exc) == "something went wrong"
