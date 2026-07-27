"""The one-line import surface: `from thermowave.components import X`.

Two things worth pinning down. First, that every public class really is
re-exported — it's easy to add a component and forget the __init__ line, and
the only symptom is an ImportError for the person following the docs. Second,
that importing these packages does not require the optional extras: the
convenience of a flat namespace would not be worth making `import
thermowave.components` fail on a machine without cantera or CoolProp.
"""

import builtins
import importlib
import inspect
import pkgutil
import sys

import pytest

import thermowave.components
import thermowave.core
import thermowave.fluids
import thermowave.maps

_PACKAGES = (
    thermowave.components,
    thermowave.fluids,
    thermowave.maps,
)
_OPTIONAL_DEPENDENCIES = ("cantera", "CoolProp", "matplotlib")


def _public_classes(package):
    """Every public class defined in a module of `package`, as
    (module_name, class_name) pairs."""
    found = []
    for info in pkgutil.iter_modules(package.__path__):
        module = importlib.import_module(f"{package.__name__}.{info.name}")
        for name, obj in vars(module).items():
            if (
                inspect.isclass(obj)
                and obj.__module__ == module.__name__
                and not name.startswith("_")
            ):
                found.append((info.name, name))
    return found


@pytest.mark.parametrize("package", _PACKAGES, ids=lambda p: p.__name__)
def test_every_public_class_is_re_exported(package):
    missing = [
        f"{module}.{name}"
        for module, name in _public_classes(package)
        if not hasattr(package, name)
    ]
    assert not missing, (
        f"{package.__name__} does not re-export: {missing}. Add them to its "
        f"__init__.py so `from {package.__name__} import <Class>` works."
    )


@pytest.mark.parametrize("package", _PACKAGES, ids=lambda p: p.__name__)
def test_all_matches_what_is_actually_exported(package):
    for name in package.__all__:
        assert hasattr(package, name), f"{package.__name__}.__all__ lists absent {name!r}"
    assert package.__all__ == sorted(package.__all__), (
        f"{package.__name__}.__all__ should stay sorted so additions are easy to spot"
    )


def test_core_re_exports_the_names_scripts_actually_reach_for():
    for name in ("Network", "NetworkState", "SolveResult", "TransientResult",
                 "ConvergenceError", "NetworkTopologyError"):
        assert hasattr(thermowave.core, name)


def test_packages_import_without_any_optional_extra_installed(monkeypatch):
    # The optional dependencies are imported inside __init__//functions, never
    # at module level, which is what lets a flat re-export stay safe. Simulate
    # all three being absent and re-import from scratch.
    real_import = builtins.__import__

    def guarded(name, *args, **kwargs):
        if name.split(".")[0] in _OPTIONAL_DEPENDENCIES:
            raise ImportError(f"No module named {name!r} (simulated)")
        return real_import(name, *args, **kwargs)

    for module in list(sys.modules):
        root = module.split(".")[0]
        if root in _OPTIONAL_DEPENDENCIES or root == "thermowave":
            monkeypatch.delitem(sys.modules, module, raising=False)
    monkeypatch.setattr(builtins, "__import__", guarded)

    for name in ("thermowave.components", "thermowave.fluids",
                 "thermowave.core", "thermowave.maps"):
        importlib.import_module(name)


def test_a_component_needing_an_extra_still_reports_it_at_construction():
    # Re-exporting Combustor must not paper over a missing cantera; the
    # error should arrive when it's actually needed, naming the extra.
    if importlib.util.find_spec("cantera") is not None:
        pytest.skip("cantera installed, so construction succeeds")
    with pytest.raises(ImportError, match=r"cantera.*extra"):
        thermowave.components.Combustor(name="comb")


def test_flat_import_gives_the_same_object_as_the_module_path():
    from thermowave.components import SimpleHeatExchanger
    from thermowave.components.simple_heat_exchanger import (
        SimpleHeatExchanger as FromModule,
    )

    assert SimpleHeatExchanger is FromModule
