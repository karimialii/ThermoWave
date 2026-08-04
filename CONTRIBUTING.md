# Contributing to ThermoWave

Thanks for considering a contribution. This document covers how to set up a
dev environment, what CI checks, and what a good pull request looks like.

## Development setup

```bash
git clone https://github.com/karimialii/ThermoWave.git
cd ThermoWave
pip install -e ".[full]"   # package + CoolProp/Cantera/plot + dev/docs tooling
```

`[full]` pulls in everything, including the optional real-fluid
(`CoolProp`) and equilibrium-chemistry (`Cantera`) backends. If you only
need the core solver, `pip install -e ".[dev]"` is enough for tests and
linting — tests that need CoolProp/Cantera are auto-skipped
(`pytest.importorskip`) when those extras aren't installed.

## Before opening a PR

Run the same checks CI runs (`.github/workflows/ci.yml`):

```bash
pytest                          # 212+ tests
ruff check src/ tests/          # lint
mypy src/thermowave             # type check (non-blocking in CI, but please don't add new errors)
```

If you touch anything under `docs/`, build the Sphinx site locally before
pushing:

```bash
pip install -e ".[docs]"
sphinx-build -b html docs docs/_build/html
```

## Making a change

- **New component?** Read
  [Writing a new component](https://github.com/karimialii/ThermoWave#writing-a-new-component)
  in the README first — it covers the residual-contribution pattern every
  component follows. Add a matching page under `docs/components/` and wire
  it into that section's `index.md`.
- **Bug fix?** Add a regression test under `tests/` that fails before your
  fix and passes after. Don't fix the symptom without a test pinning the
  behavior.
- **Breaking change?** Call it out explicitly in your PR description and
  add a note to `docs/changelog.md` and, if it's architecturally
  significant, the "Roadmap" section of the README — see existing
  "— landed." entries there for the tone/format to match.
- **Docs-only change?** Still welcome as its own PR — no need to bundle it
  with a code change.

## Code style

- Follow the conventions already in the file you're editing before
  reaching for a personal preference — this codebase favors small,
  focused components and residual functions that read like the physics
  they encode.
- `ruff` and `mypy` configuration lives in `pyproject.toml`; don't
  add per-file ignores unless there's no reasonable alternative, and
  explain why in a comment next to the ignore.
- Comments should explain *why*, not *what* — non-obvious constraints,
  workarounds, or invariants, not a restatement of the code.

## Commit messages and PRs

- Keep commits focused; a commit message's first line should say what
  changed, not just "fix bug" or "update file."
- Reference the relevant issue number if one exists.
- Small, reviewable PRs are preferred over large ones that bundle
  unrelated changes.

## Reporting bugs / requesting features

Open a [GitHub issue](https://github.com/karimialii/ThermoWave/issues)
with a minimal reproducing network (component wiring + parameters) for
bugs, or a description of the use case for feature requests.

## Security issues

Please don't open a public issue for a security vulnerability — see
[SECURITY.md](https://github.com/karimialii/ThermoWave/blob/main/SECURITY.md).

## Code of Conduct

This project follows the
[Code of Conduct](https://github.com/karimialii/ThermoWave/blob/main/CODE_OF_CONDUCT.md).
By participating, you're expected to uphold it.
