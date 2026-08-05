<!--
Thanks for the contribution. See CONTRIBUTING.md for the full guide —
this template just mirrors its checklist so nothing gets missed.
-->

## What this changes and why

<!-- What problem does this solve, or what capability does it add?
     Link the issue this addresses, if there is one. -->

Closes #

## Type of change

- [ ] Bug fix
- [ ] New feature / component
- [ ] Breaking change
- [ ] Docs only
- [ ] Other (refactor, CI, tooling, ...)

## Checklist

- [ ] `pytest` passes locally
- [ ] `ruff check src/ tests/` passes
- [ ] `mypy src/thermowave` doesn't add new errors (pre-existing ones are fine — see CI's `continue-on-error`)
- [ ] **Bug fix:** added a regression test under `tests/` that fails before this change and passes after
- [ ] **New component:** added a page under `docs/components/` and wired it into that section's `index.md`
- [ ] **Breaking change:** called out explicitly below, and noted in `docs/changelog.md` (and the README's Roadmap if architecturally significant)
- [ ] **Docs touched:** built the Sphinx site locally (`sphinx-build -b html docs docs/_build/html`) with no new warnings

## Breaking changes

<!-- If none, delete this section. Otherwise: what breaks, and what does
     a caller need to change? -->

## Test plan

<!-- How did you verify this works, beyond the test suite — a benchmark
     script, a manual repro, a before/after comparison? -->
