# Security Policy

ThermoWave is a headless numerical library (no network listeners, no
credential handling, no untrusted-input parsing beyond reading local map
files you provide it), so its attack surface is small. Still, if you find
a genuine security issue — e.g. a crafted `.cop`/`.tur` characteristic-map
file or fluid-mixture spec that causes memory corruption or arbitrary code
execution rather than just a `ValueError` — please report it privately
rather than opening a public issue.

## Reporting a vulnerability

Use GitHub's
[private vulnerability reporting](https://github.com/karimialii/ThermoWave/security/advisories/new)
for the repository. If that's unavailable, open a regular issue asking a
maintainer to reach out on a private channel, without describing the
vulnerability itself.

Please include:

- A minimal reproducing case
- The affected version (`pip show thermowave`)
- The impact you'd expect it to have

## Supported versions

ThermoWave is pre-1.0 (`0.x`); only the latest released version is
supported. Fixes land on `main` and ship in the next release rather than
being backported.
