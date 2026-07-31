# DropRebound.jl

A Julia spectral solver for drop impact and rebound (Newtonian, Oldroyd-B, and
Carreau-Yasuda shear-thinning rheology).

## CAS Derivations

Every rheology model in this repo is required (see `CLAUDE.md`) to have a
from-scratch computer-algebra derivation, with assertions, that cross-checks
against the live running solver. The pages under **CAS Derivations** are that
same derivation source -- rendered for reading, not a separate writeup. The
math you see here is exactly what CI executes and gates on; nothing is
transcribed by hand into a different, driftable document.

See the [API Reference](@ref) for the production solver's own docstrings.
