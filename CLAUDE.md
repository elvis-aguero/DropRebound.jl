Use Test-driven development. Each feature must have a test.

CI is deferred to the cloud, always.

This is a Julia-first repo. Do not introduce Python (or any other language)
tooling, dependencies, or scripts.

## Rheology model contract

A rheology model (Newtonian variants, Oldroyd-B, Carreau-Yasuda, or any
future model) may not be merged without:
1. A from-scratch CAS derivation script under `julia/derivations/` (Julia,
   using Symbolics.jl/QuadGK.jl/SpecialFunctions.jl) with assertions — not
   just a reference to a LaTeX writeup.
2. At least one assertion that cross-checks the derivation against the
   actual running Julia solver (`using DropSolver` directly — no subprocess
   bridge needed), not only against the script's own algebra.
3. Every explanatory comment block readable by a critical reader who has not
   opened the Julia source — state what's being derived and why before the
   math, restate boxed results in plain English, and explain what a failing
   assertion would mean physically.
4. The script executing clean in CI (`.github/workflows/ci.yml`,
   `derivations` job) before the model's production Julia implementation is
   considered done.
