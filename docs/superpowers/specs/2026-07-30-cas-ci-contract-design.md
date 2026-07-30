# CAS-before-code contract + CI, and an enhanced Oldroyd-B derivation

Date: 2026-07-30

## Background

`main` already contains a full Carreau shear-thinning rheology (`julia/src/st_extension.jl`,
`julia/test/test_carreau.jl`) merged from `dev/shear-thinning`, along with a symbolic
derivation notebook (`notebooks/shear_thinning_derivation.ipynb`, 19 sympy/mpmath
assertions). Oldroyd-B (`julia/src/ob_extension.jl`) has tests
(`test_ob.jl`, `test_ob_eigenvalue.jl`) and a LaTeX writeup
(`docs/section_oldroydB.tex`) but no from-scratch symbolic derivation notebook.

Nothing in the repo runs any of this automatically — there is no `.github/workflows`
at all. The intended contract ("a rheology model is not implemented until its math has
been derived and checked by a CAS, and that check runs in CI") exists only as intent,
not as an enforced process. This was surfaced while scoping an unrelated request (a
Cross-fluid feasibility study) — Cross work is explicitly deferred until this contract
is real.

## Goals

1. Stand up CI (GitHub Actions) that runs the Julia test suite and executes every
   derivation notebook headlessly, so a broken test *or* a broken mathematical
   assertion fails the build.
2. Make the notebooks actually verify the implementation, not just their own algebra:
   each notebook gets at least one cell that shells out to the real Julia solver,
   runs a small simulation, and asserts the numeric result matches the notebook's
   symbolic/asymptotic prediction.
3. Write a substantially deeper Oldroyd-B derivation notebook — not a structural
   mirror of Carreau's, but one that derives the polymer-stress ODE and the exact
   Jacobian coefficients used in `ob_extension.jl` from the upper-convected Maxwell
   constitutive law, and checks them line-for-line against the shipped code.
3b. Strengthen the existing Carreau notebook with the same live cross-check (it
    currently stops at "Next steps: ... 3. Validate" without doing so).
4. Document the contract explicitly (CLAUDE.md) so it's enforced going forward.

Out of scope: any Cross-model work. That starts only after this lands.

## Notebook writing standard

Every derivation notebook (new or touched) must serve two readers at once, not just
one:

- **The CAS itself** — self-consistency and strict derivation: symbolic algebra,
  assertions, numeric cross-checks against the live solver. This is what CI runs.
- **A critical human reader** — someone auditing the physics who did not write the
  code. For every derivation step this means: state *what* is being computed and
  *why* it matters *before* the symbolic manipulation, restate every boxed/derived
  result in plain English immediately after it (not just LaTeX), and say explicitly
  what a passing vs. failing assertion would mean physically — not only
  algebraically. An assertion cell with no surrounding prose is a gap, even if it
  passes.

Concretely: each markdown cell should read as a standalone explanation someone could
follow without opening `ob_extension.jl` or `st_extension.jl` first, and each
assertion needs a one-line "this checks that ___, which would fail if ___" comment,
not just a bare `assert`. This standard applies to the new Oldroyd-B notebook in
full, and to the new sections being added to the Carreau notebook (component 4) —
existing Carreau sections are already reasonably narrated and don't need a rewrite
pass unless a specific section is found thin during that work.

## Components

### 1. `.github/workflows/ci.yml`

Two jobs, triggered on push and pull_request to `main`:

- **`julia-tests`**: `julia-actions/setup-julia@v2` (1.12), `julia-actions/cache@v2`,
  then `julia --project=julia -e 'using Pkg; Pkg.test()'`. Covers the full existing
  suite (Newtonian, OB, Carreau, contact, postprocessing — 148+ tests plus whatever
  Carreau added).
- **`notebooks`**: needs *both* toolchains, since notebooks now call Julia directly.
  `julia-actions/setup-julia@v2` + `astral-sh/setup-uv@v3`, `uv sync` (uses the
  existing `pyproject.toml`/`uv.lock`), then for each notebook under `notebooks/`:
  `uv run jupyter nbconvert --to notebook --execute --stdout > /dev/null`.
  A failed Python `assert` raises inside the kernel, which nbconvert surfaces as a
  non-zero exit — no special assertion-catching machinery needed.

### 2. Live Julia cross-check bridge

A small shared helper, `notebooks/_juliabridge.py`, wrapping
`subprocess.run(["julia", "--project=julia", "-e", <script>], capture_output=True)`
and JSON-decoding stdout (the inline Julia script prints a `JSON3`/manual
`println(json)` line of the requested numbers — e.g. run `solve_drop!` for a
lightly-damped, no-contact case at given `(Oh, M)` or `(Oh, De1, beta_s)` or
`(Oh, lambda_c, n)`, fit the log-decrement and frequency from the saved `A[l]`
time series, and return `{decay, freq}`). Both derivation notebooks import this
and add a final "cross-check against the shipped solver" section asserting
symbolic/asymptotic prediction ≈ live Julia output within a stated tolerance
(loose, e.g. 5-10%, consistent with the tolerances already used in
`test_ob_eigenvalue.jl`/README's "errors below 5%").

This makes each notebook a genuine regression check on the implementation: if
someone changes `ob_extension.jl` or `st_extension.jl` in a way that breaks the
physics, the notebook's own CI run fails, not just the Julia unit tests.

### 3. `notebooks/oldroyd_b_derivation.ipynb` (new, from scratch)

Mirrors the *rigor*, not the outline, of the Carreau notebook. Sections:

1. Upper-convected Maxwell polymer stress + Newtonian solvent (Oldroyd-B)
   constitutive law, linearized about the spherical base state, using the same
   Reid velocity field (`U(x)`, `V(x)`) already derived for Newtonian/Carreau.
2. Derive the polymer stress auxiliary equation
   `De₁ Ṡₙ = (1-β_s)Ȧₙ - Sₙ` from the linearized upper-convected derivative — the
   notebook currently doesn't exist, so this is presently just asserted in code
   comments (`ob_extension.jl:41-45`) with no independent derivation anywhere.
3. Derive the modified mode equation
   (`Äₙ + 2Oh(n-1)(2n+1)(β_s Ȧₙ + Sₙ) + n(n+2)(n-1)Aₙ + n Bₙ = 0`) from the
   stress balance at the free surface.
4. **Code parity assertions**: symbolically reduce the BDF-discretized Block-S
   Jacobian entries (`(ak + dt/De1)`, `-dt(1-β_s)/De1`, `D2_diag = 2·Oh·(n-1)(2n+1)`)
   and assert they equal the literal expressions in `julia/src/ob_extension.jl`
   (lines 63-66, 156-166) — a direct side-by-side check against the source, which
   nothing today does.
5. Reproduce the Oldroyd-B characteristic equation (the one
   `test_ob_eigenvalue.jl` already root-finds independently in Julia) symbolically
   and check several `(Oh, De1, beta_s)` roots against that test's expected values.
6. Newtonian limits: `De1 → 0` and `beta_s → 1` both collapse exactly to the Reid
   result (mirrors Carreau's `ε_ST=0` assertions).
7. Live cross-check (component 2): run the real OB solver, compare decay
   rate/frequency to the notebook's prediction.

### 4. Carreau notebook enhancement

Add one new final section performing the live cross-check (component 2) at a
concrete `(Oh, λ_c, n, M, l)`, closing "Next steps" item 3 in the existing summary
table.

### 5. Contract documentation

Add to `CLAUDE.md`: a rheology model may not merge without (a) a from-scratch CAS
derivation notebook with assertions, (b) at least one assertion that cross-checks
against the live running Julia code, (c) the notebook executing clean in CI.

## Testing / acceptance

- `git push` on the branch triggers both CI jobs; both must go green.
- Locally verifiable before pushing: `julia --project=julia -e 'using Pkg; Pkg.test()'`
  and `uv run jupyter nbconvert --to notebook --execute notebooks/*.ipynb`.
- No changes to `julia/src/*` production code in this branch — this is
  derivation/CI/docs only. If the parity assertions in step 3.4 reveal an actual
  mismatch against `ob_extension.jl`, that's a separate finding to report, not
  silently patch.

## Explicitly out of scope

- Any Cross-fluid-model derivation or implementation.
- Refactoring `ob_extension.jl`/`st_extension.jl` production code.
