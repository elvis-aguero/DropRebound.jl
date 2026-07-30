# CAS-before-code contract + CI Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stand up GitHub Actions CI that runs the Julia test suite and executes every
CAS derivation notebook headlessly, add a live Julia cross-check to both the Carreau
notebook and a new, deep Oldroyd-B derivation notebook, and write the CAS-before-code
contract into `CLAUDE.md`.

**Architecture:** Two GitHub Actions jobs (`julia-tests`, `notebooks`). A shared
`notebooks/_juliabridge.py` shells out to `julia --project=julia` to run real
simulations and return `{decay, freq}` as JSON, so notebooks can assert their
symbolic/asymptotic predictions against the actual running solver, not just against
their own algebra. The new Oldroyd-B notebook ports the existing frequency-domain
derivation in `docs/section_oldroydB.tex` into sympy with assertions, then derives
(for the first time anywhere in the repo) why the time-domain Block-S auxiliary
variable in `julia/src/ob_extension.jl` is the correct state-space realization of
that frequency-domain result, and pins down the De₁ convention ambiguity flagged in
the tex's own "Convention note".

**Tech Stack:** sympy, mpmath, Julia 1.12 (`DropSolver`), `uv`, GitHub Actions
(`julia-actions/setup-julia`, `astral-sh/setup-uv`).

**Spec:** `docs/superpowers/specs/2026-07-30-cas-ci-contract-design.md`

---

## Chunk 1: CI workflow + Julia bridge helper

**Files:**
- Create: `.github/workflows/ci.yml`
- Create: `notebooks/_juliabridge.py`
- Test: exercised by pushing the branch (CI itself is the test harness for this chunk)

### Task 1.0: Create the branch

- [ ] **Step 1: Create and check out the feature branch from `main`**

```bash
git status --short   # confirm no unexpected uncommitted changes first
git checkout main
git pull --ff-only origin main
git checkout -b ci/cas-derivation-contract
```

Expected: `Switched to a new branch 'ci/cas-derivation-contract'`. All subsequent
commits in this plan land on this branch.

### Task 1.1: Write the CI workflow

- [ ] **Step 1: Create `.github/workflows/ci.yml`**

```yaml
name: CI

on:
  push:
    branches: [main]
  pull_request:
    branches: [main]

jobs:
  julia-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: '1.12'
      - uses: julia-actions/cache@v2
      - uses: julia-actions/julia-buildpkg@v1
        with:
          project: julia
      - uses: julia-actions/julia-runtest@v1
        with:
          project: julia

  notebooks:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: julia-actions/setup-julia@v2
        with:
          version: '1.12'
      - uses: julia-actions/cache@v2
      - name: Precompile DropSolver
        run: julia --project=julia -e 'using Pkg; Pkg.instantiate(); using DropSolver'
      - uses: astral-sh/setup-uv@v3
      - name: Install notebook dependencies
        run: uv sync
      - name: Execute every derivation notebook
        run: |
          set -e
          for nb in notebooks/*.ipynb; do
            echo "::group::$nb"
            uv run jupyter nbconvert --to notebook --execute --stdout "$nb" > /dev/null
            echo "::endgroup::"
          done
```

- [ ] **Step 2: Verify YAML is well-formed**

Run: `python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"`
Expected: no output (parses cleanly). If `yaml` module is missing, run
`uv run python3 -c "import yaml; yaml.safe_load(open('.github/workflows/ci.yml'))"`
instead (uv's synced env has PyYAML as a jupyter transitive dependency).

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: run Julia test suite and execute CAS notebooks on push/PR"
```

### Task 1.2: Write the shared Julia bridge helper

This is the piece that lets a notebook cell run the *real* solver and compare its
output to the notebook's derived prediction — the part that's been missing from
both existing notebooks.

**Files:**
- Create: `notebooks/_juliabridge.py`

- [ ] **Step 1: Write `notebooks/_juliabridge.py`**

```python
"""Run the real Julia DropSolver from a notebook and extract decay/frequency.

Every CAS derivation notebook should use this to check its symbolic or
asymptotic prediction against the actual running solver — not just against
its own algebra. If this bridge and a notebook's symbolic result disagree,
that's a real bug (in the notebook's derivation, or in the Julia code), not
a rounding issue to explain away.
"""
import json
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
JULIA_PROJECT = REPO_ROOT / "julia"


def run_julia_decay(julia_setup_lines, l, t_periods=6.0, n_save=50, timeout=120):
    """Simulate a single free-oscillating mode and fit its decay rate/frequency.

    `julia_setup_lines` is a Julia source fragment (string) that must define
    `cfg::SimConstants`, `ob::OBParams`, and `st::STParams` (use `STParams()`
    for the Newtonian/OB-only case) in scope, using mode `l`. This function
    supplies the initial condition, run length, and the log-decrement /
    zero-crossing fit — the same method `julia/scripts/run_newtonian.jl` uses,
    so results are comparable to that script's validated output.

    Returns {"decay": float, "freq": float, "n_periods_used": int}.
    """
    script = f"""
    using Pkg; Pkg.activate(raw"{JULIA_PROJECT}")
    using DropSolver

    l = {l}
    {julia_setup_lines}

    omega_guess = sqrt(Float64(l * (l - 1) * (l + 2)))
    T_period = 2*pi / omega_guess

    init = DropState(cfg.M)
    init.A[l] = 0.02
    init.z = 2.0
    init.dt = cfg.dt_max
    init.cp = 0

    times, states = solve_drop!(cfg, ob, init;
        st = st,
        t_end = {t_periods} * T_period,
        save_every = T_period / {n_save},
        dt_init = cfg.dt_max)

    Al = [s.A[l] for s in states]
    decay = -log(abs(Al[end]) / abs(Al[1])) / (times[end] - times[1])

    sign_changes = findall(i -> Al[i] * Al[i+1] < 0, 1:length(Al)-1)
    freq = length(sign_changes) >= 4 ?
        pi / (sum(diff(times[sign_changes])) / (length(sign_changes) - 1)) : NaN

    println("JULIABRIDGE_RESULT:" * "{{\\"decay\\": $(decay), \\"freq\\": $(freq), \\"n_periods_used\\": {t_periods}}}")
    """
    result = subprocess.run(
        ["julia", "-e", script],
        capture_output=True, text=True, timeout=timeout, cwd=REPO_ROOT,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Julia bridge call failed (exit {result.returncode}):\n"
            f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
        )
    for line in result.stdout.splitlines():
        if line.startswith("JULIABRIDGE_RESULT:"):
            return json.loads(line[len("JULIABRIDGE_RESULT:"):])
    raise RuntimeError(
        f"Julia bridge produced no JULIABRIDGE_RESULT line.\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )
```

- [ ] **Step 2: Smoke-test the bridge from the command line (Newtonian case)**

Run:
```bash
uv run python3 -c "
from pathlib import Path
import sys
sys.path.insert(0, 'notebooks')
from _juliabridge import run_julia_decay

setup = '''
M = 6; Oh = 0.05; Bo = 1e-6
theta_vec = make_theta_vec(M)
precomp = precompute_integrals(NaN, M)[1]
dt_max = make_dt_max(M)
cfg = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
ob = OBParams()
st = STParams()
'''
r = run_julia_decay(setup, l=2)
print(r)
expected_decay = (2-1)*(2*2+1)*0.05
print('expected decay ~', expected_decay)
"
```
Expected: prints a dict like `{'decay': 0.248..., 'freq': 2.44..., 'n_periods_used': 6.0}`
and `expected decay ~ 0.25` (`(2-1)*(2*2+1)*0.05 = 0.25`, the Lamb formula
`(l-1)(2l+1)Oh` used identically in `run_newtonian.jl`/`test_ob_eigenvalue.jl`) —
the two `decay` numbers should agree to within a few percent.

- [ ] **Step 3: Commit**

```bash
git add notebooks/_juliabridge.py
git commit -m "feat: add live Julia cross-check bridge for CAS notebooks"
```

---

## Chunk 2: Carreau notebook — add the missing live cross-check

**Files:**
- Modify: `notebooks/shear_thinning_derivation.ipynb` (append cells after cell 50,
  the existing "Summary" cell; do not alter any existing cell except the one
  "Next steps" bullet named below — this chunk only adds new cells otherwise)

The notebook's own summary (cell 50) lists as unfinished: *"2. Validate: run the
Julia solver with small ε_ST and confirm the decay rate correction matches the
analytical prediction from §8."* This task does exactly that. (Verified by
dumping the notebook's actual cells — it has 51 cells total, indices 0-50; the
"Summary" markdown cell is index 50, not 49, and the still-open item to close is
list item 2 ("Validate..."), not item 3 (item 3, "Contact-zone accuracy...",
stays open — do not touch it.))

### Task 2.1: Add the live cross-check section

- [ ] **Step 1: Open the notebook and note current cell count**

Run: `uv run python3 -c "import json; print(len(json.load(open('notebooks/shear_thinning_derivation.ipynb'))['cells']))"`
Expected: `51` (cells 0-50; new cells append starting at index 51).

- [ ] **Step 2: Append a markdown cell introducing the check**

Content (write for a critical reader — explain *why* this specific check is the
one that matters, not just that a check exists):

```markdown
---
## 10. Cross-Check Against the Running Solver

Every assertion so far has checked that this notebook's *own* algebra is
self-consistent (Reid's ODE is satisfied, limits recover Newtonian, signs are
right). None of them have touched `julia/src/st_extension.jl` — the actual code
a user runs. This section closes that gap.

**What we're checking:** the boxed slow-amplitude equation from §8,
$\dot a = -\gamma_l^{(0)} a\,[1 - \tfrac{3}{4}\varepsilon_{ST}\Lambda^2\Gamma_l a^2]$,
predicts a *reduced* effective decay rate at finite amplitude $a$ compared to the
pure-Newtonian rate $\gamma_l^{(0)}$. We run the real solver at a small-but-finite
amplitude with $\varepsilon_{ST} > 0$, extract the decay rate from the time series
exactly as `run_newtonian.jl` does, and compare it to this formula.

**What a failure would mean:** if the live decay rate does *not* fall between the
Newtonian rate and the notebook's corrected prediction (or overshoots it by more
than the tolerance below), either `build_residual_st!`/`build_jacobian_st` has a
sign or coefficient bug, or the derivation above has one — the whole point of this
section is that we can no longer tell ourselves it's "probably fine."
```

- [ ] **Step 3: Append a code cell computing the finite-Oh Γ_l needed for the check**

Reuse the `compute_gamma_l` function already defined in cell 46 (`§6.5`) at the
same `Oh` you'll use for the live run — do not hardcode the inviscid-limit
`GAMMA2_INVISCID` constant from `test_carreau.jl`, since the live run is at finite
Oh and §6.5 exists precisely to handle that case. Example:

```python
import sys
sys.path.insert(0, '.')
from _juliabridge import run_julia_decay

Oh_check = 0.05
l_check = 2
Gamma_l_check = compute_gamma_l(Oh_check, l=l_check)  # from cell 46
print(f"Gamma_{l_check}(Oh={Oh_check}) = {Gamma_l_check}")
```

- [ ] **Step 4: Append a code cell running the live Julia simulation**

```python
lambda_c_check = 0.3
eps_ST_check = 0.3   # n = 1 - 2*eps_ST = 0.4, a fairly strong shear-thinning fluid
a_amplitude = 0.15   # must match the `init.A[l] = ...` value inside _juliabridge

sigma0_check = (l_check * (l_check - 1) * (l_check + 2)) ** 0.5
Lambda_check = lambda_c_check * sigma0_check
gamma0_check = (l_check - 1) * (2 * l_check + 1) * Oh_check

setup = f'''
M = 6; Oh = {Oh_check}; Bo = 1e-6
theta_vec = make_theta_vec(M)
precomp = precompute_integrals(NaN, M)[1]
dt_max = make_dt_max(M)
cfg = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
ob = OBParams()
Gamma = fill({Gamma_l_check}, M-1)
st = STParams({eps_ST_check}, {lambda_c_check}, Gamma)
'''
result = run_julia_decay(setup, l=l_check)
print(result)
```

- [ ] **Step 5: Append the assertion cell**

```python
predicted_gamma_eff = gamma0_check * (1 - 0.75 * eps_ST_check * Lambda_check**2 * Gamma_l_check * a_amplitude**2)

# ASSERTION 20: live solver decay rate matches the slow-amplitude-equation
# prediction within 15% (loose — this is a leading-order asymptotic result at
# finite, not infinitesimal, amplitude, so some mismatch is expected; the point
# is confirming the *direction and rough magnitude* of the effect against the
# real running code, not exact agreement).
rel_err = abs(result["decay"] - predicted_gamma_eff) / gamma0_check
assert rel_err < 0.15, (
    f"Live decay {result['decay']:.5f} vs predicted {predicted_gamma_eff:.5f} "
    f"(Newtonian baseline {gamma0_check:.5f}): {rel_err:.1%} mismatch, expected <15%"
)
# ASSERTION 20b: shear-thinning must reduce the decay rate below Newtonian —
# if this fails while ASSERTION 20 passes, both formula and code agree but
# happen to agree on the wrong sign, which the other asymptotic assertions
# (11-15a) should already have caught; keep this as a second independent check.
assert result["decay"] < gamma0_check, (
    f"Live decay {result['decay']:.5f} should be below Newtonian {gamma0_check:.5f}"
)
print(f"OK: live decay {result['decay']:.5f} vs predicted {predicted_gamma_eff:.5f} "
      f"(Newtonian baseline {gamma0_check:.5f}), {rel_err:.1%} mismatch")
```

- [ ] **Step 6: Execute the whole notebook and confirm it passes**

Run: `uv run jupyter nbconvert --to notebook --execute --stdout notebooks/shear_thinning_derivation.ipynb > /dev/null`
Expected: exits 0. If ASSERTION 20 fails, do not loosen the tolerance to make it
pass — first check by hand whether `eps_ST_check`/`a_amplitude` are actually small
enough for the asymptotic formula to apply (try halving `a_amplitude`); if the
mismatch persists at small amplitude, that indicates a real bug and should be
reported, not silently patched around.

- [ ] **Step 7: Update the "Next steps" list in the existing summary cell (cell 50)**

Change item **2** ("Validate: run the Julia solver with small ε_ST and confirm the
decay rate correction matches the analytical prediction from §8.") from future
tense to a completed-item note pointing at §10, e.g. replace that bullet with:
`2. ~~Validate against the running solver~~ — done in §10.` Leave items 1 and 3
untouched.

- [ ] **Step 8: Save notebook outputs and commit**

```bash
git add notebooks/shear_thinning_derivation.ipynb
git commit -m "feat: add live solver cross-check to Carreau derivation notebook"
```

---

## Chunk 3: New Oldroyd-B derivation notebook

**Files:**
- Create: `notebooks/oldroyd_b_derivation.ipynb`

This notebook must satisfy the dual-purpose writing standard from the spec: every
markdown cell readable by a critical human reader with no other file open, and
every assertion accompanied by a "this checks ___, would fail if ___" comment.

Structure mirrors the Carreau notebook's rigor (imports → derive → assert →
narrate), but the content is new — it ports the existing frequency-domain
derivation in `docs/section_oldroydB.tex` into sympy, then derives (for the first
time in this repo) *why* the time-domain Block-S auxiliary variable in
`ob_extension.jl` correctly realizes that frequency-domain result, which today is
simply asserted in a code comment with no derivation anywhere.

### Task 3.1: Notebook scaffold and constitutive law (§1)

- [ ] **Step 1: Create `notebooks/oldroyd_b_derivation.ipynb`** with a title
  markdown cell and an imports cell (`sympy`, `mpmath`), following the header/
  notation style of cell 0 in `shear_thinning_derivation.ipynb` (state the
  notation table: $x=r/R$, $\sigma$, $q$, $\alpha$, $U(x)$, $b_l$, Oh — plus the
  new symbols $\lambda_1,\lambda_2,\beta_s,De_1,De_2,\mu_{\rm eff}(\sigma)$).

- [ ] **Step 2: §1 — Port the constitutive law and linearisation**

  Reproduce symbolically, with sympy `Eq` objects, the chain in
  `docs/section_oldroydB.tex:62-154` (eq:OB_raw → eq:OB_beta → eq:mueff_OB /
  eq:mueff_both): define $\tau, \lambda_1, \lambda_2, \beta_s, \mu, \sigma$ as
  sympy symbols, write the frequency-domain constitutive relation
  $(1-\lambda_1\sigma)\tau = 2\mu(1-\lambda_2\sigma)e$, and solve symbolically
  for $\mu_{\rm eff}(\sigma) = \mu(1-\lambda_2\sigma)/(1-\lambda_1\sigma)$.

  - **ASSERTION 1**: `sp.simplify` shows `mu_eff.subs(lambda_2, beta_s*lambda_1)`
    equals the $\beta_s$-form in eq:mueff_both. *Checks that the two
    parameterisations the tex carries side-by-side are actually algebraically
    identical — would fail if a transcription error slipped in porting the tex.*
  - **ASSERTION 2**: `mu_eff.subs(lambda_1, lambda_2)` simplifies to exactly `mu`.
    *Checks the Newtonian limit; would fail if a sign were wrong in the ratio.*
  - **ASSERTION 3**: `mu_eff.subs(lambda_2, 0)` matches the Maxwell effective
    viscosity `mu/(1-lambda_1*sigma)`. *Checks the Maxwell limit ($\beta_s=0$).*

  Markdown before the derivation must explain, in plain language, what "effective
  viscosity" means here: that a linear, frequency-independent Newtonian relation
  $\tau=2\mu e$ becomes frequency-*dependent* once the polymer has memory, and
  that this is the entire mechanism by which Oldroyd-B differs from Newtonian at
  linear order — there is no new nonlinearity, only a modified (complex,
  frequency-dependent) viscosity plugged into otherwise-unchanged Newtonian
  equations.

### Task 3.2: q*², characteristic equation, and universal limits (§2)

- [ ] **Step 1: Port eq:lamsig → eq:qstar_OB_final** (`section_oldroydB.tex:198-337,
  461-505`): define $q^2,\alpha^2,De_1,De_2$ as sympy symbols, derive
  $q_*^2 = q^2(\alpha^2-De_1q^2)/(\alpha^2-De_2q^2)$ from $\mu_{\rm eff}$ by
  substitution (not by copying the boxed result — actually perform the sympy
  substitution and simplify to it).

  - **ASSERTION 4**: symbolically simplify $\alpha_*^2/q_*^2$ (eq:alphastar_OB)
    and confirm it equals $\alpha^2/q^2$ — the ratio-cancellation identity
    eq:ratio_cancel. *Checks the key structural fact that makes the OB
    characteristic equation's LHS identical to Reid's; would fail if the
    substitution for $q_*^2$ were wrong.*
  - **ASSERTION 5**: `q_star_sq.subs(De1, De2)` simplifies to `q**2` exactly
    (Newtonian: $q_*=q$). *Checks eq:qstar_OB's stated Newtonian limit.*
  - **ASSERTION 6**: `q_star_sq.subs(De2, 0)` matches the Maxwell $q_*^2$ from
    the "Four clean limits" table (`section_oldroydB.tex:710-740`).
  - **ASSERTION 7**: `q_star_sq.subs(q**2, alpha**2/De1)` simplifies to `0`
    (relaxation zero) and the denominator `alpha**2 - De2*q**2` vanishes at
    `q**2 = alpha**2/De2` (retardation pole) — check both algebraically, not
    just by quoting the tex's notebox.

  Markdown: explain in plain language what a "pole" and a "zero" of $q_*^2$ mean
  physically here (fully relaxed polymer vs. frozen/solid-like polymer stress) —
  a critical reader should understand *why* Oldroyd-B has qualitatively new
  behavior (the retardation pole) that neither Newtonian nor Maxwell has, before
  hitting the algebra that produces it.

- [ ] **Step 2: Reproduce the characteristic equation eq:char_OB numerically**
  using the same Bessel-ratio machinery already in `src/viscous-drop.ipynb`
  (`import_root_gabbard`, `Q`) — adapt it to accept $q_*$ in place of $q$ per
  eq:qstar_OB_final, mirroring how `test_ob_eigenvalue.jl` already
  reimplements this independently in Julia (`bessel_ratio`,
  `_bessel_ratio_miller`).

  - **ASSERTION 8**: at `De1 = De2 = 0` and `Oh = 0` exactly, the root-finder
    recovers the same $(q_a, q_b)$ as `import_root_gabbard(0, l)` in
    `src/viscous-drop.ipynb` for `l = 2, 3, 5` — this is the one regime where
    `import_root_gabbard` is exact, not a continuation seed (its own docstring
    states it is only exact at `Oh = 0`). *Checks the zero-viscosity limit of
    the OB root-finder against the already-validated inviscid Newtonian roots.*
  - **ASSERTION 8b**: at `De1 = De2 = 0` and finite `Oh` (`0.05`, `0.3`), the
    root-finder must instead be checked against a Newtonian reference that is
    actually valid at finite Oh — either (i) `continuation_coef_finder(Oh, l)`
    from `src/viscous-drop.ipynb`, which refines the `Oh=0` seed via Newton
    continuation to the true finite-Oh root, or (ii) an independent finite-Oh
    root-find reusing the Bessel-ratio machinery `test_ob_eigenvalue.jl`
    already implements in Julia (`bessel_ratio`/`_bessel_ratio_miller`), called
    at `De1=0` as a cross-language check. Do **not** compare against
    `import_root_gabbard(Oh, l)` at nonzero `Oh` — it is a continuation seed,
    not a validated finite-Oh value, and a mismatch there would not indicate a
    bug in the new code. *Checks the OB root-finder against a genuinely
    finite-Oh-valid Newtonian reference; would fail if the $q_*$ substitution
    were implemented incorrectly in the numerical solver even though the
    symbolic algebra above passed.*

### Task 3.3: Why Block-S is correct — the missing derivation (§3)

This is the section that doesn't exist anywhere today: `ob_extension.jl:41-45`
simply states the Block-S ODE in a comment. Derive it.

- [ ] **Step 1: Port the memory-kernel result** (`section_oldroydB.tex:655-707`,
  eq:integrodiff_OB–eq:OB_kernel): state the time-domain kernel
  $K_{\rm OB}(s) = \beta_s\delta(s) + \frac{1-\beta_s}{\lambda_1}e^{-s/\lambda_1}$
  whose Laplace-type transform reproduces $\mu_{\rm eff}(\sigma)/\mu$, symbolically
  (compute the transform in sympy and confirm it matches eq:mueff_both after the
  $\beta_s$ substitution).

  - **ASSERTION 9**: the transform of the kernel equals $\mu_{\rm eff}(\sigma)/\mu$
    from §1. **Sign-convention warning**: the tex's kernel transform uses the
    convention $\int_0^\infty K(s)e^{+\sigma s}\,ds = \mu_{\rm eff}(\sigma)/\mu$
    (positive exponent, matching its $e^{-\sigma t}$ *decay* convention where
    $\mathrm{Re}(\sigma)>0$ means decay — see `section_oldroydB.tex:134`, and
    the notebook's own notation table in cell 0's counterpart cell for this
    notebook). A direct `sp.laplace_transform(kernel, s, sigma)` computes the
    standard $e^{-\sigma s}$ transform and will produce a sign-flipped
    denominator ($1+\lambda_1\sigma$ instead of $1-\lambda_1\sigma$) unless you
    either substitute $\sigma \to -\sigma$ afterward or do the integral manually
    with the $+\sigma s$ exponent (handle the $\beta_s\delta(s)$ term via
    `sp.DiracDelta` and note $\int_0^\infty \delta(s)f(s)ds = f(0)/2$ vs.
    $f(0)$ conventions differ across sympy versions — verify which one
    `sp.integrate` uses before trusting the $\beta_s$ term's coefficient).
    Check the sign carefully and confirm against the known correct answer
    (eq:mueff_both) rather than accepting whatever sign sympy returns. *This is
    the crux check: it confirms the time-domain kernel and the frequency-domain
    effective viscosity are the same object, viewed two ways — would fail if
    the kernel weights, timescale, or transform sign were transcribed wrong
    from the tex.*

- [ ] **Step 2: Derive the state-space realization.** Show that defining $S(t)$
  by the ODE $\lambda_1\dot S + S = (1-\beta_s)\lambda_1\,\dot{\text{(strain rate)}}$... concretely,
  in the code's actual dimensionless form (capillary time $\tau$, mode amplitude
  $\dot A_n$): $De_1\,\dot S_n + S_n = (1-\beta_s)\dot A_n$, and that the
  polymer-branch contribution to the damping force is *exactly* $S_n$ while the
  solvent branch contributes $\beta_s\dot A_n$ directly (no auxiliary variable
  needed) — reproduces the same convolution as $K_{\rm OB}$ from Step 1. Do this
  by solving the linear ODE for $S$ via Laplace transform symbolically
  ($\tilde S(\sigma) = \frac{1-\beta_s}{1+De_1\sigma}\dot{\tilde A}(\sigma)$,
  matching sign convention to the tex's $e^{-\sigma t}$ decay convention — check
  this sign carefully, it is a common place to get backwards) and confirming
  $\beta_s + \tilde S(\sigma)/\dot{\tilde A}(\sigma)$ equals $\mu_{\rm eff}(\sigma)/\mu$
  from §1.

  - **ASSERTION 10**: the above equality holds symbolically (sympy `simplify`
    of the difference is `0`). *This is the derivation `ob_extension.jl` is
    missing: it proves the auxiliary variable $S_n$ the code integrates is not
    an ad hoc convenience but the unique minimal state-space realization of the
    polymer memory kernel. A failure here would mean the code's Block-S ODE
    computes something other than what the frequency-domain physics requires.*

  Markdown for this section must be written assuming the reader has never seen
  a state-space realization of a transfer function before: explain in one or
  two sentences why a first-order linear ODE in an auxiliary variable can encode
  an infinite-memory convolution exactly (because the kernel is a single
  exponential), and that this is precisely why $S_n$ needs no history buffer
  in the Julia code — it *is* the compressed history.

### Task 3.4: Resolve the De₁ convention question (§4)

`docs/section_oldroydB.tex:176-196` explicitly flags that its own $De_1 =
\lambda_1\sigma_{l;0}$ convention differs from Zrnić & Brenn's (referenced to the
capillary timescale). Separately, `README.md`'s dimensionless-parameters table
defines $De_1 = \lambda_1/\tau_{\rm cap}$ (capillary time, no $l$-dependence), while
`julia/src/types.jl:14` comments `OBParams.De1` as "λ₁·σ_{l;0}" (mode-dependent).
These are candidate-different conventions for the same field name; this section
must determine, by testing against the code's actual numerical behavior rather
than by re-reading the comments, which one `OBParams.De1` actually implements.

- [ ] **Step 1: State the two candidate readings explicitly** in a markdown cell:
  (a) `De1` is $\lambda_1/\tau_{\rm cap}$, a single fluid property independent of
  mode number $l$ (matches README); (b) `De1` is $\lambda_1\sigma_{l;0}$, and thus
  the *same* `OBParams.De1` value passed to the solver implies a *different*
  physical $\lambda_1$ depending which mode $l$ is dominant (matches
  `types.jl`'s comment and the tex).

- [ ] **Step 2: Read the code first — this is the primary evidence, not the live
  run.** `ob_extension.jl`'s Block-S ODE (`De1 * dS/dτ = (1-beta_s)*Adot - S`) is
  written in dimensionless capillary time $\tau$ with a single scalar `ob.De1`,
  applied identically across all modes $n=2\ldots M$ in `build_residual_ob!`
  (residual at `julia/src/ob_extension.jl:57-70`, Jacobian coupling at
  `julia/src/ob_extension.jl:154-167`) — there is no per-mode $\sigma_{l;0}$
  rescaling of `De1` anywhere in either block. Quote those exact line ranges in
  this cell. By inspection alone this already forces reading (a) (a single
  fluid property in capillary time, mode-independent, matching README's
  `De1 = λ1/τ_cap`): if `De1` meant $\lambda_1\sigma_{l;0}$ per mode, the code
  would need to multiply by a mode-dependent `sigma0[k]` somewhere in Block S,
  and it does not.

  - **ASSERTION 11**: given that reading, `types.jl`'s comment ("λ₁·σ_{l;0}") is
    misleading as written — the code does not implement a per-mode Deborah
    number. State this conclusion from the code-reading evidence above, then
    confirm it numerically as a second, independent check (not the discovery
    mechanism): pick two different modes $l=2$ and $l=5$, run the live solver
    (Chunk 1's bridge) at the *same* `OBParams(De1=0.5, beta_s=0.5)` exciting
    each mode separately, and confirm the extracted decay rates are each
    consistent with reading (a)'s characteristic equation (§2) using the
    *same* `De1=0.5` for both — reading (b) would require using
    `De1 * sigma0(2)` for one and `De1 * sigma0(5)` for the other to match,
    which should visibly fail if attempted.
  - Report the outcome honestly: if the numeric check confirms the code-reading
    conclusion, say so and recommend fixing the `types.jl` comment in a
    follow-up (do **not** fix it in this notebook-only branch — out of scope,
    see plan header). If the numeric check instead contradicts the code-reading
    conclusion, that means something in this task's own reasoning is wrong —
    stop and re-examine before writing the summary, don't report a contradiction
    as the final answer.

### Task 3.5: Code-parity assertions against `ob_extension.jl` (§5)

- [ ] **Step 1: Reduce the BDF-discretized Block-S coefficients symbolically**
  from the ODE in Task 3.3 Step 2 (`De1*Ṡ + S = (1-beta_s)*Adot`), using the
  same BDF1 discretization convention as `julia/src/bdf.jl` (reuse or reproduce
  its coefficient formulas — read that file first), and confirm the resulting
  Jacobian entries match `julia/src/ob_extension.jl` literally:
  - `(ak + dt/De1)` for $\partial R_S/\partial S$ — code at `ob_extension.jl:163`
  - `-dt*(1-beta_s)/De1` for $\partial R_S/\partial\dot A$ — code at
    `ob_extension.jl:162`
  - `D2_diag = 2*Oh*(n-1)*(2n+1)` for the coupling of $S$ into $R_2$ — code at
    `ob_extension.jl:156`, and confirm this is the *same* $D_2$ used in the
    pure-Newtonian $R_2$ (`residual.jl:75`), i.e. that Oldroyd-B does not
    introduce a second, independent damping coefficient — it redistributes the
    existing one between $\beta_s\dot A$ and $S$.

  - **ASSERTION 12, 13, 14**: one per bullet above — symbolic expression from
    the notebook's own BDF reduction equals the literal Julia expression
    (compare as strings/sympy `Eq` after substituting the same symbol names, or
    numerically evaluate both at several `(dt, De1, beta_s, Oh, n)` tuples and
    assert equality to machine precision). *Each checks one specific Jacobian
    entry in the shipped code against an independent symbolic derivation —
    would fail if `ob_extension.jl` had a transcription bug that none of the
    existing finite-difference Jacobian tests (`test_ob.jl`) happened to
    exercise (finite-difference checks can pass even when a term is right by
    the wrong reasoning, e.g. compensating sign errors elsewhere; a
    from-scratch symbolic parity check does not have that failure mode).*

### Task 3.6: Live solver cross-check and Newtonian limits (§6)

- [ ] **Step 1**: using `notebooks/_juliabridge.py` (Chunk 1), run the real OB
  solver at three concrete `(Oh, De1, beta_s, l)` points —
  `(0.05, 0.5, 0.5, 2)`, `(0.05, 1.0, 0.3, 2)`, `(0.3, 0.5, 0.5, 3)` — and compare
  the extracted decay rate/frequency to the roots of eq:char_OB found
  numerically in Task 3.2 Step 2 (same `(Oh, De1, beta_s, l)` in both). Follow
  the exact same pattern as Chunk 2 Task 2.1 Steps 3-5 (compute prediction → run
  bridge → assert within tolerance). `test_ob_eigenvalue.jl` already validates
  this *exact* comparison — real time-domain BDF simulation vs. root-found
  eigenvalue, same method this step uses — at 5%. Use that same 5% tolerance
  here; do not loosen it without first checking (by running at higher `M` or
  more periods) whether the mismatch is a genuine discrepancy or just an
  under-resolved run.
  - **ASSERTION 15**: as above, one per `(Oh, De1, beta_s, l)` point.
- [ ] **Step 2**: assert `De1 → 0` and `beta_s → 1` both reduce the live solver's
  decay/frequency to the pure-Newtonian values from `run_newtonian.jl`'s method,
  mirroring Carreau's `ε_ST=0` limit checks (ASSERTIONS 2, 3, 11 in that notebook).
  - **ASSERTION 16, 17**.

### Task 3.7: Summary table and execution

- [ ] **Step 1**: add a final markdown cell with a numbered assertion-status
  table (mirroring `shear_thinning_derivation.ipynb` cell 50, the "Summary" cell),
  plus an explicit
  "What we found" paragraph reporting the Task 3.4 De₁-convention outcome
  regardless of which way it went.

- [ ] **Step 2: Execute the full notebook**

Run: `uv run jupyter nbconvert --to notebook --execute --stdout notebooks/oldroyd_b_derivation.ipynb > /dev/null`
Expected: exits 0, all assertions pass. Do not mark any task above complete by
weakening an assertion's tolerance or deleting a hard check to get a green run —
if something doesn't hold, that's a finding for the summary, not a reason to
soften the check.

- [ ] **Step 3: Commit**

```bash
git add notebooks/oldroyd_b_derivation.ipynb
git commit -m "feat: add from-scratch Oldroyd-B CAS derivation with code-parity checks"
```

---

## Chunk 4: Contract documentation and PR

**Files:**
- Modify: `CLAUDE.md`

### Task 4.1: Write the contract into CLAUDE.md

- [ ] **Step 1: Append to `CLAUDE.md`**

```markdown

## Rheology model contract

A rheology model (Newtonian variants, Oldroyd-B, Carreau, or any future model)
may not be merged without:
1. A from-scratch CAS derivation notebook under `notebooks/` with sympy/mpmath
   assertions — not just a reference to a LaTeX writeup.
2. At least one assertion that cross-checks the derivation against the actual
   running Julia solver (see `notebooks/_juliabridge.py`), not only against the
   notebook's own algebra.
3. Every markdown cell readable by a critical reader who has not opened the
   Julia source — state what's being derived and why before the math, restate
   boxed results in plain English, and explain what a failing assertion would
   mean physically.
4. The notebook executing clean in CI (`.github/workflows/ci.yml`, `notebooks`
   job) before the model's Julia implementation is considered done.
```

- [ ] **Step 2: Commit**

```bash
git add CLAUDE.md
git commit -m "docs: write the CAS-before-code rheology contract into CLAUDE.md"
```

### Task 4.2: Push and open PR

- [ ] **Step 1**: `git push -u origin ci/cas-derivation-contract`
- [ ] **Step 2**: confirm both CI jobs go green on the pushed branch before
  opening a PR (check via `gh run list --branch ci/cas-derivation-contract` or
  the Actions tab).
- [ ] **Step 3**: open the PR with `gh pr create`, summarizing the CI setup, the
  Carreau live cross-check, the new OB notebook (call out the Task 3.4 De₁
  finding explicitly in the PR description, whichever way it went), and that
  Cross-model work remains out of scope pending this landing.
