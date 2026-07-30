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


def find_ob_eigenvalue_only(Oh, De1, beta_s, l, timeout=60):
    """Return the Oldroyd-B characteristic-equation root only (no simulation).

    Reuses `find_ob_eigenvalue` from `julia/test/test_ob_eigenvalue.jl` (see
    `run_ob_eigenvalue_check`'s docstring for why we call into that file rather
    than reimplementing it). Unlike `run_ob_eigenvalue_check`/`run_ob_sim`,
    this works correctly for any `l` — `run_ob_sim` itself hardcodes
    `init.A[2]` regardless of its own `l` keyword, so it silently excites mode
    2 while comparing against mode `l`'s eigenvalue whenever `l != 2`. That is
    a real latent bug in that test helper (never exercised because every test
    in that file calls `run_ob_sim` with the default `l=2`), noted here rather
    than fixed, since fixing test code is out of scope for this notebook.

    Returns {"gamma_exact": float, "omega_exact": float}.
    """
    test_file = JULIA_PROJECT / "test" / "test_ob_eigenvalue.jl"
    script = f"""
    using Pkg; Pkg.activate(raw"{JULIA_PROJECT}")
    using DropSolver

    src = read(raw"{test_file}", String)
    idx = findfirst("@testset", src)
    Base.include_string(Main, src[1:idx[1]-1])

    sigma_exact = find_ob_eigenvalue({Oh}, {De1}, {beta_s}, {l})
    println("JULIABRIDGE_RESULT:" * "{{\\"gamma_exact\\": $(real(sigma_exact)), \\"omega_exact\\": $(imag(sigma_exact))}}")
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


def run_ob_eigenvalue_check(Oh, De1, beta_s, l=2, timeout=120):
    """Cross-check the Oldroyd-B char. equation root against a live simulation.

    Reuses `find_ob_eigenvalue`/`run_ob_sim`/`extract_decay_freq` directly from
    `julia/test/test_ob_eigenvalue.jl` (290 lines of dedicated, already-CI-green
    tests covering the Bessel-ratio root-finder) instead of reimplementing a
    second, independent Bessel root-finder in Python. This is a deliberate
    choice: it reuses already-audited code rather than risking a second,
    redundant implementation with its own bugs. Only the top of that file (the
    function definitions, before the first `@testset`) is evaluated, so no
    tests run here — just the helper functions become callable.

    Returns {"gamma_exact": float, "omega_exact": float,
             "gamma_sim": float, "omega_sim": float}, where "_exact" is the
    characteristic-equation root and "_sim" is extracted from a real
    `solve_drop!` time series.
    """
    test_file = JULIA_PROJECT / "test" / "test_ob_eigenvalue.jl"
    script = f"""
    using Pkg; Pkg.activate(raw"{JULIA_PROJECT}")
    using DropSolver

    src = read(raw"{test_file}", String)
    idx = findfirst("@testset", src)
    Base.include_string(Main, src[1:idx[1]-1])

    times, A2, sigma_exact = run_ob_sim({Oh}, {De1}, {beta_s}; l={l})
    gamma_sim, omega_sim = extract_decay_freq(times, A2)

    println("JULIABRIDGE_RESULT:" * "{{\\"gamma_exact\\": $(real(sigma_exact)), \\"omega_exact\\": $(imag(sigma_exact)), \\"gamma_sim\\": $(gamma_sim), \\"omega_sim\\": $(omega_sim)}}")
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
