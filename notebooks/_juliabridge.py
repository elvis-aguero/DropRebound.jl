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
