# Contact-selector investigation — archived, deliberately unmerged

This branch is not part of the project. It exists so that a day of diagnostic work is
recoverable rather than lost to a temporary directory, and it is kept off `main` because
none of it belongs in the shipped solver.

**Nothing here is a fix.** Ten hypotheses were tested; two were confirmed as diagnoses,
eight were falsified, and no working replacement for the contact selector came out of it.
The value is in what was ruled out and in the ground truth that was obtained.

## What was established

**The reference implementation's rule includes a `±2` probe that this repository's solver
does not have.** Reading `get_next_step_v5.m` and `solve_motion_v2.m` (around line 288):
three candidates `cp−1, cp, cp+1` are ranked by tangency residual and the contact count
moves by at most one, but a *fourth* evaluation at `cp±2` is a timestep test. If `cp+2`
scores better than `cp+1`, the contact wants to grow faster than one node per step, which
is not licence to jump — it means `dt` cannot resolve the growth, so the step is refused
and `dt` halved. That couples contact growth to time refinement.

Adding the probe makes this solver reproduce the reference **exactly** at `K = 1, M = 90`:

| case | ours | reference |
|---|---|---|
| `We=1, Oh=0.304` | CoR 0.3129, `t_c` 2.1830 | 0.3138, 2.183 |
| `Oh=0.30, We=5` | `t_c` 1.8632, `|ζ|` 0.669 | 1.863, 0.668 |
| `Oh=0.30, We=10` | `t_c` 1.7587, `|ζ|` 0.837 | 1.759, 0.836 |

**Tangency selection is degenerate once the interior is resolved (`K > 1`).** The residual
`|gap(cp+1) − gap(cp)|` is a first difference at one node pair. With one radial trial
function per mode the surface cannot zero it without genuinely flattening; with two, a
one-node contact can leave the substrate smoothly and satisfy it spuriously. Measured
directly (`h17.jl`): at `K = 1` the tangency-optimal contact extent leads the incumbent by
one node and the contact follows it down, tracking in 104 of 120 steps. At `K = 2` the two
agree for six steps, then the optimum collapses to `cp ≈ 1` and never recovers — 9 of 64.

**The energy audit is validated and is the sharpest available regression guard.** With the
film's work accounted (`Σ_l Q_l ζ̇_l + 𝔉 v`), the identity
`ΔE = ∫P_film − ∫ξ̇ᵀCξ̇` closes to 2.4% on a working case. On a failing one the drop
manufactures roughly forty times the energy it arrived with, against an energy ceiling
computed independently from the surface stiffness.

## What was falsified

| hypothesis | verdict |
|---|---|
| high-`We` failures are the linearisation breaking down | **false** — the reference completes the same cases with `|ζ| = 1.2`, inside the energy bound |
| the failure is over-determination of the onset system | false — square, residual `1e-10` |
| feasibility scope (upper-hemisphere nodes) causes `dt` collapse | false — `mine`/`ancestor` verdicts agree at every candidate |
| the `errortan` scale mismatch at `cp = 0` is the cause | false — residuals are all the same order |
| `cp = 0` is an absorbing state | false — neither run returns to it after contact |
| the KKT solve diverges at `K = 2` | false — every candidate `:ok`, `cond(A) = 3e7` |
| an exact spectral `dh/dθ` fixes the degeneracy | false — 10/64 versus 9/64; a sharper measure of a degenerate quantity is still degenerate |
| an edge-pressure gate fixes it | false — no gain, and it regresses `Oh=0.30, We=10, K=1` from working to failing |

## Method notes worth keeping

Four traps recurred and each cost real time:

* **Downsampled export read as simulation.** The reference's `saving_frequency` is coarser
  than its timestep, so an exported contact sequence of `0, 4` is `0,1,2,3,4` compressed.
  Inferring that the reference jumps at onset was wrong, and it had already caught me once
  earlier in the same session.
* **Instrumenting only accepted steps.** A branch counter that records nothing on rejected
  steps reported `RJ-shrink = 0` for a failure mode that was 45 of 47 rejections.
* **Near-zero denominators.** Relative errors normalised by a quantity that is itself
  vanishing reported `4e11` where the numerator was machine epsilon.
* **Generalising from one state.** The edge-pressure gate came from a single measurement at
  `t = 0`, where forcing `cp ≥ 2` necessarily means pulling the drop down. Of course the
  pressure was negative there; it does not follow that it is elsewhere.

## Contents

`julia/` — hypothesis tests, in rough order: `falsify.jl`, `f3.jl` (energy audit and its
control), `h11.jl`–`h19.jl` (the selector chain), `predict.jl` (the falsified prediction).
`analyse_probe.jl`, `onset.jl`, `free_compare.jl` read the reference implementation's output.

`matlab/` — probes driving the reference implementation directly. These require
`harrislab-brown/LowWeberDropRebound` cloned alongside, and one patch to it: its
`theta_vector` uses `vpasolve`, which needs the Symbolic Math Toolbox; the same nodes come
from Golub–Welsch on the Legendre Jacobi matrix. Its `initial_height` assertion also
compares a cgs value against a dimensionless one, so pass a large height to bypass it.
Contact times in its exported data are in **seconds** and need dividing by
`sqrt(ρR³/σ)`.

## Where it left off

The recommendation is not another selector. The discrete contact problem is a linear
complementarity problem, and it should be solved as one rather than by ranking candidates
with a hand-built functional. See the note on `main`.
