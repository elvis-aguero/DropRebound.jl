# API Reference

Everything the `DropSolver` module exports, grouped by the order a run meets it:
build a configuration, integrate a trajectory, then measure the bounce. The
later groups are for readers extending the model: the damping and frequency
coefficients, the non-Newtonian closures, and the residual and Jacobian blocks
those closures modify.

To simply run a simulation you need `SimConstants`, `DropState`, one rheology
parameter struct (`OBParams`, `STParams`, or `STExactParams`), and
[`solve_drop!`](@ref); the Getting started example on the home page uses nothing
else.

Unexported helpers are collected in the Internals section at the foot of the
page. They are implementation detail and may change without notice.

## Setting up a run

`SimConstants` holds everything fixed for a run, meaning mode count, Ohnesorge and Bond
numbers, the collocation grid and the quadrature tables. `DropState` holds
everything that evolves. Both are built once, before the time loop.

```@autodocs
Modules = [DropSolver]
Pages   = ["types.jl", "integrals.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Running a simulation

`solve_drop!` is the entry point: an adaptive BDF1/BDF2 loop that searches for
the contact-point count at each step. `solve_drop_v1!` is the alternative
formulation, which tracks a continuous contact angle `θ*` instead of a discrete
count of contacting collocation points.

```@autodocs
Modules = [DropSolver]
Pages   = ["timestepper.jl", "timestepper_v1.jl", "newton.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Measuring the bounce

Reduce a completed trajectory to the quantities an experiment reports, namely contact
time, coefficient of restitution and contact radius, or reconstruct the drop
outline for plotting.

```@autodocs
Modules = [DropSolver]
Pages   = ["postprocessing.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Damping and frequency coefficients

Per-mode damping ``\lambda_l`` and squared frequency ``\omega_l^2``, either from
Lamb's small-Ohnesorge asymptotics or from the exact roots of Reid's (1960)
characteristic equation. `drop_viscous_coeffs` is what `SimConstants` calls; the
root-finding machinery and its interpolation tables are exposed because the
shear-thinning closures re-evaluate these coefficients at a shear-dependent
effective Ohnesorge number on every step.

```@autodocs
Modules = [DropSolver]
Pages   = ["reid.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Non-Newtonian rheology

Each closure adds a correction on top of the Newtonian residual and Jacobian.
Oldroyd-B carries an extra polymer-stress variable per mode; the two
Carreau-Yasuda routes differ in whether the shear-thinning law is Taylor
expanded (`build_residual_st!`, valid only while `eps_ST` is small) or evaluated
exactly with full multi-mode shear coupling (`build_residual_st_exact!`).

```@autodocs
Modules = [DropSolver]
Pages   = ["ob_extension.jl", "st_extension.jl", "st_exact_extension.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Residual, Jacobian, and contact detection

The nonlinear system solved at every timestep, the state-vector packing it
expects, and the geometric tests that decide where the drop touches the wall.

```@autodocs
Modules = [DropSolver]
Pages   = ["residual.jl", "jacobian.jl", "residual_v1.jl", "contact.jl"]
Public  = true
Private = false
Order   = [:type, :constant, :function]
```

## Internals

Not exported, and not part of the public interface: caches, root-finding
helpers, quadrature bases, and the tuning constants behind the contact and
post-processing heuristics. Documented for readers of the source; do not depend
on them.

```@autodocs
Modules = [DropSolver]
Public  = false
Private = true
Order   = [:type, :constant, :function]
```

## Variational assembly

The model of *Shear-Thinning Drops* as stated: three quadratic forms and the
Euler–Lagrange equations, with the interior retained as part of the state. The
coordinates are the interior *displacement* amplitudes on a Ritz basis, and the
surface amplitude is their boundary trace.

Assembling the coefficient matrices needs only one derivative of the velocity, and
both forms are Hessians of quadratic forms, hence symmetric, which is a
correctness test rather than a remark. With a constant viscosity the assembly
reproduces Reid's exact `λ_l` and `ω_l²`.

```@autodocs
Modules = [DropSolver]
Pages   = ["variational.jl"]
Public  = true
Private = false
Order   = [:type, :function]
```

## Impact and contact

Time integration of the variational model against a solid substrate. The drop's
centre of mass and the interior displacement amplitudes advance together under BDF2,
and the air-film pressure enters as a Legendre series closed by collocation: the gap
vanishes at the contacting nodes and the pressure vanishes at the free ones.

The collocation nodes are `θ = π` together with the zeros of ``P_M``, which cluster
near the poles and so resolve the contact where it forms.

Two closures find the contact extent, and they are interchangeable. [`simulate`](@ref) runs a
primal active set: it grows the contact while a free node has penetrated, releases while the
outermost contacting node's pressure pulls, and stops when neither holds. [`simulate_lcp`](@ref)
assembles the affine map from nodal pressure to gap and solves ``h \ge 0``, ``p \ge 0``,
``p_i h_i = 0`` by pivoting, so no node's status is assumed and the contact set need not be an
interval. See *Choosing a Solver* for what the difference is worth.

```@autodocs
Modules = [DropSolver]
Pages   = ["variational_solve.jl"]
Public  = true
Private = false
Order   = [:type, :function]
```

## Choosing a backend

The three solver paths differ in two independent choices: how the momentum equation
is assembled, and how the contact set is found. [`Backend`](@ref) names a combination
of the two, so a study can be swept over solvers without rewriting the call.

```julia
using DropSolver

for b in (Backend(formulation = :variational,    contact = :active_set),
          Backend(formulation = :variational,    contact = :lcp, forcing = :nodal),
          Backend(formulation = :nonvariational, contact = :tangency))
    r = run_impact(b; We = 0.5, Oh = 0.0373, Bo = 0.0189, M = 45, K = 3)
    println(label(b), "  CoR = ", r.cor)
end
```

[`run_impact`](@ref) returns a named tuple carrying `ok`, which is `false` when the
run terminated early rather than throwing, so a sweep over many parameter sets does
not stop at the first failure. [`drop_outline`](@ref) reconstructs the interface shape
from a stored mode vector, for animation or for eyeballing a suspect run.

```@autodocs
Modules = [DropSolver]
Pages   = ["backends.jl"]
Public  = true
Private = false
Order   = [:type, :function]
```

## Logging

The solver logs through Julia's standard `Logging`, so it is controlled the way
any other Julia package is. At the default level a run announces its parameters
and rheology, the dimensionless time of contact onset and of lift-off, and warns
if the time step approaches its floor.

Per-step diagnostics — accepted steps, step halvings, Jacobian cache hits and
misses — are debug level:

```julia
ENV["JULIA_DEBUG"] = "DropSolver"
```

To silence a run entirely, wrap it:

```julia
using Logging
with_logger(NullLogger()) do
    simulate(p)
end
```

## The bundled scripts

`scripts/` holds everything that uses the solver rather than being part of it.
All of them are run from the repository root with `--project=.`, except the ones
that draw figures, which need the plotting stack in the documentation
environment and so take `--project=docs`.

| script | what it does |
|---|---|
| `run_newtonian.jl` | excites ``l=2`` on a free drop and extracts the decay rate and frequency, against Lamb (1932) at several Ohnesorge numbers |
| `run_ob_case.jl` | the same decay measurement for Newtonian against Oldroyd-B at rising ``\mathrm{De}_1``, showing elasticity suppressing the damping |
| `run_eigenvalue_sweep.jl` | sweeps ``(\mathrm{Oh},\mathrm{De}_1,\beta_s)`` and checks the measured rate and frequency against the exact root of the Oldroyd-B characteristic equation |
| `run_impact.jl` | one Newtonian and one Oldroyd-B impact side by side, printing height, contact-node count and leading amplitude |
| `run_sweep.jl` | a Cartesian product of parameters streamed to CSV, resuming from whatever rows the file already holds |
| `run_animation.jl` | renders an impact to MP4, rasterising in pure Julia and piping raw RGB24 to `ffmpeg` |
| `compare_solvers.jl` | regenerates the table on *Choosing a Solver* |
| `validate_gabbard.jl`, `validate_shear_thinning.jl` | the two validation campaigns of *Home* |
| `callmap.jl`, `callmap_draw.jl` | a static call graph of `src/`, and its rendering |

Every script writes under `outputs/`: numbers to `outputs/csv/`, images to
`outputs/figures/`. A figure carries the name of the script that drew it, so
`figure_vorticity.jl` produces `figure_vorticity.png`.
