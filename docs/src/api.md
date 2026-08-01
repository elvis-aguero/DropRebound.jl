# API Reference

Everything the `DropSolver` module exports, grouped by the order a run meets it:
build a configuration, integrate a trajectory, then measure the bounce. The
later groups are for readers extending the model — the damping and frequency
coefficients, the non-Newtonian closures, and the residual and Jacobian blocks
those closures modify.

To simply run a simulation you need `SimConstants`, `DropState`, one rheology
parameter struct (`OBParams`, `STParams`, or `STExactParams`), and
[`solve_drop!`](@ref); the Getting started example on the home page uses nothing
else.

Unexported helpers are collected in the Internals section at the foot of the
page. They are implementation detail and may change without notice.

## Setting up a run

`SimConstants` holds everything fixed for a run — mode count, Ohnesorge and Bond
numbers, the collocation grid, the quadrature tables — and `DropState` holds
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

Reduce a completed trajectory to the quantities an experiment reports — contact
time, coefficient of restitution, contact radius — or reconstruct the drop
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
