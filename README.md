# DropRebound.jl

Julia solver for the impact and rebound of a liquid drop on a flat substrate, with support for Newtonian and viscoelastic (Oldroyd-B) rheology.

## What it does

A spherical drop falls under gravity, deforms on contact with a solid surface, spreads, then rebounds. This package solves the linearized spectral equations governing that process: the drop shape is expanded in Legendre modes, contact is tracked through a discrete set of collocation points on the drop surface, and the equations are integrated with an adaptive BDF time-stepper.

## Dimensionless parameters

Two numbers control the Newtonian problem:

| Parameter | Definition | Physical meaning |
|-----------|-----------|-----------------|
| Oh (Ohnesorge) | $\nu\sqrt{\rho/(\sigma R)}$ | Ratio of viscous dissipation to surface-tension restoring force. Oh ≪ 1: nearly inviscid, bounces with little energy loss. Oh ~ 1: heavily damped, may not rebound. |
| Fr (Froude) | $\sigma/(\rho g R^2)$ | Ratio of surface tension to gravitational body force. Large Fr: surface tension dominates, drop stays nearly spherical at rest. Small Fr: gravity flattens it. |

As a reference point: a water–glycerol drop of radius 0.2 mm falling at 9 cm/s gives Oh ≈ 0.30, Fr ≈ 54.

### Viscoelastic drops — the Oldroyd-B model

The Oldroyd-B model describes dilute polymer solutions: a Newtonian solvent mixed with long elastic polymer chains. The chains stretch during impact and store energy, which is released during rebound — polymer drops can bounce higher and ring longer than equivalent Newtonian drops. Two additional parameters:

| Parameter | Definition | Physical meaning |
|-----------|-----------|-----------------|
| De₁ | $\lambda_1 \cdot \omega_{\rm cap}$ | Deborah number: polymer relaxation time $\lambda_1$ scaled by the capillary oscillation frequency. De₁ ≫ 1 means the polymer barely relaxes during one oscillation — elastic behavior dominates. |
| β_s | $\mu_s/\mu$ | Solvent fraction. β_s = 1: purely Newtonian. β_s → 0: purely polymeric. Typical polymer solutions: β_s ∈ [0.3, 0.9]. |

Setting De₁ = 0 or β_s = 1 recovers Newtonian behavior.

## Installation

```julia
# activate the environment from the repo root
julia --project=julia
```

Or from any Julia session:

```julia
using Pkg
Pkg.activate("/path/to/km-viscous-drop/julia")
using DropSolver
```

Requires Julia 1.12+. No external packages — only `LinearAlgebra` and `Logging` from the standard library.

## Quick start

```julia
using DropSolver

M  = 20       # Legendre modes — more modes = finer shape resolution
Oh = 0.3038   # Ohnesorge number (viscosity / surface tension)
Fr = 53.9     # Froude number    (surface tension / gravity)

dt_max    = make_dt_max(M)
theta_vec = make_theta_vec(M)
precomp   = precompute_integrals(NaN, M)[1]
cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
ob        = OBParams(0.0, 1.0)   # Newtonian; use OBParams(De1, beta_s) for polymer

init    = DropState(M)
init.z  = 1.1      # drop center height (z = 1 means just touching the substrate)
init.v  = -0.281   # impact velocity (dimensionless, negative = falling)
init.dt = dt_max
init.cp = 0        # no contact initially

times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=8.0, save_every=0.05)

# Inspect contact phase
in_contact = filter(s -> s.cp > 0, states)
println("Contact frames: ", length(in_contact))
println("Peak deformation A₂: ", maximum(s.A[2] for s in in_contact))
```

## Logging

The solver uses Julia's standard `Logging` library. By default you see:

- `[ Info]  solve_drop! starting` — parameters and rheology at the start of each run
- `[ Info]  Contact onset` / `Lift-off` — with the exact dimensionless time
- `[ Warn]  dt approaching minimum` — if the solver is struggling numerically

Enable per-step diagnostics (accepted steps, dt halvings, Jacobian cache events):

```julia
ENV["JULIA_DEBUG"] = "DropSolver"
```

Suppress all output:

```julia
using Logging
with_logger(NullLogger()) do
    times, states = solve_drop!(cfg, ob, deepcopy(init); ...)
end
```

## API

### Setup

```julia
make_theta_vec(M)              # M+1 Gauss-Legendre collocation angles
make_dt_max(M)                 # CFL-stable maximum time step
precompute_integrals(NaN, M)   # precomputed Legendre integral matrices (returns tuple)
SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
OBParams(De1, beta_s)          # rheology; OBParams(0.0, 1.0) = Newtonian
DropState(M)                   # zero-initialized state; set .z, .v, .dt, .cp before use
```

### Solver

```julia
times, states = solve_drop!(cfg, ob, init;
                             t_end      = 10.0,
                             save_every = 0.05,
                             dt_init    = make_dt_max(M),
                             dt_min     = 1e-6)
```

Returns `times::Vector{Float64}` and `states::Vector{DropState}` at each saved frame. The solver adapts the time step automatically and detects contact transitions. It errors with a clear message if the step size falls below `dt_min`.

### State fields

| Field | Description |
|-------|-------------|
| `A[n]` | Shape amplitudes for Legendre modes n = 2…M (A[1] ≡ 0 by symmetry) |
| `Adot[n]` | Time derivatives Ȧ_n |
| `B[n]` | Pressure amplitudes B₀…Bₘ |
| `z` | Center-of-mass height (z = 1: drop just touching substrate) |
| `v` | Center-of-mass velocity (negative = falling) |
| `t` | Current simulation time |
| `cp` | Contact point count (0 = airborne, > 0 = in contact) |

## Showcase scripts

Run from the repo root.

### Viscous decay — Newtonian free oscillation

```bash
julia --project=julia julia/scripts/run_newtonian.jl
```

Excites the l = 2 shape mode and extracts the decay rate and oscillation frequency from the time series, comparing against the Lamb (1932) analytical result at several Ohnesorge numbers.

### Viscoelastic oscillation — effect of polymer stress

```bash
julia --project=julia julia/scripts/run_ob_case.jl
```

Compares decay rate and frequency for Newtonian vs Oldroyd-B at increasing De₁. Shows that polymer elasticity suppresses viscous damping — the drop rings longer.

### Eigenvalue validation

```bash
julia --project=julia julia/scripts/run_eigenvalue_sweep.jl
```

Sweeps (Oh, De₁, β_s) and checks the simulated decay rate and frequency against the exact root of the Oldroyd-B characteristic equation. Errors are below 5% across the parameter space.

### Drop impact — Newtonian vs Oldroyd-B

```bash
julia --project=julia julia/scripts/run_impact.jl
```

Runs both fluid types at M = 20 and prints a side-by-side time series of center-of-mass height, contact point count, and leading shape amplitude — showing how polymer stress modifies spreading and rebound.

## Numerical method

The Legendre spectral expansion is exact for the linearized problem. Key implementation choices:

- **Gauss-Legendre collocation**: angular collocation points are the zeros of $P_M$ plus the south pole. This gives a well-conditioned Vandermonde matrix (condition number ~16 at M = 20; uniform spacing gives ~10¹⁶).
- **CFL time step**: $\Delta t_{\max} = 2\pi / (\sqrt{M(M+2)(M-1)} \cdot 8)$, scaling as $M^{-3/2}$.
- **Jacobian caching**: the linearized system has a constant Jacobian for fixed (M, contact state, dt, BDF order). $J^{-1}$ is computed once and reused — each time step costs one matrix–vector product.

## Tests

```bash
julia --project=julia -e 'using Pkg; Pkg.test()'
```

125 tests cover: Legendre recursion, BDF coefficients, Newtonian free-oscillation decay and frequency (< 5% vs Lamb), Oldroyd-B eigenvalue parity (< 5% vs characteristic equation), and drop impact physics (contact detection, Newtonian vs OB divergence during contact).

## Reference

Based on the v3 linearized solver described in:

> Agüero-Vera et al., *Spectral simulation of viscous drop impact with Oldroyd-B rheology*, in preparation.
