# DropSolver

A Julia implementation of the linearized spectral solver for viscous and viscoelastic drop impact on a flat substrate. Implements the v3 linearized BDF time-stepper from Agüero-Vera et al., with Newtonian and Oldroyd-B (OB) rheology.

## Physics

A liquid drop of radius $R$ falls under gravity and impacts a flat substrate at $z = 0$. The drop shape is expanded in Legendre modes $A_n(t)$ ($n = 2, \ldots, M$) on top of a spherical base state; the substrate contact is tracked by a discrete contact-point count $c_p \in \{0, 1, \ldots, M\}$.

The governing equations (in dimensionless units scaled by $R$, $\sqrt{\rho R^3 / \sigma}$) are:

- **R1/R2** – BDF-discretized shape equations for $A_n$ and $\dot{A}_n$
- **R3** – contact condition $h(\theta_i) = 0$ at each contact angle
- **R4** – pressure-free condition $\sum_n B_n P_n(\cos\theta_i) = 0$ at free angles
- **R6/R7** – COM kinematics: $\dot{z} = v$ and $\dot{v} = -1/\mathrm{Fr} - B_1$ (v3 linearization)

For Oldroyd-B fluids the viscous damping term is replaced by an effective constitutive law parameterized by the Deborah number $\mathrm{De}_1 = \lambda_1 \sigma_{2;0}$ and solvent fraction $\beta_s = \mu_s / \mu$.

Two dimensionless groups control the Newtonian problem:

| Symbol | Definition | Physical meaning |
|--------|-----------|-----------------|
| $\mathrm{Oh} = \nu\sqrt{\rho/(\sigma R)}$ | Ohnesorge | viscosity / surface tension |
| $\mathrm{Fr} = \sigma / (\rho g R^2)$ | Froude | surface tension / gravity |

## Key design choices

**Gauss-Legendre collocation** (`make_theta_vec(M)`): the $M+1$ angular collocation points are the $M$ Gauss-Legendre nodes of $P_M(\cos\theta)$ plus the south pole $\theta = \pi$. This gives a well-conditioned Legendre-Vandermonde matrix (condition number $\approx 16$ for $M = 20$, vs $10^{16}$ for uniform spacing).

**CFL-based time step** (`make_dt_max(M)`): the maximum stable time step follows the MATLAB reference formula $\Delta t_\max = 2\pi / (\sqrt{M(M+2)(M-1)} \cdot 8)$, which scales as $M^{-3/2}$.

**Jacobian caching**: the system is linear, so $J^{-1}$ is exact and constant for fixed $(M, c_p, \Delta t, \text{order})$. It is cached on first use and reused for all subsequent steps with the same key.

## Requirements

- Julia 1.12+
- No external packages (only `LinearAlgebra` from the standard library)

## Installation

```julia
# from the julia/ directory
julia --project=.
```

Or activate in any Julia session:

```julia
using Pkg; Pkg.activate("/path/to/km-viscous-drop/julia")
using DropSolver
```

## Running the tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

All 125 tests pass, including:
- Legendre polynomial recursion
- BDF coefficient correctness
- Newtonian Lamb oscillation decay rate and frequency (< 5% error vs analytical)
- OB eigenvalue validation (< 5% error vs characteristic equation root)
- MATLAB parity: contact time, max spreading radius, and CoR at M=20 (all < 5% vs MATLAB N=90 reference)

## Showcase scripts

All scripts are run from the `julia/` directory.

### Newtonian free oscillation — Lamb limit

```bash
julia --project=. scripts/run_newtonian.jl
```

Prints a table comparing the solver's decay rate $\gamma$ and frequency $\omega$ against the Lamb (1932) analytical result for $l = 2$ oscillations at several Ohnesorge numbers. Expected output:

```
Oh      γ_lamb   γ_fit    err%     ω_lamb   ω_fit    err%
Oh=0.010  ...     ...      < 5%     ...      ...      < 2%
Oh=0.050  ...     ...      < 5%     ...      ...      < 2%
Oh=0.100  ...     ...      < 5%     ...      ...      < 2%
```

### Oldroyd-B oscillation — effect of elasticity

```bash
julia --project=. scripts/run_ob_case.jl
```

Compares the $l = 2$ decay rate and frequency for Newtonian vs OB at increasing $\mathrm{De}_1$. Demonstrates that viscoelasticity suppresses viscous damping.

### OB eigenvalue sweep — characteristic equation parity

```bash
julia --project=. scripts/run_eigenvalue_sweep.jl
```

Sweeps $(\mathrm{Oh}, \mathrm{De}_1, \beta_s)$ and compares the simulated decay rate and frequency against the root of the OB characteristic equation. All errors < 5%.

### Drop impact: MATLAB parity table

```bash
julia --project=. scripts/run_matlab_parity.jl
```

Runs the canonical impact case ($\mathrm{Oh} = 0.3038$, $\mathrm{Fr} = 53.9$, $v_0 = -0.281$) at $M = 6, 10, 20, 40, 60$ and compares against the MATLAB reference (N=90 modes). Expected output at M=20:

```
M=20   t_c=2.938  err=1.7%   r_max=0.392  err=1.3%   CoR=0.476  err=1.6%
```

### Drop impact: Newtonian vs Oldroyd-B trajectory

```bash
julia --project=. scripts/run_impact.jl
```

Runs both Newtonian and OB drops at $M = 20$ and prints a side-by-side time series of $(t, z, c_p, A_2)$, showing how viscoelasticity modifies the contact and rebound phases.

## API reference

### Configuration

```julia
make_theta_vec(M)   # M+1 Gauss-Legendre collocation angles (descending, south-pole first)
make_dt_max(M)      # CFL-stable maximum time step: 2π / (√(M(M+2)(M-1)) · 8)

precompute_integrals(NaN, M)  # precomputed Legendre integrals (returns (matrix, ...) tuple)
SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
OBParams(De1, beta_s)         # OBParams(0.0, 1.0) = Newtonian
DropState(M)                  # zero-initialized state; set .z, .v, .dt, .cp manually
```

### Running a simulation

```julia
times, states = solve_drop!(cfg, ob, init;
                             t_end      = 10.0,
                             save_every = 0.05,
                             dt_init    = make_dt_max(M),
                             dt_min     = 1e-6)
```

Returns `times::Vector{Float64}` and `states::Vector{DropState}` at each saved frame. The solver uses adaptive BDF1/BDF2 time-stepping with automatic contact detection and dt halving on failure.

### Minimal example

```julia
using DropSolver

M = 20; Oh = 0.3038; Fr = 53.9

dt_max    = make_dt_max(M)
theta_vec = make_theta_vec(M)
precomp   = precompute_integrals(NaN, M)[1]
cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
ob        = OBParams(0.0, 1.0)   # Newtonian

init   = DropState(M)
init.z = 1.1          # drop center height (must be > 1 - deformation for no contact)
init.v = -0.281       # dimensionless impact velocity
init.dt = dt_max
init.cp = 0

times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=8.0, save_every=0.05)

contact = filter(s -> s.cp > 0, states)
println("Contact frames: $(length(contact))")
println("Max cp: $(maximum(s.cp for s in contact))")
```

### State fields

| Field | Type | Description |
|-------|------|-------------|
| `A[n]` | `Vector{Float64}` | Deformation amplitudes $A_n$, $n = 1 \ldots M$ ($A_1 \equiv 0$) |
| `Adot[n]` | `Vector{Float64}` | Time derivatives $\dot{A}_n$ |
| `B[n]` | `Vector{Float64}` | Pressure amplitudes $B_0 \ldots B_M$ (length $M+1$) |
| `z` | `Float64` | Center-of-mass height |
| `v` | `Float64` | Center-of-mass velocity |
| `t` | `Float64` | Current time |
| `cp` | `Int` | Number of contact points |

## Reference

Based on the MATLAB v3 solver described in:

> Agüero-Vera et al., *Spectral simulation of viscous drop impact with Oldroyd-B rheology*, in preparation.

MATLAB reference values (N=90 modes, $\mathrm{Oh} = 0.3038$, $\mathrm{Fr} = 53.9$):
- Contact time: 2.99
- Maximum spreading radius: 0.397
- Coefficient of restitution: 0.484
