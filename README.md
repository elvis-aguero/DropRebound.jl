# DropRebound.jl

Julia solver for the impact and rebound of a liquid drop on a flat substrate, with support for Newtonian, viscoelastic (Oldroyd-B), and shear-thinning (Carreau-Yasuda) rheology.

<table>
<tr>
<td align="center" width="50%">

![Oldroyd-B drop impact — Oh = 0.30, De₁ = 0.5, β_s = 0.5](docs/impact_ob.gif)

*Oldroyd-B (Oh = 0.30, De₁ = 0.5, β_s = 0.5, We = 0.5, M = 90)*

</td>
<td align="center" width="50%">

![Carreau shear-thinning drop impact — Oh = 0.30, ε_ST = 0.15, λ_c = 0.02](docs/impact_st.gif)

*Carreau shear-thinning (Oh = 0.30, ε_ST = 0.08, λ_c = 0.02, We = 0.5, M = 90)*

</td>
</tr>
</table>

## What it does

A spherical drop falls under gravity, deforms on contact with a solid surface, spreads, then rebounds. This package solves the linearized spectral equations governing that process: the drop shape is expanded in Legendre modes around the spherical base state (valid for small deformations, $|A_n| \ll 1$), contact is tracked through a discrete set of collocation points on the drop surface, and the equations are integrated with an adaptive BDF1/BDF2 time-stepper.

The Julia module name is `DropSolver` (`using DropSolver`).

## Dimensionless parameters

All quantities are non-dimensionalized by the drop radius $R$ and the capillary time $\tau_{\rm cap} = \sqrt{\rho R^3 / \sigma}$ (the natural oscillation period of a free drop). Simulation output times are in units of $\tau_{\rm cap}$.

Three numbers control the Newtonian problem:

| Parameter | Definition | Physical meaning |
|-----------|-----------|-----------------|
| Oh (Ohnesorge) | $\nu\sqrt{\rho/(\sigma R)}$ | Ratio of viscous dissipation to surface-tension restoring force. Oh ≪ 1: nearly inviscid, bounces with little energy loss. Oh ~ 1: heavily damped, may not rebound. |
| Bo (Bond) | $\rho g R^2 / \sigma$ | Ratio of gravitational body force to surface tension. Bo ≪ 1: surface tension dominates, drop stays nearly spherical at rest. Bo ~ 1: gravity flattens the drop significantly. |
| We (Weber) | $\rho v_0^2 R / \sigma$ | Ratio of impact kinetic energy to surface tension. We ≪ 1: gentle impact, drop barely deforms. We ~ 1: vigorous spreading. In the solver, We = v₀² where v₀ is the dimensionless impact velocity set via `init.v`. |

The gravity term in the center-of-mass equation is $-\mathrm{Bo}$ (dimensionless), so `Bo = 0` gives a gravity-free oscillation.

As a reference point: a water–glycerol drop of radius 0.2 mm falling at 9 cm/s gives Oh ≈ 0.30, Bo ≈ 0.019, We ≈ 0.079.

### Viscoelastic drops — the Oldroyd-B model

The Oldroyd-B model describes dilute polymer solutions: a Newtonian solvent of viscosity $\mu_s$ mixed with long elastic polymer chains that contribute an additional viscosity $\mu_p = \mu - \mu_s$. During deformation the chains stretch and store elastic energy; when released they recoil and drive extra rebound — polymer drops can bounce higher and ring longer than equivalent Newtonian drops.

The constitutive law is a convected Maxwell model for the polymer stress combined with a Newtonian solvent. It reduces to Newtonian flow when $\lambda_1 \to 0$ or $\mu_p \to 0$.

Two additional parameters:

| Parameter | Definition | Physical meaning |
|-----------|-----------|-----------------|
| De₁ | $\lambda_1 / \tau_{\rm cap}$ | Deborah number: polymer relaxation time $\lambda_1$ scaled by the capillary time. De₁ ≪ 1: polymer relaxes fast — behavior is nearly Newtonian. De₁ ≫ 1: polymer stores elastic energy across many oscillations — elastic effects dominate. |
| β_s | $\mu_s / \mu$ | Solvent fraction. β_s = 1: purely Newtonian solvent ($\mu_p = 0$). β_s → 0: purely polymeric. Typical polymer solutions: β_s ∈ [0.3, 0.9]. |

Setting `De1 = 0` or `beta_s = 1` recovers Newtonian behavior exactly.

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

M  = 20        # Legendre modes — more modes = finer shape resolution
Oh = 0.3038    # Ohnesorge number
Bo = 1/53.9    # Bond number (≈ 0.0186; Bo = ρgR²/σ)

dt_max    = make_dt_max(M)
theta_vec = make_theta_vec(M)
precomp   = precompute_integrals(NaN, M)[1]  # Legendre integral matrix; pass NaN for default grid
cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
ob        = OBParams(0.0, 1.0)               # Newtonian; use OBParams(De1, beta_s) for polymer

init    = DropState(M)
init.z  = 1.1      # drop center height (z = 1 means just touching the substrate)
init.v  = -0.281   # impact velocity (dimensionless, negative = falling)
init.dt = dt_max
init.cp = 0        # no contact initially

times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=8.0, save_every=0.05)

# Inspect contact phase
in_contact = filter(s -> s.cp > 0, states)
println("Contact frames:    ", length(in_contact))
println("Peak deformation:  A₂ = ", maximum(s.A[2] for s in in_contact))
```

Switching to a polymer drop:

```julia
ob_poly = OBParams(0.5, 0.5)   # De₁ = 0.5, β_s = 0.5
times_p, states_p = solve_drop!(cfg, ob_poly, deepcopy(init); t_end=8.0, save_every=0.05)

# Polymer stress amplitude at end of contact
last_contact = last(filter(s -> s.cp > 0, states_p))
println("Peak polymer stress S₂ = ", last_contact.S[2])
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
make_theta_vec(M)
# → M+1 Gauss-Legendre collocation angles (θ values in (π/2, π], south pole first)

make_dt_max(M)
# → CFL-stable maximum time step: 2π / (√(M(M+2)(M-1)) · 8)

precompute_integrals(NaN, M)[1]
# → (M+1)×(M+1) matrix of ∫Pₙ(u)/u³ du integrals used in the contact pressure block.
#   Pass NaN to use the default internal grid. The [1] selects the integral matrix
#   from the returned (matrix, angles) tuple; only the matrix is needed here.

SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
OBParams(De1, beta_s)   # OBParams(0.0, 1.0) = Newtonian
DropState(M)            # zero-initialized state; set .z, .v, .dt, .cp before use
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

| Field | Type | Description |
|-------|------|-------------|
| `A[n]` | `Vector{Float64}` | Shape amplitudes for Legendre modes n = 2…M (`A[1] ≡ 0` by symmetry) |
| `Adot[n]` | `Vector{Float64}` | Time derivatives $\dot{A}_n$ |
| `S[n]` | `Vector{Float64}` | Polymer stress auxiliary variables (Oldroyd-B only; zero for Newtonian) |
| `B[n]` | `Vector{Float64}` | Pressure amplitudes B₀…Bₘ (length M+1) |
| `z` | `Float64` | Center-of-mass height (z = 1: drop just touching substrate) |
| `v` | `Float64` | Center-of-mass velocity (negative = falling) |
| `t` | `Float64` | Current simulation time (in units of $\tau_{\rm cap}$) |
| `dt` | `Float64` | Time step used to reach this state |
| `cp` | `Int` | Contact point count (0 = airborne, > 0 = in contact) |
| `theta_star` | `Float64` | Continuous contact angle θ* ∈ (π/2, π] (π = no contact; used by `solve_drop_v1!`) |

### Alternative solver

`solve_drop_v1!(cfg, init; ...)` implements the same physics with a continuous contact angle θ* instead of discrete contact points. It is provided for comparison and is not recommended for production runs.

## Showcase scripts

Run from the repo root.

### Viscous decay — Newtonian free oscillation

```bash
julia --project=julia julia/scripts/run_newtonian.jl
```

Excites the $l = 2$ shape mode and extracts the decay rate $\gamma$ and oscillation frequency $\omega$ from the time series, comparing against the Lamb (1932) analytical result at several Ohnesorge numbers.

### Viscoelastic oscillation — effect of polymer stress

```bash
julia --project=julia julia/scripts/run_ob_case.jl
```

Compares $\gamma$ and $\omega$ for Newtonian vs Oldroyd-B at increasing De₁. Shows that polymer elasticity suppresses viscous damping — the drop rings longer.

### Eigenvalue validation

```bash
julia --project=julia julia/scripts/run_eigenvalue_sweep.jl
```

Sweeps (Oh, De₁, β_s) and checks the simulated decay rate and frequency against the exact root of the Oldroyd-B characteristic equation. Errors are below 5% across the parameter space.

### Drop impact — Newtonian vs Oldroyd-B

```bash
julia --project=julia julia/scripts/run_impact.jl
```

Runs both fluid types at $M = 20$ and prints a side-by-side time series of center-of-mass height, contact point count, and leading shape amplitude — showing how polymer stress modifies spreading and rebound.

### Parameter sweep

```bash
julia --project=julia julia/scripts/run_sweep.jl [output.csv]
```

Runs a Cartesian product of (Oh, Bo, We, De₁, β_s, M) and streams results to a CSV file. Each row contains the four KPIs: contact time, coefficient of restitution, maximum spreading radius, and maximum $|A_2|$ deformation. Interrupted runs resume automatically — rows already in the CSV are skipped. Edit the grid constants at the top of the script to change the sweep range.

### Animation

```bash
julia --project=julia julia/scripts/run_animation.jl [output.mp4]
```

Renders a drop-impact simulation as an MP4 video. Frames are rasterized in pure Julia (no plotting packages) and piped as raw RGB24 to `ffmpeg`. Requires `ffmpeg` on `PATH`. The 2D output shows the meridional cross-section of the drop with the substrate at $z = 0$.

## Postprocessing API

```julia
# Reconstruct the drop surface as (xs, zs) Cartesian coordinates
xs, zs = drop_profile(state, cfg; n_theta=200)

# Dimensionless contact-patch radius at a single state
r = compute_contact_radius(state, cfg)

# Extract all KPIs from a completed simulation run
kpis = extract_kpis(times, states, cfg)
# → SweepKPIs(contact_time, cor, max_radius, max_A2)
```

`SweepKPIs` fields:

| Field | Description |
|-------|-------------|
| `contact_time` | Duration of the contact phase in capillary times |
| `cor` | Coefficient of restitution $\sqrt{|E_{\rm out}/E_{\rm in}|}$ |
| `max_radius` | Maximum dimensionless spreading radius during contact |
| `max_A2` | Maximum $|A_2|$ shape amplitude during contact |

All fields are `NaN` / `0.0` when no contact occurs.

## Numerical method

The Legendre spectral expansion is exact for the linearized problem. Key implementation choices:

- **Gauss-Legendre collocation**: angular collocation points are the zeros of $P_M$ plus the south pole. This gives a well-conditioned Vandermonde matrix (condition number ~16 at $M = 20$; uniform spacing gives ~$10^{16}$).
- **CFL time step**: $\Delta t_{\max} = 2\pi / (\sqrt{M(M+2)(M-1)} \cdot 8)$, scaling as $M^{-3/2}$.
- **Jacobian caching**: the linearized system has a constant Jacobian for fixed ($M$, contact state, $\Delta t$, BDF order). $J^{-1}$ is computed once and reused — each time step costs one matrix–vector product.

## Tests

```bash
julia --project=julia -e 'using Pkg; Pkg.test()'
```

148 tests cover: Legendre recursion, BDF coefficients, Newtonian free-oscillation decay and frequency (< 5% vs Lamb), Oldroyd-B eigenvalue parity (< 5% vs characteristic equation), drop impact physics (contact detection, Newtonian vs OB divergence during contact), MATLAB parity KPIs, and postprocessing functions (drop profile, contact radius, KPI extraction).

## License

MIT License. See [LICENSE](LICENSE) for details.

## Reference

> Agüero-Vera et al., *Spectral simulation of viscous drop impact with Oldroyd-B rheology*, in preparation.
> *(arXiv link and DOI will be added upon submission.)*
