# Julia Oldroyd-B Drop Impact Solver — Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Julia equivalent of the MATLAB `km-dropplet-solidsubstrate-v3` simulation that adds Oldroyd-B viscoelasticity via the internal-variable (auxiliary ODE) method, reproduces the Newtonian limit exactly, and validates against Zrnić & Brenn (2024) modal frequencies.

**Architecture:** The solver uses the same Legendre spectral decomposition (M harmonics, M+1 angular samples), adaptive BDF1/BDF2 time-stepping, and Newton-Raphson at each step as the MATLAB code. Oldroyd-B memory is handled by an auxiliary ODE per mode: `De₁·Ṡₙ = (1−β_s)·Ȧₙ − Sₙ`, which turns the integro-differential equation into a local ODE system, adds M−1 new unknowns (S₂…S_M), and trivially recovers the Newtonian case at De₁=0 or β_s=1. All state variables are dimensionless capillary-inertial units (length ∼ R₀, time ∼ √(ρR₀³/σ)).

**Tech Stack:** Julia 1.12, LinearAlgebra (stdlib), Test (stdlib), no external dependencies by design (matches MATLAB's self-contained approach).

---

## Equation Reference

In capillary-inertial units, the coupled ODE system for mode n (n = 2…M) is:

```
Block 1 (kinematics):     Ȧₙ  =  dAₙ/dτ
Block 2 (OB motion):      Äₙ  +  2·Oh·(n−1)(2n+1)·(β_s·Ȧₙ + Sₙ)
                               +  n(n+2)(n−1)·Aₙ
                               +  cp·n·Bₙ  =  0       [cp=1 if contact, 0 otherwise]
Block S (polymer ODE):    De₁·Ṡₙ  +  Sₙ  =  (1−β_s)·Ȧₙ
Block 3 (contact):        cos(θ_c)·(1 + ΣPₙ(θ_c)·Aₙ) + z_COM = 0   [at each contact angle]
Block 4 (pressure BC):    ΣBₙ·Pₙ(cos(θ)) = 0                        [at each free angle]
Block 6 (COM kinematics): dz_COM/dτ = v_COM
Block 7 (COM dynamics):   dv_COM/dτ = −1/Fr − cp·B₁
```

Newtonian limit: De₁=0 or β_s=1 → Block S forces Sₙ=0 instantly → Block 2 reduces to `2·Oh·(n−1)(2n+1)·Ȧₙ`, reproducing the MATLAB Newtonian equation exactly.

State vector (OB, length 4M):
```
X = [A₂…A_M,  Ȧ₂…Ȧ_M,  S₂…S_M,  B₀…B_M,  z_COM,  v_COM]
     (M−1)    (M−1)    (M−1)    (M+1)     1       1
```
Newtonian state vector (length 3M+1): same without the S block.

BDF1 coefficients: `coefs = [-1, 1]`  
BDF2 coefficients (rk = dt/dt_prev): `ak=(1+2rk)/(1+rk), bk=−(1+rk), ck=rk²/(1+rk)`, `coefs=[ck, bk, ak]`

---

## File Map

```
km-viscous-drop/julia/
├── Project.toml
├── src/
│   ├── DropSolver.jl          # module root; re-exports public symbols
│   ├── types.jl               # DropState, SolverConfig, OBParams, SimConstants
│   ├── legendre.jl            # collect_Pl, collect_dPl (recurrence-based, vectorized)
│   ├── integrals.jl           # precompute_integrals: ∫_{cos a}^{cos b} Pₙ(u)/u³ du
│   ├── bdf.jl                 # bdf_coefficients(order, dt, dt_prev) → Vector
│   ├── residual.jl            # build_residual!(R, X, hist, dt, cp, cfg) — Newtonian 7-block
│   ├── jacobian.jl            # build_jacobian(X, hist, dt, cp, cfg) → Matrix  — Newtonian
│   ├── ob_extension.jl        # extend for OB: adds Block S rows/cols to R and J
│   ├── newton.jl              # newton_solve!(X, R!, J; tol, maxiter) + Jacobian cache
│   ├── contact.jl             # contact_error(state, theta_vec, cp) → Float64
│   └── timestepper.jl         # solve_drop!(cfg; ...) → (times, states)
└── test/
    ├── runtests.jl
    ├── test_legendre.jl       # recurrence, orthogonality, derivatives
    ├── test_integrals.jl      # numerical integrals vs Wolfram Alpha reference values
    ├── test_bdf.jl            # coefficient formulas; BDF2 convergence order on a toy ODE
    ├── test_residual.jl       # each block individually, expected zeros
    ├── test_contact.jl        # contact_error = 0 at a known tangent point
    ├── test_newtonian.jl      # free l=2 oscillation: ω, γ vs Lamb formula
    └── test_ob.jl             # β_s=1 → same as Newtonian; De₁ small → same Lamb limit
```

---

## Chunk 1: Foundation

### Task 1: Project scaffold

**Files:**
- Create: `km-viscous-drop/julia/Project.toml`
- Create: `km-viscous-drop/julia/src/DropSolver.jl`
- Create: `km-viscous-drop/julia/test/runtests.jl`

- [ ] **Step 1: Create the directory structure**

```bash
cd /Users/eaguerov/Documents/Github/km-viscous-drop
mkdir -p julia/src julia/test julia/scripts
```

- [ ] **Step 2: Write Project.toml**

```toml
name = "DropSolver"
uuid = "a1b2c3d4-e5f6-7890-abcd-ef1234567890"
version = "0.1.0"

[deps]

[compat]
julia = "1.12"

[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[targets]
test = ["Test"]
```

- [ ] **Step 3: Write the module root `src/DropSolver.jl`**

```julia
module DropSolver

include("types.jl")
include("legendre.jl")
include("integrals.jl")
include("bdf.jl")
include("residual.jl")
include("jacobian.jl")
include("ob_extension.jl")
include("newton.jl")
include("contact.jl")
include("timestepper.jl")

export DropState, SolverConfig, OBParams, SimConstants
export collect_Pl, collect_dPl
export precompute_integrals
export bdf_coefficients
export build_residual!, build_jacobian
export extend_residual_ob!, extend_jacobian_ob!
export newton_solve!
export contact_error
export solve_drop!

end
```

- [ ] **Step 4: Write minimal `test/runtests.jl`**

```julia
using Test

@testset "DropSolver" begin
    include("test_legendre.jl")
    include("test_integrals.jl")
    include("test_bdf.jl")
    include("test_residual.jl")
    include("test_contact.jl")
    include("test_newtonian.jl")
    include("test_ob.jl")
end
```

- [ ] **Step 5: Verify project loads**

```bash
cd /Users/eaguerov/Documents/Github/km-viscous-drop/julia
julia --project=. -e 'using DropSolver; println("ok")'
```
Expected: `ok`

- [ ] **Step 6: Commit**

```bash
git add julia/
git commit -m "feat: scaffold Julia DropSolver project"
```

---

### Task 2: Legendre polynomial utilities

**Files:**
- Create: `km-viscous-drop/julia/src/legendre.jl`
- Create: `km-viscous-drop/julia/test/test_legendre.jl`

The MATLAB `collectPl(N, x)` returns a matrix where column n contains Pₙ(xᵢ) for each input point xᵢ. `collectdnPl(N, x)` returns first derivatives P'ₙ(xᵢ). Both use the three-term recurrence.

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_legendre.jl
using Test
using DropSolver: collect_Pl, collect_dPl

@testset "Legendre polynomials" begin
    @testset "P₀, P₁, P₂ at specific points" begin
        x = [0.0, 0.5, 1.0]
        P = collect_Pl(4, x)   # returns matrix (length(x) × N+1), 1-indexed: P[:,n+1] = Pₙ
        # P₀ = 1 everywhere
        @test P[:, 1] ≈ ones(3)
        # P₁ = x
        @test P[:, 2] ≈ x
        # P₂(x) = (3x²-1)/2
        @test P[:, 3] ≈ @. (3x^2 - 1) / 2
        # P₂(1) = 1
        @test P[3, 3] ≈ 1.0
    end

    @testset "Recurrence consistency" begin
        x = collect(range(-1, 1, 21))
        P = collect_Pl(10, x)
        # (n+1)Pₙ₊₁ = (2n+1)x Pₙ - n Pₙ₋₁
        for n in 1:9
            @test (n+1) .* P[:, n+2] ≈ (2n+1) .* x .* P[:, n+1] - n .* P[:, n]
        end
    end

    @testset "Orthogonality ∫₋₁¹ Pₙ Pₘ dx = 2/(2n+1) δₙₘ" begin
        N = 6
        x = collect(range(-1.0, 1.0, 200))
        dx = x[2] - x[1]
        P = collect_Pl(N, x)
        for n in 0:N-1, m in 0:N-1
            integral = sum(P[:, n+1] .* P[:, m+1]) * dx
            expected = n == m ? 2.0/(2n+1) : 0.0
            @test abs(integral - expected) < 0.02   # coarse quadrature tolerance
        end
    end

    @testset "Derivatives P'ₙ at specific points" begin
        x = [0.0, 0.5]
        dP = collect_dPl(4, x)
        # P'₁ = 1
        @test dP[:, 2] ≈ ones(2)
        # P'₂(x) = 3x
        @test dP[:, 3] ≈ 3 .* x
        # P'₃(x) = (15x²-3)/2
        @test dP[:, 4] ≈ @. (15x^2 - 3) / 2
    end

    @testset "Derivative recurrence (1-x²)P'ₙ = n(Pₙ₋₁ - xPₙ)" begin
        x = collect(range(-0.9, 0.9, 15))
        P  = collect_Pl(6, x)
        dP = collect_dPl(6, x)
        for n in 1:5
            lhs = @. (1 - x^2) * dP[:, n+1]
            rhs = n .* (P[:, n] .- x .* P[:, n+1])
            @test lhs ≈ rhs  atol=1e-12
        end
    end
end
```

- [ ] **Step 2: Run to verify failures**

```bash
cd /Users/eaguerov/Documents/Github/km-viscous-drop/julia
julia --project=. test/test_legendre.jl
```
Expected: `UndefVarError: collect_Pl not defined`

- [ ] **Step 3: Implement `src/legendre.jl`**

```julia
# src/legendre.jl

"""
    collect_Pl(N, x) → Matrix{Float64}  size = (length(x), N+1)

Evaluate Legendre polynomials P₀…Pₙ at each point in x using the
three-term recurrence. Column index is 1-based: column k holds P_{k-1}.
Mirrors MATLAB collectPl(N, x).
"""
function collect_Pl(N::Int, x::AbstractVector{<:Real})
    n = length(x)
    P = zeros(n, N + 1)
    P[:, 1] .= 1.0          # P₀ = 1
    if N >= 1
        P[:, 2] .= x        # P₁ = x
    end
    for k in 2:N
        # (k)·P_k = (2k-1)·x·P_{k-1} - (k-1)·P_{k-2}
        @. P[:, k+1] = ((2k - 1) * x * P[:, k] - (k - 1) * P[:, k-1]) / k
    end
    return P
end

"""
    collect_dPl(N, x) → Matrix{Float64}  size = (length(x), N+1)

Evaluate first derivatives P'₀…P'ₙ via the identity
    (1-x²)P'ₙ = n(P_{n-1} - x·Pₙ)
Column index is 1-based: column k holds P'_{k-1}.
"""
function collect_dPl(N::Int, x::AbstractVector{<:Real})
    P  = collect_Pl(N, x)
    dP = zeros(size(P))
    # P'₀ = 0, P'₁ = 1 (handled by formula for n≥1 with special care at |x|=1)
    for n in 1:N
        denom = 1.0 .- x .^ 2
        safe  = abs.(denom) .> 1e-14
        dP[safe,  n+1] .= n .* (P[safe, n] .- x[safe] .* P[safe, n+1]) ./ denom[safe]
        # At x = ±1: P'ₙ(±1) = ±n(n+1)/2
        dP[.!safe, n+1] .= sign.(x[.!safe]) .^ (n+1) .* n .* (n+1) ./ 2
    end
    return dP
end
```

- [ ] **Step 4: Run tests**

```bash
julia --project=. test/test_legendre.jl
```
Expected: all pass

- [ ] **Step 5: Commit**

```bash
git add julia/src/legendre.jl julia/test/test_legendre.jl
git commit -m "feat: Legendre polynomial utilities (collect_Pl, collect_dPl)"
```

---

### Task 3: Types

**Files:**
- Create: `km-viscous-drop/julia/src/types.jl`

No tests needed — structs are thin; they get exercised through every downstream test.

- [ ] **Step 1: Write `src/types.jl`**

```julia
# src/types.jl

"""Physical and numerical parameters for one simulation run."""
struct SimConstants
    M         :: Int       # harmonics_qtt (number of Legendre modes, A₁…A_M)
    N_angles  :: Int       # angular_sampling = M+1
    Oh        :: Float64   # Ohnesorge number
    Fr        :: Float64   # Froude number
    theta_vec :: Vector{Float64}  # M+1 angles ∈ (π/2, π], Gauss-Legendre quadrature roots
    precomp_I :: Matrix{Float64}  # precomputed_integrals output (M+1 × M+1)
    dt_max    :: Float64          # maximum allowed time step
end

"""Oldroyd-B parameters. Set De1=0 or beta_s=1 for Newtonian."""
struct OBParams
    De1    :: Float64   # λ₁·σ_{l;0}  (relaxation Deborah number)
    beta_s :: Float64   # λ₂/λ₁ = μ_s/μ  (solvent fraction, ∈ [0,1])
end
OBParams() = OBParams(0.0, 1.0)   # Newtonian default

"""State at a single time step."""
mutable struct DropState
    A    :: Vector{Float64}   # deformation amplitudes A₁…A_M (A₁ always 0)
    Adot :: Vector{Float64}   # dA/dτ  (same indexing)
    S    :: Vector{Float64}   # polymer stress aux vars S₁…S_M (S₁ always 0; OB only)
    B    :: Vector{Float64}   # pressure amplitudes B₀…B_M (length M+1)
    z    :: Float64           # center-of-mass height
    v    :: Float64           # center-of-mass velocity
    t    :: Float64           # current simulation time
    dt   :: Float64           # step size used to reach this state
    cp   :: Int               # contact_points
end

function DropState(M::Int)
    DropState(zeros(M), zeros(M), zeros(M), zeros(M+1),
              0.0, 0.0, 0.0, 0.0, 0)
end
```

- [ ] **Step 2: Verify it loads**

```bash
julia --project=. -e 'using DropSolver; s=DropState(5); println(s.A)'
```
Expected: `[0.0, 0.0, 0.0, 0.0, 0.0]`

- [ ] **Step 3: Commit**

```bash
git add julia/src/types.jl
git commit -m "feat: DropState, SimConstants, OBParams types"
```

---

### Task 4: Precomputed integrals

**Files:**
- Create: `km-viscous-drop/julia/src/integrals.jl`
- Create: `km-viscous-drop/julia/test/test_integrals.jl`

Mirrors MATLAB `precompute_integrals`. Computes `∫_{cos(a_{i+1})}^{cos(aᵢ)} Pₙ(u)/u³ du` for each interval between consecutive angles and each mode n=0…N. Uses simple Gauss-Legendre quadrature internally via `sum`+`range` (no external package).

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_integrals.jl
using Test
using DropSolver: precompute_integrals

@testset "Precomputed integrals" begin
    # Reference values from Wolfram Alpha:
    # ∫_{-1}^{cos(3π/4)} P₂(u)/u³ du
    # cos(3π/4) = -√2/2 ≈ -0.7071
    # Wolfram: ≈ -0.9375 (from the MATLAB docstring example)
    angles = [π, 3π/4]
    N = 3
    M_mat, ang_out = precompute_integrals(angles, N)
    # M_mat is (1 × 4): modes 0,1,2,3 for the single interval [π, 3π/4]
    @test size(M_mat, 1) == 1
    @test size(M_mat, 2) == N + 1
    # P₂ integral should be close to the Wolfram value
    @test abs(M_mat[1, 3] - (-0.9375)) < 0.01

    @testset "P₀ integral = ∫du/u³ analytically" begin
        # ∫_{a}^{b} u⁻³ du = -1/(2b²) + 1/(2a²)
        a = cos(π)         # = -1
        b = cos(3π/4)      # = -√2/2
        expected_P0 = -1/(2*b^2) + 1/(2*a^2)
        @test abs(M_mat[1, 1] - expected_P0) < 1e-4
    end

    @testset "Matrix shape with NaN input" begin
        M2, _ = precompute_integrals(NaN, 5)
        @test size(M2, 2) == 6   # modes 0..5
    end
end
```

- [ ] **Step 2: Run to verify failures**

```bash
julia --project=. test/test_integrals.jl
```
Expected: `UndefVarError: precompute_integrals not defined`

- [ ] **Step 3: Implement `src/integrals.jl`**

```julia
# src/integrals.jl

"""
    precompute_integrals(angles, N) → (M_mat, angles_out)

Compute the integral ∫_{cos(a_{i})}^{cos(a_{i+1})} Pₙ(u)/u³ du for each
consecutive pair of angles and each n = 0…N.

`angles`: vector of angles in (π/2, π], or NaN for N+1 uniform samples.
Returns:
- `M_mat`:     (n_intervals × N+1) matrix  [mirrors MATLAB precompute_integrals]
- `angles_out`: the processed angle vector
"""
function precompute_integrals(angles_in, N::Int)
    if isnan(angles_in) || (angles_in isa Number && isnan(angles_in))
        angles = range(π, 0; length = N + 1) |> collect
    elseif angles_in isa Number
        angles = range(π, 0; length = Int(angles_in)) |> collect
    else
        angles = collect(Float64, angles_in)
    end

    # Keep only angles in (π·0.51, π]
    angles = filter(θ -> θ > π * 0.51 && θ <= π, angles)
    # Replace with midpoints between consecutive pairs, prepend π
    if length(angles) > 1
        mids = (angles[1:end-1] .+ angles[2:end]) ./ 2
        angles = vcat(π, mids)
    end

    n_intervals = length(angles) - 1
    M_mat = zeros(n_intervals, N + 1)

    for i in 1:n_intervals
        a = cos(angles[i])
        b = cos(angles[i + 1])
        # Gaussian quadrature over [b, a] (b < a since cos is decreasing)
        n_quad = 50
        t, w = gauss_legendre_nodes(n_quad, b, a)
        Pt = collect_Pl(N, t)   # (n_quad × N+1)
        integrand = Pt ./ (t .^ 3)  # element-wise
        M_mat[i, :] = integrand' * w
    end

    return M_mat, angles
end

"""Simple Gauss-Legendre quadrature nodes and weights on [a, b] via Golub-Welch."""
function gauss_legendre_nodes(n::Int, a::Float64, b::Float64)
    # Build symmetric tridiagonal matrix
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    T = diagm(1 => β, -1 => β)
    vals, vecs = eigen(Symmetric(T))
    # Nodes and weights on [-1, 1]
    x = vals
    w = 2 .* vecs[1, :] .^ 2
    # Transform to [a, b]
    xab = @. (b - a) / 2 * x + (b + a) / 2
    wab = @. (b - a) / 2 * w
    return xab, wab
end
```

Note: `eigen`, `diagm`, `Symmetric` come from `LinearAlgebra` (stdlib). Add `using LinearAlgebra` inside the module.

- [ ] **Step 4: Add `using LinearAlgebra` to `DropSolver.jl` before the includes**

- [ ] **Step 5: Run tests**

```bash
julia --project=. test/test_integrals.jl
```
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add julia/src/integrals.jl julia/test/test_integrals.jl
git commit -m "feat: precomputed Legendre-over-u3 integrals with Gauss-Legendre quadrature"
```

---

### Task 5: BDF coefficient utility

**Files:**
- Create: `km-viscous-drop/julia/src/bdf.jl`
- Create: `km-viscous-drop/julia/test/test_bdf.jl`

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_bdf.jl
using Test
using DropSolver: bdf_coefficients

@testset "BDF coefficients" begin
    @testset "BDF1" begin
        c = bdf_coefficients(1, 0.1, NaN)
        @test length(c) == 2
        @test c ≈ [-1.0, 1.0]
    end

    @testset "BDF2 uniform step" begin
        # rk = dt/dt_prev = 1 → ak=3/2, bk=-2, ck=1/2
        c = bdf_coefficients(2, 0.1, 0.1)
        @test length(c) == 3
        @test c ≈ [0.5, -2.0, 1.5]
    end

    @testset "BDF2 formula matches MATLAB" begin
        dt = 0.05; dt_prev = 0.1
        rk = dt / dt_prev
        ak = (1 + 2rk)/(1 + rk)
        bk = -(1 + rk)
        ck = rk^2/(1 + rk)
        c = bdf_coefficients(2, dt, dt_prev)
        @test c ≈ [ck, bk, ak]
    end

    @testset "BDF2 convergence on dy/dt = -y, y(0)=1" begin
        # Exact: y(t) = exp(-t). Verify BDF2 is 2nd order.
        errors = Float64[]
        for N in [20, 40, 80]
            dt = 1.0 / N
            y = [1.0, exp(-dt)]            # exact first two steps
            for k in 2:N-1
                c = bdf_coefficients(2, dt, dt)
                # BDF2: c[1]*y[k-1] + c[2]*y[k] + c[3]*y_{k+1} = dt*(-y_{k+1})
                # (c[3] + dt)*y_{k+1} = -(c[1]*y[k-1] + c[2]*y[k])
                rhs = -(c[1]*y[end-1] + c[2]*y[end])
                push!(y, rhs / (c[3] + dt))
            end
            push!(errors, abs(y[end] - exp(-1.0)))
        end
        # error should roughly quarter each time N doubles → 2nd order
        @test errors[2]/errors[1] < 0.35
        @test errors[3]/errors[2] < 0.35
    end
end
```

- [ ] **Step 2: Run to verify failures**

```bash
julia --project=. test/test_bdf.jl
```

- [ ] **Step 3: Implement `src/bdf.jl`**

```julia
# src/bdf.jl

"""
    bdf_coefficients(order, dt, dt_prev) → Vector{Float64}

Return BDF coefficients `c` such that the BDF approximation is
    c[end]*y^k + c[end-1]*y^{k-1} + ... = dt * f(y^k)

order=1: c = [-1, 1]
order=2: c = [ck, bk, ak] with rk = dt/dt_prev
"""
function bdf_coefficients(order::Int, dt::Float64, dt_prev::Float64)
    if order == 1
        return [-1.0, 1.0]
    elseif order == 2
        rk = dt / dt_prev
        ak = (1 + 2rk) / (1 + rk)
        bk = -(1 + rk)
        ck = rk^2 / (1 + rk)
        return [ck, bk, ak]
    else
        error("Only BDF orders 1 and 2 are supported.")
    end
end
```

- [ ] **Step 4: Run tests — all pass**

- [ ] **Step 5: Commit**

```bash
git add julia/src/bdf.jl julia/test/test_bdf.jl
git commit -m "feat: BDF1/BDF2 coefficient utility with convergence test"
```

---

## Chunk 2: Solver Core

### Task 6: Residual function (Newtonian, no contact)

**Files:**
- Create: `km-viscous-drop/julia/src/residual.jl`
- Create: `km-viscous-drop/julia/test/test_residual.jl`

This is the core numerical kernel. Build one block at a time, testing each.

State-vector packing convention (Newtonian, length 3M+1):
```
idx(A)   = 1        : M-1      (A₂…A_M, 1-based offset; A₁ always 0)
idx(Adot)= M        : 2M-2
idx(B)   = 2M-1     : 3M-1     (B₀…B_M, length M+1)
idx(z)   = 3M
idx(v)   = 3M+1
```

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_residual.jl
using Test, LinearAlgebra
using DropSolver

function make_equilibrium_state(M)
    # Spherical drop, no motion, no contact — should give zero residual
    s = DropState(M)
    s.z = 1.0         # COM one radius above substrate → no contact
    s.v = 0.0
    s.cp = 0
    return s
end

@testset "Residual blocks (Newtonian, no contact)" begin
    M = 4
    Oh = 0.1; Fr = 100.0
    theta_vec, precomp = let
        x = range(π, 0; length = M+1) |> collect
        x, precompute_integrals(NaN, M)[1]
    end
    cfg = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, 0.01)
    ob  = OBParams()     # Newtonian

    @testset "Equilibrium → zero residual (BDF1, cp=0)" begin
        s0  = make_equilibrium_state(M)
        s0.dt = 0.01
        R   = zeros(3M + 1)
        build_residual!(R, s0, [s0], 0.01, 0, cfg, ob)
        # R1 block: Aₙ - Aₙ_prev - dt*Ȧₙ = 0 - 0 - 0 ✓
        @test norm(R[1:M-1]) < 1e-14
        # R2 block: Ȧₙ - Ȧₙ_prev + dt*(capillary)*Aₙ = 0 ✓
        @test norm(R[M:2M-2]) < 1e-14
        # pressure block: zero pressure + zero outside → zero
        @test norm(R[2M-1:3M-1]) < 1e-14
        # COM blocks
        @test abs(R[3M])   < 1e-14
        @test abs(R[3M+1]) < 1e-14
    end

    @testset "Nonzero A₂ gives nonzero R2 (capillary restoring force)" begin
        s0   = make_equilibrium_state(M)
        s0.dt = 0.01
        s1   = make_equilibrium_state(M)
        s1.A[2] = 0.1  # A₂ = 0.1 (mode index 2 = column 2)
        s1.dt = 0.01
        R = zeros(3M + 1)
        build_residual!(R, s1, [s0], 0.01, 0, cfg, ob)
        # R2 for n=2: dt * n(n+2)(n-1) * A₂ = 0.01 * 2*4*1 * 0.1 = 0.008
        @test abs(R[M]) - 0.008 < 1e-10
    end
end
```

- [ ] **Step 2: Run to verify failures**

```bash
julia --project=. test/test_residual.jl
```

- [ ] **Step 3: Implement `src/residual.jl`**

```julia
# src/residual.jl

"""
    pack_X(state, M) → Vector{Float64}

Pack a DropState into the Newton-iteration state vector (Newtonian, length 3M+1).
"""
function pack_X(s::DropState, M::Int)
    vcat(s.A[2:end], s.Adot[2:end], s.B, s.z, s.v)
end

"""
    unpack_X!(state, X, M)

Unpack the Newton-iteration vector X back into a DropState in-place.
"""
function unpack_X!(s::DropState, X::AbstractVector, M::Int)
    s.A[2:end]    .= X[1:M-1]
    s.Adot[2:end] .= X[M:2M-2]
    s.B           .= X[2M-1:3M-1]
    s.z            = X[3M]
    s.v            = X[3M+1]
end

"""
    build_residual!(R, state, history, dt, cp, cfg, ob)

Fill the residual vector R in-place.

- `state`:   proposed new DropState (current iterate)
- `history`: Vector of DropState (most recent last); length 1 → BDF1, length 2 → BDF2
- `dt`:      current time step
- `cp`:      contact_points (integer ≥ 0)
- `cfg`:     SimConstants
- `ob`:      OBParams (ignored if ob.De1==0 && ob.beta_s==1)
"""
function build_residual!(R::AbstractVector, state::DropState,
                         history::Vector{DropState}, dt::Float64,
                         cp::Int, cfg::SimConstants, ob::OBParams)
    M  = cfg.M
    Oh = cfg.Oh
    Fr = cfg.Fr
    θv = cfg.theta_vec
    order = length(history)

    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)

    # ── Unpack current state ──────────────────────────────────────────────────
    A    = state.A[2:end]        # length M-1
    Adot = state.Adot[2:end]
    B    = state.B               # length M+1 (B₀…B_M)
    z    = state.z
    v    = state.v
    ns   = 2:M                   # mode indices

    # History matrices (M-1 × order)
    prev_A    = hcat([h.A[2:end]    for h in history]...)
    prev_Adot = hcat([h.Adot[2:end] for h in history]...)
    prev_z    = [h.z for h in history]
    prev_v    = [h.v for h in history]

    # ── Block R1: dA/dτ = Ȧ ──────────────────────────────────────────────────
    # c[end]*A^k + (lower history) = dt * Ȧ^k
    R[1:M-1] .= dropdims(sum(c[1:order] .* prev_A + c[end] .* A', dims=2), dims=2) .- dt .* Adot

    # Actually use the standard BDF form:
    # sum_i c[i] * A^{k-order+i} = dt * f → rearranged as residual = 0
    R[1:M-1] .= c[end] .* A .+ sum(c[j] .* prev_A[:, j] for j in 1:order) .- dt .* Adot

    # ── Block R2: ODE for Ȧ (Newtonian damping + capillary + pressure) ───────
    D1 = @. ns * (ns + 2) * (ns - 1)          # capillary coefficient
    D2 = @. 2Oh * (ns - 1) * (2ns + 1)        # viscous damping

    # For OB: replace Adot with β_s*Adot + S in the damping term
    if ob.De1 > 0.0 && ob.beta_s < 1.0
        S = state.S[2:end]
        effective_damp = D2 .* (ob.beta_s .* Adot .+ S)
    else
        effective_damp = D2 .* Adot
    end

    pressure_term = cp > 0 ? ns .* B[3:end] : zeros(M-1)

    R[M:2M-2] .= c[end] .* Adot .+ sum(c[j] .* prev_Adot[:, j] for j in 1:order) .+
                 dt .* (D1 .* A .+ effective_damp .+ pressure_term)

    # ── Block R3+R4: contact + pressure BCs ─────────────────────────────────
    # (handled jointly: N_angles = M+1 equations for M+1 pressure unknowns)
    θ_contact = θv[1:cp]
    θ_free    = θv[cp+1:end]

    # R3: cos(θ_c)*(1 + ΣPₙ(θ_c)*Aₙ) + z = 0  for each contact angle
    if cp > 0
        Pc = collect_Pl(M, cos.(θ_contact))   # (cp × M+1)
        Pc_A = Pc[:, 2:end]                   # discard P₀ column (A₁=0)
        R3 = cos.(θ_contact) .* (1.0 .+ Pc_A * A) .+ z
        R[2M-1 : 2M-1+cp-1] .= R3
    end

    # R4: ΣBₙ·Pₙ(cos(θ)) = 0 for each free angle
    n_free = length(θ_free)
    if n_free > 0
        Pf = collect_Pl(M, cos.(θ_free))      # (n_free × M+1)
        Pf_all = hcat(ones(n_free), Pf)       # include P₀ column (B₀ term)
        R4 = Pf_all * B
        R[2M-1+cp : 2M-1+cp+n_free-1] .= R4
    end

    # ── Block R6: dz/dτ = v ─────────────────────────────────────────────────
    R[3M] = c[end] * z + sum(c[j] * prev_z[j] for j in 1:order) - dt * v

    # ── Block R7: dv/dτ = -1/Fr - cp*B₁ ────────────────────────────────────
    R[3M+1] = c[end] * v + sum(c[j] * prev_v[j] for j in 1:order) -
              dt * (-1.0/Fr - (cp > 0 ? B[2] : 0.0))
end
```

- [ ] **Step 4: Fix indexing if tests fail** (the Pf_all construction above needs careful index check)

- [ ] **Step 5: Run tests — all pass**

- [ ] **Step 6: Commit**

```bash
git add julia/src/residual.jl julia/test/test_residual.jl
git commit -m "feat: Newtonian residual function (7-block BDF, no-contact + contact)"
```

---

### Task 7: Analytical Jacobian

**Files:**
- Create: `km-viscous-drop/julia/src/jacobian.jl`

The Jacobian is the block-structured matrix `∂R/∂X`. Mirror the MATLAB `JacobianCalculator_v3.m` block-by-block, then verify `J ≈ FiniteDiff(R)`.

- [ ] **Step 1: Add finite-difference Jacobian to test file** (for verification only, not production)

Add to `test/test_residual.jl`:
```julia
@testset "Jacobian matches finite-difference approximation" begin
    M = 4; Oh = 0.1; Fr = 100.0
    # ... (build cfg, ob, state, history as above)
    J   = build_jacobian(s1, [s0], 0.01, 0, cfg, ob)
    eps = 1e-6
    X0  = pack_X(s1, M)
    R0  = zeros(3M+1); build_residual!(R0, s1, [s0], 0.01, 0, cfg, ob)
    J_fd = zeros(3M+1, 3M+1)
    for j in 1:3M+1
        Xp = copy(X0); Xp[j] += eps
        sp = DropState(M); unpack_X!(sp, Xp, M); sp.cp=0
        Rp = zeros(3M+1); build_residual!(Rp, sp, [s0], 0.01, 0, cfg, ob)
        J_fd[:, j] = (Rp .- R0) ./ eps
    end
    @test norm(J .- J_fd) / norm(J_fd) < 1e-4
end
```

- [ ] **Step 2: Implement `src/jacobian.jl`**

```julia
# src/jacobian.jl

"""
    build_jacobian(state, history, dt, cp, cfg, ob) → Matrix{Float64}

Analytical Jacobian ∂R/∂X for the Newtonian (and OB) residual.
Returns a (3M+1 × 3M+1) matrix for Newtonian, (4M × 4M) for OB.
"""
function build_jacobian(state::DropState, history::Vector{DropState},
                        dt::Float64, cp::Int, cfg::SimConstants, ob::OBParams)
    M  = cfg.M
    Oh = cfg.Oh
    θv = cfg.theta_vec
    Fr = cfg.Fr
    order = length(history)
    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak = c[end]   # coefficient on the current time step

    ns = 2:M
    D1 = diagm(@. ns * (ns + 2) * (ns - 1))
    D2 = diagm(@. 2Oh * (ns - 1) * (2ns + 1))

    Nm = M - 1       # number of free deformation modes
    Nb = M + 1       # number of pressure unknowns
    sz = 3M + 1      # Newtonian state vector length

    J = zeros(sz, sz)

    # ── R1 / dA and dAdot ────────────────────────────────────────────────────
    J[1:Nm, 1:Nm]     .= ak .* I(Nm)         # ∂R1/∂A
    J[1:Nm, Nm+1:2Nm] .= -dt .* I(Nm)        # ∂R1/∂Ȧ

    # ── R2 / dA, dAdot, dB ───────────────────────────────────────────────────
    J[Nm+1:2Nm, 1:Nm]     .= dt .* D1         # ∂R2/∂A
    # viscous damping: ∂R2/∂Ȧ = ak + dt*D2*(β_s for OB, 1 for Newton)
    damp_factor = ob.De1 > 0 && ob.beta_s < 1 ? ob.beta_s : 1.0
    J[Nm+1:2Nm, Nm+1:2Nm] .= (ak + dt * damp_factor) .* I(Nm) .+ dt .* (D2 .- (ak + dt)*I(Nm))
    # Actually: ∂R2/∂Ȧ = ak*I + dt*(D2_diag * beta_s)
    J[Nm+1:2Nm, Nm+1:2Nm] .= diagm(fill(ak, Nm)) .+ dt .* (damp_factor .* D2)
    # pressure: only B₂…B_M enter R2 (B₀ and B₁ do not drive mode motion, per v3)
    if cp > 0
        J[Nm+1:2Nm, 2Nm+3:3Nm+1] .= dt .* diagm(collect(Float64, ns))   # ∂R2/∂B[3:end]
    end

    # ── R3+R4: contact + pressure BCs ────────────────────────────────────────
    θ_contact = θv[1:cp]
    θ_free    = θv[cp+1:end]
    n_free    = length(θ_free)

    if cp > 0
        Pc   = collect_Pl(M, cos.(θ_contact))
        PcA  = Pc[:, 2:end]
        # ∂R3/∂A = cos(θ_c) * Pₙ(cos(θ_c))
        J[2Nm+1:2Nm+cp, 1:Nm] .= cos.(θ_contact) .* PcA
        # ∂R3/∂z = 1 for each contact row
        J[2Nm+1:2Nm+cp, 3M] .= 1.0
    end

    if n_free > 0
        Pf    = collect_Pl(M, cos.(θ_free))
        Pfall = hcat(ones(n_free), Pf)
        J[2Nm+cp+1:2Nm+cp+n_free, 2Nm+1:3Nm+1] .= Pfall
    end

    # ── R6 / dz, dv ─────────────────────────────────────────────────────────
    J[3M, 3M]   = ak
    J[3M, 3M+1] = -dt

    # ── R7 / dv, dB₁ ────────────────────────────────────────────────────────
    J[3M+1, 3M+1] = ak
    if cp > 0
        J[3M+1, 2Nm+2] = -dt    # ∂R7/∂B₁ = -dt (B₁ = X[2Nm+2])
    end

    return J
end
```

- [ ] **Step 3: Run the finite-difference test**

```bash
julia --project=. test/test_residual.jl
```
Expected: all pass including FD check.

- [ ] **Step 4: Fix any index mismatches** — iterate until FD test passes.

- [ ] **Step 5: Commit**

```bash
git add julia/src/jacobian.jl julia/test/test_residual.jl
git commit -m "feat: analytical Jacobian, verified against finite differences"
```

---

### Task 8: Newton-Raphson solver with Jacobian cache

**Files:**
- Create: `km-viscous-drop/julia/src/newton.jl`

Mirror `get_next_step_v5.m`: Newton iteration with Jacobian reuse across steps with the same (cp, dt, order) key.

- [ ] **Step 1: Implement `src/newton.jl`**

```julia
# src/newton.jl

"""Cache mapping (cp, dt, order) → precomputed J⁻¹ for the linearized model."""
const _jac_cache = Dict{Tuple{Int,Float64,Int}, Matrix{Float64}}()

"""
    newton_solve!(X, R!, J_fn; tol=1e-10, maxiter=100) → (converged, iters)

Newton-Raphson: Xₙ₊₁ = Xₙ − J(Xₙ)⁻¹ R(Xₙ).
Uses the Jacobian cache when available (same (cp, dt, order) key).
"""
function newton_solve!(X::Vector{Float64},
                       R!::Function,    # R!(buf, X) fills buf in-place
                       J_fn::Function;  # J_fn(X) → Matrix
                       cache_key::Union{Nothing,Tuple} = nothing,
                       tol::Float64 = 1e-10,
                       maxiter::Int = 100)
    buf = similar(X)
    R!(buf, X)
    best_val = norm(buf)
    best_X   = copy(X)

    for iter in 1:maxiter
        # Compute or retrieve Jacobian inverse
        if cache_key !== nothing && haskey(_jac_cache, cache_key)
            Jinv = _jac_cache[cache_key]
            δX   = Jinv * buf
        else
            J = J_fn(X)
            if cache_key !== nothing
                Jinv = inv(J)
                _jac_cache[cache_key] = Jinv
                δX = Jinv * buf
            else
                δX = J \ buf
            end
        end

        X .= X .- δX
        R!(buf, X)
        val = norm(buf)
        if val < best_val
            best_val = val
            best_X  .= X
        end
        if val < tol || norm(δX) < 1e-12
            break
        end
    end

    X .= best_X
    return best_val < tol
end

"""Clear the Jacobian cache (call when contact pattern or step size changes)."""
clear_jac_cache!() = empty!(_jac_cache)
```

- [ ] **Step 2: Add a smoke test in `test/test_residual.jl`**

```julia
@testset "Newton converges on equilibrium perturbation" begin
    # Start slightly off equilibrium (A₂ = 0.01), Newton should return to equilibrium
    # in a handful of iterations.
    # ... set up R! and J_fn closures using build_residual! and build_jacobian
    # Verify converged == true and norm(R) < tol
end
```

- [ ] **Step 3: Commit**

```bash
git add julia/src/newton.jl
git commit -m "feat: Newton-Raphson solver with Jacobian caching"
```

---

### Task 9: Contact detection

**Files:**
- Create: `km-viscous-drop/julia/src/contact.jl`
- Create: `km-viscous-drop/julia/test/test_contact.jl`

Mirror the MATLAB contact-error metric: tangent slope at the contact boundary, ∞ if penetration.

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_contact.jl
using Test
using DropSolver: contact_error, DropState

@testset "Contact detection" begin
    M = 6
    s = DropState(M)
    θv = collect(range(π, π/2; length=M+1))

    @testset "No contact, drop above substrate" begin
        s.z = 1.5
        @test contact_error(s, θv, 0) < 1e-10   # cp=0: should be fine
    end

    @testset "Penetration → infinite error" begin
        s.z = -2.0    # drop below substrate
        @test contact_error(s, θv, 1) == Inf
    end
end
```

- [ ] **Step 2: Implement `src/contact.jl`**

```julia
# src/contact.jl

"""
    drop_height(state, θ, M) → Float64

Surface height z(θ) = cos(θ)·(1 + Σ Pₙ(cos θ)·Aₙ) + z_COM.
"""
function drop_height(state::DropState, θ::Float64)
    M = length(state.A)
    P = collect_Pl(M - 1, [cos(θ)])   # (1 × M)
    Σ = sum(P[1, n+1] * state.A[n+1] for n in 1:M-1)
    return cos(θ) * (1.0 + Σ) + state.z
end

"""
    contact_error(state, theta_vec, cp) → Float64

Tangent-slope error at the contact boundary (∞ if penetration detected).
Mirror of MATLAB's errortan computation.
"""
function contact_error(state::DropState, theta_vec::Vector{Float64}, cp::Int)
    if cp <= 0
        return 0.0
    end

    # Check penetration below contact angle
    check_angles = theta_vec[theta_vec .> π/2]
    if cp > 0
        check_angles = filter(θ -> θ < theta_vec[cp], check_angles)
    end
    if any(drop_height(state, θ) < 0 for θ in check_angles)
        return Inf
    end
    if cp < 0
        return Inf
    end

    # Tangent slope at boundary
    θ_in  = theta_vec[cp]
    θ_out = theta_vec[cp + 1]
    err   = drop_height(state, θ_out) - drop_height(state, θ_in)
    return err
end
```

- [ ] **Step 3: Run tests — pass**

- [ ] **Step 4: Commit**

```bash
git add julia/src/contact.jl julia/test/test_contact.jl
git commit -m "feat: contact detection and tangent-slope error metric"
```

---

### Task 10: Adaptive BDF time-stepper

**Files:**
- Create: `km-viscous-drop/julia/src/timestepper.jl`

This is the main simulation loop. Mirror `solve_motion_v2.m`: try cp ∈ {cp_prev−1, cp_prev, cp_prev+1} at each step, accept the one with minimum |contact_error|, halve dt if Newton fails, upgrade to BDF2 after 10 steps, save states periodically.

- [ ] **Step 1: Implement `src/timestepper.jl`**

```julia
# src/timestepper.jl

"""
    solve_drop!(cfg, ob, init_state; t_end, save_every, dt_init) → (times, states)

Main time-integration loop.  Uses adaptive BDF1/BDF2 with contact search.
"""
function solve_drop!(cfg::SimConstants, ob::OBParams, init::DropState;
                     t_end::Float64 = 10.0,
                     save_every::Float64 = 0.1,
                     dt_init::Float64 = cfg.dt_max)

    history = [init]
    dt      = dt_init
    t       = init.t

    saved_times  = Float64[t]
    saved_states = DropState[deepcopy(init)]
    next_save    = t + save_every
    step_count   = 0

    while t < t_end
        order = min(length(history), 2)

        # Try cp candidates
        cp_prev = history[end].cp
        best_state = nothing
        best_err   = Inf

        for cp in max(0, cp_prev - 1) : cp_prev + 1
            X0 = pack_X(history[end], cfg.M)

            R! = (buf, X) -> begin
                s = deepcopy(history[end]); unpack_X!(s, X, cfg.M); s.cp = cp
                fill!(buf, 0.0)
                build_residual!(buf, s, history[order==2 ? [end-1,end] : [end]], dt, cp, cfg, ob)
            end

            J_fn = X -> begin
                s = deepcopy(history[end]); unpack_X!(s, X, cfg.M); s.cp = cp
                build_jacobian(s, history[order==2 ? [end-1,end] : [end]], dt, cp, cfg, ob)
            end

            X = copy(X0)
            key = (cp, round(dt; sigdigits=4), order)
            converged = newton_solve!(X, R!, J_fn; cache_key=key)

            if converged
                candidate = deepcopy(history[end])
                unpack_X!(candidate, X, cfg.M)
                candidate.t   = t + dt
                candidate.dt  = dt
                candidate.cp  = cp
                err = abs(contact_error(candidate, cfg.theta_vec, cp))
                if err < best_err
                    best_err   = err
                    best_state = candidate
                end
            end
        end

        if best_state === nothing || best_err == Inf
            # Step failed: halve dt and retry
            dt /= 2
            if dt * cfg.dt_max < 1e-12
                error("Time step underflow at t=$t")
            end
            clear_jac_cache!()
            continue
        end

        # Accept the step
        t = best_state.t
        push!(history, best_state)
        step_count += 1
        if length(history) > 2; popfirst!(history); end  # keep at most 2

        # Try to recover to dt_max
        dt = min(dt * 1.1, cfg.dt_max)

        if t >= next_save
            push!(saved_times, t)
            push!(saved_states, deepcopy(best_state))
            next_save = t + save_every
        end
    end

    return saved_times, saved_states
end
```

- [ ] **Step 2: Smoke test: free l=2 mode oscillation for 2 periods**

Add to `test/test_newtonian.jl` (see Task 11).

- [ ] **Step 3: Commit**

```bash
git add julia/src/timestepper.jl
git commit -m "feat: adaptive BDF1/BDF2 time-stepper with contact search"
```

---

## Chunk 3: OB Extension and Validation

### Task 11: Oldroyd-B polymer stress block

**Files:**
- Create: `km-viscous-drop/julia/src/ob_extension.jl`
- Create: `km-viscous-drop/julia/test/test_ob.jl`

When `ob.De1 > 0` and `ob.beta_s < 1`, the state vector gains S₂…S_M. `extend_residual_ob!` appends M−1 rows for Block S. `extend_jacobian_ob!` augments the Jacobian with the corresponding block rows/columns.

State vector (OB, length 4M):
```
X_ob = [A₂…A_M, Ȧ₂…Ȧ_M, S₂…S_M, B₀…B_M, z_COM, v_COM]
```

- [ ] **Step 1: Write the failing tests**

```julia
# test/test_ob.jl
using Test, LinearAlgebra
using DropSolver

@testset "Oldroyd-B extension" begin
    M = 4; Oh = 0.05; Fr = 200.0

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, 0.005)
    ob_newtonian = OBParams(0.0, 1.0)     # Newtonian OB
    ob_ob        = OBParams(0.5, 0.5)     # actual OB fluid

    @testset "OB state vector length is 4M" begin
        s = DropState(M)
        X = pack_X_ob(s, M)
        @test length(X) == 4M
    end

    @testset "OB residual, zero state → zero" begin
        s0 = DropState(M); s0.z = 1.0; s0.dt = 0.005
        R  = zeros(4M)
        build_residual_ob!(R, s0, [s0], 0.005, 0, cfg, ob_ob)
        @test norm(R) < 1e-12
    end

    @testset "Newtonian OB ≡ Newtonian (β_s=1, De₁=0)" begin
        s0 = DropState(M); s0.A[2] = 0.05; s0.z = 1.0; s0.dt = 0.005
        # Newtonian residual
        R_N  = zeros(3M+1)
        build_residual!(R_N, s0, [s0], 0.005, 0, cfg, ob_newtonian)
        # OB with Newtonian params
        R_OB = zeros(4M)
        build_residual_ob!(R_OB, s0, [s0], 0.005, 0, cfg, ob_newtonian)
        # The first 3M+1 entries of R_OB should match R_N
        @test R_OB[1:3M+1] ≈ R_N  atol=1e-14
        # The last M-1 entries (S block) should be zero (S=0, De₁=0 → Sₙ=0)
        @test norm(R_OB[3M+2:end]) < 1e-14
    end

    @testset "OB Jacobian matches finite differences" begin
        s0 = DropState(M); s0.A[2] = 0.02; s0.z = 1.0; s0.dt = 0.005
        J   = build_jacobian_ob(s0, [s0], 0.005, 0, cfg, ob_ob)
        X0  = pack_X_ob(s0, M)
        buf = zeros(4M)
        R0  = zeros(4M); build_residual_ob!(R0, s0, [s0], 0.005, 0, cfg, ob_ob)
        eps = 1e-6
        J_fd = zeros(4M, 4M)
        for j in 1:4M
            Xp = copy(X0); Xp[j] += eps
            sp = deepcopy(s0); unpack_X_ob!(sp, Xp, M)
            Rp = zeros(4M); build_residual_ob!(Rp, sp, [s0], 0.005, 0, cfg, ob_ob)
            J_fd[:, j] = (Rp .- R0) ./ eps
        end
        @test norm(J .- J_fd) / (norm(J_fd) + 1e-14) < 1e-4
    end
end
```

- [ ] **Step 2: Implement `src/ob_extension.jl`**

```julia
# src/ob_extension.jl

function pack_X_ob(s::DropState, M::Int)
    vcat(s.A[2:end], s.Adot[2:end], s.S[2:end], s.B, s.z, s.v)
end

function unpack_X_ob!(s::DropState, X::AbstractVector, M::Int)
    Nm = M - 1
    s.A[2:end]    .= X[1:Nm]
    s.Adot[2:end] .= X[Nm+1:2Nm]
    s.S[2:end]    .= X[2Nm+1:3Nm]
    s.B           .= X[3Nm+1:3Nm+M+1]
    s.z            = X[end-1]
    s.v            = X[end]
end

"""
    build_residual_ob!(R, state, history, dt, cp, cfg, ob)

Fill the 4M-length OB residual. Calls `build_residual!` for the first 3M+1
entries, then appends the polymer Block S (M-1 rows).
"""
function build_residual_ob!(R::AbstractVector, state::DropState,
                             history::Vector{DropState}, dt::Float64,
                             cp::Int, cfg::SimConstants, ob::OBParams)
    M     = cfg.M
    Nm    = M - 1
    order = length(history)
    c     = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak    = c[end]

    # Newtonian part (first 3M+1 entries)
    build_residual!(view(R, 1:3M+1), state, history, dt, cp, cfg, ob)

    # Block S: De₁·dSₙ/dτ = (1-β_s)·Ȧₙ - Sₙ
    # BDF: ak*Sₙ^k + Σ c[j]*Sₙ^{j} = dt*((1-β_s)*Ȧₙ^k - Sₙ^k) / De₁
    # Rearranged: ak*Sₙ^k + Σ c[j]*Sₙ^{j} - dt*(1-β_s)/De₁*Ȧₙ^k + dt/De₁*Sₙ^k = 0
    S    = state.S[2:end]
    Adot = state.Adot[2:end]
    prev_S = hcat([h.S[2:end] for h in history]...)

    if ob.De1 > 0
        R[3M+2:4M] .= ak .* S .+ sum(c[j] .* prev_S[:, j] for j in 1:order) .+
                      dt ./ ob.De1 .* S .- dt .* (1 - ob.beta_s) ./ ob.De1 .* Adot
    else
        # De₁=0: Sₙ forced to 0 instantly (no dynamics)
        R[3M+2:4M] .= S
    end
end

"""
    build_jacobian_ob(state, history, dt, cp, cfg, ob) → Matrix

Analytical 4M × 4M Jacobian for the OB system.
"""
function build_jacobian_ob(state::DropState, history::Vector{DropState},
                            dt::Float64, cp::Int, cfg::SimConstants, ob::OBParams)
    M     = cfg.M
    Nm    = M - 1
    order = length(history)
    c     = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak    = c[end]

    J = zeros(4M, 4M)
    # Top-left (3M+1 × 3M+1): Newtonian Jacobian
    J_N = build_jacobian(state, history, dt, cp, cfg, ob)
    J[1:3M+1, 1:3M+1] .= J_N

    # Block S row (rows 3M+2:4M):
    # ∂R_S/∂Adot = -dt*(1-β_s)/De₁  (columns Nm+1:2Nm)
    # ∂R_S/∂S    = ak + dt/De₁       (columns 2Nm+1:3Nm)
    if ob.De1 > 0
        J[3M+2:4M, Nm+1:2Nm]  .= -dt * (1 - ob.beta_s) / ob.De1 .* I(Nm)
        J[3M+2:4M, 2Nm+1:3Nm] .= (ak + dt / ob.De1) .* I(Nm)
    else
        J[3M+2:4M, 2Nm+1:3Nm] .= I(Nm)   # De₁=0: Sₙ=0
    end

    # Update Newtonian R2 block: ∂R2/∂S = dt*D2 (from damping term β_s*Ȧ + S)
    if ob.De1 > 0 && ob.beta_s < 1
        ns = 2:M
        D2_diag = @. 2 * cfg.Oh * (ns - 1) * (2ns + 1)
        J[Nm+1:2Nm, 2Nm+1:3Nm] .= dt .* diagm(D2_diag)
    end

    return J
end
```

- [ ] **Step 3: Run the OB tests — all pass**

- [ ] **Step 4: Commit**

```bash
git add julia/src/ob_extension.jl julia/test/test_ob.jl
git commit -m "feat: Oldroyd-B polymer stress auxiliary ODE (Block S, 4M state vector)"
```

---

### Task 12: Newtonian validation

**Files:**
- Create: `km-viscous-drop/julia/test/test_newtonian.jl`
- Create: `km-viscous-drop/julia/scripts/run_newtonian.jl`

Validate the solver reproduces the Lamb small-viscosity prediction for a free l=2 oscillation (no substrate): frequency ω ≈ √(l(l−1)(l+2)) and decay γ ≈ (l−1)(2l+1)·Oh with l=2.

- [ ] **Step 1: Write the test**

```julia
# test/test_newtonian.jl
using Test, LinearAlgebra
using DropSolver

@testset "Newtonian validation: Lamb l=2 limit" begin
    M  = 6; Oh = 0.02; Fr = 1e6   # nearly zero gravity
    l  = 2
    ω_lamb = sqrt(l*(l-1)*(l+2))         # = sqrt(8) ≈ 2.828
    γ_lamb = (l-1)*(2l+1)*Oh             # = 5 * 0.02 = 0.1

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(M*(M+2)*(M-1)) * 8)
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
    ob        = OBParams()

    # Initial condition: small A₂ = 0.05, no velocity, above substrate
    init = DropState(M)
    init.A[2]  = 0.05
    init.z     = 2.0       # far above substrate
    init.v     = 0.0
    init.dt    = dt_max
    init.cp    = 0

    # Integrate 4 full periods
    T_period  = 2π / ω_lamb
    times, states = solve_drop!(cfg, ob, init;
                                t_end    = 4T_period,
                                save_every = T_period/50,
                                dt_init  = dt_max)

    # Extract A₂(t)
    A2 = [s.A[2] for s in states]

    # Fit exponential envelope: |A₂(t)| ≈ A₀ * exp(-γ*t)
    # Use first and last saved point
    t1, t2 = times[1], times[end]
    γ_fit  = -log(abs(A2[end]) / abs(A2[1])) / (t2 - t1)
    @test abs(γ_fit - γ_lamb) / γ_lamb < 0.15   # 15% tolerance

    # Estimate oscillation frequency from zero crossings
    sign_changes = findall(i -> A2[i]*A2[i+1] < 0, 1:length(A2)-1)
    if length(sign_changes) >= 2
        # Each pair of sign changes = one half period
        half_periods = diff(times[sign_changes])
        ω_fit = π / mean(half_periods)
        @test abs(ω_fit - ω_lamb) / ω_lamb < 0.05  # 5% tolerance
    end
end
```

- [ ] **Step 2: Run the test** (this exercises the full stack for the first time)

```bash
julia --project=. test/test_newtonian.jl
```

- [ ] **Step 3: Fix any issues in residual/Jacobian index alignment** (iterate)

- [ ] **Step 4: Write `scripts/run_newtonian.jl`** as a runnable example that prints a convergence table comparing γ_fit and ω_fit to Lamb for Oh = 0.01, 0.05, 0.1.

- [ ] **Step 5: Commit**

```bash
git add julia/test/test_newtonian.jl julia/scripts/run_newtonian.jl
git commit -m "test: Newtonian l=2 oscillation validated against Lamb limit"
```

---

### Task 13: OB validation against Zrnić data

**Files:**
- Create: `km-viscous-drop/julia/scripts/run_ob_case.jl`

Zrnić & Brenn (2024) Table 1 gives dimensionless frequencies and decay rates for the l=2 mode at specific (Oh, De₁, β_s) values. Use these as reference.

Zrnić reference: for l=2, Oh=0.01, De₁=0.5, β_s=0.5 (UCM-like):
- Mode frequency shifts upward relative to Newtonian.
- Decay rate decreases (viscoelastic drop rings longer).

- [ ] **Step 1: Add OB mode test to `test/test_ob.jl`**

```julia
@testset "OB l=2: decay rate decreases with De₁" begin
    M  = 6; Oh = 0.02; Fr = 1e6
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(M*(M+2)*(M-1)) * 8)
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.A[2] = 0.05; init.z = 2.0; init.dt = dt_max

    # Newtonian run
    ob_N    = OBParams(0.0, 1.0)
    _, st_N = solve_drop!(cfg, ob_N, init; t_end=30.0, save_every=0.1)
    A2_N    = [s.A[2] for s in st_N]

    # OB run with De₁=0.5, β_s=0.5
    ob_OB    = OBParams(0.5, 0.5)
    _, st_OB = solve_drop!(cfg, ob_OB, deepcopy(init); t_end=30.0, save_every=0.1)
    A2_OB    = [s.A[2] for s in st_OB]

    # OB decay should be slower (drop rings longer)
    γ_N  = -log(abs(A2_N[end]) / abs(A2_N[1])) / (st_N[end].t - st_N[1].t)
    γ_OB = -log(abs(A2_OB[end]) / abs(A2_OB[1])) / (st_OB[end].t - st_OB[1].t)
    @test γ_OB < γ_N    # viscoelastic drop rings longer
end
```

- [ ] **Step 2: Run and verify**

- [ ] **Step 3: Write `scripts/run_ob_case.jl`** that produces a table of (De₁, β_s) → (ω, γ) and compare to Zrnić when data is available.

- [ ] **Step 4: Final commit**

```bash
git add julia/test/test_ob.jl julia/scripts/run_ob_case.jl
git commit -m "test: OB l=2 validation — viscoelastic drop decays slower than Newtonian"
```

---

## Summary of deliverables

| File | Purpose |
|------|---------|
| `src/types.jl` | Data structures |
| `src/legendre.jl` | Legendre polynomial evaluations |
| `src/integrals.jl` | Precomputed ∫Pₙ/u³ integrals |
| `src/bdf.jl` | BDF1/BDF2 coefficient formulas |
| `src/residual.jl` | 7-block Newtonian residual |
| `src/jacobian.jl` | Analytical 7-block Jacobian |
| `src/ob_extension.jl` | OB polymer Block S, 4M system |
| `src/newton.jl` | Newton-Raphson with Jacobian caching |
| `src/contact.jl` | Contact detection and error metric |
| `src/timestepper.jl` | Adaptive BDF time-stepper |
| `test/test_newtonian.jl` | Lamb-limit frequency/decay validation |
| `test/test_ob.jl` | OB Newtonian limit + Zrnić qualitative check |

**Validation criterion:** All tests pass; Newtonian l=2 frequency matches Lamb to within 5% and decay to within 15%; OB decay rate is strictly lower than Newtonian for De₁>0, β_s<1.
