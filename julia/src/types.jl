"""Physical and numerical parameters for one simulation run."""
struct SimConstants
    M         :: Int       # harmonics_qtt (number of Legendre modes, A₁…A_M)
    N_angles  :: Int       # angular_sampling = M+1
    Oh        :: Float64   # Ohnesorge number
    Bo        :: Float64   # Bond number (ρgR²/σ): gravity / surface tension
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

"""Carreau shear-thinning parameters. Set eps_ST=0 for Newtonian."""
struct STParams
    eps_ST   :: Float64             # (1-n)/2 ≥ 0; zero = Newtonian
    lambda_c :: Float64             # Carreau relaxation time (non-dimensional)
    Gamma    :: Vector{Float64}     # Γ_l for modes 2..M (length M-1)
                                    # Compute from shear_thinning_derivation.ipynb
end

STParams() = STParams(0.0, 0.0, Float64[])

"""State at a single time step."""
mutable struct DropState
    A          :: Vector{Float64}   # deformation amplitudes A₁…A_M (A₁ always 0)
    Adot       :: Vector{Float64}   # dA/dτ  (same indexing)
    S          :: Vector{Float64}   # polymer stress aux vars S₁…S_M (S₁ always 0; OB only)
    B          :: Vector{Float64}   # pressure amplitudes B₀…B_M (length M+1)
    z          :: Float64           # center-of-mass height
    v          :: Float64           # center-of-mass velocity
    t          :: Float64           # current simulation time
    dt         :: Float64           # step size used to reach this state
    cp         :: Int               # contact_points (discrete solver; unused in v1)
    theta_star :: Float64           # continuous contact angle θ* ∈ (π/2, π]; π = no contact
end

function DropState(M::Int)
    DropState(zeros(M), zeros(M), zeros(M), zeros(M+1),
              0.0, 0.0, 0.0, 0.0, 0, π)
end

"""
    make_theta_vec(M) → Vector{Float64}

Build the M+1-angle Gauss-Legendre collocation grid used by the solver.
Returns θ values sorted descending (south pole π first), computed as
θ = π prepended to acos(roots of P_M) via the Jacobi tridiagonal matrix.
Matches MATLAB: theta_vector = [pi; acos(vpasolve(legendreP(M, x)))].
"""
function make_theta_vec(M::Int)
    bvec  = [sqrt(Float64(k)^2 / (4*Float64(k)^2 - 1)) for k in 1:M-1]
    J     = diagm(1 => bvec, -1 => bvec)
    nodes = sort(eigvals(J))   # M Gauss-Legendre nodes in (-1, 1)
    sort([Float64(π); acos.(nodes)], rev=true)
end

"""
    make_dt_max(M) → Float64

CFL-like stability limit for M Legendre modes: dt ≤ 2π / (√(M(M+2)(M-1)) · 8).
Matches MATLAB: max_dt = 2*pi / (sqrt(M*(M+2)*(M-1)) * 8).
"""
make_dt_max(M::Int) = 2π / (sqrt(M * (M + 2) * (M - 1)) * 8)
