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
