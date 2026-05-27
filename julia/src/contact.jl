"""
    drop_height(state, θ) → Float64

Surface height at polar angle θ: z(θ) = cos(θ)*(1 + Σₙ Pₙ(cosθ)*Aₙ) + z_COM.
Note: modes in state.A are A[1]=A₁ (always 0), A[2]=A₂, ..., A[M]=A_M.
Summation is over n=1..M-1 (A[2..M] = A₂..A_M).
"""
function drop_height(state::DropState, θ::Float64)
    M = length(state.A)     # = M
    P = collect_Pl(M - 1, [cos(θ)])   # (1 × M), col k = P_{k-1}
    # Sum: Σₙ₌₁^{M-1} P_n(cosθ) * A_{n+1} ... but A[1]=0, A[2]=A₂...
    # collect_Pl gives P₀..P_{M-1}. We need P₁..P_{M-1} times A₂..A_M.
    # P[:,2] = P₁, P[:,3] = P₂, ..., P[:,M] = P_{M-1}
    Σ = sum(P[1, n+1] * state.A[n+1] for n in 1:M-1)
    return cos(θ) * (1.0 + Σ) + state.z
end

"""
    contact_error(state, theta_vec, cp) → Float64

Tangent-slope error at the contact boundary. Returns:
- 0.0 if cp ≤ 0 (no contact)
- Inf if penetration detected below the contact angle
- drop_height(θ_out) - drop_height(θ_in): should be ≥ 0 for valid contact
"""
function contact_error(state::DropState, theta_vec::Vector{Float64}, cp::Int)
    if cp <= 0
        return 0.0
    end

    # Check penetration: any angle more equatorial than θ_cp should have height ≥ 0
    for θ in theta_vec[cp+1:end]
        if drop_height(state, θ) < 0
            return Inf
        end
    end

    # Also check if the contact point itself is penetrating
    if drop_height(state, theta_vec[cp]) < -1e-10
        return Inf
    end

    # Tangent slope at the boundary between contact and free zone
    θ_in  = theta_vec[cp]
    θ_out = theta_vec[min(cp + 1, length(theta_vec))]
    return drop_height(state, θ_out) - drop_height(state, θ_in)
end
