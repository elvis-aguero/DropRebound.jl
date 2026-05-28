"""
    drop_height(state, θ) → Float64

Surface height at polar angle θ: z(θ) = cos(θ)*(1 + Σₙ Pₙ(cosθ)*Aₙ) + z_COM.
Note: modes in state.A are A[1]=A₁ (always 0), A[2]=A₂, ..., A[M]=A_M.
Summation uses Pₙ(cosθ)*Aₙ for n=2..M (mode index matches Legendre index).
"""
function drop_height(state::DropState, θ::Float64)
    M = length(state.A)     # = M
    P = collect_Pl(M, [cos(θ)])   # (1 × M+1), col k = P_{k-1}
    # Sum: Σₙ₌₂^M Pₙ(cosθ) * Aₙ
    # collect_Pl(M,...) gives columns P₀..P_M.
    # P[:,n+1] = P_n(cosθ). Mode Aₙ stored at state.A[n] (1-based).
    Σ = sum(P[1, n+1] * state.A[n] for n in 2:M)
    return cos(θ) * (1.0 + Σ) + state.z
end

"""
    update_theta_star!(state) → Bool

Scalar Newton solve to find θ* such that h(θ*) = drop_height(state, θ*) = 0.
Updates state.theta_star in-place.  Returns true if contact status changed.

Rules:
- If h(π) > 0 (south pole above substrate): no contact → θ* = π
- If h(π) ≤ 0: contact is active; find θ* ∈ (π/2, π) where h(θ*) = 0
  using bisection (robust) followed by one Newton polish step.
"""
function update_theta_star!(state::DropState)
    PI = Float64(π)
    was_contact = state.theta_star < PI - 1e-10

    if drop_height(state, PI) > 0.0
        state.theta_star = PI
        return was_contact
    end

    # Contact active: bisect on [π/2+ε, π] for the zero of h
    lo::Float64 = PI/2 + 1e-6
    hi::Float64 = PI
    h_lo = drop_height(state, lo)

    if h_lo <= 0.0
        state.theta_star = lo
        return !was_contact
    end

    for _ in 1:50
        mid = (lo + hi) / 2
        h_mid = drop_height(state, mid)
        abs(h_mid) < 1e-12 && break
        if h_mid > 0.0
            lo = mid
        else
            hi = mid
        end
        hi - lo < 1e-14 && break
    end
    θ_new = (lo + hi) / 2

    # One Newton polish
    h_val = drop_height(state, θ_new)
    M = length(state.A)
    P  = collect_Pl(M, [cos(θ_new)])
    Σ  = sum(P[1, n+1] * state.A[n] for n in 2:M)
    dP = sum(-sin(θ_new) * P[1, n+1] * state.A[n] for n in 2:M)
    dh = -sin(θ_new) * (1.0 + Σ) + cos(θ_new) * dP
    if abs(dh) > 1e-15
        θ_new = clamp(θ_new - h_val / dh, PI/2 + 1e-6, PI - 1e-10)
    end

    state.theta_star = θ_new
    return !was_contact
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
