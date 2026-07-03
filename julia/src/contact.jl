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

"""
    cyl_radius(state, θ) → Float64

Cylindrical (in-plane) radius of the interface, x(θ) = sin(θ)·(1 + Σₙ Pₙ(cosθ)·Aₙ),
mode index matching Legendre index (n = 2..M). Companion to `drop_height`, which
returns the height Z(θ) = cos(θ)·(1 + Σ) + z.
"""
function cyl_radius(state::DropState, θ::Float64)
    M = length(state.A)
    P = collect_Pl(M, [cos(θ)])
    Σ = sum(P[1, n+1] * state.A[n] for n in 2:M)
    return sin(θ) * (1.0 + Σ)
end

"""
    contact_angle(state, θ) → Float64

Apparent (dynamic) contact angle θ_d that the free surface makes with the substrate at
polar angle θ, measured through the liquid (θ measured from the north pole; contact pole
at θ=π). Derived in docs/DropRebound_ContactLine.tex eq (thetad), transcribed to the
south-pole-at-π convention: θ_d = π − atan2(−Z_θ, −x_θ), with

    ζ    = 1 + Σₙ Aₙ Pₙ(cosθ),        ζ_θ = −sinθ · Σₙ Aₙ P'ₙ(cosθ),
    x_θ  = ζ_θ sinθ + ζ cosθ,          Z_θ = ζ_θ cosθ − ζ sinθ.

For an undeformed sphere (Aₙ≡0) this gives θ_d = θ (Lemma "undeformed reduction").
"""
function contact_angle(state::DropState, θ::Float64)
    M  = length(state.A)
    c  = cos(θ); s = sin(θ)
    P  = collect_Pl(M, [c])
    dP = collect_dPl(M, [c])
    ζ  = 1.0 + sum(P[1,  n+1] * state.A[n] for n in 2:M)
    ζθ = -s  * sum(dP[1, n+1] * state.A[n] for n in 2:M)
    xθ = ζθ * s + ζ * c
    Zθ = ζθ * c - ζ * s
    return π - atan(-Zθ, -xθ)
end

"""
    contact_angle_error(state, theta_vec, cp, theta_e) → Float64

Generalization of `contact_error` to a finite equilibrium contact angle θ_e (the
Milestone-1 contact-set selection residual; see docs/DropRebound_ContactLine.tex eq
(disc_angle)). Returns 0 if cp≤0, Inf on penetration, else the discrete contact-angle
residual across the contact edge,

    | Δx · sin(θ_e) + ΔZ · cos(θ_e) |,

with ΔZ, Δx the differences of height and cylindrical radius between the first free node
and the last contact node. At θ_e = π (sin=0, cos=−1) this reduces to |ΔZ| =
|`contact_error`|, i.e. GA's tangency residual, term for term. Linear in the amplitudes
and bounded (no 1/x_θ), so it stays well defined near maximal spreading.
"""
function contact_angle_error(state::DropState, theta_vec::Vector{Float64},
                             cp::Int, theta_e::Real)
    if cp <= 0
        return 0.0
    end
    # Penetration checks: identical to contact_error
    for θ in theta_vec[cp+1:end]
        if drop_height(state, θ) < 0
            return Inf
        end
    end
    if drop_height(state, theta_vec[cp]) < -1e-10
        return Inf
    end

    θ_in  = theta_vec[cp]
    θ_out = theta_vec[min(cp + 1, length(theta_vec))]
    ΔZ = drop_height(state, θ_out) - drop_height(state, θ_in)
    Δx = cyl_radius(state,  θ_out) - cyl_radius(state,  θ_in)
    return abs(Δx * sin(theta_e) + ΔZ * cos(theta_e))
end

"""
    contact_edge_radius(state, theta_vec) → Float64

Cylindrical radius of the contact line for the given state, evaluated at the midpoint
angle between the last contact node and the first free node (GA's boundary convention).
Returns 0.0 when cp ≤ 0.
"""
function contact_edge_radius(state::DropState, theta_vec::Vector{Float64})
    cp = state.cp
    cp <= 0 && return 0.0
    θ_edge = (theta_vec[cp] + theta_vec[min(cp + 1, length(theta_vec))]) / 2
    return cyl_radius(state, θ_edge)
end

