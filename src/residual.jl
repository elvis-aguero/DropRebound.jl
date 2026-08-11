"""
    pack_X(state, M) → Vector{Float64}

Pack a DropState into the Newton state vector (Newtonian, length 3M+1).

Layout:
  X[1      : M-1]   = A₂…A_M      (deformation amplitudes)
  X[M      : 2M-2]  = Ȧ₂…Ȧ_M     (velocities)
  X[2M-1   : 3M-1]  = B₀…B_M     (pressure coefficients, length M+1)
  X[3M]             = z           (COM height)
  X[3M+1]           = v           (COM velocity)
"""
function pack_X(s::DropState, M::Int)
    vcat(s.A[2:end], s.Adot[2:end], s.B, s.z, s.v)
end

"""
    unpack_X!(state, X, M)

Unpack Newton vector X back into a DropState in-place.
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

Fill residual vector R (length 3M+1) in-place. Seven blocks:
  R1: kinematics dA/dτ = Ȧ
  R2: equation of motion Äₙ + damping + capillary + pressure = 0
  R3: contact BCs (cp equations)
  R4: free-surface pressure BCs (M+1-cp equations)
  R6: dz/dτ = v
  R7: dv/dτ = -Bo - cp*B₁

- history: Vector of DropState, most recent last. length=1 → BDF1, length=2 → BDF2.
- cp: number of contact angles (0 = no contact)
"""
function build_residual!(R::AbstractVector, state::DropState,
                         history::Vector{DropState}, dt::Float64,
                         cp::Int, cfg::SimConstants, ob::OBParams)
    M  = cfg.M
    Bo = cfg.Bo
    θv = cfg.theta_vec
    order = length(history)

    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)

    # Unpack current state (modes 2…M, 1-based Julia index 2:M)
    A    = state.A[2:end]        # length M-1
    Adot = state.Adot[2:end]     # length M-1
    B    = state.B               # length M+1 (B₀…B_M)
    z    = state.z
    v    = state.v
    ns   = 2:M                   # mode numbers

    # History arrays: each is (M-1 × order)
    prev_A    = hcat([h.A[2:end]    for h in history]...)
    prev_Adot = hcat([h.Adot[2:end] for h in history]...)
    prev_z    = [h.z for h in history]
    prev_v    = [h.v for h in history]

    # ── Block R1: BDF(dA/dτ = Ȧ) ────────────────────────────────────────────
    # BDF form: c[end]*A^k + Σ c[j]*A^{k-order+j} - dt*Ȧ^k = 0
    R[1:M-1] .= c[end] .* A .+ sum(c[j] .* prev_A[:, j] for j in 1:order) .- dt .* Adot

    # ── Block R2: BDF(Äₙ + 2λₙ(Oh)*Ȧₙ + ωₙ²(Oh)*Aₙ + cp*n*Bₙ = 0) ──────────
    # λₙ, ωₙ² are Lamb's small-Oh asymptotics by default, or Reid's exact
    # finite-Oh values when cfg.viscous == :reid (see src/reid.jl).
    D1 = cfg.drop_omega2   # restoring/frequency-squared
    D2 = 2 .* cfg.drop_lambda   # viscous damping

    # For OB: replace Adot term with (β_s*Ȧₙ + Sₙ)
    if ob.De1 > 0.0 && ob.beta_s < 1.0
        S_modes = state.S[2:end]
        effective_damp = D2 .* (ob.beta_s .* Adot .+ S_modes)
    else
        effective_damp = D2 .* Adot
    end

    # Pressure term: only modes 2…M, using B[3:end] (B₂…B_M)
    pressure_term = cp > 0 ? Float64.(ns) .* B[3:end] : zeros(M-1)

    R[M:2M-2] .= c[end] .* Adot .+ sum(c[j] .* prev_Adot[:, j] for j in 1:order) .+
                 dt .* (D1 .* A .+ effective_damp .+ pressure_term)

    # ── Blocks R3+R4: M+1 equations for M+1 unknowns B₀…B_M ────────────────
    # R3 (contact): cos(θ_c)*(1 + ΣₙPₙ(cosθ_c)*Aₙ) + z = 0  for cp contact angles
    # R4 (free):    ΣₙBₙ*Pₙ(cosθ) = 0  for (M+1-cp) free angles

    # Starting row for the pressure block: 2M-1 (1-based)
    row_p = 2M - 1

    if cp > 0
        cos_c = cos.(θv[1:cp])
        Pc    = collect_Pl(M, cos_c)    # (cp × M+1), column k = P_{k-1}
        PcA   = Pc[:, 3:end]            # (cp × Nm), P₂…P_M evaluated at contact points
        R3    = cos_c .* (1.0 .+ PcA * A) .+ z
        R[row_p : row_p + cp - 1] .= R3
    end

    n_free = M + 1 - cp
    if n_free > 0
        cos_f = cos.(θv[cp+1:end])
        Pf    = collect_Pl(M, cos_f)         # (n_free × M+1), column k = P_{k-1}
        # ΣBₙ*Pₙ: B[1]=B₀ (P₀=1), B[2]=B₁ (P₁=cosθ), ..., B[M+1]=B_M (P_M)
        R4    = Pf * B   # Pf is (n_free × M+1), B is length M+1
        R[row_p + cp : row_p + M] .= R4
    end

    # ── Block R6: BDF(dz/dτ = v) ────────────────────────────────────────────
    R[3M] = c[end] * z + sum(c[j] * prev_z[j] for j in 1:order) - dt * v

    # ── Block R7: BDF(dv/dτ = -Bo - cp*B₁) ──────────────────────────────────
    B1_force = cp > 0 ? B[2] : 0.0
    R[3M+1] = c[end] * v + sum(c[j] * prev_v[j] for j in 1:order) -
              dt * (-Bo - B1_force)
end
