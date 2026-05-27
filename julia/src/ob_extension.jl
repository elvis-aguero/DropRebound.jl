"""
    pack_X_ob(s, M) → Vector{Float64}

Pack a DropState into the OB extended Newton state vector (length 4M).

Layout:
  X[1      : M-1]   = A₂…A_M      (Nm = M-1 deformation amplitudes)
  X[M      : 2M-2]  = Ȧ₂…Ȧ_M
  X[2M-1   : 3M-3]  = S₂…S_M     (Nm polymer stress aux vars)
  X[3M-2   : 4M-2]  = B₀…B_M     (M+1 pressure coefficients)
  X[4M-1]           = z
  X[4M]             = v
"""
function pack_X_ob(s::DropState, M::Int)
    vcat(s.A[2:end], s.Adot[2:end], s.S[2:end], s.B, s.z, s.v)
end

"""
    unpack_X_ob!(state, X, M)

Unpack OB extended vector X (length 4M) back into a DropState in-place.
"""
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

Fill the 4M-length OB residual in-place.

First 3M+1 entries: Newtonian part (calls build_residual!).
Last M-1 entries (rows 3M+2:4M): Block S polymer ODE.

Block S: De₁·dSₙ/dτ = (1-β_s)·Ȧₙ - Sₙ
BDF form: ak·Sₙ^k + Σcⱼ·Sₙ^{k-order+j} + (dt/De₁)·Sₙ^k - dt(1-β_s)/De₁·Ȧₙ^k = 0
i.e.      (ak + dt/De₁)·Sₙ^k + Σcⱼ·Sₙ^{k-order+j} - dt(1-β_s)/De₁·Ȧₙ^k = 0

When De₁=0: polymer stress is instantaneously zero, residual is just Sₙ.
"""
function build_residual_ob!(R::AbstractVector, state::DropState,
                             history::Vector{DropState}, dt::Float64,
                             cp::Int, cfg::SimConstants, ob::OBParams)
    M     = cfg.M
    order = length(history)
    c     = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)

    # Newtonian part (first 3M+1 entries)
    build_residual!(view(R, 1:3M+1), state, history, dt, cp, cfg, ob)

    # Block S (rows 3M+2:4M, length M-1)
    S    = state.S[2:end]
    Adot = state.Adot[2:end]
    prev_S = hcat([h.S[2:end] for h in history]...)   # (M-1 × order)

    if ob.De1 > 0.0
        R[3M+2:4M] .= (c[end] + dt/ob.De1) .* S .+
                      sum(c[j] .* prev_S[:, j] for j in 1:order) .-
                      dt * (1 - ob.beta_s) / ob.De1 .* Adot
    else
        # De₁=0: polymer stress forced to zero instantly
        R[3M+2:4M] .= S
    end
end

"""
    build_jacobian_ob(state, history, dt, cp, cfg, ob) → Matrix{Float64} (4M × 4M)

Analytical Jacobian for the OB extended system.

State-vector layout (length 4M):
  X[1:Nm]         = A₂…A_M      (Nm = M-1)
  X[Nm+1:2Nm]     = Ȧ₂…Ȧ_M
  X[2Nm+1:3Nm]    = S₂…S_M
  X[3Nm+1:3Nm+M+1]= B₀…B_M
  X[4M-1]         = z
  X[4M]           = v
"""
function build_jacobian_ob(state::DropState, history::Vector{DropState},
                            dt::Float64, cp::Int, cfg::SimConstants, ob::OBParams)
    M     = cfg.M
    Nm    = M - 1
    order = length(history)
    c     = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak    = c[end]

    J = zeros(4M, 4M)

    # ── OB state vector X_ob layout (length 4M) ─────────────────────────────
    # X_ob[1:Nm]            = A₂…A_M
    # X_ob[Nm+1:2Nm]        = Ȧ₂…Ȧ_M
    # X_ob[2Nm+1:3Nm]       = S₂…S_M
    # X_ob[3Nm+1:3Nm+M+1]   = B₀…B_M   (length M+1)
    # X_ob[4M-1]            = z
    # X_ob[4M]              = v
    #
    # OB residual R_ob layout (length 4M):
    # R_ob[1:3M+1]   = Newtonian residual R_N (rows unchanged from Newtonian layout)
    # R_ob[3M+2:4M]  = Block S polymer ODE
    #
    # IMPORTANT: R_ob rows 1:3M+1 follow the Newtonian *residual* row layout:
    #   rows 1:Nm          = R1  (A kinematics)
    #   rows Nm+1:2Nm      = R2  (Adot equation of motion)
    #   rows 2Nm+1:3Nm+2   = R3/R4 (pressure BCs, same as Newtonian rows 2M-1:3M-1)
    #   row  3M            = R6  (z kinematics)
    #   row  3M+1          = R7  (v equation of motion)
    #
    # Residual rows are UNCHANGED (no row permutation).
    # Only the columns are permuted: S block inserted between Adot and B.
    # Column mapping: Newtonian → OB
    #   N cols 1:Nm       → OB cols 1:Nm          (A, unchanged)
    #   N cols Nm+1:2Nm   → OB cols Nm+1:2Nm      (Adot, unchanged)
    #   N cols 2Nm+1:3Nm+2 → OB cols 3Nm+1:4Nm+2 (B, shifted by Nm)
    #   N col  3Nm+3      → OB col  4Nm+3 = 4M-1  (z, shifted by Nm)
    #   N col  3Nm+4      → OB col  4Nm+4 = 4M    (v, shifted by Nm)

    J_N = build_jacobian(state, history, dt, cp, cfg, ob)

    # Newtonian column groups (in Newtonian X_N)
    N_A_cols    = 1:Nm
    N_Adot_cols = Nm+1:2Nm
    N_B_cols    = 2Nm+1:3Nm+2    # 2M-1:3M-1
    N_z_col     = 3Nm+3          # 3M
    N_v_col     = 3Nm+4          # 3M+1

    # Corresponding OB column groups
    OB_A_cols    = 1:Nm
    OB_Adot_cols = Nm+1:2Nm
    OB_B_cols    = 3Nm+1:4Nm+2   # = 3Nm+1:3Nm+M+1
    OB_z_col     = 4M-1
    OB_v_col     = 4M

    # Build column permutation vectors
    N_col_perm  = [collect(N_A_cols);    collect(N_Adot_cols);    collect(N_B_cols);    N_z_col;    N_v_col]
    OB_col_perm = [collect(OB_A_cols); collect(OB_Adot_cols); collect(OB_B_cols); OB_z_col; OB_v_col]

    # Rows are unchanged (Newtonian residual row i → OB residual row i for i=1:3M+1)
    # Copy Newtonian Jacobian into OB Jacobian with column remapping only
    for (nc, oc) in zip(N_col_perm, OB_col_perm)
        J[1:3M+1, oc] = J_N[1:3M+1, nc]
    end

    # ── Additional OB coupling: S affects R2 ─────────────────────────────────
    # R2 (rows Nm+1:2Nm) uses effective_damp = D2*(β_s*Ȧ + S) when De1>0 && beta_s<1
    # → ∂R2/∂S = dt * D2  (rows Nm+1:2Nm, OB S-cols 2Nm+1:3Nm)
    # Note: the Newtonian Jacobian already has ∂R2/∂Ȧ with damp_factor = beta_s
    # so we only need to add the S contribution here.
    if ob.De1 > 0.0 && ob.beta_s < 1.0
        ns = collect(Float64, 2:M)
        D2_diag = @. 2 * cfg.Oh * (ns - 1) * (2ns + 1)
        J[Nm+1:2Nm, 2Nm+1:3Nm] .= dt .* diagm(D2_diag)

        # Block S rows (3M+2:4M):
        # ∂R_S/∂Ȧ  = -dt*(1-β_s)/De₁   (OB cols Nm+1:2Nm)
        # ∂R_S/∂S  = ak + dt/De₁        (OB cols 2Nm+1:3Nm)
        J[3M+2:4M, Nm+1:2Nm]  .= -dt * (1 - ob.beta_s) / ob.De1 * I(Nm)
        J[3M+2:4M, 2Nm+1:3Nm] .= (ak + dt / ob.De1) * I(Nm)
    else
        # De₁=0: ∂R_S/∂S = I
        J[3M+2:4M, 2Nm+1:3Nm] .= I(Nm)
    end

    return J
end
