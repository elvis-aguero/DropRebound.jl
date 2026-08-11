"""
    build_jacobian(state, history, dt, cp, cfg, ob) → Matrix{Float64}

Compute the analytical Jacobian ∂R/∂X of the residual built by `build_residual!`.

State-vector layout (length 3M+1):
  X[1:Nm]       = A₂…A_M      (Nm = M-1 deformation amplitudes)
  X[Nm+1:2Nm]   = Ȧ₂…Ȧ_M
  X[2M-1:3M-1]  = B₀…B_M     (M+1 pressure coefficients)
  X[3M]         = z
  X[3M+1]       = v
"""
function build_jacobian(state::DropState, history::Vector{DropState},
                        dt::Float64, cp::Int, cfg::SimConstants, ob::OBParams)
    M  = cfg.M
    θv = cfg.theta_vec
    order = length(history)
    c  = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak = c[end]

    ns = collect(Float64, 2:M)   # mode numbers 2..M as Float64
    Nm = M - 1
    sz = 3M + 1

    # λₙ, ωₙ² are Lamb's small-Oh asymptotics by default, or Reid's exact
    # finite-Oh values when cfg.viscous == :reid (see src/reid.jl).
    D1 = Diagonal(cfg.drop_omega2)
    D2 = Diagonal(2 .* cfg.drop_lambda)
    damp_factor = (ob.De1 > 0.0 && ob.beta_s < 1.0) ? ob.beta_s : 1.0

    J = zeros(sz, sz)

    # ── ∂R1/∂A, ∂R1/∂Ȧ ─────────────────────────────────────────────────────
    # R1: c[end]*A + Σc[j]*A_prev[j] - dt*Ȧ = 0
    J[1:Nm, 1:Nm]     .= ak .* I(Nm)
    J[1:Nm, Nm+1:2Nm] .= -dt .* I(Nm)

    # ── ∂R2/∂A, ∂R2/∂Ȧ, ∂R2/∂B ─────────────────────────────────────────────
    # R2: c[end]*Ȧ + ... + dt*(D1*A + damp_factor*D2*Ȧ + ns.*B[3:end]) = 0
    J[Nm+1:2Nm, 1:Nm]     .= dt .* Matrix(D1)
    J[Nm+1:2Nm, Nm+1:2Nm] .= diagm(fill(ak, Nm) .+ dt .* damp_factor .* diag(D2))
    if cp > 0
        # Pressure term: B₂…B_M at X[2M+1:3M-1], coefficient dt*n for mode n=2..M
        J[Nm+1:2Nm, 2M+1:3M-1] .= dt .* diagm(ns)
    end

    # ── ∂R3/∂A, ∂R3/∂z (contact rows, only when cp > 0) ────────────────────
    # R3_i = cos(θ_c_i)*(1 + Σₙ₌₂ᴹ Pₙ(cosθ_c_i)*Aₙ) + z = 0
    if cp > 0
        cos_c = cos.(θv[1:cp])
        Pc    = collect_Pl(M, cos_c)    # (cp × M+1), column k = P_{k-1}
        PcA   = Pc[:, 3:end]            # P₂…P_M for modes A₂…A_M (Nm columns)
        # ∂R3/∂Aₙ = cos(θ) * Pₙ(cosθ) for n=2..M
        J[2M-1:2M-1+cp-1, 1:Nm]  .= cos_c .* PcA
        # ∂R3/∂z = 1
        J[2M-1:2M-1+cp-1, 3M]    .= 1.0
    end

    # ── ∂R4/∂B (free-surface rows) ──────────────────────────────────────────
    # R4 = Pf * B = 0 where B at X[2M-1:3M-1]
    n_free = M + 1 - cp
    if n_free > 0
        cos_f = cos.(θv[cp+1:end])
        Pf    = collect_Pl(M, cos_f)    # (n_free × M+1)
        J[2M-1+cp:3M-1, 2M-1:3M-1] .= Pf
    end

    # ── ∂R6/∂z, ∂R6/∂v ─────────────────────────────────────────────────────
    # R6: c[end]*z + ... - dt*v = 0
    J[3M,   3M]   = ak
    J[3M,   3M+1] = -dt

    # ── ∂R7/∂v, ∂R7/∂B₁ ────────────────────────────────────────────────────
    # R7: c[end]*v + ... - dt*(-Bo - cp*B₁) = 0
    J[3M+1, 3M+1] = ak
    if cp > 0
        # B₁ = X[2M], coefficient from -dt*(-B₁) = +dt*B₁
        J[3M+1, 2M] = dt
    end

    return J
end
