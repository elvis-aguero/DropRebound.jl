"""
    build_residual_v1!(R, state, history, dt, cfg)

v1 (MATLAB-parity) residual for the continuous-θ* Newtonian solver.

Row layout (length 3M+1):
  R1  [1:M-1]       BDF kinematics dA/dτ = Ȧ
  R2  [M:2M-2]      BDF mode dynamics
  P-block [2M-1:3M-1]   M+1 rows for pressure B:
    No contact:  M+1 GL nodes on [-1,1], p=0 everywhere
    Contact:     row 2M-1 = contact-height h(θ*)=0;  rows 2M:3M-1 = M GL free-surface nodes
  R6  [3M]          BDF dz/dτ = v
  R7  [3M+1]        BDF dv/dτ = -1/Fr + F_contact

R7 contact force uses v1 formula: (3/2)·z_prev²·ΣBₙ·∫Pₙ(u)/u³du  (z_prev frozen → linear in B)
"""
function build_residual_v1!(R::AbstractVector, state::DropState,
                             history::Vector{DropState}, dt::Float64,
                             cfg::SimConstants)
    M    = cfg.M
    Oh   = cfg.Oh
    Fr   = cfg.Fr
    θs   = state.theta_star
    contact = θs < π - 1e-10

    order = length(history)
    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)

    A    = state.A[2:end]
    Adot = state.Adot[2:end]
    B    = state.B
    z    = state.z
    v    = state.v
    ns   = 2:M

    prev_A    = hcat([h.A[2:end]    for h in history]...)
    prev_Adot = hcat([h.Adot[2:end] for h in history]...)
    prev_z    = [h.z for h in history]
    prev_v    = [h.v for h in history]

    # ── R1 ──────────────────────────────────────────────────────────────────
    R[1:M-1] .= c[end] .* A .+ sum(c[j] .* prev_A[:, j] for j in 1:order) .- dt .* Adot

    # ── R2 ──────────────────────────────────────────────────────────────────
    D1 = @. Float64(ns) * (Float64(ns) + 2) * (Float64(ns) - 1)
    D2 = @. 2Oh * (Float64(ns) - 1) * (2*Float64(ns) + 1)
    pressure_term = contact ? Float64.(ns) .* B[3:end] : zeros(M-1)

    R[M:2M-2] .= c[end] .* Adot .+ sum(c[j] .* prev_Adot[:, j] for j in 1:order) .+
                 dt .* (D1 .* A .+ D2 .* Adot .+ pressure_term)

    # ── Pressure block [2M-1 : 3M-1] (M+1 rows) ────────────────────────────
    row_p = 2M - 1

    if contact
        # R3 (row_p): contact height h(θ*) = 0
        cθ = cos(θs)
        Pc = collect_Pl(M, [cθ])
        PA = Pc[1, 3:end]
        R[row_p] = cθ * (1.0 + dot(PA, A)) + z

        # R4 (rows row_p+1:row_p+M): p=0 at M GL nodes on free surface [cθ, 1]
        u_free, _ = gauss_legendre_nodes(M, cθ, 1.0)
        Pf = collect_Pl(M, u_free)
        R[row_p+1 : row_p+M] .= Pf * B
    else
        # No contact: M+1 GL nodes on [-1, 1], enforce p=0 everywhere
        u_all, _ = gauss_legendre_nodes(M+1, -1.0, 1.0)
        Pall = collect_Pl(M, u_all)
        R[row_p : row_p+M] .= Pall * B
    end

    # ── R6 ──────────────────────────────────────────────────────────────────
    R[3M] = c[end] * z + sum(c[j] * prev_z[j] for j in 1:order) - dt * v

    # ── R7 ──────────────────────────────────────────────────────────────────
    if contact
        z_prev = history[end].z
        Iv = integral_at_theta_star(θs, M)
        F_contact = (3/2) * z_prev^2 * dot(B, Iv)
    else
        F_contact = 0.0
    end

    R[3M+1] = c[end] * v + sum(c[j] * prev_v[j] for j in 1:order) -
              dt * (-1.0/Fr + F_contact)
end

"""
    build_jacobian_v1(state, history, dt, cfg) → Matrix{Float64}

Linear coefficient matrix for the v1 continuous-θ* solver.
Matches the row/column layout of build_residual_v1!.
"""
function build_jacobian_v1(state::DropState, history::Vector{DropState},
                            dt::Float64, cfg::SimConstants)
    M    = cfg.M
    Oh   = cfg.Oh
    θs   = state.theta_star
    contact = θs < π - 1e-10

    order = length(history)
    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    ak = c[end]

    ns = collect(Float64, 2:M)
    Nm = M - 1
    sz = 3M + 1

    D1 = Diagonal(ns .* (ns .+ 2) .* (ns .- 1))
    D2 = Diagonal(2Oh .* (ns .- 1) .* (2 .* ns .+ 1))

    J = zeros(sz, sz)

    # ── R1 ──────────────────────────────────────────────────────────────────
    J[1:Nm, 1:Nm]     .= ak .* I(Nm)
    J[1:Nm, Nm+1:2Nm] .= -dt .* I(Nm)

    # ── R2 ──────────────────────────────────────────────────────────────────
    J[Nm+1:2Nm, 1:Nm]     .= dt .* Matrix(D1)
    J[Nm+1:2Nm, Nm+1:2Nm] .= diagm(fill(ak, Nm) .+ dt .* diag(D2))
    if contact
        J[Nm+1:2Nm, 2M+1:3M-1] .= dt .* diagm(ns)
    end

    # ── Pressure block [row_p : row_p+M] ────────────────────────────────────
    row_p = 2M - 1

    if contact
        cθ = cos(θs)
        # R3: ∂h/∂A = cθ·Pₙ(cθ), ∂h/∂z = 1
        Pc  = collect_Pl(M, [cθ])
        PcA = Pc[1, 3:end]
        J[row_p, 1:Nm]  .= cθ .* PcA
        J[row_p, 3M]     = 1.0

        # R4: ∂(Pf·B)/∂B = Pf
        u_free, _ = gauss_legendre_nodes(M, cθ, 1.0)
        Pf = collect_Pl(M, u_free)
        J[row_p+1:row_p+M, 2M-1:3M-1] .= Pf
    else
        # R4 (all): ∂(Pall·B)/∂B = Pall
        u_all, _ = gauss_legendre_nodes(M+1, -1.0, 1.0)
        Pall = collect_Pl(M, u_all)
        J[row_p:row_p+M, 2M-1:3M-1] .= Pall
    end

    # ── R6 ──────────────────────────────────────────────────────────────────
    J[3M,   3M]   = ak
    J[3M,   3M+1] = -dt

    # ── R7 ──────────────────────────────────────────────────────────────────
    J[3M+1, 3M+1] = ak
    if contact
        z_prev = history[end].z
        Iv = integral_at_theta_star(θs, M)
        J[3M+1, 2M-1:3M-1] .= -dt .* (3/2) .* z_prev^2 .* Iv
    end

    return J
end
