using Test, LinearAlgebra
using DropSolver

# Exact Γ₂ in the inviscid limit (symbolic result, shear_thinning_derivation.ipynb cell 30)
const GAMMA2_INVISCID = 1783566 / 385   # ≈ 4632.6

@testset "Carreau unit tests" begin
    M = 4; Oh = 0.05; Bo = 1e-12
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, 0.005)
    ob        = OBParams()

    @testset "Newtonian limit: eps_ST=0 → identical to build_residual!" begin
        st_zero = STParams(0.0, 0.0, zeros(M-1))
        s0 = DropState(M); s0.A[2] = 0.05; s0.z = 1.0; s0.dt = 0.005
        R_N  = zeros(3M+1); build_residual!(R_N, s0, [s0], 0.005, 0, cfg, ob)
        R_ST = zeros(3M+1); build_residual_st!(R_ST, s0, [s0], 0.005, 0, cfg, ob, st_zero)
        @test R_ST ≈ R_N atol=1e-14
    end

    @testset "R2 correction formula: shear-thinning reduces Adot term" begin
        eps_ST = 0.1; lambda_c = 0.1
        Gamma  = fill(GAMMA2_INVISCID, M-1)
        st = STParams(eps_ST, lambda_c, Gamma)

        s_prev = DropState(M); s_prev.Adot[2] = 0.3; s_prev.z = 1.0; s_prev.dt = 0.005
        s_curr = DropState(M); s_curr.A[2] = 0.05; s_curr.Adot[2] = 0.25; s_curr.z = 1.0

        dt = 0.005
        R_N  = zeros(3M+1); build_residual!(R_N, s_curr, [s_prev], dt, 0, cfg, ob)
        R_ST = zeros(3M+1); build_residual_st!(R_ST, s_curr, [s_prev], dt, 0, cfg, ob, st)

        ns     = collect(Float64, 2:M)
        D2     = @. 2Oh * (ns - 1) * (2ns + 1)
        sigma0 = @. sqrt(ns * (ns - 1) * (ns + 2))
        Gamma_eff    = Gamma .* (lambda_c .* sigma0) .^ 2
        shear_sq_lag = sum(Gamma_eff[k] * s_prev.Adot[k+1]^2 for k in 1:M-1)

        expected_correction = dt .* D2 .* s_curr.Adot[2:end] .* (eps_ST * shear_sq_lag)
        @test R_N[M:2M-2] - R_ST[M:2M-2] ≈ expected_correction atol=1e-14

        # Non-R2 blocks unchanged
        @test R_ST[1:M-1]    ≈ R_N[1:M-1]    atol=1e-14
        @test R_ST[2M-1:end] ≈ R_N[2M-1:end] atol=1e-14
    end

    @testset "Jacobian finite-difference check" begin
        eps_ST = 0.1; lambda_c = 0.1
        st = STParams(eps_ST, lambda_c, fill(GAMMA2_INVISCID, M-1))

        s_prev = DropState(M); s_prev.Adot[2] = 0.3; s_prev.z = 1.0; s_prev.dt = 0.005
        s_curr = DropState(M); s_curr.A[2] = 0.05; s_curr.Adot[2] = 0.25; s_curr.z = 1.0

        dt = 0.005
        J   = build_jacobian_st(s_curr, [s_prev], dt, 0, cfg, ob, st)
        X0  = pack_X(s_curr, M)
        R0  = zeros(3M+1); build_residual_st!(R0, s_curr, [s_prev], dt, 0, cfg, ob, st)
        eps_fd = 1e-6
        J_fd = zeros(3M+1, 3M+1)
        for j in 1:3M+1
            Xp = copy(X0); Xp[j] += eps_fd
            sp = deepcopy(s_curr); unpack_X!(sp, Xp, M)
            Rp = zeros(3M+1); build_residual_st!(Rp, sp, [s_prev], dt, 0, cfg, ob, st)
            J_fd[:, j] = (Rp .- R0) ./ eps_fd
        end
        @test norm(J .- J_fd) / (norm(J_fd) + 1e-14) < 1e-4
    end
end

@testset "Carreau physics: shear-thinning reduces decay rate" begin
    M = 2; Oh = 0.02; Bo = 1e-6
    l = 2
    omega_l = sqrt(Float64(l*(l-1)*(l+2)))
    gamma_0 = Float64((l-1)*(2l+1)) * Oh

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams()

    a0 = 0.05; eps_ST = 0.1; lambda_c = 0.1
    st = STParams(eps_ST, lambda_c, [GAMMA2_INVISCID])

    init_N  = DropState(M); init_N.A[2]  = a0; init_N.z = 2.0; init_N.dt = dt_max
    init_ST = DropState(M); init_ST.A[2] = a0; init_ST.z = 2.0; init_ST.dt = dt_max

    T_period = 2π / omega_l
    t_end = 8 * T_period

    t_N,  st_N  = solve_drop!(cfg, ob, init_N;  t_end=t_end, save_every=T_period/50)
    t_ST, st_ST = solve_drop!(cfg, ob, init_ST; st=st, t_end=t_end, save_every=T_period/50)

    A2_N  = [s.A[2] for s in st_N]
    A2_ST = [s.A[2] for s in st_ST]

    γ_N  = -log(abs(A2_N[end])  / abs(A2_N[1]))  / (t_N[end]  - t_N[1])
    γ_ST = -log(abs(A2_ST[end]) / abs(A2_ST[1])) / (t_ST[end] - t_ST[1])

    @test γ_ST < γ_N
    @test γ_N - γ_ST > 0.001
end

@testset "Carreau quantitative: decay rate monotone in eps_ST" begin
    # The slow amplitude equation predicts da/dt = -γ₀a[1 - (3/4)ε_ST Λ²Γ_l a²ω²].
    # Regardless of whether the perturbation regime holds, larger eps_ST should
    # always yield a smaller effective decay rate (more shear-thinning = less damping).
    # This test verifies monotonicity: γ(eps_ST=0) > γ(eps_ST=0.05) > γ(eps_ST=0.10).
    M = 2; Oh = 0.02; Bo = 1e-6
    omega_l = sqrt(Float64(2*(2-1)*(2+2)))

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams()

    a0 = 0.05; lambda_c = 0.1
    T_period = 2π / omega_l
    t_end = 8 * T_period

    function run_decay(eps_ST_val)
        st_run = STParams(eps_ST_val, lambda_c, [GAMMA2_INVISCID])
        init   = DropState(M); init.A[2] = a0; init.z = 2.0; init.dt = dt_max
        t_r, s_r = solve_drop!(cfg, ob, init; st=st_run, t_end=t_end, save_every=T_period/10)
        A2 = [s.A[2] for s in s_r]
        -log(abs(A2[end]) / abs(A2[1])) / (t_r[end] - t_r[1])
    end

    γ0   = run_decay(0.00)
    γ05  = run_decay(0.05)
    γ10  = run_decay(0.10)

    # More shear-thinning → lower effective decay rate
    @test γ0 > γ05
    @test γ05 > γ10
    # The effect is non-negligible (at least 0.5% between steps)
    @test γ0 - γ05 > 0.0005
    @test γ05 - γ10 > 0.0005
end
