using Test, LinearAlgebra
using DropSolver

# Exact Γ₂^(a) in the inviscid limit for several shape exponents a
# (shear_thinning_derivation.ipynb, Section 6, general-a Γ_l^(a) computation).
const GAMMA2_INVISCID_A = Dict(
    1.0 => 439.1737,
    1.5 => 1416.7890,
    2.0 => 1783566 / 385,   # exact rational, a=2 (standard Carreau)
    3.0 => 51233.2220,
)

@testset "Carreau-Yasuda: backward compatibility" begin
    @testset "3-arg STParams defaults to a=2.0 (standard Carreau)" begin
        st = STParams(0.1, 0.1, [GAMMA2_INVISCID_A[2.0]])
        @test st.a == 2.0
    end

    @testset "STParams() defaults to a=2.0" begin
        @test STParams().a == 2.0
    end

    @testset "4-arg constructor with a=2.0 is identical to the 3-arg (default) constructor" begin
        M = 4; Oh = 0.05; Bo = 1e-12
        theta_vec = collect(range(π, 0; length=M+1))
        precomp   = precompute_integrals(NaN, M)[1]
        cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, 0.005)
        ob        = OBParams()

        Gamma = fill(GAMMA2_INVISCID_A[2.0], M-1)
        st_3arg = STParams(0.1, 0.1, Gamma)
        st_4arg = STParams(0.1, 0.1, Gamma, 2.0)

        s_prev = DropState(M); s_prev.Adot[2] = 0.3; s_prev.z = 1.0; s_prev.dt = 0.005
        s_curr = DropState(M); s_curr.A[2] = 0.05; s_curr.Adot[2] = 0.25; s_curr.z = 1.0
        dt = 0.005

        R_3 = zeros(3M+1); build_residual_st!(R_3, s_curr, [s_prev], dt, 0, cfg, ob, st_3arg)
        R_4 = zeros(3M+1); build_residual_st!(R_4, s_curr, [s_prev], dt, 0, cfg, ob, st_4arg)
        @test R_3 ≈ R_4 atol=1e-14
    end
end

@testset "Carreau-Yasuda unit tests (general shape exponent a)" begin
    M = 4; Oh = 0.05; Bo = 1e-12
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, 0.005)
    ob        = OBParams()

    @testset "Newtonian limit: eps_ST=0 -> identical to build_residual!, for any a" begin
        for a_val in (1.0, 1.5, 3.0)
            st_zero = STParams(0.0, 0.0, zeros(M-1), a_val)
            s0 = DropState(M); s0.A[2] = 0.05; s0.z = 1.0; s0.dt = 0.005
            R_N  = zeros(3M+1); build_residual!(R_N, s0, [s0], 0.005, 0, cfg, ob)
            R_ST = zeros(3M+1); build_residual_st!(R_ST, s0, [s0], 0.005, 0, cfg, ob, st_zero)
            @test R_ST ≈ R_N atol=1e-14
        end
    end

    @testset "R2 correction formula (general a): matches hand-computed |Adot|^a expression" begin
        for a_val in (1.0, 1.5, 3.0)
            eps_ST = 0.02; lambda_c = 0.1
            Gamma  = fill(GAMMA2_INVISCID_A[a_val], M-1)
            st = STParams(eps_ST, lambda_c, Gamma, a_val)

            s_prev = DropState(M); s_prev.Adot[2] = 0.3; s_prev.z = 1.0; s_prev.dt = 0.005
            s_curr = DropState(M); s_curr.A[2] = 0.05; s_curr.Adot[2] = 0.25; s_curr.z = 1.0
            dt = 0.005

            R_N  = zeros(3M+1); build_residual!(R_N, s_curr, [s_prev], dt, 0, cfg, ob)
            R_ST = zeros(3M+1); build_residual_st!(R_ST, s_curr, [s_prev], dt, 0, cfg, ob, st)

            ns     = collect(Float64, 2:M)
            D2     = @. 2Oh * (ns - 1) * (2ns + 1)
            sigma0 = @. sqrt(ns * (ns - 1) * (ns + 2))
            Gamma_eff     = Gamma .* (lambda_c .* sigma0) .^ a_val
            shear_pow_lag = sum(Gamma_eff[k] * abs(s_prev.Adot[k+1])^a_val for k in 1:M-1)

            expected_correction = dt .* D2 .* s_curr.Adot[2:end] .* (eps_ST * shear_pow_lag)
            @test R_N[M:2M-2] - R_ST[M:2M-2] ≈ expected_correction atol=1e-14

            # Non-R2 blocks unchanged
            @test R_ST[1:M-1]    ≈ R_N[1:M-1]    atol=1e-14
            @test R_ST[2M-1:end] ≈ R_N[2M-1:end] atol=1e-14
        end
    end

    @testset "Jacobian finite-difference check (general a)" begin
        for a_val in (1.0, 1.5, 3.0)
            eps_ST = 0.02; lambda_c = 0.1
            st = STParams(eps_ST, lambda_c, fill(GAMMA2_INVISCID_A[a_val], M-1), a_val)

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
end

@testset "Carreau-Yasuda physics: shear-thinning reduces decay rate for a != 2" begin
    M = 2; Oh = 0.02; Bo = 1e-6
    l = 2
    omega_l = sqrt(Float64(l*(l-1)*(l+2)))

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams()

    a0 = 0.05
    T_period = 2π / omega_l
    t_end = 8 * T_period

    for shape_a in (1.5, 3.0)
        eps_ST = 0.02; lambda_c = 0.1
        st = STParams(eps_ST, lambda_c, [GAMMA2_INVISCID_A[shape_a]], shape_a)

        init_N  = DropState(M); init_N.A[2]  = a0; init_N.z = 2.0; init_N.dt = dt_max
        init_ST = DropState(M); init_ST.A[2] = a0; init_ST.z = 2.0; init_ST.dt = dt_max

        t_N,  s_N  = solve_drop!(cfg, ob, init_N;  t_end=t_end, save_every=T_period/50)
        t_ST, s_ST = solve_drop!(cfg, ob, init_ST; st=st, t_end=t_end, save_every=T_period/50)

        A2_N  = [s.A[2] for s in s_N]
        A2_ST = [s.A[2] for s in s_ST]

        γ_N  = -log(abs(A2_N[end])  / abs(A2_N[1]))  / (t_N[end]  - t_N[1])
        γ_ST = -log(abs(A2_ST[end]) / abs(A2_ST[1])) / (t_ST[end] - t_ST[1])

        @test γ_ST < γ_N
    end
end
