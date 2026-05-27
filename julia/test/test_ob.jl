using Test, LinearAlgebra
using DropSolver

@testset "Oldroyd-B extension" begin
    M = 4; Oh = 0.05; Fr = 1e12   # gravity negligible

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, 0.005)
    ob_newtonian = OBParams(0.0, 1.0)
    ob_ob        = OBParams(0.5, 0.5)

    @testset "OB state vector length is 4M" begin
        s = DropState(M)
        X = pack_X_ob(s, M)
        @test length(X) == 4M
    end

    @testset "OB residual, zero state → zero" begin
        s0 = DropState(M); s0.z = 1.0; s0.dt = 0.005
        R  = zeros(4M)
        build_residual_ob!(R, s0, [s0], 0.005, 0, cfg, ob_ob)
        @test norm(R) < 1e-12
    end

    @testset "Newtonian OB ≡ Newtonian (β_s=1, De₁=0)" begin
        s0 = DropState(M); s0.A[2] = 0.05; s0.z = 1.0; s0.dt = 0.005
        R_N  = zeros(3M+1)
        build_residual!(R_N, s0, [s0], 0.005, 0, cfg, ob_newtonian)
        R_OB = zeros(4M)
        build_residual_ob!(R_OB, s0, [s0], 0.005, 0, cfg, ob_newtonian)
        @test R_OB[1:3M+1] ≈ R_N  atol=1e-14
        @test norm(R_OB[3M+2:end]) < 1e-14
    end

    @testset "OB Jacobian matches finite differences" begin
        s0 = DropState(M); s0.A[2] = 0.02; s0.z = 1.0; s0.dt = 0.005
        J   = build_jacobian_ob(s0, [s0], 0.005, 0, cfg, ob_ob)
        X0  = pack_X_ob(s0, M)
        R0  = zeros(4M); build_residual_ob!(R0, s0, [s0], 0.005, 0, cfg, ob_ob)
        eps_fd = 1e-6
        J_fd = zeros(4M, 4M)
        for j in 1:4M
            Xp = copy(X0); Xp[j] += eps_fd
            sp = deepcopy(s0); unpack_X_ob!(sp, Xp, M)
            Rp = zeros(4M); build_residual_ob!(Rp, sp, [s0], 0.005, 0, cfg, ob_ob)
            J_fd[:, j] = (Rp .- R0) ./ eps_fd
        end
        @test norm(J .- J_fd) / (norm(J_fd) + 1e-14) < 1e-4
    end
end

@testset "OB l=2: decay rate decreases with De₁" begin
    M  = 6; Oh = 0.02; Fr = 1e6
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.A[2] = 0.05; init.z = 2.0; init.v = 0.0
    init.dt = dt_max; init.cp = 0

    t_end = 30.0

    # Newtonian run
    ob_N    = OBParams(0.0, 1.0)
    t_N, st_N = solve_drop!(cfg, ob_N, deepcopy(init); t_end=t_end, save_every=0.1)
    A2_N    = [s.A[2] for s in st_N]
    γ_N     = -log(abs(A2_N[end]) / abs(A2_N[1])) / (t_N[end] - t_N[1])

    # OB run with De₁=0.5, β_s=0.5
    ob_OB    = OBParams(0.5, 0.5)
    t_OB, st_OB = solve_drop!(cfg, ob_OB, deepcopy(init); t_end=t_end, save_every=0.1)
    A2_OB    = [s.A[2] for s in st_OB]
    γ_OB     = -log(abs(A2_OB[end]) / abs(A2_OB[1])) / (t_OB[end] - t_OB[1])

    @test γ_OB < γ_N    # viscoelastic drop rings longer
end
