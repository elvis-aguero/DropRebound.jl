using Test, LinearAlgebra
using DropSolver

function make_eq_state(M)
    s = DropState(M)
    s.z = 1.0
    s.v = 0.0
    s.cp = 0
    return s
end

@testset "Residual (Newtonian, no contact)" begin
    M = 4; Oh = 0.1; Fr = 1e12
    θv = collect(range(π, π/2 + 0.01; length = M+1))
    precomp = precompute_integrals(NaN, M)[1]
    cfg = SimConstants(M, M+1, Oh, Fr, θv, precomp, 0.01)
    ob  = OBParams()

    @testset "Equilibrium BDF1 → zero residual" begin
        s0 = make_eq_state(M); s0.dt = 0.01
        R  = zeros(3M + 1)
        build_residual!(R, s0, [s0], 0.01, 0, cfg, ob)
        @test norm(R) < 1e-13
    end

    @testset "Nonzero A₂ → nonzero capillary term in R2" begin
        s0 = make_eq_state(M); s0.dt = 0.01
        s1 = make_eq_state(M); s1.dt = 0.01
        s1.A[2] = 0.1
        R = zeros(3M + 1)
        build_residual!(R, s1, [s0], 0.01, 0, cfg, ob)
        # R2 for n=2: dt * n(n+2)(n-1) * A₂ = 0.01 * 2*4*1 * 0.1 = 0.008
        @test abs(R[M] - 0.008) < 1e-10
    end

    @testset "pack/unpack round-trip" begin
        s = make_eq_state(M)
        s.A[2] = 0.5; s.Adot[3] = -0.3; s.B[1] = 1.2; s.z = 2.0; s.v = -0.1
        X = pack_X(s, M)
        s2 = DropState(M)
        unpack_X!(s2, X, M)
        @test s2.A[2:end]    ≈ s.A[2:end]
        @test s2.Adot[2:end] ≈ s.Adot[2:end]
        @test s2.B           ≈ s.B
        @test s2.z           ≈ s.z
        @test s2.v           ≈ s.v
    end
end
