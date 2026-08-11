@testset "Contact" begin
    using DropSolver: contact_error, drop_height, DropState

    @testset "Contact detection" begin
        M = 6
        s = DropState(M)
        θv = collect(range(π, π/2; length=M+1))

        @testset "cp=0 always returns 0" begin
            s.z = 1.5
            @test contact_error(s, θv, 0) == 0.0
            s.z = -2.0
            @test contact_error(s, θv, 0) == 0.0
        end

        @testset "Penetration → Inf" begin
            s.z = -2.0   # drop center far below substrate
            @test contact_error(s, θv, 1) == Inf
        end

        @testset "Drop above substrate, cp=1 → finite error" begin
            s2 = DropState(M)
            s2.z = 1.0   # COM well above substrate
            err = contact_error(s2, θv, 1)
            @test isfinite(err)
        end
    end
end
