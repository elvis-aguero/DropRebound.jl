using Test
using DropSolver: collect_Pl, collect_dPl

@testset "Legendre polynomials" begin
    @testset "P₀, P₁, P₂ at specific points" begin
        x = [0.0, 0.5, 1.0]
        P = collect_Pl(4, x)   # returns matrix (length(x) × N+1), 1-indexed: P[:,n+1] = Pₙ
        # P₀ = 1 everywhere
        @test P[:, 1] ≈ ones(3)
        # P₁ = x
        @test P[:, 2] ≈ x
        # P₂(x) = (3x²-1)/2
        @test P[:, 3] ≈ @. (3x^2 - 1) / 2
        # P₂(1) = 1
        @test P[3, 3] ≈ 1.0
    end

    @testset "Recurrence consistency" begin
        x = collect(range(-1, 1, 21))
        P = collect_Pl(10, x)
        # (n+1)Pₙ₊₁ = (2n+1)x Pₙ - n Pₙ₋₁
        for n in 1:9
            @test (n+1) .* P[:, n+2] ≈ (2n+1) .* x .* P[:, n+1] - n .* P[:, n]
        end
    end

    @testset "Orthogonality ∫₋₁¹ Pₙ Pₘ dx = 2/(2n+1) δₙₘ" begin
        N = 6
        x = collect(range(-1.0, 1.0, 200))
        dx = x[2] - x[1]
        P = collect_Pl(N, x)
        for n in 0:N-1, m in 0:N-1
            integral = sum(P[:, n+1] .* P[:, m+1]) * dx
            expected = n == m ? 2.0/(2n+1) : 0.0
            @test abs(integral - expected) < 0.02   # coarse quadrature tolerance
        end
    end

    @testset "Derivatives P'ₙ at specific points" begin
        x = [0.0, 0.5]
        dP = collect_dPl(4, x)
        # P'₁ = 1
        @test dP[:, 2] ≈ ones(2)
        # P'₂(x) = 3x
        @test dP[:, 3] ≈ 3 .* x
        # P'₃(x) = (15x²-3)/2
        @test dP[:, 4] ≈ @. (15x^2 - 3) / 2
    end

    @testset "Derivative recurrence (1-x²)P'ₙ = n(Pₙ₋₁ - xPₙ)" begin
        x = collect(range(-0.9, 0.9, 15))
        P  = collect_Pl(6, x)
        dP = collect_dPl(6, x)
        for n in 1:5
            lhs = @. (1 - x^2) * dP[:, n+1]
            rhs = n .* (P[:, n] .- x .* P[:, n+1])
            @test lhs ≈ rhs  atol=1e-12
        end
    end
end
