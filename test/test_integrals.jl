using Test
using DropSolver: precompute_integrals

@testset "Precomputed integrals" begin
    angles = [π, 3π/4]
    N = 3
    M_mat, ang_out = precompute_integrals(angles, N)
    # Single interval [π, 3π/4]: one row, N+1 columns
    @test size(M_mat, 1) == 1
    @test size(M_mat, 2) == N + 1

    @testset "P₀ integral = ∫du/u³ analytically" begin
        # ∫_{cos(π)}^{cos(3π/4)} u⁻³ du = [-1/(2u²)]_{-1}^{-√2/2}
        a = cos(π)       # = -1.0
        b = cos(3π/4)    # = -√2/2 ≈ -0.7071
        expected_P0 = -1/(2*b^2) + 1/(2*a^2)
        @test abs(M_mat[1, 1] - expected_P0) < 1e-4
    end

    @testset "NaN input gives shape" begin
        M2, _ = precompute_integrals(NaN, 5)
        @test size(M2, 2) == 6   # modes 0..5
    end
end
