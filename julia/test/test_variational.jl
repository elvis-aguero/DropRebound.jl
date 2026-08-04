# The variational assembly, against the Newtonian theory it must reproduce.
#
# This is the first thing the solver has to get right: with a constant viscosity,
# the variational statement of the model -- three quadratic forms and the
# Euler-Lagrange equations, with the interior retained as part of the state -- must
# return Reid's exact coefficients. `reid_lambda_omega2` computes them by an
# entirely separate route (Newton on the characteristic equation, via a Bessel
# ratio), so agreement is evidence about the assembly and not about itself.
#
# It also ties the one loose end the derivation page could not: that the surface
# equation obtained by varying the surface energy agrees with the traction balance.
# If it did not, these numbers would not match.

using Test
using DropSolver
using LinearAlgebra

@testset "variational assembly" begin

    @testset "the forms are Hessians, so symmetric" begin
        for (l, K, Oh) in ((2, 4, 0.05), (5, 3, 0.3), (8, 5, 0.006))
            F = assemble(RitzBasis(l, K), Oh)
            @test F.M ≈ F.M'
            @test F.C ≈ F.C'
            @test isposdef(Symmetric(F.M))          # a kinetic energy
            @test all(eigvals(Symmetric(F.C)) .> -1e-12)  # a dissipation rate
        end
    end

    @testset "the surface energy is rank one" begin
        # It reaches the interior only through the boundary trace, so it cannot
        # resolve interior structure. A higher rank would mean the stiffness had
        # acquired a dependence on the profile shape, which surface tension has no
        # way to supply.
        for K in 2:6
            @test rank(assemble(RitzBasis(3, K), 0.1).G) == 1
        end
    end

    @testset "one trial function is the inviscid limit" begin
        # phi_1 = x^(l+1) is potential flow, so K = 1 must give Rayleigh's
        # frequency exactly and Lamb's damping exactly -- the two classical results
        # the variational route reproduces without a Bessel function in sight.
        for l in 2:6
            F = assemble(RitzBasis(l, 1), 1.0)
            ω² = F.G[1, 1] / F.M[1, 1]
            @test ω² ≈ l * (l - 1) * (l + 2)  rtol=1e-10
            λ_per_Oh = F.C[1, 1] / (2 * F.M[1, 1])
            @test λ_per_Oh ≈ (l - 1) * (2l + 1)  rtol=1e-10
        end
    end

    @testset "it reproduces Reid at finite Ohnesorge" begin
        # Across three decades of Oh and both dynamical regimes. At Oh = 1, l = 2 the
        # root pair has gone REAL (above the critical Ohnesorge), which is why the
        # coefficients are Vieta on the two slowest roots rather than on a conjugate
        # pair -- reading them as a conjugate pair is wrong by 64% there.
        for (Oh, l, K, tol) in ((0.006, 2, 8, 1e-5),
                                (0.05,  2, 8, 1e-10),
                                (0.3,   2, 8, 1e-10),
                                (1.0,   2, 8, 1e-10),   # overdamped: real roots
                                (0.05,  4, 8, 1e-6),
                                (0.05,  8, 8, 1e-5),
                                (0.3,   6, 8, 1e-7))
            λ, ω² = dominant_pair(RitzBasis(l, K), Oh)
            rλ, rω², _ = reid_lambda_omega2(Oh, l)
            @test λ  ≈ rλ   rtol=tol
            @test ω² ≈ rω²  rtol=tol
        end
    end

    @testset "refining the basis converges" begin
        # Quadratic convergence is the property that makes a handful of trial
        # functions worth as much as a fine radial grid. Check the error actually
        # falls rather than just being small at one K.
        Oh, l = 0.05, 2
        rλ, _, _ = reid_lambda_omega2(Oh, l)
        errs = [abs(dominant_pair(RitzBasis(l, K), Oh)[1] - rλ) / rλ for K in 2:6]
        @test errs[end] < errs[1] / 1e4
        @test issorted(errs; rev=true) || errs[end] < 1e-9
    end

    @testset "the basis has a conditioning limit, and it is visible" begin
        # The monomial trial functions are Vandermonde-like: cond(M) grows by roughly
        # three decades per two functions added. Whitening by M's Cholesky factor
        # buys usable room, but the limit is real and a caller must be able to see
        # it rather than discover it as a silently wrong answer -- which is what
        # happened before the whitening, with omega^2 collapsing to 1e-9 at K = 8.
        cs = [decay_rates(RitzBasis(2, K), 0.05).cond_M for K in 2:8]
        @test issorted(cs)
        @test cs[1] < 1e3
        @test cs[end] > 1e9      # so an unbounded K is not safe, by construction
    end
end
