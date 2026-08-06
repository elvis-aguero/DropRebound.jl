# The cached coupled assembly, and the limits of the radial basis.
#
# With a variable viscosity the dissipation operator is rebuilt on every Picard sweep, and
# profiling put the entire cost of a shear-thinning run there -- 101 ms per sweep against
# 0.1 ms for the whole contact solve. Almost none of it depended on the state: `strain_at` is
# geometry, the mass matrix carries no viscosity at all, and only the scalar `eta` at each
# quadrature point changes between sweeps. Caching both the pairwise contractions and the
# strain basis itself took a full shear-thinning impact from 145 s to 1.4 s.
#
# A cache that is 100x faster and 1e-14 different is not a cache, it is a second
# implementation. So the first thing checked here is that it is EXACT, and exact in the strong
# sense: the same floating-point operations in the same order, not merely agreement to
# tolerance. That was verified in a scratch script when the cache was written, which is the
# same pattern that let an artefact survive weeks of numerical checking, so it lives here now.

using Test
using DropSolver
using LinearAlgebra
using Random: Xoshiro

"""The pre-cache assembly, kept verbatim as the reference the cache must reproduce."""
function reference_assembly(b::ModalBasis, Oh::Real, eta, nx::Int, nmu::Int)
    N = DropSolver.ndof(b)
    xs, wxs = DropSolver.gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = DropSolver.gauss_legendre_nodes(nmu, -1.0, 1.0)
    Mm, Cm = zeros(N, N), zeros(N, N)
    for (x, wx) in zip(xs, wxs), (mu, wmu) in zip(mus, wmus)
        F = DropSolver.strain_at(b, x, mu)
        w = wx * wmu * 2pi * x^2
        ev = eta(x, mu)
        for p in 1:N, q in p:N
            Mm[p, q] += w * (F[p,1]*F[q,1] + F[p,2]*F[q,2])
            Cm[p, q] += w * 2 * ev * DropSolver.ddot_strain(view(F, p, :), view(F, q, :))
        end
    end
    for p in 1:N, q in 1:p-1
        Mm[p, q] = Mm[q, p]; Cm[p, q] = Cm[q, p]
    end
    (M = Mm, C = Cm .* Oh)
end

@testset "cached coupled assembly" begin

    @testset "the cache is bitwise identical to the direct loop" begin
        # Not `isapprox`. If the cached path ever differs by so much as an ulp, the two are
        # no longer the same computation and the reason needs finding rather than tolerating.
        for (ls, K) in ((2:8, 2), (2:14, 2), (2:20, 1), (2:6, 3))
            b = ModalBasis(collect(ls), K)
            for eta in ((x, mu) -> 1.0,                         # constant: decouples
                        (x, mu) -> 0.6 + 0.3*mu,                # angular: couples modes
                        (x, mu) -> 0.4 + 0.5*x^2,               # radial only
                        (x, mu) -> 0.4 + 0.5*x^2*(1 - 0.3*mu))  # both
                ref = reference_assembly(b, 0.3, eta, 40, 48)
                got = assemble_coupled(b, 0.3; eta = eta, nx = 40, nmu = 48)
                @test got.M == ref.M
                @test got.C == ref.C
            end
        end
    end

    @testset "the shear-rate fast path is the same computation" begin
        # `eta_rate` + `state` evaluates the viscosity from the cached strain basis instead of
        # rebuilding it. That is the half of the assembly the first version of this cache
        # missed: the contractions were cached and the closure rebuilt `strain_at` on the way
        # in, so a shear-thinning run saw almost no speedup.
        for (ls, K) in ((2:14, 2), (2:8, 3))
            b = ModalBasis(collect(ls), K); N = DropSolver.ndof(b)
            a = zeros(N); a[1] = 0.35; a[min(4,N)] = -0.12; a[min(7,N)] = 0.08
            visc(gd) = carreau(gd; lambda_c = 12.0, a = 2.0, n = 0.5, eta_inf_ratio = 0.01)
            slow = assemble_coupled(b, 0.3; eta = (x, mu) -> visc(shear_rate(b, a, x, mu)))
            fast = assemble_coupled(b, 0.3; eta_rate = visc, state = a)
            @test fast.C == slow.C
            @test fast.M == slow.M
        end
    end

    @testset "shear_rate_at agrees with shear_rate" begin
        b = ModalBasis(collect(2:10), 2); N = DropSolver.ndof(b)
        a = randn(Xoshiro(1), N) .* 0.1
        g = DropSolver.coupled_geometry(b, 40, 48)
        for q in (1, 17, 191, 900, length(g.w))
            @test DropSolver.shear_rate_at(g, a, q) ≈
                  shear_rate(b, a, g.xs[q], g.mus[q]) rtol=1e-12
        end
    end

    @testset "the memory budget is respected" begin
        # The footprint grows as ndof^2, so it has to be bounded or a large truncation would
        # silently allocate gigabytes. Above the budget the direct loop is used instead --
        # slower, but it runs.
        @test DropSolver.coupled_cache_bytes(26, 40, 48) < 10_000_000
        @test DropSolver.coupled_cache_bytes(178, 40, 48) < DropSolver.COUPLED_CACHE_BUDGET
        @test DropSolver.coupled_cache_bytes(258, 40, 48) > DropSolver.COUPLED_CACHE_BUDGET
        ## monotone in every argument, since a bound that is not monotone is not a bound
        @test DropSolver.coupled_cache_bytes(50, 40, 48) >
              DropSolver.coupled_cache_bytes(40, 40, 48)
        @test DropSolver.coupled_cache_bytes(50, 60, 48) >
              DropSolver.coupled_cache_bytes(50, 40, 48)
    end
end

# The radial basis has a ceiling, and it is the basis rather than the physics.
#
# The Ritz functions are x^(l+1), x^(l+3), ... -- the Taylor terms of j_l(qx) about x = 0. The
# span is right and K = 1 is potential flow exactly, but as a BASIS it is Vandermonde: the
# per-mode mass matrix loses roughly one and a half digits per added function at l = 2 and
# two and a half at l = 90. Double precision then runs out, and the failure is not graceful --
# at M = 45, K = 6 the solver returns a restitution of 102, and on the shear-thinning fluid
# K >= 4 returns NaN or a plausible-looking 0.998.
#
# These tests pin the window rather than the failure, so that raising K past what the basis
# can carry is caught here instead of in a results table.
@testset "the radial basis conditioning ceiling" begin
    @testset "conditioning grows geometrically in K, faster at high l" begin
        for l in (2, 20, 90)
            cs = [cond(assemble_newtonian(ModalBasis([l], K), 0.05).M) for K in 1:6]
            @test issorted(cs)                      # monotone
            @test cs[1] ≈ 1.0 atol=1e-8             # one function is perfectly conditioned
            ## each added function costs at least a factor of ten
            @test all(cs[k+1] / cs[k] > 10 for k in 1:5)
        end
        ## and the growth is worse at high l, which is why large M reaches the wall first
        c20 = cond(assemble_newtonian(ModalBasis([20], 5), 0.05).M)
        c90 = cond(assemble_newtonian(ModalBasis([90], 5), 0.05).M)
        @test c90 > c20
    end

    @testset "the usable window, measured" begin
        # Where the assembly still carries information. Beyond this the answers are not
        # merely inaccurate, they are arbitrary.
        @test cond(assemble_newtonian(ModalBasis([45], 5), 0.05).M) < 1e14   # usable
        @test cond(assemble_newtonian(ModalBasis([45], 6), 0.05).M) > 1e15   # not
        @test cond(assemble_newtonian(ModalBasis([90], 4), 0.05).M) < 1e13   # usable
        @test cond(assemble_newtonian(ModalBasis([90], 5), 0.05).M) > 1e14   # not
    end

    @testset "the span is fine; it is the basis that is ill-conditioned" begin
        # Orthogonalising the same span in the mass inner product, keeping x^(l+1) first so
        # the potential-flow limit is untouched, drops the conditioning by ten orders of
        # magnitude and leaves the damping unchanged. That is the evidence that the ceiling is
        # a basis choice rather than a limit of the model, and it is what a replacement radial
        # family would have to preserve.
        for l in (2, 20)
            K = 6
            F = assemble_newtonian(ModalBasis([l], K), 0.05)
            T = Matrix{Float64}(I, K, K)
            for j in 2:K
                v = zeros(K); v[j] = 1.0
                for i in 1:j-1
                    u = T[:, i]; v -= (dot(u, F.M, v)/dot(u, F.M, u)) * u
                end
                T[:, j] = v / sqrt(dot(v, F.M, v))
            end
            Mo = T' * F.M * T
            @test cond(Mo) < 1e4                       # measured 1.26 at l=2, 65 at l=20
            @test cond(Mo) < cond(F.M) / 1e3           # and vastly better than the monomials
        end
    end
end
