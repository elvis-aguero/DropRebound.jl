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
# These tests pin the MONOMIAL window, which is why `:legendre` is now the default: same
# span, same K=1 potential-flow limit, conditioning better by up to ten orders of magnitude at
# low l. They are kept because the monomial basis is still selectable and because the ceiling
# is the reason the default changed.
@testset "the radial basis conditioning ceiling (:monomial, the historical basis)" begin
    @testset "conditioning grows geometrically in K, faster at high l" begin
        for l in (2, 20, 90)
            cs = [cond(assemble(RitzBasis(l, K, :monomial), 0.05).M) for K in 1:6]
            @test issorted(cs)                      # monotone
            @test cs[1] ≈ 1.0 atol=1e-8             # one function is perfectly conditioned
            ## each added function costs at least a factor of ten
            @test all(cs[k+1] / cs[k] > 10 for k in 1:5)
        end
        ## and the growth is worse at high l, which is why large M reaches the wall first
        c20 = cond(assemble(RitzBasis(20, 5, :monomial), 0.05).M)
        c90 = cond(assemble(RitzBasis(90, 5, :monomial), 0.05).M)
        @test c90 > c20
    end

    @testset "the usable window, measured" begin
        # Where the assembly still carries information. Beyond this the answers are not
        # merely inaccurate, they are arbitrary.
        @test cond(assemble(RitzBasis(45, 5, :monomial), 0.05).M) < 1e14   # usable
        @test cond(assemble(RitzBasis(45, 6, :monomial), 0.05).M) > 1e15   # not
        @test cond(assemble(RitzBasis(90, 4, :monomial), 0.05).M) < 1e13   # usable
        @test cond(assemble(RitzBasis(90, 5, :monomial), 0.05).M) > 1e14   # not
        ## the reconditioned basis clears the same points by orders of magnitude, which is
        ## the whole reason it is the default
        @test cond(assemble(RitzBasis(45, 6, :legendre), 0.05).M) <
              cond(assemble(RitzBasis(45, 6, :monomial), 0.05).M) / 100
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

# The empirical radial window, K(l).
#
# Rather than recomputing the conditioning ceiling for every mode at every run, it is fitted
# once from data and encoded. These tests are what makes that safe: they check the law against
# the thing it was fitted to, so a change to the basis that moves the ceiling fails here
# instead of silently producing a singular block.
@testset "the radial window law" begin
    @testset "the law never exceeds the measured ceiling" begin
        # Measured on the :legendre basis at cond < 1e10. The law must sit at or below this
        # at every mode -- overshooting means proposing a K that cannot be solved.
        measured = Dict(2=>40, 3=>29, 4=>21, 5=>16, 6=>14, 8=>11, 10=>9, 12=>8,
                        16=>7, 20=>6, 25=>5, 30=>5, 40=>4, 50=>4, 60=>4, 75=>4, 90=>3)
        for (l, kmax) in measured
            @test radial_window(l) <= kmax
        end
        ## and it must not be uselessly conservative either
        @test radial_window(2)  >= 20
        @test radial_window(8)  >= 8
        @test radial_window(90) >= 2
    end

    @testset "the law actually caps the conditioning" begin
        # The point of the law, checked directly rather than through the fit.
        nrmcond(A) = (d = sqrt.(abs.(diag(A))); cond(A ./ (d*d')))
        for l in (2, 5, 10, 20, 45, 90)
            K = radial_window(l)
            @test K >= 1
            c = nrmcond(assemble(RitzBasis(l, K, :legendre), 0.05).M)
            @test c < 1e11                      # the fitted threshold, with slack for the fit
        end
    end

    @testset "it is monotone and bounded" begin
        ws = [radial_window(l) for l in 2:90]
        @test issorted(ws; rev = true)          # never grants a higher mode more functions
        @test all(w -> 1 <= w <= 40, ws)
        ## the exponent, read back off the law itself
        @test radial_window(4)/radial_window(2) < 0.8      # decays, and not trivially
        @test radial_window(90) < radial_window(8) < radial_window(2)
    end
end

# The cache key must include the basis kind.
#
# It did not, and the consequence was worse than a wrong number: a run with `:monomial`
# populated the cache, a subsequent `:legendre` run reused that geometry, and the two bases
# returned IDENTICAL results. That looks exactly like a confirmation that they span the same
# space -- which they do, so the false positive was indistinguishable from the thing being
# checked. It also made the reconditioned basis appear to fail at K = 4 exactly where the
# monomial one does, which reversed the conclusion of a comparison against experiment.
@testset "the geometry cache distinguishes basis kinds" begin
    ls = collect(2:8); K = 3
    gm = DropSolver.coupled_geometry(ModalBasis(ls, K, :monomial), 40, 48)
    gl = DropSolver.coupled_geometry(ModalBasis(ls, K, :legendre), 40, 48)
    ## different bases, genuinely different geometry
    @test gm.D !== gl.D
    @test maximum(abs, gm.D .- gl.D) > 1e-6 * maximum(abs, gm.D)
    @test maximum(abs, gm.M .- gl.M) > 1e-6 * maximum(abs, gm.M)
    ## and asking twice for the same one returns the cached object, not a rebuild
    @test DropSolver.coupled_geometry(ModalBasis(ls, K, :legendre), 40, 48) === gl

    ## the property the bug destroyed: same span, so the SPECTRUM agrees even though the
    ## matrices do not. Checked after the geometry has been built both ways, which is the
    ## order that used to give a false pass.
    for l in (2, 8), Kk in (1, 2, 3)
        a = DropSolver.dominant_pair(RitzBasis(l, Kk, :monomial), 0.05)[1]
        b = DropSolver.dominant_pair(RitzBasis(l, Kk, :legendre), 0.05)[1]
        @test isapprox(a, b; rtol = 1e-8)
    end
end

# The basis kind must reach every assembly path.
#
# It did not. `assemble_newtonian` built `RitzBasis(l, b.K)` without `b.kind`, so every
# Newtonian run used the default basis whatever was asked for -- and a comparison of the two
# bases then compared one basis with itself, agreed perfectly, and read as confirmation. That
# is the third place today the same omission appeared (the geometry cache key, the results
# store key, and here), and in all three the symptom was two configurations returning identical
# numbers, which is exactly what the correctness property predicts.
#
# So the thing to assert first is that the arms DIFFER. Agreement is only evidence once that is
# established.
@testset "the basis kind reaches the assembly" begin
    b_mono = ModalBasis(collect(2:20), 6, :monomial)
    b_leg  = ModalBasis(collect(2:20), 6, :legendre)
    Fm = assemble_newtonian(b_mono, 0.05)
    Fl = assemble_newtonian(b_leg, 0.05)
    ## the arms are genuinely different matrices
    @test maximum(abs, Fm.M .- Fl.M) > 1e-6 * maximum(abs, Fm.M)
    @test cond(Fm.M) > 100 * cond(Fl.M)
    ## and only then is their agreement on the physics meaningful
    for l in (2, 8, 20), K in (1, 2, 3)
        a = DropSolver.dominant_pair(RitzBasis(l, K, :monomial), 0.05)[1]
        c = DropSolver.dominant_pair(RitzBasis(l, K, :legendre), 0.05)[1]
        @test isapprox(a, c; rtol = 1e-7)
    end
    ## The failure this fixed: at M = 45, K = 6 the monomial basis manufactures a surface
    ## amplitude of 23 -- a drop twenty-three times its own radius -- and a restitution of 102.
    ##
    ## Detected on the AMPLITUDE after a SHORT march rather than the restitution after a full
    ## one. The divergence is fully present by t = 2, so marching to t_max costs four times as
    ## much and establishes nothing further. Written the expensive way this testset ran for
    ## nearly four minutes, which is not a price a unit test should charge.
    short = (We = 0.5, Bo = 0.0188, Oh = 0.0373, M = 45, K = 6, t_max = 2.0)
    amp(p) = (r = simulate_lcp(p);
              maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a))
    am = try amp(ImpactParams(; short..., basis_kind = :monomial)) catch; Inf end
    al = try amp(ImpactParams(; short..., basis_kind = :legendre)) catch; Inf end
    @test am > 5.0                     # monomials: nonphysical, measured 23
    @test al < 1.0                     # reconditioned: measured 0.436
    @test al < am / 20
end
