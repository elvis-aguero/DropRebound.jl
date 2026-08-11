# The conjugate (nodal-load) contact forcing, `force_mode = :nodal`.
#
# The production forcing carries the film pressure as a degree-M Legendre field and loads the
# drop with `Q_l = -(4 pi/(2l+1)) p_{c,l}`. That is the exact virtual work of a smooth radial
# pressure, but it is not the transpose of the gap Jacobian, so the contact compliance is
# asymmetric and the complementarity problem is not a convex programme.
#
# `:nodal` takes the contact unknown to be the vertical LOAD at each node instead. The forcing
# is then `H'` with no freedom left in it, because `H'` is the derivative of the constraint,
# and the compliance becomes `H A^-1 H' + (1/(m beta^2)) 1 1'` -- symmetric positive definite.
#
# IT IS NOT THE PRODUCTION FORCING, and these tests pin why rather than leaving it to a page.
# It reproduces the validated reference materially worse, and the two are NOT two
# discretisations of one model: each converges to its own limit and the gap does not close.
# What it is good for is the case it was derived for -- a curved or deformable substrate, where
# the contact problem must be assembled additively from two bodies and only the conjugate form
# composes.

using Test
using DropSolver
using LinearAlgebra

@testset "nodal (conjugate) contact forcing" begin

    @testset "the compliance is symmetric positive definite" begin
        # The whole point of the formulation, and the one thing the production forcing cannot
        # offer. Symmetry here is structural: the same matrix appears on both sides of
        # `H A^-1 H'`, and the centre-of-mass term is a positive multiple of an outer square.
        for (M, K) in ((20, 2), (30, 2), (45, 2), (61, 1), (90, 1))
            p = ImpactParams(We = 0.5, Bo = 0.0189, Oh = 0.0373, M = M, K = K,
                             force_mode = :nodal)
            F0 = assemble_newtonian(DropSolver.basis(p), p.Oh)
            Vf = lu(DropSolver.legendre_vandermonde(p))
            s0 = DropSolver.initial_state(p)
            W, _, idx, _ = DropSolver.contact_lcp(p, s0, s0, p.dt0; F0 = F0, Vfac = Vf)
            @test maximum(abs, W - W') / maximum(abs, W) < 1e-12
            @test minimum(eigvals(Symmetric(0.5 .* (W .+ W')))) > 0
            @test length(idx) > 4
        end
    end

    @testset "the net contact force is the weight at rest" begin
        # The replacement for the hydrostatic anchor. With loads as the unknowns the statement
        # is simply that they sum to the weight, which is cleaner than the spectral version
        # (`p_{c,1} -> -Bo`) because no harmonic is involved. `equivalent_pressure` stores
        # `-sum(lam)/m` in the `l = 1` slot, so this reads it back out.
        for Bo in (0.05, 0.2)
            p = ImpactParams(We = 1e-8, Bo = Bo, Oh = 0.5, M = 30, K = 3, t_max = 60.0,
                             force_mode = :nodal)
            r = simulate_lcp(p)
            tail = max(1, length(r.pc) - length(r.pc) ÷ 5)
            sums = [-(4pi/3) * r.pc[i][2] for i in tail:length(r.pc)]
            @test isapprox(sum(sums)/length(sums), (4pi/3)*Bo; rtol = 0.05)
        end
    end

    @testset "it integrates, and complementarity holds exactly" begin
        for (We, Oh, M, K) in ((1.0, 0.3038, 30, 2), (0.5, 0.0373, 30, 2), (2.0, 0.05, 20, 2))
            p = ImpactParams(We = We, Bo = 0.0189, Oh = Oh, M = M, K = K, t_max = 25.0,
                             force_mode = :nodal)
            r = simulate_lcp(p)
            @test r.lcp_resid_max < 1e-8
            @test r.rejects == 0
            @test all(isfinite, r.z)
            m = proximity_metrics(p, r)
            @test 0 < m.cor <= 1                       # no energy created
            @test 0 < m.tc < 20
        end
    end

    @testset "it is NOT the production forcing, and this records why" begin
        # A negative result, pinned so that it cannot be quietly undone by someone flipping the
        # default. Two independent statements.
        #
        # FIRST, against the validated MATLAB reference, CoR = 0.3138 at We = 1, Bo = 0.0189,
        # Oh = 0.303767, M = 90, K = 1. The Legendre forcing reproduces it to 0.29 per cent;
        # the nodal forcing misses by 8.8.
        ref = 0.3138
        pl = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 90, K = 1, t_max = 25.0)
        pn = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 90, K = 1, t_max = 25.0,
                          force_mode = :nodal)
        rl = simulate_lcp(pl); rn = simulate_lcp(pn)
        @test abs(rl.cor - ref)/ref < 0.01             # measured 0.0029
        @test abs(rn.cor - ref)/ref > 0.05             # measured 0.088
        @test abs(rl.cor - ref) < abs(rn.cor - ref)    # and the ordering is the point

        # SECOND, they are not two discretisations of one model. If they were, refining M would
        # close the gap. It does not: each converges to its own limit and the difference sits at
        # about nine per cent in restitution. Measured 7.6, 10.0, 9.0, 9.3, 9.3 per cent at
        # M = 14, 20, 30, 45, 60.
        gaps = Float64[]
        for M in (20, 45)
            a = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = M, K = 3, t_max = 25.0)
            b = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = M, K = 3, t_max = 25.0,
                             force_mode = :nodal)
            ma = proximity_metrics(a, simulate_lcp(a)); mb = proximity_metrics(b, simulate_lcp(b))
            push!(gaps, abs(mb.cor - ma.cor)/ma.cor)
        end
        @test all(g -> g > 0.05, gaps)                 # the gap is real at both truncations
        @test abs(gaps[2] - gaps[1]) < 0.05            # and it does not close with refinement
    end

    @testset "shear thinning runs through the conjugate forcing too" begin
        # The Picard closure on eta is independent of which forcing is used, so this checks that
        # the two are genuinely orthogonal choices rather than entangled.
        eta = gd -> carreau(gd; lambda_c = 10.0, a = 2.0, n = 0.5, eta_inf_ratio = 0.01)
        p = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 30, K = 3, eta = eta,
                         t_max = 25.0, force_mode = :nodal)
        r = simulate_lcp(p)
        @test r.lcp_resid_max < 1e-8
        @test 0 < r.cor <= 1
        ## and the Newtonian limit of the shear-thinning path reproduces the constant-viscosity
        ## one, through the conjugate forcing as through the other
        e1 = gd -> carreau(gd; lambda_c = 1e-8, a = 2.0, n = 0.5, eta_inf_ratio = 0.0)
        pv = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 30, K = 3, eta = e1,
                          t_max = 25.0, force_mode = :nodal)
        pc = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 30, K = 3,
                          t_max = 25.0, force_mode = :nodal)
        @test isapprox(simulate_lcp(pv).cor, simulate_lcp(pc).cor; rtol = 1e-4)
    end
end
