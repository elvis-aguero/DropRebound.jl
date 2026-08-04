# The variational contact solver, against the validated Newtonian implementation.
#
# The reference numbers come from running Gabbard et al. (2025) directly -- the MATLAB
# code of the paper, at its own production settings (M = 90 harmonics, Legendre-root
# collocation, its own fixed timestep) and at the parameters of its own sweep:
#
#   sigma = 20.5, rho = 0.96, R = 0.0203 cm, nu = 0.2 St, Bo = 0.0189, We = 1
#   -> Oh = 0.303767
#   -> CoR = 0.3138, contact time = 2.183, max contact 25 of 91 nodes, z_min = 0.63178
#
# Reproducing those at K = 1 is the check that matters, because K = 1 IS their model:
# one radial trial function is potential flow, and its modal damping is Lamb's. With
# K >= 2 the interior is resolved and the damping becomes Reid's, so the numbers
# SHOULD move -- and the direction is predictable in advance, because Lamb over-damps.

using Test
using DropSolver

const REF = (We = 1.0, Bo = 0.0189, Oh = 0.303767)
const REF_COR  = 0.3138
const REF_TC   = 2.183
const REF_CP   = 25
const REF_ZMIN = 0.63178

@testset "variational contact solver vs the validated Newtonian code" begin

    @testset "the collocation nodes are Legendre roots, clustered at the pole" begin
        # Not a cosmetic choice. Uniform nodes make the first contact a wide wedge,
        # which takes an impulsive pressure to hold and returns more energy than the
        # drop arrived with. The clustering is what makes a two-node contact small.
        p = ImpactParams(; REF..., M = 90, K = 1)
        @test length(p.nodes) == 91              # square: one node per pressure harmonic
        @test p.nodes[1] ≈ pi                    # a node exactly at the pole
        @test issorted(p.nodes; rev = true)
        @test rad2deg(pi - p.nodes[2]) < 2.0     # second node within two degrees
        # and they really are the zeros of P_M, not merely clustered
        @test abs(DropSolver.legendre_angular(90, cos(p.nodes[2])).P) < 1e-9
    end

    @testset "K = 1 reproduces the reference run" begin
        p = ImpactParams(; REF..., M = 90, K = 1, t_max = 3.0)
        r = simulate(p)
        @test maximum(r.cp) == REF_CP                     # identical contact extent
        @test isapprox(r.cor, REF_COR; atol = 0.01)
        @test isapprox(r.tc, REF_TC; rtol = 0.08)
        @test isapprox(minimum(r.z), REF_ZMIN; rtol = 0.01)
        @test r.rejects < 10                              # not limping through on dt
    end

    @testset "the interior changes the answer, in the predictable direction" begin
        # Lamb over-damps every mode at this Ohnesorge -- 37 per cent at l = 2, 143
        # per cent at l = 8 -- so resolving the interior must REDUCE the damping, and
        # less damping must mean more rebound and deeper penetration. A failure here
        # would mean the interior was feeding energy the wrong way.
        r1 = simulate(ImpactParams(; REF..., M = 45, K = 1, t_max = 3.0))
        r2 = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 3.0))
        @test r2.cor > r1.cor                             # more rebound
        @test minimum(r2.z) < minimum(r1.z)               # deeper penetration
        @test r2.cor / r1.cor > 1.05                      # and not a small effect
    end

    @testset "two radial functions suffice" begin
        # The convergence that makes the extension cheap: the exact interior profile is
        # a spherical Bessel function, these are its Taylor terms, and it truncates fast.
        cors = [simulate(ImpactParams(; REF..., M = 45, K = K, t_max = 3.0)).cor
                for K in 2:4]
        @test maximum(cors) - minimum(cors) < 2e-3
    end

    @testset "the answer is insensitive to the harmonic truncation" begin
        # Measured: 0.36206, 0.36098, 0.36150, 0.36183 at M = 20, 30, 45, 90 -- a
        # spread of 0.3 per cent with no trend, so CoR is converged from M = 20 and the
        # truncation is not a tuning knob. It is not asserted that larger M is closer
        # to M = 90, because they are all equally close; an earlier version of this
        # test demanded monotone improvement and failed for that reason.
        #
        # Contact time is NOT converged the same way -- 2.199 to 2.258 over the same
        # range, drifting upward -- because it depends on when the outermost node
        # releases, and how finely the outer contact is resolved is exactly what M
        # controls. That is a real M-dependence, so it is bounded, not denied.
        Ms = (20, 30, 45, 90)
        cors = [simulate(ImpactParams(; REF..., M = M, K = 2, t_max = 3.0)).cor
                for M in Ms]
        @test maximum(cors) - minimum(cors) < 0.01 * minimum(cors)
        tcs = [simulate(ImpactParams(; REF..., M = M, K = 2, t_max = 3.0)).tc
               for M in Ms]
        @test maximum(tcs) - minimum(tcs) < 0.05 * minimum(tcs)
    end

    @testset "the active-set iteration is stable where residual ranking was not" begin
        # The regime that broke the ranked search: Oh = 0.023, almost no damping, so
        # the high-l modes ring and two candidate contact counts score almost equally.
        # Ranking chattered -- the film pressure oscillated in sign with order-one
        # amplitude and dt collapsed -- and it did so SPORADICALLY in M, failing at
        # M = 35 and 45 while coming through at 30 and 40. Under the active-set
        # iteration every one of these runs completes with no rejected step and CoR
        # converges monotonically, which is the evidence that the chatter was the
        # tie-breaking and not the physics.
        low = (We = 0.9899, Bo = 0.05263, Oh = 0.0233)
        rs = [simulate(ImpactParams(; low..., M = M, K = 2, t_max = 6.0))
              for M in (16, 30, 35, 45)]
        @test all(r -> r.rejects == 0, rs)
        @test all(r -> isfinite(r.cor) && 0 < r.cor < 1, rs)
        @test all(r -> r.tc > 2.0, rs)
        cs = [r.cor for r in rs]
        @test maximum(cs) - minimum(cs) < 0.02 * minimum(cs)
    end

    @testset "the film pressure pushes, almost everywhere, without being told to" begin
        # Signorini's inequality is NOT imposed anywhere in the solver. The validated
        # code does not impose it either, and an earlier version here that did rejected
        # admissible steps at onset and stalled. So its holding is a RESULT rather than
        # a constraint: through contact the net vertical force should be upward, which
        # for a pressure concentrated near mu = -1 means p_c,1 < 0.
        #
        # It holds for the bulk of contact and FAILS BRIEFLY at both transitions. This
        # test pins the measured extent of that failure rather than asserting it away,
        # because the failure is the model's, not the test's: nothing forbids the film
        # from pulling, so a few steps either side of onset and release come out
        # adhesive. What makes it tolerable is the magnitude -- suction of order 0.01
        # against a peak push of 3.4, so about a third of a per cent -- and what would
        # make it intolerable is that number growing. If this test starts failing
        # because the count or the magnitude rose, the film needs a real release
        # criterion, not a looser tolerance.
        r = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 3.0))
        inc = findfirst(>(0), r.cp); lastc = findlast(>(0), r.cp)
        incontact = [i for i in inc:lastc if r.cp[i] > 2]
        @test !isempty(incontact)
        suction = [i for i in incontact if r.pc1[i] >= 0]
        @test length(suction) < 0.05 * length(incontact)          # measured: 1.7%
        peak_push = maximum(abs(r.pc1[i]) for i in incontact)
        @test all(i -> r.pc1[i] < 0.02 * peak_push, incontact)     # measured: 0.3%
    end

    @testset "the contact evolves as one connected episode" begin
        r = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 3.0))
        # ONE contact, not a sequence of micro-bounces -- but the set does flicker
        # empty briefly right after first touch. Measured: two episodes, of two steps
        # and one step, at t = 0.013 and t = 0.023 out of a contact lasting to t = 2.27,
        # so 0.34 per cent of the contact, and both while the drop is still approaching
        # at v = -0.99 rather than leaving. That is onset transient, not separation, and
        # the honest assertion bounds it instead of denying it. What would matter is a
        # gap late in contact or one long enough to let the drop move: neither occurs.
        inc = findfirst(>(0), r.cp); lastc = findlast(>(0), r.cp)
        span = r.t[lastc] - r.t[inc]
        empt = [i for i in inc:lastc if r.cp[i] == 0]
        @test length(empt) < 0.02 * (lastc - inc)
        @test all(i -> (r.t[i] - r.t[inc]) < 0.05 * span, empt)   # all near onset
        @test all(i -> r.v[i] < 0, empt)                          # still approaching
        # It may move SEVERAL nodes in a step, because the active-set iteration runs to
        # convergence within the step rather than being capped at one move -- measured
        # maximum 4. The old one-node-per-step limit was a property of the ranked
        # search, and asserting it here would be asserting the design that chattered.
        # What must stay bounded is the jump relative to the contact itself.
        @test maximum(abs(r.cp[i+1] - r.cp[i]) for i in 1:(length(r.cp)-1)) <=
              max(4, maximum(r.cp) ÷ 4)
        # It grows to a single peak and retreats, with occasional reversals -- 7 in ~880
        # steps. Bound them rather than forbid them.
        pk = argmax(r.cp)
        rev = count(i -> r.cp[i+1] < r.cp[i], 1:(pk-1)) +
              count(i -> r.cp[i+1] > r.cp[i], pk:(length(r.cp)-1))
        @test rev < 0.05 * length(r.cp)
    end
end
