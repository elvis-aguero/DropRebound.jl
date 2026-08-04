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

    @testset "refining the harmonics converges" begin
        cors = [simulate(ImpactParams(; REF..., M = M, K = 2, t_max = 3.0)).cor
                for M in (20, 45, 90)]
        @test abs(cors[2] - cors[3]) < 0.01 * cors[3]     # M = 45 within one per cent
        @test abs(cors[2] - cors[3]) < abs(cors[1] - cors[3])
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

    @testset "the contact set moves one node at a time" begin
        r = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 3.0))
        @test all(abs(r.cp[i+1] - r.cp[i]) <= 1 for i in 1:(length(r.cp)-1))
        # One connected contact episode, not a chatter of many: the set is nonempty
        # from first touch to release with no gap in between.
        inc = findfirst(>(0), r.cp); lastc = findlast(>(0), r.cp)
        @test all(>(0), r.cp[inc:lastc])
        # It grows to a single peak and retreats, but NOT monotonically -- an active
        # set that may move one node per step, ranked by an edge residual, reverses
        # occasionally. Measured: 7 reversals while growing, 1 while retreating, out
        # of ~900 steps. Bound the reversals rather than forbid them.
        pk = argmax(r.cp)
        rev = count(i -> r.cp[i+1] < r.cp[i], 1:(pk-1)) +
              count(i -> r.cp[i+1] > r.cp[i], pk:(length(r.cp)-1))
        @test rev < 0.05 * length(r.cp)
    end
end
