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
        p = ImpactParams(; REF..., M = 90, K = 1, t_max = 25.0)
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
        r1 = simulate(ImpactParams(; REF..., M = 45, K = 1, t_max = 25.0))
        r2 = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 25.0))
        @test r2.cor > r1.cor                             # more rebound
        @test minimum(r2.z) < minimum(r1.z)               # deeper penetration
        @test r2.cor / r1.cor > 1.05                      # and not a small effect
    end

    @testset "two radial functions suffice" begin
        # The convergence that makes the extension cheap: the exact interior profile is
        # a spherical Bessel function, these are its Taylor terms, and it truncates fast.
        cors = [simulate(ImpactParams(; REF..., M = 45, K = K, t_max = 25.0)).cor
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
        cors = [simulate(ImpactParams(; REF..., M = M, K = 2, t_max = 25.0)).cor
                for M in Ms]
        @test maximum(cors) - minimum(cors) < 0.01 * minimum(cors)
        tcs = [simulate(ImpactParams(; REF..., M = M, K = 2, t_max = 25.0)).tc
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
        rs = [simulate(ImpactParams(; low..., M = M, K = 2, t_max = 25.0))
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
        r = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 25.0))
        inc = findfirst(>(0), r.cp); lastc = findlast(>(0), r.cp)
        incontact = [i for i in inc:lastc if r.cp[i] > 2]
        @test !isempty(incontact)
        suction = [i for i in incontact if r.pc1[i] >= 0]
        @test length(suction) < 0.05 * length(incontact)          # measured: 1.7%
        peak_push = maximum(abs(r.pc1[i]) for i in incontact)
        @test all(i -> r.pc1[i] < 0.02 * peak_push, incontact)     # measured: 0.3%
    end

    @testset "the contact evolves as one connected episode" begin
        r = simulate(ImpactParams(; REF..., M = 45, K = 2, t_max = 25.0))
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

# Shear thinning, in the time stepper rather than in the assembly.
#
# The assembly is already covered by test_variational.jl, where the Gaunt selection
# rule is verified exactly and Carreau-Yasuda is checked on real states. What is new
# here is that a state-dependent viscosity survives being integrated through a
# contact, and that its effect on the KPIs has the sign and the ordering the rheology
# demands. A shear-thinning drop dissipates less than a Newtonian one with the same
# zero-shear viscosity, because the viscosity falls exactly where the shear is
# strongest, so it must rebound HARDER -- and more so the stronger the thinning.
@testset "shear thinning through a contact" begin
    # A modest truncation, because a non-constant viscosity forces the full coupled
    # reassembly -- O(ndof^2) quadratures -- instead of the cached block-diagonal
    # operator a constant viscosity allows. That is a cost, not an approximation.
    base = (We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 30, K = 3, t_max = 25.0)

    rn = simulate(ImpactParams(; base...))
    @test rn.rejects == 0
    @test 0 < rn.cor < 1

    carr(lam, n) = ImpactParams(; base...,
        eta = gd -> carreau(gd; lambda_c = lam, a = 2.0, n = n, eta_inf_ratio = 0.01))

    @testset "a shear-thinning drop rebounds harder than a Newtonian one" begin
        p = carr(10.0, 0.5)
        @test !p.eta_const                       # the probe really sees a variable eta
        r = simulate(p)
        ## A rejected step is no longer a symptom of limping: the viscosity iteration is
        ## monitored, and a step whose iteration will not converge is rejected so that dt
        ## halves, which is the only lever that restores the contraction. So what must
        ## hold is that the iteration DID converge on every accepted step, and that it
        ## did not need the whole sweep budget to get there on most of them.
        @test r.eta_resid_max <= p.eta_tol       # every accepted step converged
        @test r.rejects < 0.05 * length(r.t)     # and rejection stayed occasional
        @test r.cor > rn.cor
        @test 0 < r.cor < 1
    end

    @testset "and harder the stronger the thinning" begin
        # Two knobs, each with a direction fixed by the rheology and not by fitting:
        # raising lambda_c thins at a lower shear rate, lowering n thins more steeply.
        # Both must increase the rebound, monotonically.
        c_weak   = simulate(carr(1.0,  0.5)).cor
        c_strong = simulate(carr(10.0, 0.5)).cor
        c_steep  = simulate(carr(10.0, 0.3)).cor
        @test rn.cor < c_weak < c_strong < c_steep
        @test c_steep / rn.cor > 1.5             # measured: +91 per cent
    end

    @testset "the Newtonian limit is recovered through the stepper" begin
        # lambda_c -> 0 sends eta -> 1 pointwise, so the whole shear-thinning path --
        # a coupled reassembly every sweep -- must reproduce the cached constant
        # viscosity operator it bypasses. This is the check that the two code paths
        # agree, not merely that each is self-consistent.
        r = simulate(carr(1e-8, 0.3))
        @test isapprox(r.cor, rn.cor; rtol = 1e-4)
        @test isapprox(r.tc,  rn.tc;  rtol = 1e-4)
    end
end

# The 3000 ppm shear-thinning fluid, and the one thing it settles.
#
# This fluid's ZERO-SHEAR Ohnesorge is 57. A Newtonian drop that viscous does not
# rebound at all -- it arrives, spreads, and stays down. The experiments nonetheless
# measure restitution coefficients between 0.57 and 0.85. So shear thinning is not a
# correction to the rebound of this drop; it is the reason there is a rebound to
# describe. That is a qualitative prediction with nothing fitted to it, and it is the
# sharpest test the rheology gets: a model that merely adjusted the damping a little
# would rebound at both viscosities and would be wrong about the physics while looking
# right on a plot.
@testset "the 3000 ppm fluid rebounds only because it thins" begin
    ## Cross fit of the fluid, relabelled to Carreau-Yasuda by lambda_c = K, a = m,
    ## n = 1 - m -- the correspondence derived on "Cross-Model Fluids".
    eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
    K_cross, m_cross = 18.48081673111359, 0.7430524574330837
    R, sigma, rho = 0.0003, 0.0728, 1000.0
    t_cap = sqrt(rho * R^3 / sigma)
    Oh_0 = eta_0 / sqrt(rho * sigma * R)
    @test Oh_0 > 50                                   # genuinely, extremely viscous

    common = (We = 0.19, Bo = 0.012, Oh = Oh_0, M = 30, K = 3, t_max = 25.0)

    ## At the zero-shear viscosity the drop must NOT come back.
    rn = simulate(ImpactParams(; common...))
    @test isnan(rn.cor) || rn.cor < 0.05              # no rebound to speak of

    ## With the fluid's own thinning it must, and near the measured value: the
    ## experiments at this Weber number average CoR ~ 0.80.
    rst = simulate(ImpactParams(; common...,
        eta = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                            n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0)))
    ## NOT `eta_resid_max <= eta_tol`. That asks whether the residual is below the
    ## number we told it to reach, which is true by construction whenever the loop
    ## exits normally, and it fails the moment the tolerance is changed even if
    ## nothing physical moved. It is a statement about the machinery, not the science.
    ##
    ## The question worth asking is whether the default settings change the answer.
    ## Tighten the tolerance by three decades and the restitution and contact time
    ## must not move: if they do, the default is buying accuracy the model needs, and
    ## if they do not, the extra sweeps are being spent for nothing.
    ##
    ## Physical meaning of a failure: the reported bounce depends on a convergence
    ## threshold, so it is a property of the iteration rather than of the fluid.
    rst_tight = simulate(ImpactParams(; common..., eta_tol = 1e-9,
        eta = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                            n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0)))
    @test isapprox(rst.cor, rst_tight.cor; rtol = 2e-3)
    @test isapprox(rst.tc,  rst_tight.tc;  rtol = 2e-3)
    @test 0.6 < rst.cor < 0.95
    @test isapprox(rst.cor, 0.80; rtol = 0.20)
    ## and the contact time within the same band as the measurement, ~2.7
    @test isapprox(rst.tc, 2.7; rtol = 0.30)
end

@testset "a viscosity failure does not collapse the step size" begin
    # The march meets two different viscosity-iteration failures and they need opposite
    # responses. Early on the iteration DIVERGES, residual of order 1e20, and halving dt
    # rescues it. Later it reaches a FLOOR just above the tolerance, and halving raises
    # that floor, because `adot_star = beta*a + hv_a` with `beta = c0/dt` is a difference
    # of large quantities, so the increment between sweeps is cancellation noise scaled by
    # beta. Measured at M = 90: 3.0e-6 at the natural step, 2.0e-2 at dt = 5.6e-8.
    #
    # Three versions of the handler were wrong before this one. Stopping on any viscosity
    # failure killed every run at t = 0, since the first step always trips. Stopping when
    # the residual worsens fired on the divergent cluster, which is not monotone. Both
    # passed the suite, because every impact test then ran at M = 14, K = 2, where the
    # first step does not trip at all. They run at M = 30 now for that reason.
    #
    # So this runs at a resolution that DOES exercise the path. A failure means either the
    # early divergence is no longer being rescued by step reduction, or the floor is no
    # longer being detected and the march is spiralling to dt_min again.
    eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
    K_cross, m_cross = 18.48081673111359, 0.7430524574330837
    R, sigma, rho = 0.0003, 0.0728, 1000.0
    t_cap = sqrt(rho * R^3 / sigma)
    Oh_0 = eta_0 / sqrt(rho * sigma * R)
    eta_fn = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                           n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0)
    p = ImpactParams(We = 0.3643, Bo = 0.012, Oh = Oh_0, M = 30, K = 3,
                     eta = eta_fn, t_max = 25.0)
    r = simulate_lcp(p)
    m = proximity_metrics(p, r)
    @test !isempty(r.t)
    @test r.t[end] > 2.0                       # the bounce completed, not a truncated march
    @test 0.70 < m.cor < 0.80                  # and gives the converged value, ~0.7535
end

@testset "the defaults themselves produce a bounce" begin
    # Every other test names its own M and K, so none of them exercises what a user
    # actually gets. These pass neither: whatever DEFAULT_M and DEFAULT_K are, this is
    # the configuration the package ships, and it has to work.
    #
    # It is not hypothetical. DEFAULT_M and DEFAULT_K were set independently, and the
    # pair M = 90, K = 3 could not complete a shear-thinning march at all until the step
    # controller was fixed: the viscosity iteration hit a floor the tolerance was below,
    # dt halved twenty-three times and the run died having computed a third of a bounce
    # while reporting a restitution of 0.15. Nothing caught it, because every impact test
    # then ran at M = 14.
    #
    # Physical meaning of a failure: the package's own defaults do not produce a drop
    # that bounces, whatever the tests at chosen resolutions say.
    @test DropSolver.DEFAULT_M == 90
    @test DropSolver.DEFAULT_K == 3

    @testset "Newtonian" begin
        p = ImpactParams(We = 0.5, Bo = 0.019, Oh = 0.0373)
        r = simulate_lcp(p)
        m = proximity_metrics(p, r)
        @test r.t[end] > 2.0                    # completed, not truncated
        @test 0.70 < m.cor < 0.80               # converged value is 0.7556
        @test 2.0 < m.tc < 3.0
    end

    @testset "shear thinning on the 3000 ppm fluid" begin
        eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
        K_cross, m_cross = 18.48081673111359, 0.7430524574330837
        R, sigma, rho = 0.0003, 0.0728, 1000.0
        t_cap = sqrt(rho * R^3 / sigma)
        Oh_0 = eta_0 / sqrt(rho * sigma * R)
        p = ImpactParams(We = 0.3643, Bo = 0.012, Oh = Oh_0,
                         eta = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                                             n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0))
        r = simulate_lcp(p)
        m = proximity_metrics(p, r)
        @test r.t[end] > 2.0                    # this is the case that used to die at 0.890
        @test 0.70 < m.cor < 0.80               # measured 0.753965 at the defaults
    end
end
