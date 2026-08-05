# Physical invariants, and the cross-validation that neither solver can fake alone.
#
# Four things nothing else in the suite covered. Two are consequences of the physics with an
# ANALYTIC target, so they cannot be satisfied by a plausible-looking wrong answer. One pins
# the agreement between two independent contact algorithms. And one covers the metrics
# themselves, which every experimental comparison now depends on and which had no test at all
# an hour after they were written.

using Test
using DropSolver
using LinearAlgebra

mean_(v) = sum(v) / length(v)

@testset "physical invariants" begin

    @testset "hydrostatic balance: the film carries the weight" begin
        # THE ONLY TEST WITH AN EXACT TARGET FOR THE CONTACT FORCE MAGNITUDE. Release the
        # drop from rest and let it settle: at equilibrium the film must support exactly the
        # weight. The net vertical force is F = -(4 pi/3) p_c,1 and the mass is 4 pi/3, so
        #
        #     p_c,1  ->  -Bo
        #
        # exactly, independent of Oh, M, K and of which contact algorithm is used. Every
        # other pressure test in this suite constrains only a sign or a relative magnitude;
        # this one has a number.
        for (Bo, Oh) in ((0.05, 0.5), (0.20, 0.5), (0.05, 1.0))
            for solver in (simulate, simulate_lcp)
                p = ImpactParams(We = 1e-8, Bo = Bo, Oh = Oh, M = 20, K = 2, t_max = 60.0)
                r = solver(p)
                ## the tail of the run, once the transient has damped
                tail = max(1, length(r.pc1) - length(r.pc1) ÷ 5)
                settled = r.pc1[tail:end]
                @test !isempty(settled)
                @test isapprox(mean_(settled), -Bo; rtol = 0.15)
            end
        end
    end

    @testset "dissipation is ordered by the viscosity, as a matrix inequality" begin
        # A THEOREM, not a trend. If eta(gammadot) <= eta_0 pointwise then the dissipation
        # form is ordered the same way, so C(eta_0) - C(eta) is positive semi-definite --
        # every eigenvalue non-negative, not merely a smaller trace. This is why a
        # shear-thinning drop must rebound harder, and it is the strongest available check on
        # the eta-weighted quadrature: a sign or factor error there breaks the ordering while
        # leaving every scalar KPI plausible.
        b = ModalBasis(2:8, 2)
        st = [0.5/i * (k == 1 ? 1.0 : 0.3) for i in 1:length(b.ls) for k in 1:b.K]
        CN = assemble_coupled(b, 0.1).C                        # eta == 1 everywhere
        for (lam, n) in ((3.0, 0.5), (10.0, 0.5), (10.0, 0.3))
            etaf = (x, mu) -> carreau(shear_rate(b, st, x, mu);
                                      lambda_c = lam, a = 2.0, n = n, eta_inf_ratio = 0.02)
            ## the premise: the viscosity really is below one everywhere it is sampled
            @test all(etaf(x, mu) <= 1.0
                      for x in (0.2,0.5,0.8,1.0), mu in (-0.95,-0.4,0.0,0.6,0.95))
            CS = assemble_coupled(b, 0.1; eta = etaf).C
            D = Symmetric(0.5 .* ((CN .- CS) .+ (CN .- CS)'))
            ev = eigvals(D)
            @test minimum(ev) > -1e-10 * maximum(abs, ev)      # PSD: the ordering holds
            @test maximum(ev) > 0                              # and is not vacuous
        end
    end

    @testset "free flight: the shape does not know the drop is falling" begin
        # Before contact the only coupling between shape and centre of mass is the film, so
        # with no film the two must be exactly independent: the shape trajectory cannot depend
        # on the impact speed, and the centre of mass must follow exact ballistics. This is
        # what makes the contact force the ONLY thing that couples them, which the whole
        # variational structure asserts and nothing else here checks.
        Bo = 0.05
        function preflight(We)
            p = ImpactParams(We = We, Bo = Bo, Oh = 0.3, M = 20, K = 2, t_max = 25.0)
            r = simulate(p)
            i = findfirst(>(0), r.cp)
            i === nothing && return (nothing, nothing, nothing, p)
            k = max(1, i - 2)                                  # strictly before contact
            (r.t[1:k], [surface_amplitudes(p, a) for a in r.a[1:k]], r.z[1:k], p)
        end
        t1, s1, z1, p1 = preflight(1.0)
        @test t1 !== nothing
        ## the drop starts undeformed, so in free flight it must STAY undeformed: any shape
        ## motion before contact would be forcing from nowhere
        @test maximum(maximum(abs, s) for s in s1) < 1e-10
        ## and the centre of mass follows z = z0 + v0 t - Bo t^2/2 exactly
        v0 = -sqrt(1.0)
        worst = maximum(abs(z1[k] - (z1[1] + v0*t1[k] - 0.5*Bo*t1[k]^2)) for k in eachindex(t1))
        @test worst < 1e-8
    end

    @testset "the two contact algorithms agree where both are valid" begin
        # THE STRONGEST CROSS-VALIDATION IN THIS PROJECT, and it was not pinned. A ranked
        # candidate search and a complementarity solve share the assembly but nothing about
        # how they choose the contact set, so agreement is evidence about the assembly and
        # about both closures at once. Measured at four decimal places; if either drifts, this
        # is what catches it.
        for (We, Oh) in ((1.0, 0.303767), (5.0, 0.30), (10.0, 0.30))
            p = ImpactParams(We = We, Bo = 0.0189, Oh = Oh, M = 45, K = 1, t_max = 25.0)
            rs = simulate(p); rl = simulate_lcp(p)
            @test isfinite(rs.cor) && isfinite(rl.cor)
            @test isapprox(rs.cor, rl.cor; rtol = 1e-3)
            @test isapprox(rs.tc,  rl.tc;  rtol = 1e-3)
            ## and under the proximity metric too, which is what the comparisons use
            ms = proximity_metrics(p, rs); ml = proximity_metrics(p, rl)
            @test isapprox(ms.cor, ml.cor; rtol = 1e-3)
            @test isapprox(ms.tc,  ml.tc;  rtol = 1e-3)
        end
    end

    @testset "the proximity metrics mean what they say" begin
        # These define every experimental comparison and had no test. Four properties, each
        # falsifiable: contact is detected when the surface comes within the threshold, a
        # larger threshold can only lengthen the contact, restitution is the speed ratio at
        # the detected endpoints, and the drop is descending at the first detection and
        # rising at the last.
        p = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 45, K = 2, t_max = 25.0)
        r = simulate(p)
        m = proximity_metrics(p, r; h_thresh = 0.02)
        @test m.n_detected > 0
        @test m.i_first < m.i_last
        ## the definition, recomputed from the trajectory
        @test m.tc ≈ r.t[m.i_last] - r.t[m.i_first]
        @test m.cor ≈ abs(r.v[m.i_last] / r.v[m.i_first])
        ## approaching at the start, leaving at the end
        @test m.v_in < 0
        @test m.v_out > 0
        ## a looser threshold detects contact no later and loses it no sooner
        loose = proximity_metrics(p, r; h_thresh = 0.05)
        tight = proximity_metrics(p, r; h_thresh = 0.005)
        @test loose.tc >= m.tc >= tight.tc
        @test loose.i_first <= m.i_first <= tight.i_first
        ## and the minimum gap actually reached is below every threshold that detected contact
        @test m.min_gap < 0.005
    end
end

mean_(v) = sum(v) / length(v)
