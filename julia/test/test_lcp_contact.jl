# The contact set solved for rather than searched over.
#
# `simulate_lcp` determines which nodes are in contact by solving the complementarity
# problem `h = A_c p + b`, `h >= 0`, `p >= 0`, `p'h = 0` at every step. The tests below
# are about the two things that makes possible and the searching solver cannot do: it
# satisfies the pressure inequality, and it is not restricted to a disc.

using Test
using DropSolver
using LinearAlgebra

const LREF = (We = 1.0, Bo = 0.0189, Oh = 0.303767)

@testset "LCP contact" begin

    @testset "the complementarity conditions are satisfied, not approximated" begin
        # The searching solver enforces h >= 0 and h*p = 0 by construction but NOT
        # p >= 0 -- measured at 0.015 per cent of node-steps in a case that works and 41
        # per cent in one that does not. This is the condition the LCP adds, so it is the
        # thing to check, and the complementarity residual is a certificate rather than a
        # tolerance: it vanishes exactly when all three hold.
        p = ImpactParams(; LREF..., M = 45, K = 2, t_max = 25.0)
        r = simulate_lcp(p)
        @test r.lcp_resid_max < 1e-8
        @test r.rejects == 0                       # the solve never had to be retried
        @test r.lcp_sweeps_max < 100               # and converged cheaply
        ## the pressure is non-negative at every node of every step, by construction
        V = DropSolver.legendre_vandermonde(p)
        worst = 0.0
        for i in eachindex(r.t)
            worst = min(worst, minimum(V * r.pc[i]))
        end
        @test worst > -1e-6 * maximum(maximum(abs, V * r.pc[i]) for i in eachindex(r.t))
    end

    @testset "restitution agrees with two independent methods" begin
        # CoR at K = 1, M = 90 is 0.3129 from the LCP, 0.3129 from the tangency rule, and
        # 0.3138 from the reference MATLAB implementation. Three routes to the same
        # number, which is what makes it evidence rather than coincidence.
        r = simulate_lcp(ImpactParams(; LREF..., M = 90, K = 1, t_max = 25.0))
        @test isapprox(r.cor, 0.3138; atol = 0.01)
    end

    @testset "the cases that defeated the searching solver now run" begin
        # Oh = 0.03 at We = 5 and 10 manufactured forty times the drop's energy under the
        # active-set search and starved under tangency selection. The reference
        # implementation completes them with |zeta| = 0.99 and 1.22.
        for (We, Oh, zref) in ((5.0, 0.03, 0.99), (10.0, 0.03, 1.22), (10.0, 0.30, 0.84))
            p = ImpactParams(We = We, Bo = 0.0189, Oh = Oh, M = 90, K = 1, t_max = 30.0)
            r = simulate_lcp(p)
            @test r.rejects == 0
            @test r.lcp_resid_max < 1e-6
            amp = maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a)
            @test isapprox(amp, zref; rtol = 0.15)      # and the amplitude is physical
            ## the energy ceiling from the surface stiffness, which the old runs blew
            ## through by a factor of five
            F = assemble_newtonian(DropSolver.basis(p), Oh)
            N = length(r.a[1])
            cheap = minimum(begin e = zeros(N); e[i] = 1.0; 0.5*dot(e, F.G, e) end
                            for i in 1:N if begin e = zeros(N); e[i] = 1.0
                                                 dot(e, F.G, e) end > 0)
            @test amp <= sqrt(0.5*(4pi/3)*We / cheap)
        end
    end

    @testset "the contact is ANNULAR, and it converges in M" begin
        # The axiom on the model page -- "contact occupies a single connected patch about
        # the pole" -- is false. Solved without that restriction, the pressure is exactly
        # zero at the pole and order ten on a ring a few nodes out: a dimple with air
        # trapped under it, which is a known feature of drop impact and which the interval
        # parametrisation cannot represent at all.
        #
        # That it is not a discretisation artefact is the point of this test. An artefact
        # would span a fixed NUMBER of nodes and its angular extent would shrink with the
        # node spacing. Measured instead: the free arc at the pole converges to about 11
        # degrees (10.90, 11.17, 11.44 at M = 45, 60, 90) while the nodes spanning it grow
        # (3, 4, 6). The angle is the invariant, so the feature is resolved rather than
        # created by the grid.
        #
        # What this does NOT establish is that a dimple survives NONLINEAR kinematics.
        # These amplitudes reach |zeta| ~ 0.4, where a linear shape expansion is stretched,
        # and settling that needs a nonlinear reference rather than a finer grid.
        function free_arc(M)
            p = ImpactParams(; LREF..., M = M, K = 2, t_max = 25.0)
            F0 = assemble_newtonian(DropSolver.basis(p), p.Oh)
            Vf = lu(DropSolver.legendre_vandermonde(p))
            prev = curr = DropSolver.initial_state(p); dt = p.dt0
            for _ in 1:200000
                curr.z <= 0.95 && break
                st, nxt, _ = DropSolver.try_step_lcp(p, prev, curr, dt; F0=F0, Vfac=Vf)
                st === :ok || (dt /= 2; dt < p.dt_min && break; continue)
                prev, curr = curr, nxt; dt = min(2dt, p.dt0)
            end
            Ac, bv, idx, _ = DropSolver.contact_lcp(p, prev, curr, dt; F0=F0, Vfac=Vf)
            pv, _, _ = DropSolver.lcp_pgs(Ac, bv)
            act = [i for i in eachindex(pv) if pv[i] > 1e-10*max(maximum(pv), 1)]
            isempty(act) && return (NaN, 0)
            (180 - rad2deg(p.nodes[idx[first(act)]]), first(act) - 1)
        end
        a45, n45 = free_arc(45)
        a60, n60 = free_arc(60)
        a90, n90 = free_arc(90)
        @test all(isfinite, (a45, a60, a90))
        @test n45 > 0 && n90 > n45                    # more nodes resolve the same arc
        @test isapprox(a90, a60; rtol = 0.10)          # and the ARC is what converges
        @test 5.0 < a90 < 25.0                         # a real angular extent, not a node
    end
end

# Shear thinning through the complementarity solve.
#
# A variable viscosity makes the problem NONLINEAR: eta depends on the strain rate, the
# dissipation operator on eta, and the compliance on that operator -- so A_c is a function
# of the velocity the solve produces. It is closed by the same Picard iteration the
# searching solver uses, freezing eta at an extrapolated strain rate and re-evaluating at
# the answer. Each iterate is a genuine LCP with a symmetric PSD compliance, so
# complementarity holds exactly at every sweep and only eta lags.
@testset "shear thinning through the LCP" begin
    ## the 3000 ppm fluid, from its own Cross fit -- nothing fitted to the impact data
    eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
    K_cross, m_cross = 18.48081673111359, 0.7430524574330837
    R, sigma, rho = 0.0003, 0.0728, 1000.0
    t_cap = sqrt(rho * R^3 / sigma)
    Oh_0 = eta_0 / sqrt(rho * sigma * R)
    etaf = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                         n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0)

    @testset "it integrates, and complementarity still holds exactly" begin
        p = ImpactParams(We = 0.1912, Bo = 0.012, Oh = Oh_0, M = 14, K = 2,
                         eta = etaf, t_max = 25.0)
        @test !p.eta_const                        # the probe really sees a variable eta
        r = simulate_lcp(p)
        ## the LCP residual is unaffected by the nonlinearity: every Picard iterate is a
        ## genuine complementarity solve, so this stays at machine level
        @test r.lcp_resid_max < 1e-10
        ## and the viscosity iteration converged on every accepted step
        @test r.eta_resid_max <= p.eta_tol
        @test r.rejects < 0.05 * length(r.t)
        @test 0 < r.cor < 1
    end

    @testset "and it agrees with experiment at least as well as the search" begin
        # Measured at We = 0.1912, where the experiments average CoR 0.804 and tc 2.693:
        #   LCP    CoR 0.7872 (2.1%)   tc 2.4597 (8.7%)
        #   search CoR 0.7673 (4.6%)   tc 2.7180 (0.9%)
        # So the LCP is better on restitution and worse on contact time -- the same split
        # as the Newtonian comparison, which is itself evidence that the difference is in
        # the contact model rather than in the rheology.
        p = ImpactParams(We = 0.1912, Bo = 0.012, Oh = Oh_0, M = 14, K = 2,
                         eta = etaf, t_max = 25.0)
        rl = simulate_lcp(p); rs = simulate(p)
        @test abs(rl.cor - 0.804)/0.804 < 0.10        # within ten per cent on restitution
        @test abs(rl.cor - 0.804) < abs(rs.cor - 0.804)   # and closer than the search
        @test abs(rl.tc - 2.693)/2.693 < 0.15         # contact time within fifteen
    end

    @testset "the Newtonian limit is recovered through the LCP path" begin
        # lambda_c -> 0 sends eta -> 1 pointwise, so the expensive path -- a coupled
        # reassembly and a fresh compliance every sweep -- must reproduce the cached
        # constant-viscosity one. This checks the two code paths against each other rather
        # than each against itself.
        base = (We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 14, K = 2, t_max = 25.0)
        rn = simulate_lcp(ImpactParams(; base...))
        rt = simulate_lcp(ImpactParams(; base...,
                eta = gd -> carreau(gd; lambda_c = 1e-10, a = 2.0, n = 0.4,
                                    eta_inf_ratio = 0.01)))
        @test isapprox(rt.cor, rn.cor; rtol = 1e-3)
        @test isapprox(rt.tc,  rn.tc;  rtol = 1e-3)
    end
end
