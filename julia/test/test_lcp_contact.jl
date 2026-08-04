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
