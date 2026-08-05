# The nonvariational model with a complementarity contact closure.
#
# This is the fourth cell of a two-by-two: {variational, nonvariational} formulation crossed
# with {ranked search, complementarity} contact closure. The point of having all four live is
# that a disagreement then localises itself. If the two formulations agree but the two closures
# do not, the contact treatment is responsible; if the two closures agree but the formulations
# do not, the physics of the interior is. With only one of each, every disagreement is
# confounded and every agreement is uninformative.
#
# WHY THIS CELL IS EXACT RATHER THAN APPROXIMATE. In the Newtonian case every block of the
# nonvariational residual -- kinematics, equation of motion, the pressure conditions, and the
# centre-of-mass pair -- is affine in the unknowns. So the map from contact pressure to gap is
# affine too, the complementarity problem is the exact step condition rather than a
# linearisation of it, and its solution is the step. The first test below is that premise:
# if the Jacobian ever stops being state-independent, everything after it is void.

using Test
using DropSolver
using LinearAlgebra

const NV_M  = 8
const NV_Oh = 0.3038
const NV_Bo = 1/53.9

function nv_cfg(M = NV_M; Oh = NV_Oh, Bo = NV_Bo)
    SimConstants(M, M+1, Oh, Bo, make_theta_vec(M),
                 precompute_integrals(NaN, M)[1], make_dt_max(M))
end
nv_init(M = NV_M; z = 1.02, v = -0.4) = begin
    s = DropState(M); s.z = z; s.v = v; s.dt = make_dt_max(M); s.cp = 0; s
end

@testset "nonvariational LCP" begin

    @testset "the premise: the Newtonian system really is affine" begin
        # The whole construction rests on this. Build the Jacobian at two unrelated states and
        # at the same step size: if the model is affine they are the same matrix to rounding.
        # A failure here would mean the compliance map is only a local linearisation, the
        # complementarity solution is no longer the exact step, and the LCP would need an
        # outer Newton iteration wrapped around it.
        cfg = nv_cfg(); ob = OBParams(0.0, 1.0); dt = cfg.dt_max
        s1 = nv_init(); s2 = nv_init()
        s2.A[2] = 0.13; s2.A[3] = -0.07; s2.Adot[2] = 0.4; s2.z = 0.91; s2.v = 0.2
        s2.B .= collect(range(-0.3, 0.5; length = length(s2.B)))
        for cp in (0, 3, NV_M+1)
            J1 = build_jacobian(s1, [deepcopy(s1)], dt, cp, cfg, ob)
            J2 = build_jacobian(s2, [deepcopy(s1)], dt, cp, cfg, ob)
            @test maximum(abs, J1 .- J2) < 1e-12 * max(1.0, maximum(abs, J1))
        end
    end

    @testset "it runs, and the drop bounces" begin
        cfg = nv_cfg(; Oh = 0.05, Bo = 0.01)
        _, states = solve_drop_lcp!(cfg, OBParams(), nv_init(NV_M; z = 1.05, v = -0.4);
                                    t_end = 15.0, save_every = 0.05)
        @test all(isfinite(s.z) for s in states)
        @test all(all(isfinite, s.A) for s in states)
        @test any(s.cp > 0 for s in states)                 # it touched
        i = findlast(s -> s.cp > 0, states)
        @test i < length(states)                            # and it let go
        @test states[end].v > 0                             # moving up at the end
    end

    @testset "the compliance is positive definite but NOT symmetric" begin
        # This is the structural difference between the two formulations, and it decides which
        # contact solver is legitimate. The variational compliance is a Hessian of an energy,
        # so it is symmetric and the complementarity problem is a convex quadratic program that
        # a projected Gauss-Seidel sweep solves. Here there is no such energy -- the model is
        # built from per-mode damping coefficients -- and the compliance comes out asymmetric
        # by more than half its own size, so that sweep is not applicable.
        #
        # What survives is positive definiteness, which is all the existence argument needs: a
        # positive definite matrix is a P-matrix, and a P-matrix LCP has exactly one solution
        # for every right-hand side. So the problem is well posed and only the algorithm had to
        # change, which is why this is asserted rather than left as a remark.
        cfg = nv_cfg(); ob = OBParams(0.0, 1.0)
        _, _, _, lower = nv_compliance(cfg, ob, [nv_init()], cfg.dt_max)
        A_c, _, _, _ = nv_compliance(cfg, ob, [nv_init()], cfg.dt_max)
        @test lower == collect(1:length(lower))            # contiguous, from the south pole
        @test all(cos(cfg.theta_vec[i]) < 0 for i in lower)
        As = A_c[lower, lower]
        @test maximum(abs, As - As') / maximum(abs, As) > 0.1     # measured 0.56: asymmetric
        @test minimum(eigvals(Symmetric(0.5*(As + As')))) > 0     # yet positive definite
        ## and the full matrix is NOT usable: its diagonal changes sign at the equator,
        ## because pushing on the crown of the drop opens no gap underneath it
        @test any(<(0), diag(A_c))
        @test all(>(0), diag(A_c)[lower])
    end

    @testset "complementarity actually holds at every step" begin
        # The defining property, and the reason this closure needs no candidate ranking:
        # at each accepted step the pressure is non-negative, the gap is non-negative, and
        # they are never both positive at the same node. A search-based closure satisfies the
        # first two by construction and the third only if its chosen contact count happens to
        # be right.
        #
        # Checked with a TWO-SIDED measure. `max_i min(h_i, p_i)` is small whenever one of the
        # pair is small at each node, so it certifies states in which the drop has penetrated
        # the substrate -- it passed a run whose gap reached -1.7e-3 while reporting success.
        cfg = nv_cfg(); ob = OBParams(0.0, 1.0)
        _, _, _, lower = nv_compliance(cfg, ob, [nv_init()], cfg.dt_max)
        _, states = solve_drop_lcp!(cfg, ob, nv_init(); t_end = 6.0, save_every = 0.05)
        @test length(states) > 50
        for s in states
            p = nv_nodal_pressure(cfg, s)[lower]
            g = nv_gap(cfg, s)[lower]
            @test minimum(p) > -1e-8                        # the wall pushes, never pulls
            @test minimum(g) > -1e-8                        # and the drop never enters it
            @test maximum(abs.(p .* g)) < 1e-8              # pressure only where contact is
        end
    end

    @testset "it agrees with the ranked search where the search is reliable" begin
        # Two closures over the same formulation. They choose the contact set by unrelated
        # criteria -- one ranks candidates by an edge residual, the other solves an
        # inequality-constrained problem with no candidates at all -- so agreement is evidence
        # about the formulation rather than about either closure.
        cfg = nv_cfg(; Oh = 0.05, Bo = 0.01); ob = OBParams()
        _, ss = solve_drop!(    cfg, ob, nv_init(NV_M; z = 1.05, v = -0.4); t_end = 15.0, save_every = 0.05)
        _, sl = solve_drop_lcp!(cfg, ob, nv_init(NV_M; z = 1.05, v = -0.4); t_end = 15.0, save_every = 0.05)
        ## contact interval, measured the same way on both
        span(states, times) = begin
            i = findfirst(s -> s.cp > 0, states); j = findlast(s -> s.cp > 0, states)
            (i === nothing || j === nothing) ? NaN : times[j] - times[i]
        end
        ts, _ = solve_drop!(    cfg, ob, nv_init(NV_M; z=1.05, v=-0.4); t_end=15.0, save_every=0.05)
        tl, _ = solve_drop_lcp!(cfg, ob, nv_init(NV_M; z=1.05, v=-0.4); t_end=15.0, save_every=0.05)
        @test isapprox(span(ss, ts), span(sl, tl); rtol = 0.15)
        ## and the rebound speed
        vs = ss[end].v; vl = sl[end].v
        @test vs > 0 && vl > 0
        @test isapprox(vs, vl; rtol = 0.15)
    end
end
