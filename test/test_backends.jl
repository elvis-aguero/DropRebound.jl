# The backend API.
#
# `Backend` names a formulation, a contact closure and a forcing; `run_impact` runs it and
# returns the same NamedTuple whichever was named. That uniformity is the whole point -- the
# comparison figure and the animations are written against it and never ask which solver
# produced a frame -- so it is what these tests check, rather than the physics, which is
# covered where each solver is tested.

using Test
using DropSolver
using LinearAlgebra

## Every cell, named explicitly. Writing one of these as `Backend()` couples the
## enumeration to whichever cell is currently the default, so changing the default
## silently turns two entries into duplicates instead of failing honestly.
const BK_ALL = [Backend(contact = :lcp),
                Backend(contact = :lcp,        forcing = :nodal),
                Backend(contact = :active_set),
                Backend(contact = :active_set, forcing = :nodal),
                Backend(formulation = :nonvariational, contact = :tangency)]

@testset "backend API" begin

    @testset "impossible combinations are refused, not substituted" begin
        # Silently falling back to something that does run is how a comparison ends up
        # comparing one solver with itself. Each of these names a cell that does not exist.
        @test_throws ErrorException Backend(formulation = :nonvariational, contact = :lcp)
        @test_throws ErrorException Backend(formulation = :nonvariational, forcing = :nodal)
        @test_throws ErrorException Backend(contact = :tangency)          # variational has none
        @test_throws ErrorException Backend(formulation = :quantum)
        @test_throws ErrorException Backend(contact = :magic)
        @test_throws ErrorException Backend(forcing = :magic)
        ## and the labels are distinct, since the store is keyed on them
        @test length(unique(label.(BK_ALL))) == length(BK_ALL)
    end

    @testset "every backend returns the same fields" begin
        # The property the figure and the animation both rely on.
        want = (:t, :z, :v, :zeta, :ls, :cp, :cor, :tc, :wall, :ok, :backend, :diag)
        wantdiag = (:min_gap, :n_detected, :rejects, :lcp_resid, :eta_sweeps,
                    :cor_internal, :tc_internal)
        for b in BK_ALL
            r = run_impact(b; We = 0.5, Bo = 0.0188, Oh = 0.0373, M = 30, K = 3, t_max = 25.0)
            @test all(k -> haskey(r, k), want)
            @test all(k -> haskey(r.diag, k), wantdiag)
            @test r.backend == label(b)
            if r.ok
                @test length(r.t) == length(r.z) == length(r.v) == length(r.zeta) == length(r.cp)
                @test length(r.zeta[1]) == length(r.ls)
                @test isfinite(r.cor) && isfinite(r.tc)
            end
            @test r.wall > 0
        end
    end

    @testset "the backends agree with the functions they wrap" begin
        # `run_impact` must not be a second implementation. Same parameters, same numbers.
        kw = (We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 30, K = 3, t_max = 25.0)
        p = ImpactParams(; kw...)
        @test run_impact(Backend(contact = :lcp); kw...).cor ≈
              proximity_metrics(p, simulate_lcp(p)).cor
        @test run_impact(Backend(contact = :active_set); kw...).cor ≈
              proximity_metrics(p, simulate(p)).cor
        pn = ImpactParams(; kw..., force_mode = :nodal)
        @test run_impact(Backend(contact = :lcp, forcing = :nodal); kw...).cor ≈
              proximity_metrics(pn, simulate_lcp(pn)).cor
        ## and the default resolves to a cell that exists
        @test label(Backend()) in label.(BK_ALL)
    end

    @testset "a solver that gives up reports it instead of throwing" begin
        # A sweep has to be able to record which cases a backend cannot do, so a run that
        # fails must come back as data rather than as an exception. That contract is what
        # is tested here, in two halves: `check_converged` decides, and `run_impact`
        # reports the decision without throwing.
        #
        # WHY THE VERDICT IS NOT PINNED TO A CASE. This test used to assert that the
        # nonvariational search fails at We = 1, Oh = 0.3038, M = 30. It does, on some
        # machines. On CI it failed on 10 August and completed a bounce on 11 August with
        # byte-identical sources, the same Julia 1.12.6 and the same runner image; on macOS
        # it returns NaN outright. The march there is unstable, so whether it survives turns
        # on floating-point details that are not reproducible between runs, and a test that
        # asserts the outcome of an unstable march tests the machine rather than the code.
        #
        # So the decision logic is tested directly, where it is deterministic.
        @test DropSolver.check_converged(NaN, 1.0, 25.0, 0.5, "t") == false      # non-finite
        @test DropSolver.check_converged(0.5, 24.0, 25.0, 0.5, "t") == false     # never released
        @test DropSolver.check_converged(1e-18, 2.0, 25.0, 0.5, "t") == false    # no rebound
        @test DropSolver.check_converged(1.5, 2.0, 25.0, 0.5, "t") == false      # cor > 1
        @test DropSolver.check_converged(0.5, 2.0, 25.0, 0.0, "t") == false      # hit the wall
        @test DropSolver.check_converged(0.5, 2.0, 25.0, 0.5, "t"; released = false) == false
        @test DropSolver.check_converged(0.5, 2.0, 25.0, 0.5, "t"; dt_final = 1e-9,
                                         dt_min = 1e-9) == false                 # step collapsed
        @test DropSolver.check_converged(0.5, 2.0, 25.0, 0.5, "t") == true       # a real bounce

        # And the reporting contract is tested on the marginal case, whichever way it goes:
        # it must return, it must say which backend ran, and its metrics must agree with its
        # own verdict -- usable when `ok`, NaN when not. Both branches are legitimate here.
        r = run_impact(Backend(formulation = :nonvariational, contact = :tangency);
                       We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 30, K = 3, t_max = 25.0)
        @test r.backend == "nonvar/tangency"
        @test r.ok isa Bool
        @test r.ok ? isfinite(r.cor) : isnan(r.cor)

        ## misuse still throws, because it is a mistake rather than a result
        @test_throws ErrorException run_impact(Backend(); We = 1.0, Bo = 0.02, Oh = 0.3,
                                               M = 12, K = 1, eta_nonvar = 1)
    end

    @testset "the outline is the surface the solver produced" begin
        # `drop_outline` is what both the figure and the animation draw, and it is written
        # against the common form rather than against any solver, so it has to reproduce the
        # surface from `zeta` alone.
        r = run_impact(Backend(); We = 0.5, Bo = 0.0188, Oh = 0.0373, M = 30, K = 3,
                       t_max = 25.0)
        @test r.ok
        i = argmax([maximum(abs, z) for z in r.zeta])
        x, y = drop_outline(r.zeta[i], r.ls, r.z[i]; nth = 200)
        @test length(x) == length(y) == 200
        @test all(isfinite, x) && all(isfinite, y)
        ## The outline runs theta = 0 to pi, so its LAST point is the south pole and must sit
        ## at the centre of mass minus the polar radius. Checked at the pole rather than at the
        ## lowest point of the outline: at peak deformation the drop is flattened and its
        ## lowest point is a ring, not the pole -- an earlier version of this test assumed
        ## otherwise and failed for that reason rather than for a fault in the code.
        rpole = 1.0 + sum(r.zeta[i][k] * DropSolver.legendre_angular(l, -1.0).P
                          for (k, l) in enumerate(r.ls))
        @test isapprox(y[end], r.z[i] - rpole; atol = 1e-8)
        @test isapprox(x[end], 0.0; atol = 1e-12)
        ## and the flattening is real: the lowest point of the outline is at or below the pole
        @test minimum(y) <= y[end] + 1e-12
        ## an undeformed drop is a unit circle about z
        x0, y0 = drop_outline(zeros(length(r.ls)), r.ls, 2.0; nth = 200)
        @test maximum(abs, sqrt.(x0.^2 .+ (y0 .- 2.0).^2) .- 1.0) < 1e-12
    end
end
