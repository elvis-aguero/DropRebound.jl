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

const BK_ALL = [Backend(),
                Backend(forcing = :nodal),
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
            r = run_impact(b; We = 0.5, Bo = 0.0188, Oh = 0.0373, M = 14, K = 2, t_max = 25.0)
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
        kw = (We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 20, K = 2, t_max = 25.0)
        p = ImpactParams(; kw...)
        @test run_impact(Backend(); kw...).cor ≈ proximity_metrics(p, simulate_lcp(p)).cor
        @test run_impact(Backend(contact = :active_set); kw...).cor ≈
              proximity_metrics(p, simulate(p)).cor
        pn = ImpactParams(; kw..., force_mode = :nodal)
        @test run_impact(Backend(forcing = :nodal); kw...).cor ≈
              proximity_metrics(pn, simulate_lcp(pn)).cor
    end

    @testset "a solver that gives up reports it instead of throwing" begin
        # A sweep has to be able to record which cases a backend cannot do. The nonvariational
        # search fails at the canonical reference point We = 1, Oh = 0.3038 -- which is worth
        # knowing on its own account, and is why the animations use a different case.
        r = run_impact(Backend(formulation = :nonvariational, contact = :tangency);
                       We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 20, K = 2, t_max = 25.0)
        @test r.ok == false
        @test isnan(r.cor)
        @test r.backend == "nonvar/tangency"
        ## misuse still throws, because it is a mistake rather than a result
        @test_throws ErrorException run_impact(Backend(); We = 1.0, Bo = 0.02, Oh = 0.3,
                                               M = 12, K = 1, eta_nonvar = 1)
    end

    @testset "the outline is the surface the solver produced" begin
        # `drop_outline` is what both the figure and the animation draw, and it is written
        # against the common form rather than against any solver, so it has to reproduce the
        # surface from `zeta` alone.
        r = run_impact(Backend(); We = 0.5, Bo = 0.0188, Oh = 0.0373, M = 14, K = 2,
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
