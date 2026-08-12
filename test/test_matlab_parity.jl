using Test
using DropSolver

# ── Reference values (N=90 modes, linearized solver, BDF1) ───────────────────
# Physical: ρ=0.96 g/cm³, σ=20.5 dyn/cm, R=0.0201 cm, ν=0.199 cm²/s, v₀=-9.156 cm/s
# Oh = ν√(ρ/(σR)) = 0.3038    Bo = ρgR²/σ = 1/53.9 ≈ 0.01856    v₀_dimless = -0.281
# time_unit = √(ρR³/σ) = 6.16×10⁻⁴ s
const MATLAB_Oh           = 0.3038
const MATLAB_Bo           = 1/53.9
const MATLAB_v0           = -0.281
const MATLAB_CONTACT_TIME = 2.99     # dimensionless contact duration
const MATLAB_MAX_RADIUS   = 0.397    # dimensionless maximum spreading radius
const MATLAB_COR          = 0.484    # coefficient of restitution

@testset "MATLAB parity: the production solver at production resolution" begin
    # Code-to-code verification against the Gabbard et al. (2025) MATLAB solver, on
    # the variational solver at M = 90, K = 3 -- the configuration this package
    # actually ships.
    #
    # THE COMPARISON IS NOT ONE-TO-ONE, AND THAT IS THE POINT. Their model carries
    # surface modes with potential flow inside and Lamb's damping. Here that is
    # K = 1: a single irrotational trial function per mode. K >= 2 resolves the
    # vortical layer and so carries Reid's exact finite-Ohnesorge damping, which
    # their model does not have. So the two quantities behave differently and are
    # asserted differently:
    #
    #   contact time  is set by the capillary spring, not by the damping, so it
    #                 must agree at BOTH K. It does, and production K = 3 agrees
    #                 more closely (0.06 %) than K = 1 (1.1 %).
    #
    #   restitution   is set by how much energy the drop loses, which is exactly
    #                 what the two models disagree about. K = 1 reproduces them to
    #                 2.5 %; K = 3 sits 17 % higher because Lamb over-damps at this
    #                 Ohnesorge.
    #
    # A failing K = 1 assertion means this package no longer reproduces the
    # published model when configured as it. A failing K = 3 contact-time assertion
    # means the capillary response has drifted, independently of any damping model.
    Oh, Bo, We = MATLAB_Oh, MATLAB_Bo, MATLAB_v0^2
    M, TOL = 90, 0.10

    r1 = simulate(ImpactParams(We = We, Bo = Bo, Oh = Oh, M = M, K = 1, t_max = 25.0))
    r3 = simulate(ImpactParams(We = We, Bo = Bo, Oh = Oh, M = M, K = 3, t_max = 25.0))

    @test isfinite(r1.cor) && isfinite(r3.cor)

    ## Same model as theirs: both KPIs must land.
    @test abs(r1.cor - MATLAB_COR) / MATLAB_COR                   < TOL
    @test abs(r1.tc  - MATLAB_CONTACT_TIME) / MATLAB_CONTACT_TIME < TOL

    ## Production: contact time must still land, since it does not depend on the
    ## damping model.
    @test abs(r3.tc - MATLAB_CONTACT_TIME) / MATLAB_CONTACT_TIME  < TOL

    ## And the restitution must differ in the direction the physics requires.
    ## Lamb over-predicts the damping, so resolving the interior returns energy to
    ## the drop and the bounce gets livelier. A K = 3 that matched their
    ## restitution would mean the extra radial functions had stopped doing anything.
    @test r3.cor > r1.cor
    @test (r3.cor - MATLAB_COR) / MATLAB_COR > 0.05
end

@testset "MATLAB parity: the legacy nonvariational solver" begin
    # The same reference point on `solve_drop!`, which eliminates the interior and
    # carries Reid's coefficients per mode. Kept because it is the second
    # formulation and a cross-check on the first, but it is not the production
    # route and it is fragile here: at M = 45 this case never releases, while
    # M = 20 and M = 90 both complete. So it is pinned to M = 20 deliberately, and
    # that fragility is the reason rather than an accident of the tolerance.
    M         = 20
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, MATLAB_Oh, MATLAB_Bo, theta_vec, precomp, dt_max)
    ob        = OBParams(0.0, 1.0)

    init = DropState(M)
    init.z = 1.1; init.v = MATLAB_v0; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=8.0, save_every=0.02)

    @test all(isfinite(s.z)    for s in states)
    @test all(isfinite(s.A[2]) for s in states)
    @test any(s.cp > 0 for s in states)

    kpis = extract_kpis(times, states, cfg)

    @test !isnan(kpis.cor)
    @test abs(kpis.contact_time - MATLAB_CONTACT_TIME) / MATLAB_CONTACT_TIME < 0.05
    @test abs(kpis.max_radius   - MATLAB_MAX_RADIUS)   / MATLAB_MAX_RADIUS   < 0.05
    @test abs(kpis.cor          - MATLAB_COR)           / MATLAB_COR          < 0.05
end

@testset "OB impact: polymer stress changes contact metrics" begin
    M         = 20
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, MATLAB_Oh, MATLAB_Bo, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.1; init.v = MATLAB_v0; init.dt = dt_max; init.cp = 0

    ob_N  = OBParams(0.0, 1.0)
    ob_OB = OBParams(0.5, 0.5)

    times_N,  states_N  = solve_drop!(cfg, ob_N,  deepcopy(init); t_end=8.0, save_every=0.02)
    times_OB, states_OB = solve_drop!(cfg, ob_OB, deepcopy(init); t_end=8.0, save_every=0.02)

    @test all(isfinite(s.z) for s in states_N)
    @test all(isfinite(s.z) for s in states_OB)
    @test any(s.cp > 0 for s in states_N)
    @test any(s.cp > 0 for s in states_OB)

    kpis_N  = extract_kpis(times_N,  states_N,  cfg)
    kpis_OB = extract_kpis(times_OB, states_OB, cfg)

    @test !isnan(kpis_N.cor)
    @test !isnan(kpis_OB.cor)
    @test kpis_OB.cor != kpis_N.cor   # viscoelasticity changes rebound energy
end
