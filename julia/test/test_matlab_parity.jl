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

@testset "MATLAB parity: Newtonian KPIs at M=20 (<5% error)" begin
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
