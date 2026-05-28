using Test
using DropSolver

# Shared minimal setup: M=6, Oh=0.1, Bo=1e-6 (gravity-free), free oscillation
function _pp_setup(; M=6, Oh=0.1, Bo=1e-6)
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    return cfg, dt_max
end

@testset "drop_profile: shape and normalization" begin
    cfg, dt_max = _pp_setup()

    # Spherical drop (all A=0): r(θ)=1, so xs=sinθ, zs=z+cosθ
    state = DropState(cfg.M)
    state.z = 2.0
    xs, zs = drop_profile(state, cfg; n_theta=100)

    @test length(xs) == 100
    @test length(zs) == 100
    @test all(isfinite, xs)
    @test all(isfinite, zs)
    # North pole (θ=0): xs[1]≈0, zs[1]≈z+1
    @test abs(xs[1])        < 1e-12
    @test abs(zs[1] - 3.0) < 1e-10
    # South pole (θ=π): xs[end]≈0, zs[end]≈z-1
    @test abs(xs[end])       < 1e-12
    @test abs(zs[end] - 1.0) < 1e-10
    # Max radius near equator (θ≈π/2): xs < 1 but close (discrete sampling)
    @test maximum(xs) > 0.99

    # With A₂ deformation: profile should deviate from sphere
    state2 = DropState(cfg.M)
    state2.z   = 2.0
    state2.A[2] = 0.1
    xs2, zs2 = drop_profile(state2, cfg; n_theta=100)
    @test !all(xs2 .≈ xs)
end

@testset "compute_contact_radius: zero when airborne, positive in contact" begin
    cfg, _ = _pp_setup()

    state_air = DropState(cfg.M)
    state_air.cp = 0
    @test compute_contact_radius(state_air, cfg) == 0.0

    # Manually set cp to first collocation point
    state_c = DropState(cfg.M)
    state_c.cp = 1
    r = compute_contact_radius(state_c, cfg)
    @test r > 0.0
    @test r < 1.0   # radius must be sub-unitary for small contact angles
end

@testset "extract_kpis: no contact → NaN / zeros" begin
    cfg, _ = _pp_setup()
    # Airborne trajectory: drop high up, short t_end so it never touches
    init = DropState(cfg.M)
    init.z  = 5.0
    init.v  = 0.0
    init.dt = make_dt_max(cfg.M)
    init.cp = 0

    times, states = solve_drop!(cfg, OBParams(0.0, 1.0), deepcopy(init);
                                 t_end=0.5, save_every=0.1)

    kpis = extract_kpis(times, states, cfg)
    @test isnan(kpis.contact_time)
    @test isnan(kpis.cor)
    @test kpis.max_radius == 0.0
    @test kpis.max_A2     == 0.0
end

@testset "extract_kpis: impact run produces valid KPIs" begin
    M         = 20
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, 0.3038, 1/53.9, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.1; init.v = -0.281; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, OBParams(0.0, 1.0), deepcopy(init);
                                 t_end=8.0, save_every=0.05)

    kpis = extract_kpis(times, states, cfg)

    @test kpis.contact_time > 0.0
    @test isfinite(kpis.contact_time)
    @test isfinite(kpis.cor)
    @test kpis.cor > 0.0
    @test kpis.max_radius > 0.0
    @test kpis.max_A2     > 0.0
end
