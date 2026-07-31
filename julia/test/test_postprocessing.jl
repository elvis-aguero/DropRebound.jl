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

@testset "extract_kpis: multi-bounce trajectory uses ONLY the first contact segment" begin
    # Regression: a drop that bounces, falls back under gravity, and
    # re-contacts within the simulated window (routine at low We with a
    # shear-thinning fluid -- a small rebound velocity plus gravity easily
    # produces a second impact before t_end) must have cor/contact_time
    # computed from the FIRST bounce only. Using contact_idx[1]:contact_idx[end]
    # (spanning every contact anywhere in the trajectory) silently mixes in
    # the free-fall speed regained between bounces and can even report
    # cor > 1 -- not an energy-conservation bug, a measurement-extraction one.
    cfg, dt_max = _pp_setup()
    M = cfg.M

    mk(cp, v, z) = begin
        s = DropState(M)
        s.cp = cp; s.v = v; s.z = z
        s
    end
    # idx: 1        2       3(c)    4(c)    5        6        7(c)     8
    #      airborne airborne first bounce   airborne airborne 2nd bounce airborne
    states = [
        mk(0, -1.0, 1.5),   # 1: pre-contact (v_in for bounce 1)
        mk(0, -1.0, 1.2),   # 2
        mk(1,  0.0, 1.0),   # 3: first contact segment starts
        mk(1,  0.3, 1.0),   # 4: first contact segment ends
        mk(0,  0.5, 1.05),  # 5: airborne after bounce 1 (v_out for bounce 1)
        mk(0, -0.9, 1.02),  # 6: falling back (gravity), re-approaching
        mk(1,  0.0, 1.0),   # 7: SECOND, separate contact segment
        mk(0,  0.4, 1.03),  # 8: airborne after bounce 2
    ]
    # Non-uniform on purpose: bounce 1 (idx3->idx4) must span at least
    # MIN_BOUNCE_STEPS*dt_max to be accepted as a real bounce, and the
    # airborne gap before bounce 2 (idx4->idx7) must exceed dt_max so it is
    # a genuine separation, not a bridged gap.
    @assert dt_max > 0.02
    times = [0.0, 0.4, 0.8, 0.8 + 10.5 * dt_max, 2.2, 2.8, 2.8 + 2 * dt_max, 4.0]

    kpi = extract_kpis(times, states, cfg)

    @test kpi.contact_time ≈ times[4] - times[3]   # first segment only (idx 3-4), not idx 3-7
    v_in, z_in = states[2].v, states[2].z
    v_out, z_out = states[5].v, states[5].z
    expected_cor = sqrt(abs((0.5 * v_out^2 + cfg.Bo * (z_out - z_in)) / (0.5 * v_in^2)))
    @test kpi.cor ≈ expected_cor
    @test kpi.cor < 1.0   # this hand-built case is physically a genuine loss
end

@testset "extract_kpis: brief sub-dt_max cp=0 gap is bridged, not a real liftoff" begin
    # Regression: with the coupled multi-mode model, the very first
    # geometric touch of a still nearly-spherical drop is a degenerate
    # tangency -- a tiny elastic "kiss" that recedes for roughly one dt_max
    # (deformation still at the machine-zero floor, momentum essentially
    # unchanged) before the real, dynamically significant bounce begins a
    # few frames later. Treating that brief recession as the end of the
    # first bounce collapsed contact_time to near-zero (median tc error 62%
    # in the 20-sample validation sweep). A genuine separation-and-reimpact
    # requires the whole drop to retreat and fall back under gravity, which
    # cannot happen within one dt_max, so a cp=0 gap no longer than dt_max
    # must be bridged; a gap longer than dt_max still ends the segment
    # (exercised by the multi-bounce test above, whose gap is 0.3 >> dt_max
    # for this M=6 setup).
    cfg, dt_max = _pp_setup()
    M = cfg.M

    mk(cp, v, z) = begin
        s = DropState(M)
        s.cp = cp; s.v = v; s.z = z
        s
    end
    # idx: 1        2       3(c)  4(gap) 5(c)  6(c)   7        8
    states = [
        mk(0, -1.0, 1.5),   # 1: pre-contact (v_in)
        mk(0, -1.0, 1.2),   # 2
        mk(1,  0.0, 1.0),   # 3: contact begins
        mk(0,  0.05, 1.0),  # 4: sub-dt_max cp=0 gap -- must be bridged
        mk(2,  0.1, 1.0),   # 5: contact resumes
        mk(1,  0.3, 1.0),   # 6: contact segment ends here
        mk(0,  0.5, 1.05),  # 7: airborne after bounce (v_out)
        mk(0,  0.45, 1.06), # 8
    ]
    # Non-uniform on purpose: the idx3->idx5 gap must stay well under
    # dt_max to be bridged, while the whole segment (idx3->idx6) must span
    # at least MIN_BOUNCE_STEPS*dt_max to be accepted as a real bounce.
    @assert dt_max > 0.02
    times = [0.0, 0.4, 0.8, 0.8 + 0.4 * dt_max, 0.8 + 0.8 * dt_max, 0.8 + 10.5 * dt_max, 2.2, 2.8]

    kpi = extract_kpis(times, states, cfg)

    @test kpi.contact_time ≈ times[6] - times[3]   # bridges idx4, ends at idx6
    v_in, v_out = states[2].v, states[7].v
    z_in, z_out = states[2].z, states[7].z
    expected_cor = sqrt(abs((0.5 * v_out^2 + cfg.Bo * (z_out - z_in)) / (0.5 * v_in^2)))
    @test kpi.cor ≈ expected_cor
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
