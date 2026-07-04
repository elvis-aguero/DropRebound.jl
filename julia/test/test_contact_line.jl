using Test
using DropSolver

# Milestone 1: finite equilibrium contact angle enters ONLY through the contact-set
# selection residual. See docs/DropRebound_ContactLine.tex.

@testset "Contact line: C1 undeformed sphere θ_d = θ_c" begin
    # For zero deformation, the apparent contact angle equals the polar contact angle.
    s = DropState(6)   # all amplitudes zero
    for θ in range(π/2 + 0.05, π - 0.05; length = 9)
        @test isapprox(contact_angle(s, θ), θ; atol = 1e-10)
    end
end

@testset "Contact line: C1 with deformation is a small perturbation of θ_c" begin
    s = DropState(6)
    s.A[2] = 1e-4
    θ = 3π/4
    @test isapprox(contact_angle(s, θ), θ; atol = 5e-3)   # O(A) correction
end

@testset "Contact line: C8 residual reduces to GA tangency at θ_e = π" begin
    # contact_angle_error(...,π) must equal |contact_error(...)| for arbitrary states.
    M = 6
    θv = make_theta_vec(M)
    for trial in 1:20
        s = DropState(M)
        # deterministic pseudo-amplitudes (no RNG dependence)
        for n in 2:M
            s.A[n] = 0.02 * sin(1.7 * n + 0.3 * trial)
        end
        s.z = 0.9 + 0.01 * trial
        for cp in 1:M
            ga  = abs(contact_error(s, θv, cp))
            cl  = contact_angle_error(s, θv, cp, π)
            @test isapprox(cl, ga; atol = 1e-12, rtol = 1e-10)
        end
        @test contact_angle_error(s, θv, 0, π) == 0.0   # cp≤0 short-circuit
    end
end

@testset "Contact line: exact reduction of solve_drop! at θ_e=π (CLParams default)" begin
    # A full impact trajectory with the default CLParams() must be bit-for-bit identical
    # to the baseline (no cl kwarg): θ_e=π ⇒ is_cl=false ⇒ contact_error path.
    M = 6; Oh = 0.3038; Bo = 1/53.9
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.1; init.v = -0.281; init.dt = dt_max; init.cp = 0

    _, base = solve_drop!(cfg, OBParams(), deepcopy(init); t_end=8.0, save_every=0.05)
    _, cl0  = solve_drop!(cfg, OBParams(), deepcopy(init); cl=CLParams(), t_end=8.0, save_every=0.05)

    @test length(base) == length(cl0)
    for (a, b) in zip(base, cl0)
        @test a.z  == b.z
        @test a.v  == b.v
        @test a.cp == b.cp
        @test a.A  == b.A
        @test a.B  == b.B
    end
end

@testset "Contact line: finite θ_e runs and can change contact behaviour" begin
    # A more wetting substrate (θ_e < π) should still produce a valid impact, and the
    # contact history need not be identical to the non-wetting case.
    M = 6; Oh = 0.3038; Bo = 1/53.9
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.1; init.v = -0.281; init.dt = dt_max; init.cp = 0

    _, wet = solve_drop!(cfg, OBParams(), deepcopy(init);
                         cl=CLParams(0.85π, 0.0), t_end=8.0, save_every=0.05)

    @test all(isfinite(s.z)    for s in wet)
    @test all(isfinite(s.A[2]) for s in wet)
    @test any(s.cp > 0 for s in wet)                 # contact still occurs
end

# ── Milestone 2: contact-line mobility ξ (continuous tracker) ────────────────────
@testset "Contact line: M2 mobility limits and monotonicity" begin
    M = 12; Oh = 0.1; Bo = 0.01
    cfg = SimConstants(M, M+1, Oh, Bo, make_theta_vec(M),
                       precompute_integrals(NaN, M)[1], make_dt_max(M))
    base() = (i = DropState(M); i.z = 1.05; i.v = -0.5; i.dt = make_dt_max(M); i.cp = 0; i)
    cor(θe, ξ) = begin
        t, s = solve_drop!(cfg, OBParams(), base(); cl = CLParams(θe, ξ),
                           t_end = 12.0, save_every = 0.05)
        extract_kpis(t, s, cfg).cor
    end

    cor_nonwet = cor(π,     0.0)     # perfectly non-wetting reference
    cor_qs     = cor(0.90π, 0.0)     # wetting, quasi-static (Milestone 1)
    cor_eps    = cor(0.90π, 1e-8)    # wetting, near-zero friction
    cor_pin    = cor(0.90π, 50.0)    # wetting, strong pinning

    # All physical: no energy injection into the ballistic COM motion.
    for c in (cor_nonwet, cor_qs, cor_eps, cor_pin)
        @test isfinite(c) && c <= 1.0 + 1e-6
    end
    # ξ→0 limit: the mobility tracker reduces to the quasi-static selection.
    @test isapprox(cor_eps, cor_qs; atol = 1e-3)
    # Wetting lowers restitution relative to non-wetting.
    @test cor_qs < cor_nonwet - 1e-2
    # ξ→∞ limit: strong pinning suppresses wetting capture, recovering non-wetting COR.
    @test isapprox(cor_pin, cor_nonwet; atol = 2e-2)
    # Monotone effect of friction: pinning raises restitution above the mobile case.
    @test cor_pin > cor_qs + 1e-2
end

# ── Event-located time stepping (opt-in; docs §Discrete) ─────────────────────────
@testset "Event location: backwards compatible + agrees with legacy" begin
    M = 14; Oh = 0.03; Bo = 0.02
    cfg = SimConstants(M, M+1, Oh, Bo, make_theta_vec(M),
                       precompute_integrals(NaN, M)[1], make_dt_max(M))
    mk(We) = (i = DropState(M); i.z = 1.0; i.v = -sqrt(We);
              i.dt = make_dt_max(M); i.cp = 0; i)

    # Bit-for-bit: event_location=false must reproduce the default (legacy) path.
    # Clear the module-level Jacobian cache before each run so a warm-cache inv(A)*b
    # vs cold A\b (machine-precision) difference doesn't masquerade as a discrepancy.
    clear_jac_cache!()
    t0, s0 = solve_drop!(cfg, OBParams(), mk(0.253); t_end = 8.0, save_every = 0.05)
    clear_jac_cache!()
    t1, s1 = solve_drop!(cfg, OBParams(), mk(0.253); t_end = 8.0, save_every = 0.05,
                         event_location = false)
    @test t0 == t1
    @test all(a.z == b.z && a.cp == b.cp && a.A == b.A for (a, b) in zip(s0, s1))

    # Event-located stepping agrees with legacy on the robust KPIs (contact time, COR).
    # (Spreading time is validated at production resolution in scripts/run_event_validation.jl;
    # its broad, flat maximum makes it too mesh-sensitive to gate at test resolution.)
    for We in (0.253, 2.269)
        o = extract_kpis(solve_drop!(cfg, OBParams(), mk(We);
                                     t_end = 12.0, save_every = 0.01)..., cfg)
        n = extract_kpis(solve_drop!(cfg, OBParams(), mk(We);
                                     t_end = 12.0, save_every = 0.01,
                                     event_location = true)..., cfg)
        @test isfinite(n.contact_time) && isfinite(n.cor)
        @test abs(n.contact_time - o.contact_time) / o.contact_time < 0.04
        @test abs(n.cor - o.cor) / o.cor < 0.04
        @test n.cor <= 1.0 + 1e-6
    end

    # Guard: event_location with finite friction (ξ>0) must fall back to legacy
    # stepping (the geometric event is not the controlling transition), producing the
    # same result as legacy — not a dt-floor crash.
    clear_jac_cache!()
    tf, sf = solve_drop!(cfg, OBParams(), mk(0.253); cl = CLParams(0.9π, 1.0),
                         t_end = 10.0, save_every = 0.05)
    clear_jac_cache!()
    te, se = solve_drop!(cfg, OBParams(), mk(0.253); cl = CLParams(0.9π, 1.0),
                         t_end = 10.0, save_every = 0.05, event_location = true)
    @test tf == te
    @test all(a.z == b.z && a.cp == b.cp for (a, b) in zip(sf, se))
end
