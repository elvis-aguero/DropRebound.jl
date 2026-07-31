using Test, LinearAlgebra
using DropSolver

@testset "Non-perturbative Carreau-Yasuda (STExactParams)" begin

    @testset "characteristic_shear_K matches the exact-rational derivation values" begin
        # julia/derivations/carreau_yasuda_nonperturbative_derivation.jl Section 2:
        # K_l^2 = 3, 4, 9/2, 24/5, 5, 36/7, 21/4, 16/3, 27/5 for l=2..10.
        expected = Dict(2 => 3.0, 3 => 4.0, 4 => 4.5, 5 => 4.8, 6 => 5.0,
            7 => 36 / 7, 8 => 5.25, 9 => 16 / 3, 10 => 5.4)
        for (l, K2) in expected
            @test characteristic_shear_K(l)^2 ≈ K2 atol=1e-10
        end
    end

    @testset "lambda_c=0 -> identical to build_residual!/build_jacobian" begin
        M = 4
        Oh0 = 0.05
        Bo = 1e-12
        theta_vec = collect(range(pi, 0; length=M+1))
        precomp = precompute_integrals(NaN, M)[1]
        cfg = SimConstants(M, M + 1, Oh0, Bo, theta_vec, precomp, 0.005)
        ob = OBParams()

        stx = STExactParams(M, Oh0, 0.0, 2.0, 0.5; viscous=:lamb)  # lambda_c=0 -> gammadot_char irrelevant

        s0 = DropState(M)
        s0.A[2] = 0.05
        s0.z = 1.0
        s0.dt = 0.005

        R_N = zeros(3M + 1)
        build_residual!(R_N, s0, [s0], 0.005, 0, cfg, ob)
        R_X = zeros(3M + 1)
        build_residual_st_exact!(R_X, s0, [s0], 0.005, 0, cfg, ob, stx)
        @test R_X ≈ R_N atol=1e-12

        J_N = build_jacobian(s0, [s0], 0.005, 0, cfg, ob)
        J_X = build_jacobian_st_exact(s0, [s0], 0.005, 0, cfg, ob, stx)
        @test J_X ≈ J_N atol=1e-12
    end

    @testset "Jacobian matches finite differences (:lamb)" begin
        M = 4
        Oh0 = 0.3
        Bo = 1e-6
        theta_vec = collect(range(pi, 0; length=M+1))
        precomp = precompute_integrals(NaN, M)[1]
        cfg = SimConstants(M, M + 1, Oh0, Bo, theta_vec, precomp, 0.01)
        ob = OBParams()
        stx = STExactParams(M, Oh0, 5.0, 1.5, 0.3; viscous=:lamb)

        s_prev = DropState(M)
        s_prev.Adot[2] = 0.2
        s_prev.z = 1.0
        s_prev.dt = 0.01
        history = [s_prev]
        dt = 0.01
        cp = 0

        s_curr = DropState(M)
        s_curr.A[2] = 0.04
        s_curr.Adot[2] = 0.18
        s_curr.z = 1.02
        s_curr.v = -0.1
        s_curr.B .= 0.01

        X0 = pack_X(s_curr, M)
        function R!(buf, Xv)
            s = deepcopy(s_curr)
            unpack_X!(s, Xv, M)
            build_residual_st_exact!(buf, s, history, dt, cp, cfg, ob, stx)
        end
        n = length(X0)
        R0 = zeros(n)
        R!(R0, X0)
        J_analytic = build_jacobian_st_exact(s_curr, history, dt, cp, cfg, ob, stx)

        J_fd = zeros(n, n)
        h = 1e-6
        for j in 1:n
            Xp = copy(X0)
            Xp[j] += h
            Rp = zeros(n)
            R!(Rp, Xp)
            J_fd[:, j] .= (Rp .- R0) ./ h
        end
        @test J_analytic ≈ J_fd rtol=1e-4 atol=1e-6
    end

    @testset "STExactParams(:reid) construction" begin
        M = 6
        Oh0 = 57.4
        stx = STExactParams(M, Oh0, 30507.0, 0.7431, 0.99956; viscous=:reid)
        @test stx.viscous === :reid
        @test stx.reid_cache !== nothing
        @test length(stx.reid_cache) == M - 1
        @test length(stx.e_rr) == M - 1
    end

    @testset "oh_eff_all_coupled: single active mode is close to characteristic_shear_K" begin
        # These are two DIFFERENT calculations, not the same one two ways: the
        # single-mode formula evaluates mu_eff ONCE at a single characteristic
        # (volume-RMS) shear rate, while oh_eff_all_coupled evaluates mu_eff
        # POINTWISE and averages the RATIO over the whole volume. For a
        # nonlinear (here concave/saturating) mu_eff, these need not agree
        # (Jensen's inequality) except where the mode's own shear field is
        # spatially uniform enough that "the RMS value" is representative.
        # Measured: exact for l=2 (ratio=1.0000), and a growing gap for higher
        # l (l=6: ratio=0.819) as the mode's own field becomes more spatially
        # oscillatory -- expected, not a bug: this is exactly why the coupled,
        # volume-averaged calculation is the more accurate one to use in
        # production, not merely an equivalent reformulation.
        # Constructed with headroom (M=10, only exciting l up to 6) so that
        # dealiasing_cutoff(10)=9 doesn't filter out any mode this test
        # actually exercises -- see the dedicated dealiasing test below for
        # why that filter exists and what M=6 alone would have excluded.
        M = 10
        L_MAX = 6
        Oh0 = 57.4
        lambda_c = 30507.0
        a = 0.7431
        eps_ST = 0.99956
        stx = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid)

        @test isapprox(begin
            Adot_vec = zeros(M - 1); Adot_vec[1] = 0.03
            oh_eff_all_coupled(stx, Oh0, Adot_vec)[1]
        end, begin
            Kl = characteristic_shear_K(2)
            Oh0 * (1 + (lambda_c * Kl * 0.03)^a)^(-eps_ST)
        end; rtol=1e-6)   # l=2 matches essentially exactly

        prev_ratio = 1.0
        for l in 3:L_MAX
            k = l - 1
            Adot_l = 0.03
            Kl = characteristic_shear_K(l)
            Oh_eff_expected = Oh0 * (1 + (lambda_c * Kl * Adot_l)^a)^(-eps_ST)
            Adot_vec = zeros(M - 1)
            Adot_vec[k] = Adot_l
            Oh_eff = oh_eff_all_coupled(stx, Oh0, Adot_vec)[k]
            ratio = Oh_eff / Oh_eff_expected
            @test 0.5 < ratio < 1.0   # coupled result is smaller, but same order of magnitude
            @test ratio < prev_ratio  # gap grows monotonically with l, as expected
            prev_ratio = ratio
        end
    end

    @testset "oh_eff_all_coupled: Newtonian limit is Oh0 exactly regardless of coupling" begin
        M = 6
        Oh0 = 3.0
        stx = STExactParams(M, Oh0, 0.0, 2.0, 0.5; viscous=:lamb)  # lambda_c=0 -> no thinning anywhere
        for Adot_vec in ([0.3, 0.0, 0.0, 0.0, 0.0], [0.3, 0.0, 0.6, -0.2, 0.0])
            Oh_eff = oh_eff_all_coupled(stx, Oh0, Adot_vec)
            @test all(x -> isapprox(x, Oh0; atol=1e-8), Oh_eff)
        end
    end

    @testset "oh_eff_all_coupled: exciting a second mode enhances thinning of the first" begin
        # The actual point of the multi-mode fix: mode 2's effective Oh must
        # drop FURTHER once another mode (e.g. l=5, as contact would excite)
        # becomes active alongside it, versus mode 2 acting alone.
        M = 6
        Oh0 = 57.4
        stx = STExactParams(M, Oh0, 30507.0, 0.7431, 0.99956; viscous=:reid)

        Adot2 = 0.02
        base = zeros(M - 1)
        base[1] = Adot2   # mode 2 alone
        Oh_eff_self_only = oh_eff_all_coupled(stx, Oh0, base)[1]

        prev = Oh_eff_self_only
        for Adot5 in (0.02, 0.05, 0.1)
            vec = zeros(M - 1)
            vec[1] = Adot2
            vec[4] = Adot5   # mode 5 (index 4 = l-1 for l=5)
            Oh_eff_2 = oh_eff_all_coupled(stx, Oh0, vec)[1]
            @test Oh_eff_2 < prev   # thinning enhanced monotonically with the second mode's amplitude
            prev = Oh_eff_2
        end
        @test prev < Oh_eff_self_only
    end

    @testset "oh_eff_all_coupled: eta_inf_ratio floors Oh_eff, defaults to no floor" begin
        # Regression: with eta_inf_ratio=0.0 (the default), extreme shear can
        # drive Oh_eff arbitrarily close to 0. A real Cross-model fluid's
        # infinite-shear viscosity is never exactly zero, so Oh_eff must never
        # drop below Oh0*eta_inf_ratio once that ratio is supplied.
        M = 6
        Oh0 = 57.4
        lambda_c = 30507.0
        a = 0.7431
        eps_ST = 0.99956
        huge_shear = fill(1e6, M - 1)   # enormous Adot -> mu_eff/mu_0 saturates near its floor

        stx_no_floor = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid)
        Oh_eff_no_floor = oh_eff_all_coupled(stx_no_floor, Oh0, huge_shear)
        @test all(x -> x < 1e-6 * Oh0, Oh_eff_no_floor)   # no floor -> collapses toward 0

        eta_inf_ratio = 0.00044
        stx_floor = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid, eta_inf_ratio=eta_inf_ratio)
        Oh_eff_floor = oh_eff_all_coupled(stx_floor, Oh0, huge_shear)
        @test all(x -> x >= Oh0 * eta_inf_ratio - 1e-8, Oh_eff_floor)
        @test all(x -> isapprox(x, Oh0 * eta_inf_ratio; rtol=1e-3), Oh_eff_floor)   # saturates AT the floor

        @test_throws ArgumentError STExactParams(M, Oh0, lambda_c, a, eps_ST; eta_inf_ratio=1.0)
        @test_throws ArgumentError STExactParams(M, Oh0, lambda_c, a, eps_ST; eta_inf_ratio=-0.1)
    end

    @testset "_ringing_outlier_mask: flags only a genuine >factor amplitude outlier" begin
        # A single mode 100x+ every other active mode -> flagged (the
        # truncation-ringing signature). A single active mode with everything
        # else exactly zero -> NOT flagged (the vacuous-outlier guard).
        # Comparable-magnitude modes -> none flagged.
        mask1 = DropSolver._ringing_outlier_mask([1e-15, 1e-15, 1e-15, 0.673])
        @test mask1 == [false, false, false, true]

        mask2 = DropSolver._ringing_outlier_mask([0.03, 0.0, 0.0, 0.0])
        @test mask2 == [false, false, false, false]

        mask3 = DropSolver._ringing_outlier_mask([-0.0159, 0.00719, -0.0137, 0.0163])
        @test mask3 == [false, false, false, false]
    end

    @testset "oh_eff_all_coupled: amplitude-outlier ringing does not contaminate lower modes, but comparable-magnitude high modes are kept" begin
        # Regression for the actual bug found in the 20-sample validation
        # sweep: a live solver run at contact onset showed every mode except
        # the single highest retained one (l=M) sitting at floating-point
        # noise (~1e-15), while l=M itself carried a large, spurious value
        # (order 1e-1) -- classic Gibbs-type ringing from representing the
        # contact-onset kink in a finite Legendre truncation, not real
        # physics. Unfiltered, that alone drove Oh_eff for EVERY mode
        # (including mode 2, which had not moved at all) from Oh0=57.4 down
        # to ~0.014 -- a ~4000x, entirely spurious thinning event -- which
        # collapsed predicted contact times to near-zero (median tc error
        # 62% in the full validation sweep).
        #
        # A blanket "exclude the top ~10% of modes by INDEX" fix (an earlier
        # version of this filter) fixed that case but broke a DIFFERENT one:
        # a live low-We run showed mode M itself reaching genuine, large
        # amplitude COMPARABLE to mode 2's (ratio ~1.0), not a >100x outlier
        # -- excluding it by index alone threw away real physics and cut
        # predicted CoR roughly in half. The amplitude-RATIO mask (checked
        # here) correctly keeps that case while still catching the ringing.
        M = 12
        Oh0 = 57.4
        lambda_c = 30507.0
        a = 0.7431
        eps_ST = 0.99956
        stx = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid)

        Adot_vec = fill(1e-15, M - 1)   # every low/mid mode at the noise floor
        Adot_vec[end] = 0.286           # l=M alone carries a large, spurious value
        Oh_eff = oh_eff_all_coupled(stx, Oh0, Adot_vec)
        @test isapprox(Oh_eff[1], Oh0; rtol=1e-6)   # mode 2: unaffected by the filtered-out ringing

        # A genuinely large value in a mid mode (e.g. l=8) must still thin
        # mode 2 -- confirming the mask targets only true amplitude outliers,
        # not real higher-mode coupling (the actual point of this extension).
        Adot_vec2 = zeros(M - 1)
        Adot_vec2[1] = 1e-15
        Adot_vec2[7] = 0.286   # l=8
        Oh_eff2 = oh_eff_all_coupled(stx, Oh0, Adot_vec2)
        @test Oh_eff2[1] < 0.5 * Oh0

        # The actual regression: l=M at genuine, comparable-to-mode-2
        # amplitude (live low-We data, We=0.0158, i=60) must NOT be excluded.
        Adot_vec3 = [-0.0159, 0.00719, -0.0137, 0.0135, -0.00666, 0.00287,
                     -0.00241, 0.002, -0.00129, 0.000523, 0.0163]
        Oh_eff3_with_l12 = oh_eff_all_coupled(stx, Oh0, Adot_vec3)
        Adot_vec3_no_l12 = copy(Adot_vec3)
        Adot_vec3_no_l12[end] = 0.0
        Oh_eff3_without_l12 = oh_eff_all_coupled(stx, Oh0, Adot_vec3_no_l12)
        @test !isapprox(Oh_eff3_with_l12[1], Oh_eff3_without_l12[1]; rtol=1e-3)
    end

    @testset "Real validation fluid: bounded, decaying free oscillation (live Newton solve)" begin
        # The parameters where st_extension.jl's perturbative correction is
        # already unphysical (negative effective damping) -- see
        # julia/derivations/carreau_yasuda_nonperturbative_derivation.jl.
        M = 2
        Oh0 = 57.4
        lambda_c = 30507.0
        a = 0.7431
        eps_ST = 0.99956
        Bo = 1e-6

        stx_lamb = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:lamb)
        stx_reid = STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid)

        for stx in (stx_lamb, stx_reid)
            theta_vec = make_theta_vec(M)
            precomp = precompute_integrals(NaN, M)[1]
            sigma0 = sqrt(Float64(M * (M - 1) * (M + 2)))
            dt = 2 * pi / (sigma0 * 40)
            cfg = SimConstants(M, M + 1, Oh0, Bo, theta_vec, precomp, dt)
            ob = OBParams()

            s = DropState(M)
            s.A[2] = 0.05
            s.z = 2.0
            s.dt = dt
            history = [deepcopy(s)]

            A2_hist = Float64[s.A[2]]
            for _ in 1:400
                X0 = pack_X(history[end], M)
                R!(buf, Xv) = begin
                    st = deepcopy(history[end])
                    unpack_X!(st, Xv, M)
                    build_residual_st_exact!(buf, st, history, dt, 0, cfg, ob, stx)
                end
                J_fn(Xv) = begin
                    st = deepcopy(history[end])
                    unpack_X!(st, Xv, M)
                    build_jacobian_st_exact(st, history, dt, 0, cfg, ob, stx)
                end
                X = copy(X0)
                converged = newton_solve!(X, R!, J_fn)
                @test converged
                new_state = deepcopy(history[end])
                unpack_X!(new_state, X, M)
                new_state.t += dt
                new_state.dt = dt
                push!(history, new_state)
                length(history) > 2 && popfirst!(history)
                push!(A2_hist, new_state.A[2])
            end

            @test all(isfinite, A2_hist)
            @test maximum(abs.(A2_hist)) < 1.0
            @test abs(A2_hist[end]) < abs(A2_hist[1])
        end
    end
end
