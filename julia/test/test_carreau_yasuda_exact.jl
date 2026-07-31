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

    @testset "STExactParams(:reid) construction and lookup accuracy" begin
        M = 6
        Oh0 = 57.4
        stx = STExactParams(M, Oh0, 30507.0, 0.7431, 0.99956; viscous=:reid)
        @test stx.viscous === :reid
        @test length(stx.K) == M - 1
        @test stx.reid_cache !== nothing
        @test length(stx.reid_cache) == M - 1

        theta_vec = collect(range(pi, 0; length=M+1))
        precomp = precompute_integrals(NaN, M)[1]
        cfg = SimConstants(M, M + 1, Oh0, 1e-6, theta_vec, precomp, 0.001)

        # Spot-check a few (mode, Adot) pairs against the exact (uncached) solve.
        for k in 1:M-1, Adot_char in (0.0, 1e-3, 0.5)
            lam_fast, om2_fast = oh_eff_lambda_omega2(stx, cfg, k, Adot_char)
            l = k + 1
            gammadot_char = stx.K[k] * abs(Adot_char)
            Oh_eff = Oh0 * (1 + (stx.lambda_c * gammadot_char)^stx.a)^(-stx.eps_ST)
            lam_exact, om2_exact, resid = reid_lambda_omega2(Oh_eff, l)
            if resid < 1e-6
                @test isapprox(lam_fast, lam_exact; rtol=0.02)
                @test isapprox(om2_fast, om2_exact; rtol=0.02)
            end
        end
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
