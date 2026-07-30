using Test
using DropSolver: reid_char, dominant_root, second_root, reid_lambda_omega2,
                   drop_viscous_coeffs, _dominant_root_direct

@testset "Reid (arbitrary-Oh drop viscous coefficients)" begin

    @testset "dominant_root satisfies reid_char = 0" begin
        for Oh in (0.001, 0.006, 0.05, 0.1, 0.3, 1.85, 57.4)
            for l in (2, 4, 8, 16)
                q = dominant_root(Oh, l)
                @test abs(reid_char(q, Oh, l)) < 1e-8
            end
        end
    end

    @testset "second_root satisfies reid_char = 0" begin
        for Oh in (0.05, 0.5, 5.0, 57.4)
            for l in (2, 4, 8)
                q1 = dominant_root(Oh, l)
                q2 = second_root(Oh, l, q1)
                @test abs(reid_char(q2, Oh, l)) < 1e-8
            end
        end
    end

    @testset "reduces to Lamb as Oh -> 0" begin
        # Lamb is the small-viscosity asymptotics of exactly this equation, so the
        # relative discrepancy must vanish monotonically as Oh -> 0.
        prev_err = Inf
        for Oh in (1e-2, 1e-3, 1e-4)
            lam, _, resid = reid_lambda_omega2(Oh, 2)
            @test resid < 1e-8
            lamb = Oh * (2 - 1) * (2 * 2 + 1)
            err = abs(lamb - lam) / lamb
            @test err < prev_err
            prev_err = err
        end
        @test prev_err < 0.01
    end

    @testset "reduces to Molacek & Bush's high-Oh D_l limit" begin
        # D_l = omega2*lambda / (l^2*Oh*A_l), A_l = omega0^2/omega2 (see the
        # derivation script's Section 3 for why this recasting is valid).
        for l in (2, 3, 4, 5, 10)
            om0sq = Float64(l) * (l - 1) * (l + 2)
            target_D = (l - 1) * (2l^2 + 4l + 3) / (l^2 * (2l + 1))
            lam, om2, resid = reid_lambda_omega2(1000.0, l)
            @test resid < 1e-8
            A = om0sq / om2
            D = lam * A / (l^2 * 1000.0)
            @test abs(D - target_D) / target_D < 1e-3
        end
    end

    @testset "wrong-branch regression: l=16, Oh=0.3" begin
        # A direct (non-continued) Lamb-seeded Newton solve converges to a real,
        # more strongly damped, non-dominant higher overtone here (reproduced
        # from SpectralKM.jl's own documented failure) -- confirmed by an
        # independent real-axis scan of reid_char to be a genuine root, just not
        # the dominant one. dominant_root's continuation finds the genuinely
        # less-damped (per Reid's own smallest-Re(sigma) ordering) dominant root
        # instead -- which is complex here, since Reid's exact equation is still
        # underdamped at a point where Lamb's own cruder coefficient (148.5) has
        # already exceeded omega_{l,0} (65.7) and would spuriously predict
        # overdamping. So the right check is NOT "close to Lamb's c" (known to be
        # a poor approximation at this extreme (l,Oh), per SpectralKM's own
        # measured error table) but "continuation finds the less-damped root".
        l, Oh = 16, 0.3
        q_direct = _dominant_root_direct(Oh, l)
        q_tracked = dominant_root(Oh, l)
        @test abs(reid_char(q_direct, Oh, l)) < 1e-8
        @test abs(reid_char(q_tracked, Oh, l)) < 1e-8
        sigma_direct = q_direct^2 * Oh
        sigma_tracked = q_tracked^2 * Oh
        @test abs(sigma_tracked - sigma_direct) / abs(sigma_direct) > 0.5
        @test real(sigma_tracked) < real(sigma_direct)

        lam, _, resid = reid_lambda_omega2(Oh, l)
        @test resid < 1e-8
        @test lam ≈ real(sigma_tracked)   # reid_lambda_omega2 uses the tracked (dominant) root
    end

    @testset "omega2 stays near the inviscid omega_{l,0}^2 (viscosity damps, doesn't restore)" begin
        for l in (2, 4, 8), Oh in (0.8, 1.0, 3.0, 10.0)
            om0sq = Float64(l) * (l - 1) * (l + 2)
            _, om2, resid = reid_lambda_omega2(Oh, l)
            @test resid < 1e-8
            @test 0.2 * om0sq < om2 < 5 * om0sq
        end
    end

    @testset "drop_viscous_coeffs(:lamb) is bit-for-bit the closed-form Lamb formula" begin
        M = 12
        lambda, omega2 = drop_viscous_coeffs(M, 0.03, :lamb)
        @test length(lambda) == M - 1 == length(omega2)
        for (i, n) in enumerate(2:M)
            @test lambda[i] == 0.03 * (n - 1) * (2n + 1)
            @test omega2[i] == Float64(n) * (n - 1) * (n + 2)
        end
    end

    @testset "drop_viscous_coeffs(:reid) matches reid_lambda_omega2 per mode" begin
        M = 6
        Oh = 0.05
        lambda, omega2 = drop_viscous_coeffs(M, Oh, :reid)
        for (i, n) in enumerate(2:M)
            lam, om2, _ = reid_lambda_omega2(Oh, n)
            @test lambda[i] ≈ lam
            @test omega2[i] ≈ om2
        end
    end

    @testset "invalid model symbol throws" begin
        @test_throws ArgumentError drop_viscous_coeffs(6, 0.05, :bogus)
    end
end
