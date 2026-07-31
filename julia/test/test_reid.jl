using Test
using DropSolver: reid_char, dominant_root, second_root, reid_lambda_omega2,
                   drop_viscous_coeffs, _dominant_root_direct,
                   ReidTable, build_reid_table, reid_lambda_omega2_fast, build_reid_cache

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

    @testset "dominant_root's adaptive step count avoids the wrong-branch jump" begin
        # Regression: a FIXED step count (this function's earlier design) gave a
        # discontinuous jump to a more strongly damped, non-dominant root for
        # l=10 somewhere past Oh~900 (e.g. lambda 61483 -> 273551 between
        # Oh=958 and Oh=1014, both satisfying reid_char to near machine
        # precision -- two genuinely different roots, not non-convergence).
        # Adaptive step count (constant per-step Oh ratio) must track the
        # dominant branch continuously through this range instead.
        l = 10
        prev_lam = nothing
        for Oh in exp.(range(log(500.0), log(2000.0); length=30))
            lam, om2, resid = reid_lambda_omega2(Oh, l)
            @test resid < 1e-8
            if prev_lam !== nothing
                @test 0.5 < lam / prev_lam < 2.0   # no discontinuous jump between adjacent points
            end
            prev_lam = lam
        end
    end

    @testset "ReidTable / build_reid_table: interpolation matches the exact solve" begin
        for l in (2, 5, 10, 16, 20)
            table = build_reid_table(l; Oh_min=1e-4, Oh_max=1e4, n=150)
            @test table.l == l
            @test length(table.log_Oh) == length(table.lambda) == length(table.omega2) == 150
            max_err_lam = 0.0
            max_err_om2 = 0.0
            for Oh in exp.(range(log(2e-4), log(5e3); length=60))
                lam_exact, om2_exact, resid = reid_lambda_omega2(Oh, l)
                resid > 1e-6 && continue
                lam_fast, om2_fast = reid_lambda_omega2_fast(table, Oh)
                max_err_lam = max(max_err_lam, abs(lam_fast - lam_exact) / max(abs(lam_exact), 1e-10))
                max_err_om2 = max(max_err_om2, abs(om2_fast - om2_exact) / max(abs(om2_exact), 1e-10))
            end
            @test max_err_lam < 0.01
            @test max_err_om2 < 0.01
        end
    end

    @testset "reid_lambda_omega2_fast clamps outside the table range" begin
        table = build_reid_table(2; Oh_min=1e-2, Oh_max=1e2, n=50)
        lam_lo, om2_lo = reid_lambda_omega2_fast(table, 1e-6)
        @test lam_lo == table.lambda[1]
        @test om2_lo == table.omega2[1]
        lam_hi, om2_hi = reid_lambda_omega2_fast(table, 1e6)
        @test lam_hi == table.lambda[end]
        @test om2_hi == table.omega2[end]
    end

    @testset "build_reid_cache builds one table per mode l=2..M" begin
        M = 8
        cache = build_reid_cache(M; Oh_min=1e-3, Oh_max=1e3, n=40)
        @test length(cache) == M - 1
        for (i, l) in enumerate(2:M)
            @test cache[i].l == l
        end
    end

    @testset "fast lookup is much cheaper than the exact solve" begin
        table = build_reid_table(5; Oh_min=1e-3, Oh_max=1e3, n=60)
        reid_lambda_omega2_fast(table, 0.4)  # warm up
        reid_lambda_omega2(0.4, 5)
        t0 = time_ns()
        for _ in 1:1000
            reid_lambda_omega2_fast(table, 0.4)
        end
        t_fast = time_ns() - t0
        t0 = time_ns()
        for _ in 1:10
            reid_lambda_omega2(0.4, 5)
        end
        t_exact = (time_ns() - t0) * 100   # normalize to 1000 calls
        @test t_fast < t_exact / 100   # at least 100x faster
    end
end
