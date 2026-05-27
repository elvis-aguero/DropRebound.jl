using Test
using DropSolver: bdf_coefficients

@testset "BDF coefficients" begin
    @testset "BDF1" begin
        c = bdf_coefficients(1, 0.1, NaN)
        @test length(c) == 2
        @test c ≈ [-1.0, 1.0]
    end

    @testset "BDF2 uniform step" begin
        # rk = dt/dt_prev = 1 → ak=3/2, bk=-2, ck=1/2
        c = bdf_coefficients(2, 0.1, 0.1)
        @test length(c) == 3
        @test c ≈ [0.5, -2.0, 1.5]
    end

    @testset "BDF2 formula matches MATLAB" begin
        dt = 0.05; dt_prev = 0.1
        rk = dt / dt_prev
        ak = (1 + 2rk)/(1 + rk)
        bk = -(1 + rk)
        ck = rk^2/(1 + rk)
        c = bdf_coefficients(2, dt, dt_prev)
        @test c ≈ [ck, bk, ak]
    end

    @testset "BDF2 convergence on dy/dt = -y, y(0)=1" begin
        # Exact: y(t) = exp(-t). Verify BDF2 is 2nd order.
        errors = Float64[]
        for N in [20, 40, 80]
            dt = 1.0 / N
            y = [1.0, exp(-dt)]            # exact first two steps
            for k in 2:N-1
                c = bdf_coefficients(2, dt, dt)
                # BDF2: c[1]*y[k-1] + c[2]*y[k] + c[3]*y_{k+1} = dt*(-y_{k+1})
                # (c[3] + dt)*y_{k+1} = -(c[1]*y[k-1] + c[2]*y[k])
                rhs = -(c[1]*y[end-1] + c[2]*y[end])
                push!(y, rhs / (c[3] + dt))
            end
            push!(errors, abs(y[end] - exp(-1.0)))
        end
        # error should roughly halve each time N doubles (at least first-order)
        @test errors[2]/errors[1] < 0.6
        @test errors[3]/errors[2] < 0.6
    end
end
