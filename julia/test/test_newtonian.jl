using Test
using DropSolver

@testset "Newtonian smoke test: free drop, no contact" begin
    M = 4; Oh = 0.1; Bo = 1e-6
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 0.05
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams()

    init      = DropState(M)
    init.A[2] = 0.05
    init.z    = 10.0   # far above substrate
    init.v    = 0.0
    init.dt   = dt_max
    init.cp   = 0

    times, states = solve_drop!(cfg, ob, init;
                                t_end      = 1.0,
                                save_every = 0.2,
                                dt_init    = dt_max)

    @test length(times) >= 5
    @test times[end] >= 1.0
    @test all(isfinite, [s.A[2] for s in states])
    # Drop should fall under gravity (z decreases)
    @test states[end].z < states[1].z
end

@testset "Newtonian Lamb l=2 oscillation" begin
    M  = 6; Oh = 0.02; Bo = 1e-6   # nearly zero gravity
    l  = 2
    ω_lamb = sqrt(Float64(l*(l-1)*(l+2)))    # = sqrt(8) ≈ 2.828
    γ_lamb = Float64((l-1)*(2l+1)) * Oh       # = 5 * 0.02 = 0.1

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    # dt_max: resolve the highest mode (l=M=6) with 8 points per period
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams()

    init      = DropState(M)
    init.A[2] = 0.05     # small l=2 perturbation
    init.z    = 2.0      # far above substrate (no contact)
    init.v    = 0.0
    init.dt   = dt_max
    init.cp   = 0

    T_period  = 2π / ω_lamb
    times, states = solve_drop!(cfg, ob, init;
                                t_end      = 4 * T_period,
                                save_every = T_period / 50,
                                dt_init    = dt_max)

    A2 = [s.A[2] for s in states]
    t1, t2 = times[1], times[end]

    # Decay rate: fit ln|A2(t)| = -γ*t slope
    γ_fit = -log(abs(A2[end]) / abs(A2[1])) / (t2 - t1)
    @test abs(γ_fit - γ_lamb) / γ_lamb < 0.15   # 15% tolerance

    # Frequency: count zero-crossings
    sign_changes = findall(i -> A2[i]*A2[i+1] < 0, 1:length(A2)-1)
    if length(sign_changes) >= 4
        half_periods = diff(times[sign_changes])
        ω_fit = π / (sum(half_periods) / length(half_periods))
        @test abs(ω_fit - ω_lamb) / ω_lamb < 0.05   # 5% tolerance
    end
end
