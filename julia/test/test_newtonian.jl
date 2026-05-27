using Test
using DropSolver

@testset "Newtonian smoke test: free drop, no contact" begin
    M = 4; Oh = 0.1; Fr = 1e6
    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 0.05
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
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
