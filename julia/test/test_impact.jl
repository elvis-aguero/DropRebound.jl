using Test
using DropSolver

# Use MATLAB canonical parameters (Oh=0.3038, Fr=53.9) at M=6 for speed.
# Quantitative KPI parity at M=20 is covered by test_matlab_parity.jl.
const _IMP_M   = 6
const _IMP_Oh  = 0.3038
const _IMP_Fr  = 53.9
const _IMP_v0  = -0.281

@testset "Drop impact: Newtonian" begin
    dt_max    = make_dt_max(_IMP_M)
    theta_vec = make_theta_vec(_IMP_M)
    precomp   = precompute_integrals(NaN, _IMP_M)[1]
    cfg       = SimConstants(_IMP_M, _IMP_M+1, _IMP_Oh, _IMP_Fr, theta_vec, precomp, dt_max)

    init = DropState(_IMP_M)
    init.z = 1.1; init.v = _IMP_v0; init.A[2] = 0.0; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, OBParams(0.0, 1.0), deepcopy(init); t_end=8.0, save_every=0.05)

    @test all(isfinite(s.z)    for s in states)
    @test all(isfinite(s.A[2]) for s in states)
    @test any(s.cp > 0 for s in states)

    contact_idx = findfirst(s -> s.cp > 0, states)
    @test states[contact_idx].z < states[1].z
end

@testset "Drop impact: OB vs Newtonian differ during contact" begin
    dt_max    = make_dt_max(_IMP_M)
    theta_vec = make_theta_vec(_IMP_M)
    precomp   = precompute_integrals(NaN, _IMP_M)[1]
    cfg       = SimConstants(_IMP_M, _IMP_M+1, _IMP_Oh, _IMP_Fr, theta_vec, precomp, dt_max)

    init = DropState(_IMP_M)
    init.z = 1.1; init.v = _IMP_v0; init.A[2] = 0.0; init.dt = dt_max; init.cp = 0

    ob_N  = OBParams(0.0, 1.0)
    ob_OB = OBParams(0.5, 0.5)

    _, states_N  = solve_drop!(cfg, ob_N,  deepcopy(init); t_end=8.0, save_every=0.05)
    _, states_OB = solve_drop!(cfg, ob_OB, deepcopy(init); t_end=8.0, save_every=0.05)

    @test any(s.cp > 0 for s in states_N)
    @test any(s.cp > 0 for s in states_OB)

    max_A2_N  = maximum(abs(s.A[2]) for s in states_N  if s.cp > 0; init=0.0)
    max_A2_OB = maximum(abs(s.A[2]) for s in states_OB if s.cp > 0; init=0.0)

    @test max_A2_N > 0 || max_A2_OB > 0
    @test abs(max_A2_N - max_A2_OB) > 0
end
