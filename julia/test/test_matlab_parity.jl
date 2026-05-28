using Test
using DropSolver
using DropSolver: collect_Pl

# ── MATLAB reference values (N=90 modes, v3 linearized solver, BDF1) ──────────
# Physical: ρ=0.96 g/cm³, σ=20.5 dyn/cm, R=0.0201 cm, ν=0.199 cm²/s, v₀=-9.156 cm/s
# Oh = ν√(ρ/(σR)) = 0.3038    Fr = σ/(ρgR²) = 53.9    v₀_dimless = -0.281
# time_unit = √(ρR³/σ) = 6.16×10⁻⁴ s
const MATLAB_Oh           = 0.3038
const MATLAB_Fr           = 53.9
const MATLAB_v0           = -0.281
const MATLAB_CONTACT_TIME = 2.99     # dimensionless contact duration
const MATLAB_MAX_RADIUS   = 0.397    # dimensionless maximum spreading radius
const MATLAB_COR          = 0.484    # coefficient of restitution

function compute_contact_radius(state::DropState, cfg::SimConstants)
    cp = state.cp
    cp == 0 && return 0.0
    θ = cfg.theta_vec[cp]
    M = cfg.M
    P = collect_Pl(M, [cos(θ)])
    deform = sum(P[1, n+1] * state.A[n] for n in 2:M)
    return sin(θ) * (1.0 + deform)
end

function compute_cor(states, Fr)
    contact_idx = findall(s -> s.cp > 0, states)
    isempty(contact_idx) && return NaN
    first_c = contact_idx[1]
    last_c  = contact_idx[end]
    last_c == length(states) && return NaN
    v_in  = first_c > 1 ? states[first_c - 1].v : states[first_c].v
    z_in  = first_c > 1 ? states[first_c - 1].z : states[first_c].z
    v_out = states[last_c + 1].v
    z_out = states[last_c + 1].z
    E_in  = 0.5 * v_in^2
    E_out = 0.5 * v_out^2 + (z_out - z_in) / Fr
    E_in ≈ 0 && return NaN
    return sqrt(abs(E_out / E_in))
end

@testset "MATLAB parity: Newtonian KPIs at M=20 (<5% error)" begin
    M         = 20
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, MATLAB_Oh, MATLAB_Fr, theta_vec, precomp, dt_max)
    ob        = OBParams(0.0, 1.0)

    init = DropState(M)
    init.z = 1.1; init.v = MATLAB_v0; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=8.0, save_every=0.02)

    @test all(isfinite(s.z)    for s in states)
    @test all(isfinite(s.A[2]) for s in states)
    @test any(s.cp > 0 for s in states)

    contact_idx = findall(s -> s.cp > 0, states)
    t_contact = times[contact_idx[end]] - times[contact_idx[1]]
    max_r     = maximum(compute_contact_radius(s, cfg) for s in states if s.cp > 0)
    cor       = compute_cor(states, MATLAB_Fr)

    @test !isnan(cor)
    @test abs(t_contact - MATLAB_CONTACT_TIME) / MATLAB_CONTACT_TIME < 0.05
    @test abs(max_r     - MATLAB_MAX_RADIUS)   / MATLAB_MAX_RADIUS   < 0.05
    @test abs(cor       - MATLAB_COR)           / MATLAB_COR          < 0.05
end

@testset "OB impact: polymer stress changes contact metrics" begin
    M         = 20
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, MATLAB_Oh, MATLAB_Fr, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.1; init.v = MATLAB_v0; init.dt = dt_max; init.cp = 0

    ob_N  = OBParams(0.0, 1.0)
    ob_OB = OBParams(0.5, 0.5)

    _, states_N  = solve_drop!(cfg, ob_N,  deepcopy(init); t_end=8.0, save_every=0.02)
    _, states_OB = solve_drop!(cfg, ob_OB, deepcopy(init); t_end=8.0, save_every=0.02)

    @test all(isfinite(s.z)    for s in states_N)
    @test all(isfinite(s.z)    for s in states_OB)
    @test any(s.cp > 0 for s in states_N)
    @test any(s.cp > 0 for s in states_OB)

    cor_N  = compute_cor(states_N,  MATLAB_Fr)
    cor_OB = compute_cor(states_OB, MATLAB_Fr)

    @test !isnan(cor_N)
    @test !isnan(cor_OB)
    @test cor_OB != cor_N   # viscoelasticity changes rebound energy
end
