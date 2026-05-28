#!/usr/bin/env julia
# MATLAB parity: Newtonian KPI convergence vs MATLAB N=90 reference
# Uses GL-node theta_vec and CFL-based dt_max (matching MATLAB exactly)
using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using DropSolver: collect_Pl
using Printf

const MATLAB_Oh           = 0.3038
const MATLAB_Fr           = 53.9
const MATLAB_v0           = -0.281
const MATLAB_CONTACT_TIME = 2.99
const MATLAB_MAX_RADIUS   = 0.397
const MATLAB_COR          = 0.484

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

function run_case(M; t_end=8.0, save_every=0.02)
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, MATLAB_Oh, MATLAB_Fr, theta_vec, precomp, dt_max)
    ob        = OBParams(0.0, 1.0)

    init = DropState(M); init.z = 1.1; init.v = MATLAB_v0; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, ob, deepcopy(init); t_end=t_end, save_every=save_every)

    contact_idx = findall(s -> s.cp > 0, states)
    isempty(contact_idx) && return NaN, NaN, NaN

    t_contact = times[contact_idx[end]] - times[contact_idx[1]]
    max_r     = maximum(compute_contact_radius(s, cfg) for s in states if s.cp > 0)
    cor       = compute_cor(states, MATLAB_Fr)
    return t_contact, max_r, cor
end

println("MATLAB parity: Newtonian drop impact  (Oh=$(MATLAB_Oh), Fr=$(MATLAB_Fr), v0=$(MATLAB_v0))")
println("MATLAB ref (N=90): contact_time=$(MATLAB_CONTACT_TIME)  max_radius=$(MATLAB_MAX_RADIUS)  CoR=$(MATLAB_COR)")
println()
@printf("%-6s  %8s  %8s  %8s  %8s  %8s  %8s\n",
        "M", "t_c", "err%", "r_max", "err%", "CoR", "err%")
println("-"^66)

for M in [6, 10, 20, 40, 60]
    @printf("M=%-4d  ", M)
    flush(stdout)
    try
        tc, mr, cor = run_case(M)
        if isnan(tc)
            @printf("%8s  %8s  %8s  %8s  %8s  %8s\n",
                    "no_contact", "-", "-", "-", "-", "-")
        else
            e_tc  = 100 * abs(tc  - MATLAB_CONTACT_TIME) / MATLAB_CONTACT_TIME
            e_mr  = 100 * abs(mr  - MATLAB_MAX_RADIUS)   / MATLAB_MAX_RADIUS
            e_cor = 100 * abs(cor - MATLAB_COR)           / MATLAB_COR
            @printf("%8.3f  %7.1f%%  %8.4f  %7.1f%%  %8.4f  %7.1f%%\n",
                    tc, e_tc, mr, e_mr, cor, e_cor)
        end
    catch e
        @printf("  ERROR: %s\n", e)
    end
    flush(stdout)
end

println()
println("MATLAB ref:     $(MATLAB_CONTACT_TIME)            $(MATLAB_MAX_RADIUS)            $(MATLAB_COR)")
