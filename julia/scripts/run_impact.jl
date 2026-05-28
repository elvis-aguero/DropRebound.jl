#!/usr/bin/env julia
# Drop impact showcase: Newtonian vs Oldroyd-B side-by-side trajectory table.
# Parameters: Oh=0.1 (moderate viscosity), Bo=2.0 (strong gravity), M=20 modes.
# Uses GL-node theta_vec and CFL-based dt_max.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

M = 20; Oh = 0.1; Bo = 2.0

dt_max    = make_dt_max(M)
theta_vec = make_theta_vec(M)
precomp   = precompute_integrals(NaN, M)[1]
cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)

init = DropState(M)
init.z    = 1.1
init.v    = -0.5
init.A[2] = 0.0
init.dt   = dt_max
init.cp   = 0

println("Drop impact: Oh=$(Oh), Bo=$(Bo), M=$(M), z0=$(init.z), v0=$(init.v)")
println("dt_max=$(round(dt_max; sigdigits=3))")
println()

println("Running Newtonian (De1=0, beta_s=1)...")
times_N, states_N = solve_drop!(cfg, OBParams(0.0, 1.0), deepcopy(init);
                                 t_end=8.0, save_every=0.05)

println("Running Oldroyd-B (De1=0.5, beta_s=0.5)...")
times_OB, states_OB = solve_drop!(cfg, OBParams(0.5, 0.5), deepcopy(init);
                                   t_end=8.0, save_every=0.05)

println()
@printf("%-8s  %-8s  %-6s  %-8s  |  %-8s  %-6s  %-8s\n",
        "t", "z_N", "cp_N", "A2_N", "z_OB", "cp_OB", "A2_OB")
println("-"^70)

stride = 5
n_frames = min(length(times_N), length(times_OB))
for i in 1:stride:n_frames
    sN  = states_N[i]
    sOB = states_OB[i]
    @printf("%-8.3f  %-8.4f  %-6d  %-8.4f  |  %-8.4f  %-6d  %-8.4f\n",
            times_N[i],
            sN.z,  sN.cp,  sN.A[2],
            sOB.z, sOB.cp, sOB.A[2])
end
