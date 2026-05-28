#!/usr/bin/env julia
# OB validation: compare (ω, γ) for Newtonian vs OB at several (De₁, β_s)

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver, Printf

function run_case(Oh, De1, beta_s; M=6, Fr=1e6, t_end=30.0)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = make_dt_max(M)
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
    ob        = OBParams(De1, beta_s)

    init      = DropState(M)
    init.A[2] = 0.05; init.z = 2.0; init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, ob, deepcopy(init);
                                t_end=t_end, save_every=0.1, dt_init=dt_max)
    A2 = [s.A[2] for s in states]
    γ  = -log(abs(A2[end]) / abs(A2[1])) / (times[end] - times[1])

    sign_changes = findall(i -> A2[i]*A2[i+1] < 0, 1:length(A2)-1)
    if length(sign_changes) >= 4
        hp = diff(times[sign_changes])
        ω  = π / (sum(hp) / length(hp))
    else
        ω = NaN
    end
    return γ, ω
end

Oh = 0.02
println("Oh=$Oh. Comparing Newtonian vs Oldroyd-B l=2 oscillations")
println()
println("De₁    β_s    γ         ω")
println("-"^40)

for (De1, bs) in [(0.0,1.0), (0.2,0.7), (0.5,0.5), (1.0,0.3)]
    γ, ω = run_case(Oh, De1, bs)
    @printf("%.1f    %.1f    %.4f    %.4f\n", De1, bs, γ, isnan(ω) ? 0.0 : ω)
end
