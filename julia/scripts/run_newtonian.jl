#!/usr/bin/env julia
# Run Newtonian validation: print convergence table for Oh = 0.01, 0.05, 0.1

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

println("Oh      γ_lamb   γ_fit    err%     ω_lamb   ω_fit    err%")
println("-"^65)

for Oh in [0.01, 0.05, 0.1]
    M = 6; Fr = 1e6; l = 2
    ω_lamb = sqrt(Float64(l*(l-1)*(l+2)))
    γ_lamb = Float64((l-1)*(2l+1)) * Oh

    theta_vec = collect(range(π, 0; length=M+1))
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = 2π / (sqrt(Float64(M*(M+2)*(M-1))) * 8)
    cfg       = SimConstants(M, M+1, Oh, Fr, theta_vec, precomp, dt_max)
    ob        = OBParams()

    init      = DropState(M)
    init.A[2] = 0.05; init.z = 2.0; init.dt = dt_max; init.cp = 0

    T_period  = 2π / ω_lamb
    times, states = solve_drop!(cfg, ob, init;
                                t_end      = 6 * T_period,
                                save_every = T_period / 50,
                                dt_init    = dt_max)

    A2 = [s.A[2] for s in states]
    γ_fit = -log(abs(A2[end]) / abs(A2[1])) / (times[end] - times[1])

    sign_changes = findall(i -> A2[i]*A2[i+1] < 0, 1:length(A2)-1)
    if length(sign_changes) >= 4
        hp = diff(times[sign_changes])
        ω_fit = π / (sum(hp) / length(hp))
        @printf("%.3f   %.4f   %.4f   %5.1f%%   %.4f   %.4f   %5.1f%%\n",
                Oh, γ_lamb, γ_fit, 100*abs(γ_fit-γ_lamb)/γ_lamb,
                ω_lamb, ω_fit, 100*abs(ω_fit-ω_lamb)/ω_lamb)
    else
        @printf("%.3f   %.4f   %.4f   %5.1f%%   (insufficient zero crossings)\n",
                Oh, γ_lamb, γ_fit, 100*abs(γ_fit-γ_lamb)/γ_lamb)
    end
end
