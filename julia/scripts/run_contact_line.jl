#!/usr/bin/env julia
# Milestone 1 + 2 validation at M=20: sweep equilibrium contact angle θ_e and
# contact-line friction ξ, reporting rebound metrics.
#   θ_e = π, ξ = 0  reproduces the perfectly non-wetting GA model exactly.
# The method is best-conditioned near non-wetting (edge-load error grows as O(sinθ_d)
# away from θ_e=π; see docs/DropRebound_ContactLine.tex Remark "edge-load"), so the
# quantitatively trustworthy band is θ_e ≳ 0.9π. COR ≤ 1 is the energy-injection gate.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

const M   = 20
const Oh  = 0.1
const Bo  = 0.01
const V0  = -0.5

cfg = SimConstants(M, M+1, Oh, Bo, make_theta_vec(M),
                   precompute_integrals(NaN, M)[1], make_dt_max(M))

function run_case(theta_e, xi)
    init = DropState(M)
    init.z = 1.05; init.v = V0; init.dt = make_dt_max(M); init.cp = 0
    t, s = solve_drop!(cfg, OBParams(), init;
                       cl = CLParams(theta_e, xi), t_end = 15.0, save_every = 0.02)
    extract_kpis(t, s, cfg), maximum(x.cp for x in s)
end

println("== θ_e sweep (ξ = 0, quasi-static) ==  [Oh=$Oh Bo=$Bo We=$(V0^2) M=$M]")
println("θ_e/π   COR      t_contact   max_r_c   max_cp")
println("-"^52)
for frac in [1.00, 0.98, 0.96, 0.94, 0.92, 0.90]
    k, mcp = run_case(frac * π, 0.0)
    @printf("%.2f    %.4f   %-9.4f   %.4f    %d\n",
            frac, k.cor, k.contact_time, k.max_radius, mcp)
end

println("\n== ξ sweep (θ_e = 0.90π) ==   COR should rise with ξ (pinning resists capture)")
println("ξ        COR      t_contact   max_r_c   max_cp")
println("-"^52)
for xi in [0.0, 0.5, 1.0, 2.0, 5.0, 20.0]
    k, mcp = run_case(0.90π, xi)
    @printf("%-6.1f   %.4f   %-9.4f   %.4f    %d\n",
            xi, k.cor, k.contact_time, k.max_radius, mcp)
end
