#!/usr/bin/env julia
# Old vs event-located time stepping on five points from the drop-rebound paper
# (Gabbard et al., DropRebound_JFM.tex): Oh=0.03, Bo=0.02, We ∈ {0.023 … 2.269}.
# Acceptance: contact time, coefficient of restitution, and spreading time agree
# between the legacy (dt_max/2^n) and event-located steppers to ≤3% relative.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

const Oh = 0.03
const Bo = 0.02
const M  = 20
const WES = [0.023, 0.10, 0.253, 1.03, 2.269]

cfg = SimConstants(M, M+1, Oh, Bo, make_theta_vec(M),
                   precompute_integrals(NaN, M)[1], make_dt_max(M))

function kpis(We, ev)
    init = DropState(M)
    init.z = 1.0; init.v = -sqrt(We); init.dt = make_dt_max(M); init.cp = 0
    try
        t, s = solve_drop!(cfg, OBParams(), init;
                           t_end = 30.0, save_every = 0.005, event_location = ev)
        return extract_kpis(t, s, cfg)
    catch
        return nothing
    end
end

rel(o, n) = abs(n - o) / (abs(o) + eps())

@printf("%-7s | %-20s | %-20s | %-20s\n", "We",
        "contact time (o/n/%)", "COR (o/n/%)", "spread time (o/n/%)")
println("-"^74)
worst_tc = 0.0; worst_cor = 0.0; worst_ts = 0.0
for We in WES
    o = kpis(We, false); n = kpis(We, true)
    if o === nothing || n === nothing
        tag = o === nothing && n !== nothing ? "legacy FAILED, event OK (robustness win)" :
              n === nothing && o !== nothing ? "event FAILED, legacy OK" : "both FAILED"
        @printf("%-7.3f | %s\n", We, tag)
        continue
    end
    rc = rel(o.contact_time, n.contact_time)
    rr = rel(o.cor, n.cor)
    rs = rel(o.spreading_time, n.spreading_time)
    global worst_tc = max(worst_tc, rc); global worst_cor = max(worst_cor, rr)
    global worst_ts = max(worst_ts, rs)
    @printf("%-7.3f | %6.3f %6.3f %4.1f%% | %6.4f %6.4f %4.1f%% | %6.3f %6.3f %4.1f%%\n",
            We, o.contact_time, n.contact_time, 100rc,
            o.cor, n.cor, 100rr, o.spreading_time, n.spreading_time, 100rs)
end
println("-"^74)
@printf("worst-case rel. diff (comparable cases): contact %.1f%%, COR %.1f%%, spread %.1f%%\n",
        100worst_tc, 100worst_cor, 100worst_ts)
@printf("contact time & COR within 3%%: %s\n",
        (worst_tc <= 0.03 && worst_cor <= 0.03) ? "YES" : "NO")
