#!/usr/bin/env julia
# The three solvers, side by side: {variational, nonvariational} formulation and
# {ranked search, complementarity} contact closure.
#
# Read the table along its axes rather than cell by cell. `simulate` and `simulate_lcp` share a
# formulation and differ only in the closure, so a difference between them is the closure's.
# `simulate` and `solve_drop!` share a closure and differ in the formulation, so a difference
# between them is the formulation's. One solver alone cannot separate the two.
#
# All three are measured with the SAME contact definition -- the surface is in contact when any
# point of it lies below 0.02R, contact time runs from the first such instant to the last, and
# restitution compares centre-of-mass speed at those two instants. Anything else would compare
# the definitions instead of the solvers.
#
# Usage: julia --project=julia julia/scripts/compare_solvers.jl

using Printf, LinearAlgebra, Statistics
using DropSolver

const H_THRESH = 0.02
const NANG     = 240

"""Proximity contact metrics for a nonvariational run, matching `proximity_metrics`."""
function nv_proximity(times, states; h_thresh = H_THRESH)
    ths = range(pi/2, pi; length = NANG)
    touch = [minimum(drop_height(s, th) for th in ths) < h_thresh for s in states]
    i = findfirst(touch); j = findlast(touch)
    (i === nothing || j === nothing || i == j) && return (tc = NaN, cor = NaN, ok = false)
    (tc = times[j] - times[i], cor = abs(states[j].v / states[i].v), ok = true)
end

const CELLS = ["var/search", "var/lcp", "nonvar/search"]

"""One case in all three cells."""
function three(We, Oh, Bo; Mv = 30, Kv = 2, Mn = 20, t_end = 25.0)
    out = Dict{String,Any}()
    p = ImpactParams(We = We, Bo = Bo, Oh = Oh, M = Mv, K = Kv, t_max = t_end)
    for (name, f) in (("var/search", simulate), ("var/lcp", simulate_lcp))
        t0 = time()
        out[name] = try
            r = f(p); m = proximity_metrics(p, r; h_thresh = H_THRESH)
            (tc = m.tc, cor = m.cor, wall = time()-t0, ok = isfinite(m.cor))
        catch
            (tc = NaN, cor = NaN, wall = time()-t0, ok = false)
        end
    end
    ## Reid's exact per-mode coefficients: this formulation's closure, not Lamb's asymptotics
    cfg = SimConstants(Mn, Mn+1, Oh, Bo, make_theta_vec(Mn),
                       precompute_integrals(NaN, Mn)[1], make_dt_max(Mn); viscous = :reid)
    t0 = time()
    out["nonvar/search"] = try
        s = DropState(Mn); s.z = 1.0 + 4*H_THRESH; s.v = -sqrt(We)
        s.dt = make_dt_max(Mn); s.cp = 0
        ts, st = solve_drop!(cfg, OBParams(), s; t_end = t_end, save_every = 0.005)
        m = nv_proximity(ts, st)
        (tc = m.tc, cor = m.cor, wall = time()-t0, ok = m.ok)
    catch
        (tc = NaN, cor = NaN, wall = time()-t0, ok = false)
    end
    out
end

const CASES = [(0.2, 0.0373, 0.02), (0.5, 0.0373, 0.02), (1.0, 0.0767, 0.02),
               (1.0, 0.3038, 0.0189), (2.0, 0.3038, 0.0189)]

println("Contact and restitution at h = 0.02R, first touch to last release.\n")
@printf("%-6s %-8s | %-25s | %-25s | %s\n", "We", "Oh", "contact time", "restitution", "wall (s)")
@printf("%-6s %-8s | %-8s %-8s %-8s| %-8s %-8s %-8s| %s\n", "", "",
        "v/srch","v/lcp","nv/srch","v/srch","v/lcp","nv/srch","")
rows = []
for (We, Oh, Bo) in CASES
    o = three(We, Oh, Bo); push!(rows, (We, Oh, o))
    f(x) = isfinite(x) ? @sprintf("%.4f", x) : "  --    "
    @printf("%-6.3g %-8.4g | %-8s %-8s %-8s| %-8s %-8s %-8s| %s\n", We, Oh,
            f(o["var/search"].tc), f(o["var/lcp"].tc), f(o["nonvar/search"].tc),
            f(o["var/search"].cor), f(o["var/lcp"].cor), f(o["nonvar/search"].cor),
            join([@sprintf("%.1f", o[c].wall) for c in CELLS], "/"))
end

println("\nWhat each axis is worth (median relative difference over the cases above):")
function meddiff(a, b, field)
    d = [abs(getproperty(o[a],field) - getproperty(o[b],field)) / abs(getproperty(o[a],field))
         for (_,_,o) in rows
         if isfinite(getproperty(o[a],field)) && isfinite(getproperty(o[b],field)) &&
            getproperty(o[a],field) != 0]
    isempty(d) ? NaN : median(d)
end
for field in (:tc, :cor)
    @printf("  %-4s closure     (var/search   vs var/lcp)      : %5.1f%%\n",
            field, 100*meddiff("var/search","var/lcp",field))
    @printf("  %-4s formulation (var/search   vs nonvar/search): %5.1f%%\n",
            field, 100*meddiff("var/search","nonvar/search",field))
    @printf("  %-4s both        (var/lcp      vs nonvar/search): %5.1f%%\n",
            field, 100*meddiff("var/lcp","nonvar/search",field))
end
@printf("\ncases completed: %s\n",
        join([@sprintf("%s %d/%d", c, count(r -> r[3][c].ok, rows), length(rows)) for c in CELLS], "  "))
