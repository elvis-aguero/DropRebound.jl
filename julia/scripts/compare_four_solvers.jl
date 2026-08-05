#!/usr/bin/env julia
# The two-by-two: {variational, nonvariational} x {ranked search, complementarity}.
#
# The reason for keeping all four is that a disagreement between any two of them localises
# itself. Read the table by rows and columns rather than cell by cell: agreement down a column
# means the formulation is not what is driving a difference, agreement across a row means the
# contact closure is not. With one formulation and one closure, every disagreement is
# confounded with every other and none of them can be attributed.
#
# All four are measured with the SAME contact definition -- the surface is in contact when any
# point of it lies below 0.02R, contact time runs from the first such instant to the last, and
# restitution compares the centre-of-mass speed at those two instants. Anything else would
# compare the definitions instead of the solvers.
#
# Usage: julia --project=julia julia/scripts/compare_four_solvers.jl

using Printf, LinearAlgebra, Statistics
using DropSolver

const H_THRESH = 0.02
const NANG     = 240

"""Proximity contact metrics for a nonvariational run, matching `proximity_metrics`."""
function nv_proximity(cfg, times, states; h_thresh = H_THRESH)
    ths = range(pi/2, pi; length = NANG)
    touch = falses(length(states))
    for (k, s) in enumerate(states)
        touch[k] = minimum(drop_height(s, th) for th in ths) < h_thresh
    end
    i = findfirst(touch); j = findlast(touch)
    (i === nothing || j === nothing || i == j) &&
        return (tc = NaN, cor = NaN, ok = false)
    (tc = times[j] - times[i], cor = abs(states[j].v / states[i].v), ok = true)
end

"""One case in all four cells."""
function four(We, Oh, Bo; Mv = 30, Kv = 2, Mn = 20, t_end = 25.0)
    v0 = -sqrt(We)
    out = Dict{String,Any}()

    ## variational, both closures
    p = ImpactParams(We = We, Bo = Bo, Oh = Oh, M = Mv, K = Kv, t_max = t_end)
    for (name, f) in (("var/search", simulate), ("var/lcp", simulate_lcp))
        t0 = time()
        try
            r = f(p); m = proximity_metrics(p, r; h_thresh = H_THRESH)
            out[name] = (tc = m.tc, cor = m.cor, wall = time()-t0, ok = isfinite(m.cor))
        catch e
            out[name] = (tc = NaN, cor = NaN, wall = time()-t0, ok = false)
        end
    end

    ## nonvariational, both closures. Reid's exact per-mode coefficients, which is this
    ## formulation's closure -- not Lamb's small-Oh asymptotics.
    cfg = SimConstants(Mn, Mn+1, Oh, Bo, make_theta_vec(Mn),
                       precompute_integrals(NaN, Mn)[1], make_dt_max(Mn);
                       viscous = :reid)
    mkinit() = begin
        s = DropState(Mn); s.z = 1.0 + 4*H_THRESH; s.v = v0
        s.dt = make_dt_max(Mn); s.cp = 0; s
    end
    for (name, f) in (("nonvar/search", solve_drop!), ("nonvar/lcp", solve_drop_lcp!))
        t0 = time()
        try
            ts, st = f(cfg, OBParams(), mkinit(); t_end = t_end, save_every = 0.005)
            m = nv_proximity(cfg, ts, st)
            out[name] = (tc = m.tc, cor = m.cor, wall = time()-t0, ok = m.ok)
        catch e
            out[name] = (tc = NaN, cor = NaN, wall = time()-t0, ok = false)
        end
    end
    out
end

const CELLS = ["var/search", "var/lcp", "nonvar/search", "nonvar/lcp"]
const CASES = [(0.2, 0.0373, 0.02), (0.5, 0.0373, 0.02), (1.0, 0.0767, 0.02),
               (1.0, 0.3038, 0.0189), (2.0, 0.3038, 0.0189)]

println("Contact and restitution at h = 0.02R, first touch to last release.\n")
@printf("%-6s %-8s | %-32s | %-32s | %s\n", "We", "Oh",
        "contact time", "restitution", "wallclock (s)")
@printf("%-6s %-8s | %-7s %-7s %-8s %-7s| %-7s %-7s %-8s %-7s| %s\n", "", "",
        "v/srch","v/lcp","nv/srch","nv/lcp","v/srch","v/lcp","nv/srch","nv/lcp","")
rows = []
for (We, Oh, Bo) in CASES
    o = four(We, Oh, Bo); push!(rows, (We, Oh, o))
    f(x) = isfinite(x) ? @sprintf("%.4f", x) : "  --   "
    @printf("%-6.3g %-8.4g | %-7s %-7s %-8s %-7s| %-7s %-7s %-8s %-7s| %s\n", We, Oh,
            f(o["var/search"].tc), f(o["var/lcp"].tc),
            f(o["nonvar/search"].tc), f(o["nonvar/lcp"].tc),
            f(o["var/search"].cor), f(o["var/lcp"].cor),
            f(o["nonvar/search"].cor), f(o["nonvar/lcp"].cor),
            join([@sprintf("%.1f", o[c].wall) for c in CELLS], "/"))
end

println("\nHow much each axis matters (median relative difference over the cases above):")
function meddiff(a, b, field)
    d = Float64[]
    for (_, _, o) in rows
        x, y = getproperty(o[a], field), getproperty(o[b], field)
        (isfinite(x) && isfinite(y) && x != 0) && push!(d, abs(x-y)/abs(x))
    end
    isempty(d) ? NaN : median(d)
end
for field in (:tc, :cor)
    @printf("  %-4s  closure, within variational   : %5.1f%%\n", field,
            100*meddiff("var/search","var/lcp",field))
    @printf("  %-4s  closure, within nonvariational: %5.1f%%\n", field,
            100*meddiff("nonvar/search","nonvar/lcp",field))
    @printf("  %-4s  formulation, under search     : %5.1f%%\n", field,
            100*meddiff("var/search","nonvar/search",field))
    @printf("  %-4s  formulation, under LCP        : %5.1f%%\n", field,
            100*meddiff("var/lcp","nonvar/lcp",field))
end
@printf("\ncases each cell completed: %s\n",
        join([@sprintf("%s %d/%d", c, count(r -> r[3][c].ok, rows), length(rows)) for c in CELLS], "  "))
