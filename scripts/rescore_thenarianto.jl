# Re-draw the Thenarianto comparison from runs already on disk, without simulating.
#
# The Weber points come from whatever the last sweep wrote to
# outputs/csv/figure_thenarianto_model.csv; each is re-scored from the cached
# trajectory under the current postprocessing. Nothing is integrated. This exists to
# answer "did the metric change fix the curve" in seconds rather than in an hour, and
# to keep a half-finished sweep from being the only way to look.
#
# It is not a substitute for the sweep: the Weber points are wherever the previous
# adaptive pass happened to land, so coverage is whatever it was then. A point whose
# run is not cached is skipped and counted rather than simulated.
#
#   julia --project=docs scripts/rescore_thenarianto.jl
#
# Writes outputs/figures/figure_thenarianto_rescored_{glaco,bsi}.png.

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
include(joinpath(@__DIR__, "_runcache.jl"))

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")

## The settings the cached sweep ran at. These must match, or every lookup misses.
const M_RUN, K_RUN, T_MAX = 90, 3, 25.0
const R_SPLIT = 0.8e-3

const SURFACES = [(key = "glaco", csv = "glaco_restitution.csv", title = "Glaco coating"),
                  (key = "bsi",   csv = "bsi_restitution.csv",   title = "Black silicon")]

function read_data(csv)
    We, C, R = Float64[], Float64[], Float64[]
    for ln in eachline(joinpath(ROOT, "data", csv))
        (startswith(ln, "#") || startswith(ln, "We")) && continue
        v = parse.(Float64, split(strip(ln), ','))
        push!(We, v[1]); push!(C, v[2]); push!(R, v[5])
    end
    (; We, C, R)
end

function main()
    rows = NamedTuple[]
    for ln in eachline(joinpath(OUT, "figure_thenarianto_model.csv"))
        startswith(ln, "surface") && continue
        f = split(strip(ln), ',')
        length(f) >= 6 || continue
        push!(rows, (key = f[1], big = f[2] == "large", R = parse(Float64, f[3]),
                     Oh = parse(Float64, f[4]), Bo = parse(Float64, f[5]),
                     We = parse(Float64, f[6])))
    end
    @printf("%d stored points\n", length(rows))

    hit, miss = 0, 0
    scored = NamedTuple[]
    for r in rows
        p = ImpactParams(We = r.We, Bo = r.Bo, Oh = r.Oh, M = M_RUN, K = K_RUN, t_max = T_MAX)
        s = cache_load(cache_path(p))
        if s === nothing; miss += 1; continue; end
        hit += 1
        m = score(p, s; h_thresh = 0.02)
        push!(scored, (r..., cor = (m.released && isfinite(m.cor)) ? m.cor : 0.0))
    end
    @printf("cache: %d hits, %d misses (misses are skipped, not simulated)\n\n", hit, miss)

    for sf in SURFACES
        d = read_data(sf.csv)
        big, sml = d.R .>= R_SPLIT, d.R .< R_SPLIT
        lo, hi = minimum(d.We), maximum(d.We)
        tks = [10.0^k for k in floor(Int, log10(lo)):ceil(Int, log10(hi))]
        plt = plot(xscale = :log10, xlabel = "Weber number  We", ylabel = "restitution  ε",
                   xticks = tks, xlims = (10.0^(log10(lo)-0.25), 10.0^(log10(hi)+0.25)),
                   title = @sprintf("%s  —  re-scored from cache (M = %d, K = %d)", sf.title, M_RUN, K_RUN),
                   titlefontsize = 11, size = (960, 640), dpi = 200, framestyle = :axes,
                   grid = false, legend = :bottomright, legendfontsize = 9,
                   guidefontsize = 13, tickfontsize = 12, ylims = (0, 1),
                   foreground_color_axis = :gray40, foreground_color_border = :gray40,
                   left_margin = 8Plots.mm, bottom_margin = 8Plots.mm,
                   right_margin = 5Plots.mm, top_margin = 5Plots.mm)

        for isbig in (true, false), Rv in sort(unique([x.R for x in scored
                                                       if x.key == sf.key && x.big == isbig]))
            pts = sort([x for x in scored if x.key == sf.key && x.big == isbig && x.R == Rv],
                       by = x -> x.We)
            length(pts) >= 2 || continue
            plot!(plt, [x.We for x in pts], [x.cor for x in pts];
                  c = isbig ? :indianred : :steelblue, lw = 2.2, marker = :circle, ms = 2.5,
                  alpha = 0.9, label = "")
        end
        scatter!(plt, d.We[sml], d.C[sml]; mc = :white, msc = :steelblue, msw = 1.3, ms = 5,
                 label = @sprintf("R = %.3f–%.3f mm (n=%d)", 1000minimum(d.R[sml]),
                                  1000maximum(d.R[sml]), sum(sml)))
        scatter!(plt, d.We[big], d.C[big]; mc = :indianred, msc = :indianred, msw = 0, ms = 5,
                 alpha = 0.85,
                 label = @sprintf("R = %.2f–%.2f mm (n=%d)", 1000minimum(d.R[big]),
                                  1000maximum(d.R[big]), sum(big)))
        out = joinpath(FIGS, "figure_thenarianto_rescored_$(sf.key).png")
        savefig(plt, out); println("wrote ", out)
    end

    ## The question this script exists to answer.
    println("\nlargest jump between neighbouring stored points, per curve:")
    for sf in SURFACES, isbig in (true, false),
        Rv in sort(unique([x.R for x in scored if x.key == sf.key && x.big == isbig]))
        pts = sort([x for x in scored if x.key == sf.key && x.big == isbig && x.R == Rv],
                   by = x -> x.We)
        length(pts) >= 2 || continue
        d = maximum(abs(pts[i+1].cor - pts[i].cor) for i in 1:length(pts)-1)
        @printf("  %-6s R = %.3f mm  n=%3d   max |Δε| = %.4f\n", sf.key, 1000Rv, length(pts), d)
    end
end

main()
