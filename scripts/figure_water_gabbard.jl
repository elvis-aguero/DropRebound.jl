# Our own water measurements against the model, overlaid with the least viscous
# silicone oil in Gabbard et al. (2025) -- an independent apparatus, an independent
# fluid, at a comparable but not identical Ohnesorge.
#
# WHY THIS IS A SECOND CHECK AND NOT A REPEAT OF THE WATER CURVE. Our water drops and
# Gabbard et al.'s least-viscous oil sit at similar Oh (0.0068 against roughly 0.021 to
# 0.030) but are otherwise nothing alike: different fluid, different drop generator,
# different lab, different measurement of restitution. The model is handed the same
# equations, the same discretisation and the same contact treatment for both, with
# nothing tuned per dataset. If the model tracked our water and missed theirs, or the
# reverse, that would say the agreement was a property of one apparatus rather than of
# the physics.
#
# WHICH ROWS ARE "1 cSt". Gabbard et al.'s own Table 1 reports only the range of
# viscosities used, not a per-fluid label, so the least-viscous fluid has to be
# identified from the data itself -- and a plain threshold on Oh is not enough. Sorting
# all 935 rows by Oh shows a clean gap: 284 rows sit in [0.0214, 0.0300] with Bo varying
# smoothly and continuously from 0.016 to 0.075 (one fluid run across several nozzle
# sizes, the way the paper's own experiment varies R), then NOTHING until a separate
# 40-row cluster at Oh in [0.0139, 0.0169] with Bo in [0.19, 0.42] -- a different,
# unrelated population, not a continuation of the first. An Oh threshold alone (the
# first draft of this script used Oh < 0.025) cuts through the first group and pulls in
# the second, which is how two visually distinct clusters of colour ended up on one
# scatter with one label. The 284-row group is used here; the other is excluded.
#
#   julia --project=docs -t 6 scripts/figure_water_gabbard.jl
#
# Writes outputs/csv/figure_water_gabbard.csv and
# outputs/figures/figure_water_gabbard.png.

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
BLAS.set_num_threads(1)          # one BLAS thread per worker; see parallel_curves
include(joinpath(@__DIR__, "_curve.jl"))
include(joinpath(@__DIR__, "_runcache.jl"))

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
const DATA = joinpath(ROOT, "data")
mkpath(OUT); mkpath(FIGS)

## `--preview` renders the same figure at throwaway resolution, to judge the layout
## without waiting on a production sweep; see figure_thenarianto_model.jl.
const PREVIEW = "--preview" in ARGS
const M_RUN, K_RUN, T_MAX = PREVIEW ? (30, 2, 15.0) : (90, 3, 25.0)
const TH, WE_N0, WE_MAX_PTS = 0.03, 8, PREVIEW ? 16 : 90
const WE_X_TH, WE_X_MAX = 0.02, PREVIEW ? 0.5 : 0.10

## Our own water drop (matches figure_water_contamination.jl's clean case).
const RHO, GRAV, R_DROP, GAMMA, MU = 997.0, 9.81, 3.0e-4, 0.0728, 1.0e-3
const OH_WATER = MU / sqrt(RHO * GAMMA * R_DROP)
const BO_WATER = RHO * GRAV * R_DROP^2 / GAMMA

## The single self-consistent low-Oh cluster identified above -- see the header.
const GAB_OH_MIN, GAB_OH_MAX = 0.021, 0.030

## One colour scale for Bond, shared by both experimental datasets and both model
## curves -- see figure_thenarianto_model.jl for why: population identity is
## carried by marker shape (circle for ours, square for Gabbard) rather than by
## colour, so a shared scale can answer "what is this point's Bond number" for
## everything on the page instead of colour meaning two different things at once.
const BO_CMAP = :viridis
bo_color(bo, lo, hi) = get(cgrad(BO_CMAP), clamp((log10(bo) - lo) / (hi - lo), 0.0, 1.0))

"""
Robust linear size scale for Ohnesorge: the 5th-95th percentile of every Oh value
on the figure (both datasets pooled) maps linearly to marker size 2-10, clamped
beyond that range -- see figure_thenarianto_model.jl, which uses the same scale
(and the same reason for 2-10 rather than a literal 20-150: Plots.jl's `ms` is
roughly a point diameter, not matplotlib's `s`).
"""
function oh_size(oh_all, oh; ms_lo = 2.0, ms_hi = 10.0)
    lo, hi = quantile(oh_all, 0.05), quantile(oh_all, 0.95)
    t = lo == hi ? 0.5 : clamp((oh - lo) / (hi - lo), 0.0, 1.0)
    ms_lo + (ms_hi - ms_lo) * t
end

"""Indices that sort `oh` from largest to smallest, so a scatter plotted in that
order draws its biggest markers first and its smallest last -- see
figure_thenarianto_model.jl for why."""
biggest_first(oh) = sortperm(oh; rev = true)

"""Water restitution measured by our experiment: Weber and epsilon."""
function read_water()
    lines = readlines(joinpath(DATA, "metrics_Water.csv"))
    hdr = split(strip(lines[1]), ';')
    function idx(pat)
        i = findfirst(h -> occursin(pat, lowercase(h)), hdr)
        i === nothing && error("no column $pat in metrics_Water.csv")
        i
    end
    iw, ie = idx("we"), idx("epsilon")
    out = NTuple{2,Float64}[]
    for ln in lines[2:end]
        f = split(strip(ln), ';'); length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f, ',' => '.')); any(isnothing, v) && continue
        push!(out, (v[iw], v[ie]))
    end
    sort!(out, by = first)
end

"""Gabbard et al. (2025) rows in the single low-Oh cluster: We, cor, Oh, Bo."""
function read_gabbard()
    out = NTuple{4,Float64}[]
    for ln in eachline(joinpath(DATA, "gabbard2025_restitution.csv"))
        (startswith(ln, "#") || startswith(ln, "We")) && continue
        v = tryparse.(Float64, split(strip(ln), ','))
        any(isnothing, v) && continue
        we, cor, oh, bo = v[1], v[3], v[5], v[6]
        GAB_OH_MIN <= oh <= GAB_OH_MAX || continue
        push!(out, (we, cor, oh, bo))
    end
    sort!(out, by = first)
end

"""One model curve at fixed (Oh, Bo), over [lo, hi] in We."""
function curve(Oh, Bo, lo, hi, tag)
    ev = function (w)
        p = ImpactParams(We = w, Bo = Bo, Oh = Oh, M = M_RUN, K = K_RUN, t_max = T_MAX)
        is_measurable(p, 0.02) || return 0.0
        s = cached_series(p)
        m = score(p, s; h_thresh = 0.02)
        (m.released && isfinite(m.cor)) ? m.cor : 0.0
    end
    adaptive_curve(ev, lo, hi; th = TH, n0 = WE_N0, maxpts = WE_MAX_PTS,
                   x_th = WE_X_TH, x_max = WE_X_MAX, tag = tag)
end

function main()
    water = read_water()
    gab = read_gabbard()
    @printf("water:   %d points, Oh = %.4f, Bo = %.4f\n", length(water), OH_WATER, BO_WATER)
    @printf("gabbard: %d points, Oh in [%.4f, %.4f], Bo in [%.4f, %.4f]\n",
            length(gab), extrema(x -> x[3], gab)..., extrema(x -> x[4], gab)...)

    med_bo = median(last.(gab))
    gab_med = argmin(x -> abs(x[4] - med_bo), gab)

    we_water = first.(water)
    we_gab = first.(gab)

    ## Each curve's own roll-off, `We = 2*Bo*h`, so the left edge shows the rise
    ## rather than clipping it -- see figure_thenarianto_model.jl for why this
    ## matters when the data floor sits above the model's own cutoff.
    we_roll_water = 2 * BO_WATER * 0.02
    we_roll_gab = 2 * gab_med[4] * 0.02
    lo = min(minimum(we_water), minimum(we_gab), 0.5 * we_roll_water, 0.5 * we_roll_gab)
    hi = 1.15 * max(maximum(we_water), maximum(we_gab))
    tks = [10.0^k for k in floor(Int, log10(lo)):ceil(Int, log10(hi))]

    progress_reset()
    jobs = [(key = "water", Oh = OH_WATER, Bo = BO_WATER),
            (key = "gab",   Oh = gab_med[3], Bo = gab_med[4])]
    curves = parallel_curves(jobs) do j
        xs, ys = curve(j.Oh, j.Bo, lo, hi, j.key)
        (; j.key, xs, ys)
    end

    PREVIEW || open(joinpath(OUT, "figure_water_gabbard.csv"), "w") do io
        println(io, "series,Oh,Bo,We,cor")
        for (j, c) in zip(jobs, curves), (x, y) in zip(c.xs, c.ys)
            @printf(io, "%s,%.6g,%.6g,%.8g,%.8g\n", j.key, j.Oh, j.Bo, x, y)
        end
    end

    plt = plot(xscale = :log10, xlabel = "We", ylabel = "restitution  ε",
               xticks = tks, xlims = (10.0^(log10(lo) - 0.2), 10.0^(log10(hi) + 0.1)),
               size = (960, 660), dpi = 200,
               framestyle = :axes, grid = false, legend = :bottomright,
               legendfontsize = 9, guidefontsize = 13, tickfontsize = 12,
               ylims = (0, 1), foreground_color_axis = :gray40,
               foreground_color_border = :gray40, left_margin = 8Plots.mm,
               bottom_margin = 8Plots.mm, right_margin = 12Plots.mm, top_margin = 5Plots.mm)

    ## One shared colour scale for Bond, and one shared size scale for Ohnesorge,
    ## across both datasets and both model curves.
    bo_lo, bo_hi = extrema(log10.(vcat(last.(gab), BO_WATER)))
    oh_all = vcat(fill(OH_WATER, length(water)), (x -> x[3]).(gab))

    ## Measurements first, model on top -- the comparison is "does the line track
    ## the cloud", which reads backwards if the cloud sits over the line. Between
    ## the two datasets, Gabbard is drawn first and ours last: ours sits at the
    ## bottom of the pooled Oh range and is therefore the smallest marker on the
    ## whole figure, and the smallest marker has to be in the LAST layer drawn or
    ## it is just buried under whichever larger marker happens to land on it,
    ## no matter what order it is in within its own series. The colorbar is
    ## requested on this last call for the same reason it has to be: it is a
    ## plot-level attribute, not a per-series one, so an earlier
    ## `colorbar = true` would be silently overridden by a later `colorbar = false`.
    gab_oh = (x -> x[3]).(gab)
    ig = biggest_first(gab_oh)   # biggest markers first within this series too
    scatter!(plt, we_gab[ig], (x -> x[2]).(gab)[ig];
             marker_z = log10.(last.(gab)[ig]), color = BO_CMAP, clims = (bo_lo, bo_hi),
             marker = :rect, ms = oh_size.(Ref(oh_all), gab_oh[ig]),
             msc = :gray30, msw = 0.4, colorbar = false, label = "Gabbard 2025")
    scatter!(plt, we_water, last.(water);
             marker_z = fill(log10(BO_WATER), length(water)), color = BO_CMAP,
             clims = (bo_lo, bo_hi), marker = :circle,
             ms = oh_size.(Ref(oh_all), fill(OH_WATER, length(water))),
             msc = :gray30, msw = 0.5,
             colorbar = true, colorbar_title = "\nlog₁₀ Bo", colorbar_titlefontsize = 11,
             label = "Ours")

    for (j, c) in zip(jobs, curves)
        keep = isfinite.(c.ys)
        sum(keep) >= 2 || continue
        col = bo_color(j.Bo, bo_lo, bo_hi)
        plot!(plt, c.xs[keep], c.ys[keep]; c = col, lw = 3, ls = :solid, label = "")
    end
    plot!(plt, Float64[], Float64[]; c = bo_color(BO_WATER, bo_lo, bo_hi), lw = 3,
          label = "model, our water")
    plot!(plt, Float64[], Float64[]; c = bo_color(gab_med[4], bo_lo, bo_hi), lw = 3,
          label = "model, Gabbard")

    out = joinpath(FIGS, "figure_water_gabbard" * (PREVIEW ? "_preview" : "") * ".png")
    savefig(plt, out); println("\nwrote ", out)

    println()
    for (j, c) in zip(jobs, curves)
        @printf("%-6s Oh = %.4f  Bo = %.4f   %s\n", j.key, j.Oh, j.Bo,
                curve_report(c.xs, c.ys, TH))
    end
end

main()
