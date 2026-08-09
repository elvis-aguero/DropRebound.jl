# Vorticity through a bounce: where rotation is created, and when.
#
# The drop arrives irrotational. Nothing in a free capillary oscillation makes
# vorticity if it starts without any, so every bit of rotation on this figure was
# generated during the impact and then diffused. Watching it appear is watching the
# contact do work on the interior, which is the part of the bounce the restitution
# number summarises away.
#
# WHAT IS PLOTTED
#   (a) a space-time map of |omega| just under the free surface, so the vertical axis
#       is position on the drop and the horizontal axis is time
#   (b) enstrophy, the volume integral of omega^2, which is the single number saying
#       how much rotation the drop is carrying at that instant
#
# The drop is against the substrate for the whole record: it starts in contact and the
# march stops when it releases, so the time axis IS the contact interval. Release is
# marked; there is no free flight to compare against on this plot.

using Printf, LinearAlgebra
using Plots
using DropSolver

gr()

const ASSET = joinpath(@__DIR__, "..", "..", "docs", "src", "assets")
const OUT   = joinpath(@__DIR__, "..", "..", "results")
mkpath(ASSET); mkpath(OUT)

## The vorticity formula, its derivation, the K = 1 self-check and the palettes all
## live in the shared helper, so this figure and the animation cannot drift apart.
include(joinpath(@__DIR__, "_vorticity.jl"))

# ---------------------------------------------------------------------------
# Run one impact and evaluate the field on it.
# ---------------------------------------------------------------------------
const WE, BO, OH = 2.0, 0.019, 0.0373
const M_RUN, K_RUN = DropSolver.DEFAULT_M, DropSolver.DEFAULT_K

p = ImpactParams(We = WE, Bo = BO, Oh = OH, M = M_RUN, K = K_RUN, t_max = 25.0)

assert_irrotational(We = WE, Bo = BO, Oh = OH)

println("running We = $WE, Oh = $OH at M = $M_RUN, K = $K_RUN ..."); flush(stdout)
r = DropSolver.simulate_lcp(p)
@printf("  %d frames, restitution %.4f, contact time %.3f\n", length(r.t), r.cor, r.tc)

## The contact interval, first touch to last release, as the solver defines it.
touch   = findfirst(>(0), r.cp)
release = findlast(>(0), r.cp)
t_in, t_out = r.t[touch], r.t[release]

## Subsample frames for the map: a few hundred columns is past what the image can show.
const NCOL = 400
cols = unique(round.(Int, range(1, length(r.t); length = min(NCOL, length(r.t)))))
ts   = r.t[cols]

## Evaluate just under the surface. Exactly at x = 1 the trial functions are all equal
## to one and the map is dominated by that normalisation rather than by the flow.
## The vertical axis is measured FROM THE CONTACT POLE, so 0 is the point that touches
## first and 180 is the crown. That puts the interesting end of the drop at the bottom
## of the map and makes the axis ascending, which is what the heatmap wants.
##
## The two poles themselves are left out. Vorticity vanishes on the axis by symmetry
## (it is proportional to dP/dtheta, which is zero there), so nothing is lost, and the
## angular routine divides by sin(theta) and returns NaN if handed mu = +-1 exactly.
const X_EVAL = 0.97
const NTH    = 220
angs = range(1.0, 179.0; length = NTH)      # degrees from the contact pole
mus  = -cos.(deg2rad.(angs))                # mu = -1 at the contact pole

W = [vorticity(p, r.adot[c], X_EVAL, mu) for mu in mus, c in cols]
@assert all(isfinite, W) "the surface field has non-finite entries"
ens = [enstrophy(p, r.adot[c]) for c in cols]

@printf("  peak |omega| = %.3g, peak enstrophy = %.3g\n", maximum(abs, W), maximum(ens))
@printf("  enstrophy at first touch %.3g, at release %.3g\n",
        ens[findmin(abs.(ts .- t_in))[2]], ens[findmin(abs.(ts .- t_out))[2]])

# ---------------------------------------------------------------------------
# The figure. Minimal: two panels, one shared time axis, no gridlines competing with
# the map, and type large enough to survive being shrunk into a column.
# ---------------------------------------------------------------------------
default(fontfamily = "sans-serif", framestyle = :box,
        guidefontsize = 15, tickfontsize = 13, titlefontsize = 16, legendfontsize = 13)

pa = heatmap(ts, angs, abs.(W), c = ABYSS, colorbar_title = "  |ω|",
             xlabel = "", ylabel = "angle from contact pole  (deg)",
             title = @sprintf("Vorticity through a bounce   (We = %.1f, Oh = %.3f)", WE, OH),
             yticks = ([0, 45, 90, 135, 180], ["0", "45", "90", "135", "180"]),
             grid = false, ylims = (0, 180))
## The drop is touching for the whole record: it starts against the substrate and the
## march stops as soon as it releases. So shading the contact interval would shade the
## entire axis and say nothing. Only the release is marked.
vline!(pa, [t_out], c = :white, lw = 2, ls = :dash, label = "")
annotate!(pa, t_out, 168, text("release  ", 13, :white, :right))

pb = plot(ts, ens, lw = 3, c = DEEP, label = "",
          xlabel = "time  (capillary times)", ylabel = "enstrophy  ∫ω² dV",
          grid = true, gridalpha = 0.15)
vline!(pb, [t_out], c = DEEP, lw = 2, ls = :dash, label = "")

plt = plot(pa, pb; layout = grid(2, 1, heights = [0.62, 0.38]),
           size = (1100, 780), dpi = 150,
           left_margin = 9Plots.mm, bottom_margin = 7Plots.mm, right_margin = 5Plots.mm)

savefig(plt, joinpath(ASSET, "vorticity_bounce.png"))
println("wrote ", joinpath(ASSET, "vorticity_bounce.png"))

## The underlying series, so the figure is reproducible without rerunning the impact.
open(joinpath(OUT, "vorticity_bounce.csv"), "w") do io
    println(io, "t,enstrophy,max_abs_omega_surface")
    for (j, t) in enumerate(ts)
        println(io, join((t, ens[j], maximum(abs, view(W, :, j))), ","))
    end
end
println("wrote ", joinpath(OUT, "vorticity_bounce.csv"))
