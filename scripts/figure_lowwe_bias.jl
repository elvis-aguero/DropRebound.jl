# Where the Newtonian model agrees with experiment, and where it stops.
#
# Two independent campaigns, plotted as the signed error of the model against
# the same quantity: (eps_model - eps_exp)/eps_exp, per Weber group.
#
#   Gabbard et al. (2025)   935 impacts,  Oh in [0.014, 0.79],  Bo in [0.016, 0.071]
#   Glaco surface           174 impacts,  Oh in [0.0031, 0.024], Bo in [1e-4, 0.29]
#
# The two barely overlap: the Glaco drops are either much larger (1.5 mm, so
# Bo = 0.29) or much smaller (32-45 um, so Bo ~ 1e-4) than anything in Gabbard.
# Between them they cover a factor of 250 in Ohnesorge and 3000 in Bond, which is
# why disagreeing about the same thing is worth something.
#
# WHAT THE FIGURE SHOWS. For We above about 0.1 the model sits within roughly ten
# per cent of both campaigns, with no clear sign. Below that it runs high, and the
# excess grows as the impact gets gentler -- to +21 per cent against Gabbard and
# much more against the smallest Glaco drops.
#
# WHY THAT IS EXPECTED. The model dissipates energy only in the bulk of the drop.
# It has no work of adhesion, no contact-angle hysteresis and no dissipation in the
# air film, and those losses do not vanish with the impact energy. As We falls, the
# kinetic energy available falls with it while those losses do not, so they take a
# growing share, and a model without them predicts a bounce that is too elastic.
# In the limit We -> 0 the model returns eps -> 0.95 where a real drop settles.
#
# This is a statement about where the model applies, not a correction to it: no
# fitted loss is introduced here.
#
#   julia --project=docs scripts/figure_lowwe_bias.jl
#
# Writes outputs/figures/figure_lowwe_bias.png.

using Printf, Statistics, DelimitedFiles
using Plots
gr()

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(FIGS)

function load(file, wcol, ecol, scol)
    d, h = readdlm(joinpath(OUT, file), ','; header = true)
    h = vec(string.(h))
    g(n) = [x isa Number ? Float64(x) : NaN for x in d[:, findfirst(==(n), h)]]
    We, ex, sm = g(wcol), g(ecol), g(scol)
    keep = isfinite.(We) .& isfinite.(ex) .& isfinite.(sm) .& (ex .> 0)
    We[keep], (sm[keep] .- ex[keep]) ./ ex[keep]
end

Wg, Eg = load("gabbard_validation.csv", "We", "cor_exp", "cor_K2")
Wl, El = load("validate_glaco.csv",     "We", "cor_exp", "cor_sim")

plt = plot(xscale = :log10, xlabel = "Weber number  We",
           ylabel = "model error   (ε_model − ε_exp) / ε_exp",
           legend = :topright, size = (1000, 640), dpi = 160,
           framestyle = :axes, grid = false, guidefontsize = 15, tickfontsize = 13,
           legendfontsize = 12, ylims = (-0.6, 2.35),
           foreground_color_axis = :gray40, foreground_color_border = :gray40,
           left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,
           right_margin = 6Plots.mm, top_margin = 8Plots.mm)

## The band the model is actually good to, and the line it should sit on.
hspan!(plt, [-0.10, 0.10]; c = :gray, alpha = 0.13, label = "±10 %")
hline!(plt, [0.0]; c = :gray40, lw = 1.5, ls = :dash, label = "")

scatter!(plt, Wg, Eg; c = :steelblue, ms = 7, msw = 0, alpha = 0.8,
         label = "Gabbard et al. (2025),  Oh 0.014–0.79")
scatter!(plt, Wl, El; c = :indianred, ms = 7, msw = 0, marker = :diamond, alpha = 0.8,
         label = "Glaco surface,  Oh 0.003–0.024")

annotate!(plt, 0.09, 1.75, text("model too elastic:\nlosses it does not carry", 10, :center, :gray30))
annotate!(plt, 1.2, -0.45, text("agreement to ~10 %", 10, :center, :gray30))

out = joinpath(FIGS, "figure_lowwe_bias.png")
savefig(plt, out)
println("wrote ", out)

@printf("\n%-14s %6s %10s %10s\n", "We", "n", "Gabbard", "Glaco")
for (lo, hi) in ((0.0, 0.03), (0.03, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 20.0))
    mg = (Wg .>= lo) .& (Wg .< hi); ml = (Wl .>= lo) .& (Wl .< hi)
    @printf("[%.2f,%-6.2f) %6d %9s %10s\n", lo, hi, sum(mg) + sum(ml),
            sum(mg) == 0 ? "--" : @sprintf("%+.1f%%", 100median(Eg[mg])),
            sum(ml) == 0 ? "--" : @sprintf("%+.1f%%", 100median(El[ml])))
end
