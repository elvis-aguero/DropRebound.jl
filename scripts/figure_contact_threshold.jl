# What the 0.02R contact line does to the restitution it defines.
#
# Contact here is a measurement convention, not a solver setting: a run is post-
# processed by `proximity_metrics`, which calls the drop "in contact" from the first
# instant a surface point comes within `h_thresh` of the substrate to the last, and
# reports restitution as the ratio of centre-of-mass speeds at those two instants.
# The published number uses h_thresh = 0.02R, chosen to match what a camera can
# resolve.
#
# WHY THAT CHOICE MIGHT NOT BE INNOCENT. The drop crosses the line before it touches
# and re-crosses it after it leaves, and it is accelerating under gravity the whole
# time. The speed picked up over that 0.02R is a fixed quantity, sqrt(2 Bo * 0.02)
# in these units, while the impact speed itself is sqrt(We). So the correction is
# negligible at large Weber number and grows without bound as We falls -- which is
# exactly the regime where the model refuses to roll off, and exactly the regime
# being compared against experiment.
#
# The test is free. `proximity_metrics` is post-processing, so one sweep can be
# scored at every threshold rather than one sweep per threshold.
#
# A fixed logarithmic grid is used rather than the adaptive sampler of `_curve.jl`,
# because every point is then independent and the whole sweep runs in parallel; the
# adaptive one has to see each result before choosing the next.
#
#   julia --project=docs -t 6 scripts/figure_contact_threshold.jl
#
# Writes outputs/csv/figure_contact_threshold.csv and
# outputs/figures/figure_contact_threshold.png.

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
BLAS.set_num_threads(1)
include(joinpath(@__DIR__, "_curve.jl"))

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

const M_RUN, K_RUN, T_MAX = 90, 3, 25.0

## The median small drop of the Glaco set: R = 0.049 mm. Median Bond number of the
## population that shows no roll-off, which is what is under suspicion.
const R_DROP = 4.9e-5
const RHO, SIGMA, MU, GRAV = 997.0, 0.072, 1.0e-3, 9.81
const OH = MU / sqrt(RHO * SIGMA * R_DROP)
const BO = RHO * GRAV * R_DROP^2 / SIGMA

const THRESH = [0.02, 0.01, 0.005, 0.002, 0.0]
const WE_LO, WE_HI, WE_PER_DEC = 5.6e-5, 5.0, 10

function main()
    n = round(Int, (log10(WE_HI) - log10(WE_LO)) * WE_PER_DEC) + 1
    Wes = exp10.(range(log10(WE_LO), log10(WE_HI); length = n))
    @printf("R = %.3f mm  Oh = %.4f  Bo = %.3e   %d Weber points, threads = %d\n",
            1000R_DROP, OH, BO, n, Threads.nthreads())
    @printf("speed gained falling through 0.02R under gravity: sqrt(2*Bo*0.02) = %.3e\n",
            sqrt(2 * BO * 0.02))
    @printf("  ... equal to the impact speed sqrt(We) at We = %.2e\n\n", 2 * BO * 0.02)

    progress_reset()
    rows = parallel_curves(Wes) do We
        p = ImpactParams(We = We, Bo = BO, Oh = OH, M = M_RUN, K = K_RUN, t_max = T_MAX)
        t0 = time()
        r = simulate(p)
        cors = [proximity_metrics(p, r; h_thresh = h).cor for h in THRESH]
        tcs  = [proximity_metrics(p, r; h_thresh = h).tc  for h in THRESH]
        _progress("thresh", 0, We, cors[1], time() - t0)
        (We = We, cors = cors, tcs = tcs)
    end

    open(joinpath(OUT, "figure_contact_threshold.csv"), "w") do io
        println(io, "We," * join(["cor_h$(h)" for h in THRESH], ",") * "," *
                    join(["tc_h$(h)" for h in THRESH], ","))
        for r in rows
            @printf(io, "%.8g,%s,%s\n", r.We,
                    join((@sprintf("%.8g", c) for c in r.cors), ","),
                    join((@sprintf("%.8g", t) for t in r.tcs), ","))
        end
    end

    cols = [:black, :steelblue, :seagreen, :goldenrod, :indianred]
    plt = plot(xscale = :log10, xlabel = "Weber number  We", ylabel = "restitution  ε",
               title = @sprintf("Effect of the contact threshold on ε   (R = %.3f mm, Oh = %.4f, Bo = %.1e)",
                                1000R_DROP, OH, BO),
               titlefontsize = 11, size = (960, 640), dpi = 200,
               framestyle = :axes, grid = false, legend = :bottomright,
               legendfontsize = 10, guidefontsize = 13, tickfontsize = 12,
               ylims = (0, 1.05), foreground_color_axis = :gray40,
               foreground_color_border = :gray40, left_margin = 8Plots.mm,
               bottom_margin = 8Plots.mm, right_margin = 5Plots.mm, top_margin = 5Plots.mm)
    for (j, h) in enumerate(THRESH)
        y = [r.cors[j] for r in rows]
        keep = isfinite.(y)
        plot!(plt, [r.We for r in rows][keep], y[keep]; c = cols[j], lw = 2.5,
              label = h == 0.0 ? "h = 0 (touching)" : @sprintf("h = %.3fR", h))
    end
    out = joinpath(FIGS, "figure_contact_threshold.png")
    savefig(plt, out)
    println("\nwrote ", out)

    @printf("\n%-10s %s\n", "We", join([@sprintf("%10s", h == 0.0 ? "h=0" : "h=$h") for h in THRESH], ""))
    for r in rows[1:max(1, length(rows) ÷ 12):end]
        @printf("%-10.2e %s\n", r.We,
                join([@sprintf("%10s", isfinite(c) ? @sprintf("%.4f", c) : "--") for c in r.cors], ""))
    end
end

main()
