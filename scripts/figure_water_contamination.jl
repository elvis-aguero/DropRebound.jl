# Can contamination account for the water discrepancy?
#
# Water is the one fluid of the five whose restitution the model gets wrong, and it is
# wrong in one direction: the model bounces too well, worst at low Weber number. The
# collaborators' hypothesis is that the water was not clean. This figure tests that
# hypothesis quantitatively rather than by sensitivity: it puts the measured water
# restitution on one axis and, next to it, the model run as though the liquid had been
# progressively more contaminated.
#
# THE SUBTLETY THAT DECIDES WHETHER THE FIGURE MEANS ANYTHING. An experimentalist does
# not measure a Weber number. They measure a radius and an impact speed, and they compute
#
#     We_reported = rho R v0^2 / gamma_clean
#
# using a handbook surface tension for clean water. If the liquid was in fact
# contaminated, that reported Weber number is not the Weber number the drop experienced --
# but it IS the abscissa the data sits at, and it cannot be recomputed after the fact
# because the true gamma is exactly what is unknown. So each model curve here is
# evaluated at the TRUE dimensionless groups of its liquid,
#
#     We_true = rho R v0^2 / gamma ,  Bo_true = rho g R^2 / gamma ,  Oh_true = mu / sqrt(rho gamma R)
#
# and then plotted at the abscissa the experiment would have reported for the same drop.
# Plotting each liquid at its own We instead -- which is the right thing to do when the
# question is "what does gamma do to the curve", and is what
# `figure_surface_tension_viscosity.jl` does -- would slide the curves sideways relative
# to the data and make a worse fit look like a better one.
#
# A PREDICTION THAT FALLS OUT OF THIS, AND IS ASSERTED BELOW. The drop stops rebounding
# when it can no longer climb back to the 0.02R measurement line, We_true = 2 Bo_true h.
# Both sides carry 1/gamma, so gamma cancels and the condition is v0^2 = 2 g R h: a fixed
# impact SPEED. Contamination therefore cannot move where the roll-off sits in reported
# Weber number. If the measured roll-off were somewhere else, no amount of surface
# tension would explain it, and the hypothesis under test here would be dead on arrival.
#
# WHAT THIS CANNOT REPRESENT. A surfactant does more than lower gamma and raise mu. It
# makes the interface itself dissipative -- Marangoni stresses from a non-uniform
# coverage, and an intrinsic surface viscosity -- and this model has neither. So a
# contaminated liquid that still misses the data does not clear contamination as the
# cause; it only shows that the part of contamination expressible as bulk properties is
# not enough.
#
#   julia --project=docs -t 6 scripts/figure_water_contamination.jl
#
# Writes outputs/csv/figure_water_contamination.csv and
# outputs/figures/figure_water_contamination.png.

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

const M_RUN, K_RUN, T_MAX = 90, 3, 25.0

## The experiment's own liquid and geometry. GAMMA_REPORTED is the value the
## experimentalist divided by, so it defines the abscissa and must not be varied.
const RHO             = 997.0
const GRAV            = 9.81
const R_DROP          = 3.0e-4        # m
const GAMMA_REPORTED  = 0.0728        # N/m, handbook clean water at 20 C
const MU_CLEAN        = 1.0e-3        # Pa s

## Contamination ladder. Heavier than plausible on purpose: a milder sweep already
## exists (figure_surface_tension_viscosity) and moves restitution by a few per cent,
## which is smaller than the discrepancy, so the question here is whether ANY bulk-
## property contamination reaches the data. 30 mN/m is about where a common surfactant
## saturates the interface, so the last two rows are close to the physical limit of what
## "dirty water" can mean without changing the liquid.
const CASES = [
    (key = "clean", gamma = 0.0728, mu = 1.0e-3, c = :black,      dash = :solid),
    (key = "g50",   gamma = 0.0500, mu = 1.0e-3, c = :steelblue,  dash = :solid),
    (key = "g40",   gamma = 0.0400, mu = 1.0e-3, c = :seagreen,   dash = :solid),
    (key = "g30",   gamma = 0.0300, mu = 1.0e-3, c = :goldenrod,  dash = :solid),
    (key = "g30m2", gamma = 0.0300, mu = 2.0e-3, c = :indianred,  dash = :dash),
]

label_of(cs) = cs.mu == MU_CLEAN ?
    @sprintf("γ = %.0f mN/m", 1000cs.gamma) :
    @sprintf("γ = %.0f mN/m, μ = %.0f mPa·s", 1000cs.gamma, 1000cs.mu)

"""Dimensionless groups actually experienced by a drop whose reported Weber is `we_rep`."""
function groups(cs, we_rep)
    ## The reported Weber fixes the impact speed, because R and rho are measured:
    ##   we_rep = rho R v0^2 / gamma_reported  =>  v0^2 = we_rep * gamma_reported/(rho R)
    v0sq = we_rep * GAMMA_REPORTED / (RHO * R_DROP)
    (We = RHO * R_DROP * v0sq / cs.gamma,
     Bo = RHO * GRAV * R_DROP^2 / cs.gamma,
     Oh = cs.mu / sqrt(RHO * cs.gamma * R_DROP))
end

"""Water restitution measured by the experiment: reported Weber, and epsilon."""
function read_water()
    lines = readlines(joinpath(DATA, "metrics_Water.csv"))
    hdr = split(strip(lines[1]), ';')
    ## `something(findfirst(...), error(...))` would call `error` eagerly, on every
    ## lookup including the successful ones -- Julia evaluates both arguments before
    ## `something` runs. The branch has to be explicit.
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

function main()
    exp_pts = read_water()
    we_lo_data, we_hi_data = extrema(first, exp_pts)

    ## The roll-off abscissa, and the claim that it does not depend on gamma. Asserted
    ## rather than asserted-in-prose: if a future change to the metric breaks the
    ## cancellation, this is where it should be caught.
    ##
    ## A failure here would mean the reported Weber number at which the drop can no
    ## longer climb back to the measurement line depends on the surface tension of the
    ## liquid -- which would make the roll-off position a second, independent handle on
    ## contamination, and would invalidate the reading of this figure.
    we_roll = 0.04 * (RHO * GRAV * R_DROP^2 / GAMMA_REPORTED)
    for cs in CASES
        g = groups(cs, we_roll)
        @assert isapprox(g.We, 2 * g.Bo * 0.02; rtol = 1e-12) "roll-off moved with gamma"
    end

    ## Start below the roll-off so the curve reaches zero inside the frame, and stop
    ## just past the last measurement.
    we_lo, we_hi = 0.5 * we_roll, 1.25 * we_hi_data
    @printf("clean-water Bo = %.5f   roll-off at reported We = %.3e (gamma-independent)\n",
            RHO * GRAV * R_DROP^2 / GAMMA_REPORTED, we_roll)
    @printf("experiments: %d points, reported We %.4g to %.4g\n", length(exp_pts),
            we_lo_data, we_hi_data)
    @printf("sweeping reported We %.3e to %.3g, %d cases, threads = %d\n\n",
            we_lo, we_hi, length(CASES), Threads.nthreads())

    progress_reset()
    curves = parallel_curves(CASES) do cs
        ev = function (we_rep)
            g = groups(cs, we_rep)
            p = ImpactParams(We = g.We, Bo = g.Bo, Oh = g.Oh,
                             M = M_RUN, K = K_RUN, t_max = T_MAX)
            ## Below the geometric threshold the drop cannot reach the line at all.
            ## That is a definition, not a simulation, and running it would return NaN
            ## and be read as a solver failure.
            is_measurable(p, 0.02) || return 0.0
            s = cached_series(p)
            m = score(p, s; h_thresh = 0.02)
            (m.released && isfinite(m.cor)) ? m.cor : 0.0
        end
        xs, ys = adaptive_curve(ev, we_lo, we_hi;
                                th = 0.02, x_max = 0.10, x_th = 0.02,
                                n0 = 8, maxpts = 90, tag = cs.key)
        (cs = cs, xs = xs, ys = ys)
    end

    open(joinpath(OUT, "figure_water_contamination.csv"), "w") do io
        println(io, "case,gamma_N_per_m,mu_Pa_s,we_reported,we_true,Oh,Bo,cor")
        for cu in curves, (x, y) in zip(cu.xs, cu.ys)
            g = groups(cu.cs, x)
            @printf(io, "%s,%.6g,%.6g,%.8g,%.8g,%.6g,%.6g,%.8g\n",
                    cu.cs.key, cu.cs.gamma, cu.cs.mu, x, g.We, g.Oh, g.Bo, y)
        end
    end

    tks = [10.0^k for k in floor(Int, log10(we_lo)):ceil(Int, log10(we_hi))]
    plt = plot(xscale = :log10,
               xlabel = "Weber number reported by the experiment",
               ylabel = "restitution  ε",
               size = (1000, 700), dpi = 200,
               xticks = tks, xlims = (we_lo / 1.4, we_hi * 1.15), ylims = (0, 1.0),
               framestyle = :axes, grid = false, legend = :bottomright,
               legendfontsize = 11, guidefontsize = 15, tickfontsize = 13,
               foreground_color_axis = :gray40, foreground_color_border = :gray40,
               left_margin = 9Plots.mm, bottom_margin = 9Plots.mm,
               right_margin = 6Plots.mm, top_margin = 6Plots.mm)

    scatter!(plt, [x[1] for x in exp_pts], [x[2] for x in exp_pts];
             mc = :gray55, msc = :gray35, msw = 0.7, ms = 5, alpha = 0.75,
             label = @sprintf("experiment, water (n = %d)", length(exp_pts)))
    for cu in curves
        plot!(plt, cu.xs, cu.ys; c = cu.cs.c, lw = 3, ls = cu.cs.dash,
              label = label_of(cu.cs))
    end
    vline!(plt, [we_roll]; c = :gray70, lw = 1.2, ls = :dot, label = "")
    annotate!(plt, we_roll * 1.15, 0.93,
              text("roll-off, set by v₀ alone", 9, :left, :gray45))

    out = joinpath(FIGS, "figure_water_contamination.png")
    savefig(plt, out); println("\nwrote ", out)

    ## Which case actually sits closest to the data, over the Weber range the data
    ## covers. Median relative residual against the nearest experiments in reported We.
    @printf("\n%-24s %8s %8s %10s\n", "case", "Oh_true", "Bo_true", "med |Δε|/ε")
    for cu in curves
        g = groups(cu.cs, 1.0)
        ds = Float64[]
        for (x, y) in zip(cu.xs, cu.ys)
            (we_lo_data <= x <= we_hi_data) || continue
            near = [e[2] for e in exp_pts if 1/1.3 <= e[1] / x <= 1.3]
            isempty(near) && continue
            push!(ds, abs(y - median(near)) / median(near))
        end
        isempty(ds) && continue
        @printf("%-24s %8.4f %8.4f %9.1f%%\n", label_of(cu.cs), g.Oh, g.Bo,
                100 * median(ds))
    end

    ## Continuity of each curve, in the sampler's own words: whether the ordinate
    ## tolerance was met, or the abscissa ran out of room first at a genuine step.
    println()
    for cu in curves
        @printf("%-6s %s\n", cu.cs.key, curve_report(cu.xs, cu.ys, 0.02))
    end
end

main()
