# The restitution curve, and what surface tension and viscosity each do to it.
#
# Two panels, each varying one liquid property and holding the other fixed, so
# the two effects can be read off separately:
#
#   left    gamma = 72.8 down to 28 mN/m at mu = 1.0 mPa s
#   right   mu    = 1.0 up to 8.0 mPa s   at gamma = 72.8 mN/m
#
# The ranges are deliberately wider than plausible contamination. A milder sweep
# was run first and moved restitution by under two per cent, which says nothing
# about whether the property matters -- only that the range was too narrow to see.
#
# WHAT IS HELD FIXED, AND WHY IT MATTERS. The drop radius is fixed at 0.3 mm,
# because that is what an experiment sets. The horizontal axis is the Weber number
# of the liquid actually being simulated, so each curve is that liquid's own
# restitution curve rather than a rescaling of one curve.
#
# The two properties do not enter symmetrically, and the figure is arranged to
# show that. Viscosity appears in one group only,
#
#     Oh = mu / sqrt(rho gamma R) ,
#
# so raising mu at fixed We slides the curve down and does nothing else. Surface
# tension appears in all three,
#
#     We = rho U^2 R / gamma ,   Oh = mu / sqrt(rho gamma R) ,   Bo = rho g R^2 / gamma ,
#
# so a curve at lower gamma is a curve at higher Oh AND higher Bo, and it also
# reaches a given We at a lower impact speed. Plotted against We the surface
# tension family therefore separates only through Oh and Bo, which is the honest
# way to show it: gamma has no effect on restitution that is not carried by a
# dimensionless group.
#
#   julia --project=docs scripts/figure_surface_tension_viscosity.jl
#
# Writes outputs/csv/figure_surface_tension_viscosity.csv and
# outputs/figures/figure_surface_tension_viscosity.png.

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
BLAS.set_num_threads(1)   # one BLAS thread per worker thread; see parallel_curves
include(joinpath(@__DIR__, "_curve.jl"))

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

const M_RUN, K_RUN, T_MAX = 90, 3, 25.0

const RHO    = 997.0
const GRAV   = 9.81
const R_DROP = 3.0e-4        # m
const GAMMA0 = 0.0728        # N/m, clean water
const MU0    = 1.0e-3        # Pa s, clean water

## The Weber range of the experiments, not a round number: the interesting part of
## the curve is the roll-off as We -> 0, and a grid that starts at 0.05 cannot show
## whether the model reproduces it.
const WE_LO, WE_HI = 0.004, 3.0
## One threshold sets both guarantees of `_curve.jl`: no consecutive pair of
## samples differs by more than TH in restitution, and the sweep is followed to
## smaller We until restitution itself drops below TH.
const TH         = 0.03
const WE_N0      = 6
const WE_MAX_PTS = 55
const WE_FLOOR   = 1.0e-8

"""Our own water impacts: `We;tc;epsilon`, semicolon separated with comma decimals."""
function read_water()
    We, eps = Float64[], Float64[]
    for (i, ln) in enumerate(eachline(joinpath(ROOT, "data", "metrics_Water.csv")))
        i == 1 && continue
        f = split(strip(ln), ';')
        length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f, ',' => '.'))
        any(isnothing, v) && continue
        push!(We, v[1]); push!(eps, v[3])
    end
    We, eps
end

groups(gamma, mu) = (Oh = mu / sqrt(RHO * gamma * R_DROP),
                     Bo = RHO * GRAV * R_DROP^2 / gamma)

"""Restitution curve for a liquid, on the adaptive Weber grid of `_curve.jl`."""
function curve(gamma, mu)
    g = groups(gamma, mu)
    p(We) = ImpactParams(We = We, Bo = g.Bo, Oh = g.Oh, M = M_RUN, K = K_RUN, t_max = T_MAX)
    ev(We) = cor_or_settle(() -> simulate(p(We)), () -> simulate_lcp(p(We)), T_MAX)
    xs, ys = adaptive_curve(ev, WE_LO, WE_HI; th = TH, n0 = WE_N0,
                            maxpts = WE_MAX_PTS, xfloor = WE_FLOOR,
                            tag = @sprintf("γ=%.1f μ=%.1f", 1000gamma, 1000mu))
    xs, ys, g
end

## Wide on purpose. 72.8 is clean water; 30 is about where a saturated surfactant
## takes it, and nothing milder produces a visible effect.
const GAMMAS = [0.0728, 0.050, 0.038, 0.028]
## A full decade, so the lever is visible rather than inferred.
const MUS    = [1.0e-3, 2.0e-3, 4.0e-3, 8.0e-3]
const CG     = [:black, :steelblue, :goldenrod, :indianred]

function main()
    common = (xscale = :log10, xlabel = "Weber number  We", framestyle = :box,
              grid = true, gridalpha = 0.15, tickfontsize = 11, guidefontsize = 13,
              legendfontsize = 10, titlefontsize = 12, ylims = (0.0, 1.0),
              legend = :bottomleft)

    pγ = plot(; ylabel = "restitution  ε",
              title = "Surface tension   (μ = 1.0 mPa·s)", common...)
    pμ = plot(; ylabel = "restitution  ε",
              title = "Viscosity   (γ = 72.8 mN/m)", common...)

    ## Our own water impacts, on both panels: the question the figure exists to
    ## answer is whether either property, pushed to a value an impure liquid could
    ## plausibly reach, moves the model onto them.
    We_e, eps_e = read_water()
    for p in (pγ, pμ)
        scatter!(p, We_e, eps_e; c = :gray35, ms = 4.5, msw = 0, alpha = 0.5,
                 label = "our water experiments (n = $(length(We_e)))")
    end
    @printf("our water experiments: n = %d, We in [%.4g, %.4g], median ε = %.3f\n\n",
            length(We_e), minimum(We_e), maximum(We_e), median(eps_e))

    rows = NTuple{6,Float64}[]

    @printf("Water drop, R = 0.3 mm.  Adaptive Weber grid over [%.3g, %.1f], threads = %d.\n\n",
            WE_LO, WE_HI, Threads.nthreads())

    ## Every curve is independent, so all of them are computed at once and only
    ## then plotted, in order. See `parallel_curves` for what this does and does
    ## not buy.
    jobs = vcat([(:γ, i, g, MU0) for (i, g) in enumerate(GAMMAS)],
                [(:μ, i, m, GAMMA0) for (i, m) in enumerate(MUS)])
    progress_reset()
    done = parallel_curves(j -> j[1] === :γ ? curve(j[3], j[4]) : curve(j[4], j[3]), jobs)

    println("--- surface tension, at mu = 1.0 mPa s ---")
    for (i, g) in enumerate(GAMMAS)
        c_we, c, d = done[i]
        @printf("  γ = %4.1f mN/m -> Oh = %.5f, Bo = %.5f   ε(We=0.5) = %.4f  We_min = %.2g  ->  %s\n",
                1000g, d.Oh, d.Bo, c[argmin(abs.(c_we .- 0.5))], c_we[1], curve_report(c_we, c, TH))
        plot!(pγ, c_we, c; c = CG[i], lw = 2.5, marker = :circle, ms = 3.5,
              label = @sprintf("γ = %.1f mN/m  (Oh = %.4f)", 1000g, d.Oh))
        for (j, We) in enumerate(c_we)
            push!(rows, (1000g, 1000MU0, We, d.Oh, d.Bo, c[j]))
        end
        flush(stdout)
    end

    println("\n--- viscosity, at gamma = 72.8 mN/m ---")
    for (i, m) in enumerate(MUS)
        c_we, c, d = done[length(GAMMAS) + i]
        @printf("  μ = %4.1f mPa·s -> Oh = %.5f, Bo = %.5f   ε(We=0.5) = %.4f  We_min = %.2g  ->  %s\n",
                1000m, d.Oh, d.Bo, c[argmin(abs.(c_we .- 0.5))], c_we[1], curve_report(c_we, c, TH))
        plot!(pμ, c_we, c; c = CG[i], lw = 2.5, marker = :circle, ms = 3.5,
              label = @sprintf("μ = %.1f mPa·s  (Oh = %.4f)", 1000m, d.Oh))
        for (j, We) in enumerate(c_we)
            push!(rows, (1000GAMMA0, 1000m, We, d.Oh, d.Bo, c[j]))
        end
        flush(stdout)
    end

    fig = plot(pγ, pμ; layout = (1, 2), size = (1150, 480), dpi = 200,
               left_margin = 7Plots.mm, bottom_margin = 7Plots.mm, top_margin = 5Plots.mm,
               plot_title = "Water drop, R = 0.3 mm, M = 90, K = 3: can either property move the model onto our water data?",
               plot_titlefontsize = 12)
    out = joinpath(FIGS, "figure_surface_tension_viscosity.png")
    savefig(fig, out)

    open(joinpath(OUT, "figure_surface_tension_viscosity.csv"), "w") do io
        println(io, "gamma_mN_m,mu_mPa_s,We,Oh,Bo,cor")
        for r in rows
            @printf(io, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", r...)
        end
    end
    println("\nwrote ", out)
    println("wrote ", joinpath(OUT, "figure_surface_tension_viscosity.csv"))
end

main()
