# How far does restitution move if the liquid is not quite water?
#
# THE QUESTION. A water drop that rebounds lower than predicted may be a failure
# of the model, or it may not be water. Surfactant picked up from glassware, from
# the air, or from the syringe lowers surface tension; dissolved or suspended
# material raises viscosity. This figure puts a number on both, so that a measured
# discrepancy can be compared against the contamination that would explain it.
#
# WHY THE AXES ARE PHYSICAL AND NOT DIMENSIONLESS. An experiment fixes the drop
# radius and the impact speed; it does not fix Weber or Ohnesorge. Surface tension
# appears in all three groups the model takes,
#
#     We = rho U^2 R / gamma ,   Oh = mu / sqrt(rho gamma R) ,   Bo = rho g R^2 / gamma ,
#
# so lowering gamma at fixed (R, U) raises all three at once, and the change in
# restitution is the net of three effects rather than any one of them. Plotting
# against gamma keeps that coupling; plotting against We would hide it.
#
# WHAT THIS DOES NOT COVER. The model carries surface tension as a constant. A real
# surfactant also brings Marangoni stresses from gradients in surface coverage, and
# a surface viscosity, neither of which is here. Both oppose the flow and so remove
# further energy. The curves below are therefore a LOWER BOUND on what a given
# contamination does to restitution, not an estimate of it.
#
#   julia --project=docs scripts/figure_contamination.jl
#
# Writes outputs/csv/figure_contamination.csv and outputs/figures/figure_contamination.png.

using Printf, Statistics
using Plots
using DropSolver
gr()

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

const M_RUN, K_RUN, T_MAX = 90, 3, 25.0

## Clean water, and the drop the concentration series uses.
const RHO    = 997.0
const GRAV   = 9.81
const R_DROP = 3.0e-4        # m
const GAMMA0 = 0.0728        # N/m
const MU0    = 1.0e-3        # Pa s

"""Dimensionless groups for a drop of radius `R` hitting at speed `U`."""
dims(gamma, mu, U; R = R_DROP) =
    (We = RHO * U^2 * R / gamma, Oh = mu / sqrt(RHO * gamma * R), Bo = RHO * GRAV * R^2 / gamma)

"""Impact speed that gives Weber number `We` in clean water."""
speed_for(We) = sqrt(We * GAMMA0 / (RHO * R_DROP))

run_cor(gamma, mu, U) = let d = dims(gamma, mu, U)
    r = simulate(ImpactParams(We = d.We, Bo = d.Bo, Oh = d.Oh,
                              M = M_RUN, K = K_RUN, t_max = T_MAX))
    r.cor
end

const WE_REF = [0.2, 0.5, 1.0]
const SPEEDS = speed_for.(WE_REF)
const COLS   = [:steelblue, :goldenrod, :indianred]

## Surfactant contamination: pure water is 72 mN/m; a saturated common surfactant
## takes it to about 35. Viscosity: the range a dissolved contaminant plausibly spans.
const GAMMAS = collect(range(GAMMA0, 0.035; length = 10))
const MUS    = collect(range(MU0, 2.2e-3; length = 10))

function main()
    rows = NTuple{7,Float64}[]

    pγ = plot(xlabel = "surface tension  γ  (mN/m)", ylabel = "restitution  ε",
              legend = :bottomright, framestyle = :box, grid = true, gridalpha = 0.15,
              guidefontsize = 12, tickfontsize = 11, legendfontsize = 10,
              title = "Contamination lowers γ", titlefontsize = 12, xflip = true)
    pμ = plot(xlabel = "viscosity  μ  (mPa·s)", ylabel = "restitution  ε",
              legend = :bottomleft, framestyle = :box, grid = true, gridalpha = 0.15,
              guidefontsize = 12, tickfontsize = 11, legendfontsize = 10,
              title = "Contamination raises μ", titlefontsize = 12)

    println("Drop radius $(1000R_DROP) mm.  Clean water: γ = $(1000GAMMA0) mN/m, μ = $(1000MU0) mPa·s.\n")

    for (i, (We_ref, U)) in enumerate(zip(WE_REF, SPEEDS))
        @printf("=== We = %.1f in clean water  (U = %.1f cm/s) ===\n", We_ref, 100U)

        cγ = Float64[]
        for g in GAMMAS
            c = run_cor(g, MU0, U); push!(cγ, c)
            d = dims(g, MU0, U)
            push!(rows, (We_ref, 1000g, 1000MU0, d.We, d.Oh, d.Bo, c))
        end
        cμ = Float64[]
        for m in MUS
            c = run_cor(GAMMA0, m, U); push!(cμ, c)
            d = dims(GAMMA0, m, U)
            push!(rows, (We_ref, 1000GAMMA0, 1000m, d.We, d.Oh, d.Bo, c))
        end

        plot!(pγ, 1000 .* GAMMAS, cγ; c = COLS[i], lw = 2.5, marker = :circle, ms = 3,
              label = @sprintf("We = %.1f", We_ref))
        plot!(pμ, 1000 .* MUS, cμ; c = COLS[i], lw = 2.5, marker = :circle, ms = 3,
              label = @sprintf("We = %.1f", We_ref))

        ## Logarithmic sensitivity at the clean-water point: the fractional change in
        ## restitution per fractional change in the property. This is the number to
        ## carry away -- it converts an observed shortfall into a required contamination.
        sγ = (log(cγ[2]) - log(cγ[1])) / (log(GAMMAS[2]) - log(GAMMAS[1]))
        sμ = (log(cμ[2]) - log(cμ[1])) / (log(MUS[2]) - log(MUS[1]))
        @printf("  clean ε = %.4f\n", cγ[1])
        @printf("  γ 72 -> 50 mN/m (-31%%): ε %.4f -> %.4f  (%+.1f%%)\n",
                cγ[1], cγ[findmin(abs.(GAMMAS .- 0.050))[2]],
                100 * (cγ[findmin(abs.(GAMMAS .- 0.050))[2]] / cγ[1] - 1))
        @printf("  μ 1.0 -> 2.0 mPa·s (+100%%): ε %.4f -> %.4f  (%+.1f%%)\n",
                cμ[1], cμ[findmin(abs.(MUS .- 2.0e-3))[2]],
                100 * (cμ[findmin(abs.(MUS .- 2.0e-3))[2]] / cμ[1] - 1))
        @printf("  dlnε/dlnγ = %+.3f     dlnε/dlnμ = %+.3f\n\n", sγ, sμ)
        flush(stdout)
    end

    fig = plot(pγ, pμ; layout = (1, 2), size = (1100, 460), dpi = 200,
               left_margin = 6Plots.mm, bottom_margin = 6Plots.mm, top_margin = 4Plots.mm,
               plot_title = @sprintf("Water drop, R = %.1f mm: how far ε moves if the liquid is not clean water", 1000R_DROP),
               plot_titlefontsize = 11)
    out = joinpath(FIGS, "figure_contamination.png")
    savefig(fig, out)

    open(joinpath(OUT, "figure_contamination.csv"), "w") do io
        println(io, "We_clean,gamma_mN_m,mu_mPa_s,We,Oh,Bo,cor")
        for r in rows
            @printf(io, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g\n", r...)
        end
    end
    println("wrote ", out)
    println("wrote ", joinpath(OUT, "figure_contamination.csv"))
end

main()
