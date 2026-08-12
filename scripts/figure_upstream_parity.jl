# This solver against the Gabbard et al. (2025) MATLAB solver, code to code.
#
# Not against their experiments -- against their program, run verbatim from the
# published repository at conditions this package also runs. Two independent
# implementations of the same model should agree to several digits, and where
# they do not, the difference is a statement about one of them rather than about
# the physics.
#
# WHICH OF OUR SOLVERS IS THE COMPARABLE ONE. Their model carries surface modes
# only, with potential flow inside and Lamb's damping. That is exactly K = 1
# here: one radial trial function per mode, x^(l+1), which is irrotational. So
# K = 1 is the parity test and it should match closely. K = 3 resolves the
# vortical layer and carries Reid's exact damping, which their model does not
# have, so it is expected to differ -- the gap between the two columns is the
# physics this package adds, not an error.
#
# THEIR SIDE OF THE TABLE was produced by scripts/upstream_water_matlab.m, run
# against their repository, and is stored in outputs/csv/water_upstream.csv. It is
# not recomputed here: MATLAB is not a dependency of this repository.
#
#   julia --project=docs scripts/figure_upstream_parity.jl
#
# Writes outputs/csv/figure_upstream_parity.csv and
# outputs/figures/figure_upstream_parity.png.

using Printf, Statistics
using Plots
using DropSolver
gr()

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

## The conditions the MATLAB driver was run at: water, R = 0.3 mm.
const OH   = 0.006766649524509585
const BO   = 0.0121
const M_RUN, T_MAX = 90, 25.0

function read_upstream()
    We, cor = Float64[], Float64[]
    for ln in eachline(joinpath(OUT, "water_upstream.csv"))
        f = split(strip(ln), ',')
        length(f) < 2 && continue
        v = tryparse.(Float64, f)
        any(isnothing, v) && continue
        push!(We, v[1]); push!(cor, v[2])
    end
    We, cor
end

function main()
    We, up = read_upstream()
    @printf("Gabbard et al. (2025) MATLAB solver vs this package\n")
    @printf("water: Oh = %.6f, Bo = %.4f, M = %d\n\n", OH, BO, M_RUN)
    @printf("%-9s %-12s %-12s %-10s %-12s %-10s\n",
            "We", "upstream", "ours K=1", "rel diff", "ours K=3", "rel diff")

    k1, k3, d1, d3 = Float64[], Float64[], Float64[], Float64[]
    for (i, w) in enumerate(We)
        r1 = simulate(ImpactParams(We = w, Bo = BO, Oh = OH, M = M_RUN, K = 1, t_max = T_MAX))
        r3 = simulate(ImpactParams(We = w, Bo = BO, Oh = OH, M = M_RUN, K = 3, t_max = T_MAX))
        push!(k1, r1.cor); push!(k3, r3.cor)
        e1 = (r1.cor - up[i]) / up[i]
        e3 = (r3.cor - up[i]) / up[i]
        push!(d1, e1); push!(d3, e3)
        @printf("%-9.4g %-12.6f %-12.6f %+9.3f%%  %-12.6f %+9.3f%%\n",
                w, up[i], r1.cor, 100e1, r3.cor, 100e3)
        flush(stdout)
    end

    ## Their solver is validated over Oh in [0.0139, 0.79]. Water is below that, and
    ## above We ~ 1 their run collapses to cor ~ 0.33 while every mode-count check
    ## says it is converged. Parity is therefore quoted where both solvers are inside
    ## their stated range, and the collapse is shown rather than hidden.
    good = We .<= 0.7
    @printf("\nPARITY over We <= 0.7 (n = %d):\n", sum(good))
    @printf("  K = 1 : median |diff| %.4f %%   worst %.4f %%\n",
            100median(abs.(d1[good])), 100maximum(abs.(d1[good])))
    @printf("  K = 3 : median |diff| %.3f %%    worst %.3f %%\n",
            100median(abs.(d3[good])), 100maximum(abs.(d3[good])))
    nsig = -log10(median(abs.(d1[good])))
    @printf("  agreement at K = 1 is %.1f significant digits\n", nsig)

    open(joinpath(OUT, "figure_upstream_parity.csv"), "w") do io
        println(io, "We,cor_upstream,cor_ours_K1,cor_ours_K3,reldiff_K1,reldiff_K3")
        for i in eachindex(We)
            @printf(io, "%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n",
                    We[i], up[i], k1[i], k3[i], d1[i], d3[i])
        end
    end

    common = (xscale = :log10, xlabel = "Weber number  We", framestyle = :box,
              grid = true, gridalpha = 0.15, tickfontsize = 10, guidefontsize = 12,
              legendfontsize = 9, titlefontsize = 11)

    p1 = plot(; ylabel = "restitution  ε", title = "Water, Oh = 0.0068", common...)
    plot!(p1, We, up; c = :black, lw = 3, marker = :circle, ms = 5,
          label = "Gabbard et al. (2025) solver")
    plot!(p1, We, k1; c = :crimson, lw = 2, ls = :dash, marker = :diamond, ms = 5,
          label = "this package, K = 1 (same model)")
    plot!(p1, We, k3; c = :steelblue, lw = 2, marker = :utriangle, ms = 5,
          label = "this package, K = 3 (resolved interior)")
    vline!(p1, [0.7]; c = :gray60, ls = :dot, lw = 1.5, label = "")
    annotate!(p1, 1.6, 0.55, text("their run collapses\nhere, ours does not", 8, :center, :gray40))

    p2 = plot(; ylabel = "difference from their solver", yscale = :log10,
              title = "K = 1 is the same model, and matches it", common...)
    plot!(p2, We, max.(abs.(d1), 1e-9); c = :crimson, lw = 2, marker = :diamond, ms = 5,
          label = "K = 1")
    plot!(p2, We, max.(abs.(d3), 1e-9); c = :steelblue, lw = 2, marker = :utriangle, ms = 5,
          label = "K = 3")
    hline!(p2, [1e-2]; c = :gray50, ls = :dash, lw = 1.5, label = "1 % (2 significant digits)")
    vline!(p2, [0.7]; c = :gray60, ls = :dot, lw = 1.5, label = "")

    fig = plot(p1, p2; layout = (1, 2), size = (1150, 470), dpi = 200,
               left_margin = 7Plots.mm, bottom_margin = 7Plots.mm, top_margin = 4Plots.mm)
    out = joinpath(FIGS, "figure_upstream_parity.png")
    savefig(fig, out)
    println("\nwrote ", out)
    println("wrote ", joinpath(OUT, "figure_upstream_parity.csv"))
end

main()
