# The solver against water rebounding from a Glaco superhydrophobic surface.
#
# WHY THIS DATASET. The Gabbard et al. (2025) campaign, which is the package's
# main Newtonian check, spans Oh in [0.014, 0.79]. Water at the drop sizes used
# in the concentration series sits at Oh = 0.0068, below that floor, so every
# statement the documentation makes about water is an extrapolation. This dataset
# reaches Oh = 0.0031, and 119 of its 174 impacts lie below Gabbard's floor, so it
# tests exactly the regime the other campaign does not.
#
# HOW THE COMPARISON IS SET UP. Drop radius varies across the dataset by a factor
# of sixty, and both dimensionless groups are tied to it: Oh ~ R^-1/2 and
# Bo ~ R^2. Binning by Ohnesorge therefore bins radius, and with it Bond, so each
# band is a single physical drop size swept in impact speed. That is the
# comparison the model can make: one (Oh, Bo), a sweep in We.
#
# Within a band the experiments are grouped into Weber bins of at least MIN_REPS
# repetitions, and the group median is compared against a run at the band's median
# Oh and Bo. Comparing against individual shots would report the experimental
# scatter as model error.
#
#   julia --project=docs scripts/validate_glaco.jl
#
# Writes outputs/csv/validate_glaco.csv and outputs/figures/validate_glaco.png.

using Printf, Statistics
using Plots
using DropSolver
gr()

const ROOT = joinpath(@__DIR__, "..")
const DATA = joinpath(ROOT, "data", "glaco_restitution.csv")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

const M_RUN   = 90        # production resolution
const K_RUN   = 3
const T_MAX   = 25.0
const MIN_REPS = 4
const WE_MAX  = 4.0       # above this the linearised model is out of its range

## Bands in Ohnesorge. The second band of the raw data (7 points, all at We > 10)
## is not represented here: it has no overlap with the Weber range the model can
## be asked about, so a band would be drawn through nothing.
const OH_EDGES = [(0.0030, 0.0055), (0.0085, 0.0130), (0.0130, 0.0190), (0.0190, 0.0250)]

function read_glaco()
    We, C, Oh, Bo = Float64[], Float64[], Float64[], Float64[]
    for ln in eachline(DATA)
        (startswith(ln, "#") || startswith(ln, "We")) && continue
        f = parse.(Float64, split(strip(ln), ','))
        ## cor = 0 means no rebound was recorded; those are not restitution
        ## measurements and are excluded rather than averaged in as zeros.
        f[2] > 0 || continue
        push!(We, f[1]); push!(C, f[2]); push!(Oh, f[3]); push!(Bo, f[4])
    end
    We, C, Oh, Bo
end

"""Groups of at least `MIN_REPS` experiments sharing a logarithmic Weber bin."""
function wegroups(We, C, idx, nbins)
    isempty(idx) && return NTuple{4,Float64}[]
    lo, hi = log10(minimum(We[idx])), log10(maximum(We[idx]))
    hi > lo || return NTuple{4,Float64}[]
    edges = range(lo, hi; length = nbins + 1)
    out = NTuple{4,Float64}[]
    for b in 1:nbins
        g = [i for i in idx if edges[b] <= log10(We[i]) <= edges[b+1] + (b == nbins ? 1e-9 : 0)]
        length(g) >= MIN_REPS || continue
        push!(out, (median(We[g]), median(C[g]), length(g) > 1 ? std(C[g]) : 0.0, length(g)))
    end
    out
end

function main()
    We, C, Oh, Bo = read_glaco()
    @printf("%d usable impacts, Oh in [%.4f, %.4f]\n\n", length(We), minimum(Oh), maximum(Oh))

    plt = plot(xscale = :log10, xlabel = "Weber number  We",
               ylabel = "restitution  ε", legend = :bottomleft,
               size = (1000, 700), dpi = 160, framestyle = :axes, grid = false,
               guidefontsize = 17, tickfontsize = 15, legendfontsize = 12,
               ylims = (0.0, 1.0), foreground_color_axis = :gray40,
               foreground_color_border = :gray40, left_margin = 10Plots.mm,
               bottom_margin = 10Plots.mm, right_margin = 6Plots.mm, top_margin = 6Plots.mm)
    cols = [:steelblue, :seagreen, :goldenrod, :indianred]

    rows = NTuple{7,Float64}[]
    errs = Float64[]
    for (bi, (lo, hi)) in enumerate(OH_EDGES)
        idx = [i for i in eachindex(We) if lo <= Oh[i] < hi && We[i] <= WE_MAX]
        length(idx) >= MIN_REPS || continue
        oh_b, bo_b = median(Oh[idx]), median(Bo[idx])
        groups = wegroups(We, C, idx, 6)
        isempty(groups) && continue

        @printf("Oh band [%.4f,%.4f)  model at Oh = %.4f, Bo = %.4f   (n = %d)\n",
                lo, hi, oh_b, bo_b, length(idx))
        @printf("  %-10s %-6s %-18s %-10s %s\n", "We", "n", "exp (med±sd)", "model", "err")

        sim_We, sim_C = Float64[], Float64[]
        for (we_g, c_g, sd_g, n_g) in groups
            r = simulate(ImpactParams(We = we_g, Bo = bo_b, Oh = oh_b,
                                      M = M_RUN, K = K_RUN, t_max = T_MAX))
            e = isfinite(r.cor) ? (r.cor - c_g) / c_g : NaN
            isfinite(e) && push!(errs, abs(e))
            isfinite(r.cor) && (push!(sim_We, we_g); push!(sim_C, r.cor))
            push!(rows, (oh_b, bo_b, we_g, n_g, c_g, sd_g, isfinite(r.cor) ? r.cor : NaN))
            @printf("  %-10.4g %-6d %.4f ± %.4f     %-10.4f %+.1f%%\n",
                    we_g, n_g, c_g, sd_g, r.cor, 100e)
            flush(stdout)
        end

        lbl = @sprintf("Oh = %.4f", oh_b)
        scatter!(plt, We[idx], C[idx]; c = cols[bi], ms = 5, msw = 0, alpha = 0.35, label = "")
        !isempty(sim_We) && plot!(plt, sim_We, sim_C; c = cols[bi], lw = 4, label = lbl)
        println()
    end

    open(joinpath(OUT, "validate_glaco.csv"), "w") do io
        println(io, "Oh,Bo,We,n_exp,cor_exp,cor_exp_sd,cor_sim")
        for r in rows
            @printf(io, "%.8g,%.8g,%.8g,%d,%.8g,%.8g,%.8g\n",
                    r[1], r[2], r[3], Int(r[4]), r[5], r[6], r[7])
        end
    end

    out = joinpath(FIGS, "validate_glaco.png")
    savefig(plt, out)
    println("wrote ", out)
    println("wrote ", joinpath(OUT, "validate_glaco.csv"))
    @printf("\nGLACO SUMMARY (M = %d, K = %d)\n", M_RUN, K_RUN)
    @printf("  CoR median |err| %.1f%%   mean %.1f%%   worst %.1f%%   (n = %d Weber groups)\n",
            100median(errs), 100mean(errs), 100maximum(errs), length(errs))
end

main()
