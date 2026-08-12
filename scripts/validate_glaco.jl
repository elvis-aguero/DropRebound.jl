# The solver against Thenarianto et al. (2023): water on a Glaco superhydrophobic
# surface.
#
# WHY THIS DATASET. Water is already checked against our own 53 impacts, which the
# concentration series reproduces to 8.6 per cent. What it was not checked against
# was a second, independent campaign at the same Ohnesorge: the Gabbard et al.
# (2025) data spans Oh in [0.014, 0.79], and water at the drop sizes used here sits
# at Oh = 0.0068, below that floor. Thenarianto et al. reach Oh = 0.0031, and 119
# of their 174 impacts lie below Gabbard's floor, so they cover the regime the
# other campaign does not.
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

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
BLAS.set_num_threads(1)   # one BLAS thread per worker thread; see parallel_curves
include(joinpath(@__DIR__, "_curve.jl"))

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

## Adaptive Weber sampling; see scripts/_curve.jl. One threshold sets both
## guarantees: no consecutive pair of samples differs by more than TH in
## restitution, and each curve is followed down until restitution falls below TH.
const TH         = 0.03
const WE_N0      = 6
const WE_MAX_PTS = 55
const WE_FLOOR   = 1.0e-8

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

"""Model results already computed, keyed by band Ohnesorge. Empty if none stored."""
function cached()
    path = joinpath(OUT, "validate_glaco.csv")
    isfile(path) || return Dict{Float64,Vector{NTuple{2,Float64}}}()
    d = Dict{Float64,Vector{NTuple{2,Float64}}}()
    for ln in eachline(path)
        startswith(ln, "Oh") && continue
        f = tryparse.(Float64, split(strip(ln), ','))
        (length(f) < 7 || any(isnothing, f[[1, 3, 7]])) && continue
        isfinite(f[7]) || continue
        push!(get!(d, f[1], NTuple{2,Float64}[]), (f[3], f[7]))
    end
    d
end

function main()
    ## `--replot` redraws from outputs/csv/validate_glaco.csv without simulating, so
    ## that a change to the figure does not cost a fresh sweep at M = 90.
    replot = "--replot" in ARGS
    store = replot ? cached() : Dict{Float64,Vector{NTuple{2,Float64}}}()
    replot && isempty(store) && error("--replot given but outputs/csv/validate_glaco.csv has no usable rows")

    We, C, Oh, Bo = read_glaco()
    @printf("%d usable impacts, Oh in [%.4f, %.4f]%s\n\n", length(We), minimum(Oh), maximum(Oh),
            replot ? "  (replot: no simulations run)" : "")

    plt = plot(xscale = :log10, xlabel = "Weber number  We",
               ylabel = "restitution  ε", legend = :bottomleft,
               title = "Thenarianto et al. (2023), water on Glaco  vs  this solver (M = 90, K = 3)\npoints: experiment      lines: model",
               titlefontsize = 11,
               size = (1000, 700), dpi = 160, framestyle = :axes, grid = false,
               guidefontsize = 14, tickfontsize = 13, legendfontsize = 11,
               ylims = (0.0, 1.0), foreground_color_axis = :gray40,
               foreground_color_border = :gray40, left_margin = 10Plots.mm,
               bottom_margin = 10Plots.mm, right_margin = 6Plots.mm, top_margin = 6Plots.mm)
    cols = [:steelblue, :seagreen, :goldenrod, :indianred]

    ## Each band is one drop size and they do not interact, so all of them are swept
    ## at once and only then printed and plotted, in order.
    bands = NamedTuple[]
    for (bi, (lo, hi)) in enumerate(OH_EDGES)
        idx = [i for i in eachindex(We) if lo <= Oh[i] < hi && We[i] <= WE_MAX]
        length(idx) >= MIN_REPS || continue
        groups = wegroups(We, C, idx, 6)
        isempty(groups) && continue
        push!(bands, (bi = bi, lo = lo, hi = hi, idx = idx,
                      oh = median(Oh[idx]), bo = median(Bo[idx]), groups = groups))
    end
    @printf("%d bands, threads = %d\n\n", length(bands), Threads.nthreads())

    progress_reset()
    swept = parallel_curves(bands) do b
        ## The stored Ohnesorge is the band median written to eight digits, so it does
        ## not compare equal to the one recomputed here. Match the nearest key instead.
        hit = if replot && !isempty(store)
            k = collect(keys(store))
            store[k[argmin(abs.(k .- b.oh))]]
        else
            NTuple{2,Float64}[]
        end
        ## The model is drawn as a curve, not as a polyline through the handful of
        ## Weber groups the experiments happen to form: at one to six points a band,
        ## that is a sketch of the model rather than the model. So the band's own
        ## Weber span is swept adaptively, with the group medians forced in so the
        ## error column below is still an exact comparison rather than an
        ## interpolation.
        ev(w) = if replot
            j = findfirst(t -> isapprox(t[1], w; rtol = 1e-6), hit)
            j === nothing ? NaN : hit[j][2]
        else
            cor_or_settle(simulate(ImpactParams(We = w, Bo = b.bo, Oh = b.oh,
                                                M = M_RUN, K = K_RUN, t_max = T_MAX)), T_MAX)
        end
        if replot
            xs = sort([t[1] for t in hit])
            (xs, [ev(x) for x in xs])
        else
            adaptive_curve(ev, minimum(We[b.idx]), min(maximum(We[b.idx]), WE_MAX);
                           th = TH, forced = [g[1] for g in b.groups],
                           n0 = WE_N0, maxpts = WE_MAX_PTS, xfloor = WE_FLOOR,
                           tag = @sprintf("Bo=%.4f", b.bo))
        end
    end

    rows = NTuple{7,Float64}[]
    errs = Float64[]
    for (bn, b) in enumerate(bands)
        bi, idx, oh_b, bo_b, groups = b.bi, b.idx, b.oh, b.bo, b.groups
        sim_We, sim_C = swept[bn]

        @printf("Oh band [%.4f,%.4f)  model at Oh = %.4f, Bo = %.4f   (n = %d)\n",
                b.lo, b.hi, oh_b, bo_b, length(idx))
        @printf("  %-10s %-6s %-18s %-10s %s\n", "We", "n", "exp (med±sd)", "model", "err")

        for (we_g, c_g, sd_g, n_g) in groups
            j = findfirst(x -> isapprox(x, we_g; rtol = 1e-6), sim_We)
            cor = j === nothing ? NaN : sim_C[j]
            e = isfinite(cor) ? (cor - c_g) / c_g : NaN
            isfinite(e) && push!(errs, abs(e))
            @printf("  %-10.4g %-6d %.4f ± %.4f     %-10.4f %+.1f%%\n",
                    we_g, n_g, c_g, sd_g, cor, 100e)
            flush(stdout)
        end
        for (x, y) in zip(sim_We, sim_C)
            k = findfirst(g -> isapprox(g[1], x; rtol = 1e-6), groups)
            push!(rows, (oh_b, bo_b, x, k === nothing ? 0 : groups[k][4],
                         k === nothing ? NaN : groups[k][2],
                         k === nothing ? NaN : groups[k][3], y))
        end
        println("  sampler: ", curve_report(sim_We, sim_C, TH))
        keep = isfinite.(sim_C)
        sim_We, sim_C = sim_We[keep], sim_C[keep]

        ## One legend entry per band rather than two: the colour carries the band, and
        ## the marker-versus-line distinction is stated once in the axis label. A band
        ## whose experiments support only one Weber group gets a marker, since a line
        ## through one point draws nothing at all.
        lbl = @sprintf("Bo = %.4f  (Oh = %.4f, R = %.0f µm, n = %d)",
                       bo_b, oh_b, 1e6 * sqrt(bo_b * 0.0728 / (997 * 9.81)), length(idx))
        scatter!(plt, We[idx], C[idx]; c = cols[bi], ms = 5, msw = 0, alpha = 0.32, label = "")
        if length(sim_We) >= 2
            plot!(plt, sim_We, sim_C; c = cols[bi], lw = 4, label = lbl)
        elseif !isempty(sim_We)
            scatter!(plt, sim_We, sim_C; c = cols[bi], ms = 11, marker = :star5,
                     msw = 0, label = lbl)
        end
        println()
    end

    ## Never rewritten in replot mode: the stored results are the input there, and
    ## an earlier version of this path overwrote them with NaN.
    replot || open(joinpath(OUT, "validate_glaco.csv"), "w") do io
        println(io, "Oh,Bo,We,n_exp,cor_exp,cor_exp_sd,cor_sim")
        ## Rows with n_exp = 0 are curve points with no experimental group at that
        ## Weber number; their cor_exp columns are NaN by construction.
        for r in rows
            @printf(io, "%.8g,%.8g,%.8g,%d,%.8g,%.8g,%.8g\n",
                    r[1], r[2], r[3], Int(r[4]), r[5], r[6], r[7])
        end
    end

    out = joinpath(FIGS, "validate_glaco.png")
    savefig(plt, out)
    println("wrote ", out)
    replot || println("wrote ", joinpath(OUT, "validate_glaco.csv"))
    @printf("\nGLACO SUMMARY (M = %d, K = %d)\n", M_RUN, K_RUN)
    @printf("  CoR median |err| %.1f%%   mean %.1f%%   worst %.1f%%   (n = %d Weber groups)\n",
            100median(errs), 100mean(errs), 100maximum(errs), length(errs))
end

main()
