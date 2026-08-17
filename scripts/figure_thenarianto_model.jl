# This solver against Thenarianto et al. (2023), one figure per surface.
#
# WHAT THE MODEL IS AND IS NOT TOLD. It is given We, Oh and Bo. It is told nothing
# about the surface -- there is no contact angle here, no hysteresis, no work of
# adhesion -- so identical (We, Oh, Bo) returns an identical curve on Glaco and on
# black silicon. Any difference between the two datasets is therefore a lower bound
# on the physics this model does not carry, and the paper says what that physics is:
# contact-line pinning, with e = sqrt(1 - 6*alpha/We) below the inertial-capillary
# regime and alpha a measured property of the surface, 2.9e-3 for Glaco against
# 1.1e-4 for black silicon.
#
# WHY THE SMALL DROPS GET A BAND AND THE LARGE ONES A LINE. Each surface carries two
# droplet populations. The dispensed millimetre drops span a factor of 1.1 in
# Ohnesorge and 1.3 in Bond, which is a single point as far as the model is
# concerned, so one curve represents them honestly. The sprayed submillimetre drops
# span a factor of four in Ohnesorge and four hundred in Bond, and Bond is already
# known to change the answer, so one curve there would be a fiction. They get three:
# the smallest, median and largest radius in the group, drawn as a band.
#
# That band is also the test of the paper's central claim. They report that
# restitution collapses onto a single function of We for a given surface, whatever
# the drop size. If the three curves lie on top of each other the model reproduces
# that collapse; if they fan out, the model predicts a size dependence the
# experiments deny. Either way it is a result rather than a nuisance.
#
#   julia --project=docs -t 6 scripts/figure_thenarianto_model.jl
#
# Writes outputs/csv/figure_thenarianto_model.csv and
# outputs/figures/figure_thenarianto_model_{glaco,bsi}.png.

using Printf, Statistics, LinearAlgebra
using Plots
using DropSolver
gr()
BLAS.set_num_threads(1)   # one BLAS thread per worker thread; see parallel_curves
include(joinpath(@__DIR__, "_curve.jl"))
include(joinpath(@__DIR__, "_runcache.jl"))

const ROOT = joinpath(@__DIR__, "..")
const OUT  = joinpath(ROOT, "outputs", "csv")
const FIGS = joinpath(ROOT, "outputs", "figures")
mkpath(OUT); mkpath(FIGS)

## `--preview` renders the same figure at throwaway resolution, to judge the layout
## without waiting on a production sweep. It is never the published artefact: the
## output is written beside the real one with a _preview suffix and says so on its face.
const PREVIEW = "--preview" in ARGS
const M_RUN, K_RUN, T_MAX = PREVIEW ? (30, 2, 15.0) : (90, 3, 25.0)

## One threshold, two guarantees; see scripts/_curve.jl.
const TH         = 0.03
const WE_N0      = 6
const WE_MAX_PTS = PREVIEW ? 14 : 150
const WE_X_TH    = 0.02   # decades; the x sensitivity, and a stopping rule
const WE_X_MAX   = PREVIEW ? 0.5 : 0.10   # decades; longest segment drawn
const WE_CEIL    = 5.0    # the linearised model is out of its range above this

const R_SPLIT = 0.8e-3    # m; the gap between the two populations is wide

## Water, as the workbooks assume.
const RHO, SIGMA, MU, GRAV = 997.0, 0.072, 1.0e-3, 9.81
oh_of(R) = MU / sqrt(RHO * SIGMA * R)
bo_of(R) = RHO * GRAV * R^2 / SIGMA

"""
Per-point marker size and colour for one population's experimental scatter,
encoding Bond and Ohnesorge directly rather than folding every drop in a
population into one size and one colour.

Size follows `sqrt(Bo)` -- proportional to `R` itself, since rho, g and sigma are
fixed across one surface's water drops -- normalised within the population's OWN
range. The two populations are not put on one shared scale: the header already
notes the sprayed population's Bond spans a factor of four hundred and the
dispensed one only 1.3, and a shared scale would flatten the interesting spread
in the first population to show the (already known, already why there are two
populations) gap between them instead. Colour follows Oh the same way, walked
along a gradient anchored on one hue per population, so a data point can be read
against the smallest/median/largest-R model curves it should sit nearest to.
"""
function bo_oh_style(R, grad; ms_lo = 3.5, ms_hi = 8.5)
    Bo, Oh = bo_of.(R), oh_of.(R)
    rt(x, lo, hi) = lo == hi ? 0.5 : clamp((x - lo) / (hi - lo), 0.0, 1.0)
    slo, shi = extrema(sqrt.(Bo))
    olo, ohi = extrema(Oh)
    ms = [ms_lo + (ms_hi - ms_lo) * rt(sqrt(b), slo, shi) for b in Bo]
    mc = [get(grad, rt(o, olo, ohi)) for o in Oh]
    (; ms, mc)
end

const GRAD_SML = cgrad([:lightskyblue1, :steelblue, :midnightblue])
const GRAD_BIG = cgrad([:navajowhite, :indianred, :darkred])

const SURFACES = [
    (key = "glaco", csv = "glaco_restitution.csv", title = "Glaco coating", alpha = 2.9e-3),
    (key = "bsi",   csv = "bsi_restitution.csv",   title = "Black silicon", alpha = 1.1e-4),
]

function read_data(csv)
    We, C, R = Float64[], Float64[], Float64[]
    for ln in eachline(joinpath(ROOT, "data", csv))
        (startswith(ln, "#") || startswith(ln, "We")) && continue
        v = parse.(Float64, split(strip(ln), ','))
        push!(We, v[1]); push!(C, v[2]); push!(R, v[5])
    end
    (; We, C, R)
end

"""
One model curve at the radius `R`, scored two ways from the same runs.

  theirs   contact where the drop actually touches, energy ratio.
           Thenarianto et al. measure `e = U_R/U_I = sqrt(h_n/h_{n-1})`, referenced to
           apex heights and so gravity-corrected; there is no camera-pixel threshold
           in their definition.

  ours-2   contact at 0.02R, energy ratio. What this package now reports by default,
           and what Gabbard et al. publish as `coef_restitution_exp`.

Restitution is post-processing, so both come from one simulation. The refinement
follows `theirs`, since that is the curve being compared against their data.
"""
function curve(R, lo, hi, tag)
    Oh, Bo = oh_of(R), bo_of(R)
    mk(w) = ImpactParams(We = w, Bo = Bo, Oh = Oh, M = M_RUN, K = K_RUN, t_max = T_MAX)
    ours = Dict{Float64,NTuple{2,Float64}}()   # We -> (our cor, the We it is measured at)
    lk = ReentrantLock()

    function ev(w)
        p = mk(w)
        s = cached_series(p)                   # simulated once, ever; see _runcache.jl
        m = score(p, s; h_thresh = 0.02)       # ours = ours-2, the only metric now

        ## A drop that never releases has no restitution to report, and NaN would drop
        ## the point from the curve. It is drawn at zero instead, which is what the drop
        ## did: it stayed down. That is the roll-off, not a gap in the data.
        ##
        ## The 0.02R curve is also plotted against the Weber number that criterion
        ## actually measures, and is withheld entirely where the ancestor would refuse
        ## the case -- below `We = 2 Bo h`, no free-flight state exists at the line.
        c = (m.released && isfinite(m.cor)) ? m.cor : 0.0
        lock(lk) do
            ours[w] = (c, is_measurable(p, 0.02) ? we_measured(p, 0.02) : NaN)
        end
        c
    end

    xs, ys = adaptive_curve(ev, lo, min(hi, WE_CEIL); th = TH, n0 = WE_N0,
                            maxpts = WE_MAX_PTS, x_th = WE_X_TH, x_max = WE_X_MAX, tag = tag)

    ## WE_CEIL is a declared boundary of trust, not a measured one -- see the header.
    ## When the data reach past it, the ceiling becomes the last sampled point, and if
    ## the model is already breaking down there, that sample is the one place a real
    ## discontinuity could be mistaken for a resolved result. Observed on the black
    ## silicon large-drop curve: a smooth decline out to We = 4.88, then a step to
    ## We = 5.0 nine times steeper than the trend it interrupts. Trust the disclaimer
    ## over the sample.
    hi > WE_CEIL && ((xs, ys) = trim_ceiling_jump(xs, ys, WE_CEIL, TH; x_th = WE_X_TH))

    got = [get(ours, x, (NaN, NaN)) for x in xs]
    (; R, Oh, Bo, xs, ys, ys2 = first.(got), xs2 = last.(got))
end

"""
The radii to run for one population: one for the tight group, three for the wide one.

The lower end of the sweep is NOT simply the smallest Weber number measured. The model
cuts off at `We = 2 Bo h`, below which the rebound cannot carry the drop back to the
measurement line and restitution is zero by construction, and for the smaller drops that
sits well below anything the experiments reach -- 3.1e-6 against a smallest measurement of
5.6e-5 at R = 0.024 mm. A sweep that starts at the data therefore never sees its own cutoff
and reports a curve that is flat all the way down. Each curve is taken below its own.
"""
function plan(d, big)
    m = big ? (d.R .>= R_SPLIT) : (d.R .< R_SPLIT)
    R, We = d.R[m], d.We[m]
    cut = 0.04 * bo_of(minimum(R))          # the SMALLEST cutoff in the group: smallest R,
                                            # smallest Bo. Going below that one puts every
                                            # curve of the group past its own.
    lo, hi = min(minimum(We), 0.3 * cut), maximum(We)
    radii = big ? [median(R)] : [minimum(R), median(R), maximum(R)]
    (; radii, lo, hi, n = sum(m), rlo = minimum(R), rhi = maximum(R))
end

function main()
    jobs, meta = Any[], Any[]
    for s in SURFACES, big in (true, false)
        p = plan(read_data(s.csv), big)
        for (i, R) in enumerate(p.radii)
            push!(jobs, (R = R, lo = p.lo, hi = p.hi,
                         tag = @sprintf("%s/%s/R=%.3fmm", s.key, big ? "mm" : "um", 1000R)))
            push!(meta, (key = s.key, big = big, idx = i, nr = length(p.radii)))
        end
    end
    @printf("%d model curves, M = %d, K = %d, threads = %d\n\n",
            length(jobs), M_RUN, K_RUN, Threads.nthreads())

    progress_reset()
    curves = parallel_curves(j -> curve(j.R, j.lo, j.hi, j.tag), jobs)

    PREVIEW || open(joinpath(OUT, "figure_thenarianto_model.csv"), "w") do io
        println(io, "surface,population,R_m,Oh,Bo,We,cor_theirs,We_ours,cor_ours")
        for (m, c) in zip(meta, curves), (x, y, x2, y2) in zip(c.xs, c.ys, c.xs2, c.ys2)
            @printf(io, "%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n",
                    m.key, m.big ? "large" : "small", c.R, c.Oh, c.Bo, x, y, x2, y2)
        end
    end

    for s in SURFACES
        d = read_data(s.csv)
        big, sml = d.R .>= R_SPLIT, d.R .< R_SPLIT

        ## The plotted range has to hold the whole roll-off of the median-R curve, the
        ## solid line the comparison is actually made against, not just the range the
        ## experiment covers. That curve's own cutoff, `2*Bo*h`, can sit below the
        ## smallest measured We -- the median BSI drop rolls off at We = 2.1e-5 against
        ## a smallest measurement of 4.2e-5 -- and clipping the x-axis at the data
        ## floor there cuts the roll-off itself out of the figure: the solid curve
        ## would appear to start already at its plateau, with no rise to show for it.
        we_roll_med = 2 * bo_of(median(d.R[sml])) * 0.02
        lo, hi = min(minimum(d.We), 0.5 * we_roll_med), maximum(d.We)
        tks = [10.0^k for k in floor(Int, log10(lo)):ceil(Int, log10(hi))]

        plt = plot(xscale = :log10, xlabel = "Weber number  We", ylabel = "restitution  ε",
                   xticks = tks,
                   xlims = (10.0^(log10(lo) - 0.25), 10.0^(log10(hi) + 0.25)),
                   title = @sprintf("%s  —  experiment vs this solver (M = %d, K = %d)%s",
                                    s.title, M_RUN, K_RUN, PREVIEW ? "   [PREVIEW: low resolution]" : ""),
                   titlefontsize = 11, size = (960, 660), dpi = 200,
                   framestyle = :axes, grid = false, legend = :bottomright,
                   legendfontsize = 9, guidefontsize = 13, tickfontsize = 12,
                   ylims = (0, 1), foreground_color_axis = :gray40,
                   foreground_color_border = :gray40, left_margin = 8Plots.mm,
                   bottom_margin = 8Plots.mm, right_margin = 5Plots.mm, top_margin = 5Plots.mm)

        ## Model first, so the measurements sit on top of it.
        for (m, c) in zip(meta, curves)
            m.key == s.key || continue
            keep = isfinite.(c.ys)
            sum(keep) >= 2 || continue
            col = m.big ? :indianred : :steelblue
            main = m.nr == 1 || m.idx == 2
            ## Solid: their metric, which is the comparison being made.
            plot!(plt, c.xs[keep], c.ys[keep];
                  c = col, lw = main ? 3 : 1.5, ls = main ? :solid : :dot,
                  label = "", alpha = main ? 1.0 : 0.7)
            ## Dash-dot, median curve only: the SAME runs scored at 0.02R instead. The
            ## gap between the two lines is what the contact convention alone is worth,
            ## with the physics held fixed.
        end

        sml_style = bo_oh_style(d.R[sml], GRAD_SML)
        big_style = bo_oh_style(d.R[big], GRAD_BIG; ms_lo = 4.5, ms_hi = 7.5)

        scatter!(plt, d.We[sml], d.C[sml]; mc = sml_style.mc, ms = sml_style.ms,
                 msc = :gray30, msw = 0.5,
                 label = @sprintf("R = %.3f–%.3f mm  (n = %d)",
                                  1000minimum(d.R[sml]), 1000maximum(d.R[sml]), sum(sml)))
        scatter!(plt, d.We[big], d.C[big]; mc = big_style.mc, ms = big_style.ms,
                 msc = :gray30, msw = 0.3, alpha = 0.85,
                 label = @sprintf("R = %.2f–%.2f mm  (n = %d)",
                                  1000minimum(d.R[big]), 1000maximum(d.R[big]), sum(big)))
        plot!(plt, Float64[], Float64[]; c = :gray30, lw = 3, label = "model, median R")
        plot!(plt, Float64[], Float64[]; c = :gray30, lw = 1.5, ls = :dot,
              label = "model, smallest / largest R")
        annotate!(plt, 10.0^(log10(lo) + 0.3), 0.12,
                  text("marker size ∝ √Bo (∝ R), shade ∝ Oh -- both within their own population",
                       8, :left, :gray45))

        out = joinpath(FIGS, "figure_thenarianto_model_" * s.key *
                             (PREVIEW ? "_preview" : "") * ".png")
        savefig(plt, out)
        println("wrote ", out)
    end

    println()
    for (m, c) in zip(meta, curves)
        @printf("%-6s %-6s R = %.3f mm  Oh = %.4f  Bo = %.2e   %s\n",
                m.key, m.big ? "large" : "small", 1000c.R, c.Oh, c.Bo,
                curve_report(c.xs, c.ys, TH))
        for (a, b, dy) in jumps(c.xs, c.ys, TH)
            @printf("         step of %.3f in restitution between We = %.4g and %.4g\n", dy, a, b)
        end
    end
end

main()
