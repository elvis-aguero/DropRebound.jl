# The 3000 ppm shear-thinning experiments, overlaid with three solver backends.
#
# One fluid, nothing fitted: Oh, Bo and the Carreau-Yasuda parameters all come from the fluid's
# own Cross characterisation. The three curves differ ONLY in how the drop is solved, so the
# spread between them is the numerical method's contribution and the distance from the points
# is the model's.
#
#   nonvar/tangency   surface amplitudes, Reid coefficients, contact ranked by tangency error
#   var/lcp           interior in the state, contact by complementarity, spectral film pressure
#   var/lcp/nodal     the same, with the conjugate nodal-load forcing
#
# Usage: julia --project=docs julia/scripts/figure_backends_shear_thinning.jl

using Printf, Plots
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))
gr()

const DATA = joinpath(@__DIR__, "..", "derivations", "data", "metrics_3000ppm.csv")

eta_0, eta_inf   = 8.433817577956766, 0.0037320997942061666
K_cross, m_cross = 18.48081673111359, 0.7430524574330837
R, SIGMA, G, BO  = 0.0003, 0.0728, 9.81, 0.012
RHO   = BO * SIGMA / (G * R^2)
T_CAP = sqrt(RHO * R^3 / SIGMA)
OH_0  = eta_0 / sqrt(RHO * SIGMA * R)
LAM, A_CY, N_CY = K_cross / T_CAP, m_cross, 1 - m_cross
ETA_RATIO = eta_inf / eta_0
EPS_ST    = (1 - N_CY) / A_CY

## the export is semicolon separated with comma decimals, and tc is in seconds
rows = NTuple{3,Float64}[]
for (i, ln) in enumerate(eachline(DATA))
    i == 1 && continue
    f = split(strip(ln), ';'); length(f) < 3 && continue
    v = tryparse.(Float64, replace.(f[1:3], ',' => '.'))
    any(isnothing, v) && continue
    push!(rows, (v[1], v[2], v[3] / T_CAP))
end
We_e = [r[1] for r in rows]; cor_e = [r[2] for r in rows]; tc_e = [r[3] for r in rows]
@printf("%d experiments, Oh_0 = %.1f, Bo = %.3f\n", length(rows), OH_0, BO)

eta_st = gd -> carreau(gd; lambda_c = LAM, a = A_CY, n = N_CY, eta_inf_ratio = ETA_RATIO)

## RESOLUTION. The variational curves run at M = 30, K = 3, which is where the shear-thinning
## restitution converges: K = 2 is 0.8 per cent low and K = 4 and 5 agree with K = 3 to four
## decimals. Earlier versions of this figure used M = 14, K = 2 because a run cost minutes;
## caching the assembly geometry made it seconds, so there is no longer a reason to publish the
## coarser number. The nonvariational curve keeps M = 20 -- its closure is diagonal and it has
## no K at all.
const CURVES = [
    (b = Backend(formulation = :nonvariational, contact = :tangency),
     M = 20, K = 2, colour = :darkorange, marker = :diamond, dash = :dash),
    (b = Backend(),                 M = 30, K = 3, colour = :crimson,  marker = :circle,    dash = :solid),
    (b = Backend(forcing = :nodal), M = 30, K = 3, colour = :seagreen, marker = :utriangle, dash = :dot),
]

## The grid spans the 5th to 95th percentile of the measured Weber numbers, NOT their raw
## extremes. The lowest measured point is We = 0.0045, and there the nonvariational solver
## drives its step size into the floor and grinds for an hour without finishing -- a gentle
## impact means a long contact, and its adaptive stepping cannot resolve the release. Plotting
## a curve over a range one of the backends cannot run is not a comparison.
grid = exp.(range(log(quantile(We_e, 0.05)), log(quantile(We_e, 0.95)); length = 8))
@printf("Weber grid: %.4g to %.4g, %d points\n", first(grid), last(grid), length(grid))
results = Dict{String,Any}()
for c in CURVES
    nm = label(c.b); cor_m = Float64[]; tc_m = Float64[]; wes = Float64[]; wall = 0.0
    @printf("\n%s (M = %d, K = %d)\n", nm, c.M, c.K); flush(stdout)
    for We in grid
        r = run_impact(c.b; We = We, Bo = BO, Oh = OH_0, M = c.M, K = c.K, t_max = 25.0,
                       eta = c.b.formulation === :variational ? eta_st : nothing,
                       eta_nonvar = c.b.formulation === :variational ? nothing :
                           STExactParams(c.M, OH_0, LAM, A_CY, EPS_ST;
                                         viscous = :reid, eta_inf_ratio = ETA_RATIO))
        wall += r.wall
        if r.ok
            push!(wes, We); push!(cor_m, r.cor); push!(tc_m, r.tc)
            @printf("  We=%-9.4g CoR=%.4f  tc=%.3f  (%.0fs)\n", We, r.cor, r.tc, r.wall); flush(stdout)
        else
            @printf("  We=%-9.4g did not complete            (%.0fs)\n", We, r.wall); flush(stdout)
        end
    end
    results[nm] = (we = wes, cor = cor_m, tc = tc_m, wall = wall, c = c)
    ## cached, so that changing how this is DRAWN never costs another sweep
    open(joinpath(@__DIR__, "..", "..", "results", "backends_shear_thinning.csv"),
         nm == label(CURVES[1].b) ? "w" : "a") do io
        nm == label(CURVES[1].b) && println(io, "backend,We,cor,tc")
        for k in eachindex(wes)
            @printf(io, "%s,%.10g,%.10g,%.10g\n", nm, wes[k], cor_m[k], tc_m[k])
        end
    end
    @printf("  -> %d of %d completed, %.0f s total\n", length(wes), length(grid), wall)
end

common = (xscale = :log10, xlabel = "Weber number  We", framestyle = :box,
          grid = true, gridalpha = 0.15, tickfontsize = 8, guidefontsize = 9,
          legendfontsize = 7, titlefontsize = 10)

p1 = scatter(We_e, cor_e; label = "experiment (n=$(length(rows)))", mc = :steelblue,
             ms = 3.0, msw = 0.3, ma = 0.55, ylabel = "coefficient of restitution",
             title = "Restitution", legend = :bottomleft, common...)
## The contact-time axis is bounded to the measured range. At the lowest Weber numbers the
## nonvariational backend does not release at all and reports 20 to 25 capillary times against
## a measured 4; drawn to scale that one curve flattens everything else into a line. Those
## points are clipped and called out in the caption rather than quietly dropped.
p2 = scatter(We_e, tc_e; label = "experiment", mc = :steelblue, ms = 3.0, msw = 0.3,
             ma = 0.55, ylabel = "contact time  /  (ρR³/σ)^{1/2}",
             ylim = (0.0, 1.35*maximum(tc_e)),
             title = "Contact time", legend = :topright, common...)
for c in CURVES
    r = results[label(c.b)]
    isempty(r.we) && continue
    plot!(p1, r.we, r.cor; label = label(c.b), lc = c.colour, lw = 2, ls = c.dash,
          marker = c.marker, ms = 3.5, mc = c.colour, msw = 0)
    plot!(p2, r.we, r.tc;  label = label(c.b), lc = c.colour, lw = 2, ls = c.dash,
          marker = c.marker, ms = 3.5, mc = c.colour, msw = 0)
end

fig = plot(p1, p2; layout = (1, 2), size = (1000, 420), dpi = 200,
           left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)
out = joinpath(@__DIR__, "..", "..", "docs", "figures", "backends_shear_thinning")
savefig(fig, out * ".png"); savefig(fig, out * ".svg")
@printf("\nwrote %s.png\n", out)

println("\nwallclock over the ten-point sweep:")
for c in CURVES
    r = results[label(c.b)]
    @printf("  %-18s %6.0f s   (%d of %d points)\n", label(c.b), r.wall, length(r.we), length(grid))
end
