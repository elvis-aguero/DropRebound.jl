# Overlay of the 3000 ppm shear-thinning experiments and the variational solver, for
# both KPIs. Ten simulations spanning the measured Weber range, one fluid, nothing
# fitted: Oh, Bo and the Carreau-Yasuda parameters all come from the fluid's own
# Cross characterisation.
using Printf, Plots
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))
gr()

const DATA = joinpath(@__DIR__, "..", "derivations", "data", "metrics_3000ppm.csv")

eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
K_cross, m_cross = 18.48081673111359, 0.7430524574330837
R, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
T_CAP = sqrt(RHO * R^3 / SIGMA)
OH_0  = eta_0 / sqrt(RHO * SIGMA * R)
LAM, A_CY, N_CY = K_cross / T_CAP, m_cross, 1 - m_cross
ETA_RATIO = eta_inf / eta_0

## the export is semicolon separated with comma decimals, and tc is in seconds
rows = NTuple{3,Float64}[]
for (i, ln) in enumerate(eachline(DATA))
    i == 1 && continue
    f = split(strip(ln), ';'); length(f) < 3 && continue
    v = tryparse.(Float64, replace.(f[1:3], ',' => '.'))
    any(isnothing, v) && continue
    push!(rows, (v[1], v[2], v[3] / T_CAP))
end
We_e  = [r[1] for r in rows]; cor_e = [r[2] for r in rows]; tc_e = [r[3] for r in rows]
@printf("%d experiments, Oh_0 = %.1f\n", length(rows), OH_0)

eta_st = gd -> carreau(gd; lambda_c = LAM, a = A_CY, n = N_CY, eta_inf_ratio = ETA_RATIO)
grid = exp.(range(log(minimum(We_e)), log(maximum(We_e)); length = 10))
cor_m = Float64[]; tc_m = Float64[]
for We in grid
    r = simulate(ImpactParams(We = We, Bo = BO, Oh = OH_0, M = 14, K = 2,
                              eta = eta_st, t_max = 25.0))
    push!(cor_m, r.cor); push!(tc_m, r.tc)
    @printf("  We=%-8.4g CoR=%.4f tc=%.3f\n", We, r.cor, r.tc)
end

## The Newtonian drop at the same zero-shear viscosity, for contrast. It does not
## rebound at all, so it has no CoR to plot -- which is the whole point.
rn = simulate(ImpactParams(We = 0.19, Bo = BO, Oh = OH_0, M = 14, K = 2, t_max = 25.0))
@printf("Newtonian at Oh_0 = %.1f : CoR = %s (never releases)\n", OH_0, string(rn.cor))

common = (xscale = :log10, xlabel = "Weber number  We",
          framestyle = :box, grid = true, gridalpha = 0.15, tickfontsize = 8,
          guidefontsize = 9, legendfontsize = 7, titlefontsize = 10)

p1 = scatter(We_e, cor_e; label = "experiment (n=$(length(rows)))",
             mc = :steelblue, ms = 3.5, msw = 0, ma = 0.55,
             ylabel = "coefficient of restitution", title = "Restitution",
             ylims = (0.0, 1.0), legend = :bottomleft, common...)
plot!(p1, grid, cor_m; label = "model, shear-thinning", lc = :crimson, lw = 2,
      marker = :circle, ms = 4, mc = :crimson)
## Stated rather than plotted: a Newtonian drop at this viscosity never releases, so
## it contributes no point, and a legend entry for an empty series is noise.
annotate!(p1, 0.055, 0.30,
          text("Newtonian at the same Oh₀:\nno rebound at any We", 7, :center,
               RGB(0.4, 0.4, 0.4)))

p2 = scatter(We_e, tc_e; label = "experiment", mc = :steelblue, ms = 3.5,
             msw = 0, ma = 0.55, ylabel = "contact time  t_c / √(ρR³/σ)",
             title = "Contact time", legend = :topright, common...)
plot!(p2, grid, tc_m; label = "model, shear-thinning", lc = :crimson, lw = 2,
      marker = :circle, ms = 4, mc = :crimson)

fig = plot(p1, p2; layout = (1, 2), size = (980, 420), dpi = 200,
           left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, top_margin = 3Plots.mm,
           plot_title = "3000 ppm drop,  Oh₀ = $(round(OH_0, digits=1))",
           plot_titlefontsize = 10)
out = joinpath(@__DIR__, "..", "..", "docs", "figures", "shear_thinning_overlay")
savefig(fig, out * ".png"); savefig(fig, out * ".svg")
println("wrote $(out).png and .svg")
