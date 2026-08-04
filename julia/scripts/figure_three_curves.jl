# CoR against Weber number for the SAME fluid treated three ways: Newtonian at its
# zero-shear viscosity, Newtonian at its infinite-shear plateau, and the actual
# Carreau-Yasuda law.
#
# The two Newtonian curves are not alternatives to the third -- they BRACKET it, and
# the bracketing is a theorem rather than an expectation. Since
# eta_inf <= eta(gammadot) <= eta_0 pointwise, the dissipation form is ordered the same
# way, so the shear-thinning drop must dissipate less than the eta_0 drop and more than
# the eta_inf one. Its restitution must therefore lie between the two curves. If it does
# not, the eta-weighted quadrature has a sign or factor error that no comparison with a
# single measured number would reveal.
using Printf, Statistics, Plots
using DropSolver
gr()

const DATA = joinpath(@__DIR__, "..", "derivations", "data", "metrics_3000ppm.csv")
eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
K_cross, m_cross = 18.48081673111359, 0.7430524574330837
R, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
T_CAP = sqrt(RHO * R^3 / SIGMA)
OH_0   = eta_0   / sqrt(RHO * SIGMA * R)
OH_INF = eta_inf / sqrt(RHO * SIGMA * R)
LAM, A_CY, N_CY = K_cross / T_CAP, m_cross, 1 - m_cross
RATIO = eta_inf / eta_0
@printf("Oh_0 = %.3g   Oh_inf = %.4g   ratio = %.3g   lambda_c = %.3g\n",
        OH_0, OH_INF, RATIO, LAM)

rows = NTuple{3,Float64}[]
for (i, ln) in enumerate(eachline(DATA))
    i == 1 && continue
    f = split(strip(ln), ';'); length(f) < 3 && continue
    v = tryparse.(Float64, replace.(f[1:3], ',' => '.'))
    any(isnothing, v) && continue
    push!(rows, (v[1], v[2], v[3] / T_CAP))
end
We_e = [r[1] for r in rows]; cor_e = [r[2] for r in rows]; tc_e = [r[3] for r in rows]

etaf = gd -> carreau(gd; lambda_c = LAM, a = A_CY, n = N_CY, eta_inf_ratio = RATIO)
grid = exp.(range(log(minimum(We_e)), log(maximum(We_e)); length = 10))

const M_RUN, K_RUN = 14, 2
run1(We, Oh; eta = nothing) = simulate(
    eta === nothing ? ImpactParams(We = We, Bo = BO, Oh = Oh, M = M_RUN, K = K_RUN,
                                   t_max = 25.0)
                    : ImpactParams(We = We, Bo = BO, Oh = Oh, M = M_RUN, K = K_RUN,
                                   eta = eta, t_max = 25.0))

c0 = Float64[]; t0v = Float64[]
ci = Float64[]; tiv = Float64[]
cs = Float64[]; tsv = Float64[]
@printf("\n%-9s | %-16s | %-16s | %s\n", "We", "Newt eta_0", "Newt eta_inf", "Carreau-Yasuda")
for We in grid
    a = run1(We, OH_0);   b = run1(We, OH_INF);   c = run1(We, OH_0; eta = etaf)
    push!(c0, a.cor); push!(t0v, a.tc)
    push!(ci, b.cor); push!(tiv, b.tc)
    push!(cs, c.cor); push!(tsv, c.tc)
    @printf("%-9.4g | CoR %-7s tc %5.2f | CoR %.4f tc %5.2f | CoR %.4f tc %5.2f\n",
            We, isfinite(a.cor) ? @sprintf("%.4f", a.cor) : "none", a.tc,
            b.cor, b.tc, c.cor, c.tc)
end

## the bracketing, checked
viol = count(i -> isfinite(cs[i]) && isfinite(ci[i]) && cs[i] > ci[i] + 1e-6, eachindex(cs))
@printf("\nbracketing: shear-thinning above the eta_inf bound at %d of %d points\n",
        viol, length(cs))

common = (xscale = :log10, xlabel = "Weber number  We", framestyle = :box,
          grid = true, gridalpha = 0.15, tickfontsize = 8, guidefontsize = 9,
          legendfontsize = 7, titlefontsize = 10)

p1 = scatter(We_e, cor_e; label = "experiment (n=$(length(rows)))", mc = :steelblue,
             ms = 3.5, msw = 0, ma = 0.5, ylabel = "coefficient of restitution",
             title = "Restitution", ylims = (0.0, 1.0), legend = :bottomleft, common...)
plot!(p1, grid, ci; label = "Newtonian at η∞  (upper bound)", lc = :seagreen,
      lw = 2, ls = :dash, marker = :utriangle, ms = 3.5, mc = :seagreen)
plot!(p1, grid, cs; label = "Carreau–Yasuda", lc = :crimson, lw = 2.5,
      marker = :circle, ms = 4, mc = :crimson)
annotate!(p1, 0.055, 0.10,
          text("Newtonian at η₀ : no rebound at any We", 7, :center, RGB(0.35,0.35,0.35)))

p2 = scatter(We_e, tc_e; label = "experiment", mc = :steelblue, ms = 3.5, msw = 0,
             ma = 0.5, ylabel = "contact time  t_c / √(ρR³/σ)", title = "Contact time",
             legend = :topright, common...)
plot!(p2, grid, tiv; label = "Newtonian at η∞", lc = :seagreen, lw = 2, ls = :dash,
      marker = :utriangle, ms = 3.5, mc = :seagreen)
plot!(p2, grid, tsv; label = "Carreau–Yasuda", lc = :crimson, lw = 2.5,
      marker = :circle, ms = 4, mc = :crimson)

fig = plot(p1, p2; layout = (1, 2), size = (980, 420), dpi = 200,
           left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, top_margin = 3Plots.mm,
           plot_title = "3000 ppm fluid, three rheologies  (Oh₀ = $(round(OH_0,digits=1)), Oh∞ = $(round(OH_INF,digits=4)))",
           plot_titlefontsize = 9)
out = joinpath(@__DIR__, "..", "..", "three_rheologies")
savefig(fig, out * ".png"); savefig(fig, out * ".svg")
println("wrote $(out).png")
