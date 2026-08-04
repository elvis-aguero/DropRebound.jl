# CoR and contact time against Weber number for one fluid treated three ways:
# Newtonian at its zero-shear viscosity, Newtonian at its infinite-shear plateau, and
# the actual Carreau-Yasuda law in between.
#
# The two Newtonian curves are not alternatives to the third, they BRACKET it, and the
# bracketing is a theorem rather than an expectation: since
# eta_inf <= eta(gammadot) <= eta_0 pointwise, the dissipation form is ordered the same
# way, so the shear-thinning drop dissipates less than the eta_0 drop and more than the
# eta_inf one. Its restitution must lie between the two curves. A crossing would mean
# the eta-weighted quadrature has a sign or factor error -- something no comparison
# against a single measured number would reveal.
#
# The rheology here is ILLUSTRATIVE, chosen so both bounds actually rebound. At the
# 3000 ppm experimental fluid's own values the eta_0 bound is degenerate: Oh_0 = 57
# never releases, so it contributes no curve at all and the bracket is vacuous on that
# side. Oh_0 = 0.3 with a plateau ratio of 0.1 keeps both bounds live AND inside the
# Ohnesorge range the experiments cover -- which matters for contact time: the measured
# tc is flat in We above We ~ 0.1, and how steeply the model's tc falls depends on Oh.
# At Oh = 0.5, outside that range, it falls about half again too fast.
using Printf, Plots
using DropSolver
gr()

const BO = 0.0189
const OH_0, RATIO = 0.3, 0.1
const OH_INF = OH_0 * RATIO
const LAM, A_CY, N_CY = 10.0, 2.0, 0.5
const M_RUN, K_RUN = 14, 2

@printf("Oh_0 = %.3g   Oh_inf = %.3g   lambda_c = %.3g   a = %.1f   n = %.1f\n",
        OH_0, OH_INF, LAM, A_CY, N_CY)
etaf = gd -> carreau(gd; lambda_c = LAM, a = A_CY, n = N_CY, eta_inf_ratio = RATIO)
@printf("eta/eta_0 at gammadot = 1 is %.3f, so Oh_eff there is about %.3f\n\n",
        etaf(1.0), OH_0 * etaf(1.0))

grid = [0.02 * (3.0/0.02)^(k/8) for k in 0:10]   # the original 9 points, plus two more
run1(We, Oh; eta = nothing) = eta === nothing ?
    simulate(ImpactParams(We = We, Bo = BO, Oh = Oh, M = M_RUN, K = K_RUN, t_max = 25.0)) :
    simulate(ImpactParams(We = We, Bo = BO, Oh = Oh, M = M_RUN, K = K_RUN,
                          eta = eta, t_max = 25.0))

c0 = Float64[]; t0v = Float64[]; ci = Float64[]; tiv = Float64[]
cs = Float64[]; tsv = Float64[]
z0 = Float64[]; zi = Float64[]; zs = Float64[]
amp(p, r) = maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a)
@printf("%-8s | %-22s | %-22s | %s\n", "We", "Newtonian eta_0", "Newtonian eta_inf",
        "Carreau-Yasuda")
for We in grid
    a = run1(We, OH_0); b = run1(We, OH_INF); c = run1(We, OH_0; eta = etaf)
    pa = ImpactParams(We=We, Bo=BO, Oh=OH_0, M=M_RUN, K=K_RUN)
    push!(z0, amp(pa, a)); push!(zi, amp(pa, b)); push!(zs, amp(pa, c))
    push!(c0, a.cor); push!(t0v, a.tc)
    push!(ci, b.cor); push!(tiv, b.tc)
    push!(cs, c.cor); push!(tsv, c.tc)
    @printf("%-8.4g | CoR %-7s tc %5.3f |z|%.2f | CoR %-7s tc %5.3f |z|%.2f | CoR %-7s tc %5.3f |z|%.2f\n",
            We, isfinite(a.cor) ? @sprintf("%.4f",a.cor) : "NONE", a.tc, z0[end],
                isfinite(b.cor) ? @sprintf("%.4f",b.cor) : "NONE", b.tc, zi[end],
                isfinite(c.cor) ? @sprintf("%.4f",c.cor) : "NONE", c.tc, zs[end])
end

# --- the bracketing, as a check rather than a claim -------------------------------
lo = count(i -> cs[i] < c0[i] - 1e-6, eachindex(cs))     # below the eta_0 bound
hi = count(i -> cs[i] > ci[i] + 1e-6, eachindex(cs))     # above the eta_inf bound
@printf("\nbracketing violations: %d below the eta_0 curve, %d above the eta_inf curve",
        lo, hi)
@printf("  (of %d points)\n", length(cs))
@printf("shear-thinning CoR spans %.3f..%.3f, between %.3f..%.3f and %.3f..%.3f\n",
        minimum(cs), maximum(cs), minimum(c0), maximum(c0), minimum(ci), maximum(ci))

## VALIDITY MASK. The shape expansion is linear in the amplitude, so a point with
## |zeta| approaching one is outside the model regardless of whether the solve
## converged -- and a run that never released has no restitution at all. Both are
## excluded rather than drawn, because a plotted point implies a claim.
ok(c, z) = [isfinite(c[i]) && z[i] < 0.5 ? c[i] : NaN for i in eachindex(c)]
c0p, cip, csp = ok(c0, z0), ok(ci, zi), ok(cs, zs)
t0p, tip, tsp = ok(t0v, z0), ok(tiv, zi), ok(tsv, zs)
for (nm, cc) in (("eta_0", c0p), ("eta_inf", cip), ("Carreau", csp))
    n_ok = count(!isnan, cc)
    @printf("%-8s usable at %d of %d We points; drops out above We = %.2f\n", nm, n_ok,
            length(cc), n_ok == length(cc) ? Inf : grid[n_ok])
end

common = (xscale = :log10, xlabel = "We", framestyle = :box, grid = true,
          gridalpha = 0.15, tickfontsize = 8, guidefontsize = 9, legendfontsize = 7,
          titlefontsize = 10)

p1 = plot(grid, cip; label = "η∞  (Oh = $OH_INF)", lc = :seagreen, lw = 2,
          ls = :dash, marker = :utriangle, ms = 3.5, mc = :seagreen,
          ylabel = "coefficient of restitution", title = "Restitution",
          ylims = (0.0, 1.0), legend = :bottomleft, common...)
plot!(p1, grid, csp; label = "Carreau–Yasuda", lc = :crimson, lw = 2.5,
      marker = :circle, ms = 4.5, mc = :crimson)
plot!(p1, grid, c0p; label = "η₀  (Oh = $OH_0)", lc = :darkorange, lw = 2,
      ls = :dashdot, marker = :dtriangle, ms = 3.5, mc = :darkorange)

p2 = plot(grid, tip; label = "η∞", lc = :seagreen, lw = 2, ls = :dash,
          marker = :utriangle, ms = 3.5, mc = :seagreen,
          ylabel = "contact time  t_c / √(ρR³/σ)", title = "Contact time",
          legend = :topright, common...)
plot!(p2, grid, tsp; label = "Carreau–Yasuda", lc = :crimson, lw = 2.5,
      marker = :circle, ms = 4.5, mc = :crimson)
plot!(p2, grid, t0p; label = "η₀", lc = :darkorange, lw = 2, ls = :dashdot,
      marker = :dtriangle, ms = 3.5, mc = :darkorange)

fig = plot(p1, p2; layout = (1, 2), size = (980, 420), dpi = 200,
           left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, top_margin = 3Plots.mm,
           plot_title = "Oh₀ = $OH_0,  Oh∞ = $OH_INF,  λ_c = $LAM",
           plot_titlefontsize = 9)
out = joinpath(@__DIR__, "..", "..", "docs", "figures", "three_rheologies")
savefig(fig, out * ".png"); savefig(fig, out * ".svg")
println("wrote $(out).png")
