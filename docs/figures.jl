# Generates the figures used by the documentation into `docs/src/assets/`.
#
# Figures are produced here rather than inside the derivation pages so that the
# plotting code does not appear in the rendered text. The pages reference the
# results with ordinary Markdown image syntax.

using Plots, Printf
using DropSolver

gr()
default(fontfamily="sans-serif", framestyle=:box, grid=true, gridalpha=0.18,
        legendfontsize=8, guidefontsize=9, titlefontsize=10, tickfontsize=8,
        size=(760, 430), dpi=140)

const ASSETS = joinpath(@__DIR__, "src", "assets")
mkpath(ASSETS)
# PNG rather than SVG: GR's SVG backend writes content outside its declared
# viewBox on multi-block plots, which clips the right-hand side.
save(name, plt) = savefig(plt, joinpath(ASSETS, name))

# ---------------------------------------------------------------------------
# 1. Root locus: how the two dominant roots move as viscosity increases, and
#    where the oscillatory branch becomes aperiodic, beside what that means for
#    the surface. Referenced from *The Free Viscous Drop*.
# ---------------------------------------------------------------------------
function root_locus()
    l  = 2
    s0 = sqrt(l * (l - 1) * (l + 2))
    ohs = exp10.(range(-2, 0.45; length = 400))
    re1, im1, re2, im2 = Float64[], Float64[], Float64[], Float64[]
    for Oh in ohs
        q1 = dominant_root(Oh, l); q2 = second_root(Oh, l, q1)
        s1 = q1^2 * Oh / s0;       s2 = q2^2 * Oh / s0
        push!(re1, real(s1)); push!(im1, imag(s1))
        push!(re2, real(s2)); push!(im2, imag(s2))
    end
    icrit  = findfirst(i -> abs(im1[i]) < 1e-9, eachindex(im1))
    ohcrit = ohs[icrit]

    # LEFT: where the roots live.
    p1 = plot(xlabel = "decay rate  Re σ / σ₀", ylabel = "frequency  Im σ / σ₀",
              title = "Roots move as viscosity rises", legend = :topright,
              xlims = (0, 2.6), ylims = (-1.3, 1.3),
              guidefontsize = 11, titlefontsize = 12, tickfontsize = 10,
              legendfontsize = 10)
    hline!(p1, [0.0], c = :gray70, lw = 1, ls = :dot, label = "")
    osc = 1:icrit
    plot!(p1, re2[osc], im2[osc], lw = 3, c = :steelblue, label = "oscillatory")
    plot!(p1, re1[osc], im1[osc], lw = 3, c = :steelblue, label = "")
    plot!(p1, re1[icrit:end], im1[icrit:end], lw = 3, c = :indianred, label = "overdamped")
    plot!(p1, re2[icrit:end], im2[icrit:end], lw = 3, c = :indianred, label = "")
    scatter!(p1, [re1[icrit]], [0.0], m = :circle, ms = 8, c = :black, msw = 0,
             label = @sprintf("critical Oh = %.2f", ohcrit))
    # Past the critical point the pair splits into two real roots that move in
    # opposite directions. Without saying so, the red branch reads as one line
    # drawn through the critical point rather than as two.
    annotate!(p1, 2.05, 0.20, text("fast root", 10, :indianred, :center))
    annotate!(p1, 0.42, 0.20, text("slow creep", 10, :indianred, :center))

    # RIGHT: what that means for the surface, which is the point of the left panel.
    ts = range(0, 12; length = 600)
    p2 = plot(xlabel = "time  (capillary times)", ylabel = "surface amplitude  ζ₂",
              title = "…and what the drop does", legend = :topright,
              guidefontsize = 11, titlefontsize = 12, tickfontsize = 10,
              legendfontsize = 10)
    for (Oh, c, lab) in ((0.05, :steelblue, "Oh = 0.05  (rings)"),
                         (1.20, :indianred, "Oh = 1.20  (creeps back)"))
        q = dominant_root(Oh, l); s = q^2 * Oh / s0
        y = real.(exp.(-s .* ts))
        plot!(p2, ts, y ./ maximum(abs, y), lw = 3, c = c, label = lab)
    end
    hline!(p2, [0.0], c = :gray70, lw = 1, ls = :dot, label = "")

    plt = plot(p1, p2; layout = (1, 2), size = (1040, 420), dpi = 150,
               left_margin = 6Plots.mm, bottom_margin = 6Plots.mm)
    save("reid_root_locus.png", plt)
end

# ---------------------------------------------------------------------------
# 2. Where the classical small-viscosity formula stops being usable.
# ---------------------------------------------------------------------------
function reid_vs_lamb()
    ohs = exp10.(range(-3, 1; length=140))
    plt = plot(xscale=:log10, yscale=:log10, xlabel="Ohnesorge number  Oh",
               ylabel="damping rate  λₗ", legend=:topleft,
               title="Exact damping vs the small-viscosity approximation")
    for (l, c) in zip((2, 4, 8), (:steelblue, :seagreen, :indianred))
        exact = [reid_lambda_omega2(Oh, l)[1] for Oh in ohs]
        lamb  = [Oh * (l - 1) * (2l + 1) for Oh in ohs]
        plot!(plt, ohs, lamb,  lw=1.6, ls=:dash, c=c, label="Lamb, l = $l")
        plot!(plt, ohs, exact, lw=2.4, c=c, label="Reid (exact), l = $l")
    end
    save("reid_vs_lamb.png", plt)
end

# ---------------------------------------------------------------------------
# 3. The measured flow curve of the validation fluid, and the part of it the
#    drop actually samples during an impact.
# ---------------------------------------------------------------------------
function viscosity_curve()
    ETA_INF = 0.0037320997942061666; ETA_0 = 8.433817577956766
    K = 18.48081673111359; m = 0.7430524574330837
    ratio = ETA_INF / ETA_0
    gd = exp10.(range(-4, 6; length=400))
    eta = @. ratio + (1 - ratio) / (1 + (K * gd)^m)

    plt = plot(gd, eta, xscale=:log10, yscale=:log10, lw=2.6, c=:steelblue,
               xlabel="shear rate  γ̇   (s⁻¹)", ylabel="η(γ̇) / η₀",
               title="Flow curve of the validation fluid", label="Cross model")
    hline!(plt, [1.0], c=:gray, ls=:dash, lw=1.2, label="rest viscosity η₀")
    hline!(plt, [ratio], c=:gray, ls=:dot, lw=1.2, label="infinite-shear floor η∞")
    vspan!(plt, [1e1, 1e4], c=:orange, alpha=0.13, label="range sampled during impact")
    save("cy_flow_curve.png", plt)
end

# ---------------------------------------------------------------------------
# 4. What each modelling choice does to the mode-coupling matrix.
# ---------------------------------------------------------------------------
function coupling_structure()
    M, gap = 26, 7
    W = 3M + 2gap
    A = fill(NaN, M, W)
    for i in 1:M, j in 1:M
        ## Parity: G^{l'}_{l l''} vanishes unless l + l' + l'' is even, so for any
        ## single viscosity harmonic roughly half the entries are zero even in the
        ## fully coupled case. Drawing the right panel solid would contradict the
        ## selection rule the figure illustrates.
        even = iseven(i + j)
        A[i, j]             = (i == j)                  ? 1.0 : 0.0
        A[i, j + M + gap]   = (abs(i - j) <= 4 && even) ? 1.0 : 0.0
        A[i, j + 2M + 2gap] = even                      ? 1.0 : 0.0
    end
    plt = heatmap(A, c=cgrad([:white, :steelblue]), cbar=false,
                  ticks=false, framestyle=:none, grid=false, size=(940, 360),
                  ylims=(-4.5, M + 0.5))
    for (k, ttl) in enumerate(("η constant\nmodes stay independent",
                               "η varies slowly\nnearby modes couple",
                               "η varies sharply\nevery parity-allowed pair couples"))
        cx = (k - 1) * (M + gap) + M / 2 + 0.5
        annotate!(plt, cx, -2.3, text(ttl, 10, :black, :center))
    end
    save("coupling_structure.png", plt)
end

function build_all()
    root_locus(); reid_vs_lamb(); viscosity_curve(); coupling_structure()
    @info "figures written" dir=ASSETS
end
