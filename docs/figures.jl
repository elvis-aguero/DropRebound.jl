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
#    where the oscillatory branch becomes aperiodic. This is the picture that
#    Section 8 of the Reid derivation describes in words.
# ---------------------------------------------------------------------------
function root_locus()
    l = 2
    s0 = sqrt(l * (l - 1) * (l + 2))          # inviscid frequency
    ohs = exp10.(range(-2, 0.45; length=400))
    re1, im1, re2, im2 = Float64[], Float64[], Float64[], Float64[]
    for Oh in ohs
        q1 = dominant_root(Oh, l)
        q2 = second_root(Oh, l, q1)
        s1 = q1^2 * Oh / s0
        s2 = q2^2 * Oh / s0
        push!(re1, real(s1)); push!(im1, imag(s1))
        push!(re2, real(s2)); push!(im2, imag(s2))
    end
    icrit = findfirst(i -> abs(im1[i]) < 1e-9, eachindex(im1))
    ohcrit = ohs[icrit]

    plt = plot(xlabel="decay rate   Re(σ) / σ₀",
               ylabel="oscillation frequency   Im(σ) / σ₀",
               title="How a viscous drop stops ringing   (mode l = 2)",
               legend=:topright, xlims=(0, 2.6), ylims=(-1.25, 1.25))
    vspan!(plt, [0, real(ohcrit) * 0 + re1[icrit]], c=:steelblue, alpha=0.06, label="")
    hline!(plt, [0.0], c=:gray, lw=1, ls=:dot, label="")

    plot!(plt, re1, im1, lw=2.6, c=:steelblue, label="")
    plot!(plt, re2, im2, lw=2.6, c=:steelblue, label="conjugate pair of roots")

    # ticks along the locus, on the upper branch only, so labels do not collide
    for Oh in (0.05, 0.2, 0.5)
        i = argmin(abs.(ohs .- Oh))
        scatter!(plt, [re2[i]], [im2[i]], m=:circle, ms=4, c=:steelblue, msw=0, label="")
        annotate!(plt, re2[i] + 0.07, im2[i] + 0.09,
                  text("Oh = $(Oh)", 7, :steelblue, :left))
    end

    scatter!(plt, [re1[icrit]], [0.0], m=:circle, ms=7, c=:black,
             label=@sprintf("critical point, Oh = %.3f", ohcrit))

    annotate!(plt, 0.30,  0.55, text("oscillatory:\nconjugate pair", 8, :steelblue, :center))
    annotate!(plt, 1.75,  0.62, text("aperiodic: two real roots", 8, :indianred, :center))
    annotate!(plt, 1.95,  0.14, text("fast-decaying root →", 7, :indianred, :center))
    annotate!(plt, 0.42, -0.16, text("← slow creep root", 7, :indianred, :center))
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
