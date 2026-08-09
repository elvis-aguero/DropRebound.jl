# Vorticity through a bounce: where rotation is created, and when.
#
# The drop arrives irrotational. Nothing in a free capillary oscillation makes
# vorticity if it starts without any, so every bit of rotation on this figure was
# generated during the impact and then diffused. Watching it appear is watching the
# contact do work on the interior, which is the part of the bounce the restitution
# number summarises away.
#
# WHAT IS PLOTTED
#   (a) a space-time map of |omega| just under the free surface, so the vertical axis
#       is position on the drop and the horizontal axis is time
#   (b) enstrophy, the volume integral of omega^2, which is the single number saying
#       how much rotation the drop is carrying at that instant
#
# The drop is against the substrate for the whole record: it starts in contact and the
# march stops when it releases, so the time axis IS the contact interval. Release is
# marked; there is no free flight to compare against on this plot.

using Printf, LinearAlgebra
using Plots
using DropSolver

gr()

const ASSET = joinpath(@__DIR__, "..", "..", "docs", "src", "assets")
const OUT   = joinpath(@__DIR__, "..", "..", "results")
mkpath(ASSET); mkpath(OUT)

# ---------------------------------------------------------------------------
# THE VORTICITY OF A RITZ MODE
#
# The interior velocity is a Stokes stream function. For one surface mode `l` with
# radial trial function `f(x)`,
#
#     Psi = f(x) C_l(theta),        C_l = -sin(theta) P_l'(theta) / L,   L = l(l+1)
#
# which is what makes the code's field read u_r = f/x^2 P_l and u_th = f'/(x L) P_l'.
# Axisymmetric flow has one vorticity component, omega_phi = -(1/(x sin th)) E^2 Psi,
# with E^2 = d_xx + (sin th / x^2) d_th( (1/sin th) d_th ).
#
# The angular operator acting on C_l returns -L C_l -- that is the Gegenbauer equation,
# and it is why this basis is the natural one -- so E^2 Psi = (f'' - L f/x^2) C_l and
#
#     omega_phi = (1/(x L)) (f'' - L f / x^2) dP_l/dtheta .
#
# The sin(theta) cancels, so this is finite on the axis and can be evaluated at the
# pole directly. Note the bracket: a trial function with f'' = L f/x^2 carries no
# vorticity at all. The first radial function is f = x^(l+1), for which
# f'' = l(l+1) x^(l-1) = L f/x^2 exactly, so K = 1 is irrotational and every bit of
# rotation on this figure comes from k >= 2. That is asserted below rather than
# assumed, because if it failed the figure would be plotting a basis artefact.
#
# The amplitudes used are `adot`, not `a`: the mass matrix is built from the field
# components as a kinetic energy, so it is the RATE that multiplies the field to give
# the velocity. Using `a` would plot the vorticity of the displacement, which is not
# a physical quantity.
# ---------------------------------------------------------------------------
"""
    vorticity(p, adot, x, mu) -> Float64

The azimuthal vorticity of the interior flow at radius `x` and `mu = cos(theta)`.
"""
function vorticity(p::ImpactParams, adot::AbstractVector, x::Real, mu::Real)
    b = DropSolver.basis(p)
    w = 0.0
    for (i, l) in enumerate(b.ls)
        L  = l * (l + 1)
        A  = DropSolver.legendre_angular(l, mu)
        rb = DropSolver.RitzBasis(l, p.K, p.basis_kind)
        s  = 0.0
        for k in 1:p.K
            f   = DropSolver.phi(rb, k, x)
            d2f = DropSolver.d2phi(rb, k, x)
            s  += adot[DropSolver.dofindex(b, i, k)] * (d2f - L * f / x^2)
        end
        w += s * A.dPdth / (x * L)
    end
    w
end

"""Enstrophy, the volume integral of `omega^2`, by Gauss-Legendre in both directions."""
function enstrophy(p::ImpactParams, adot::AbstractVector; nx = 24, nmu = 32)
    xs, wxs   = DropSolver.gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = DropSolver.gauss_legendre_nodes(nmu, -1.0, 1.0)
    s = 0.0
    for (x, wx) in zip(xs, wxs), (mu, wmu) in zip(mus, wmus)
        s += wx * wmu * 2pi * x^2 * vorticity(p, adot, x, mu)^2   # dV = 2 pi x^2 dx dmu
    end
    s
end

# ---------------------------------------------------------------------------
# A colour palette in the manner of MATLAB's `abyss` and `sky`.
#
# Neither ships with Julia and neither can be copied out of MATLAB here, so these are
# hand-built gradients with the same character: `abyss` runs from near-black navy up
# through blue to a pale sky, `sky` is the light half of that range alone. They are
# sequential, which suits |omega|; the sign is carried by panel (b) instead.
# ---------------------------------------------------------------------------
const ABYSS = cgrad([RGB(0.000, 0.012, 0.086), RGB(0.008, 0.106, 0.278),
                     RGB(0.016, 0.239, 0.478), RGB(0.055, 0.400, 0.667),
                     RGB(0.302, 0.612, 0.816), RGB(0.749, 0.882, 0.945)])
const SKY   = cgrad([RGB(0.055, 0.220, 0.400), RGB(0.145, 0.420, 0.667),
                     RGB(0.353, 0.639, 0.831), RGB(0.643, 0.831, 0.925),
                     RGB(0.878, 0.949, 0.980)])

const PALETTE = ABYSS
const DEEP    = RGB(0.016, 0.239, 0.478)     # a single line colour from the same family

# ---------------------------------------------------------------------------
# Run one impact and evaluate the field on it.
# ---------------------------------------------------------------------------
const WE, BO, OH = 2.0, 0.019, 0.0373
const M_RUN, K_RUN = DropSolver.DEFAULT_M, DropSolver.DEFAULT_K

p = ImpactParams(We = WE, Bo = BO, Oh = OH, M = M_RUN, K = K_RUN, t_max = 25.0)

## SELF-CHECK before anything is plotted: the irrotational mode must carry no
## vorticity. If this fires, the radial derivative or the mode indexing is wrong and
## every number below is decoration.
let p1 = ImpactParams(We = WE, Bo = BO, Oh = OH, M = 8, K = 1, t_max = 1.0)
    v = [vorticity(p1, ones(DropSolver.ndof(DropSolver.basis(p1))), x, mu)
         for x in (0.3, 0.7, 0.95), mu in (-0.9, -0.2, 0.5)]
    @assert maximum(abs, v) < 1e-10 "K = 1 is not irrotational: max |omega| = $(maximum(abs, v))"
    println("self-check: K = 1 carries no vorticity (max |omega| = ",
            @sprintf("%.1e", maximum(abs, v)), ")")
end

println("running We = $WE, Oh = $OH at M = $M_RUN, K = $K_RUN ..."); flush(stdout)
r = DropSolver.simulate_lcp(p)
@printf("  %d frames, restitution %.4f, contact time %.3f\n", length(r.t), r.cor, r.tc)

## The contact interval, first touch to last release, as the solver defines it.
touch   = findfirst(>(0), r.cp)
release = findlast(>(0), r.cp)
t_in, t_out = r.t[touch], r.t[release]

## Subsample frames for the map: a few hundred columns is past what the image can show.
const NCOL = 400
cols = unique(round.(Int, range(1, length(r.t); length = min(NCOL, length(r.t)))))
ts   = r.t[cols]

## Evaluate just under the surface. Exactly at x = 1 the trial functions are all equal
## to one and the map is dominated by that normalisation rather than by the flow.
## The vertical axis is measured FROM THE CONTACT POLE, so 0 is the point that touches
## first and 180 is the crown. That puts the interesting end of the drop at the bottom
## of the map and makes the axis ascending, which is what the heatmap wants.
##
## The two poles themselves are left out. Vorticity vanishes on the axis by symmetry
## (it is proportional to dP/dtheta, which is zero there), so nothing is lost, and the
## angular routine divides by sin(theta) and returns NaN if handed mu = +-1 exactly.
const X_EVAL = 0.97
const NTH    = 220
angs = range(1.0, 179.0; length = NTH)      # degrees from the contact pole
mus  = -cos.(deg2rad.(angs))                # mu = -1 at the contact pole

W = [vorticity(p, r.adot[c], X_EVAL, mu) for mu in mus, c in cols]
@assert all(isfinite, W) "the surface field has non-finite entries"
ens = [enstrophy(p, r.adot[c]) for c in cols]

@printf("  peak |omega| = %.3g, peak enstrophy = %.3g\n", maximum(abs, W), maximum(ens))
@printf("  enstrophy at first touch %.3g, at release %.3g\n",
        ens[findmin(abs.(ts .- t_in))[2]], ens[findmin(abs.(ts .- t_out))[2]])

# ---------------------------------------------------------------------------
# The figure. Minimal: two panels, one shared time axis, no gridlines competing with
# the map, and type large enough to survive being shrunk into a column.
# ---------------------------------------------------------------------------
default(fontfamily = "sans-serif", framestyle = :box,
        guidefontsize = 15, tickfontsize = 13, titlefontsize = 16, legendfontsize = 13)

pa = heatmap(ts, angs, abs.(W), c = PALETTE, colorbar_title = "  |ω|",
             xlabel = "", ylabel = "angle from contact pole  (deg)",
             title = @sprintf("Vorticity through a bounce   (We = %.1f, Oh = %.3f)", WE, OH),
             yticks = ([0, 45, 90, 135, 180], ["0", "45", "90", "135", "180"]),
             grid = false, ylims = (0, 180))
## The drop is touching for the whole record: it starts against the substrate and the
## march stops as soon as it releases. So shading the contact interval would shade the
## entire axis and say nothing. Only the release is marked.
vline!(pa, [t_out], c = :white, lw = 2, ls = :dash, label = "")
annotate!(pa, t_out, 168, text("release  ", 13, :white, :right))

pb = plot(ts, ens, lw = 3, c = DEEP, label = "",
          xlabel = "time  (capillary times)", ylabel = "enstrophy  ∫ω² dV",
          grid = true, gridalpha = 0.15)
vline!(pb, [t_out], c = DEEP, lw = 2, ls = :dash, label = "")

plt = plot(pa, pb; layout = grid(2, 1, heights = [0.62, 0.38]),
           size = (1100, 780), dpi = 150,
           left_margin = 9Plots.mm, bottom_margin = 7Plots.mm, right_margin = 5Plots.mm)

savefig(plt, joinpath(ASSET, "vorticity_bounce.png"))
println("wrote ", joinpath(ASSET, "vorticity_bounce.png"))

## The underlying series, so the figure is reproducible without rerunning the impact.
open(joinpath(OUT, "vorticity_bounce.csv"), "w") do io
    println(io, "t,enstrophy,max_abs_omega_surface")
    for (j, t) in enumerate(ts)
        println(io, join((t, ens[j], maximum(abs, view(W, :, j))), ","))
    end
end
println("wrote ", joinpath(OUT, "vorticity_bounce.csv"))
