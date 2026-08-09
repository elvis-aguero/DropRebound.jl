# The impact animation, with the interior coloured by |omega|.
#
# `animate_backends.jl` draws the same case as a flat silhouette, which shows the shape
# and nothing else. The shape is only half the story: the drop arrives irrotational, and
# everything that happens to the interior during contact is rotation being created at
# the substrate and carried into the bulk. Colouring the inside shows that directly.
#
# Usage:
#   julia --project=docs julia/scripts/animate_vorticity.jl              # Newtonian
#   julia --project=docs julia/scripts/animate_vorticity.jl --shear      # the 3000 ppm fluid
#   julia --project=docs julia/scripts/animate_vorticity.jl --preview    # three stills, no gif
#
# The march is run with `stop_on_release = false`, so it carries on past the bounce and
# the video shows the drop actually leaving. The default stops at release, which is all
# the impact metrics need but makes for a video that ends mid-air.
#
# HOW THE INTERIOR IS DRAWN
#
# The theory is linear, so the velocity field lives on the REFERENCE unit sphere while
# the outline is the deformed surface. Drawing the field on the deformed shape means
# choosing a map between them, and the one used here is the obvious radial stretch: a
# point at angle theta and physical distance s from the centre is assigned the
# reference radius x = s / R(theta), so x = 1 is the free surface wherever it happens
# to be. To the order the theory is valid, this agrees with any other reasonable
# choice; it is a plotting convention, not a result.

using Printf, LinearAlgebra
using Plots
using DropSolver
gr()

include(joinpath(@__DIR__, "_vorticity.jl"))

const SHEAR   = "--shear" in ARGS
const PREVIEW = "--preview" in ARGS
const FIGS  = joinpath(@__DIR__, "..", "..", "docs", "figures")
mkpath(FIGS)

## The 3000 ppm fluid, from its own Cross characterisation, exactly as animate_backends
## sets it up so the two videos are of the same physical case.
eta_0, eta_inf   = 8.433817577956766, 0.0037320997942061666
K_cross, m_cross = 18.48081673111359, 0.7430524574330837
R_D, SIGMA, G    = 0.0003, 0.0728, 9.81
BO_ST = 0.012
RHO   = BO_ST * SIGMA / (G * R_D^2)
T_CAP = sqrt(RHO * R_D^3 / SIGMA)
OH_ST = eta_0 / sqrt(RHO * SIGMA * R_D)
eta_st = gd -> carreau(gd; lambda_c = K_cross / T_CAP, a = m_cross,
                       n = 1 - m_cross, eta_inf_ratio = eta_inf / eta_0)

const CASE = SHEAR ?
    (We = 0.3643, Bo = BO_ST,  Oh = OH_ST,  tag = "shear",     eta = eta_st) :
    (We = 0.5,    Bo = 0.0188, Oh = 0.0373, tag = "newtonian", eta = nothing)

## K = 3 is not optional here. K = 1 is irrotational by construction, so a video built
## on it would be uniformly zero; K is the number of radial functions and the vorticity
## lives entirely in the ones past the first.
const M_RUN, K_RUN = DropSolver.DEFAULT_M, DropSolver.DEFAULT_K
const N_FRAMES = 300
const FPS      = 30
const NPIX     = 300           # raster resolution of the interior
## Enough flight after release to see the drop leave, and no more: the vertical window has
## to hold every frame, so a long tail of rising drop shrinks the contact phase to nothing.
const T_END    = 4.0

assert_irrotational()

@printf("case %s: We=%.4g Bo=%.4g Oh=%.4g at M=%d K=%d\n",
        CASE.tag, CASE.We, CASE.Bo, CASE.Oh, M_RUN, K_RUN)
## `eta` has to be absent for a Newtonian run, not `nothing`: the constructor probes it
## as a callable to decide whether the viscosity is constant.
p = SHEAR ?
    ImpactParams(We = CASE.We, Bo = CASE.Bo, Oh = CASE.Oh, M = M_RUN, K = K_RUN,
                 t_max = T_END, eta = CASE.eta, stop_on_release = false) :
    ImpactParams(We = CASE.We, Bo = CASE.Bo, Oh = CASE.Oh, M = M_RUN, K = K_RUN,
                 t_max = T_END, stop_on_release = false)
r = DropSolver.simulate_lcp(p)
@printf("  %d frames, restitution %.4f, contact time %.3f\n", length(r.t), r.cor, r.tc)

keep = unique(round.(Int, range(1, length(r.t); length = min(N_FRAMES, length(r.t)))))

# ---------------------------------------------------------------------------
# Kernels. Both are state-independent, so they are built once and every frame is a
# contraction against them.
# ---------------------------------------------------------------------------
const NMU = 240
## The poles are excluded: vorticity vanishes on the axis anyway (it carries a factor
## dP/dtheta) and the angular routine divides by sin(theta).
mus  = collect(range(-0.99999, 0.99999; length = NMU))
xs, Phi = radial_kernel(p)
Pang, Dang = angular_kernel(p, mus)

surf_amp(i) = DropSolver.surface_amplitudes(p, r.a[i])
radius_of(zeta) = 1.0 .+ vec(transpose(zeta) * Pang)     # R at every mu on the grid

## COLOUR SCALE. Vorticity here is not spread evenly through the drop: it is made in a
## thin layer at the contact and diffuses inward, so the field spans orders of magnitude
## across one frame -- on the Newtonian case the median is 0.06 against a peak of 25. A
## linear scale set by the peak draws everything except that layer black, which pictures
## the colour scale rather than the flow.
##
## So the quantity drawn is log10|omega|, clipped at both ends. The ceiling is a high
## quantile of the actual interior values rather than the peak, so the brief transient at
## first touch clips instead of setting the scale for the whole video; the floor is three
## decades below it, which covers the bulk without spending colour on numerical dust.
## Unlike a power law, the numbers on the colourbar ARE the quantity plotted.
const DECADES = 3.0
const LOG_FLOOR = 1e-12
lo, hi = let vals = Float64[]
    for i in keep[1:4:end]
        append!(vals, vec(vorticity_grid(p, r.adot[i], Phi, Dang)))
    end
    sort!(vals)
    q  = vals[max(1, round(Int, 0.995 * length(vals)))]
    hi = log10(q)
    @printf("  |omega|: peak %.3g, median %.3g, ceiling %.3g -> log10 range [%.2f, %.2f]\n",
            vals[end], vals[max(1, length(vals) ÷ 2)], q, hi - DECADES, hi)
    (hi - DECADES, hi)
end

## A common window, from the extreme outline over all kept frames.
ymax = 1.15 * maximum(maximum(radius_of(surf_amp(i))) + r.z[i] for i in keep)
const XLIM = (-1.75, 1.75)
const YLIM = (-0.10, ymax)

# ---------------------------------------------------------------------------
# One frame: rasterise |omega| inside the outline, NaN outside so the background shows.
# ---------------------------------------------------------------------------
gx = collect(range(XLIM[1], XLIM[2]; length = NPIX))
gy = collect(range(YLIM[1], YLIM[2]; length = NPIX))

function frame_image(i)
    zeta = surf_amp(i)
    Rmu  = radius_of(zeta)                       # deformed radius on the mu grid
    Wg   = vorticity_grid(p, r.adot[i], Phi, Dang)   # |omega| on (x, mu)
    z    = r.z[i]
    img  = fill(NaN, length(gy), length(gx))
    nx, nmu = length(xs), length(mus)
    for (cx, X) in enumerate(gx), (cy, Y) in enumerate(gy)
        dy = Y - z
        s  = hypot(X, dy)
        s < 1e-9 && continue
        mu = dy / s
        ## interpolate the deformed radius and the field at this direction
        fj = (mu - mus[1]) / (mus[end] - mus[1]) * (nmu - 1) + 1
        j  = clamp(floor(Int, fj), 1, nmu - 1); tj = clamp(fj - j, 0.0, 1.0)
        Rh = (1 - tj) * Rmu[j] + tj * Rmu[j+1]
        s > Rh && continue                        # outside the drop
        x  = s / Rh                               # radial stretch onto the outline
        fi = (x - xs[1]) / (xs[end] - xs[1]) * (nx - 1) + 1
        ii = clamp(floor(Int, fi), 1, nx - 1); ti = clamp(fi - ii, 0.0, 1.0)
        w = (1 - ti) * ((1 - tj) * Wg[ii,   j] + tj * Wg[ii,   j+1]) +
                 ti  * ((1 - tj) * Wg[ii+1, j] + tj * Wg[ii+1, j+1])
        img[cy, cx] = clamp(log10(max(w, LOG_FLOOR)), lo, hi)   # clipped at both ends
    end
    img
end

default(fontfamily = "sans-serif", framestyle = :box, titlefontsize = 11,
        guidefontsize = 11, tickfontsize = 10)

render(i) = begin
    heatmap(gx, gy, frame_image(i); c = ABYSS, clims = (lo, hi),
            aspect_ratio = 1, xlim = XLIM, ylim = YLIM,
            ## GR places the colourbar title hard against its tick labels, so the gap has
            ## to be bought with margin rather than asked for.
            colorbar = true, colorbar_title = "log₁₀ |ω|", colorbar_titlefontsize = 11,
            legend = false, grid = false,
            size = (700, 520), dpi = 130,
            right_margin = 22Plots.mm, left_margin = 3Plots.mm, bottom_margin = 3Plots.mm,
            xlabel = "x / R", ylabel = "y / R",
            title = r.cp[i] > 0 ?
                @sprintf("t = %5.2f    contact nodes = %d", r.t[i], r.cp[i]) :
                @sprintf("t = %5.2f    free flight", r.t[i]))
    ## the outline itself, drawn over the fill, and the substrate
    xo, yo = drop_outline(surf_amp(i), DropSolver.basis(p).ls, r.z[i])
    plot!(vcat(-reverse(xo), xo), vcat(reverse(yo), yo);
          lc = RGB(0.75, 0.88, 0.95), lw = 1.4, label = "")
    hline!([0.0]; lc = :black, lw = 2.5, label = "")
end

if PREVIEW
    ## three stills spanning the record, to check layout without paying for 300 frames
    sel = [keep[6], keep[length(keep) ÷ 3], keep[end - 12]]
    plt = plot((render(i) for i in sel)...; layout = (1, 3), size = (1740, 520))
    out = joinpath(@__DIR__, "..", "..", "results", "preview_vorticity_" * CASE.tag * ".png")
    savefig(plt, out); println("wrote ", out)
    exit()
end

println("rendering $(length(keep)) frames ..."); flush(stdout)
anim = @animate for (n, i) in enumerate(keep)
    render(i)
    n % 60 == 0 && (@printf("  %d / %d\n", n, length(keep)); flush(stdout))
end

out = joinpath(FIGS, "impact_vorticity_" * CASE.tag * ".gif")
gif(anim, out; fps = FPS, show_msg = false)
println("wrote ", out)
