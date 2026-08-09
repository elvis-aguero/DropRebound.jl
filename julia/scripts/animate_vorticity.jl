# The impact animation, with the interior coloured by |omega|.
#
# `animate_backends.jl` draws the same case as a flat silhouette, which shows the shape
# and nothing else. The shape is only half the story: the drop arrives irrotational, and
# everything that happens to the interior during contact is rotation being created at
# the substrate and carried into the bulk. Colouring the inside shows that directly.
#
# Usage:
#   julia --project=docs julia/scripts/animate_vorticity.jl            # Newtonian
#   julia --project=docs julia/scripts/animate_vorticity.jl --shear    # the 3000 ppm fluid
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

const SHEAR = "--shear" in ARGS
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
const N_FRAMES = 220
const FPS      = 30
const NPIX     = 300           # raster resolution of the interior

assert_irrotational()

@printf("case %s: We=%.4g Bo=%.4g Oh=%.4g at M=%d K=%d\n",
        CASE.tag, CASE.We, CASE.Bo, CASE.Oh, M_RUN, K_RUN)
## `eta` has to be absent for a Newtonian run, not `nothing`: the constructor probes it
## as a callable to decide whether the viscosity is constant.
p = SHEAR ?
    ImpactParams(We = CASE.We, Bo = CASE.Bo, Oh = CASE.Oh, M = M_RUN, K = K_RUN,
                 t_max = 25.0, eta = CASE.eta) :
    ImpactParams(We = CASE.We, Bo = CASE.Bo, Oh = CASE.Oh, M = M_RUN, K = K_RUN,
                 t_max = 25.0)
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
## across one frame. A linear scale set by the peak renders everything except that layer
## as black, which is a picture of the colour scale rather than of the flow.
##
## So the map is a power law, |omega|^GAMMA, with the colourbar ticked in |omega| so the
## numbers a reader takes off it are still vorticity. The ceiling is a high quantile of
## the actual interior values rather than the peak, and the brief transient at first
## touch clips against it.
const GAMMA = 0.4
vmax = let vals = Float64[]
    for i in keep[1:4:end]
        g = vorticity_grid(p, r.adot[i], Phi, Dang)
        append!(vals, vec(g))
    end
    sort!(vals)
    q = vals[max(1, round(Int, 0.995 * length(vals)))]
    @printf("  |omega|: peak %.3g, median %.3g, colour ceiling %.3g (gamma = %.2f)\n",
            vals[end], vals[max(1, length(vals) ÷ 2)], q, GAMMA)
    q
end

## GR will not take custom colourbar ticks, so the bar is labelled for the quantity it
## actually shows, |omega|^GAMMA, rather than carrying |omega| numbers it is not drawing.
## The peak and median above give the absolute scale.
const CB_LABEL = @sprintf("  |ω|^%.1f", GAMMA)

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
        img[cy, cx] = w^GAMMA          # drawn on the power scale, ticked in |omega|
    end
    img
end

default(fontfamily = "sans-serif", framestyle = :box, titlefontsize = 11,
        guidefontsize = 11, tickfontsize = 10)

println("rendering $(length(keep)) frames ..."); flush(stdout)
anim = @animate for (n, i) in enumerate(keep)
    heatmap(gx, gy, frame_image(i); c = ABYSS, clims = (0.0, vmax^GAMMA),
            aspect_ratio = 1, xlim = XLIM, ylim = YLIM,
            colorbar = true, colorbar_title = CB_LABEL,
            legend = false, grid = false,
            size = (600, 520), dpi = 130, right_margin = 6Plots.mm,
            xlabel = "", ylabel = "",
            title = @sprintf("t = %5.2f    contact nodes = %d", r.t[i], r.cp[i]))
    ## the outline itself, drawn over the fill, and the substrate
    xo, yo = drop_outline(surf_amp(i), DropSolver.basis(p).ls, r.z[i])
    plot!(vcat(-reverse(xo), xo), vcat(reverse(yo), yo);
          lc = RGB(0.75, 0.88, 0.95), lw = 1.4, label = "")
    hline!([0.0]; lc = :black, lw = 2.5, label = "")
    n % 40 == 0 && (@printf("  %d / %d\n", n, length(keep)); flush(stdout))
end

out = joinpath(FIGS, "impact_vorticity_" * CASE.tag * ".gif")
gif(anim, out; fps = FPS, show_msg = false)
println("wrote ", out)
