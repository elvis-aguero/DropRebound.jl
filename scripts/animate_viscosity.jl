# The impact animation with the interior coloured by the LOCAL VISCOSITY.
#
# This is the picture behind the whole effective-viscosity question. A generalized
# Newtonian drop has no single viscosity: eta follows the local shear rate, so the drop
# is a moving map of stiff and soft regions. Every model that replaces the rheology with
# one number -- an effective Ohnesorge, a characteristic shear rate -- is averaging over
# exactly this map, and the video says how much is being thrown away.
#
# The drop arrives with a uniform velocity and no strain anywhere, so it starts at the
# rest viscosity eta_0 throughout and thins from the contact outward. Watching where it
# stays thick is watching where an effective viscosity is wrong.
#
# Usage:
#   julia --project=docs scripts/animate_viscosity.jl            # 3000 ppm
#   julia --project=docs scripts/animate_viscosity.jl --preview  # three stills

using Printf, LinearAlgebra
using Plots
using DropSolver
gr()

include(joinpath(@__DIR__, "_vorticity.jl"))

const PREVIEW = "--preview" in ARGS
const FIGS = joinpath(@__DIR__, "..", "outputs", "figures")
mkpath(FIGS)

## 3000 ppm: the most strongly thinning fluid in the series, so the map has the widest
## range to show. Its viscosity spans three and a half decades between rest and plateau.
const R_D, SIGMA, G = 0.0003, 0.0728, 9.81
const BO    = 0.012
const RHO   = BO * SIGMA / (G * R_D^2)
const T_CAP = sqrt(RHO * R_D^3 / SIGMA)
const ETA_INF, ETA_0 = 0.0037320997942061666, 8.433817577956766
const K_CROSS, M_CROSS = 18.48081673111359, 0.7430524574330837
const OH_0  = ETA_0 / sqrt(RHO * SIGMA * R_D)
const RATIO = ETA_INF / ETA_0
## `carreau` takes the shear rate in SOLVER units, which is why lambda_c carries the
## capillary time: lambda_c * gammadot_solver = k * gammadot_physical.
eta_st(gd) = carreau(gd; lambda_c = K_CROSS / T_CAP, a = M_CROSS,
                     n = 1 - M_CROSS, eta_inf_ratio = RATIO)

const M_RUN, K_RUN = DropSolver.DEFAULT_M, DropSolver.DEFAULT_K
const T_END    = 4.0
const N_FRAMES = 300
const FPS      = 30
const NPIX     = 300

@printf("3000 ppm: Oh_0 = %.4g, eta_inf/eta_0 = %.4g (%.1f decades of thinning)\n",
        OH_0, RATIO, -log10(RATIO))
p = ImpactParams(We = 0.3643, Bo = BO, Oh = OH_0, M = M_RUN, K = K_RUN,
                 eta = eta_st, t_max = T_END, stop_on_release = false)
r = DropSolver.simulate_lcp(p)
@printf("  %d frames, restitution %.4f, contact time %.3f\n", length(r.t), r.cor, r.tc)

keep = unique(round.(Int, range(1, length(r.t); length = min(N_FRAMES, length(r.t)))))

# ---------------------------------------------------------------------------
# Kernels: the strain basis and the angular basis, both state-independent.
# ---------------------------------------------------------------------------
const NX, NMU = 50, 120
xs  = collect(range(0.02, 1.0; length = NX))
mus = collect(range(-0.99999, 0.99999; length = NMU))
Kern = strain_kernel(p, xs, mus)
Pang, _ = angular_kernel(p, mus)
@printf("  strain kernel %d x %d (%.0f MB)\n", size(Kern, 1), size(Kern, 2),
        length(Kern) * 8 / 1e6)

surf_amp(i) = DropSolver.surface_amplitudes(p, r.a[i])
radius_of(zeta) = 1.0 .+ vec(transpose(zeta) * Pang)

## Drawn relative to the PLATEAU, not to rest. Scaling to rest gives a picture that is
## uniformly dark, because this fluid is thinned by two to three decades everywhere from
## the first instant -- which is the finding, but it is not a legible image. Against the
## plateau, zero means "as thin as this fluid can get" and positive means "still thicker
## than the plateau", which is exactly the deviation an effective viscosity would have to
## account for.
"""log10 of the local viscosity above the infinite-shear plateau, on the (x, mu) grid."""
function visc_grid(i)
    gd = shear_rate_grid(Kern, r.adot[i], NX, NMU)
    log10.(clamp.(eta_st.(gd), RATIO, 1.0) ./ RATIO)
end

## Ceiling from the field itself, so the scale spans what actually occurs.
const LO = 0.0
const HI = let v = Float64[]
    for i in keep[1:8:end]; append!(v, vec(visc_grid(i))); end
    sort!(v); q = v[max(1, round(Int, 0.90 * length(v)))]
    @printf("  eta/eta_inf: median %.3g, 90th %.3g, rest is %.0fx the plateau\n",
            10^v[length(v) ÷ 2], 10^q, 1 / RATIO)
    q
end

ymax = 1.15 * maximum(maximum(radius_of(surf_amp(i))) + r.z[i] for i in keep)
const XLIM, YLIM = (-1.75, 1.75), (-0.10, ymax)
gx = collect(range(XLIM[1], XLIM[2]; length = NPIX))
gy = collect(range(YLIM[1], YLIM[2]; length = NPIX))

function frame_image(i)
    zeta = surf_amp(i); Rmu = radius_of(zeta); Wg = visc_grid(i); z = r.z[i]
    img = fill(NaN, length(gy), length(gx))
    for (cx, X) in enumerate(gx), (cy, Y) in enumerate(gy)
        dy = Y - z; s = hypot(X, dy)
        s < 1e-9 && continue
        mu = dy / s
        fj = (mu - mus[1]) / (mus[end] - mus[1]) * (NMU - 1) + 1
        j = clamp(floor(Int, fj), 1, NMU - 1); tj = clamp(fj - j, 0.0, 1.0)
        Rh = (1 - tj) * Rmu[j] + tj * Rmu[j+1]
        s > Rh && continue
        x = s / Rh                                   # radial stretch onto the outline
        fi = (x - xs[1]) / (xs[end] - xs[1]) * (NX - 1) + 1
        ii = clamp(floor(Int, fi), 1, NX - 1); ti = clamp(fi - ii, 0.0, 1.0)
        img[cy, cx] = (1 - ti) * ((1 - tj) * Wg[ii, j] + tj * Wg[ii, j+1]) +
                           ti  * ((1 - tj) * Wg[ii+1, j] + tj * Wg[ii+1, j+1])
    end
    img
end

## The one-number summary the video argues against: the volume-median viscosity, which
## is what an effective-viscosity model would be trying to use.
median_visc(i) = (v = sort(vec(visc_grid(i))); v[max(1, length(v) ÷ 2)])

default(fontfamily = "sans-serif", framestyle = :box, titlefontsize = 11,
        guidefontsize = 11, tickfontsize = 10)

render(i) = begin
    heatmap(gx, gy, frame_image(i); c = ABYSS, clims = (LO, HI),
            aspect_ratio = 1, xlim = XLIM, ylim = YLIM,
            colorbar = true, colorbar_title = "log₁₀ η/η∞", colorbar_titlefontsize = 11,
            legend = false, grid = false, size = (700, 520), dpi = 130,
            right_margin = 22Plots.mm, left_margin = 3Plots.mm, bottom_margin = 3Plots.mm,
            xlabel = "x / R", ylabel = "y / R",
            title = @sprintf("t = %5.2f    median η/η∞ = %.3g", r.t[i], 10^median_visc(i)))
    xo, yo = drop_outline(surf_amp(i), DropSolver.basis(p).ls, r.z[i])
    plot!(vcat(-reverse(xo), xo), vcat(reverse(yo), yo);
          lc = RGB(0.75, 0.88, 0.95), lw = 1.4, label = "")
    hline!([0.0]; lc = :black, lw = 2.5, label = "")
end

if PREVIEW
    sel = [keep[4], keep[length(keep) ÷ 3], keep[end - 12]]
    plt = plot((render(i) for i in sel)...; layout = (1, 3), size = (2100, 520))
    out = joinpath(@__DIR__, "..", "outputs", "figures", "animate_viscosity_preview.png")
    savefig(plt, out); println("wrote ", out); exit()
end

println("rendering $(length(keep)) frames ..."); flush(stdout)
anim = @animate for (n, i) in enumerate(keep)
    render(i)
    n % 60 == 0 && (@printf("  %d / %d\n", n, length(keep)); flush(stdout))
end
out = joinpath(FIGS, "animate_viscosity_shear.gif")
gif(anim, out; fps = FPS, show_msg = false)
println("wrote ", out)

## The summary curve, so the video's claim is also a number: how much of the drop is
## thinned, and by how much, through the bounce.
open(joinpath(@__DIR__, "..", "outputs", "csv", "viscosity_field.csv"), "w") do io
    println(io, "t,median_eta_over_etainf,max_eta_over_etainf,frac_within_2x_plateau")
    for i in keep
        v = vec(visc_grid(i))
        println(io, join((r.t[i], 10^(sort(v)[max(1, length(v) ÷ 2)]),
                          10^maximum(v), count(<(log10(2.0)), v) / length(v)), ","))
    end
end
println("wrote outputs/csv/viscosity_field.csv")
