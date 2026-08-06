# One animation per backend, on the same case, for looking at.
#
# The KPIs are two numbers per run and they hide everything that is not restitution or contact
# time. A drop that passes through the substrate, rings at the truncation, or develops a
# corner still reports a plausible CoR -- the annular-contact artefact survived several rounds
# of numerical checking and was obvious the moment the surface was drawn. So: draw the surface.
#
# Every backend is rendered from the SAME data structure, because `run_impact` returns the
# surface history in a common form. Nothing here knows which solver produced a frame.
#
# Usage:
#   julia --project=docs julia/scripts/animate_backends.jl            # Newtonian reference case
#   julia --project=docs julia/scripts/animate_backends.jl --shear    # the 3000 ppm fluid

using Printf, Plots
using DropSolver
gr()

const SHEAR = "--shear" in ARGS

## The shear-thinning fluid, from its own Cross characterisation
eta_0, eta_inf   = 8.433817577956766, 0.0037320997942061666
K_cross, m_cross = 18.48081673111359, 0.7430524574330837
R, SIGMA, G      = 0.0003, 0.0728, 9.81
BO_ST = 0.012
RHO   = BO_ST * SIGMA / (G * R^2)
T_CAP = sqrt(RHO * R^3 / SIGMA)
OH_ST = eta_0 / sqrt(RHO * SIGMA * R)
LAM, A_CY, N_CY = K_cross / T_CAP, m_cross, 1 - m_cross
ETA_RATIO = eta_inf / eta_0
EPS_ST    = (1 - N_CY) / A_CY
eta_st = gd -> carreau(gd; lambda_c = LAM, a = A_CY, n = N_CY, eta_inf_ratio = ETA_RATIO)

## The Newtonian case is chosen so that ALL THREE backends complete it. The obvious reference
## point, We = 1 at Oh = 0.3038, is one of the two in the comparison table where the
## nonvariational search gives up, and a missing video is the one outcome that teaches nothing.
const CASE = SHEAR ?
    (We = 0.3643, Bo = BO_ST,  Oh = OH_ST,  Mv = 14, Mn = 20, tag = "shear") :
    (We = 0.5,    Bo = 0.0188, Oh = 0.0373, Mv = 30, Mn = 30, tag = "newtonian")

const BACKENDS = [
    Backend(formulation = :nonvariational, contact = :tangency),
    Backend(),
    Backend(forcing = :nodal),
]

const FPS       = 30
const N_FRAMES  = 200          # subsampled from the march, so all backends share a cadence

"""Render one run to a gif. Returns the output path."""
function animate(r, tag::String, nm::String)
    keep = round.(Int, range(1, length(r.t); length = min(N_FRAMES, length(r.t))))
    ## a common vertical window, so the three videos are comparable at a glance
    ymax = 1.15 * maximum(maximum(abs, drop_outline(r.zeta[i], r.ls, r.z[i])[2]) for i in keep)
    fname = joinpath(@__DIR__, "..", "..", "docs", "figures",
                     "impact_" * tag * "_" * replace(nm, "/" => "-") * ".gif")
    anim = @animate for i in keep
        x, y = drop_outline(r.zeta[i], r.ls, r.z[i])
        plot(vcat(-reverse(x), x), vcat(reverse(y), y);
             seriestype = :shape, fillalpha = 0.28, fc = :steelblue, lc = :steelblue, lw = 1.6,
             aspect_ratio = 1, xlim = (-1.9, 1.9), ylim = (-0.12, ymax),
             legend = false, framestyle = :box, grid = false,
             size = (460, 420), dpi = 130,
             title = @sprintf("%s   t = %5.2f   cp = %d", nm, r.t[i], r.cp[i]),
             titlefontsize = 9)
        ## the substrate, and the 0.02R contact threshold the metrics use
        hline!([0.0]; lc = :black, lw = 2.5, label = "")
        hline!([0.02]; lc = :grey, lw = 0.8, ls = :dot, label = "")
    end
    gif(anim, fname; fps = FPS, show_msg = false)
    fname
end

@printf("case: %s  We=%.4g Bo=%.4g Oh=%.4g\n", CASE.tag, CASE.We, CASE.Bo, CASE.Oh)
for b in BACKENDS
    nm = label(b)
    isvar = b.formulation === :variational
    M = isvar ? CASE.Mv : CASE.Mn
    r = run_impact(b; We = CASE.We, Bo = CASE.Bo, Oh = CASE.Oh, M = M, K = 2, t_max = 25.0,
                   save_every = 0.01,
                   eta = (isvar && SHEAR) ? eta_st : nothing,
                   eta_nonvar = (!isvar && SHEAR) ?
                       STExactParams(M, OH_ST, LAM, A_CY, EPS_ST;
                                     viscous = :reid, eta_inf_ratio = ETA_RATIO) : nothing)
    if !r.ok
        @printf("  %-18s did not complete -- no animation\n", nm)
        continue
    end
    ## a sanity number alongside the picture: how far the surface went below the substrate
    worst = minimum(minimum(drop_outline(r.zeta[i], r.ls, r.z[i])[2]) for i in eachindex(r.t))
    f = animate(r, CASE.tag, nm)
    @printf("  %-18s CoR=%.4f tc=%.3f  frames=%d  worst surface height=%+.2e  -> %s\n",
            nm, r.cor, r.tc, length(r.t), worst, basename(f))
end
