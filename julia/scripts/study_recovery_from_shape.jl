# CAN THE SURFACE ALONE MEASURE THE RHEOLOGY?
#
# The companion study inverted two numbers per impact -- restitution and contact time --
# and got error bars too wide to call an instrument. But a high-speed camera does not
# return two numbers. It returns the whole outline, at every frame, for the entire
# bounce and the ringing afterwards. That is hundreds of measurements per drop, and the
# question here is whether they close the gap.
#
# WHAT IS ASSUMED OBSERVABLE. The surface is decomposed into Legendre amplitudes
# zeta_l(t) for l = 2..8 -- roughly what a camera can actually resolve on a millimetric
# drop -- plus the centre-of-mass height z(t), sampled on a fixed time grid spanning the
# contact AND the free flight after it. The free-flight ringing is included deliberately:
# it is small-amplitude, so it probes lower shear rates than the impact does, and if the
# thinning knee is visible anywhere it should be visible there.
#
# NOISE. Edge detection on a high-speed frame locates the interface to well under a per
# cent of the radius. The noise here is 0.5 per cent of R on the lowest mode, growing
# linearly with l because fine angular structure is progressively harder to resolve.
#
# THE HONEST CAVEAT, stated before the answer. Successive video frames do not carry
# independent noise: edge-detection error is correlated in time, and the modes are not
# independent of each other either. Treating every sample as independent OVERSTATES the
# information, so the error bars below are a best case. The time grid is deliberately
# coarse -- about a tenth of a capillary time between samples -- to keep that overstatement
# as small as is reasonable. A real inversion would need the measured noise covariance.

using Printf, LinearAlgebra
using DropSolver

flushln(a...) = (println(a...); flush(stdout))
const OUT = joinpath(@__DIR__, "..", "..", "results")
mkpath(OUT)

const R_D, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
const T_CAP = sqrt(RHO * R_D^3 / SIGMA)

## 1000 ppm, same fluid as the companion study so the two answers compare directly
const ETA_INF0 = 0.002747884967209792
const ETA_0    = 0.2601531742575956
const K_TRUE   = 1.076060672638565
const M_TRUE   = 0.6724627161400011
const OH_0     = ETA_0 / sqrt(RHO * SIGMA * R_D)

const M_RUN, K_RUN = 45, 3
const WES   = [0.3, 1.0, 2.5]
const LOBS  = 2:8                       # modes a camera can resolve
const TGRID = collect(0.1:0.1:5.0)      # capillary times, through contact and flight
const SIG_L = [0.005 * (l / 2) for l in LOBS]    # fraction of R, worse for fine structure
const SIG_Z = 0.005

eta_fn(k, m, ratio) = gd -> carreau(gd; lambda_c = k / T_CAP, a = m, n = 1 - m,
                                    eta_inf_ratio = ratio)

"""Linear interpolation of a time series onto the fixed observation grid."""
function resample(ts, ys, grid)
    out = similar(grid)
    for (i, t) in enumerate(grid)
        if t <= ts[1];   out[i] = ys[1];   continue; end
        if t >= ts[end]; out[i] = ys[end]; continue; end
        j = searchsortedlast(ts, t)
        w = (t - ts[j]) / (ts[j+1] - ts[j])
        out[i] = (1 - w) * ys[j] + w * ys[j+1]
    end
    out
end

"""
The full observable vector: every retained mode amplitude and the centre of mass, on the
observation grid, at every Weber number.
"""
function observe_shape(theta)
    k, m, ratio = 10^theta[1], theta[2], 10^theta[3]
    f = eta_fn(k, m, ratio)
    obs = Float64[]
    for We in WES
        ## stop_on_release = false so the record covers the ringing after the bounce
        p = ImpactParams(We = We, Bo = BO, Oh = OH_0, M = M_RUN, K = K_RUN,
                         eta = f, t_max = 5.2, stop_on_release = false)
        r = DropSolver.simulate_lcp(p)
        length(r.t) < 20 && return nothing
        ls = DropSolver.basis(p).ls
        amps = [DropSolver.surface_amplitudes(p, a) for a in r.a]
        for l in LOBS
            j = findfirst(==(l), ls)
            j === nothing && return nothing
            append!(obs, resample(r.t, [a[j] for a in amps], TGRID))
        end
        append!(obs, resample(r.t, r.z, TGRID))
    end
    obs
end

const THETA0 = [log10(K_TRUE), M_TRUE, log10(ETA_INF0 / ETA_0)]
const PNAMES = ["log10 k", "m", "log10 eta_inf/eta_0"]
const DTHETA = [0.04, 0.03, 0.06]

## Noise vector matching the observable layout: modes then centre of mass, per We
const SIGMA_VEC = let s = Float64[]
    for _ in WES
        for (i, _) in enumerate(LOBS); append!(s, fill(SIG_L[i], length(TGRID))); end
        append!(s, fill(SIG_Z, length(TGRID)))
    end
    s
end

flushln("== can the surface alone measure the rheology? ==")
@printf("fluid: 1000 ppm, Oh_0 = %.4g, knee 1/k = %.3g /s\n", OH_0, 1 / K_TRUE)
@printf("observables: zeta_l for l = %s, plus z, on %d times to t = %.1f, at We = %s\n",
        collect(LOBS), length(TGRID), TGRID[end], WES)
@printf("that is %d numbers per parameter evaluation\n", length(SIGMA_VEC))
@printf("noise: %.3f R on l = 2 growing to %.3f R on l = 8; %.3f R on the centre of mass\n\n",
        SIG_L[1], SIG_L[end], SIG_Z)

t0 = time()
base = observe_shape(THETA0)
base === nothing && error("baseline failed")
flushln("baseline computed ($(round(time()-t0)) s)")

J = zeros(length(base), 3)
for j in 1:3
    tp = copy(THETA0); tp[j] += DTHETA[j]
    tm = copy(THETA0); tm[j] -= DTHETA[j]
    op, om = observe_shape(tp), observe_shape(tm)
    (op === nothing || om === nothing) && error("perturbation $j failed")
    J[:, j] = (op .- om) ./ (2 * DTHETA[j])
    @printf("  sensitivity to %-22s done (%.0f s)\n", PNAMES[j], time() - t0)
end

Jt = J ./ SIGMA_VEC
C  = inv(transpose(Jt) * Jt)
sd = sqrt.(diag(C))

flushln("\n== error bars from the surface history ==")
@printf("%-22s %10s %10s %11s   %s\n", "parameter", "true", "sigma", "rel. err",
        "(scalars gave)")
prev = [1.1359, 0.2889, 0.1347]          # the companion study, restitution + contact time
for j in 1:3
    rel = j == 2 ? sd[j] / abs(THETA0[j]) : (10^sd[j] - 1)
    @printf("%-22s %10.4f %10.4f %10.1f%%   %10.4f  (x%.1f better)\n",
            PNAMES[j], THETA0[j], sd[j], 100 * rel, prev[j], prev[j] / sd[j])
end

sv = svdvals(Jt)
@printf("\nsingular values: %s\ncondition number: %.3g\n",
        join((@sprintf("%.3g", s) for s in sv), ", "), sv[1] / sv[end])

## Which part of the record carries the information: the impact, or the ringing after it?
nT = length(TGRID); nblk = (length(LOBS) + 1) * nT
contact_mask = falses(length(base))
for (w, _) in enumerate(WES), b in 0:length(LOBS)
    off = (w - 1) * nblk + b * nT
    for i in 1:nT
        TGRID[i] <= 2.4 && (contact_mask[off + i] = true)      # contact lasts about 2.4
    end
end
for (lbl, mask) in (("during contact", contact_mask), ("free ringing only", .!contact_mask))
    Js = Jt[mask, :]
    s = try sqrt.(abs.(diag(inv(transpose(Js) * Js)))) catch; fill(NaN, 3) end
    @printf("%-20s sigma(log10 k) %9.4f  sigma(m) %8.4f  sigma(log10 ratio) %8.4f\n",
            lbl, s[1], s[2], s[3])
end

## And the physical reason, measured: what shear rates does each phase probe?
flushln("\n== shear rates probed, by phase ==")
let p = ImpactParams(We = 1.0, Bo = BO, Oh = OH_0, M = M_RUN, K = K_RUN,
                     eta = eta_fn(K_TRUE, M_TRUE, ETA_INF0 / ETA_0),
                     t_max = 5.2, stop_on_release = false)
    r = DropSolver.simulate_lcp(p); b = DropSolver.basis(p)
    inc = findlast(>(0), r.cp)
    for (lbl, rng) in (("contact", 2:inc), ("free flight", (inc+1):length(r.t)))
        gds = Float64[]
        for i in round.(Int, range(first(rng), last(rng); length = 40)),
            x in (0.3, 0.6, 0.85, 0.98), mu in (-0.95, -0.5, 0.0, 0.5)
            push!(gds, DropSolver.shear_rate(b, r.adot[i], x, mu) / T_CAP)
        end
        sort!(gds); q(a) = gds[max(1, round(Int, a * length(gds)))]
        @printf("  %-12s 10th %.3g   median %.3g   90th %.3g /s\n", lbl, q(0.1), q(0.5), q(0.9))
    end
    @printf("  knee at %.3g /s\n", 1 / K_TRUE)
end

@printf("\ntotal wall time %.1f min\n", (time() - t0) / 60)
