# CAN A BOUNCE MEASURE A FLOW CURVE?
#
# The proposal is to use drop rebound as a high-shear-rate rheometer: impact a drop of
# unknown complex fluid, measure what the bounce does, and invert for the viscosity
# curve. The forward model is this solver. The question this script answers is whether
# the inverse problem is WELL POSED -- whether the data actually constrain the
# parameters, or whether wildly different fluids bounce indistinguishably.
#
# This is deliberately the cheapest decisive test. Rather than run a full nonlinear
# inversion, it linearises: build the sensitivity of every observable to every
# parameter, weight by realistic experimental noise, and read the achievable error bars
# straight off the Fisher information. If the parameters come back with useless error
# bars here, no clever optimiser will rescue them, and the idea dies for the price of
# thirty-five impacts.
#
# THE FLUID is the 1000 ppm solution, mid-range in the validation series. Its Cross fit
# is treated as ground truth, and the question is whether bounces alone could have
# produced it.
#
# THE PARAMETERS are the three that describe the SHAPE of the flow curve:
#
#   log10 k          the thinning onset, in seconds
#   m                the thinning exponent
#   log10(eta_inf/eta_0)   the high-shear plateau, relative to rest
#
# The zero-shear viscosity eta_0 is taken as known, because a conventional rheometer
# measures it easily. The whole point of the instrument is the part rheometers cannot
# reach: what happens above about 1000 per second.
#
# THE OBSERVABLES are what a high-speed camera actually returns: the restitution and the
# contact time, at several impact speeds.
#
# WHAT COMES OUT is (a) the error bar on each parameter, (b) which combinations of
# parameters the bounce constrains and which it cannot separate, and (c) the resulting
# uncertainty on the flow curve itself at each shear rate -- which is the number a
# rheologist would ask for.

using Printf, LinearAlgebra
using DropSolver

flushln(a...) = (println(a...); flush(stdout))
const OUT = joinpath(@__DIR__, "..", "outputs", "csv")
mkpath(OUT)

# ---------------------------------------------------------------------------
# The fluid and the experiment
# ---------------------------------------------------------------------------
const R_D, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
const T_CAP = sqrt(RHO * R_D^3 / SIGMA)

## 1000 ppm, from its own Cross characterisation
const ETA_INF0 = 0.002747884967209792
const ETA_0    = 0.2601531742575956
const K_TRUE   = 1.076060672638565        # thinning onset, s
const M_TRUE   = 0.6724627161400011       # thinning exponent
const OH_0     = ETA_0 / sqrt(RHO * SIGMA * R_D)

## M = 45 rather than the production 90: it differs by 0.17 per cent in restitution on
## this fluid, which is an order of magnitude below the experimental noise the study
## assumes, and it makes a 35-impact Jacobian affordable.
const M_RUN, K_RUN = 45, 3

const WES = [0.2, 0.5, 1.0, 2.0, 3.0]

## EXPERIMENTAL NOISE, and the whole answer scales with it. Restitution from high-speed
## imaging repeats to a couple of points in the second decimal; contact time to a couple
## of per cent. Both are absolute standard deviations on a single measurement.
const SIG_COR = 0.02
const SIG_TC_REL = 0.02

"""Carreau viscosity in the solver's convention, from Cross parameters."""
eta_fn(k, m, ratio) = gd -> carreau(gd; lambda_c = k / T_CAP, a = m, n = 1 - m,
                                    eta_inf_ratio = ratio)

"""The observable vector: restitution and contact time at every Weber number."""
function observe(theta)
    k, m, ratio = 10^theta[1], theta[2], 10^theta[3]
    f = eta_fn(k, m, ratio)
    cors = Float64[]; tcs = Float64[]
    for We in WES
        r = run_impact(Backend(contact = :lcp); We = We, Bo = BO, Oh = OH_0,
                       M = M_RUN, K = K_RUN, eta = f, t_max = 25.0)
        r.ok || return nothing
        push!(cors, r.cor); push!(tcs, r.tc)
    end
    vcat(cors, tcs)
end

const THETA0 = [log10(K_TRUE), M_TRUE, log10(ETA_INF0 / ETA_0)]
const PNAMES = ["log10 k", "m", "log10 eta_inf/eta_0"]
## Step sizes: large enough to clear the solver's own tolerance, small enough to stay
## linear. Roughly 10 per cent in k, 5 per cent in m, 15 per cent in the plateau.
const DTHETA = [0.04, 0.03, 0.06]

flushln("== drop-impact rheometry: is the inverse problem well posed? ==")
@printf("fluid: 1000 ppm, Oh_0 = %.4g, k = %.4g s, m = %.4g, eta_inf/eta_0 = %.4g\n",
        OH_0, K_TRUE, M_TRUE, ETA_INF0 / ETA_0)
@printf("observables: restitution and contact time at We = %s\n", WES)
@printf("assumed noise: sigma(cor) = %.3f absolute, sigma(tc)/tc = %.3f\n\n",
        SIG_COR, SIG_TC_REL)

t0 = time()
base = observe(THETA0)
base === nothing && error("the baseline fluid did not bounce at every We")
ncor = length(WES)
flushln("baseline:")
for (i, We) in enumerate(WES)
    @printf("  We %4.1f   cor %.4f   tc %.4f\n", We, base[i], base[ncor + i])
end

## Noise vector, in the same order as the observables
sig = vcat(fill(SIG_COR, ncor), SIG_TC_REL .* base[ncor+1:end])

# ---------------------------------------------------------------------------
# Sensitivity by central differences, in units of the noise.
#
# Jt[i, j] is "how many error bars does observable i move when parameter j changes by
# one step". A column of small numbers is a parameter the bounce cannot see.
# ---------------------------------------------------------------------------
flushln("\ncomputing sensitivities ...")
J = zeros(length(base), length(THETA0))
for j in eachindex(THETA0)
    tp = copy(THETA0); tp[j] += DTHETA[j]
    tm = copy(THETA0); tm[j] -= DTHETA[j]
    op, om = observe(tp), observe(tm)
    (op === nothing || om === nothing) && error("perturbed fluid $j failed to bounce")
    J[:, j] = (op .- om) ./ (2 * DTHETA[j])
    @printf("  %-22s done (%.0f s elapsed)\n", PNAMES[j], time() - t0)
end

Jt = J ./ sig                                  # sensitivity in noise units
F  = transpose(Jt) * Jt                        # Fisher information
C  = inv(F)                                    # parameter covariance
sd = sqrt.(diag(C))

flushln("\n== how well is each parameter determined by ONE bounce per We? ==")
@printf("%-22s %12s %12s %12s\n", "parameter", "true", "sigma", "rel. err")
for j in eachindex(THETA0)
    rel = j == 2 ? sd[j] / abs(THETA0[j]) : (10^sd[j] - 1)   # log params -> per cent
    @printf("%-22s %12.4f %12.4f %11.1f%%\n", PNAMES[j], THETA0[j], sd[j], 100 * rel)
end

sv = svdvals(Jt)
@printf("\nsingular values: %s\n", join((@sprintf("%.3g", s) for s in sv), ", "))
@printf("condition number: %.3g\n", sv[1] / sv[end])
U, S, V = svd(Jt)
flushln("\nthe worst-determined combination of parameters (last right singular vector):")
for j in eachindex(THETA0)
    @printf("   %-22s %+.3f\n", PNAMES[j], V[j, end])
end

## Does contact time add anything over restitution alone? This is the claim that the
## bounce carries more information than a single spreading number.
for (lbl, rows) in (("restitution only", 1:ncor),
                    ("contact time only", (ncor+1):length(base)),
                    ("both", 1:length(base)))
    Js = Jt[rows, :]
    Cs = try inv(transpose(Js) * Js) catch; fill(NaN, 3, 3) end
    s  = sqrt.(abs.(diag(Cs)))
    @printf("\n%-18s sigma(log10 k) %8.4f   sigma(m) %8.4f   sigma(log10 ratio) %8.4f",
            lbl, s[1], s[2], s[3])
end
flushln()

# ---------------------------------------------------------------------------
# The number a rheologist would ask for: the error bar on the FLOW CURVE.
#
# Propagate the parameter covariance through eta(gammadot) at each shear rate. Where the
# band is tight the instrument measures the viscosity; where it is wide it does not.
# ---------------------------------------------------------------------------
flushln("\n== implied uncertainty on the flow curve ==")
@printf("%12s %14s %12s\n", "gammadot", "eta/eta_0", "1-sigma")
open(joinpath(OUT, "recovery_flow_curve_band.csv"), "w") do io
    println(io, "gammadot,eta_over_eta0,sigma_log10_eta")
    for gd in exp10.(range(-1, 5; length = 25))
        g(th) = log10(carreau(gd; lambda_c = 10^th[1] / T_CAP, a = th[2], n = 1 - th[2],
                              eta_inf_ratio = 10^th[3]))
        grad = [(g(THETA0 .+ [i == j ? 1e-4 : 0.0 for j in 1:3]) -
                 g(THETA0 .- [i == j ? 1e-4 : 0.0 for j in 1:3])) / 2e-4 for i in 1:3]
        sg = sqrt(max(dot(grad, C * grad), 0.0))
        val = carreau(gd; lambda_c = K_TRUE / T_CAP, a = M_TRUE, n = 1 - M_TRUE,
                      eta_inf_ratio = ETA_INF0 / ETA_0)
        println(io, join((gd, val, sg), ","))
        gd in exp10.(range(-1, 5; length = 25))[1:3:end] &&
            @printf("%12.3g %14.5f %11.1f%%\n", gd, val, 100 * (10^sg - 1))
    end
end

# ---------------------------------------------------------------------------
# Which shear rates does the drop actually sample? The instrument can only report the
# viscosity where the flow puts the fluid, so this bounds what is measurable in
# principle, independently of noise.
# ---------------------------------------------------------------------------
flushln("\n== shear rates the drop actually probes ==")
let p = ImpactParams(We = 1.0, Bo = BO, Oh = OH_0, M = M_RUN, K = K_RUN,
                     eta = eta_fn(K_TRUE, M_TRUE, ETA_INF0 / ETA_0), t_max = 25.0)
    r = DropSolver.simulate_lcp(p)
    b = DropSolver.basis(p)
    gds = Float64[]
    for i in round.(Int, range(2, length(r.t); length = 60))
        for x in (0.3, 0.6, 0.85, 0.98), mu in (-0.95, -0.6, -0.2, 0.3, 0.8)
            push!(gds, DropSolver.shear_rate(b, r.adot[i], x, mu) / T_CAP)   # per second
        end
    end
    sort!(gds)
    q(f) = gds[max(1, round(Int, f * length(gds)))]
    @printf("  median %.3g /s, 90th %.3g /s, 99th %.3g /s, max %.3g /s\n",
            q(0.5), q(0.9), q(0.99), gds[end])
    @printf("  the Cross onset 1/k is %.3g /s, so the drop probes %s the knee\n",
            1 / K_TRUE, q(0.9) > 1 / K_TRUE ? "past" : "below")
end

@printf("\ntotal wall time %.1f min\n", (time() - t0) / 60)
