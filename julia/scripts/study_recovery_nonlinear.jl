# NONLINEAR RECOVERY: does the inversion actually return the right fluid?
#
# The two previous studies were LINEARISED. They ask how far the observables move for a
# small parameter change and read error bars off that, which assumes the likelihood is a
# well-behaved paraboloid. It answers "how precise, if the inversion works". It cannot
# see a second minimum somewhere else in parameter space, a curved degeneracy valley, or
# a posterior so skewed that the error bar is meaningless.
#
# So this maps the whole thing. Rather than run an optimiser and report where it landed
# -- which hides exactly the pathologies worth finding -- it evaluates the forward model
# on a grid covering the plausible region, then for each of many synthetic noisy
# datasets computes chi^2 everywhere and looks at the entire surface.
#
# The grid is the expensive part and it is computed once. Every noise realisation after
# that is free, so the recovery statistics come from hundreds of synthetic experiments
# rather than the handful an optimiser-per-dataset could afford.
#
# ONE Weber number, not three, to keep the grid affordable. That costs about a factor of
# sqrt(3) in precision against the previous study, and is stated in the output rather
# than hidden.
#
# Stage 1 (this script, slow) writes the grid to results/recovery_grid.csv.
# Stage 2 (study_recovery_posterior.jl, seconds) reads it and does the statistics.

using Printf, LinearAlgebra
using DropSolver

flushln(a...) = (println(a...); flush(stdout))
const OUT = joinpath(@__DIR__, "..", "..", "results")
mkpath(OUT)

const R_D, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
const T_CAP = sqrt(RHO * R_D^3 / SIGMA)

const ETA_INF0 = 0.002747884967209792
const ETA_0    = 0.2601531742575956
const K_TRUE   = 1.076060672638565
const M_TRUE   = 0.6724627161400011
const OH_0     = ETA_0 / sqrt(RHO * SIGMA * R_D)

const M_RUN, K_RUN = 45, 3
const WE    = 1.0
const LOBS  = 2:8
const TGRID = collect(0.1:0.1:5.0)

eta_fn(k, m, ratio) = gd -> carreau(gd; lambda_c = k / T_CAP, a = m, n = 1 - m,
                                    eta_inf_ratio = ratio)

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

"""The surface observables: mode amplitudes and centre of mass on the fixed time grid."""
function observe(theta)
    k, m, ratio = 10^theta[1], theta[2], 10^theta[3]
    p = ImpactParams(We = WE, Bo = BO, Oh = OH_0, M = M_RUN, K = K_RUN,
                     eta = eta_fn(k, m, ratio), t_max = 5.2, stop_on_release = false)
    r = try DropSolver.simulate_lcp(p) catch; return nothing end
    length(r.t) < 20 && return nothing
    ls = DropSolver.basis(p).ls
    amps = [DropSolver.surface_amplitudes(p, a) for a in r.a]
    obs = Float64[]
    for l in LOBS
        j = findfirst(==(l), ls); j === nothing && return nothing
        append!(obs, resample(r.t, [a[j] for a in amps], TGRID))
    end
    append!(obs, resample(r.t, r.z, TGRID))
    all(isfinite, obs) ? obs : nothing
end

## The box. Centred on the truth, and wide enough to contain any parameter set a real
## experiment might plausibly confuse with it: roughly two standard deviations from the
## three-Weber study, scaled up for the single Weber number used here.
const KS  = collect(range(-0.90, 1.00; length = 5))     # log10 k
const MS  = collect(range( 0.45, 0.90; length = 5))     # m
const RS  = collect(range(-2.10, -1.85; length = 5))    # log10 eta_inf/eta_0

flushln("== nonlinear recovery, stage 1: the forward grid ==")
@printf("truth: log10 k = %.4f, m = %.4f, log10 ratio = %.4f\n",
        log10(K_TRUE), M_TRUE, log10(ETA_INF0 / ETA_0))
@printf("grid: %d x %d x %d = %d forward solves at We = %.1f, M = %d\n",
        length(KS), length(MS), length(RS), length(KS)*length(MS)*length(RS), WE, M_RUN)

t0 = time(); n = 0; nfail = 0
open(joinpath(OUT, "recovery_grid.csv"), "w") do io
    ## `do` opens a closure, so the counters have to be declared global or Julia treats
    ## them as fresh locals and the first read throws. This cost one silent 40-minute
    ## wait on a job that had already died.
    global n, nfail
    ## header: the parameter triple, then every observable
    println(io, "log10k,m,log10ratio," *
            join(("o$(i)" for i in 1:(length(LOBS)+1)*length(TGRID)), ","))
    for kk in KS, mm in MS, rr in RS
        n += 1
        o = observe([kk, mm, rr])
        if o === nothing
            nfail += 1
            @printf("  [%3d] log10k %+.2f m %.2f ratio %+.2f  FAILED\n", n, kk, mm, rr)
            flush(stdout); continue
        end
        println(io, join(vcat(kk, mm, rr, o), ",")); flush(io)
        n % 10 == 0 && (@printf("  %3d / %d done (%.1f min)\n", n,
                                length(KS)*length(MS)*length(RS), (time()-t0)/60);
                        flush(stdout))
    end
end
@printf("\n%d solves, %d failed, %.1f min\n", n, nfail, (time() - t0) / 60)
flushln("wrote ", joinpath(OUT, "recovery_grid.csv"))
