#!/usr/bin/env julia
# dt-convergence study at M=100 (finer than the paper's L=90).
# For each case, refine dt_max by 1, 1/2, 1/4 and report:
#   (a) old-vs-new COM-trajectory deviation at each dt level — shrinking ⇒ the
#       old/new difference is pure discretization (both converge to one answer);
#   (b) self-convergence of each stepper (successive dt levels) — Cauchy check.
# Deviations are max|Δz| over the overlapping time span, as a fraction of R.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

const M   = 100
const Oh  = 0.1
const Bo  = 0.01
const V0  = -0.5
const TE  = 9.0
const SE  = 0.02
const FACTORS = [1.0, 0.5, 0.25]

const DTM0 = make_dt_max(M)
const TV   = make_theta_vec(M)
const PC   = precompute_integrals(NaN, M)[1]

function interp(xs, ys, q)
    out = similar(q)
    for (k, x) in enumerate(q)
        if x <= xs[1]; out[k] = ys[1]
        elseif x >= xs[end]; out[k] = ys[end]
        else
            i = searchsortedlast(xs, x); t = (x - xs[i]) / (xs[i+1] - xs[i])
            out[k] = ys[i] + t * (ys[i+1] - ys[i])
        end
    end
    out
end

function traj(θe, ξ, ev, dtfac)
    cfg = SimConstants(M, M+1, Oh, Bo, TV, PC, DTM0 * dtfac)
    i = DropState(M); i.z = 1.05; i.v = V0; i.dt = DTM0 * dtfac; i.cp = 0
    try
        t, s = solve_drop!(cfg, OBParams(), i; cl = CLParams(θe, ξ),
                           t_end = TE, save_every = SE, event_location = ev)
        return (t, [x.z for x in s])
    catch
        return nothing
    end
end

devmax(a, b) = begin
    (a === nothing || b === nothing) && return NaN
    lo = max(a[1][1], b[1][1]); hi = min(a[1][end], b[1][end])
    g  = collect(range(lo, hi; length = 800))
    maximum(abs.(interp(b[1], b[2], g) .- interp(a[1], a[2], g)))
end

for (θe, ξ) in [(0.92π, 0.0), (0.92π, 1.0)]
    @printf("\n=== case θ_e=%.2fπ, ξ=%.1f  (M=%d) ===\n", θe/π, ξ, M)
    olds = Dict(); news = Dict()
    for f in FACTORS
        olds[f] = traj(θe, ξ, false, f)
        news[f] = traj(θe, ξ, true,  f)
        @printf("  dt_max×%.2f : old-vs-new max|Δz|/R = %.4f\n", f, devmax(olds[f], news[f]))
    end
    @printf("  self-convergence (max|Δz|/R between successive dt levels):\n")
    for (fa, fb) in [(1.0, 0.5), (0.5, 0.25)]
        @printf("    legacy  dt×%.2f→%.2f : %.4f\n", fa, fb, devmax(olds[fa], olds[fb]))
        @printf("    event   dt×%.2f→%.2f : %.4f\n", fa, fb, devmax(news[fa], news[fb]))
    end
end
