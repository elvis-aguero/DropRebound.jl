# Adaptive sampling of a curve that is expensive to evaluate.
#
# Shared by the figure scripts that draw restitution against Weber number, where a
# uniform grid is wrong twice over: it wastes points on the flat middle, and it
# under-resolves the roll-off at small Weber number, which is the part being argued
# about.
#
# THE SAMPLER GUARANTEES TWO THINGS, and they are the reason it exists.
#
#   1. RESOLUTION.  No two consecutive samples differ by more than `th` in y. A
#      plotted segment is then never a straight line standing in for a bend the
#      sampler never looked at.
#
#   2. REACH.  The sweep is extended to smaller x until y itself falls below `th`.
#      A curve that rolls off has to be followed until it has actually rolled off;
#      stopping while y is still large reports the beginning of a decline as though
#      it were the whole of one.
#
# The same `th` sets both, which is the point: the smallest feature the curve is
# resolved to is also the smallest value it is followed down to.
#
# WHY NOT SIMPLY REFINE THE STEEPEST INTERVAL. That was the first version and it
# failed silently. On this curve the steepest descent in y is at large x, so the
# budget drained there and the low-x roll-off kept whatever the initial grid gave
# it -- three or four points across two decades, with the shape between them
# interpolated rather than computed. Bounding |dy| everywhere fixes it by
# construction, because an unresolved bend IS a large |dy|.

"""
    adaptive_curve(ev, lo, hi; th, forced, n0, maxpts, xfloor) -> (xs, ys)

Sample `ev` over `[lo, hi]`, logarithmically, until no two consecutive samples
differ by more than `th`, extending below `lo` until `ev` itself drops below `th`
or `xfloor` is reached.

`forced` points are always included, so a caller needing values at particular
abscissae gets them without a second pass. `ev` may return `NaN` for a point that
did not converge; intervals with a `NaN` endpoint are never selected for
refinement, so one failure cannot absorb the whole budget, and the descent stops
rather than chasing a run that no longer completes.

Returns the samples in increasing `x`. `maxpts` is a backstop, not a target: if it
binds, the two guarantees above have NOT been met, and the caller is told so by
`curve_report`.
"""
## Progress goes to stderr, line by line, under a lock and explicitly flushed.
##
## Redirecting stdout to a file makes it block-buffered, so a long sweep that
## prints only when it finishes looks identical to one that has hung -- and once
## the curves were computed concurrently, ahead of all the printing, there was no
## output at all until the end. Reporting each evaluation as it lands is what makes
## "how far along is it" answerable, so it is not optional decoration.
const _PROGRESS_LOCK = ReentrantLock()
const _PROGRESS_T0 = Ref(0.0)

progress_reset() = (_PROGRESS_T0[] = time())

function _progress(tag, n, x, y, dt)
    lock(_PROGRESS_LOCK) do
        el = _PROGRESS_T0[] == 0.0 ? 0.0 : time() - _PROGRESS_T0[]
        println(stderr, "    [", lpad(round(Int, el), 5), "s] ", tag,
                " pt ", lpad(n, 2), "  We = ", rpad(round(x, sigdigits = 3), 9),
                " cor = ", rpad(isfinite(y) ? round(y, digits = 4) : "NaN", 8),
                " (", round(dt, digits = 1), "s)")
        flush(stderr)
    end
end

function adaptive_curve(ev, lo, hi;
                        th::Float64 = 0.01,
                        forced::AbstractVector{<:Real} = Float64[],
                        n0::Int = 6, maxpts::Int = 60, xfloor::Float64 = 1e-7,
                        tag::AbstractString = "")
    raw, count = ev, Ref(0)
    ev = function (x)
        t0 = time()
        y = raw(x)
        _progress(tag, (count[] += 1), x, y, time() - t0)
        y
    end
    xs = collect(exp10.(range(log10(lo), log10(hi); length = n0)))
    for f in forced
        lo <= f <= hi || continue
        any(x -> isapprox(x, f; rtol = 1e-6), xs) || push!(xs, float(f))
    end
    sort!(xs)
    ys = Float64[ev(x) for x in xs]

    ## REACH: walk down in x until the curve has actually collapsed below `th`.
    ## A decade per step is too coarse to keep |dy| bounded on its own, but the
    ## refinement pass below cleans up whatever this leaves behind.
    while length(xs) < maxpts && isfinite(ys[1]) && ys[1] > th && xs[1] > xfloor
        xn = max(xs[1] / 3.0, xfloor)
        pushfirst!(xs, xn)
        pushfirst!(ys, ev(xn))
    end

    ## RESOLUTION: bisect the widest jump in y until none exceeds `th`.
    while length(xs) < maxpts
        worst, at = -Inf, 0
        for i in 1:length(xs)-1
            d = abs(ys[i+1] - ys[i])
            isfinite(d) || continue
            d > worst && (worst = d; at = i)
        end
        (at == 0 || worst <= th) && break
        xm = sqrt(xs[at] * xs[at+1])        # geometric midpoint: halves the log interval
        insert!(xs, at + 1, xm)
        insert!(ys, at + 1, ev(xm))
    end
    xs, ys
end

"""
    parallel_curves(f, items) -> Vector

Evaluate `f` on each item concurrently, returning results in the order given.

Curves are independent, so this is the free axis of parallelism: one thread per
curve, no interaction. Within a curve the sampler is inherently sequential -- each
new abscissa is chosen from the results so far -- so that is left alone.

Measured on this solver: about 2.4x on six physical cores, not six. The limit is
not scheduling. One impact allocates on the order of a gigabyte, so several at once
saturate memory bandwidth; separate processes scale *worse* (1.9x) because they
lose the shared cache without fixing the allocation. Reducing the solver's
allocation is the change that would make this scale, not a better scheduler.

Run with `julia -t auto` (or `-t 6`); with one thread this is an ordinary map.
"""
function parallel_curves(f, items)
    out = Vector{Any}(undef, length(items))
    Threads.@threads :dynamic for i in eachindex(items)
        out[i] = f(items[i])
    end
    out
end

"""
    cor_or_settle(r, t_max) -> Float64

Restitution of a finished run, with a drop that never left the substrate recorded
as zero rather than as a failure.

This distinction is the whole low-Weber story. `simulate` reports a run that stays
in contact for the whole march as `cor = NaN`, which is right when the question is
"did this run converge" and wrong when the question is "how well does the drop
bounce": a drop that settles has a restitution, and it is zero. Left as NaN it
would be dropped from the curve, and a family of curves would appear to stop for
numerical reasons exactly where the physical answer became most interesting.

A genuine numerical failure -- non-finite metrics with a contact time that is not
simply the whole march -- still returns NaN.
"""
function cor_or_settle(r, t_max)
    isfinite(r.cor) && return r.cor
    isfinite(r.tc) && r.tc >= 0.9 * t_max && return 0.0
    NaN
end

"""
    curve_report(xs, ys, th) -> String

One line stating whether the two guarantees hold, for printing next to the curve.
A sampler that quietly failed to meet them would make the figure look better than
the computation behind it.
"""
function curve_report(xs, ys, th)
    fin = isfinite.(ys)
    gaps = [abs(ys[i+1] - ys[i]) for i in 1:length(ys)-1 if fin[i] && fin[i+1]]
    worst = isempty(gaps) ? NaN : maximum(gaps)
    reached = any(fin) && minimum(ys[fin]) <= th
    string(length(xs), " pts, max |Δ| = ", isnan(worst) ? "n/a" : string(round(worst, digits = 4)),
           worst <= th ? " (resolved)" : " (NOT resolved)",
           reached ? ", rolled off" : ", NOT rolled off")
end
