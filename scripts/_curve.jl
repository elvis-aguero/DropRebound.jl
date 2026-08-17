# Adaptive sampling of a curve that is expensive to evaluate.
#
# Shared by the figure scripts that draw restitution against Weber number, where a
# uniform grid is wrong twice over: it wastes points on the flat middle, and it
# under-resolves the roll-off at small Weber number, which is the part being argued
# about.
#
# THERE ARE THREE TOLERANCES AND NONE OF THEM IS A TARGET.
#
#   th     resolution in y. No two consecutive samples differ by more than this, so
#          a plotted segment is never a straight line standing in for a bend the
#          sampler never looked at.
#
#   x_max  coarseness in x, in decades. No plotted segment spans more than this,
#          however flat it is. Without it, |dy| alone is a linear test on a
#          logarithmic axis and a flat stretch is drawn as one straight line across
#          decades that were never computed.
#
#   x_th   the floor. An interval already this narrow is left alone however large
#          its jump: at that separation the jump is a real step in the curve, and
#          bisecting further buys nothing but simulations.
#
# Refinement ends when every interval is within `th` in y AND within `x_max` in x,
# or has reached `x_th`. None of the three is a target to be driven to.
#
# WHY NOT SIMPLY REFINE THE STEEPEST INTERVAL. That was the first version and it
# failed silently. On this curve the steepest descent in y is at large x, so the
# budget drained there and the low-x roll-off kept whatever the initial grid gave
# it -- three or four points across two decades, with the shape between them
# interpolated rather than computed. Bounding |dy| everywhere fixes it by
# construction, because an unresolved bend IS a large |dy|.
#
# WHY NOT DRIVE y DOWN TO th AS WELL. The second version did, extending the sweep
# to smaller x until restitution itself fell below the tolerance. That is only
# meaningful for a curve that rolls off, and these do not always: at Bond ~ 1e-4 a
# drop is effectively weightless and keeps rebounding however gently it lands. The
# sweep then runs decades past the smallest measurement chasing a floor that is not
# there -- 58 of 166 points, on one run, spent outside the data entirely.

"""
    adaptive_curve(ev, lo, hi; th, x_max, x_th, forced, n0, maxpts) -> (xs, ys)

Sample `ev` over `[lo, hi]`, logarithmically, refining any interval that spans more
than `th` in y OR more than `x_max` decades in x, and stopping once an interval is
narrower than `x_th` decades.

`forced` points are always included, so a caller needing values at particular
abscissae gets them without a second pass. `ev` may return `NaN` for a point that
did not converge; intervals with a `NaN` endpoint are never selected for
refinement, so one failure cannot absorb the whole budget.

Returns the samples in increasing `x`. `maxpts` is a backstop: if it binds, neither
tolerance was satisfied and `curve_report` says so.
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
                        n0::Int = 6, maxpts::Int = 60,
                        x_max::Float64 = 0.10, x_th::Float64 = 0.02,
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

    ## An interval needs refining if it is coarse in EITHER axis, and may be refined
    ## only while it is still wide enough in x to be worth it.
    ##
    ##   too coarse in y   |dy| > th
    ##   too coarse in x   dlog10(x) > x_max
    ##   refinable         dlog10(x) > x_th
    ##
    ## The x_max term is not decoration. Judging an interval by |dy| alone is a
    ## linear test on a logarithmic axis, and the two do not compose: where the curve
    ## is flat -- 0.98 to 0.95 across three decades of Weber -- the y test is
    ## satisfied at once and the segment is drawn as a straight line spanning decades
    ## with nothing computed inside it. Flat in y is not the same as resolved.
    ##
    ## Priority goes to whichever interval is worst relative to its own tolerance, so
    ## the two criteria compete on equal terms instead of one always winning.
    while length(xs) < maxpts
        worst, at = -Inf, 0
        for i in 1:length(xs)-1
            dx = log10(xs[i+1] / xs[i])
            dx > x_th || continue
            dy = abs(ys[i+1] - ys[i])
            isfinite(dy) || (dy = 0.0)
            score = max(dy / th, dx / x_max)
            score > worst && (worst = score; at = i)
        end
        (at == 0 || worst <= 1.0) && break
        xm = sqrt(xs[at] * xs[at+1])        # geometric midpoint: halves the log interval
        insert!(xs, at + 1, xm)
        insert!(ys, at + 1, ev(xm))
    end
    xs, ys
end

"""
    jumps(xs, ys, th; x_th) -> Vector{Tuple{Float64,Float64,Float64}}

Intervals that are narrow in `x` and still span more than `th` in `y`: the places
the curve genuinely steps rather than bends. Returned as `(x_lo, x_hi, dy)`.
"""
function jumps(xs, ys, th; x_th::Float64 = 0.02)
    out = Tuple{Float64,Float64,Float64}[]
    for i in 1:length(xs)-1
        d = abs(ys[i+1] - ys[i])
        (isfinite(d) && d > th && log10(xs[i+1] / xs[i]) <= x_th) || continue
        push!(out, (xs[i], xs[i+1], d))
    end
    out
end

"""
    trim_ceiling_jump(xs, ys, ceiling, th; x_th = 0.02) -> (xs, ys)

Drop the last sample if it sits at `ceiling` and is a genuine jump from its
neighbour, by the same test `jumps` uses: narrower than `x_th` decades in `x`,
larger than `th` in `y`.

A caller passes `ceiling` because it does not trust the curve above that point --
`adaptive_curve` is called with `hi = ceiling` for exactly that reason. But the
last sample is still evaluated AT the ceiling, and if the underlying model is
already breaking down there, that sample can land as a step rather than a bend:
too narrow in `x` for the sampler to have refined further, too large in `y` to
be a smooth continuation of the curve up to that point. Plotting it anyway would
show a real discontinuity as though it were a resolved result, in a regime the
caller has already disclaimed. Leaves `xs`, `ys` untouched if the last point is
not at `ceiling`, or is a smooth continuation rather than a step.
"""
function trim_ceiling_jump(xs, ys, ceiling, th; x_th::Float64 = 0.02)
    length(xs) >= 2 && isapprox(xs[end], ceiling; rtol = 1e-9) || return xs, ys
    js = jumps(xs, ys, th; x_th = x_th)
    (!isempty(js) && isapprox(js[end][2], xs[end]; rtol = 1e-9)) || return xs, ys
    xs[1:end-1], ys[1:end-1]
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
    cor_or_settle(primary, fallback, t_max) -> Float64

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
function cor_or_settle(primary, fallback, t_max)
    r = primary()
    isfinite(r.cor) && return r.cor
    ## Not a settle at all -- the march died some other way.
    (isfinite(r.tc) && r.tc >= 0.9 * t_max) || return NaN

    ## A march that runs to t_max still in contact means one of two very different
    ## things, and they must not be conflated. At small Weber number the drop has
    ## genuinely settled and its restitution is zero. At large Weber number the
    ## active-set closure is documented to fail in exactly this way -- it walks to
    ## the contact set from the previous step and stalls when the two are far apart,
    ## which is what a fast, weakly damped impact produces -- and calling that a
    ## settled drop would put a hard zero in the middle of a healthy curve.
    ##
    ## So the closure that does not have that failure mode is asked before the answer
    ## is called physics. If it also never releases, the drop really does stay down.
    r2 = fallback()
    isfinite(r2.cor) && return r2.cor
    (isfinite(r2.tc) && r2.tc >= 0.9 * t_max) && return 0.0
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
    nj = length(jumps(xs, ys, th))
    string(length(xs), " pts, max |Δ| = ", isnan(worst) ? "n/a" : string(round(worst, digits = 4)),
           worst <= th ? " (resolved in y)" :
               nj > 0 ? " (resolved in x;  step" * (nj == 1 ? ")" : "s)") : " (NOT resolved)")
end
