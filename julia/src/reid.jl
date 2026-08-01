"""
Reid (1960) exact linear viscous drop oscillation, valid at arbitrary Ohnesorge
number, replacing Lamb's small-viscosity asymptotics used everywhere else in
this package. Derivation and cross-checks against two independently-published
limits: `julia/derivations/reid_finite_oh_derivation.jl`.

The per-mode ODE this package integrates (build_residual!/build_jacobian!) is
the unit-mass oscillator `Addot_l + 2*lambda_l*Adot_l + omega_l^2*A_l = forcing`.
Lamb's formula is `lambda_l = Oh*(l-1)*(2l+1)`, `omega_l^2 = l*(l-1)*(l+2)` --
the exact Oh -> 0 limit of Reid's characteristic equation, with error growing
in both Oh and l. `drop_viscous_coeffs` computes either that closed form
(`model = :lamb`, bit-for-bit identical to what this package always did) or
Reid's exact finite-Oh values (`model = :reid`) via Vieta's formulas on the two
dominant eigenvalues of Reid's own characteristic equation.
"""

"""
    sph_bessel_ratio(l, q) -> Complex

`Q_l(q) = j_{l+1}(q)/j_l(q)` by downward recurrence from the three-term
spherical Bessel recurrence `j_{l-1}+j_{l+1} = ((2l+1)/q) j_l`, divided by `j_l`
and inverted. Never evaluates a Bessel function directly -- evaluating
`j_l(q)` itself overflows at small Oh, where `q = sqrt(sigma/Oh)` has large
`|Im(q)|`, while the ratio stays O(1). A direct `besselj`-based ratio overflows
at ordinary parameters -- for example l=16, Oh=0.3.
"""
function sph_bessel_ratio(l::Integer, q::Number)
    # Capped defensively: the recurrence's working length scales with abs(q),
    # so an ill-behaved caller passing a huge |q| (e.g. a Newton step that
    # overshot before step-capping was added) would otherwise turn into a
    # silent multi-minute hang rather than a fast, inaccurate answer that a
    # residual check downstream can reject.
    pad = min(max(60, l ÷ 2 + ceil(Int, abs(q))), 5000)
    n0 = l + pad + min(ceil(Int, abs(q)), 5000)
    Q = q / (2 * n0 + 3)
    for n in n0:-1:(l+1)
        Q = 1 / ((2n + 1) / q - Q)
    end
    return Q
end

"""
    reid_char(q, Oh, l)

Residual of Reid's exact characteristic equation for mode `l`:
`alpha^4/q^4 + 1 = (2(l-1)/q^2) * [l + (l+1)*(q - 2l*Q)/(q - 2Q)]`, with
`alpha^2 = sqrt(l(l-1)(l+2))/Oh` and `Q(q) = j_{l+1}(q)/j_l(q)`. Zero of this
function iff `q` is an eigenvalue wavenumber of the linearized viscous drop
problem at this `(Oh, l)`.
"""
function reid_char(q, Oh, l)
    alpha2 = sqrt(l * (l - 1) * (l + 2)) / Oh
    Q = sph_bessel_ratio(l, q)
    lhs = alpha2^2 / q^4 + 1
    rhs = 2 * (l - 1) / q^2 * (l + (l + 1) * (q - 2l * Q) / (q - 2Q))
    lhs - rhs
end

function _safe_eval(f, q)
    try
        v = f(q)
        return isfinite(abs(v)) ? v : nothing
    catch
        return nothing
    end
end

"""
    _newton_complex(f, q0; ...)

Damped, step-capped Newton with backtracking: each step is capped at a
fraction of `|q|` and halved until the residual decreases. An UNDAMPED Newton
step can overshoot into an argument where `sph_bessel_ratio`'s downward
recurrence needs an enormous number of terms to converge (its working length
scales with `abs(q)`), which produces no exception -- only a multi-minute hang.
The step cap keeps `q` away from such arguments.
"""
function _newton_complex(f, q0; maxiter=300, tol=1e-13, step_cap_frac=0.5)
    q = q0
    F = _safe_eval(f, q)
    F === nothing && return q0
    for _ in 1:maxiter
        h = max(1e-8, abs(q) * 1e-8)
        Fh = _safe_eval(f, q + h * im)
        Fh === nothing && break
        dF = (Fh - F) / (h * im)
        abs(dF) < 1e-300 && break
        step = F / dF
        cap = step_cap_frac * max(abs(q), 1e-3)
        abs(step) > cap && (step *= cap / abs(step))
        accepted = false
        t = 1.0
        for _ in 1:40
            cand = q - t * step
            Fc = _safe_eval(f, cand)
            if Fc !== nothing && abs(Fc) < abs(F)
                q, F, accepted = cand, Fc, true
                break
            end
            t /= 2
        end
        accepted || break
        abs(step) < tol * max(abs(q), 1.0) && break
    end
    q
end

_dominant_root_direct(Oh, l) = begin
    sigma0 = sqrt(l * (l - 1) * (l + 2))
    gamma0 = (l - 1) * (2l + 1) * Oh
    q0 = sqrt(complex(gamma0 / Oh, -sigma0 / Oh))
    imag(q0) > 0 && (q0 = -q0)
    q = _newton_complex(qv -> reid_char(qv, Oh, l), q0)
    imag(q) > 0 ? conj(q) : q
end

"""
    dominant_root(Oh, l; oh_start=1e-4, max_step_ratio=1.15)

The dominant (least-damped) root of Reid's characteristic equation, found by
continuation in Oh: start at an Oh small enough that Lamb's asymptotic guess
is genuinely accurate, then walk geometrically up to the target, seeding each
step's Newton solve with the previous step's converged root. A single
Lamb-seeded Newton solve at the target Oh is NOT reliable in general -- it can
converge to a more strongly damped, non-dominant root once Oh grows (confirmed
at l=16, Oh=0.3 in the "Finite-Ohnesorge Coefficients" derivation).

The number of steps is chosen ADAPTIVELY so each step's Oh ratio stays below
`max_step_ratio`, rather than a fixed step count. A fixed count gives a
per-step ratio too coarse to track the dominant branch at higher `l`: with 24
steps from Oh=1e-4, the continuation jumps to a more strongly damped root past
Oh ~ 900 for l=10, and at smaller Oh for larger `l` (for instance lambda
61483 -> 273551 between Oh=958 and Oh=1014, both satisfying the characteristic
equation to near machine precision -- two genuinely different roots, not a
convergence failure). A ratio of ~1.15 is stable to within 15% step-to-step
for `l` up to 40 across Oh 1e-4 to 1e4.
"""
function dominant_root(Oh, l; oh_start=1e-4, max_step_ratio=1.15)
    Oh <= oh_start && return _dominant_root_direct(Oh, l)
    log_ratio = log(Oh / oh_start)
    nsteps = max(1, ceil(Int, log_ratio / log(max_step_ratio)))
    q = _dominant_root_direct(oh_start, l)
    for ohv in exp.(range(log(oh_start), log(Oh); length=nsteps + 1))[2:end]
        q = _newton_complex(qv -> reid_char(qv, ohv, l), q)
    end
    q
end

"""
    second_root(Oh, l, q1)

The second dominant root, given the first (`q1`, from `dominant_root`).
Underdamped (`q1` complex): the true conjugate pair, no separate search.
Overdamped (`q1` real): a genuinely separate, much smaller real root (a
"creep" mode) -- NOT obtainable by continuing `conj(q1)`'s branch through the
critical Oh, which collapses onto the same root as `q1` instead (see the
"Finite-Ohnesorge Coefficients" derivation, Section 2). Found from an analytic small-q asymptotic
guess (derived by substituting the small-argument ratio `Q(q) ~ q/(2l+3)` into
`reid_char` and balancing its two singular terms), refined by Newton.
"""
function second_root(Oh, l, q1)
    if abs(imag(q1)) > 1e-9 * abs(q1)
        return conj(q1)
    end
    alpha2 = sqrt(l * (l - 1) * (l + 2)) / Oh
    q0 = alpha2 * sqrt((2l + 1) / (2 * (l - 1) * (2l^2 + 4l + 3)))
    _newton_complex(qv -> reid_char(qv, Oh, l), complex(q0); maxiter=200, tol=1e-13)
end

"""
    reid_lambda_omega2(Oh, l) -> (lambda, omega2, resid)

Exact finite-Oh coefficients of the unit-mass oscillator whose eigenvalues are
Reid's own two dominant roots, via Vieta's formulas:
`lambda = (sigma_1+sigma_2)/2`, `omega2 = sigma_1*sigma_2`
(`sigma_i = q_i^2*Oh`, this package's time-unit convention). `resid`
is the achieved characteristic-equation residual at the dominant root, for the
caller to assert on.
"""
function reid_lambda_omega2(Oh, l)
    q1 = dominant_root(Oh, l)
    q2 = second_root(Oh, l, q1)
    s1 = q1^2 * Oh
    s2 = q2^2 * Oh
    lambda = real(s1 + s2) / 2
    omega2 = real(s1 * s2)
    resid = abs(reid_char(q1, Oh, l))
    lambda, omega2, resid
end

"""
    drop_viscous_coeffs(M, Oh, model::Symbol) -> (lambda::Vector, omega2::Vector)

Per-mode damping and squared frequency for modes `l = 2..M`, returned as
length-`(M-1)` vectors indexed the same way as `A[2:end]`/`Adot[2:end]`
throughout this package (`lambda[k]`, `omega2[k]` correspond to mode `l=k+1`).

- `model = :lamb` -- Lamb's small-Ohnesorge asymptotics,
  `lambda_l = Oh*(l-1)*(2l+1)`, `omega_l^2 = l*(l-1)*(l+2)`; the exact `Oh -> 0`
  limit of the `:reid` result, with error growing in both `Oh` and `l`.
- `model = :reid` -- exact roots of Reid's characteristic equation. Falls back
  to Lamb for any mode whose root fails to converge (residual >= 1e-6),
  warning which modes and why, rather than silently using a wrong coefficient.
"""
function drop_viscous_coeffs(M::Int, Oh::Float64, model::Symbol)
    model in (:lamb, :reid) ||
        throw(ArgumentError("viscous model must be :lamb or :reid, got $model"))
    ns = 2:M
    lambda_lamb = [Oh * (n - 1) * (2n + 1) for n in ns]
    omega2_lamb = [Float64(n) * (n - 1) * (n + 2) for n in ns]
    model === :lamb && return lambda_lamb, omega2_lamb

    lambda = copy(lambda_lamb)
    omega2 = copy(omega2_lamb)
    fallback = Int[]
    for (i, n) in enumerate(ns)
        lam, om2, resid = reid_lambda_omega2(Oh, n)
        if !isfinite(lam) || !isfinite(om2) || resid >= 1e-6 || lam <= 0
            push!(fallback, n)
            continue
        end
        lambda[i], omega2[i] = lam, om2
    end
    if !isempty(fallback)
        @warn """
            Reid continuation failed for some drop modes; those modes fall back to
            Lamb's small-viscosity formula and are not using the arbitrary-Oh model.
            """ modes = fallback Oh = Oh
    end
    lambda, omega2
end

# ------------------------------------------------------------------------------
# Fast, cached lookup: reid_lambda_omega2 costs ~2ms/call (a ~24-step Oh
# continuation, each step a damped Newton solve). That's fine once per
# SimConstants, but the non-perturbative Carreau-Yasuda extension needs
# lambda_l(Oh_eff_l(t)), omega_l^2(Oh_eff_l(t)) at a DIFFERENT Oh every
# residual/Jacobian evaluation -- thousands of calls per simulation, which
# would cost tens of seconds to minutes per run at ~2ms each. A precomputed
# table + interpolation reduces each lookup to a few array accesses and a
# linear interpolation, at a one-time setup cost of (table size) x ~2ms.
# ------------------------------------------------------------------------------

"""
Precomputed `lambda_l(Oh)`, `omega_l^2(Oh)` on a log-spaced Oh grid for one
mode `l`, for fast interpolated lookup. Built once (`build_reid_table`), then
queried many times per simulation (`reid_lambda_omega2_fast`).
"""
struct ReidTable
    l      :: Int
    log_Oh :: Vector{Float64}
    lambda :: Vector{Float64}
    omega2 :: Vector{Float64}
end

"""
    build_reid_table(l; Oh_min=1e-6, Oh_max=1e6, n=100) -> ReidTable

Tabulate `reid_lambda_omega2(Oh, l)` at `n` log-spaced points from `Oh_min` to
`Oh_max`. The default range spans ~12 orders of magnitude -- far wider than
any physical rest-to-fully-thinned Ohnesorge swing this package's shear-
thinning models are expected to produce (the real validation fluid this was
built for spans Oh ~ 0.025 to 57) -- so lookups landing outside the table
should not occur in practice; `reid_lambda_omega2_fast` clamps to the nearest
edge if one ever does, rather than extrapolating into a regime with no
guarantee of accuracy.
"""
function build_reid_table(l::Int; Oh_min::Float64=1e-6, Oh_max::Float64=1e6, n::Int=100)
    log_Oh = collect(range(log(Oh_min), log(Oh_max); length=n))
    lambda = Vector{Float64}(undef, n)
    omega2 = Vector{Float64}(undef, n)

    # Track the dominant root by SEQUENTIAL continuation across the table's
    # own (fine) grid, rather than calling reid_lambda_omega2/dominant_root
    # independently at each point. This matters: dominant_root's internal
    # continuation uses a FIXED 24-step ladder from oh_start=1e-4 to whatever
    # target Oh it's given, which can be too coarse (and jump to a more
    # strongly damped, non-dominant branch) at higher l and Oh -- confirmed
    # during development: reid_lambda_omega2(Oh, 10) jumps discontinuously
    # (e.g. lambda 61483 -> 273551 between Oh=958 and Oh=1014, both with
    # tiny char-eq residuals -- i.e. two DIFFERENT genuine roots, not a
    # convergence failure) somewhere past Oh~900 for l=10, and at
    # correspondingly smaller Oh for higher l. Seeding each grid point's
    # Newton solve from the PREVIOUS (adjacent, already-tracked) grid point
    # keeps every step small relative to the table's own fine spacing,
    # avoiding this failure mode regardless of l or how wide a range is
    # requested.
    q1 = _dominant_root_direct(exp(log_Oh[1]), l)
    prev_lambda = NaN
    for (i, lo) in enumerate(log_Oh)
        Oh = exp(lo)
        q1 = _newton_complex(qv -> reid_char(qv, Oh, l), q1)
        q2 = second_root(Oh, l, q1)
        s1 = q1^2 * Oh
        s2 = q2^2 * Oh
        lam = real(s1 + s2) / 2
        om2 = real(s1 * s2)
        resid = abs(reid_char(q1, Oh, l))
        if !isfinite(lam) || !isfinite(om2) || resid >= 1e-6 || lam <= 0
            lam = Oh * (l - 1) * (2l + 1)              # Lamb fallback at this grid point
            om2 = Float64(l * (l - 1) * (l + 2))
            q1 = _dominant_root_direct(Oh, l)            # resync the tracked branch
        elseif isfinite(prev_lambda) && prev_lambda > 0 &&
               !(0.2 < lam / prev_lambda < 5.0)
            # A >5x jump between adjacent (fine-grid) points despite a tiny
            # residual means continuation still landed on a different root
            # than its neighbor -- fall back to Lamb here too rather than
            # bake a discontinuity into the table silently.
            lam = Oh * (l - 1) * (2l + 1)
            om2 = Float64(l * (l - 1) * (l + 2))
            q1 = _dominant_root_direct(Oh, l)
        end
        lambda[i], omega2[i] = lam, om2
        prev_lambda = lam
    end
    ReidTable(l, log_Oh, lambda, omega2)
end

"""
    reid_lambda_omega2_fast(table, Oh) -> (lambda, omega2)

Linear interpolation (in log(Oh)) of a `ReidTable` built by `build_reid_table`.
Clamps to the table's nearest edge for `Oh` outside its range.
"""
function reid_lambda_omega2_fast(table::ReidTable, Oh::Float64)
    logOh = log(Oh)
    g = table.log_Oh
    if logOh <= g[1]
        return table.lambda[1], table.omega2[1]
    elseif logOh >= g[end]
        return table.lambda[end], table.omega2[end]
    end
    idx = clamp(searchsortedlast(g, logOh), 1, length(g) - 1)
    t = (logOh - g[idx]) / (g[idx+1] - g[idx])
    lambda = table.lambda[idx] + t * (table.lambda[idx+1] - table.lambda[idx])
    omega2 = table.omega2[idx] + t * (table.omega2[idx+1] - table.omega2[idx])
    lambda, omega2
end

"""
    build_reid_cache(M; kwargs...) -> Vector{ReidTable}

`ReidTable`s for modes `l = 2..M`, indexed the same way as
`drop_viscous_coeffs`'s output (`cache[k]` corresponds to mode `l=k+1`).
Build once per simulation (or reuse across many, since the table depends
only on `M` and the Oh range, not on any particular fluid's shear-thinning
parameters) and pass to `reid_lambda_omega2_fast`.
"""
build_reid_cache(M::Int; kwargs...) = [build_reid_table(l; kwargs...) for l in 2:M]
