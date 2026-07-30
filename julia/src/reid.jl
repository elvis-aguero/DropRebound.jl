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
`|Im(q)|`, while the RATIO stays O(1). (Confirmed during development: a direct
`besselj`-based ratio threw an overflow exception at l=16, Oh=0.3, a perfectly
ordinary parameter pair.)
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
scales with `abs(q)`), which is not an exception but a silent multi-minute
hang -- confirmed during development. The step cap prevents `q` from ever
reaching such an argument in the first place.
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
    dominant_root(Oh, l; oh_start=1e-4, nsteps=24)

The dominant (least-damped) root of Reid's characteristic equation, found by
continuation in Oh: start at an Oh small enough that Lamb's asymptotic guess
is genuinely accurate, then walk geometrically up to the target, seeding each
step's Newton solve with the previous step's converged root. A single
Lamb-seeded Newton solve at the target Oh is NOT reliable in general -- it can
converge to a more strongly damped, non-dominant root once Oh grows (confirmed
at l=16, Oh=0.3 in `julia/derivations/reid_finite_oh_derivation.jl`).
"""
function dominant_root(Oh, l; oh_start=1e-4, nsteps=24)
    Oh <= oh_start && return _dominant_root_direct(Oh, l)
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
derivation script's Section 2). Found from an analytic small-q asymptotic
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
(`sigma_i = q_i^2*Oh`, this package's existing time-unit convention). `resid`
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

- `model = :lamb` -- `lambda_l = Oh*(l-1)*(2l+1)`, `omega_l^2 = l*(l-1)*(l+2)`.
  Bit-for-bit identical to what `build_residual!`/`build_jacobian!` have always
  computed inline.
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
            Reid continuation failed for some drop modes; fell back to Lamb there. Those
            modes are NOT using the arbitrary-Oh model -- worth investigating.
            """ modes = fallback Oh = Oh
    end
    lambda, omega2
end
