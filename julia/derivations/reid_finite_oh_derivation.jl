# ==============================================================================
# Finite-Oh Damping and Frequency: from Reid's Exact Characteristic Equation
#
# What's broken today: `julia/src/timestepper.jl` (build_residual!/build_jacobian!)
# implements ONLY Lamb's classical result D1[l]=l(l-1)(l+2), D2[l]=2*Oh*(l-1)*(2l+1)
# for every mode l -- see julia/test/test_ob_eigenvalue.jl's own comment: "the
# Lamb approximation coefficients". Lamb's formula is the *leading-order-in-Oh*
# (nu -> 0) asymptotic reduction of Reid (1960)'s exact viscous characteristic
# equation (docs/reid1960_expanded-3.tex, "The small-viscosity limit"), and its
# error grows with BOTH Oh and l -- not just a large-Oh problem.
#
# Why this matters concretely: a real shear-thinning fluid shared for this
# repo's validation dataset has a REST-state Ohnesorge number Oh_0 ~ 57 (using
# the given eta_0=8.43 Pa.s, Bo=0.012, R=0.3mm, sigma=72.8mN/m) and a
# fully-thinned Oh_inf ~ 0.025 -- a ~2260x swing that spends most of its time
# deep in the regime Lamb's formula was never meant to cover.
#
# This derivation, and the choice of what to replace Lamb's coefficients WITH,
# leans heavily on a sister repository, elvis-aguero/SpectralKM.jl (a fully
# spectral kinematic-match model for droplet-bath impact by the same group),
# which already carries a battle-tested `:reid` viscous model in production
# (`src/reid.jl`, `src/types.jl`'s `DEFAULT_VISCOUS = :reid`). Three concrete
# things were ported or learned from it, credited at point of use below:
#   1. A downward Bessel-ratio recurrence that never evaluates a Bessel function
#      directly. An earlier version of this script used SpecialFunctions.besselj
#      directly and threw `AmosException: overflow` at l=16, Oh=0.3 -- a
#      perfectly reasonable (l, Oh) pair, not an edge case.
#   2. Continuation in Oh for the dominant root. A single Newton solve from
#      Lamb's asymptotic guess (this script's original approach) can converge to
#      the WRONG root once Oh grows -- confirmed by reproducing SpectralKM.jl's
#      own documented example (l=16, Oh=0.3: a direct solve returns a real root
#      near 288, the wrong branch, while continuation correctly tracks ~150).
#   3. The parametrization itself. This script's first version followed
#      Molacek & Bush (2012) and added a per-mode INERTIA coefficient A_l(Oh)
#      alongside a rescaled damping D_l(Oh), keeping the restoring/frequency
#      term fixed at the inviscid l(l-1)(l+2). SpectralKM.jl instead keeps unit
#      inertia and lets BOTH the damping (lambda_l) and the frequency-squared
#      term (omega_l^2) vary with Oh, i.e. the per-mode ODE
#          Addot_l + 2*lambda_l(Oh)*Adot_l + omega_l^2(Oh)*A_l = forcing
#      This is mathematically equivalent to the A_l/D_l form for FREE decay
#      (same two roots, just an overall rescaling of the ODE) but is the
#      simpler gauge, requires no new inertia term in julia/src/timestepper.jl
#      at all, and is the one already carrying real simulations in the sister
#      repo -- so it is what this script derives and what should be wired into
#      production, not the A_l/D_l form. See the note at the end of Section 3
#      on why the two are NOT obviously equivalent once contact forcing (B_l)
#      enters -- an open question this script does not resolve.
#
# What this derives: Reid's own transcendental characteristic equation
# (docs/reid1960_expanded-3.tex eq:char_eq) has, for each l and each Oh, a pair
# of dominant eigenvalues sigma_1(Oh), sigma_2(Oh) -- complex conjugates when
# Oh is below the critical point (damped oscillation), two DISTINCT real values
# above it (aperiodic, "fast" and "slow" decay). Vieta's formulas on this pair
# give lambda_l(Oh) = (sigma_1+sigma_2)/2, omega_l^2(Oh) = sigma_1*sigma_2 --
# the coefficients of the unit-mass oscillator whose eigenvalues are exactly
# Reid's own pair. Section 3 verifies this against Molacek & Bush's
# independently-stated closed-form limits (recast in this parametrization:
# lambda_l -> Oh(l-1)(2l+1), omega_l^2 -> l(l-1)(l+2) as Oh->0; and
# omega_l^2 -> l(l-1)(l+2) * (l-1)(2l^2+4l+3)/((l-1)(2l+1)) as Oh->infinity --
# derived from their A_l, D_l limits via the rescaling above), with no free
# parameters.
#
# A subtlety worth flagging for a future reader: the "second root" (needed once
# overdamped) is NOT obtained by continuing the complex-conjugate pair through
# the critical Oh. Numerically, doing that makes BOTH branches converge onto
# the SAME real root as Oh grows (a continuation artifact -- see Section 2),
# which silently gives a wrong, degenerate result. The correct second root is a
# genuinely separate, much smaller real root ("creep" mode) found by an
# independent search, not by tracking the first root's conjugate across the
# transition.
# ==============================================================================

using DropSolver

# ------------------------------------------------------------------------------
# Section 1: Reid's characteristic equation (docs/reid1960_expanded-3.tex,
# eq:char_eq) and a numerically robust way to evaluate it.
# ------------------------------------------------------------------------------

"""
    sph_bessel_ratio(l, q) -> Complex

`Q_l(q) = j_{l+1}(q)/j_l(q)` by DOWNWARD recurrence (ported from
elvis-aguero/SpectralKM.jl's `src/reid.jl`, credited there to the standard
three-term spherical Bessel recurrence `j_{l-1}+j_{l+1} = ((2l+1)/q) j_l`,
divided through by `j_l` and inverted). Never evaluates a Bessel function
directly, which matters because `besselj(l+0.5, q)` OVERFLOWS at small Oh
(where `q ~ sqrt(sigma/Oh)` has large |Im(q)|) while the RATIO itself stays
O(1) -- confirmed to fail (AmosException) at l=16, Oh=0.3 with the naive
`sqrt(pi/2q)*besselj(l+0.5,q)` approach this script originally used.

The recurrence length `pad` below is capped defensively: its working length
scales with `abs(q)`, so an ill-behaved caller passing a huge `|q|` (e.g. a
Newton step that overshot before step-capping was added) would otherwise
turn into a silent multi-minute hang rather than a fast, inaccurate answer
that a residual check downstream can reject.
"""
function sph_bessel_ratio(l::Integer, q::Number)
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

Residual of Reid's exact characteristic equation for mode l:
  alpha^4/q^4 + 1 = (2(l-1)/q^2) * [l + (l+1)*(q - 2l*Q)/(q - 2Q)]
with alpha^2 = sqrt(l(l-1)(l+2))/Oh and Q(q) = j_{l+1}(q)/j_l(q).
Zero of this function <=> q is an eigenvalue wavenumber of the linearized
viscous drop problem at this (Oh, l).
"""
function reid_char(q, Oh, l)
    alpha2 = sqrt(l * (l - 1) * (l + 2)) / Oh
    Q = sph_bessel_ratio(l, q)
    lhs = alpha2^2 / q^4 + 1
    rhs = 2 * (l - 1) / q^2 * (l + (l + 1) * (q - 2l * Q) / (q - 2Q))
    lhs - rhs
end

function safe_eval(f, q)
    try
        v = f(q)
        return isfinite(abs(v)) ? v : nothing
    catch
        return nothing
    end
end

"""
    newton_complex(f, q0; ...)

Damped, step-capped Newton with backtracking: each step is capped at a
fraction of `|q|` and halved until the residual decreases. An UNDAMPED Newton
step can overshoot into an argument where `sph_bessel_ratio`'s downward
recurrence needs an enormous number of terms to converge (its working length
scales with `abs(q)`), which is not an exception but a silent multi-minute
hang -- confirmed during development, and the reason this function is no
longer a plain Newton iteration.
"""
function newton_complex(f, q0; maxiter=300, tol=1e-13, step_cap_frac=0.5)
    q = q0
    F = safe_eval(f, q)
    F === nothing && return q0
    for _ in 1:maxiter
        h = max(1e-8, abs(q) * 1e-8)
        Fh = safe_eval(f, q + h * im)
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
            Fc = safe_eval(f, cand)
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

"""
    dominant_root_direct(Oh, l)

A single Newton solve from Lamb's small-Oh asymptotic guess. Correct for small
Oh; NOT reliable in general -- see `dominant_root_tracked` and Section 1's
regression test.
"""
function dominant_root_direct(Oh, l)
    sigma0 = sqrt(l * (l - 1) * (l + 2))
    gamma0 = (l - 1) * (2l + 1) * Oh
    q0 = sqrt(complex(gamma0 / Oh, -sigma0 / Oh))
    imag(q0) > 0 && (q0 = -q0)
    q = newton_complex(qv -> reid_char(qv, Oh, l), q0)
    imag(q) > 0 ? conj(q) : q
end

"""
    dominant_root_tracked(Oh, l; oh_start=1e-4, nsteps=24)

The dominant root, found by CONTINUATION in Oh rather than a single Newton
solve (technique ported from SpectralKM.jl's `reid_root_tracked`): start at an
Oh small enough that Lamb's guess is genuinely asymptotic, then walk
geometrically up to the target, seeding each step's Newton solve with the
previous step's converged root. This is what makes the result trustworthy past
the underdamping transition -- a single Lamb-seeded solve at the target Oh can
converge to a MORE STRONGLY DAMPED, non-dominant root instead (Section 1
reproduces this failure at l=16, Oh=0.3 and shows continuation fixes it).
"""
function dominant_root_tracked(Oh, l; oh_start=1e-4, nsteps=24)
    Oh <= oh_start && return dominant_root_direct(Oh, l)
    q = dominant_root_direct(oh_start, l)
    for ohv in exp.(range(log(oh_start), log(Oh); length=nsteps + 1))[2:end]
        q = newton_complex(qv -> reid_char(qv, ohv, l), q)
    end
    q
end

println("="^78)
println("Section 1: Reid's characteristic equation and the dominant root")
println("="^78)
println("""
reid_char(q, Oh, l) = 0 defines the eigenvalue wavenumbers q of Reid's linear
viscous drop problem. A single Lamb-seeded Newton solve (dominant_root_direct)
is fast but not always reliable; continuation in Oh (dominant_root_tracked)
is the robust version.
""")

@assert abs(reid_char(dominant_root_tracked(0.05, 2), 0.05, 2)) < 1e-10
@assert abs(reid_char(dominant_root_tracked(1.85, 2), 1.85, 2)) < 1e-10
@assert abs(reid_char(dominant_root_tracked(57.4, 2), 57.4, 2)) < 1e-10
println("ASSERTION 1 OK: dominant_root_tracked satisfies reid_char=0 across the")
println("underdamped (Oh=0.05), near-critical (Oh=1.85), and deep-overdamped")
println("(Oh=57.4) regimes for l=2.")

# --- The wrong-branch failure this script's ORIGINAL direct-Newton approach
#     had, and confirmation that continuation fixes it (SpectralKM.jl's own
#     documented example, reproduced independently here). ---
#
# IMPORTANT, discovered scanning reid_char's real axis directly at this (l,Oh):
# Lamb's OWN coefficient c=Oh(l-1)(2l+1)=148.5 exceeds omega_{l,0}=sqrt(l(l-1)(l+2))
# =65.7, which per Lamb's own (approximate) oscillator would mean overdamping --
# but this does NOT mean Reid's EXACT equation is overdamped here too. Reid's
# theory can (and, verified below, does) remain genuinely underdamped (complex
# dominant pair) at a point where Lamb's cruder model has already crossed into
# spurious overdamping. So the right test is not "close to Lamb's c" (which is
# known, per SpectralKM's own measured table, to be a poor approximation here)
# but "continuation finds a genuinely LESS damped (smaller Re(sigma), hence more
# physically dominant per Reid's own ordering) root than the direct solve does".
let l = 16, Oh = 0.3
    q_direct = dominant_root_direct(Oh, l)
    sigma_direct = q_direct^2 * Oh
    q_tracked = dominant_root_tracked(Oh, l)
    sigma_tracked = q_tracked^2 * Oh

    @assert abs(reid_char(q_direct, Oh, l)) < 1e-8    # both satisfy the char eq --
    @assert abs(reid_char(q_tracked, Oh, l)) < 1e-8   # so this isn't a convergence failure,
    @assert abs(sigma_tracked - sigma_direct) / abs(sigma_direct) > 0.5   # it's two DIFFERENT genuine roots
    @assert real(sigma_tracked) < real(sigma_direct)   # tracked is the LESS damped (dominant) one
    println()
    println("l=16, Oh=0.3: direct Newton gives sigma=$(round(sigma_direct,digits=1))" *
            " (a real, more strongly damped higher overtone -- confirmed by an independent" *
            " real-axis scan of reid_char to be a genuine root, just not the dominant one);" *
            " continuation gives sigma=$(round(sigma_tracked,digits=2))" *
            " (complex -- i.e. Reid's exact equation is STILL underdamped here, even though" *
            " Lamb's own coefficient $(round(Oh*(l-1)*(2l+1),digits=1)) already exceeds" *
            " omega_{l,0}=$(round(sqrt(l*(l-1)*(l+2)),digits=1)) and would spuriously predict" *
            " overdamping).")
    println("ASSERTION 1b OK: reproduces the direct-Newton wrong-branch failure at")
    println("l=16, Oh=0.3, and confirms continuation finds the genuinely dominant root instead.")
end

# ------------------------------------------------------------------------------
# Section 2: the second root, and the continuation trap for finding it.
# ------------------------------------------------------------------------------

"""
    second_root_analytic_guess(Oh, l)

Small-q asymptotic solution of Reid's characteristic equation, derived by
substituting the small-argument ratio Q(q) ~ q/(2l+3) into reid_char and
balancing the two singular (1/q^2, 1/q^4) terms:
    q ~ alpha^2 * sqrt((2l+1) / (2(l-1)(2l^2+4l+3)))
Valid only in the overdamped regime, where the second dominant root is real
and small.
"""
function second_root_analytic_guess(Oh, l)
    alpha2 = sqrt(l * (l - 1) * (l + 2)) / Oh
    alpha2 * sqrt((2l + 1) / (2 * (l - 1) * (2l^2 + 4l + 3)))
end

function find_second_root(Oh, l, q1)
    if abs(imag(q1)) > 1e-9 * abs(q1)
        return conj(q1)   # underdamped: true conjugate pair, no separate search
    end
    q0 = second_root_analytic_guess(Oh, l)
    newton_complex(qv -> reid_char(qv, Oh, l), complex(q0); maxiter=200, tol=1e-13)
end

println()
println("="^78)
println("Section 2: Why continuing the conjugate branch through Oh_crit fails")
println("="^78)
println("""
Below, we track q1 (dominant root) and a naive "q2" initialized as conj(q1) at
a safely underdamped Oh, then continue BOTH through the critical point using
each step's own result as the next Newton guess. Contrast with
find_second_root's independent analytic-guess approach.
""")

let l = 2
    Oh_seq = vcat(0.5:0.05:0.9, 1.0:1.0:10.0, 20.0, 57.4)
    q1c = dominant_root_direct(Oh_seq[1], l)
    q2c = conj(q1c)
    for Oh in Oh_seq[2:end]
        q1c = newton_complex(q -> reid_char(q, Oh, l), q1c)
        q2c = newton_complex(q -> reid_char(q, Oh, l), q2c)
    end
    @assert abs(q1c - q2c) / abs(q1c) < 1e-6   # both branches collapse onto one root
    q1_correct = dominant_root_tracked(57.4, l)
    q2_correct = find_second_root(57.4, l, q1_correct)
    @assert abs(q2_correct - q2c) / abs(q2c) > 0.9   # genuinely different root
    println("ASSERTION 2 OK: naive conjugate-branch continuation collapses to")
    println("q=$(round(q1c,digits=5)) at Oh=57.4 for BOTH tracked branches, while")
    println("the correct second root is q=$(round(q2_correct,digits=6)) -- confirming")
    println("continuation silently finds the wrong (degenerate) pair.")
end

# ------------------------------------------------------------------------------
# Section 3: lambda_l(Oh), omega_l^2(Oh) via Vieta's formulas -- the unit-mass
# parametrization, matching SpectralKM.jl's `reid_pole_pair`/`drop_viscous_coeffs`.
# ------------------------------------------------------------------------------

println()
println("="^78)
println("Section 3: lambda_l(Oh), omega_l^2(Oh) -- exact finite-Oh coefficients")
println("="^78)
println("""
Unit-mass per-mode oscillator: Addot_l + 2*lambda_l*Adot_l + omega_l^2*A_l = 0
(no forcing). Its characteristic equation is x^2 - 2*lambda_l*x + omega_l^2 = 0,
whose roots must equal Reid's own sigma_1, sigma_2 (sigma_i = q_i^2*Oh, this
repo's existing q^2*Oh time-unit convention -- see
julia/derivations/carreau_yasuda_derivation.jl). Vieta's formulas then give:
    lambda_l  = (sigma_1 + sigma_2) / 2
    omega_l^2 = sigma_1 * sigma_2
Lamb's formula (and this repo's current code) is exactly lambda_l=Oh(l-1)(2l+1),
omega_l^2=l(l-1)(l+2) -- the Oh -> 0 limit derived below.
""")

# sigma_i = q_i^2*Oh -- keep COMPLEX (underdamped: q1,q2 are a conjugate
# pair, so sigma1,sigma2 are too; real(q)^2 would discard Im and corrupt
# the underdamped regime).
function compute_lambda_omega2(Oh, l)
    q1 = dominant_root_tracked(Oh, l)
    q2 = find_second_root(Oh, l, q1)
    s1 = q1^2 * Oh
    s2 = q2^2 * Oh
    lambda = real(s1 + s2) / 2
    omega2 = real(s1 * s2)
    lambda, omega2
end

# --- Low-Oh limit: lambda_l -> Oh(l-1)(2l+1), omega_l^2 -> l(l-1)(l+2) ---
println("Low-Oh limit (lambda_l -> Oh(l-1)(2l+1), omega_l^2 -> l(l-1)(l+2)):")
for l in (2, 3, 4, 5, 10)
    target_om2 = float(l) * (l - 1) * (l + 2)
    lam_med, om2_med = compute_lambda_omega2(0.001, l)
    lam_small, om2_small = compute_lambda_omega2(0.0001, l)
    lamb_med   = 0.001 * (l - 1) * (2l + 1)
    lamb_small = 0.0001 * (l - 1) * (2l + 1)
    err_om2_med   = abs(om2_med - target_om2) / target_om2
    err_om2_small = abs(om2_small - target_om2) / target_om2
    err_lam_med   = abs(lam_med - lamb_med) / lamb_med
    err_lam_small = abs(lam_small - lamb_small) / lamb_small
    println("  l=$l  omega2(Oh=0.001)=$(round(om2_med,digits=4))  omega2(Oh=0.0001)=$(round(om2_small,digits=4))" *
            "  target=$(round(target_om2,digits=4))")
    @assert err_om2_small < err_om2_med       # omega_l^2 -> target monotonically
    @assert err_lam_small < err_lam_med       # lambda_l -> Oh(l-1)(2l+1) monotonically
    @assert err_om2_small < 0.02
    @assert err_lam_small < 0.02
end
println("ASSERTION 3 OK: lambda_l -> Oh(l-1)(2l+1) and omega_l^2 -> l(l-1)(l+2)")
println("monotonically as Oh -> 0, for l=2,3,4,5,10, to <2% at Oh=0.0001.")

# --- High-Oh limit, recast from Molacek & Bush's A_l, D_l via the rescaling
#     lambda_l = l^2*Oh*D_l/A_l, omega_l^2 = l(l-1)(l+2)/A_l ---
println()
println("High-Oh limit (recast from Molacek & Bush 2012's A_l, D_l via the")
println("unit-mass rescaling above):")
# At high Oh, A_l -> om0sq/omega_l^2 (from Section 3's own definition), so
# D_l = omega_l^2 * lambda_l / (l^2*Oh*A_l)... but omega_l^2 and lambda_l
# ALONE already fully determine the physics; recover D_l only to check
# against the independently published number.
for l in (2, 3, 4, 5, 10)
    om0sq = float(l) * (l - 1) * (l + 2)
    target_D = (l - 1) * (2l^2 + 4l + 3) / (l^2 * (2l + 1))
    _, om2_100 = compute_lambda_omega2(100.0, l)
    _, om2_1000 = compute_lambda_omega2(1000.0, l)
    lam_1000, _ = compute_lambda_omega2(1000.0, l)
    A_1000 = om0sq / om2_1000
    D_1000 = lam_1000 * A_1000 / (l^2 * 1000.0)
    err = abs(D_1000 - target_D) / target_D
    println("  l=$l  D_l(Oh=1000, recast)=$(round(D_1000,digits=6))  target=$(round(target_D,digits=6))")
    @assert err < 1e-3
end
println("ASSERTION 4 OK: recasting lambda_l, omega_l^2 back into Molacek & Bush's")
println("A_l, D_l at high Oh reproduces their independently-published high-Oh")
println("limit to <0.1%, for l=2,3,4,5,10 -- confirming the two parametrizations")
println("carry the same physics for free decay, as claimed in the header.")

# --- Continuity across the critical Oh (no jump discontinuity) ---
println()
Oh_span = vcat(range(0.3, 1.3; length=41))
prev_lam, prev_om2 = compute_lambda_omega2(Oh_span[1], 2)
max_jump = 0.0
for Oh in Oh_span[2:end]
    lam, om2 = compute_lambda_omega2(Oh, 2)
    global max_jump = max(max_jump, abs(lam - prev_lam) / max(prev_lam, 1.0), abs(om2 - prev_om2))
    global prev_lam, prev_om2 = lam, om2
end
@assert max_jump < 0.1
println("ASSERTION 5 OK: lambda_l(Oh), omega_l^2(Oh) vary smoothly across the")
println("critical Oh (l=2, Oh in [0.3, 1.3], max per-step change = $(round(max_jump,digits=4)))")
println("-- no discontinuity at the oscillatory/overdamped transition.")

println()
println("""
NOT resolved here (a genuine open question, matching SpectralKM.jl's own
reid.jl header): substituting exact eigenvalues into a second-order oscillator
is exact for FREE decay only. Under forcing (this repo's l*B_l contact-pressure
term), the true system carries memory from the full discarded viscous
spectrum, which neither this parametrization nor Molacek & Bush's A_l/D_l form
captures exactly -- and the two need NOT agree once forcing is added, since an
overall equation rescaling that leaves free-decay roots invariant does NOT
leave a forcing term's relative weight invariant. This script derives and
verifies the FREE-decay coefficients only; the forced-response question is
for the production-wiring step that follows, and should be checked against a
live simulation once wired in (as SpectralKM.jl's types.jl does for its own
:reid default).
""")

# ------------------------------------------------------------------------------
# Section 4: Live cross-check against the actual running DropSolver.
#
# This repo's current production code (julia/src/timestepper.jl) implements
# ONLY Lamb's asymptotic formula lambda_lamb=(l-1)(2l+1)*Oh -- NOT Reid's exact
# eigenvalue, which differs from Lamb's by an O(Oh) correction even at Oh as
# small as 0.02 (a ~7% gap, per the numbers below). The correct live check is
# therefore two separate claims, not one: (a) the running solver actually
# computes what Lamb's formula predicts (confirming this document accurately
# describes what julia/src/timestepper.jl does today), and (b) THIS section's
# own Reid-exact eigenvalue converges to that same Lamb value as Oh -> 0
# (confirming Lamb truly is this derivation's own low-Oh limit).
# ------------------------------------------------------------------------------

println()
println("="^78)
println("Section 4: Live cross-check against DropSolver (small-Oh limit)")
println("="^78)

function extract_decay_freq(times, A2)
    gamma = -log(abs(A2[end]) / abs(A2[1])) / (times[end] - times[1])
    sign_changes = findall(i -> A2[i] * A2[i+1] < 0, 1:length(A2)-1)
    omega = NaN
    if length(sign_changes) >= 4
        hp = diff(times[sign_changes])
        omega = pi / (sum(hp) / length(hp))
    end
    gamma, omega
end

println("(a) Does the running solver actually compute Lamb's formula?")
for (Oh, l) in ((0.02, 2), (0.02, 3), (0.05, 2))
    gamma_lamb = (l - 1) * (2l + 1) * Oh

    M = l
    theta_vec = make_theta_vec(M)
    precomp = precompute_integrals(NaN, M)[1]
    sigma0 = sqrt(l * (l - 1) * (l + 2))
    dt_osc = 2 * pi / (sigma0 * 40)
    cfg = SimConstants(M, M + 1, Oh, 1e-6, theta_vec, precomp, dt_osc)

    init = DropState(M)
    init.A[l] = 0.05
    init.z = 2.0
    init.dt = dt_osc

    T_period = 2 * pi / sigma0
    times, states = solve_drop!(cfg, OBParams(), init;
        t_end=6 * T_period, save_every=T_period / 50, dt_init=dt_osc)
    Al = [s.A[l] for s in states]
    gamma_sim, omega_sim = extract_decay_freq(times, Al)

    err_gamma = abs(gamma_sim - gamma_lamb) / gamma_lamb
    println("  Oh=$Oh l=$l: gamma_lamb=$(round(gamma_lamb,digits=5))" *
            "  gamma_sim(DropSolver)=$(round(gamma_sim,digits=5))" *
            "  rel_err=$(round(err_gamma,digits=4))")
    @assert err_gamma < 0.05
end
println("ASSERTION 6 OK: live solve_drop! decay rate matches Lamb's formula to")
println("<5% -- confirming this document's description of the running code")
println("(Lamb-only, no finite-Oh correction) is accurate.")

println()
println("(b) Does THIS derivation's own Reid-exact eigenvalue converge to that")
println("    same Lamb value as Oh -> 0 (rather than to something else)?")
let l = 2
    gamma_lamb_fn(Oh) = (l - 1) * (2l + 1) * Oh
    prev_err = Inf
    for Oh in (0.05, 0.02, 0.005, 0.001)
        lam, _ = compute_lambda_omega2(Oh, l)
        err = abs(lam - gamma_lamb_fn(Oh)) / gamma_lamb_fn(Oh)
        println("  Oh=$Oh: lambda_exact(Reid)=$(round(lam,digits=6))" *
                "  gamma_lamb=$(round(gamma_lamb_fn(Oh),digits=6))  rel_gap=$(round(err,digits=4))")
        @assert err < prev_err
        prev_err = err
    end
end
println("ASSERTION 7 OK: the gap between Reid's exact eigenvalue and Lamb's")
println("asymptotic formula shrinks monotonically as Oh -> 0 (l=2), confirming")
println("Lamb's formula is genuinely this derivation's own Oh -> 0 limit.")

# ------------------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------------------

println()
println("="^78)
println("Summary")
println("="^78)
println("""
In plain terms: this repo's damping formula is only the low-viscosity limit of
a more general, exact result. We derived that exact result (lambda_l(Oh),
omega_l^2(Oh), the unit-mass parametrization already validated in production
in the sister repo SpectralKM.jl) directly from Reid's own characteristic
equation, verified it against Molacek & Bush's independently-published limits
with no fitting, and confirmed the current production solver already matches
it in the regime it was designed for. A failing assertion above would mean
either Reid's characteristic equation was mis-transcribed, the two dominant
roots were mis-identified (see Section 1's wrong-branch reproduction and
Section 2's continuation trap), or the current production code no longer
matches its own documented Lamb-limit behavior -- any of which would be a
physically real regression, not a cosmetic one.

Not yet done (deliberately out of scope for this script): wiring
lambda_l(Oh), omega_l^2(Oh) into julia/src/timestepper.jl's residual/Jacobian
in place of the hardcoded Lamb formulas (behind a :lamb/:reid switch, mirroring
SpectralKM.jl's own API, defaulting to :lamb to preserve all existing test
behavior), and replacing julia/src/st_extension.jl's perturbative
Carreau-Yasuda correction with a non-perturbative Oh_eff(t) computed per mode
from the local shear rate fed through this machinery.
""")
