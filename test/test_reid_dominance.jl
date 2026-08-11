# Is the root `dominant_root` tracks actually the dominant one?
#
# Every other Reid test asserts on the characteristic-equation residual, and
# that check is structurally incapable of answering this: reid_char has many
# roots, and ALL of them drive the residual to ~1e-14. A continuation that
# slipped onto a different branch would look perfect. These tests use an
# independent criterion -- find the roots by multi-start Newton from a net over
# the complex q-plane, then compare damping rates.
#
# `dominant_root` was previously only exercised to l = 10 while production runs
# reach M = 90.

"""
Distinct roots of `reid_char(., Oh, l)` found by multi-start Newton, keeping
only those with `Re(sigma) > 0` (`sigma = q^2*Oh`, decaying) and, if
`complex_only`, only those with a genuine oscillation.
"""
function _rd_roots(Oh::Float64, l::Int; complex_only::Bool = true,
                   nrad::Int = 18, nang::Int = 28)
    scale = sqrt(sqrt(l * (l - 1) * (l + 2)) / Oh)
    roots = ComplexF64[]
    for rr in exp.(range(log(0.05 * max(scale, 1.0)), log(30 * max(scale, 1.0)); length = nrad)),
        th in range(-pi, pi; length = nang)
        q = try
            DropSolver._newton_complex(qv -> DropSolver.reid_char(qv, Oh, l),
                                       rr * cis(th); maxiter = 300, tol = 1e-14)
        catch
            continue
        end
        (isfinite(abs(q)) && abs(q) > 1e-8) || continue
        abs(DropSolver.reid_char(q, Oh, l)) < 1e-8 || continue
        q = imag(q) > 0 ? conj(q) : q
        s = q^2 * Oh
        real(s) > 1e-9 || continue
        complex_only && !(abs(imag(s)) > 1e-6 * abs(s)) && continue
        any(abs(q - r) < 1e-6 * max(1.0, abs(q)) for r in roots) || push!(roots, q)
    end
    roots
end

"""
Largest ratio between the tracked damping at adjacent points of a fine
log-spaced `Oh` sweep. The tracked branch is a smooth function of `Oh`, so on a
grid this fine the ratio stays near 1; a continuation that slipped onto a
different root shows up as a step change.
"""
function _rd_worst_jump(l::Int, max_step_ratio::Float64;
                        Oh_lo = 1e-3, Oh_hi = 1e4, n = 140)
    prev, worst = NaN, 0.0
    for Oh in exp.(range(log(Oh_lo), log(Oh_hi); length = n))
        s = real(DropSolver.dominant_root(Oh, l; max_step_ratio)^2 * Oh) / Oh
        if isfinite(prev) && prev > 0 && s > 0
            worst = max(worst, max(s / prev, prev / s))
        end
        prev = s
    end
    worst
end

@testset "Reid root dominance" begin

    @testset "Tracked branch is continuous in Oh" begin
        # The dominance test below is confined to the underdamped regime, so it
        # says nothing about Oh past critical -- which is where the branch jump
        # this package documents was actually observed (l = 10, lambda
        # 61483 -> 273551 across Oh 958 -> 1014). This one sweeps the whole
        # range and needs no notion of which root is dominant: the tracked
        # branch is a smooth function of Oh, so a jump is a jump.
        #
        # Correct tracking holds the worst adjacent ratio to <= 1.21 across
        # l = 2..120 over six decades of Oh; breaking the step ratio to 1e6 (a
        # single Lamb-seeded Newton solve at the target) gives 8x to 62x. The
        # threshold sits in that gap.
        for l in (2, 10, 30, 90, 120)
            @test _rd_worst_jump(l, 1.15) < 2.0
        end
    end

    @testset "Tracked root is the least-damped complex root" begin
        # In the underdamped regime the pair is {q1, conj(q1)}, so "q1 is
        # dominant" means exactly: no complex root decays more slowly. Spans
        # l = 2..120, past the M = 90 production runs use.
        #
        # Which (l, Oh) are underdamped is decided from the independently-found
        # roots, NOT from the tracked root: keying the filter on the code under
        # test lets a mutation silently skip cases instead of failing them.
        for l in (2, 10, 30, 60, 90, 120), Oh in (0.01, 0.05, 0.2)
            rs = _rd_roots(Oh, l)
            isempty(rs) && continue                          # genuinely overdamped
            s1 = DropSolver.dominant_root(Oh, l)^2 * Oh
            @test abs(imag(s1)) > 1e-6 * abs(s1)             # so must the tracked root be
            @test real(s1) <= minimum(real(q^2 * Oh) for q in rs) * (1 + 1e-8)
        end
    end

    @testset "Overdamped pair is {continuation branch, creep mode}" begin
        # Past the critical Oh the pair is deliberately NOT the two
        # slowest-decaying roots: reid_char then has a dense set of overdamped
        # roots, and Reid's reduction keeps the analytic continuation of the
        # Lamb branch together with the slow creep mode. Pinned here because
        # the two-slowest reading is the natural wrong guess, and a probe built
        # on it reports the textbook l=2 case as broken.
        #
        # The creep root has a closed form: substituting the small-argument
        # ratio Q(q) ~ q/(2l+3) into reid_char and balancing its two singular
        # terms gives sigma_creep = l(l+2)(2l+1) / (2*Oh*(2l^2+4l+3)).
        # That form is a LARGE-Oh asymptotic, so it is held to a tight
        # tolerance only well past the critical Oh, and the approach itself is
        # asserted separately: the relative error must shrink as Oh grows.
        # (At Oh = 2, just past critical for l = 2, it is still 3.9% off.)
        creep_err(l, Oh) = begin
            q1 = DropSolver.dominant_root(Oh, l)
            @test abs(imag(q1^2 * Oh)) < 1e-6 * abs(q1^2 * Oh)   # is overdamped
            s2 = real(DropSolver.second_root(Oh, l, q1)^2 * Oh)
            @test s2 < real(q1^2 * Oh)          # genuinely the slow one of the pair
            want = l * (l + 2) * (2l + 1) / (2 * Oh * (2l^2 + 4l + 3))
            abs(s2 - want) / want
        end
        for l in (2, 5, 10)
            errs = [creep_err(l, Oh) for Oh in (2.0, 10.0, 57.0)]
            @test errs[3] < 1e-3                       # converged far past critical
            @test errs[1] > errs[2] > errs[3]          # and approaching monotonically
        end
    end

    @testset "Residual alone cannot detect a wrong branch" begin
        # The reason these tests exist. A non-dominant root satisfies the
        # characteristic equation just as well as the dominant one, so any
        # residual threshold accepts both.
        l, Oh = 10, 10.0
        rs = _rd_roots(Oh, l; complex_only = false)
        @test length(rs) >= 3                       # several genuine roots
        @test all(abs(DropSolver.reid_char(q, Oh, l)) < 1e-8 for q in rs)
        damping = sort!([real(q^2 * Oh) for q in rs])
        @test damping[end] > 10 * damping[1]        # wildly different physics
    end
end
