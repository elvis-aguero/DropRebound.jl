# The closure comparison of `audit_contact_page.jl`, repeated on the 3000 ppm
# shear-thinning fluid at the PRODUCTION truncation M = 90.
#
# M = 45 is not a fair test for this fluid: its viscosity is rebuilt inside every
# Picard sweep, and at M = 45 two of four Weber numbers exhaust the step-halving
# budget before releasing. The published claim is about the solver as it is run,
# so it is measured at the truncation it is run with.

using Printf
using DropSolver

flushln(a...) = (println(a...); flush(stdout))

const R_DROP = 0.0003
const SIGMA  = 0.0728
const RHO    = 1000.0
const BO     = 0.012
const T_CAP  = sqrt(RHO * R_DROP^3 / SIGMA)

eta_inf, eta_0, k, n = 0.0037320997942061666, 8.433817577956766,
                       18.48081673111359, 0.7430524574330837
Oh_0   = eta_0 / sqrt(RHO * SIGMA * R_DROP)
eta_fn = gd -> carreau(gd; lambda_c = k / T_CAP, a = n, n = 1 - n,
                       eta_inf_ratio = eta_inf / eta_0)

flushln(@sprintf("3000 ppm, Oh_0 = %.4g, M = 90", Oh_0))
dcor = Float64[]; dtc = Float64[]
for We in (0.2, 0.5, 1.0, 2.0)
    kw = (We = We, Bo = BO, Oh = Oh_0, M = 90, K = 3, eta = eta_fn, t_max = 25.0)
    a = run_impact(Backend(contact = :lcp); kw...)
    b = run_impact(Backend(contact = :active_set); kw...)
    if !(a.ok && b.ok)
        flushln(@sprintf("We %4.2f | SKIPPED (lcp ok=%s, as ok=%s)", We, a.ok, b.ok)); continue
    end
    push!(dcor, abs(a.cor - b.cor)); push!(dtc, abs(a.tc - b.tc))
    flushln(@sprintf("We %4.2f | dcor %.3e  dtc %.3e  (%.0fs / %.0fs)",
                     We, abs(a.cor - b.cor), abs(a.tc - b.tc), a.wall, b.wall))
end
flushln(@sprintf("compared %d of 4; worst |dcor| = %.3e, worst |dtc| = %.3e",
                 length(dcor), isempty(dcor) ? NaN : maximum(dcor),
                 isempty(dtc) ? NaN : maximum(dtc)))
