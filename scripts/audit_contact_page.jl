# Measures every quantitative claim made by `docs/src/contact.md`.
#
# The page asserts things a reader cannot check by reading it: how large the
# complementarity problem actually is, how closely the two closures agree, and
# whether the contact region the solver returns is a single patch containing the
# pole. Each is measured here so the published number has a source.
#
#   1. LCP size against M+1, which is NOT the same number: the complementarity
#      problem is posed on the lower hemisphere, and at odd M one node lands on
#      the equator with an identically zero clearance row.
#   2. Closure agreement over the grid the page quotes.
#   3. Contact-set structure: is the pole ever released while a ring still
#      presses, and is the set contiguous.

using Printf, LinearAlgebra
using DropSolver

flushln(a...) = (println(a...); flush(stdout))
const OUT = joinpath(@__DIR__, "..", "outputs", "csv")

# ---------------------------------------------------------------------------
# 1. How big is the complementarity problem, really?
# ---------------------------------------------------------------------------
flushln("== LCP size vs collocation count ==")
flushln(@sprintf("%5s %8s %8s %10s", "M", "nodes", "lcp", "equator?"))
for M in (20, 21, 30, 31, 45, 60, 61, 90, 91)
    q  = ImpactParams(We = 0.5, Bo = 0.019, Oh = 0.0373, M = M, K = 3)
    F0 = DropSolver.assemble_newtonian(DropSolver.basis(q), q.Oh)
    Vf = lu(DropSolver.legendre_vandermonde(q))
    s0 = DropSolver.initial_state(q)
    _, _, idx, _ = DropSolver.contact_lcp(q, s0, s0, q.dt0; F0 = F0, Vfac = Vf)
    nn  = length(q.nodes)
    eq  = any(i -> abs(cos(q.nodes[i])) < 1e-8, 1:nn)
    flushln(@sprintf("%5d %8d %8d %10s", M, nn, length(idx), eq ? "yes" : "no"))
end

# ---------------------------------------------------------------------------
# 3. Does the contact set contain the pole, and is it one patch?
#
# `r.pc` is the Legendre pressure field per frame, so the nodal pressures are
# `V * pc` and a node is in contact exactly where that is positive. This is the
# LCP unknown itself, recovered without instrumenting the solver.
# ---------------------------------------------------------------------------
function set_structure(We, Oh, M)
    ## simulate_lcp rather than run_impact: only the raw solver return carries the
    ## per-frame pressure field, which is what the contact set has to be read from.
    q = ImpactParams(We = We, Bo = 0.019, Oh = Oh, M = M, K = 3, t_max = 25.0)
    r = DropSolver.simulate_lcp(q)
    V = DropSolver.legendre_vandermonde(q)
    lower = [i for i in 1:length(q.nodes) if cos(q.nodes[i]) < -1e-8]
    frames = pole_free = noncontig = 0
    for pc in r.pc
        pn  = (V * pc)[lower]                     # nodal pressures, pole first
        act = findall(>(1e-12), pn)
        isempty(act) && continue
        frames += 1
        1 in act || (pole_free += 1)              # a ring pressing, pole released
        (act == first(act):last(act)) || (noncontig += 1)
    end
    (frames, pole_free, noncontig, r.cor, r.noncontiguous_steps)
end

flushln("\n== contact-set structure ==")
flushln(@sprintf("%5s %5s %7s | %8s %10s %12s", "M", "We", "Oh", "frames", "pole free", "non-contig"))
for M in (45, 90), (We, Oh) in ((0.5, 0.0373), (2.0, 0.0373))
    f, pf, nc, cor, ncs = set_structure(We, Oh, M)
    flushln(@sprintf("%5d %5.1f %7.4f | %8d %10d %12d   (cor %.4f, solver says %d)",
                     M, We, Oh, f, pf, nc, cor, ncs))
end

# ---------------------------------------------------------------------------
# 2. Do the two closures agree, over the grid the page quotes?
# ---------------------------------------------------------------------------
flushln("\n== closure agreement, M = 45 ==")
ohs = exp10.(range(log10(0.02), log10(0.7); length = 5))
wes = range(0.1, 3.0; length = 5)
dcor = Float64[]; dtc = Float64[]
open(joinpath(OUT, "audit_contact_page_closures.csv"), "w") do io
    println(io, "We,Oh,cor_lcp,cor_as,dcor,tc_lcp,tc_as,wall_lcp,wall_as")
    for Oh in ohs, We in wes
        kw = (We = We, Bo = 0.019, Oh = Oh, M = 45, K = 3, t_max = 25.0)
        a = run_impact(Backend(contact = :lcp); kw...)
        b = run_impact(Backend(contact = :active_set); kw...)
        if !(a.ok && b.ok)
            flushln(@sprintf("We %4.2f Oh %6.4f | SKIPPED (lcp ok=%s, as ok=%s)", We, Oh, a.ok, b.ok))
            continue
        end
        d = abs(a.cor - b.cor); push!(dcor, d); push!(dtc, abs(a.tc - b.tc))
        println(io, join((We, Oh, a.cor, b.cor, d, a.tc, b.tc, a.wall, b.wall), ","))
        flush(io)
        flushln(@sprintf("We %4.2f Oh %6.4f | dcor %.3e  dtc %.3e  (%.0fs / %.0fs)",
                         We, Oh, d, abs(a.tc - b.tc), a.wall, b.wall))
    end
end
## The same comparison on the 3000 ppm shear-thinning fluid, whose viscosity is
## rebuilt inside the step and so exercises a path the Newtonian grid never touches.
flushln("\n== closure agreement, 3000 ppm shear thinning, M = 45 ==")
let R_DROP = 0.0003, SIGMA = 0.0728, RHO = 1000.0, BO = 0.012
    T_CAP = sqrt(RHO * R_DROP^3 / SIGMA)
    eta_inf, eta_0, k, n = 0.0037320997942061666, 8.433817577956766,
                           18.48081673111359, 0.7430524574330837
    Oh_0 = eta_0 / sqrt(RHO * SIGMA * R_DROP)
    eta_fn = gd -> carreau(gd; lambda_c = k / T_CAP, a = n, n = 1 - n,
                           eta_inf_ratio = eta_inf / eta_0)
    st = Float64[]
    for We in (0.2, 0.5, 1.0, 2.0)
        kw = (We = We, Bo = BO, Oh = Oh_0, M = 45, K = 3, eta = eta_fn, t_max = 25.0)
        a = run_impact(Backend(contact = :lcp); kw...)
        b = run_impact(Backend(contact = :active_set); kw...)
        if !(a.ok && b.ok)
            flushln(@sprintf("We %4.2f | SKIPPED (lcp ok=%s, as ok=%s)", We, a.ok, b.ok)); continue
        end
        push!(st, abs(a.cor - b.cor))
        flushln(@sprintf("We %4.2f | dcor %.3e  dtc %.3e  (%.0fs / %.0fs)",
                         We, abs(a.cor - b.cor), abs(a.tc - b.tc), a.wall, b.wall))
    end
    flushln(@sprintf("shear-thinning worst |dcor| = %.3e over %d cases",
                     isempty(st) ? NaN : maximum(st), length(st)))
end

flushln(@sprintf("\ncompared %d of 25 impacts; worst |dcor| = %.3e, worst |dtc| = %.3e",
                 length(dcor), isempty(dcor) ? NaN : maximum(dcor),
                 isempty(dtc) ? NaN : maximum(dtc)))
flushln("done")
