# How well is the CONTACT PATCH resolved, as distinct from how converged the KPIs are.
#
# Restitution converging in M says the integrated bounce is right. It says nothing
# about whether the contact region is represented during the impact, and those are
# different questions. Contact is imposed at M+1 collocation nodes clustered at the
# pole, so the patch is only ever a whole number of nodes wide: at M = 45 a patch
# covering a tenth of the lower hemisphere is a couple of nodes across, and it has to
# grow through that handful in the few steps a high-Weber impact allows.
#
# Three numbers per run:
#   maxcp     how many nodes are in contact at the widest, so how many distinct
#             contact radii the discretisation can even express
#   jumps>1   accepted steps where the patch changed by more than one node, which is
#             the patch moving faster than the grid can follow
#   theta_c   the half-angle of the widest patch, so the physical size being resolved

using Printf
using DropSolver

med(v) = (w = sort(v); n = length(w); n == 0 ? NaN :
          iseven(n) ? (w[n÷2] + w[n÷2+1])/2 : w[(n+1)÷2])
flushln(a...) = (println(a...); flush(stdout))

const OUT = joinpath(@__DIR__, "..", "outputs", "csv")

open(joinpath(OUT, "audit_contact_resolution.csv"), "w") do io
    println(io, "backend,M,K,We,Oh,ok,maxcp,nodes,frac_nodes,theta_c_deg,jumps_gt1,steps,wall_s,cor")
    flushln(@sprintf("%5s %5s %6s | %6s %7s %9s %8s %7s %8s",
                     "M", "We", "Oh", "maxcp", "%nodes", "theta_c", "jumps>1", "wall_s", "cor"))
    for (We, Oh) in ((0.5, 0.0373), (2.0, 0.0373), (3.0, 0.0200)), M in (30, 45, 60, 90)
        b = Backend(contact = :lcp)
        r = run_impact(b; We = We, Bo = 0.019, Oh = Oh, M = M, K = 3, t_max = 25.0)
        if isempty(r.cp)
            flushln(@sprintf("%5d %5.1f %6.3f | run produced no frames", M, We, Oh)); continue
        end
        cp = r.cp; nodes = M + 1
        mx = maximum(cp)
        d  = abs.(diff(cp)); nz = filter(>(0), d)
        ## the collocation nodes are theta = pi and the zeros of P_M, descending from
        ## the pole, so the mx-th node is the edge of the widest patch
        p  = ImpactParams(We = We, Bo = 0.019, Oh = Oh, M = M, K = 3)
        thc = mx >= 1 && mx <= length(p.nodes) ? rad2deg(pi - p.nodes[mx]) : NaN
        println(io, join((label(b), M, 3, We, Oh, r.ok, mx, nodes,
                          round(mx/nodes, digits=4), round(thc, digits=2),
                          count(>(1), d), length(cp), round(r.wall, digits=2), r.cor), ","))
        flushln(@sprintf("%5d %5.1f %6.3f | %6d %6.1f%% %8.2f° %8d %7.1f %8.4f",
                         M, We, Oh, mx, 100*mx/nodes, thc, count(>(1), d), r.wall, r.cor))
        flush(io)
    end
end
flushln("done")
