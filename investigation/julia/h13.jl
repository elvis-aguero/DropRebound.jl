# Does the tangency residual have a clean minimum in cp?  The ancestor's rule assumes it
# does -- it only ever looks two nodes either side. Compare K=1 (works) with K=2 (fails)
# at the same M, from an identical early-contact state.
using DropSolver, Printf
const D = DropSolver

function residual_profile(K; M = 90, nwarm = 6)
    p  = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.303767, M = M, K = K, t_max = 25.0)
    F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
    function et_of(prev, curr, dt, cp)
        (0 <= cp <= nn) || return (Inf, curr)
        st, nx = D.try_step(p, prev, curr, dt, cp; F0 = F0)
        st === :diverge && return (Inf, curr)
        g = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
        all(i -> g[i] >= 0, (cp+1):nn) || return (Inf, nx)
        cp == 0 && return (0.0, nx)
        cp >= nn && return (Inf, nx)
        (abs(g[cp+1] - g[cp]), nx)
    end
    ## warm up with the plain ancestor rule so both K land in a comparable state
    prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
    for _ in 1:nwarm
        e2,n2 = et_of(prev,curr,dt,curr.cp-1); e3,n3 = et_of(prev,curr,dt,curr.cp)
        e4,n4 = et_of(prev,curr,dt,curr.cp+1)
        if e3 > e4 || e3 > e2
            e5,_ = et_of(prev,curr,dt,curr.cp+2)
            if e4 <= e2 && e4 < e5
                prev, curr = curr, n4
            else
                dt /= 2
            end
        elseif isfinite(e3)
            prev, curr = curr, n3
        else
            dt /= 2
        end
    end
    ets = [et_of(prev, curr, dt, c)[1] for c in 0:10]
    (p, curr, dt, ets)
end

for K in (1, 2)
    p, curr, dt, ets = residual_profile(K)
    @printf("\nK=%d  (state: t=%.5f z=%.6f cp=%d, dt=%.3e)\n", K, curr.t, curr.z, curr.cp, dt)
    @printf("  cp:  ")
    for c in 0:10; @printf("%-11d", c); end
    @printf("\n  et:  ")
    for e in ets; @printf("%-11.4g", e); end
    fin = [(c, ets[c+1]) for c in 0:10 if isfinite(ets[c+1]) && c > 0]
    if !isempty(fin)
        mn = argmin([f[2] for f in fin])
        @printf("\n  finite minimum at cp=%d (et=%.4g);  monotone increasing after it? %s\n",
                fin[mn][1], fin[mn][2],
                issorted([f[2] for f in fin[mn:end]]) ? "yes" : "NO -- rule's assumption fails")
    else
        println("\n  no finite residual anywhere")
    end
end
