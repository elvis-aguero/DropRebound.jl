# H15: cp=0 carries a SENTINEL, not a residual. Exclude it from every tangency
# comparison; it is enterable only when cp=1 is infeasible (genuine release).
using DropSolver, Printf
const D = DropSolver

function march(p; excl0::Bool, nst = 60000)
    F0 = p.eta_const ? assemble_newtonian(D.basis(p), p.Oh) : nothing
    nn = length(p.nodes)
    function et_of(prev, curr, dt, cp)
        (0 <= cp <= nn) || return (Inf, curr, false)
        st, nx = D.try_step(p, prev, curr, dt, cp; F0=F0)
        st === :diverge && return (Inf, curr, false)
        g = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
        feas = all(i -> g[i] >= 0, (cp+1):nn)
        feas || return (Inf, nx, false)
        cp == 0 && return (excl0 ? Inf : 0.0, nx, true)   # sentinel, not a number
        cp >= nn && return (Inf, nx, true)
        (abs(g[cp+1] - g[cp]), nx, true)
    end
    prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
    ts=[0.0]; zs=[curr.z]; vs=[curr.v]; cps=[0]; as=[copy(curr.a)]; rej=0
    for _ in 1:nst
        curr.t >= p.t_max && break
        e2,n2,f2 = et_of(prev,curr,dt,curr.cp-1)
        e3,n3,f3 = et_of(prev,curr,dt,curr.cp)
        e4,n4,_  = et_of(prev,curr,dt,curr.cp+1)
        nxt = nothing; ncp = curr.cp
        if e3 > e4 || e3 > e2
            if e4 <= e2
                e5,_,_ = et_of(prev,curr,dt,curr.cp+2)
                e4 < e5 && (nxt = n4; ncp = curr.cp+1)
            else
                e1,_,_ = et_of(prev,curr,dt,curr.cp-2)
                e2 < e1 && (nxt = n2; ncp = curr.cp-1)
            end
        elseif isfinite(e3)
            nxt = n3
        end
        ## RELEASE: with cp=0 excluded from comparison, it is reachable only when the
        ## smallest contact is itself inadmissible -- which is what release means.
        if nxt === nothing && excl0 && curr.cp == 1 && !f3
            e0,n0,f0 = et_of(prev,curr,dt,0)
            f0 && (nxt = n0; ncp = 0)
        end
        if nxt === nothing
            rej += 1; dt /= 2; dt < p.dt_min && break; continue
        end
        prev, curr = curr, nxt
        push!(ts,curr.t); push!(zs,curr.z); push!(vs,curr.v); push!(cps,ncp)
        push!(as,copy(curr.a)); dt = p.dt0
        ncp == 0 && curr.v > 0 && curr.z > 1.0 && any(>(0), cps) && break
    end
    amp = maximum(maximum(abs, surface_amplitudes(p,a)) for a in as)
    (cor=D.restitution(vs,cps,p.We), tc=D.contact_time(ts,cps), maxcp=maximum(cps),
     amp=amp, rej=rej, n=length(ts))
end

cases = [("REFERENCE  M=45 K=2", ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=45,K=2,t_max=25.0), 0.3138, 2.183),
         ("REFERENCE  M=90 K=2", ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=90,K=2,t_max=25.0), 0.3138, 2.183),
         ("Oh=.03 We=5  M=90 K=1", ImpactParams(We=5.0,Bo=0.0189,Oh=0.03,M=90,K=1,t_max=30.0), NaN, 2.047),
         ("Oh=.03 We=10 M=90 K=1", ImpactParams(We=10.0,Bo=0.0189,Oh=0.03,M=90,K=1,t_max=30.0), NaN, 1.955),
         ("Oh=.03 We=10 M=90 K=2", ImpactParams(We=10.0,Bo=0.0189,Oh=0.03,M=90,K=2,t_max=30.0), NaN, 1.955)]
for ex in (false, true)
    @printf("\n===== cp=0 excluded from comparisons: %s =====\n", ex ? "YES (H15)" : "no")
    @printf("%-24s %-9s %-8s %-6s %-8s %-6s %s\n","case","CoR","tc","maxcp","max|z|","rej","target tc")
    for (lab,p,cA,tA) in cases
        r = march(p; excl0 = ex)
        @printf("%-24s %-9s %-8.4f %-6d %-8.3f %-6d %.3f\n", lab,
                isfinite(r.cor) ? @sprintf("%.4f",r.cor) : "NONE", r.tc, r.maxcp, r.amp, r.rej, tA)
    end
end
