# THE PREDICTION. Gate = no penetration outside contact AND p_c >= 0 at the contact
# EDGE. Rank = finite-difference tangency. Probe = +-2 as the dt-too-big test.
# cp = 0 is a sentinel, entered only when cp = 1 fails the gate.
#
# Falsifiable both ways: K=1 must be UNCHANGED (the gate should rarely bind where
# tangency is already non-degenerate), and K=2 must now run.
using DropSolver, Printf
const D = DropSolver

function march(p; gate::Bool, nst = 200000)
    F0 = p.eta_const ? assemble_newtonian(D.basis(p), p.Oh) : nothing
    nn = length(p.nodes)
    function ev(prev, curr, dt, cp)
        (0 <= cp <= nn) && return _ev(prev, curr, dt, cp)
        (Inf, curr, false)
    end
    function _ev(prev, curr, dt, cp)
        st, nx = D.try_step(p, prev, curr, dt, cp; F0 = F0)
        st === :diverge && return (Inf, curr, false)
        g = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
        all(i -> g[i] >= 0, (cp+1):nn) || return (Inf, nx, false)
        ## THE GATE: the film may not pull at the contact edge.
        if gate && cp > 0 && D.pc_at(p, nx.pc, p.nodes[cp]) < 0
            return (Inf, nx, false)
        end
        cp == 0 && return (0.0, nx, true)
        cp >= nn && return (Inf, nx, true)
        (abs(g[cp+1] - g[cp]), nx, true)
    end
    prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
    ts=[0.0]; zs=[curr.z]; vs=[curr.v]; cps=[0]; as=[copy(curr.a)]; rej=0
    for _ in 1:nst
        curr.t >= p.t_max && break
        e2,n2,_ = ev(prev,curr,dt,curr.cp-1)
        e3,n3,f3 = ev(prev,curr,dt,curr.cp)
        e4,n4,_ = ev(prev,curr,dt,curr.cp+1)
        nxt = nothing; ncp = curr.cp
        if e3 > e4 || e3 > e2
            if e4 <= e2
                e5,_,_ = ev(prev,curr,dt,curr.cp+2)
                if e4 < e5; nxt = n4; ncp = curr.cp+1; end
            else
                e1,_,_ = ev(prev,curr,dt,curr.cp-2)
                if e2 < e1; nxt = n2; ncp = curr.cp-1; end
            end
        elseif isfinite(e3)
            nxt = n3
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
    jm = count(i -> abs(cps[i+1]-cps[i]) > 1, 1:length(cps)-1)
    (cor=D.restitution(vs,cps,p.We), tc=D.contact_time(ts,cps), maxcp=maximum(cps),
     amp=amp, rej=rej, n=length(ts), jumps=jm)
end

cases = [
 ("REFERENCE M=90 K=1", ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=90,K=1,t_max=25.0), 0.3138, 2.183),
 ("REFERENCE M=90 K=2", ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=90,K=2,t_max=25.0), 0.3138, 2.183),
 ("REFERENCE M=45 K=2", ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=45,K=2,t_max=25.0), 0.3138, 2.183),
 ("Oh=.03 We=5  K=1",   ImpactParams(We=5.0, Bo=0.0189,Oh=0.03,M=90,K=1,t_max=30.0), NaN, 2.047),
 ("Oh=.03 We=10 K=1",   ImpactParams(We=10.0,Bo=0.0189,Oh=0.03,M=90,K=1,t_max=30.0), NaN, 1.955),
 ("Oh=.30 We=5  K=1",   ImpactParams(We=5.0, Bo=0.0189,Oh=0.30,M=90,K=1,t_max=30.0), NaN, 1.863),
 ("Oh=.30 We=10 K=1",   ImpactParams(We=10.0,Bo=0.0189,Oh=0.30,M=90,K=1,t_max=30.0), NaN, 1.759),
 ("Oh=.03 We=10 K=2",   ImpactParams(We=10.0,Bo=0.0189,Oh=0.03,M=90,K=2,t_max=30.0), NaN, 1.955),
 ("Oh=.30 We=10 K=2",   ImpactParams(We=10.0,Bo=0.0189,Oh=0.30,M=90,K=2,t_max=30.0), NaN, 1.759),
]
for g in (false, true)
    @printf("\n================ edge pressure gate: %s ================\n", g ? "ON" : "off")
    @printf("%-20s %-9s %-9s %-6s %-8s %-6s %-4s %s\n",
            "case","CoR","tc","maxcp","max|z|","rej","j>1","target tc")
    for (lab,p,cA,tA) in cases
        r = march(p; gate = g)
        ok = isfinite(r.cor) && r.tc > 1.0 && r.amp < 2.0
        @printf("%-20s %-9s %-9.4f %-6d %-8.3f %-6d %-4d %-8.3f %s\n", lab,
                isfinite(r.cor) ? @sprintf("%.4f",r.cor) : "NONE",
                r.tc, r.maxcp, r.amp, r.rej, r.jumps, tA, ok ? "" : "<-- FAILS")
    end
end
