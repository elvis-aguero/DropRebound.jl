# Which branch actually fires, counted -- including on REJECTED steps, which every
# earlier trace of mine failed to record.
using DropSolver, Printf
const D = DropSolver
p  = ImpactParams(We=1.0, Bo=0.0189, Oh=0.303767, M=90, K=2, t_max=25.0)
F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
function et_of(prev, curr, dt, cp)
    (0 <= cp <= nn) || return (Inf, curr)
    st, nx = D.try_step(p, prev, curr, dt, cp; F0=F0)
    st === :diverge && return (Inf, curr)
    g = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
    all(i -> g[i] >= 0, (cp+1):nn) || return (Inf, nx)
    cp == 0 && return (0.0, nx)
    cp >= nn && return (Inf, nx)
    (abs(g[cp+1] - g[cp]), nx)
end
counts = Dict("grow"=>0, "hold"=>0, "shrink"=>0, "RJ-grow"=>0, "RJ-shrink"=>0, "RJ-Inf"=>0)
prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
first_shrink_e1 = nothing
for step in 1:400
    global prev, curr, dt, first_shrink_e1
    curr.t >= p.t_max && break
    e2,n2 = et_of(prev,curr,dt,curr.cp-1); e3,n3 = et_of(prev,curr,dt,curr.cp)
    e4,n4 = et_of(prev,curr,dt,curr.cp+1)
    nxt = nothing; br = ""
    if e3 > e4 || e3 > e2
        if e4 <= e2
            e5,_ = et_of(prev,curr,dt,curr.cp+2)
            if e4 < e5; nxt = n4; br = "grow" else br = "RJ-grow" end
        else
            e1,_ = et_of(prev,curr,dt,curr.cp-2)
            if first_shrink_e1 === nothing
                first_shrink_e1 = (step, curr.cp, e1, e2, e3, e4)
            end
            if e2 < e1; nxt = n2; br = "shrink" else br = "RJ-shrink" end
        end
    elseif isfinite(e3); nxt = n3; br = "hold"
    else br = "RJ-Inf" end
    counts[br] += 1
    if nxt === nothing
        dt /= 2; dt < p.dt_min && break
    else
        prev, curr = curr, nxt; dt = p.dt0
    end
end
println("branch counts over the whole march (rejections INCLUDED):")
for k in ("grow","hold","shrink","RJ-grow","RJ-shrink","RJ-Inf")
    @printf("  %-10s %d\n", k, counts[k])
end
if first_shrink_e1 !== nothing
    s,cp,e1,e2,e3,e4 = first_shrink_e1
    @printf("\nfirst time the shrink probe was consulted: step %d, cp=%d\n", s, cp)
    @printf("  et(cp-2)=%.4g  et(cp-1)=%.4g  et(cp)=%.4g  et(cp+1)=%.4g\n", e1,e2,e3,e4)
    @printf("  cp-2 = %d  -> is the probe target cp=0?  %s\n", cp-2, cp-2 == 0 ? "YES" : "NO")
else
    println("\nTHE SHRINK BRANCH IS NEVER ENTERED -- et(0) is never consulted at all.")
end
