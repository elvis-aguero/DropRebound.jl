# H19: where tangency is degenerate at K=2, does the EDGE PRESSURE discriminate?
# At the correct extent the film pressure should approach zero at the edge; at a
# spuriously small contact it should not.
using DropSolver, Printf
const D = DropSolver
p  = ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=90,K=2,t_max=25.0)
F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
function step_at(prev,curr,dt,cp)
    st,nx = D.try_step(p,prev,curr,dt,cp;F0=F0)
    (st, nx)
end
# march 6 steps with the finite-difference rule to reach the state where it collapses
prev,curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
function fd(prev,curr,dt,cp)
    (0<=cp<=nn) || return (Inf,curr)
    st,nx = step_at(prev,curr,dt,cp); st===:diverge && return (Inf,curr)
    g=[D.gap(p,nx.a,nx.z,p.nodes[i]) for i in 1:nn]
    all(i->g[i]>=0,(cp+1):nn) || return (Inf,nx)
    cp==0 && return (0.0,nx); cp>=nn && return (Inf,nx)
    (abs(g[cp+1]-g[cp]),nx)
end
for _ in 1:6
    global prev,curr,dt
    e2,n2=fd(prev,curr,dt,curr.cp-1); e3,n3=fd(prev,curr,dt,curr.cp); e4,n4=fd(prev,curr,dt,curr.cp+1)
    if e3>e4 || e3>e2
        if e4<=e2
            e5,_=fd(prev,curr,dt,curr.cp+2); e4<e5 ? (prev,curr=curr,n4) : (dt/=2)
        else
            e1,_=fd(prev,curr,dt,curr.cp-2); e2<e1 ? (prev,curr=curr,n2) : (dt/=2)
        end
    elseif isfinite(e3); prev,curr=curr,n3
    else dt/=2 end
end
@printf("state: t=%.5f z=%.6f cp=%d dt=%.3e\n\n", curr.t, curr.z, curr.cp, dt)
@printf("%-4s %-10s %-12s %-12s %-12s %s\n","cp","status","tangency","p_c at edge","p_c at edge-1","p_c peak on patch")
for cp in 1:8
    st,nx = step_at(prev,curr,dt,cp)
    g=[D.gap(p,nx.a,nx.z,p.nodes[i]) for i in 1:nn]
    feas = all(i->g[i]>=0,(cp+1):nn)
    tg = abs(g[cp+1]-g[cp])
    pe  = D.pc_at(p, nx.pc, p.nodes[cp])
    pe1 = cp>1 ? D.pc_at(p, nx.pc, p.nodes[cp-1]) : NaN
    pk  = maximum(abs(D.pc_at(p,nx.pc,p.nodes[i])) for i in 1:cp)
    @printf("%-4d %-10s %-12.4g %-12.4g %-12.4g %.4g\n", cp,
            feas ? string(st) : "infeas", tg, pe, pe1, pk)
end
println("\nIf p_c at the edge is near zero for exactly one cp, it is the missing selector.")
