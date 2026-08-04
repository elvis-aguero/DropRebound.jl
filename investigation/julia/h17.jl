# H17: does the tangency-optimal cp track the descent at K=1 but not K=2?
using DropSolver, Printf
const D = DropSolver
function probe_run(K; M=90, nst=120)
    p  = ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=M,K=K,t_max=25.0)
    F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
    et_of = (prev,curr,dt,cp) -> begin
        (0 <= cp <= nn) || return (Inf, curr)
        st,nx = D.try_step(p,prev,curr,dt,cp;F0=F0)
        st === :diverge && return (Inf,curr)
        g = [D.gap(p,nx.a,nx.z,p.nodes[i]) for i in 1:nn]
        all(i->g[i]>=0,(cp+1):nn) || return (Inf,nx)
        cp==0 && return (0.0,nx); cp>=nn && return (Inf,nx)
        (abs(g[cp+1]-g[cp]), nx)
    end
    prev,curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
    rows = []
    for step in 1:nst
        curr.t >= p.t_max && break
        # where is the minimum over cp >= 1, right now?
        ets = [et_of(prev,curr,dt,c)[1] for c in 1:min(30,nn)]
        fin = findall(isfinite, ets)
        argm = isempty(fin) ? 0 : fin[argmin(ets[fin])]
        push!(rows, (step, curr.t, curr.z, curr.cp, argm))
        e2,n2 = et_of(prev,curr,dt,curr.cp-1); e3,n3 = et_of(prev,curr,dt,curr.cp)
        e4,n4 = et_of(prev,curr,dt,curr.cp+1)
        nxt = nothing
        if e3 > e4 || e3 > e2
            if e4 <= e2
                e5,_ = et_of(prev,curr,dt,curr.cp+2); e4 < e5 && (nxt = n4)
            else
                e1,_ = et_of(prev,curr,dt,curr.cp-2); e2 < e1 && (nxt = n2)
            end
        elseif isfinite(e3); nxt = n3 end
        if nxt === nothing
            dt /= 2; dt < p.dt_min && break
        else
            prev,curr = curr,nxt; dt = p.dt0
        end
    end
    rows
end
for K in (1,2)
    rows = probe_run(K)
    @printf("\nK=%d : %d steps recorded\n", K, length(rows))
    @printf("  %-6s %-9s %-10s %-5s %s\n","step","t","z","cp","argmin et over cp>=1")
    for r in rows[1:min(12,end)]
        @printf("  %-6d %-9.5f %-10.6f %-5d %d\n", r[1],r[2],r[3],r[4],r[5])
    end
    if length(rows) > 12
        @printf("  ... last: step %d t=%.5f z=%.6f cp=%d argmin=%d\n",
                rows[end][1],rows[end][2],rows[end][3],rows[end][4],rows[end][5])
    end
    tracks = count(r -> r[5] == r[4], rows)
    @printf("  argmin coincides with the incumbent cp in %d of %d steps\n", tracks, length(rows))
end
