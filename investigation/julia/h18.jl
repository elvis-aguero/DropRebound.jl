# H18: an EXACT spectral tangency residual -- |dh/dtheta| at the contact edge, computed
# from the coefficients rather than by a two-node difference -- tracks the descent at
# K=2 where the finite difference collapses.
using DropSolver, Printf
const D = DropSolver

"""Exact dh/dtheta at one angle, from the spectral coefficients."""
function dgap_dtheta(p, a, th)
    z = surface_amplitudes(p, a); ls = p.ls
    ang = [D.legendre_angular(l, cos(th)) for l in ls]
    r   = 1 + sum(z[i]*ang[i].P for i in eachindex(ls))
    drdth = sum(z[i]*ang[i].dPdth for i in eachindex(ls))
    -sin(th)*r + cos(th)*drdth
end

function profiles(K; M=90, nst=90)
    p  = ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=M,K=K,t_max=25.0)
    F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
    ## returns (finite-difference et, exact-slope et, next state)
    function both(prev,curr,dt,cp)
        (0 <= cp <= nn) || return (Inf,Inf,curr)
        st,nx = D.try_step(p,prev,curr,dt,cp;F0=F0)
        st === :diverge && return (Inf,Inf,curr)
        g = [D.gap(p,nx.a,nx.z,p.nodes[i]) for i in 1:nn]
        all(i->g[i]>=0,(cp+1):nn) || return (Inf,Inf,nx)
        cp==0 && return (0.0,Inf,nx)
        cp>=nn && return (Inf,Inf,nx)
        fd = abs(g[cp+1]-g[cp])
        ex = abs(dgap_dtheta(p, nx.a, p.nodes[cp]))
        (fd, ex, nx)
    end
    prev,curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
    rows = []
    for step in 1:nst
        curr.t >= p.t_max && break
        cmax = min(28, nn)
        fds = Float64[]; exs = Float64[]
        for c in 1:cmax
            f,e,_ = both(prev,curr,dt,c); push!(fds,f); push!(exs,e)
        end
        ff = findall(isfinite,fds); fe = findall(isfinite,exs)
        amin_fd = isempty(ff) ? 0 : ff[argmin(fds[ff])]
        amin_ex = isempty(fe) ? 0 : fe[argmin(exs[fe])]
        push!(rows,(step,curr.z,curr.cp,amin_fd,amin_ex))
        ## advance with the FINITE-DIFFERENCE rule (as now) so both are measured on the
        ## same trajectory -- this isolates the residual, not the dynamics
        f2,_,n2 = both(prev,curr,dt,curr.cp-1); f3,_,n3 = both(prev,curr,dt,curr.cp)
        f4,_,n4 = both(prev,curr,dt,curr.cp+1)
        nxt = nothing
        if f3 > f4 || f3 > f2
            if f4 <= f2
                f5,_,_ = both(prev,curr,dt,curr.cp+2); f4 < f5 && (nxt = n4)
            else
                f1,_,_ = both(prev,curr,dt,curr.cp-2); f2 < f1 && (nxt = n2)
            end
        elseif isfinite(f3); nxt = n3 end
        if nxt === nothing
    dt /= 2
    dt < p.dt_min && break
else
    prev, curr = curr, nxt
    dt = p.dt0
end
    end
    rows
end
for K in (1,2)
    rows = profiles(K)
    @printf("\nK=%d\n", K)
@printf("  %-5s %-10s %-4s %-14s %s\n","step","z","cp","argmin FIN.DIFF","argmin EXACT SLOPE")
    for r in rows[1:min(10,end)]
        @printf("  %-5d %-10.6f %-4d %-14d %d\n", r[1],r[2],r[3],r[4],r[5])
    end
    length(rows)>10 && @printf("  ... last: step %d z=%.6f cp=%d  fd=%d  exact=%d\n",
        rows[end][1],rows[end][2],rows[end][3],rows[end][4],rows[end][5])
    @printf("  tracks cp:  finite-difference %d/%d,  EXACT SLOPE %d/%d\n",
            count(r->r[4]==r[3],rows), length(rows), count(r->r[5]==r[3],rows), length(rows))
end
