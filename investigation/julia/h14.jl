# H14: at the state where K=2 stalls, is the failure :diverge (solve blew up) or
# :penetrate (geometry)?  And what is the KKT conditioning?
using DropSolver, Printf, LinearAlgebra
const D = DropSolver
p  = ImpactParams(We=1.0, Bo=0.0189, Oh=0.303767, M=90, K=2, t_max=25.0)
F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)
function et_of(prev, curr, dt, cp)
    (0 <= cp <= nn) || return (:oob, Inf, curr)
    st, nx = D.try_step(p, prev, curr, dt, cp; F0=F0)
    st === :diverge && return (:diverge, Inf, nx)
    g = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
    all(i -> g[i] >= 0, (cp+1):nn) || return (:penetrate, Inf, nx)
    cp == 0 && return (:ok, 0.0, nx)
    (:ok, abs(g[cp+1] - g[cp]), nx)
end
# march to the stall
prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
for _ in 1:17
    global prev, curr, dt
    _,e2,n2 = et_of(prev,curr,dt,curr.cp-1); _,e3,n3 = et_of(prev,curr,dt,curr.cp)
    _,e4,n4 = et_of(prev,curr,dt,curr.cp+1)
    if e3 > e4 || e3 > e2
        _,e5,_ = et_of(prev,curr,dt,curr.cp+2)
        if e4 <= e2 && e4 < e5; prev,curr = curr,n4
        else
            _,e1,_ = et_of(prev,curr,dt,curr.cp-2)
            e2 < e1 ? (prev,curr = curr,n2) : (dt /= 2)
        end
    elseif isfinite(e3); prev,curr = curr,n3
    else dt /= 2 end
end
@printf("stalled at t=%.5f z=%.6f cp=%d dt=%.3e\n", curr.t, curr.z, curr.cp, dt)
@printf("max|zeta| here = %.4g\n", maximum(abs, surface_amplitudes(p, curr.a)))
println("\ncandidate statuses at this state, and at successively halved dt:")
for (lab, d) in (("dt", dt), ("dt/4", dt/4), ("dt/64", dt/64), ("dt/4096", dt/4096))
    @printf("  %-8s ", lab)
    for c in (curr.cp-1, curr.cp, curr.cp+1, curr.cp+2)
        s, e, _ = et_of(prev, curr, d, c)
        @printf("c=%d:%-10s et=%-10.3g ", c, string(s), e)
    end
    println()
end
# conditioning of the KKT block at this state
b = D.basis(p); N = ndof(b); npc = D.pc_len(p)
β = 1.5/dt
A = β^2*F0.M + β*F0.C + F0.G
Qm = hcat((D.force_column(p, j) for j in 1:npc)...)
Hm = zeros(npc, N); Zm = zeros(npc, npc)
for i in 1:npc
    if i <= curr.cp
        row, _ = D.gap_row(p, p.nodes[i]); Hm[i,:] = row; Zm[i,2] += -1/β^2
    else
        for j in 1:npc; Zm[i,j] = D.legendre_angular(D.pc_l(j), cos(p.nodes[i])).P; end
    end
end
KKT = [A -Qm; Hm Zm]
@printf("\ncond(A)=%.3g   cond(KKT)=%.3g   cond(M)=%.3g\n",
        cond(A), cond(KKT), cond(F0.M))
