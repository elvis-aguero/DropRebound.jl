# All five tangency residuals the ancestor's rule consults, so the branch it takes is
# visible rather than inferred.
using DropSolver, Printf
const D = DropSolver
p  = ImpactParams(We=1.0, Bo=0.0189, Oh=0.303767, M=45, K=2, t_max=25.0)
F0 = assemble_newtonian(D.basis(p), p.Oh); nn = length(p.nodes)

function et_of(prev, curr, dt, cp)
    cp < 0 && return (Inf, curr)               # the ancestor: errortan = Inf for cp < 0
    cp > nn && return (Inf, curr)
    st, nx = D.try_step(p, prev, curr, dt, cp; F0=F0)
    st === :diverge && return (Inf, curr)
    gaps = [D.gap(p, nx.a, nx.z, p.nodes[i]) for i in 1:nn]
    all(i -> gaps[i] >= 0, (cp+1):nn) || return (Inf, nx)
    cp == 0 && return (0.0, nx)
    cp >= nn && return (Inf, nx)
    (abs(gaps[cp+1] - gaps[cp]), nx)
end

prev, curr = D.initial_state(p), D.initial_state(p); dt = p.dt0
@printf("step cp |   et(cp-2)    et(cp-1)    et(cp)      et(cp+1)    et(cp+2)  | branch\n")
for step in 1:12
    global prev, curr, dt
    e1,_  = et_of(prev, curr, dt, curr.cp-2)
    e2,n2 = et_of(prev, curr, dt, curr.cp-1)
    e3,n3 = et_of(prev, curr, dt, curr.cp)
    e4,n4 = et_of(prev, curr, dt, curr.cp+1)
    e5,_  = et_of(prev, curr, dt, curr.cp+2)
    br, nxt = if e3 > e4 || e3 > e2
        if e4 <= e2
            e4 < e5 ? ("grow", n4) : ("RECALC: probe says dt too big", nothing)
        else
            e2 < e1 ? ("shrink", n2) : ("RECALC: shrink blocked by et(cp-2)", nothing)
        end
    elseif !isfinite(e3)
        ("RECALC: incumbent Inf", nothing)
    else
        ("hold", n3)
    end
    @printf("%4d %2d | %-11.4g %-11.4g %-11.4g %-11.4g %-11.4g | %s\n",
            step, curr.cp, e1, e2, e3, e4, e5, br)
    if nxt === nothing
        dt /= 2; dt < p.dt_min && (println("  dt floor"); break)
    else
        prev, curr = curr, nxt
    end
end
