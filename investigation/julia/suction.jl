using DropSolver, Printf
const D = DropSolver
p = ImpactParams(We=5.0, Bo=0.0189, Oh=0.03, M=20, K=2, t_max=6.0)
r = simulate(p)
@printf("run: CoR=%s tc=%.3f maxcp=%d steps=%d\n",
        isfinite(r.cor) ? @sprintf("%.3f", r.cor) : "NONE", r.tc, maximum(r.cp), length(r.t))
worst = 0.0; wstep = 1; wnode = 0; nbad = 0
for i in eachindex(r.t)
    global worst, wstep, wnode, nbad
    c = r.cp[i]
    c < 2 && continue
    vals = [D.pc_at(p, r.pc[i], p.nodes[j]) for j in 1:c]
    mn = minimum(vals); j = argmin(vals)
    if mn < 0
        nbad += 1
        if mn < worst
            worst = mn; wstep = i; wnode = j
        end
    end
end
@printf("steps with negative pressure somewhere in the patch: %d of %d\n", nbad, length(r.t))
if wnode > 0
    @printf("worst %.4g at step %d, node %d of %d in contact, t=%.4f, max|zeta|=%.3f\n",
            worst, wstep, wnode, r.cp[wstep], r.t[wstep],
            maximum(abs, surface_amplitudes(p, r.a[wstep])))
    @printf("interior node (not the edge)? %s\n", wnode < r.cp[wstep] ? "YES" : "no")
end
