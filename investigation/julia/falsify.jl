using DropSolver, Printf, LinearAlgebra
const D = DropSolver

# ---- H3: is the SELECTOR the causal difference?  My code at K=1, M=90 matches the
#          ancestor in every respect except how the contact count is chosen.
println("=== H3: my code at K=1, M=90 (ancestor-matched except the selector) ===")
for (Oh, We) in ((0.03, 10.0), (0.30, 10.0))
    p = ImpactParams(We = We, Bo = 0.0189, Oh = Oh, M = 90, K = 1, t_max = 30.0)
    r = simulate(p)
    z = maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a)
    @printf("  Oh=%.2f We=%.1f K=1 M=90 : CoR=%-8s tc=%-6.3f maxcp=%-3d max|zeta|=%.3f\n",
            Oh, We, isfinite(r.cor) ? @sprintf("%.4f", r.cor) : "NONE", r.tc,
            maximum(r.cp), z)
end
println("  ancestor, same cases      : released, tc=1.955 |z|=1.216 | tc=1.759 |z|=0.836")

# ---- H1: is the energy gain accounted for by the FILM's work, or does the budget
#          fail to close at all?  The second would mean a numerical leak, not a pulling
#          film, and would kill H1 outright.
println("\n=== H1: energy audit THROUGH contact (Oh=0.03, We=5, M=45, K=2) ===")
p = ImpactParams(We = 5.0, Bo = 0.0189, Oh = 0.03, M = 45, K = 2, t_max = 8.0)
r = simulate(p)
b = basis(p); F = assemble_newtonian(b, p.Oh); mass = 4pi/3
Emech(a, ad, z, v) = 0.5*dot(ad, F.M, ad) + 0.5*dot(a, F.G, a) +
                     0.5*mass*v^2 + mass*p.Bo*z
function film_power(i)
    Q = zeros(ndof(b))
    for j in eachindex(r.pc[i]); Q .+= r.pc[i][j] .* D.force_column(p, j); end
    dot(Q, r.adot[i]) + (-(mass) * r.pc[i][2]) * r.v[i]   # shape work + CoM work
end
E0 = Emech(r.a[1], r.adot[1], r.z[1], r.v[1])
Wf = 0.0; Dis = 0.0
for i in 1:length(r.t)-1
    dt = r.t[i+1] - r.t[i]
    global Wf  += 0.5*dt*(film_power(i) + film_power(i+1))
    global Dis += 0.5*dt*(dot(r.adot[i], F.C, r.adot[i]) + dot(r.adot[i+1], F.C, r.adot[i+1]))
end
E1 = Emech(r.a[end], r.adot[end], r.z[end], r.v[end])
@printf("  E(0)=%.4g   E(end)=%.4g   dE=%+.4g\n", E0, E1, E1-E0)
@printf("  film work   = %+.4g   dissipated = %.4g\n", Wf, Dis)
@printf("  predicted dE = Wf - Dis = %+.4g   ACTUAL dE = %+.4g\n", Wf-Dis, E1-E0)
res = abs((E1-E0) - (Wf-Dis)) / max(abs(E1-E0), 1.0)
@printf("  budget closes to %.3g relative  ->  %s\n", res,
        res < 0.05 ? "CLOSES: the film really is the source (H1 stands)" :
                     "FAILS TO CLOSE: energy from outside the equations (H1 DEAD)")
@printf("  did the film do NET POSITIVE work? %s\n", Wf > 0 ? "YES" : "NO -- H1 dead")

# ---- H2: are the tangency residuals for cp=0 and cp>0 even on the same scale?
println("\n=== H2: errortan scale across candidates, early contact ===")
p2 = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 45, K = 2, t_max = 25.0)
s0 = D.initial_state(p2); prev, curr = s0, s0
for _ in 1:6
    st, nx = D.try_step(p2, prev, curr, p2.dt0, curr.cp;
                        F0 = assemble_newtonian(basis(p2), p2.Oh))
    st === :ok || break
    global prev, curr = curr, nx
end
@printf("  at t=%.4f, cp=%d:\n", curr.t, curr.cp)
for cand in 0:3
    st, nx = D.try_step(p2, prev, curr, p2.dt0, cand;
                        F0 = assemble_newtonian(basis(p2), p2.Oh))
    et = if cand == 0
        abs(D.gap(p2, nx.a, nx.z, p2.nodes[1]))
    else
        abs(D.gap(p2, nx.a, nx.z, p2.nodes[cand+1]) - D.gap(p2, nx.a, nx.z, p2.nodes[cand]))
    end
    @printf("    cand=%d  status=%-12s errortan=%.4e\n", cand, st, et)
end
println("  if cp=0 is orders of magnitude smaller, it always wins -> H2 stands")
