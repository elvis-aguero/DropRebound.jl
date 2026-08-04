using DropSolver, Printf, LinearAlgebra
const D = DropSolver
BAD = (We = 10.0, Bo = 0.0189, Oh = 0.03, M = 90, K = 1)

# ===== H6 first: is the energy formula trustworthy at all? =====================
# The control I omitted last time. Run the identical audit on a case that WORKS. If the
# budget fails to close there too, the formula is wrong and the H1 verdict was noise.
function audit(p, label)
    r = simulate(p)
    b = D.basis(p); F = assemble_newtonian(b, p.Oh); mass = 4pi/3
    Em(a, ad, z, v) = 0.5*dot(ad, F.M, ad) + 0.5*dot(a, F.G, a) + 0.5*mass*v^2 + mass*p.Bo*z
    function pw(i)
        Q = zeros(ndof(b))
        for j in eachindex(r.pc[i]); Q .+= r.pc[i][j] .* D.force_column(p, j); end
        dot(Q, r.adot[i]) - mass * r.pc[i][2] * r.v[i]
    end
    Wf = 0.0; Dis = 0.0
    for i in 1:length(r.t)-1
        dt = r.t[i+1] - r.t[i]
        Wf  += 0.5*dt*(pw(i) + pw(i+1))
        Dis += 0.5*dt*(dot(r.adot[i],F.C,r.adot[i]) + dot(r.adot[i+1],F.C,r.adot[i+1]))
    end
    dE = Em(r.a[end],r.adot[end],r.z[end],r.v[end]) - Em(r.a[1],r.adot[1],r.z[1],r.v[1])
    rel = abs(dE - (Wf - Dis)) / max(abs(dE), 1.0)
    @printf("  %-22s dE=%+11.4g  Wf-Dis=%+11.4g  closes to %8.3g  %s\n",
            label, dE, Wf-Dis, rel, rel < 0.05 ? "CLOSES" : "FAILS")
    (r, rel)
end
println("=== H6: does the audit close on a case that WORKS? ===")
rw, relw = audit(ImpactParams(We=1.0,Bo=0.0189,Oh=0.303767,M=45,K=2,t_max=25.0), "Oh=0.30 We=1 (works)")
rb, relb = audit(ImpactParams(; BAD..., t_max=30.0), "Oh=0.03 We=10 (breaks)")
println(relw < 0.05 ? "  -> formula VALIDATED; the failure at Oh=0.03 is real" :
                      "  -> formula WRONG; H1's test was uninterpretable (H6 confirmed)")

# ===== H5: does the contact set ever release? =================================
println("\n=== H5: is the contact monotonically growing? ===")
cp = rb.cp
inc = count(i -> cp[i+1] > cp[i], 1:length(cp)-1)
dec = count(i -> cp[i+1] < cp[i], 1:length(cp)-1)
@printf("  cp: max=%d of %d nodes;  increases %d times, DECREASES %d times\n",
        maximum(cp), length(ImpactParams(; BAD...).nodes), inc, dec)
@printf("  first 30 cp values: %s\n", string(cp[1:min(30,end)]))
println(dec == 0 ? "  -> NEVER releases: H5 CONFIRMED" :
                   "  -> it does release, so growth is not one-way: H5 DEAD")

# ===== H4: is dt control implicated? ==========================================
# dt_min = dt0 makes halving impossible: the first rejection ends the run. If the
# amplitude has already blown up by then, dt control cannot be the cause.
println("\n=== H4: dt control, with halving disabled ===")
p4 = ImpactParams(; BAD..., t_max = 30.0, dt_min = ImpactParams(; BAD...).dt0)
r4 = simulate(p4)
z4 = maximum(maximum(abs, surface_amplitudes(p4, a)) for a in r4.a)
@printf("  no halving allowed: steps=%d rejects=%d maxcp=%d max|zeta|=%.3f\n",
        length(r4.t), r4.rejects, maximum(r4.cp), z4)
println(z4 > 2 ? "  -> blows up WITHOUT any dt change: H4 DEAD, dt control exonerated" :
                 "  -> stays bounded, so dt control is implicated: H4 stands")
