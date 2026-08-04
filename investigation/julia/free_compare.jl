# FREE OSCILLATION, Oh = 0, NO CONTACT -- the maxdev-code check with the contact model
# switched off. Their code integrates the modal ODE; ours integrates the variational
# assembly. Both must agree with each other and with Rayleigh's analytic result, and
# any disagreement here would mean the contact debugging had been chasing a symptom.
using DelimitedFiles, Printf, LinearAlgebra
using DropSolver

dir = joinpath(@__DIR__, "lwdr_free")
sigma = 20.5; rho = 0.96; Ro = 0.0203
tu = sqrt(rho * Ro^3 / sigma)

scal = readdlm(joinpath(dir, "scalars.csv"), ',')
defo = readdlm(joinpath(dir, "deformation.csv"), ',')
cpv  = Int.(vec(readdlm(joinpath(dir, "contact.csv"), ',')))

t  = scal[:,1] ./ tu
Z2 = defo[:,2] ./ Ro            # zeta_2(t), non-dimensional
dtv = scal[:,4] ./ tu

@printf("their run: %d states, max contact points = %d (must be 0), dt = %.4g\n",
        length(t), maximum(cpv), dtv[end])
@printf("zeta_2(0) = %.6f   t spans [%.3f, %.3f]\n", Z2[1], t[1], t[end])

# --- analytic: Rayleigh, l = 2, Oh = 0 -------------------------------------------
om = sqrt(2 * 1 * 4)                        # l(l-1)(l+2) = 8
ana(tt) = Z2[1] * cos(om * tt)
err_them = maximum(abs, Z2 .- ana.(t)) / Z2[1]
@printf("\nRayleigh omega_2 = %.8f  (period %.5f)\n", om, 2pi/om)
@printf("THEIR zeta_2 vs analytic:  max rel error = %.3e\n", err_them)

# --- ours: BDF2 on the variational assembly, same dt, K = 1 ----------------------
# K = 1 is the right comparison: their initial condition is a pure surface amplitude
# with no interior structure beyond potential flow, which is exactly what the single
# trial function x^(l+1) represents.
function free_march(l, K, Oh, z0, dt, tend)
    b = ModalBasis([l], K)
    F = assemble_newtonian(b, Oh)
    N = ndof(b)
    # zeta = sum_k a_k, so a pure zeta_2 with no interior structure is a_1 = z0
    a0 = zeros(N); a0[1] = z0
    ts = [0.0]; zs = [z0]
    prev_a = copy(a0); prev_ad = zeros(N)
    curr_a = copy(a0); curr_ad = zeros(N)
    prev_dt = dt
    tt = 0.0
    tr = [1.0 for _ in 1:K]                 # every trial function is 1 at x = 1
    while tt < tend
        r = dt / prev_dt
        c0 = (1 + 2r)/(1 + r); c1 = -(1 + r); c2 = r^2/(1 + r)
        beta = c0/dt
        hv_a  = (c1*curr_a  + c2*prev_a) / dt
        hv_ad = (c1*curr_ad + c2*prev_ad) / dt
        A = beta^2 * F.M + beta * F.C + F.G
        rhs = -F.M * (beta*hv_a + hv_ad) - F.C * hv_a
        a_new = A \ rhs
        ad_new = beta*a_new + hv_a
        prev_a, prev_ad = curr_a, curr_ad
        curr_a, curr_ad = a_new, ad_new
        prev_dt = dt
        tt += dt
        push!(ts, tt); push!(zs, dot(tr, curr_a))
    end
    (ts, zs)
end

for K in (1, 2, 3)
    ts, zs = free_march(2, K, 0.0, Z2[1], dtv[end], t[end])
    e = maximum(abs, zs .- ana.(ts)) / Z2[1]
    @printf("OURS  K=%d  zeta_2 vs analytic:  max rel error = %.3e\n", K, e)
end

# --- direct maxdev-code comparison at matched times -------------------------------
ts, zs = free_march(2, 1, 0.0, Z2[1], dtv[end], t[end])
const DEV = Ref(0.0)
for (i, tt) in enumerate(t)
    j = searchsortedfirst(ts, tt)
    j > length(ts) && break
    DEV[] = max(DEV[], abs(zs[j] - Z2[i]))
end
@printf("\nOURS vs THEIRS directly, K=1: max |dzeta_2| = %.3e  (%.3e relative)\n",
        DEV[], DEV[] / Z2[1])

# --- and with viscosity, against Reid rather than Rayleigh -----------------------
println("\n--- damped, against Reid (independent route) ---")
for Oh in (0.01, 0.05)
    for K in (1, 4)
        b = ModalBasis([2], K)
        lam, om2 = dominant_pair(RitzBasis(2, K), Oh)
        rl, ro, _ = reid_lambda_omega2(Oh, 2)
        @printf("  Oh=%.3f K=%d : lambda %.6f vs Reid %.6f   omega2 %.4f vs %.4f\n",
                Oh, K, lam, rl, om2, ro)
    end
end
