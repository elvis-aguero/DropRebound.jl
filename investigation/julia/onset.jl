# How the contact set is BORN in a working simulation -- the question my seeding
# argument turned on, and which the downsampled output could not answer.
using DelimitedFiles, Printf

dir = joinpath(@__DIR__, "lwdr_probe")
sigma = 20.5; rho = 0.96; Ro = 0.0203
tu = sqrt(rho * Ro^3 / sigma); vu = Ro / tu; pu = rho * vu^2

scal = readdlm(joinpath(dir, "scalars.csv"), ',')
pres = readdlm(joinpath(dir, "pressure.csv"), ',')
th   = vec(readdlm(joinpath(dir, "theta.csv"), ','))

t = scal[:,1] ./ tu; z = scal[:,2] ./ Ro; v = scal[:,3] ./ vu
cp = Int.(scal[:,4]); dt = scal[:,5] ./ tu
P = pres ./ pu; M = size(pres,2) - 1

function legP(l, x)
    p0, p1 = 1.0, x
    l == 0 && return p0
    for k in 1:(l-1); p0, p1 = p1, ((2k+1)*x*p1 - k*p0)/(k+1); end
    p1
end
pcf(r, mu) = sum(P[r,j] * legP(j-1, mu) for j in 1:(M+1))

@printf("states=%d   dt in [%.3g, %.3g] time units   saved every ~%.2e\n",
        length(t), minimum(dt), maximum(dt), t[2]-t[1])
@printf("node spacing near the pole: theta_1=%.2f deg, theta_2=%.2f, theta_3=%.2f\n",
        rad2deg(th[1]), rad2deg(th[2]), rad2deg(th[3]))

println("\n--- FIRST 16 SAVED STATES ---")
for i in 1:min(16, length(t))
    nc = cp[i]
    if nc > 0
        f = [pcf(i, cos(th[k])) for k in 1:nc]
        @printf("  i=%2d t=%.5f dt=%.2e z=%.6f v=%8.5f cp=%2d p_1=%9.4g  p_c %8.3g..%8.3g\n",
                i, t[i], dt[i], z[i], v[i], nc, P[i,2], minimum(f), maximum(f))
    else
        @printf("  i=%2d t=%.5f dt=%.2e z=%.6f v=%8.5f cp= 0\n",
                i, t[i], dt[i], z[i], v[i])
    end
end

inc = findfirst(>(0), cp); lastc = findlast(>(0), cp)
nz = cp[cp .> 0]
@printf("\nfirst contact: state %d, t=%.5f, cp=%d\n", inc, t[inc], cp[inc])
@printf("SMALLEST nonzero contact count over the whole run: %d\n", minimum(nz))
@printf("contact counts seen in the first 40 states: %s\n", string(cp[1:min(40,end)]))
@printf("t_c=%.4f   CoR=%.4f\n", t[lastc]-t[inc], abs(v[min(lastc+5,length(v))]/v[max(inc-1,1)]))

# does the pressure peak at onset?
pmax = [cp[i] > 0 ? maximum(pcf(i, cos(th[k])) for k in 1:cp[i]) : 0.0 for i in 1:length(t)]
@printf("\nmax patch pressure over time: at onset=%.3g, global max=%.3g at t=%.4f (cp=%d)\n",
        pmax[inc], maximum(pmax), t[argmax(pmax)], cp[argmax(pmax)])
