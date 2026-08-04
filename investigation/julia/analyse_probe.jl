# Read the Gabbard et al. probe run and answer the question the guessing was about:
# what does the film pressure actually look like at contact onset?
using DelimitedFiles, Printf

dir = joinpath(@__DIR__, "lwdr_probe")
meta = Dict{String,Float64}()
for ln in eachline(joinpath(dir, "meta.txt"))
    k, v = split(strip(ln))
    meta[k] = parse(Float64, v)
end
M = Int(meta["M"]); Bo = meta["Bo"]; Oh = meta["Oh"]; We = meta["We"]

sigma = 20.5; rho = 0.96; Ro = 0.0203
tu = sqrt(rho * Ro^3 / sigma)
vu = Ro / tu
pu = rho * vu^2                       # = sigma/Ro, the capillary pressure scale

scal = readdlm(joinpath(dir, "scalars.csv"), ',')
pres = readdlm(joinpath(dir, "pressure.csv"), ',')
defo = readdlm(joinpath(dir, "deformation.csv"), ',')
th   = vec(readdlm(joinpath(dir, "theta.csv"), ','))

t  = scal[:, 1] ./ tu
z  = scal[:, 2] ./ Ro
v  = scal[:, 3] ./ vu
cp = Int.(scal[:, 4])
P  = pres ./ pu                       # column j = harmonic l = j-1
Z  = defo ./ Ro                       # column j = harmonic l = j

@printf("We=%.3f  Oh=%.4f  Bo=%.4f  M=%d  angles=%d\n", We, Oh, Bo, M, length(th))
@printf("states=%d   t in [%.3f, %.3f]  (units of sqrt(rho R^3/sigma))\n",
        length(t), t[1], t[end])

# --- KPIs ------------------------------------------------------------------------
inc = findfirst(>(0), cp)
if inc === nothing
    println("NO CONTACT in this run")
else
    lastc = findlast(>(0), cp)
    @printf("contact: step %d (t=%.4f) -> step %d (t=%.4f)   t_c = %.4f\n",
            inc, t[inc], lastc, t[lastc], t[lastc] - t[inc])
    @printf("max contact points = %d  of %d angles (theta_c = %.2f deg)\n",
            maximum(cp), length(th), rad2deg(th[maximum(cp)]))
    vin  = v[inc-1]
    vout = v[min(lastc + 3, length(v))]
    @printf("v_in = %.4f   v_out = %.4f   CoR = %.4f\n", vin, vout, abs(vout/vin))
end

# --- the pressure spectrum at onset and at peak ----------------------------------
legP(l, x) = begin
    p0, p1 = 1.0, x
    l == 0 && return p0
    for k in 1:(l-1); p0, p1 = p1, ((2k+1)*x*p1 - k*p0)/(k+1); end
    p1
end
pc_field(row, mu) = sum(P[row, j] * legP(j-1, mu) for j in 1:(M+1))

if inc !== nothing
    peak = inc - 1 + argmax(cp[inc:findlast(>(0), cp)])
    for (label, row) in (("onset  ", inc), ("onset+1", inc+1), ("peak   ", peak))
        row > length(t) && continue
        @printf("\n--- %s : step %d  t=%.4f  cp=%d  z=%.4f  v=%.4f\n",
                label, row, t[row], cp[row], z[row], v[row])
        sp = P[row, :]
        @printf("  p_l for l=0..8 : %s\n",
                join((@sprintf("%.4g", x) for x in sp[1:9]), "  "))
        @printf("  |p_l| max=%.4g at l=%d ; |p_90|/|p_max| = %.2e\n",
                maximum(abs, sp), argmax(abs.(sp)) - 1,
                abs(sp[end]) / maximum(abs, sp))
        # reconstructed field on and off the patch
        nc = max(cp[row], 1)
        mus_in  = [cos(th[i]) for i in 1:nc]
        f_in    = [pc_field(row, m) for m in mus_in]
        @printf("  p_c on the patch: min=%.4g  max=%.4g   (positive everywhere: %s)\n",
                minimum(f_in), maximum(f_in), all(>(0), f_in))
        outi = (nc+1):min(nc+6, length(th))
        @printf("  p_c just outside: %s\n",
                join((@sprintf("%.2g", pc_field(row, cos(th[i]))) for i in outi), " "))
    end
end

# --- trajectory ------------------------------------------------------------------
println("\n--- trajectory ---")
step = max(1, length(t) ÷ 22)
for i in 1:step:length(t)
    @printf("  t=%7.4f  z=%8.5f  v=%8.4f  cp=%3d  p_1=%10.4g\n",
            t[i], z[i], v[i], cp[i], P[i, 2])
end
