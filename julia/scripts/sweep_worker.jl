# One shard of the concentration sweep. The driver splits the (fluid, Weber) work list
# across processes; this runs the shard whose index it is given.
#
# WHY PROCESSES RATHER THAN THREADS. The coupled assembly streams the whole geometry
# cache once per Picard sweep and is memory-bandwidth bound, so threading inside one
# solve saturates early: eight threads buy about 2x. Eight processes each doing a
# different Weber number buy close to 8x, because each is a separate stream and the
# machine has bandwidth to spare. The cost is one 574 MB cache per process, which is
# why the driver sizes the pool against RAM as well as against cores.

using Printf, Dates
using DropSolver

const SHARD  = parse(Int, ARGS[1])
const NSHARD = parse(Int, ARGS[2])
const OUTCSV = ARGS[3]

const R_DROP, SIGMA, RHO, BO = 0.0003, 0.0728, 1000.0, 0.012
const T_CAP = sqrt(RHO * R_DROP^3 / SIGMA)
const M_RUN, K_RUN = DropSolver.DEFAULT_M, DropSolver.DEFAULT_K

const FLUIDS = [
 (ppm =  300, file = "metrics_300ppm.csv",   eta_inf = 0.002044613704709438,
  eta_0 = 0.07434592880102645, k = 0.861266870703513,  n = 0.6688965094281969),
 (ppm = 1000, file = "metrics_1000ppm.csv",  eta_inf = 0.002747884967209792,
  eta_0 = 0.2601531742575956,  k = 1.076060672638565,  n = 0.6724627161400011),
 (ppm = 2000, file = "metrics_2000ppm.csv",  eta_inf = 0.00501454294688455,
  eta_0 = 1.466021351677301,   k = 2.8234142918640557, n = 0.7640580520083537),
 (ppm = 3000, file = "metrics_3000ppm_2.csv", eta_inf = 0.0037320997942061666,
  eta_0 = 8.433817577956766,   k = 18.48081673111359,  n = 0.7430524574330837),
 (ppm =    0, file = "metrics_Water.csv",     eta_inf = 1.0e-3,
  eta_0 = 1.0e-3,              k = 0.0,               n = 1.0),
]

function read_metrics(path)
    lines = readlines(path); hdr = split(strip(lines[1]), ';')
    function idx(pat)
        i = findfirst(h -> occursin(pat, lowercase(h)), hdr)
        i === nothing && error("no column matching $pat in $(basename(path)): $hdr")
        i
    end
    iw, ie, it = idx("we"), idx("epsilon"), idx("tc")
    rows = NTuple{3,Float64}[]
    for ln in lines[2:end]
        f = split(strip(ln), ';'); length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f, ',' => '.')); any(isnothing, v) && continue
        push!(rows, (v[iw], v[ie], v[it] / T_CAP))
    end
    rows
end

## the full work list, identical in every shard so the partition is consistent
const DATA = joinpath(@__DIR__, "..", "data")
work = NamedTuple[]
for fl in FLUIDS
    rows = read_metrics(joinpath(DATA, fl.file))
    Wes = exp10.(range(log10(minimum(r[1] for r in rows)),
                       log10(maximum(r[1] for r in rows)); length = 9))
    for We in Wes
        push!(work, (fl = fl, We = We, Oh_0 = fl.eta_0/sqrt(RHO*SIGMA*R_DROP)))
    end
end

open(OUTCSV, "w") do io
    for (i, w) in enumerate(work)
        (i - 1) % NSHARD == SHARD || continue
        fl = w.fl; newt = fl.k == 0.0
        eta_fn = newt ? nothing :
                 gd -> carreau(gd; lambda_c = fl.k/T_CAP, a = fl.n, n = 1 - fl.n,
                               eta_inf_ratio = fl.eta_inf/fl.eta_0)
        t0 = time()
        local m, ok
        try
            p = newt ? ImpactParams(We = w.We, Bo = BO, Oh = w.Oh_0, M = M_RUN, K = K_RUN,
                                    t_max = 25.0) :
                       ImpactParams(We = w.We, Bo = BO, Oh = w.Oh_0, M = M_RUN, K = K_RUN,
                                    eta = eta_fn, t_max = 25.0)
            r = simulate_lcp(p); m = proximity_metrics(p, r)
            ok = !isempty(r.t) && r.t[end] > 2.0 && isfinite(m.cor)
        catch
            m = (cor = NaN, tc = NaN); ok = false
        end
        println(io, join((fl.ppm, "simulation", w.We, m.cor, m.tc, w.Oh_0, ok,
                          round(time()-t0, digits=1)), ","))
        flush(io)
        @printf("[shard %d] ppm=%-4d We=%.4f cor=%.4f %s %.0fs\n",
                SHARD, fl.ppm, w.We, m.cor, ok ? "" : "REJECTED", time()-t0); flush(stdout)
    end
end
