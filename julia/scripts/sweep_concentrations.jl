# The shear-thinning series against experiment, one sweep per concentration.
#
# Each fluid is a Cross fit, relabelled to Carreau-Yasuda by lambda_c = k, a = m,
# n_CY = 1 - m, the correspondence derived on "Cross-Model Fluids". Only the fluid
# changes between sweeps: R, sigma, rho and Bo are shared, so the concentration enters
# through the zero-shear Ohnesorge number and the shape of eta(gammadot).
#
# The CSVs are a European-locale export: semicolons, comma decimals, tc in seconds.
# Their column ORDER is not stable -- the 2026 files are `We;tc;epsilon` and the older
# 3000 ppm file is `We;epsilon;tc` -- so columns are located by header name. Reading by
# position would silently swap restitution and contact time for one of the two.

using Printf, Dates
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const DATA = joinpath(@__DIR__, "..", "data")
const OUT  = joinpath(@__DIR__, "..", "..", "results")

## shared drop and substrate
const R_DROP = 0.0003        # m
const SIGMA  = 0.0728        # N/m
const RHO    = 1000.0        # kg/m^3
const BO     = 0.012
const T_CAP  = sqrt(RHO * R_DROP^3 / SIGMA)

"""Cross fits. `n` here is the Cross exponent `m`, so the Carreau-Yasuda pair is
`a = n`, `n_CY = 1 - n`."""
const FLUIDS = [
 (ppm =  300, file = "metrics_300ppm.csv",   eta_inf = 0.002044613704709438,
  eta_0 = 0.07434592880102645, k = 0.861266870703513,  n = 0.6688965094281969),
 (ppm = 1000, file = "metrics_1000ppm.csv",  eta_inf = 0.002747884967209792,
  eta_0 = 0.2601531742575956,  k = 1.076060672638565,  n = 0.6724627161400011),
 (ppm = 2000, file = "metrics_2000ppm.csv",  eta_inf = 0.00501454294688455,
  eta_0 = 1.466021351677301,   k = 2.8234142918640557, n = 0.7640580520083537),
 (ppm = 3000, file = "metrics_3000ppm_2.csv", eta_inf = 0.0037320997942061666,
  eta_0 = 8.433817577956766,   k = 18.48081673111359,  n = 0.7430524574330837),
 ## Water is the Newtonian end of the series: no thinning, so no Cross fit. At 20 C
 ## its dynamic viscosity is 1.0e-3 Pa s, which gives Oh = 6.8e-3 here. It is the
 ## control the shear-thinning fluids are read against.
 (ppm =    0, file = "metrics_Water.csv",     eta_inf = 1.0e-3,
  eta_0 = 1.0e-3,              k = 0.0,               n = 1.0),
]

"""Read a locale export, locating `We`, `epsilon` and `tc` by HEADER."""
function read_metrics(path)
    lines = readlines(path)
    hdr = split(strip(lines[1]), ';')
    ## `something(x, error(...))` would throw unconditionally: `something` is a
    ## function, so its second argument is evaluated whether or not the first is nothing.
    function idx(pat)
        i = findfirst(h -> occursin(pat, lowercase(h)), hdr)
        i === nothing && error("no column matching $pat in $(basename(path)): $hdr")
        i
    end
    iw, ie, it = idx("we"), idx("epsilon"), idx("tc")
    rows = NTuple{3,Float64}[]
    for ln in lines[2:end]
        f = split(strip(ln), ';'); length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f, ',' => '.'))
        any(isnothing, v) && continue
        push!(rows, (v[iw], v[ie], v[it] / T_CAP))     # We, restitution, tc in capillary times
    end
    rows
end

const M_RUN = parse(Int, get(ENV, "SWEEP_M", string(DropSolver.DEFAULT_M)))
const K_RUN = DropSolver.DEFAULT_K

## MEMORY. The coupled geometry cache is ndof^2 * nq * 8 bytes and is built once per
## basis, then shared read-only by every sweep and every thread. At the production
## truncation that is 574 MB for a 2 GB budget, so it is cached rather than falling
## back to full quadrature. Report it before spending an hour on the assumption.
let b = ModalBasis(2:M_RUN, K_RUN, :legendre)
    by = DropSolver.coupled_cache_bytes(DropSolver.ndof(b), 40, 48)
    @printf("M = %d, K = %d, ndof = %d, threads = %d\n",
            M_RUN, K_RUN, DropSolver.ndof(b), Threads.nthreads())
    @printf("geometry cache %.0f MB of a %d MB budget -> %s\n", by/1e6,
            DropSolver.COUPLED_CACHE_BUDGET ÷ 1_000_000,
            by <= DropSolver.COUPLED_CACHE_BUDGET ? "cached" : "NOT cached, slow path")
    flush(stdout)
end

open(joinpath(OUT, "sweep_concentrations.csv"), "w") do io
    println(io, "ppm,kind,We,cor,tc,Oh_0,ok,wall_s")
    for fl in FLUIDS
        exp_rows = read_metrics(joinpath(DATA, fl.file))
        Oh_0     = fl.eta_0 / sqrt(RHO * SIGMA * R_DROP)
        lambda_c = fl.k / T_CAP
        newtonian = fl.k == 0.0
        eta_fn   = newtonian ? nothing :
                   gd -> carreau(gd; lambda_c = lambda_c, a = fl.n,
                                 n = 1 - fl.n, eta_inf_ratio = fl.eta_inf / fl.eta_0)
        @printf("%s: %d experiments, Oh_0 = %.3g%s\n",
                fl.ppm == 0 ? "water" : "$(fl.ppm) ppm", length(exp_rows), Oh_0,
                newtonian ? " (Newtonian)" :
                @sprintf(", lambda_c = %.4g, a = %.4g", lambda_c, fl.n)); flush(stdout)
        for (We, eps, tc) in exp_rows
            println(io, join((fl.ppm, "experiment", We, eps, tc, Oh_0, true, 0.0), ","))
        end
        ## the simulation curve, spanning the measured Weber range
        Wes = exp10.(range(log10(minimum(r[1] for r in exp_rows)),
                           log10(maximum(r[1] for r in exp_rows)); length = 9))
        for We in Wes
            t0 = time()
            p = newtonian ?
                ImpactParams(We = We, Bo = BO, Oh = Oh_0, M = M_RUN, K = K_RUN, t_max = 25.0) :
                ImpactParams(We = We, Bo = BO, Oh = Oh_0, M = M_RUN, K = K_RUN,
                             eta = eta_fn, t_max = 25.0)
            local m, ok
            try
                r = simulate_lcp(p); m = proximity_metrics(p, r)
                ok = !isempty(r.t) && r.t[end] > 2.0 && isfinite(m.cor)
            catch
                m = (cor = NaN, tc = NaN); ok = false
            end
            w = time() - t0
            println(io, join((fl.ppm, "simulation", We, m.cor, m.tc, Oh_0, ok,
                              round(w, digits=1)), ","))
            @printf("   We=%.4f  cor=%.4f  tc=%.3f  %s  %.0fs\n",
                    We, m.cor, m.tc, ok ? "" : "REJECTED", w); flush(stdout)
            flush(io)
        end
    end
end
println("wrote results/sweep_concentrations.csv")
