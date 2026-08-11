# The variational solver against the shear-thinning experiments.
#
# The fluid is the 3000 ppm polymer solution of derivations/data, characterised
# by a Cross fit and converted to Carreau-Yasuda by the relabelling derived on
# "Cross-Model Fluids": lambda_c = K, a = m, n = 1 - m.
#
# The measured file is in the locale it was exported in -- semicolon separated, comma
# decimal separator -- and its contact time is in SECONDS, so it has to be divided by
# the capillary time sqrt(rho R^3 / sigma) before it can be compared with anything.
# The Weber number and the restitution are already dimensionless.

using Printf
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const DATA = joinpath(@__DIR__, "..", "derivations", "data", "metrics_3000ppm.csv")

## Cross fit of the fluid, from scripts/validate_shear_thinning.jl
const ETA_INF = 0.0037320997942061666   # Pa s
const ETA_0   = 8.433817577956766       # Pa s
const K_CROSS = 18.48081673111359       # s
const M_CROSS = 0.7430524574330837
const BO      = 0.012
const R       = 0.0003                  # m
const SIGMA   = 0.0728                  # N/m
const RHO     = 1000.0                  # kg/m^3, aqueous solution

const T_CAP = sqrt(RHO * R^3 / SIGMA)             # capillary time, s
const OH    = ETA_0 / sqrt(RHO * SIGMA * R)       # zero-shear Ohnesorge
## Cross -> Carreau-Yasuda relabelling
const LAMBDA_C = K_CROSS / T_CAP                  # nondimensionalised by capillary time
const A_CY     = M_CROSS
const N_CY     = 1 - M_CROSS
const ETA_RATIO = ETA_INF / ETA_0

@printf("fluid: Oh_0 = %.4g   Bo = %.4g   t_cap = %.4g s\n", OH, BO, T_CAP)
@printf("Carreau-Yasuda: lambda_c = %.4g (nondim)  a = %.4g  n = %.4g  eta_inf/eta_0 = %.3g\n",
        LAMBDA_C, A_CY, N_CY, ETA_RATIO)

"""Read the European-locale export: semicolons, comma decimals, tc in seconds."""
function read_metrics(path)
    rows = NTuple{3,Float64}[]
    for (i, ln) in enumerate(eachline(path))
        i == 1 && continue                              # header
        f = split(strip(ln), ';')
        length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f[1:3], ',' => '.'))
        any(isnothing, v) && continue
        push!(rows, (v[1], v[2], v[3] / T_CAP))         # We, CoR, tc nondimensional
    end
    rows
end

data = read_metrics(DATA)
@printf("%d experiments; We in [%.4g, %.4g], CoR in [%.3f, %.3f], tc/t_cap in [%.2f, %.2f]\n",
        length(data), minimum(d[1] for d in data), maximum(d[1] for d in data),
        minimum(d[2] for d in data), maximum(d[2] for d in data),
        minimum(d[3] for d in data), maximum(d[3] for d in data))

const M_RUN = 14      # a variable viscosity forces the full coupled reassembly per step
const K_RUN = 2

eta_st = gd -> carreau(gd; lambda_c = LAMBDA_C, a = A_CY, n = N_CY,
                       eta_inf_ratio = ETA_RATIO)

## A Weber sweep rather than one run per experiment: the fluid is fixed, so Oh and Bo
## are the same for every point and only We varies.
Wes = sort(unique(round.([d[1] for d in data]; digits = 3)))
grid = exp.(range(log(minimum(Wes)), log(maximum(Wes)); length = 8))

@printf("\n%-8s | %-26s | %-26s | %s\n", "We", "shear-thinning", "Newtonian (Oh_0)", "experiment")
res = Any[]
for We in grid
    rst = simulate(ImpactParams(We = We, Bo = BO, Oh = OH, M = M_RUN, K = K_RUN,
                                eta = eta_st, t_max = 25.0))
    rnw = simulate(ImpactParams(We = We, Bo = BO, Oh = OH, M = M_RUN, K = K_RUN,
                                t_max = 25.0))
    near = [d for d in data if abs(log(d[1]) - log(We)) < 0.30]
    ec = isempty(near) ? NaN : mean(d[2] for d in near)
    et = isempty(near) ? NaN : mean(d[3] for d in near)
    push!(res, (We, rst.cor, rst.tc, rnw.cor, rnw.tc, ec, et, length(near)))
    @printf("%-8.4g | CoR %.4f tc %8.3f | CoR %.4f tc %8.3f | CoR %.4f tc %6.3f (n=%d)\n",
            We, rst.cor, rst.tc, rnw.cor, rnw.tc, ec, et, length(near))
end

println("\n================ SUMMARY ================")
for (name, ic, it) in (("shear-thinning", 2, 3), ("Newtonian", 4, 5))
    okc = [r for r in res if isfinite(r[6]) && isfinite(r[ic])]
    okt = [r for r in res if isfinite(r[7]) && isfinite(r[it])]
    if isempty(okc)
    @printf("%-15s no run completed -- see below\n", name)
    continue
end
rc = [abs(r[ic] - r[6])/r[6] for r in okc]
    rt = [abs(r[it] - r[7])/r[7] for r in okt]
    sc = [(r[ic] - r[6])/r[6] for r in okc]
    @printf("%-15s CoR median |err| %5.1f%%  signed %+6.1f%%  (n=%d)   tc median |err| %5.1f%%  (n=%d)\n",
            name, 100*median(rc), 100*median(sc), length(okc),
            isempty(rt) ? NaN : 100*median(rt), length(okt))
end
