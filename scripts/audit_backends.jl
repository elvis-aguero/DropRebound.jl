# A measured account of the three backends: trustworthiness, stability, cost.
#
# Everything here is written to CSV under `outputs/csv/` so that a claim about a
# backend can be traced to the run that produced it rather than to memory. Four
# measurements, each with its own file:
#
#   1. stability   -- completion over a (We, Oh) grid, with the failure mode
#   2. cost        -- wallclock against M, Newtonian and shear-thinning
#   3. newtonian   -- against the 935 impacts of Gabbard et al. (2025)
#   4. thinning    -- against the 72 impacts on the 3000 ppm fluid
#
# Run as: julia --project=. scripts/audit_backends.jl [stage]

using Printf, Dates
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const OUT  = joinpath(@__DIR__, "..", "outputs", "csv")
const DATA = joinpath(@__DIR__, "..", "data")
const STAMP = Dates.format(now(), "yyyy-mm-ddTHH:MM:SS")

const BACKENDS = [
    Backend(formulation = :variational,    contact = :active_set),
    Backend(formulation = :variational,    contact = :lcp),
    Backend(formulation = :nonvariational, contact = :tangency),
]

writerow(io, xs...) = println(io, join(xs, ","))

# ---------------------------------------------------------------- 1. stability
# Does the backend finish, and if not, how does it fail? A backend that is
# accurate on the cases it completes but completes half of them is not usable
# for a sweep, and nothing in the suite measures that.
function stage_stability()
    open(joinpath(OUT, "audit_stability.csv"), "w") do io
        writerow(io, "backend","We","Oh","Bo","M","K","ok","cor","tc","wall_s","rejects","mode")
        for b in BACKENDS, We in (0.1, 0.3, 1.0, 2.0, 3.0), Oh in (0.02, 0.05, 0.1, 0.3, 0.7)
            r = run_impact(b; We = We, Bo = 0.019, Oh = Oh, M = 45, K = 3, t_max = 25.0)
            mode = r.ok ? "ok" :
                   (isempty(r.t)                 ? "threw"      :
                    !isfinite(r.cor)             ? "no_metrics" :
                    r.tc >= 25.0                 ? "no_release" : "other")
            writerow(io, label(b), We, Oh, 0.019, 45, 3, r.ok,
                     r.cor, r.tc, round(r.wall, digits=3), r.diag.rejects, mode)
            @printf("  %-16s We=%.1f Oh=%.2f  %-11s %6.1fs\n",
                    label(b), We, Oh, mode, r.wall)
            flush(io)
        end
    end
end

# -------------------------------------------------------------------- 2. cost
function stage_cost()
    ## The 3000 ppm fluid, Cross fit relabelled as Carreau-Yasuda; same constants as
    ## `validate_shear_thinning_variational.jl`.
    T_CAP = sqrt(1000.0 * 0.0003^3 / 0.0728)
    eta_st = gd -> carreau(gd; lambda_c = 18.48081673111359 / T_CAP,
                           a = 0.7430524574330837, n = 1 - 0.7430524574330837,
                           eta_inf_ratio = 0.0037320997942061666 / 8.433817577956766)
    open(joinpath(OUT, "audit_cost.csv"), "w") do io
        writerow(io, "backend","rheology","M","K","ok","wall_s","cor")
        for b in BACKENDS, M in (14, 20, 30, 45, 60, 90)
            r = run_impact(b; We = 0.5, Bo = 0.019, Oh = 0.0373, M = M, K = 3, t_max = 25.0)
            writerow(io, label(b), "newtonian", M, 3, r.ok, round(r.wall, digits=3), r.cor)
            @printf("  %-16s newtonian M=%2d  %7.1fs  ok=%s\n", label(b), M, r.wall, r.ok)
            flush(io)
        end
        ## Shear thinning is the expensive path; only the variational backends can
        ## evaluate a viscosity field, so the nonvariational one is absent by design.
        for b in BACKENDS[1:2], M in (14, 20, 30, 45)
            r = run_impact(b; We = 0.3643, Bo = 0.012,
                           Oh = 8.433817577956766/sqrt(1000.0*0.0728*0.0003), M = M, K = 3,
                           t_max = 25.0, eta = eta_st)
            writerow(io, label(b), "shear_thinning", M, 3, r.ok, round(r.wall, digits=3), r.cor)
            @printf("  %-16s thinning  M=%2d  %7.1fs  ok=%s\n", label(b), M, r.wall, r.ok)
            flush(io)
        end
    end
end

stage = length(ARGS) >= 1 ? ARGS[1] : "all"
@printf("audit stage=%s  stamp=%s\n", stage, STAMP)
stage in ("all","stability") && (println("[1] stability"); stage_stability())
stage in ("all","cost")      && (println("[2] cost");      stage_cost())
println("done")
