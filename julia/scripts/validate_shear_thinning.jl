#!/usr/bin/env julia
# Validates the non-perturbative Carreau-Yasuda extension (STExactParams,
# julia/src/st_exact_extension.jl) against real shear-thinning experimental
# data: julia/derivations/data/metrics_3000ppm.csv.
#
# Target (user-specified): sample 20 RANDOM experiments (not cherry-picked),
# median relative error on both CoR and contact time should be <= 20%.
#
# Usage:
#   julia --project=.. scripts/validate_shear_thinning.jl [n_samples] [seed]

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))
using Random

using Logging
global_logger(ConsoleLogger(stderr, Logging.Warn))

# ---------------------------------------------------------------------------
# Fluid characterization (Cross model fit, converted to Carreau-Yasuda via
# lambda_c=K, a=m, n=1-m -- see julia/derivations/cross_fluid_derivation.jl).
#
# `eps_ST` is the Carreau-Yasuda SHAPE EXPONENT, (1-n)/a -- see the law in
# `STExactParams`'s docstring:
#     mu_eff/mu_0 = eta_inf_ratio + (1-eta_inf_ratio)*[1+(lambda_c*gammadot)^a]^(-eps_ST)
# The thinning AMPLITUDE Delta=(eta_0-eta_inf)/eta_0 is already carried by the
# separate `(1-eta_inf_ratio)` factor, so it must NOT also be folded into the
# exponent. This script previously passed `Delta` here, conflating the two;
# harmless only because Delta ~ 0.9996 for this particular fluid, but for a
# fluid with Delta = 0.5 the exponent would have been wrong by a factor of two.
# Derived rather than hardcoded so the next fluid inherits the fix. See
# `julia/derivations/generalized_newtonian_hierarchy_derivation.jl` Step 9.
# ---------------------------------------------------------------------------
const ETA_INF = 0.0037320997942061666   # Pa s
const ETA_0   = 8.433817577956766       # Pa s
const K_CROSS = 18.48081673111359       # s
const M_CROSS = 0.7430524574330837
const BO      = 0.012
const R       = 0.0003                  # m
const SIGMA   = 0.0728                  # N/m
const G       = 9.81                    # m/s^2

const RHO      = BO * SIGMA / (G * R^2)
const T_SIGMA  = sqrt(RHO * R^3 / SIGMA)
const OH0      = ETA_0 / sqrt(RHO * SIGMA * R)
const LAMBDA_C = K_CROSS / T_SIGMA        # dimensionless
const A_SHAPE  = M_CROSS
const N_CY     = 1 - M_CROSS              # Cross -> CY power-law index
const EPS_ST   = (1 - N_CY) / A_SHAPE     # == 1.0 exactly under this mapping
const ETA_INF_RATIO = ETA_INF / ETA_0     # physical floor on mu_eff/mu_0 (and hence Oh_eff/Oh0)

@assert isapprox(EPS_ST, 1.0; atol=1e-12) "Cross->CY mapping must give exponent 1"

@info "Fluid parameters" RHO T_SIGMA OH0 LAMBDA_C A_SHAPE EPS_ST ETA_INF_RATIO

# ---------------------------------------------------------------------------
# Load the CSV (semicolon-separated, comma-decimal: "0,481266206;0,818908412;0,0017")
# ---------------------------------------------------------------------------
function load_experiments(path)
    rows = NamedTuple{(:We, :epsilon, :tc), Tuple{Float64,Float64,Float64}}[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue   # header
        isempty(strip(line)) && continue
        parts = split(line, ';')
        length(parts) < 3 && continue
        We  = parse(Float64, replace(parts[1], ',' => '.'))
        eps = parse(Float64, replace(parts[2], ',' => '.'))
        tc  = parse(Float64, replace(parts[3], ',' => '.'))
        push!(rows, (We=We, epsilon=eps, tc=tc))
    end
    rows
end

const DATA_PATH = joinpath(@__DIR__, "..", "derivations", "data", "metrics_3000ppm.csv")
const EXPERIMENTS = load_experiments(DATA_PATH)
@info "Loaded experiments" n=length(EXPERIMENTS)

# ---------------------------------------------------------------------------
# Run one experiment through the solver
# ---------------------------------------------------------------------------
function predict(We::Float64, M::Int, stx::STExactParams)
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    ## `viscous=:reid` to match the `stx` built in main(). These coefficients are
    ## dead on the shear-thinning path (build_residual_st_exact! overwrites every
    ## entry that reads them), but leaving the default `:lamb` here is a latent
    ## trap: anything later added to the residual that reads cfg.drop_lambda
    ## outside the R2 block would silently mix Lamb and Reid coefficients.
    cfg       = SimConstants(M, M + 1, OH0, BO, theta_vec, precomp, dt_max; viscous=:reid)

    init = DropState(M)
    init.z = 1.05
    init.v = -sqrt(We)
    init.dt = dt_max
    init.cp = 0

    times, states = solve_drop!(cfg, OBParams(), init; stx=stx,
        t_end=15.0, save_every=dt_max / 4)
    kpi = extract_kpis(times, states, cfg)
    kpi.cor, kpi.contact_time * T_SIGMA
end

function main(n_samples::Int, seed::Int; M::Int=12)
    rng = MersenneTwister(seed)
    idx = randperm(rng, length(EXPERIMENTS))[1:n_samples]
    sample = EXPERIMENTS[idx]

    stx = STExactParams(M, OH0, LAMBDA_C, A_SHAPE, EPS_ST; viscous=:reid, eta_inf_ratio=ETA_INF_RATIO)

    cor_errs = Float64[]
    tc_errs  = Float64[]
    println()
    println(rpad("We", 10), rpad("eps_meas", 10), rpad("eps_pred", 10),
             rpad("tc_meas(s)", 12), rpad("tc_pred(s)", 12), "cor_err  tc_err")
    for row in sample
        cor_pred, tc_pred = try
            predict(row.We, M, stx)
        catch e
            @warn "solve failed" We = row.We exception = e
            (NaN, NaN)
        end
        cor_err = isfinite(cor_pred) ? abs(cor_pred - row.epsilon) / row.epsilon : NaN
        tc_err  = isfinite(tc_pred)  ? abs(tc_pred - row.tc) / row.tc : NaN
        isfinite(cor_err) && push!(cor_errs, cor_err)
        isfinite(tc_err) && push!(tc_errs, tc_err)
        println(rpad(round(row.We, digits=4), 10), rpad(round(row.epsilon, digits=4), 10),
                 rpad(round(cor_pred, digits=4), 10), rpad(round(row.tc, digits=5), 12),
                 rpad(round(tc_pred, digits=5), 12),
                 rpad(round(cor_err, digits=3), 9), round(tc_err, digits=3))
    end

    println()
    println("median relative error, CoR: ", round(median(cor_errs), digits=4),
            "  (", length(cor_errs), "/", n_samples, " solved)")
    println("median relative error, tc:  ", round(median(tc_errs), digits=4),
            "  (", length(tc_errs), "/", n_samples, " solved)")
end

n_samples = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 20
seed      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
M_arg     = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 12
main(n_samples, seed; M=M_arg)
