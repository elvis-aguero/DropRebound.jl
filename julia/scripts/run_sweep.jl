#!/usr/bin/env julia
# Parameter sweep over (Oh, Bo, We, De1, beta_s, M).
# Streams results to a CSV file; skips rows already present (resume-safe).
#
# Usage:
#   julia --project=.. scripts/run_sweep.jl [output.csv]
#
# Customise the parameter grids below, then run.  Each row in the CSV is one
# (Oh, Bo, We, De1, beta_s, M) combination together with the KPIs extracted
# by extract_kpis: contact_time, cor, max_radius, max_A2.

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf
using Logging

# ---------------------------------------------------------------------------
# Parameter grids — edit these
# ---------------------------------------------------------------------------
const Oh_grid    = [0.05, 0.10, 0.20, 0.30]
const Bo_grid    = [1/53.9]           # ≈ 0.0186 (water-glycerol, R=0.2 mm)
const We_grid    = [0.05, 0.10, 0.20] # impact KE: We = v₀²
const De1_grid   = [0.0, 0.5]
const beta_s_grid = [1.0, 0.5]        # 1.0 = Newtonian
const M_grid     = [20]

# ---------------------------------------------------------------------------
# Output file
# ---------------------------------------------------------------------------
const CSV_PATH = length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "..", "outputs", "csv", "sweep.csv")
const CSV_HEADER = "Oh,Bo,We,De1,beta_s,M,contact_time,cor,max_radius,max_A2\n"

# ---------------------------------------------------------------------------
# Build the set of already-completed rows so we can resume
# ---------------------------------------------------------------------------
function load_completed(path)
    completed = Set{NTuple{6,Float64}}()
    isfile(path) || return completed
    for line in eachline(path)
        startswith(line, "Oh") && continue   # header
        parts = split(line, ',')
        length(parts) < 6 && continue
        try
            key = ntuple(i -> parse(Float64, parts[i]), 6)
            push!(completed, key)
        catch
        end
    end
    return completed
end

# ---------------------------------------------------------------------------
# Run one simulation and return a CSV line
# ---------------------------------------------------------------------------
function run_one(Oh, Bo, We, De1, beta_s, M)
    v0 = -sqrt(We)   # falling, so negative
    M_int = Int(M)

    dt_max    = make_dt_max(M_int)
    theta_vec = make_theta_vec(M_int)
    precomp   = precompute_integrals(NaN, M_int)[1]
    cfg       = SimConstants(M_int, M_int+1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams(De1, beta_s)

    init = DropState(M_int)
    init.z  = 1.1
    init.v  = v0
    init.dt = dt_max
    init.cp = 0

    t_end = max(10.0, 5.0 / (Oh + 0.01))   # run longer for low-Oh drops

    times, states = with_logger(NullLogger()) do
        solve_drop!(cfg, ob, deepcopy(init); t_end=t_end, save_every=0.05)
    end

    kpis = extract_kpis(times, states, cfg)
    return kpis
end

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
mkpath(dirname(CSV_PATH))

# Open in append mode; write header only if file is new
file_exists = isfile(CSV_PATH)
completed   = load_completed(CSV_PATH)

io = open(CSV_PATH, "a")
file_exists || write(io, CSV_HEADER)

total   = length(Oh_grid) * length(Bo_grid) * length(We_grid) *
          length(De1_grid) * length(beta_s_grid) * length(M_grid)
done    = 0
skipped = 0

@printf("Sweep: %d combinations → %s\n", total, CSV_PATH)
println("-"^60)

for Oh in Oh_grid, Bo in Bo_grid, We in We_grid,
    De1 in De1_grid, beta_s in beta_s_grid, M in M_grid

    # skip Newtonian-equivalent duplicates: De1>0 + beta_s=1 is just Newtonian
    De1 > 0 && beta_s == 1.0 && continue

    key = (Oh, Bo, We, De1, beta_s, Float64(M))
    if key in completed
        skipped += 1
        continue
    end

    done += 1
    @printf("[%4d/%d] Oh=%.3f  Bo=%.4f  We=%.3f  De1=%.2f  β_s=%.2f  M=%d  ",
            done, total - skipped, Oh, Bo, We, De1, beta_s, M)
    flush(stdout)

    t_start = time()
    try
        kpis = run_one(Oh, Bo, We, De1, beta_s, M)
        elapsed = time() - t_start
        @printf("→  τ_c=%.3f  COR=%.3f  r_max=%.3f  [%.1fs]\n",
                kpis.contact_time, kpis.cor, kpis.max_radius, elapsed)
        @printf(io, "%.6g,%.6g,%.6g,%.6g,%.6g,%d,%.6g,%.6g,%.6g,%.6g\n",
                Oh, Bo, We, De1, beta_s, M,
                kpis.contact_time, kpis.cor, kpis.max_radius, kpis.max_A2)
        flush(io)
    catch e
        elapsed = time() - t_start
        @printf("→  ERROR: %s  [%.1fs]\n", e, elapsed)
        @printf(io, "%.6g,%.6g,%.6g,%.6g,%.6g,%d,NaN,NaN,NaN,NaN\n",
                Oh, Bo, We, De1, beta_s, M)
        flush(io)
    end
end

close(io)
@printf("\nDone. %d computed, %d skipped (already in CSV).\n", done, skipped)
