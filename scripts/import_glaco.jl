# Convert `data/Glaco_summary.xlsx` to plain CSV.
#
# The repository keeps experimental data as CSV, not as spreadsheets: a CSV diffs,
# a workbook does not, and the numbers here are read by figure scripts that should
# not depend on a spreadsheet reader. This is the same route the Gabbard data took.
#
# An .xlsx is a zip of XML, so the conversion needs no package: `unzip` puts the
# sheet on disk and the cell values are read off with a regular expression. Run it
# once after replacing the workbook; the CSV it writes is what everything else uses.
#
#   julia --project=. scripts/import_glaco.jl
#
# The sheet carries five columns. Two are raw measurements (impact speed, drop
# diameter), and the other three are quantities Excel derived from them; those are
# recomputed here rather than trusted, and the script stops if they disagree.

using Printf

const ROOT = joinpath(@__DIR__, "..")
const XLSX = joinpath(ROOT, "data", "Glaco_summary.xlsx")
const OUT  = joinpath(ROOT, "data", "glaco_restitution.csv")

## Water, at the values the workbook's own U_cap formula uses: SQRT(72/(997*D/2)).
const SIGMA = 0.072     # N/m
const RHO   = 997.0     # kg/m^3
const MU    = 1.0e-3    # Pa s
const GRAV  = 9.81      # m/s^2

tmp = mktempdir()
run(pipeline(`unzip -o -q $XLSX -d $tmp`; stdout = devnull))
sheet = read(joinpath(tmp, "xl", "worksheets", "sheet1.xml"), String)

## Columns, from row 1 of the workbook: A U_impact (cm/s), B D (mm),
## C U_cap (cm/s), D U_norm, E epsilon. Blank cells simply do not appear.
rows = Dict{Int,Dict{Char,Float64}}()
for m in eachmatch(r"<c r=\"([A-E])(\d+)\"[^>]*>(?:<f>.*?</f>)?<v>([^<]*)</v></c>", sheet)
    col, row, val = m.captures[1][1], parse(Int, m.captures[2]), tryparse(Float64, m.captures[3])
    row == 1 && continue                       # the header, which is a string index
    val === nothing && continue
    get!(rows, row, Dict{Char,Float64}())[col] = val
end

## Wrapped in a function: a top-level loop is soft scope, so `dropped += 1` there
## would silently become a fresh local and throw on read.
function convert_rows(rows)
    kept, dropped = NTuple{6,Float64}[], 0
    for r in sort(collect(keys(rows)))
        c = rows[r]
        all(haskey(c, k) for k in ('A', 'B', 'E')) || (dropped += 1; continue)
        U, D_m, eps = c['A'] / 100, c['B'] / 1000, c['E']    # SI: m/s, m, dimensionless
        R = D_m / 2
        (isfinite(U) && isfinite(R) && isfinite(eps) && R > 0 && U > 0) ||
            (dropped += 1; continue)

        ## Recomputed, not read: the workbook's U_cap and U_norm are checks, not inputs.
        U_cap = sqrt(SIGMA / (RHO * R))
        haskey(c, 'C') && @assert isapprox(U_cap * 100, c['C']; rtol = 1e-6) "U_cap disagrees at row $r"
        haskey(c, 'D') && @assert isapprox(U / U_cap, c['D']; rtol = 1e-6) "U_norm disagrees at row $r"

        We = RHO * U^2 * R / SIGMA
        Oh = MU / sqrt(RHO * SIGMA * R)
        Bo = RHO * GRAV * R^2 / SIGMA
        push!(kept, (We, eps, Oh, Bo, R, U))
    end
    sort!(kept, by = first)
    kept, dropped
end

kept, dropped = convert_rows(rows)
open(OUT, "w") do io
    println(io, "# Drop rebound on a Glaco superhydrophobic surface, water.")
    println(io, "# Converted from data/Glaco_summary.xlsx by scripts/import_glaco.jl.")
    println(io, "#")
    println(io, "# The workbook carries impact speed and drop diameter; We, Oh and Bo are")
    println(io, "# computed here from sigma = $SIGMA N/m, rho = $RHO kg/m^3, mu = $MU Pa s,")
    println(io, "# the same water properties its own U_cap column assumes.")
    println(io, "We,cor,Oh,Bo,R_m,U_impact_m_s")
    for (We, eps, Oh, Bo, R, U) in kept
        @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", We, eps, Oh, Bo, R, U)
    end
end

@printf("wrote %s\n  %d rows kept, %d skipped\n", OUT, length(kept), dropped)
@printf("  We  in [%.4g, %.4g]\n", minimum(x -> x[1], kept), maximum(x -> x[1], kept))
@printf("  cor in [%.4g, %.4g]\n", minimum(x -> x[2], kept), maximum(x -> x[2], kept))
@printf("  Oh  in [%.4g, %.4g]\n", minimum(x -> x[3], kept), maximum(x -> x[3], kept))
@printf("  Bo  in [%.4g, %.4g]\n", minimum(x -> x[4], kept), maximum(x -> x[4], kept))
@printf("  R   in [%.4g, %.4g] m\n", minimum(x -> x[5], kept), maximum(x -> x[5], kept))
