# Convert the Thenarianto et al. (2023) workbooks to plain CSV.
#
#   data/Glaco_summary.xlsx  ->  data/glaco_restitution.csv   (Glaco coating)
#   data/BSi_summary.xlsx    ->  data/bsi_restitution.csv     (black silicon)
#
# Both are water on a superhydrophobic surface, and both carry the same five
# columns, so one converter serves them. The surface is the only thing that
# differs, which is exactly what makes the pair worth having: same liquid, same
# drop sizes, different substrate.
#
# The repository keeps experimental data as CSV, not as spreadsheets: a CSV diffs,
# a workbook does not, and the numbers here are read by figure scripts that should
# not depend on a spreadsheet reader. This is the same route the Gabbard data took.
#
# An .xlsx is a zip of XML, so the conversion needs no package: `unzip` puts the
# sheet on disk and the cell values are read off with a regular expression. Run it
# once after replacing a workbook; the CSVs it writes are what everything else uses.
#
#   julia --project=. scripts/import_thenarianto.jl
#
# Each sheet carries five columns. Two are raw measurements (impact speed, drop
# diameter), and the other three are quantities Excel derived from them; those are
# recomputed here rather than trusted, and the script stops if they disagree.

using Printf

const ROOT = joinpath(@__DIR__, "..")

## Water, at the values the workbooks' own U_cap formula uses: SQRT(72/(997*D/2)).
const SIGMA = 0.072     # N/m
const RHO   = 997.0     # kg/m^3
const MU    = 1.0e-3    # Pa s
const GRAV  = 9.81      # m/s^2

const SOURCES = [
    (book = "Glaco_summary.xlsx", csv = "glaco_restitution.csv", surface = "a Glaco coating"),
    (book = "BSi_summary.xlsx",   csv = "bsi_restitution.csv",   surface = "black silicon (BSi)"),
]

"""Cell values of the first worksheet, as `row => Dict(column => value)`."""
function read_sheet(xlsx)
    tmp = mktempdir()
    run(pipeline(`unzip -o -q $xlsx -d $tmp`; stdout = devnull))
    sheet = read(joinpath(tmp, "xl", "worksheets", "sheet1.xml"), String)
    rows = Dict{Int,Dict{Char,Float64}}()
    for m in eachmatch(r"<c r=\"([A-E])(\d+)\"[^>]*>(?:<f>.*?</f>)?<v>([^<]*)</v></c>", sheet)
        col, row = m.captures[1][1], parse(Int, m.captures[2])
        val = tryparse(Float64, m.captures[3])
        (row == 1 || val === nothing) && continue     # row 1 is the header, as string indices
        get!(rows, row, Dict{Char,Float64}())[col] = val
    end
    rows
end

## THE "D (mm)" COLUMN IS A DIAMETER THROUGHOUT, AND PART OF ONE SHEET DISAGREES.
##
## In BSi_summary the U_cap column is a live formula down to row 134 and hardcoded
## from row 135 on, and the hardcoded values are the formula evaluated with the full
## D rather than D/2 -- smaller by exactly sqrt(2). Since U_norm is U/U_cap and the
## Weber number is U_norm squared, those 36 rows carry a Weber number twice the
## correct one, and a radius twice the correct one if U_cap is inverted for it.
##
## The radius is therefore taken from the D column, halved, everywhere. That is what
## makes the drop sizes agree with the experiments as performed: the last block
## becomes R = 1.33 mm, inside the 1.3-1.5 mm the authors report, where inverting
## their U_cap would give 2.65 mm -- a drop near the capillary length, which is not
## what was run.
##
## Their U_cap is still checked, and a mismatch is reported rather than adopted. An
## earlier version of this script inverted U_cap to get the radius, which sounds
## conservative and is not: it reproduces the sheet exactly, including its error.
##
## Wrapped in a function: a top-level loop is soft scope, so `dropped += 1` there
## would silently become a fresh local and throw on read.
function convert_rows(rows)
    kept, dropped, mismatched = NTuple{6,Float64}[], 0, Int[]
    for r in sort(collect(keys(rows)))
        c = rows[r]
        ## Columns, from row 1 of both workbooks: A U_impact (cm/s), B D (mm),
        ## C U_cap (cm/s), D U_norm, E epsilon. Blank cells simply do not appear.
        all(haskey(c, k) for k in ('A', 'B', 'E')) || (dropped += 1; continue)
        U, D_m, eps = c['A'] / 100, c['B'] / 1000, c['E']    # SI: m/s, m, dimensionless
        R = D_m / 2
        (isfinite(U) && isfinite(R) && isfinite(eps) && R > 0 && U > 0) ||
            (dropped += 1; continue)

        ## Recomputed, not read. Where the sheet disagrees it is the sheet that is
        ## wrong, so the row is kept with the recomputed value and the row noted.
        U_cap = sqrt(SIGMA / (RHO * R))
        haskey(c, 'C') && !isapprox(U_cap * 100, c['C']; rtol = 1e-5) && push!(mismatched, r)

        We = RHO * U^2 * R / SIGMA
        Oh = MU / sqrt(RHO * SIGMA * R)
        Bo = RHO * GRAV * R^2 / SIGMA
        push!(kept, (We, eps, Oh, Bo, R, U))
    end
    sort!(kept, by = first)
    kept, dropped, mismatched
end

function convert(src)
    xlsx = joinpath(ROOT, "data", src.book)
    isfile(xlsx) || (@warn "workbook not found, skipping" book = src.book; return)
    kept, dropped, mismatched = convert_rows(read_sheet(xlsx))
    out = joinpath(ROOT, "data", src.csv)

    open(out, "w") do io
        println(io, "# Drop rebound on $(src.surface), water.")
        println(io, "# Thenarianto et al. (2023).")
        println(io, "# Converted from data/$(src.book) by scripts/import_thenarianto.jl.")
        println(io, "#")
        println(io, "# The workbook carries impact speed and drop diameter; We, Oh and Bo are")
        println(io, "# computed here from sigma = $SIGMA N/m, rho = $RHO kg/m^3, mu = $MU Pa s,")
        println(io, "# the same water properties its own U_cap column assumes. The radius is")
        println(io, "# taken from that U_cap column rather than from the \"D (mm)\" column, which")
        println(io, "# changes meaning partway through BSi_summary; see the script for why.")
        println(io, "We,cor,Oh,Bo,R_m,U_impact_m_s")
        for (We, eps, Oh, Bo, R, U) in kept
            @printf(io, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n", We, eps, Oh, Bo, R, U)
        end
    end

    @printf("wrote %s\n  %d rows kept, %d skipped\n", out, length(kept), dropped)
    ## Not fatal, but not silent: the recomputed value is used, and the rows whose
    ## stored U_cap disagrees are named so the discrepancy can be taken up with the
    ## authors rather than inherited.
    if !isempty(mismatched)
        @warn """$(src.book): $(length(mismatched)) rows carry a stored U_cap that \
                 disagrees with sqrt(sigma/(rho R)) at R = D/2 (rows \
                 $(first(mismatched))-$(last(mismatched))). Those values are the \
                 formula evaluated with D instead of D/2, low by sqrt(2), which \
                 would put their Weber number a factor of two high. Recomputed here."""
    end
    @printf("  We  in [%.4g, %.4g]   cor in [%.4g, %.4g]\n",
            minimum(x -> x[1], kept), maximum(x -> x[1], kept),
            minimum(x -> x[2], kept), maximum(x -> x[2], kept))
    @printf("  Oh  in [%.4g, %.4g]   Bo  in [%.4g, %.4g]   R in [%.3g, %.3g] mm\n\n",
            minimum(x -> x[3], kept), maximum(x -> x[3], kept),
            minimum(x -> x[4], kept), maximum(x -> x[4], kept),
            1000minimum(x -> x[5], kept), 1000maximum(x -> x[5], kept))
end

foreach(convert, SOURCES)
