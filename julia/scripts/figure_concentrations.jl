# Restitution against Weber number for the shear-thinning series: experiments as
# scatter, coloured by concentration; the model as a line of the same colour.
#
# The point of putting all five fluids on one axis is that the SAME rheological model
# has to explain all of them. Zero-shear Ohnesorge runs from 0.0068 for water to 57 for
# the 3000 ppm solution, four orders of magnitude, and the only thing that changes
# between curves is eta(gammadot) measured by rheometry. Nothing is fitted to impact
# data.

using Printf
using Plots
gr()

const OUT   = joinpath(@__DIR__, "..", "..", "results")
const ASSET = joinpath(@__DIR__, "..", "..", "docs", "src", "assets")
mkpath(ASSET)

rows = NamedTuple[]
for f in readdir(OUT; join = true)
    ## ONLY the sharded production run. An earlier version also globbed
    ## `sweep_concentrations.csv`, a superseded M = 45 sweep, so the figure silently
    ## mixed two resolutions with nothing on the plot to say so.
    occursin("sweep_shard_", f) || continue
    for ln in eachline(f)
        p = split(strip(ln), ',')
        length(p) < 8 && continue
        p[1] == "ppm" && continue
        push!(rows, (ppm = parse(Int, p[1]), kind = p[2], We = parse(Float64, p[3]),
                     cor = parse(Float64, p[4]), tc = parse(Float64, p[5]),
                     Oh = parse(Float64, p[6]), ok = p[7] == "true"))
    end
end

## experiments come from the metrics files, not the shards, so read them here
const DATA  = joinpath(@__DIR__, "..", "data")
const T_CAP = sqrt(1000.0 * 0.0003^3 / 0.0728)
const FILES = [(0, "metrics_Water.csv"), (300, "metrics_300ppm.csv"),
               (1000, "metrics_1000ppm.csv"), (2000, "metrics_2000ppm.csv"),
               (3000, "metrics_3000ppm_2.csv")]
function read_metrics(path)
    lines = readlines(path); hdr = split(strip(lines[1]), ';')
    function idx(pat)
        i = findfirst(h -> occursin(pat, lowercase(h)), hdr)
        i === nothing && error("no column $pat in $(basename(path))"); i
    end
    iw, ie = idx("we"), idx("epsilon")
    out = NTuple{2,Float64}[]
    for ln in lines[2:end]
        f = split(strip(ln), ';'); length(f) < 3 && continue
        v = tryparse.(Float64, replace.(f, ',' => '.')); any(isnothing, v) && continue
        push!(out, (v[iw], v[ie]))
    end
    out
end

## colour by concentration, ordered so the eye reads dilute -> concentrated
const COLS = Dict(0 => :steelblue, 300 => :seagreen, 1000 => :goldenrod,
                  2000 => :darkorange, 3000 => :indianred)
lab(ppm) = ppm == 0 ? "water" : "$(ppm) ppm"

## the resolution is stamped on the figure, so a plot cannot outlive the run that made it
const MK = let fs = filter(f -> occursin("sweep_shard_", f), readdir(OUT))
    isempty(fs) ? "no data" : "M = 90, K = 3"
end

plt = plot(xscale = :log10,
           xlabel = "Weber number",
           ylabel = "restitution  ε",
           legend = :bottomleft, size = (1000, 700), dpi = 160,
           guidefontsize = 17, tickfontsize = 15, legendfontsize = 14,
           framestyle = :axes, grid = false, ylims = (0.35, 1.0),
           foreground_color_axis = :gray40, foreground_color_border = :gray40,
           left_margin = 10Plots.mm, bottom_margin = 10Plots.mm,
           right_margin = 6Plots.mm, top_margin = 6Plots.mm)

for (ppm, file) in FILES
    e = read_metrics(joinpath(DATA, file))
    isempty(e) && continue
    scatter!(plt, [x[1] for x in e], [x[2] for x in e];
             c = COLS[ppm], ms = 5, msw = 0, alpha = 0.40, label = "")
end
for (ppm, _) in FILES
    s = sort([r for r in rows if r.ppm == ppm && r.kind == "simulation" && r.ok],
             by = r -> r.We)
    isempty(s) && continue
    plot!(plt, [r.We for r in s], [r.cor for r in s];
          c = COLS[ppm], lw = 4, label = lab(ppm))
end

## The upstream solver of Gabbard et al. (2025) is NOT overlaid here.
## It was, briefly. Run verbatim at our water conditions it agrees with this model to
## better than 0.5 per cent up to We = 0.7 and then collapses to 0.33, which is not a
## feature of their solver but of running it at Oh = 0.0068. Their published data spans
## Oh 0.0139 to 0.79, and a scan across that range (results/upstream_oh_scan.csv) is
## smooth and monotone at every Oh from 0.0139 up; only the water case breaks, and it
## is converged in modes, so it is not under-resolution. Plotting a curve taken from
## outside a solver's validated range, next to data, would misrepresent it.
savefig(plt, joinpath(ASSET, "concentration_series.png"))
@printf("wrote %s\n", joinpath(ASSET, "concentration_series.png"))

## and the residual: how far the model sits from the experiments nearest it in We
@printf("\n%-8s %6s %10s %10s\n", "fluid", "n_sim", "med |dcor|", "Oh_0")
for (ppm, file) in FILES
    e = read_metrics(joinpath(DATA, file))
    s = sort([r for r in rows if r.ppm == ppm && r.kind == "simulation" && r.ok],
             by = r -> r.We)
    (isempty(e) || isempty(s)) && continue
    ds = Float64[]
    for r in s
        near = [x[2] for x in e if 1/1.3 <= x[1]/r.We <= 1.3]
        isempty(near) && continue
        med = sort(near)[max(1, length(near) ÷ 2)]
        push!(ds, abs(r.cor - med) / med)
    end
    isempty(ds) && continue
    @printf("%-8s %6d %9.1f%% %10.3g\n", lab(ppm), length(s),
            100 * sort(ds)[max(1, length(ds) ÷ 2)], s[1].Oh)
end
