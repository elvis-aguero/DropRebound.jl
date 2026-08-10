# Draws the call map produced by `callmap.jl`.
#
# Layout is the reading order: level 0 at the top depends on nothing else in the
# package, and each row below it uses only what is above. Every box lists the names that
# file defines, so the granularity is function-level without the graph becoming a
# hairball of 136 nodes.
#
# Deliberately plain: no colour scale, no legend beyond two swatches, and type large
# enough to read at a glance. The picture has one job, which is to tell a reader where to
# start and what is unreachable.

using Printf
using Plots
gr()

const ROOT = joinpath(@__DIR__, "..", "..")
const RES  = joinpath(ROOT, "results")

nodes = Dict{String,NamedTuple}()
for (i, ln) in enumerate(eachline(joinpath(RES, "callmap_nodes.csv")))
    i == 1 && continue
    p = split(ln, ','; limit = 4)
    nodes[p[1]] = (level = parse(Int, p[2]), n = parse(Int, p[3]),
                   defs = isempty(p[4]) ? String[] : split(p[4], ' '))
end
edges = Tuple{String,String}[]
for (i, ln) in enumerate(eachline(joinpath(RES, "callmap.csv")))
    i == 1 && continue
    p = split(strip(ln), ',')
    length(p) == 2 && push!(edges, (p[1], p[2]))
end

## `DropSolver.jl` includes and re-exports everything, so drawing its edges would
## connect it to all twenty boxes and say nothing. It is shown, without them.
const HUB = "DropSolver.jl"
edges = filter(e -> e[1] != HUB, edges)

## Files nothing else uses: entry points, or dead. Marked, not hidden.
others  = filter(!=(HUB), collect(keys(nodes)))
unused  = Set(f for f in others if !any(e -> e[2] == f && e[1] != HUB, edges))

# ---------------------------------------------------------------------------
# Geometry. One unit is one line of text.
# ---------------------------------------------------------------------------
const BW    = 40.0     # box width
const GAPX  = 6.0      # between boxes in a row
const GAPY  = 7.0      # between rows
const HEAD  = 3.4      # title lines inside a box
const LINE  = 1.45     # vertical pitch of one definition name

levels = sort(unique(v.level for v in values(nodes)))
rows   = Dict(l => sort([f for (f, v) in nodes if v.level == l]) for l in levels)
boxh(f) = HEAD + max(nodes[f].n, 1) * LINE + 1.2

## Row y positions, level 0 at the top. In a function: a top-level `for` opens soft
## scope, and an accumulator assigned inside it becomes a fresh local.
function rowlayout(levels, rows, boxh)
    rowtop = Dict{Int,Float64}(); y = 0.0
    for l in levels
        rowtop[l] = y
        y -= (maximum(boxh(f) for f in rows[l]) + GAPY)
    end
    (rowtop, y)
end
rowtop, ybot = rowlayout(levels, rows, boxh)
## box x positions, rows centred on each other
maxw = maximum(length(rows[l]) * (BW + GAPX) for l in levels)
pos = Dict{String,NTuple{2,Float64}}()
for l in levels
    w  = length(rows[l]) * (BW + GAPX) - GAPX
    x0 = (maxw - GAPX - w) / 2
    for (i, f) in enumerate(rows[l])
        pos[f] = (x0 + (i - 1) * (BW + GAPX), rowtop[l])
    end
end

const INK   = RGB(0.10, 0.12, 0.16)
const FILL  = RGB(0.945, 0.955, 0.97)
const EDGE  = RGB(0.62, 0.67, 0.74)
const DEAD  = RGB(0.99, 0.93, 0.86)
const DEADL = RGB(0.85, 0.45, 0.15)
const ARROW = RGB(0.68, 0.72, 0.78)

plt = plot(; size = (9000, 9000), dpi = 100, legend = false, grid = false,
           framestyle = :none, background_color = :white,
           xlims = (-4, maxw + 2), ylims = (ybot - 5, 16))

## edges first, so boxes sit on top of them
for (a, b) in edges
    haskey(pos, a) && haskey(pos, b) || continue
    xa, ya = pos[a]; xb, yb = pos[b]
    x1, y1 = xa + BW / 2, ya - boxh(a)          # bottom of the user
    x2, y2 = xb + BW / 2, yb                    # top of the dependency
    ## a slight bow so parallel edges do not overprint
    xm = (x1 + x2) / 2 + 0.06 * (y1 - y2)
    plot!(plt, [x1, xm, x2], [y1, (y1 + y2) / 2, y2];
          lc = ARROW, lw = 4.0, alpha = 0.85, label = "")
end

for (f, v) in nodes
    x, yt = pos[f]; h = boxh(f)
    isdead = f in unused
    plot!(plt, Shape([x, x + BW, x + BW, x], [yt, yt, yt - h, yt - h]);
          fc = isdead ? DEAD : FILL, lc = isdead ? DEADL : EDGE,
          lw = isdead ? 5 : 3, label = "")
    annotate!(plt, x + 1.2, yt - 1.7, text(f, 50, INK, :left, "Courier"))
    annotate!(plt, x + BW - 1.2, yt - 1.7,
              text("L$(v.level)", 40, EDGE, :right, "Courier"))
    for (i, d) in enumerate(v.defs)
        annotate!(plt, x + 1.8, yt - HEAD - (i - 1) * LINE - 0.5,
                  text(d, 34, RGB(0.24, 0.28, 0.34), :left, "Courier"))
    end
    isempty(v.defs) && annotate!(plt, x + 1.8, yt - HEAD - 0.5,
                                 text("(module: includes and exports)", 32, EDGE, :left))
end

annotate!(plt, 0, 12.0, text("DropRebound.jl call map", 104, INK, :left))
annotate!(plt, 0, 6.6,
          text("read top to bottom: level 0 depends on nothing else in the package; " *
               "each row uses only the rows above it", 46, RGB(0.42,0.47,0.54), :left))
annotate!(plt, 0, 2.6,
          text("arrows point from a file to a file whose names it uses   |   " *
               "orange outline: used by nothing else in the package", 46, DEADL, :left))

out = joinpath(ROOT, "docs", "figures", "callmap.png")
mkpath(dirname(out)); savefig(plt, out)
@printf("wrote %s  (%d boxes, %d edges, %d marked unused)\n",
        out, length(nodes), length(edges), length(unused))
