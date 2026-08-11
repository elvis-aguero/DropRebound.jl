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

const ROOT = joinpath(@__DIR__, "..")
const RES  = joinpath(ROOT, "outputs", "csv")

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
           xlims = (-4, maxw + 2), ylims = (ybot - 5, 26))

## HUBS: the files most depended upon. Those are the ones to understand first, so
## their incoming edges are coloured and everything else stays grey. Colouring all
## sixteen targets would be a rainbow; colouring the six that matter is a map.
indeg = Dict{String,Int}()
for (_, b) in edges; indeg[b] = get(indeg, b, 0) + 1; end
hubs = first(sort(collect(keys(indeg)); by = k -> -indeg[k]), 6)
const PALETTE = [RGB(0.09,0.50,0.52), RGB(0.83,0.52,0.09), RGB(0.28,0.33,0.71),
                 RGB(0.76,0.26,0.41), RGB(0.28,0.56,0.27), RGB(0.51,0.31,0.64)]
hubcol = Dict(h => PALETTE[i] for (i, h) in enumerate(hubs))
ecol(b) = get(hubcol, b, ARROW)

## Edges run from the TOP of the user up to the BOTTOM of the dependency, which sits on
## an earlier row. Anchoring them to the far sides of each box, as a previous version
## did, drew every line straight through the boxes it connected.
for (a, b) in edges
    haskey(pos, a) && haskey(pos, b) || continue
    xa, ya = pos[a]; xb, yb = pos[b]
    x1, y1 = xa + BW / 2, ya                    # top of the user
    x2, y2 = xb + BW / 2, yb - boxh(b)          # bottom of the dependency
    c = ecol(b); ishub = haskey(hubcol, b)
    ## a slight bow so parallel edges do not overprint
    xm = (x1 + x2) / 2 + 0.05 * (y2 - y1)
    plot!(plt, [x1, xm, x2], [y1, (y1 + y2) / 2, y2];
          lc = c, lw = ishub ? 5.0 : 3.0, alpha = ishub ? 0.85 : 0.5, label = "")
    ## arrowhead at the dependency end, pointing at it
    dx, dy = x2 - xm, y2 - (y1 + y2) / 2
    n = hypot(dx, dy); n < 1e-9 && continue
    ux, uy = dx / n, dy / n; px, py = -uy, ux
    hl, hw = 1.6, 0.62
    plot!(plt, Shape([x2, x2 - hl*ux + hw*px, x2 - hl*ux - hw*px],
                     [y2, y2 - hl*uy + hw*py, y2 - hl*uy - hw*py]);
          fc = c, lc = c, alpha = ishub ? 0.9 : 0.55, label = "")
end

for (f, v) in nodes
    x, yt = pos[f]; h = boxh(f)
    isdead = f in unused
    bc = isdead ? DEADL : get(hubcol, f, EDGE)
    plot!(plt, Shape([x, x + BW, x + BW, x], [yt, yt, yt - h, yt - h]);
          fc = isdead ? DEAD : FILL, lc = bc,
          lw = (isdead || haskey(hubcol, f)) ? 6 : 3, label = "")
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

annotate!(plt, 0, 21.5, text("DropRebound.jl call map", 104, INK, :left))
annotate!(plt, 0, 15.5,
          text("read top to bottom: level 0 depends on nothing else in the package; " *
               "each row uses only the rows above it", 46, RGB(0.42,0.47,0.54), :left))
annotate!(plt, 0, 11.0,
          text("arrows point from a file to a file whose names it uses   |   " *
               "orange outline: used by nothing else in the package", 46, DEADL, :left))

## the hub key, drawn in the hubs' own colours
let x = 0.0
    annotate!(plt, x, 6.0, text("most depended upon:", 40, RGB(0.42,0.47,0.54), :left))
    x += 27
    for h in hubs
        annotate!(plt, x, 6.0, text(h, 40, hubcol[h], :left, "Courier"))
        x += 1.05 * length(h) + 5
    end
end

out = joinpath(ROOT, "outputs", "figures", "callmap_draw.png")
mkpath(dirname(out)); savefig(plt, out)
@printf("wrote %s  (%d boxes, %d edges, %d marked unused)\n",
        out, length(nodes), length(edges), length(unused))
