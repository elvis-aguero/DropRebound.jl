#!/usr/bin/env julia
# Model against experiment, read from the provenance store rather than resimulated.
#
# Contact and restitution are the PROXIMITY definitions: contact exists whenever any point
# of the surface is below 0.02R, contact time runs from the first such instant to the last,
# and restitution compares the centre-of-mass speed at those two instants. That is what the
# measurement does, and unlike a criterion based on the solver's contact set it does not
# fragment when the drop separates by a thousandth of a radius.
#
# Every comparison point needs at least five experimental replicates, so the target is a mean
# with a spread. The verdict is whether the model lies INSIDE that spread. Two models whose
# median errors differ by a point are not ranked by that difference if both sit inside one
# standard deviation -- and if both sit outside it, the difference between them is not the
# interesting quantity.

using Printf
include(joinpath(@__DIR__, "_stats.jl"))

const ROOT = joinpath(@__DIR__, "..", "..")
const STORE = joinpath(ROOT, "results", "runs.csv")
const MIN_REPS = 5

function read_csv(path; sep = ',', dec = '.')
    hdr = String[]; rows = Vector{String}[]
    for ln in eachline(path)
        startswith(ln, "#") && continue
        f = String.(split(chomp(ln), sep))
        if isempty(hdr); hdr = f; continue; end
        length(f) == length(hdr) || continue
        push!(rows, dec == '.' ? f : replace.(f, dec => '.'))
    end
    (hdr, rows)
end

col(hdr, rows, name) = begin
    j = findfirst(==(name), hdr)
    [tryparse(Float64, r[j]) for r in rows]
end
scol(hdr, rows, name) = begin
    j = findfirst(==(name), hdr)
    [r[j] for r in rows]
end

# ---- the store -------------------------------------------------------------------
shdr, srows = read_csv(STORE)
S = (solver = scol(shdr, srows, "solver"), rheo = scol(shdr, srows, "rheology"),
     We = col(shdr, srows, "We"), Oh = col(shdr, srows, "Oh"),
     cor = col(shdr, srows, "cor"), tc = col(shdr, srows, "tc"),
     M = col(shdr, srows, "M"), K = col(shdr, srows, "K"))
@printf("store: %d runs\n", length(S.We))

"""Log-linear interpolation of a stored model curve at one Weber number."""
function model_at(solver, rheo, Oh, We, field)
    sel = [i for i in eachindex(S.We) if S.solver[i] == solver && S.rheo[i] == rheo &&
           isapprox(S.Oh[i], Oh; rtol = 1e-6) && isfinite(getproperty(S, field)[i])]
    length(sel) < 2 && return NaN
    ws = S.We[sel]; vs = getproperty(S, field)[sel]
    o = sortperm(ws); ws = ws[o]; vs = vs[o]
    (We < ws[1] || We > ws[end]) && return NaN        # never extrapolate
    k = searchsortedlast(ws, We); k = clamp(k, 1, length(ws)-1)
    t = (log(We) - log(ws[k])) / (log(ws[k+1]) - log(ws[k]))
    vs[k] + t*(vs[k+1] - vs[k])
end

"""Groups of at least MIN_REPS experiments sharing a Weber bin."""
function wegroups(we, y, sel, ng)
    edges = exp.(range(log(quantile(we[sel],0.02)), log(quantile(we[sel],0.98)); length=ng+1))
    out = []
    for k in 1:ng
        g = [i for i in sel if edges[k] <= we[i] < edges[k+1]]
        length(g) >= MIN_REPS || continue
        push!(out, (We = median(we[g]), m = mean(y[g]), sd = std(y[g]), n = length(g)))
    end
    out
end

function report(label, groups_cor, groups_tc, solvers, rheo, Oh)
    @printf("\n%s\n", label)
    @printf("  %-9s %-4s | %-18s %-10s %-10s | %-18s %-10s %s\n",
            "We","n","CoR exp (mean+-sd)","search","LCP","tc exp (mean+-sd)","search","LCP")
    acc = Dict(s => (c=Float64[], t=Float64[], ic=0, it=0, n=0) for s in solvers)
    stats = Dict(s => Dict(:c=>Float64[], :t=>Float64[], :ic=>[0], :it=>[0]) for s in solvers)
    for g in groups_cor
        gt = findfirst(x -> isapprox(x.We, g.We; rtol=0.35), groups_tc)
        mk(v,m,s) = !isfinite(v) ? "   --    " : (abs(v-m) <= s ? @sprintf("%.4f*",v) : @sprintf("%.4f ",v))
        vals = Dict(s => (model_at(s, rheo, Oh, g.We, :cor),
                          model_at(s, rheo, Oh, g.We, :tc)) for s in solvers)
        tgt = gt === nothing ? nothing : groups_tc[gt]
        @printf("  %-9.4g %-4d | %.4f+-%.4f     %-10s %-10s | %s %-10s %s\n",
                g.We, g.n, g.m, g.sd,
                mk(vals["search"][1],g.m,g.sd), mk(vals["lcp"][1],g.m,g.sd),
                tgt === nothing ? "    too few      " : @sprintf("%.3f+-%.3f (n=%2d)", tgt.m, tgt.sd, tgt.n),
                tgt === nothing ? "   --    " : mk(vals["search"][2],tgt.m,tgt.sd),
                tgt === nothing ? "   --    " : mk(vals["lcp"][2],tgt.m,tgt.sd))
        for s in solvers
            c, t = vals[s]
            if isfinite(c)
                push!(stats[s][:c], abs(c-g.m)/g.m)
                abs(c-g.m) <= g.sd && (stats[s][:ic][1] += 1)
            end
            if tgt !== nothing && isfinite(t)
                push!(stats[s][:t], abs(t-tgt.m)/tgt.m)
                abs(t-tgt.m) <= tgt.sd && (stats[s][:it][1] += 1)
            end
        end
    end
    stats
end

println("\n" * "="^104)
println("NEWTONIAN -- Gabbard et al. (2025).  Proximity metrics at h = 0.02R.")
println("="^104)
chdr, crows = read_csv(joinpath(ROOT,"julia","data","gabbard2025_restitution.csv"))
thdr, trows = read_csv(joinpath(ROOT,"julia","data","gabbard2025_contact_time.csv"))
cWe = col(chdr,crows,"We"); cY = col(chdr,crows,"cor"); cOh = col(chdr,crows,"Oh")
tWe = col(thdr,trows,"We"); tY = col(thdr,trows,"tc_over_tsigma"); tOh = col(thdr,trows,"Oh")
BANDS = [(0.018,0.030,0.0233),(0.030,0.050,0.0373),(0.050,0.100,0.0767),
         (0.20,0.35,0.2889),(0.50,0.80,0.6849)]
tot = Dict("search"=>(Float64[],Float64[],[0],[0],[0]), "lcp"=>(Float64[],Float64[],[0],[0],[0]))
for (lo,hi,Oh) in BANDS
    sc = [i for i in eachindex(cWe) if lo <= cOh[i] < hi]
    st = [i for i in eachindex(tWe) if lo <= tOh[i] < hi]
    (length(sc) < MIN_REPS || length(st) < MIN_REPS) && continue
    gc = wegroups(cWe, cY, sc, 6); gt = wegroups(tWe, tY, st, 6)
    isempty(gc) && continue
    s = report(@sprintf("Oh in [%.3f,%.3f), model at Oh=%.4f", lo, hi, Oh),
               gc, gt, ["search","lcp"], "newtonian", Oh)
    for k in ("search","lcp")
        append!(tot[k][1], s[k][:c]); append!(tot[k][2], s[k][:t])
        tot[k][3][1] += s[k][:ic][1]; tot[k][4][1] += s[k][:it][1]
        tot[k][5][1] += length(s[k][:c])
    end
end
println("\n  * = inside one experimental standard deviation")
println("\nNEWTONIAN SUMMARY")
for k in ("search","lcp")
    c,t,ic,it,n = tot[k]
    @printf("  %-7s CoR med |err| %5.1f%%  inside 1sd %2d/%2d | tc med |err| %5.1f%%  inside 1sd %2d/%2d\n",
            k, isempty(c) ? NaN : 100*median(c), ic[1], n[1],
            isempty(t) ? NaN : 100*median(t), it[1], length(t))
end

println("\n" * "="^104)
println("SHEAR THINNING -- 3000 ppm fluid.  Proximity metrics at h = 0.02R.")
println("="^104)
R,sg,rh = 0.0003,0.0728,1000.0; TCAP = sqrt(rh*R^3/sg)
OH_ST = 8.433817577956766/sqrt(rh*sg*R)
rows = NTuple{3,Float64}[]
for (i,ln) in enumerate(eachline(joinpath(ROOT,"julia","derivations","data","metrics_3000ppm.csv")))
    i == 1 && continue
    f = split(strip(ln),';'); length(f) < 3 && continue
    v = tryparse.(Float64, replace.(f[1:3], ',' => '.')); any(isnothing,v) && continue
    push!(rows, (v[1], v[2], v[3]/TCAP))
end
sWe = [r[1] for r in rows]; sC = [r[2] for r in rows]; sT = [r[3] for r in rows]
gc = wegroups(sWe, sC, collect(eachindex(sWe)), 7)
gt = wegroups(sWe, sT, collect(eachindex(sWe)), 7)
s = report("3000 ppm, one fluid", gc, gt, ["search","lcp"], "carreau_3000ppm", OH_ST)
println("\n  * = inside one experimental standard deviation")
println("\nSHEAR-THINNING SUMMARY")
for k in ("search","lcp")
    @printf("  %-7s CoR med |err| %5.1f%%  inside 1sd %2d/%2d | tc med |err| %5.1f%%  inside 1sd %2d/%2d\n",
            k, isempty(s[k][:c]) ? NaN : 100*median(s[k][:c]), s[k][:ic][1], length(s[k][:c]),
            isempty(s[k][:t]) ? NaN : 100*median(s[k][:t]), s[k][:it][1], length(s[k][:t]))
end
