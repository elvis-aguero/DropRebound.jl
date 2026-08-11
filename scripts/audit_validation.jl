# The backends against the impact experiments of Gabbard et al. (2025).
#
# 935 measurements spanning Oh in [0.014, 0.79]. Rather than one run per point,
# the data are binned by Ohnesorge and a Weber sweep is run at each bin's median
# Oh and Bo, which is the comparison the paper itself makes and separates the two
# dependences instead of scattering them.
#
# Each backend is run at the resolution it can actually be trusted at, which is
# not the same for all three: see `outputs/csv/audit_cost.csv`.

using Printf, Statistics
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const DATA = joinpath(@__DIR__, "..", "data")
const OUT  = joinpath(@__DIR__, "..", "outputs", "csv")

function read_exp(name)
    rows = Vector{Vector{Float64}}(); hdr = String[]
    for ln in eachline(joinpath(DATA, name))
        startswith(ln, "#") && continue
        f = split(strip(ln), ',')
        if isempty(hdr); hdr = String.(f); continue; end
        v = tryparse.(Float64, f); any(isnothing, v) && continue
        push!(rows, Float64[v...])
    end
    reduce(vcat, (r' for r in rows))
end

cor_d = read_exp("gabbard2025_restitution.csv")
tc_d  = read_exp("gabbard2025_contact_time.csv")
iW, iY, iO, iB = 1, 3, 5, 6

const OH_EDGES = [0.0, 0.025, 0.05, 0.10, 0.20, 0.45, 1.0]
bin_of(oh) = findfirst(i -> OH_EDGES[i] <= oh < OH_EDGES[i+1], 1:length(OH_EDGES)-1)

const RUNS = [(Backend(formulation=:variational, contact=:active_set), 45),
              (Backend(formulation=:variational, contact=:lcp),        45),
              (Backend(formulation=:nonvariational, contact=:tangency), 30)]

open(joinpath(OUT, "audit_validation.csv"), "w") do io
    println(io, "backend,M,Oh_bin,Oh,Bo,We,n_cor,n_tc,cor_sim,cor_exp,tc_sim,tc_exp,ok")
    for (b, M) in RUNS, bi in 1:(length(OH_EDGES)-1)
        selc = [i for i in axes(cor_d,1) if bin_of(cor_d[i,iO]) == bi]
        selt = [i for i in axes(tc_d,1)  if bin_of(tc_d[i,iO])  == bi]
        length(selc) < 12 && continue
        Oh = median(cor_d[selc, iO]); Bo = median(cor_d[selc, iB])
        Wes = exp10.(range(log10(quantile(cor_d[selc,iW], 0.1)),
                           log10(quantile(cor_d[selc,iW], 0.9)); length = 6))
        for We in Wes
            r = run_impact(b; We=We, Bo=Bo, Oh=Oh, M=M, K=3, t_max=25.0)
            ## nearest experimental points in We, within a factor of 1.3
            nc = [i for i in selc if 1/1.3 <= cor_d[i,iW]/We <= 1.3]
            nt = [i for i in selt if 1/1.3 <= tc_d[i,iW]/We  <= 1.3]
            ce = isempty(nc) ? NaN : median(cor_d[nc, iY])
            te = isempty(nt) ? NaN : median(tc_d[nt, iY])
            ## A run that never releases reports tc = t_max and CoR = 0; it is not a
            ## measurement of anything and must not enter a median as if it were.
            good = r.ok && r.tc < 0.9*25.0 && r.cor > 1e-6
            println(io, join((label(b), M, bi, round(Oh,digits=4), round(Bo,digits=4),
                              round(We,digits=4), length(nc), length(nt),
                              r.cor, ce, r.tc, te, good), ","))
            @printf("  %-16s bin%d Oh=%.3f We=%.3f  cor %.3f/%.3f  tc %.2f/%.2f %s\n",
                    label(b), bi, Oh, We, r.cor, ce, r.tc, te, good ? "" : "  REJECTED")
            flush(io)
        end
    end
end
println("done")
