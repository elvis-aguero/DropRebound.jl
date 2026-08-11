# The solver against the experiments of Gabbard et al. (2025).
#
# The two KPIs are the coefficient of restitution and the contact time. There are 935
# experimental points spanning Oh from 0.014 to 0.79, so rather than one run per point
# the data are binned by Ohnesorge and a Weber sweep is run at each bin's median Oh
# and Bo. That is the comparison the paper itself makes -- CoR against We along
# curves of fixed Oh -- and it separates the two dependences instead of scattering
# them.
#
# Both K = 1 and K = 2 are run at every point, because the difference between them is
# the whole question: K = 1 reproduces the published Newtonian model, with Lamb's
# modal damping, and K = 2 resolves the interior and so carries Reid's. Whether that
# 16 per cent shift in CoR at We = 1 moves the model toward the measurements or away
# from them is the thing worth knowing, and it cannot be argued, only run.

using Printf
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const DATA = joinpath(@__DIR__, "..", "data")

"""Read one of the extracted experiment files, skipping its comment header."""
function read_exp(name)
    rows = Vector{Vector{Float64}}()
    hdr = String[]
    for ln in eachline(joinpath(DATA, name))
        startswith(ln, "#") && continue
        f = split(strip(ln), ',')
        if isempty(hdr)
            hdr = String.(f); continue
        end
        v = tryparse.(Float64, f)
        any(isnothing, v) && continue
        push!(rows, Float64[v...])
    end
    (hdr, reduce(vcat, (r' for r in rows)))
end

hdr_c, cor_d = read_exp("gabbard2025_restitution.csv")
hdr_t, tc_d  = read_exp("gabbard2025_contact_time.csv")
@printf("restitution: %d points, columns %s\n", size(cor_d, 1), join(hdr_c, ","))
@printf("contact time: %d points, columns %s\n", size(tc_d, 1), join(hdr_t, ","))

iW, iY, iO, iB = 1, 3, 5, 6                   # We, value, Oh, Bo in both files
@printf("Oh spans [%.4f, %.4f]   We spans [%.4g, %.4g]\n",
        minimum(cor_d[:, iO]), maximum(cor_d[:, iO]),
        minimum(cor_d[:, iW]), maximum(cor_d[:, iW]))

# --- bin by Ohnesorge ------------------------------------------------------------
const OH_EDGES = [0.0, 0.025, 0.05, 0.10, 0.20, 0.45, 1.0]

function bin_of(oh)
    for i in 1:(length(OH_EDGES)-1)
        OH_EDGES[i] <= oh < OH_EDGES[i+1] && return i
    end
    nothing
end

const M_RUN = 45      # converged to 1% in CoR against M = 90, and 15x cheaper
# A cap, not a schedule: the march stops on its own once the drop has left the
# substrate and is rising, so this only bounds runs that never release. It has to
# comfortably exceed the longest contact in the data, which reaches 5.4 at the lowest
# Weber, plus the flight back through z = 1.
const T_MAX = 25.0

results = Any[]
for bi in 1:(length(OH_EDGES)-1)
    sel = [i for i in 1:size(cor_d,1) if bin_of(cor_d[i, iO]) == bi]
    length(sel) < 12 && continue
    Oh = median(cor_d[sel, iO]); Bo = median(cor_d[sel, iB])
    Wes = cor_d[sel, iW]
    # a Weber sweep spanning this bin, logarithmic because the data are
    grid = exp.(range(log(quantile(Wes, 0.05)), log(quantile(Wes, 0.95)); length = 7))
    @printf("\n=== Oh bin %d: [%.3f, %.3f)  n=%d  median Oh=%.4f Bo=%.5f\n",
            bi, OH_EDGES[bi], OH_EDGES[bi+1], length(sel), Oh, Bo)
    @printf("%-9s | %-17s | %-17s | %s\n", "We", "CoR K=1", "CoR K=2", "experiment (binned)")
    for We in grid
        row = Any[We, Oh, Bo]
        for K in (1, 2)
            r = simulate(ImpactParams(; We = We, Bo = Bo, Oh = Oh,
                                      M = M_RUN, K = K, t_max = T_MAX))
            push!(row, r.cor, r.tc, maximum(r.cp))
        end
        # nearest experimental points in We within this Oh bin
        near = [i for i in sel if abs(log(cor_d[i, iW]) - log(We)) < 0.25]
        exp_cor = isempty(near) ? NaN : mean(cor_d[near, iY])
        exp_n = length(near)
        tnear = [i for i in 1:size(tc_d,1)
                 if bin_of(tc_d[i, iO]) == bi && abs(log(tc_d[i, iW]) - log(We)) < 0.25]
        exp_tc = isempty(tnear) ? NaN : mean(tc_d[tnear, iY])
        push!(row, exp_cor, exp_tc, exp_n)
        push!(results, row)
        @printf("%-9.4g | CoR %.4f tc %.3f | CoR %.4f tc %.3f | CoR %.4f tc %.3f (n=%d)\n",
                We, row[4], row[5], row[7], row[8], exp_cor, exp_tc, exp_n)
    end
end

# --- overall agreement -----------------------------------------------------------
println("\n================ SUMMARY ================")
for (K, ic, it) in ((1, 4, 5), (2, 7, 8))
    # A model NaN is a run that did not complete, not an agreement of zero, so both
    # sides have to be finite before a residual means anything -- and how many were
    # dropped is reported, because a silent drop reads as coverage that was not there.
    okc = [i for i in eachindex(results)
           if isfinite(results[i][10]) && isfinite(results[i][ic])]
    okt = [i for i in eachindex(results)
           if isfinite(results[i][11]) && isfinite(results[i][it])]
    rc = [abs(results[i][ic] - results[i][10])/results[i][10] for i in okc]
    rt = [abs(results[i][it] - results[i][11])/results[i][11] for i in okt]
    sc = [(results[i][ic] - results[i][10])/results[i][10] for i in okc]
    st = [(results[i][it] - results[i][11])/results[i][11] for i in okt]
    @printf("K=%d : CoR  median |err| %.1f%%  mean %.1f%%  median SIGNED %+.1f%%  (n=%d of %d)\n",
            K, 100*median(rc), 100*mean(rc), 100*median(sc), length(okc), length(results))
    @printf("      tc   median |err| %.1f%%  mean %.1f%%  median SIGNED %+.1f%%  (n=%d of %d)\n",
            100*median(rt), 100*mean(rt), 100*median(st), length(okt), length(results))
end

open(joinpath(@__DIR__, "..", "outputs", "csv", "gabbard_validation.csv"), "w") do io
    println(io, "We,Oh,Bo,cor_K1,tc_K1,cp_K1,cor_K2,tc_K2,cp_K2,cor_exp,tc_exp,n_exp")
    for r in results
        println(io, join(r, ","))
    end
end
println("\nwrote gabbard_validation.csv")
