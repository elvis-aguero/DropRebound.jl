# NONLINEAR RECOVERY, stage 2: the posterior.
#
# Reads the forward grid written by `study_recovery_nonlinear.jl` and asks the question
# the linearised studies could not: given a noisy synthetic experiment, where does the
# likelihood actually put the fluid, and is that anywhere near the truth?
#
# The grid point nearest the true parameters supplies the noiseless data. Noise is added,
# chi^2 is evaluated at EVERY grid point, and the posterior weight exp(-chi^2/2) is
# normalised over the grid. That gives marginal distributions directly, without assuming
# they are Gaussian and without an optimiser that could stop in a local minimum and
# report success.
#
# Hundreds of realisations, because the grid is already paid for and every extra
# synthetic experiment after it costs nothing.
#
# WHAT TO LOOK FOR
#   bias      the posterior mean minus the truth: does the inversion sit on the answer
#   spread    the posterior standard deviation: does it agree with the linearised bars
#   shape     whether the weight is concentrated or smeared along a valley
#
# A recovery that is unbiased with a spread matching the Fisher prediction means the
# linearised error bars were honest. A spread much larger, or a posterior with weight in
# two places, means they were not.

using Printf, LinearAlgebra, Random

const OUT = joinpath(@__DIR__, "..", "..", "results")
const LOBS = 2:8
const NT   = 50
const SIG_L = [0.005 * (l / 2) for l in LOBS]
const SIG_Z = 0.005
const K_TRUE, M_TRUE = 1.076060672638565, 0.6724627161400011
const RATIO_TRUE = 0.002747884967209792 / 0.2601531742575956
const NREAL = 400

## noise vector matching the observable layout: seven modes then the centre of mass
const SIGMA_VEC = let s = Float64[]
    for i in eachindex(SIG_L); append!(s, fill(SIG_L[i], NT)); end
    append!(s, fill(SIG_Z, NT)); s
end

rows = Vector{Vector{Float64}}()
open(joinpath(OUT, "recovery_grid.csv")) do io
    readline(io)
    for ln in eachline(io)
        v = parse.(Float64, split(strip(ln), ','))
        length(v) == 3 + length(SIGMA_VEC) && push!(rows, v)
    end
end
isempty(rows) && error("no usable rows in recovery_grid.csv")

P = [r[1:3] for r in rows]                  # parameter triples
O = [r[4:end] for r in rows]                # observable vectors
@printf("== nonlinear recovery, stage 2 ==\n%d grid points, %d observables each\n",
        length(P), length(SIGMA_VEC))

truth = [log10(K_TRUE), M_TRUE, log10(RATIO_TRUE)]
itrue = argmin([sum(abs2, (p .- truth) ./ [0.5, 0.12, 0.06]) for p in P])
@printf("truth        log10 k %+.4f  m %.4f  log10 ratio %+.4f\n", truth...)
@printf("nearest node log10 k %+.4f  m %.4f  log10 ratio %+.4f\n", P[itrue]...)
data0 = O[itrue]

rng = Random.Xoshiro(20260809)
names = ["log10 k", "m", "log10 ratio"]
means = [Float64[] for _ in 1:3]
maps  = [Float64[] for _ in 1:3]
nmulti = 0

for _ in 1:NREAL
    d = data0 .+ SIGMA_VEC .* randn(rng, length(SIGMA_VEC))
    chi2 = [sum(abs2, (d .- o) ./ SIGMA_VEC) for o in O]
    w = exp.(-(chi2 .- minimum(chi2)) ./ 2)
    w ./= sum(w)
    for j in 1:3
        push!(means[j], sum(w[i] * P[i][j] for i in eachindex(P)))
        push!(maps[j], P[argmin(chi2)][j])
    end
    ## a posterior with weight in two separated places is the pathology worth catching
    heavy = findall(>(0.01), w)
    length(heavy) > 1 && maximum(norm(P[a] .- P[b]) for a in heavy, b in heavy) > 0.9 &&
        (global nmulti += 1)
end

sd(v) = (m = sum(v)/length(v); sqrt(sum(abs2, v .- m)/max(1, length(v)-1)))

flushed = @sprintf("\n%d synthetic experiments\n", NREAL); print(flushed)
@printf("%-14s %10s %10s %10s %10s\n", "parameter", "truth", "post.mean", "bias", "post.sd")
for j in 1:3
    m = sum(means[j]) / length(means[j])
    @printf("%-14s %10.4f %10.4f %+10.4f %10.4f\n", names[j], P[itrue][j], m,
            m - P[itrue][j], sd(means[j]))
end

@printf("\nrealisations with posterior weight in two separated regions: %d of %d\n",
        nmulti, NREAL)
@printf("grid spacing: log10 k %.3f, m %.3f, log10 ratio %.3f\n",
        P[2][1] == P[1][1] ? NaN : 0.475, 0.1125, 0.0625)
println("""
Read the spreads against the grid spacing: a posterior standard deviation smaller than
one grid step is resolution-limited here and should be treated as an upper bound, not a
measurement. The linearised study (three Weber numbers) gave 0.367, 0.083 and 0.017;
this run uses one Weber number, so its bars should be about sqrt(3) wider before the two
are compared.""")
