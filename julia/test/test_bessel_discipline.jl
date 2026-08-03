# Spherical Bessel discipline.
#
# `julia/src/reid.jl` documents the trap: `j_l(q)` OVERFLOWS at small Ohnesorge,
# where `q = sqrt(sigma/Oh)` has large `|Im q|`, while the ratio
# `Q = j_{l+1}/j_l` stays O(1). `sph_bessel_ratio` exists for exactly that
# reason and never evaluates a Bessel function directly.
#
# There is a second failure, subtler and worse because it is silent:
# `sqrt(pi/(2z)) * besselj(l + 1/2, z)` does not agree with the entire function
# `j_l` for complex `z` with negative imaginary part -- the branch cut of the
# prefactor and the half-integer order do not cancel there. A check built on it
# can return a small residual for the WRONG reason. That happened: a
# well-posedness test reported agreement to 1e-12 with a hand-rolled evaluator
# and disagreement at 4e-1 once it used `sph_bessel_ratio`, and the hand-rolled
# answer was the wrong one.
#
# The sister repository SpectralKM states the same rule from the other side:
# "hand-rolling a robust Bessel function is materially riskier than hand-rolling
# Legendre ... this is a correctness-critical dependency rather than a style one."
#
# THE RULE. Do not write your own spherical Bessel evaluation. Use
# `sph_bessel_ratio` where a ratio suffices -- it almost always does, because the
# characteristic equation and every boundary condition are homogeneous in `j_l`.
# Where a value is genuinely needed at complex argument, use the SCALED
# `besseljx` and form the quantity as a ratio so the exponential factors cancel,
# as `eta_spectrum_derivation.jl` does.

using Test

const SRC_DIRS = [joinpath(@__DIR__, "..", "src"),
                  joinpath(@__DIR__, "..", "derivations")]

# The forbidden construction: an unscaled half-integer besselj turned into a
# spherical Bessel by a sqrt(pi/2z) prefactor.
const FORBIDDEN = r"sqrt\(\s*pi\s*/\s*\(?\s*2\s*\*?\s*z"i

# Known debt, with the reason. These predate the rule, live on UNPUBLISHED pages,
# and are called at COMPLEX argument -- so they carry the defect above, not an
# exemption from it. Listed so the rule can be enforced now and the debt paid
# separately. Do not add to this list to make a new check pass.
const KNOWN_DEBT = Set([
    "carreau_yasuda_derivation.jl",             # its own reid_char, complex qv
    "carreauYasuda_firstprinciples_derivation.jl",  # sph_jl at complex q0
])

@testset "spherical Bessel discipline" begin
    offenders = String[]
    for d in SRC_DIRS
        isdir(d) || continue
        for f in readdir(d)
            endswith(f, ".jl") || continue
            body = read(joinpath(d, f), String)
            # only flag code, not the prose that explains the rule
            code = join([l for l in split(body, "\n")
                         if !startswith(strip(l), "#") || occursin("#src", l)], "\n")
            occursin(FORBIDDEN, code) && push!(offenders, f)
        end
    end

    new_offenders = setdiff(Set(offenders), KNOWN_DEBT)
    @test isempty(new_offenders) ||
          error("hand-rolled spherical Bessel evaluation in: " *
                join(sort(collect(new_offenders)), ", ") *
                "\nUse sph_bessel_ratio for ratios, or scaled besseljx formed as a " *
                "ratio for values. See this file's header for why.")

    # The debt list must not rot: if a file is fixed, drop it from the list, so
    # the list can only shrink.
    stale = setdiff(KNOWN_DEBT, Set(offenders))
    @test isempty(stale) ||
          error("KNOWN_DEBT lists files that no longer offend: " *
                join(sort(collect(stale)), ", ") * ". Remove them.")

    @info "bessel discipline" offenders = sort(offenders) debt = length(KNOWN_DEBT)
end
