# References to chapters that no longer exist under the name given.
#
# The site refers to its own chapters by italicised title, as in *Variational
# Mechanics*. Documenter validates `[text](@ref)` links and file links; it has no
# notion of an italicised title, so every such reference is unchecked. Renaming or
# archiving a chapter therefore breaks them silently.
#
# Both have happened. Archiving `carreau_yasuda_multimode_derivation.jl` left two
# live pages pointing at *Carreau-Yasuda: Multi-Mode*. Renaming Reid's chapter left
# *The Viscous Drop: Reid (1960)* dangling in the closures chapter.
#
# WHY THIS IS A REGISTRY AND NOT A HEURISTIC. The obvious test is to flag every
# italicised Title Case phrase that does not match a live chapter. That produces
# about seventy hits on this corpus, nearly all of them ordinary emphasis: *Proof*,
# *Assumption*, *In plain English*. An allowlist that large stops being read.
#
# An earlier version of this file checked live pages against the H1 headings of
# archived scripts. It passed, and it was vacuous: prose refers to chapters by
# short names, and no short name equals an H1. It was caught by reinstating a known
# dangling reference and watching the test go on passing.
#
# So the list below is maintained by hand, and it must be: adding to it is part of
# renaming or archiving a chapter, which is a deliberate act by a person.

using Test

const DOCS_ROOT = normpath(joinpath(@__DIR__, "..", "..", "docs"))
const DERIV_DIR = normpath(joinpath(@__DIR__, "..", "derivations"))

"Chapter titles that once appeared on the site and no longer do."
const RETIRED_TITLES = [
    # archived: the model evaluates eta pointwise, not per mode
    "Carreau-Yasuda: Multi-Mode",
    "Carreau-Yasuda: Single-Mode",
    # renamed when Part II was reframed around the eigenvalue problem
    "The Viscous Drop: Reid (1960)",
    # archived alongside the weakly nonlinear route
    "Why Amplitude Expansion Fails",
    "The Differential Formulation",
]

function published_scripts()
    mk = read(joinpath(DOCS_ROOT, "make.jl"), String)
    m = match(r"const PUBLISHED = \[(.*?)\]"s, mk)
    m === nothing ? String[] :
        [c.captures[1] for c in eachmatch(r"\"([a-zA-Z0-9_]+\.jl)\"", m.captures[1])]
end

@testset "no live page cites a retired chapter title" begin
    pubs = published_scripts()
    @test !isempty(pubs)

    live = [joinpath(DOCS_ROOT, "src", f)
            for f in readdir(joinpath(DOCS_ROOT, "src")) if endswith(f, ".md")]
    append!(live, [joinpath(DERIV_DIR, s) for s in pubs])
    push!(live, joinpath(DOCS_ROOT, "make.jl"))
    @test length(live) > 5

    offences = Tuple{String,Int,String}[]
    for path in live
        isfile(path) || continue
        for (i, line) in enumerate(eachline(path))
            for t in RETIRED_TITLES
                occursin(t, line) && push!(offences, (basename(path), i, t))
            end
        end
    end

    # Physical meaning of a failure: a reader is being sent to a chapter that is
    # not on the site under that name. Repoint the reference, or restore the
    # chapter and drop it from RETIRED_TITLES.
    for (f, l, t) in offences
        @info "live page cites a retired chapter" file = f line = l chapter = t
    end
    @test isempty(offences)
end
