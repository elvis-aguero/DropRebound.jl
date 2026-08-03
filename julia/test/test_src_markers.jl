# Literate `#src` hygiene.
#
# `docs/make.jl` already fails the build when an `@example` block renders with
# nothing but comments in it -- which happens when a `##` comment sits in a block
# whose every executable line is marked `#src`. That guard is correct but it only
# fires after a full Literate + Documenter build, which is the slowest job in CI.
#
# This is the same rule as a text scan, so the failure surfaces in milliseconds
# instead of seventeen minutes. It has caught the identical mistake twice:
# `@variables` lines that leaked into the rendered page, and `##` comments added
# to an all-`#src` assertion block by a hurried edit.
#
# THE RULE. In a derivation script, every `##`-prefixed line -- Literate's
# "comment that stays inside the code block" -- must carry a trailing `#src`.
# Prose belongs in single-`#` lines, which Literate renders as Markdown; `##`
# lines exist only to annotate code, and all the code in these files is `#src`.

using Test

const DERIV_DIR = joinpath(@__DIR__, "..", "derivations")

@testset "Literate #src markers" begin
    offenders = Tuple{String,Int,String}[]
    for f in sort(readdir(DERIV_DIR))
        endswith(f, ".jl") || continue
        for (i, line) in enumerate(readlines(joinpath(DERIV_DIR, f)))
            s = strip(line)
            startswith(s, "##") || continue
            occursin("#src", line) && continue
            push!(offenders, (f, i, first(s, 70)))
        end
    end

    @test isempty(offenders) || error(
        "`##` comment lines without a trailing `#src`:\n" *
        join(["  $f:$i  $t" for (f, i, t) in offenders], "\n") *
        "\n\nThese render inside the @example block. If every executable line in " *
        "that block is `#src`, the block renders as comments only and " *
        "docs/make.jl fails the build. Add `#src`, or move the text to a " *
        "single-`#` prose line.")

    # Guard against the rule being satisfied vacuously by a file with no `##`
    # comments at all -- the annotation style is used heavily and should stay.
    total = sum(count(l -> startswith(strip(l), "##"),
                      readlines(joinpath(DERIV_DIR, f)))
                for f in readdir(DERIV_DIR) if endswith(f, ".jl"))
    @test total > 50
    @info "src markers" annotated_lines = total
end
