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

# The rule above catches unmarked COMMENTS. It does not catch unmarked CODE, and
# that is the other half of the same failure: a `@printf` or `@variables` line
# added to an otherwise-hidden assertion block leaks into the rendered page and
# runs in a module where `Printf` was never imported. Both have happened; the
# second broke the build twice.
#
# Heuristic, and it is only a heuristic: inside a `let ... end` block that is
# overwhelmingly `#src`, a line WITHOUT `#src` is almost certainly an oversight.
# The threshold is deliberately high so that genuinely mixed blocks -- ones that
# intend to show code on the page -- are not flagged.
@testset "no unmarked code inside hidden blocks" begin
    # LOCAL CONTEXT, not block parsing. A first version tracked `let ... end`
    # depth, but these blocks contain `for` loops, so the first inner `end` closed
    # the block early and the check missed the very line that had just broken the
    # build. Nesting-aware parsing of Julia by regex is not worth attempting.
    #
    # The rule instead: a code line lacking `#src` whose nearest code neighbours on
    # BOTH sides carry `#src` is an oversight. A genuine on-page code block has
    # several consecutive unmarked lines, so its interior is never flagged; only an
    # isolated unmarked line among marked ones is.
    suspects = Tuple{String,Int,String}[]
    iscode(l) = !isempty(strip(l)) && !startswith(strip(l), "#")
    for f in sort(readdir(DERIV_DIR))
        endswith(f, ".jl") || continue
        lines = readlines(joinpath(DERIV_DIR, f))
        code = [i for i in eachindex(lines) if iscode(lines[i])]
        for (n, i) in enumerate(code)
            occursin("#src", lines[i]) && continue
            (n == 1 || n == length(code)) && continue   # && binds tighter than ||
            prev, nxt = code[n-1], code[n+1]
            # CONTIGUOUS neighbours only. Skipping over prose to find the nearest
            # code line flags legitimate on-page snippets -- a one-line definition
            # displayed between two paragraphs has `#src` code somewhere above and
            # below it, and is not a leak. Requiring adjacency separates "isolated
            # unmarked line inside a hidden block" from "deliberately visible code".
            (i - prev <= 2 && nxt - i <= 2) || continue
            if occursin("#src", lines[prev]) && occursin("#src", lines[nxt])
                push!(suspects, (f, i, first(strip(lines[i]), 60)))
            end
        end
    end
    @test isempty(suspects) || error(
        "code without `#src` inside an otherwise-hidden block:\n" *
        join(["  $f:$i  $t" for (f, i, t) in suspects], "\n") *
        "\n\nThese lines render into the @example block and execute in a module " *
        "that has none of the derivation's imports. Add `#src`.")
end
