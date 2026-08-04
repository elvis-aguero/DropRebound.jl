using Documenter
using Literate
using DropSolver

include(joinpath(@__DIR__, "prettylatex.jl"))

# Figures are generated before the pages are rendered, so that the plotting
# code stays out of the documentation and the pages can reference the results
# with ordinary Markdown image syntax.
include(joinpath(@__DIR__, "figures.jl"))
build_all()

const DERIVATIONS_SRC = joinpath(@__DIR__, "..", "julia", "derivations")
const DERIVATIONS_OUT = joinpath(@__DIR__, "src", "derivations")

isdir(DERIVATIONS_OUT) && rm(DERIVATIONS_OUT; recursive=true)
mkpath(DERIVATIONS_OUT)

# Every derivation script is rendered. Titles and grouping live in `pages`
# below, so this list only has to stay in sync with the directory -- the
# assertion enforces that, so adding a script cannot silently skip the site.
const PUBLISHED = [
    "reid1960_full_derivation.jl",
    "reid_finite_oh_derivation.jl",
    "generalized_newtonian_hierarchy_derivation.jl",
    "shear_thinning_closures_derivation.jl",
    "carreau_yasuda_multimode_derivation.jl",
    "cross_fluid_derivation.jl",
    "eta_spectrum_derivation.jl",
    "oldroyd_b_derivation.jl",
]

# Superseded derivations stay in the repository and keep running in CI, but are
# not part of the documentation: they describe code paths the solver no longer
# uses, and publishing them alongside the current theory would only mislead.
const UNPUBLISHED = [
    "carreau_yasuda_derivation.jl",
    # The single-mode closure, superseded twice over: first by the multi-mode
    # version below, and then by the model itself, which evaluates the viscosity
    # pointwise rather than through any per-mode characteristic shear rate.
    "carreau_yasuda_nonperturbative_derivation.jl",
    # Its result -- that no small-amplitude expansion of this problem exists,
    # because eps^a is non-analytic at eps=0 for non-integer a -- is now derived
    # in "Shear-Thinning Drops" where it is load-bearing. A chapter whose title
    # is that something does not work is not a chapter.
    "carreauYasuda_firstprinciples_derivation.jl",
    # It opens "this is the model the validation pipeline runs", and it is not: it
    # builds each mode's viscosity from that mode's own amplitude through a scalar
    # characteristic shear rate. The model computes eta pointwise from the shear rate
    # of the full summed strain field, and that the invariant does NOT superpose over
    # modes is the reason the coupling exists at all. A page that presents the
    # superseded closure as the live one is the most misleading kind of stale.
    "carreau_yasuda_multimode_derivation.jl",
]

all_scripts = Set(filter(f -> endswith(f, ".jl"), readdir(DERIVATIONS_SRC)))
@assert Set(vcat(PUBLISHED, UNPUBLISHED)) == all_scripts "docs/make.jl is out of sync with julia/derivations/*.jl -- add the new/renamed script to PUBLISHED or UNPUBLISHED"

# script name -> path of its rendered page, for use in `pages`
const PAGE = Dict{String,String}()
for script in PUBLISHED
    Literate.markdown(joinpath(DERIVATIONS_SRC, script), DERIVATIONS_OUT; documenter=true)
    PAGE[script] = joinpath("derivations", replace(script, ".jl" => ".md"))
end

# A `##` comment inside a block whose every executable line is `#src` makes
# Literate emit an @example box containing only comments -- a grey rectangle
# with nothing in it. This has slipped through repeatedly, so fail the build.
for md in readdir(DERIVATIONS_OUT; join=true)
    endswith(md, ".md") || continue
    body = read(md, String)
    for m in eachmatch(r"````@example[^\n]*\n(.*?)````"s, body)
        code = filter(l -> !isempty(strip(l)) && !startswith(strip(l), "#"),
                      split(m.captures[1], "\n"))
        isempty(code) && error("""
            Empty code block in $(basename(md)): an @example block contains only
            comments. Mark the offending `##` lines with a trailing `#src`.
            Block began: $(first(split(strip(m.captures[1]), "\n")))""")
    end
end

makedocs(
    sitename = "DropRebound.jl",
    modules = [DropSolver],
    repo = Remotes.GitHub("elvis-aguero", "DropRebound.jl"),
    pages = [
        "Home" => "index.md",
        # Grouped by the question a reader arrives with, not by filename or by
        # the order things happened to be written. Each page is a derivation
        # that also runs in CI, so the physics on the site is executable.
        #
        # (The hand-authored `reid1960/` and `carreau_yasuda_fp/` chapters were
        # removed: they restated, in a second voice, what these scripts derive.)
        "Newtonian Theory" => [
            "The Viscous Drop: Reid (1960)"    => PAGE["reid1960_full_derivation.jl"],
            "Finite-Ohnesorge Coefficients"    => PAGE["reid_finite_oh_derivation.jl"],
        ],
        "Shear-Thinning Fluids" => [
            # Read in this order: the model first, stated completely and
            # without approximation; then the descent from it to something
            # runnable; then the constitutive law a real fluid is characterised
            # by; and last the measurement that prices the mode-coupling
            # concessions the descent makes.
            "Shear-Thinning Drops"             => PAGE["generalized_newtonian_hierarchy_derivation.jl"],
            "Shear-Thinning Drops: Closures"   => PAGE["shear_thinning_closures_derivation.jl"],
            "Cross-Model Fluids"               => PAGE["cross_fluid_derivation.jl"],
            "Angular Bandwidth of Viscosity"   => PAGE["eta_spectrum_derivation.jl"],
        ],
        "Viscoelastic Fluids" => [
            "Oldroyd-B"                        => PAGE["oldroyd_b_derivation.jl"],
        ],
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        inventory_version = "0.1.0",
        edit_link = "main",
        canonical = "https://elvis-aguero.github.io/DropRebound.jl",
    ),
)

deploydocs(
    repo = "github.com/elvis-aguero/DropRebound.jl.git",
    devbranch = "main",
    push_preview = true,
)
