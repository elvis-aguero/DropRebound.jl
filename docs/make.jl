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
    "carreauYasuda_firstprinciples_derivation.jl",
    "carreau_yasuda_nonperturbative_derivation.jl",
    "carreau_yasuda_multimode_derivation.jl",
    "cross_fluid_derivation.jl",
    "oldroyd_b_derivation.jl",
]

# Superseded derivations stay in the repository and keep running in CI, but are
# not part of the documentation: they describe code paths the solver no longer
# uses, and publishing them alongside the current theory would only mislead.
const UNPUBLISHED = [
    "carreau_yasuda_derivation.jl",
]

all_scripts = Set(filter(f -> endswith(f, ".jl"), readdir(DERIVATIONS_SRC)))
@assert Set(vcat(PUBLISHED, UNPUBLISHED)) == all_scripts "docs/make.jl is out of sync with julia/derivations/*.jl -- add the new/renamed script to PUBLISHED or UNPUBLISHED"

# script name -> path of its rendered page, for use in `pages`
const PAGE = Dict{String,String}()
for script in PUBLISHED
    Literate.markdown(joinpath(DERIVATIONS_SRC, script), DERIVATIONS_OUT; documenter=true)
    PAGE[script] = joinpath("derivations", replace(script, ".jl" => ".md"))
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
            # Read in this order: the map first, then why the obvious route is
            # closed, then the two closures that are actually implemented.
            "A Hierarchy of Models"            => PAGE["generalized_newtonian_hierarchy_derivation.jl"],
            "Why Amplitude Expansion Fails"    => PAGE["carreauYasuda_firstprinciples_derivation.jl"],
            "Carreau-Yasuda: Single-Mode"      => PAGE["carreau_yasuda_nonperturbative_derivation.jl"],
            "Carreau-Yasuda: Multi-Mode"       => PAGE["carreau_yasuda_multimode_derivation.jl"],
            "Cross-Model Fluids"               => PAGE["cross_fluid_derivation.jl"],
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
