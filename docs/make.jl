using Documenter
using Literate
using DropSolver

include(joinpath(@__DIR__, "prettylatex.jl"))

const DERIVATIONS_SRC = joinpath(@__DIR__, "..", "julia", "derivations")
const DERIVATIONS_OUT = joinpath(@__DIR__, "src", "derivations")

isdir(DERIVATIONS_OUT) && rm(DERIVATIONS_OUT; recursive=true)
mkpath(DERIVATIONS_OUT)

# Explicit reading order (not alphabetical): the physics foundation first,
# then its production-wiring continuation, then each rheology extension in
# the order it builds on what came before. `carreau_yasuda_derivation.jl`
# is EXPLICITLY marked legacy/superseded in its title -- see
# `julia/derivations/carreau_yasuda_multimode_derivation.jl`'s header for
# why (it derives `julia/src/st_extension.jl`, which the real validation
# pipeline no longer uses).
const READING_ORDER = [
    "reid1960_full_derivation.jl" => "Part 1 -- Reid (1960): the viscous drop, from scratch",
    "reid_finite_oh_derivation.jl" => "Part 2 -- Finite-Oh damping and frequency (production wiring)",
    "oldroyd_b_derivation.jl" => "Part 3 -- Oldroyd-B viscoelasticity",
    "carreau_yasuda_nonperturbative_derivation.jl" => "Part 4 -- Carreau-Yasuda, exact single-mode",
    "carreau_yasuda_multimode_derivation.jl" => "Part 5 -- Carreau-Yasuda, multi-mode coupling (current model)",
    "cross_fluid_derivation.jl" => "Part 6 -- Cross-model fluids",
    "carreau_yasuda_derivation.jl" => "Part 7 -- Carreau-Yasuda, weakly-nonlinear (legacy, superseded by Part 5)",
    "carreauYasuda_firstprinciples_derivation.jl" => "Part 8 -- Carreau-Yasuda, why amplitude perturbation theory fails",
]

all_scripts = Set(filter(f -> endswith(f, ".jl"), readdir(DERIVATIONS_SRC)))
ordered_scripts = [p for p in READING_ORDER if p.first in all_scripts]
@assert Set(first.(ordered_scripts)) == all_scripts "READING_ORDER in docs/make.jl is out of sync with julia/derivations/*.jl -- add the new/renamed script to READING_ORDER"

derivation_pages = Pair{String,String}[]
for (script, title) in ordered_scripts
    src = joinpath(DERIVATIONS_SRC, script)
    Literate.markdown(src, DERIVATIONS_OUT; documenter=true)
    md_name = replace(script, ".jl" => ".md")
    push!(derivation_pages, title => joinpath("derivations", md_name))
end

makedocs(
    sitename = "DropRebound.jl",
    modules = [DropSolver],
    repo = Remotes.GitHub("elvis-aguero", "DropRebound.jl"),
    pages = [
        "Home" => "index.md",
        "Reid (1960) [pilot: hand-authored prose]" => [
            "Introduction" => "reid1960/01-introduction.md",
            "Mathematical Preliminaries" => "reid1960/02-preliminaries.md",
            "Problem Setup" => "reid1960/03-problem-setup.md",
            "Linearized Governing Equations" => "reid1960/04-linearized-equations.md",
            "The Pressure Field" => "reid1960/05-pressure-field.md",
            "The Velocity ODE" => "reid1960/06-velocity-ode.md",
            "Boundary Conditions" => "reid1960/07-boundary-conditions.md",
            "The Characteristic Equation" => "reid1960/08-characteristic-equation.md",
            "Structure of the Solutions" => "reid1960/09-structure-of-solutions.md",
            "Molaček & Bush Connection + Summary" => "reid1960/10-molacek-bush-summary.md",
        ],
        "Carreau-Yasuda [pilot: hand-authored prose]" => [
            "Why Amplitude Perturbation Theory Fails" => "carreau_yasuda_fp/01-why-perturbation-fails.md",
            "What Survives" => "carreau_yasuda_fp/02-what-survives.md",
        ],
        "CAS Derivations" => derivation_pages,
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
