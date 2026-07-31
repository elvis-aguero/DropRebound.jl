using Documenter
using Literate
using DropSolver

const DERIVATIONS_SRC = joinpath(@__DIR__, "..", "julia", "derivations")
const DERIVATIONS_OUT = joinpath(@__DIR__, "src", "derivations")

isdir(DERIVATIONS_OUT) && rm(DERIVATIONS_OUT; recursive=true)
mkpath(DERIVATIONS_OUT)

derivation_scripts = sort(filter(f -> endswith(f, ".jl"), readdir(DERIVATIONS_SRC)))

derivation_pages = Pair{String,String}[]
for script in derivation_scripts
    src = joinpath(DERIVATIONS_SRC, script)
    Literate.markdown(src, DERIVATIONS_OUT; documenter=true)
    md_name = replace(script, ".jl" => ".md")
    title = replace(splitext(script)[1], "_" => " ")
    push!(derivation_pages, title => joinpath("derivations", md_name))
end

makedocs(
    sitename = "DropRebound.jl",
    modules = [DropSolver],
    pages = [
        "Home" => "index.md",
        "CAS Derivations" => derivation_pages,
        "API Reference" => "api.md",
    ],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        inventory_version = "0.1.0",
        edit_link = "main",
    ),
)
