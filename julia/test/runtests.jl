using Test
using DropSolver

@testset "DropSolver" begin
    ## the live variational model: assembly, the two contact algorithms, and the physics
    include("test_variational.jl")
    include("test_variational_impact.jl")
    include("test_conservation.jl")
    include("test_lcp_contact.jl")
    include("test_physics_invariants.jl")
    include("test_experimental_expectations.jl")
    ## Reid's dispersion relation, which the variational limits are checked against
    include("test_reid.jl")
    include("test_reid_dominance.jl")
    ## hygiene over the published derivations
    include("test_bessel_discipline.jl")
    include("test_src_markers.jl")
    include("test_claim_ledger.jl")
end
