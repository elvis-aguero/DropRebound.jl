using Test
using DropSolver

@testset "DropSolver" begin
    include("test_legendre.jl")
    include("test_integrals.jl")
    include("test_bdf.jl")
    include("test_residual.jl")
    include("test_contact.jl")
    include("test_newtonian.jl")
    include("test_ob.jl")
end
