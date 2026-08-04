using Test
using DropSolver

@testset "DropSolver" begin
    include("test_legendre.jl")
    include("test_integrals.jl")
    include("test_bdf.jl")
    include("test_reid.jl")
    include("test_reid_dominance.jl")
    include("test_residual.jl")
    include("test_contact.jl")
    include("test_timestepper.jl")
    include("test_newtonian.jl")
    include("test_ob.jl")
    include("test_ob_eigenvalue.jl")
    include("test_impact.jl")
    include("test_matlab_parity.jl")
    include("test_postprocessing.jl")
    include("test_carreau.jl")
    include("test_carreau_yasuda.jl")
    include("test_carreau_yasuda_exact.jl")
    include("test_oh_eff_quadrature.jl")
    include("test_convergence_order.jl")
    include("test_claim_ledger.jl")
    include("test_bessel_discipline.jl")
    include("test_src_markers.jl")
    include("test_variational.jl")
    include("test_variational_impact.jl")
end
