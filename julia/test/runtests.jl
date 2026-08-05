using Test
using DropSolver

@testset "DropSolver" begin
    ## ---- the variational formulation ------------------------------------------------
    include("test_variational.jl")            # assembly structure, Reid/Lamb limits, Gaunt band
    include("test_variational_impact.jl")     # the searching closure against the ancestor
    include("test_conservation.jl")           # energy, impulse, volume on the trajectory
    include("test_lcp_contact.jl")            # the complementarity closure
    include("test_physics_invariants.jl")     # analytic targets and two-solver agreement
    include("test_experimental_expectations.jl")  # what the high-speed video shows

    ## ---- the nonvariational formulation ---------------------------------------------
    ## Surface amplitudes with Reid's exact per-mode coefficients. Kept live as the second
    ## leg of the two-by-two in docs/src/solvers.md: it shares no contact logic and no state
    ## layout with the variational solvers, so agreement between them is real evidence.
    include("test_legendre.jl")
    include("test_integrals.jl")
    include("test_bdf.jl")
    include("test_residual.jl")
    include("test_contact.jl")
    include("test_timestepper.jl")
    include("test_newtonian.jl")
    include("test_impact.jl")
    include("test_matlab_parity.jl")
    include("test_postprocessing.jl")
    include("test_convergence_order.jl")
    include("test_nonvariational_lcp.jl")     # the fourth cell, and why it needs its own solver
    ## viscoelastic and shear-thinning extensions of the nonvariational model
    include("test_ob.jl")
    include("test_ob_eigenvalue.jl")
    include("test_carreau.jl")
    include("test_carreau_yasuda.jl")
    include("test_carreau_yasuda_exact.jl")
    include("test_oh_eff_quadrature.jl")

    ## ---- shared foundations and hygiene ---------------------------------------------
    include("test_reid.jl")
    include("test_reid_dominance.jl")
    include("test_bessel_discipline.jl")
    include("test_src_markers.jl")
    include("test_claim_ledger.jl")
end
