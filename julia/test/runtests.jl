using Test
using DropSolver

@testset "DropSolver" begin
    ## ---- the variational formulation, both contact closures --------------------------
    include("test_geometry_cache.jl")      # the cached assembly, and the radial basis ceiling
    include("test_variational.jl")            # assembly structure, Reid/Lamb limits, Gaunt band
    include("test_variational_impact.jl")     # the searching closure against the ancestor
    include("test_conservation.jl")           # energy, impulse, volume on the trajectory
    include("test_lcp_contact.jl")            # the complementarity closure
    include("test_backends.jl")              # the unified solver entry point
    include("test_nodal_forcing.jl")          # the conjugate forcing, and why it is not production
    include("test_physics_invariants.jl")     # analytic targets and two-solver agreement
    include("test_experimental_expectations.jl")  # what the high-speed video shows

    ## ---- the nonvariational formulation, searching closure --------------------------
    ## Surface amplitudes advanced with Reid's exact per-mode coefficients, no interior
    ## state. Live as the second formulation in docs/src/solvers.md: it shares neither the
    ## state layout nor the contact logic of the variational solvers, so agreement between
    ## them is evidence rather than a self-check. It also carries the Oldroyd-B extension.
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
