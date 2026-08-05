module DropSolver

using LinearAlgebra
using Logging

# Included in order of dependency
include("reid.jl")
include("types.jl")

# Export types
export SimConstants, OBParams, STParams, DropState, make_theta_vec, make_dt_max

# Export Reid finite-Oh viscous model
export reid_char, dominant_root, second_root, reid_lambda_omega2, drop_viscous_coeffs
export ReidTable, build_reid_table, reid_lambda_omega2_fast, build_reid_cache

# Export residual and Jacobian functions
export pack_X, unpack_X!, build_residual!, build_jacobian

# Export OB extension functions
export pack_X_ob, unpack_X_ob!, build_residual_ob!, build_jacobian_ob

# Export Carreau (shear-thinning) extension functions
export build_residual_st!, build_jacobian_st

# Export non-perturbative Carreau-Yasuda extension functions
export characteristic_shear_K, STExactParams, oh_eff_all_coupled, lambda_omega2_from_oh_eff,
       build_residual_st_exact!, build_jacobian_st_exact, OUTLIER_FACTOR

# Export Newton solver
export newton_solve!, clear_jac_cache!

# Export utilities used in tests
export precompute_integrals

# Export contact detection functions
export drop_height, contact_error, update_theta_star!

# Export time-steppers
export solve_drop!, solve_drop_v1!
export solve_drop_lcp!, nv_compliance, nv_gap, nv_nodal_pressure
export lcp_active_set, lcp_residual

# Export v1 (continuous θ*) residual/Jacobian
export build_residual_v1!, build_jacobian_v1, integral_at_theta_star

include("legendre.jl")
include("integrals.jl")
include("bdf.jl")
include("residual.jl")
include("jacobian.jl")
include("ob_extension.jl")
include("st_extension.jl")
include("st_exact_extension.jl")
include("newton.jl")
include("contact.jl")
include("timestepper.jl")
include("residual_v1.jl")
include("timestepper_v1.jl")
include("postprocessing.jl")
include("nonvariational_lcp.jl")

# The variational assembly of the shear-thinning model (see the derivation page
# "Shear-Thinning Drops"). Kept separate from the eliminated-interior solver above:
# this one retains the interior as part of the state.
include("variational.jl")
export RitzBasis, assemble, decay_rates, dominant_pair
export ModalBasis, ndof, assemble_coupled, assemble_newtonian, strain_at, shear_rate, block_norm, carreau
include("variational_solve.jl")
export ImpactParams, ImpactState, simulate, simulate_lcp, proximity_metrics, try_step_lcp, lcp_pgs, initial_state, surface_amplitudes
# Export postprocessing
export SweepKPIs, extract_kpis, compute_contact_radius, drop_profile

end # module
