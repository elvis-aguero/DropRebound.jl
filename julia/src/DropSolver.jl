module DropSolver

using LinearAlgebra

# Included in order of dependency
include("types.jl")

# Export types
export SimConstants, OBParams, DropState

# Export residual and Jacobian functions
export pack_X, unpack_X!, build_residual!, build_jacobian

# Export Newton solver
export newton_solve!, clear_jac_cache!

# Export utilities used in tests
export precompute_integrals

# Export contact detection functions
export drop_height, contact_error
include("legendre.jl")
include("integrals.jl")
include("bdf.jl")
include("residual.jl")
include("jacobian.jl")
include("ob_extension.jl")
include("newton.jl")
include("contact.jl")
include("timestepper.jl")

end # module
