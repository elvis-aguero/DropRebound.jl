module DropSolver

using LinearAlgebra

# Included in order of dependency
include("types.jl")

# Export types
export SimConstants, OBParams, DropState

# Export residual functions
export pack_X, unpack_X!, build_residual!

# Export utilities used in tests
export precompute_integrals
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
