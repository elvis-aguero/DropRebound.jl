module DropSolver

# using LinearAlgebra  # will be used in actual implementations

# Included in order of dependency
include("types.jl")
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
