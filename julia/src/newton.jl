"""Cache mapping (cp, dt, order[, is_ob]) → precomputed J⁻¹."""
const _jac_cache = Dict{Any, Matrix{Float64}}()

"""
    newton_solve!(X, R!, J_fn; cache_key, tol, maxiter) → Bool

Newton-Raphson solver. Returns true if converged (norm(R) < tol).
X is updated in-place to the best iterate found.
"""
function newton_solve!(X::Vector{Float64},
                       R!::Function,
                       J_fn::Function;
                       cache_key::Union{Nothing,Tuple} = nothing,
                       tol::Float64 = 1e-10,
                       maxiter::Int = 100)
    buf = similar(X)
    R!(buf, X)
    best_val = norm(buf)
    best_X   = copy(X)

    for iter in 1:maxiter
        if cache_key !== nothing && haskey(_jac_cache, cache_key)
            Jinv = _jac_cache[cache_key]
            δX   = Jinv * buf
        else
            J = J_fn(X)
            if cache_key !== nothing
                Jinv = inv(J)
                _jac_cache[cache_key] = Jinv
                @debug "Jacobian cached" key=cache_key
                δX = Jinv * buf
            else
                δX = J \ buf
            end
        end

        X .= X .- δX
        R!(buf, X)
        val = norm(buf)
        if val < best_val
            best_val = val
            best_X  .= X
        end
        if val < tol || norm(δX) < 1e-12
            break
        end
    end

    X .= best_X
    return best_val < tol
end

"""Clear the Jacobian cache (call when cp or dt changes)."""
clear_jac_cache!() = empty!(_jac_cache)
