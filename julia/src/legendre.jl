"""
    collect_Pl(N, x) → Matrix{Float64}  size = (length(x), N+1)

Evaluate Legendre polynomials P₀…P_N at each point in x using the
three-term recurrence. Column index is 1-based: column k holds P_{k-1}.
Mirrors MATLAB collectPl(N, x).
"""
function collect_Pl(N::Int, x::AbstractVector{<:Real})
    n = length(x)
    P = zeros(n, N + 1)
    P[:, 1] .= 1.0          # P₀ = 1
    if N >= 1
        P[:, 2] .= x        # P₁ = x
    end
    for k in 2:N
        # (k)·P_k = (2k-1)·x·P_{k-1} - (k-1)·P_{k-2}
        @. P[:, k+1] = ((2k - 1) * x * P[:, k] - (k - 1) * P[:, k-1]) / k
    end
    return P
end

"""
    collect_dPl(N, x) → Matrix{Float64}  size = (length(x), N+1)

Evaluate first derivatives P'₀…P'_N via the identity
    (1-x²)P'ₙ = n(P_{n-1} - x·Pₙ)
Column index is 1-based: column k holds P'_{k-1}.
"""
function collect_dPl(N::Int, x::AbstractVector{<:Real})
    P  = collect_Pl(N, x)
    dP = zeros(size(P))
    for n in 1:N
        denom = 1.0 .- x .^ 2
        safe  = abs.(denom) .> 1e-14
        dP[safe,  n+1] .= n .* (P[safe, n] .- x[safe] .* P[safe, n+1]) ./ denom[safe]
        # At x = ±1: P'ₙ(±1) = ±n(n+1)/2
        dP[.!safe, n+1] .= sign.(x[.!safe]) .^ (n+1) .* n .* (n+1) ./ 2
    end
    return dP
end
