"""
    precompute_integrals(angles_in, N) → (M_mat, angles_out)

Compute the integral ∫_{cos(aᵢ)}^{cos(a_{i+1})} Pₙ(u)/u³ du for each
consecutive pair of angles and each n = 0…N.

`angles_in`: vector of angles in (π/2, π], or NaN for N+1 uniform samples.
Returns:
- `M_mat`:      (n_intervals × N+1) matrix
- `angles_out`: the processed angle vector
"""
function precompute_integrals(angles_in, N::Int)
    if angles_in isa Number && isnan(angles_in)
        angles = range(π, 0; length = N + 2) |> collect
    elseif angles_in isa Number
        angles = range(π, 0; length = Int(angles_in)) |> collect
    else
        angles = collect(Float64, angles_in)
    end

    # Keep only angles in (π*0.51, π]
    angles = filter(θ -> θ > π * 0.51 && θ <= π, angles)

    n_intervals = length(angles) - 1
    M_mat = zeros(n_intervals, N + 1)

    for i in 1:n_intervals
        a = cos(angles[i])
        b = cos(angles[i + 1])
        # Gauss-Legendre over [a, b] (a < b since cos is decreasing on (0,π))
        n_quad = 50
        t, w = gauss_legendre_nodes(n_quad, a, b)
        Pt = collect_Pl(N, t)   # (n_quad × N+1)
        integrand = Pt ./ (t .^ 3)
        M_mat[i, :] = integrand' * w
    end

    return M_mat, angles
end

"""
    integral_at_theta_star(theta_star, M; n_quad=60) → Vector{Float64}

Compute ∫_{-1}^{cos(θ*)} Pₙ(u)/u³ du for n = 0…M.
Returns a length-(M+1) vector.  Valid for θ* ∈ (π/2, π) so the upper
limit cos(θ*) ∈ (-1, 0) — the singularity at u=0 is never reached.
"""
function integral_at_theta_star(theta_star::Float64, M::Int; n_quad::Int=60)
    a = -1.0
    b = cos(theta_star)
    t, w = gauss_legendre_nodes(n_quad, a, b)
    Pt = collect_Pl(M, t)          # (n_quad × M+1)
    integrand = Pt ./ (t .^ 3)    # (n_quad × M+1)
    return integrand' * w          # length M+1
end

"""Gauss-Legendre nodes and weights on [a, b] via Golub-Welsch."""
function gauss_legendre_nodes(n::Int, a::Float64, b::Float64)
    β = [k / sqrt(4k^2 - 1) for k in 1:n-1]
    T = diagm(1 => β, -1 => β)
    vals, vecs = eigen(Symmetric(T))
    x = vals
    w = 2 .* vecs[1, :] .^ 2
    xab = @. (b - a) / 2 * x + (b + a) / 2
    wab = @. (b - a) / 2 * w
    return xab, wab
end
