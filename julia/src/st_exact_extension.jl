"""
Non-perturbative Carreau-Yasuda shear-thinning extension. Derivation:
`julia/derivations/carreau_yasuda_nonperturbative_derivation.jl`.

Unlike `st_extension.jl` (which Taylor-expands the Carreau-Yasuda viscosity
law to first order in `eps_ST` -- invalid once `eps_ST` is not small, which
is the case for the real shear-thinning validation fluid this was built for,
`eps_ST ~ 0.9996`), this computes an effective Ohnesorge number per mode,
`Oh_eff_l(t)`, from mode `l`'s own instantaneous shear rate through the EXACT
(non-Taylor-expanded) Carreau-Yasuda law, then replaces BOTH the restoring
(`D1`) and damping (`D2`) coefficients via `julia/src/reid.jl`'s exact
finite-Oh relations -- not just `D2`, since `omega_l^2` only equals the
inviscid value in the small-Oh limit and this fluid's rest Oh is ~57.
"""

"""
    characteristic_shear_K(l) -> Float64

Geometric constant relating mode `l`'s characteristic (volume-RMS) shear rate
to its velocity amplitude: `gammadot_char_l = characteristic_shear_K(l) *
abs(Adot_l)`. From the exact single-mode potential-flow strain-rate field
(same "inviscid mode shape" already used by `st_extension.jl`'s `Gamma`);
see the derivation script for the exact-rational cross-check (`K_2^2=3`,
`K_3^2=4`, `K_4^2=9/2`, ... reproduced here in Float64).
"""
function characteristic_shear_K(l::Int)
    function legendre_P_dP(x::Float64)
        l == 0 && return 1.0, 0.0
        Pm1, P = 1.0, x
        for n in 1:l-1
            P, Pm1 = ((2n + 1) * x * P - n * Pm1) / (n + 1), P
        end
        Plm1 = l == 1 ? 1.0 : begin
            Qm1, Q = 1.0, x
            for n in 1:l-2
                Q, Qm1 = ((2n + 1) * x * Q - n * Qm1) / (n + 1), Q
            end
            Q
        end
        dP = l * (Plm1 - x * P) / (1 - x^2)
        P, dP
    end
    function g_l(x::Float64)
        X, Xp = legendre_P_dP(x)
        e_rr = l * (l - 1) * X
        e_thth = x * Xp - l^2 * X
        e_phph = l * X - x * Xp
        e_rth = -(l - 1) * sqrt(1 - x^2) * Xp
        2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2)
    end
    # Gauss-Legendre quadrature (40 points, exact for the polynomial integrand
    # up to the degrees arising here) rather than the exact-rational polynomial
    # integration used in the derivation script -- simpler in production, and
    # cross-checked against the exact values there (K_2^2=3, K_3^2=4, ...).
    nodes, weights = gauss_legendre_nodes(40, -1.0, 1.0)
    I_l = sum(w * g_l(x) for (x, w) in zip(nodes, weights))
    sqrt((3 / (2 * (2l - 1))) * I_l / l^2)
end

"""
Non-perturbative Carreau-Yasuda shear-thinning parameters.

- `lambda_c`, `a`, `eps_ST`: same meaning as `STParams` (`eps_ST=(1-n)/a`,
  `mu_eff/mu_0 = [1+(lambda_c*gammadot)^a]^(-eps_ST)`, EXACT -- no Taylor
  expansion, unlike `STParams`/`st_extension.jl`).
- `K`: `characteristic_shear_K(l)` for `l=2..M`, length M-1.
- `viscous`: `:lamb` (cheap, closed-form lambda_l/omega_l^2 at the exact
  Oh_eff_l) or `:reid` (exact finite-Oh via a cached/interpolated
  `ReidTable` per mode -- accurate at any Oh_eff_l, ~17000x cheaper per
  lookup than re-solving Reid's characteristic equation every step).
- `reid_cache`: `nothing` for `:lamb`; a `Vector{ReidTable}` (length M-1,
  built once via `build_reid_cache`) for `:reid`.
"""
struct STExactParams
    lambda_c   :: Float64
    a          :: Float64
    eps_ST     :: Float64
    K          :: Vector{Float64}
    viscous    :: Symbol
    reid_cache :: Union{Nothing,Vector{ReidTable}}
end

"""
    STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid, Oh_min=1e-4, Oh_max=nothing, n_table=150)

Build `STExactParams` for modes `l=2..M` at rest Ohnesorge number `Oh0`.
`Oh_max` defaults to `10*Oh0` (comfortable headroom above the rest state,
since shear-thinning only ever REDUCES Oh_eff below Oh0 for a fluid with
`eps_ST > 0`); pass explicitly if a wider margin is needed.
"""
function STExactParams(M::Int, Oh0::Float64, lambda_c::Float64, a::Float64,
    eps_ST::Float64; viscous::Symbol=:reid, Oh_min::Float64=1e-4,
    Oh_max::Union{Nothing,Float64}=nothing, n_table::Int=150)
    viscous in (:lamb, :reid) ||
        throw(ArgumentError("viscous must be :lamb or :reid, got $viscous"))
    K = [characteristic_shear_K(l) for l in 2:M]
    Oh_hi = Oh_max === nothing ? 10Oh0 : Oh_max
    reid_cache = if viscous === :reid
        build_reid_cache(M; Oh_min=Oh_min, Oh_max=Oh_hi, n=n_table)
    else
        nothing
    end
    STExactParams(lambda_c, a, eps_ST, K, viscous, reid_cache)
end

"""
    oh_eff_lambda_omega2(stx, cfg, k, Adot_char) -> (lambda, omega2)

Effective `Oh_eff_l` for mode index `k` (mode `l=k+1`) from the characteristic
shear rate `stx.K[k]*abs(Adot_char)`, through the exact Carreau-Yasuda law,
then `(lambda, omega2)` via Lamb's closed form or a cached Reid lookup.
"""
function oh_eff_lambda_omega2(stx::STExactParams, cfg::SimConstants, k::Int, Adot_char::Float64)
    l = k + 1
    gammadot_char = stx.K[k] * abs(Adot_char)
    Oh_eff = cfg.Oh * (1 + (stx.lambda_c * gammadot_char)^stx.a)^(-stx.eps_ST)
    if stx.viscous === :reid
        reid_lambda_omega2_fast(stx.reid_cache[k], Oh_eff)
    else
        Oh_eff * (l - 1) * (2l + 1), Float64(l * (l - 1) * (l + 2))
    end
end

"""
    build_residual_st_exact!(R, state, history, dt, cp, cfg, ob, stx)

Like `build_residual!`, but with the R2 block's D1 (restoring) and D2
(damping) replaced per-mode by `oh_eff_lambda_omega2` -- evaluated at the
LAGGED (`history[end]`) Ȧ, matching `st_extension.jl`'s own lagging
convention, so the Jacobian stays exact within a Newton step.
"""
function build_residual_st_exact!(R::AbstractVector, state::DropState,
    history::Vector{DropState}, dt::Float64, cp::Int, cfg::SimConstants,
    ob::OBParams, stx::STExactParams)
    M = cfg.M
    build_residual!(R, state, history, dt, cp, cfg, ob)

    order = length(history)
    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    Adot_prev_lag = history[end].Adot[2:end]

    ns = 2:M
    D1 = Vector{Float64}(undef, M - 1)
    D2 = Vector{Float64}(undef, M - 1)
    for k in 1:M-1
        lam, om2 = oh_eff_lambda_omega2(stx, cfg, k, Adot_prev_lag[k])
        D1[k] = om2
        D2[k] = 2lam
    end

    Adot = state.Adot[2:end]
    A = state.A[2:end]
    prev_A = hcat([h.A[2:end] for h in history]...)
    prev_Adot = hcat([h.Adot[2:end] for h in history]...)

    if ob.De1 > 0.0 && ob.beta_s < 1.0
        S_modes = state.S[2:end]
        effective_damp = D2 .* (ob.beta_s .* Adot .+ S_modes)
    else
        effective_damp = D2 .* Adot
    end
    B = state.B
    pressure_term = cp > 0 ? Float64.(ns) .* B[3:end] : zeros(M - 1)

    R[M:2M-2] .= c[end] .* Adot .+ sum(c[j] .* prev_Adot[:, j] for j in 1:order) .+
                 dt .* (D1 .* A .+ effective_damp .+ pressure_term)
end

"""
    build_jacobian_st_exact(state, history, dt, cp, cfg, ob, stx) → Matrix{Float64}

Analytical Jacobian for `build_residual_st_exact!`. D1, D2 are evaluated at
the same lagged Ȧ as the residual, so they're constant within a Newton step
and the Jacobian structure is unchanged from `build_jacobian`'s -- only the
D1, D2 diagonal VALUES differ (per-mode, from `oh_eff_lambda_omega2`).
"""
function build_jacobian_st_exact(state::DropState, history::Vector{DropState},
    dt::Float64, cp::Int, cfg::SimConstants, ob::OBParams, stx::STExactParams)
    J = build_jacobian(state, history, dt, cp, cfg, ob)

    M = cfg.M
    Nm = M - 1
    ak = bdf_coefficients(length(history), dt,
        length(history) == 2 ? history[end-1].dt : NaN)[end]
    Adot_prev_lag = history[end].Adot[2:end]
    damp_factor = (ob.De1 > 0.0 && ob.beta_s < 1.0) ? ob.beta_s : 1.0

    for k in 1:Nm
        lam, om2 = oh_eff_lambda_omega2(stx, cfg, k, Adot_prev_lag[k])
        J[Nm+k, k] = dt * om2
        J[Nm+k, Nm+k] = ak + dt * damp_factor * 2lam
    end
    J
end
