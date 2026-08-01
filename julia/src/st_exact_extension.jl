"""
Non-perturbative, multi-mode-coupled Carreau-Yasuda shear-thinning extension.
Derivation: `julia/derivations/carreau_yasuda_nonperturbative_derivation.jl`
(single-mode groundwork) and `julia/derivations/carreau_yasuda_multimode_derivation.jl`
(the coupling this file actually implements).

Two simpler routes are wrong here. `st_extension.jl` Taylor-expands
Carreau-Yasuda to first order in `eps_ST`, which is invalid once `eps_ST` is
not small -- the validation fluid this was built for has `eps_ST ~ 0.9996`.
Computing each mode's effective Oh from that mode's OWN velocity alone is also
wrong. This file instead computes an effective Ohnesorge number per
mode from the TRUE, fully-superposed local shear-rate field of every
currently active mode -- shear rate is a kinematic field quantity, so a
mode excited by contact still contributes to the shear (and hence thinning)
felt by every other active mode at the same point in the drop, not just its
own. `Oh_eff_l(t)` then replaces BOTH the restoring (`D1`) and damping
(`D2`) coefficients via `julia/src/reid.jl`'s exact finite-Oh relations.
"""

function _legendre_P_dP(l::Int, x::Float64)
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

"""
    _strain_basis(l, x) -> (e_rr, e_thth, e_phph, e_rth)

Strain-rate tensor components for mode `l`'s own potential-flow field, per
unit `Adot_l` and per unit `r^(l-2)` (the r-dependence is a single power per
mode, applied separately when building the full multi-mode field). See
the "Carreau-Yasuda: Multi-Mode" derivation, Section 1.
"""
function _strain_basis(l::Int, x::Float64)
    X, Xp = _legendre_P_dP(l, x)
    e_rr = (l - 1) * X
    e_thth = (x * Xp - l^2 * X) / l
    e_phph = (l * X - x * Xp) / l
    e_rth = -(l - 1) / l * sqrt(1 - x^2) * Xp
    e_rr, e_thth, e_phph, e_rth
end

"""
    characteristic_shear_K(l) -> Float64

Geometric constant relating mode `l`'s OWN characteristic (volume-RMS) shear
rate to its velocity amplitude in ISOLATION (`gammadot_char_l =
characteristic_shear_K(l) * abs(Adot_l)`), from mode `l`'s single-mode
potential-flow field. This is a diagnostic: it is the value the multi-mode
coupled `oh_eff_all_coupled` must reduce to when only one mode is active.
`build_residual_st_exact!`/`build_jacobian_st_exact` never call it — they
always use the full coupling.
"""
function characteristic_shear_K(l::Int)
    g_l(x) = begin
        e_rr, e_thth, e_phph, e_rth = _strain_basis(l, x)
        2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2)
    end
    nodes, weights = gauss_legendre_nodes(40, -1.0, 1.0)
    I_l = sum(w * g_l(x) for (x, w) in zip(nodes, weights))
    # _strain_basis is already normalized per unit Adot_l (needed for the
    # multi-mode sum in oh_eff_all_coupled), so no further division by l^2
    # is needed here (an earlier version of this function used a per-unit-
    # Phi_l basis and divided by l^2 at this step; that basis and this one
    # are NOT interchangeable -- caught by this function's own test
    # (julia/test/test_carreau_yasuda_exact.jl) suddenly returning K_l^2/l^2
    # instead of K_l^2 after unifying the two bases).
    sqrt((3 / (2 * (2l - 1))) * I_l)
end

"""
Non-perturbative, multi-mode-coupled Carreau-Yasuda shear-thinning parameters.

- `lambda_c`, `a`, `eps_ST`: same meaning as `STParams` (`eps_ST=(1-n)/a`,
  `mu_eff/mu_0 = eta_inf_ratio + (1-eta_inf_ratio)*[1+(lambda_c*gammadot)^a]^(-eps_ST)`,
  EXACT -- no Taylor expansion).
- `eta_inf_ratio`: `eta_inf/eta_0`, the fluid's infinite-shear viscosity
  floor relative to its zero-shear viscosity. Defaults to `0.0`, the plain
  Carreau-Yasuda law, in which `mu_eff/mu_0 -> 0` as `gammadot -> infinity`.
  For a REAL Cross-model fluid the ratio is not exactly zero (`eta_inf` is
  always finite, if small), and leaving it at zero lets `Oh_eff` fall below
  `Oh0*eta_inf/eta_0`, the lowest Ohnesorge number the fluid can physically
  reach at any shear rate. On the validation fluid at `We = 0.7649` that floor
  is ~0.0254, and the unfloored model reached `Oh_eff` of 0.0074 and 0.0044.
- `r_nodes,r_wts,x_nodes,x_wts`: Gauss-Legendre quadrature for the volume
  integral (`r in [0,1]`, `x=cos(theta) in [-1,1]`).
- `e_rr,e_thth,e_phph,e_rth`: precomputed, STATE-INDEPENDENT strain-rate
  basis values at every quadrature node, for modes `l=2..M`
  (`e_rr[k][i,j]` = mode `l=k+1`'s `e_rr` at `(r_nodes[i], x_nodes[j])`,
  including the `r^(l-2)` factor). Built once; reused every residual/
  Jacobian evaluation, since only the CURRENT `Adot` vector (not this
  geometry) changes step to step.
- `viscous`: `:lamb` (cheap closed form at the exact `Oh_eff_l`) or `:reid`
  (cached/interpolated `ReidTable` lookup, accurate at any `Oh_eff_l`).
- `reid_cache`: `nothing` for `:lamb`; `Vector{ReidTable}` (length M-1) for `:reid`.
"""
struct STExactParams
    lambda_c :: Float64
    a        :: Float64
    eps_ST   :: Float64
    eta_inf_ratio :: Float64
    r_nodes  :: Vector{Float64}
    r_wts    :: Vector{Float64}
    x_nodes  :: Vector{Float64}
    x_wts    :: Vector{Float64}
    e_rr     :: Vector{Matrix{Float64}}
    e_thth   :: Vector{Matrix{Float64}}
    e_phph   :: Vector{Matrix{Float64}}
    e_rth    :: Vector{Matrix{Float64}}
    viscous    :: Symbol
    reid_cache :: Union{Nothing,Vector{ReidTable}}
end

"""
    STExactParams(M, Oh0, lambda_c, a, eps_ST; viscous=:reid, eta_inf_ratio=0.0,
                  Oh_min=1e-4, Oh_max=nothing, n_table=150, n_r=20, n_x=30)

Build `STExactParams` for modes `l=2..M` at rest Ohnesorge number `Oh0`.
`Oh_max` defaults to `10*Oh0` (shear-thinning only ever reduces `Oh_eff`
below `Oh0` for `eps_ST > 0`); pass explicitly if a wider margin is needed.
`n_r`, `n_x` control the volume-quadrature resolution for the coupling
integral -- 20x30 (600 points) was sufficient in testing; increase if
`Oh_eff` estimates look under-resolved for a very large `M`.
"""
function STExactParams(M::Int, Oh0::Float64, lambda_c::Float64, a::Float64,
    eps_ST::Float64; viscous::Symbol=:reid, eta_inf_ratio::Float64=0.0,
    Oh_min::Float64=1e-4,
    Oh_max::Union{Nothing,Float64}=nothing, n_table::Int=150,
    n_r::Int=20, n_x::Int=30)
    viscous in (:lamb, :reid) ||
        throw(ArgumentError("viscous must be :lamb or :reid, got $viscous"))
    (0.0 <= eta_inf_ratio < 1.0) ||
        throw(ArgumentError("eta_inf_ratio must be in [0,1), got $eta_inf_ratio"))

    r_nodes, r_wts = gauss_legendre_nodes(n_r, 0.0, 1.0)
    x_nodes, x_wts = gauss_legendre_nodes(n_x, -1.0, 1.0)
    nl = M - 1
    e_rr = Vector{Matrix{Float64}}(undef, nl)
    e_thth = Vector{Matrix{Float64}}(undef, nl)
    e_phph = Vector{Matrix{Float64}}(undef, nl)
    e_rth = Vector{Matrix{Float64}}(undef, nl)
    for (idx, l) in enumerate(2:M)
        Mrr = zeros(n_r, n_x)
        Mthth = zeros(n_r, n_x)
        Mphph = zeros(n_r, n_x)
        Mrth = zeros(n_r, n_x)
        for (i, r) in enumerate(r_nodes), (j, x) in enumerate(x_nodes)
            b = _strain_basis(l, x)
            rp = r^(l - 2)
            Mrr[i, j] = b[1] * rp
            Mthth[i, j] = b[2] * rp
            Mphph[i, j] = b[3] * rp
            Mrth[i, j] = b[4] * rp
        end
        e_rr[idx], e_thth[idx], e_phph[idx], e_rth[idx] = Mrr, Mthth, Mphph, Mrth
    end

    Oh_hi = Oh_max === nothing ? 10Oh0 : Oh_max
    reid_cache = viscous === :reid ? build_reid_cache(M; Oh_min=Oh_min, Oh_max=Oh_hi, n=n_table) : nothing

    STExactParams(lambda_c, a, eps_ST, eta_inf_ratio, r_nodes, r_wts, x_nodes, x_wts,
        e_rr, e_thth, e_phph, e_rth, viscous, reid_cache)
end

"""
How much larger than every trusted (low- and mid-index) mode a top-of-spectrum
mode's velocity must be before it is treated as truncation ringing and dropped
from the coupled shear field (value: 20). Ringing is the un-representable
leftover of a non-smooth boundary condition — the contact-onset kink — in a
finite Legendre expansion; counting it as real shear would thin the fluid on a
signal that is a discretization artifact, not motion.
"""
const OUTLIER_FACTOR = 20.0

# Floor on the trusted-mode reference amplitude, so a fully quiescent trusted
# set (~1e-15 of real floating-point noise, from steps that have not yet moved)
# does not make every nonzero candidate look infinitely large by comparison.
const RINGING_NOISE_FLOOR = 1e-9

# _ringing_outlier_mask: `true` at index k means mode k (l=k+1) is excluded from
# the coupled local shear field built by oh_eff_all_coupled.
#
# Two signals are combined, because neither alone works:
#
# - INDEX restricts which modes can even be CANDIDATES for exclusion, to roughly
#   the top `candidate_fraction` (minimum 1) nearest the truncation boundary --
#   the only place the un-representable leftover can concentrate. Amplitude ratio
#   ALONE cannot make this restriction: a candidate's amplitude must be compared
#   against something, and comparing against ALL other modes (including other
#   high, possibly-also-large ones) breaks the moment one of those "other" modes
#   is itself just floating-point noise (~1e-15, not exactly 0.0) sitting next to
#   a genuinely large low-index mode. An amplitude-only version of this check
#   flagged real mid-mode physics (l=8) as ringing purely because a nearby low
#   mode happened to be at the noise floor instead of exactly zero.
#
# - AMPLITUDE decides whether a candidate is ACTUALLY ringing: excluded only if
#   its magnitude exceeds `factor` times the largest magnitude among the trusted
#   modes, floored at RINGING_NOISE_FLOOR. This distinguishes observed ringing (a
#   single top mode at ~1e-1 while every trusted mode sits at ~1e-15) from
#   genuine higher-mode physics (a low-We run where mode M reached ~0.0163,
#   comparable to mode 2's ~0.0159 -- neither is noise, so mode M is kept).
#
# The two cases this was tuned against are in
# julia/derivations/carreau_yasuda_multimode_derivation.jl, Sections 4 and 4b.
function _ringing_outlier_mask(Adot_vec::AbstractVector{Float64};
    candidate_fraction::Float64=0.1, factor::Float64=OUTLIER_FACTOR)
    nl = length(Adot_vec)
    M = nl + 1
    n_candidates = max(1, ceil(Int, candidate_fraction * M))
    n_trusted = nl - n_candidates
    absvals = abs.(Adot_vec)
    trusted_max = n_trusted <= 0 ? 0.0 : maximum(@view absvals[1:n_trusted])
    reference = max(trusted_max, RINGING_NOISE_FLOOR)

    mask = falses(nl)
    for k in (n_trusted + 1):nl
        absvals[k] > factor * reference && (mask[k] = true)
    end
    mask
end

"""
    oh_eff_all_coupled(stx, Oh0, Adot_vec) -> Vector{Float64}

`Oh_eff_l` for every mode `l=2..M` (length M-1, indexed like `Adot_vec`),
given the CURRENT full mode-velocity vector. Each is a dissipation-weighted
average of the true local viscosity ratio `mu_eff(S_total)/mu_0` over the
whole drop, where `S_total` is the shear-rate invariant of the FULLY
SUPERPOSED strain field from every active mode in `Adot_vec` that is NOT
flagged by `_ringing_outlier_mask` (not mode `l` alone), weighted by mode
`l`'s own sensitivity `dS^2/dAdot_l` to that field. Every mode (including any
masked-out one) still receives an `Oh_eff` value -- the mask only controls
what counts as "real" shear when building the field itself.
See the "Carreau-Yasuda: Multi-Mode" derivation, Section 3.
"""
function oh_eff_all_coupled(stx::STExactParams, Oh0::Float64, Adot_vec::AbstractVector{Float64})
    nl = length(Adot_vec)
    n_r, n_x = length(stx.r_nodes), length(stx.x_nodes)
    exclude = _ringing_outlier_mask(Adot_vec)

    Trr = zeros(n_r, n_x)
    Tthth = zeros(n_r, n_x)
    Tphph = zeros(n_r, n_x)
    Trth = zeros(n_r, n_x)
    for k in 1:nl
        (Adot_vec[k] == 0.0 || exclude[k]) && continue
        @. Trr += Adot_vec[k] * stx.e_rr[k]
        @. Tthth += Adot_vec[k] * stx.e_thth[k]
        @. Tphph += Adot_vec[k] * stx.e_phph[k]
        @. Trth += Adot_vec[k] * stx.e_rth[k]
    end

    muratio = similar(Trr)
    for i in 1:n_r, j in 1:n_x
        S = sqrt(2 * (Trr[i, j]^2 + Tthth[i, j]^2 + Tphph[i, j]^2 + 2 * Trth[i, j]^2))
        muratio[i, j] = stx.eta_inf_ratio +
                        (1 - stx.eta_inf_ratio) * (1 + (stx.lambda_c * S)^stx.a)^(-stx.eps_ST)
    end

    Oh_eff = Vector{Float64}(undef, nl)
    for k in 1:nl
        num = 0.0
        den = 0.0
        for i in 1:n_r, j in 1:n_x
            wl = 2 * (Trr[i, j] * stx.e_rr[k][i, j] + Tthth[i, j] * stx.e_thth[k][i, j] +
                      Tphph[i, j] * stx.e_phph[k][i, j] + 2 * Trth[i, j] * stx.e_rth[k][i, j])
            wl == 0.0 && continue
            weight = abs(wl) * stx.r_wts[i] * stx.x_wts[j] * stx.r_nodes[i]^2
            num += weight * muratio[i, j]
            den += weight
        end
        Oh_eff[k] = den == 0.0 ? Oh0 : Oh0 * num / den
    end
    Oh_eff
end

"""
    lambda_omega2_from_oh_eff(stx, Oh_eff) -> (lambda::Vector, omega2::Vector)

Per-mode `(lambda_l, omega_l^2)` from `Oh_eff_l` (one entry per mode,
`Oh_eff[k]` for mode `l=k+1`), via Lamb's closed form or a cached Reid lookup.
"""
function lambda_omega2_from_oh_eff(stx::STExactParams, Oh_eff::AbstractVector{Float64})
    nl = length(Oh_eff)
    lambda = Vector{Float64}(undef, nl)
    omega2 = Vector{Float64}(undef, nl)
    for k in 1:nl
        l = k + 1
        if stx.viscous === :reid
            lambda[k], omega2[k] = reid_lambda_omega2_fast(stx.reid_cache[k], Oh_eff[k])
        else
            lambda[k] = Oh_eff[k] * (l - 1) * (2l + 1)
            omega2[k] = Float64(l * (l - 1) * (l + 2))
        end
    end
    lambda, omega2
end

"""
    build_residual_st_exact!(R, state, history, dt, cp, cfg, ob, stx)

Like `build_residual!`, but with the R2 block's D1 (restoring) and D2
(damping) replaced per-mode by `oh_eff_all_coupled`/`lambda_omega2_from_oh_eff`
-- evaluated at the LAGGED (`history[end]`), FULL mode-velocity vector,
matching `st_extension.jl`'s own lagging convention, so the Jacobian stays
diagonal and exact within a Newton step even though `Oh_eff_l` for every
mode now depends on every other mode's (lagged) velocity.
"""
function build_residual_st_exact!(R::AbstractVector, state::DropState,
    history::Vector{DropState}, dt::Float64, cp::Int, cfg::SimConstants,
    ob::OBParams, stx::STExactParams)
    M = cfg.M
    build_residual!(R, state, history, dt, cp, cfg, ob)

    order = length(history)
    c = bdf_coefficients(order, dt, order == 2 ? history[end-1].dt : NaN)
    Adot_prev_lag = history[end].Adot[2:end]

    Oh_eff = oh_eff_all_coupled(stx, cfg.Oh, Adot_prev_lag)
    lambda, omega2 = lambda_omega2_from_oh_eff(stx, Oh_eff)
    D1 = omega2
    D2 = 2 .* lambda

    ns = 2:M
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

Analytical Jacobian for `build_residual_st_exact!`. `Oh_eff_l` (hence D1, D2)
is evaluated at the same lagged, full Adot vector as the residual, so it is
constant within a Newton step and the Jacobian structure is unchanged from
`build_jacobian`'s (still diagonal in the D1/D2 block) -- only the diagonal
VALUES differ, from the coupled `Oh_eff`.
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

    Oh_eff = oh_eff_all_coupled(stx, cfg.Oh, Adot_prev_lag)
    lambda, omega2 = lambda_omega2_from_oh_eff(stx, Oh_eff)

    for k in 1:Nm
        J[Nm+k, k] = dt * omega2[k]
        J[Nm+k, Nm+k] = ak + dt * damp_factor * 2lambda[k]
    end
    J
end
