"""
    build_residual_st!(R, state, history, dt, cp, cfg, ob, st)

Fill residual R (length 3M+1) with Carreau-Yasuda shear-thinning correction to R2.

Calls build_residual! for the Newtonian part, then corrects the R2 block:
  R2[k] -= dt * D2[k] * Ȧ_curr[k] * ε_ST * shear_pow_lag

where shear_pow_lag = Σ_k Gamma_eff[k] * |Ȧ_prev[k]|^a is lagged (history[end]).
`st.a` is the Carreau-Yasuda shape exponent (a=2 recovers the original,
always-quadratic Carreau correction exactly — see
julia/derivations/carreau_yasuda_derivation.jl §7, ASSERTION 17 for the
exact reduction).
When st.eps_ST == 0 the result is identical to build_residual!.
"""
function build_residual_st!(R::AbstractVector, state::DropState,
                             history::Vector{DropState}, dt::Float64,
                             cp::Int, cfg::SimConstants, ob::OBParams,
                             st::STParams)
    M = cfg.M
    build_residual!(R, state, history, dt, cp, cfg, ob)

    st.eps_ST == 0.0 && return

    ns        = collect(Float64, 2:M)
    D2        = @. 2cfg.Oh * (ns - 1) * (2ns + 1)
    sigma0    = @. sqrt(ns * (ns - 1) * (ns + 2))
    Gamma_eff = st.Gamma .* (st.lambda_c .* sigma0) .^ st.a

    Adot_prev     = history[end].Adot[2:end]
    shear_pow_lag = sum(Gamma_eff[k] * abs(Adot_prev[k])^st.a for k in 1:M-1)

    Adot_curr = state.Adot[2:end]
    R[M:2M-2] .-= dt .* D2 .* Adot_curr .* (st.eps_ST * shear_pow_lag)
end

"""
    build_jacobian_st(state, history, dt, cp, cfg, ob, st) → Matrix{Float64}

Analytical Jacobian for the Carreau-Yasuda system. Identical to build_jacobian
except the ∂R2/∂Ȧ diagonal is multiplied by (1 - ε_ST * shear_pow_lag).
shear_pow_lag is lagged (uses history[end].Adot), so the Jacobian is constant
across Newton iterations within a step — unaffected by `a` other than through
this lagged scalar, since the residual correction is still linear in the
current Ȧ.
"""
function build_jacobian_st(state::DropState, history::Vector{DropState},
                            dt::Float64, cp::Int, cfg::SimConstants,
                            ob::OBParams, st::STParams)
    J = build_jacobian(state, history, dt, cp, cfg, ob)

    st.eps_ST == 0.0 && return J

    M  = cfg.M
    Nm = M - 1
    ns = collect(Float64, 2:M)
    D2 = @. 2cfg.Oh * (ns - 1) * (2ns + 1)
    sigma0    = @. sqrt(ns * (ns - 1) * (ns + 2))
    Gamma_eff = st.Gamma .* (st.lambda_c .* sigma0) .^ st.a

    Adot_prev     = history[end].Adot[2:end]
    shear_pow_lag = sum(Gamma_eff[k] * abs(Adot_prev[k])^st.a for k in 1:M-1)

    for k in 1:Nm
        J[Nm+k, Nm+k] -= dt * D2[k] * st.eps_ST * shear_pow_lag
    end
    return J
end
