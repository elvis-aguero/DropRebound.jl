"""
Key performance indicators from a single drop impact simulation.

Fields:
  contact_time  – duration of the contact phase in capillary times (NaN if no contact
                  occurred, or if the drop was still in contact at t_end)
  cor           – coefficient of restitution √(|E_out/E_in|) (NaN if invalid)
  max_radius    – maximum dimensionless contact radius during contact
  max_A2        – maximum |A₂| deformation amplitude during contact
"""
struct SweepKPIs
    contact_time :: Float64
    cor          :: Float64
    max_radius   :: Float64
    max_A2       :: Float64
end

"""
    compute_contact_radius(state, cfg) → Float64

Dimensionless radius of the contact patch: sin(θ_cp) × (1 + shape deformation at θ_cp).
Returns 0.0 when cp = 0.
"""
function compute_contact_radius(state::DropState, cfg::SimConstants)
    state.cp == 0 && return 0.0
    θ = cfg.theta_vec[state.cp]
    P = collect_Pl(cfg.M, [cos(θ)])
    deform = sum(P[1, n+1] * state.A[n] for n in 2:cfg.M)
    return sin(θ) * (1.0 + deform)
end

"""
    drop_profile(state, cfg; n_theta=200) → (xs, zs)

Reconstruct the meridional cross-section of the drop surface as two vectors of
Cartesian coordinates (x, z), sampled at `n_theta` angles from θ=0 (top) to
θ=π (south pole). The profile is the right half of the axisymmetric outline;
reflect x → -x for the full silhouette.

Coordinates are dimensionless (scaled by R): x = r sinθ, z = z_com + r cosθ.
"""
function drop_profile(state::DropState, cfg::SimConstants; n_theta::Int=200)
    θ_vals = range(0.0, π; length=n_theta)
    cosθ   = cos.(θ_vals)
    P      = collect_Pl(cfg.M, cosθ)
    r      = ones(n_theta)
    for n in 2:cfg.M
        r .+= state.A[n] .* P[:, n+1]
    end
    xs = r .* sin.(θ_vals)
    zs = state.z .+ r .* cosθ
    return xs, zs
end

"""
    com_energy(state, cfg) → Float64

Centre-of-mass mechanical energy per unit mass, ½v² + Bo·z. Conserved in free flight
and strictly decreasing across a bounce (equivalently COR ≤ 1); used as the robust
energy-injection acceptance gate. (The full internal energy budget, which requires the
liquid–gas cap-area accounting, is not evaluated here — see the coefficient of
restitution in `extract_kpis`.)
"""
com_energy(state::DropState, cfg::SimConstants) = 0.5 * state.v^2 + cfg.Bo * state.z

"""
    extract_kpis(times, states, cfg) → SweepKPIs

Extract key performance indicators from a completed `solve_drop!` run.
All KPIs are NaN / 0.0 when no contact occurs in the simulation.
"""
function extract_kpis(times::Vector{Float64},
                      states::Vector{DropState},
                      cfg::SimConstants)
    contact_idx = findall(s -> s.cp > 0, states)

    isempty(contact_idx) && return SweepKPIs(NaN, NaN, 0.0, 0.0)

    first_c = contact_idx[1]
    last_c  = contact_idx[end]

    contact_time = times[last_c] - times[first_c]

    cor = if last_c < length(states) && first_c > 1
        v_in  = states[first_c - 1].v
        z_in  = states[first_c - 1].z
        v_out = states[last_c + 1].v
        z_out = states[last_c + 1].z
        E_in  = 0.5 * v_in^2
        E_out = 0.5 * v_out^2 + cfg.Bo * (z_out - z_in)
        (E_in ≈ 0.0) ? NaN : sqrt(abs(E_out / E_in))
    else
        NaN
    end

    max_radius = maximum(compute_contact_radius(s, cfg) for s in states if s.cp > 0)
    max_A2     = maximum(abs(s.A[2]) for s in states if s.cp > 0)

    return SweepKPIs(contact_time, cor, max_radius, max_A2)
end
