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
Minimum duration, in multiples of `cfg.dt_max`, that a contact segment must
span before `extract_kpis` accepts it as a bounce (value: 5). A contact patch
grows as `radius ~ sqrt(t - t_onset)` (Gabbard et al. 2025), so a segment
lasting only a step or two is still a near-singular point contact -- a
transient elastic "kiss" -- rather than an established macroscopic contact.
"""
const MIN_BOUNCE_STEPS = 5

"""
Largest coefficient of restitution `extract_kpis` will accept from a contact
segment before rejecting it and moving on to the next candidate (value: 1.02).
A passive, dissipative impact cannot gain energy net of the gravitational term
already carried in the `cor` formula, so `cor > 1` beyond quadrature and
discretization noise marks a numerical breakdown rather than a fast bounce --
which duration alone cannot detect.
"""
const MAX_ALLOWED_COR = 1.02

# Provenance for the two constants above, kept out of the docstrings so it does
# not render in the manual.
#
# MIN_BOUNCE_STEPS: duration alone does not separate a spurious micro-kiss from
# a real but fast bounce -- both appeared at 3-7 steps across the 20-sample
# validation sweep (julia/derivations/carreau_yasuda_multimode_derivation.jl,
# Section 4). 5 sits above the micro-kiss's consistently observed 3-step
# duration and below the shortest genuine bounce observed, at 6 steps.
# MAX_ALLOWED_COR is the primary filter for segments that clear this floor but
# are still numerically invalid.
#
# MAX_ALLOWED_COR: the breakdown this rejects was traced to a shape grown large
# enough, after an earlier violent contact cycle, to leave the small-deformation
# regime the linearized model assumes; it shows up as a discontinuous velocity
# jump of several times the impact speed at re-contact. A real short bounce and
# a broken segment can share a duration, but only the broken one violates energy
# conservation.

"""
    _tentative_cor(states, cfg, fc, lc) -> Float64

Coefficient of restitution for the candidate segment `[fc, lc]`, or `NaN` if
it cannot be evaluated (segment touches either end of the trajectory). Used
by `extract_kpis` to screen candidate segments against `MAX_ALLOWED_COR`
before accepting one.
"""
function _tentative_cor(states::Vector{DropState}, cfg::SimConstants, fc::Int, lc::Int)
    (fc <= 1 || lc >= length(states)) && return NaN
    v_in  = states[fc - 1].v
    z_in  = states[fc - 1].z
    v_out = states[lc + 1].v
    z_out = states[lc + 1].z
    E_in  = 0.5 * v_in^2
    E_in ≈ 0.0 && return NaN
    E_out = 0.5 * v_out^2 + cfg.Bo * (z_out - z_in)
    sqrt(abs(E_out / E_in))
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
    extract_kpis(times, states, cfg) → SweepKPIs

Extract key performance indicators from a completed `solve_drop!` run.
All KPIs are NaN / 0.0 when no contact occurs in the simulation.

`contact_time` and `cor` describe the FIRST bounce only (first contact onset
to the end of that same, contiguous contact segment) -- matching how both
are measured experimentally (a single impact-and-rebound event) and how
Gabbard et al. (2025) define them. A drop that bounces, falls back under
gravity, and re-contacts within the simulated window (routine at low We with
a shear-thinning fluid, since a small rebound velocity plus gravity easily
produces a second, third, ... impact before `t_end`) must NOT have its
first-bounce CoR computed against the LAST such contact instead of the
first -- that would silently mix in the free-fall speed gained between
bounces, which is real physics but not what CoR means for a single impact.

A `cp == 0` gap between two contact frames is bridged (treated as still part
of the same bounce) when its duration is no longer than `cfg.dt_max`, the
solver's own base adaptive timestep. That is the same physical scale that
justifies the solver's `cp_prev-1:cp_prev+1` contact search: `cp` cannot jump
discontinuously within one timestep of continuous evolution, and a genuine
separation-and-reimpact requires the whole drop to retreat and fall back under
gravity, which cannot complete within one `dt_max`. A gap longer than `dt_max`
DOES end the segment — that is what separates the true first bounce from the
later, gravity-driven re-impacts of a multi-bounce trajectory.

A bridged segment is accepted as "the first bounce" only if it lasts at least
`MIN_BOUNCE_STEPS * cfg.dt_max`. At a near-point contact patch (`cp` at its
smallest values) the coupled Carreau-Yasuda model can predict a real but
extremely brief elastic "kiss": local shear at a near-singular point contact
is enormous, crashing the local effective viscosity and momentarily removing
essentially all damping, so a tiny dimple springs back within a handful of
timesteps, before the growing patch and mounting contact force produce the
real, sustained bounce a moment later. That kiss is dynamically real rather
than numerical noise, but it is orders of magnitude shorter than any
resolvable experimental measurement (a shape-camera frame interval, or even
the capillary oscillation period) and is not what "contact time" means either
experimentally or in Gabbard et al. (2025). Shorter segments are skipped in
favor of the next candidate.

A segment that passes the duration floor is ADDITIONALLY rejected if its
tentative `cor` would exceed `MAX_ALLOWED_COR`. Duration alone cannot
distinguish a real-but-fast bounce from a segment beginning at a numerical
breakdown — a shape that has grown too large, after an earlier violent contact
cycle, for this linearized model's small-deformation assumption, which shows up
as a discontinuous velocity jump of several times the impact speed at
re-contact. Only the energy-conservation check tells the two apart.

`max_radius`/`max_A2` still reflect the WHOLE trajectory (all bounces).
"""
function extract_kpis(times::Vector{Float64},
                      states::Vector{DropState},
                      cfg::SimConstants)
    contact_idx = findall(s -> s.cp > 0, states)

    isempty(contact_idx) && return SweepKPIs(NaN, NaN, 0.0, 0.0)

    first_c = nothing
    last_c  = nothing
    ci = 1
    while ci <= length(contact_idx)
        fc = contact_idx[ci]
        lc = fc
        i = fc
        while i < length(states)
            if states[i+1].cp > 0
                lc = i + 1
                i += 1
            else
                j = i + 1
                while j <= length(states) && states[j].cp == 0
                    j += 1
                end
                if j <= length(states) && (times[j] - times[i]) <= cfg.dt_max
                    lc = j
                    i = j
                else
                    break
                end
            end
        end

        if times[lc] - times[fc] >= MIN_BOUNCE_STEPS * cfg.dt_max
            tentative = _tentative_cor(states, cfg, fc, lc)
            if isnan(tentative) || tentative <= MAX_ALLOWED_COR
                first_c, last_c = fc, lc
                break
            end
        end
        while ci <= length(contact_idx) && contact_idx[ci] <= lc
            ci += 1
        end
    end

    first_c === nothing && return SweepKPIs(NaN, NaN, 0.0, 0.0)

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
