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
Minimum number of `dt_max` steps a bridged contact segment must span to be
accepted as a real bounce in `extract_kpis`, rather than a transient
elastic "kiss" at a near-point contact patch. The patch-growth law already
assumed at contact onset elsewhere in this repo (`radius ~ sqrt(t-t_onset)`,
Gabbard et al. 2025) needs several steps before a patch is a meaningfully
established, macroscopic contact rather than a single near-singular point.
Duration ALONE cannot fully distinguish a spurious micro-kiss from a real
but fast bounce (both were observed at 3-7 steps across the 20-sample
validation sweep, see `julia/derivations/carreau_yasuda_multimode_derivation.jl`
Section 4) -- 5 was chosen as comfortably above the spurious micro-kiss's
consistently-observed 3-step duration while remaining below the shortest
genuine bounce observed (6 steps); `MAX_ALLOWED_COR` below is the primary
filter for segments that pass this duration floor but are still numerically
invalid.
"""
const MIN_BOUNCE_STEPS = 5

"""
Largest coefficient of restitution `extract_kpis` will accept from a
candidate contact segment before rejecting it in favor of the next
candidate. A passive, dissipative fluid impact can never gain kinetic
energy net of the gravitational PE change already accounted for in the
`cor` formula, so `cor > 1` (beyond a small allowance for quadrature/
discretization noise) is not a physically valid bounce -- it is a signature
of a genuine numerical breakdown (traced live: a shape that has grown large
enough, after an earlier violent contact cycle, to leave the small-
deformation regime this linearized model assumes shows a discontinuous
velocity jump of several times the impact speed at re-contact). Rejecting
such segments and continuing to the next candidate is what actually
distinguishes them from a real, fast bounce -- duration alone cannot (see
`MIN_BOUNCE_STEPS`): a real bounce that happens to be short and a numerically
broken one can have the same duration, but only the broken one violates
energy conservation.
"""
const MAX_ALLOWED_COR = 1.02

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
solver's own base adaptive timestep. This is the same physical scale already
used to justify contact continuity elsewhere in this repo (the cp-search
window is `cp_prev-1:cp_prev+1` because cp cannot jump discontinuously
within one timestep of continuous evolution): a genuine separation-and-
reimpact requires the whole drop to retreat and fall back under gravity,
which cannot complete within one `dt_max`. A gap longer than `dt_max` DOES
end the segment, since that duration is long enough to be a real separation
(this is exactly what distinguishes the true first bounce from later,
gravity-driven re-impacts in a multi-bounce trajectory).

A bridged segment is only accepted as "the first bounce" if its duration is
at least `MIN_BOUNCE_STEPS * cfg.dt_max` (see that constant's docstring):
right at a near-point contact patch (`cp` at its smallest values), the
coupled Carreau-Yasuda model can predict a real but extremely brief elastic
"kiss" -- local shear right at a near-singular point contact is enormous,
crashing the local effective viscosity and momentarily removing essentially
all damping, so a tiny dimple springs back within a handful of timesteps,
before the growing patch and mounting contact force produce the real,
sustained bounce a short time later. That kiss is dynamically real (not
numerical noise -- see `julia/derivations/carreau_yasuda_multimode_derivation.jl`
Section 4), but it is many orders of magnitude shorter than any resolvable
experimental measurement (a shape-camera frame interval, or even the
capillary oscillation period itself) and is not what "contact time" means
either experimentally or in Gabbard et al. (2025). Segments shorter than
this threshold are skipped in favor of the next candidate segment.

A bridged segment that passes the duration floor is ADDITIONALLY rejected
if its tentative `cor` would exceed `MAX_ALLOWED_COR` (see that constant's
docstring): duration alone cannot distinguish a real-but-fast bounce from a
segment starting at a genuine numerical breakdown (a shape that has grown
too large, after an earlier violent contact cycle, for this linearized
model's small-deformation assumption -- observed live as a discontinuous
velocity jump of several times the impact speed right at re-contact). Both
failure modes were observed at similar (6-13 step) durations in the
20-sample validation sweep, so only the energy-conservation check reliably
tells them apart.

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
