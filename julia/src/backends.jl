# One entry point for every solver in the package.
#
# There are two formulations and three ways of closing the contact, and until now each was
# reached through a different function with a different signature, a different parameter object
# and a different result type. That is fine when you are working on one of them and hopeless
# when you want to compare them, which is the only reason to keep more than one.
#
# `Backend` names a choice on each axis and `run_impact` returns the SAME NamedTuple whichever
# is chosen, including the surface history in a common form so that a plot or an animation does
# not need to know which solver produced it.
#
#   formulation  :variational     interior displacement amplitudes, Euler-Lagrange
#                :nonvariational  surface amplitudes with Reid's per-mode coefficients
#
#   contact      :lcp             complementarity; the contact set is an output
#                :active_set      primal active-set iteration on the same two inequalities,
#                                 restricted to a pole-anchored interval
#                :tangency        candidate contact counts ranked by tangency error
#
#   forcing      :legendre        film pressure as a degree-M Legendre field  (variational)
#                :nodal           vertical loads at the collocation nodes     (variational)
#
# Not every combination exists. `:tangency` is the nonvariational solver's closure and
# `:active_set` is the variational one's; `:forcing` is meaningless without an interior. The
# constructor rejects the combinations that do not, rather than silently substituting.

"""
    Backend(; formulation = :variational, contact = :lcp, forcing = :legendre)

A choice of solver. See [`run_impact`](@ref).

    Backend()                                             # variational, complementarity, spectral pressure
    Backend(forcing = :nodal)                             # ... with conjugate nodal loads
    Backend(contact = :active_set)                        # variational, primal search
    Backend(formulation = :nonvariational, contact = :tangency)
"""
struct Backend
    formulation::Symbol
    contact::Symbol
    forcing::Symbol
end

function Backend(; formulation::Symbol = :variational, contact::Symbol = :lcp,
                 forcing::Symbol = :legendre)
    formulation in (:variational, :nonvariational) ||
        error("formulation must be :variational or :nonvariational, got $formulation")
    if formulation === :variational
        contact in (:lcp, :active_set) ||
            error("the variational formulation offers :lcp or :active_set, got $contact")
        forcing in (:legendre, :nodal) ||
            error("forcing must be :legendre or :nodal, got $forcing")
    else
        contact === :tangency ||
            error("the nonvariational formulation offers :tangency only, got $contact")
        forcing === :legendre ||
            error("forcing is a variational choice; the nonvariational solver has no interior")
    end
    Backend(formulation, contact, forcing)
end

"""Short human-readable name, stable enough to key a results store on."""
function label(b::Backend)
    b.formulation === :nonvariational && return "nonvar/tangency"
    f = b.forcing === :nodal ? "/nodal" : ""
    (b.contact === :lcp ? "var/lcp" : "var/active-set") * f
end

Base.show(io::IO, b::Backend) = print(io, "Backend(", label(b), ")")

"""
    run_impact(b::Backend; We, Bo, Oh, M, K, t_max, eta, eta_nonvar, save_every)

Run one impact and return a NamedTuple with the same fields whatever the backend:

  * `t`, `z`, `v`       -- time, centre-of-mass height and velocity
  * `zeta`              -- surface amplitudes per frame, a vector of vectors
  * `ls`                -- the harmonic degrees `zeta` is indexed by
  * `cp`                -- contact node count per frame
  * `cor`, `tc`         -- restitution and contact time under the PROXIMITY definition,
                           identical for every backend: contact whenever any surface point is
                           below `0.02R`, first touch to last release
  * `wall`              -- seconds
  * `ok`                -- whether the run produced usable metrics

`eta` is the dimensionless viscosity function for the variational backends. The nonvariational
solver cannot take one -- it has no interior field to evaluate it on -- so shear thinning is
passed to it as `eta_nonvar`, an `STExactParams`, and it is an error to give one without the
other.
"""
function run_impact(b::Backend; kw...)
    ## A solver that gives up is a RESULT, not an exception: a sweep must be able to record
    ## which cases a backend cannot do. Only genuine misuse -- an impossible Backend, or a
    ## viscosity handed to the formulation that cannot evaluate it -- is allowed to throw, and
    ## those are raised before any solving starts.
    t0 = time()
    try
        return _run_impact(b; kw...)
    catch e
        e isa ErrorException && occursin("eta", e.msg) && rethrow()
        return (t = Float64[], z = Float64[], v = Float64[], zeta = Vector{Float64}[],
                ls = Int[], cp = Int[], cor = NaN, tc = NaN, wall = time() - t0,
                ok = false, backend = label(b),
                diag = (min_gap = NaN, n_detected = 0, rejects = 0, lcp_resid = NaN,
                        eta_sweeps = 0, cor_internal = NaN, tc_internal = NaN))
    end
end

function _run_impact(b::Backend; We::Real, Bo::Real, Oh::Real, M::Int = 30, K::Int = 2,
                    t_max::Real = 25.0, eta = nothing, eta_nonvar = nothing,
                    save_every::Real = 0.005, h_thresh::Real = 0.02)
    t0 = time()
    if b.formulation === :variational
        eta_nonvar === nothing ||
            error("eta_nonvar is for the nonvariational backend; pass `eta` instead")
        p = ImpactParams(We = We, Bo = Bo, Oh = Oh, M = M, K = K, t_max = t_max,
                         force_mode = b.forcing,
                         eta = eta === nothing ? (gd -> 1.0) : eta)
        r = b.contact === :lcp ? simulate_lcp(p) : simulate(p)
        m = proximity_metrics(p, r; h_thresh = h_thresh)
        return (t = r.t, z = r.z, v = r.v,
                zeta = [surface_amplitudes(p, a) for a in r.a],
                ls = collect(p.ls), cp = r.cp,
                cor = m.cor, tc = m.tc, wall = time() - t0,
                ok = isfinite(m.cor) && isfinite(m.tc) && 0 < m.tc < t_max,
                backend = label(b),
                ## solver-specific diagnostics, kept so a stored row is still reproducible
                diag = (min_gap = m.min_gap, n_detected = m.n_detected,
                        rejects = r.rejects,
                        lcp_resid = hasproperty(r, :lcp_resid_max) ? r.lcp_resid_max : NaN,
                        eta_sweeps = hasproperty(r, :eta_sweeps_max) ? r.eta_sweeps_max : 0,
                        cor_internal = r.cor, tc_internal = r.tc))
    end
    ## nonvariational
    eta === nothing ||
        error("the nonvariational solver has no interior field; pass `eta_nonvar` instead")
    cfg = SimConstants(M, M + 1, Oh, Bo, make_theta_vec(M),
                       precompute_integrals(NaN, M)[1], make_dt_max(M); viscous = :reid)
    s = DropState(M); s.z = 1.0 + 4*h_thresh; s.v = -sqrt(We)
    s.dt = make_dt_max(M); s.cp = 0
    ts, st = eta_nonvar === nothing ?
        solve_drop!(cfg, OBParams(), s; t_end = t_max, save_every = save_every) :
        solve_drop!(cfg, OBParams(), s; stx = eta_nonvar, t_end = t_max, save_every = save_every)
    ths = range(pi/2, pi; length = 240)
    touch = [minimum(drop_height(x, th) for th in ths) < h_thresh for x in st]
    i = findfirst(touch); j = findlast(touch)
    ok = !(i === nothing || j === nothing || i == j)
    (t = ts, z = [x.z for x in st], v = [x.v for x in st],
     zeta = [x.A[2:end] for x in st], ls = collect(2:M), cp = [x.cp for x in st],
     cor = ok ? abs(st[j].v / st[i].v) : NaN,
     tc  = ok ? ts[j] - ts[i] : NaN,
     wall = time() - t0, ok = ok, backend = label(b),
     diag = (min_gap = minimum(minimum(drop_height(x, th) for th in ths) for x in st),
             n_detected = count(touch), rejects = 0, lcp_resid = NaN, eta_sweeps = 0,
             cor_internal = NaN, tc_internal = NaN))
end

"""
    drop_outline(zeta, ls, z; nth = 400) -> (x, y)

Cross-section of the drop for plotting: `r(theta) = 1 + sum zeta_l P_l(cos theta)` lifted by the
centre-of-mass height. Backend-independent, because `run_impact` returns `zeta` in a common
form -- which is the point of having it do so.
"""
function drop_outline(zeta::AbstractVector, ls::AbstractVector, z::Real; nth::Int = 400)
    th = range(0, pi; length = nth)
    r = similar(collect(th))
    for (k, t) in enumerate(th)
        mu = cos(t)
        s = 1.0
        for (i, l) in enumerate(ls)
            s += zeta[i] * legendre_angular(l, mu).P
        end
        r[k] = s
    end
    (r .* sin.(th), r .* cos.(th) .+ z)
end
