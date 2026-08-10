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
    check_converged(cor, tc, t_max, minz, label) -> Bool

Whether a finished march produced a bounce, as opposed to finishing.

A run that never releases still returns numbers. It reports `tc` a hair under
`t_max` and a restitution of order `1e-15`, and an `ok` test of the form
`0 < tc < t_max` waves it through. Every one of twenty-five nonvariational runs
passed that test while the drop sat on the substrate, and the audit that found it
had already drawn the wrong conclusion from it once.

So the test is on the physics, not on the absence of an exception, and a failure is
warned about rather than returned quietly: a caller who does not inspect `ok` should
still hear about it.
"""
function check_converged(cor, tc, t_max, minz, lbl; t_final = NaN, dt_final = NaN,
                         dt_min = NaN, released = true)
    reasons = String[]
    (isfinite(cor) && isfinite(tc))     || push!(reasons, "non-finite metrics")
    tc >= 0.9 * t_max                   && push!(reasons, "never released (tc = $(round(tc, digits=2)) of t_max = $t_max)")
    isfinite(cor) && cor <= 1e-6        && push!(reasons, "restitution is zero to machine precision")
    isfinite(cor) && cor > 1.0          && push!(reasons, "restitution exceeds 1, which is unphysical")
    isfinite(minz) && minz <= 1e-9      && push!(reasons, "centre of mass reached the substrate")
    ## A march can also stop EARLY. The tight-tolerance run at M = 60 died at
    ## t = 0.235 of a requested 25 with the drop still on the substrate, and every
    ## test above passes for it: tc is small, the restitution is finite and under
    ## one, the drop never reached z = 0. What marks it is that the integrator gave
    ## up rather than the drop leaving, so the tests are on the integrator.
    released || push!(reasons, "the march ended with the drop still in contact")
    isfinite(dt_final) && isfinite(dt_min) && dt_final <= 1.01 * dt_min &&
        push!(reasons, "step size collapsed to its floor (dt = $dt_final)")
    isempty(reasons) && return true
    @warn "run did not converge to a bounce; treat its metrics as meaningless" backend=lbl reasons
    false
end

"""
    run_impact(b::Backend; We, Bo, Oh, M, K, t_max, eta, eta_nonvar, eta_tol, ob, save_every)

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

`ob` is an [`OBParams`](@ref) for a viscoelastic drop, and likewise only the nonvariational
backend accepts it: the polymer stress needs a state the variational formulation does not
carry. Omitting it runs a Newtonian drop.
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
        ## A bad keyword or a missing method is a caller mistake. Swallowing it as
        ## "the solver could not do this case" is how a typo becomes a data point.
        (e isa MethodError || e isa UndefKeywordError) && rethrow()
        ## Handing a backend a parameter it structurally cannot use is misuse, not a case
        ## the solver "could not do". Both such guards name the offending parameter.
        e isa ErrorException &&
            (occursin("eta", e.msg) || occursin("Oldroyd-B", e.msg)) && rethrow()
        return (t = Float64[], z = Float64[], v = Float64[], zeta = Vector{Float64}[],
                ls = Int[], cp = Int[], cor = NaN, tc = NaN, wall = time() - t0,
                ok = false, backend = label(b),
                diag = (min_gap = NaN, n_detected = 0, rejects = 0, lcp_resid = NaN,
                        eta_sweeps = 0, cor_internal = NaN, tc_internal = NaN))
    end
end

function _run_impact(b::Backend; We::Real, Bo::Real, Oh::Real,
                    M::Int = DEFAULT_M, K::Int = DEFAULT_K,
                    t_max::Real = 25.0, eta = nothing, eta_nonvar = nothing,
                    save_every::Real = 0.005, h_thresh::Real = 0.02,
                    eta_tol = nothing, ob = nothing)
    t0 = time()
    if b.formulation === :variational
        eta_nonvar === nothing ||
            error("eta_nonvar is for the nonvariational backend; pass `eta` instead")
        ## Oldroyd-B needs a polymer-stress state, which the variational formulation does
        ## not carry. Refusing here is better than silently integrating a Newtonian drop.
        ob === nothing ||
            error("Oldroyd-B needs the nonvariational formulation; " *
                  "use Backend(formulation = :nonvariational, contact = :tangency)")
        p = ImpactParams(We = We, Bo = Bo, Oh = Oh, M = M, K = K, t_max = t_max,
                         force_mode = b.forcing, eta_tol = eta_tol,
                         eta = eta === nothing ? (gd -> 1.0) : eta)
        r = b.contact === :lcp ? simulate_lcp(p) : simulate(p)
        m = proximity_metrics(p, r; h_thresh = h_thresh)
        return (t = r.t, z = r.z, v = r.v,
                zeta = [surface_amplitudes(p, a) for a in r.a],
                ls = collect(p.ls), cp = r.cp,
                cor = m.cor, tc = m.tc, wall = time() - t0,
                ## A completed bounce ends with the drop moving UP. A march that the
                ## integrator abandoned ends wherever it gave up, which for the M = 60
                ## tight-tolerance case was mid-impact at t = 0.235. Comparing times or
                ## frame indices does not separate these, because the march terminates
                ## at release, so a good run has no frames after last contact either.
                ok = check_converged(m.cor, m.tc, t_max,
                                     isempty(r.z) ? NaN : minimum(r.z), label(b);
                                     released = !isempty(r.v) && r.v[end] > 0),
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
    ## `OBParams()` is the Newtonian default; a caller who wants a viscoelastic drop
    ## supplies one. Before this keyword existed there was NO route from `run_impact` to
    ## Oldroyd-B at all, while the solver-choice page told readers to use it here.
    obp = ob === nothing ? OBParams() : ob
    ts, st = eta_nonvar === nothing ?
        solve_drop!(cfg, obp, s; t_end = t_max, save_every = save_every) :
        solve_drop!(cfg, obp, s; stx = eta_nonvar, t_end = t_max, save_every = save_every)
    ths = range(pi/2, pi; length = 240)
    touch = [minimum(drop_height(x, th) for th in ths) < h_thresh for x in st]
    i = findfirst(touch); j = findlast(touch)
    touched = !(i === nothing || j === nothing || i == j)
    cor_ = touched ? abs(st[j].v / st[i].v) : NaN
    tc_  = touched ? ts[j] - ts[i] : NaN
    zs   = [x.z for x in st]
    ok   = touched && check_converged(cor_, tc_, t_max, isempty(zs) ? NaN : minimum(zs),
                                      label(b))
    (t = ts, z = zs, v = [x.v for x in st],
     zeta = [x.A[2:end] for x in st], ls = collect(2:M), cp = [x.cp for x in st],
     cor = cor_, tc = tc_,
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
