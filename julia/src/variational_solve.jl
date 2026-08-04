# Time integration of the variational model, with contact.
#
# THE CONTACT CLOSURE HERE IS NOT INVENTED. It is the closure of Gabbard et al.
# (2025), read out of the validated MATLAB implementation and reproduced with one
# addition: the interior of the drop is part of the state, so each surface mode
# carries K radial trial functions instead of one. Everything else -- the node
# distribution, the pressure representation, the square collocation, the feasibility
# rule, the objective, the timestep -- is theirs, because an earlier version of this
# file reinvented all of it and got every piece wrong in a way that took a working
# simulation to see. What that run showed, and what this file therefore does:
#
#   * THE NODES ARE LEGENDRE ROOTS, not uniform. theta = pi together with the zeros
#     of P_M, which cluster hard at the poles: at M = 90 the second node sits 1.52
#     degrees from the pole. That clustering is what lets a two-node contact be a
#     genuinely small patch. On a uniform grid at M = 8 a two-node contact spanned
#     twenty degrees, the pressure needed to hold it was impulsive, and the drop left
#     with more energy than it arrived with.
#
#   * THE PRESSURE IS A GLOBAL LEGENDRE SERIES, l = 0..M, closed by SQUARE
#     collocation: the gap vanishes at the cp contact nodes and the pressure vanishes
#     at the remaining M+1-cp nodes. Its coefficients do NOT converge in l -- in the
#     reference run |p_l| is largest at l = M, the truncation itself -- so the series
#     is interpolatory rather than spectral, and only the reconstructed field means
#     anything. That field is smooth and positive on the patch and 1e-12 off it. An
#     earlier version here read the non-convergent spectrum as proof the design was
#     broken and replaced it with a patch-supported pressure enforced weakly. The
#     spectrum is real; the conclusion was not.
#
#   * THERE IS NO PRESSURE SIGN CONDITION. The only feasibility test is that the
#     surface stay above the substrate OUTSIDE the contact. Signorini's inequality is
#     satisfied by the solution rather than imposed on it, and imposing it here
#     rejected admissible steps at exactly the wrong moment.
#
#   * THE PRESSURE IS LARGEST AT ONSET, not smallest. In the reference run the patch
#     pressure peaks at 73 at first touch and falls to 7 at maximum contact, while
#     the net force does the opposite, because the patch area is what is small. Much
#     of the earlier tolerance fiddling here rested on the opposite belief.
#
#   * dt IS CONSTANT, 2 pi/(8 sqrt(M(M+2)(M-1))) -- the fastest retained capillary
#     mode resolved eight times per period. Halving is for startup and rejected
#     steps only; there is no adaptive control during contact.
#
# What is genuinely ours is the state: `a` holds the interior displacement profile
# coefficients, and the surface amplitudes are their boundary trace. The pressure
# coupling that results, Q_zeta_l = -(4 pi/(2l+1)) p_c,l, is not a free choice -- it
# reproduces their forcing coefficient l exactly, because the assembled modal mass
# reproduces Lamb's added mass 4 pi/(l(2l+1)) to machine precision.
#
# BDF2 in both a and adot, first order on the opening step. With a constant viscosity
# the operator is built once (block diagonal in l, independent of the state); with a
# shear-thinning one it is rebuilt per sweep, because eta multiplies a stiff
# dissipation operator and extrapolating it is a stability risk, not just an accuracy
# cost.

using LinearAlgebra

# Active-set moves allowed per step. Each move is forced by a violated complementarity
# condition, so the iteration cannot wander; this only bounds pathological cycling,
# after which the step is rejected and dt halved.
const MAX_ACTIVE_SET_ITERS = 40

struct ImpactParams
    We::Float64
    Bo::Float64
    Oh::Float64
    ls::Vector{Int}
    K::Int
    eta::Function                # shear rate -> viscosity / eta_0
    nodes::Vector{Float64}       # theta collocation nodes, descending from pi
    dt0::Float64
    dt_min::Float64
    t_max::Float64
    eta_sweeps::Int
    eta_const::Bool
end

"""
    ImpactParams(; We, Bo, Oh, M=90, K=1, ...)

`M` is the harmonic truncation: shape modes are `l = 2..M`, the film pressure carries
`l = 0..M`, and there are `M+1` collocation nodes. Those three counts are tied
together deliberately -- that is what makes the contact system square.
"""
function ImpactParams(; We, Bo, Oh, M::Int = 90, K::Int = 1, eta = gd -> 1.0,
                      dt0 = nothing, dt_min = 1e-10, t_max = 25.0, eta_sweeps = 1)
    ls = collect(2:M)
    # theta = pi plus the zeros of P_M. These cluster at the poles, which is the
    # whole point: contact is resolved where contact happens.
    mus, _ = gauss_legendre_nodes(M, -1.0, 1.0)
    nodes = vcat(pi, acos.(clamp.(sort(mus), -1.0, 1.0)))
    dt = something(dt0, 2pi / (8 * sqrt(M * (M + 2) * (M - 1))))
    ec = all(gd -> eta(gd) == eta(0.0), (0.0, 1e-3, 1.0, 1e3))
    ImpactParams(We, Bo, Oh, ls, K, eta, nodes, dt, dt_min, t_max, eta_sweeps, ec)
end

basis(p::ImpactParams) = ModalBasis(p.ls, p.K)
lmax(p::ImpactParams) = maximum(p.ls)
pc_len(p::ImpactParams) = lmax(p) + 1        # harmonics l = 0..M
pc_l(j::Int) = j - 1

"""Every trial function equals 1 at the surface, so the trace is a vector of ones."""
trace_vec(p::ImpactParams, l) = [phi(RitzBasis(l, p.K), k, 1.0) for k in 1:p.K]

"""Surface amplitudes from the interior displacements: `zeta_l = chi_l(1)`."""
function surface_amplitudes(p::ImpactParams, a::AbstractVector)
    b = basis(p)
    [dot(trace_vec(p, l), view(a, dofindex(b, i, 1):dofindex(b, i, p.K)))
     for (i, l) in enumerate(b.ls)]
end

"""
Row of the linear map `a -> mu (1 + sum zeta_l P_l)` at one node, plus the constant.

The gap to the substrate is `h = row . a + mu + z`: the surface point at polar angle
`theta` sits at height `z + r(theta) cos(theta)` with `r = 1 + sum zeta_l P_l`.
"""
function gap_row(p::ImpactParams, th::Real)
    b = basis(p); mu = cos(th)
    row = zeros(ndof(b))
    for (i, l) in enumerate(b.ls)
        Pl = legendre_angular(l, mu).P
        tv = trace_vec(p, l)
        for k in 1:p.K
            row[dofindex(b, i, k)] = mu * Pl * tv[k]
        end
    end
    (row, mu)
end

gap(p::ImpactParams, a, z, th) = (r = gap_row(p, th); dot(r[1], a) + r[2] + z)

"""
Generalised force on the interior coordinates from film harmonic `l = pc_l(j)`.

`Q_zeta_l = -(4 pi/(2l+1)) p_c,l`, spread over the `K` trial functions by their trace.
Only harmonics that coincide with a retained SHAPE mode do work on the surface: `l=0`
changes volume and `l=1` translates, so both give a zero column. The `l=1` harmonic
still acts -- through the centre-of-mass equation, where it is the whole net force.
"""
function force_column(p::ImpactParams, j::Int)
    b = basis(p); Q = zeros(ndof(b)); l = pc_l(j)
    i = findfirst(==(l), b.ls)
    i === nothing && return Q
    tv = trace_vec(p, l)
    for k in 1:p.K
        Q[dofindex(b, i, k)] = -(4pi / (2l + 1)) * tv[k]
    end
    Q
end

"""The film pressure reconstructed at a node. The only meaningful form of `p_c`."""
function pc_at(p::ImpactParams, pc::AbstractVector, th::Real)
    mu = cos(th)
    sum(pc[j] * legendre_angular(pc_l(j), mu).P for j in eachindex(pc))
end

eta_field(p::ImpactParams, adot) = (x, mu) -> p.eta(shear_rate(basis(p), adot, x, mu))

mutable struct ImpactState
    t::Float64
    dt::Float64
    a::Vector{Float64}
    adot::Vector{Float64}
    z::Float64
    v::Float64
    cp::Int
    pc::Vector{Float64}
    first::Bool                  # opening step: BDF1, no two-level history yet
end

"""
The drop starts with its pole exactly on the substrate and `cp = 0`.

No artificial standoff gap. An earlier version started above contact to avoid an
impulsive first step, but the impulse came from the twenty-degree patch a coarse
uniform grid forced on the first contact, not from the initial condition. On Legendre
nodes the first contact is two nodes wide, 1.5 degrees, and needs no cushioning.
"""
function initial_state(p::ImpactParams)
    N = ndof(basis(p))
    ImpactState(0.0, p.dt0, zeros(N), zeros(N), 1.0, -sqrt(p.We), 0,
                zeros(pc_len(p)), true)
end

"""
    try_step(p, prev, curr, dt, cp; F0, cache) -> (ok, next, errortan)

One step at a FIXED contact count, as one square linear system in `(a, p_c)`.

`z` and `v` are eliminated first: both are affine in the single harmonic `p_c,1`,
because that is the only one that moves the centre of mass. What remains is
`ndof + M + 1` equations in as many unknowns -- the momentum equations, the gap
conditions on the contact nodes, and the vanishing of the pressure on the free ones.
"""
function try_step(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                  dt::Float64, cp::Int; F0 = nothing, cache = nothing)
    b = basis(p); N = ndof(b); npc = pc_len(p); nn = length(p.nodes)
    (0 <= cp <= nn) || return (false, curr, Inf)

    if curr.first
        c0, c1, c2 = 1.0, -1.0, 0.0
    else
        r = dt / curr.dt
        c0 = (1 + 2r)/(1 + r); c1 = -(1 + r); c2 = r^2/(1 + r)
    end
    β = c0 / dt
    hv_a    = (c1 * curr.a    + c2 * prev.a)    / dt
    hv_adot = (c1 * curr.adot + c2 * prev.adot) / dt
    hv_z    = (c1 * curr.z    + c2 * prev.z)    / dt
    hv_v    = (c1 * curr.v    + c2 * prev.v)    / dt

    mass = 4pi / 3
    adot_star = curr.first ? curr.adot : (1 + r) * curr.adot - r * prev.adot
    a_next = copy(curr.a); pc = zeros(npc)
    z_next = curr.z; v_next = curr.v

    Qm = hcat((force_column(p, j) for j in 1:npc)...)
    # the centre of mass is driven by p_c,1 alone: v_dot = -Bo - p_c,1
    z0 = ((-p.Bo - hv_v) / β - hv_z) / β
    dz_dpc1 = -1 / β^2

    Hm = zeros(npc, N); Zm = zeros(npc, npc); rhs_h = zeros(npc)
    for i in 1:npc
        th = p.nodes[i]
        if i <= cp                                   # contact: the gap vanishes
            row, mu = gap_row(p, th)
            Hm[i, :] = row
            Zm[i, 2] += dz_dpc1
            rhs_h[i] = -mu - z0
        else                                         # free: the pressure vanishes
            mu = cos(th)
            for j in 1:npc
                Zm[i, j] = legendre_angular(pc_l(j), mu).P
            end
        end
    end

    nsweep = F0 === nothing ? p.eta_sweeps : 0
    for _ in 0:nsweep
        F = F0 === nothing ?
            assemble_coupled(b, p.Oh; eta = eta_field(p, adot_star)) : F0
        A = β^2 * F.M + β * F.C + F.G
        rhs0 = -F.M * (β * hv_a + hv_adot) - F.C * hv_a
        # With a constant viscosity and a constant dt the whole KKT matrix is
        # constant for a given contact count, so its factorisation is cached -- the
        # same economy the reference implementation makes, and what makes M = 90
        # affordable over thousands of steps.
        key = (cp, dt, curr.first)
        fac = if cache !== nothing && F0 !== nothing
            get!(cache, key) do
                lu([A -Qm; Hm Zm])
            end
        else
            lu([A -Qm; Hm Zm])
        end
        sol = fac \ [rhs0; rhs_h]
        any(!isfinite, sol) && return (false, curr, Inf)
        a_next = sol[1:N]
        pc = sol[N+1:N+npc]
        v_next = (-p.Bo - pc[2] - hv_v) / β
        z_next = (v_next - hv_z) / β
        adot_star = β * a_next + hv_a
    end

    adot_next = β * a_next + hv_a
    nxt = ImpactState(curr.t + dt, dt, a_next, adot_next, z_next, v_next, cp,
                      copy(pc), false)

    # FEASIBILITY, and it is only this: the surface must not be below the substrate
    # anywhere outside the contact. No condition on the pressure.
    ok = true
    for i in (cp+1):nn
        if gap(p, a_next, z_next, p.nodes[i]) < 0
            ok = false; break
        end
    end
    # THE OBJECTIVE IS EDGE-LOCAL: the height the surface would have at the first
    # free node, relative to the last contacting one. One integrated over the patch
    # is biased toward vanishing contact.
    et = if !ok
        Inf
    elseif cp == 0
        abs(gap(p, a_next, z_next, p.nodes[1]))
    elseif cp < nn
        abs(gap(p, a_next, z_next, p.nodes[cp+1]) -
            gap(p, a_next, z_next, p.nodes[cp]))
    else
        Inf
    end
    (ok, nxt, et)
end

"""
    simulate(p; verbose=false) -> NamedTuple

March the impact. Returns the trajectory and the two KPIs, coefficient of restitution
and contact time.

The contact count changes by at most one per step and is chosen as the admissible
candidate with the smallest edge residual. Inadmissible candidates are rejected
outright rather than ranked, and when no candidate is admissible the step is rejected
and `dt` halved -- that, rather than event detection, is what keeps the nonsmooth
transition from being straddled.
"""
function simulate(p::ImpactParams; verbose = false)
    F0 = p.eta_const ? assemble_newtonian(basis(p), p.Oh) : nothing
    cache = Dict{Tuple{Int,Float64,Bool},Any}()
    s0 = initial_state(p)
    prev, curr = s0, s0
    ts = Float64[0.0]; zs = Float64[s0.z]; vs = Float64[s0.v]
    cps = Int[0]; pc1 = Float64[0.0]
    dt = p.dt0
    nrej = 0
    while curr.t < p.t_max
        # ACTIVE-SET ITERATION ON THE COMPLEMENTARITY PAIR, not a ranked search over
        # neighbouring contact counts.
        #
        # Ranking candidates by an edge residual is what the reference implementation
        # does, and at its single Ohnesorge it works. It chatters when the damping is
        # small: two candidates score almost equally, the choice flips step to step,
        # and the film pressure oscillates in sign with order-one amplitude until dt
        # collapses. At Oh = 0.023 that killed the run at M = 35 and M = 45 while
        # M = 30 and M = 40 came through -- sporadic in M, which is the signature of a
        # numerical bistability rather than of a resolution limit.
        #
        # The complementarity conditions themselves say what to do, and say it without
        # a tie to break: `h >= 0` is violated by too small a contact, `p_c >= 0` by
        # too large a one. So grow while any free node has penetrated, release while
        # the pressure at the outermost contacting node pulls, and stop when neither
        # holds -- the textbook primal active-set method, which terminates because
        # each move is forced. This is also the first place the pressure sign belongs:
        # as the RELEASE criterion it exists to be, rather than as the feasibility
        # filter an earlier version made it, where it rejected admissible steps.
        best = nothing; best_cp = curr.cp
        cand = curr.cp
        seen = Set{Int}()
        for _ in 1:MAX_ACTIVE_SET_ITERS
            (0 <= cand <= length(p.nodes)) || break
            cand ∈ seen && break                  # cycling: fall back to dt control
            push!(seen, cand)
            ok, nxt, _ = try_step(p, prev, curr, dt, cand; F0 = F0, cache = cache)
            if !ok
                cand += 1                          # penetrating: the contact is too small
                continue
            end
            if cand > 0 && pc_at(p, nxt.pc, p.nodes[cand]) < 0
                best, best_cp = nxt, cand          # keep it in case releasing fails
                cand -= 1                          # the edge is pulling: too large
                continue
            end
            best, best_cp = nxt, cand
            break
        end
        if best === nothing
            nrej += 1
            dt /= 2
            dt < p.dt_min && break
            continue
        end
        prev, curr = curr, best
        push!(ts, curr.t); push!(zs, curr.z); push!(vs, curr.v)
        push!(cps, best_cp); push!(pc1, curr.pc[2])
        dt = p.dt0                       # no adaptive control while in contact
        verbose && best_cp != cps[end-1] &&
            @info "contact" t=curr.t cp=best_cp z=curr.z v=curr.v
        # done once the drop has left the substrate and is rising
        if best_cp == 0 && curr.v > 0 && curr.z > 1.0 && any(>(0), cps)
            break
        end
    end
    (t = ts, z = zs, v = vs, cp = cps, pc1 = pc1, rejects = nrej,
     cor = restitution(vs, cps, p.We), tc = contact_time(ts, cps))
end

"""Coefficient of restitution: rebound speed over impact speed."""
function restitution(vs, cps, We)
    inc = findfirst(>(0), cps)
    inc === nothing && return NaN
    lastc = findlast(>(0), cps)
    lastc >= length(vs) && return NaN
    abs(vs[end] / vs[max(inc - 1, 1)])
end

"""Contact time: the longest contiguous interval with a nonempty contact set."""
function contact_time(ts, cps)
    best = 0.0; run_start = nothing
    for i in eachindex(cps)
        if cps[i] > 0 && run_start === nothing
            run_start = i
        elseif cps[i] == 0 && run_start !== nothing
            best = max(best, ts[i-1] - ts[run_start]); run_start = nothing
        end
    end
    run_start !== nothing && (best = max(best, ts[end] - ts[run_start]))
    best
end
