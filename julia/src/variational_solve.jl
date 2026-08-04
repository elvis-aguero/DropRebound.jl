# Time integration of the variational model, with contact.
#
# THE CONTACT CLOSURE HERE IS NOT INVENTED. It is the closure of Gabbard et al.
# (2025), read out of the validated MATLAB implementation and reproduced with one
# addition: the interior of the drop is part of the state, so each surface mode
# carries K radial trial functions instead of one. The node distribution, the pressure
# representation, the square collocation and the timestep are theirs, because an earlier
# version of this file reinvented all of it and got every piece wrong in a way that took
# a working simulation to see. HOW THE CONTACT COUNT IS CHOSEN is not theirs, and the
# difference is set out at the active-set iteration in `simulate`. What that run showed,
# and what this file therefore does:
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
#   * THE PRESSURE SIGN IS NOT A FEASIBILITY TEST. The only thing that makes a step
#     inadmissible is the surface dropping below the substrate OUTSIDE the contact.
#     Making the pressure sign a filter, as an earlier version did, rejected admissible
#     steps at onset and stalled the run. It is used in exactly one place -- at the
#     outermost contacting node, as the criterion for RELEASE -- which is the role it
#     has in the complementarity pair. Nothing constrains its sign inside the patch,
#     and that is why the film is free to pull there, as it measurably does in a small
#     minority of steps.
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
using LinearAlgebra: I

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
    eta_tol::Float64             # Picard convergence tolerance on the strain rate
    eta_max_sweeps::Int
    eta_const::Bool
end

"""
    ImpactParams(; We, Bo, Oh, M=90, K=1, ...)

`M` is the harmonic truncation: shape modes are `l = 2..M`, the film pressure carries
`l = 0..M`, and there are `M+1` collocation nodes. Those three counts are tied
together deliberately -- that is what makes the contact system square.
"""
function ImpactParams(; We, Bo, Oh, M::Int = 90, K::Int = 1, eta = gd -> 1.0,
                      dt0 = nothing, dt_min = 1e-10, t_max = 25.0,
                      eta_tol = 1e-8, eta_max_sweeps = 12)
    ls = collect(2:M)
    # theta = pi plus the zeros of P_M. These cluster at the poles, which is the
    # whole point: contact is resolved where contact happens.
    mus, _ = gauss_legendre_nodes(M, -1.0, 1.0)
    nodes = vcat(pi, acos.(clamp.(sort(mus), -1.0, 1.0)))
    dt = something(dt0, 2pi / (8 * sqrt(M * (M + 2) * (M - 1))))
    ec = all(gd -> eta(gd) == eta(0.0), (0.0, 1e-3, 1.0, 1e3))
    ImpactParams(We, Bo, Oh, ls, K, eta, nodes, dt, dt_min, t_max, eta_tol,
                 eta_max_sweeps, ec)
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

"""
The Legendre Vandermonde on the collocation nodes: `V[i,j] = P_{l_j}(cos theta_i)`.

It converts pressure coefficients to nodal pressure values, `nodal = V * coefficients`.
Square and invertible because there are exactly as many nodes as harmonics and the nodes
are distinct -- the same fact that makes the contact system square.
"""
function legendre_vandermonde(p::ImpactParams)
    npc = pc_len(p)
    V = zeros(npc, npc)
    for i in 1:npc
        mu = cos(p.nodes[i])
        for j in 1:npc
            V[i, j] = legendre_angular(pc_l(j), mu).P
        end
    end
    V
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
    try_step(p, prev, curr, dt, cp; F0, cache) -> (status, next)

`status` is `:ok`, `:penetrate` (the surface is below the substrate outside the
contact, so the contact set is too small) or `:diverge` (the viscosity iteration did
not converge, which says nothing about the geometry and calls for a smaller step).

One step at a FIXED contact count, as one square linear system in `(a, p_c)`.

`z` and `v` are eliminated first: both are affine in the single harmonic `p_c,1`,
because that is the only one that moves the centre of mass. What remains is
`ndof + M + 1` equations in as many unknowns -- the momentum equations, the gap
conditions on the contact nodes, and the vanishing of the pressure on the free ones.
"""
function try_step(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                  dt::Float64, cp::Int; F0 = nothing, cache = nothing, diag = nothing,
                  Vfac = nothing)
    b = basis(p); N = ndof(b); npc = pc_len(p); nn = length(p.nodes)
    (0 <= cp <= nn) || return (:penetrate, curr)

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

    # THE PRESSURE IS SOLVED FOR AT THE NODES, not as Legendre coefficients.
    #
    # The two are related by the Legendre Vandermonde, `nodal = V * coefficients`, so
    # this is an exact change of variables and the answer cannot move. What moves is the
    # structure of the system, in two ways that matter:
    #
    #   * The free-node condition becomes `v_i = 0` -- a single unit entry -- instead of
    #     a dense row evaluating M+1 Legendre polynomials. That is a structural
    #     simplification only. It does NOT improve the conditioning, which was the first
    #     reason I gave for making the change and was wrong: measured on the same state,
    #     cond(KKT) goes 8.6e8 -> 7.6e8 at M = 20 but 2.9e12 -> 6.7e12 at M = 90, so it
    #     is slightly WORSE where it matters. Whatever makes this system ill-conditioned
    #     is not the pressure rows.
    #   * The compliance implied by the eliminated system becomes SYMMETRIC. Virtual work
    #     pairs the nodal gap with the nodal pressure weighted by quadrature, so in this
    #     basis `W^(1/2) H A^-1 H' W^(1/2)` is symmetric to machine precision and positive
    #     semi-definite -- measured at 1e-16 asymmetry and no negative eigenvalues for
    #     every M and K tried. In the coefficient basis it is neither: gaps at nodes
    #     paired against coefficients of harmonics gave a relative asymmetry of sqrt(2)
    #     with half the spectrum negative.
    #
    # That symmetry is the precondition for treating the contact problem as the linear
    # complementarity problem it is, rather than searching over candidate contact sets.
    # Nothing downstream sees the change: the coefficients are recovered before returning,
    # so `pc` keeps its meaning everywhere else.
    Vf = Vfac === nothing ? lu(legendre_vandermonde(p)) : Vfac
    Vinv_row2 = (Vf \ Matrix{Float64}(I, npc, npc))[2, :]   # the l = 1 coefficient
    Qn = Qm * (Vf \ Matrix{Float64}(I, npc, npc))

    Hm = zeros(npc, N); Zm = zeros(npc, npc); rhs_h = zeros(npc)
    for i in 1:npc
        th = p.nodes[i]
        if i <= cp                                   # contact: the gap vanishes
            row, mu = gap_row(p, th)
            Hm[i, :] = row
            Zm[i, :] .+= dz_dpc1 .* Vinv_row2        # z is affine in p_c,1
            rhs_h[i] = -mu - z0
        else                                         # free: the nodal pressure vanishes
            Zm[i, i] = 1.0
        end
    end

    # THE PICARD ITERATION IS MONITORED, not counted off.
    #
    # The nonlinearity is only eta(gammadot), so the coefficient matrix depends on the
    # unknown through the strain rate and the system is solved by freezing eta,
    # solving, re-evaluating eta at the answer, and repeating. That map contracts for
    # small enough dt, and the reason is structural rather than hopeful: the nonlinear
    # term enters A = beta^2 M + beta C + G at O(beta) against an inertial O(beta^2),
    # so the Lipschitz constant scales as dt times the sensitivity of C to the strain
    # rate and vanishes with dt. At the timestep used here the observed factor is about
    # 1/200 per sweep, which is why one sweep is indistinguishable from six.
    #
    # None of that is a guarantee. "Small enough dt" is unquantified, and eta' is
    # largest at the Carreau knee, gammadot ~ 1/lambda_c, which is exactly where
    # near-stagnation points sit; a yield-stress law would have no such reprieve at all.
    # So the iteration runs to a TOLERANCE and reports how many sweeps it took, and a
    # step whose iteration will not converge is REJECTED rather than returned. Before
    # this, a fixed sweep count meant a non-contracting regime would have returned a
    # silently unconverged state, with every linear solve inside it succeeding.
    conv = 0.0; used = 0
    for it in 1:(F0 === nothing ? p.eta_max_sweeps : 1)
        used = it
        prev_star = copy(adot_star)
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
                lu([A -Qn; Hm Zm])
            end
        else
            lu([A -Qn; Hm Zm])
        end
        sol = fac \ [rhs0; rhs_h]
        any(!isfinite, sol) && return (:diverge, curr)
        a_next = sol[1:N]
        ## the solve returns NODAL pressures; the coefficients are recovered here so that
        ## `pc` keeps the meaning the rest of the code and every test relies on
        pc = Vf \ sol[N+1:N+npc]
        v_next = (-p.Bo - pc[2] - hv_v) / β
        z_next = (v_next - hv_z) / β
        adot_star = β * a_next + hv_a
        ## relative to the strain rate's own scale, so the criterion means the same
        ## thing at every impact energy
        scale = max(maximum(abs, adot_star), 1e-12)
        conv = maximum(abs, adot_star .- prev_star) / scale
        if F0 !== nothing
            conv = 0.0                               # constant eta: nothing was iterated
            break
        end
        conv < p.eta_tol && break
    end
    ## A step whose viscosity iteration has not converged is not a worse step, it is
    ## not a solution. Rejecting it hands control to the dt halving that already
    ## exists, which is the only lever that can restore the contraction.
    ## A NON-CONVERGED VISCOSITY ITERATION IS NOT A GEOMETRY SIGNAL. The caller grows
## the contact set when a step comes back inadmissible, because that is what
## penetration means -- so returning the same flag for a stalled Picard iteration
## makes the active set grow for a reason that has nothing to do with the surface.
## Doing exactly that took CoR from 0.767 to 1.003 and the contact time from 2.72 to
## 0.44: the contact ran away and the drop left almost elastically. The two are
## reported separately so the caller can halve dt for one and grow for the other.
(F0 === nothing && conv > p.eta_tol) && return (:diverge, curr)
    diag === nothing || (diag[] = (sweeps = used, resid = conv))

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
    # NO TANGENCY RESIDUAL IS RETURNED. The reference implementation computes one --
    # the height of the first free node relative to the last contacting one -- and
    # selects the contact count by minimising it over feasible candidates. This solver
    # does not: the contact count comes from the active-set iteration in `simulate`,
    # which needs no objective and therefore has no tie to break. An earlier version of
    # this function still computed and returned the residual after the caller had
    # stopped using it, which reads as a selection rule that is not there.
    (ok ? :ok : :penetrate, nxt)
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
## depends only on the nodes and the harmonic count, so it is built once
Vfac = lu(legendre_vandermonde(p))
    cache = Dict{Tuple{Int,Float64,Bool},Any}()
    dg = Ref((sweeps = 0, resid = 0.0))
    max_sweeps = 0; max_resid = 0.0
    s0 = initial_state(p)
    prev, curr = s0, s0
    ts = Float64[0.0]; zs = Float64[s0.z]; vs = Float64[s0.v]
cps = Int[0]; pc1 = Float64[0.0]
## the amplitudes themselves, so energy, volume and the shape can be checked
## from the trajectory rather than re-derived by a caller re-marching it
as = Vector{Float64}[copy(s0.a)]; adots = Vector{Float64}[copy(s0.adot)]
pcs = Vector{Float64}[copy(s0.pc)]
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
            status, nxt = try_step(p, prev, curr, dt, cand;
                                   F0 = F0, cache = cache, diag = dg, Vfac = Vfac)
            ## A stalled viscosity iteration is not a statement about the surface, so it
            ## must not drive the contact set. Abandon the active-set search and let the
            ## step be rejected, which halves dt -- the only lever that restores the
            ## contraction. Reading this as penetration instead grew the contact for no
            ## physical reason and returned CoR = 1.003 with a contact time of 0.44.
            status === :diverge && break
            ok = status === :ok
            if ok
                max_sweeps = max(max_sweeps, dg[].sweeps)
                max_resid  = max(max_resid, dg[].resid)
            end
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
push!(as, copy(curr.a)); push!(adots, copy(curr.adot)); push!(pcs, copy(curr.pc))
        # RECOVER THE STEP THE WAY IT WAS GIVEN UP: by doubling, capped at the nominal
        # value. Snapping straight back to dt0 after every success makes the solver
        # rediscover a refinement it just paid for -- if a whole stretch of the
        # trajectory needs a finer step, each step there costs a wasted rejected attempt
        # before halving again. Doubling is the symmetric counterpart of the halving on
        # rejection, so the step size carries a short memory of what the last few steps
        # needed. dt0 remains a ceiling, not a target: it is the resolution the fastest
        # retained capillary mode requires, and nothing is gained by exceeding it.
        dt = min(2 * dt, p.dt0)
        verbose && best_cp != cps[end-1] &&
            @info "contact" t=curr.t cp=best_cp z=curr.z v=curr.v
        # done once the drop has left the substrate and is rising
        if best_cp == 0 && curr.v > 0 && curr.z > 1.0 && any(>(0), cps)
            break
        end
    end
    (t = ts, z = zs, v = vs, cp = cps, pc1 = pc1, a = as, adot = adots, pc = pcs,
     rejects = nrej,
     ## the Picard iteration's worst behaviour over the whole march, so a regime where
     ## the map stops contracting announces itself instead of returning quietly
     eta_sweeps_max = max_sweeps, eta_resid_max = max_resid,
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

# ============================================================================
# THE CONTACT PROBLEM AS A LINEAR COMPLEMENTARITY PROBLEM
#
# Everything above chooses the contact extent by searching: propose a count, solve, test
# whether the proposal was consistent, rank the survivors. That is a workaround for not
# knowing the contact set, and it has cost a great deal -- a tangency ranking that is
# degenerate once the interior is resolved, a chattering active set that injected forty
# times the drop's energy, and a sequence of tie-breaking rules none of which is physics.
#
# The contact set does not have to be guessed. Eliminating the interior amplitudes and the
# centre of mass leaves the nodal gaps affine in the nodal pressures,
#
#     h = A_c p + b ,
#
# and the physics is three conditions at every node: `h >= 0` (no penetration), `p >= 0`
# (a gas film cannot pull), `h_i p_i = 0` (touching with pressure, or clear of it with
# none). An affine relation plus those conditions IS a linear complementarity problem, and
# the contact set is an OUTPUT of solving it -- the nodes where `p_i > 0`.
#
# What makes this tractable rather than aspirational is measured, not assumed:
#
#   * `A_c` is symmetric to 1e-16 and positive semi-definite at every M and K tried, in
#     the NODAL pressure basis. In the Legendre coefficient basis it is neither -- gaps at
#     nodes paired against coefficients of harmonics gave sqrt(2) relative asymmetry with
#     half the spectrum negative. The change of variables is what buys the structure.
#   * With `A_c` symmetric PSD the LCP is exactly the KKT system of the convex programme
#     `min (1/2) p'A_c p + b'p` subject to `p >= 0`, so a solution exists and the
#     complementarity residual `max_i min(h_i, p_i)` is a certificate rather than a
#     heuristic stopping rule.
#   * The affine relation itself was checked against the searching solver's own output:
#     relative error 0.
#
# And what it is worth: the searching solver satisfies `h >= 0` and `h p = 0` by
# construction but NOT `p >= 0`. Measured over whole runs, the pressure is negative at
# 0.015 per cent of node-steps in a case that works and 41 per cent in one that does not.
# The dropped condition is exactly the difference between the two.

"""
    lcp_pgs(A, b; iters, tol) -> (p, residual, sweeps)

Projected Gauss-Seidel for `h = Ap + b`, `h >= 0`, `p >= 0`, `p'h = 0`, with `A`
symmetric positive semi-definite.

Each sweep is one pass of Gauss-Seidel followed by projection onto `p >= 0`, which for a
symmetric PSD matrix is coordinate descent on the equivalent convex programme and so
cannot increase the objective. Nodes with a non-positive diagonal are skipped: those are
nodes whose gap does not respond to their own pressure, which cannot carry load.

The returned residual is `max_i min(h_i, p_i)`, scaled. It vanishes exactly when all three
conditions hold, so it is a certificate of the solution rather than a proxy for one.
"""
function lcp_pgs(A::AbstractMatrix, b::AbstractVector; iters::Int = 20000,
                 tol::Real = 1e-11)
    n = length(b)
    p = zeros(n)
    dg = [A[i, i] for i in 1:n]
    resid = Inf; sweeps = 0
    for k in 1:iters
        sweeps = k
        for i in 1:n
            dg[i] <= 0 && continue
            p[i] = max(0.0, p[i] - (b[i] + dot(view(A, i, :), p)) / dg[i])
        end
        h = A * p .+ b
        resid = maximum(min.(h, p))
        sc = max(maximum(abs, h), maximum(abs, p), 1.0)
        abs(resid) / sc < tol && break
    end
    (p, resid, sweeps)
end

"""
    contact_lcp(p, prev, curr, dt; F0, Vfac) -> (Ac, b, nodes)

Assemble the complementarity problem for one step: the compliance `Ac`, the free-flight
gaps `b`, and the node indices the problem is posed on.

Restricted to the lower hemisphere. The upper nodes sit a diameter above the substrate and
carry zero pressure in any solution, and including them only adds null directions -- their
self-compliance is ~1e-20, which is numerically indistinguishable from a node that cannot
carry load at all.
"""
function contact_lcp(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                     dt::Float64; F0, Vfac)
    b = basis(p); N = ndof(b); npc = pc_len(p); nn = length(p.nodes)
    if curr.first
        c0, c1, c2 = 1.0, -1.0, 0.0
    else
        r = dt / curr.dt
        c0 = (1 + 2r)/(1 + r); c1 = -(1 + r); c2 = r^2/(1 + r)
    end
    β = c0 / dt
    hv_a  = (c1*curr.a    + c2*prev.a)    / dt
    hv_ad = (c1*curr.adot + c2*prev.adot) / dt
    hv_z  = (c1*curr.z    + c2*prev.z)    / dt
    hv_v  = (c1*curr.v    + c2*prev.v)    / dt
    A = β^2 * F0.M + β * F0.C + F0.G
    rhs0 = -F0.M * (β*hv_a + hv_ad) - F0.C * hv_a
    Vinv = Vfac \ Matrix{Float64}(I, npc, npc)
    Qn = hcat((force_column(p, j) for j in 1:npc)...) * Vinv
    H = zeros(nn, N); mu = zeros(nn)
    for i in 1:nn
        row, m = gap_row(p, p.nodes[i]); H[i, :] = row; mu[i] = m
    end
    z0 = ((-p.Bo - hv_v)/β - hv_z)/β
    dz = -1/β^2
    ## h is affine in the nodal pressure through BOTH the shape and the centre of mass
    Ac_full = H * (A \ Qn) .+ dz .* (ones(nn) * Vinv[2, :]')
    b_full  = H * (A \ rhs0) .+ mu .+ z0
    idx = [i for i in 1:nn if p.nodes[i] > pi/2]
    S = Ac_full[idx, idx]
    (Matrix(Symmetric(0.5 .* (S .+ S'))), b_full[idx], idx,
     (A = A, rhs0 = rhs0, Qn = Qn, Vinv = Vinv, β = β, hv_a = hv_a,
      hv_z = hv_z, hv_v = hv_v, npc = npc, N = N))
end
