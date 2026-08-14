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
"""
Default angular truncation. `M+1` collocation nodes carry the contact constraint, so
this is not only a spectral resolution: it is how many distinct contact radii the
discretisation can represent at all. At `M = 45` a patch covering a tenth of the
lower hemisphere is a couple of nodes wide, which is too coarse to follow contact
growth at high Weber number even where the restitution has stopped moving.
"""
const DEFAULT_M = 90

"""
Default radial truncation. `K = 1` is potential flow, which over-damps by up to
170 per cent at experimental Ohnesorge; `K = 2` is the bare minimum and `K = 3`
reproduces Reid's exact damping to under two per cent across the validated range.
"""
const DEFAULT_K = 3

"""
    default_eta_tol(M) -> Float64

Default Picard tolerance on the interior strain rate, relative and in the infinity norm.
It scales with the truncation because the reachable residual does.

The iterate is `adot_star = beta*a + hv_a` with `beta = c0/dt` and `dt ~ M^(-3/2)`, so
the sweep-to-sweep increment is a difference of large quantities and its floor rises with
`M`. Measured on the 3000 ppm fluid, the floor is below 1e-6 at M = 45 and M = 60 and
3.0e-6 at M = 90, which a cubic in `M` tracks: `4e-12 * M^3` passes through all three.

The default sits about three decades above that fit rather than on it. A tolerance below
the floor is not merely slow, it is fatal: the step is rejected, `dt` halves, the floor
rises as `1/dt`, and the march dies having computed part of a bounce.
"""
default_eta_tol(M::Integer) = 1e-8 * float(M)^3

"""
How many consecutive step reductions a viscosity failure is allowed before the march
gives up.

Calibrated against measurement, and the margin is narrow. Rescued episodes on the
3000 ppm fluid at M = 45 take one to seven halvings; the 2000 ppm fluid at M = 90 takes
FOURTEEN on its very first step, with the residual falling four decades, jumping back up
to 1e46, and only then converging. The unrescuable episode that motivated this bound ran
twenty-three.

So the window between "still worth shrinking" and "never going to work" is 14 to 23, and
a bound of 12 sat inside the legitimate range: it stopped the 2000 ppm sweep dead at
t = 0 without a single step accepted. Eighteen leaves room above the largest observed
rescue while still cutting the spiral well before dt_min.

If a fluid ever needs more than this, the symptom is a march that reports no contact at
all, `cor = 1` and `tc = 0`, rather than a wrong answer.
"""
const MAX_ETA_HALVINGS = 18

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
    basis_kind::Symbol           # :legendre (well conditioned) or :monomial (historical)
    force_mode::Symbol           # :legendre (radial spectral field) or :nodal (conjugate loads)
    stop_on_release::Bool        # end the march when the drop leaves, rather than at t_max
end

"""
    ImpactParams(; We, Bo, Oh, M=DEFAULT_M, K=DEFAULT_K, ...)

`M` is the harmonic truncation: shape modes are `l = 2..M`, the film pressure carries
`l = 0..M`, and there are `M+1` collocation nodes. Those three counts are tied
together deliberately -- that is what makes the contact system square.

`eta_max_sweeps` bounds the Picard iteration on the viscosity, and it is headroom rather
than a tuning knob. How many sweeps a step actually needs depends strongly on resolution:
on the 3000 ppm fluid the worst step takes 17 sweeps at `M = 30` and 45 at `M = 60`, which
is why the cap is generous. The iteration is what a shear-thinning run spends its time on,
around 80 per cent of a step at `M = 60`, with the contact solve under one per cent.

`stop_on_release` ends the march once the drop has left the substrate and is rising, which
is all the impact metrics need and is why it defaults to true. Set it false to keep marching
to `t_max` and record the free flight after the bounce -- what the animations use, since a
video that stops at release does not show the drop leaving.
"""
function ImpactParams(; We, Bo, Oh, M::Int = DEFAULT_M, K::Int = DEFAULT_K, eta = gd -> 1.0,
                      dt0 = nothing, dt_min = 1e-10, t_max = 25.0,
                      eta_tol = nothing, eta_max_sweeps = 100, force_mode::Symbol = :legendre,
                      basis_kind::Symbol = :legendre, stop_on_release::Bool = true)
    ls = collect(2:M)
    tol = something(eta_tol, default_eta_tol(M))
    # theta = pi plus the zeros of P_M. These cluster at the poles, which is the
    # whole point: contact is resolved where contact happens.
    mus, _ = gauss_legendre_nodes(M, -1.0, 1.0)
    nodes = vcat(pi, acos.(clamp.(sort(mus), -1.0, 1.0)))
    dt = something(dt0, 2pi / (8 * sqrt(M * (M + 2) * (M - 1))))
    ec = all(gd -> eta(gd) == eta(0.0), (0.0, 1e-3, 1.0, 1e3))
    force_mode in (:legendre, :nodal) ||
        error("force_mode must be :legendre or :nodal, got $force_mode")
    basis_kind in (:legendre, :monomial) ||
        error("basis_kind must be :legendre or :monomial, got $basis_kind")
    ImpactParams(We, Bo, Oh, ls, K, eta, nodes, dt, dt_min, t_max, tol,
                 eta_max_sweeps, ec, basis_kind, force_mode, stop_on_release)
end

basis(p::ImpactParams) = ModalBasis(p.ls, p.K, p.basis_kind)
lmax(p::ImpactParams) = maximum(p.ls)
pc_len(p::ImpactParams) = lmax(p) + 1        # harmonics l = 0..M
pc_l(j::Int) = j - 1

"""Every trial function equals 1 at the surface, so the trace is a vector of ones."""
trace_vec(p::ImpactParams, l) = [phi(RitzBasis(l, p.K, p.basis_kind), k, 1.0) for k in 1:p.K]

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
            assemble_coupled(b, p.Oh; eta_rate = p.eta, state = adot_star) : F0
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

The contact set is found by a primal active set on the complementarity pair: grow while
a free node has penetrated, release while the outermost contacting node's pressure
pulls, stop when neither holds. Each move is forced by a violated condition, so there is
no candidate to rank and no tie to break.

This replaced a ranked search over neighbouring contact counts scored by an edge
residual, which is what the reference implementation does and what the nonvariational
solver still does. That search chatters when the damping is small: two candidates score
almost equally, the choice flips step to step, and the film pressure oscillates in sign
until `dt` collapses, which made it fail sporadically in `M`. The active set also
matches the reference case more closely, giving a contact time of 2.1830 against its
2.183 where the search gave 2.300.
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
        if p.stop_on_release && best_cp == 0 && curr.v > 0 && curr.z > 1.0 && any(>(0), cps)
            break
        end
    end
    traj = (t = ts, z = zs, v = vs, cp = cps, pc1 = pc1, a = as, adot = adots, pc = pcs,
            rejects = nrej,
            ## the Picard iteration's worst behaviour over the whole march, so a regime
            ## where the map stops contracting announces itself instead of returning quietly
            eta_sweeps_max = max_sweeps, eta_resid_max = max_resid)
    merge(traj, _kpis(p, traj))
end

"""
Coefficient of restitution: rebound speed over impact speed.

Measured at the FRAME AFTER THE LAST CONTACT, not at the end of the record. The two are
the same whenever the march stops at release, which is the default and is why this went
unnoticed: the break fires on the very step that clears the substrate, so the last frame
IS the release frame. They stop being the same the moment a run keeps integrating -- with
`stop_on_release = false` the drop is still rising and still decelerating under gravity,
and reading `vs[end]` after five more capillary times reported 0.602 for a bounce whose
restitution is 0.750. Anchoring to the release frame gives the same number either way.
"""
function restitution(vs, cps, We)
    inc = findfirst(>(0), cps)
    inc === nothing && return NaN
    lastc = findlast(>(0), cps)
    lastc >= length(vs) && return NaN
    abs(vs[lastc + 1] / vs[max(inc - 1, 1)])
end

"""
Contact time: from first touch to final separation.

The SPAN, not the longest contiguous run of nonempty contact. Those differ, and the
difference is not cosmetic. When the contact set is solved for rather than assumed to be a
disc it comes out annular, and its node count dips briefly to zero as the dimple evolves.
Nine such dips in a run of 2516 steps split one physical contact into ten runs and reported
a contact time ten per cent short -- which I then attributed to a ring releasing sooner
than a disc. It was arithmetic, not physics: measured as a span, that same run gives 2.1830
against the reference implementation's 2.183.

The span is also what an experiment measures. A camera sees the drop down from first touch
to final separation and cannot resolve a separation lasting a few timesteps, so merging
them is the faithful comparison rather than a convenience. `contact_gap_fraction` is
reported alongside so a span that has quietly merged two genuine bounces cannot pass as one
contact.
"""
function contact_time(ts, cps)
    inc = findfirst(>(0), cps)
    inc === nothing && return 0.0
    ## FIRST TOUCH TO LAST RELEASE, unconditionally. Contact is contact: as long as one
    ## point is touching, the drop is in contact, and intervening separations do not start
    ## a new measurement. An earlier version of mine reported only the first episode,
    ## merging separations below a threshold -- which is a definition I invented, and it
    ## understated the contact time by up to a factor of three on runs where the drop
    ## genuinely bounces mid-contact.
    ts[findlast(>(0), cps)] - ts[inc]
end

"""
Fraction of the contact span during which the contact set was momentarily empty.

Small means one contact with brief dimple transients. Large means the span has merged
events that are genuinely separate, and should not be read as a single contact time.
"""
function contact_gap_fraction(ts, cps)
    inc = findfirst(>(0), cps)
    inc === nothing && return 0.0
    lastc = findlast(>(0), cps)
    span = ts[lastc] - ts[inc]
    span <= 0 && return 0.0
    empty_t = 0.0
    for i in inc:(lastc-1)
        cps[i] == 0 && (empty_t += ts[i+1] - ts[i])
    end
    empty_t / span
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
#   * The affine relation itself was checked against the searching solver's own output:
#     relative error 0.
#   * `A_c` is NOT symmetric. It is asymmetric by about forty per cent at every M and K
#     tried, so the LCP is not the KKT system of a convex programme and a projected
#     Gauss-Seidel sweep does not solve it. `lcp_active_set` does, without assuming
#     symmetry, and the residual it reports is measured against this map rather than a
#     symmetrised copy of it.
#
# TWO EARLIER CLAIMS HERE WERE WRONG, and they cost real results, so they are recorded
# rather than deleted:
#
#   * "`A_c` is symmetric to 1e-16 in the NODAL pressure basis; the change of variables is
#     what buys the structure." The measurement behind this was of `H A^-1 H'`, which is
#     symmetric for ANY `H` whenever `A` is symmetric -- it could not have come out
#     otherwise, so it tested nothing. The matrix actually assembled is `H A^-1 Q_n`, gap
#     Jacobian against force operator, and those are not transposes of each other. Moving to
#     nodal variables WAS necessary -- `p >= 0` is not a sign condition on Legendre
#     coefficients -- but necessary is not sufficient, and symmetry needs `Q_n = -H' W`.
#   * "the residual `max_i min(h_i, p_i)` is a certificate." It is small whenever one of
#     each pair is small, so it certifies penetration. Combined with the symmetrised
#     surrogate it passed states with the drop 4.7 per cent of a radius inside the wall, and
#     those states are where the annular contact and the closure-dependent contact time both
#     came from. Both retracted; see `test/test_lcp_contact.jl`.
#
# And what it is worth: the searching solver satisfies `h >= 0` and `h p = 0` by
# construction but NOT `p >= 0`. Measured over whole runs, the pressure is negative at
# 0.015 per cent of node-steps in a case that works and 41 per cent in one that does not.
# The dropped condition is exactly the difference between the two.

"""
    lcp_pgs(A, b; iters, tol) -> (p, residual, sweeps)

Projected Gauss-Seidel for `h = Ap + b`, `h >= 0`, `p >= 0`, `p'h = 0`, requiring `A`
symmetric positive semi-definite.

NOT USED BY THE VARIATIONAL CONTACT SOLVER, whose compliance is asymmetric -- see
`contact_lcp`. Retained because the requirement is worth stating explicitly and because a
symmetric problem is cheaper to sweep than to pivot. Use [`lcp_active_set`](@ref) unless you
have checked that your `A` is symmetric.

Each sweep is one pass of Gauss-Seidel followed by projection onto `p >= 0`, which for a
symmetric PSD matrix is coordinate descent on the equivalent convex programme and so
cannot increase the objective. Nodes with a non-positive diagonal are skipped: those are
nodes whose gap does not respond to their own pressure, which cannot carry load.

The returned residual is TWO-SIDED: it measures violation of `h >= 0` and of `p >= 0`
separately as well as complementarity, so it vanishes only for an actual solution. See
[`lcp_residual`](@ref).

THIS USED TO BE `max_i min(h_i, p_i)`, AND THAT WAS WRONG. `min(h_i, p_i)` is small whenever
EITHER entry is small, so a node with a gap of -0.012 and a pressure of 3 scores -0.012 and
the maximum over nodes returns something negligible. The quantity measures "are both positive
at once" and never "is each one positive at all". Used as the stopping test it ended the
iteration after a single sweep on every step of every run -- the one-sided residual was
already at 1e-13 by then -- and the accepted states had the drop penetrating the substrate by
up to 1.2 per cent of a radius, which is comparable to the 0.02R threshold the contact metrics
are measured at. The iteration was never converging; it was being told it had.

Iterating properly is free. The complementarity solve is 0.1 to 0.6 per cent of the cost of a
step at M from 20 to 90, where assembling the compliance is the rest, so there was no
performance reason for the early exit either.
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
        resid = lcp_residual(A, b, p)
        resid < tol && break
    end
    (p, resid, sweeps)
end

"""
    lcp_active_set(A, b; tol, maxpivots) -> (p, residual, pivots)

Solve `h = Ap + b`, `h >= 0`, `p >= 0`, `p_i h_i = 0` by active-set pivoting, WITHOUT assuming
`A` is symmetric.

This is the solver the variational contact problem actually needs. [`lcp_pgs`](@ref) is
Gauss-Seidel on the convex programme `min ½pᵀAp + bᵀp`, which is equivalent to the
complementarity problem only when `A` is symmetric; the assembled compliance is asymmetric by
about forty per cent, for the reason given in `contact_lcp`. Sweeping it anyway converges to
something that is not a solution, and the symmetrised surrogate it was previously handed
admitted states with the drop inside the substrate.

Guess a contact set, solve the equality problem on it, then repair the guess: drop any node
that came out pulling, add any node the surface has passed through. Each move is forced, so
the iteration cannot cycle between two sets for the same reason the searching closure's
cannot.
"""
function lcp_active_set(A::AbstractMatrix, b::AbstractVector;
                        tol::Real = 1e-12, maxpivots::Int = 400)
    n = length(b)
    active = falses(n)
    p = zeros(n)
    pivots = 0
    for _ in 1:maxpivots
        pivots += 1
        p .= 0.0
        idx = findall(active)
        if !isempty(idx)
            sub = A[idx, idx]
            ps = try
                sub \ (-b[idx])
            catch
                pinv(Matrix(sub)) * (-b[idx])       # a rank-deficient set is not fatal
            end
            p[idx] .= ps
        end
        sc = max(maximum(abs, b), maximum(abs, p), 1.0)
        ## a node that came out pulling cannot be carrying load -- release the worst
        if !isempty(idx)
            j = argmin(p[idx])
            if p[idx[j]] < -tol*sc
                active[idx[j]] = false
                continue
            end
        end
        h = A*p .+ b
        ## a free node the surface has passed through must be in contact -- seize the worst
        free = findall(.!active)
        if !isempty(free)
            j = argmin(h[free])
            if h[free[j]] < -tol*sc
                active[free[j]] = true
                continue
            end
        end
        break
    end
    (p, lcp_residual(A, b, p), pivots)
end

"""
    lcp_residual(A, b, p) -> Float64

How badly `p` fails to solve `h = Ap + b`, `h >= 0`, `p >= 0`, `p_i h_i = 0`, scaled and
counting all three conditions. Zero only for a solution.

Each condition is measured on its own, which is the whole point: a measure that combines them
can be small while one of them is badly violated, and the combination `max_i min(h_i, p_i)`
does exactly that.
"""
function lcp_residual(A::AbstractMatrix, b::AbstractVector, p::AbstractVector)
    h = A*p .+ b
    sc = max(maximum(abs, h), maximum(abs, p), 1.0)
    max(maximum(max.(-h, 0.0)),          # the drop must not enter the substrate
        maximum(max.(-p, 0.0)),          # the substrate must not pull
        maximum(abs.(h .* p))) / sc      # pressure only where there is contact
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
    ## h is affine in the contact unknown through BOTH the shape and the centre of mass
    Ac_full = H * (A \ Qn) .+ dz .* (ones(nn) * Vinv[2, :]')
    b_full  = H * (A \ rhs0) .+ mu .+ z0
    ## Lower hemisphere only, selected with a TOLERANCE rather than by the sign of cos(theta).
    ##
    ## At odd M the polynomial P_M has a root at mu = 0, so one collocation node lands exactly
    ## on the equator -- and there the gap row is identically zero, because every entry carries
    ## a factor cos(theta). Such a node cannot be in contact and its pressure does nothing, but
    ## if it is admitted it puts a zero row and column into the compliance and the matrix is
    ## singular. Whether `cos(theta) < 0` admitted it was decided by rounding: at M = 45 and
    ## M = 61 the root evaluates to about -3e-16 and the node was kept, at M = 21, 31 and 91 it
    ## evaluates positive and the node was dropped. M = 45 is the truncation the validation
    ## sweeps use, which is where the otherwise unexplained 1e-35 smallest eigenvalue came from.
    EQUATOR_TOL = 1e-8
    idx = [i for i in 1:nn if cos(p.nodes[i]) < -EQUATOR_TOL]
    ## THE COMPLIANCE IS RETURNED AS ASSEMBLED, NOT SYMMETRISED.
    ##
    ## It is asymmetric by about forty per cent, and this used to be forced symmetric with
    ## `Symmetric(0.5(S + S'))` before being handed to the projected Gauss-Seidel sweep. That
    ## made the sweep solve a DIFFERENT problem from the one this map defines, and the residual
    ## -- computed against the symmetrised surrogate -- reported 1e-11 while the state it
    ## returned had the drop 4.7 per cent of a radius inside the substrate.
    ##
    ## WHY IT IS ASYMMETRIC, since the obvious guesses are wrong. It is not a missing cos(theta)
    ## between a vertical gap and a radial pressure: that would predict perfect alignment at the
    ## pole, where the two directions coincide, and the misalignment is WORST there (cosine
    ## similarity between the force column and its own gap row is 0.23 at the pole and 0.99 at
    ## the equator). The cause is that the film pressure is carried as a degree-M Legendre FIELD
    ## with forcing Q_l = -(4pi/(2l+1)) p_{c,l}, so a unit nodal pressure enters as the Galerkin
    ## force of the polynomial interpolating a delta at that node -- an oscillation spread over
    ## the whole sphere rather than a load at one place. A complementarity problem needs the
    ## constraint and its multiplier to be dual, and interpolation does not make them dual.
    ##
    ## So the contact problem here is genuinely NOT the convex programme a symmetric compliance
    ## would give. It is still a well-posed complementarity problem, and `lcp_active_set` solves
    ## it without assuming symmetry.
    ##
    ## THE CONJUGATE ALTERNATIVE, selected by `force_mode = :nodal`. Take the contact unknown to
    ## be the vertical LOAD at each node rather than a pressure field. Differentiating the
    ## constraint then gives the forcing with no freedom left in it,
    ##
    ##     A xi = f + H' lam,      m zddot = -m Bo + 1' lam,
    ##
    ## because H' and 1' ARE the derivatives of h with respect to xi and z. Eliminating both,
    ##
    ##     h = W lam + b,   W = H A^-1 H' + (1/(m beta^2)) 1 1',
    ##
    ## and every term is symmetric positive semi-definite: the first because A is, the second
    ## because it is a positive multiple of an outer square. The problem is then exactly the KKT
    ## system of min (1/2) lam' W lam + b' lam over lam >= 0 -- a convex programme, with a
    ## solution always and a unique one when W is definite, which it is once the degenerate
    ## equator node is excluded.
    ##
    ## The price is that the multiplier is a load and not a pointwise pressure. The pressure is
    ## recoverable only as a diagnostic, and not at the pole at all, because the quadrature
    ## weight there is zero to machine precision.
    if p.force_mode === :nodal
        Hl = H[idx, :]
        m_drop = 4pi/3
        W = Hl * (A \ transpose(Hl)) .+ (1/(m_drop*β^2)) .* (ones(length(idx)) * ones(length(idx))')
        return (W, b_full[idx], idx,
                (A = A, rhs0 = rhs0, Qn = Qn, Vinv = Vinv, β = β, hv_a = hv_a,
                 hv_z = hv_z, hv_v = hv_v, npc = npc, N = N, H = H, Hl = Hl,
                 m_drop = m_drop))
    end
    (Ac_full[idx, idx], b_full[idx], idx,
     (A = A, rhs0 = rhs0, Qn = Qn, Vinv = Vinv, β = β, hv_a = hv_a,
      hv_z = hv_z, hv_v = hv_v, npc = npc, N = N, H = H, Hl = H[idx, :],
      m_drop = 4pi/3))
end

"""
    shape_from_contact(p, aux, idx, x) -> Vector

Interior amplitudes implied by a contact solution.

In `:legendre` mode `x` is the nodal pressure and enters through `Q_n`. In `:nodal` mode it is
the vector of vertical loads and enters through `H'` -- the transpose of the very constraint it
is dual to, which is what makes the compliance symmetric.
"""
function shape_from_contact(p::ImpactParams, aux, idx, x::AbstractVector)
    if p.force_mode === :nodal
        return aux.A \ (aux.rhs0 + transpose(aux.Hl) * x)
    end
    pfull = zeros(aux.npc); pfull[idx] = x
    aux.A \ (aux.rhs0 + aux.Qn * pfull)
end

"""
    equivalent_pressure(p, aux, idx, x) -> Vector

Film-pressure harmonics `l = 0..M`, for reporting.

In `:legendre` mode these are the solved unknowns. In `:nodal` mode they are a DIAGNOSTIC and
nothing more: the unknowns there are loads, and what is reported is the harmonic content that
would produce the same generalised force,

    p_{c,l} = -(2l+1)/(4 pi) * (Q . tv_l) / (tv_l . tv_l),

together with `p_{c,1} = -sum(lam)/m`, which makes the net vertical force exact rather than
fitted. A pointwise pressure is not recoverable at the pole in this mode at all, because the
quadrature weight there is zero to machine precision -- the load is finite and the area it acts
over is not resolved.
"""
function equivalent_pressure(p::ImpactParams, aux, idx, x::AbstractVector)
    if p.force_mode !== :nodal
        pfull = zeros(aux.npc); pfull[idx] = x
        return aux.Vinv * pfull
    end
    b = basis(p)
    Q = transpose(aux.Hl) * x
    pc = zeros(aux.npc)
    for (i, l) in enumerate(b.ls)
        tv = trace_vec(p, l)
        rng = dofindex(b, i, 1):dofindex(b, i, p.K)
        pc[l + 1] = -(2l + 1)/(4pi) * dot(view(Q, rng), tv) / dot(tv, tv)
    end
    pc[2] = -sum(x) / aux.m_drop            # exact net vertical force
    pc
end

"""
    try_step_lcp(p, prev, curr, dt; F0, Vfac) -> (status, next, diag)

One step with the contact set determined by complementarity rather than searched for.

There is no contact-count argument, and that is the point: the set is an output. The
pressure comes from the LCP, the interior amplitudes and the centre of mass follow from it
by back-substitution, and no candidate is ever proposed, ranked or rejected.

`status` is `:ok`, or `:diverge` if the complementarity residual does not clear tolerance
-- in which case the caller should halve `dt`, as for any nonlinear solve that fails.
"""
function try_step_lcp(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                      dt::Float64; F0 = nothing, Vfac, tol::Real = 1e-9)
    ## SHEAR THINNING MAKES THE COMPLEMENTARITY PROBLEM NONLINEAR.
    ##
    ## eta depends on the strain rate, the dissipation operator depends on eta, and the
    ## compliance A_c = H A^-1 Q depends on that operator -- so with a variable viscosity
    ## A_c is a function of the very velocity the solve produces. The LCP is no longer
    ## linear, and it is closed the same way the searching solver closes its own
    ## nonlinearity: freeze eta at an extrapolated strain rate, solve the resulting LINEAR
    ## complementarity problem, re-evaluate eta at the answer, repeat. Each iterate is a
    ## genuine LCP -- asymmetric, as in the constant-viscosity case -- so complementarity
    ## holds exactly at
    ## every sweep and only eta lags.
    ##
    ## The cost is real: a variable viscosity needs a fresh coupled assembly AND a fresh
    ## compliance -- M+1 back-substitutions -- on every sweep, where the constant-viscosity
    ## path reuses both.
    adot_star = curr.first ? curr.adot : (1 + dt/curr.dt) * curr.adot - (dt/curr.dt) * prev.adot
    conv = 0.0; used = 0
    star_back = Float64[]                      # iterate before last, for the secant step
    Ac = zeros(0,0); bv = Float64[]; idx = Int[]; aux = nothing
    pact = Float64[]; resid = Inf; sweeps = 0
    for it in 1:(F0 === nothing ? p.eta_max_sweeps : 1)
        used = it
        prev_star = copy(adot_star)
        F = F0 === nothing ?
            assemble_coupled(basis(p), p.Oh; eta_rate = p.eta, state = adot_star) : F0
        Ac, bv, idx, aux = contact_lcp(p, prev, curr, dt; F0 = F, Vfac = Vfac)
        pact, resid, sweeps = lcp_active_set(Ac, bv)
        ## `resid` is already the two-sided measure, so this rejects a step whose gap or
        ## pressure has the wrong sign as well as one that is not complementary
        resid > tol && return (:diverge, curr,
            (resid = resid, sweeps = sweeps, nact = 0, contiguous = true,
             eta_sweeps = used, eta_resid = conv))
        a_try = shape_from_contact(p, aux, idx, pact)
        adot_star = aux.β * a_try + aux.hv_a
        if F0 !== nothing
            conv = 0.0
            break
        end
        ## The residual is measured on the RAW fixed-point step, before any
        ## acceleration, so it stays the honest distance between successive iterates.
        scv = max(maximum(abs, adot_star), 1e-12)
        conv = maximum(abs, adot_star .- prev_star) / scv
        conv < p.eta_tol && break
        ## Irons-Tuck (Anderson with depth one) on the NEXT guess. Plain Picard on this
        ## fixed point converges linearly, and the rate worsens with resolution: the
        ## worst step needs 20 sweeps at M = 45, 16 at M = 60 and 63 at M = 90, where it
        ## approaches the cap and the march dies. A secant step along the last two
        ## increments costs three vector operations on `ndof` and does not change the
        ## fixed point, only the path taken to it.
        if !isempty(star_back)
            d1 = prev_star .- star_back
            d2 = adot_star .- prev_star
            dd = d2 .- d1
            ndd = dot(dd, dd)
            if ndd > 0
                adot_star = adot_star .- (dot(d2, dd) / ndd) .* d2
            end
        end
        star_back = prev_star
    end
    (F0 === nothing && conv > p.eta_tol) && return (:diverge, curr,
        (resid = resid, sweeps = sweeps, nact = 0, contiguous = true,
         eta_sweeps = used, eta_resid = conv))

    β = aux.β
    a_next = shape_from_contact(p, aux, idx, pact)
    pc = equivalent_pressure(p, aux, idx, pact)
    ## Centre-of-mass acceleration is -Bo plus the film force per unit mass, which is -p_{c,1}
    ## for a spectral pressure and sum(lam)/m for nodal loads. `equivalent_pressure` puts the
    ## second into the same slot, so this line is identical in both modes.
    v_next = (-p.Bo - pc[2] - aux.hv_v) / β
    z_next = (v_next - aux.hv_z) / β
    adot_next = β * a_next + aux.hv_a
    ## The active set is whatever the solve says it is. It is NOT constrained to be an
    ## interval, so whether it comes out contiguous is a RESULT rather than an assumption.
    act = [i for i in eachindex(pact) if pact[i] > 1e-10 * max(maximum(pact), 1.0)]
    nact = length(act)
    contiguous = isempty(act) || (act == collect(first(act):last(act)) && first(act) == 1)
    nxt = ImpactState(curr.t + dt, dt, a_next, adot_next, z_next, v_next, nact, pc, false)
    (:ok, nxt, (resid = resid, sweeps = sweeps, nact = nact, contiguous = contiguous,
                eta_sweeps = used, eta_resid = conv))
end

"""
    simulate_lcp(p) -> NamedTuple

March the impact with the contact set determined by complementarity at every step.

No candidate search, no tangency ranking, no `±2` probe, no contact-count limit. `dt`
halves only when the complementarity solve itself fails to converge.
"""
function simulate_lcp(p::ImpactParams)
    ## A constant viscosity lets the operator be built once; a variable one is rebuilt
    ## inside the Picard loop, which is what makes the shear-thinning path expensive.
    F0 = p.eta_const ? assemble_newtonian(basis(p), p.Oh) : nothing
    Vfac = lu(legendre_vandermonde(p))
    s0 = initial_state(p)
    prev, curr = s0, s0
    ts = Float64[0.0]; zs = Float64[s0.z]; vs = Float64[s0.v]
    cps = Int[0]; pc1 = Float64[0.0]
    as = Vector{Float64}[copy(s0.a)]; adots = Vector{Float64}[copy(s0.adot)]
    pcs = Vector{Float64}[copy(s0.pc)]
    dt = p.dt0; nrej = 0; worst_resid = 0.0; max_sweeps = 0; noncontig = 0
    eta_sw = 0; eta_rs = 0.0
    eta_halvings = 0       # consecutive step reductions spent on a viscosity failure
    while curr.t < p.t_max
        st, nxt, dg = try_step_lcp(p, prev, curr, dt; F0 = F0, Vfac = Vfac)
        if st !== :ok
            nrej += 1
            ## HALVING dt RESCUES SOME VISCOSITY FAILURES AND IS FUTILE FOR OTHERS,
            ## and the residual does not say which. Measured on the 3000 ppm fluid,
            ## episodes that halving DID rescue contain residuals of 1.2e-4, 4.4e-4,
            ## 8.2e-1, 3.8e0 and 1.4e13; the episode it could not rescue sat at 3.0e-6.
            ## The same magnitudes appear on both sides, so every threshold on
            ## `eta_resid` classifies one of them wrongly. Four versions of this test
            ## did exactly that.
            ##
            ## What separates them is how MANY halvings it takes. Rescued episodes end
            ## after one to seven; the fatal one ran twenty-three consecutive halvings
            ## with the residual climbing as 1/dt the whole way, because
            ## `adot_star = beta*a + hv_a` with `beta = c0/dt` is a difference of large
            ## quantities and shrinking the step raises its own noise floor.
            ##
            ## So the step is allowed a bounded number of reductions. Beyond that the
            ## iteration is not going to be rescued by a smaller step and the march
            ## says so, instead of grinding to dt_min and reporting part of a bounce as
            ## if it were a result.
            if dg.eta_resid > p.eta_tol && dg.resid <= 1e-8
                eta_halvings += 1
                if eta_halvings > MAX_ETA_HALVINGS
                    @warn "the viscosity iteration is not recovering as the step" *
                          " shrinks; no smaller step will help, so the march stops" *
                          " here (raise eta_tol, or lower M)" t=curr.t dt=dt*1.0 *
                          1.0 eta_resid=dg.eta_resid halvings=eta_halvings
                    break
                end
            end
            dt /= 2
            dt < p.dt_min && break
            continue
        end
        eta_halvings = 0
        worst_resid = max(worst_resid, abs(dg.resid))
        max_sweeps = max(max_sweeps, dg.sweeps)
        eta_sw = max(eta_sw, dg.eta_sweeps); eta_rs = max(eta_rs, dg.eta_resid)
        dg.contiguous || (noncontig += 1)
        prev, curr = curr, nxt
        push!(ts, curr.t); push!(zs, curr.z); push!(vs, curr.v)
        push!(cps, dg.nact); push!(pc1, curr.pc[2])
        push!(as, copy(curr.a)); push!(adots, copy(curr.adot)); push!(pcs, copy(curr.pc))
        dt = min(2*dt, p.dt0)
        if p.stop_on_release && dg.nact == 0 && curr.v > 0 && curr.z > 1.0 && any(>(0), cps)
            break
        end
    end
    traj = (t = ts, z = zs, v = vs, cp = cps, pc1 = pc1, a = as, adot = adots, pc = pcs,
            rejects = nrej, lcp_resid_max = worst_resid, lcp_sweeps_max = max_sweeps,
            eta_sweeps_max = eta_sw, eta_resid_max = eta_rs,
            noncontiguous_steps = noncontig)
    merge(traj, _kpis(p, traj))
end

# ============================================================================
# EXPERIMENT-LIKE CONTACT METRICS
#
# A camera does not see "the pressure is positive at a collocation node". It sees the drop
# arrive at the wall within the resolution of the image, and it sees it leave. So contact is
# defined here the way the measurement defines it: contact exists whenever ANY point of the
# surface is below a threshold height above the substrate, contact time runs from the first
# such frame to the last, and restitution compares the centre-of-mass speed at those two
# instants.
#
# The threshold matters for more than fidelity. Judging contact by the pressure being
# nonzero makes a separation of 0.001R end the contact, which fragments one physical
# impact into several and made a contact time read three times short. At a threshold of
# 0.02R those separations are simply not resolved -- which is exactly what a camera does --
# and the metric stops depending on the solver's internal notion of a contact set.

"""
    _kpis(p, traj) -> (cor, tc, gap_fraction)

The KPIs of a finished march, from the one definition of contact this package has.

`proximity_metrics` is that definition. It is called here rather than duplicated, so that
`simulate`, `simulate_lcp` and every script reading `r.cor` report the same quantity as the
figures and the validation. Two implementations existed before this and they disagreed: the
solver returned a bare velocity ratio keyed on the contact-node count, while the metrics
used an energy ratio keyed on a 0.02R proximity line. Numbers taken from one and compared
against the other differed by ten per cent at Bond numbers where gravity matters, and by
much more near the rebound threshold.
"""
function _kpis(p::ImpactParams, traj)
    m = proximity_metrics(p, traj)
    (cor = m.cor, tc = m.tc, gap_fraction = contact_gap_fraction(traj.t, traj.cp))
end

"""
    min_gap_series(p, r; nang = 240) -> Vector

Smallest gap between the drop surface and the substrate at each stored step.

Separated out because it is the only part of `proximity_metrics` that touches the
shape, it is by far the most expensive part, and **it does not depend on the
threshold**. Caching it lets a run be re-scored at any `h_thresh`, and under any
later change of restitution convention, without simulating anything again.
"""
function min_gap_series(p::ImpactParams, r; nang::Int = 240)
    ths = range(pi, pi/2; length = nang)          # the lower hemisphere, where contact is
    mus = cos.(ths)
    ## precompute the Legendre matrix once: the per-step cost is then one matrix-vector
    ## product rather than tens of thousands of polynomial evaluations
    P = [legendre_angular(l, mu).P for mu in mus, l in p.ls]
    minh = fill(Inf, length(r.t))
    for i in eachindex(r.t)
        z = surface_amplitudes(p, r.a[i])
        rad = 1 .+ P * z                          # surface radius at each sampled angle
        h = r.z[i] .+ mus .* rad                  # height above the substrate
        minh[i] = minimum(h)
    end
    minh
end

"""
    proximity_metrics(p, r; h_thresh = 0.02, nang = 240) -> NamedTuple

Contact time and restitution measured by proximity, the way an experiment measures them.

Contact exists at a step when the minimum gap over the whole surface is below `h_thresh`
(in units of the drop radius). `tc` is the interval between the first and last such step.

RESTITUTION IS AN ENERGY RATIO, AND THE GRAVITY TERM IS NOT OPTIONAL.

    cor  = sqrt(|E_out / E_in|),   E_in  = v_in^2 / 2,
                                   E_out = v_out^2 / 2 + (z_out - z_in) * Bo

`cor_vel = |v_out / v_in|` is also returned, and is what this function used to report.
The two agree whenever the centre of mass leaves at the height it arrived at, and the
second is wrong whenever it does not, because the work gravity does on the drop while it
is squashed is then counted as restitution.

How wrong is not academic. At `Bo = 0.29` -- a millimetre drop, the largest in the
Thenarianto et al. (2023) set -- and `We = 3e-3`, the velocity ratio returns 3.9. A
restitution above one is energy from nowhere: the drop is detected on the way in while it
is barely moving, and gravity accelerates it over the remaining `h_thresh` of fall.

The energy form is what the ancestor of this package publishes (Gabbard et al. 2025,
`coef_restitution_exp`), and it is what the experiment measures: their protocol fits
parabolas to the free-flight centre-of-mass track and evaluates the velocity at the
contact instants, which is gravity-corrected by construction.

WHERE CONTACT IS DECLARED IS A SEPARATE CHOICE, and the energy form does not rescue a bad
one. At `Bo = 0.29` the speed picked up falling through `0.02R` is `sqrt(2 Bo h_thresh)`,
twice the impact speed at that Weber number, and the energy ratio still returns 3.8. Pass
`h_thresh = 0` to declare contact where the drop actually touches, which is what
Thenarianto et al. measure and what stays finite there.

The surface is sampled on `nang` angles rather than at the collocation nodes, because the
criterion is about the surface and not about the discretisation -- "a single point below the
line", not "a node below the line".
"""
function proximity_metrics(p::ImpactParams, r; h_thresh::Real = 0.02, nang::Int = 240,
                           minh::Union{Nothing,AbstractVector} = nothing)
    minh = minh === nothing ? min_gap_series(p, r; nang = nang) : minh
    detected = minh .< h_thresh
    i1 = findfirst(detected)
    if i1 === nothing
        return (tc = 0.0, cor = NaN, cor_vel = NaN, i_first = 0, i_last = 0,
                v_in = NaN, v_out = NaN, min_gap = minimum(minh), n_detected = 0,
                released = false, h_thresh = h_thresh)
    end

    ## LIFTOFF IS AN EVENT, NOT THE END OF THE MARCH.
    ##
    ## The outgoing state is read where the drop rises back above the line, which is
    ## the last index at which it is still detected. A drop that never leaves has no
    ## such index and no outgoing state: there is no restitution to report, and the
    ## caller is told `released = false` rather than handed a number.
    ##
    ## Taking `findlast(detected)` instead evaluates the outgoing energy at whatever
    ## step the march happened to stop on. For a settled drop that is `t_max`, by which
    ## time the centre of mass has sunk well below where it touched, so the
    ## gravitational term makes E_out NEGATIVE and `sqrt(abs(...))` turns a drop that
    ## never bounced into a restitution of six. This is what the ancestor avoids by
    ## computing `Eout` only inside the branch that detects liftoff.
    ## Two ways the drop can stop being detected, and only one of them is a liftoff.
    ## `stop_on_release` ends the march at release, so a bouncing drop has no
    ## un-detected step at all -- the trajectory simply stops. A settled drop also has
    ## none, because it is still down when the clock runs out. The two are told apart
    ## by whether the clock ran out, not by the absence of a later step.
    ## Liftoff ends the FIRST contact: the first stored step at which the surface is back
    ## above the line. Not the last step in contact anywhere in the march -- a drop that
    ## bounces, falls back under gravity and lands again would then have its outgoing
    ## state read from the wrong impact, and a march long enough to catch the second
    ## landing would report the first bounce as a drop that never left.
    ##
    ## The gap does oscillate while the drop is down, by about 1e-3 of a radius, so "first
    ## step out of contact" would be unreliable against a threshold of zero. Against 0.02R
    ## it is not: the whole oscillation sits far below the line, and the first crossing is
    ## the departure. This is the same reason the threshold is 0.02R and not "the pressure
    ## is nonzero".
    ##
    ## No crossing at all means the drop never got back above the line, which is a drop
    ## that stayed down: no outgoing state, and no restitution to report.
    ## Contact ends at the last step the drop is still down. After that it is in free
    ## flight, and the energy referenced to the measurement line is conserved there, so it
    ## does not matter whether the march was integrated all the way back up to the line or
    ## stopped at release: `E_out` is the same either way. Whether the drop CLEARS the line
    ## is then not a separate question -- it clears it exactly when that energy is
    ## positive, since reaching the line is what `E_out = 0` means.
    ##
    ## This is the same back-extrapolation the ancestor applies on the way in, and it is
    ## what lets a truncated trajectory be scored without integrating the flight.
    ## Contact ends at the last step of the FIRST contact episode. Not the last contact in
    ## the march: with a long enough clock the drop bounces, falls back and lands again,
    ## and `findlast` would read the outgoing state off the wrong impact.
    ##
    ## Three ways the episode can end, and they must not be confused:
    ##   a step out of contact exists   -> that is the liftoff
    ##   none, and the clock ran out    -> the drop settled and never left
    ##   none, and the clock did not    -> `stop_on_release` ended the march AT liftoff
    j = findfirst(k -> !detected[k], (i1 + 1):length(detected))
    i2 = if j !== nothing
        i1 + j - 1
    elseif r.t[end] >= 0.9 * p.t_max
        0
    else
        length(detected)
    end
    if i2 == 0
        return (tc = r.t[end] - r.t[i1], cor = NaN, cor_vel = NaN,
                i_first = i1, i_last = 0, v_in = r.v[i1], v_out = NaN,
                min_gap = minimum(minh), n_detected = count(detected),
                released = false, h_thresh = h_thresh)
    end

    ## The drop leaves from a different height than it arrived at, and gravity did work
    ## over that displacement. Counting it as restitution is what makes the velocity
    ## ratio exceed one.
    ## BOTH ENERGIES ARE REFERENCED TO THE MEASUREMENT LINE, NOT TO THE CONTACT INSTANT.
    ##
    ## The line sits `h_thresh` above the substrate. The drop crosses it on the way in,
    ## keeps falling to the substrate, and crosses it again on the way out, so the natural
    ## reference for an energy budget is the line itself -- and that is what the ancestor
    ## uses (`Vin_exp = Vin + t0*g`, `CM_in_exp = (1+pixel_adim)*CM_in`).
    ##
    ##   v at the line   v_line^2 = v_contact^2 - 2 Bo h        (free fall over h)
    ##   height          z_ref    = z_contact (1 + h)
    ##
    ## Referencing to the contact instant instead -- which this function did -- inflates
    ## E_in and, far worse, removes a cancellation in E_out. Near the roll-off E_out is a
    ## difference of two comparable terms, the outgoing kinetic energy against the
    ## potential energy still owed on the climb back to the line, and it collapses toward
    ## zero. That collapse IS the roll-off: it is what takes restitution from 0.53 to 0
    ## over a factor of 1.5 in Weber number. Measured from z_contact the two terms never
    ## cancel, restitution sits flat near 0.9, and the roll-off does not exist at all.
    ##
    ## Below `v_contact^2 = 2 Bo h` there is no crossing to reference to: the drop never
    ## reached the line in free flight, `E_in` is negative, and the case is one the
    ## ancestor skips outright. `is_measurable` reports that.
    v_line_sq = r.v[i1]^2 - 2 * p.Bo * h_thresh
    z_ref     = r.z[i1] * (1 + h_thresh)
    E_in  = v_line_sq / 2
    E_out = r.v[i2]^2 / 2 + (r.z[i2] - z_ref) * p.Bo
    ## The drop reaches the line only if the energy referenced to it is positive.
    released = E_out > 0

    ## CONTACT TIME RUNS TO THE LINE, NOT TO WHEREVER THE MARCH STOPPED.
    ##
    ## Energy is conserved in free flight, so `cor` does not care whether the trajectory
    ## was integrated back up to the line or cut off at release. Time is not conserved, and
    ## `tc` does: a march ended by `stop_on_release` stops while the surface is still below
    ## the line, and reports a contact time short by the climb. That made `tc` depend on a
    ## solver option, which it must not.
    ##
    ## The remaining climb is free flight, so it is closed-form. Taking the drop as rigid
    ## over it, the gap rises with the centre of mass, and the first crossing is
    ##
    ##     dz = h_thresh - gap(i2),   v t - Bo t^2 / 2 = dz
    ##     t  = (v - sqrt(v^2 - 2 Bo dz)) / Bo        (the smaller root)
    ##
    ## Added only when the march was truncated; when the trajectory already contains the
    ## crossing, `i2` is that step and there is nothing to add.
    t_climb = 0.0
    if i2 == length(detected) && released
        dz = h_thresh - minh[i2]
        v2 = r.v[i2]
        if dz > 0 && v2 > 0
            disc = v2^2 - 2 * p.Bo * dz
            t_climb = disc > 0 ? (p.Bo > 0 ? (v2 - sqrt(disc)) / p.Bo : dz / v2) : 0.0
        end
    end

    (tc = r.t[i2] - r.t[i1] + t_climb,
     cor = (E_in > 0 && released) ? sqrt(E_out / E_in) : NaN,
     e_in = E_in, e_out = E_out,
     cor_vel = abs(r.v[i2] / r.v[i1]),
     i_first = i1, i_last = i2, v_in = r.v[i1], v_out = r.v[i2],
     min_gap = minimum(minh), n_detected = count(detected),
     released = released, h_thresh = h_thresh,
     measurable = is_measurable(p, h_thresh), we_measured = we_measured(p, h_thresh))
end

"""
    is_measurable(p, h_thresh = 0.02) -> Bool

Whether a drop launched at this Weber number could have reached the measurement line
from free flight at all.

The line sits `h_thresh` above the substrate, and a drop falling to it under gravity
arrives with `v^2 = We + 2*Bo*h_thresh`. Run backwards, a drop that arrives at `sqrt(We)`
was moving at `sqrt(We - 2*Bo*h_thresh)` when it crossed the line -- which is imaginary
once `We < 2*Bo*h_thresh`. There is then no free-flight state at the line to reference the
impact to, and the nominal Weber number has stopped describing the experiment: the drop is
being delivered by gravity, not by its launch.

The ancestor refuses these outright (`if Vn^2 < 2*g*pixel_adim*Ro ... continue`), and so
should any comparison drawn against it. At `Bo = 0.03` the boundary is `We = 1.3e-3`, well
inside the range the low-Weber experiments cover.
"""
is_measurable(p::ImpactParams, h_thresh::Real = 0.02) = p.We > 2 * p.Bo * h_thresh

"""
    we_measured(p, h_thresh = 0.02) -> Float64

The Weber number referenced to the measurement line rather than to the launch.

`We_measured = We - 2*Bo*h_thresh`, which is what the ancestor reports as
`Westar*Vin_exp^2/Vn^2`. The two agree to a part in a thousand once `We` exceeds
`2*Bo*h_thresh` by a decade, and diverge without bound below that -- so a figure plotted
against the nominal Weber number puts its lowest points in the wrong place.
"""
we_measured(p::ImpactParams, h_thresh::Real = 0.02) =
    max(p.We - 2 * p.Bo * h_thresh, 0.0)
