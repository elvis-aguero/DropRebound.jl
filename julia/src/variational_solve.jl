# Time integration of the variational model, with contact.
#
# !!! STATUS: NOT WORKING. Do not use, and do not read a result out of it. !!!
#
# The assembly this builds on IS verified -- see test_variational.jl, where the
# Newtonian limit reproduces Reid to machine precision and the shear-thinning
# coupling obeys the Gaunt selection rule exactly. This file is the time stepper on
# top of it, and it stalls at contact onset. It is committed with the diagnosis
# recorded rather than deleted, because the diagnosis is most of the remaining work.
#
# WHERE IT STOPS. The drop free-falls, closes the initial gap, and then every
# candidate contact count is rejected, dt halves to the floor, and the run ends at
# t ~ 0.03 having never entered contact. So there is no coefficient of restitution
# and no contact time out of this yet.
#
# WHAT IS ALREADY RULED OUT, and each cost a fix:
#   * p_c indexed by shape mode. It must live on l = 0..Lmax, because p_{c,1} is
#     what drives the centre of mass and l = 0, 1 are not shape modes. Fixed.
#   * activating only the lowest cp pressure harmonics. l = 0 does no work on any
#     shape mode, so a one-node contact could not enforce its own gap condition and
#     the KKT system was singular. Fixed by following the validated ancestor: ALL
#     coefficients are unknowns, closed by gap = 0 at contact nodes and p_c = 0 at
#     free nodes -- a square nodal system.
#   * testing the sign of the pressure COEFFICIENTS. p_c is concentrated near
#     mu = -1, so p_{c,1} is negative for a perfectly admissible push; the condition
#     is on the reconstructed FIELD. Fixed.
#   * starting with the pole exactly at h = 0. That makes the first contact step
#     impulsive: the shape must flatten within one dt and the solve answers with
#     |zeta| ~ 10 and pressure coefficients in the hundreds. Fixed by starting above
#     contact so the dt control can resolve the onset.
#   * absolute feasibility tolerances. At first touch the pressure is essentially
#     zero and its sign is numerically indeterminate -- the physical solution sits ON
#     the boundary, which is what makes onset degenerate. Made relative. NOT
#     sufficient: the stall survives this.
#
# WHAT I WOULD LOOK AT NEXT, in order. The surviving suspicion is that the onset
# step is still over-constrained: a single contact node pins the gap while the
# square nodal system simultaneously forces p_c to vanish at the other Lmax nodes,
# which may leave no admissible pressure at all for a nearly-tangent surface.
# LowWeberDropRebound seeds the first contact from the geometric crossing rather
# than solving for it, and SpectralKM's provenance records the tangency residual as
# identically degenerate at onset -- so the first step probably needs its own rule,
# which is exactly what both ancestors do and what this file does not.
#
# The contact machinery otherwise follows LowWeberDropRebound (Gabbard et al. 2025):
#
#   * the contact set is an INTEGER -- the number of angular nodes in contact --
#     and it changes by at most one per step;
#   * candidates leaving the surface below the substrate outside the contact region
#     are INFEASIBLE and rejected, not ranked;
#   * the objective is local to the contact edge, because one integrated over the
#     patch is biased toward vanishing contact;
#   * when no admissible neighbouring set exists the step is rejected and dt halved.
#     That, not event detection, is what stops the nonsmooth transition being
#     straddled.
#
# BDF2 in both a and adot, with a fixed-point sweep on eta: eta multiplies a stiff
# dissipation operator, so extrapolating it is a stability risk and not only an
# accuracy cost.

using LinearAlgebra

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
    gap0::Float64
end

function ImpactParams(; We, Bo, Oh, ls = 2:16, K = 2, eta = gd -> 1.0,
                      n_nodes = nothing, dt0 = nothing, dt_min = 1e-8,
                      t_max = 30.0, eta_sweeps = 1, gap0 = 0.02)
    lsv = collect(ls); Mmax = maximum(lsv)
    # THE COLLOCATION GRID IS INDEPENDENT OF THE MODE LIST. Tying it to Lmax + 1 was
# a leftover from the global-pressure design, where one node per pressure
# coefficient was what made the system square. With the pressure on the patch
# there is no such constraint, and the grid has a job of its own: it has to
# resolve how far the contact advances in one step. At Lmax = 8 the nodes were
# 10 degrees apart -- a patch far coarser than the drop descends per step -- and
# the shape simply dimpled to keep the neighbouring nodes above the substrate
# instead of building pressure.
nn = something(n_nodes, 90)
    nodes = [pi - (i - 1) * (pi / 2) / nn for i in 1:nn]
    dt = something(dt0, 2pi / (16 * sqrt(Mmax * (Mmax - 1) * (Mmax + 2))))
    ImpactParams(We, Bo, Oh, lsv, K, eta, nodes, dt, dt_min, t_max, eta_sweeps, gap0)
end

basis(p::ImpactParams) = ModalBasis(p.ls, p.K)
trace_vec(p::ImpactParams, l) = [phi(RitzBasis(l, p.K), k, 1.0) for k in 1:p.K]

"""Surface amplitudes from the interior displacements: `zeta_l = chi_l(1)`."""
function surface_amplitudes(p::ImpactParams, a::AbstractVector)
    b = basis(p)
    [dot(trace_vec(p, l), view(a, dofindex(b, i, 1):dofindex(b, i, p.K)))
     for (i, l) in enumerate(b.ls)]
end

"""Row of the linear map `a -> mu (1 + sum zeta_l P_l)` at one node, and the constant."""
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
    (row, mu)          # h = row . a + mu + z
end

gap(p::ImpactParams, a, z, th) = (r = gap_row(p, th); dot(r[1], a) + r[2] + z)

# THE FILM PRESSURE LIVES ON THE CONTACT PATCH, not on the whole sphere.
#
# This is the one design decision in the file that had to be got right, and getting
# it wrong is instructive. The natural-looking alternative is to carry a global
# pressure `p_c = sum_l p_{c,l} P_l(mu)` with as many coefficients as angular nodes,
# and to close the system by demanding the gap vanish at the contact nodes and the
# PRESSURE vanish at the free ones. That system is square and perfectly consistent
# -- its relative residual is 1e-10 -- and it is physically wrong. A degree-Lmax
# polynomial forced through zero at every free node is a near-delta spike at the
# pole: its nodal VALUE is small, but its `l = 1` moment, which is the entire net
# vertical force, is not. So a one-node contact launched the drop at CoR = 12.
#
# A contact of vanishing extent must exert vanishing net force. That is guaranteed
# only if the pressure's SUPPORT shrinks with the patch, so the pressure is expanded
# in shifted Legendre polynomials on `mu in [-1, mu_c]` with one coefficient per
# contact node -- exactly as many unknowns as gap conditions, and no free-node
# conditions at all, because the support is imposed rather than inferred. SpectralKM
# makes the same choice for the same reason.

"""Contact-edge cosine: the first FREE node bounds the patch from above."""
function patch_mu_c(p::ImpactParams, cp::Int)
    cp <= 0 && return -1.0
    th_last = p.nodes[cp]
    th_next = cp < length(p.nodes) ? p.nodes[cp + 1] : 2 * p.nodes[cp] - p.nodes[cp - 1]
    cos((th_last + th_next) / 2)     # the contact edge lies between the two
end

"""Shifted Legendre basis function `j = 1..cp` on the patch, as a function of mu."""
function patch_basis(mu_c::Float64, j::Int, mu::Real)
    mu > mu_c && return 0.0
    s = 2 * (mu - (-1.0)) / (mu_c + 1.0) - 1.0        # patch -> [-1, 1]
    legendre_angular(j - 1, s).P
end

"""
Coupling of patch pressure mode `j` to the spherical harmonic `l`.

`Q_{zeta_l} = -(4 pi/(2l+1)) p_{c,l}` and `p_{c,l} = ((2l+1)/2) int Q_j P_l dmu`, so
the `2l+1` cancels and every harmonic -- including `l = 1`, which carries the net
vertical force `F = -(4 pi/3) p_{c,1}` -- gets the same expression. One formula for
the shape forcing and the centre-of-mass forcing is not a coincidence: both are the
work the film pressure does on a surface motion.
"""
function patch_moment(p::ImpactParams, mu_c::Float64, j::Int, l::Int; nq = 24)
    xq, wq = gauss_legendre_nodes(nq, -1.0, mu_c)
    -2pi * sum(wq[q] * patch_basis(mu_c, j, xq[q]) * legendre_angular(l, xq[q]).P
               for q in eachindex(xq))
end

"""Generalised force on the interior coordinates from patch pressure mode `j`."""
function force_column(p::ImpactParams, mu_c::Float64, j::Int)
    b = basis(p); Q = zeros(ndof(b))
    for (i, l) in enumerate(b.ls)
        c = patch_moment(p, mu_c, j, l)
        tv = trace_vec(p, l)
        for k in 1:p.K
            Q[dofindex(b, i, k)] = c * tv[k]
        end
    end
    Q
end

"""The film pressure reconstructed at a node -- what the sign condition applies to."""
function pc_at(p::ImpactParams, pc::AbstractVector, mu_c::Float64, th::Real)
    isempty(pc) && return 0.0
    mu = cos(th)
    sum(pc[j] * patch_basis(mu_c, j, mu) for j in eachindex(pc))
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
end

function initial_state(p::ImpactParams)
    N = ndof(basis(p))
    # Start slightly ABOVE contact and let free flight bring the drop in. Starting
# with the pole exactly at h = 0 makes the first contact step impulsive: the
# shape has to flatten within one dt, and the solve answers with |zeta| ~ 10 and
# pressure coefficients in the hundreds. With a gap to close, the dt-halving
# control resolves the onset instead.
ImpactState(0.0, p.dt0, zeros(N), zeros(N), 1.0 + p.gap0, -sqrt(p.We), 0, Float64[])
end

"""
    try_step(p, prev, curr, dt, cp) -> (ok, next, edge)

One BDF2 step at a FIXED contact count. `cp` pinned gaps supply `cp` conditions on
the active film-pressure coefficients; the rest are zero, so `p_c` is supported on
the contact region by construction rather than by a penalty.
"""
function try_step(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                  dt::Float64, cp::Int)
    b = basis(p); N = ndof(b)
    r = dt / curr.dt
    c0 = (1 + 2r) / (1 + r); c1 = -(1 + r); c2 = r^2 / (1 + r)
    β = c0 / dt
    hv_a    = (c1 * curr.a    + c2 * prev.a)    / dt     # adot   = β a + hv_a
    hv_adot = (c1 * curr.adot + c2 * prev.adot) / dt     # addot  = β adot + hv_adot
    hv_z    = (c1 * curr.z    + c2 * prev.z)    / dt
    hv_v    = (c1 * curr.v    + c2 * prev.v)    / dt

    adot_star = (1 + r) * curr.adot - r * prev.adot
    # A FEW PRESSURE MODES, ENFORCED WEAKLY OVER THE WHOLE PATCH.
    #
    # This is the design, and both halves of it are load-bearing. The count of
    # pressure modes is small and FIXED -- it does not grow with the contact -- while
    # the gap condition is imposed at EVERY contacting node, so the constraint rows
    # outnumber the pressure unknowns and the gap is satisfied in least squares rather
    # than pointwise. Each half fixes a failure the other does not:
    #
    #   * one shifted-Legendre mode per contact node is square but over-flexible. The
    #     true film pressure has a square-root edge singularity, a polynomial chasing
    #     it oscillates, and the reconstructed pressure ran -263, -4690, -27579 as the
    #     contact grew to two, three, four nodes.
    #   * a single amplitude pinning a single node is stable but too weak to hold the
    #     patch flat. The contact energy went into a local inertial dimple instead of
    #     decelerating the centre of mass: at z = 0.962 the drop had dimpled PAST the
    #     substrate, every gap was positive, and the contact released with the drop
    #     four per cent below the wall and still travelling at its impact speed.
    #
    # Three modes against tens of nodes is the ancestor's ratio, and for the same
    # reason: the pressure needs enough freedom to carry a flattening profile, and
    # least squares over a fine grid is what keeps that freedom from being spent on
    # the singular edge.
    np = min(cp, 3)
    mass = 4pi / 3
    mu_c = patch_mu_c(p, cp)
    fz = [patch_moment(p, mu_c, j, 1) for j in 1:np]   # net vertical force per mode
    Qm = np == 0 ? zeros(N, 0) : hcat((force_column(p, mu_c, j) for j in 1:np)...)

    a_next = copy(curr.a); pc = zeros(np); z_next = curr.z; v_next = curr.v

    for _ in 0:p.eta_sweeps
        F = assemble_coupled(b, p.Oh; eta = eta_field(p, adot_star))
        A = β^2 * F.M + β * F.C + F.G
        z0 = ((-p.Bo - hv_v) / β - hv_z) / β        # where the pole falls unforced
        rhs0 = -F.M * (β * hv_a + hv_adot) - F.C * hv_a

        if np == 0
            a_next = A \ rhs0
            pc = Float64[]
            v_next = (-p.Bo - hv_v) / β
        else
            Hm = zeros(cp, N); Zm = zeros(cp, np); rhs_h = zeros(cp)
            for i in 1:cp
                row, mu = gap_row(p, p.nodes[i])
                Hm[i, :] = row
                for j in 1:np
                    Zm[i, j] = fz[j] / (mass * β^2)   # z is affine in the pressure
                end
                rhs_h[i] = -mu - z0
            end
            # Eliminate the interior coordinates against the momentum equations, then
            # solve the OVERDETERMINED gap conditions for the few pressure amplitudes.
            # a = A^-1 (rhs0 + Qm pc), so the gap residual is affine in pc alone.
            W = A \ hcat(rhs0, Qm)
            S = Hm * W[:, 2:end] + Zm
            rvec = Hm * W[:, 1] - rhs_h
            # TRUNCATED least squares, and the truncation is physics rather than
            # hygiene. Over a narrow patch the pressure modes act on the gap almost
            # collinearly: whatever their profile, a small patch influences the
            # surface mainly through its net force, so the higher moments change the
            # gap in nearly the same way and only about one combination of amplitudes
            # is actually determined. Solving anyway gives the undetermined
            # combinations enormous values that cancel -- the reconstructed pressure
            # ran to -1.9e9 while the gap stayed satisfied to 1e-5. Dropping the
            # directions the gap cannot see sets them to zero instead, which is the
            # right answer for a quantity no equation constrains: the scheme then
            # degrades gracefully to a uniform pressure on a narrow patch and takes
            # up profile structure only once the patch is wide enough to demand it.
            U, sv, V = svd(S)
            keep = sv .> 1e-8 * sv[1]
            pc = -V[:, keep] * ((U[:, keep]' * rvec) ./ sv[keep])
            a_next = W[:, 1] + W[:, 2:end] * pc
            v_next = (-p.Bo + dot(fz, pc) / mass - hv_v) / β
        end
        z_next = (v_next - hv_z) / β
        adot_star = β * a_next + hv_a
    end
    adot_next = β * a_next + hv_a
    nxt = ImpactState(curr.t + dt, dt, a_next, adot_next, z_next, v_next, cp, copy(pc))

    # NON-PENETRATION, in absolute units of the drop radius -- a fixed length, and not
    # the quantity under test. An earlier version scaled every tolerance by the maximum
    # of the field it was testing, which at onset compared a 1e-4 pressure against a
    # threshold derived from that same 1e-4: an absolute test at exactly the point
    # where the physical solution sits ON the boundary.
    gaps = [gap(p, a_next, z_next, th) for th in p.nodes]
    ok = all(i -> gaps[i] > -1e-6, (cp+1):length(p.nodes))
    # The patch interior is held to a looser tolerance because it is held WEAKLY: the
    # least-squares gap condition leaves a residual there by construction, and
    # demanding it vanish would be demanding the pressure space be as large as the
    # contact -- the over-flexible design this replaced.
    cp > 0 && (ok &= minimum(view(gaps, 1:cp)) > -1e-2)
    # THE PRESSURE MUST PUSH, tested on the reconstructed field rather than on the
    # coefficients: p_c is concentrated near mu = -1, so its l = 1 coefficient is
    # negative for a perfectly admissible push.
    # ACROSS THE WHOLE PATCH, not only at the contact nodes. The nodes all sit in the
    # lower part of the patch -- the patch edge lies beyond the last of them -- so
    # sampling the nodes alone leaves the outer patch unwatched. With more than one
    # pressure mode the profile then goes positive at every node and strongly negative
    # outside them, and the net l = 1 moment comes out downward: the contact
    # ACCELERATED the drop into the wall, from -1.0 to -1.83 in one step, with every
    # nodal pressure admissible. Signorini is a condition on the field, everywhere the
    # field is supported.
    if cp > 0
        mus = range(-1.0, mu_c; length = 16)
        pcs = [pc_at(p, pc, mu_c, acos(clamp(m, -1, 1))) for m in mus]
        ok &= minimum(pcs) > -1e-6 * max(maximum(pcs), 1.0)
    end
    edge = if !ok
        Inf
    elseif cp == 0
        abs(gap(p, a_next, z_next, p.nodes[1]))
    elseif cp < length(p.nodes)
        abs(gap(p, a_next, z_next, p.nodes[cp+1]))
    else
        Inf
    end
    (ok, nxt, edge)
end

"""
    simulate(p) -> NamedTuple

March the impact. Returns the trajectory and the two KPIs: coefficient of
restitution and contact time.
"""
function simulate(p::ImpactParams; verbose = false)
    s0 = initial_state(p)
    prev, curr = s0, s0
    ts = Float64[0.0]; zs = Float64[s0.z]; vs = Float64[s0.v]; cps = Int[0]
    dt = p.dt0
    nrej = 0
    while curr.t < p.t_max
        moved = false
        # ACTIVE-SET ORDER: stay, then grow, then release. Not smallest-first. The
        # two feasibility conditions are one-sided, so the order they are tried in is
        # part of the rule and not an implementation detail. A candidate with no
        # contact passes the pressure test vacuously, so smallest-first releases the
        # drop the instant the pole gap is marginally non-negative, and the contact
        # can then never grow past one node: the run flip-flops between cp = 0 and
        # cp = 1 and stalls with the drop already 0.06 below the substrate. Trying the
        # incumbent set first, growing only against penetration and releasing only
        # against a suction pressure, is the standard active-set iteration -- and it
        # is what the complementarity pair actually says.
        best = nothing; best_cp = curr.cp
        # Growth of more than one node per step is allowed, because the contact radius
        # goes like sqrt(penetration) and therefore sweeps arbitrarily fast at onset:
        # a strict one-node limit there is not a stability safeguard, it just collapses
        # dt trying to resolve a square-root corner. Release stays one node at a time,
        # where no such corner exists.
        # THE CONTACT EXTENT IS SEEDED FROM THE GEOMETRIC CROSSING, not crawled up one
        # node at a time. This is the single most important thing the stepper gets from
        # LowWeberDropRebound, and doing it the other way fails in a way that looks
        # like a solver problem and is not. Growing the set incrementally means the
        # first contact step engages one node, a one-node contact barely decelerates
        # anything, so the drop sinks another 0.016 before the next step -- and the
        # crawl never catches the crossing, which advances like sqrt(penetration). Four
        # steps in, the drop sat four per cent BELOW the substrate with two nodes in
        # contact, and the pressure needed to lift the surface back out in a single
        # step was 1e6. Nothing about that is a conditioning failure; the state was
        # already wrong. Seeding from where the surface actually crosses the wall keeps
        # the extent right from first touch, and the crossing then advances on its own.
        zfree = curr.z + dt * curr.v
        gx = count(th -> gap(p, curr.a, zfree, th) <= 0, p.nodes)
        cands = Int[]
        for c in (gx, gx + 1, gx - 1, gx + 2, gx - 2, curr.cp, 0)
            0 <= c <= length(p.nodes) && c ∉ cands && push!(cands, c)
        end
        for cand in cands
            ok, nxt, _ = try_step(p, prev, curr, dt, cand)
            if ok
                best, best_cp = nxt, cand
                break
            end
        end
        if best === nothing
            dt /= 2
            nrej += 1
            dt < p.dt_min && break
            continue
        end
        prev, curr = curr, best
        push!(ts, curr.t); push!(zs, curr.z); push!(vs, curr.v); push!(cps, curr.cp)
        moved = true
        # leave once the drop is clear of the wall and rising
        if curr.cp == 0 && curr.v > 0 && curr.z > 1.0 + p.gap0 && any(>(0), cps)
            break
        end
        dt = min(dt * 1.05, p.dt0)
    end
    (t = ts, z = zs, v = vs, cp = cps,
     cor = restitution(vs, cps, p), tc = contact_time(ts, cps), rejects = nrej)
end

"""Coefficient of restitution: rebound speed over impact speed."""
function restitution(vs, cps, p::ImpactParams)
    idx = findall(>(0), cps)
    isempty(idx) && return NaN
    j = last(idx)
    j >= length(vs) && return NaN
    vout = maximum(view(vs, j:length(vs)))
    vout <= 0 ? NaN : vout / sqrt(p.We)
end

"""
Contact time as the LONGEST CONTIGUOUS interval, not the first-to-last span.

A marginally resolved run chatters -- brief spurious re-entries just after the
physical rebound -- and the span then overstates the contact time while the FIRST
interval understates it badly on exactly those runs. This is recorded in
SpectralKM's provenance as a defect that survived a dozen call sites.
"""
function contact_time(ts, cps)
    best = 0.0; run_start = -1.0
    for i in eachindex(cps)
        if cps[i] > 0 && run_start < 0
            run_start = ts[i]
        elseif cps[i] == 0 && run_start >= 0
            best = max(best, ts[i] - run_start); run_start = -1.0
        end
    end
    run_start >= 0 && (best = max(best, ts[end] - run_start))
    best
end
