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
    nn = something(n_nodes, Mmax + 1)   # square: one node per p_c coefficient
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

"""
Film-pressure coefficients live on `l = 0..Lmax`, NOT on the shape modes: `p_{c,1}`
is what drives the centre of mass and `l = 0, 1` are absent from the shape list.
Indexed `1..Lmax+1` for `l = 0..Lmax`.
"""
pc_len(p::ImpactParams) = maximum(p.ls) + 1
pc_l(j::Int) = j - 1

"""Generalised force on the interior coordinates from the film harmonic `l = pc_l(j)`.

`Q_{zeta_l} = -(4 pi/(2l+1)) p_{c,l}`, and only harmonics that coincide with a
retained SHAPE mode do work on the surface -- `l = 0` changes volume and `l = 1`
translates, neither of which is a shape mode, so both give a zero column here. The
`l = 1` harmonic still acts, through the centre-of-mass equation.
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

"""The film pressure reconstructed at a node -- what the sign condition applies to."""
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
end

function initial_state(p::ImpactParams)
    N = ndof(basis(p))
    # Start slightly ABOVE contact and let free flight bring the drop in. Starting
# with the pole exactly at h = 0 makes the first contact step impulsive: the
# shape has to flatten within one dt, and the solve answers with |zeta| ~ 10 and
# pressure coefficients in the hundreds. With a gap to close, the dt-halving
# control resolves the onset instead.
ImpactState(0.0, p.dt0, zeros(N), zeros(N), 1.0 + p.gap0, -sqrt(p.We), 0, zeros(pc_len(p)))
end

"""
    try_step(p, prev, curr, dt, cp) -> (ok, next, edge)

One BDF2 step at a FIXED contact count. `cp` pinned gaps supply `cp` conditions on
the active film-pressure coefficients; the rest are zero, so `p_c` is supported on
the contact region by construction rather than by a penalty.
"""
function try_step(p::ImpactParams, prev::ImpactState, curr::ImpactState,
                  dt::Float64, cp::Int)
    b = basis(p); N = ndof(b); npc = pc_len(p)
    r = dt / curr.dt
    c0 = (1 + 2r) / (1 + r); c1 = -(1 + r); c2 = r^2 / (1 + r)
    β = c0 / dt
    hv_a    = (c1 * curr.a    + c2 * prev.a)    / dt     # adot   = β a + hv_a
    hv_adot = (c1 * curr.adot + c2 * prev.adot) / dt     # addot  = β adot + hv_adot
    hv_z    = (c1 * curr.z    + c2 * prev.z)    / dt
    hv_v    = (c1 * curr.v    + c2 * prev.v)    / dt

    adot_star = (1 + r) * curr.adot - r * prev.adot
    a_next = copy(curr.a); pc = zeros(npc); z_next = curr.z; v_next = curr.v

    for _ in 0:p.eta_sweeps
        F = assemble_coupled(b, p.Oh; eta = eta_field(p, adot_star))
        A = β^2 * F.M + β * F.C + F.G
        rhs0 = -F.M * (β * hv_a + hv_adot) - F.C * hv_a
        # ALL film-pressure coefficients are unknowns, closed by a square set of
        # nodal conditions: the gap vanishes at the cp contact nodes, and the
        # PRESSURE vanishes at the remaining free nodes. That is the design of the
        # validated ancestor, and it is what makes the system non-singular -- an
        # earlier version activated only the lowest cp harmonics, but l = 0 does no
        # work on any shape mode and l = 1 acts only through the centre of mass, so
        # a one- or two-node contact had no way to enforce its own gap condition.
        Qm = hcat((force_column(p, jj) for jj in 1:npc)...)
        z0 = ((-p.Bo - hv_v) / β - hv_z) / β
        dz_dpc1 = -1 / β^2                    # z is affine in the l = 1 harmonic
        Hm = zeros(npc, N); Zm = zeros(npc, npc); rhs_h = zeros(npc)
        for i in 1:npc
            th = p.nodes[i]
            if i <= cp                        # contact: h = 0
                row, mu = gap_row(p, th)
                Hm[i, :] = row
                Zm[i, 2] += dz_dpc1
                rhs_h[i] = -mu - z0
            else                              # free: p_c = 0
                mu = cos(th)
                for jj in 1:npc
                    Zm[i, jj] = legendre_angular(pc_l(jj), mu).P
                end
            end
        end
        sol = [A -Qm; Hm Zm] \ [rhs0; rhs_h]
        a_next = sol[1:N]
        pc = sol[N+1:N+npc]
        v_next = (-p.Bo - pc[2] - hv_v) / β
        z_next = (v_next - hv_z) / β
        adot_star = β * a_next + hv_a
    end
    adot_next = β * a_next + hv_a
    nxt = ImpactState(curr.t + dt, dt, a_next, adot_next, z_next, v_next, cp, copy(pc))

    # Tolerances are RELATIVE to the scale each quantity actually has. At first
# touch the pressure is essentially zero and its sign is numerically
# indeterminate -- the physical solution sits ON the feasibility boundary, which
# is what makes onset degenerate. An absolute 1e-8 rejects it and the run stalls
# at t = 0.026 with every candidate infeasible.
pcs = [pc_at(p, pc, th) for th in p.nodes]
pc_scale = max(maximum(abs, pcs), 1e-12)
gap_scale = max(maximum(abs, [gap(p, a_next, z_next, th) for th in p.nodes]), 1e-12)
ok = all(i -> gap(p, a_next, z_next, p.nodes[i]) > -1e-6 * gap_scale,
         (cp+1):length(p.nodes))
    # the sign condition is on the reconstructed FIELD at the contact nodes, not on
# the coefficients: p_c is concentrated near mu = -1, so p_{c,1} is NEGATIVE for a
# perfectly admissible push, and testing coefficients rejects every real solution.
ok &= all(i -> pcs[i] > -1e-4 * pc_scale, 1:cp)
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
        best = nothing; best_edge = Inf; best_cp = curr.cp
        for cand in max(0, curr.cp - 1):min(curr.cp + 1, length(p.nodes))
            ok, nxt, edge = try_step(p, prev, curr, dt, cand)
            if ok && edge < best_edge
                best, best_edge, best_cp = nxt, edge, cand
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
