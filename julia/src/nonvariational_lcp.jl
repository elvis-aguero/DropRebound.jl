# A complementarity contact closure for the nonvariational model.
#
# This completes a two-by-two. The formulation can be variational (interior Ritz amplitudes,
# Euler-Lagrange) or nonvariational (surface amplitudes with Reid's exact per-mode damping as
# the closure); the contact can be found by ranking candidate contact counts or by solving a
# complementarity problem. All four combinations run, and that is the point: when two of them
# disagree, the cell that differs identifies the cause. With one formulation and one closure,
# every disagreement is confounded with every other.
#
# WHAT MAKES THIS EXACT. For the Newtonian case every block of `build_residual!` is affine in
# the unknowns -- the kinematics and the centre-of-mass pair are BDF differences, the equation
# of motion is linear in amplitude, velocity and pressure, and both pressure conditions are
# linear. So writing the residual as
#
#     R(X) = J X + R(0)
#
# is not an approximation, it is an identity, and J is the same matrix at every state. Split X
# into the state part and the pressure part, eliminate the state, and the gap becomes an affine
# function of the contact pressure:
#
#     h = A_c p + b
#
# The step condition is then exactly a linear complementarity problem: h >= 0 (the drop does
# not enter the wall), p >= 0 (the wall pushes, never pulls), and p_i h_i = 0 at every node
# (pressure only where there is contact). No candidate contact count is proposed, ranked or
# rejected, because the contact set is an output.
#
# Two practical notes. The construction uses cp = M+1 when building J: that switches the
# pressure terms on -- they are gated behind `cp > 0` in the residual, and the gate is
# vacuous here because a zero pressure vector contributes nothing anyway -- and it makes all
# M+1 pressure-block rows gap rows, which is exactly the map needed. And the complementarity
# is imposed on NODAL pressures rather than on the Legendre coefficients B, because the
# constraint "the wall pushes" is a statement about the pressure at a point; in the coefficient
# basis it is not even sign-definite.
#
# Restricted to the Newtonian model on purpose. Oldroyd-B carries polymer stress as extra
# state and the shear-thinning closures make the viscosity depend on the solution, so in
# neither case is the system affine and neither would give an exact LCP.

"""
    lcp_active_set(A, b; tol) -> (p, resid, pivots)

Solve `h = A p + b`, `h >= 0`, `p >= 0`, `p_i h_i = 0` by active-set pivoting.

WHY NOT THE PROJECTED GAUSS-SEIDEL SWEEP USED ON THE VARIATIONAL SIDE. That sweep is
Gauss-Seidel on the convex quadratic program `min ½ pᵀAp + bᵀp`, and it is only equivalent to
the complementarity problem when `A` is symmetric. The variational formulation delivers that
for free -- its compliance is a Hessian of an energy, so symmetry is structural. This
formulation does not: the compliance here is asymmetric by more than fifty per cent relative,
because the model is assembled from per-mode damping coefficients rather than from a single
quadratic functional, and there is no energy whose second derivative it could be. Sweeping it
anyway converges to something that is not a solution.

It is still positive definite in the sense that `xᵀAx > 0`, hence a P-matrix, hence the
problem has exactly one solution for every `b` -- so the difficulty is the algorithm, not the
formulation. Active-set pivoting is the standard remedy and terminates for P-matrices: guess a
contact set, solve the equality problem on it, then repair the guess by dropping any node that
came out pulling and adding any node the surface has passed through.

The returned residual is TWO-SIDED, unlike `max_i min(h_i, p_i)`: it measures violation of
`h >= 0` and of `p >= 0` separately as well as complementarity. The one-sided measure is small
whenever one of the two is small at each node, so it reports success on a state where the drop
has penetrated the wall -- which is exactly the failure it hid here.
"""
function lcp_active_set(A::AbstractMatrix, b::AbstractVector;
                        tol::Real = 1e-12, maxpivots::Int = 200)
    n = length(b)
    active = falses(n)
    p = zeros(n)
    pivots = 0
    for _ in 1:maxpivots
        pivots += 1
        p .= 0.0
        idx = findall(active)
        if !isempty(idx)
            ps = A[idx, idx] \ (-b[idx])
            p[idx] .= ps
        end
        sc = max(maximum(abs, b), maximum(abs, p), 1.0)
        ## a node that came out pulling cannot be in contact -- drop the worst
        if !isempty(idx)
            j = argmin(p[idx])
            if p[idx[j]] < -tol*sc
                active[idx[j]] = false
                continue
            end
        end
        h = A*p .+ b
        ## a free node the surface has passed through must be in contact -- add the worst
        free = findall(.!active)
        if !isempty(free)
            j = argmin(h[free])
            if h[free[j]] < -tol*sc
                active[free[j]] = true
                continue
            end
        end
        return (p, lcp_residual(A, b, p), pivots)
    end
    (p, lcp_residual(A, b, p), pivots)
end

"""
    lcp_residual(A, b, p) -> Float64

How badly `p` fails to solve the complementarity problem, counting both inequalities and the
complementarity itself, scaled. Zero only for an actual solution.
"""
function lcp_residual(A::AbstractMatrix, b::AbstractVector, p::AbstractVector)
    h = A*p .+ b
    sc = max(maximum(abs, h), maximum(abs, p), 1.0)
    max(maximum(max.(-h, 0.0)), maximum(max.(-p, 0.0)),
        maximum(abs.(h .* p))) / sc
end

"""
    nv_nodal_pressure(cfg, state) -> Vector

Contact pressure evaluated at the M+1 collocation angles. The Legendre coefficients `B`
carried in the state are the pressure's spectrum; this is the pressure itself, which is the
quantity that has to be non-negative.
"""
nv_nodal_pressure(cfg::SimConstants, s::DropState) =
    collect_Pl(cfg.M, cos.(cfg.theta_vec)) * s.B

"""
    nv_gap(cfg, state) -> Vector

Height of the surface above the substrate at the M+1 collocation angles,
`h(θ) = cos θ (1 + Σ_n P_n(cos θ) A_n) + z`. Zero where the drop touches.
"""
function nv_gap(cfg::SimConstants, s::DropState)
    c = cos.(cfg.theta_vec)
    P = collect_Pl(cfg.M, c)
    c .* (1.0 .+ P[:, 3:end] * s.A[2:end]) .+ s.z
end

"""
    nv_compliance(cfg, ob, history, dt) -> (A_c, b, recover, lower)

The affine map from nodal contact pressure to gap for one BDF step, plus a closure that
recovers the full state from a pressure vector. Exact for the Newtonian model.

`A_c[i,j]` is how much the gap at node `i` opens per unit pressure at node `j` -- a compliance,
so it should be symmetric positive semi-definite, and it is that structure the projected
Gauss-Seidel sweep relies on.
"""
function nv_compliance(cfg::SimConstants, ob::OBParams,
                       history::Vector{DropState}, dt::Float64)
    M = cfg.M; sz = 3M + 1; nB = M + 1
    idx_B = (2M-1):(3M-1)                 # pressure coefficients B_0 .. B_M
    idx_s = vcat(1:(2M-2), 3M, 3M+1)      # A, Adot, z, v -- everything else

    ## cp = M+1: pressure terms active, and every pressure-block row is a gap row
    probe = deepcopy(history[end])
    J  = build_jacobian(probe, history, dt, M+1, cfg, ob)
    zst = deepcopy(probe); unpack_X!(zst, zeros(sz), M)
    r0 = zeros(sz); build_residual!(r0, zst, history, dt, M+1, cfg, ob)

    Jss = J[idx_s, idx_s]; JsB = J[idx_s, idx_B]; rd = r0[idx_s]
    Gs  = J[idx_B, idx_s]; g0 = r0[idx_B]
    F = lu(Jss)

    ## state as a function of the pressure coefficients: Jss s + JsB B + rd = 0
    Sb = -(F \ Matrix(JsB))               # ds/dB
    s0 = -(F \ rd)                        # state with no contact pressure
    ## gap = Gs s + g0
    Cb = Gs * Sb
    b  = Gs * s0 .+ g0
    V  = collect_Pl(M, cos.(cfg.theta_vec))   # nodal pressure = V B
    A_c = Cb / V                              # gap = A_c * (nodal p) + b

    ## WHICH NODES MAY CARRY PRESSURE. The collocation angles span the whole sphere, not just
    ## the wetted side, so most of them sit on top of the drop where the expression above is
    ## the height of the crown rather than a distance to the substrate. Imposing
    ## complementarity there is meaningless, and it announces itself: the diagonal of A_c
    ## changes sign exactly at the equator, because pushing down on the crown opens no gap
    ## underneath. Contact is therefore posed only on the lower hemisphere, and the remaining
    ## nodes keep the free condition p = 0 -- which is what the searching closure does as
    ## well, since it always takes its contact angles from the south pole upward.
    lower = findall(<(0.0), cos.(cfg.theta_vec))

    function recover(p::Vector{Float64})
        B = V \ p
        s = s0 .+ Sb * B
        X = zeros(sz)
        X[idx_s] .= s
        X[idx_B] .= B
        X
    end
    (A_c, b, recover, lower)
end

"""
    solve_drop_lcp!(cfg, ob, init; t_end, save_every, dt_init, dt_min)
        -> (times, states)

March the nonvariational model with a complementarity contact closure. Same signature and
same return as [`solve_drop!`](@ref), so the two are interchangeable.

There is no contact search. At each step the gap-versus-pressure map is assembled and the
complementarity problem solved; the contact set falls out as the nodes carrying pressure. A
step is rejected only if the solve fails to certify complementarity or the state stops being
finite, and `dt` then halves -- never because no candidate contact count was admissible,
which is the failure mode the searching closure has.
"""
function solve_drop_lcp!(cfg::SimConstants, ob::OBParams, init::DropState;
                         t_end::Float64      = 10.0,
                         save_every::Float64 = 0.1,
                         dt_init::Float64    = cfg.dt_max,
                         dt_min::Float64     = cfg.dt_max * 1e-6)
    (ob.De1 > 0.0 && ob.beta_s < 1.0) &&
        error("solve_drop_lcp! is Newtonian-only: Oldroyd-B carries extra state, so the " *
              "system is not affine and the complementarity problem is not exact.")
    M = cfg.M
    t = 0.0; dt = dt_init
    curr = deepcopy(init); curr.dt = dt
    hist = DropState[deepcopy(curr)]
    times = Float64[0.0]; out = DropState[deepcopy(curr)]
    next_save = save_every
    while t < t_end
        dt = min(dt, t_end - t)
        dt <= 0 && break
        ok = false
        while dt >= dt_min
            A_c, b, recover, lower = nv_compliance(cfg, ob, hist, dt)
            ps, resid, _ = lcp_active_set(A_c[lower, lower], b[lower])
            p = zeros(length(b)); p[lower] .= ps
            sc = max(maximum(abs, b), maximum(abs, p), 1.0)
            if resid < 1e-8
                X = recover(p)
                if all(isfinite, X)
                    cand = deepcopy(curr)
                    unpack_X!(cand, X, M)
                    cand.dt = dt
                    cand.cp = count(>(1e-10 * sc), p)
                    if all(isfinite, cand.A) && isfinite(cand.z)
                        push!(hist, cand); length(hist) > 2 && popfirst!(hist)
                        curr = cand; t += dt; ok = true
                        break
                    end
                end
            end
            dt /= 2
        end
        ok || break
        if t >= next_save - 1e-12
            push!(times, t); push!(out, deepcopy(curr))
            next_save += save_every
        end
        dt = min(dt * 1.1, cfg.dt_max)      # ramp back, as the searching solver does
    end
    (times, out)
end
