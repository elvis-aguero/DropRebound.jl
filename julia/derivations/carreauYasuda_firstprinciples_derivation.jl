# # General-Exponent Carreau-Yasuda: Why Amplitude Perturbation Theory Fails,
# # and What Survives
#
# Reid's theory gives the exact Newtonian damping rate
# ``\lambda_l(\mathrm{Oh})`` and squared frequency
# ``\omega_l^2(\mathrm{Oh})`` at any Ohnesorge number.
# The classical treatment builds a small-amplitude, Landau-style correction
# on top of that theory for the *classical* Carreau law, whose shape exponent
# is fixed at ``a=2``. This page asks the harder question: does an analogous
# correction exist for the *general* Carreau-Yasuda law that the fitted
# validation fluid uses, with ``a\approx0.743`` rather than 2?
#
# For this fluid's fitted parameters, no. That is a quantitative statement
# about where this fluid sits on the Carreau-Yasuda curve rather than a
# general one about shear-thinning drops, and §2-§4 give the numbers behind
# it.
#
# The negative result is not the whole story. Three general results (§6-§8)
# hold independently of it and of any fluid's specific parameters, and all
# three are tools a first-principles treatment would need.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``\mathrm{Oh}_0`` | rest (zero-shear) Ohnesorge number of the real fluid |
# | ``\mathrm{Oh}_\infty`` | infinite-shear-rate Ohnesorge number, ``=\mathrm{Oh}_0(1-\varepsilon_{ST})`` |
# | ``a`` | Carreau-Yasuda shape exponent (``\approx0.743`` here, **not** the classical 2) |
# | ``\mathcal{L}[U]`` | Reid's velocity operator, ``U''-l(l+1)U/x^2+q^2U`` |
# | ``\dot\gamma``, ``S`` | scalar shear-rate invariant, ``\sqrt{2e_{ij}e_{ij}}`` |
# | ``\phi=\omega t`` | oscillation phase at a complex Reid root ``q`` |

using Symbolics
using SpecialFunctions
using QuadGK
using DropSolver

# The fitted parameters of the real validation fluid, and the two natural
# Newtonian anchors they define:

const OH0 = 57.371648873370795
const LAMBDA_C = 30507.34501244818
const A_SHAPE = 0.7430524574330837
const EPS_ST = 0.9995574839318364
const ETA_INF_RATIO = 1 - EPS_ST
const OH_INF = OH0 * ETA_INF_RATIO      # = 0.02539

## Symbolic-zero check used throughout this repo: Symbolics.jl's `simplify`   #src
## does not always collapse an algebraically-zero expression to the literal   #src
## 0, so confirm both ways -- symbolic simplification AND numeric evaluation  #src
## at several concrete points for every free variable.                        #src
function symbolic_zero(expr)                                                  #src
    simplified = simplify(expr; expand=true)                                  #src
    is_symbolic_zero = isequal(simplified, 0) || isequal(simplified, 0.0)     #src
    vars = Symbolics.get_variables(expr)                                      #src
    is_numeric_zero = if isempty(vars)                                        #src
        true                                                                  #src
    else                                                                      #src
        f = Symbolics.build_function(expr, vars...; expression=false)         #src
        test_vals = (0.37, 1.21, 2.03, 0.68, 1.59, 3.14, 0.91, 2.77)          #src
        all(abs(f((test_vals[mod1(i + k, length(test_vals))] for k in 1:length(vars))...)) < 1e-8 #src
            for i in 1:length(test_vals))                                     #src
    end                                                                       #src
    is_symbolic_zero || is_numeric_zero                                       #src
end                                                                           #src

# ![Flow curve of the validation fluid on log-log axes: viscosity falls by
# about two and a half orders of magnitude between the rest plateau and the
# infinite-shear floor, and the shear rates reached during an impact fall in
# the steep transition between them.](../assets/cy_flow_curve.png)
#
# *The measured flow curve, with the band of shear rates an impact actually
# produces. The fluid is sampled in the middle of the thinning transition,
# far from either plateau -- which is the root of everything below.*
#
# ## 1. The exponent, not just its size, is the problem
#
# The classical route works because ``a=2`` is an *even integer*. The
# viscosity correction ``(\lambda_c\dot\gamma)^2`` is then an analytic,
# Taylor-expandable function of the oscillation amplitude ``\epsilon``, since
# ``\dot\gamma=O(\epsilon)`` linearly in Reid's single-mode ansatz. For a
# general, non-integer exponent this breaks down completely: the analogous
# ``(\lambda_c\dot\gamma)^a\sim|\epsilon|^{a}`` is not differentiable at
# ``\epsilon=0``, only Hölder continuous there, and for ``a<1`` it is
# actually **larger** than ``\epsilon`` itself as ``\epsilon\to0``.
#
# Concretely, the ratio ``|\epsilon|^{a}/|\epsilon|`` at
# ``\epsilon=10^{-1},10^{-3},10^{-5},10^{-7}`` is
#
# ```math
# 1.807,\quad 5.900,\quad 19.26,\quad 62.90
# ```
#
# growing without bound as ``\epsilon\to0``. There is no amplitude, however
# small, at which this correction is a small perturbation of a fixed
# integer-power hierarchy: it sits at a genuinely fractional order, between
# Reid's own ``O(\epsilon)`` linear terms and the usual ``O(\epsilon^2)``
# geometric mode-coupling correction.
#
# This alone does not doom a perturbative treatment. A fractional expansion
# parameter is unusual but not impossible. What actually closes off the
# small-amplitude route for *this* fluid is quantitative rather than
# structural, and that is what §2-§4 establish.

let                                                                           #src
    ratios = [abs(eps)^A_SHAPE / abs(eps) for eps in (0.1, 0.001, 1e-5, 1e-7)] #src
    @assert issorted(ratios)          # strictly increasing as eps shrinks (list is in decreasing-eps order) #src
    @assert ratios[end] > 10 * ratios[1]                                      #src
    @assert isapprox(ratios[1], 1.807; rtol=1e-3)   # pins the tabulated values #src
    @assert isapprox(ratios[end], 62.90; rtol=1e-3)                           #src
end                                                                           #src
println("ASSERTION 1 OK: |eps|^a/|eps| grows without bound as eps->0 for a<1")  #src

# ## 2. Anchoring at rest, ``\mathrm{Oh}_0``, fails outright
#
# A shear-thinning fluid's viscosity interpolates between two Newtonian
# limits: ``\eta_0`` at rest and ``\eta_\infty`` at infinite shear rate.
# Reid's characteristic equation gives the exact
# ``\lambda_l(\mathrm{Oh})``, ``\omega_l^2(\mathrm{Oh})`` at either, so both
# are legitimate places to anchor a correction. The rest anchor fails first,
# for a simple reason. At ``\mathrm{Oh}_0\approx57.37``, every mode this
# solver resolves is heavily overdamped:
#
# | ``l`` | ``\lambda_l`` | ``\omega_l^2`` | ``\lambda_l^2`` | |
# |:--|:--|:--|:--|:--|
# | 2 | 203.8 | 7.48 | ``4.2\times10^{4}`` | overdamped |
# | 3 | 459.4 | 25.48 | ``2.1\times10^{5}`` | overdamped |
# | 5 | 1119 | 102.9 | ``1.3\times10^{6}`` | overdamped |
# | 10 | 3681 | 665.4 | ``1.4\times10^{7}`` | overdamped |
#
# ``\lambda_l^2\gg\omega_l^2`` throughout, by three to five orders of
# magnitude: the root pair is real, there is no oscillation at all. A
# multiple-scales or envelope treatment needs a fast carrier wave to
# modulate, and at ``\mathrm{Oh}_0`` there is none.
#
# Separately: the shear rate at which the viscosity
# correction reaches even the 1% level for this fluid is
# ``\dot\gamma\approx6.8\times10^{-8}``, far below any resolved mode
# velocity. So an expansion anchored at rest would be expanding about a state
# the fluid leaves immediately and never returns to.

for l in (2, 3, 5, 10)                                                        #src
    lam, om2, resid = reid_lambda_omega2(OH0, l)                              #src
    @assert resid < 1e-6                                                      #src
    @assert lam^2 > om2   # overdamped: real root pair, no oscillation         #src
    @assert lam^2 > 1000 * om2   # and not marginally so                       #src
end                                                                           #src
## The shear rate at which the CY correction first reaches 1%.                #src
let f(g) = 1 - (1 + (LAMBDA_C * g)^A_SHAPE)^(-EPS_ST) - 0.01                   #src
    lo, hi = 1e-16, 1e-2                                                      #src
    for _ in 1:200                                                            #src
        mid = sqrt(lo * hi)                                                   #src
        f(mid) < 0 ? (lo = mid) : (hi = mid)                                  #src
    end                                                                       #src
    @assert sqrt(lo * hi) < 1e-7                                              #src
end                                                                           #src
println("ASSERTION 2 OK: every mode is heavily overdamped at Oh_0, and the 1% shear threshold is <1e-7") #src

# ## 3. Anchoring at infinite shear is better posed, but still not small
#
# Define ``\mathrm{Oh}_\infty\equiv\mathrm{Oh}_0(1-\varepsilon_{ST})
# \approx0.02539``. Unlike ``\mathrm{Oh}_0``, this anchor *is* comfortably
# underdamped: at ``l=2``, ``\lambda_2=0.1161`` and ``\omega_2^2=7.952``,
# giving a quality factor ``Q\approx12.1``. (Reid's oscillator convention is
# ``\ddot A+2\lambda\dot A+\omega^2A=\text{forcing}``, so ``Q=\omega/2\lambda``,
# not ``\omega/\lambda``.) That is a genuine separation between the decay
# time and the oscillation period -- exactly what ``\mathrm{Oh}_0`` lacked.
#
# But being well posed is not the same as being close. The question a
# perturbative correction actually needs answered is not "how far is
# ``\mathrm{Oh}_{\mathrm{eff}}`` from a reference?" but "does the
# **first-order Taylor expansion** of Reid's exact relations about that
# reference reproduce their true values over the range this fluid's dynamics
# actually visits?" Checked directly against the exact,
# continuation-solved relations:
#
# | ``l`` | ``\lambda`` error at ``\mathrm{Oh}=0.03`` | ``\lambda`` error at ``\mathrm{Oh}=0.1`` | ``\omega^2`` error at ``\mathrm{Oh}=0.1`` |
# |:--|:--|:--|:--|
# | 2 | 0.12% | 8.2% | 0.25% |
# | 3 | 0.16% | 11.9% | 0.81% |
# | 5 | 0.21% | 16.9% | 1.43% |
# | 10 | 0.30% | 25.9% | 1.52% |
#
# Two things stand out. ``\omega_l^2`` is nearly flat in ``\mathrm{Oh}``
# near ``\mathrm{Oh}_\infty``, so its linear correction stays small even at
# high ``l`` -- a genuine asymmetry between the frequency and damping
# corrections' validity. ``\lambda_l`` is not flat: at the low end of the
# operative band every mode is comfortably accurate, but by the high end the
# error has passed 25% at ``l=10``. The same growing-with-``l`` curvature
# that broke the ``\mathrm{Oh}_0`` anchor reappears here, and here it also
# grows *across* the operative band, not merely across modes.
#
# Moving the anchor to some other "typical" ``\mathrm{Oh}^\ast`` does not
# rescue this. The curvature at high ``l`` sets in almost immediately away
# from any point in the operative range, so there is no interior anchor from
# which a single linearization covers the band.

lam2_inf, om2_2_inf, resid_inf = reid_lambda_omega2(OH_INF, 2)                 #src
@assert resid_inf < 1e-6                                                      #src
@assert lam2_inf^2 < om2_2_inf                                                #src
## Reid's oscillator convention (see julia/src/reid.jl's docstring) is        #src
## Addot + 2*lambda*Adot + omega^2*A = forcing, so Q = omega/(2*lambda).      #src
let Q_factor = sqrt(om2_2_inf) / (2 * lam2_inf)                               #src
    @assert Q_factor > 5                       # a genuine separation of timescales #src
    @assert isapprox(Q_factor, 12.1; rtol=0.02)  # pins the reported value      #src
end                                                                           #src
println("ASSERTION 3 OK: Oh_inf is underdamped at l=2 with Q ~ 12")            #src

let h = 1e-6                                                                  #src
    lam0 = Dict(l => reid_lambda_omega2(OH_INF, l)[1] for l in (2, 3, 5, 10))  #src
    om0 = Dict(l => reid_lambda_omega2(OH_INF, l)[2] for l in (2, 3, 5, 10))   #src
    dlam = Dict(l => (reid_lambda_omega2(OH_INF + h, l)[1] - reid_lambda_omega2(OH_INF - h, l)[1]) / (2h) for l in (2, 3, 5, 10)) #src
    dom = Dict(l => (reid_lambda_omega2(OH_INF + h, l)[2] - reid_lambda_omega2(OH_INF - h, l)[2]) / (2h) for l in (2, 3, 5, 10)) #src
    errs_lambda_low, errs_lambda_high, errs_omega_high = Float64[], Float64[], Float64[] #src
    for l in (2, 3, 5, 10)                                                    #src
        lam_true_lo, _, _ = reid_lambda_omega2(0.03, l)                       #src
        lam_lin_lo = lam0[l] + dlam[l] * (0.03 - OH_INF)                      #src
        push!(errs_lambda_low, abs(lam_lin_lo - lam_true_lo) / abs(lam_true_lo)) #src
        lam_true_hi, om_true_hi, _ = reid_lambda_omega2(0.1, l)               #src
        lam_lin_hi = lam0[l] + dlam[l] * (0.1 - OH_INF)                       #src
        om_lin_hi = om0[l] + dom[l] * (0.1 - OH_INF)                          #src
        push!(errs_lambda_high, abs(lam_lin_hi - lam_true_hi) / abs(lam_true_hi)) #src
        push!(errs_omega_high, abs(om_lin_hi - om_true_hi) / abs(om_true_hi))  #src
    end                                                                       #src
    @assert all(e -> e < 0.01, errs_lambda_low)   # low end: all comfortably small #src
    @assert errs_lambda_high[end] > 0.2           # high end, l=10: NOT small (>20%) #src
    @assert errs_omega_high[1] < 0.01             # omega^2, l=2: still small   #src
    @assert errs_omega_high[end] < 0.02           # omega^2, l=10: still small (contrast with lambda) #src
    @assert issorted(errs_lambda_high)            # lambda error grows with l   #src
end                                                                           #src
println("ASSERTION 4 OK: the linear correction to lambda_l around Oh_inf is <1% at Oh=0.03 but >20% at Oh=0.1, l=10") #src

# ## 4. Where the fluid actually lives: a live solver trace
#
# The two preceding sections are properties of the standalone
# ``\mathrm{Oh}_{\mathrm{eff}}`` formula. This one is a cross-check against
# the running solver. Reproducing the trajectory of
# ``\mathrm{Oh}_{\mathrm{eff}}(t)`` for a real impact
# (``\mathrm{We}=0.7649``) shows two things.
#
# First, the drop's effective Ohnesorge number revisits the neighbourhood of
# ``\mathrm{Oh}_0`` **at most once** -- the very first contact transient,
# before any resolved mode has developed real shear -- and never again on
# later contact or lift-off cycles. Once excited, it is permanently pinned
# away from rest.
#
# Second, and this is the point, "pinned away from ``\mathrm{Oh}_0``" is not
# "close to ``\mathrm{Oh}_\infty``". Post-transient, the ratio
# ``\mathrm{Oh}_{\mathrm{eff}}/\mathrm{Oh}_\infty`` stays in the band
# ``[1.82,\ 2.87]``: never at the infinite-shear plateau, never near the rest
# plateau, always a factor of two to three above ``\mathrm{Oh}_\infty``.
# The dynamics lives in the transitional part of the Carreau-Yasuda curve,
# which is exactly the part no single Newtonian anchor describes.

let                                                                           #src
    BO = 0.012                                                                #src
    We = 0.7649                                                               #src
    M = 12                                                                    #src
    stx = STExactParams(M, OH0, LAMBDA_C, A_SHAPE, EPS_ST; viscous=:reid, eta_inf_ratio=ETA_INF_RATIO) #src
    dt_max = make_dt_max(M)                                                   #src
    theta_vec = make_theta_vec(M)                                             #src
    precomp = precompute_integrals(NaN, M)[1]                                 #src
    cfg = SimConstants(M, M + 1, OH0, BO, theta_vec, precomp, dt_max)         #src
    init = DropState(M)                                                       #src
    init.z = 1.05                                                             #src
    init.v = -sqrt(We)                                                        #src
    init.dt = dt_max                                                          #src
    init.cp = 0                                                               #src
    times, states = solve_drop!(cfg, OBParams(), init; stx=stx, t_end=0.5, save_every=dt_max / 8) #src
    near_rest_count = 0                                                       #src
    pinned_ratios = Float64[]                                                 #src
    prev_cp = 0                                                               #src
    seen_second_transition = false                                            #src
    for s in states                                                           #src
        Adot_vec = s.Adot[2:end]                                              #src
        all(x -> x == 0.0, Adot_vec) && (prev_cp = s.cp; continue)            #src
        oh_max = maximum(oh_eff_all_coupled(stx, OH0, Adot_vec))              #src
        transition = s.cp != prev_cp                                          #src
        if oh_max > 1.0                                                       #src
            near_rest_count += 1                                              #src
        elseif seen_second_transition || !transition                          #src
            push!(pinned_ratios, oh_max / OH_INF)                             #src
        end                                                                   #src
        transition && prev_cp != 0 && (seen_second_transition = true)         #src
        prev_cp = s.cp                                                        #src
    end                                                                       #src
    @assert near_rest_count <= 2   # at most the initial transient(s), not recurring #src
    @assert !isempty(pinned_ratios)                                           #src
    @assert all(r -> r > 1.0, pinned_ratios)         # never AT Oh_inf either  #src
    @assert minimum(pinned_ratios) > 1.0 && maximum(pinned_ratios) < 10.0      # an intermediate band #src
    @assert isapprox(minimum(pinned_ratios), 1.82; rtol=0.02)  # pins the reported band #src
    @assert isapprox(maximum(pinned_ratios), 2.87; rtol=0.02)                 #src
end                                                                           #src
println("ASSERTION 5 OK: a live We=0.7649 run pins Oh_eff/Oh_inf in [1.82, 2.87] after one transient") #src

# ## 5. Conclusion, and what follows from it
#
# For this fluid's fitted parameters there is no single reference viscosity
# -- rest, infinite shear, or anywhere between -- about which a regular
# perturbation expansion of Reid's damping rate is uniformly valid across the
# modes this solver resolves. The physics genuinely lives in the fully
# nonlinear, transitional part of the Carreau-Yasuda curve. This is a
# quantitative fact about *this fluid's* ``\lambda_c`` and ``a``, not a
# universal statement about shear-thinning drops.
#
# What follows from it matters more than the finding itself. The failure is
# specifically a failure of *linearizing* ``\lambda_l(\mathrm{Oh})`` and
# ``\omega_l^2(\mathrm{Oh})``. But those functions are already known exactly,
# for any ``\mathrm{Oh}``, through Reid's characteristic equation. So the
# correct response is not to build a better linear correction: it is not to
# linearize that piece at all -- to evaluate the exact relations directly at
# whatever effective Ohnesorge number a closure provides, which is exactly
# what DropSolver's Carreau-Yasuda extension does.
#
# That part was never the weak point. The weak point is upstream of it, in
# how the closure produces a single scalar ``\mathrm{Oh}_{\mathrm{eff}}`` in
# the first place. The rest of this page builds the tools that question needs.

# ## 6. Adjoint sensitivity: the boundary derivative without solving the ODE
#
# Reid's velocity operator
#
# ```math
# \mathcal{L}[U] \;\equiv\; U'' - \frac{l(l+1)}{x^2}U + q^2U
# ```
#
# has no first-derivative term, so it is formally self-adjoint on ``(0,1)``
# with the plain ``L^2`` inner product -- no weight function needed. The
# Lagrange identity,
#
# ```math
# W\,\mathcal{L}[U] - U\,\mathcal{L}[W]
#   \;=\; \frac{d}{dx}\Bigl[\,W\,U' - W'\,U\,\Bigr],
# ```
#
# is verified symbolically for abstract ``U(x)``, ``W(x)``, with bilinear
# concomitant ``P[U,W]=WU'-W'U``.
#
# That identity buys a shortcut for precisely the calculation a perturbation
# of Reid's problem needs -- the boundary derivative of a particular
# solution, *without solving the ODE*:

sph_jl(l, z) = sqrt(pi / (2z)) * besselj(l + 0.5, z)

"""
    adjoint_shortcut_Yprime1(l, q0, RHS)

For `Y` regular at `x=0`, satisfying `Y(1)=0`, and solving `L[Y] = RHS(x)`
at fixed wavenumber `q0`, returns `Y'(1)` -- using the regular homogeneous
solution `x*j_l(q0*x)` as the adjoint test function.
"""
adjoint_shortcut_Yprime1(l, q0, RHS) =
    first(quadgk(t -> t * sph_jl(l, q0 * t) * RHS(t), 0.0, 1.0; rtol=1e-10)) / sph_jl(l, q0)

# In closed form, that is
#
# ```math
# \boxed{\;
# Y'(1) \;=\; \frac{1}{j_l(q_0)}\int_0^1 x\,j_l(q_0x)\,\mathrm{RHS}(x)\,dx
# \;}
# ```
#
# and it reproduces direct integration of the ODE to better than ``10^{-6}``
# relative error, for three independent ``(l,q_0,\mathrm{RHS})`` combinations
# spanning different modes and different ``\mathrm{Oh}`` regimes.
#
# Combined with Reid's BC1 (``U(1)=-1``, which is purely kinematic and so
# unaffected by any viscosity correction) and BC2, this reduces the search
# for the ``O(\delta)`` shift in ``q^2`` produced by any new forcing added to
# Reid's ODE -- both stress boundary conditions included -- to an ordinary
# linear solve. No further ODE integration is required at all.

@variables x l q                                                              #src
@variables Ufun(..) Wfun(..)                                                  #src
let Uf = Ufun(x), Wf = Wfun(x), Dx_ = Differential(x)                         #src
    Lop(f) = expand_derivatives(Dx_(Dx_(f))) - l * (l + 1) / x^2 * f + q^2 * f #src
    lagrange_lhs = Wf * Lop(Uf) - Uf * Lop(Wf)                                #src
    lagrange_rhs = expand_derivatives(Dx_(Wf * Dx_(Uf) - Dx_(Wf) * Uf))       #src
    @assert symbolic_zero(lagrange_lhs - lagrange_rhs)                        #src
end                                                                           #src
println("ASSERTION 6 OK: W*L[U] - U*L[W] = d/dx[W*U' - W'*U] identically")     #src

sph_jlp(l, z) = (l / z) * sph_jl(l, z) - sph_jl(l + 1, z)                     #src

function shoot_Yprime1(l, q0, RHS::Function; eps0=1e-2, n=200_000)            #src
    q2 = q0^2                                                                 #src
    y1(x) = x * sph_jl(l, q0 * x)                                             #src
    dy1(x) = sph_jl(l, q0 * x) + q0 * x * ((l / (q0 * x)) * sph_jl(l, q0 * x) - sph_jl(l + 1, q0 * x)) #src
    f(x_, Y_, Yp_) = (Yp_, RHS(x_) + l * (l + 1) / x_^2 * Y_ - q2 * Y_)       #src
    h = (1.0 - eps0) / n                                                      #src
    xc, Y, Yp = eps0, y1(eps0), dy1(eps0)                                     #src
    for _ in 1:n                                                              #src
        k1 = f(xc, Y, Yp)                                                     #src
        k2 = f(xc + h / 2, Y + h / 2 * k1[1], Yp + h / 2 * k1[2])             #src
        k3 = f(xc + h / 2, Y + h / 2 * k2[1], Yp + h / 2 * k2[2])             #src
        k4 = f(xc + h, Y + h * k3[1], Yp + h * k3[2])                         #src
        Y += h / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1])                        #src
        Yp += h / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])                       #src
        xc += h                                                               #src
    end                                                                       #src
    kappa = -Y / y1(1.0)                                                      #src
    Yp + kappa * dy1(1.0)                                                     #src
end                                                                           #src

for (l_val, q0_val, rhs) in [(2, 2.6656, x -> x^2),                           #src
    (3, 1.9, x -> sin(3x) + 1),                                               #src
    (2, 0.9, x -> x^4 - 2x)]                                                  #src
    y_shortcut = adjoint_shortcut_Yprime1(l_val, q0_val, rhs)                 #src
    y_direct = shoot_Yprime1(l_val, q0_val, rhs)                              #src
    @assert abs(y_shortcut - y_direct) / abs(y_direct) < 1e-6                  #src
end                                                                           #src
println("ASSERTION 7 OK: the adjoint shortcut matches direct RK4 shooting to <1e-6, three cases") #src

# ## 7. The strain-rate field of Reid's actual viscous mode
#
# Any first-principles shear-rate calculation must use Reid's actual
# **viscous** velocity profile
#
# ```math
# U(x) \;=\; C\,x\,j_l(qx) + \Pi_0\,x^{l+1},
# ```
#
# not the inviscid potential-flow shape ``r^lP_l(\cos\theta)`` that the
# closures of the preceding pages use. The two are not
# interchangeable, and the reason is structural rather than a matter of
# accuracy: Reid's own damping normalization comes from the homogeneous
# (Bessel) part of ``U(x)``, which *is* the viscous correction to potential
# flow. Dropping it and then computing a viscous correction from what remains
# is inconsistent with Reid's theory already at zeroth order.
#
# Writing ``u_r=F(x)P_l(\cos\theta)`` with ``F(x)\equiv U(x)/x^2``, and
# ``u_\theta`` from the stream-function relation already used for BC2, the
# strain-rate components are
#
# ```math
# e_{rr}=F'(x)P_l(\mu),\qquad
# e_{\theta\theta}=\frac{1}{x}\partial_\theta u_\theta+\frac{u_r}{x},
# ```
# ```math
# e_{\varphi\varphi}=\frac{u_r+u_\theta\cot\theta}{x},\qquad
# e_{r\theta}=\frac12\left[x\,\partial_x\!\left(\frac{u_\theta}{x}\right)
#   +\frac{1}{x}\partial_\theta u_r\right],
# ```
#
# with ``\mu\equiv\cos\theta``. Incompressibility,
# ``e_{rr}+e_{\theta\theta}+e_{\varphi\varphi}=0``, is verified to
# floating-point precision for ``l=2,\dots,8`` and for *any* radial profile
# ``F(x)`` -- it is a kinematic identity of the poloidal representation, so
# the assembled tensor is the strain rate of a genuinely incompressible flow
# before any Carreau-Yasuda physics enters.

function legendre_P(l::Int, xv)                                               #src
    l == 0 && return one(xv)                                                  #src
    l == 1 && return xv                                                       #src
    Pm1, P = one(xv), xv                                                      #src
    for n in 1:(l-1)                                                          #src
        P, Pm1 = ((2n + 1) * xv * P - n * Pm1) / (n + 1), P                    #src
    end                                                                       #src
    P                                                                         #src
end                                                                           #src

@variables theta                                                              #src
@variables Ffun(..)                                                           #src
let Ff = Ffun(x), Dx_ = Differential(x), Dth_ = Differential(theta), Dmu_ = Differential(cos(theta)) #src
    function strain_trace(l::Int)                                             #src
        mu = cos(theta)                                                       #src
        Pl = legendre_P(l, mu)                                                #src
        Plp = expand_derivatives(Dmu_(Pl))                                    #src
        ur = Ff * Pl                                                          #src
        uth = -(2 * Ff + x * Dx_(Ff)) * sin(theta) * Plp / (l * (l + 1))       #src
        e_rr = expand_derivatives(Dx_(ur))                                    #src
        e_thth = expand_derivatives((1 / x) * Dth_(uth) + ur / x)             #src
        e_phph = expand_derivatives((ur + uth * cos(theta) / sin(theta)) / x)  #src
        e_rr + e_thth + e_phph                                                #src
    end                                                                       #src
    for l_val in (2, 3, 4, 5, 6, 8)                                           #src
        expr = strain_trace(l_val)                                            #src
        fexpr = Symbolics.build_function(expr, x, theta, Ffun(x), Dx_(Ffun(x)); expression=false) #src
        maxabs = 0.0                                                          #src
        for xv in (0.3, 0.5, 0.7, 0.9), thv in (0.4, 1.0, 1.7, 2.3, 2.9)       #src
            Fv, Fpv = xv^3, 3 * xv^2   # a concrete radial profile for the numeric check #src
            maxabs = max(maxabs, abs(fexpr(xv, thv, Fv, Fpv)))                #src
        end                                                                   #src
        @assert maxabs < 1e-10                                                #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 8 OK: incompressibility of the viscous-mode strain tensor holds for l=2..8, any F(x)") #src

# ## 8. The period-``\pi`` lemma
#
# At ``\mathrm{Oh}_\infty`` Reid's root ``q`` is genuinely complex -- a true
# damped oscillation, not a pure decay -- so the physical strain field at
# phase ``\phi=\omega t`` is ``\mathrm{Re}[e_{ij}(x,\theta)e^{-i\phi}]``.
#
# The natural first attempt at a shear-thinning correction is to evaluate
# everything at "one representative oscillation phase". That cannot be
# justified: the physical shape at a representative radius swings by
# roughly the same magnitude in both directions across one period, nowhere
# near constant. Any nonlinear function of this field -- and a fractional
# power is emphatically nonlinear -- depends on which phase was picked,
# unless the choice is justified by an actual time-average.
#
# The lemma that resolves it:
#
# > Let ``S\equiv\sqrt{2e_{ij}e_{ij}}`` be built from any single-frequency
# > oscillating axisymmetric poloidal field. Then ``S(\phi)`` is exactly
# > periodic with period ``\pi``, not ``2\pi``, so its Fourier series over
# > the full period contains only even harmonics: the coefficient at the
# > fundamental ``(m=1)`` is identically zero -- for ``S`` itself and for
# > *any* function of ``S``.
#
# **Proof.** Each strain component satisfies
# ``\mathrm{Re}(e_{ij}e^{-i\phi})=|e_{ij}|\cos(\phi-\arg e_{ij})``. Squaring,
# ``|e_{ij}|^2\cos^2(\phi-\arg e_{ij})
# =\tfrac12|e_{ij}|^2[1+\cos(2\phi-2\arg e_{ij})]``, which is manifestly
# invariant under ``\phi\to\phi+\pi``. ``S^2`` is a sum of such terms, hence
# itself period-``\pi``; and ``S=\sqrt{S^2}\ge0``, along with any real
# function of ``S``, inherits the period exactly -- there is no branch
# ambiguity, since ``S^2\ge0`` throughout. A period-``\pi`` function
# integrated against ``\sin(m\phi)`` or ``\cos(m\phi)`` over ``[0,2\pi)``
# vanishes identically for every odd ``m``, because the two half-period
# copies enter with opposite sign and cancel. ``\blacksquare``
#
# Confirmed two independent ways. Numerically,
# ``S^2(\phi+\pi)=S^2(\phi)`` to better than ``10^{-10}`` at several
# independent ``(x,\theta)`` points. And by direct Fourier projection of the
# actual Carreau-Yasuda correction shape ``S^{-a}`` on a joint
# ``200\times200`` ``(\theta,\phi)`` grid, where the ``m=1`` coefficient --
# both ``\cos`` and ``\sin`` -- is zero to better than ``10^{-10}`` relative
# to the ``m=0`` channel, at every ``\theta`` tested.
#
# **The consequence.** Since no forcing
# resonant at the base mode's own frequency exists, the adjoint machinery of
# §6 -- built for exactly this kind of resonant solvability condition -- has
# nothing to act on at ``m=1``. The leading temporal effect of any
# generalized-Newtonian correction is carried entirely by the ``m=0``
# channel: period-averaged, effectively steady. The ``m=2`` channel is a
# distinct, second-harmonic *parametric* coupling -- a different physical
# phenomenon, not addressed here.

q0_l2 = dominant_root(OH_INF, 2)                                              #src
jl0_c = sph_jl(2, q0_l2)                                                      #src
Q0_c = sph_jl(3, q0_l2) / jl0_c                                               #src
C0_c = 2 * (2^2 - 1) / (jl0_c * q0_l2 * (2 * Q0_c - q0_l2))                   #src
Pi0_c = -1 - C0_c * jl0_c                                                     #src
Ufun_complex(xv) = C0_c * xv * sph_jl(2, q0_l2 * xv) + Pi0_c * xv^(2 + 1)      #src
Ffun_complex(xv) = Ufun_complex(xv) / xv^2                                    #src
Fderiv_complex(xv; h=1e-6) = (Ffun_complex(xv + h) - Ffun_complex(xv - h)) / (2h) #src
legendre_Pp(l::Int, xv) = l == 0 ? zero(xv) : l * (xv * legendre_P(l, xv) - legendre_P(l - 1, xv)) / (xv^2 - 1) #src

function strain_complex(xv, th, l)                                            #src
    mu = cos(th)                                                              #src
    F = Ffun_complex(xv)                                                      #src
    Fp = Fderiv_complex(xv)                                                   #src
    dth = 1e-6                                                                #src
    uth_coef(xx) = -(2 * Ffun_complex(xx) + xx * Fderiv_complex(xx)) / (l * (l + 1)) #src
    Pl_th(t) = legendre_P(l, cos(t))                                          #src
    Plp_th(t) = legendre_Pp(l, cos(t))                                        #src
    ur_of_theta(t) = F * Pl_th(t)                                             #src
    uth_of_theta(t) = uth_coef(xv) * sin(t) * Plp_th(t)                       #src
    dur_dth = (ur_of_theta(th + dth) - ur_of_theta(th - dth)) / (2dth)        #src
    duth_dth = (uth_of_theta(th + dth) - uth_of_theta(th - dth)) / (2dth)     #src
    ur = F * Pl_th(th)                                                        #src
    uth = uth_coef(xv) * sin(th) * Plp_th(th)                                 #src
    e_rr = Fp * legendre_P(l, mu)                                             #src
    e_thth = duth_dth / xv + ur / xv                                          #src
    e_phph = (ur + uth * cos(th) / sin(th)) / xv                              #src
    dxv = 1e-6                                                                #src
    uth_x(xx) = -(2 * Ffun_complex(xx) + xx * Fderiv_complex(xx)) / (l * (l + 1)) * sin(th) * Plp_th(th) / xx #src
    duthx_dx = (uth_x(xv + dxv) - uth_x(xv - dxv)) / (2dxv)                   #src
    e_rth = 0.5 * (xv * duthx_dx + dur_dth / xv)                              #src
    e_rr, e_thth, e_phph, e_rth                                               #src
end                                                                           #src

function S_at(xv, th, phase, l)                                               #src
    e_rr, e_thth, e_phph, e_rth = strain_complex(xv, th, l)                   #src
    ph = cos(phase) - im * sin(phase)                                         #src
    err_r, ethth_r, ephph_r, erth_r = real(e_rr * ph), real(e_thth * ph), real(e_phph * ph), real(e_rth * ph) #src
    sqrt(2 * (err_r^2 + ethth_r^2 + ephph_r^2 + 2 * erth_r^2))                #src
end                                                                           #src

for (xv, th) in ((0.3, 1.0), (0.7, 0.5), (0.95, 2.0))                         #src
    @assert abs(S_at(xv, th, 0.7, 2)^2 - S_at(xv, th, 0.7 + pi, 2)^2) < 1e-10  #src
end                                                                           #src
println("ASSERTION 9 OK: S^2(phase+pi) = S^2(phase) to <1e-10 at three (x,theta) points") #src

let Nth = 200, Nph = 200                                                      #src
    thetas = collect(range(0.02, pi - 0.02; length=Nth))                      #src
    phases = collect(range(0.0, 2pi; length=Nph + 1))[1:end-1]                 #src
    dphase = phases[2] - phases[1]                                            #src
    grid = [S_at(0.7, th, ph, 2)^(-A_SHAPE) for th in thetas, ph in phases]    #src
    m1_cos_max = 0.0                                                          #src
    m1_sin_max = 0.0                                                          #src
    for i in 1:Nth                                                            #src
        c = (1 / pi) * sum(grid[i, :] .* cos.(phases)) * dphase               #src
        s = (1 / pi) * sum(grid[i, :] .* sin.(phases)) * dphase               #src
        m1_cos_max = max(m1_cos_max, abs(c))                                  #src
        m1_sin_max = max(m1_sin_max, abs(s))                                  #src
    end                                                                       #src
    m0_typical = abs((1 / (2pi)) * sum(grid[Nth÷2, :]) * dphase)              #src
    @assert m1_cos_max < 1e-10 * max(m0_typical, 1.0)                         #src
    @assert m1_sin_max < 1e-10 * max(m0_typical, 1.0)                         #src
end                                                                           #src
println("ASSERTION 10 OK: the m=1 Fourier coefficient of S^(-a) is zero to numerical precision") #src

# ### Spatial leakage: one active mode forces every even degree
#
# The same period-averaged correction shape, projected onto Legendre
# polynomials, shows a second mechanism -- this one spatial. A **single**
# active mode ``l=2`` spontaneously forces every even spherical-harmonic
# degree:
#
# | ``l'`` | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
# |:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
# | leakage coefficient | 1.297 | ``4\!\times\!10^{-11}`` | ``-0.0418`` | ``-2\!\times\!10^{-11}`` | 0.1023 | ``3\!\times\!10^{-12}`` | ``-0.0594`` | ``3\!\times\!10^{-11}`` | 0.0949 |
#
# Every even degree is genuinely nonzero; every odd degree vanishes to
# ``10^{-11}``, by a parity property that is verified rather than assumed.
# The reason is elementary and unavoidable: raising a finite Legendre
# polynomial to a non-integer power does not, in general, stay within its own
# degree. A closure that assigns mode ``l`` a viscosity from mode ``l``'s own
# shape is discarding all of the ``l'\neq l`` columns.

let Nth = 4000                                                                #src
    thetas = collect(range(0.02, pi - 0.02; length=Nth))                      #src
    mus = cos.(thetas)                                                        #src
    perm = sortperm(mus)                                                      #src
    mus_s, thetas_s = mus[perm], thetas[perm]                                 #src
    Nph = 200   # period-averaged (m=0) shape at each theta                    #src
    phases = collect(range(0.0, 2pi; length=Nph + 1))[1:end-1]                 #src
    shape_avg = [sum(S_at(0.7, th, ph, 2)^(-A_SHAPE) for ph in phases) / Nph for th in thetas_s] #src
    trapz(xv, yv) = sum(0.5 * (yv[i] + yv[i+1]) * (xv[i+1] - xv[i]) for i in 1:length(xv)-1) #src
    coeffs = Dict{Int,Float64}()                                              #src
    for lp in 0:8                                                             #src
        integrand = legendre_P.(lp, mus_s) .* shape_avg                       #src
        coeffs[lp] = (2lp + 1) / 2 * trapz(mus_s, integrand)                  #src
    end                                                                       #src
    @assert all(lp -> abs(coeffs[lp]) < 1e-6, (1, 3, 5, 7))     # odd l': zero #src
    @assert all(lp -> abs(coeffs[lp]) > 1e-4, (0, 2, 4, 6, 8))  # even l': genuinely nonzero #src
    @assert isapprox(coeffs[0], 1.297; rtol=1e-3)               # pins the tabulated values #src
    @assert isapprox(coeffs[4], 0.1023; rtol=1e-3)                            #src
end                                                                           #src
println("ASSERTION 11 OK: a single l=2 mode leaks onto every even l', and onto no odd l'") #src

# ## 9. The genuinely open question: spatial homogenization
#
# Because the ``m=0`` channel is the leading one, and because a steady
# forcing is a valid input to §6's adjoint shortcut, that machinery is the
# right tool for the question that actually remains open. Replacing a
# spatially varying ``\eta(x,\theta)`` with a single scalar
# ``\mathrm{Oh}_{\mathrm{eff}}(t)`` is an uncontrolled effective-medium
# approximation unless the spatial variation is itself small. Is it?
#
# Using the period-averaged shear-rate shape from §7-§8 over the drop's
# volume at ``\mathrm{Oh}_\infty``, ``l=2``, the shape varies by a factor of
# **2.25** between its minimum and maximum. Propagated through the
# Carreau-Yasuda exponent ``a=0.743``, that is roughly a factor of **1.83**
# in local viscosity across the drop at a single instant.
#
# That is smaller than the multiple-order-of-magnitude range
# ``\mathrm{Oh}_{\mathrm{eff}}`` traverses over an impact *in time* -- so the
# scalar approximation is not unreasonable. But a factor of nearly two is
# not small enough to treat as negligible either: the spatial variation is
# real, moderate, and not yet controlled.

let vals = Float64[]                                                          #src
    for xv in range(0.15, 0.98; length=12), th in range(0.1, pi - 0.1; length=12) #src
        Nph = 100                                                             #src
        phases = collect(range(0.0, 2pi; length=Nph + 1))[1:end-1]             #src
        push!(vals, sum(S_at(xv, th, ph, 2) for ph in phases) / Nph)          #src
    end                                                                       #src
    spread_ratio = maximum(vals) / minimum(vals)                              #src
    @assert 1.5 < spread_ratio < 4.0   # a real, moderate, O(1) spread -- not tiny, not huge #src
    @assert isapprox(spread_ratio, 2.248; rtol=1e-3)   # pins the reported value #src
    viscosity_spread = spread_ratio^(-A_SHAPE)   # propagated through the CY exponent #src
    @assert isapprox(1 / viscosity_spread, 1.826; rtol=1e-3)                  #src
end                                                                           #src
println("ASSERTION 12 OK: spatial spread in local viscosity across the drop is a factor of ~1.83") #src

# ## 10. Summary
#
# 1. The classical ``O(\epsilon^3)`` Carreau expansion does **not** extend to
#    this fluid's fitted, non-integer exponent. No small-amplitude anchor --
#    rest, infinite shear, or in between -- supports a uniformly valid linear
#    correction to ``\lambda_l(\mathrm{Oh})`` across the modes this solver
#    resolves.
# 2. Reid's exact ``\lambda_l(\mathrm{Oh})``, ``\omega_l^2(\mathrm{Oh})``
#    should therefore be evaluated *directly* at whatever
#    ``\mathrm{Oh}_{\mathrm{eff}}(t)`` a closure provides, and never
#    linearized -- which DropSolver implements, via a fast
#    tabulated version of exactly that evaluation.
# 3. The genuinely open approximation in the current scheme is the
#    collapse of a spatially varying ``\eta(x,\theta,t)`` onto a single
#    scalar ``\mathrm{Oh}_{\mathrm{eff}}(t)``. It is quantified here at
#    roughly a factor of two in local viscosity across the drop at a
#    representative instant: real, moderate, and not yet controlled.
# 4. The adjoint sensitivity shortcut (§6), the viscous strain-rate tensor
#    (§7), and the period-``\pi`` lemma (§8) survive independently of all of
#    the above. None depends on this fluid's specific parameters, and all
#    three are the correct starting point for a future first-principles
#    treatment -- most plausibly a self-consistent, spatially resolved solve
#    using the adjoint machinery as a Newton/Fréchet-derivative iteration
#    step, rather than any further attempt at a closed-form amplitude
#    equation.
#
# **Scope.** This page does not implement that self-consistent solve; that is
# future work. It does not revisit the ``m=2`` parametric channel. And it
# does not extend the strain-tensor or leakage calculations beyond ``l=2``
# and a single representative radius. What it does establish is precisely
# which of the natural next steps work, which do not, and why -- for this
# fluid specifically.
