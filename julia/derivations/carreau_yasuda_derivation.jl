# # Weakly Nonlinear Carreau-Yasuda Drop (legacy: superseded by the multi-mode model)
#
# !!! warning "Status"
#     This page derives `julia/src/st_extension.jl`, the *perturbative*
#     shear-thinning correction. The real validation pipeline no longer uses
#     it: for the fitted validation fluid its small-parameter assumptions
#     fail outright, and the damping multiplier can go negative --
#     see `carreau_yasuda_nonperturbative_derivation.jl` for the
#     quantification and `carreau_yasuda_multimode_derivation.jl` for the
#     model that replaced it. What remains valuable here is the machinery:
#     the geometric integral ``\Gamma_l^{(a)}``, the secular factor
#     ``C(a)``, and their finite-``\mathrm{Oh}`` versions, all of which the
#     later work reuses.
#
# The derivation itself: a perturbation expansion for a shear-thinning drop
# oscillating on a flat substrate, using the Carreau-Yasuda constitutive law
# -- Carreau with a free shape exponent ``a``, rather than ``a`` hardcoded to
# 2. It follows the same non-dimensionalisation as the Newtonian (Reid 1960)
# and Oldroyd-B analyses.
#
# **Why Carreau-Yasuda rather than plain Carreau?** The standard Carreau
# model's nonlinearity is always built from ``(\lambda_c\dot\gamma)^2``, so
# the *shape* of the transition between the zero-shear plateau and the
# power-law region is fixed. Carreau-Yasuda adds one parameter, ``a``,
# controlling that shape; ``a=2`` recovers standard Carreau exactly. Real
# experimental shear-thinning fits -- including Cross-model
# characterizations, which convert directly to an equivalent Carreau-Yasuda
# ``a`` -- commonly need ``a\neq2`` to fit well.
#
# Every section derives one piece of the chain and then checks the result
# against an independent computation. Two checks recur and are worth naming
# up front. **Reduction at ``a=2``:** every generalized result must collapse
# onto the older, separately validated Carreau-only result, which is a strong
# constraint because the two are computed by completely different routes.
# **Newtonian limits:** setting ``n=1``, ``\lambda_c=0``, or
# ``\varepsilon_{ST}=0`` must return the unmodified theory for *any* ``a``.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``x=r/R`` | dimensionless radial coordinate |
# | ``\epsilon`` | oscillation amplitude, the small parameter |
# | ``\sigma=-\gamma+i\omega`` | complex modal frequency (``\gamma>0`` = decay) |
# | ``q``, ``\sigma=q^2\mathrm{Oh}`` | viscous wavenumber |
# | ``\alpha^2=\sigma_{l;0}/\mathrm{Oh}``, ``\sigma_{l;0}=\sqrt{l(l-1)(l+2)}`` | |
# | ``\mathrm{Oh}=\mu/\sqrt{\rho T_1 R}`` | Ohnesorge number |
# | ``U(x)`` | radial velocity eigenfunction (Reid) |
# | ``b_l(t)`` | dimensionless amplitude of mode ``l`` (the solver's ``A_l``) |
# | ``P_l(\cos\theta)`` | Legendre polynomial of degree ``l`` |
# | ``a`` | the Carreau-Yasuda **shape exponent**, a model parameter -- distinct from the mode amplitude ``a(t)`` of §7-§8; the shape exponent is always `a_shape` in code |

using Symbolics
using QuadGK
using SpecialFunctions
using LinearAlgebra
using DropSolver

nth_derivative_at_zero(expr, var, n) = begin                                  #src
    d = expr                                                                  #src
    for _ in 1:n                                                              #src
        d = Symbolics.derivative(d, var)                                      #src
    end                                                                       #src
    simplify(substitute(d, Dict(var => 0)))                                   #src
end                                                                           #src

## Robust symbolic-equality check: Symbolics.jl's `simplify` does not always  #src
## collapse an algebraically-zero difference to the literal 0 (a known        #src
## limitation, weaker than sympy's simplifier here).  Rather than trust       #src
## `isequal(simplify(a-b), 0)` alone, evaluate both expressions at several    #src
## concrete numeric points for every free variable -- an independent, purely  #src
## numeric check that is immune to this limitation.                           #src
const _FIXED_TEST_VALUES = (0.31, 0.57, 1.13, 1.94, 2.71)                     #src

function numerically_equal(expr1, expr2, test_points::Dict=Dict())            #src
    vars = collect(union(Symbolics.get_variables(expr1), Symbolics.get_variables(expr2))) #src
    isempty(vars) && return isequal(simplify(expr1 - expr2), 0)               #src
    f1 = Symbolics.build_function(expr1, vars...; expression=false)           #src
    f2 = Symbolics.build_function(expr2, vars...; expression=false)           #src
    for trial in 1:length(_FIXED_TEST_VALUES)                                 #src
        vals = [haskey(test_points, v) ? test_points[v] :                     #src
                _FIXED_TEST_VALUES[mod1(trial + i, length(_FIXED_TEST_VALUES))] #src
                for (i, v) in enumerate(vars)]                                #src
        isapprox(f1(vals...), f2(vals...); atol=1e-8, rtol=1e-6) || return false #src
    end                                                                       #src
    true                                                                      #src
end                                                                           #src

# ## 1. The constitutive law and its leading correction
#
# The Carreau-Yasuda effective viscosity is
#
# ```math
# \mu_{\mathrm{eff}}(\dot\gamma)
#   = \mu_0\left[1+(\lambda_c\dot\gamma)^{a}\right]^{(n-1)/a},
# ```
#
# with ``n\in(0,1]`` the power-law index (``n=1`` is Newtonian) and ``a>0``
# the shape exponent (``a=2`` is standard Carreau). Here ``\mu_\infty=0``,
# matching `julia/src/st_extension.jl`, which has never carried a separate
# infinite-shear viscosity. That is not a loss of generality for this
# linearized solver: a nonzero ``\mu_\infty`` from a real fluid
# characterization folds entirely into a redefined ``\varepsilon_{ST}``.
#
# Writing ``\dot\gamma=\epsilon\hat\gamma`` and expanding for small
# ``\epsilon``, the leading correction is
#
# ```math
# \boxed{\;
# \mu_{\mathrm{eff}} \;\approx\; \mu_0\left[1
#   + \frac{n-1}{a}\,(\lambda_c\hat\gamma)^{a}\,\epsilon^{a}\right]
# \;}
# ```
#
# so the correction enters at fractional order ``\epsilon^{a}`` with
# coefficient ``p=(n-1)/a``.
#
# That expansion is not asserted -- it is verified by repeated symbolic
# differentiation, Taylor's theorem applied directly, at concrete integer
# ``a=2,3,4`` (Symbolics.jl has no `series()` the way sympy does). Two things
# are checked at each ``a``: that the ``\epsilon^{a}`` coefficient matches
# the hand-built form, and that *every* intermediate power ``\epsilon^{1}``
# through ``\epsilon^{a-1}`` vanishes identically. The second is the load
# bearing one -- a surviving lower-order term would mean the correction
# enters earlier than claimed and the whole ordering of the expansion is
# wrong.
#
# The limits then follow: ``n=1`` and ``\lambda_c=0`` each return ``\mu_0``
# exactly, for any ``a``. And at ``a=2`` the correction reduces to
# ``\tfrac12\mu_0(n-1)(\lambda_c\hat\gamma\epsilon)^2``, i.e.
# ``\varepsilon_{ST}=(1-n)/2`` -- exactly the pre-generalization Carreau
# result.

@variables mu_0 lambda_c eps gammadot_hat a_shape n_idx                       #src

x_small = (lambda_c * gammadot_hat)^a_shape * eps^a_shape                     #src
p_exp = (n_idx - 1) / a_shape                                                 #src
mu_CY_leading = mu_0 * (1 + p_exp * x_small)                                  #src

for a_val in (2, 3, 4)                                                        #src
    mu_CY_concrete = mu_0 * (1 + (lambda_c * eps * gammadot_hat)^a_val)^((n_idx - 1) / a_val) #src
    coeff_exact = nth_derivative_at_zero(mu_CY_concrete, eps, a_val) / factorial(a_val) #src
    hand_built = substitute(mu_CY_leading, Dict(a_shape => a_val))            #src
    coeff_hand = Symbolics.coeff(expand(hand_built - mu_0), eps^a_val)        #src
    @assert numerically_equal(coeff_exact, coeff_hand)                        #src
    for lower_power in 1:(a_val-1)                                            #src
        lower = nth_derivative_at_zero(mu_CY_concrete, eps, lower_power)      #src
        @assert isequal(lower, 0) "a=$a_val: unexpected nonzero O(eps^$lower_power) term: $lower" #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 1 OK: leading-order expansion verified by repeated differentiation at a=2,3,4") #src

@assert numerically_equal(substitute(mu_CY_leading, Dict(n_idx => 1)), mu_0)   #src
println("ASSERTION 2 OK: n=1 recovers Newtonian (zero correction) for any a")  #src
@assert numerically_equal(substitute(mu_CY_leading, Dict(lambda_c => 0)), mu_0) #src
println("ASSERTION 3 OK: lambda_c=0 recovers Newtonian for any a")             #src

eps_ST_CY = (1 - n_idx) / a_shape                                             #src
correction_a2 = simplify(substitute(mu_CY_leading - mu_0, Dict(a_shape => 2))) #src
target_original = eps .^ 2 .* gammadot_hat .^ 2 .* lambda_c .^ 2 .* mu_0 .* (n_idx - 1) ./ 2 #src
@assert numerically_equal(correction_a2, target_original)                     #src
println("ASSERTION 4 OK: at a=2 the correction reduces to eps_ST = (1-n)/2")   #src

# ## 2. Reid's velocity field, and the boundary conditions that fix it
#
# Everything in §2-§5 is Newtonian base-flow machinery. It does not depend on
# the shear-thinning law at all -- these are exact properties of the
# Newtonian base flow, identical whether the higher-order correction is
# Carreau, Carreau-Yasuda, Cross, or absent.
#
# Reid's radial velocity eigenfunction solves
#
# ```math
# \left[\frac{d^2}{dx^2}-\frac{l(l+1)}{x^2}+q^2\right]U(x)
#   \;=\; q^2\,\Pi_0\,x^{l+1},
# \qquad
# U(x) = C\,x\,j_l(qx) + \Pi_0\,x^{l+1}.
# ```
#
# Spherical Bessel functions of integer order reduce to elementary
# functions -- for ``l=2``,
# ``j_2(z)=(3/z^3-1/z)\sin z-(3/z^2)\cos z`` -- so the ODE can be verified
# with ordinary calculus, no special-function differentiation rules required.
# The residual collapses to symbolic zero.
#
# The two boundary conditions fix ``C`` and ``\Pi_0`` through a ``2\times2``
# linear system. Solving it reproduces Reid's own equation (38),
#
# ```math
# C = \frac{2(l-1)(l+1)}{(2l-q^2)\,j_l - 2q\,j_l'},
# ```
#
# and BC1 in the form ``C j_l+\Pi_0=-1``, i.e. ``U(1)=-1``, holds
# identically. If either failed, the eigenfunction being differentiated for
# the strain rate would not be a solution of Reid's problem at all.

@variables x_sym q_sym C_sym Pi0_sym                                          #src
j2(z) = (3 / z^3 - 1 / z) * sin(z) - (3 / z^2) * cos(z)                       #src

l_test = 2                                                                    #src
U_expr = C_sym * x_sym * j2(q_sym * x_sym) + Pi0_sym * x_sym^(l_test + 1)      #src
U_expr_p = Symbolics.derivative(U_expr, x_sym)                                #src
U_expr_pp = Symbolics.derivative(U_expr_p, x_sym)                             #src
ODE_lhs = U_expr_pp - l_test * (l_test + 1) / x_sym^2 * U_expr + q_sym^2 * U_expr #src
ODE_rhs = q_sym^2 * Pi0_sym * x_sym^(l_test + 1)                              #src
residual = simplify(ODE_lhs - ODE_rhs)                                        #src

## Verify both symbolically AND numerically: a concrete-value check catches   #src
## cases where simplify() cannot fully collapse a transcendental expression   #src
## that IS actually zero.                                                     #src
symbolic_zero = isequal(residual, 0)                                          #src
numeric_zero = let                                                            #src
    f = Symbolics.build_function(residual, x_sym, q_sym, C_sym, Pi0_sym; expression=false) #src
    all(abs(f(xv, qv, 1.3, -0.7)) < 1e-9 for xv in (0.2, 0.5, 0.9), qv in (1.1, 2.7, 4.3)) #src
end                                                                           #src
@assert symbolic_zero || numeric_zero "Reid ODE not satisfied: $residual"     #src
println("ASSERTION 5 OK: U(x) = C*x*j_2(qx) + Pi0*x^3 satisfies Reid's ODE")   #src

@variables jl_sym jlp_sym l_s                                                 #src
A_mat = [jl_sym 1;                                                            #src
    -q_sym^2*jl_sym+2*(l_s^2 + l_s - 1)*jl_sym-2*q_sym*jlp_sym 2*(l_s^2 - 1)]  #src
b_vec = [-1, 0]                                                               #src
sol = Symbolics.simplify.(A_mat \ b_vec)                                      #src
C_sol, Pi0_sol = sol[1], sol[2]                                               #src

@assert isequal(simplify(C_sol * jl_sym + Pi0_sol + 1), 0)                     #src
println("ASSERTION 6 OK: C*j_l + Pi0 = -1 (BC1 verified for Pi0)")             #src
C_reid = 2 * (l_s - 1) * (l_s + 1) / ((2 * l_s - q_sym^2) * jl_sym - 2 * q_sym * jlp_sym) #src
@assert numerically_equal(C_sol, C_reid)                                      #src
println("ASSERTION 7 OK: the BC solution for C matches Reid eq. (38)")         #src
@assert numerically_equal(C_sol * jl_sym + Pi0_sol + 1, 0)                     #src
println("ASSERTION 8 OK: BC1 satisfied, U(1) = -1")                           #src

# ## 3. Incompressibility of the velocity field
#
# With ``u_r=U(x)P_l`` and ``u_\theta=V(x)\,dP_l/d\theta``, continuity fixes
#
# ```math
# V(x) = \frac{(x^2U)'}{l(l+1)\,x},
# ```
#
# and the divergence ``(x^2U)'/x^2-l(l+1)V/x`` then vanishes identically for
# an arbitrary ``U(x)`` -- a kinematic identity, not a property of Reid's
# particular solution.

@variables theta_sym                                                          #src
l_num = 2                                                                     #src
V_from_U(Uexpr) = Symbolics.derivative(x_sym^2 * Uexpr, x_sym) / (l_num * (l_num + 1) * x_sym) #src
@variables Utest(x_sym)                                                       #src
Vtest = V_from_U(Utest)                                                       #src
div_check = Symbolics.derivative(x_sym^2 * Utest, x_sym) / x_sym^2 - l_num * (l_num + 1) * Vtest / x_sym #src
@assert numerically_equal(div_check, 0)                                       #src
println("ASSERTION 9 OK: u_r = U*P_l, u_theta = V*(dP_l/dtheta) is incompressible") #src

# ## 4. The strain-rate tensor
#
# In the axisymmetric spherical representation, with ``f\equiv U`` and
# ``g\equiv V``,
#
# ```math
# e_{rr}=f'P_l,\qquad
# e_{r\theta}=\tfrac12\!\left(g'-\frac{g}{x}+\frac{f}{x}\right)\frac{dP_l}{d\theta},
# ```
# ```math
# e_{\theta\theta}=\frac{f}{x}P_l+\frac{g}{x}\frac{d^2P_l}{d\theta^2},\qquad
# e_{\varphi\varphi}=\frac{f}{x}P_l+\frac{g\cot\theta}{x}\frac{dP_l}{d\theta},
# ```
#
# and the scalar shear-rate invariant is
# ``\dot\gamma^2=2(e_{rr}^2+e_{\theta\theta}^2+e_{\varphi\varphi}^2+2e_{r\theta}^2)``.

Pl_expr = (3 * cos(theta_sym)^2 - 1) / 2                                      #src
dPl_dth = Symbolics.derivative(Pl_expr, theta_sym)                            #src
d2Pl_dth = Symbolics.derivative(dPl_dth, theta_sym)                           #src

@variables f_sym(x_sym) g_sym(x_sym)                                          #src
fxp = Symbolics.derivative(f_sym, x_sym)                                      #src
gxp = Symbolics.derivative(g_sym, x_sym)                                      #src
e_rr = fxp * Pl_expr                                                          #src
e_rth = (1 // 2) * (gxp - g_sym / x_sym + f_sym / x_sym) * dPl_dth            #src
e_thth = (f_sym / x_sym) * Pl_expr + (g_sym / x_sym) * d2Pl_dth               #src
e_phph = (f_sym / x_sym) * Pl_expr + (g_sym * cos(theta_sym) / sin(theta_sym) / x_sym) * dPl_dth #src
gdot_sq = 2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2)                    #src

# ## 5. Angular integration
#
# The two Legendre-orthogonality integrals the projection needs are
#
# ```math
# \int_0^\pi P_l^2\sin\theta\,d\theta=\frac{2}{2l+1},
# \qquad
# \int_0^\pi\left(\frac{dP_l}{d\theta}\right)^2\sin\theta\,d\theta
#   =\frac{2l(l+1)}{2l+1},
# ```
#
# which at ``l=2`` give ``0.4`` and ``2.4``. Since Symbolics.jl has no
# general symbolic `integrate`, both are confirmed to ``10^{-10}`` by QuadGK
# quadrature of the explicit ``l=2`` polynomials -- a fully rigorous numeric
# check of known closed forms.

Pl_f = Symbolics.build_function(Pl_expr, theta_sym; expression=false)         #src
dPl_f = Symbolics.build_function(dPl_dth, theta_sym; expression=false)        #src
I_Pl2, _ = quadgk(th -> Pl_f(th)^2 * sin(th), 0, pi)                          #src
I_dPl2, _ = quadgk(th -> dPl_f(th)^2 * sin(th), 0, pi)                        #src
l_check = 2                                                                   #src
@assert abs(I_Pl2 - 2 / (2 * l_check + 1)) < 1e-10                            #src
println("ASSERTION 10 OK: int P_2^2 sin(theta) dtheta = 2/(2l+1) = 0.4")       #src
@assert abs(I_dPl2 - 2 * l_check * (l_check + 1) / (2 * l_check + 1)) < 1e-10  #src
println("ASSERTION 11 OK: int (dP_2/dtheta)^2 sin(theta) dtheta = 2l(l+1)/(2l+1) = 2.4") #src

# ## 6. The correction integral ``\Gamma_l^{(a)}``
#
# The geometric factor that carries the shear-thinning correction into the
# mode equation is
#
# ```math
# \Gamma_l^{(a)}
#  = \frac{\displaystyle\int_0^1\!\!\left(\int_0^\pi
#      \dot\gamma^{\,a+2}\sin\theta\,d\theta\right)x^2\,dx}
#         {N_l^{(a+2)/2}} .
# ```
#
# The normalization generalizes the old ``a=2``-hardcoded ``N_l^2``: since
# ``\dot\gamma^{a+2}\sim U^{a+2}`` while ``N_l\sim U^2``, the exponent
# ``(a+2)/2`` on ``N_l`` is what keeps ``\Gamma_l^{(a)}`` dimensionless in
# powers of ``U`` for every ``a``. At ``a=2`` it is the old convention
# unchanged.
#
# ### The inviscid limit factorizes exactly
#
# In the inviscid limit ``U(x)=\Pi_0x^{l+1}=-x^{l+1}``, every strain
# component is proportional to ``x^{l}``, so
#
# ```math
# \dot\gamma^2(x,\theta) = x^{2l}\,H(\theta) \quad\text{exactly},
# \qquad
# H(\theta) = 3\cos^4\theta + 11\cos^2\theta + 13
# \quad (l=2),
# ```
#
# and -- this is the part that makes a general exponent tractable -- the
# factorization survives being raised to *any* power, not just squared. So
# the radial and angular integrals separate for every ``a``, with
# ``N_l=\int_0^1U^2x^2dx=1/9`` at ``l=2``.

l_val = 2                                                                     #src
Pi0_inviscid = -1                                                             #src
f_inviscid = Pi0_inviscid * x_sym^(l_val + 1)                                 #src
fp_inviscid = Symbolics.derivative(f_inviscid, x_sym)                         #src
@assert isequal(substitute(f_inviscid, Dict(x_sym => 1)), -1)                  #src
println("ASSERTION 12 OK: BC1 satisfied in the inviscid limit, U(1) = -1")     #src

g_inv_l2 = simplify(Symbolics.derivative(x_sym^2 * f_inviscid, x_sym) / (l_val * (l_val + 1) * x_sym)) #src
gp_inv_l2 = Symbolics.derivative(g_inv_l2, x_sym)                             #src
e_rr_inv = fp_inviscid * Pl_expr                                              #src
e_rth_inv = (1 // 2) * (gp_inv_l2 - g_inv_l2 / x_sym + f_inviscid / x_sym) * dPl_dth #src
e_thth_inv = (f_inviscid / x_sym) * Pl_expr + (g_inv_l2 / x_sym) * d2Pl_dth    #src
e_phph_inv = (f_inviscid / x_sym) * Pl_expr + (g_inv_l2 * cos(theta_sym) / (sin(theta_sym) * x_sym)) * dPl_dth #src
gdot_sq_inv = 2 * (e_rr_inv^2 + e_thth_inv^2 + e_phph_inv^2 + 2 * e_rth_inv^2) #src
H_theta_sym = simplify(substitute(gdot_sq_inv, Dict(x_sym => 1)))             #src

factorization_residual = simplify(gdot_sq_inv - x_sym^4 * H_theta_sym)        #src
symbolic_zero_13 = isequal(factorization_residual, 0)                         #src
numeric_zero_13 = let                                                         #src
    f = Symbolics.build_function(factorization_residual, x_sym, theta_sym; expression=false) #src
    all(abs(f(xv, thv)) < 1e-9 for xv in (0.1, 0.4, 0.7, 0.99), thv in (0.3, 1.1, 2.0, 2.9)) #src
end                                                                           #src
@assert symbolic_zero_13 || numeric_zero_13 "Factorization failed: $factorization_residual" #src
println("ASSERTION 13 OK: gdot^2 = x^(2l) * H(theta) exactly (inviscid, l=2)")  #src

# The angular shape and the resulting integral, as actually evaluated:

H_theta(th) = 3 * cos(th)^4 + 11 * cos(th)^2 + 13

"""
    Gamma_2_inviscid(a)

``\\Gamma_2^{(a)}`` in the inviscid limit: the angular integral of
``H(\\theta)^{(a+2)/2}`` divided by the radial factor ``l(a+2)+3`` and the
normalization ``N_2^{(a+2)/2}=(1/9)^{(a+2)/2}``.
"""
function Gamma_2_inviscid(a)
    p = (a + 2) / 2
    ang, _ = quadgk(th -> H_theta(th)^p * sin(th), 0, pi)
    (ang / (2 * (a + 2) + 3)) / (1 / 9)^p
end

# The symbolic ``H(\theta)`` from the strain tensor agrees with that closed
# polynomial to ``10^{-10}`` across the whole range of ``\theta``, and at
# ``a=2`` the integral reproduces the old Carreau-only notebook's
# sympy-exact value ``\Gamma_2=1783566/385=4632.6390`` to nine digits.
# Growth with ``a`` is steep, which is the reason the exponent cannot be
# treated as a minor detail:
#
# | ``a`` | 0.5 | 1 | 1.5 | 2 | 3 | 4 |
# |:--|:--|:--|:--|:--|:--|:--|
# | ``\Gamma_2^{(a)}`` | 138.3 | 439.2 | 1417 | 4632.6 | 51233 | 588151 |

H_f = Symbolics.build_function(H_theta_sym, theta_sym; expression=false)      #src
@assert all(abs(H_f(th) - H_theta(th)) < 1e-10 for th in range(0.01, pi - 0.01, length=20)) #src
println("ASSERTION 14 OK: H(theta) = 3cos^4 + 11cos^2 + 13 confirmed")         #src

Gamma2_exact_a2 = 1783566 / 385                                               #src
@assert abs(Gamma_2_inviscid(2.0) - Gamma2_exact_a2) / Gamma2_exact_a2 < 1e-9  #src
println("ASSERTION 15 OK: Gamma_2^(a=2) reproduces the exact 1783566/385")     #src
let vals = Dict(a => Gamma_2_inviscid(a) for a in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)) #src
    @assert isapprox(vals[0.5], 138.309; rtol=1e-4)   # pins the tabulated row  #src
    @assert isapprox(vals[4.0], 588150.79; rtol=1e-4)                         #src
    @assert issorted([vals[a] for a in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)])        #src
end                                                                           #src

# ### General ``l``, still inviscid
#
# Repeating the construction for ``l=2,3,4,5`` (with ``N_l=1/(2l+5)``, and
# the ``\theta``-derivatives of ``P_l`` taken by finite differences):
#
# | ``\Gamma_l^{(a)}`` | ``a=1`` | ``a=2`` | ``a=3`` |
# |:--|:--|:--|:--|
# | ``l=2`` | 439.2 | 4632.6 | 51233 |
# | ``l=3`` | 659.4 | 9043.6 | 134397 |
# | ``l=4`` | 906.1 | 15349 | 290166 |
# | ``l=5`` | 1187.6 | 24031 | 557087 |
#
# The general-``l`` routine and the exact ``l=2`` symbolic result agree to a
# relative ``3.7\times10^{-8}`` at ``a=2`` -- limited by the finite-difference
# ``\theta``-derivatives, not by the construction -- and ``\Gamma_l^{(a)}>0``
# throughout, which it must be for a quantity that is an integral of a
# positive quantity.

function Gamma_l_a_general_l(l_val_, a_val)                                   #src
    fp_c(xv) = -(l_val_ + 1) * xv^l_val_                                      #src
    f_c(xv) = -xv^(l_val_ + 1)                                                #src
    Plfun(u) = begin   # Legendre polynomial via Bonnet's recursion (u = cos(theta)) #src
        P0, P1 = one(u), u                                                    #src
        l_val_ == 0 && return P0                                              #src
        l_val_ == 1 && return P1                                              #src
        Pkm1, Pk = P0, P1                                                     #src
        for k in 1:(l_val_-1)                                                 #src
            Pkp1 = ((2k + 1) * u * Pk - k * Pkm1) / (k + 1)                    #src
            Pkm1, Pk = Pk, Pkp1                                               #src
        end                                                                   #src
        Pk                                                                    #src
    end                                                                       #src
    function H_at_theta(thv)                                                  #src
        xv = 1.0                                                              #src
        u = cos(thv)                                                          #src
        s = sin(thv)                                                          #src
        h = 1e-6   # numeric differentiation of the Legendre function in theta  #src
        Pl_here = Plfun(cos(thv))                                             #src
        dPl_here = (Plfun(cos(thv + h)) - Plfun(cos(thv - h))) / (2h)          #src
        d2Pl_here = (Plfun(cos(thv + h)) - 2 * Plfun(cos(thv)) + Plfun(cos(thv - h))) / h^2 #src
        g_c = (2 * f_c(xv) + xv * fp_c(xv)) / (l_val_ * (l_val_ + 1))          #src
        fpp_c = -(l_val_ + 1) * l_val_ * xv^(l_val_ - 1)                       #src
        gp_c = (3 * fp_c(xv) + xv * fpp_c) / (l_val_ * (l_val_ + 1))           #src
        e_rr = fp_c(xv) * Pl_here                                             #src
        e_rth = 0.5 * (gp_c - g_c / xv + f_c(xv) / xv) * dPl_here             #src
        e_thth = (f_c(xv) / xv) * Pl_here + (g_c / xv) * d2Pl_here            #src
        e_phph = (f_c(xv) / xv) * Pl_here + (g_c * u / (s * xv)) * dPl_here    #src
        2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2)                       #src
    end                                                                       #src
    exponent = (a_val + 2) / 2                                                #src
    ang, _ = quadgk(th -> H_at_theta(th)^exponent * sin(th), 1e-6, pi - 1e-6)  #src
    N_l_val = 1 / (2 * l_val_ + 5)                                            #src
    (ang / (l_val_ * (a_val + 2) + 3)) / N_l_val^exponent                     #src
end                                                                           #src

## Numeric theta-derivatives (finite difference) limit precision here vs the  #src
## exact symbolic l=2 result, hence 1e-4 rather than machine precision.       #src
@assert abs(Gamma_l_a_general_l(2, 2.0) - Gamma2_exact_a2) / Gamma2_exact_a2 < 1e-4 #src
println("ASSERTION 16 OK: the general-l routine matches the exact l=2, a=2 value") #src
for l_val_ in (2, 3, 4, 5), a_val in (0.5, 1.0, 1.5, 2.0, 3.0)                #src
    @assert Gamma_l_a_general_l(l_val_, a_val) > 0                            #src
end                                                                           #src
@assert isapprox(Gamma_l_a_general_l(5, 3.0), 557086.63; rtol=1e-4)   # pins the table #src
println("ASSERTION 17 OK: Gamma_l^(a) > 0 for l=2..5 and a in {0.5,1,1.5,2,3}") #src

# ### Finite ``\mathrm{Oh}``
#
# At finite ``\mathrm{Oh}`` the eigenvalue ``q`` acquires a negative
# imaginary part and the eigenfunction becomes complex,
# ``U=U_R+iU_I``. The old ``a=2``-hardcoded notebook computed
# ``\langle\dot\gamma^4\rangle`` through a closed-form algebraic identity
# specific to the exponent 4. That identity does not generalize, so here
# ``\langle[\cdot]^{(a+2)/2}\rangle`` is computed by direct numerical
# time-averaging over one period -- fully general, and reducible to the old
# route only at ``a=2``, which is the one case where both are available to
# compare.
#
# One normalization subtlety carries over. The old ``\Gamma_l`` is a purely
# spatial integral with no time average baked in (§7-§8 apply the temporal
# average separately), and the old code divides its raw time-averaged
# quantity by ``\langle\sin^4\rangle=3/8`` "to match the inviscid
# convention", in its own words. Generalized, that divisor is
#
# ```math
# \langle\sin^{a+2}\rangle
#  = \frac{\Gamma\!\left(\frac{a+3}{2}\right)}
#         {\sqrt\pi\,\Gamma\!\left(\frac{a+4}{2}\right)},
# ```
#
# which equals ``3/8`` exactly at ``a=2``.
#
# The two routes agree at ``a=2`` to relative error ``\le3\times10^{-16}`` --
# machine precision -- at every ``\mathrm{Oh}`` tested. The sweep:
#
# | ``\mathrm{Oh}`` | 0.001 | 0.01 | 0.05 | 0.1 | 0.2 | 0.3 | 0.5 |
# |:--|:--|:--|:--|:--|:--|:--|:--|
# | ``\Gamma_2^{(1)}`` | 417.1 | 367.5 | 269.7 | 201.4 | 147.4 | 131.9 | 123.0 |
# | ``\Gamma_2^{(2)}`` | 4283.7 | 3540.9 | 2270.4 | 1505.6 | 972.7 | 833.0 | 756.0 |
# | ``\Gamma_2^{(3)}`` | 46002 | 35477 | 19796 | 11619 | 6596 | 5400 | 4766 |
#
# ``\Gamma_l^{(a)}(\mathrm{Oh})`` is positive and strictly decreasing in
# ``\mathrm{Oh}`` for every exponent -- viscosity smooths the mode shape and
# reduces the shear it generates, so an increasing column would be
# physically wrong, not merely surprising.
#
# Two sanity anchors sit at the ends of the sweep. At ``\mathrm{Oh}=0.001``
# the viscous correction relative to the inviscid ``\Gamma_2`` is ``7.5\%``,
# inside the ``4\%``-``12\%`` band the old notebook established -- so even
# a nearly inviscid drop carries a several-percent correction. And the decay
# rate at ``\mathrm{Oh}=0.3`` is ``\gamma=1.097`` against Lamb's
# ``(l-1)(2l+1)\mathrm{Oh}=1.5``, i.e. 27% below: at that Ohnesorge number
# Lamb's thin-boundary-layer formula is already a poor approximation, which
# is exactly why the finite-``\mathrm{Oh}`` treatment exists.

function sph_bessel_j(l, z)                                                   #src
    sqrt(pi / (2z)) * besselj(l + 0.5, z)                                     #src
end                                                                           #src
bessel_ratio(l, z) = sph_bessel_j(l + 1, z) / sph_bessel_j(l, z)              #src

function reid_char(qv, Oh, l)                                                 #src
    alpha2 = sqrt(l * (l - 1) * (l + 2)) / Oh                                 #src
    Qv = bessel_ratio(l, qv)                                                  #src
    lhs = alpha2^2 / qv^4 + 1                                                 #src
    rhs = 2 * (l - 1) / qv^2 * (l + (l + 1) * (qv - 2 * l * Qv) / (qv - 2 * Qv)) #src
    lhs - rhs                                                                 #src
end                                                                           #src

function find_eigenvalue(Oh_val, l=2; maxiter=200, tol=1e-13)                 #src
    sigma0 = sqrt(l * (l - 1) * (l + 2))                                      #src
    gamma0 = (l - 1) * (2l + 1) * Oh_val                                      #src
    q = sqrt(complex(gamma0 / Oh_val, -sigma0 / Oh_val))                       #src
    imag(q) > 0 && (q = -q)                                                   #src
    for _ in 1:maxiter                                                        #src
        F = reid_char(q, Oh_val, l)                                           #src
        h = max(1e-7, abs(q) * 1e-7)                                          #src
        dF = (reid_char(q + h * im, Oh_val, l) - F) / (h * im)                #src
        q -= F / dF                                                           #src
        abs(F) < tol && break                                                 #src
    end                                                                       #src
    imag(q) > 0 ? conj(q) : q                                                 #src
end                                                                           #src

function eigenfunction(x_arr, qv, l=2)                                        #src
    jlq = sph_bessel_j(l, qv)                                                 #src
    jl1q = sph_bessel_j(l + 1, qv)                                            #src
    Cv = 2 * (l^2 - 1) / (qv * (2 * jl1q - qv * jlq))                          #src
    Pi0v = -1 - Cv * jlq                                                      #src
    n = length(x_arr)                                                         #src
    U = zeros(ComplexF64, n)                                                  #src
    dU = zeros(ComplexF64, n)                                                 #src
    d2U = zeros(ComplexF64, n)                                                #src
    for (i, xv) in enumerate(x_arr)                                           #src
        xv < 1e-12 && continue                                                #src
        jlv = sph_bessel_j(l, qv * xv)                                        #src
        jl1v = sph_bessel_j(l + 1, qv * xv)                                   #src
        jl2v = sph_bessel_j(l + 2, qv * xv)                                   #src
        U[i] = Cv * xv * jlv + Pi0v * xv^(l + 1)                               #src
        dU[i] = Cv * ((l + 1) * jlv - qv * xv * jl1v) + Pi0v * (l + 1) * xv^l  #src
        d2U[i] = Cv * (l * (l + 1) / xv * jlv - (2l + 3) * qv * jl1v + qv^2 * xv * jl2v) + Pi0v * l * (l + 1) * xv^(l - 1) #src
    end                                                                       #src
    U, dU, d2U                                                                #src
end                                                                           #src

function legendre_arrays(l, u_pts)                                            #src
    n = length(u_pts)                                                         #src
    Pl = zeros(n)                                                             #src
    Pl1 = zeros(n)                                                            #src
    for (i, u) in enumerate(u_pts)                                            #src
        P0, P1 = 1.0, u                                                       #src
        if l == 0                                                             #src
            Pl[i] = P0                                                        #src
        elseif l == 1                                                         #src
            Pl[i] = P1                                                        #src
        else                                                                  #src
            pkm1, pk = P0, P1                                                 #src
            for k in 1:(l-1)                                                  #src
                pkp1 = ((2k + 1) * u * pk - k * pkm1) / (k + 1)                #src
                pkm1, pk = pk, pkp1                                           #src
            end                                                               #src
            Pl[i] = pk                                                        #src
        end                                                                   #src
        if l - 1 == 0   # P_{l-1}                                             #src
            Pl1[i] = 1.0                                                      #src
        else                                                                  #src
            pkm1, pk = 1.0, u                                                 #src
            for k in 1:(l-2)                                                  #src
                pkp1 = ((2k + 1) * u * pk - k * pkm1) / (k + 1)                #src
                pkm1, pk = pk, pkp1                                           #src
            end                                                               #src
            Pl1[i] = pk                                                       #src
        end                                                                   #src
    end                                                                       #src
    dPl = similar(Pl)                                                         #src
    d2Pl = similar(Pl)                                                        #src
    sin_th = similar(Pl)                                                      #src
    for i in 1:n                                                              #src
        u = u_pts[i]                                                          #src
        sin_th[i] = sqrt(max(1 - u^2, 0.0))                                   #src
        dPl[i] = abs(u) < 1 - 1e-14 ? l * (Pl1[i] - u * Pl[i]) / (1 - u^2) : 0.0 #src
    end                                                                       #src
    for i in 1:n                                                              #src
        u = u_pts[i]                                                          #src
        d2Pl[i] = abs(u) < 1 - 1e-14 ? (2 * u * dPl[i] - l * (l + 1) * Pl[i]) / (1 - u^2) : 0.0 #src
    end                                                                       #src
    Pl, dPl, d2Pl, sin_th                                                     #src
end                                                                           #src

function gdot_sq_matrix(f_arr, df_arr, x_arr, l, u_pts, d2f_arr)              #src
    nx, nth = length(x_arr), length(u_pts)                                    #src
    g_arr = (2 .* f_arr .+ x_arr .* df_arr) ./ (l * (l + 1))                   #src
    dg_arr = (3 .* df_arr .+ x_arr .* d2f_arr) ./ (l * (l + 1))                #src
    Pl, dPl, d2Pl, sin_th = legendre_arrays(l, u_pts)                         #src
    H = zeros(nx, nth)                                                        #src
    for i in 1:nx                                                             #src
        xv = x_arr[i]                                                         #src
        xv < 1e-12 && continue                                                #src
        f_, fp, g_, gp = f_arr[i], df_arr[i], g_arr[i], dg_arr[i]              #src
        for j in 1:nth                                                        #src
            e_rr = fp * Pl[j]                                                 #src
            e_thth = (f_ / xv) * Pl[j] - (g_ / xv) * (u_pts[j] * dPl[j] - sin_th[j]^2 * d2Pl[j]) #src
            e_phph = (f_ / xv) * Pl[j] - (g_ / xv) * (u_pts[j] * dPl[j])       #src
            e_rth = 0.5 * sin_th[j] * dPl[j] * (g_ / xv - gp - f_ / xv)        #src
            H[i, j] = 2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2)         #src
        end                                                                   #src
    end                                                                       #src
    H                                                                         #src
end                                                                           #src

function gdot_cross_matrix(f1, df1, d2f1, f2, df2, d2f2, x_arr, l, u_pts)     #src
    nx, nth = length(x_arr), length(u_pts)                                    #src
    g1 = (2 .* f1 .+ x_arr .* df1) ./ (l * (l + 1))                            #src
    dg1 = (3 .* df1 .+ x_arr .* d2f1) ./ (l * (l + 1))                         #src
    g2 = (2 .* f2 .+ x_arr .* df2) ./ (l * (l + 1))                            #src
    dg2 = (3 .* df2 .+ x_arr .* d2f2) ./ (l * (l + 1))                         #src
    Pl, dPl, d2Pl, sin_th = legendre_arrays(l, u_pts)                         #src
    Hc = zeros(nx, nth)                                                       #src
    for i in 1:nx                                                             #src
        xv = x_arr[i]                                                         #src
        xv < 1e-12 && continue                                                #src
        f1_, fp1, g1_, gp1 = f1[i], df1[i], g1[i], dg1[i]                      #src
        f2_, fp2, g2_, gp2 = f2[i], df2[i], g2[i], dg2[i]                      #src
        for j in 1:nth                                                        #src
            e_rr1 = fp1 * Pl[j]                                               #src
            e_thth1 = (f1_ / xv) * Pl[j] - (g1_ / xv) * (u_pts[j] * dPl[j] - sin_th[j]^2 * d2Pl[j]) #src
            e_phph1 = (f1_ / xv) * Pl[j] - (g1_ / xv) * (u_pts[j] * dPl[j])    #src
            e_rth1 = 0.5 * sin_th[j] * dPl[j] * (g1_ / xv - gp1 - f1_ / xv)    #src
            e_rr2 = fp2 * Pl[j]                                               #src
            e_thth2 = (f2_ / xv) * Pl[j] - (g2_ / xv) * (u_pts[j] * dPl[j] - sin_th[j]^2 * d2Pl[j]) #src
            e_phph2 = (f2_ / xv) * Pl[j] - (g2_ / xv) * (u_pts[j] * dPl[j])    #src
            e_rth2 = 0.5 * sin_th[j] * dPl[j] * (g2_ / xv - gp2 - f2_ / xv)    #src
            Hc[i, j] = 2 * (e_rr1 * e_rr2 + e_thth1 * e_thth2 + e_phph1 * e_phph2 + 2 * e_rth1 * e_rth2) #src
        end                                                                   #src
    end                                                                       #src
    Hc                                                                        #src
end                                                                           #src

## Gauss-Legendre nodes/weights via Jacobi-matrix eigendecomposition -- no    #src
## external quadrature-rule package needed, matching julia/src/types.jl's own #src
## make_theta_vec approach.                                                   #src
function gauss_legendre(n)                                                    #src
    bvec = [k / sqrt(4 * k^2 - 1) for k in 1:n-1]                              #src
    J = diagm(1 => bvec, -1 => bvec)                                          #src
    F = eigen(J)                                                              #src
    F.values, 2 .* (F.vectors[1, :]) .^ 2                                     #src
end                                                                           #src

norm_p(p) = gamma(p + 0.5) / (sqrt(pi) * gamma(p + 1))                        #src

function trapz(xs, ys)                                                        #src
    s = 0.0                                                                   #src
    for i in 1:length(xs)-1                                                   #src
        s += (xs[i+1] - xs[i]) * (ys[i+1] + ys[i]) / 2                         #src
    end                                                                       #src
    s                                                                         #src
end                                                                           #src

function compute_gamma_l_a(Oh_val, a_val, l=2; n_x=300, n_theta=100, n_phi=128) #src
    qv = find_eigenvalue(Oh_val, l)                                           #src
    sigma = qv^2 * Oh_val                                                     #src
    gamma_, omega_ = real(sigma), -imag(sigma)                                #src
    x_arr = range(1e-4, 1.0, length=n_x)                                      #src
    U, dU, d2U = eigenfunction(x_arr, qv, l)                                  #src
    u_pts, wts = gauss_legendre(n_theta)                                      #src
    fR, dfR, d2fR = real.(U), real.(dU), real.(d2U)                            #src
    fI, dfI, d2fI = imag.(U), imag.(dU), imag.(d2U)                            #src
    H_R = gdot_sq_matrix(fR, dfR, x_arr, l, u_pts, d2fR)                      #src
    H_I = gdot_sq_matrix(fI, dfI, x_arr, l, u_pts, d2fI)                      #src
    H_cross = gdot_cross_matrix(fR, dfR, d2fR, fI, dfI, d2fI, x_arr, l, u_pts) #src
    p = (a_val + 2) / 2                                                       #src
    time_avg = zeros(n_x, n_theta)                                            #src
    for k in 1:n_phi                                                          #src
        phv = 2pi * (k - 1) / n_phi                                           #src
        sinp, cosp = sin(phv), cos(phv)                                       #src
        @. time_avg += abs(H_R * sinp^2 + H_I * cosp^2 + 2 * H_cross * sinp * cosp)^p #src
    end                                                                       #src
    time_avg ./= n_phi                                                        #src
    time_avg ./= norm_p(p)                                                    #src
    numerator = trapz(x_arr, (time_avg * wts) .* x_arr .^ 2)                   #src
    N_l = trapz(x_arr, (abs.(U) .^ 2) .* x_arr .^ 2)                           #src
    numerator / N_l^p, gamma_, omega_                                         #src
end                                                                           #src

## The old (a=2-only) closed-form route, kept solely as the independent       #src
## comparison for the generalized time-averaging above.                       #src
function compute_gamma_l_original(Oh_val, l=2; n_x=300, n_theta=100)          #src
    qv = find_eigenvalue(Oh_val, l)                                           #src
    sigma = qv^2 * Oh_val                                                     #src
    gamma_, omega_ = real(sigma), -imag(sigma)                                #src
    x_arr = range(1e-4, 1.0, length=n_x)                                      #src
    U, dU, d2U = eigenfunction(x_arr, qv, l)                                  #src
    u_pts, wts = gauss_legendre(n_theta)                                      #src
    H_R = gdot_sq_matrix(real.(U), real.(dU), x_arr, l, u_pts, real.(d2U))    #src
    H_I = gdot_sq_matrix(imag.(U), imag.(dU), x_arr, l, u_pts, imag.(d2U))    #src
    H_cross = gdot_cross_matrix(real.(U), real.(dU), real.(d2U), imag.(U), imag.(dU), imag.(d2U), x_arr, l, u_pts) #src
    H_sq = @. H_R^2 + H_I^2 + (2 / 3) * H_R * H_I + (4 / 3) * H_cross^2        #src
    numerator = trapz(x_arr, (H_sq * wts) .* x_arr .^ 2)                       #src
    N_l = trapz(x_arr, (abs.(U) .^ 2) .* x_arr .^ 2)                           #src
    numerator / N_l^2, gamma_, omega_                                         #src
end                                                                           #src

let mismatches = Float64[]                                                    #src
    for Oh_v in (0.001, 0.05, 0.1, 0.3)                                       #src
        G_gen, _, _ = compute_gamma_l_a(Oh_v, 2.0, 2)                          #src
        G_orig, _, _ = compute_gamma_l_original(Oh_v, 2)                       #src
        abs(G_gen - G_orig) / G_orig > 1e-6 && push!(mismatches, Oh_v)         #src
    end                                                                       #src
    @assert isempty(mismatches)                                               #src
end                                                                           #src
println("ASSERTION 18 OK: the generalized finite-Oh Gamma matches the old closed form at a=2") #src
@assert abs(norm_p(2.0) - 0.375) < 1e-12                                      #src
println("ASSERTION 19 OK: <sin^(a+2)> = 3/8 exactly at a=2, matching the old convention") #src

const _OH_SWEEP = (0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)                     #src
sweep_a2 = Tuple{Float64,Float64,Float64,Float64}[]                           #src
for Oh_v in _OH_SWEEP                                                         #src
    G, gam, om = compute_gamma_l_a(Oh_v, 2.0, 2)                               #src
    push!(sweep_a2, (Oh_v, G, gam, om))                                       #src
end                                                                           #src
for a_val in (1.0, 2.0, 3.0)                                                  #src
    vals = [compute_gamma_l_a(Oh_v, a_val, 2)[1] for Oh_v in _OH_SWEEP]        #src
    @assert all(v -> v > 0, vals)                                             #src
    @assert all(i -> vals[i] > vals[i+1], 1:length(vals)-1)                    #src
end                                                                           #src
@assert isapprox(sweep_a2[1][2], 4283.69; rtol=1e-4)   # pins the tabulated sweep #src
@assert isapprox(sweep_a2[end][2], 756.03; rtol=1e-4)                         #src
println("ASSERTION 20 OK: Gamma_l^(a)(Oh) positive and strictly decreasing in Oh, for a in {1,2,3}") #src

let rel_corr = (Gamma2_exact_a2 - sweep_a2[1][2]) / Gamma2_exact_a2           #src
    @assert 0.04 < rel_corr < 0.12                                            #src
end                                                                           #src
println("ASSERTION 21 OK: the Oh=0.001 viscous correction (a=2) sits in the 4%-12% band") #src

let gamma_03 = first(r[3] for r in sweep_a2 if abs(r[1] - 0.3) < 1e-9)         #src
    gamma_Lamb = (2 - 1) * (2 * 2 + 1) * 0.3                                  #src
    @assert abs(gamma_03 - gamma_Lamb) / gamma_Lamb < 0.30                     #src
    @assert isapprox(gamma_03, 1.0971; rtol=1e-3)   # pins the reported value   #src
end                                                                           #src
println("ASSERTION 22 OK: gamma(Oh=0.3, a=2) = 1.097, 27% below Lamb's 1.5")   #src

# ## 7. The modified mode equation and its effective damping
#
# The shear-thinning correction enters the mode dynamics through a
# generalized dissipation function ``\delta\Phi=-(\text{coeff})|\dot
# b|^{a+2}``, giving a dissipative force
#
# ```math
# \frac{\partial(\delta\Phi)}{\partial\dot b}
#   \;=\; -(a+2)\,|\dot b|^{a}\,\dot b .
# ```
#
# That differentiation step is checked numerically at
# ``a\in\{0.5,1,1.5,2,3\}`` rather than trusted, because it is the exact step
# an earlier attempt at the same generalization for the Cross model got wrong
# -- writing ``|\dot b|^{a-1}`` and thereby predicting a divergence that does
# not exist.
#
# The effective damping coefficient is therefore
#
# ```math
# D_{\mathrm{eff}} \;=\; D_l^{(0)}\left[1
#   - \varepsilon_{ST}\,\Lambda^{a}\,\Gamma_l^{(a)}\,|\dot b|^{a}\right],
# \qquad
# D_l^{(0)} = 2\,\mathrm{Oh}\,(l-1)(2l+1),
# ```
#
# which at ``l=2`` is ``10\,\mathrm{Oh}``. At ``a=2`` this collapses exactly
# onto the old notebook's ``D_{\mathrm{eff}}`` -- comparing on equal footing
# by substituting ``|\dot b|\to\sqrt{\dot b^2}``, since the two are identical
# for real ``\dot b`` -- and at ``\varepsilon_{ST}=0`` it returns the
# Newtonian damping for any ``a``.

function check_force_order(a_val; x0=0.37, h=1e-6)                            #src
    f(xv) = -abs(xv)^(a_val + 2)                                              #src
    dfdx = (f(x0 + h) - f(x0 - h)) / (2h)                                     #src
    claimed = -(a_val + 2) * abs(x0)^a_val * x0                               #src
    dfdx, claimed                                                             #src
end                                                                           #src
let mismatches = Float64[]                                                    #src
    for a_val in (0.5, 1.0, 1.5, 2.0, 3.0)                                    #src
        d, c = check_force_order(a_val)                                       #src
        abs(d - c) < 1e-4 || push!(mismatches, a_val)                         #src
    end                                                                       #src
    @assert isempty(mismatches)                                               #src
end                                                                           #src
println("ASSERTION 23 OK: the generalized force ~ |bdot|^a * bdot, verified numerically") #src

@variables Oh_sym Lambda_sym eps_ST_sym Gamma_l_sym bdot_sym                  #src
D0_newtonian = 2 * Oh_sym * (l_val - 1) * (2 * l_val + 1)                     #src
D_eff = D0_newtonian * (1 - eps_ST_sym * Lambda_sym^a_shape * Gamma_l_sym * abs(bdot_sym)^a_shape) #src
D_eff_a2 = simplify(substitute(D_eff, Dict(a_shape => 2, abs(bdot_sym) => sqrt(bdot_sym^2)))) #src
D_eff_original = D0_newtonian * (1 - eps_ST_sym * Lambda_sym^2 * Gamma_l_sym * bdot_sym^2) #src
@assert numerically_equal(D_eff_a2, D_eff_original)                           #src
println("ASSERTION 24 OK: the generalized D_eff reduces exactly to the old D_eff at a=2") #src
@assert numerically_equal(substitute(D_eff, Dict(eps_ST_sym => 0)), D0_newtonian) #src
println("ASSERTION 25 OK: eps_ST=0 recovers the Newtonian damping for any a")  #src

# ## 8. Secular terms and the slow-amplitude equation
#
# The old ``a=2``-hardcoded derivation's cubic damping term projects onto
# ``\sin\omega t`` through ``\langle\sin^4\rangle=3/8``, producing the
# familiar ``3/4`` secular factor. For the generalized
# ``|\dot b|^{a}\dot b`` damping -- non-polynomial for any non-even ``a`` --
# the analogous projection is a classical Wallis-type integral, closed-form
# through the Gamma function for any real ``a``:

C_of_a(a_val) = (2 / sqrt(pi)) * gamma((a_val + 3) / 2) / gamma((a_val + 4) / 2)

# Against direct quadrature of ``\frac{2}{\pi}\int_0^\pi\sin^{a+2}\theta\,
# d\theta``, this is exact to better than ``10^{-11}`` at every exponent
# tested:
#
# | ``a`` | 0.5 | 1 | 1.5 | 2 | 3 | 4 |
# |:--|:--|:--|:--|:--|:--|:--|
# | ``C(a)`` | 0.9153 | 0.8488 | 0.7949 | 0.7500 | 0.6791 | 0.6250 |
#
# ``C(2)=3/4`` exactly, which is the entry that has to come out right, and
# the direct ``\langle\sin^4\rangle`` period integral independently returns
# ``3/4`` as well.
#
# Assembling everything gives the generalized slow-amplitude equation:
#
# ```math
# \boxed{\;
# \frac{da}{dt} = -\gamma_l^{(0)}\,a
#   \left[1 - C(a_{\text{shape}})\,\varepsilon_{ST}\,\Lambda^{a_{\text{shape}}}\,
#   \Gamma_l^{(a_{\text{shape}})}\,a^{a_{\text{shape}}}\right]
# \;}
# ```
#
# At ``a=2`` it collapses exactly onto the old notebook's boxed equation
# ``da/dt=-\gamma_0a(1-\tfrac34\varepsilon_{ST}\Lambda^2\Gamma_la^2)``,
# verified at several concrete parameter points rather than symbolically --
# an equally rigorous check that avoids fragile symbolic simplification of
# `gamma()`. At ``\varepsilon_{ST}=0`` it returns ``da/dt=-\gamma_0a`` for
# any ``a``. And the correction vanishes as the amplitude ``a\to0`` for every
# ``a_{\text{shape}}>0``, so the equation is well posed at small amplitude --
# which is precisely the property the mistaken ``a^{m-1}`` order counting
# would have destroyed.

let mismatches = Float64[]                                                    #src
    for a_val in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)                               #src
        closed = C_of_a(a_val)                                                #src
        numeric, _ = quadgk(th -> sin(th)^(a_val + 2), 0, pi)                  #src
        numeric *= 2 / pi                                                     #src
        abs(closed - numeric) > 1e-9 && push!(mismatches, a_val)               #src
    end                                                                       #src
    @assert isempty(mismatches)                                               #src
end                                                                           #src
println("ASSERTION 26 OK: the closed-form C(a) matches direct quadrature for every a tested") #src
@assert abs(C_of_a(2.0) - 0.75) < 1e-12                                       #src
println("ASSERTION 27 OK: C(2) = 3/4, matching the old validated secular factor") #src
let coeff_34 = first(quadgk(t -> sin(t)^4, 0, 2pi)) / pi                      #src
    @assert abs(coeff_34 - 0.75) < 1e-10                                      #src
end                                                                           #src
println("ASSERTION 28 OK: the direct <sin^4> period integral still gives exactly 3/4") #src

function da_dt_general(gamma0, a_amp, eps_ST, Lambda, Gamma_l, a_val)         #src
    -gamma0 * a_amp * (1 - C_of_a(a_val) * eps_ST * Lambda^a_val * Gamma_l * a_amp^a_val) #src
end                                                                           #src
function da_dt_original(gamma0, a_amp, eps_ST, Lambda, Gamma_l)               #src
    -gamma0 * a_amp * (1 - 0.75 * eps_ST * Lambda^2 * Gamma_l * a_amp^2)       #src
end                                                                           #src

let mismatches = []                                                           #src
    for (gamma0, a_amp, eps_ST, Lambda, Gamma_l) in                           #src
        [(1.0, 0.3, 0.1, 0.2, 5.0), (0.5, 0.1, 0.05, 0.15, 100.0), (2.0, 0.4, 0.02, 0.3, 4632.639)] #src
        gen = da_dt_general(gamma0, a_amp, eps_ST, Lambda, Gamma_l, 2.0)       #src
        orig = da_dt_original(gamma0, a_amp, eps_ST, Lambda, Gamma_l)         #src
        abs(gen - orig) > 1e-9 && push!(mismatches, (gamma0, a_amp, eps_ST, Lambda, Gamma_l)) #src
    end                                                                       #src
    @assert isempty(mismatches) "mismatches: $mismatches"                      #src
end                                                                           #src
println("ASSERTION 29 OK: the generalized slow-amplitude equation collapses exactly to the old one at a=2") #src

for (gamma0, a_amp, Lambda, Gamma_l) in [(1.0, 0.3, 0.2, 5.0), (2.0, 0.1, 0.3, 100.0)] #src
    @assert abs(da_dt_general(gamma0, a_amp, 0.0, Lambda, Gamma_l, 2.5) - (-gamma0 * a_amp)) < 1e-12 #src
end                                                                           #src
## Well-posedness: the bracketed correction must vanish as the amplitude does. #src
for a_shape_val in (0.5, 1.0, 2.0, 3.0)                                       #src
    corr = [C_of_a(a_shape_val) * 0.1 * 0.2^a_shape_val * 5.0 * amp^a_shape_val #src
            for amp in (1e-1, 1e-3, 1e-6, 1e-12)]                             #src
    @assert all(i -> corr[i] > corr[i+1], 1:length(corr)-1)                    #src
    @assert corr[end] < 1e-6                                                  #src
end                                                                           #src
println("ASSERTION 30 OK: eps_ST=0 gives da/dt = -gamma0*a, and the correction vanishes as a -> 0 for every a") #src

# ## What this produced, and why it was replaced
#
# This is a strict generalization of the Carreau-only derivation: at every
# point where the two can be compared at ``a=2``, the match is exact or at
# machine precision. What it put into the code was small --
# `STParams` gained an `a::Float64` field defaulting to 2.0, so existing
# callers were unaffected, and the two places `julia/src/st_extension.jl`
# hardcoded the exponent 2 (`Gamma_eff = Gamma .* (lambda_c.*sigma0).^2` and
# `shear_sq_lag = sum(Gamma_eff.*Adot_prev.^2)`) became `.^a`, with an
# `abs()` added on the second since ``a`` need not be an even integer.
#
# The conversion recipe for real experimental data: given Cross-model fit
# parameters ``(K,m,\mu_0,\mu_\infty)``, take ``\lambda_c=K``, ``a=m``,
# ``n=1-m``, ``\varepsilon_{ST}=[(\mu_0-\mu_\infty)/\mu_0](1-n)/a``, and
# compute ``\Gamma_l^{(a)}(\mathrm{Oh})`` with §6's finite-``\mathrm{Oh}``
# machinery.
#
# **And then the real fluid arrived.** Its fitted parameters --
# ``\varepsilon_{ST}\approx0.9996``, dimensionless
# ``\lambda_c\approx3\times10^{4}`` -- sit far outside the double-small-
# parameter regime this entire expansion assumes. The bracketed factor
# ``1-\varepsilon_{ST}\Lambda^a\Gamma_l^{(a)}|\dot b|^a`` has no lower bound
# once ``\varepsilon_{ST}\approx1``, so at realistic shear rates it turns
# negative and the damping term starts injecting energy. That is not a
# tolerance problem, and it is why
# `carreau_yasuda_nonperturbative_derivation.jl` abandons the expansion for
# the exact law, and `carreau_yasuda_multimode_derivation.jl` goes further
# and abandons the per-mode shear rate as well.
#
# What survives and is still used: ``\Gamma_l^{(a)}``, ``C(a)``, the
# ``H(\theta)`` factorization, and the finite-``\mathrm{Oh}`` eigenfunction
# machinery of §6.
