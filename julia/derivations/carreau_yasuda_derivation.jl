# # Weakly Nonlinear Shear-Thinning Drop: Carreau-Yasuda Derivation
#
# Derives the perturbation expansion for a shear-thinning drop oscillating on
# a flat substrate, using the Carreau-Yasuda constitutive law -- Carreau with
# a free shape exponent ``a``, rather than ``a`` hardcoded to 2. Follows the
# same non-dimensionalisation as the Newtonian (Reid 1960) and Oldroyd-B
# analyses in the companion LaTeX documents (`docs/reid1960_expanded-3.tex`,
# `docs/section_oldroydB.tex`).
#
# **Why Carreau-Yasuda, not plain Carreau?** The standard Carreau model's
# nonlinearity is always built from ``(\lambda_c\dot\gamma)^2`` -- the shape
# of the transition between the zero-shear plateau and the power-law region
# is fixed. Carreau-Yasuda adds one parameter, ``a``, controlling that
# transition shape; ``a=2`` recovers standard Carreau exactly. Real
# experimental shear-thinning fits (including Cross-model characterizations,
# which convert directly to an equivalent Carreau-Yasuda ``a``) commonly
# need ``a \neq 2`` to fit well.
#
# Each section derives one piece of the chain symbolically or numerically,
# then asserts the result against an independent check. **A failing
# assertion means the derivation broke** -- fix the math, not the assertion.
#
# ## Notation
# (matches `reid1960_expanded-3.tex`)
#
# | symbol | meaning |
# |:--|:--|
# | ``x = r/R`` | dimensionless radial coordinate |
# | ``\varepsilon`` | oscillation amplitude, small parameter |
# | ``\sigma = -\gamma + i\omega`` | complex modal frequency (``\gamma>0`` = decay) |
# | ``q``, ``\sigma = q^2 \mathrm{Oh}`` | viscous wavenumber |
# | ``\alpha^2 = \sigma_{l;0}/\mathrm{Oh}``, ``\sigma_{l;0} = \sqrt{l(l-1)(l+2)}`` | |
# | ``\mathrm{Oh} = \mu/\sqrt{\rho T_1 R}`` | Ohnesorge number |
# | ``U(x)`` | radial velocity eigenfunction (Reid) |
# | ``b_l(t)`` | dimensionless amplitude of mode ``l`` (solver's ``A_l``) |
# | ``P_l(\cos\theta)`` | Legendre polynomial of degree ``l`` |
# | ``a`` | the Carreau-Yasuda SHAPE EXPONENT (a model parameter) -- distinct from the mode AMPLITUDE ``a(t)`` used in Sections 7-8; the shape exponent is always `a_shape` in code |

using Symbolics
using QuadGK
using SpecialFunctions
using DropSolver

nth_derivative_at_zero(expr, var, n) = begin
    d = expr
    for _ in 1:n
        d = Symbolics.derivative(d, var)
    end
    simplify(substitute(d, Dict(var => 0)))
end

"""
Robust symbolic-equality check: Symbolics.jl's `simplify` does not always
collapse an algebraically-zero difference to the literal `0` (a known
limitation, weaker than sympy's simplifier here). Rather than trust
`isequal(simplify(a-b), 0)` alone, evaluate both expressions at several
concrete numeric points for every free variable and require they match to
near machine precision -- an independent, purely numeric check that is
immune to this limitation, used throughout this script wherever a symbolic
identity needs verifying.
"""
const _FIXED_TEST_VALUES = (0.31, 0.57, 1.13, 1.94, 2.71)

function numerically_equal(expr1, expr2, test_points::Dict=Dict())
    vars = collect(union(Symbolics.get_variables(expr1), Symbolics.get_variables(expr2)))
    isempty(vars) && return isequal(simplify(expr1 - expr2), 0)
    f1 = Symbolics.build_function(expr1, vars...; expression=false)
    f2 = Symbolics.build_function(expr2, vars...; expression=false)
    for trial in 1:length(_FIXED_TEST_VALUES)
        vals = [haskey(test_points, v) ? test_points[v] :
                _FIXED_TEST_VALUES[mod1(trial + i, length(_FIXED_TEST_VALUES))]
                for (i, v) in enumerate(vars)]
        v1 = f1(vals...)
        v2 = f2(vals...)
        isapprox(v1, v2; atol=1e-8, rtol=1e-6) || return false
    end
    true
end

# ## Section 1: Carreau-Yasuda Constitutive Law and Perturbation Expansion
#
# The Carreau-Yasuda effective viscosity is
# ```math
# \mu_{\mathrm{eff}}(\dot\gamma) = \mu_0 \left[1 + (\lambda_c\dot\gamma)^{a}\right]^{(n-1)/a}
# ```
# ``n \in (0,1]``: power-law index (``n=1``: Newtonian). ``a>0``: shape
# exponent (``a=2``: standard Carreau). ``\mu_\infty=0`` (matches
# `julia/src/st_extension.jl`, which has never had a separate infinite-shear
# viscosity) -- a nonzero ``\mu_\infty`` from a real fluid characterization
# folds entirely into a redefined ``\varepsilon_{ST}`` below (see the note in
# Section 7), so this is not a loss of generality for this repo's linearized
# solver.

@variables mu_0 lambda_c eps gammadot_hat a_shape n_idx

x_small = (lambda_c * gammadot_hat)^a_shape * eps^a_shape
p_exp = (n_idx - 1) / a_shape
mu_CY_leading = mu_0 * (1 + p_exp * x_small)
println("Leading-order expansion: mu_CY ~ mu_0*(1 + p*x), p=(n-1)/a, x=(lambda_c*gammadot)^a")
println("mu_CY_leading = ", mu_CY_leading)

#md # ```@eval
#md # using Symbolics, Markdown
#md # @variables mu_0 lambda_c eps gammadot_hat a_shape n_idx
#md # x_small = (lambda_c * gammadot_hat)^a_shape * eps^a_shape
#md # p_exp = (n_idx - 1) / a_shape
#md # mu_CY_leading = mu_0 * (1 + p_exp * x_small)
#md # Markdown.parse("**Leading-order expansion** (as derived, not transcribed):\n```math\n" *
#md #     Main.pretty_latex(mu_CY_leading) * "\n```")
#md # ```

# ASSERTION 1: the leading-order expansion is the correct O(eps^a) Taylor term
# of (1+x)^p -- verified directly via repeated symbolic differentiation
# (Taylor's theorem) at CONCRETE integer a (Symbolics.jl has no built-in
# series() the way sympy does), confirming all intermediate powers vanish.
for a_val in (2, 3, 4)
    mu_CY_concrete = mu_0 * (1 + (lambda_c*eps*gammadot_hat)^a_val)^((n_idx-1)/a_val)
    coeff_exact = nth_derivative_at_zero(mu_CY_concrete, eps, a_val) / factorial(a_val)
    hand_built  = substitute(mu_CY_leading, Dict(a_shape => a_val))
    coeff_hand  = Symbolics.coeff(expand(hand_built - mu_0), eps^a_val)
    match = numerically_equal(coeff_exact, coeff_hand)
    println("a=$a_val: exact coeff=", simplify(coeff_exact), "  hand-built=", simplify(coeff_hand), "  match=$match")
    @assert match
    for lower_power in 1:(a_val-1)
        lower = nth_derivative_at_zero(mu_CY_concrete, eps, lower_power)
        @assert isequal(lower, 0) "a=$a_val: unexpected nonzero O(eps^$lower_power) term: $lower"
    end
end
println("ASSERTION 1 OK: hand-built leading-order expansion verified via repeated")
println("differentiation at a=2,3,4; all intermediate powers confirmed to vanish")

@assert numerically_equal(substitute(mu_CY_leading, Dict(n_idx => 1)), mu_0)
println("ASSERTION 2 OK: n=1 recovers Newtonian (zero correction) for any a")

@assert numerically_equal(substitute(mu_CY_leading, Dict(lambda_c => 0)), mu_0)
println("ASSERTION 3 OK: lambda_c=0 recovers Newtonian for any a")

# ASSERTION 4: exact reduction to the ORIGINAL (pre-generalization) Carreau
# result at a=2: eps_ST = (1-n)/2.
eps_ST_CY = (1 - n_idx) / a_shape
correction_a2 = simplify(substitute(mu_CY_leading - mu_0, Dict(a_shape => 2)))
target_original = eps.^2 .* gammadot_hat.^2 .* lambda_c.^2 .* mu_0 .* (n_idx - 1) ./ 2
@assert numerically_equal(correction_a2, target_original)
println("ASSERTION 4 OK: exact match at a=2 (eps_ST -> ", simplify(substitute(eps_ST_CY, Dict(a_shape=>2))),
        ", matching (1-n)/2)")

println()
println("="^78)
println("Sections 2-5: Newtonian base-flow machinery (rheology-independent)")
println("="^78)
println("""
The Reid velocity field, strain-rate tensor, and angular-integration
machinery below do not depend on the shear-thinning law at all -- they are
exact properties of the Newtonian base flow, and are identical regardless of
which shear-thinning correction (Carreau, Carreau-Yasuda, Cross, or none) is
applied at higher order.
""")

println("--- Section 2: Reid velocity field, l=2 ---")
println("""
U(x) = C*x*j_l(qx) + Pi0*x^(l+1) solves Reid's ODE:
  [d^2/dx^2 - l(l+1)/x^2 + q^2] U(x) = q^2 * Pi0 * x^(l+1)

Spherical Bessel functions of integer order reduce to elementary functions;
for l=2, j_2(z) = (3/z^3 - 1/z)*sin(z) - (3/z^2)*cos(z). We use this explicit
elementary form so Symbolics.jl's ordinary calculus (no special-function
differentiation rules needed) verifies the ODE directly.
""")

@variables x_sym q_sym C_sym Pi0_sym
j2(z) = (3/z^3 - 1/z)*sin(z) - (3/z^2)*cos(z)

l_test = 2
U_expr = C_sym * x_sym * j2(q_sym * x_sym) + Pi0_sym * x_sym^(l_test + 1)
U_expr_p  = Symbolics.derivative(U_expr, x_sym)
U_expr_pp = Symbolics.derivative(U_expr_p, x_sym)
ODE_lhs = U_expr_pp - l_test*(l_test+1)/x_sym^2 * U_expr + q_sym^2 * U_expr
ODE_rhs = q_sym^2 * Pi0_sym * x_sym^(l_test + 1)
residual = simplify(ODE_lhs - ODE_rhs)
println("ODE residual (should be 0): ", residual)

# ASSERTION 5: verify residual == 0 both symbolically AND numerically (a
# concrete-value check catches cases where simplify() can't fully collapse
# a transcendental expression that IS actually zero).
symbolic_zero = isequal(residual, 0)
numeric_zero = let
    f = Symbolics.build_function(residual, x_sym, q_sym, C_sym, Pi0_sym; expression=false)
    all(abs(f(xv, qv, 1.3, -0.7)) < 1e-9 for xv in (0.2, 0.5, 0.9), qv in (1.1, 2.7, 4.3))
end
@assert symbolic_zero || numeric_zero "Reid ODE not satisfied: $residual"
println("ASSERTION 5 OK: U(x) = C*x*j_2(qx) + Pi0*x^3 satisfies Reid's ODE",
        symbolic_zero ? " (symbolically)" : " (verified numerically at several (x,q,C,Pi0))")

println()
println("--- Boundary conditions: solve the 2x2 system for C, Pi0 ---")
@variables jl_sym jlp_sym l_s
A_mat = [jl_sym                                             1;
         -q_sym^2*jl_sym + 2*(l_s^2+l_s-1)*jl_sym - 2*q_sym*jlp_sym    2*(l_s^2-1)]
b_vec = [-1, 0]
sol = Symbolics.simplify.(A_mat \ b_vec)
C_sol, Pi0_sol = sol[1], sol[2]

bc1_pi0 = simplify(C_sol * jl_sym + Pi0_sol + 1)
@assert isequal(bc1_pi0, 0)
println("ASSERTION 6 OK: C*j_l + Pi0 = -1 (BC1 verified for Pi0)")

C_reid = 2*(l_s-1)*(l_s+1) / ((2*l_s - q_sym^2)*jl_sym - 2*q_sym*jlp_sym)
@assert numerically_equal(C_sol, C_reid)
println("ASSERTION 7 OK: BC solution for C matches Reid eq. (38)")

@assert numerically_equal(C_sol * jl_sym + Pi0_sol + 1, 0)
println("ASSERTION 8 OK: BC1 satisfied: U(1) = -1")

println()
println("--- Section 3: velocity field incompressibility ---")
@variables theta_sym

# V(x) = (x^2 U)' / (l(l+1) x), verify continuity: (x^2 U)'/x^2 - l(l+1)*V/x = 0
l_num = 2
V_from_U(Uexpr) = Symbolics.derivative(x_sym^2 * Uexpr, x_sym) / (l_num*(l_num+1)*x_sym)

@variables Utest(x_sym)
Vtest = V_from_U(Utest)
div_check = Symbolics.derivative(x_sym^2*Utest, x_sym)/x_sym^2 - l_num*(l_num+1)*Vtest/x_sym
@assert numerically_equal(div_check, 0)
println("ASSERTION 9 OK: velocity field u_r=U*P_l, u_theta=V*(dP_l/dtheta) is incompressible")
println("V(x) = ", Vtest)

println()
println("--- Section 4: strain rate tensor (l=2) ---")
Pl_expr = (3*cos(theta_sym)^2 - 1)/2
dPl_dth = Symbolics.derivative(Pl_expr, theta_sym)
d2Pl_dth = Symbolics.derivative(dPl_dth, theta_sym)

@variables f_sym(x_sym) g_sym(x_sym)
fxp = Symbolics.derivative(f_sym, x_sym)
gxp = Symbolics.derivative(g_sym, x_sym)

e_rr = fxp * Pl_expr
e_rth = (1//2) * (gxp - g_sym/x_sym + f_sym/x_sym) * dPl_dth
e_thth = (f_sym/x_sym)*Pl_expr + (g_sym/x_sym)*d2Pl_dth
e_phph = (f_sym/x_sym)*Pl_expr + (g_sym*cos(theta_sym)/sin(theta_sym)/x_sym)*dPl_dth
gdot_sq = 2*(e_rr^2 + e_thth^2 + e_phph^2 + 2*e_rth^2)
println("Strain rate components defined (l=2). Proceeding to angular integration...")

println()
println("--- Section 5: angular integration (Legendre orthogonality), via QuadGK ---")
# Numeric definite integrals over theta in [0,pi] -- Symbolics.jl has no
# general symbolic integrate(); QuadGK.jl gives a fully rigorous numeric
# check of the known closed-form Legendre-orthogonality identities.
Pl_f  = Symbolics.build_function(Pl_expr, theta_sym; expression=false)
dPl_f = Symbolics.build_function(dPl_dth, theta_sym; expression=false)

I_Pl2, _  = quadgk(th -> Pl_f(th)^2 * sin(th), 0, pi)
I_dPl2, _ = quadgk(th -> dPl_f(th)^2 * sin(th), 0, pi)

l_check = 2
@assert abs(I_Pl2 - 2/(2*l_check+1)) < 1e-10
println("ASSERTION 10 OK: int P_2^2 sin(theta) dtheta = 2/(2l+1) = ", 2/(2*l_check+1), " (quadgk: $I_Pl2)")
@assert abs(I_dPl2 - 2*l_check*(l_check+1)/(2*l_check+1)) < 1e-10
println("ASSERTION 11 OK: int (dP_2/dtheta)^2 sin(theta) dtheta = 2l(l+1)/(2l+1) = ",
        2*l_check*(l_check+1)/(2*l_check+1), " (quadgk: $I_dPl2)")

println()
println("="^78)
println("Section 6: The Correction Integral Gamma_l^(a) -- generalized")
println("="^78)
println("""
Gamma_l^(a) = [int_0^1 (int_0^pi gammadot^(a+2) sin(theta) dtheta) x^2 dx] / N_l^((a+2)/2)

generalizing the old (a=2-hardcoded) Gamma_l's N_l^2 normalization, since
gammadot^(a+2) ~ U^(a+2) while N_l ~ U^2, so N_l^((a+2)/2) ~ U^(a+2) -- the
same power, keeping Gamma_l^(a) dimensionless in powers of U for any a.

Inviscid limit: U(x) = Pi0*x^(l+1) = -x^(l+1) makes every strain component
proportional to x^l, so gammadot^2(x,theta) = x^(2l)*H(theta) EXACTLY -- this
factorization holds for any power you then raise it to, not just squared.
""")

l_val = 2
Pi0_inviscid = -1
f_inviscid = Pi0_inviscid * x_sym^(l_val+1)
fp_inviscid = Symbolics.derivative(f_inviscid, x_sym)

@assert isequal(substitute(f_inviscid, Dict(x_sym => 1)), -1)
println("ASSERTION 12 OK: BC1 satisfied in inviscid limit: U(1) = -1")

N_l_inviscid = simplify((f_inviscid^2 * x_sym^2 |> expr -> begin
    1//9   # int_0^1 x^8 dx = 1/9 (f_inviscid^2 = x^6, times x^2 Jacobian = x^8)
end))
println("N_l (inviscid) = int U^2 x^2 dx = ", N_l_inviscid)

g_inv_l2 = simplify(Symbolics.derivative(x_sym^2 * f_inviscid, x_sym) / (l_val*(l_val+1)*x_sym))
gp_inv_l2 = Symbolics.derivative(g_inv_l2, x_sym)
e_rr_inv   = fp_inviscid * Pl_expr
e_rth_inv  = (1//2)*(gp_inv_l2 - g_inv_l2/x_sym + f_inviscid/x_sym)*dPl_dth
e_thth_inv = (f_inviscid/x_sym)*Pl_expr + (g_inv_l2/x_sym)*d2Pl_dth
e_phph_inv = (f_inviscid/x_sym)*Pl_expr + (g_inv_l2*cos(theta_sym)/(sin(theta_sym)*x_sym))*dPl_dth
gdot_sq_inv = 2*(e_rr_inv^2 + e_thth_inv^2 + e_phph_inv^2 + 2*e_rth_inv^2)

H_theta = simplify(substitute(gdot_sq_inv, Dict(x_sym => 1)))
println("H(theta) = gdot_sq / x^(2l) at x=1 = ", H_theta)

factorization_residual = simplify(gdot_sq_inv - x_sym^4 * H_theta)
symbolic_zero_13 = isequal(factorization_residual, 0)
numeric_zero_13 = let
    f = Symbolics.build_function(factorization_residual, x_sym, theta_sym; expression=false)
    all(abs(f(xv, thv)) < 1e-9 for xv in (0.1,0.4,0.7,0.99), thv in (0.3,1.1,2.0,2.9))
end
@assert symbolic_zero_13 || numeric_zero_13 "Factorization failed: $factorization_residual"
println("ASSERTION 13 OK: gdot_sq = x^(2l) * H(theta) exactly (inviscid l=2)")

H_f = Symbolics.build_function(H_theta, theta_sym; expression=false)
H_poly_ref(th) = 3*cos(th)^4 + 11*cos(th)^2 + 13
@assert all(abs(H_f(th) - H_poly_ref(th)) < 1e-10 for th in range(0.01, pi-0.01, length=20))
println("ASSERTION 14 OK: H(theta) = 3cos^4(theta) + 11cos^2(theta) + 13 confirmed")

function Gamma_l_a_inviscid_l2(a_val)
    exponent = (a_val + 2) / 2
    ang, _ = quadgk(th -> H_poly_ref(th)^exponent * sin(th), 0, pi)
    radial_denom = l_val*(a_val + 2) + 3
    (ang / radial_denom) / (1/9)^exponent
end

G2_a2 = Gamma_l_a_inviscid_l2(2.0)
target = 1783566 / 385
println("Gamma_2^(a=2) computed = ", G2_a2, "   old notebook's exact value = ", target)
@assert abs(G2_a2 - target)/target < 1e-9
println("ASSERTION 15 OK: Gamma_l^(a) reduces EXACTLY to the old (Python) notebook's",
        " sympy-exact Gamma_2 = 1783566/385 at a=2")

println()
println("Gamma_2^(a) for several shape exponents a (inviscid limit):")
for a_val in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
    println("  a=$a_val: Gamma_2^(a) = ", round(Gamma_l_a_inviscid_l2(a_val), digits=4))
end

println()
println("="^78)
println("Section 7: The Modified Mode Equation, and Its Effective Damping")
println("="^78)
println("""
Generalized dissipation function: delta_Phi = -(coeff)*|bdot|^(a+2), giving a
generalized dissipative force d(delta_Phi)/d(bdot) ~ |bdot|^a * bdot --
verified below by direct NUMERICAL differentiation (the exact step an
earlier, unrelated attempt at this same generalization for the Cross model
got wrong; checked carefully here for the same reason).
""")

function check_force_order(a_val; x0=0.37, h=1e-6)
    f(xv) = -abs(xv)^(a_val+2)
    dfdx = (f(x0+h) - f(x0-h)) / (2h)
    claimed = -(a_val+2)*abs(x0)^a_val*x0
    dfdx, claimed
end

mismatches = Float64[]
for a_val in (0.5, 1.0, 1.5, 2.0, 3.0)
    d, c = check_force_order(a_val)
    matchv = abs(d-c) < 1e-4
    println("a=$a_val: numeric d(delta_Phi)/d(bdot)=$d  claimed=-(a+2)|bdot|^a*bdot=$c  match=$matchv")
    matchv || push!(mismatches, a_val)
end
@assert isempty(mismatches)
println("ASSERTION 16 OK: generalized force ~ |bdot|^a * bdot verified numerically")

@variables Oh_sym Lambda_sym eps_ST_sym Gamma_l_sym bdot_sym
D0_newtonian = 2*Oh_sym*(l_val-1)*(2*l_val+1)
D_eff = D0_newtonian * (1 - eps_ST_sym * Lambda_sym^a_shape * Gamma_l_sym * abs(bdot_sym)^a_shape)
println("Newtonian damping D_l^(0) [l=2] = ", D0_newtonian)
println("Effective damping D_eff = ", D_eff)

# ASSERTION 17: at a=2, D_eff MUST reduce exactly to the old notebook's D_eff:
# D0*(1 - eps_ST*Lambda^2*Gamma_l*bdot^2). Since |bdot|^2 == bdot^2 for real
# bdot, substitute abs(bdot)->sqrt(bdot^2) to compare on equal footing.
D_eff_a2 = simplify(substitute(D_eff, Dict(a_shape => 2, abs(bdot_sym) => sqrt(bdot_sym^2))))
D_eff_original = D0_newtonian * (1 - eps_ST_sym * Lambda_sym^2 * Gamma_l_sym * bdot_sym^2)
match17 = numerically_equal(D_eff_a2, D_eff_original)
println("D_eff at a=2: ", D_eff_a2)
println("Original D_eff: ", D_eff_original)
@assert match17
println("ASSERTION 17 OK: generalized D_eff reduces EXACTLY to the old notebook's D_eff at a=2")

@assert numerically_equal(substitute(D_eff, Dict(eps_ST_sym => 0)), D0_newtonian)
println("ASSERTION 18 OK: eps_ST=0 recovers Newtonian damping for any a")

println()
println("="^78)
println("Section 8: Secular Terms and the Generalized Slow-Amplitude Equation")
println("="^78)
println("""
The old (a=2-hardcoded) notebook's cubic damping term projects onto
sin(omega*t) via <sin^4>=3/8, giving the 3/4 secular factor. For the
generalized |bdot|^a*bdot damping (non-polynomial for non-even a), the
analogous projection is a classical Wallis-type integral, closed-form via
the Gamma function for any real a:
  C(a) = (2/sqrt(pi)) * Gamma((a+3)/2) / Gamma((a+4)/2)
verified below against QuadGK direct quadrature (not just cited).
""")

C_of_a(a_val) = (2/sqrt(pi)) * gamma((a_val+3)/2) / gamma((a_val+4)/2)

mismatches2 = Float64[]
for a_val in (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)
    closed = C_of_a(a_val)
    numeric, _ = quadgk(th -> sin(th)^(a_val+2), 0, pi)
    numeric *= 2/pi
    err = abs(closed - numeric)
    println("a=$a_val: closed-form=$closed  numeric-quadrature=$numeric  |diff|=$err")
    err > 1e-9 && push!(mismatches2, a_val)
end
@assert isempty(mismatches2)
println("ASSERTION 19 OK: closed-form C(a) matches direct quadrature for all a tested")

@assert abs(C_of_a(2.0) - 0.75) < 1e-12
println("ASSERTION 20 OK: C(2) = 3/4, exactly matching the old notebook's validated secular factor")

# Cross-check against a direct symbolic period integral (mirrors the old
# notebook's explicit sin^4 integral exactly, now via QuadGK since Symbolics
# has no integrate()):
coeff_34, _ = quadgk(t -> sin(t)^4, 0, 2pi)
coeff_34 /= pi
@assert abs(coeff_34 - 0.75) < 1e-10
println("ASSERTION 21 OK: direct <sin^4> period integral (quadgk) still gives exactly 3/4")

println()
println("Assembling the generalized slow-amplitude equation:")
println("da/dt = -gamma0*a*[1 - C(a)*eps_ST*Lambda^a*Gamma_l^(a)*a^a]")

# ASSERTION 22: at a=2, MUST collapse to the old boxed equation
# da/dt = -gamma0*a*(1 - (3/4)*eps_ST*Lambda^2*Gamma_l*a^2). Verify at several
# concrete numeric points rather than symbolically (avoids fragile gamma()
# symbolic simplification), which is an equally rigorous check.
function da_dt_general(gamma0, a_amp, eps_ST, Lambda, Gamma_l, a_val)
    -gamma0 * a_amp * (1 - C_of_a(a_val) * eps_ST * Lambda^a_val * Gamma_l * a_amp^a_val)
end
function da_dt_original(gamma0, a_amp, eps_ST, Lambda, Gamma_l)
    -gamma0 * a_amp * (1 - 0.75 * eps_ST * Lambda^2 * Gamma_l * a_amp^2)
end

mismatches3 = []
for (gamma0, a_amp, eps_ST, Lambda, Gamma_l) in
        [(1.0, 0.3, 0.1, 0.2, 5.0), (0.5, 0.1, 0.05, 0.15, 100.0), (2.0, 0.4, 0.02, 0.3, 4632.639)]
    gen = da_dt_general(gamma0, a_amp, eps_ST, Lambda, Gamma_l, 2.0)
    orig = da_dt_original(gamma0, a_amp, eps_ST, Lambda, Gamma_l)
    err = abs(gen - orig)
    err > 1e-9 && push!(mismatches3, (gamma0, a_amp, eps_ST, Lambda, Gamma_l, err))
end
@assert isempty(mismatches3) "mismatches: $mismatches3"
println("ASSERTION 22 OK: generalized slow-amplitude equation collapses EXACTLY to the")
println("old notebook's boxed equation at a=2, verified at several concrete parameter points")

for (gamma0, a_amp, Lambda, Gamma_l) in [(1.0, 0.3, 0.2, 5.0), (2.0, 0.1, 0.3, 100.0)]
    @assert abs(da_dt_general(gamma0, a_amp, 0.0, Lambda, Gamma_l, 2.5) - (-gamma0*a_amp)) < 1e-12
end
println("ASSERTION 23 OK: eps_ST=0 gives da/dt = -gamma0*a (Newtonian) for any a")

println()
println("Boxed result:")
println("  da/dt = -gamma_l^(0) * a * [1 - C(a)*eps_ST*Lambda^a*Gamma_l^(a)*a^a]")
println("  C(a) = (2/sqrt(pi)) * Gamma((a+3)/2)/Gamma((a+4)/2),  C(2) = 3/4")
println("well-posed (correction -> 0 as amplitude -> 0) for every a>0, and reduces")
println("to exactly the old Carreau-only boxed equation at a=2.")

println()
println("="^78)
println("Section 6 (continued): Gamma_l^(a) for l=2,3,4,5, inviscid limit")
println("="^78)

# f_inv = -x^(l+1); build strain components the same way as l=2 above, for general l.
function Gamma_l_a_general_l(l_val_, a_val)
    fp_c(xv) = -(l_val_+1)*xv^l_val_
    f_c(xv)  = -xv^(l_val_+1)
    Plfun(u) = begin   # Legendre polynomial via Bonnet's recursion (u = cos(theta))
        P0, P1 = one(u), u
        l_val_ == 0 && return P0
        l_val_ == 1 && return P1
        Pkm1, Pk = P0, P1
        for k in 1:(l_val_-1)
            Pkp1 = ((2k+1)*u*Pk - k*Pkm1) / (k+1)
            Pkm1, Pk = Pk, Pkp1
        end
        Pk
    end
    function H_at_theta(thv)
        xv = 1.0
        u = cos(thv); s = sin(thv)
        h = 1e-6   # numeric differentiation of the Legendre function in theta (small h)
        Pl_here = Plfun(cos(thv))
        dPl_here = (Plfun(cos(thv+h)) - Plfun(cos(thv-h))) / (2h)
        d2Pl_here = (Plfun(cos(thv+h)) - 2*Plfun(cos(thv)) + Plfun(cos(thv-h))) / h^2
        g_c  = (2*f_c(xv) + xv*fp_c(xv)) / (l_val_*(l_val_+1))
        fpp_c = -(l_val_+1)*l_val_*xv^(l_val_-1)
        gp_c = (3*fp_c(xv) + xv*fpp_c) / (l_val_*(l_val_+1))
        e_rr   = fp_c(xv) * Pl_here
        e_rth  = 0.5*(gp_c - g_c/xv + f_c(xv)/xv)*dPl_here
        e_thth = (f_c(xv)/xv)*Pl_here + (g_c/xv)*d2Pl_here
        e_phph = (f_c(xv)/xv)*Pl_here + (g_c*u/(s*xv))*dPl_here
        2*(e_rr^2 + e_thth^2 + e_phph^2 + 2*e_rth^2)
    end
    exponent = (a_val+2)/2
    ang, _ = quadgk(th -> H_at_theta(th)^exponent * sin(th), 1e-6, pi-1e-6)
    radial_denom = l_val_*(a_val + 2) + 3
    N_l_val = 1/(2*l_val_ + 5)
    (ang/radial_denom) / N_l_val^exponent
end

println("Gamma_l^(a) (inviscid limit) for l=2,3,4,5, a=1,2,3:")
for l_val_ in (2,3,4,5)
    row = [Gamma_l_a_general_l(l_val_, a_val) for a_val in (1.0, 2.0, 3.0)]
    println("  l=$l_val_:  a=1: ", round(row[1],digits=2), "  a=2: ", round(row[2],digits=2),
            "  a=3: ", round(row[3],digits=2))
end

G_l2_a2 = Gamma_l_a_general_l(2, 2.0)
rel_err_24 = abs(G_l2_a2 - target)/target
println()
println("Cross-check: general-l function at l=2,a=2 -> $G_l2_a2  vs exact $target  rel_err=$rel_err_24")
@assert rel_err_24 < 1e-4   # numeric theta-derivatives (finite difference) limit precision vs the exact symbolic l=2 result
println("ASSERTION 24 OK: general-l Gamma_l^(a) matches the exact l=2,a=2 value")

for l_val_ in (2,3,4,5), a_val in (0.5, 1.0, 1.5, 2.0, 3.0)
    @assert Gamma_l_a_general_l(l_val_, a_val) > 0
end
println("ASSERTION 25 OK: Gamma_l^(a) > 0 for l=2..5, a in {0.5,1,1.5,2,3} (inviscid limit)")

println()
println("="^78)
println("Section 6.5: Finite-Oh Correction Integral Gamma_l^(a)(Oh)")
println("="^78)
println("""
At finite Oh, the eigenvalue q acquires a negative imaginary part and the
eigenfunction becomes complex, U=U_R+i*U_I. The old (a=2-hardcoded)
notebook computed <gammadot^4> via a closed-form algebraic identity specific
to the exponent 4; for general a we compute <[...]^((a+2)/2)> by DIRECT
NUMERICAL time-averaging over one period -- fully general, verified below to
reproduce the old closed-form result to high precision at a=2 (the only
case where both methods are available to compare).

Normalization: the old Gamma_l is defined WITHOUT any time-average baked
in (a purely spatial integral; §7-8 apply the temporal average separately).
The old code divides its raw time-averaged quantity by 3/8=<sin^4> "to match
the inviscid convention" (its own words). We generalize that normalization
to <sin^(a+2)> = Gamma((a+3)/2) / [sqrt(pi)*Gamma((a+4)/2)] -- verified below
to equal 3/8 exactly at a=2.
""")

function sph_bessel_j(l, z)
    sqrt(pi/(2z)) * besselj(l + 0.5, z)
end

function bessel_ratio(l, z)
    sph_bessel_j(l+1, z) / sph_bessel_j(l, z)
end

function reid_char(qv, Oh, l)
    alpha2 = sqrt(l*(l-1)*(l+2)) / Oh
    Qv = bessel_ratio(l, qv)
    lhs = alpha2^2/qv^4 + 1
    rhs = 2*(l-1)/qv^2 * (l + (l+1)*(qv - 2*l*Qv)/(qv - 2*Qv))
    lhs - rhs
end

function find_eigenvalue(Oh_val, l=2; maxiter=200, tol=1e-13)
    sigma0 = sqrt(l*(l-1)*(l+2))
    gamma0 = (l-1)*(2l+1) * Oh_val
    q = sqrt(complex(gamma0/Oh_val, -sigma0/Oh_val))
    imag(q) > 0 && (q = -q)
    for _ in 1:maxiter
        F = reid_char(q, Oh_val, l)
        h = max(1e-7, abs(q)*1e-7)
        dF = (reid_char(q + h*im, Oh_val, l) - F) / (h*im)
        step = F/dF
        q -= step
        abs(F) < tol && break
    end
    imag(q) > 0 ? conj(q) : q
end

function eigenfunction(x_arr, qv, l=2)
    jlq  = sph_bessel_j(l, qv)
    jl1q = sph_bessel_j(l+1, qv)
    Cv   = 2*(l^2-1) / (qv*(2*jl1q - qv*jlq))
    Pi0v = -1 - Cv*jlq
    n = length(x_arr)
    U = zeros(ComplexF64, n); dU = zeros(ComplexF64, n); d2U = zeros(ComplexF64, n)
    for (i, xv) in enumerate(x_arr)
        xv < 1e-12 && continue
        jlv  = sph_bessel_j(l, qv*xv)
        jl1v = sph_bessel_j(l+1, qv*xv)
        jl2v = sph_bessel_j(l+2, qv*xv)
        U[i]   = Cv*xv*jlv + Pi0v*xv^(l+1)
        dU[i]  = Cv*((l+1)*jlv - qv*xv*jl1v) + Pi0v*(l+1)*xv^l
        d2U[i] = Cv*(l*(l+1)/xv*jlv - (2l+3)*qv*jl1v + qv^2*xv*jl2v) + Pi0v*l*(l+1)*xv^(l-1)
    end
    U, dU, d2U
end

function legendre_arrays(l, u_pts)
    n = length(u_pts)
    Pl = zeros(n); Pl1 = zeros(n)
    for (i,u) in enumerate(u_pts)
        P0, P1 = 1.0, u
        if l == 0
            Pl[i] = P0
        elseif l == 1
            Pl[i] = P1
        else
            pkm1, pk = P0, P1
            for k in 1:(l-1)
                pkp1 = ((2k+1)*u*pk - k*pkm1)/(k+1)
                pkm1, pk = pk, pkp1
            end
            Pl[i] = pk
        end
        if l-1 == 0   # P_{l-1}
            Pl1[i] = 1.0
        else
            pkm1, pk = 1.0, u
            for k in 1:(l-2)
                pkp1 = ((2k+1)*u*pk - k*pkm1)/(k+1)
                pkm1, pk = pk, pkp1
            end
            Pl1[i] = pk
        end
    end
    dPl = similar(Pl); d2Pl = similar(Pl); sin_th = similar(Pl)
    for i in 1:n
        u = u_pts[i]
        sin_th[i] = sqrt(max(1-u^2, 0.0))
        if abs(u) < 1 - 1e-14
            dPl[i] = l*(Pl1[i] - u*Pl[i])/(1-u^2)
        else
            dPl[i] = 0.0
        end
    end
    for i in 1:n
        u = u_pts[i]
        if abs(u) < 1 - 1e-14
            d2Pl[i] = (2*u*dPl[i] - l*(l+1)*Pl[i])/(1-u^2)
        else
            d2Pl[i] = 0.0
        end
    end
    Pl, dPl, d2Pl, sin_th
end

function gdot_sq_matrix(f_arr, df_arr, x_arr, l, u_pts, d2f_arr)
    nx, nth = length(x_arr), length(u_pts)
    g_arr = (2 .* f_arr .+ x_arr .* df_arr) ./ (l*(l+1))
    dg_arr = (3 .* df_arr .+ x_arr .* d2f_arr) ./ (l*(l+1))
    Pl, dPl, d2Pl, sin_th = legendre_arrays(l, u_pts)
    H = zeros(nx, nth)
    for i in 1:nx
        xv = x_arr[i]
        xv < 1e-12 && continue
        f_, fp, g_, gp = f_arr[i], df_arr[i], g_arr[i], dg_arr[i]
        for j in 1:nth
            e_rr   = fp * Pl[j]
            e_thth = (f_/xv)*Pl[j] - (g_/xv)*(u_pts[j]*dPl[j] - sin_th[j]^2*d2Pl[j])
            e_phph = (f_/xv)*Pl[j] - (g_/xv)*(u_pts[j]*dPl[j])
            e_rth  = 0.5*sin_th[j]*dPl[j]*(g_/xv - gp - f_/xv)
            H[i,j] = 2*(e_rr^2 + e_thth^2 + e_phph^2 + 2*e_rth^2)
        end
    end
    H
end

function gdot_cross_matrix(f1, df1, d2f1, f2, df2, d2f2, x_arr, l, u_pts)
    nx, nth = length(x_arr), length(u_pts)
    g1 = (2 .* f1 .+ x_arr .* df1) ./ (l*(l+1)); dg1 = (3 .* df1 .+ x_arr .* d2f1) ./ (l*(l+1))
    g2 = (2 .* f2 .+ x_arr .* df2) ./ (l*(l+1)); dg2 = (3 .* df2 .+ x_arr .* d2f2) ./ (l*(l+1))
    Pl, dPl, d2Pl, sin_th = legendre_arrays(l, u_pts)
    Hc = zeros(nx, nth)
    for i in 1:nx
        xv = x_arr[i]
        xv < 1e-12 && continue
        f1_,fp1,g1_,gp1 = f1[i],df1[i],g1[i],dg1[i]
        f2_,fp2,g2_,gp2 = f2[i],df2[i],g2[i],dg2[i]
        for j in 1:nth
            e_rr1   = fp1 * Pl[j]
            e_thth1 = (f1_/xv)*Pl[j] - (g1_/xv)*(u_pts[j]*dPl[j] - sin_th[j]^2*d2Pl[j])
            e_phph1 = (f1_/xv)*Pl[j] - (g1_/xv)*(u_pts[j]*dPl[j])
            e_rth1  = 0.5*sin_th[j]*dPl[j]*(g1_/xv - gp1 - f1_/xv)
            e_rr2   = fp2 * Pl[j]
            e_thth2 = (f2_/xv)*Pl[j] - (g2_/xv)*(u_pts[j]*dPl[j] - sin_th[j]^2*d2Pl[j])
            e_phph2 = (f2_/xv)*Pl[j] - (g2_/xv)*(u_pts[j]*dPl[j])
            e_rth2  = 0.5*sin_th[j]*dPl[j]*(g2_/xv - gp2 - f2_/xv)
            Hc[i,j] = 2*(e_rr1*e_rr2 + e_thth1*e_thth2 + e_phph1*e_phph2 + 2*e_rth1*e_rth2)
        end
    end
    Hc
end

# Gauss-Legendre nodes/weights via Jacobi matrix eigendecomposition (no
# external quadrature-rule package needed -- matches julia/src/types.jl's
# own make_theta_vec approach).
using LinearAlgebra
function gauss_legendre(n)
    bvec = [k / sqrt(4*k^2 - 1) for k in 1:n-1]
    J = diagm(1 => bvec, -1 => bvec)
    F = eigen(J)
    nodes = F.values
    weights = 2 .* (F.vectors[1, :]).^2
    nodes, weights
end

norm_p(p) = gamma(p+0.5) / (sqrt(pi)*gamma(p+1))

function compute_gamma_l_a(Oh_val, a_val, l=2; n_x=300, n_theta=100, n_phi=128)
    qv = find_eigenvalue(Oh_val, l)
    sigma = qv^2 * Oh_val
    gamma_, omega_ = real(sigma), -imag(sigma)

    x_arr = range(1e-4, 1.0, length=n_x)
    U, dU, d2U = eigenfunction(x_arr, qv, l)
    u_pts, wts = gauss_legendre(n_theta)
    fR, dfR, d2fR = real.(U), real.(dU), real.(d2U)
    fI, dfI, d2fI = imag.(U), imag.(dU), imag.(d2U)
    H_R = gdot_sq_matrix(fR, dfR, x_arr, l, u_pts, d2fR)
    H_I = gdot_sq_matrix(fI, dfI, x_arr, l, u_pts, d2fI)
    H_cross = gdot_cross_matrix(fR, dfR, d2fR, fI, dfI, d2fI, x_arr, l, u_pts)

    p = (a_val + 2)/2
    time_avg = zeros(n_x, n_theta)
    for k in 1:n_phi
        phv = 2pi*(k-1)/n_phi
        sinp, cosp = sin(phv), cos(phv)
        @. time_avg += abs(H_R*sinp^2 + H_I*cosp^2 + 2*H_cross*sinp*cosp)^p
    end
    time_avg ./= n_phi
    time_avg ./= norm_p(p)

    numerator = trapz(x_arr, (time_avg * wts) .* x_arr.^2)
    N_l = trapz(x_arr, (abs.(U).^2) .* x_arr.^2)
    numerator / N_l^p, gamma_, omega_
end

function trapz(xs, ys)
    s = 0.0
    for i in 1:length(xs)-1
        s += (xs[i+1]-xs[i]) * (ys[i+1]+ys[i]) / 2
    end
    s
end

function compute_gamma_l_original(Oh_val, l=2; n_x=300, n_theta=100)
    qv = find_eigenvalue(Oh_val, l)
    sigma = qv^2 * Oh_val
    gamma_, omega_ = real(sigma), -imag(sigma)
    x_arr = range(1e-4, 1.0, length=n_x)
    U, dU, d2U = eigenfunction(x_arr, qv, l)
    u_pts, wts = gauss_legendre(n_theta)
    H_R = gdot_sq_matrix(real.(U), real.(dU), x_arr, l, u_pts, real.(d2U))
    H_I = gdot_sq_matrix(imag.(U), imag.(dU), x_arr, l, u_pts, imag.(d2U))
    H_cross = gdot_cross_matrix(real.(U),real.(dU),real.(d2U), imag.(U),imag.(dU),imag.(d2U), x_arr, l, u_pts)
    H_sq = @. H_R^2 + H_I^2 + (2/3)*H_R*H_I + (4/3)*H_cross^2
    numerator = trapz(x_arr, (H_sq * wts) .* x_arr.^2)
    N_l = trapz(x_arr, (abs.(U).^2) .* x_arr.^2)
    numerator / N_l^2, gamma_, omega_
end

println("Comparing generalized (numerical time-average) vs. old (closed-form) at a=2:")
mismatches4 = Float64[]
for Oh_v in (0.001, 0.05, 0.1, 0.3)
    G_gen, gam_gen, om_gen = compute_gamma_l_a(Oh_v, 2.0, 2)
    G_orig, gam_orig, om_orig = compute_gamma_l_original(Oh_v, 2)
    rel_err = abs(G_gen - G_orig)/G_orig
    println("Oh=$Oh_v: generalized=$(round(G_gen,digits=6))  original=$(round(G_orig,digits=6))  rel_err=$rel_err")
    rel_err > 1e-6 && push!(mismatches4, Oh_v)
end
@assert isempty(mismatches4)
println("ASSERTION 26 OK: generalized finite-Oh Gamma_l^(a) matches the old notebook's")
println("compute_gamma_l to high precision at a=2, for every Oh tested")

@assert abs(norm_p(2.0) - 0.375) < 1e-12
println("ASSERTION 27 OK: norm_p(2) = 3/8 exactly, matching the old notebook's convention")

Gamma_inv_ref = 1783566/385
println()
println("Oh sweep, Gamma_l^(a)(Oh) for a=1,2,3:")
sweep_results_a2 = Tuple{Float64,Float64,Float64,Float64}[]
for Oh_v in (0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)
    row = Float64[]
    for a_val in (1.0, 2.0, 3.0)
        G, gam, om = compute_gamma_l_a(Oh_v, a_val, 2)
        push!(row, G)
        a_val == 2.0 && push!(sweep_results_a2, (Oh_v, G, gam, om))
    end
    println("  Oh=$Oh_v:  a=1: $(round(row[1],digits=2))  a=2: $(round(row[2],digits=2))  a=3: $(round(row[3],digits=2))")
end

for a_val in (1.0, 2.0, 3.0)
    vals = [compute_gamma_l_a(Oh_v, a_val, 2)[1] for Oh_v in (0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.5)]
    @assert all(v -> v > 0, vals)
    @assert all(i -> vals[i] > vals[i+1], 1:length(vals)-1)
end
println("ASSERTION 28 OK: Gamma_l^(a)(Oh) positive and strictly decreasing with Oh, for a in {1,2,3}")

G_001 = sweep_results_a2[1][2]
rel_corr_001 = (Gamma_inv_ref - G_001)/Gamma_inv_ref
@assert 0.04 < rel_corr_001 < 0.12
println("ASSERTION 29 OK: Oh=0.001 viscous correction (a=2) = $rel_corr_001, in the")
println("  (4%-12%) band the old notebook established")

gamma_03 = first(r[3] for r in sweep_results_a2 if abs(r[1]-0.3) < 1e-9)
gamma_Lamb = (2-1)*(2*2+1)*0.3
rel_diff_Lamb = abs(gamma_03 - gamma_Lamb)/gamma_Lamb
@assert rel_diff_Lamb < 0.30
println("ASSERTION 30 OK: gamma(Oh=0.3, a=2) = $gamma_03 within 30% of Lamb ($gamma_Lamb)")

println()
println("="^78)
println("SUMMARY: 30/30 assertions passed. This is a strict generalization of the")
println("Carreau-only derivation -- every place the old derivation could be checked")
println("against this one at a=2, the match is exact or near-machine-precision.")
println()
println("Implemented in julia/src/st_extension.jl and julia/src/types.jl:")
println("  STParams gains an a::Float64 field (default 2.0, backward-compatible).")
println("  The two places the code hardcoded the exponent 2")
println("  (Gamma_eff = Gamma .* (lambda_c.*sigma0).^2 and")
println("   shear_sq_lag = sum(Gamma_eff.*Adot_prev.^2)) now use .^a, with abs()")
println("  added on the second since a need not be an even integer.")
println()
println("Ready for real experimental data: convert Cross-model fit parameters")
println("(K, m, mu0, mu_infty) via lambda_c=K, a=m, n=1-m,")
println("eps_ST=[(mu0-mu_infty)/mu0]*(1-n)/a, compute Gamma_l^(a)(Oh) via this")
println("script's Section 6.5 machinery, and STParams/build_residual_st! are ready.")
println("="^78)
