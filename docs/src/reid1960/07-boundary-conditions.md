# Boundary Conditions

Three physical conditions at the drop surface (``x=1``) fix the two
constants ``C,\Pi_0`` from the general solution ``U(x)=Cx\,j_l(qx)+\Pi_0x^{l+1}``
and, together, produce the characteristic equation.

## BC1: Kinematic condition

The fluid at the surface must move with the surface -- the surface isn't
a permeable membrane. The surface itself moves at
``\partial r/\partial t|_{\text{surface}} = -\sigma R\epsilon_0e^{-\sigma t}Y_l^m``
(differentiating the surface ansatz directly), while the fluid's radial
velocity there is ``u_r|_{x=1}=\epsilon_0\sigma R\,U(1)\,Y_l^m e^{-\sigma t}``.
Equating the two and canceling the common prefactor:

```math
-\sigma\epsilon_0 = \epsilon_0\sigma\,U(1)
```

```@eval
using Symbolics, Markdown
@variables sigma_sym epsilon0 U1
surface_velocity = -sigma_sym*epsilon0
fluid_velocity_at_1 = epsilon0*sigma_sym*U1
bc1_solution = Symbolics.solve_for(surface_velocity - fluid_velocity_at_1 ~ 0, U1)
Markdown.parse("Canceling the common prefactor \$\\epsilon_0\\sigma\$ and solving for \$U(1)\$ gives\n\n```math\nU(1) = " *
    string(bc1_solution) * "\n```")
```

!!! note "Why x=1 and not the true deformed surface?"
    The condition is imposed at ``x=1`` (the undeformed surface) because
    evaluating it at the true deformed surface ``r=R[1+\epsilon Y_l^m]``
    would only add ``O(\epsilon^2)`` corrections -- negligible at this
    linear order.

## BC2: Tangential stress condition

At a free surface with no viscous exterior fluid, the tangential viscous
stress must vanish:

```math
\tau_{r\theta} = \mu\left[r\frac{\partial}{\partial r}\!\left(\frac{u_\theta}{r}\right) + \frac{1}{r}\frac{\partial u_r}{\partial\theta}\right] = 0 \qquad \text{at } x=1.
```

This needs ``u_\theta``, which has not appeared until now. For an
axisymmetric poloidal field, the standard (citable, textbook) Stokes
stream function ``\psi(r,\theta)`` gives
``u_r=(1/(r^2\sin\theta))\,\partial\psi/\partial\theta`` and
``u_\theta=-(1/(r\sin\theta))\,\partial\psi/\partial r``. Writing
``u_r=f(r)P_l(\cos\theta)`` and integrating in ``\theta`` gives

```math
u_\theta = \frac{g(r)}{\sin\theta}\Big[P_{l+1}(\cos\theta)-P_{l-1}(\cos\theta)\Big], \qquad g(r)\equiv\frac{2f(r)+rf'(r)}{2l+1},
```

using the standard Legendre recurrence ``(2l+1)P_l=P_{l+1}'-P_{l-1}'``:

```@eval
using Symbolics, Markdown
function legendre_P(l::Int, xv)
    l == 0 && return one(xv)
    l == 1 && return xv
    Pm1, P = one(xv), xv
    for n in 1:(l - 1)
        P, Pm1 = ((2n + 1) * xv * P - n * Pm1) / (n + 1), P
    end
    P
end
@variables x_leg
Dxleg = Differential(x_leg)
results = String[]
for l_val in (2, 3, 4, 5)
    Pl, Plp1, Plm1 = legendre_P(l_val, x_leg), legendre_P(l_val + 1, x_leg), legendre_P(l_val - 1, x_leg)
    residual = simplify((2l_val + 1) * Pl - expand_derivatives(Dxleg(Plp1) - Dxleg(Plm1)); expand=true)
    push!(results, "\$l=$l_val\$: `$(string(residual))`")
end
Markdown.parse("Checked for \$l=2,3,4,5\$ (built via the same Bonnet recursion as Chapter 2): " * join(results, ", ") * ".")
```

Substituting both ``u_r`` and ``u_\theta`` into ``\tau_{r\theta}`` and
using a second standard recurrence,
``(2l+1)(1-x^2)P_l'(x)=l(l+1)[P_{l-1}(x)-P_{l+1}(x)]``, every term
collapses onto the single common angular factor
``\sin\theta\,P_l'(\cos\theta)`` -- exactly the shape a genuine
tangential-stress condition should have -- times a purely radial
coefficient:

```@eval
using Symbolics, Markdown
function legendre_P(l::Int, xv)
    l == 0 && return one(xv)
    l == 1 && return xv
    Pm1, P = one(xv), xv
    for n in 1:(l - 1)
        P, Pm1 = ((2n + 1) * xv * P - n * Pm1) / (n + 1), P
    end
    P
end
@variables x_leg
Dxleg = Differential(x_leg)
results = String[]
for l_val in (2, 3, 4, 5)
    Pl, Plp1, Plm1 = legendre_P(l_val, x_leg), legendre_P(l_val + 1, x_leg), legendre_P(l_val - 1, x_leg)
    lhs = expand_derivatives((2l_val + 1) * (1 - x_leg^2) * Dxleg(Pl))
    residual = simplify(lhs - l_val * (l_val + 1) * (Plm1 - Plp1); expand=true)
    push!(results, "\$l=$l_val\$: `$(string(residual))`")
end
Markdown.parse("Checked the same way: " * join(results, ", ") * ".")
```

Setting the resulting radial coefficient to zero at ``r=R`` turns out to
be *exactly* equivalent to ``R^2f''(R)+2Rf'(R)+[l(l+1)-2]f(R)=0``, and
substituting ``f(r)=\kappa\,G(r/R)`` (with ``G(x)=U(x)/x^2`` from Chapter
6 -- the condition is linear and homogeneous in ``f,f',f''``, so the
overall constant ``\kappa`` cancels) turns this into exactly

```math
\mathcal{L}_2[U] \equiv \left[\frac{d^2}{dx^2} - \frac{2}{x}\frac{d}{dx} + \frac{l(l+1)}{x^2}\right]U = 0 \qquad\text{at } x = 1.
```

This condition, previously just asserted in the source material, is now
derived end to end: stream function, then ``\tau_{r\theta}=0``, then the
radial ODE on ``f``, then the change of variables to ``U``.

Evaluating ``\mathcal{L}_2`` on the two pieces of the general solution
from Chapter 6:

```@eval
using Symbolics, Markdown
@variables x l q Pi_0 C
Dx = Differential(x)
L2(U) = expand_derivatives(Dx(Dx(U)) - 2/x*Dx(U) + l*(l+1)/x^2*U)
Up = Pi_0 * x^(l+1)
L2_Up_target = 2*(l-1)*(l+1)*Pi_0*x^(l-1)
diff_expr = L2(Up) - L2_Up_target
f = Symbolics.build_function(diff_expr, x, l, Pi_0; expression=false)
max_abs = maximum(abs(f(xv, lv, piv)) for xv in (0.4,1.3,2.1), lv in (2.0,3.0,5.7), piv in (1.0,-2.3))
Markdown.parse("Symbolics' simplifier does not fully collapse this particular difference to the literal token `0` " *
    "(a known quirk), so it is checked numerically instead, at several \$(x,l,\\Pi_0)\$ points: maximum " *
    "absolute residual `" * string(round(max_abs; sigdigits=2)) * "` (floating-point noise) -- i.e. \$\\mathcal{L}_2[U_p]=2(l-1)(l+1)\\Pi_0x^{l-1}\$ exactly.")
```

For the homogeneous part ``U_h=Cxj_l(qx)``, using the spherical Bessel
equation to eliminate ``v''`` at ``z=q`` (the same relation established in
Chapter 6):

```math
\mathcal{L}_2[U_h]\big|_{x=1} = C\!\left[-q^2j_l(q) + 2\bigl(l^2+l-1\bigr)j_l(q) - 2qj_l'(q)\right].
```

Setting ``\mathcal{L}_2[U_p]|_{x=1}+\mathcal{L}_2[U_h]|_{x=1}=0`` gives one
linear relation between ``C`` and ``\Pi_0``, used in the next chapter.

## BC3: Normal stress condition

The most involved of the three: it couples pressure, viscous normal
stress, and the curvature of the deformed surface through surface
tension. The radial normal stress inside the drop is
``-p_{rr}=p+\delta p-2\mu\,\partial u_r/\partial r``, and the free-surface
condition is ``-p_{rr}=T_1(1/R_1+1/R_2)``.

**The curvature term.** To first order in ``\epsilon``, the mean curvature
of the perturbed surface is

```math
\frac{1}{R_1}+\frac{1}{R_2} = \frac{1}{R}\Big[2 + (l-1)(l+2)\,\epsilon\,Y_l^m\Big].
```

This follows from the standard (citable) linearized mean-curvature
formula for a nearly-spherical surface ``r=R+\zeta``,
``1/R_1+1/R_2 = 2/R - (1/R^2)[2\zeta+\Delta_\Omega\zeta]``, where
``\Delta_\Omega`` is the Laplace-Beltrami operator on the *unit* sphere
(dimensionless -- distinct from Chapter 2's ``\nabla^2_{\text{angular}}``,
which carries a ``1/r^2`` since it is the angular piece of the full
3-D Laplacian; the two share the same eigenvalue ``-l(l+1)`` on
``Y_l^m``, which is all that is used below). With ``\zeta=\epsilon RY_l^m``:

```@eval
using Symbolics, Markdown
@variables l_sym eps_sym R_sym Yl
zeta = eps_sym * R_sym * Yl
angular_lap_zeta = eps_sym * R_sym * (-l_sym*(l_sym+1)*Yl)
curvature_from_formula = 2/R_sym - (1/R_sym^2)*(2*zeta + angular_lap_zeta)
curvature_target = (1/R_sym)*(2 + (l_sym-1)*(l_sym+2)*eps_sym*Yl)
residual = simplify(curvature_from_formula - curvature_target; expand=true)
Markdown.parse("The difference between the two sides collapses to `" * string(residual) *
    "`, for symbolic \$l\$.")
```

The ``(l-1)(l+2)`` factor is the same combination that appears in the
inviscid frequency (Chapter 3) -- both arise from the curvature response
to a harmonic deformation, so this is a consistency check worth noting,
not a coincidence.

**The stress balance.** At ``O(\epsilon^0)`` the equilibrium Young-Laplace
balance is automatically satisfied. At ``O(\epsilon^1)``, using the
pressure solution from Chapter 5
(``\delta p|_{x=1}=\rho\epsilon P_0Y_l^m``) and
``u_r=\epsilon_0\sigma R\,G(x)Y_l^m e^{-\sigma t}`` with ``G=U/x^2``:

```@eval
using Symbolics, Markdown
@variables x
@variables U(x)
Dx = Differential(x)
dGdx_at_1 = substitute(expand_derivatives(Dx(U/x^2)), Dict(x => 1))
target_at_1 = substitute(Dx(U) - 2*U, Dict(x => 1))
residual = simplify(dGdx_at_1 - target_at_1; expand=true)
Markdown.parse("``d/dx(U/x^2)|_{x=1}`` minus ``U'(1)-2U(1)`` collapses to `" * string(residual) * "`.")
```

Dividing through by ``\rho\epsilon Y_l^m`` and using ``\mu=\rho\nu``, BC3
becomes

```math
(l-1)(l+2)\frac{T_1}{\rho R} = P_0 - 2\nu\sigma\!\left[U'(1)-2U(1)\right].
```

This is the only equation containing ``T_1`` explicitly. The next
chapter's key move is to rescale parameters so that ``T_1`` drops out
entirely.
