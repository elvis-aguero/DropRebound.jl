# The Velocity Field and its Governing ODE

This is the technical heart of the whole problem: one ODE that packages up
the entire viscous, incompressible flow field consistent with the
linearized momentum equation.

## Poloidal representation

Since the pressure gradient is purely poloidal, the momentum equation
forces the velocity field to be purely poloidal too -- nothing drives a
toroidal part. We write the radial component as

```math
u_r = \epsilon_0\,\sigma R\,\frac{U(x)}{x^2}\,Y_l^m\,e^{-\sigma t},
```

where ``U(x)`` is a dimensionless scalar function. The ``1/x^2`` scaling
is conventional; it is exactly what makes the ODE for ``U`` below take the
clean, Bessel-equation form from Chapter 2 -- that is the entire point of
choosing it.

!!! note "Where the Newtonian assumption enters"
    The momentum equation this chapter's ODE comes from,
    ``-\sigma\bm u = -\nabla(\delta p/\rho) + \nu\nabla^2\bm u``, carries a
    single viscous term, ``\nu\nabla^2\bm u`` -- a scalar, constant ``\nu``
    multiplying the Laplacian of the velocity. That is the Newtonian
    constitutive law: stress proportional to strain rate with a
    viscosity that does not depend on the flow itself. A shear-thinning
    fluid (Carreau-Yasuda) replaces that scalar with a viscosity that
    depends on the local strain rate, which in general no longer factors
    out of the Laplacian this cleanly. Everything in this chapter and the
    next -- the ODE for ``U``, both boundary conditions, and the
    characteristic equation -- depends on this one substitution and does
    not carry over unmodified.

## Deriving the ODE for $U(x)$

Writing ``G(x) = U(x)/x^2``, so that ``u_r \propto G(x)``, and using the
scalar-Laplacian eigenvalue identity from Chapter 2, the ``r``-component of
the vector Laplacian for an axisymmetric poloidal field is

```math
\left[\nabla^2\bm u\right]_r = \nabla^2 u_r + \frac{2}{r^2}u_r + \frac{2}{r}\frac{\partial u_r}{\partial r}.
```

!!! note "A fact we cite rather than re-derive"
    Reaching this formula from the raw vector-Laplacian components requires
    the incompressibility constraint to eliminate the angular (``u_\theta``)
    part in favor of radial derivatives of ``u_r`` -- a genuine
    vector-calculus derivation (the poloidal stream-function
    representation), but a different one from the ODE algebra that
    follows, and not where a silent sign error would actually hide. We
    take it as given here, at the same level the source material does.
    Everything from here on **is** verified directly.

Substituting ``u_r = \sigma R\,G(x)\,Y_l^m`` (the ``\epsilon_0 e^{-\sigma t}``
prefactor cancels throughout), the three pieces combine into a single
operator on ``G``:

```@eval
using Symbolics, Markdown
@variables x l Gfun(..)
G = Gfun(x)
Dx = Differential(x)
term_scalar_lap = Dx(Dx(G)) + 2/x*Dx(G) - l*(l+1)/x^2*G
term_2_over_r2 = 2/x^2*G
term_2_over_r_deriv = 2/x*Dx(G)
combined = expand_derivatives(term_scalar_lap + term_2_over_r2 + term_2_over_r_deriv)
target = expand_derivatives(Dx(Dx(G)) + 4/x*Dx(G) + (2 - l*(l+1))/x^2*G)
Markdown.parse("Summing the three pieces, the difference from ``G'' + (4/x)G' + (2-l(l+1))/x^2\\,G`` collapses to `" *
    string(simplify(combined - target; expand=true)) * "`: the operator on the right is exactly the combined vector Laplacian.")
```

The ``r``-component of the momentum equation (after canceling
``Y_l^m e^{-\sigma t}`` and using ``x=r/R``, ``q^2=\sigma R^2/\nu``) is
then

```math
-q^2 G = -\frac{P_0\,l\,x^{l-1}}{\sigma\nu} + G'' + \frac{4}{x}G' + \frac{2-l(l+1)}{x^2}G.
```

Substituting ``G=U/x^2`` -- this is where the ``1/x^2`` scaling earns its
keep --

```@eval
using Symbolics, Markdown
@variables x l Ufun(..)
Dx = Differential(x)
Usym = Ufun(x)
Gsub = Usym / x^2
lhs = expand_derivatives(Dx(Dx(Gsub)) + 4/x*Dx(Gsub) + (2 - l*(l+1))/x^2*Gsub)
rhs = expand_derivatives((Dx(Dx(Usym)) - l*(l+1)/x^2*Usym)/x^2)
Markdown.parse("the difference collapses to `" * string(simplify(lhs - rhs; expand=true)) *
    "`: substituting ``G=U/x^2`` turns the equation above into one purely in ``U``.")
```

Multiplying through by ``x^2`` and using ``P_0 = \sigma^2 R^2\Pi_0/l``
(previous chapter) to write ``P_0 l/(\sigma\nu) = q^2\Pi_0``, the momentum
equation becomes exactly Reid's Eq. 9 -- the single equation that packages
up the entire linearized, viscous, incompressible flow problem:

```@eval
using Symbolics, Markdown
@variables x l q Pi_0
rhs = q^2 * Pi_0 * x^(l + 1)
Markdown.parse("```math\n\\left[\\frac{d^2}{dx^2} - \\frac{l(l+1)}{x^2} + q^2\\right] U(x) = " *
    Main.pretty_latex(rhs) * "\n```")
```

## Solving the ODE

This is a linear, second-order, inhomogeneous ODE. Its general solution is
a particular solution plus the general homogeneous solution.

**Particular solution.** Try ``U_p = \Pi_0\,x^{l+1}`` directly:

```@eval
using Symbolics, Markdown
@variables x l q Pi_0
Dx = Differential(x)
Up = Pi_0 * x^(l + 1)
residual = simplify(expand_derivatives(Dx(Dx(Up)) - l*(l+1)/x^2*Up + q^2*Up - q^2*Pi_0*x^(l+1)); expand=true)
Markdown.parse("Substituting ``U_p`` into Reid's Eq. 9, the left-hand side minus the right-hand side collapses to `" *
    string(residual) * "` -- the ``l(l+1)/x^2`` and ``q^2`` terms of ``U_p''`` exactly cancel the corresponding terms\n" *
    "on the left, leaving the ``q^2\\Pi_0 x^{l+1}`` forcing.")
```

**Homogeneous solution.** The homogeneous equation
``U''-l(l+1)U/x^2+q^2U=0`` is the Chapter 2 Bessel-substitution ODE with
``x\to qx``: writing ``U_h = x\,v(qx)`` and using the spherical Bessel
equation for ``v`` at argument ``z=qx`` collapses the whole thing to zero
identically -- so ``U_h = C\,x\,j_l(qx)`` for any constant ``C``. The other
homogeneous solution, built from the second-kind spherical Bessel function
``n_l``, diverges as ``z\to 0`` and is discarded on the same regularity
grounds as the pressure field's ``B``-term in the previous chapter.

**General solution.**

```@eval
using Symbolics, Markdown
@variables x l q C Pi_0
Markdown.parse("```math\nU(x) = " * Main.pretty_latex(C*x) * "\\,j_l(qx) + " *
    Main.pretty_latex(Pi_0*x^(l+1)) * "\n```")
```

where ``C`` and ``\Pi_0`` are fixed by the three boundary conditions
derived next.
