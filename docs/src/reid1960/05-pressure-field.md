# The Pressure Field

## Pressure satisfies Laplace's equation

Take the divergence of the momentum equation from the previous chapter:

```math
-\sigma\underbrace{(\nabla\cdot\bm u)}_{=\,0}
  = -\nabla^2\!\left(\frac{\delta p}{\rho}\right)
  + \nu\underbrace{\nabla\cdot(\nabla^2\bm u)}_{\nabla^2(\nabla\cdot\bm u)\,=\,0},
```

where ``\nabla\cdot\nabla^2\bm u = \nabla^2(\nabla\cdot\bm u) = 0`` because
divergence and Laplacian commute. Both marked terms vanish because
``\nabla\cdot\bm u = 0``, leaving

```math
\nabla^2\!\left(\frac{\delta p}{\rho}\right) = 0.
```

The pressure perturbation satisfies Laplace's equation. We need the
solution that (a) has angular dependence ``Y_l^m``, and (b) is regular at
the origin ``r=0``. Writing ``\delta p/\rho = f(r)\,Y_l^m``, the scalar
Laplacian identity from Chapter 2 turns this into the radial ODE

```math
f'' + \frac{2}{r}f' - \frac{l(l+1)}{r^2}f = 0,
```

whose general solution is ``f(r) = Ar^l + Br^{-(l+1)}``:

```@eval
using Symbolics, Markdown
@variables x l A B
Dx = Differential(x)
f = A * x^l + B * x^(-(l + 1))
residual = simplify(expand_derivatives(Dx(Dx(f)) + 2 / x * Dx(f) - l * (l + 1) / x^2 * f); expand=true)
Markdown.parse("Substituting ``f = Ax^l + Bx^{-(l+1)}`` into the ODE, for symbolic ``l``, the left-hand side collapses to `" *
    string(residual) * "`.")
```

Regularity at ``r=0`` requires ``B=0``. Writing ``x=r/R``, we therefore
have

```math
\frac{\delta p}{\rho} = \epsilon\,P_0\,x^l\,Y_l^m,
```

where ``P_0`` is a constant of integration (dimensions of velocity
squared) fixed later by the boundary conditions.

## Scalar potential for the pressure gradient

Since ``\nabla(\delta p/\rho)`` is a purely poloidal, divergence-free
vector field with angular structure ``Y_l^m``, it can be represented
through a single scalar function ``\Pi(x)`` via its radial component:

```math
\left(\nabla\frac{\delta p}{\rho}\right)_r
  = \epsilon_0\,\sigma^2 R\,\frac{\Pi(x)}{x^2}\,Y_l^m\,e^{-\sigma t}.
```

Computing the radial component of ``\nabla(\epsilon P_0 x^l Y_l^m)``
directly,

```math
\left(\nabla\frac{\delta p}{\rho}\right)_r
  = \frac{\partial}{\partial r}\!\left(\epsilon P_0 x^l Y_l^m\right)
  = \epsilon_0 e^{-\sigma t} P_0\,\frac{l}{R}\,x^{l-1}\,Y_l^m,
```

and comparing the two gives ``\Pi(x) = \Pi_0\,x^{l+1}`` with

```math
\Pi_0 = \frac{l}{\sigma^2 R^2}\,P_0.
```

This is the second of the three constants (``P_0``, or equivalently
``\Pi_0``) that the boundary conditions will pin down; the third, ``C``,
appears once the velocity field enters in the next chapter.
