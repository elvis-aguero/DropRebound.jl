# Linearized Governing Equations

The departures from equilibrium are governed by the linearized
incompressible Navier-Stokes equations,

```math
\frac{\partial\bm u}{\partial t} = -\nabla\frac{\delta p}{\rho}
                                    - \nu\,\nabla\times\nabla\times\bm u,
\qquad
\nabla\cdot\bm u = 0.
```

!!! note "An equivalent, more familiar form"
    For incompressible flow, ``\nabla\cdot\bm u = 0`` implies
    ``\nabla^2\bm u = -\nabla\times\nabla\times\bm u``, so this is
    equivalent to ``\partial_t\bm u = -\nabla(\delta p/\rho) + \nu\nabla^2\bm u``.
    Both forms appear in the literature; we use the second from here on.

With the time dependence ``\epsilon(t)=\epsilon_0 e^{-\sigma t}``, all
fields vary as ``e^{-\sigma t}``, so ``\partial_t \to -\sigma``. The
momentum equation becomes

```math
-\sigma\bm u = -\nabla\frac{\delta p}{\rho} + \nu\nabla^2\bm u.
```

This is the equation everything downstream is built from: the pressure
field (next chapter) satisfies Laplace's equation as a direct consequence
of its divergence, and the velocity field satisfies a single scalar ODE
(two chapters from now) as a consequence of its radial component.
