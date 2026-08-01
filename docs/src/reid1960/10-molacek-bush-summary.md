# Connection to Molaček & Bush (2012), and Summary

## Connection to Molaček & Bush (2012)

In their quasi-static model of drop impact, Molaček & Bush parameterize
the kinetic energy and viscous dissipation associated with each surface
harmonic mode ``m`` (their notation) via two coefficients
``A_m(\mathrm{Oh}_m)`` and ``D_m(\mathrm{Oh}_m)``, which appear in their
Eq. 31:

```math
\mathrm{K.E.} = \pi\rho R_0^5 \dot{B}^2 \sum_m A_m \frac{2b_m^2}{m(2m+1)},
\qquad
D = 8\pi\mu R_0^3 \dot{B}^2 \sum_m D_m \frac{m}{2m+1}b_m^2.
```

These coefficients are defined so that the two roots of the quadratic

```math
A_m b^2 - 2a D_m b + 1 = 0, \qquad a \equiv \mathrm{Oh}\sqrt{m(m-1)(m+2)},
```

reproduce the two dominant roots of Reid's characteristic equation for
harmonic order ``l=m``. The variable identification connecting this
quadratic to Reid's transcendental equation is
``q^2 = (b/a)\cdot\alpha^2`` (equivalently ``b = a\cdot q^2/\alpha^2``),
under which ``W(b/a) = Q_{l+1/2}(q)/q``. In other words, ``A_m`` and
``D_m`` encode the *result* of solving Reid's problem at each order
``m``, compressed into two numbers per mode via a quadratic fit to the
dominant eigenvalue pair -- this section describes a mapping between two
representations of the same physics, not a separate derivation.

**Workflow.**
1. For each mode ``m`` and Ohnesorge number ``\mathrm{Oh}``, evaluate
   ``\alpha^2=\mathrm{Oh}^{-1}\sqrt{m(m-1)(m+2)}`` and ``q^2=\sigma R^2/\nu``.
2. Numerically solve the characteristic equation for the two roots
   ``b_{1,2}`` with the smallest ``\mathrm{Re}(\sigma)``.
3. Read off ``A_m`` and ``D_m`` by matching the quadratic: sum of roots
   ``=2aD_m/A_m``; product of roots ``=1/A_m``.
4. Tabulate or fit ``A_m(\mathrm{Oh}_m)`` and ``D_m(\mathrm{Oh}_m)`` for
   use in the Lagrangian equation of motion.

**Why discard the higher roots?** The quasi-static assumption restricts
the drop shape to evolve through a one-parameter family of equilibrium
shapes. Within that family, only the fundamental surface mode at each
harmonic order ``l`` is representable -- higher radial overtones have
internal nodal surfaces absent from the quasi-static shape family, and
they decay on timescales much shorter than ``R^2/\nu``, well separated
from the impact timescale.

## Summary of key equations

| Equation | Physical content |
|:--|:--|
| ``\sigma_{l;0}^2 = l(l-1)(l+2)\,T_1/(\rho R^3)`` | Inviscid capillary frequency |
| ``q^2 = \sigma R^2/\nu`` | Complex decay rate, scaled |
| ``\alpha^2 = \sigma_{l;0}R^2/\nu`` | Inviscid frequency, scaled |
| ``Q_{l+1/2}(q) = J_{l+3/2}(q)/J_{l+1/2}(q)`` | Bessel function ratio |
| ``\alpha^4/q^4 + 1 = (2(l-1)/q^2)[\cdots]`` | Reid's characteristic equation |
| ``\sigma = (l-1)(2l+1)\nu/R^2 \pm i\,\sigma_{l;0}`` | Small-``\nu`` (Lamb) limit |

This closes the Reid (1960) derivation. The poloidal decomposition, the
pressure-field derivation, and BC1 (kinematic) depend only on
incompressibility and geometry -- they are rheology-agnostic and carry
over unchanged to every other model in this repo. The Newtonian
constitutive assumption flagged in Chapter 6 enters in more than the one
place that chapter names, though: ``\nu\nabla^2\bm u`` in the momentum
equation is one instance of it, but BC2 and BC3 both also assume a
constant ``\mu`` multiplying a linear strain rate (``\tau_{r\theta}=\mu[\ldots]``,
``-2\mu\,\partial u_r/\partial r``). A shear-thinning model changes the
momentum equation AND both stress boundary conditions, not just the
former -- Oldroyd-B and Carreau-Yasuda each have to redo that work, not
merely swap a coefficient.
