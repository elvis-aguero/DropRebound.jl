# Problem Setup

## Equilibrium state

A liquid drop of density ``\rho``, viscosity ``\mu = \rho\nu``, and surface
tension ``T_1`` is in equilibrium as a sphere of radius ``R``. With the
external pressure set to zero, the Young-Laplace equation gives the
internal equilibrium pressure:

```math
p_0 = \frac{2T_1}{R}.
```

!!! note "Young-Laplace equation"
    For an interface with principal radii of curvature ``R_1`` and
    ``R_2``, the pressure jump across the interface is
    ``\Delta p = T_1(1/R_1 + 1/R_2)``. For a sphere, ``R_1=R_2=R``, giving
    ``\Delta p = 2T_1/R``.

## Perturbed surface

We consider a small deformation of the spherical surface of the form

```math
r = R\left[1 + \epsilon\,Y_l^m(\theta,\varphi)\right], \qquad \epsilon \ll 1,
```

with time dependence

```math
\epsilon(t) = \epsilon_0\,e^{-\sigma t},
```

where ``\sigma`` is a complex number to be determined. The convention is
chosen so that ``\mathrm{Re}(\sigma) > 0`` corresponds to decay.

!!! tip "Why $e^{-\sigma t}$ rather than $e^{+i\omega t}$?"
    This follows Chandrasekhar's convention, where ``\sigma`` is complex
    and a damped oscillation corresponds to ``\sigma = \gamma + i\omega_d``
    with ``\gamma>0`` (decay) and ``\omega_d>0`` (oscillation frequency).
    In the inviscid limit, ``\gamma\to 0`` and ``\sigma\to\pm i\omega_d``.

## Inviscid frequencies

In the total absence of viscosity, Rayleigh and Lamb's classical result
gives the oscillation frequency of a spherical-harmonic mode of order
``l`` as

```math
\sigma_{l;0}^2 = l(l-1)(l+2)\,\frac{T_1}{\rho R^3}.
```

This is real (undamped) and independent of the azimuthal order ``m``. The
``l=2`` mode, for example, is the oblate-prolate oscillation familiar from
a wobbling drop of water. This formula requires a separate calculation (the
inviscid energy functional) that we do not re-derive here -- we take it as
given, exactly as Reid does, and use it only as a bookkeeping device:
viscosity enters everything that follows only through this one number.

Two consequences of the formula are checked directly rather than taken on
faith. It should be dimensionally a rate squared:
``[T_1/(\rho R^3)] = [\mathrm{N/m}]/[\mathrm{kg/m^3\cdot m^3}]
= [\mathrm{kg/s^2}]/[\mathrm{kg}] = [\mathrm{s^{-2}}]`` -- consistent.
And the ``l=1`` mode, a rigid translation of the whole drop, has no
restoring force (translating a sphere changes neither its shape nor its
surface energy), so ``\sigma_{1;0}^2`` must vanish identically:

```@eval
using Symbolics, Markdown
@variables l_sym T1 rho R
sigma_l0_sq = l_sym * (l_sym - 1) * (l_sym + 2) * T1 / (rho * R^3)
at_l1 = simplify(substitute(sigma_l0_sq, Dict(l_sym => 1)); expand=true)
Markdown.parse("```math\n\\sigma_{l;0}^2\\Big|_{l=1} = " * Main.pretty_latex(at_l1) *
    "\n```")
```
