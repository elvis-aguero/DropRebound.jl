# DropRebound.jl

## The problem

A liquid drop of radius ``R``, density ``\rho``, surface tension ``\sigma`` and dynamic
viscosity ``\eta`` falls onto a rigid, non-wetting surface at speed ``v_0``. It flattens
against the surface, spreads, recoils, and in most cases leaves again. Two numbers summarise
the encounter: how long the drop stays in contact, and what fraction of its downward speed it
recovers.

This package predicts both, for Newtonian, shear-thinning and viscoelastic liquids, by solving
the linearised equations of motion in a spectral basis. Nothing is fitted to impact data. The
fluid is characterised by rheometry, the geometry by three dimensionless groups, and the
outcome follows.

![A 3000 ppm shear-thinning drop: ten simulations against 72 experiments](assets/shear_thinning_overlay.png)

## Governing equations

Inside the drop the flow obeys the incompressible Navier–Stokes equations. For the small
deformations this package treats, the nonlinear advection term is negligible against the
unsteady term, and what remains is

```math
\rho\,\partial_t \bm u \;=\; -\nabla p + \eta\,\nabla^2\bm u , \qquad \nabla\cdot\bm u = 0 .
```

The drop is bounded by a free surface at ``r = R + \zeta(\theta,t)``, and the conditions there
carry the physics:

```math
\partial_t\zeta = u_r , \qquad
\bm n\cdot\bm\Sigma\cdot\bm t = 0 , \qquad
\bm n\cdot\bm\Sigma\cdot\bm n = \sigma\,\kappa ,
```

with ``\kappa`` the mean curvature, ``\bm n`` and ``\bm t`` the normal and tangent to the
surface, and

```math
\bm\Sigma = -p\,\bm I + 2\eta\,\bm e , \qquad
\bm e = \tfrac12\left(\nabla\bm u + \nabla\bm u^{\mathsf T}\right)
```

the stress tensor and the strain-rate tensor. Everything below is written in terms of ``\bm e``,
so it is worth noting now that it is symmetric and, for an incompressible flow, traceless.

The first condition is kinematic: the surface moves with the fluid. The second says a free
surface cannot support shear. The third balances the normal traction against capillarity.

The tangential condition is the one that shapes everything downstream. A potential flow cannot
satisfy it, so a viscous drop must carry vorticity near its surface, and any method that
forbids that vorticity will over-predict the damping.

While the drop is close to the substrate, a thin film of air separates the two. Its pressure
``p_c(\theta,t)`` acts on the lower surface and is the only agent by which the wall influences
the drop. The drop never touches the solid.

## Dimensionless groups

Lengths are scaled by ``R`` and time by the capillary time
``\tau_\sigma = \sqrt{\rho R^3/\sigma}``, the natural oscillation period of a free drop. Three
groups remain:

| group | definition | what it controls |
|---|---|---|
| Weber ``\mathrm{We}`` | ``\rho v_0^2 R/\sigma`` | impact energy against surface tension |
| Ohnesorge ``\mathrm{Oh}`` | ``\eta/\sqrt{\rho\sigma R}`` | viscous dissipation against surface tension |
| Bond ``\mathrm{Bo}`` | ``\rho g R^2/\sigma`` | weight against surface tension |

A water–glycerol drop of radius 0.2 mm falling at 9 cm/s has ``\mathrm{Oh} \approx 0.30``,
``\mathrm{Bo} \approx 0.019`` and ``\mathrm{We} \approx 0.079``.

## The one physical approximation

Surface deformations are assumed small against the radius, so the equations above are
linearised about a sphere. The approximation is controlled by ``\mathrm{We}``, and the
amplitude it produces can be checked after the fact: at ``\mathrm{We} = 1`` the largest surface
mode reaches ``|\zeta_2| \approx 0.4R``. That is not small, and results at higher Weber number
should be read with that in mind.

No further approximation enters. Everything after this point is a discretisation whose error
can be reduced by refinement, and *Resolution and Convergence* reports how far it has been
reduced.

## From the equations to a variational statement

Solving the system above directly requires the stress divergence, the pressure field, and
explicit imposition of two stress conditions on a moving boundary. There is a formulation that
avoids all three.

Take a test velocity field ``\bm v`` that is incompressible and compatible with the kinematic
condition. Contract the momentum equation with it and integrate over the drop:

```math
\int \rho\,\partial_t\bm u\cdot\bm v \,dV
\;=\; \int \left(-\nabla p + \eta\nabla^2\bm u\right)\cdot\bm v \,dV .
```

Both terms on the right integrate by parts once. For the pressure,

```math
\int -\nabla p\cdot\bm v\,dV \;=\; -\oint p\,(\bm v\cdot\bm n)\,dS + \int p\,\nabla\cdot\bm v\,dV ,
```

and the volume term vanishes because ``\bm v`` is divergence-free. For the viscous term, use
``\nabla^2\bm u = 2\nabla\cdot\bm e`` (which holds when ``\nabla\cdot\bm u = 0``) and the
symmetry of ``\bm e``, which turns ``\bm e\!:\!\nabla\bm v`` into ``\bm e\!:\!\bm e(\bm v)``.
The two surface terms combine into the traction, and what remains is

```math
\int \rho\,\partial_t\bm u\cdot\bm v\,dV
\;+\; \int 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm v)\,dV
\;=\; \oint \left(\bm\Sigma\cdot\bm n\right)\cdot\bm v \,dS .
```

The pressure has not disappeared. It has left the volume integral, where it was a constraint,
and survives inside the traction, where the normal stress condition will remove it.

The two stress conditions now sit on the right-hand side rather than acting as constraints. The
tangential condition makes the shear part of the traction vanish. The normal condition replaces
what is left by the capillary traction ``\sigma\kappa\,\bm n``, and that traction is the first
variation of the surface energy

```math
V \;=\; \sigma\left(\text{area of the deformed surface} - 4\pi R^2\right) ,
```

so the whole right-hand side is a derivative of an energy. The free-surface conditions have
become **natural**: any solution of the weak form satisfies them, and neither is ever imposed.

### From a field to coordinates

Nothing so far has reduced the problem, because ``\bm u`` is still a field. Introduce
generalised coordinates ``\bm\xi = (\xi_1,\ldots,\xi_N)`` describing the fluid's displacement
from the sphere, so that the velocity is linear in their rates,

```math
\bm u(\bm r,t) \;=\; \sum_{a} \dot\xi_a(t)\,\bm u^{(a)}(\bm r) , \qquad
\bm u^{(a)} = \frac{\partial\bm u}{\partial\dot\xi_a} ,
```

with each ``\bm u^{(a)}`` divergence-free and fixed in time. *Variational Mechanics* constructs
them; here only their existence matters.

The weak form holds for every admissible ``\bm v``, so it holds in particular for each
``\bm v = \bm u^{(a)}`` in turn. That choice is what converts it into equations of motion.
Define

```math
T = \tfrac12\int\rho\,|\bm u|^2\,dV , \qquad
\mathcal R = \eta\int \bm e\!:\!\bm e \,dV ,
```

the kinetic energy and the Rayleigh dissipation function, both quadratic in ``\dot{\bm\xi}``.
Then term by term,

```math
\int\rho\,\partial_t\bm u\cdot\bm u^{(a)}dV = \frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a} ,
\qquad
\int 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm u^{(a)})dV = \frac{\partial\mathcal R}{\partial\dot\xi_a} ,
\qquad
\oint(\bm\Sigma\cdot\bm n)\cdot\bm u^{(a)}dS = -\frac{\partial V}{\partial\xi_a} + Q_a .
```

The first holds because ``\bm u^{(a)}`` does not depend on time, the second because
``\partial\bm e/\partial\dot\xi_a = \bm e(\bm u^{(a)})`` and ``\mathcal R`` is quadratic, and
the third is the statement that the capillary part of the traction is the variation of ``V``
while ``Q_a`` collects the work done by the film pressure. Assembling them,

```math
\frac{d}{dt}\frac{\partial T}{\partial\dot{\bm\xi}}
\;+\; \frac{\partial\mathcal R}{\partial\dot{\bm\xi}}
\;+\; \frac{\partial V}{\partial\bm\xi}
\;=\; \bm Q ,
```

which is Lagrange's equations for a dissipative system. Because ``T``, ``\mathcal R`` and ``V``
are quadratic, their second derivatives are constant matrices, and this reads

```math
\bm M\ddot{\bm\xi} + \bm C\dot{\bm\xi} + \bm G\bm\xi = \bm Q ,
```

with ``\bm M`` the mass matrix from ``T``, ``\bm C`` the damping matrix from ``\mathcal R``,
and ``\bm G`` the stiffness matrix from ``V``.

Three things are gained. Only one derivative of the velocity is needed, where the strong form
needs two. The pressure never has to be solved for, because incompressibility was built into
the trial fields and the normal stress condition disposed of the rest. And the boundary
conditions are automatic.

This is Onsager's variational principle [^onsager], whose modern statement and use as an
approximation tool are set out by Wang, Qian and Xu [^wqx]. The principle is usually written
for overdamped motion with the inertial term dropped, which is not our case: a bouncing drop is
an oscillator, and ``T`` is as important as ``\mathcal R``. The form above is the inertial
generalisation, discussed by Archer [^archer]. Solving such a principle with trial functions is
the Ritz method, applied to variational problems in soft matter by Wang and co-workers [^ritz].

## Getting started

```julia
using DropSolver

p = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 45, K = 3)
r = simulate(p)

r.cor          # coefficient of restitution
r.tc           # contact time, in units of sqrt(rho R^3 / sigma)
```

``M`` truncates the angular expansion and ``K`` the radial one. *Resolution and Convergence*
gives the values at which each quantity stops moving.

## Validation

Against the impact experiments of Gabbard and co-workers, restitution agrees to 8 per cent and
contact time to 13 per cent in median, across 935 measurements spanning
``\mathrm{Oh} \in [0.014, 0.79]``.

Against a 3000 ppm shear-thinning solution, restitution agrees to 6 per cent and contact time
to 8 per cent across seven Weber groups of at least five repetitions each. A Newtonian drop at
the same zero-shear viscosity does not rebound at all, so the measured rebound exists because
the fluid thins. The Carreau–Yasuda parameters come from that fluid's own rheometry.

## How to read this

The material is ordered so that each part supplies what the next one needs.

**Variational mechanics** turns the statement above into matrices: how a trial space is chosen so
that incompressibility costs nothing, and what the three operators inherit from the functionals
they differentiate.

**The free viscous drop** solves the linearised problem exactly, by separation of variables.
It supplies the target the discretisation is checked against, and its radial structure supplies
the trial functions.

**Contact** treats the air film and the unilateral constraint between gap and pressure, and the
two ways the contact region is located.

**Shear-thinning fluids** carries the viscosity through two routes: evaluated pointwise on the
full strain field, which keeps the interior in the state, and reduced to one number per mode,
which does not.

**Viscoelastic fluids** adds polymer stress as an extra state variable.

**Using it** covers solver choice, resolution, and the interface.

[^onsager]: L. Onsager, *Reciprocal relations in irreversible processes I*, Phys. Rev. **37**, 405 (1931).
[^wqx]: H. Wang, T. Qian and X. Xu, *Onsager's variational principle in active soft matter*, [arXiv:2011.10821](https://arxiv.org/abs/2011.10821).
[^archer]: A. J. Archer, *Variational formulation for the dynamics of soft matter including inertia*, [arXiv:2607.18457](https://arxiv.org/abs/2607.18457).
[^ritz]: H. Wang, B. Zou, J. Su, D. Wang and X. Xu, *Variational methods and deep Ritz method for active elastic solids*, [arXiv:2203.15376](https://arxiv.org/abs/2203.15376).
