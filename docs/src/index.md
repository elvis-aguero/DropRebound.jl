# DropRebound.jl

## The problem

A liquid drop of radius ``R``, density ``\rho``, surface tension ``\gamma`` and dynamic
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

Write ``\Omega`` for the region the liquid occupies and ``\partial\Omega`` for its free
surface, at ``r = R + \zeta(\theta,t)``. Every integral below is over one of these two. The
conditions on ``\partial\Omega`` carry the physics:

```math
\partial_t\zeta = u_r , \qquad
\bm n\cdot\bm\Sigma\cdot\bm t = 0 , \qquad
\bm n\cdot\bm\Sigma\cdot\bm n = \gamma\,\kappa - p_c ,
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
surface cannot support shear. The third balances the normal traction against capillarity and
against ``p_c``, the pressure in the thin air film that separates the drop from the substrate
during an impact. That film is the only agent by which the wall acts on the drop, and the drop
never touches the solid. Away from the substrate ``p_c = 0`` and the third condition is the
usual Young-Laplace balance.

The tangential condition is the one that shapes everything downstream. A potential flow cannot
satisfy it, so a viscous drop must carry vorticity near its surface, and any method that
forbids that vorticity will over-predict the damping.

## Dimensionless groups

Lengths are scaled by ``R`` and time by the capillary time
``\tau_c = \sqrt{\rho R^3/\gamma}``, the natural oscillation period of a free drop. Three
groups remain:

| group | definition | what it controls |
|---|---|---|
| Weber ``\mathrm{We}`` | ``\rho v_0^2 R/\gamma`` | impact energy against surface tension |
| Ohnesorge ``\mathrm{Oh}`` | ``\eta/\sqrt{\rho\gamma R}`` | viscous dissipation against surface tension |
| Bond ``\mathrm{Bo}`` | ``\rho g R^2/\gamma`` | weight against surface tension |

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
\int_\Omega \rho\,\partial_t\bm u\cdot\bm v \,dV
\;=\; \int_\Omega \left(-\nabla p + \eta\nabla^2\bm u\right)\cdot\bm v \,dV .
```

Both terms on the right integrate by parts once. For the pressure,

```math
\int_\Omega -\nabla p\cdot\bm v\,dV
\;=\; -\oint_{\partial\Omega} p\,(\bm v\cdot\bm n)\,dS + \int_\Omega p\,\nabla\cdot\bm v\,dV ,
```

and the volume term vanishes because ``\bm v`` is divergence-free. For the viscous term, use
``\nabla^2\bm u = 2\nabla\cdot\bm e`` (which holds when ``\nabla\cdot\bm u = 0``) and the
symmetry of ``\bm e``, which turns ``\bm e\!:\!\nabla\bm v`` into ``\bm e\!:\!\bm e(\bm v)``.
The two surface terms combine into the traction, and what remains is

```math
\int_\Omega \rho\,\partial_t\bm u\cdot\bm v\,dV
\;+\; \int_\Omega 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm v)\,dV
\;=\; \oint_{\partial\Omega} \left(\bm\Sigma\cdot\bm n\right)\cdot\bm v \,dS .
```

The pressure has not disappeared. It has left the volume integral, where it was a constraint,
The pressure has not disappeared. It has left the volume integral, where it was a constraint,
and survives inside the traction, where the stress conditions will now remove it.

Split the test field on the surface into its normal and tangential parts,
``\bm v = (\bm v\cdot\bm n)\,\bm n + \bm v_t``. The traction then contributes two pieces,

```math
\left(\bm\Sigma\cdot\bm n\right)\cdot\bm v
\;=\; \underbrace{\left(\bm n\cdot\bm\Sigma\cdot\bm n\right)}_{\gamma\kappa - p_c}(\bm v\cdot\bm n)
\;+\; \underbrace{\left(\bm n\cdot\bm\Sigma\cdot\bm t\right)}_{=\,0}\,(\bm v\cdot\bm t) ,
```

and each is settled by one of the two stress conditions. The tangential condition kills the
second outright. The normal condition turns the first into capillarity and film pressure, so
the right-hand side becomes

```math
\oint_{\partial\Omega}\left(\bm\Sigma\cdot\bm n\right)\cdot\bm v\,dS
\;=\; \underbrace{\oint_{\partial\Omega}\gamma\kappa\,(\bm v\cdot\bm n)\,dS}_{\text{capillary}}
\;-\; \underbrace{\oint_{\partial\Omega} p_c\,(\bm v\cdot\bm n)\,dS}_{\text{film}} .
```

The capillary term is a derivative of an energy. Deforming the surface by a normal displacement
``\delta\zeta`` changes its area by ``\oint_{\partial\Omega} \kappa\,\delta\zeta\,dS``, so with

```math
V \;=\; \gamma\left(|\partial\Omega| - 4\pi R^2\right)
```

the excess surface energy, the first variation is ``\delta V = \oint_{\partial\Omega}
\gamma\kappa\,\delta\zeta\,dS``. Since ``\bm v\cdot\bm n`` is a rate of normal displacement,
the capillary integral is exactly ``V`` differentiated along ``\bm v``.

Collecting, the weak form reads

```math
\int_\Omega \rho\,\partial_t\bm u\cdot\bm v\,dV
\;+\; \int_\Omega 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm v)\,dV
\;+\; \delta_{\bm v} V
\;=\; -\oint_{\partial\Omega} p_c\,(\bm v\cdot\bm n)\,dS ,
```

for every admissible ``\bm v``: inertia, dissipation and capillarity on the left, and the wall
on the right. The free-surface conditions have become **natural**. They were used once, to
evaluate the boundary term, and they are never imposed on the solution; any field satisfying
this statement for all ``\bm v`` satisfies them.

Three things are gained. Only one derivative of the velocity appears, where the strong form
needs two. The pressure never has to be solved for, because incompressibility was built into
the test fields and the normal stress condition disposed of the rest. And the boundary
conditions are automatic.

What remains is to turn this from a statement about fields into a finite set of equations.
*Variational Mechanics* does that, and shows that the result is Lagrange's equations for a
damped system.

This is Onsager's variational principle [^onsager], whose modern statement and use as an
approximation tool are set out by Wang, Qian and Xu [^wqx]. The principle is usually written
for overdamped motion with the inertial term dropped, which is not our case: a bouncing drop is
an oscillator, and the first term above matters as much as the second. The inertial
generalisation is discussed by Archer [^archer]. Solving such a principle with trial functions
is the Ritz method, applied to variational problems in soft matter by Wang and co-workers
[^ritz].

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
