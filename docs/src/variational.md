# Variational Mechanics of a Dissipative Fluid

The home page ended with a weak form: for every admissible test field ``\bm v``,

```math
\int_\Omega \rho\,\partial_t\bm u\cdot\bm v\,dV
\;+\; \int_\Omega 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm v)\,dV
\;+\; \delta_{\bm v} V
\;=\; -\oint_{\partial\Omega} p_c\,(\bm v\cdot\bm n)\,dS .
```

It is exact, and on its own it computes nothing, because ``\bm u`` is a field and the statement
ranges over an infinite-dimensional space of ``\bm v``.

This page reduces it to a finite system, in four steps. First the general reduction from fields
to coordinates, which needs no assumption about drops. Then the choice of coordinates for this
particular problem. Then a representation that makes the test fields incompressible for free.
Then the quadratures that produce the matrices.

## From fields to coordinates

Suppose the velocity can be written as a finite sum

```math
\bm u(\bm r,t) \;=\; \sum_{a=1}^{N} \dot\xi_a(t)\,\bm u^{(a)}(\bm r) ,
```

where the ``\xi_a`` are generalised coordinates measuring the fluid's displacement from the
sphere, and each ``\bm u^{(a)} = \partial\bm u/\partial\dot\xi_a`` is a fixed, divergence-free
field. The velocity is linear in the rates, so the ``\bm u^{(a)}`` do not depend on time.

The weak form holds for every admissible ``\bm v``, so it holds for each ``\bm u^{(a)}`` in
turn. Making that substitution is the whole reduction: ``N`` choices of test field give ``N``
scalar equations. Take them one at a time.

**Inertia.** With ``\bm u^{(a)}`` independent of time,

```math
\int_\Omega \rho\,\partial_t\bm u\cdot\bm u^{(a)}\,dV
\;=\; \frac{d}{dt}\int_\Omega \rho\,\bm u\cdot\bm u^{(a)}\,dV
\;=\; \frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a} ,
\qquad
T = \tfrac12\int_\Omega\rho\,|\bm u|^2\,dV .
```

**Dissipation.** Since ``\bm e`` is linear in the rates,
``\partial\bm e/\partial\dot\xi_a = \bm e(\bm u^{(a)})``, and so

```math
\int_\Omega 2\eta\,\bm e(\bm u)\!:\!\bm e(\bm u^{(a)})\,dV
\;=\; \frac{\partial\mathcal R}{\partial\dot\xi_a} ,
\qquad
\mathcal R = \eta\int_\Omega \bm e\!:\!\bm e\,dV .
```

``\mathcal R`` is the Rayleigh dissipation function. It is half the rate at which the fluid
turns kinetic energy into heat, which is where the factor of two in the integrand goes.

**Capillarity.** The home page identified the capillary term as the variation of the surface
energy, so ``\delta_{\bm u^{(a)}}V = \partial V/\partial\xi_a``.

**The wall.** The right-hand side defines the generalised force,

```math
Q_a \;=\; -\oint_{\partial\Omega} p_c\,\left(\bm u^{(a)}\!\cdot\bm n\right)\,dS ,
```

the rate at which the film pressure does work on the ``a``-th coordinate. Its modal form is
derived in *Contact*; here it is enough that it is a boundary integral of a known pressure
against a known field, so it can be evaluated once the pressure is known.

Assembling the four,

```math
\frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a}
\;+\; \frac{\partial\mathcal R}{\partial\dot\xi_a}
\;+\; \frac{\partial V}{\partial\xi_a}
\;=\; Q_a ,
```

which is Lagrange's equations for a system with dissipation. Nothing so far is specific to
drops. What remains is to choose the ``\xi_a``.

## Surface modes

The drop is axisymmetric, so its surface is described by one function of polar angle,
``r = 1 + \zeta(\theta,t)``. Expand that function in Legendre polynomials,

```math
\zeta(\theta,t) \;=\; \sum_{l\ge 2} \zeta_l(t)\,P_l(\cos\theta) .
```

Each term is a **surface mode**: a fixed shape ``P_l(\cos\theta)`` whose amplitude ``\zeta_l``
varies in time. Mode ``l = 2`` is the prolate-oblate oscillation, ``l = 3`` is pear-shaped, and
the shapes get finer as ``l`` grows.

Legendre polynomials are the right basis for three reasons. They are the axisymmetric spherical
harmonics, so they are regular at both poles. They are orthogonal, which is what makes the
Newtonian problem separate. And they diagonalise the surface energy, giving

```math
V \;=\; \frac{\sigma}{2}\sum_{l\ge2}\frac{4\pi}{2l+1}\,(l-1)(l+2)\,\zeta_l^2
```

for the excess area to second order in the amplitudes. The sum starts at ``l = 2`` because
``l = 0`` changes the volume, which is forbidden, and ``l = 1`` translates the drop without
deforming it. Both give zero restoring force, as the factor ``(l-1)`` records.

## Coordinates that are displacements, not velocities

A Lagrangian needs positions and velocities. A flow is naturally described by a velocity field,
which is the wrong half of that pair, so the coordinates cannot be velocities.

Take instead the interior displacement. For each mode ``l``, let ``\chi_l(x,t)`` be the radial
profile of how far the fluid has moved from its rest position, with ``x = r/R``. Three
quantities are then tied together:

```math
\underbrace{\chi_l(x,t)}_{\text{coordinate}} , \qquad
\underbrace{\psi_l = \dot\chi_l}_{\text{velocity}} , \qquad
\underbrace{\zeta_l = \chi_l(1)}_{\text{surface amplitude}} .
```

The middle equality is the definition of a velocity. The right-hand one is the kinematic
condition ``\dot\zeta = u_r`` integrated in time: the surface rides on the fluid, so its
displacement is the boundary value of the interior displacement.

That last point resolves what otherwise looks like a miscount. There is no evolution equation
for ``\zeta_l``, and there should not be, because ``\zeta_l`` is not an independent coordinate.
The interior is the only thing that moves; the surface is its trace at ``x = 1``.

So the coordinates ``\xi_a`` of the previous section are the amplitudes of ``\chi_l``, and the
fields ``\bm u^{(a)}`` are the velocity fields those amplitudes generate. The next section
constructs them.

## Making the test fields incompressible

The reduction required every ``\bm u^{(a)}`` to be divergence-free. Imposing that as a
constraint would reintroduce a pressure as its multiplier, undoing the work of the weak form.
In axisymmetry it can be had for nothing instead.

Any axisymmetric incompressible flow derives from a Stokes stream function ``\Psi``,

```math
\bm u \;=\; \nabla\times\!\left(\frac{\Psi}{r\sin\theta}\,\hat{\bm\varphi}\right) ,
```

with ``\hat{\bm\varphi}`` the azimuthal unit vector. Any field written this way satisfies
``\nabla\cdot\bm u = 0`` identically, for every ``\Psi``, because the divergence of a curl
vanishes. The constraint is discharged by the choice of representation, so nothing in the
discrete problem is a Lagrange multiplier and no saddle-point system appears.

Separate ``\Psi`` mode by mode,

```math
\Psi \;=\; f(x)\,\Gamma_l(\theta) , \qquad
\Gamma_l(\theta) = -\frac{\sin\theta\,\partial_\theta P_l(\cos\theta)}{l(l+1)} ,
```

where ``f`` is the radial profile and ``\Gamma_l`` is the Gegenbauer function, fixed by
demanding that the radial velocity come out proportional to ``P_l(\cos\theta)`` so that it
matches the surface mode it drives. Taking the curl gives

```math
u_r = \frac{f}{x^2}\,P_l , \qquad
u_\theta = \frac{f'}{x\,l(l+1)}\,\partial_\theta P_l ,
```

and differentiating those gives the strain rates,

```math
e_{rr} = \left(\frac{f'}{x^2} - \frac{2f}{x^3}\right)P_l , \qquad
e_{r\theta} = \frac{1}{2x\,l(l+1)}
\underbrace{\left(f'' - \frac{2f'}{x} + \frac{l(l+1)f}{x^2}\right)}_{\textstyle\mathcal T[f]}
\,\partial_\theta P_l ,
```

with ``e_{\theta\theta}`` and ``e_{\varphi\varphi}`` following the same way. Their trace
vanishes identically: Legendre's equation gives
``\partial_\theta^2 P_l + \cot\theta\,\partial_\theta P_l = -l(l+1)P_l``, and the surviving
terms cancel in pairs. Incompressibility is visible in the algebra rather than asserted.

The operator ``\mathcal T`` marked above is the tangential-stress operator. The free-surface
condition ``\bm n\cdot\bm\Sigma\cdot\bm t = 0`` is precisely ``\mathcal T[f] = 0`` at ``x = 1``.
Nothing in the discretisation imposes it, because it emerged as a natural condition of the weak
form. A converged solution satisfies it and an unconverged one approaches it, which is the
mechanism behind the potential-flow failure described in *Resolution and Convergence*: a single
trial function cannot annihilate ``\mathcal T``, and the leftover tangential stress is paid for
as extra dissipation.

## The matrices

One choice is left. Restrict each radial profile to a finite span,

```math
\chi_l(x,t) \;=\; \sum_{k=1}^{K} a_k(t)\,\phi_k(x) ,
```

with trial functions regular at the origin. The coordinates ``\xi_a`` are now the ``a_k``, one
set per mode, and each carries a velocity field ``\bm u^{(a)}`` obtained by putting
``f = \phi_k`` into the formulas above.

From here the problem is nondimensional: lengths in units of ``R``, time in units of
``\tau_\sigma``, and ``\eta`` the local viscosity divided by its zero-shear value, so ``\eta=1``
for a Newtonian fluid. The Ohnesorge number is then the only material parameter left.

Substituting the expansion into the three functionals makes each one an explicit quadratic form,

```math
T = \tfrac12\,\dot{\bm a}^{\mathsf T}\bm M\,\dot{\bm a} , \qquad
\mathcal R = \tfrac12\,\dot{\bm a}^{\mathsf T}\bm C\,\dot{\bm a} , \qquad
V = \tfrac12\,\bm a^{\mathsf T}\bm G\,\bm a ,
```

whose coefficients are read off directly:

```math
M_{ab} = \int_\Omega \bm u^{(a)}\!\cdot\bm u^{(b)}\,dV , \qquad
C_{ab} = \mathrm{Oh}\!\int_\Omega 2\eta\, \bm e^{(a)}\!:\!\bm e^{(b)}\,dV , \qquad
G_{ab} = \frac{4\pi}{2l+1}(l-1)(l+2)\,\phi_a(1)\phi_b(1) .
```

The first two are quadratures over the drop. The third needs none, because the surface energy
of the previous section depends on the coordinates only through
``\zeta_l = \sum_k a_k\phi_k(1)``.

Differentiating a quadratic form is immediate, so Lagrange's equations become

```math
\bm M\ddot{\bm a} + \bm C\dot{\bm a} + \bm G\bm a = \bm Q ,
```

a damped linear system. Each matrix is the Hessian of a scalar functional, for instance
``M_{ab} = \partial^2 T/\partial\dot a_a\partial\dot a_b``. A second derivative does not care
about the order it is taken in, so all three are symmetric identically, rather than through a
cancellation that could fail.

Every entry needs one derivative of the velocity. The strong form needs two, since it carries
``\nabla^2\bm u``. That is the practical payoff of the weak form, and it is why the trial
functions can be low-order polynomials without the discretisation losing accuracy.

## What the stiffness cannot see

``\bm G`` is rank one for a single mode, as the outer product ``\phi_a(1)\phi_b(1)`` shows.
Surface energy depends on the interface, the interface is the boundary trace, and the interior
enters only through that one number per mode.

The consequence is structural. The restoring force cannot distinguish two interior fields with
the same surface amplitude, so any information about interior structure has to reach the answer
through ``\bm M`` and ``\bm C``. Refining ``K`` refines the inertia and the dissipation, never
the capillarity. That also explains the asymmetry in the convergence tables: the angular
direction is where surface energy lives and it converges early, while the radial direction
feeds only the two operators that carry interior structure.

## Where the modes couple

For a Newtonian fluid ``\eta`` is constant, the angular integrals collapse by orthogonality of
the Legendre polynomials, and the modes separate. Each ``l`` is an independent damped
oscillator, which is why an exact solution exists at all, and why *The Free Viscous Drop* can
proceed one mode at a time.

Two things break that separation. Contact is one: the film pressure is a function of angle and
does not respect the modal decomposition, so ``\bm Q`` couples every mode to every other. A
state-dependent viscosity is the second, and it is worse, because ``\eta`` then sits inside the
integrand of ``C_{ab}`` and carries its own angular spectrum. The orthogonality that
diagonalised the Newtonian problem no longer applies, and the dissipation matrix becomes dense
and state-dependent. Those two are the subjects of *Contact* and *Shear-Thinning Fluids*.

The formulation does not change. What changes is that ``\bm C`` must be rebuilt as the solution
evolves, and that is the entire cost difference between a Newtonian run and a shear-thinning
one.
