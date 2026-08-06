# Variational Mechanics of a Dissipative Fluid

The home page ended with a claim: the linearised motion of a viscous drop obeys

```math
\frac{d}{dt}\frac{\partial T}{\partial\dot{\bm\xi}}
\;+\; \frac{\partial\mathcal R}{\partial\dot{\bm\xi}}
\;+\; \frac{\partial V}{\partial\bm\xi} \;=\; \bm Q ,
```

with ``T`` the kinetic energy, ``\mathcal R`` the Rayleigh dissipation function, and ``V`` the
surface energy. That statement is exact and useless on its own, because ``\bm\xi`` is a field
and the equation is an identity over an infinite-dimensional space.

This page turns it into matrices. Three things have to happen: the trial fields must be chosen
so that incompressibility holds without a pressure, the angular and radial directions must be
separated, and the three functionals must be reduced to quadratures.

## Coordinates that are displacements

A Lagrangian needs positions and velocities. The natural description of a flow is a velocity
field, which is the wrong half of that pair, so the first choice is what plays the role of
position.

Take the interior displacement. Write the fluid's displacement from the sphere as a field
``\chi_l(x,t)`` for each surface mode ``l``, with ``x = r/R``. The velocity potential is its
time derivative, and the surface amplitude is its boundary trace:

```math
\psi_l = \dot\chi_l , \qquad \zeta_l = \chi_l(1) .
```

The second follows from integrating the kinematic condition ``\dot\zeta = u_r`` in time.

This resolves a question that otherwise looks like a defect in the counting. There is no
evolution equation for ``\zeta_l``, and there should not be, because ``\zeta_l`` is not an
independent coordinate. The surface is the boundary value of the interior, and the interior is
the only thing that moves.

## A trial space with no pressure in it

The weak form removed the pressure by requiring that test fields be divergence-free. That
requirement is now a constraint on the trial space, and constraints of that kind are usually
awkward. In axisymmetry they are free.

Any axisymmetric incompressible flow derives from a Stokes stream function,

```math
\bm u \;=\; \nabla\times\!\left(\frac{\Psi}{r\sin\theta}\,\bm e_\phi\right) ,
```

and a field written this way satisfies ``\nabla\cdot\bm u = 0`` identically, for every ``\Psi``.
The constraint is discharged by the choice of representation rather than enforced afterwards.
Nothing in the discrete problem is a Lagrange multiplier, and no saddle-point system appears.

Separating ``\Psi = f(x)\,C_l(\theta)`` mode by mode gives velocity components

```math
u_r = \frac{f}{x^2}\,P_l(\cos\theta) , \qquad
u_\theta = \frac{f'}{x\,l(l+1)}\,\partial_\theta P_l(\cos\theta) ,
```

and strain-rate components

```math
e_{rr} = \left(\frac{f'}{x^2} - \frac{2f}{x^3}\right)P_l , \qquad
e_{r\theta} = \frac{1}{2xl(l+1)}\left(f'' - \frac{2f'}{x} + \frac{l(l+1)f}{x^2}\right)\partial_\theta P_l ,
```

with ``e_{\theta\theta}`` and ``e_{\phi\phi}`` following in the same way. Their trace vanishes
because Legendre's equation supplies
``\partial_\theta^2 P_l + \cot\theta\,\partial_\theta P_l = -l(l+1)P_l``, and the two surviving
terms cancel. Incompressibility is visible in the algebra, not asserted.

The combination in ``e_{r\theta}`` is worth naming. It is the tangential-stress operator, and
the free-surface condition ``\bm n\cdot\bm\Sigma\cdot\bm t = 0`` is the statement that it
vanishes at ``x = 1``. Nothing in the discretisation imposes that. It emerged as a natural
condition of the weak form, so a converged solution satisfies it and an unconverged one
approaches it. This is the mechanism behind the potential-flow failure described in *Resolution
and Convergence*: a single trial function cannot make that operator vanish, and the residual
tangential stress is paid for as extra dissipation.

## Ritz reduction

Restrict each ``\chi_l`` to a finite span,

```math
\chi_l(x,t) \;=\; \sum_{k=1}^{K} a_k(t)\,\phi_k(x) ,
```

with trial functions that are regular at the origin. The functionals become quadratic forms in
``\bm a``, and the Lagrange equations become

```math
\bm M\ddot{\bm a} + \bm C\dot{\bm a} + \bm G\bm a = \bm Q ,
```

with

```math
M_{ab} = \int \bm u^{(a)}\!\cdot\bm u^{(b)}\,dV , \qquad
C_{ab} = \mathrm{Oh}\!\int 2\eta\, \bm e^{(a)}\!:\!\bm e^{(b)}\,dV , \qquad
G_{ab} = \frac{4\pi}{2l+1}(l-1)(l+2)\,\phi_a(1)\phi_b(1) .
```

Because ``T``, ``\mathcal R`` and ``V`` are quadratic, the Ritz and Galerkin routes coincide
here, and the matrices are Hessians of scalar functionals. They are therefore symmetric by
construction rather than by cancellation, which is a property worth testing and one the suite
does test.

Every entry needs one derivative of the velocity. The strong form needs two, since it carries
``\nabla^2\bm u``. That reduction is the practical payoff of the weak form, and it is why the
trial functions can be low-order polynomials without the discretisation losing accuracy.

## What the stiffness cannot see

``\bm G`` is rank one for a single mode. Surface energy depends on the interface, the interface
is the boundary trace ``\phi_a(1)``, and the interior enters only through that number.

The consequence is structural. The restoring force cannot distinguish two interior fields with
the same surface amplitude, so any information about interior structure has to reach the answer
through ``\bm M`` and ``\bm C``. Refining ``K`` is therefore refining the inertia and the
dissipation, never the capillarity. It also explains the asymmetry in the convergence tables:
the angular direction is where surface energy lives and it converges early, while the radial
direction feeds only the two operators that do carry interior structure.

## Where the modes couple

For a Newtonian fluid ``\eta`` is constant, the angular integrals reduce by orthogonality of
the Legendre polynomials, and the modes separate. Each ``l`` is an independent damped
oscillator, which is why an exact solution exists at all and why *The Free Viscous Drop*
can proceed one mode at a time.

Two things break that separation. Contact is one: the film pressure is a function of angle and
does not respect the modal decomposition, so ``\bm Q`` couples every mode to every other. A
state-dependent viscosity is the second, and it is worse, because ``\eta`` then sits inside
the integrand of ``C_{ab}`` and carries its own angular spectrum. The orthogonality that
diagonalised the Newtonian problem no longer applies, and the dissipation matrix becomes dense
and state-dependent. Those two are the subjects of *Contact* and *Shear-Thinning Fluids*.

The formulation itself does not change. What changes is that ``\bm C`` must be rebuilt as the
solution evolves, and that is the entire cost difference between a Newtonian run and a
shear-thinning one.
