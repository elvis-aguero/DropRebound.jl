# Contact

A free drop oscillates forever in linear theory, decaying at a rate the previous chapter fixes.
A bouncing drop does something a free drop cannot: it meets a wall, stays in touch with it for a
finite time, and leaves. Everything interesting about the bounce happens in that interval, and
the interval is not known before the calculation.

That is the difficulty. The region of the surface in contact is not an input. It appears, grows,
shrinks and vanishes as part of the solution, so the problem has a boundary whose location is
one of the unknowns.

## The drop does not touch the wall

At the pressures and speeds of interest a film of air survives between the liquid and the solid,
and the drop rebounds without wetting. The film is thin against the radius, so it is not
resolved. What it does is transmit a pressure ``p_c(\theta,t)`` to the lower surface, and that
pressure is the only channel through which the wall acts on the drop.

Two statements characterise it, and both are inequalities. The film cannot pull, because a gas
layer supports compression and not tension. And the surface cannot pass through the solid. Where
the film is thin the pressure is finite and the clearance is zero; where the surface has lifted
away the clearance is positive and the pressure is not there.

Write ``h(\theta,t)`` for the clearance between the surface and the substrate. Then

```math
h \ge 0 , \qquad p_c \ge 0 , \qquad h\,p_c = 0 .
```

Together these are the Signorini conditions, the standard statement of unilateral contact. The
third is the interesting one. It does not say where contact occurs. It says that at every point,
one of the two quantities vanishes, and which one is left open. The contact region is defined as
the set where the pressure is the nonzero member of the pair, so it is read off the solution
rather than supplied to it.

Nothing here is specific to drops. The same three conditions describe a rod resting on a table.
What is specific is what fills in ``h``.

## The clearance in terms of the state

A surface point at polar angle ``\theta`` sits at radius ``r(\theta) = 1 + \sum_l \zeta_l
P_l(\cos\theta)``, and the drop's centre of mass sits at height ``z``. The height of that point
above the substrate is therefore

```math
h(\theta) \;=\; z \;+\; \cos\theta\left(1 + \sum_l \zeta_l P_l(\cos\theta)\right) .
```

Two features of this expression matter later. It is affine in the state, so the constraint is
linear even though the contact problem is not. And it carries a factor ``\cos\theta``, because
the constraint is on vertical clearance while ``\zeta_l`` is a radial displacement. Near the
south pole the two nearly coincide; away from it they do not.

The film pressure is carried as a Legendre series ``p_c = \sum_l p_{c,l} P_l(\cos\theta)``, and
the generalised force it exerts on mode ``l`` is

```math
Q_l \;=\; -\frac{4\pi}{2l+1}\,p_{c,l} .
```

The ``l = 0`` and ``l = 1`` harmonics do no work on the shape. The first changes volume, which is
forbidden, and the second translates the drop, which is the centre-of-mass equation rather than a
shape mode. The ``l = 1`` harmonic is not idle: it is the entire net force lifting the drop.

## Collocation

Both the constraint and the pressure are continuous in ``\theta`` and both have to be reduced to
finitely many numbers. The pressure is truncated at ``l = M``, giving ``M+1`` coefficients, and
the contact conditions are imposed at ``M+1`` angles. The counts are tied because the system must
be square.

The nodes are ``\theta = \pi`` together with the zeros of ``P_M``. Legendre zeros cluster toward
the ends of the interval, which places nodes densely near the pole, which is where contact forms
and where the resolution is wanted. The pole itself is included because it is the first point to
touch and the last to leave.

At each node one of two equations is imposed: the clearance vanishes if the node is in contact,
or the pressure vanishes if it is not. Which of the two is applied at which node is precisely the
unknown that the Signorini conditions decline to specify.

## Two ways to close the system

Suppose the contact set were known. Then the choice at each node is known, the system is square
and linear, and one solve advances the step. Both available closures exploit this, and they
differ in how they arrive at the set.

**Ranked search** treats the number of contacting nodes as the quantity to determine. It proposes
candidate counts, discards any candidate whose solution would put part of the surface below the
substrate, and accepts the survivor whose edge residual is smallest. The count is allowed to move
by at most one node per step, which encodes the physical expectation that a contact patch grows
and retreats continuously. If no candidate is admissible the step is rejected and the step size
halves.

**Complementarity** proposes nothing. Eliminating the interior state from the step equations
leaves the nodal clearances affine in the nodal pressures,

```math
\bm h \;=\; \bm A_c\,\bm p + \bm b , \qquad
\bm A_c = \bm H \bm A^{-1} \bm Q_n ,
```

where ``\bm A = \beta^2\bm M + \beta\bm C + \bm G`` is the step operator, ``\bm H`` is the
constraint Jacobian, and ``\bm Q_n`` maps nodal pressure to generalised force. The entries of
``\bm A_c`` are influence coefficients: how much the clearance at node ``i`` opens per unit
pressure at node ``j``. Solving

```math
\bm h \ge 0 , \qquad \bm p \ge 0 , \qquad p_i h_i = 0
```

for ``\bm p`` is a linear complementarity problem, and its solution determines the contact set as
a by-product. There is no candidate to reject, so a step is never refused for lack of an
admissible guess, and the contact region is not required to be a single patch.

The two agree. Over 35 cases spanning ``\mathrm{Oh}\in[0.023,\,0.685]`` and
``\mathrm{We}\in[0.05,\,3]``, contact time is identical in all 35 and restitution agrees to
``2.4\times10^{-4}`` at worst. On a shear-thinning fluid they agree to ``3.5\times10^{-6}``.

The agreement is informative rather than trivial. Complementarity is free to return an annular
contact, with the pole released while a ring still presses, and it does not. At ``M = 90`` the
free arc at the pole is exactly zero. Non-contiguous sets appear only as brief transients at the
release edge, in a minority of steps that shrinks with resolution: 12 of 530 accepted steps at
``M = 30``, 4 of 967 at ``M = 45``, none at ``M = 90``. The single-patch assumption built into
the search is a consequence of the dynamics rather than an approximation imposed on it.

## Is the contact problem convex?

This question has a practical edge. If ``\bm A_c`` were symmetric and positive semidefinite, the
complementarity problem would be the optimality condition of

```math
\min_{\bm p \ge 0} \;\; \tfrac12\,\bm p^{\mathsf T}\bm A_c\,\bm p + \bm b^{\mathsf T}\bm p ,
```

a convex quadratic programme. Existence would follow, uniqueness would follow from definiteness,
and a projected Gauss–Seidel sweep would converge to the answer.

Symmetry is not an accident of the discretisation. It is a statement about conjugacy. Suppose the
constraint reads ``\bm h = \bm H\bm\xi + \bm b`` and the reaction force is the one that does no
virtual work on motions respecting the constraint. That force is

```math
\bm Q = -\bm H^{\mathsf T}\bm\lambda ,
```

with ``\bm\lambda`` the multiplier, and then

```math
\bm A_c = \bm H\bm A^{-1}\bm Q_n = -\bm H\bm A^{-1}\bm H^{\mathsf T}
```

is symmetric for any ``\bm H``, because ``\bm A`` is. Symmetry of the contact operator is
therefore equivalent to the pressure being the multiplier conjugate to the clearance.

Here it is not, and the two expressions above show why. The constraint Jacobian has entries

```math
H_{il} \;=\; \cos\theta_i\,P_l(\cos\theta_i)
```

read off the clearance, while the forcing has entries ``-\tfrac{4\pi}{2l+1}\delta_{jl}`` read off
the generalised force. These differ in two independent ways. The clearance carries
``\cos\theta_i`` and the forcing does not, because the constraint is vertical while the pressure
acts along the radius. And the clearance is sampled at nodes while the pressure is expanded in
harmonics, so even the index sets do not match.

The consequence is measurable: ``\bm A_c`` departs from symmetry by about forty per cent, and the
shipped forcing is 95 per cent orthogonal to any conjugate one. This is not a small correction to
be absorbed by a tolerance.

Nothing is broken. The problem remains a well-posed linear complementarity problem with a
solution, and active-set pivoting finds it exactly, which is why the two closures agree to the
precision quoted above. What is given up is convexity, and with it the cheap sweep and the
uniqueness argument. Recovering it would mean constraining the radial clearance at the nodes and
forcing with nodal pressures, which changes the model rather than the code.
