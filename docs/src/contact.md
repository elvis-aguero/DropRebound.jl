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

**A primal active set** reads the answer off the two inequalities. A free node that has penetrated
says the contact is too small, and a contacting node whose pressure pulls says it is too large.
So grow while any free node lies below the substrate, release while the outermost contacting node
pulls, and stop when neither holds. Each move is forced by a violated condition rather than chosen
by scoring candidates, so there is no tie to break and the iteration terminates.

**Complementarity** proposes nothing at all. Eliminating the interior state from the step equations
leaves the nodal clearances affine in the nodal pressures,

```math
\bm h \;=\; \bm A_c\,\bm p + \bm b , \qquad
\bm A_c = \bm H \bm A^{-1} \bm Q_n ,
```

where ``\bm A = \beta^2\bm M + \beta\bm C + \bm G`` is the step operator, ``\bm H`` is the
constraint Jacobian ``\partial\bm h/\partial\bm\xi``, and ``\bm Q_n`` maps nodal pressure
to generalised force. Here ``\bm M``, ``\bm C`` and ``\bm G`` are the mass, damping and
stiffness matrices of *Variational Mechanics*, and ``\beta`` is the coefficient the time
discretisation puts in front of ``d/dt``, so that ``\dot{\bm\xi} \to \beta\bm\xi +
(\text{known history})`` over one step. For backward Euler ``\beta = 1/\Delta t``. Boldface
distinguishes these matrices from the harmonic truncation ``M`` used above. The entries of
``\bm A_c`` are influence coefficients: how much the clearance at node ``i`` opens per unit
pressure at node ``j``. Solving

```math
\bm h \ge 0 , \qquad \bm p \ge 0 , \qquad p_i h_i = 0
```

for ``\bm p`` is a linear complementarity problem, solved by pivoting. Both closures therefore
determine the contact set from the same two inequalities. They differ in that the active set walks
to it one node at a time from the previous step's answer, while the complementarity solve treats
every node's status as an independent unknown and is free to return a contact region that is not a
single patch.

The two agree. Over 35 cases spanning ``\mathrm{Oh}\in[0.023,\,0.685]`` and
``\mathrm{We}\in[0.05,\,3]``, contact time is identical in all 35 and restitution agrees to
``2.4\times10^{-4}`` at worst. On a shear-thinning fluid they agree to ``3.5\times10^{-6}``.

The agreement carries information. Complementarity is free to return an annular contact, with the
pole released while a ring still presses, and it does not: at ``M = 90`` the free arc at the pole
is exactly zero, and non-contiguous sets survive only as brief transients at the release edge,
vanishing as resolution rises. A contact patch that stays a patch is therefore a result of the
dynamics rather than an assumption imposed on them.

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

Here it is not, and the reason rewards a moment's care, because the obvious explanation is wrong.

The constraint is read at nodes: ``H_{il} = \cos\theta_i\,P_l(\cos\theta_i)``. The forcing is
read from harmonics: ``Q_l = -\tfrac{4\pi}{2l+1}p_{c,l}``. One might guess the mismatch is the
factor ``\cos\theta_i``, present in the clearance because the constraint is vertical and absent
from the forcing because the pressure acts along the radius. That guess makes a prediction. At the
south pole the vertical and radial directions coincide, so the two should align there and diverge
toward the equator.

The measurement runs the other way. Take the cosine between a node's constraint row
``\bm H_{i\cdot}`` and the generalised force that a unit pressure at that same node produces. At
``M = 45`` it is 0.02 at the pole and 0.28 at the far end of the film. Alignment is worst exactly
where the geometric argument says it should be perfect, so that is not the cause.

The cause is interpolation. The film pressure is carried as a degree-``M`` Legendre field, so a
unit pressure at one node does not enter the equations as a load at that node. It enters as the
Galerkin force of the polynomial that interpolates a spike there, and that polynomial oscillates
over the whole sphere. The constraint is local and its supposed multiplier is global. Duality
fails not because the directions disagree but because the pressure variable is not attached to a
place.

Measured, ``\bm A_c`` departs from symmetry by 0.43 at ``M = 20`` and 0.37 at ``M = 45``, in the
norm ``\|\bm A_c - \bm A_c^{\mathsf T}\|/\|\bm A_c\|``.

Nothing is broken. The problem remains a well-posed linear complementarity problem with a
solution, and active-set pivoting finds it exactly, which is why the two closures agree to the
precision quoted above. What is given up is convexity, and with it the cheap sweep and the
uniqueness argument.

## Recovering convexity

Since the failure is that the multiplier is not attached to a place, the repair is to make it so.
Take the contact unknown to be the vertical **load** at each node rather than a pressure field.
Differentiating the constraint then leaves no freedom in the forcing, because ``\bm H^{\mathsf T}``
and ``\bm 1^{\mathsf T}`` are by definition the derivatives of the clearance with respect to the
shape and to the centre of mass:

```math
\bm A\bm\xi = \bm f + \bm H^{\mathsf T}\bm\lambda , \qquad
m\ddot z = -m\,\mathrm{Bo} + \bm 1^{\mathsf T}\bm\lambda ,
```

with ``m = 4\pi/3`` the drop's mass in these units, ``\bm\lambda`` the vector of nodal loads,
and ``\bm 1`` a column of ones, so that ``\bm 1^{\mathsf T}\bm\lambda`` is their sum.
Eliminating both gives

```math
\bm h = \bm W\bm\lambda + \bm b , \qquad
\bm W = \bm H\bm A^{-1}\bm H^{\mathsf T} + \frac{1}{m\beta^2}\,\bm 1\bm 1^{\mathsf T} ,
```

and every term is symmetric positive semidefinite: the first because ``\bm A`` is, the second
because it is a positive multiple of an outer square. The problem is then exactly the convex
programme written above, with a solution always and a unique one when ``\bm W`` is definite.

This route is implemented and selectable as `force_mode = :nodal`. Its price is that the
multiplier is a load rather than a pointwise pressure, so pressure is recoverable only as a
diagnostic, and not at the pole at all, because the quadrature weight there vanishes to machine
precision. The default remains the Legendre field, which gives the pressure directly and solves
without needing symmetry.
