# Contact

A free drop rings down at a rate the previous chapter fixes, and touches nothing while it does.
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
the clearance is sampled at ``M+1`` collocation angles. Tying the two counts makes the map from
pressures to clearances square, which is what the elimination below needs.

The nodes are ``\theta = \pi`` together with the zeros of ``P_M``. Legendre zeros cluster toward
the ends of ``\mu = \cos\theta \in [-1,1]``, so nodes are dense near both poles, and in
particular near the south pole, where contact forms and where the resolution is wanted. The pole
itself is included because it is the first point to touch and the last to leave.

Not all of these nodes carry an unknown. A node on the upper hemisphere is a great distance from
the substrate and can never be in contact, so its pressure is zero throughout and it is dropped.
The equator has to be dropped for a sharper reason: the clearance row of a node at
``\theta = \pi/2`` is identically zero, because every entry of it carries the factor
``\cos\theta`` that converts radial displacement to vertical clearance. Such a node has no
influence on the clearance anywhere and no clearance of its own to constrain, and admitting it
would put a zero row and column into the system. At odd ``M`` the polynomial ``P_M`` has a root
at ``\mu = 0`` exactly, so this node is present and must be removed by a tolerance rather than by
the sign of ``\cos\theta``, which at ``M = 45`` evaluates to ``-3\times10^{-16}`` and would
otherwise be admitted.

What remains is the south pole together with the zeros of ``P_M`` in ``\mu < 0``, so the contact
problem carries

```math
\left\lfloor M/2 \right\rfloor + 1
```

unknowns: 23 at ``M = 45`` and 46 at ``M = 90``, against ``M+1`` collocation nodes and ``M+1``
pressure coefficients. The truncation ``M`` therefore buys contact resolution at half the rate
the mode count suggests.

At each retained node one of two equations is imposed: the clearance vanishes if the node is in
contact, or the pressure vanishes if it is not. Which of the two is applied at which node is
precisely the unknown that the Signorini conditions decline to specify.

## Two ways to close the system

Suppose the contact set were known. Then the choice at each node is known, the system is square
and linear, and one solve advances the step. Both closures exploit this, and they differ in how
they arrive at the set.

**A primal active set** reads the answer off the two inequalities. A free node that has penetrated
says the contact is too small, and a contacting node whose pressure pulls says it is too large.
So grow while any free node lies below the substrate, release while the outermost contacting node
pulls, and stop when neither holds. Each move is forced by a violated condition rather than chosen
by scoring candidates, so there is no tie to break and the iteration terminates.

**Complementarity** proposes nothing at all. It asks instead how the clearance responds to
the pressure, and reads the contact set off the answer.

Over one step the discretisation replaces the time derivatives, so the equations of motion
become an algebraic system in the end-of-step state. Collect the shape amplitudes ``\zeta_l``
into a vector ``\bm\xi``, write ``\bm p`` for the pressures at the retained nodes, and let
``\beta`` be the coefficient the scheme puts in front of ``d/dt``. The shape equation is then

```math
\underbrace{\left(\beta^2\bm M + \beta\bm C + \bm G\right)}_{\textstyle \bm A}\bm\xi
\;=\; \bm f \;+\; \bm Q_n\,\bm p ,
```

where ``\bm M``, ``\bm C`` and ``\bm G`` are the mass, damping and stiffness matrices assembled
in the previous chapter. The boldface ``\bm M`` is that matrix and not the truncation degree
``M``, which is an unhappy collision this page inherits and lives with. Here ``\bm f`` collects
everything known from previous steps, and ``\bm Q_n`` maps nodal pressures to generalised forces
by interpolating them onto the Legendre basis and applying
``Q_l = -\tfrac{4\pi}{2l+1}p_{c,l}``. The scheme is BDF2, so ``\beta = c_0/\Delta t`` with
``c_0 = (1+2r)/(1+r)`` and ``r`` the ratio of successive steps; on the opening step it falls
back to BDF1 and ``c_0 = 1``.

The clearance is affine in the state, ``\bm h = \bm H\bm\xi + \bm 1 z + \bm b_0``, so it responds
to the pressure through two channels rather than one. The shape channel is obtained by
substituting ``\bm\xi = \bm A^{-1}(\bm f + \bm Q_n\bm p)``. The second is the centre of mass: the
``l = 1`` harmonic is the net force, and by the same discretisation
``m\beta^2 z = \cdots - m\,p_{c,1}``, so raising ``p_{c,1}`` lifts the whole drop. Writing
``\bm v_1^{\mathsf T}\bm p = p_{c,1}`` for the row that reads that coefficient off the nodal
pressures, both eliminations together give

```math
\boxed{\;\bm h \;=\; \bm A_c\,\bm p + \bm b ,\qquad
\bm A_c \;=\; \bm H\bm A^{-1}\bm Q_n \;-\; \frac{1}{\beta^2}\,\bm 1\,\bm v_1^{\mathsf T} ,\qquad
\bm b \;=\; \bm H\bm A^{-1}\bm f + \bm b_0 \;}
```

The minus sign is not a lift downwards. Contact sits at ``\theta = \pi``, where
``P_1(\cos\theta) = -1``, so a pressure bump against the substrate has ``p_{c,1} < 0`` and the
term raises the drop.

Every entry of ``\bm A_c`` is a derivative,

```math
(A_c)_{ij} \;=\; \frac{\partial h_i}{\partial p_j} ,
```

the clearance opened at node ``i`` by unit pressure at node ``j``, and ``\bm b`` is the
clearance the step would produce with no contact force at all. The contact conditions

```math
\bm h \ge 0 , \qquad \bm p \ge 0 , \qquad p_i h_i = 0
```

then close the system: a linear complementarity problem in ``\bm p`` alone, solved by pivoting.

Both closures therefore read the contact set from the same two inequalities. They differ in that
the active set walks to it one node at a time from the previous step's answer, while the
complementarity solve treats every node's status as an independent unknown and may return a
contact region that is not a single patch.

The two agree, and where they do not is as informative as where they do. Over a grid of 25
impacts spanning ``\mathrm{Oh}\in[0.02,\,0.7]`` and ``\mathrm{We}\in[0.1,\,3]`` at ``M = 45``,
restitution agrees to ``2.9\times10^{-5}`` at worst, with the largest discrepancies along the
low-viscosity edge ``\mathrm{Oh} = 0.02``. Contact time is identical in every one of them, to the
last bit. That is not a general guarantee. It holds because no step was rejected anywhere on this
grid, so both closures ran the same step sequence, and ``t_c`` is a difference of step times.
Where steps are rejected the two histories diverge and the contact times separate a little. On
the 3000 ppm shear-thinning fluid at the production truncation ``M = 90``, where the viscosity is
rebuilt inside every sweep and steps are rejected, restitution agrees to ``7.1\times10^{-7}`` and
contact time to ``2.9\times10^{-5}``.

One of the 25 is missing from that comparison. At ``\mathrm{We} = 3`` and
``\mathrm{Oh} = 0.02``, the fastest and least viscous corner, the active set never releases: it
reaches the end of the march with the drop still on the substrate, while complementarity bounces.
Reading the contact set node by node from the previous step's answer is the weaker of the two
methods at the edge of the parameter range.

Complementarity is also free to return an annular contact, with the pole released while a ring
still presses, and it very nearly does not. Counting the accepted steps whose contact set is not
a single run of adjacent nodes gives 0.7% at ``\mathrm{We} = 0.5`` and 6.4% at
``\mathrm{We} = 2`` when ``M = 45``, falling to 0.4% and 1.7% at ``M = 90``. These are brief
transients at the release edge, and they thin with resolution without reaching zero. A contact
patch that stays a patch is therefore a result of the dynamics rather than an assumption imposed
on them, but it is a result that holds to a few parts in a thousand and not exactly.

## Is the contact problem convex?

This question has a practical edge. If ``\bm A_c`` were symmetric and positive semidefinite, the
complementarity problem would be the optimality condition of

```math
\min_{\bm p \ge 0} \;\; \tfrac12\,\bm p^{\mathsf T}\bm A_c\,\bm p + \bm b^{\mathsf T}\bm p ,
```

a convex quadratic programme. Existence would follow, uniqueness would follow from definiteness,
and a projected Gauss–Seidel sweep would converge to the answer.

Symmetry is not an accident of the discretisation. It is a statement about conjugacy. Enforce the
constraint ``\bm h = \bm H\bm\xi + \bm b_0 \ge 0`` with a multiplier ``\bm\lambda \ge 0``, which
adds ``-\bm\lambda^{\mathsf T}\bm h`` to the potential. The force it contributes is minus the
gradient of that term,

```math
\bm Q \;=\; -\frac{\partial}{\partial\bm\xi}\left(-\bm\lambda^{\mathsf T}\bm h\right)
      \;=\; \bm H^{\mathsf T}\bm\lambda ,
```

so the load a constraint exerts is read off the constraint's own derivative. Then

```math
\bm A_c = \bm H\bm A^{-1}\bm Q_n = \bm H\bm A^{-1}\bm H^{\mathsf T}
```

is symmetric for any ``\bm H``, because ``\bm A`` is, and positive semidefinite with it.
Symmetry of the contact operator is therefore equivalent to the pressure being the multiplier
conjugate to the clearance.

Here it is not, and the reason is interpolation rather than geometry.

The constraint is read at nodes, ``H_{il} = \cos\theta_i\,P_l(\cos\theta_i)``, while the forcing
is read from harmonics, ``Q_l = -\tfrac{4\pi}{2l+1}\,p_{c,l}``. A unit pressure at one node
therefore does not enter as a load at that node. The film pressure is carried as a degree-``M``
Legendre field, so that unit pressure enters as the Galerkin force of the polynomial interpolating
a spike there, and that polynomial oscillates over the whole sphere. The constraint is local and
its supposed multiplier is global, so the two are not dual.

The measurement that identifies interpolation as the cause, rather than the missing
``\cos\theta_i``, is the alignment between a node's constraint row and the force its own pressure
produces:

```math
\cos\angle\left(\bm H_{i\cdot},\, \bm Q_n\bm e_i\right)
\;=\; \frac{\bm H_{i\cdot}\cdot\bm Q_n\bm e_i}{\|\bm H_{i\cdot}\|\,\|\bm Q_n\bm e_i\|} .
```

A missing ``\cos\theta_i`` would make this unity at the south pole, where vertical and radial
coincide. Measured at ``M = 45`` it is 0.02 at the pole and 0.28 at the far end of the film: worst
exactly where that explanation predicts it should be perfect.

The resulting asymmetry is large. In the norm
``\|\bm A_c - \bm A_c^{\mathsf T}\|/\|\bm A_c\|`` it is 0.43 at ``M = 20`` and 0.37 at
``M = 45``.

The problem is still a well-posed linear complementarity problem with a solution, and pivoting
finds it exactly, which is why the two closures agree to the precision above. What is given up is
convexity, and with it the projected sweep and the uniqueness argument.

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
