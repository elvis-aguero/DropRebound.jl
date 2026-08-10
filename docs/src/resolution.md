# Resolution and Convergence

Two integers set the size of the discrete problem. ``M`` truncates the angular expansion, so
the retained shape modes are ``l = 2\ldots M``, the film pressure carries harmonics
``l = 0\ldots M``, and contact is imposed at ``M+1`` collocation angles. ``K`` truncates the
radial expansion within each mode, and fixes how much interior structure that mode can hold.

The two are not interchangeable, and the natural instinct is the wrong one.

## What each buys

Restitution for a Newtonian drop at ``\mathrm{Oh} = 0.0373``, ``\mathrm{We} = 0.5``:

| refinement | change in CoR |
|---|---|
| ``M``: 30 to 45, at ``K = 2`` | 0.09 % |
| ``M``: 30 to 45, at ``K = 3`` | 0.10 % |
| ``K``: 2 to 3, at ``M = 30`` | 0.98 % |

A fifty per cent increase in angular resolution moves the answer by a tenth of a per cent. The
same increase in radial resolution moves it by one per cent. Angular convergence arrives early
and radial convergence does not, so a run at ``M = 90, K = 2`` spends 178 degrees of freedom
refining a direction that had already converged.

## Why the radial direction is the expensive one

The exact radial profile of a freely oscillating viscous drop, derived in *The Free Viscous
Drop*, is

```math
\chi_l(x) \;=\; A\,x^{l+1} \;+\; B\,x\,j_l(qx) , \qquad q^2 = \sigma_l/\mathrm{Oh} ,
```

an irrotational part plus a vortical one, with ``A`` and ``B`` fixed by the two free-surface
conditions. Here ``j_l`` is the spherical Bessel function of the first kind, and ``\sigma_l``
is the complex decay rate of mode ``l``, defined by writing the free oscillation as
``\chi_l \propto e^{-\sigma_l t}``, so that ``\mathrm{Re}\,\sigma_l > 0`` is decay and
``\mathrm{Im}\,\sigma_l`` is the oscillation frequency. The wavenumber follows from
substituting that form into the unsteady Stokes equation, which turns it into a Helmholtz
equation ``\nabla^2\bm u + (\sigma_l/\mathrm{Oh})\bm u = 0``. Setting ``K = 1`` keeps only the first. That is a
potential flow, which still dissipates, because an irrotational field has nonzero strain. It is
Lamb's calculation, and it over-predicts damping at every Ohnesorge number tested, from
+5 per cent at ``\mathrm{Oh} = 0.01`` to +170 per cent at ``\mathrm{Oh} = 1``. The error is
one-signed because a free surface requires zero tangential stress, potential flow cannot supply
it, and a field denied its vorticity layer must strain more to move the same surface.

The radial trial functions are the Taylor terms of that Bessel function, so how many are needed
follows ``|q|``:

| ``\mathrm{Oh}`` | 1.0 | 0.3038 | 0.05 | 0.01 |
|---|---|---|---|---|
| ``\|q\|`` at ``l = 2`` | 1.1 | 3.0 | 7.5 | 16.8 |
| ``K`` for damping within 1 % | 2 | 2 | 3 | 4 |

Empirically ``K \gtrsim |q|/1.5``, with ``|q| \approx \sqrt{\omega_l/\mathrm{Oh}}``, where
``\omega_l = \sqrt{l(l-1)(l+2)}`` is Rayleigh's inviscid capillary frequency and so grows like
``l^{3/2}``. That estimate replaces ``\sigma_l`` by its inviscid limit, which is why it
reproduces the table above except in the first column: at ``\mathrm{Oh} = 1`` mode 2 is close
to critically damped, ``\sigma_l`` is no longer close to ``i\omega_l``, and the entry there is
taken from the exact root instead. The requirement therefore grows with mode number and with falling
viscosity.

## The ceiling on ``K``

Every radial trial function carries a factor ``x^{l+1}``. At large ``l`` that confines all of
them to a layer of thickness ``\sim 1/l`` beneath the surface, where they become increasingly
alike, and the mode's mass matrix loses conditioning. The rate is about one and a half decimal
digits per added function at ``l = 2`` and two and a half at ``l = 90``.

The trial functions are therefore taken as ``x^{l+1}P_{k-1}(2x^2-1)`` rather than the raw
powers ``x^{l+1}, x^{l+3}, \ldots``. Both families span the same space and both give ``K = 1``
as potential flow exactly, and where both are numerically sound they agree to eight decimal
places. The difference is conditioning, which improves by up to ten orders of magnitude at low
``l``. See [`RitzBasis`](@ref).

[`radial_window`](@ref) reports how many functions a mode may carry,

```math
K(l) \;=\; \left\lfloor 0.9 \times 47.2\,l^{-0.642} \right\rfloor ,
```

a power law fitted to the measured ceiling at a normalised condition number of ``10^{10}``, then
shaved by ten per cent so that it sits below that ceiling rather than on it. The window bounds the
conditioning of the mass matrix, which is a necessary condition for a run to be trustworthy and
not a sufficient one: a pair ``(M, K)`` inside the window can still fail for reasons that have
nothing to do with conditioning. The window
closes as ``l^{-0.64}`` while the requirement grows as ``l^{3/4}``, so the two cross. Above the
crossing the high modes are under-resolved by construction, which is tolerable because they
carry little: the highest retained mode holds about ``3\times10^{-5}`` of the surface energy in
a representative run.

## Recommended settings

The defaults are ``M = 90``, ``K = 3``. Two different questions decide them, and they do
not have the same answer.

**Restitution converges early.** It moves 0.09 per cent between ``M = 45`` and ``M = 90``,
and ``K = 4`` and ``K = 5`` agree with ``K = 3`` to four decimal places. On that evidence
alone ``M = 45`` would do.

**The contact patch does not.** Contact is imposed at ``M+1`` collocation nodes, so the
contact radius can only take as many distinct values as there are nodes inside the patch.
Measured on a Newtonian drop, the patch covers a converged physical angle but a
resolution-dependent number of nodes:

| ``\mathrm{We}`` | ``M`` | nodes in contact | half-angle | steps where the patch jumped |
|---|---|---|---|---|
| 0.5 | 30 | 9 | 45.8° | 6 |
| 0.5 | 45 | 12 | 42.5° | 6 |
| 0.5 | 90 | 24 | 45.2° | 8 |
| 3.0 | 45 | 16 | 58.4° | 27 |
| 3.0 | 90 | 31 | 59.2° | 48 |

The half-angle is converged by ``M = 30``. The number of nodes spanning it is not: at
``M = 45`` the whole growth and retreat of contact is described by a dozen or so discrete
radii. A quantity integrated over the bounce, like restitution, is insensitive to that. A
quantity that depends on the contact history is not, and this package exists to compute
the second kind.

The last column is the one to watch. It counts accepted steps where the patch moved by more
than one node, and it *rises* with resolution rather than falling. At high Weber number the
contact edge outruns the grid whatever the grid, so refining buys a finer description of
the patch without buying continuity of its motion.

**Cost.** A Newtonian run at ``M = 90, K = 3`` takes 18 s with the active-set closure and
50 s with complementarity, against 1.3 s and 6.8 s at ``M = 45``.

## Practical limits

Not every pair ``(M, K)`` runs. On the shear-thinning fluid, ``M = 45, K = 4`` and
``M = 60, K = 3`` both stop before completing, while ``M = 45, K = 3`` completes. In each case
the contact solve is exact, with a complementarity residual near ``10^{-15}``, and it is the
fixed-point iteration on the viscosity that fails to reach its tolerance. Reducing the time step
does not help; it has been taken to ``10^{-10}`` without effect.

### The viscosity tolerance, and why it depends on the truncation

`eta_tol` is the Picard tolerance on the interior strain rate, relative and in the infinity
norm. It is **not a fixed number**: the default is

```math
\texttt{eta\_tol}(M) \;=\; 10^{-8}\,M^{3},
```

which is ``9.1\times10^{-4}`` at ``M = 45`` and ``7.3\times10^{-3}`` at ``M = 90``.

It scales with ``M`` because the reachable residual does. The Picard iterate is
``\bm{\dot a}^{*} = \beta\bm a + \bm h``, with ``\beta = c_0/\Delta t`` and
``\Delta t \sim M^{-3/2}``, so the sweep-to-sweep increment is a difference of two large
quantities and its floor rises with the truncation. Measured on the 3000 ppm fluid the floor is
below ``10^{-6}`` at ``M = 45`` and ``M = 60``, and ``3.0\times10^{-6}`` at ``M = 90``; a cubic
passes through all three. The default sits about three decades above that fit.

Setting `eta_tol` below the floor is not slow, it is fatal, and the mechanism is worth knowing
because it is counter-intuitive. A step whose Picard iteration misses the tolerance is rejected
and ``\Delta t`` is halved. But halving ``\Delta t`` doubles ``\beta``, and the residual floor
rises as ``1/\Delta t`` — so the retry is further from the tolerance than the attempt that
failed. Left alone the march spirals down to `dt_min` and dies having computed part of a bounce.
The solver bounds the number of consecutive halvings spent this way and warns when it gives up:

```
┌ Warning: the viscosity iteration is not recovering as the step shrinks;
│ no smaller step will help, so the march stops here (raise eta_tol, or lower M)
```

If you see that warning, raise `eta_tol` or lower ``M``. Do not lower `dt_min`.
