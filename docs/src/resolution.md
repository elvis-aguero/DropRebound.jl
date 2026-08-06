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
\chi_l(x) \;=\; A\,x^{l+1} \;+\; B\,x\,j_l(qx) , \qquad q^2 = -\sigma_l/\mathrm{Oh} ,
```

an irrotational part plus a vortical one. Setting ``K = 1`` keeps only the first. That is a
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

Empirically ``K \gtrsim |q|/1.5``, with ``|q| \approx \sqrt{\omega_l/\mathrm{Oh}}`` and
``\omega_l \sim l^{3/2}``. The requirement therefore grows with mode number and with falling
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

**Newtonian.** ``M = 45`` to ``90``, ``K = 3``. Restitution is converged in ``M`` by 45 and in
``K`` by 3; ``K = 4`` and ``K = 5`` agree with ``K = 3`` to four decimal places.

**Shear thinning.** ``M = 45``, ``K = 3``. On the 3000 ppm fluid at ``\mathrm{We} = 0.3643``:

| ``M`` at ``K = 3`` | 14 | 30 | 45 | 60 |
|---|---|---|---|---|
| CoR | 0.74305 | 0.75165 | 0.752992 | 0.75315 |

| ``K`` at ``M = 14`` | 2 | 3 | 4 | 5 |
|---|---|---|---|---|
| CoR | 0.73711 | 0.74305 | 0.74233 | 0.74235 |

Refining ``M`` from 45 to 60 moves restitution by 0.02 per cent, and ``K`` from 3 to 4 by
0.10 per cent.

A shear-thinning run costs far more than a Newtonian one at the same truncation, because a
viscosity that depends on the solution rebuilds the dissipation operator at every iteration of
the fixed-point closure. Cost grows roughly as ``M^3``: 9 s at ``M = 14``, 68 s at ``M = 30``,
250 s at ``M = 45``. The pairwise strain contractions are geometry and are cached per basis,
which is what brings those figures within reach.

## Practical limits

Not every pair ``(M, K)`` runs. On the shear-thinning fluid, ``M = 45, K = 4`` and
``M = 60, K = 3`` both stop before completing, while ``M = 45, K = 3`` completes. In each case
the contact solve is exact, with a complementarity residual near ``10^{-15}``, and it is the
fixed-point iteration on the viscosity that fails to reach its tolerance. Reducing the time step
does not help; it has been taken to ``10^{-10}`` without effect.

For ``M = 60, K = 3`` the iteration stalls just above ``10^{-8}``, so raising `eta_tol` to
``10^{-7}`` admits the step, and a control at ``M = 45, K = 3`` returns the same restitution to
six decimal places at either tolerance. That is a usable workaround rather than an account of
what sets the floor. ``M = 45, K = 4`` fails at ``10^{-6}`` as well, and is a separate matter.

The recommended settings sit below these limits, and the convergence tables above are measured
within them.
