# Oscillations of a Viscous Liquid Drop

*A self-contained pedagogical expansion of Reid (1960), adapted for this repo.*

This is the physics underneath every rheology model in this repo. Every other
derivation here (Oldroyd-B, Carreau-Yasuda) is a correction layered on top of
what gets derived across this chapter; read this one first.

The full CAS derivation and its live cross-check against the running solver
live in `julia/derivations/reid1960_full_derivation.jl`. Starting from the
next chapter, this text's load-bearing results are typeset directly from
that script's own verified symbolic objects, not transcribed by hand -- the
handful of formulas cited in this introductory chapter (the inviscid
frequency below, Lamb's correction) are exactly that, citations, the same
way Reid's own paper takes them as given rather than re-deriving them.

## What problem are we solving?

Consider a small liquid droplet (say, a water drop suspended in air) that has
been slightly deformed from its equilibrium spherical shape. Surface tension
acts as a restoring force, and the drop will oscillate. Viscosity damps those
oscillations. The question is: *at what frequency does it oscillate, and how
fast does it decay?*

This couples three distinct pieces of physics:

1. The inviscid capillary (surface-tension-driven) oscillation modes,
   characterized by the frequency ``\sigma_{l;0}`` for a deformation of
   spherical-harmonic order ``l``.
2. Viscous dissipation inside the bulk of the drop.
3. The coupling between surface deformation, internal flow, and the pressure
   field, all of which must be consistent at the boundary.

Reid's paper resolves all three simultaneously for *arbitrary* viscosity, by
reducing the problem to a single transcendental characteristic equation. Its
central insight is a rescaling: surface tension enters only through the
inviscid frequency ``\sigma_{l;0}``. After non-dimensionalization, the
characteristic equation is identical to Chandrasekhar's equation for a
self-gravitating viscous liquid globe -- a problem about surface tension turns
out to be, mathematically, the same problem as one about gravity.

## Historical context and relation to Molaček & Bush (2012)

Three analytical milestones precede Reid's result:

- **Rayleigh (1879) / Lamb (1932):** the inviscid frequencies (below).
- **Lamb (1881):** the small-viscosity correction,
  ``\sigma_{l;\nu} = (l-1)(2l+1)\nu/R^2 \pm i\,\sigma_{l;0}``.
- **Chandrasekhar (1959):** the full arbitrary-viscosity solution for a
  *self-gravitating* globe.

Reid (1960) shows the surface-tension problem is identical to Chandrasekhar's,
thereby completing the picture for surface-tension-driven drops.

The result is used directly in Molaček & Bush (2012), where the coefficients
``A_m`` and ``D_m`` governing kinetic energy and viscous dissipation for each
surface mode are extracted from Reid's characteristic equation (derived in
full later in this chapter).

## Notation

| symbol | meaning |
|:--|:--|
| ``R`` | equilibrium (undeformed) drop radius |
| ``\rho`` | fluid density |
| ``\mu = \rho\nu`` | dynamic / kinematic viscosity |
| ``T_1`` | surface tension |
| ``x = r/R`` | dimensionless radial coordinate |
| ``\theta,\varphi`` | polar, azimuthal angle |
| ``l`` | spherical-harmonic degree of the surface deformation |
| ``Y_l^m(\theta,\varphi)`` | spherical harmonic (angular shape of the deformation) |
| ``\epsilon(t) = \epsilon_0 e^{-\sigma t}`` | dimensionless amplitude of the deformation |
| ``\sigma`` | complex decay rate; ``\mathrm{Re}(\sigma)>0`` is decay, ``\mathrm{Im}(\sigma)`` is the oscillation frequency |
| ``\sigma_{l;0}`` | inviscid oscillation frequency of mode ``l`` (real, no viscosity) |
| ``U(x)`` | dimensionless radial-velocity eigenfunction (all of the unknown physics lives in this one function) |
| ``q^2 = \sigma R^2/\nu`` | viscous wavenumber |
| ``\alpha^2 = \sigma_{l;0}R^2/\nu`` | inviscid frequency, in the same units |
| ``j_l(z)`` | spherical Bessel function of the first kind, order ``l`` |
| ``Q_{l+1/2}(q) = j_{l+1}(q)/j_l(q)`` | the one Bessel-function combination the whole problem reduces to |
