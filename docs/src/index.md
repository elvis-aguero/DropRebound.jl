# DropRebound.jl

A Julia spectral solver for the impact and rebound of a liquid drop on a solid surface, for
Newtonian, shear-thinning (Carreau–Yasuda, Cross) and viscoelastic (Oldroyd-B) fluids.

A drop falls under gravity, flattens against the surface, spreads, and leaves. The package
solves the linearised spectral equations for that sequence and reports what an experiment
measures: how long the drop stays in contact, and how much of its speed it gets back.

## The model

The model is stated **variationally**. Three quadratic functionals — kinetic energy, viscous
dissipation and surface energy — determine the motion through Lagrange's equations:

```math
T[\dot{\bm\xi}] = \tfrac12\int |\bm u|^2\,dV, \qquad
\Phi[\dot{\bm\xi}] = \int 2\eta\, \bm e\!:\!\bm e \,dV, \qquad
V[\bm\xi] = \tfrac12\sum_l \frac{4\pi}{2l+1}(l-1)(l+2)\,\zeta_l^2,
```

giving ``\bm M\ddot{\bm\xi} + \bm C\dot{\bm\xi} + \bm G\bm\xi = \bm Q``.

What distinguishes this from the usual modal treatment is that **the interior flow is part of
the state**. The generalised coordinates are amplitudes of an interior *displacement* stream
function — Gegenbauer functions in angle, and the Taylor powers ``x^{l+1}, x^{l+3}, \dots`` of
the exact spherical-Bessel profile in radius. The surface amplitudes are its boundary trace;
nothing is eliminated.

Keeping the interior is what makes the damping exact rather than a small-viscosity estimate. At
the Ohnesorge number of the reference experiments, Lamb's formula over-damps by 37 per cent at
``l = 2`` and 143 per cent at ``l = 8``, and correcting it moves the coefficient of restitution
by 16 per cent. One radial function per mode (``K = 1``) is potential flow and reproduces
Lamb; ``K = 2`` resolves the interior and recovers Reid's exact rates.

The same structure carries **shear thinning** with no further approximation. The viscosity is
evaluated pointwise from the shear-rate invariant of the full strain field, and because that
invariant does not superpose over modes, ``\eta`` acquires angular structure and couples the
shape modes within a Gaunt band.

Contact is a **unilateral constraint**: the gap and the film pressure cannot both be non-zero,
neither may be negative, and the contact extent is an unknown rather than an input.

## The three solvers

Two formulations crossed with two ways of finding the contact set:

|                    | ranked search         | complementarity        |
| ------------------ | --------------------- | ---------------------- |
| **variational**    | [`simulate`](@ref)    | [`simulate_lcp`](@ref) |
| **nonvariational** | [`solve_drop!`](@ref) | —                      |

All three are maintained so that a disagreement between any two of them can be attributed to
one cause rather than to two at once. *Choosing a Solver* gives the measured differences and a
recommendation per use case.

## Getting started

```julia
using DropSolver

# Newtonian, at the parameters of the Gabbard et al. (2025) production sweep.
# M sets the shape modes l = 2..M, the film-pressure harmonics l = 0..M, and the
# M+1 collocation angles that make the contact system square. K is the number of
# radial trial functions per mode.
p = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 45, K = 2)
r = simulate(p)

r.cor          # coefficient of restitution
r.tc           # contact time, in units of √(ρR³/σ)
r.cp           # contact extent (node count) at each step
r.pc           # film-pressure coefficients
```

Shear thinning passes ``\eta(\dot\gamma)/\eta_0``, and `Oh` is then the zero-shear Ohnesorge
number:

```julia
eta = gd -> carreau(gd; lambda_c = 10.0, a = 2.0, n = 0.5, eta_inf_ratio = 0.01)
rst = simulate(ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 14, K = 2, eta = eta))
```

A variable viscosity forces the coupled operators to be rebuilt every step, so shear-thinning
runs use a smaller `M` than Newtonian ones, where the operator is block diagonal and built
once.

## Validation

Against Gabbard et al. (2025): coefficient of restitution to 8 per cent median and contact time
to 13 per cent, across 935 experiments spanning ``\mathrm{Oh} \in [0.014, 0.79]``.

Against a 3000 ppm shear-thinning solution whose zero-shear Ohnesorge number is 57 — a
Newtonian drop that viscous does not rebound at all, so the rebound the experiments measure
exists *because* the fluid thins. Restitution agrees to 7 per cent and contact time to 10 per
cent, with nothing fitted to the impact data.

## Documentation

- **Choosing a Solver** — the two axes, what each is worth in measured terms, and which of the
  three to use.
- **Newtonian Theory** — Reid's exact finite-Ohnesorge damping and frequency, and how they are
  wired into the solver.
- **Shear-Thinning Fluids** — begin with *Shear-Thinning Drops*, which states the model
  completely; then the descent from it to something runnable, the constitutive law a real fluid
  is characterised by, and the measurement that prices the mode-coupling concessions.
- **Viscoelastic Fluids** — the Oldroyd-B extension, which runs on the nonvariational solver.
- **API Reference** — the exported interface, grouped in the order a run meets it, with
  unexported internals kept separate at the end.

Each derivation page is executed in continuous integration, so the expressions and numbers
shown are computed rather than transcribed.

## References

- Lamb, H. *Hydrodynamics*, 6th ed., Cambridge University Press (1932).
- Chandrasekhar, S. "The oscillations of a viscous liquid globe."
  *Proc. London Math. Soc.* **9**, 141–149 (1959).
- Reid, W. H. "The oscillations of a viscous liquid drop."
  *Quart. Appl. Math.* **18**, 86–89 (1960).
- Richard, D., Clanet, C. & Quéré, D. "Contact time of a bouncing drop."
  *Nature* **417**, 811 (2002).
- Molaček, J. & Bush, J. W. M. "A quasi-static model of drop impact."
  *Phys. Fluids* **24**, 127103 (2012).
- Gabbard, C. T. *et al.* "Drop rebound at low Weber number." *J. Fluid Mech.* (2025).
