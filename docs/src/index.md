# DropRebound.jl

A Julia spectral solver for the impact and rebound of a liquid drop on a solid
surface, for Newtonian, viscoelastic (Oldroyd-B) and shear-thinning
(Carreau-Yasuda, Cross) fluids.

## Method

The drop's shape is expanded in Legendre modes,
``\zeta(\theta,t) = R\sum_{l\ge2} a_l(t)P_l(\cos\theta)``, and each mode
evolves as a damped oscillator

```math
\ddot a_l + 2\lambda_l\,\dot a_l + \omega_l^2\,a_l = \text{forcing},
```

with contact against the surface entering as a pressure forcing.

The damping rate and frequency come from the exact roots of Reid's (1960)
characteristic equation, which is valid at arbitrary Ohnesorge number. Passing
`viscous = :lamb` selects the classical small-viscosity limit
``\lambda_l = \mathrm{Oh}\,(l-1)(2l+1)``, ``\omega_l^2 = l(l-1)(l+2)`` instead;
the two agree as ``\mathrm{Oh}\to0`` and diverge as ``\mathrm{Oh}`` and ``l``
grow.

![Damping rate against Ohnesorge number, comparing Lamb's small-viscosity
formula with Reid's exact result for three mode numbers.](assets/reid_vs_lamb.png)

*Why the exact relation is worth the trouble: Lamb's formula (dashed) and
Reid's exact damping rate (solid) agree only in the small-viscosity limit.*

## Getting started

```julia
using DropSolver

M  = 6                        # number of Legendre modes
Oh = 0.05                     # Ohnesorge number
Bo = 1e-6                     # Bond number

# Pressure-integral table ∫ Pₙ(u)/u³ du over each interval of the near-pole
# angle grid. `precompute_integrals` returns (table, angles); passing NaN asks
# it to build its own default grid, and `[1]` keeps the table.
pressure_integrals = precompute_integrals(NaN, M)[1]

cfg = SimConstants(M, M+1, Oh, Bo,
                   make_theta_vec(M),   # Gauss-Legendre collocation angles
                   pressure_integrals,
                   make_dt_max(M);      # CFL-like step ceiling for M modes
                   viscous = :reid)     # exact finite-Oh damping / frequency

init      = DropState(M)
init.A[2] = 0.05              # excite mode l = 2
init.z    = 1.05              # centre of mass, in drop radii above the wall
init.v    = -0.5              # impact velocity (downward)
init.dt   = make_dt_max(M)

times, states = solve_drop!(cfg, OBParams(), init; t_end = 20.0)
kpis = extract_kpis(times, states, cfg)
```

`kpis` reports the first bounce: contact time ≈ 2.74 capillary times,
coefficient of restitution ≈ 0.76, maximum contact radius ≈ 0.37 drop radii.

`OBParams()` is the Newtonian default. Viscoelastic runs pass a non-trivial
`OBParams`; shear-thinning runs pass `STExactParams` through `solve_drop!`'s
`stx` keyword.

## Documentation

- **Newtonian Theory** — Reid's derivation of the exact finite-Ohnesorge
  damping and frequency, and how they are wired into the solver.
- **Shear-Thinning Fluids** — the Carreau-Yasuda and Cross models. Begin with
  *A Hierarchy of Models*, which sets out the assumptions each closure rests
  on and the range over which they hold.
- **Viscoelastic Fluids** — the Oldroyd-B extension.
- **API Reference** — the exported solver interface, grouped by the order a run
  meets it, with unexported internals kept separate at the end.

Each derivation page is executed in continuous integration, so the expressions
and numbers shown are computed rather than transcribed.

## References

- Lamb, H. *Hydrodynamics*, 6th ed., Cambridge University Press (1932).
- Chandrasekhar, S. "The oscillations of a viscous liquid globe."
  *Proc. London Math. Soc.* **9**, 141–149 (1959).
- Reid, W. H. "The oscillations of a viscous liquid drop."
  *Quart. Appl. Math.* **18**, 86–89 (1960).
- Molaček, J. & Bush, J. W. M. "A quasi-static model of drop impact."
  *Phys. Fluids* **24**, 127103 (2012).
