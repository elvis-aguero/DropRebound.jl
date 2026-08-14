# DropRebound.jl

Julia solver for the impact and rebound of a liquid drop on a flat substrate, with
Newtonian, shear-thinning (Carreau–Yasuda) and viscoelastic (Oldroyd-B) rheology.

**[Documentation](https://elvis-aguero.github.io/DropRebound.jl/dev/)** — the model,
its derivation, the solvers, and the full API.

<img src="outputs/figures/impact_ob.gif" width="48%"> <img src="outputs/figures/impact_st.gif" width="48%">

*Left: Oldroyd-B (Oh = 0.30, De₁ = 0.5, β_s = 0.5). Right: Carreau shear-thinning
(Oh = 0.30, λ_c = 0.02). Both at We = 0.5, M = 90.*

## What it does

A drop falls, deforms on contact, spreads, and rebounds. The model is stated
**variationally** — kinetic energy, viscous dissipation and surface energy, with
Lagrange's equations giving the evolution — and what distinguishes it from the usual
modal treatment is that **the interior flow is part of the state**. The generalised
coordinates are amplitudes of an interior displacement field; the surface is their
boundary trace, and nothing is eliminated. That is what makes the damping exact
rather than Lamb's small-viscosity estimate, which over-damps by 37 % at $l = 2$ at
the Ohnesorge of the reference experiments and moves restitution by 16 %.

The same structure carries shear thinning without further approximation: $\eta$ is
evaluated pointwise from the shear-rate invariant of the full strain field, and since
that invariant does not superpose over modes, the viscosity acquires angular
structure and couples the shape modes.

Contact is a unilateral constraint — gap and film pressure cannot both be non-zero,
neither may be negative, and the contact extent is an unknown rather than an input.

**Validated** against Gabbard et al. (2025), spanning $\mathrm{Oh} \in [0.014, 0.79]$,
and against a 3000 ppm shear-thinning solution whose zero-shear Ohnesorge is 57 — so a
Newtonian drop that viscous would not rebound at all, and the measured rebound exists
*because* the fluid thins. Nothing is fitted to the impact data. Current residuals are
in the figures and in `outputs/csv/`, not quoted here: they move with every re-sweep,
and a number in prose does not.

![3000 ppm shear-thinning drop: ten simulations against 72 experiments](outputs/figures/figure_shear_thinning.png)

## Quick start

Julia 1.12+, standard library only. From the repository root:

```julia
using DropSolver

# M sets shape modes l = 2..M, film-pressure harmonics l = 0..M, and M+1
# collocation angles at once, which is what makes the contact system square.
# K is the number of radial trial functions per mode; K = 1 is potential flow.
p = ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 45, K = 2)
r = simulate(p)

r.cor, r.tc          # restitution, and contact time in units of √(ρR³/γ)

# Shear thinning: pass η(γ̇)/η₀. Oh is then the ZERO-SHEAR Ohnesorge.
eta = gd -> carreau(gd; lambda_c = 10.0, a = 2.0, n = 0.5, eta_inf_ratio = 0.01)
rst = simulate(ImpactParams(We = 1.0, Bo = 0.0189, Oh = 0.3038, M = 14, K = 2, eta = eta))
```

`simulate_lcp` takes identical arguments and returns the same type, resolving contact
by complementarity rather than by an active set; `solve_drop!` is a second formulation
that eliminates the interior and carries the Oldroyd-B extension. *Choosing a Solver*
in the documentation says which to reach for.

## Repository layout

The package sits at the root in the ordinary Julia arrangement, so `--project=.` is the
environment for the solver, the tests and the scripts alike.

| | |
|---|---|
| `src/`, `test/` | the `DropSolver` module, and its test suite |
| `scripts/` | everything that *uses* the solver: sweeps, audits, figures, animations |
| `derivations/` | the from-scratch CAS derivations, which are also the theory chapters of the documentation |
| `data/` | measured input — experimental metrics and digitised published curves |
| `outputs/` | **everything any script writes**: `csv/` and `figures/` are versioned, `jld/` and `logs/` are not |

**A figure carries the name of the script that draws it** — `figure_vorticity.jl`
writes `figure_vorticity.png`, and where one script draws several the script name is
the stem. `impact_ob.gif` and `impact_st.gif` above predate the rule and have no
generating script in the repository.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## License and reference

MIT, see [LICENSE](LICENSE).

> Agüero-Vera et al., *Spectral simulation of viscous drop impact with Oldroyd-B
> rheology*, in preparation.
