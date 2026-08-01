# DropRebound.jl

A Julia spectral solver for the impact and rebound of a liquid drop on a
solid surface — for Newtonian, viscoelastic (Oldroyd-B), and shear-thinning
(Carreau-Yasuda, Cross) fluids.

## How the solver works

The drop's shape is expanded in Legendre modes,
``\zeta(\theta,t) = R\sum_{l\ge2} a_l(t)P_l(\cos\theta)``, and each mode obeys
a damped-oscillator equation

```math
\ddot a_l + 2\lambda_l\,\dot a_l + \omega_l^2\,a_l = \text{forcing},
```

with contact against the surface entering as a pressure forcing. The two
coefficients are what the theory pages are about. Rather than Lamb's
small-viscosity limit ``\lambda_l = \mathrm{Oh}\,(l-1)(2l+1)``, the solver uses
the exact roots of Reid's (1960) characteristic equation, valid at *any*
Ohnesorge number — already a 13% difference in ``\lambda_2`` at
``\mathrm{Oh} = 0.05``, growing with both ``\mathrm{Oh}`` and ``l``.

## Getting started

```julia
using DropSolver

M  = 6                        # number of Legendre modes
Oh = 0.05                     # Ohnesorge number
Bo = 1e-6                     # Bond number

cfg = SimConstants(M, M+1, Oh, Bo,
                   make_theta_vec(M),
                   precompute_integrals(NaN, M)[1],
                   make_dt_max(M);
                   viscous = :reid)     # :reid (exact) or :lamb (asymptotic)

init      = DropState(M)
init.A[2] = 0.05              # excite mode l = 2
init.z    = 2.0               # initial height
init.dt   = make_dt_max(M)

times, states = solve_drop!(cfg, OBParams(), init; t_end = 20.0)
```

Viscoelastic and shear-thinning runs swap in `OBParams` / `STExactParams`; the
theory pages explain what those parameters mean and where they are valid.

## How to read these pages

Start with **Newtonian Theory** — everything else is a correction on top of it.

For shear-thinning fluids, read **A Hierarchy of Models** first. It lays out
every rung from the exact free-surface problem down to Lamb's limit, and at
each step records the assumption made, what it discards, and how to undo it.
The other pages in that section are individual rungs; reading them without the
map is how earlier work in this repo went wrong.

**Superseded** holds one page kept for provenance. It contains a known error
and says so at the top; it does not describe the current model.

## These pages are executable

Every page is a from-scratch computer-algebra derivation with assertions,
cross-checked against the running solver — the contract is in `CLAUDE.md`. The
page you are reading *is* that script, rendered. Its verification machinery is
hidden so the page reads as a chapter, but CI executes every assertion on every
push and a failure blocks the build. Numbers quoted in the prose are pinned by
those assertions, so they cannot silently drift from what the code does.

That cuts both ways, and is worth stating plainly: a passing assertion means
the algebra is self-consistent, **not** that the physics is right. Where a
derivation is checked against something genuinely independent — Lamb's
``\mathrm{Oh}\to0`` limit, Molaček & Bush's high-``\mathrm{Oh}`` limit,
Chandrasekhar's critical values, or experimental data — the page says so
explicitly.

## Status

| model | standing |
|:--|:--|
| Newtonian (Reid 1960) | validated against two independently published limits |
| Oldroyd-B | derived and CAS-verified |
| Carreau-Yasuda / Cross | **derived, not validated against experiment** |

The shear-thinning extension reproduces measured restitution coefficient and
contact time to a median of roughly 18% and 13%, but that median hides a heavy
tail: about a quarter of sampled impacts are wrong by large factors at the
extremes of the Weber range. Its closure is a postulate rather than a
derivation, and the open questions are catalogued in *A Hierarchy of Models*.
Treat its predictions accordingly.

See the **API Reference** for the production solver's docstrings.
