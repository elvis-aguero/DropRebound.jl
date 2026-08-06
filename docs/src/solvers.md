# Choosing a Solver

Three solvers are available. They differ along two axes, namely how the model is formulated and
how the contact set is found:

|                    | active set            | complementarity        | ranked search         |
| ------------------ | --------------------- | ---------------------- | --------------------- |
| **variational**    | [`simulate`](@ref)    | [`simulate_lcp`](@ref) | not carried           |
| **nonvariational** | not carried           | not carried            | [`solve_drop!`](@ref) |

[`simulate`](@ref) and [`simulate_lcp`](@ref) take identical arguments and return identical
types, so switching closures is a one-word change. [`Backend`](@ref) names a combination of the
two axes if a study needs to sweep over solvers.

Three are kept rather than one so that a disagreement can be narrowed down. `simulate` and
`simulate_lcp` share a formulation and differ only in the closure, so any difference between them
belongs to the closure alone. `solve_drop!` differs from both in formulation and in closure at
once, so a difference against it is bounded by the two together rather than assigned to either.
The first comparison is clean; the second is a cross-check, not an attribution.

## The formulation axis

The **variational** solvers carry interior displacement amplitudes and read the surface off as a
boundary trace, as *Variational Mechanics* sets out. Two consequences matter here. The interior
flow is resolved, so a viscosity that varies through the drop can be evaluated from first
principles. And every assembled operator is a Hessian of a scalar functional, hence symmetric,
which the contact solve relies on.

The **nonvariational** solver carries surface amplitudes only. Each mode advances as a damped
oscillator whose coefficients are Reid's exact finite-Ohnesorge values, so the linear damping is
right by construction with no radial resolution to choose. There is no interior state, which is
also the limitation: there is nothing against which to evaluate a shear-dependent viscosity.

## What the axes are worth

All three are measured against one contact definition. Contact exists whenever any surface point
lies within ``0.02R`` of the substrate, contact time runs from the first such instant to the last, and restitution
compares centre-of-mass speed at those two instants.

| ``\mathrm{We}`` | ``\mathrm{Oh}`` | ``t_c`` var/act | var/lcp | nonvar | CoR var/act | var/lcp | nonvar |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.2 | 0.0373 | 2.7113 | 2.7113 | 2.8003 | 0.8227 | 0.8227 | 0.8179 |
| 0.5 | 0.0373 | 2.4901 | 2.4901 | 2.5512 | 0.7462 | 0.7462 | 0.7450 |
| 1.0 | 0.0767 | 2.3677 | 2.3677 | n/a | 0.5777 | 0.5777 | n/a |
| 1.0 | 0.3038 | 2.2735 | 2.2735 | n/a | 0.3610 | 0.3610 | n/a |
| 2.0 | 0.3038 | 2.1323 | 2.1323 | 2.2164 | 0.3026 | 0.3026 | 0.2978 |

The two variational columns are indistinguishable, for the reasons *Contact* gives: both closures
read the contact set off the same two inequalities. Contact time agrees exactly because ``dt`` is
fixed and ``t_c`` is a difference of step times, so two closures that flag the same steps return
the same number.

The nonvariational column differs by 3.3 per cent in contact time and 0.6 per cent in restitution.
That gap is the formulation and its closure together, since `solve_drop!` differs in both, and
nothing here separates the two contributions.

Robustness separates them where accuracy does not. `solve_drop!` completed three of the five
cases above and both variational cells completed all five. `simulate_lcp` also runs cases at which
`simulate` fails outright: the active set walks to the contact set from the previous step's
answer, so it can stall when that answer is far from the new one, while the pivoting solve starts
from no assumption.

Wallclock over these cases is 0.2 to 2.6 s for the variational solvers. The nonvariational
search ranges from 0.0 to 7.5 s, erratic because rejected steps dominate its cost.

## What the complementarity solve costs

Per step the compliance ``\bm A_c = \bm H\bm A^{-1}\bm Q_n`` costs ``M+1`` back-substitutions
against ``\bm A``, and that assembly rather than the complementarity solve is where the time
goes:

| ``M`` | degrees of freedom | LCP size | compliance | LCP solve | whole step |
| --- | --- | --- | --- | --- | --- |
| 20 | 38 | 11 | 0.20 ms | 0.007 ms | 0.27 ms |
| 45 | 88 | 23 | 0.91 ms | 0.005 ms | 1.05 ms |
| 90 | 178 | 46 | 6.41 ms | 0.025 ms | 7.19 ms |

The solve is 0.3 to 2.6 per cent of a step and converges in a handful of pivots, so
complementarity is not what makes a run expensive. A variable viscosity is. It forces a fresh
coupled assembly and a fresh compliance on every inner sweep, where the constant-viscosity path
reuses both.

## Which to use

**A shear-dependent viscosity** needs a variational solver. The nonvariational formulation has
no interior state to evaluate ``\eta(\dot\gamma)`` against.

**Oldroyd-B** needs the nonvariational solver, which is where the polymer-stress state lives.

**Exact linear damping with nothing to tune** is the nonvariational solver, which has Reid's
coefficients built in. The variational solver reaches the same damping given enough radial
functions, with ``K = 2`` the practical minimum and ``K = 1`` over-damping.

**A wide parameter sweep** favours `simulate_lcp`, which does not walk to the contact set and so
has no previous answer to be stranded by.

**Cross-checking a result** means running `simulate_lcp` against `solve_drop!`, since they share
neither formulation nor closure.
