# Choosing a Solver

Three solvers are available. They differ along two axes — how the model is formulated, and how
the contact set is found:

|                    | ranked search         | complementarity        |
| ------------------ | --------------------- | ---------------------- |
| **variational**    | [`simulate`](@ref)    | [`simulate_lcp`](@ref) |
| **nonvariational** | [`solve_drop!`](@ref) | —                      |

`simulate` and `simulate_lcp` take identical arguments and return identical types, so
switching closures there is a one-word change.

All three are maintained because a disagreement between any two of them can be attributed.
`simulate` and `simulate_lcp` share a formulation and differ only in the closure, so a
difference between them is the closure's. `simulate` and `solve_drop!` share a closure and
differ in the formulation, so a difference between them is the formulation's. A single solver
cannot distinguish the two.

## The two axes

### Formulation: what is in the state

The **variational** solvers carry interior displacement amplitudes and obtain the surface as a
boundary trace. The equations of motion follow from one Lagrangian, so the mass, damping and
stiffness operators are second derivatives of kinetic energy, dissipation and surface energy
respectively. Two consequences matter in practice: the interior flow is resolved, so a
viscosity that varies through the drop is representable from first principles; and every
assembled operator is symmetric, which the contact solver exploits.

The **nonvariational** solver carries surface amplitudes only, advancing each mode by a
damped-oscillator equation whose coefficients are Reid's exact finite-Ohnesorge values. There
is no interior state. The linear damping is therefore correct by construction with no
convergence parameter to choose — but there is also nothing against which to evaluate a
shear-dependent viscosity.

### Closure: how the contact set is chosen

The **ranked search** proposes candidate contact counts, discards any that would let the
surface pass through the substrate, and accepts the survivor with the smallest edge residual.
Contact grows or shrinks by at most one collocation point per accepted step; if no candidate is
admissible the step is rejected and the step size halves.

The **complementarity** closure proposes nothing. It assembles the affine map from contact
pressure to gap,

```math
h = A_c\,p + b,
```

and solves for the ``p`` satisfying

```math
h \ge 0, \qquad p \ge 0, \qquad p_i h_i = 0 .
```

The three conditions say that the drop does not enter the substrate, that the substrate pushes
rather than pulls, and that pressure appears only where there is contact. The contact set is an
output, so the rejection-and-halving mode has no analogue, and the contact set is not required
to be contiguous.

## Measured behaviour

All three are measured with one contact definition: contact exists whenever any point of the
surface lies below ``0.02R``, contact time runs from the first such instant to the last, and
restitution compares centre-of-mass speed at those two instants.

| ``We`` | ``Oh`` | ``t_c`` var/srch | var/lcp | nonvar/srch | CoR var/srch | var/lcp | nonvar/srch |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.2 | 0.0373 | 2.7113 | 2.7113 | 2.8003 | 0.8227 | 0.8227 | 0.8179 |
| 0.5 | 0.0373 | 2.4901 | 2.4901 | 2.5512 | 0.7462 | 0.7462 | 0.7450 |
| 1.0 | 0.0767 | 2.3677 | 2.3677 | — | 0.5777 | 0.5777 | — |
| 1.0 | 0.3038 | 2.2735 | 2.2735 | — | 0.3610 | 0.3610 | — |
| 2.0 | 0.3038 | 2.1323 | 2.1323 | 2.2164 | 0.3026 | 0.3026 | 0.2978 |

Reading it along the axes:

- **The closure makes no difference.** Over a wider grid — 35 cases spanning
  ``\mathrm{Oh}\in[0.023,\,0.685]`` and ``We\in[0.05,\,3]`` — contact time is bit-identical in
  all 35 and restitution agrees to 2.4e-4 in the worst case, growing with ``We`` at low
  ``\mathrm{Oh}`` and vanishing elsewhere. Contact time agrees exactly because ``dt`` is
  constant and ``t_c`` is a difference of step times, so two closures that flag the same steps
  return the same number. On the 3000 ppm shear-thinning fluid the two agree to 3.5e-6.
- **So the formulation is the only axis that moves anything**: 3.3 per cent in contact time and
  0.6 per cent in restitution, under a shared closure.
- **The contact set comes out a single patch**, even though complementarity does not require
  one. At ``M = 90`` the free arc at the pole is exactly zero. Non-contiguity survives only as a
  transient at the release edge, in a minority of steps that shrinks with resolution — 12 of 530
  accepted steps at ``M = 30``, 4 of 967 at ``M = 45``, none at ``M = 90``.
- **The searching closure is the less robust of the two.** `solve_drop!` completed three of the
  five cases above; both variational cells completed all five, and `simulate_lcp` also runs
  cases at which `simulate` fails outright.

Wallclock over these cases: the variational solvers take 0.2–2.6 s; the nonvariational search
takes 0.0–7.5 s, erratic because its cost is dominated by rejected steps.

!!! note "An earlier version of this page said otherwise"
    It reported the closure changing contact time by 7.5 to 10.2 per cent, and an annular
    contact grid-converged at an arc of about 11 degrees. Both were artefacts. `contact_lcp`
    was symmetrising a compliance that is asymmetric by about forty per cent before handing it
    to a sweep that assumes symmetry, so the pressures were not a solution of the contact
    problem and the accepted states had the drop up to 4.7 per cent of a radius inside the
    substrate. An unresolved penetration looks like a dimple from the outside. Solved against
    the true compliance, the difference is gone.

### Cost of the complementarity solve

Per step the compliance ``A_c = H A^{-1} Q_n`` costs ``M+1`` back-substitutions against
``A = \beta^2\bm M + \beta\bm C + \bm G``, and that assembly, not the complementarity solve,
is where the time goes:

| ``M`` | degrees of freedom | LCP size | compliance | LCP solve | whole step |
| --- | --- | --- | --- | --- | --- |
| 20 | 38 | 11 | 0.20 ms | 0.007 ms | 0.27 ms |
| 45 | 88 | 24 | 0.91 ms | 0.005 ms | 1.05 ms |
| 90 | 178 | 46 | 6.41 ms | 0.025 ms | 7.19 ms |

The solve is 0.3 to 2.6 per cent of a step and converges in a handful of sweeps, so
complementarity is not what makes a run expensive. A variable viscosity is: it forces a fresh
coupled assembly *and* a fresh compliance on every inner sweep, where the constant-viscosity
path reuses both.

## Why there is no nonvariational complementarity solver

The combination works and was built, but it is not carried, and the reason is worth stating
because it also explains the structure of the two that are.

Measured against its own searching closure it agreed to 0.0 per cent on both contact time and
restitution — the two closures select the same contact set at every step, so the cell never
disagreed with `solve_drop!` and bought no information.

The work is preserved on the `archive/nonvariational-lcp` branch.

## Neither compliance is symmetric, and what that costs

It would be convenient if the gap-versus-pressure map were symmetric, because the
complementarity problem would then be exactly the KKT system of the convex programme

```math
\min_{p \ge 0} \; \tfrac{1}{2} p^{\mathsf T} A_c\, p + b^{\mathsf T} p,
```

with a guaranteed solution, uniqueness when the matrix is definite, and a projected
Gauss–Seidel sweep to find it. It is not symmetric in either formulation — by about forty per
cent in the variational one and more than fifty in the nonvariational — so both use active-set
pivoting, which assumes nothing about symmetry.

The reason is the same in both cases and it is structural rather than accidental. Symmetry of
``A_c = H A^{-1} Q_n`` requires the forcing to be the transpose of the constraint Jacobian,
``Q_n = -H^{\mathsf T} W`` for some positive diagonal — that is what makes a constraint and its
multiplier conjugate. The solver constrains the **vertical** gap at **collocation nodes** but
forces the drop with a **radial** pressure projected in the **Legendre** basis, which departs
from conjugacy in two independent ways at once. Measured, the shipped forcing is 95 per cent
orthogonal to any conjugate one, so this is not a small correction.

Nothing is wrong as a result: the problem is still a well-posed complementarity problem with a
unique solution, and active-set pivoting solves it exactly, which is why the two closures agree
to the precision reported above. What is given up is the convexity, and with it the option of
the cheaper sweep and a uniqueness proof.

The derivation of what conjugacy would require, and what adopting it would cost the model, is
in `julia/derivations/contact_conjugacy_derivation.jl`.

## Which to use

- **Shear thinning, or any shear-dependent viscosity** — variational only. The nonvariational
  formulation has no interior state to evaluate ``\eta(\dot\gamma)`` against.
- **Oldroyd-B** — nonvariational only. The extension has not been ported to the variational
  formulation.
- **Exact linear damping with nothing to tune** — nonvariational; Reid's coefficients are built
  in. The variational solver reaches the same damping but needs enough radial basis functions:
  ``K = 1`` is potential flow and over-damps at experimental Ohnesorge, and ``K = 2`` is the
  practical minimum.
- **Robustness across a parameter sweep** — `simulate_lcp`. It cannot reject a step for want of
  an admissible candidate.
- **Cross-checking a result** — run `simulate_lcp` against `solve_drop!`. They share neither
  formulation nor closure.
