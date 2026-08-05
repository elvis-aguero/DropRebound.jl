# Choosing a Solver

Four solvers are available. They differ along two independent axes, and any of the four
combinations can be run:

|                    | ranked search      | complementarity      |
| ------------------ | ------------------ | -------------------- |
| **variational**    | [`simulate`](@ref) | [`simulate_lcp`](@ref) |
| **nonvariational** | [`solve_drop!`](@ref) | [`solve_drop_lcp!`](@ref) |

Within each row the two solvers take identical arguments and return identical types, so
switching closures is a one-word change.

All four are maintained because a disagreement between any two of them can then be
attributed. Agreement down a column shows the formulation is not responsible for a
difference; agreement across a row shows the contact closure is not. With a single
formulation and a single closure, every disagreement is confounded with every other.

## The two axes

### Formulation: what is in the state

The **variational** solvers carry interior displacement amplitudes and obtain the surface as
a boundary trace. The equations of motion follow from a single Lagrangian, so the mass,
damping and stiffness operators are the second derivatives of kinetic energy, dissipation and
surface energy respectively. Two consequences matter in practice: the interior flow is
resolved, so a viscosity that varies through the drop is representable from first principles;
and every assembled operator is symmetric, which the contact solver can exploit.

The **nonvariational** solvers carry surface amplitudes only, with each mode advanced by a
damped-oscillator equation whose coefficients are Reid's exact finite-Ohnesorge values. There
is no interior state. The damping is therefore correct in the linear regime by construction,
with no convergence parameter to choose — but there is also nothing to evaluate a
shear-dependent viscosity against.

### Closure: how the contact set is chosen

The **ranked search** proposes candidate contact counts, discards any that would let the
surface pass through the substrate, and accepts the survivor with the smallest edge residual.
Contact grows or shrinks by at most one collocation point per accepted step; if no candidate
is admissible the step is rejected and the step size halves.

The **complementarity** closure proposes nothing. It assembles the affine map from contact
pressure to gap,

```math
h = A_c\,p + b,
```

and solves for the unique ``p`` satisfying

```math
h \ge 0, \qquad p \ge 0, \qquad p_i h_i = 0 .
```

The three conditions say that the drop does not enter the substrate, that the substrate
pushes rather than pulls, and that pressure appears only where there is contact. The contact
set is an output. No candidate is proposed, ranked or rejected, so the rejection-and-halving
mode has no analogue, and the contact set is not required to be contiguous.

Complementarity is posed on the lower hemisphere only. The collocation angles span the whole
sphere, and on the upper half the gap expression measures the height of the crown rather than
a distance to the substrate — pushing down there opens no gap underneath, and the diagonal of
``A_c`` changes sign at the equator to say so.

## Measured behaviour

Contact and restitution below use one definition for all four solvers: contact exists whenever
any point of the surface lies below ``0.02R``, contact time runs from the first such instant to
the last, and restitution compares centre-of-mass speed at those two instants.

| ``We`` | ``Oh`` | ``t_c`` var/srch | var/lcp | nonvar/srch | nonvar/lcp | CoR var/srch | var/lcp | nonvar/srch | nonvar/lcp |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.2 | 0.0373 | 2.7113 | 2.5089 | 2.8003 | 2.8003 | 0.8227 | 0.8237 | 0.8179 | 0.8179 |
| 0.5 | 0.0373 | 2.4901 | 2.2406 | 2.5512 | 2.5512 | 0.7462 | 0.7406 | 0.7450 | 0.7450 |
| 1.0 | 0.0767 | 2.3677 | 2.0853 | — | 2.4309 | 0.5777 | 0.5687 | — | 0.5731 |
| 1.0 | 0.3038 | 2.2735 | 2.2735 | — | 2.4739 | 0.3610 | 0.3611 | — | 0.3556 |
| 2.0 | 0.3038 | 2.1323 | 2.1323 | 2.2164 | 2.3107 | 0.3026 | 0.3027 | 0.2978 | 0.2993 |

Reading the table along its axes:

- **Restitution is insensitive to both choices.** Every cell agrees within 0.8 per cent. Two
  formulations and two contact algorithms, sharing no contact logic, land on the same rebound
  speed — which is the strongest evidence available that the rebound speed is a property of
  the physics rather than of an implementation.
- **Contact time is sensitive to the closure, but only in the variational formulation** — 7.5
  per cent there against 0.0 per cent in the nonvariational model, where the two closures
  select the same contact set at every step. So the shorter contact time the variational
  complementarity solver reports is specific to that formulation and not a generic property of
  complementarity. It is consistent with contact detaching at the centre and leaving an
  annulus, a shape a surface-only formulation has no way to represent.
- **The searching closure is the least robust cell.** It completed three of five cases; the
  other three cells completed all five. Both failures are in the nonvariational model at
  moderate Ohnesorge, and its complementarity counterpart runs them.

Wallclock over the same cases: the variational solvers take 0.2–1.5 s, the nonvariational
complementarity solver a steady 1.2–1.3 s, and the nonvariational search 0.0–6.4 s — erratic
because its cost is dominated by rejected steps.

## Structural difference worth knowing

The compliance matrix ``A_c`` is symmetric in the variational formulation and asymmetric in
the nonvariational one, by more than half its own magnitude. Symmetry is structural in the
first case: ``A_c`` is a Hessian of an energy. The second formulation is assembled from
per-mode damping coefficients and is not the second derivative of anything, so nothing forces
symmetry.

This decides the algorithm. A symmetric ``A_c`` makes the complementarity problem equivalent
to the convex quadratic program

```math
\min_{p \ge 0} \; \tfrac{1}{2} p^{\mathsf T} A_c\, p + b^{\mathsf T} p,
```

which a projected Gauss–Seidel sweep solves. Without symmetry that equivalence fails and the
sweep converges to something that is not a solution. What survives is positive definiteness,
which is all that existence requires — a positive definite matrix is a P-matrix, and a
P-matrix complementarity problem has exactly one solution for every right-hand side. The
nonvariational solver therefore uses active-set pivoting, which terminates for P-matrices.

Solutions are verified with a two-sided residual, measuring violation of ``h \ge 0`` and of
``p \ge 0`` separately as well as complementarity. The one-sided quantity
``\max_i \min(h_i, p_i)`` is small whenever one of each pair is small, so it certifies states
in which the drop has passed through the substrate.

## Which to use

- **Shear-thinning or any shear-dependent viscosity** — variational only. The nonvariational
  formulation has no interior state to evaluate a viscosity against.
- **Oldroyd-B** — nonvariational search only. Polymer stress is extra state, so the system is
  not affine and the complementarity problem is not exact; `solve_drop_lcp!` rejects these
  parameters rather than solving an approximation silently.
- **Linear damping accuracy at moderate Ohnesorge without tuning** — nonvariational. Reid's
  exact coefficients are built in. The variational solver reaches the same damping but needs
  enough radial basis functions: ``K = 1`` is potential flow and over-damps at experimental
  Ohnesorge, and ``K = 2`` is the practical minimum.
- **Robustness across a parameter sweep** — either complementarity solver. Neither rejects
  steps for want of an admissible candidate.
- **Cross-checking a result** — run the diagonal. Agreement between `simulate` and
  `solve_drop_lcp!` shares neither formulation nor closure.
