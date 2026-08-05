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
| 0.2 | 0.0373 | 2.7113 | 2.4336 | 2.8003 | 0.8227 | 0.8602 | 0.8179 |
| 0.5 | 0.0373 | 2.4901 | 2.0523 | 2.5512 | 0.7462 | 0.7394 | 0.7450 |
| 1.0 | 0.0767 | 2.3677 | 1.8217 | — | 0.5777 | 0.5162 | — |
| 1.0 | 0.3038 | 2.2735 | 2.2735 | — | 0.3610 | 0.3602 | — |
| 2.0 | 0.3038 | 2.1323 | 2.1323 | 2.2164 | 0.3026 | 0.3023 | 0.2978 |

Reading it along the axes:

- **Restitution is largely insensitive to both choices.** The closure moves it by 0.9 per cent
  and the formulation by 0.6 per cent in the median. Two formulations and two contact
  algorithms that share no contact logic land on nearly the same rebound speed, which is the
  strongest available evidence that the rebound speed is a property of the physics rather than
  of an implementation. The exception is instructive: at the lowest Weber number the two
  closures differ by 4.5 per cent, and at ``\mathrm{Oh} = 0.3038`` they agree to four decimals.
  The disagreement lives at low Ohnesorge and low Weber, not everywhere.
- **Contact time is sensitive to the closure**, by 10.2 per cent within the variational
  formulation, with `simulate_lcp` consistently the shorter and the gap widening as Ohnesorge
  falls. Complementarity permits a contact set that is not an interval, and the runs use one:
  contact detaches at the pole and leaves an annulus, grid-converged at an arc of about 11
  degrees. A ranked search cannot represent that, because it only ever proposes a contiguous
  patch grown from the pole, and a patch that stays attached at the centre releases later.
- **The formulation matters least of the three effects**, 3.3 per cent in contact time under a
  shared closure.
- **The searching closure is the less robust of the two.** `solve_drop!` completed three of the
  five cases; both variational cells completed all five, and `simulate_lcp` also runs cases at
  which `simulate` fails outright.

Wallclock over these cases: the variational solvers take 0.2–1.7 s; the nonvariational search
takes 0.0–7.3 s, erratic because its cost is dominated by rejected steps.

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

It also needed a *different* complementarity algorithm. The variational compliance matrix is
symmetric, structurally: it is a Hessian of an energy. That makes the problem equivalent to the
convex quadratic program

```math
\min_{p \ge 0} \; \tfrac{1}{2} p^{\mathsf T} A_c\, p + b^{\mathsf T} p,
```

which a projected Gauss–Seidel sweep solves. The nonvariational model is assembled from
per-mode damping coefficients and is the second derivative of nothing, so its compliance comes
out asymmetric — by more than half its own magnitude — and that equivalence fails. Existence
survives, since positive definiteness makes the matrix a P-matrix and a P-matrix
complementarity problem has exactly one solution for every right-hand side, but the sweep must
be replaced by active-set pivoting. Two closures whose answers never differ are not worth two
algorithms.

The work is preserved on the `archive/nonvariational-lcp` branch.

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
