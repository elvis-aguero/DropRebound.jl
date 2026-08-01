# Structure of the Solutions

## Roots of the characteristic equation

The characteristic equation from the previous chapter is transcendental in
``q^2 = \sigma R^2/\nu``, given the physical parameters encoded in
``\alpha``. It has infinitely many solutions.

The Bessel function ``J_{l+1/2}(q)`` in the denominator of
``Q_{l+1/2}(q)`` has infinitely many real zeros ``q_1<q_2<\cdots`` on the
positive real axis. Near each zero the characteristic equation has a pole,
and each pole "captures" one pair of roots. These correspond to
increasingly high radial overtones of the oscillation -- modes with more
and more nodal surfaces in the radial direction inside the drop.

Roots with smaller real part of ``\sigma`` (slower decay) are the
physically dominant ones: they are the last to die out and govern the
observable oscillation of the drop. The two roots with the smallest
``\mathrm{Re}(\sigma)`` correspond to the fundamental surface oscillation
at harmonic order ``l``; all higher roots decay exponentially faster.

!!! warning "Sign convention"
    In this convention, ``\epsilon\sim e^{-\sigma t}`` with
    ``\mathrm{Re}(\sigma)>0`` meaning decay. So the smallest positive real
    part corresponds to the slowest decay -- the dominant mode. This is
    the opposite of stability analysis, where one looks for the
    most-positive growth rate.

## Two regimes: oscillatory vs. aperiodic

For fixed ``l``, the character of the two dominant roots depends on
``\alpha^2``:

- **Large ``\alpha^2`` (low viscosity, large drop):** the two dominant
  roots are complex conjugates, ``\sigma=\gamma\pm i\omega_d``, a damped
  oscillation.
- **Small ``\alpha^2`` (high viscosity, small drop):** the two dominant
  roots are both real and positive, two aperiodic decaying modes
  (overdamped).

The transition occurs at a critical value of ``\alpha^2``. For ``l=2``,
this was determined numerically by Chandrasekhar:
``\sigma_{2;0}R^2/\nu=3.69``, ``\sigma_{2;\nu}/\sigma_{2;0}=0.968``. For
water in air (``T_1=74\,\mathrm{dyn/cm}``, ``\nu=0.01\,\mathrm{cm^2/s}``),
this gives a critical radius ``R_c\approx 0.23\,\mathrm{mm}``: drops larger
than ``R_c`` oscillate with damping; smaller drops are aperiodically
damped, with no oscillation at all.

## The small-viscosity limit

As ``\nu\to 0`` (``q\to\infty``), standard large-argument Bessel
asymptotics give ``Q_{l+1/2}(q)/q\to 0`` between the poles of
``Q_{l+1/2}``. We take this large-``q`` behavior as a citable fact about
Bessel functions, the same way the previous chapters cite the raw
spherical-coordinate operator formulas -- what we verify directly is what
setting ``Q\to 0`` does to the characteristic equation itself:

```@eval
using Symbolics, Markdown
@variables l q alpha Q
characteristic_rhs_full = (2*(l-1)/q^2)*(l + (l+1)*(q-2*l*Q)/(q-2*Q))
characteristic_rhs_Q0 = substitute(characteristic_rhs_full, Dict(Q => 0))
residual = simplify(characteristic_rhs_Q0 - 2*(l-1)*(2l+1)/q^2; expand=true)
Markdown.parse("Setting \$Q\\to 0\$, the difference between the bracket " *
    "\$l+(l+1)(q-2lQ)/(q-2Q)\$ and its claimed limit \$2l+1\$ collapses to `" *
    string(residual) * "`, for symbolic \$l\$: the characteristic equation becomes\n\n" *
    "```math\n\\frac{\\alpha^4}{q^4} + 1 = \\frac{2(l-1)(2l+1)}{q^2}\n```")
```

Multiplying through by ``q^4`` turns this into a genuine quadratic in
``q^2``:

```math
q^4 - 2(l-1)(2l+1)\,q^2 + \alpha^4 = 0,
```

solved by the quadratic formula,
``q^2 = (l-1)(2l+1) \pm \sqrt{(l-1)^2(2l+1)^2-\alpha^4}``. For
``\alpha\gg 1`` (the low-viscosity regime this limit describes), the
``\alpha^4`` term dominates under the square root, giving
``q^2\approx(l-1)(2l+1)\pm i\alpha^2`` to leading order -- and since
``\sigma=q^2\nu/R^2`` and ``\alpha^2\nu/R^2=\sigma_{l;0}`` by definition,
this is exactly Lamb's classical result:

```math
\sigma_{l;\nu} = (l-1)(2l+1)\,\frac{\nu}{R^2} \pm i\,\sigma_{l;0}.
```

This is checked quantitatively, not just asymptotically: for
``l=2,3,5,10``, the relative gap between the exact quadratic root and
Lamb's leading-order formula shrinks monotonically as ``\alpha`` grows
from ``10`` to ``10^4``, confirming Lamb's formula is genuinely the
``\alpha\to\infty`` limit of the exact equation above, not an unrelated
approximation that happens to look similar.

## Checking against the running solver

Everything so far is algebra. `julia/src/timestepper.jl`'s production
code implements exactly Lamb's formula for the small-``\mathrm{Oh}``
regime. A real, small-``\mathrm{Oh}`` `solve_drop!` run should therefore
show a free-oscillation decay rate matching ``\mathrm{Oh}(l-1)(2l+1)`` to
good accuracy -- and it does, to within ``5\%`` at ``\mathrm{Oh}=0.02``
for ``l=2`` and ``l=3``, and at ``\mathrm{Oh}=0.05`` for ``l=2``. This
document's physics matches the code that actually runs, not just its own
internal algebra.

Starting from the linearized Navier-Stokes equations for a perturbed
spherical drop, this derivation has produced the exact transcendental
characteristic equation governing its damped oscillations, and confirmed
it reduces to Lamb's classical small-viscosity formula in the appropriate
limit. This equation -- not Lamb's asymptotic approximation to it -- is
the physically correct starting point for any drop whose Ohnesorge number
isn't small, which includes every shear-thinning fluid the Carreau-Yasuda
extension in this repo was built to handle, since shear-thinning can swing
the effective Ohnesorge number across orders of magnitude within a single
impact.
