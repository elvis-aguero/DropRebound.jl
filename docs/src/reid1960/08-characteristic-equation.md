# The Characteristic Equation

## Eliminating the surface tension

Define ``\alpha^2 \equiv \sigma_{l;0}R^2/\nu`` (the inviscid frequency,
measured in units of the viscous diffusion rate ``\nu/R^2``) alongside
``q^2=\sigma R^2/\nu`` from Chapter 6. Using
``\sigma_{l;0}^2 = l(l-1)(l+2)T_1/(\rho R^3)`` (Chapter 3) to rewrite BC3's
left side, and ``P_0=\sigma^2R^2\Pi_0/l`` (Chapter 5) to rewrite its right
side, then multiplying BC3 through by ``lR^2/\nu^2``, the left side becomes
``\alpha^4 = \sigma_{l;0}^2R^4/\nu^2`` and the right side becomes purely a
function of ``q`` and the boundary-condition unknowns:

```math
\alpha^4 = q^2\left\{q^2\Pi_0 - 2l\left[U'(1)-2U(1)\right]\right\}.
```

``T_1`` and ``\rho`` have cancelled entirely. The problem in terms of
``(\alpha,q)`` is universal -- this is Reid's key rescaling, and the reason
his result carries over unchanged to any other spherically-symmetric
restoring force (surface tension here, self-gravity in Chandrasekhar's
version of the same problem).

## Solving BC1 + BC2 for ``C, \Pi_0``

BC1 (``U(1)=Cj_l(q)+\Pi_0=-1``) and BC2 (Chapter 7, evaluated on the
general solution) are two linear equations in the two unknowns
``C,\Pi_0``. Solving them exactly and using the Bessel recurrence
``qj_l'/j_l = l-qQ_{l+1/2}(q)`` (Chapter 2) to eliminate ``j_l'`` in favor
of ``Q_{l+1/2}(q)=j_{l+1}(q)/j_l(q)``:

```math
C = \frac{2(l-1)(l+1)}{j_l(q)\,q\,\left(2Q_{l+1/2}(q)-q\right)}.
```

## The characteristic equation itself

Substituting ``C`` and ``\Pi_0`` into ``U'(1)-2U(1)`` (using ``U(1)=-1``
from BC1), then into the ``\alpha^4`` relation above, is the single
largest algebraic reduction in the whole derivation. This is exactly the
step the source material calls "lengthy" and does not show line by line --
here, it is not skipped, only carried out by machine instead of by hand:

```@eval
using Symbolics, Markdown
@variables l q Pi_0
@variables jl jlp Q C
bc1_eq = C*jl + Pi_0 ~ -1
bc2_eq = C*(-q^2*jl + 2*(l^2+l-1)*jl - 2*q*jlp) + 2*(l^2-1)*Pi_0 ~ 0
bc_solution = Symbolics.solve_for([bc1_eq, bc2_eq], [C, Pi_0])
Csol = simplify(substitute(bc_solution[1], Dict(jlp => jl*(l-q*Q)/q)); expand=true)
Pi0sol = simplify(substitute(bc_solution[2], Dict(jlp => jl*(l-q*Q)/q)); expand=true)
Uprime1_minus_2U1 = substitute(C*(jl + q*jlp) + Pi_0*(l+1) + 2, Dict(C => bc_solution[1], Pi_0 => bc_solution[2]))
Uprime1_minus_2U1 = simplify(substitute(Uprime1_minus_2U1, Dict(jlp => jl*(l-q*Q)/q)); expand=true)
alpha4_derived = simplify(q^2*(q^2*Pi0sol - 2*l*Uprime1_minus_2U1); expand=true)
characteristic_eq_rhs = q^4*((2*(l-1)/q^2)*(l + (l+1)*(q-2*l*Q)/(q-2*Q)) - 1)
residual = simplify(alpha4_derived - characteristic_eq_rhs; expand=true)
Markdown.parse("Substituting the boundary-condition solutions into the \$T_1\$-eliminated \$\\alpha^4\$ relation, the difference from Reid's closed-form right-hand side collapses to `" *
    string(residual) * "`, for symbolic \$l\$ and \$q\$:\n\n```math\n\\frac{\\alpha^4}{q^4} + 1 = \\frac{2(l-1)}{q^2}\\left[l + (l+1)\\frac{q-2lQ_{l+1/2}(q)}{q-2Q_{l+1/2}(q)}\\right]\n```")
```

where ``q^2=\sigma R^2/\nu``, ``\alpha^2=\sigma_{l;0}R^2/\nu``, and
``Q_{l+1/2}(q)=j_{l+1}(q)/j_l(q)``.

Every physically meaningful outcome of this whole derivation -- the
oscillation frequency, the damping rate, whether a drop of a given size
and viscosity oscillates or just squashes back down -- follows from
exactly this one equation.

!!! note "Universality"
    This is identical to Chandrasekhar's characteristic equation for a
    viscous self-gravitating liquid globe, with the self-gravitational
    parameter identified with ``\alpha^2``. The physical restoring force
    (surface tension vs. self-gravity) enters only through
    ``\sigma_{l;0}``, which defines ``\alpha``. The equation above holds
    for any spherically symmetric restoring force.
