# What Survives: Three General Results, and the Open Question

The previous chapter's negative finding does not mean the investigation
produced nothing. Three genuinely new, general results survive
independently of any fluid's specific parameters -- they hold for *any*
Carreau-Yasuda fluid, any Oh, any mode -- and are exactly the tools a
future, genuinely self-consistent treatment would need.

## Machinery I: adjoint sensitivity for Reid's boundary-value problem

Reid's velocity operator ``\mathcal{L}[U]\equiv U''-l(l+1)U/x^2+q^2U`` has
no first-derivative term, so it is formally self-adjoint on ``(0,1)`` with
the plain ``L^2`` inner product -- no weight function needed:

```math
W\,\mathcal{L}[U] - U\,\mathcal{L}[W] = \frac{d}{dx}\Bigl[\,W\,U' - W'\,U\,\Bigr],
```

verified symbolically for abstract ``U(x)``, ``W(x)``. This gives a
shortcut for exactly the calculation a perturbation of Reid's problem
needs -- the boundary derivative of a particular solution, without solving
the ODE:

!!! tip "Adjoint shortcut"
    Let ``Y(x)`` be regular at ``x=0``, satisfy ``Y(1)=0``, and solve
    ``\mathcal{L}[Y]=\mathrm{RHS}(x)`` at a fixed wavenumber ``q_0``. Then

    ```math
    Y'(1) = \frac{1}{j_l(q_0)}\int_0^1 x\,j_l(q_0x)\,\mathrm{RHS}(x)\,dx,
    ```

    using the regular homogeneous solution ``x\,j_l(q_0x)`` as the adjoint
    test function.

Checked against direct numerical shooting (regular-solution seeding near
``x=0``, RK4 integration, homogeneous-solution correction to enforce
``Y(1)=0``) for three independent ``(l,q_0,\mathrm{RHS})`` combinations:
agreement better than ``10^{-6}`` relative error. Combined with Reid's BC1
(``U(1)=-1``, unaffected by any viscosity correction since it is purely
kinematic) and BC2, this reduces the search for the ``O(\delta)`` shift in
``q^2`` from any new forcing added to Reid's ODE and both stress boundary
conditions to an ordinary linear solve -- no further ODE integration is
required.

## Machinery II: the strain-rate field of Reid's actual viscous mode

Any first-principles shear-rate calculation must use Reid's actual
*viscous* velocity profile ``U(x)=Cx\,j_l(qx)+\Pi_0x^{l+1}``, not the
inviscid potential-flow shape ``r^lP_l(\cos\theta)`` the existing heuristic
Carreau-Yasuda scripts use. The two are not interchangeable: Reid's own
damping normalization comes from the homogeneous (Bessel) part of
``U(x)``, which *is* the viscous correction to potential flow.

Writing ``u_r=F(x)P_l(\cos\theta)``, ``F(x)\equiv U(x)/x^2``, and
``u_\theta`` from the stream-function relation already used for BC2,

```math
e_{rr}=F'(x)P_l(\mu),\quad
e_{\theta\theta}=\frac{1}{x}\partial_\theta u_\theta+\frac{u_r}{x},\quad
e_{\varphi\varphi}=\frac{u_r+u_\theta\cot\theta}{x},\quad
e_{r\theta}=\frac12\!\left[x\partial_x\!\left(\frac{u_\theta}{x}\right)+\frac{1}{x}\partial_\theta u_r\right],
```

with ``\mu\equiv\cos\theta``. Incompressibility (``e_{rr}+e_{\theta\theta}+e_{\varphi\varphi}=0``)
was verified to floating-point precision for ``l=2,\ldots,8``, for *any*
radial profile ``F(x)`` -- a kinematic identity of the poloidal
representation, confirming the construction before any Carreau-Yasuda
physics enters.

## Machinery III: the period-$\pi$ lemma

At ``\mathrm{Oh}_\infty``, Reid's root ``q`` is genuinely complex (a true
damped oscillation, not a pure decay). A natural first attempt at a
shear-thinning correction -- evaluate everything at one representative
oscillation phase -- turns out to be indefensible: the real physical shape
at a representative radius swings by roughly the same magnitude in both
directions across one period, nowhere near constant. Any nonlinear
(especially fractional-power) function of this field depends on which
phase was chosen, unless the choice is justified by an actual time-average.

!!! note "The shear-rate invariant cannot resonate at the mode's own frequency"
    Let ``S\equiv\sqrt{2e_{ij}e_{ij}}`` be built from any single-frequency
    oscillating axisymmetric poloidal field. Then ``S(\phi)`` (``\phi=\omega t``)
    is exactly periodic with period ``\pi``, not ``2\pi``, so its Fourier
    series over the full period contains only even harmonics: the
    coefficient at the fundamental (``m=1``) is identically zero, for
    ``S`` itself and for *any* function of ``S``.

```@eval
using Markdown
Markdown.parse("**Proof.** Each strain component satisfies " *
    "\$\\mathrm{Re}(e_{ij}e^{-i\\phi})=|e_{ij}|\\cos(\\phi-\\arg e_{ij})\$. Squaring, " *
    "\$|e_{ij}|^2\\cos^2(\\phi-\\arg e_{ij})=\\tfrac12|e_{ij}|^2[1+\\cos(2\\phi-2\\arg e_{ij})]\$, " *
    "manifestly invariant under \$\\phi\\to\\phi+\\pi\$. \$S^2\$ is a sum of such terms, hence " *
    "itself period-\$\\pi\$; \$S=\\sqrt{S^2}\\ge0\$ and any real function of \$S\$ inherit the " *
    "period exactly (no branch ambiguity, since \$S^2\\ge0\$ throughout). A period-\$\\pi\$ " *
    "function integrated against \$\\sin(m\\phi)\$ or \$\\cos(m\\phi)\$ over \$[0,2\\pi)\$ vanishes " *
    "identically for every odd \$m\$, since the two half-period copies enter with opposite sign " *
    "and cancel.")
```

Confirmed two ways: numerically, ``S^2(\phi+\pi)=S^2(\phi)`` to
``10^{-15}``--``10^{-16}`` at multiple independent ``(x,\theta)`` points; and
via direct Fourier projection on a joint ``800\times800`` ``(\theta,\phi)``
grid, where the ``m=1`` coefficient (both ``\cos`` and ``\sin``) is zero to
``10^{-17}``--``10^{-18}`` at every spherical-harmonic degree tested.

This has an immediate consequence: since no forcing resonant at the base
mode's own frequency exists, the adjoint machinery above -- built for
exactly this kind of resonant solvability condition -- has nothing to act
on at ``m=1``. The leading temporal effect of any generalized-Newtonian
correction is carried entirely by the ``m=0`` (period-averaged, effectively
steady) channel; the ``m=2`` channel is a distinct, second-harmonic
*parametric* coupling, a different physical phenomenon not addressed here.

## The genuinely open question: spatial homogenization

Because the ``m=0`` channel is the leading one, and because a steady
forcing is a valid input to the adjoint shortcut above, that machinery is
the right tool for the question that actually remains open: replacing a
spatially varying ``\eta(x,\theta)`` with a single scalar
``\mathrm{Oh}_{\mathrm{eff}}(t)`` is an uncontrolled effective-medium
approximation unless the spatial variation is itself small.

```@eval
using DropSolver, SpecialFunctions, Markdown
OH0 = 57.371648873370795
LAMBDA_C = 30507.34501244818
A_SHAPE = 0.7430524574330837
EPS_ST = 0.9995574839318364
OH_INF = OH0*(1-EPS_ST)
l = 2
q0 = dominant_root(OH_INF, l)
sph_jl(l,z) = sqrt(pi/(2z)) * besselj(l+0.5, z)
Q0 = sph_jl(l+1,q0)/sph_jl(l,q0)
jl0 = sph_jl(l,q0)
C0 = 2*(l^2-1) / (jl0*q0*(2*Q0-q0))
Pi0 = -1 - C0*jl0
Ufun_c(x) = C0*x*sph_jl(l, q0*x) + Pi0*x^(l+1)
Ffun_c(x) = Ufun_c(x)/x^2
Fderiv_c(x; h=1e-6) = (Ffun_c(x+h)-Ffun_c(x-h))/(2h)
function legendre_P(l::Int, xv)
    l == 0 && return one(xv); l == 1 && return xv
    Pm1, P = one(xv), xv
    for n in 1:(l - 1); P, Pm1 = ((2n + 1) * xv * P - n * Pm1) / (n + 1), P; end
    P
end
legendre_Pp(l::Int, xv) = l==0 ? zero(xv) : l*(xv*legendre_P(l,xv) - legendre_P(l-1,xv))/(xv^2-1)
function strain_complex(xv, th, l)
    mu = cos(th); F = Ffun_c(xv); Fp = Fderiv_c(xv); dth = 1e-6
    uth_coef(xx) = -(2*Ffun_c(xx) + xx*Fderiv_c(xx)) / (l*(l+1))
    Pl_th(t)  = legendre_P(l, cos(t)); Plp_th(t) = legendre_Pp(l, cos(t))
    ur_of_theta(t)  = F*Pl_th(t)
    uth_of_theta(t) = uth_coef(xv)*sin(t)*Plp_th(t)
    dur_dth  = (ur_of_theta(th+dth)  - ur_of_theta(th-dth))  / (2dth)
    duth_dth = (uth_of_theta(th+dth) - uth_of_theta(th-dth)) / (2dth)
    ur  = F*Pl_th(th); uth = uth_coef(xv)*sin(th)*Plp_th(th)
    e_rr   = Fp*legendre_P(l,mu)
    e_thth = duth_dth/xv + ur/xv
    e_phph = (ur + uth*cos(th)/sin(th))/xv
    dxv = 1e-6
    uth_x(xx) = -(2*Ffun_c(xx)+xx*Fderiv_c(xx))/(l*(l+1)) * sin(th)*Plp_th(th) / xx
    duthx_dx = (uth_x(xv+dxv)-uth_x(xv-dxv))/(2dxv)
    e_rth  = 0.5*(xv*duthx_dx + dur_dth/xv)
    e_rr, e_thth, e_phph, e_rth
end
function S_at(xv, th, phase, l)
    e_rr, e_thth, e_phph, e_rth = strain_complex(xv, th, l)
    ph = cos(phase) - im*sin(phase)
    err_r,ethth_r,ephph_r,erth_r = real(e_rr*ph),real(e_thth*ph),real(e_phph*ph),real(e_rth*ph)
    sqrt(2*(err_r^2+ethth_r^2+ephph_r^2+2*erth_r^2))
end
vals = Float64[]
for xv in range(0.15,0.98;length=8), th in range(0.1,pi-0.1;length=8)
    Nph = 60
    phases = collect(range(0.0,2pi;length=Nph+1))[1:end-1]
    push!(vals, sum(S_at(xv,th,ph,l) for ph in phases)/Nph)
end
spread = maximum(vals)/minimum(vals)
visc_spread = spread^(-A_SHAPE)
Markdown.parse("Using the period-averaged shear-rate shape from Machinery II-III over the " *
    "drop's volume at \$\\mathrm{Oh}_\\infty\$, \$l=2\$: the shape varies by a factor of " *
    "\$\\approx $(round(spread,digits=2))\$ between its minimum and maximum. Propagated through " *
    "the Carreau-Yasuda exponent \$a=$(round(A_SHAPE,digits=3))\$, this is roughly a factor of " *
    "\$$(round(1/visc_spread,digits=2))\$ spread in local viscosity across the drop at a single " *
    "instant -- smaller than the multiple-order-of-magnitude range \$\\mathrm{Oh}_{\\mathrm{eff}}\$ " *
    "traverses over an impact *in time*, but not small enough to treat the scalar-\$\\mathrm{Oh}_" *
    "{\\mathrm{eff}}\$ approximation as an uncontrolled-but-negligible simplification.")
```

Legendre-projecting the period-averaged correction shape confirms a further
mechanism: a *single* active mode (``l=2``) spontaneously forces every even
spherical-harmonic degree (``l'=0,4,6,8`` all genuinely nonzero, odd ``l'``
zero to ``10^{-11}`` by a verified parity argument), because raising a
finite Legendre polynomial to a non-integer power does not, in general,
stay within its own degree.

## Summary

1. The classical ``O(\epsilon^3)`` Carreau expansion does not extend to
   this fluid's fitted, non-integer exponent: no small-amplitude anchor
   supports a uniformly valid linear correction to ``\lambda_l(\mathrm{Oh})``
   across the modes this solver resolves.
2. Reid's exact ``\lambda_l(\mathrm{Oh})``, ``\omega_l^2(\mathrm{Oh})``
   should be evaluated directly at whatever ``\mathrm{Oh}_{\mathrm{eff}}(t)``
   a closure provides -- never linearized -- since this repo already
   implements a fast, tabulated version of exactly that evaluation.
3. The genuinely open approximation in the existing heuristic scheme is
   the collapse of a spatially varying ``\eta(x,\theta,t)`` onto a single
   scalar ``\mathrm{Oh}_{\mathrm{eff}}(t)``, quantified here at roughly a
   factor of two in local viscosity across the drop at a representative
   instant -- real, moderate, but not yet controlled.
4. The adjoint sensitivity shortcut, the viscous strain-rate tensor, and
   the period-``\pi`` lemma survive independently of any of the above.
   None depend on this fluid's specific parameters, and all three are the
   correct starting point for a future, genuinely first-principles
   treatment -- most plausibly a self-consistent, spatially resolved solve
   using the adjoint machinery as a Newton/Fréchet-derivative iteration
   step, rather than any further attempt at a closed-form amplitude
   equation.

This chapter does not implement that self-consistent solve -- that is
future work. It also does not revisit the ``m=2`` parametric channel, nor
extend the strain-tensor/leakage calculation beyond ``l=2`` and a single
representative radius. What it establishes, to the same CAS-verified
standard as the rest of this repo's documentation -- every claim above has
a corresponding assertion in
`julia/derivations/carreauYasuda_firstprinciples_derivation.jl`, including
a live cross-check against a running `solve_drop!` trace -- is precisely
which of the natural next steps work, which do not, and why, for this
fluid specifically.
