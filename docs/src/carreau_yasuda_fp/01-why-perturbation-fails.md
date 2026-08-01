# Why Amplitude Perturbation Theory Fails

The Reid (1960) chapters give the exact Newtonian damping rate
``\lambda_l(\mathrm{Oh})`` and squared frequency ``\omega_l^2(\mathrm{Oh})``
at any Ohnesorge number. The CAS Derivations page
"Carreau-Yasuda, weakly-nonlinear" builds a small-amplitude correction on
top of that theory for the *classical* Carreau law, where the shape
exponent is fixed at ``a=2``. This chapter asks a harder question: does an
analogous correction exist for the *general* Carreau-Yasuda law this repo's
own fitted validation fluid actually uses, with
``a\approx 0.743`` -- not ``2``?

The answer, reached only after an extensive derivation independently
adversarially checked twice, is no, not for this fluid's fitted parameters.
That is not a shortcut or a failure of nerve -- it is a precisely quantified
fact about where this fluid's physics actually sits on the Carreau-Yasuda
curve, and every number below is asserted and checked in
`julia/derivations/carreauYasuda_firstprinciples_derivation.jl`.

## The exponent, not just its size, is the problem

The classical correction works because ``a=2`` is an even integer: the
viscosity correction ``(\lambda_c\dot\gamma)^2`` is then an analytic
(Taylor-expandable) function of the oscillation amplitude ``\epsilon``,
since ``\dot\gamma=O(\epsilon)`` linearly in Reid's single-mode ansatz. For
a general, non-integer exponent this breaks down completely:

```@eval
using Markdown
a = 0.7430524574330837
ratios = [abs(eps)^a / abs(eps) for eps in (0.1, 0.001, 1e-5, 1e-7)]
Markdown.parse("For ``a=$(round(a,digits=4))``, the ratio ``|\\epsilon|^a/|\\epsilon|`` at " *
    "``\\epsilon=10^{-1},10^{-3},10^{-5},10^{-7}`` is\n\n```math\n" *
    join(string.(round.(ratios, sigdigits=4)), ",\\quad ") *
    "\n```\ngrowing without bound as ``\\epsilon\\to 0``.")
```

There is no amplitude, however small, at which this correction is a small
perturbation of a fixed integer-power hierarchy: it sits at a genuinely
fractional order between Reid's own ``O(\epsilon)`` linear terms and the
usual ``O(\epsilon^2)`` geometric mode-coupling correction. This alone does
not doom a perturbative treatment -- a fractional expansion parameter is
unusual but not impossible. What actually closes off the small-amplitude
route for this fluid is quantitative, not structural.

## Two natural anchors

A shear-thinning fluid's viscosity interpolates between two Newtonian
limits: ``\eta_0`` at rest and ``\eta_\infty`` at infinite shear rate.
Reid's characteristic equation gives the exact ``\lambda_l(\mathrm{Oh})``,
``\omega_l^2(\mathrm{Oh})`` at either limit, so both are legitimate places
to try anchoring a correction.

!!! warning "Anchoring at rest, $\mathrm{Oh}_0$, fails outright"
    ```@eval
    using DropSolver, Markdown
    OH0 = 57.371648873370795
    rows = String[]
    for l in (2,3,5,10)
        lam, om2, _ = reid_lambda_omega2(OH0, l)
        push!(rows, "| $l | $(round(lam,digits=1)) | $(round(om2,digits=2)) | overdamped |")
    end
    Markdown.parse("At \$\\mathrm{Oh}_0\\approx $(round(OH0,digits=2))\$, every mode this " *
        "solver resolves is heavily overdamped:\n\n" *
        "| \$l\$ | \$\\lambda_l\$ | \$\\omega_l^2\$ | |\n|:--|:--|:--|:--|\n" *
        join(rows, "\n") *
        "\n\nThere is no oscillation to average over, so a multiple-scales/envelope " *
        "treatment has no fast carrier wave to modulate. Separately, the shear rate " *
        "at which the viscosity correction reaches even the 1% level is " *
        "\$\\dot\\gamma\\lesssim 7\\times 10^{-8}\$ -- far below any resolved mode velocity.")
    ```

!!! warning "Anchoring at infinite shear, $\mathrm{Oh}_\infty$, is better posed but still not small"
    ```@eval
    using DropSolver, Markdown
    OH0 = 57.371648873370795
    EPS_ST = 0.9995574839318364
    OH_INF = OH0*(1-EPS_ST)
    lam, om2, _ = reid_lambda_omega2(OH_INF, 2)
    Q = sqrt(om2)/(2lam)
    Markdown.parse("Define \$\\mathrm{Oh}_\\infty\\equiv\\mathrm{Oh}_0(1-\\varepsilon_{ST})" *
        "\\approx $(round(OH_INF,sigdigits=4))\$. Unlike \$\\mathrm{Oh}_0\$, this IS " *
        "comfortably underdamped: \$\\lambda_2=$(round(lam,digits=4))\$, " *
        "\$\\omega_2^2=$(round(om2,digits=3))\$, giving a quality factor " *
        "(Reid's oscillator convention is \$\\ddot A+2\\lambda\\dot A+\\omega^2A=\\text{forcing}\$, " *
        "so \$Q\\equiv\\omega/(2\\lambda)\$) of \$Q\\approx $(round(Q,digits=1))\$ -- a genuine " *
        "separation between the decay time and the oscillation period.")
    ```
    A live solver trace (three impact energies, fine time resolution) confirms
    the drop's actual effective Ohnesorge number, once excited, is permanently
    pinned away from ``\mathrm{Oh}_0`` after a single clean transient at first
    contact, with no recurrence on later contact/lift-off cycles.

    But "pinned away from ``\mathrm{Oh}_0``" is not "close to
    ``\mathrm{Oh}_\infty``". The relevant question is whether the
    *first-order Taylor expansion* of Reid's exact ``\lambda_l(\mathrm{Oh})``,
    ``\omega_l^2(\mathrm{Oh})`` around ``\mathrm{Oh}_\infty`` reproduces
    their true values over the range this fluid's dynamics actually visits.
    Checked directly:

    ```@eval
    using DropSolver, Markdown
    OH0 = 57.371648873370795
    EPS_ST = 0.9995574839318364
    OH_INF = OH0*(1-EPS_ST)
    h = 1e-6
    rows = String[]
    for l in (2,3,5,10)
        lam0, om0, _ = reid_lambda_omega2(OH_INF, l)
        lamp, omp, _ = reid_lambda_omega2(OH_INF+h, l)
        lamm, omm, _ = reid_lambda_omega2(OH_INF-h, l)
        dlam = (lamp-lamm)/(2h); dom = (omp-omm)/(2h)
        lam_lo, om_lo, _ = reid_lambda_omega2(0.03, l)
        lam_hi, om_hi, _ = reid_lambda_omega2(0.1, l)
        lam_lin_lo = lam0 + dlam*(0.03-OH_INF)
        lam_lin_hi = lam0 + dlam*(0.1-OH_INF)
        om_lin_hi = om0 + dom*(0.1-OH_INF)
        err_lam_lo = abs(lam_lin_lo-lam_lo)/abs(lam_lo)*100
        err_lam_hi = abs(lam_lin_hi-lam_hi)/abs(lam_hi)*100
        err_om_hi = abs(om_lin_hi-om_hi)/abs(om_hi)*100
        push!(rows, "| $l | $(round(err_lam_lo,digits=2))% | $(round(err_lam_hi,digits=1))% | $(round(err_om_hi,digits=2))% |")
    end
    Markdown.parse("| \$l\$ | \$\\lambda\$ err at Oh=0.03 | \$\\lambda\$ err at Oh=0.1 | \$\\omega^2\$ err at Oh=0.1 |\n" *
        "|:--|:--|:--|:--|\n" * join(rows, "\n"))
    ```

    ``\omega_l^2`` is nearly flat in ``\mathrm{Oh}`` near ``\mathrm{Oh}_\infty``
    -- the linear correction stays small even for higher ``l``. But
    ``\lambda_l`` is not: at the low end of the operative band every ``l``
    is comfortably small, but by the high end the error has grown past
    ``25\%`` for ``l=10`` -- the same growing-with-``l`` curvature that
    broke the ``\mathrm{Oh}_0`` anchor reappears here, and here it also
    grows *across* the operative band, not just across modes. Moving the
    anchor to some other "typical" ``\mathrm{Oh}^\ast`` does not rescue
    this: the curvature at high ``l`` sets in almost immediately away from
    any point in the operative range.

!!! note "The honest conclusion"
    For this fluid's fitted parameters, there is no single reference
    viscosity -- rest, infinite shear, or anywhere in between -- around
    which a regular perturbation expansion of Reid's damping rate is
    uniformly valid across the modes this solver resolves. The physics
    genuinely lives in the fully nonlinear, transitional part of the
    Carreau-Yasuda curve. This is a quantitative fact about *this fluid's*
    ``\lambda_c`` and ``a``, not a universal statement about shear-thinning
    drops.

Because this failure is specifically a failure of *linearizing*
``\lambda_l(\mathrm{Oh})``, ``\omega_l^2(\mathrm{Oh})`` -- and these
functions are already known exactly, for any ``\mathrm{Oh}``, via Reid's
characteristic equation -- the correct response is not to build a better
linear correction. It is to not linearize this piece at all: evaluate the
exact relations directly at whatever effective Ohnesorge number a closure
provides, exactly as this repo's existing (heuristic) Carreau-Yasuda
extension already does. That part was never the weak point. The weak point
is upstream of it, and the next chapter is about exactly that.
