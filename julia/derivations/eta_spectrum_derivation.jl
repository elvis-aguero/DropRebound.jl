# # The Angular Bandwidth of a Shear-Thinning Viscosity Field
#
# A generalised-Newtonian drop carries a viscosity that is a *field*,
# ``\eta(r,\theta,t)``, because the local shear rate is a field. Expanding that
# field in Legendre polynomials,
#
# ```math
# \eta(r,\mu,t) \;=\; \sum_{l'\ge 0} \eta_{l'}(r,t)\,P_{l'}(\mu),
# \qquad \mu=\cos\theta,
# ```
#
# turns the extra momentum term ``2(\nabla\eta)\cdot\bm e`` into a matrix that
# couples the drop's shape modes to one another. The selection rule for
# Legendre triple products means the harmonic ``\eta_{l'}`` connects shape modes
# ``l`` and ``l''`` only when ``|l-l''|\le l'\le l+l''``. So if ``\eta_{l'}``
# were negligible above some cutoff ``L_\eta``, the coupling matrix would be
# **banded** with half-bandwidth ``L_\eta``, and applying it would cost
# ``O(M L_\eta)`` instead of ``O(M^2)`` for a drop truncated at ``M`` shape
# modes.
#
# This page measures ``L_\eta``. The answer decides whether banding is worth
# implementing, and it also prices the cruder closure that throws away *all*
# angular structure and keeps only ``\eta_0(r)``. Those two closures are Steps 6
# and 7 of [A Hierarchy of Models](generalized_newtonian_hierarchy_derivation.md);
# this page is the evidence behind the numbers quoted there.
#
# The measurement has three parts, and each is a place the answer could go
# wrong, so each is verified separately:
#
# 1. **the velocity field.** ``\eta`` is a function of the strain rate, so the
#    strain-rate tensor of Reid's viscous eigenmodes has to be evaluated
#    accurately at mode numbers up to ``l\sim50``. This is harder than it looks
#    and Section 1 is entirely about doing it without catastrophic cancellation.
# 2. **the shear-rate invariant.** The eigenmodes oscillate; ``\eta`` responds
#    instantaneously but the mode-coupling matrix wanted by the solver is a
#    period-averaged object. Section 2 does that average in closed form up to a
#    single one-dimensional quadrature.
# 3. **the decomposition itself,** including the choice of error norm, which
#    turns out to change the answer by more than an order of magnitude and is
#    the single most misleading step. Section 3.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``R`` | equilibrium drop radius; ``x=r/R\in(0,1]`` |
# | ``l``, ``l''`` | degree of a **shape** mode; runs ``2\ldots M`` |
# | ``M`` | truncation of the shape expansion |
# | ``l'`` | degree of a harmonic of the **viscosity field**; runs ``0,1,2,\ldots`` |
# | ``L_\eta`` | largest ``l'`` retained: the half-bandwidth of the coupling |
# | ``\mathrm{Oh}`` | Ohnesorge number of the mode being evaluated |
# | ``q`` | Reid's radial wavenumber; ``q=q(\mathrm{Oh},l)``, complex |
# | ``U(x)`` | Reid's poloidal amplitude |
# | ``F,G`` | radial and tangential velocity amplitudes |
# | ``\bm e`` | strain-rate tensor |
# | ``S=\sqrt{2\,\bm e\!:\!\bm e}`` | shear-rate invariant |
# | ``\bar S`` | its average over one oscillation period |
# | ``\dot A_l`` | modal surface velocity; the drop's state vector |
# | ``\lambda_c,\ m`` | Cross time constant and exponent |
# | ``\eta_{l'}`` | ``l'``-th Legendre coefficient of ``\eta`` |
# | ``\mathcal T_1(L)`` | ``\sum_{l'>L}|\eta_{l'}|/|\eta_0|`` -- summed discarded coupling |
# | ``\mathcal T_2(L)`` | ``\bigl(\sum_{l'>L}\eta_{l'}^2\bigr)^{1/2}/|\eta_0|`` -- its RMS counterpart |
#
# ``l'`` is not a shape mode and never enters the *size* of the coupling matrix,
# which is always ``M\times M``. It controls only how far off the diagonal that
# matrix reaches.

using Printf, SpecialFunctions  #src
using DropSolver: dominant_root, sph_bessel_ratio, reid_char, gauss_legendre_nodes  #src

## The default run is sized for continuous integration. Setting  #src
## ETA_SPECTRUM_FULL=1 restores the production sweep documented in Section 5.2.  #src
const FULL = get(ENV, "ETA_SPECTRUM_FULL", "0") == "1"  #src

println("="^92)  #src
println("THE ANGULAR BANDWIDTH L_eta OF A SHEAR-THINNING VISCOSITY FIELD")  #src
println(FULL ? "mode: FULL production sweep (ETA_SPECTRUM_FULL=1)" :  #src
               "mode: reduced CI run (set ETA_SPECTRUM_FULL=1 for the full sweep)")  #src
println("="^92)  #src

# ## 1. The strain rate of a single Reid eigenmode
#
# ### 1.1 Why the textbook expression cannot be typed in literally
#
# Reid's mode-``l`` poloidal amplitude, derived in
# [The Free Viscous Drop](reid1960_full_derivation.md), is
#
# ```math
# U(x) \;=\; C\,x\,j_l(qx) \;+\; \Pi_0\,x^{\,l+1},
# ```
#
# where ``j_l`` is the spherical Bessel function of the first kind,
# ``q=q(\mathrm{Oh},l)`` is the dominant root of Reid's characteristic equation,
# and the two constants follow from the boundary conditions,
#
# ```math
# C \;=\; \frac{2(l^2-1)}{j_l(q)\,q\,(2Q-q)},
# \qquad Q=\frac{j_{l+1}(q)}{j_l(q)},
# \qquad C\,j_l(q)+\Pi_0=-1 .
# ```
#
# Written that way the expression is a trap. Reid's root is complex, and at
# large ``l`` and small ``\mathrm{Oh}`` it is *very* complex: at ``l=50``,
# ``\mathrm{Oh}=0.025`` one finds ``q\approx 91.6-70.8\,\mathrm{i}``, so
# ``j_l(q)\sim 2.5\times10^{25}`` while ``C\sim10^{-25}``. Their product is
# perfectly ordinary; the factors are not. Pushing to smaller ``\mathrm{Oh}``
# the failure is progressive rather than abrupt: at ``\mathrm{Oh}=10^{-3}``,
# ``|\operatorname{Im}q|\approx420`` and ``j_l(q)\sim10^{178}`` -- large enough
# that the field returned is accurate to only a handful of digits, while nothing
# in the output announces it -- and at ``\mathrm{Oh}=3\times10^{-4}`` the Bessel
# routine simply refuses.
#
# The fix is to never form ``C`` or ``j_l(q)`` separately. Define the
# **scale-free** constant and the **scale-free** Bessel ratio
#
# ```math
# C_j \;:=\; C\,j_l(q) \;=\; \frac{2(l^2-1)}{q\,(2Q-q)},
# \qquad \Pi_0=-1-C_j,
# \qquad W(x) \;:=\; \frac{j_l(qx)}{j_l(q)},
# ```
#
# so that ``W(1)=1`` and
#
# ```math
# \boxed{\;U(x)=C_j\,x\,W(x)+\Pi_0\,x^{\,l+1},\qquad U(1)=C_j+\Pi_0=-1\;}
# ```
#
# **In plain English:** every quantity in that line is of order one at the drop
# surface and decays monotonically inward. The enormous numbers were an artefact
# of splitting a well-behaved product into two badly-behaved factors, and the
# boundary condition ``U(1)=-1`` is now satisfied *identically*, by construction,
# rather than as the result of a near-cancellation.
#
# ``Q`` comes from `DropSolver.sph_bessel_ratio`, which propagates the ratio
# itself by downward recurrence rather than the two Bessel functions; the same
# idiom, one level up. ``W`` is built from the *exponentially scaled* cylindrical
# Bessel function ``\tilde J_\nu(z)=J_\nu(z)e^{-|\operatorname{Im}z|}``. Using
# ``j_l(z)=\sqrt{\pi/2z}\,J_{l+1/2}(z)`` and noting that for real ``x>0`` one has
# ``|\operatorname{Im}(qx)|=x|\operatorname{Im}q|``,
#
# ```math
# W(x)=x^{-1/2}\;
#      \frac{\tilde J_{l+1/2}(qx)}{\tilde J_{l+1/2}(q)}\;
#      e^{|\operatorname{Im}q|\,(x-1)} .
# ```
#
# The exponential factor never exceeds one. Deep inside the drop it underflows
# smoothly to zero, which is the physically correct answer there: a high-``l``
# viscous eigenmode is a boundary layer near the surface, not a bulk motion.

const XMIN = 1e-6  #src

struct ReidMode  #src
    l    :: Int  #src
    Oh   :: Float64  #src
    q    :: ComplexF64  #src
    Cj   :: ComplexF64      # C * j_l(q); finite where C and j_l(q) separately are not  #src
    Pi0  :: ComplexF64  #src
    jden :: ComplexF64      # the scaled denominator Jtilde_{l+1/2}(q)  #src
    aiq  :: Float64         # |Im q|  #src
end  #src

function reid_mode(Oh::Real, l::Integer)  #src
    q  = ComplexF64(dominant_root(Float64(Oh), Int(l)))  #src
    Q  = ComplexF64(sph_bessel_ratio(Int(l), q))  #src
    Cj = 2 * (l^2 - 1) / (q * (2Q - q))  #src
    Pi0 = -1 - Cj  #src
    jden = besseljx(l + 0.5, q)  #src
    (isfinite(abs(jden)) && abs(jden) > 0) ||  #src
        error("scaled Bessel denominator $jden unusable at Oh=$Oh, l=$l")  #src
    (isfinite(abs(Cj)) && isfinite(abs(Pi0))) ||  #src
        error("Cj/Pi0 not finite at Oh=$Oh, l=$l")  #src
    ReidMode(Int(l), Float64(Oh), q, Cj, Pi0, jden, abs(imag(q)))  #src
end  #src

## W(x) = j_l(qx)/j_l(q) and V(x) = j_{l-1}(qx)/j_l(q). Complex x is accepted  #src
## deliberately: it is what makes the Cauchy-integral derivative test possible.  #src
@inline function bessel_WV(m::ReidMode, x::Number)  #src
    z  = m.q * x  #src
    sc = exp(abs(imag(z)) - m.aiq) / sqrt(complex(x))  #src
    W  = besseljx(m.l + 0.5, z) / m.jden * sc  #src
    V  = besseljx(m.l - 0.5, z) / m.jden * sc  #src
    return W, V  #src
end  #src

# ### 1.2 Derivatives without finite differences
#
# The strain rate needs ``U'`` and ``U''``. Differencing them would reintroduce
# a truncation floor exactly where the field is stiffest, so both are obtained
# in closed form.
#
# Start from the spherical Bessel derivative recurrence
# ``\tfrac{d}{dz}j_l(z)=j_{l-1}(z)-\tfrac{l+1}{z}j_l(z)``. Divide by ``j_l(q)``,
# set ``z=qx`` and use ``d/dx=q\,d/dz``:
#
# ```math
# W'(x)=q\,V(x)-\frac{l+1}{x}\,W(x),
# \qquad V(x):=\frac{j_{l-1}(qx)}{j_l(q)} .
# ```
#
# The definition of ``V`` deserves a note, because the obvious alternative is a
# bug. Factoring ``W`` out of the recurrence would give the ratio
# ``j_{l-1}(qx)/j_l(qx)``, which is **singular at every zero of ``j_l(qx)``** --
# and those zeros are real and inside the drop whenever ``q`` is real, which
# happens at large ``\mathrm{Oh}``. Keeping the denominator pinned at ``j_l(q)``
# makes ``V`` a finite quantity everywhere, evaluated by exactly the same scaled
# route as ``W``, one order down.
#
# Substituting ``xW'=qxV-(l+1)W`` into ``U'=C_j(W+xW')+\Pi_0(l+1)x^l`` makes the
# ``(l+1)W`` terms partially cancel *analytically*:
#
# ```math
# \boxed{\;U'(x)=C_j\bigl[q\,x\,V(x)-l\,W(x)\bigr]+\Pi_0(l+1)x^{\,l}\;}
# \tag{A}
# ```
#
# For the second derivative, differentiating (A) again would require
# ``j_{l-2}``. It is cheaper and better conditioned to use the spherical Bessel
# **equation** instead, ``j_l''=-\tfrac{2}{z}j_l'-\bigl(1-\tfrac{l(l+1)}{z^2}\bigr)j_l``,
# which in the ``W`` variable says ``2W'+xW''=-\bigl(q^2x-\tfrac{l(l+1)}{x}\bigr)W``,
# hence
#
# ```math
# \boxed{\;U''(x)=C_j\Bigl[\frac{l(l+1)}{x}-q^2x\Bigr]W(x)+\Pi_0\,l(l+1)\,x^{\,l-1}\;}
# \tag{B}
# ```
#
# **In plain English:** (B) says the same thing as Reid's own radial equation
# ``U''=l(l+1)U/x^2-q^2U+q^2\Pi_0x^{l+1}`` -- expand ``U`` in (B) and the two
# agree term by term. That is not a coincidence to be admired but a check to be
# run, because it exercises the Bessel algebra above independently of the field
# it produces: if the recurrence manipulation were wrong, the two forms would
# disagree.

@inline function profile(m::ReidMode, x::Number)  #src
    l = m.l  #src
    W, V = bessel_WV(m, x)  #src
    xl = x^l  #src
    U   = m.Cj * x * W + m.Pi0 * x * xl  #src
    Up  = m.Cj * (m.q * x * V - l * W) + m.Pi0 * (l + 1) * xl  #src
    Upp = m.Cj * (l * (l + 1) / x - m.q^2 * x) * W + m.Pi0 * l * (l + 1) * xl / x  #src
    return U, Up, Upp  #src
end  #src

# ### 1.3 From the amplitude to the strain-rate tensor
#
# Reid's convention puts the radial velocity amplitude at ``F(x)=U(x)/x^2`` (the
# ``x^{-2}`` is what makes the radial equation Bessel-shaped), and
# incompressibility then fixes the tangential amplitude,
# ``G=-(2F+xF')/[l(l+1)]``, so that with ``P=P_l(\mu)`` and ``P'=dP_l/d\mu``,
#
# ```math
# u_r=F(x)\,P,\qquad u_\theta=G(x)\,\sin\theta\;P' .
# ```
#
# Two exact reductions are used before anything is evaluated. Expanding
# ``2F+xF'=2U/x^2+x(U'/x^2-2U/x^3)`` the ``U/x^2`` pieces cancel identically, so
#
# ```math
# G(x)=-\frac{U'(x)}{l(l+1)\,x},
# \qquad
# G'(x)=\frac{U'(x)-x\,U''(x)}{l(l+1)\,x^{2}} .
# ```
#
# This matters numerically, not just cosmetically. The literal form subtracts
# ``2U/x^2`` from itself, so its relative accuracy degrades like
# ``|U/x^2|\big/|U'/x|`` wherever ``U'`` happens to be small -- which is
# wherever the mode has an internal node. The reduced form contains no
# subtraction at all.
#
# The strain-rate components in spherical coordinates, for an axisymmetric field
# and with ``R=1`` as the length unit, are then
#
# ```math
# \begin{aligned}
# e_{rr}       &= F'(x)\,P, \\
# e_{\theta\theta} &= \frac{G(x)\bigl[l(l+1)P-\mu P'\bigr]}{x}+\frac{F(x)P}{x}, \\
# e_{\phi\phi} &= \frac{F(x)P+G(x)\,\mu P'}{x}, \\
# e_{r\theta}  &= \tfrac12\sin\theta\,P'\Bigl[G'(x)-\frac{G(x)}{x}-\frac{F(x)}{x}\Bigr].
# \end{aligned}
# ```
#
# The ``e_{\theta\theta}`` line has already had ``P''`` eliminated, and it is
# worth showing how rather than asserting it. From
# ``e_{\theta\theta}=r^{-1}\partial_\theta u_\theta+u_r/r`` one needs
# ``\partial_\theta[\sin\theta\,P'(\cos\theta)]=\mu P'-(1-\mu^2)P''``, and
# Legendre's equation ``(1-\mu^2)P''=2\mu P'-l(l+1)P`` turns that into
# ``l(l+1)P-\mu P'``. No second derivative of a Legendre polynomial is ever
# evaluated in the final expressions, which removes the least stable of the
# recurrences.

@inline function radial_amps(m::ReidMode, x::Number)  #src
    l = m.l  #src
    U, Up, Upp = profile(m, x)  #src
    x2 = x * x  #src
    F  = U / x2  #src
    Fp = Up / x2 - 2U / (x2 * x)  #src
    ll = l * (l + 1)  #src
    G  = -Up / (ll * x)  #src
    Gp = (Up - x * Upp) / (ll * x2)  #src
    return F, Fp, G, Gp  #src
end  #src

## P_l, dP_l/dmu, d2P_l/dmu2 by the standard recurrences. P'' is retained only  #src
## so that Legendre's equation can be checked; it is not used in `strain`.  #src
function legendre_PdP(l::Integer, mu::Real)  #src
    Pm1 = 1.0  #src
    P   = float(mu)  #src
    l == 0 && return 1.0, 0.0, 0.0  #src
    for n in 1:(l-1)  #src
        Pm1, P = P, ((2n + 1) * mu * P - n * Pm1) / (n + 1)  #src
    end  #src
    d = 1 - mu * mu  #src
    if abs(d) < 1e-13  #src
        s  = mu > 0 ? 1.0 : -1.0  #src
        return s^l, s^(l + 1) * l * (l + 1) / 2,  #src
               s^l * (l - 1) * l * (l + 1) * (l + 2) / 8  #src
    end  #src
    Pp  = l * (Pm1 - mu * P) / d  #src
    Ppp = (2 * mu * Pp - l * (l + 1) * P) / d  #src
    return P, Pp, Ppp  #src
end  #src

## The four components, evaluated LITERALLY -- no term pre-cancelled. That is  #src
## deliberate: it keeps the incompressibility residual an honest measurement of  #src
## how much relative precision F, F' and G still carry, rather than a tautology.  #src
function strain(m::ReidMode, x::Real, theta::Real)  #src
    l  = m.l  #src
    mu = cos(theta); st = sin(theta)  #src
    P, Pp, _ = legendre_PdP(l, mu)  #src
    F, Fp, G, Gp = radial_amps(m, x)  #src
    ll = l * (l + 1)  #src
    e_rr   = Fp * P  #src
    e_thth = (G * (ll * P - mu * Pp)) / x + F * P / x  #src
    e_phph = (F * P + G * mu * Pp) / x  #src
    e_rth  = 0.5 * st * Pp * (Gp - G / x - F / x)  #src
    return e_rr, e_thth, e_phph, e_rth  #src
end  #src

# ### 1.4 What is checked, and which check actually discriminates
#
# Three residuals are monitored on every mode used below. They are not
# interchangeable, and the difference between them is the whole reason this
# section exists.
#
# **Incompressibility is a weak test.** Summing the three diagonal components
# gives ``P\bigl[F'+2F/x+l(l+1)G/x\bigr]``, and ``G``'s own definition makes
# ``l(l+1)G/x=-(2F/x+F')``. The trace is therefore analytically zero *for any
# radial profile whatsoever* -- correct, wrong, or invented. Measuring it tells
# you only that nothing went non-finite and that enough relative precision
# survived inside the sum for the cancellation to happen. It cannot detect a
# field that is smoothly and completely wrong, and in practice it does not: the
# naive unscaled implementation passes this test at ``10^{-16}`` while violating
# the boundary condition below by ``5\times10^{-7}``.
#
# **The tangential-stress boundary condition is the discriminating test.**
# Substituting the reduced forms at the surface,
#
# ```math
# G'(1)-G(1)-F(1)=-\frac{U''(1)-2U'(1)+l(l+1)U(1)}{l(l+1)},
# ```
#
# and the numerator is exactly Reid's second boundary condition: zero tangential
# stress on a free surface. So ``e_{r\theta}(1,\theta)=0`` is not an identity
# available to any field -- it holds only if ``C_j``, ``\Pi_0`` and the root
# ``q`` genuinely satisfy the eigenvalue problem.
#
# **A failure here would mean** that the object being decomposed in Section 3 is
# not a Reid eigenmode, and every bandwidth on this page would be measuring the
# spectrum of some other velocity field. It is the one residual that can fail
# for a reason other than round-off.
#
# **The characteristic-equation residual** closes the loop back to the solver:
# ``q`` is returned by `DropSolver.dominant_root`, and `DropSolver.reid_char`
# evaluates the equation it is supposed to root. Checking one against the other
# verifies that this page and the production solver are talking about the same
# eigenvalue.

radial_grid(n=200) = sort(unique(vcat(  #src
    [1 - (1 - XMIN) * (1 - t) for t in range(0, 1; length=n)],  #src
    exp.(range(log(XMIN), 0.0; length=n)))))  #src
angle_grid(n=41) = collect(range(0, pi; length=n+2))[2:end-1]  #src

## (a) max|e_rr + e_thth + e_phph| / max|e_rr|. Analytically zero; see above  #src
##     for why passing it proves very little.  #src
function incompressibility(m::ReidMode; xs=radial_grid(), ths=angle_grid())  #src
    num = 0.0; den = 0.0  #src
    for x in xs, th in ths  #src
        a, b, c, _ = strain(m, x, th)  #src
        s = abs(a + b + c)  #src
        (isfinite(s) && isfinite(abs(a))) || return NaN  #src
        num = max(num, s); den = max(den, abs(a))  #src
    end  #src
    num / den  #src
end  #src

## (b) BC2. Returns (normalised max|e_rtheta(1,.)|, |L2[U](1)|, |U(1)+1|).  #src
function bc2_residual(m::ReidMode; xs=radial_grid(60), ths=angle_grid(21))  #src
    U, Up, Upp = profile(m, 1.0)  #src
    L2 = Upp - 2Up + m.l * (m.l + 1) * U  #src
    mx = maximum(abs(strain(m, 1.0, th)[4]) for th in ths)  #src
    scale = maximum(abs(strain(m, x, th)[1]) for x in xs, th in ths)  #src
    return mx / scale, abs(L2), abs(U + 1)  #src
end  #src

## (c) The two analytic routes to U'' must agree: (B) above versus Reid's own  #src
##     radial equation. Disagreement means the Bessel algebra is wrong.  #src
function check_Upp_vs_eq9(m::ReidMode; xs=radial_grid(60))  #src
    worst = 0.0  #src
    for x in xs  #src
        U, _, Upp = profile(m, x)  #src
        eq9 = m.l * (m.l + 1) * U / x^2 - m.q^2 * U + m.q^2 * m.Pi0 * x^(m.l + 1)  #src
        worst = max(worst, abs(Upp - eq9) / max(abs(Upp), abs(eq9), 1e-300))  #src
    end  #src
    worst  #src
end  #src

## (d) The reduced G, G' must equal the literal definitions -(2F+xF')/l(l+1)  #src
##     and -(3F'+xF'')/l(l+1), checked away from the core where the literal  #src
##     form is itself trustworthy.  #src
function check_reduced_forms(m::ReidMode; xs=0.3:0.02:1.0)  #src
    wG = 0.0; wGp = 0.0; ll = m.l * (m.l + 1)  #src
    for x in xs  #src
        U, Up, Upp = profile(m, x)  #src
        F   = U / x^2  #src
        Fp  = Up / x^2 - 2U / x^3  #src
        Fpp = Upp / x^2 - 4Up / x^3 + 6U / x^4  #src
        _, _, G, Gp = radial_amps(m, x)  #src
        wG  = max(wG,  abs(G  - (-(2F + x * Fp) / ll))  / max(abs(G), 1e-300))  #src
        wGp = max(wGp, abs(Gp - (-(3Fp + x * Fpp) / ll)) / max(abs(Gp), 1e-300))  #src
    end  #src
    wG, wGp  #src
end  #src

## (e) Analytic F', G' against the Cauchy integral formula on a circle of radius  #src
##     rho about x. F is analytic in x, so the trapezoid rule on a closed contour  #src
##     converges geometrically: there is NO truncation floor, unlike a real-axis  #src
##     finite difference, so this stays sharp even at l = 50 where the field  #src
##     varies on the scale 1/|q|.  #src
function cauchy_check(m::ReidMode; xs=0.3:0.1:0.9, rho=0.05, n=64)  #src
    wF = 0.0; wG = 0.0  #src
    for x in xs  #src
        sF = zero(ComplexF64); sG = zero(ComplexF64)  #src
        for k in 0:(n-1)  #src
            phi = 2pi * k / n  #src
            F, _, G, _ = radial_amps(m, x + rho * cis(phi))  #src
            sF += F * cis(-phi); sG += G * cis(-phi)  #src
        end  #src
        _, Fp, _, Gp = radial_amps(m, x)  #src
        wF = max(wF, abs(sF / (n * rho) - Fp) / max(abs(Fp), 1e-300))  #src
        wG = max(wG, abs(sG / (n * rho) - Gp) / max(abs(Gp), 1e-300))  #src
    end  #src
    wF, wG  #src
end  #src

## (f) Legendre's equation, which was used to eliminate P'' from e_thth.  #src
function check_legendre_ode()  #src
    worst = 0.0  #src
    for l in (2, 5, 10, 20, 30, 40, 50), mu in range(-0.97, 0.97; length=25)  #src
        P, Pp, Ppp = legendre_PdP(l, mu)  #src
        r = (1 - mu^2) * Ppp - 2mu * Pp + l * (l + 1) * P  #src
        worst = max(worst, abs(r) / max(abs(l * (l + 1) * P), 1e-300))  #src
    end  #src
    worst  #src
end  #src

let  #src
    lo = check_legendre_ode()  #src
    @printf("\nSECTION 1 -- the eigenmode evaluator\n")  #src
    @printf("  ASSERTION 1: Legendre ODE residual %.2e\n", lo)  #src
    @assert lo < 1e-10 "Legendre recurrences are wrong; e_thth was derived using this ODE"  #src

    ls  = FULL ? (2, 5, 10, 20, 30, 40, 50) : (2, 10, 30, 50)  #src
    Ohs = FULL ? (0.025, 0.2, 57.4) : (0.025, 0.2)  #src
    wch = 0.0; wbc = 0.0; wl2 = 0.0; winc = 0.0  #src
    weq = 0.0; wred = 0.0; wcau = 0.0; wbc1 = 0.0  #src
    @printf("  %5s %8s | %-10s %-10s %-10s %-10s %-10s\n",  #src
            "l", "Oh", "reid_char", "BC2 e_rth", "incompr", "U'' vs eq", "Cauchy")  #src
    for l in ls, Oh in Ohs  #src
        m = reid_mode(Oh, l)  #src
        ch = abs(reid_char(m.q, m.Oh, m.l))  #src
        bc, l2, bc1 = bc2_residual(m)  #src
        inc = incompressibility(m; xs=radial_grid(60), ths=angle_grid(21))  #src
        eq  = check_Upp_vs_eq9(m)  #src
        rg, rgp = check_reduced_forms(m)  #src
        cf, cg  = cauchy_check(m)  #src
        @printf("  %5d %8g | %-10.2e %-10.2e %-10.2e %-10.2e %-10.2e\n",  #src
                l, Oh, ch, bc, inc, eq, max(cf, cg))  #src
        wch = max(wch, ch); wbc = max(wbc, bc); wl2 = max(wl2, l2)  #src
        winc = max(winc, inc); weq = max(weq, eq)  #src
        wred = max(wred, rg, rgp); wcau = max(wcau, cf, cg); wbc1 = max(wbc1, bc1)  #src
    end  #src
    @printf("  ASSERTION 2: worst |reid_char(q)| over the grid = %.2e\n", wch)  #src
    @assert wch < 1e-9 "dominant_root and reid_char disagree: this page and the solver\n" *  #src
                       "would be using different eigenvalues"  #src
    @printf("  ASSERTION 3: worst BC2 (DISCRIMINATING) = %.2e, worst |L2[U](1)| = %.2e\n", wbc, wl2)  #src
    @assert wbc < 1e-12 "tangential stress does not vanish at the free surface: the field\n" *  #src
                        "being decomposed below is not Reid's eigenmode"  #src
    @printf("  ASSERTION 4: worst BC1 residual |U(1)+1| = %.2e\n", wbc1)  #src
    @assert wbc1 < 1e-14 "U(1) = -1 is supposed to hold by construction"  #src
    @printf("  ASSERTION 5: worst incompressibility (WEAK) = %.2e\n", winc)  #src
    @assert winc < 1e-8 "strain components lost too much relative precision"  #src
    @printf("  ASSERTION 6: worst |U'' (Bessel ODE) - U'' (Reid eq.)| = %.2e\n", weq)  #src
    @assert weq < 1e-9 "the Bessel-recurrence algebra behind eq. (B) is wrong"  #src
    @printf("  ASSERTION 7: worst reduced-vs-literal G, G' = %.2e\n", wred)  #src
    @assert wred < 1e-6 "the algebraic reduction of G, G' does not reproduce the definitions"  #src
    @printf("  ASSERTION 8: worst analytic-vs-Cauchy derivative = %.2e\n", wcau)  #src
    @assert wcau < 1e-8 "the analytic derivatives (A), (B) disagree with a contour integral"  #src
end  #src

# ## 2. Superposition, and the period-averaged shear rate
#
# A real drop is not in one eigenmode. Its state is a vector of modal surface
# velocities ``\dot A_l``, ``l=2\ldots M``, and the interior flow is the
# corresponding superposition. Reid's normalisation ``U(1)=-1`` means
# ``F(1)=-1``, so weighting mode ``l`` by ``-\dot A_l`` makes the surface radial
# velocity exactly ``u_r(1,\mu)=\sum_l\dot A_l P_l(\mu)``, which is the
# convention the solver's state vector already uses. All modes are given a
# common oscillation phase.
#
# Each mode is complex, so the physical strain rate at phase ``\varphi`` is
# ``\bm e=\operatorname{Re}(\bm E e^{\mathrm{i}\varphi})`` for a complex tensor
# ``\bm E``. The invariant that the constitutive law needs is
# ``S=\sqrt{2\,\bm e\!:\!\bm e}``, and this can be reduced before any numerics.
# For a single complex component ``E``,
# ``\operatorname{Re}(Ee^{\mathrm{i}\varphi})^2
# =\tfrac12|E|^2+\tfrac12\operatorname{Re}(E^2e^{2\mathrm{i}\varphi})``, so
# summing over components with the multiplicities ``(1,1,1,2)`` of
# ``(rr,\theta\theta,\phi\phi,r\theta)``,
#
# ```math
# S(\varphi)^2 = a + \operatorname{Re}\!\bigl(b\,e^{2\mathrm{i}\varphi}\bigr),
# \qquad
# a=\sum_k m_k|E_k|^2,\quad b=\sum_k m_k E_k^2 .
# ```
#
# **In plain English:** the shear rate has period ``\pi``, not ``2\pi`` -- it
# does not care which way the fluid is moving -- and its whole time dependence
# is one cosine. The period average is therefore a single one-dimensional
# integral of a smooth positive function,
#
# ```math
# \bar S=\frac{1}{2\pi}\int_0^{2\pi}\sqrt{a+|b|\cos u}\;du ,
# ```
#
# done with the periodic trapezoid rule. Cauchy-Schwarz guarantees
# ``a\ge|b|``, so the radicand never goes negative. The one place convergence
# degrades is ``|b|/a\to1``, where every tensor component shares a single complex
# phase, ``S\propto|\cos|``, and the integrand has a square-root cusp. That is
# the deep core, where the Bessel part of ``U`` has underflowed and the profile
# is a single complex constant times a real function -- and where ``S`` is
# negligible and ``\eta`` is unthinned anyway. The reduction is checked against
# a direct trapezoid rule over the raw components in Section 4.
#
# The superposed field gets its own boundary-condition check. Every mode
# satisfies ``e_{r\theta}(1,\theta)=0`` individually, so anything nonzero in the
# sum is cancellation loss between modes -- precisely what has to be ruled out
# before the sum is used as an input to a nonlinear constitutive law.

function build_field(Oh, ls, Adot, xs, mus)  #src
    nx, nm = length(xs), length(mus)  #src
    Err = zeros(ComplexF64, nx, nm); Ett = zeros(ComplexF64, nx, nm)  #src
    Epp = zeros(ComplexF64, nx, nm); Ert = zeros(ComplexF64, nx, nm)  #src
    bc2_worst = 0.0  #src
    for (k, l) in enumerate(ls)  #src
        w = -Adot[k]  #src
        w == 0 && continue  #src
        m = reid_mode(Oh, l); ll = l * (l + 1)  #src
        bc2_worst = max(bc2_worst, bc2_residual(m)[1])  #src
        Fv = Vector{ComplexF64}(undef, nx); Fpv = Vector{ComplexF64}(undef, nx)  #src
        Gv = Vector{ComplexF64}(undef, nx); Gpv = Vector{ComplexF64}(undef, nx)  #src
        for (i, x) in enumerate(xs)  #src
            Fv[i], Fpv[i], Gv[i], Gpv[i] = radial_amps(m, x)  #src
        end  #src
        Pl = Vector{Float64}(undef, nm); Ppl = Vector{Float64}(undef, nm)  #src
        for (j, mu) in enumerate(mus)  #src
            P, Pp, _ = legendre_PdP(l, mu); Pl[j] = P; Ppl[j] = Pp  #src
        end  #src
        @inbounds for j in 1:nm  #src
            mu = mus[j]; st = sqrt(max(0.0, 1 - mu * mu)); P = Pl[j]; Pp = Ppl[j]  #src
            for i in 1:nx  #src
                x = xs[i]; F = Fv[i]; Fp = Fpv[i]; G = Gv[i]; Gp = Gpv[i]  #src
                Err[i,j] += w * (Fp * P)  #src
                Ett[i,j] += w * ((G * (ll * P - mu * Pp)) / x + F * P / x)  #src
                Epp[i,j] += w * ((F * P + G * mu * Pp) / x)  #src
                Ert[i,j] += w * (0.5 * st * Pp * (Gp - G / x - F / x))  #src
            end  #src
        end  #src
    end  #src
    for A in (Err, Ett, Epp, Ert)  #src
        all(z -> isfinite(abs(z)), A) || error("non-finite strain component at Oh=$Oh")  #src
    end  #src
    Err, Ett, Epp, Ert, bc2_worst  #src
end  #src

## BC2 of the SUPERPOSED field, exactly at x = 1.  #src
function field_bc2(Oh, ls, Adot, mus, scale)  #src
    acc = zeros(ComplexF64, length(mus))  #src
    for (k, l) in enumerate(ls)  #src
        Adot[k] == 0 && continue  #src
        m = reid_mode(Oh, l); w = -Adot[k]  #src
        F, _, G, Gp = radial_amps(m, 1.0)  #src
        for (j, mu) in enumerate(mus)  #src
            _, Pp, _ = legendre_PdP(l, mu)  #src
            acc[j] += w * (0.5 * sqrt(max(0.0, 1 - mu * mu)) * Pp * (Gp - G - F))  #src
        end  #src
    end  #src
    maximum(abs, acc) / scale  #src
end  #src

## Period-averaged S via the (a, b) reduction. Also returns max|b|/a restricted  #src
## to points carrying more than 1% of max S, i.e. the worst cusp that matters.  #src
function sbar_field(Err, Ett, Epp, Ert; N=256)  #src
    nx, nm = size(Err)  #src
    S = zeros(nx, nm); rat = zeros(nx, nm)  #src
    cu = [cos(2pi * k / N) for k in 0:N-1]  #src
    @inbounds for j in 1:nm, i in 1:nx  #src
        A = Err[i,j]; B = Ett[i,j]; C = Epp[i,j]; D = Ert[i,j]  #src
        a = abs2(A) + abs2(B) + abs2(C) + 2 * abs2(D)  #src
        a == 0 && continue  #src
        b = abs(A * A + B * B + C * C + 2 * D * D)  #src
        rat[i,j] = b / a  #src
        s = 0.0  #src
        for n in 1:N  #src
            v = a + b * cu[n]; s += sqrt(v > 0 ? v : 0.0)  #src
        end  #src
        S[i,j] = s / N  #src
    end  #src
    sm = maximum(S); kmx = 0.0  #src
    for i in eachindex(S); S[i] > 0.01 * sm && (kmx = max(kmx, rat[i])); end  #src
    S, kmx  #src
end  #src

## Reference: the same average taken directly on the four real components, with  #src
## no algebraic reduction at all.  #src
function sbar_direct(Err, Ett, Epp, Ert; N=1024)  #src
    nx, nm = size(Err); S = zeros(nx, nm)  #src
    @inbounds for j in 1:nm, i in 1:nx  #src
        s = 0.0  #src
        for n in 0:N-1  #src
            z = cis(pi * n / N)  #src
            a = real(Err[i,j] * z); b = real(Ett[i,j] * z)  #src
            c = real(Epp[i,j] * z); d = real(Ert[i,j] * z)  #src
            s += sqrt(2 * (a * a + b * b + c * c + 2 * d * d))  #src
        end  #src
        S[i,j] = s / N  #src
    end  #src
    S  #src
end  #src

# ## 3. The viscosity field and its Legendre spectrum
#
# The fluid is the repository's fitted Cross validation fluid -- the 3000 ppm
# shear-thinning solution the solver's rebound predictions are checked against:
#
# ```math
# \frac{\eta(S)}{\eta_0}=r+\frac{1-r}{1+(\lambda_c S)^{m}},
# \qquad r=\frac{\eta_\infty}{\eta_0},
# ```
#
# with ``\eta_0=8.434``, ``\eta_\infty=3.732\times10^{-3}`` Pa s,
# ``k=18.48`` s and ``m=0.7431``, on a drop of radius
# ``R=0.3`` mm at ``\mathrm{Bo}=0.012``. Non-dimensionalising by the
# inertio-capillary time ``T_\sigma=\sqrt{\rho R^3/\sigma}`` gives
# ``\lambda_c=k/T_\sigma\approx3.05\times10^{4}``, so the knee of the flow curve
# sits at ``S\approx3.3\times10^{-5}`` in units of ``T_\sigma^{-1}``: this fluid
# is thinned by extremely gentle motion, and the interesting regime is
# ``\lambda_c S\gg1``.
#
# Feeding ``\bar S`` into the Cross law gives ``\eta(x,\mu)``, whose Legendre
# coefficients follow from orthogonality,
#
# ```math
# \eta_{l'}(x)=\frac{2l'+1}{2}\int_{-1}^{1}\eta(x,\mu)\,P_{l'}(\mu)\,d\mu ,
# ```
#
# evaluated by Gauss-Legendre in ``\mu``. Two summaries are reported: the value
# at the surface ``x=1``, and the volume average
# ``\int\eta_{l'}(x)\,x^2dx\big/\!\int x^2dx``. The volume average is the one
# quoted, because it is what enters the radial integrals of the coupling matrix;
# the surface value is a single-node diagnostic and converges more slowly.
#
# ### 3.1 The choice of error norm changes the answer
#
# This is the step where it is easiest to reach a comfortable and wrong
# conclusion, so both norms are reported.
#
# The instinctive metric is a **power fraction**: the smallest ``L`` holding, say,
# 99% of ``\sum_{l'\ge1}\eta_{l'}^2\cdot\frac{2}{2l'+1}``, the anisotropic part
# of the ``L^2`` norm. (The ``l'=0`` term has to be excluded or every answer is
# ``L=0``: ``\eta`` is overwhelmingly its own spherical mean. That variant is
# reported once below, for honesty, and then dropped.) This metric is
# **misleading here**, and it is worth being precise about why. Power is a sum of
# *squares*; the discarded coefficients enter the coupling matrix
# *additively*, each multiplying an ``O(1)`` Gaunt-type factor. A plateau of two
# hundred coefficients at ``10^{-2}`` carries a negligible fraction of the power
# and a summed magnitude of order one.
#
# The quantity that actually controls the banding error is therefore
#
# ```math
# \mathcal T_1(L)=\frac{1}{|\eta_0|}\sum_{l'>L}\bigl|\eta_{l'}\bigr| ,
# ```
#
# the relative size of the coupling thrown away by banding at ``L``. Its RMS
# counterpart ``\mathcal T_2`` is tabulated alongside it as the optimistic
# bound, and the gap between the two is itself part of the result. Define
# ``L_\eta(t)`` as the smallest ``L`` with ``\mathcal T_1(L)<t``; the headline
# number is ``t=10^{-2}``, i.e. banding that discards 1% of the coupling.
#
# One deliberate detail: the search for ``L_\eta`` stops one short of the last
# measured coefficient. ``\mathcal T_1`` at the final index is identically zero
# because the sum beyond it is empty, and accepting that would report a
# bandwidth that is an artefact of where the table was cut off rather than a
# property of the field.

## The fitted Cross validation fluid (identical constants to  #src
## julia/scripts/validate_shear_thinning.jl and docs/figures.jl).  #src
const ETA_INF = 0.0037320997942061666  #src
const ETA_0   = 8.433817577956766  #src
const K_CROSS = 18.48081673111359  #src
const M_CROSS = 0.7430524574330837  #src
const BO = 0.012; const RDROP = 0.0003; const SIGMA = 0.0728; const GRAV = 9.81  #src
const RHO       = BO * SIGMA / (GRAV * RDROP^2)  #src
const T_SIGMA   = sqrt(RHO * RDROP^3 / SIGMA)  #src
const OH0       = ETA_0 / sqrt(RHO * SIGMA * RDROP)  #src
const LAMBDA_C  = K_CROSS / T_SIGMA  #src
const ETA_RATIO = ETA_INF / ETA_0  #src

@inline eta_cross(S) = ETA_RATIO + (1 - ETA_RATIO) / (1 + (LAMBDA_C * S)^M_CROSS)  #src

## Gauss-Legendre by Newton iteration on P_n -- deliberately a different  #src
## algorithm from DropSolver's Golub-Welsch routine, so the two cross-check.  #src
function gl(n::Int, a::Float64=-1.0, b::Float64=1.0)  #src
    x = zeros(n); w = zeros(n)  #src
    for i in 1:n  #src
        z = cos(pi * (i - 0.25) / (n + 0.5)); pp = 0.0  #src
        for _ in 1:100  #src
            p0 = 1.0; p1 = 0.0  #src
            for j in 1:n  #src
                p2 = p1; p1 = p0  #src
                p0 = ((2j - 1) * z * p1 - (j - 1) * p2) / j  #src
            end  #src
            pp = n * (z * p0 - p1) / (z * z - 1)  #src
            dz = p0 / pp; z -= dz  #src
            abs(dz) < 1e-16 && break  #src
        end  #src
        x[i] = z; w[i] = 2 / ((1 - z * z) * pp * pp)  #src
    end  #src
    p = sortperm(x); xs = x[p]; ws = w[p]  #src
    @. xs = (b - a) / 2 * xs + (b + a) / 2  #src
    @. ws = (b - a) / 2 * ws  #src
    xs, ws  #src
end  #src

## Panel Gauss-Legendre in radius, graded towards the surface because high-l  #src
## modes are a thin boundary layer at x = 1.  #src
function radial_quad(; edges=[0.0,0.3,0.6,0.8,0.9,0.95,0.98,0.99,0.995,1.0], n=12)  #src
    xs = Float64[]; ws = Float64[]  #src
    for k in 1:length(edges)-1  #src
        nd, wd = gl(n, edges[k], edges[k+1]); append!(xs, nd); append!(ws, wd)  #src
    end  #src
    xs, ws  #src
end  #src

mu_quad(n) = gl(n, -1.0, 1.0)  #src

function legendre_matrix(mus, LMAX)  #src
    P = zeros(LMAX + 1, length(mus))  #src
    @inbounds for (j, mu) in enumerate(mus)  #src
        P[1,j] = 1.0  #src
        LMAX >= 1 && (P[2,j] = mu)  #src
        for n in 2:LMAX  #src
            P[n+1,j] = ((2n - 1) * mu * P[n,j] - (n - 1) * P[n-1,j]) / n  #src
        end  #src
    end  #src
    P  #src
end  #src

function eta_spectrum(S, scale, xs, wx, mus, wmu, Pmat)  #src
    LMAX = size(Pmat, 1) - 1; nx, nm = size(S)  #src
    cr = zeros(LMAX + 1, nx); etamin = Inf; etamax = -Inf  #src
    @inbounds for i in 1:nx  #src
        for j in 1:nm  #src
            et = eta_cross(scale * S[i,j])  #src
            etamin = min(etamin, et); etamax = max(etamax, et)  #src
            wj = wmu[j] * et  #src
            for lp in 0:LMAX; cr[lp+1,i] += wj * Pmat[lp+1,j]; end  #src
        end  #src
        for lp in 0:LMAX; cr[lp+1,i] *= (2lp + 1) / 2; end  #src
    end  #src
    vw = wx .* xs.^2; vw = vw ./ sum(vw)  #src
    cr * vw, cr[:, argmax(xs)], etamin, etamax  #src
end  #src

## Smallest L holding fraction f of the ANISOTROPIC power. Set lo=0 to include  #src
## the mean, which makes every answer 0 and is reported once for honesty.  #src
function Lpow(c, f; lo=1)  #src
    LMAX = length(c) - 1  #src
    p = [lp < lo ? 0.0 : c[lp+1]^2 * 2 / (2lp + 1) for lp in 0:LMAX]  #src
    tot = sum(p); tot == 0 && return -1  #src
    k = findfirst(>=(f), cumsum(p) ./ tot)  #src
    k === nothing ? -1 : k - 1  #src
end  #src

T1(c, L) = sum(abs, @view c[L+2:end]) / abs(c[1])  #src
T2(c, L) = sqrt(sum(abs2, @view c[L+2:end])) / abs(c[1])  #src

function Ltail(c, tol)  #src
    for L in 0:length(c)-2  #src
        T1(c, L) < tol && return L  #src
    end  #src
    -1  #src
end  #src
fmtL(L) = L < 0 ? "  >T" : @sprintf("%4d", L)  #src

# ## 4. Grids, and the evidence that they are fine enough
#
# Four independent discretisations enter, and a bandwidth is exactly the kind of
# quantity that a coarse grid inflates by aliasing, so each is refined and the
# answer shown to be resolved:
#
# * the ``\mu`` quadrature used for the Legendre projection,
# * the radial panel quadrature used for the volume average,
# * the number of phase points in the period average,
# * and the algebraic ``(a,b)`` reduction of Section 2, checked against a direct
#   average over the raw tensor components.
#
# A fifth check is a cross-check rather than a refinement: the Gauss-Legendre
# nodes used here are generated by Newton iteration on ``P_n``, while
# `DropSolver.gauss_legendre_nodes` uses the Golub-Welsch eigenvalue algorithm.
# Agreement to machine precision confirms that this page and the running solver
# integrate over the same points.
#
# The refinements are run on the hardest configuration available -- a real
# solver state, where ``\eta`` swings by two orders of magnitude with sharp
# transitions at the near-nodes of ``S`` -- because a smooth field would make
# the test easy and prove nothing about the case that matters.

## `family` records what a configuration is evidence FOR, so the assertions can  #src
## be scoped honestly:  #src
##   :realistic -- algebraically decaying modal spectra and real solver states,  #src
##                 i.e. what the solver actually produces. These carry the  #src
##                 headline claim L_eta >= M.  #src
##   :contrast  -- exploratory shapes included to show that L_eta tracks the  #src
##                 CONTRAST of the eta field rather than the mode count. They  #src
##                 are held only to the weaker multi-mode claim.  #src
struct Spectrum  #src
    name   :: String  #src
    ls     :: Vector{Int}  #src
    Adot   :: Vector{Float64}  #src
    family :: Symbol  #src
end  #src

single_mode(l) = Spectrum("single l=$l", [l], [1.0], :single)  #src
decay(p, M) = Spectrum("l^-$p, l=2..$M", collect(2:M),  #src
                       [float(l)^(-p) for l in 2:M], :realistic)  #src
expdecay(q, M) = Spectrum("$(q)^l, l=2..$M", collect(2:M),  #src
                          [float(q)^l for l in 2:M], :contrast)  #src
bump(l0, w, M) = Spectrum("bump at l=$l0 (w=$w), l=2..$M", collect(2:M),  #src
                          [exp(-((l - l0) / w)^2) for l in 2:M], :contrast)  #src

## Real modal velocity vectors saved from live impact solves; see the header of  #src
## each data file for the run that produced it.  #src
function load_adot(file, label)  #src
    f = joinpath(@__DIR__, "data", file)  #src
    isfile(f) || error("missing modal-velocity data file $f")  #src
    v = Float64[]  #src
    for line in eachline(f)  #src
        s = strip(line)  #src
        (isempty(s) || startswith(s, "#")) && continue  #src
        push!(v, parse(Float64, s))  #src
    end  #src
    Spectrum(label, collect(2:(length(v) + 1)), v, :realistic)  #src
end  #src

## The mu resolution is NOT reduced for CI: the Legendre projection is the one  #src
## grid whose refinement moves the answer rather than merely the error bar.  #src
const NMU    = 600  #src
const NPAN   = FULL ? 12  : 8  #src
const NPHASE = FULL ? 512 : 256  #src
const LMAX   = 140  #src
const XS, WX = radial_quad(n=NPAN)  #src
const MUS, WMU = mu_quad(NMU)  #src
const PMAT = legendre_matrix(MUS, LMAX)  #src
## The operating point at which the anisotropy ratios are quoted: deep in the  #src
## saturated shear-thinned regime, lambda_c * max(Sbar) = 1e3.  #src
const TARGET = 1e3  #src

let  #src
    @printf("\nSECTION 4 -- grids\n")  #src
    @printf("  fluid: Oh_0 = %.4f, Oh_inf = %.4f, lambda_c = %.4e, m = %.5f\n",  #src
            OH0, OH0 * ETA_RATIO, LAMBDA_C, M_CROSS)  #src
    @printf("  grids: %d radial nodes (%d per panel), %d-node GL in mu, l' <= %d, %d phases\n",  #src
            length(XS), NPAN, NMU, LMAX, NPHASE)  #src

    g1, w1 = gl(NMU)  #src
    g2, w2 = gauss_legendre_nodes(NMU, -1.0, 1.0)  #src
    dn = maximum(abs, g1 .- g2); dw = maximum(abs, w1 .- w2)  #src
    @printf("  ASSERTION 9: Newton GL vs DropSolver Golub-Welsch: |dx| %.2e, |dw| %.2e\n", dn, dw)  #src
    @assert dn < 1e-13 && dw < 1e-13 "this page integrates over different nodes than the solver"  #src

    rs = load_adot("adot_clean_M30.txt", "real M=30")  #src
    E = build_field(0.2, rs.ls, rs.Adot, XS, MUS)  #src

    ## phase-average refinement  #src
    Sref, _ = sbar_field(E[1], E[2], E[3], E[4]; N=4 * NPHASE)  #src
    dph = 0.0  #src
    for N in (div(NPHASE, 2), NPHASE)  #src
        SA, _ = sbar_field(E[1], E[2], E[3], E[4]; N=N)  #src
        d = maximum(abs, SA .- Sref) / maximum(Sref)  #src
        @printf("  phase average N=%-5d vs N=%-5d: max rel diff %.2e\n", N, 4 * NPHASE, d)  #src
        N == NPHASE && (dph = d)  #src
    end  #src
    @printf("  ASSERTION 10: phase average converged to %.2e\n", dph)  #src
    @assert dph < 1e-4 "period average is not resolved in phase"  #src

    ## the (a,b) reduction against a direct component average  #src
    S0, kmax = sbar_field(E[1], E[2], E[3], E[4]; N=NPHASE)  #src
    SD = sbar_direct(E[1], E[2], E[3], E[4])  #src
    dab = maximum(abs, S0 .- SD) / maximum(SD)  #src
    @printf("  ASSERTION 11: (a,b) reduction vs direct component average: %.2e", dab)  #src
    @printf("   (worst |b|/a where S > 1%% of max = %.4f)\n", kmax)  #src
    @assert dab < 1e-4 "the closed form S^2 = a + Re(b e^{2 i phi}) is wrong"  #src

    ## mu and radial refinement of the spectrum itself  #src
    sc = TARGET / (LAMBDA_C * maximum(S0))  #src
    cref, sref, _, _ = eta_spectrum(S0, sc, XS, WX, MUS, WMU, PMAT)  #src
    dmu = 0.0  #src
    for n in (div(NMU, 2), 2 * NMU)  #src
        mu2, wm2 = mu_quad(n); P2 = legendre_matrix(mu2, LMAX)  #src
        E2 = build_field(0.2, rs.ls, rs.Adot, XS, mu2)  #src
        S2, _ = sbar_field(E2[1], E2[2], E2[3], E2[4]; N=NPHASE)  #src
        c2, s2, _, _ = eta_spectrum(S2, sc, XS, WX, mu2, wm2, P2)  #src
        d = maximum(abs, c2 .- cref) / abs(cref[1])  #src
        @printf("  mu quadrature n=%-5d vs n=%-5d: max_l' |dc|/|c_0| = %.2e (vol) %.2e (x=1)\n",  #src
                n, NMU, d, maximum(abs, s2 .- sref) / abs(sref[1]))  #src
        n > NMU && (dmu = d)  #src
    end  #src
    @printf("  ASSERTION 12: volume-averaged spectrum converged in mu to %.2e\n", dmu)  #src
    @assert dmu < 1e-4 "the Legendre projection is not resolved in mu"  #src

    drad = 0.0  #src
    for npan in (NPAN + 6, NPAN + 12)  #src
        ## The SAME scale factor sc is used on every grid on purpose: rescaling  #src
        ## by each grid's own max S would inject an O(grid) shift into every  #src
        ## coefficient and fake a non-convergence.  #src
        xs2, wx2 = radial_quad(n=npan)  #src
        E3 = build_field(0.2, rs.ls, rs.Adot, xs2, MUS)  #src
        S3, _ = sbar_field(E3[1], E3[2], E3[3], E3[4]; N=NPHASE)  #src
        c3, _, _, _ = eta_spectrum(S3, sc, xs2, wx2, MUS, WMU, PMAT)  #src
        drad = max(drad, maximum(abs, c3 .- cref) / abs(cref[1]))  #src
        @printf("  radial quadrature %d nodes/panel (%d total) vs %d (%d): %.2e\n",  #src
                npan, length(xs2), NPAN, length(XS), drad)  #src
    end  #src
    @printf("  ASSERTION 13: volume average converged in radius to %.2e\n", drad)  #src
    @assert drad < 1e-5 "the radial volume average is not resolved"  #src

    @printf("  tail truncation: |c_%d|/|c_0| = %.2e (vol), T1(%d) = %.2e (vol)\n",  #src
            LMAX, abs(cref[end]) / abs(cref[1]), LMAX - 40, T1(cref, LMAX - 40))  #src
end  #src

# ## 5. Results
#
# ### 5.1 The reduced run
#
# The measurement below is the one this script performs by default: four modal
# spectra at ``\mathrm{Oh}=0.2``, all at the single operating point
# ``\lambda_c\max\bar S=10^{3}`` (deep in the thinned regime, where ``\eta``
# spans roughly ``6\times10^{-3}`` to ``0.7`` of ``\eta_0``), on a ``600``-node
# Gauss-Legendre grid in ``\mu`` and a ``72``-node graded radial grid. It takes
# about twenty seconds.
#
# | modal spectrum | ``M`` | ``L_\eta`` for ``\mathcal T_1<10^{-2}`` | ``\mathcal T_1(M)`` |
# |:--|--:|--:|--:|
# | single mode ``l=2`` | 2 | 8 | -- |
# | ``\propto l^{-2}``, ``M=30`` | 30 | 32 | 0.02 |
# | ``\propto l^{-1}``, ``M=30`` | 30 | 126 | 0.45 |
# | real solver state, ``M=30`` | 30 | 133 | 1.33 |
#
# Two things in that table are the actual result.
#
# **One active mode gives a narrow spectrum.** With a single ``l=2`` mode, ``S``
# has a handful of well-separated nodes, ``\eta`` is smooth, and eight harmonics
# hold it to 1%. If drops oscillated one mode at a time, banding would be a real
# saving.
#
# **A realistic multi-mode state does not.** For an algebraically decaying modal
# spectrum, and for a real solver state, ``L_\eta`` reaches or exceeds ``M``
# itself. At ``M=30`` a half-bandwidth of ``\approx133`` is needed to discard
# only 1% of the coupling -- which is to say, a dense matrix, four times wider
# than the matrix has rows. Banding buys nothing.
#
# The gap between the two error norms is worth quoting, because it is the whole
# reason the naive answer is wrong. On the real ``M=30`` state, the
# power-fraction metric of Section 3.1 reports that 90% of the anisotropic power
# lies below ``l'=1`` and 99% below ``l'=3`` -- it would sanction banding at
# ``L_\eta=3``. At ``L=3`` the summed discarded coupling is
# ``\mathcal T_1=2.96``, i.e. **three times** ``|\eta_0|``, and even the
# forgiving RMS norm is ``\mathcal T_2=0.43``. The power metric is not merely
# optimistic here; it is wrong by a factor of forty in the bandwidth.
#
# ### 5.2 The production sweep
#
# Setting `ETA_SPECTRUM_FULL=1` widens the measurement to six synthetic modal
# spectra and two real solver states, at
# ``\mathrm{Oh}\in\{0.025,0.2,2\}`` and six amplitudes spanning
# ``\lambda_c\max\bar S=10^{-2}\ldots10^{5}``, with the radial grid and phase
# average refined to the values used in the original measurement campaign. It
# takes about a minute. Four of the six rows below are inside the reduced run;
# the two marked otherwise need the full sweep:
#
# | modal spectrum | ``L_\eta`` for 1% | in the reduced run? |
# |:--|--:|:--|
# | single mode ``l=2`` | 8 | yes |
# | ``\propto l^{-2}``, ``M=30`` | 32 | yes |
# | ``\propto l^{-2}``, ``M=50`` | 50 | full sweep only |
# | ``\propto l^{-1}``, ``M=30`` | 126 | yes |
# | real solver state, ``M=30`` | 133 | yes |
# | real solver state, ``M=20`` | 139 | full sweep only |
#
# The full sweep also carries two configurations that are **not** meant to be
# realistic and are there to expose the mechanism, because the mechanism is not
# "more modes means wider". ``L_\eta`` tracks the **contrast** of the ``\eta``
# field, and mode count is only a proxy for contrast:
#
# | contrast configuration | ``\eta`` range | ``L_\eta`` for 1% |
# |:--|--:|--:|
# | ``\propto 0.85^{\,l}``, ``M=50`` | ``34\times`` | 16 |
# | Gaussian bump at ``l=10``, ``M=50`` | ``153\times`` | 79 |
#
# The exponentially decaying spectrum has 49 modes nominally active and still
# comes out *narrower than* ``M``. That is not a counterexample to the result;
# it is the result stated correctly. Its amplitudes fall off so fast that the
# superposed strain field is dominated by a handful of low modes, ``\eta`` spans
# only a factor of ``34``, and the field stays smooth. The bump, with the same
# number of modes but a hundred-and-fifty-fold swing in ``\eta``, is broad.
#
# The physical mechanism behind the contrast is **near-nodal points of the
# superposed strain field**, where ``S`` collapses and the fluid snaps back
# towards unthinned. That is a near-discontinuity in ``\mu``, and a
# near-discontinuity has a Legendre spectrum that decays only algebraically. One
# mode has few such nodes; a dozen comparable modes beating against one another
# have a dense nodal set. Real solver states have exactly that structure, which
# is why the two rows measured from live solves are the widest in the table.
#
# A sign-pattern variation was tried and is *not* an independent test, which is
# worth recording so nobody repeats it: under ``\mu\to-\mu``,
# ``P_l\to(-1)^lP_l``, so the field built from ``\dot A_l=(-1)^l l^{-p}`` is the
# exact mirror image of the one built from ``l^{-p}``. Every ``|\eta_{l'}|`` is
# bit-for-bit identical. Genuinely different *shapes* -- exponential decay, a
# bump at ``l=10``, and the real solver vectors -- are used instead.
#
# ### 5.3 What this costs the cheaper closure
#
# The same spectrum prices the closure that keeps only ``\eta_0(r)`` and throws
# away all angular structure. At the operating point, on the real ``M=30``
# state,
#
# ```math
# \frac{|\eta_1|}{|\eta_0|}=1.25,
# \qquad
# \frac{|\eta_2|}{|\eta_0|}=0.20 .
# ```
#
# **In plain English:** the angular structure of the viscosity is *comparable to
# or larger than its own mean*. Discarding it is a leading-order modelling
# error, not a perturbative correction, and any model that replaces the field by
# a single effective Ohnesorge number inherits that error.
#
# The contrast with the single-mode picture is again the point. With one ``l=2``
# mode active, ``\eta`` varies by only ``1.1``--``1.2\times`` across the drop and
# spatial homogenisation looks like a 10% effect. That estimate is real, and it
# is real only for one mode.
#
# ### 5.4 Caveat: every entry is a lower bound
#
# Coefficients above ``l'=140`` are not measured. The reported
# ``\mathcal T_1(L)`` therefore omits the tail beyond the table, so every
# ``L_\eta`` above is a **lower** bound on the bandwidth required -- equivalently,
# an upper bound on how narrow ``\eta`` is. The residual matters: on the real
# ``M=30`` state ``\mathcal T_1(100)\approx0.14``, still fourteen percent of
# ``|\eta_0|``, and ``|\eta_{140}|/|\eta_0|\approx1.2\times10^{-3}`` shows the
# spectrum has not finished decaying where the table stops. Nothing on this page
# would improve if the table were extended; the numbers would only get worse.

function run_config(sp::Spectrum, Oh; targets=(TARGET,), report=false)  #src
    Err, Ett, Epp, Ert, bc2mode = build_field(Oh, sp.ls, sp.Adot, XS, MUS)  #src
    scale_rr = maximum(abs, Err)  #src
    tr   = maximum(abs, Err .+ Ett .+ Epp) / scale_rr  #src
    bc2f = field_bc2(Oh, sp.ls, sp.Adot, MUS, scale_rr)  #src
    S, kmax = sbar_field(Err, Ett, Epp, Ert; N=NPHASE)  #src
    Smax = maximum(S)  #src
    @printf("\n  %s   Oh = %g   (%d active modes, l = %d..%d)\n",  #src
            sp.name, Oh, count(!=(0), sp.Adot), minimum(sp.ls), maximum(sp.ls))  #src
    @printf("    gates: BC2 per-mode %.2e | BC2 superposed %.2e  <-- DISCRIMINATING\n",  #src
            bc2mode, bc2f)  #src
    @printf("           incompressibility %.2e (WEAK: the trace is algebraically 0)\n", tr)  #src
    @assert bc2mode < 1e-12 "a constituent mode violates the free-surface condition"  #src
    @assert bc2f < 1e-12 "the superposition lost the free-surface condition to cancellation"  #src
    @assert tr < 1e-8 "strain components lost too much relative precision"  #src
    out = Dict{Float64,Any}()  #src
    @printf("    %-9s %-19s | %-14s | %-14s | %-14s\n", "lam*Smax", "eta/eta_0 min..max",  #src
            "Lpow 90/99/99.9", "Ltail vol 1e-2/-3", "Ltail x=1 1e-2")  #src
    for tgt in targets  #src
        sc = tgt / (LAMBDA_C * Smax)  #src
        cv, cs, emin, emax = eta_spectrum(S, sc, XS, WX, MUS, WMU, PMAT)  #src
        @printf("    %-9.0e %.3e %.3e | %4d %4d %4d | %s %s      | %s\n",  #src
                tgt, emin, emax, Lpow(cv, 0.90), Lpow(cv, 0.99), Lpow(cv, 0.999),  #src
                fmtL(Ltail(cv, 1e-2)), fmtL(Ltail(cv, 1e-3)), fmtL(Ltail(cs, 1e-2)))  #src
        out[tgt] = (cv, cs)  #src
    end  #src
    if report  #src
        cv = out[TARGET][1]  #src
        @printf("    Lpow INCLUDING l'=0 (uninformative -- eta is mostly its own mean): %d/%d/%d\n",  #src
                Lpow(cv, 0.90; lo=0), Lpow(cv, 0.99; lo=0), Lpow(cv, 0.999; lo=0))  #src
        @printf("    banding error  T1(L) / T2(L):")  #src
        for L in (0, 2, 3, 4, 8, 16, 30, 60, 100)  #src
            @printf(" L=%d:%.1e/%.1e", L, T1(cv, L), T2(cv, L))  #src
            L == 8 && @printf("\n                                 ")  #src
        end  #src
        println()  #src
    end  #src
    (S, Smax, out)  #src
end  #src

let  #src
    @printf("\nSECTION 5 -- the measurement\n")  #src
    tg = FULL ? (1e-2, 1e0, 1e2, 1e3, 1e4, 1e5) : (TARGET,)  #src

    ## (i) one active mode: the narrow case  #src
    _, _, o1 = run_config(single_mode(2), 0.2; targets=tg, report=true)  #src
    c1 = o1[TARGET][1]  #src
    L1 = Ltail(c1, 1e-2)  #src
    @printf("  ASSERTION 14: single l=2 mode gives L_eta(1%%) = %d\n", L1)  #src
    @assert 0 <= L1 <= 12 "a single active mode should give a NARROW eta spectrum;\n" *  #src
                          "if this fails the contrast that is this page's result has gone"  #src

    ## (ii) realistic multi-mode states: the broad case  #src
    configs = Any[decay(2, 30), decay(1, 30), load_adot("adot_clean_M30.txt", "real solver state M=30")]  #src
    FULL && append!(configs, Any[decay(2, 50), decay(1, 50), expdecay(0.85, 50), bump(10, 4, 50),  #src
                                load_adot("adot_clean_M20.txt", "real solver state M=20")])  #src
    creal = nothing  #src
    for sp in configs  #src
        M = maximum(sp.ls)  #src
        Ohs = FULL ? (0.025, 0.2, 2.0) : (0.2,)  #src
        for Oh in Ohs  #src
            isreal_state = startswith(sp.name, "real")  #src
            _, _, o = run_config(sp, Oh; targets=tg, report=(Oh == 0.2))  #src
            Oh == 0.2 || continue  #src
            c = o[TARGET][1]  #src
            L = Ltail(c, 1e-2)  #src
            @printf("  ASSERTION 15 [%s, %s]: L_eta(1%%) = %s vs M = %d, T1(M) = %.2f\n",  #src
                    sp.name, sp.family, fmtL(L), M, T1(c, M))  #src
            ## Weak claim, held by EVERY multi-mode configuration: turning on a  #src
            ## dozen modes makes the viscosity spectrum substantially broader  #src
            ## than the single-mode case. This is the contrast that is the  #src
            ## result; if it went away there would be nothing on this page.  #src
            @assert L < 0 || L > L1 "a multi-mode state gave a bandwidth no broader than a\n" *  #src
                                    "single active mode; the contrast this page reports has gone"  #src
            ## Strong claim, held only by the spectra a solver actually produces:  #src
            ## the bandwidth needed reaches the shape truncation itself, so  #src
            ## banding saves nothing. L < 0 means 1% was never reached anywhere  #src
            ## in the measured table, which is a stronger version of the same  #src
            ## statement.  #src
            if sp.family === :realistic  #src
                @assert L < 0 || L >= M "a realistic modal spectrum gave a bandwidth NARROWER\n" *  #src
                                        "than M; banding would then be worth implementing"  #src
            end  #src
            ## For a real solver state the statement is stronger still: banding  #src
            ## at the shape truncation leaves a discarded coupling of order  #src
            ## |eta_0| itself.  #src
            if isreal_state  #src
                @printf("               (real state: T1(M) = %.2f -- banding at L = M discards\n", T1(c, M))  #src
                @printf("                coupling of order |eta_0| itself)\n")  #src
                @assert T1(c, M) > 0.5 "banding a real solver state at L = M no longer discards\n" *  #src
                                       "an O(|eta_0|) amount of coupling"  #src
            end  #src
            isreal_state && M == 30 && (creal = c)  #src
        end  #src
    end  #src

    ## (iii) the anisotropy ratios that price the eta = eta(r) closure  #src
    @assert creal !== nothing "the real M=30 state was not measured"  #src
    r1 = abs(creal[2]) / abs(creal[1])  #src
    r2 = abs(creal[3]) / abs(creal[1])  #src
    @printf("  ASSERTION 16: |eta_1|/|eta_0| = %.3f, |eta_2|/|eta_0| = %.3f\n", r1, r2)  #src
    @assert 1.15 < r1 < 1.45 "|eta_1|/|eta_0| is not ~1.3; the pricing of the eta(r) closure\n" *  #src
                             "in Section 5.3 no longer holds"  #src
    @assert 0.15 < r2 < 0.28 "|eta_2|/|eta_0| is not ~0.2"  #src

    ## (iv) the power-fraction metric is misleading, quantified  #src
    Lp = Lpow(creal, 0.99)  #src
    @printf("  ASSERTION 17: power metric sanctions L=%d, where T1 = %.2f and T2 = %.2f\n",  #src
            Lp, T1(creal, Lp), T2(creal, Lp))  #src
    @assert Lp <= 5 "the power metric no longer gives the small answer Section 3.1 warns about"  #src
    @assert T1(creal, Lp) > 1.0 "the summed-magnitude norm no longer contradicts the power norm"  #src

    ## (v) the tail beyond the table, i.e. how much of a lower bound this is  #src
    @printf("  ASSERTION 18: unmeasured tail: |c_%d|/|c_0| = %.2e, T1(100) = %.3f\n",  #src
            LMAX, abs(creal[end]) / abs(creal[1]), T1(creal, 100))  #src
    @assert abs(creal[end]) / abs(creal[1]) > 1e-6 "the spectrum has decayed to nothing by l'=$LMAX;\n" *  #src
                                                   "the lower-bound caveat in Section 5.4 would be overstated"  #src
end  #src

# ## 6. Consequences
#
# Two rungs of the shear-thinning model hierarchy are priced by this
# measurement, and both prices are higher than they look.
#
# **Truncating ``\eta``'s spectrum at ``L_\eta``** is not a saving. It requires
# ``L_\eta\gtrsim M`` for algebraically decaying modal spectra and
# ``L_\eta\gg M`` for real solver states, so the "banded" matrix is the dense
# matrix with extra bookkeeping. The right move is to take the coupled system
# dense, or to skip past it. The exception found in Section 5.2 -- an
# exponentially decaying modal spectrum, which does band narrowly -- is not a
# reprieve: it describes a drop ringing in two or three modes, not one that has
# just been struck.
#
# **Keeping only ``\eta_0(r)``** -- a spherically symmetric viscosity -- is the
# cheapest defensible rung, because a spherically symmetric ``\eta`` commutes
# with the angular Laplacian and every mode decouples again, restoring the
# solver's one-oscillator-per-mode architecture exactly. What it costs is
# ``|\eta_1|/|\eta_0|=1.25``: not an accurate model, but an honest one, and the
# honesty is the point of measuring rather than assuming.

println("\n" * "="^92)  #src
println("ALL ASSERTIONS PASSED")  #src
println(FULL ? "(full production sweep)" :  #src
               "(reduced CI run; the l^-2 M=50 and real M=20 rows of Section 5.2, and the\n" *  #src
               " two contrast configurations, need ETA_SPECTRUM_FULL=1)")  #src
println("="^92)  #src
