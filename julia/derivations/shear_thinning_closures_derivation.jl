# # Shear-Thinning Drops: Closures
#
# The companion page *Shear-Thinning Drops* derives the coupled model without
# approximation beyond small amplitude and axisymmetry, and closes with the
# complete statement of it. That model is closed but expensive: at every
# instant it asks for the solution of a banded system of radial boundary-value
# problems, coupled to a viscosity field computed from that same solution, with
# neither available in closed form.
#
# This page is the descent from there to something a solver can run. Each
# section adopts one further assumption, says what it buys, and prices what it
# costs; the last section covers the numerical realisation, which is a
# different kind of concession and is kept separate for that reason. The
# modelling sections are ordered so that each is independent of the ones after
# it: any of them can be dropped and the remainder still stands.

using Symbolics, QuadGK, SpecialFunctions, DropSolver  #src
using Printf  #src

function symbolic_zero(expr, vars; npts::Int=12, tol::Float64=1e-10)  #src
    s = Symbolics.simplify(Symbolics.expand(expr))  #src
    if s isa Number && iszero(s)  #src
        return true, :symbolic  #src
    end  #src
    worst = 0.0  #src
    for k in 1:npts  #src
        sub = Dict(v => 0.31 + 0.17i + 0.09k + 0.011i*k for (i, v) in enumerate(vars))  #src
        val = Symbolics.value(Symbolics.substitute(s, sub))  #src
        val isa Number || return false, :unresolved  #src
        worst = max(worst, abs(Float64(val)))  #src
    end  #src
    worst < tol, :numeric  #src
end  #src

println("="^78)  #src
println("SHEAR-THINNING DROPS: CLOSURES")  #src
println("="^78)  #src

# ## The quasi-static interior
#
# **Assumption.** The interior velocity field is in equilibrium with the
# instantaneous surface motion.
#
# The model page leaves the interior problem as a parabolic evolution,
#
# ```math
# \partial_t\,\mathcal D_l[U_l]
#   = \mathrm{Oh}\sum_{l''}\mathcal R_{l l''}[U_{l''};\eta] ,
# ```
#
# which has to be marched alongside the surface amplitudes. Replacing it by the
# statement that the interior has already relaxed -- dropping ``\partial_t`` and
# solving for a single decay rate ``\sigma`` at each instant -- turns it into
# the eigenproblem
#
# ```math
# q^2\,\mathcal D_l[U_l] + \mathcal R_l[U_l;\eta] = 0,
# \qquad q^2=\sigma/\mathrm{Oh} ,
# ```
#
# one per mode, of exactly the form solved on the Newtonian pages. This is what
# makes ``\lambda_l`` and ``\omega_l^2`` meaningful as *instantaneous*
# coefficients at all: without it there are no eigenvalues, only a field.
#
# **What it costs.** The substitution is exact only if the interior relaxes
# fast compared with the surface motion, i.e. if the viscous diffusion time
# across the drop is short against the oscillation period. Those are
# ``\mathrm{Oh}^{-1}`` and ``O(1)`` respectively in these units, so the
# assumption is a **high-Ohnesorge** one, and it degrades exactly where the
# drop is least damped and rings longest. It is not a small-parameter expansion
# with a leading correction; it is a change in the type of the interior
# equation, from parabolic to algebraic, and the discarded term is
# ``\partial_t\mathcal D_l[U_l]`` in full.
#
# Every section below is stated within this assumption.

# ## The temporal closure
#
# **Assumption.** How ``\eta``'s time dependence is reduced.
#
# ``\eta`` depends on ``t`` through ``\dot\gamma(t)``, which oscillates with
# the mode. Three inequivalent choices are available, and the production code
# takes the first:
#
# | choice | what it keeps | status |
# |:--|:--|:--|
# | (a) instantaneous ``\eta(t)`` | everything | **quasi-static; unjustified** |
# | (b) period-averaged | the ``m=0`` channel | justified by the lemma below |
# | (c) Floquet | ``m=0`` and ``m=2`` | most faithful, most work |
#
# ### What the modulation actually does: the amplitude equation
#
# Take one mode in free decay and carry the time-varying damping through an
# averaging calculation. Write ``\phi=\omega t``.
# By the period-``\pi`` lemma proved below, ``\lambda`` has **only even
# harmonics** in ``\phi``:
#
# ```math
# \lambda(\phi) \;=\; \bar\lambda
#   \;+\; \sum_{k\ge1}\bigl[\alpha_k\cos 2k\phi + \beta_k\sin 2k\phi\bigr].
# ```
#
# Put ``A=a(t)\cos\theta``, ``\theta=\omega t+\psi``, with ``a`` and ``\psi``
# slow, into ``\ddot A+2\lambda(t)\dot A+\omega^2A=0``. Krylov--Bogoliubov
# averaging gives ``\dot a=-2a\,\langle\lambda\sin^2\theta\rangle``, and with
# ``\sin^2\theta=\tfrac12[1-\cos(2\phi+2\psi)]`` every harmonic ``k\ge2`` is
# orthogonal to ``\cos(2\phi+2\psi)`` and drops out. What survives is
#
# ```math
# \boxed{\;
# \frac{da}{dt} \;=\; -\,\bar\lambda\,a
#   \;+\; \frac{a}{2}\bigl(\alpha_1\cos 2\psi - \beta_1\sin 2\psi\bigr) \;}
# ```
#
# Three things follow.
#
# **The period-average is the whole secular effect.** ``\bar\lambda`` -- the
# ``m=0`` harmonic -- is the only term that survives without a phase condition.
# That is what justifies closure (b), and it justifies it *quantitatively*
# rather than by appeal to plausibility.
#
# **The ``m=2`` harmonic is a parametric term, and it sits exactly on
# resonance.** It enters multiplied by ``\cos2\psi``, so its sign depends on the
# phase -- and the period-``\pi`` lemma places it at ``2\omega``, precisely the
# principal parametric resonance of an oscillator at ``\omega``. The lemma is
# therefore not a reassurance about the second harmonic: it identifies a
# resonance rather than ruling one out.
#
# **A bound settles when that matters.** Since ``|\alpha_1\cos2\psi-\beta_1\sin2\psi|
# \le\sqrt{\alpha_1^2+\beta_1^2}``, the effective decay rate is confined to
#
# ```math
# \bar\lambda - \tfrac12\sqrt{\alpha_1^2+\beta_1^2}
# \;\le\; \lambda_{\mathrm{eff}} \;\le\;
# \bar\lambda + \tfrac12\sqrt{\alpha_1^2+\beta_1^2},
# ```
#
# so a sufficient condition for the amplitude to decay **at every phase** is
#
# ```math
# \sqrt{\alpha_1^2+\beta_1^2} \;<\; 2\bar\lambda .
# ```
#
# If that holds, discarding the second harmonic changes the decay rate by a
# bounded amount and cannot change its sign; if it fails, there are phases at
# which the mode is parametrically pumped, and closure (b) is not merely
# approximate but qualitatively wrong. The criterion is checkable for a given
# fluid and amplitude, and it is the right thing to evaluate before trusting a
# period-averaged model.
#
let  #src
    avg(f) = quadgk(ph -> f(ph), 0, 2pi; rtol=1e-12)[1] / (2pi)  #src
    worst = 0.0  #src
    for ps in (0.0, 0.7, 1.9, 3.3, 5.1)  #src
        ## the m=0 harmonic contributes 1/2, independent of phase  #src
        worst = max(worst, abs(avg(ph -> sin(ph + ps)^2) - 0.5))  #src
        ## the m=2 harmonics contribute -cos(2psi)/4 and +sin(2psi)/4  #src
        worst = max(worst, abs(avg(ph -> cos(2ph) * sin(ph + ps)^2) + cos(2ps)/4))  #src
        worst = max(worst, abs(avg(ph -> sin(2ph) * sin(ph + ps)^2) - sin(2ps)/4))  #src
        ## every higher even harmonic averages away at this order  #src
        for k in 2:5  #src
            worst = max(worst, abs(avg(ph -> cos(2k*ph) * sin(ph + ps)^2)))  #src
            worst = max(worst, abs(avg(ph -> sin(2k*ph) * sin(ph + ps)^2)))  #src
        end  #src
    end  #src
    @assert worst < 1e-10 "the averaging coefficients in the amplitude equation are wrong ($worst)"  #src
    println("  ASSERTION 8b OK: averaging da/dt = -2a<lambda sin^2 theta> gives")  #src
    println("    da/dt = -lbar*a + (a/2)(alpha_1 cos 2psi - beta_1 sin 2psi),")  #src
    println("    with every harmonic k>=2 averaging to zero (worst residual $(round(worst, sigdigits=2))).")  #src
    println("    => the period-average carries the secular decay; the m=2 harmonic is")  #src
    println("       parametric, bounded by sqrt(alpha_1^2+beta_1^2)/2.")  #src
end  #src

# ### Why (a) is not licensed
#
# ``\lambda_l`` and ``\omega_l^2`` are *eigenvalues of a boundary-value
# problem*. They are not pointwise functions of the instantaneous state.
# Substituting a time-varying ``\eta`` into them is a quasi-static
# approximation, and it needs ``\eta`` to vary slowly compared with the
# mode's own period.
#
# It does the opposite. Here is the reason, and it is exact:
#
# > **Period-``\pi`` lemma.** Let ``S=\sqrt{2\bm e\!:\!\bm e}`` be built from any
# > single-frequency oscillating axisymmetric field. Each component obeys
# > ``\mathrm{Re}(e_{ij}e^{-i\phi})^2=\tfrac12|e_{ij}|^2[1+\cos(2\phi-2\arg e_{ij})]``,
# > manifestly invariant under ``\phi\to\phi+\pi``. Hence ``S^2``, ``S``, and
# > *any function of* ``S`` -- including ``\eta`` -- have period ``\pi`` in the
# > oscillation phase, not ``2\pi``.
#
# So ``\eta`` is modulated at **exactly twice** the mode frequency: the
# fastest modulation the problem admits. The quasi-static ratio is ``O(1)``
# by construction, at every amplitude. There is no limit in which choice
# (a) becomes valid.
#
# Two corollaries worth stating:
#
# * A period-``\pi`` function has **identically zero** content at odd
#   harmonics. In particular there is no forcing at the mode's own
#   frequency -- so the leading effect is carried entirely by the ``m=0``
#   (period-averaged) channel. That is choice (b), and the lemma is its
#   justification.
# * The surviving ``m=2`` channel modulates the damping at ``2\omega``,
#   which is the *principal parametric resonance* condition. Choice (b)
#   discards it; whether that is safe is a stability question this file
#   does not answer.
#
# Fourier-decomposing ``S`` over a full ``2\pi`` period, for a strain field
# with four independent complex components, gives
#
# | harmonic | ``m=0`` | ``m=1`` | ``m=2`` |
# |:--|--:|--:|--:|
# | coefficient | ``2.2374`` | ``<5\times10^{-16}`` | ``0.2781`` |
#
# The fundamental is zero to machine precision -- not small, *absent* -- while
# the mean and the second harmonic are both ``O(1)``. There is nothing for a
# resonant solvability condition to act on.
#
# It is tempting to read the criterion off this table, and it does not work.
# The criterion is on the Fourier harmonics of ``\lambda(\phi)``; these are the
# harmonics of ``S(\phi)``. Between the two sit a saturating constitutive law
# ``S\mapsto\eta(S)`` and an eigenvalue problem ``\eta\mapsto\lambda`` which, as
# stated above, is not a pointwise function of the state. Treating the two sets
# of coefficients as proportional is precisely the substitution this step exists
# to rule out.
#
# So ``\alpha_1`` and ``\beta_1`` are **not evaluated here.** Obtaining them
# means sweeping ``\phi`` through a cycle, evaluating ``\lambda`` from the
# eigenvalue solver at ``\eta(S(\phi))`` at each phase, and Fourier-decomposing
# *that*. Until that is done the criterion stands as a criterion and the
# stability of discarding the ``m=2`` channel is open -- which is how the
# closing section records it.
#
#
let  #src
    ## An arbitrary four-component complex field. That is the right test for a  #src
    ## LEMMA -- an identity must hold for any such field -- but it is not a  #src
    ## physical state, and nothing about the model may be priced from it.  #src
    e_amp = [1.3 + 0.7im, -0.4 + 1.1im, 0.9 - 0.2im, 0.5 + 0.3im]  #src
    Sof(φ) = begin  #src
        r = [real(z*(cos(φ) - im*sin(φ))) for z in e_amp]  #src
        sqrt(2*(r[1]^2 + r[2]^2 + r[3]^2 + 2*r[4]^2))  #src
    end  #src
    worst = maximum(abs(Sof(φ) - Sof(φ + π)) for φ in range(0, 2π; length=257))  #src
    @assert worst < 1e-13 "S is not period-pi; the lemma is false as stated"  #src
    N = 4096  #src
    φs = range(0, 2π; length=N+1)[1:end-1]  #src
    c1 = sum(Sof(φ)*cos(φ) for φ in φs)*2/N  #src
    s1 = sum(Sof(φ)*sin(φ) for φ in φs)*2/N  #src
    c2 = sum(Sof(φ)*cos(2φ) for φ in φs)*2/N  #src
    c0 = sum(Sof(φ) for φ in φs)/N  #src
    @assert abs(c1) < 1e-12 && abs(s1) < 1e-12 "S has content at the mode's OWN frequency"  #src
    @assert abs(c2) > 1e-3 "the m=2 parametric channel vanished; expected it to survive"  #src
    println("  ASSERTION 8 OK: S is period-pi to $(round(worst, sigdigits=2)); its Fourier")  #src
    println("    content at m=1 is zero to machine precision, while m=0 and m=2 survive.")  #src
    println("    (A lemma check on an arbitrary field -- not a physical state.)")  #src
end  #src

# ## Truncating the viscosity spectrum
#
# **Assumption.** ``\eta``'s Legendre content above ``l'=L_\eta`` is negligible.
#
# **What it drops.** Coupling between modes further apart than ``L_\eta``. The
# matrix becomes banded and the cost of applying it falls from ``O(M^2)`` to
# ``O(M L_\eta)``.
#
# **Error.** Exactly the discarded coupling, which is a measurable number. It
# has been measured, and the measurement rules this rung out.
#
# ### Measurement, and the choice of error norm
#
# The natural-looking metric, "what fraction of ``\sum_{l'}|\eta_{l'}|^2``
# lies below ``L_\eta``", says ``L_\eta=2\ldots5`` suffices everywhere. That
# metric is **misleading**. What controls the banding error is the *summed
# magnitude* of the discarded coefficients,
# ``\mathcal T(L)=\sum_{l'>L}|\eta_{l'}|/|\eta_0|``, because those terms enter the
# coupling matrix additively, not in quadrature. The spectrum turns out to
# have a long, slowly-decaying plateau carrying little *power* but large
# *summed magnitude*.
#
# Measured on Reid's actual viscous eigenmodes (volume-averaged,
# ``\mathrm{Oh}=0.2``, in the saturated thinned regime), the smallest
# ``L_\eta`` holding ``T_1<10^{-2}``:
#
# | modal spectrum | ``L_\eta`` for 1% |
# |:--|:--|
# | single mode ``l=2`` | 8 |
# | ``\propto l^{-2}``, ``M=30`` | 32 |
# | ``\propto l^{-2}``, ``M=50`` | 50 |
# | ``\propto l^{-1}``, ``M=30`` | 126 |
# | real solver state, ``M=30`` | 133 |
# | real solver state, ``M=20`` | 139 |
#
# > **For algebraically decaying modal spectra and for real solver states,
# > ``L_\eta\gtrsim M``, reaching ``L_\eta\gg M`` in the latter. The coupling
# > matrix is effectively dense, and this rung buys nothing.**
#
# The qualifier is doing work. A spectrum that decays *exponentially* -- say
# ``\dot A_l\propto0.85^{\,l}`` -- keeps ``L_\eta\approx16`` even with fifty
# modes nominally active, because so little amplitude reaches high ``l`` that
# the viscosity field stays smooth: its contrast is only ``\approx34\times``
# rather than the ``100\times`` a real state produces. Banding would work for
# such a fluid. It does not work here because the modal spectra this solver
# actually produces decay algebraically, not exponentially.
#
# At ``M=30`` you would need a half-bandwidth of ``\approx133`` to keep 1% --
# i.e. a full matrix. Banding at ``L_\eta=8`` leaves a 40% error even under
# the forgiving RMS norm.
#
# ### Why the spectrum decays so slowly
#
# The mechanism is worth stating because it is physical, not numerical.
# ``L_\eta`` tracks the *contrast* of the ``\eta`` field: configurations
# where ``\eta`` spans ``9\times`` are narrow, those spanning ``100\times``
# are broad. The contrast comes from **near-nodal points of the superposed
# strain field**, where ``S`` collapses and the fluid snaps back toward
# unthinned. That is a near-discontinuity in ``\mu``, and a
# near-discontinuity has a Legendre spectrum that decays only
# *algebraically*. With one mode active there are few such nodes and the
# field is smooth; with a dozen modes beating against each other the nodal
# set is dense.
#
# **Caveat.** Coefficients beyond ``l'=140`` were not measured, so every
# entry above is a *lower* bound on the discarded coupling.
#
# The measurement itself, including the high-``l`` eigenfunction evaluator it
# rests on (checked to ``5\times10^{-15}`` against Reid's tangential-stress
# boundary condition), is carried out in
# [Angular Bandwidth of Viscosity](eta_spectrum_derivation.md), which reproduces
# four of the six rows above on every CI run and the remaining two under
# `ETA_SPECTRUM_FULL=1`.
#
# **How to undo it.** Raise ``L_\eta``. The structure does not change -- only
# the bandwidth. But since ``L_\eta\gtrsim M`` is required for the spectra this
# solver produces, the practical move is to skip this rung and either keep the
# coupled system dense or drop to a spherically symmetric viscosity.

# ## A spherically symmetric viscosity
#
# **Assumption.** ``L_\eta=0``, i.e. ``\eta=\eta(r,t)`` with no angular
# structure.
#
# **Why this rung is special.** A spherically symmetric ``\eta`` commutes
# with the angular Laplacian. Every mode **decouples again**:
# ``G^{0}_{l l''}=\delta_{l l''}``, the matrices return to diagonal, and the
# entire architecture of the solver -- one independent oscillator per mode --
# is recovered intact.
#
# What you give up is *only* the closed form. The radial equation now has an
# ``r``-dependent coefficient, so it is no longer Bessel's equation and must
# be solved as a numerical two-point boundary-value problem, once per mode.
# `julia/src/reid.jl`'s continuation machinery already knows how to track
# eigenvalue branches through such a solve.
#
# !!! warning "This rung is cheap, and it is not a small correction"
#     The same measurement that killed the bandwidth truncation prices this one. At the physical
#     operating point, for a real multi-mode solver state, the anisotropic
#     coefficients are ``|\eta_1|/|\eta_0|\approx1.3`` and
#     ``|\eta_2|/|\eta_0|\approx0.2``. The angular structure of ``\eta`` is
#     *comparable to or larger than its mean*. Discarding all of it is a
#     **leading-order** modeling error, not a perturbative one.
#
#     Measured with a *single* active ``l=2`` mode the picture looks far
#     milder: the viscosity then varies by only ``1.1``--``1.2\times`` across
#     the drop, which would make spatial homogenization a ``\sim10\%``
#     effect. That estimate holds for one mode and not for a realistic state:
#     with a dozen modes beating, the nodal collapse described under
#     *Truncating the viscosity spectrum*
#     makes ``\eta`` span ``100\times``. Any
#     scalar-``\mathrm{Oh}_{\mathrm{eff}}`` model inherits this error.
#
# **This is the cheapest defensible rung, not an accurate one.** It keeps the
# radial structure at zero coupling cost.
#
# ### What Reid still gives you for free here
#
# Working through the boundary conditions with ``\eta=\eta(r)`` rather than
# constant ``\mu``:
#
# * **BC1 (kinematic)** contains no viscosity at all. **Unchanged.**
# * **BC2 (tangential stress).** The free-surface condition is
#   ``\tau_{r\theta}=2\eta\,e_{r\theta}=0``. Since ``\eta\ge\eta_\infty>0``
#   everywhere, this forces ``e_{r\theta}=0`` *regardless of whether ``\eta``
#   is constant*. **Unchanged**, and so is the whole
#   ``\tau_{r\theta}=0\Rightarrow\mathcal L_2[U]=0`` reduction already derived
#   in *The Viscous Drop: Reid (1960)*.
# * **BC3 (normal stress)** carries ``\eta`` multiplicatively, so it becomes
#   the *surface* value ``\eta_s=\eta(\dot\gamma|_{r=R})``. Same form, one
#   state-dependent coefficient.
#
# So the full scope of "adapt Reid for a shear-thinning fluid" is: **one
# extra term in the momentum equation, one coefficient in BC3, and nothing
# else.**

let  #src
    println("\nSTEP 7 (A7): eta = eta(r) restores exact diagonality")  #src
    Pl(l, m) = l == 0 ? one(m) : l == 1 ? m :  #src
        begin am, b = one(m), m; for n in 1:l-1; b, am = ((2n+1)*m*b - n*am)/(n+1), b; end; b end  #src
    worst_off = 0.0  #src
    worst_diag_err = 0.0  #src
    for l in 0:10, lpp in 0:10  #src
        g = (2l+1)/2 * quadgk(m -> Pl(l,m)*Pl(0,m)*Pl(lpp,m), -1, 1; rtol=1e-12)[1]  #src
        if l == lpp  #src
            worst_diag_err = max(worst_diag_err, abs(g - 1))  #src
        else  #src
            worst_off = max(worst_off, abs(g))  #src
        end  #src
    end  #src
    @assert worst_off < 1e-11 "l'=0 produced off-diagonal coupling ($worst_off)"  #src
    @assert worst_diag_err < 1e-11 "l'=0 diagonal entry is not 1 ($worst_diag_err)"  #src
    @printf("  largest off-diagonal entry of G^0 : %.2e  (must be 0)\n", worst_off)  #src
    @printf("  largest |diagonal - 1| of G^0     : %.2e  (must be 0)\n", worst_diag_err)  #src
    println("  ASSERTION 9 OK: G^0 = identity, so a spherically symmetric viscosity")  #src
    println("    leaves every mode DECOUPLED. Legendre structure fully intact.")  #src
    println("    Physical meaning of a failure: even a purely radial viscosity profile")  #src
    println("    would scatter energy between surface modes, which would contradict the")  #src
    println("    rotational symmetry of the problem.")  #src
end  #src

# ### The radial equation
#
# With ``\eta=\eta(r)`` the angular structure is untouched, so the problem stays
# diagonal in ``l`` exactly as the coupled system predicts, and each mode again
# reduces to a
# single radial equation. Carrying a variable ``\eta`` through the
# stream-function form of the momentum equation adds terms proportional to
# ``\eta'`` and ``\eta''`` and nothing else; setting ``\eta'=0`` returns Reid's
# operator ``\mathcal D_l(\mathcal D_l+q^2)U=0`` identically, which is how the
# construction is checked.
#
# What remains is a linear two-point boundary-value eigenproblem with variable
# coefficients: no longer Bessel, but entirely standard. The continuation
# machinery that tracks Reid's eigenvalue branches applies to it unchanged.
#
@variables rr tt qq hh1 hh2  #src
let  #src
    Drr = Differential(rr); Dtt = Differential(tt)  #src
    LPn(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    LPpn(l,m) = l==0 ? zero(m) : l*(m*LPn(l,m)-LPn(l-1,m))/(m^2-1)  #src
    angl(l) = sin(tt)^2 * LPpn(l, cos(tt)) / (l*(l+1))  #src
    function curlmom(l, U, H)  #src
        psi = U * angl(l)  #src
        u_r  =  Symbolics.expand_derivatives( Dtt(psi) / (rr^2*sin(tt)) )  #src
        u_th = -Symbolics.expand_derivatives( Drr(psi) / (rr*sin(tt)) )  #src
        e_rr = Symbolics.expand_derivatives(Drr(u_r))  #src
        e_tt = Symbolics.expand_derivatives(Dtt(u_th)/rr + u_r/rr)  #src
        e_pp = u_r/rr + u_th*cos(tt)/(sin(tt)*rr)  #src
        e_rt = Symbolics.expand_derivatives((Dtt(u_r)/rr + Drr(u_th) - u_th/rr)/2)  #src
        t_rr = 2H*e_rr; t_tt = 2H*e_tt; t_pp = 2H*e_pp; t_rt = 2H*e_rt  #src
        div_r = Symbolics.expand_derivatives(Drr(rr^2*t_rr)/rr^2 + Dtt(t_rt*sin(tt))/(rr*sin(tt)) - (t_tt+t_pp)/rr)  #src
        div_t = Symbolics.expand_derivatives(Drr(rr^2*t_rt)/rr^2 + Dtt(t_tt*sin(tt))/(rr*sin(tt)) + t_rt/rr - cos(tt)/(sin(tt)*rr)*t_pp)  #src
        Symbolics.expand_derivatives((Drr(rr*(qq*u_th + div_t)) - Dtt(qq*u_r + div_r))/rr)  #src
    end  #src
    Dln(e,l) = Symbolics.expand_derivatives(Drr(Drr(e))) - l*(l+1)*e/rr^2  #src
    reidop(U,l) = Symbolics.expand_derivatives(Dln(Dln(U,l) + qq*U, l))  #src
    function ev(e, rv, tv, qv, a1=0.0, a2=0.0)  #src
        sub = Dict(rr=>rv, tt=>tv, qq=>qv, hh1=>a1, hh2=>a2)  #src
        Float64(eval(Symbolics.toexpr(Symbolics.substitute(Symbolics.expand_derivatives(e), sub))))  #src
    end  #src
    PTS = ((0.43,0.7),(0.61,1.3),(0.82,2.1),(0.37,2.6),(0.95,0.44))  #src
    for l in (2,3,4)  #src
        U = rr^(l+1) + 0.7*rr^(l+3) - 0.3*rr^(l+5) + 0.11*rr^(l+2)  #src
        C = curlmom(l, U, 1)  #src
        T = -reidop(U,l)*angl(l)/(rr*sin(tt))  #src
        rats = [ev(C,rv,tv,3.7)/ev(T,rv,tv,3.7) for (rv,tv) in PTS]  #src
        spread = (maximum(rats)-minimum(rats))/abs(rats[1])  #src
        @assert spread < 1e-10 "curl of momentum is not proportional to Reid's operator (l=$l, spread=$spread)"  #src
        @assert all(abs(x - 1) < 1e-9 for x in rats) "proportionality constant is not 1 at l=$l"  #src
    end  #src
    for l in (2,3)  #src
        U = rr^(l+1) + 0.7*rr^(l+3)  #src
        extra = Symbolics.expand_derivatives(curlmom(l, U, 1 + hh1*rr + hh2*rr^2) - curlmom(l, U, 1))  #src
        @assert abs(ev(extra,0.61,1.3,3.7,0.0,0.0)) < 1e-10 "variable-eta terms must vanish when eta is constant"  #src
        @assert abs(ev(extra,0.61,1.3,3.7,1.0,0.0)) > 1e-3 "a linear eta(r) must produce a nonzero correction"  #src
        @assert abs(ev(extra,0.61,1.3,3.7,0.0,1.0)) > 1e-3 "a quadratic eta(r) must produce a nonzero correction"  #src
    end  #src
end  #src
# ## Freezing the radial profile
#
# **Assumption.** The radial variation of ``\eta`` is small enough to replace
# by a single number ``\eta_{\mathrm{eff}}``, taken as the
# dissipation-weighted average of ``\eta`` over the drop. Equivalently, the
# whole field is summarised by one effective Ohnesorge number
# ``\mathrm{Oh}_{\mathrm{eff}}=\eta_{\mathrm{eff}}/\sqrt{\rho T_1R}``.
#
# **What it buys.** The radial equation becomes constant-coefficient again,
# i.e. Bessel's equation, and Reid's closed-form characteristic equation
# returns verbatim -- evaluated at a shifted Ohnesorge number.
#
# **This is where the current production code sits**, combined with choice (a)
# of the temporal closure.
#
# Its error is the one measured in the previous section: for a realistic
# multi-mode state the
# viscosity spans a factor of order ``100`` across the drop, and the angular
# coefficients ``|\eta_1|/|\eta_0|\approx1.3`` and
# ``|\eta_2|/|\eta_0|\approx0.2`` are comparable to or larger than the mean.
# Collapsing all of that onto one number per mode is a leading-order
# approximation, not a perturbative one. Combined with the temporal closure of
# the instantaneous temporal closure, this rung carries two
# uncontrolled approximations at once -- which is worth knowing when reading its
# predictions.
#
# **How to undo it.** Keep the radial profile: the previous section.

# ## Specialising the constitutive law
#
# Cross is not a competing constitutive law; it is Carreau-Yasuda on the slice
# ``p=(n-1)/a=-1``, and it is how this repository's validation fluid is actually
# characterised. Establishing that takes a short detour through which exponents
# make ``\eta`` well behaved -- which is also where a persistent confusion about
# the amplitude non-analyticity can be cleared up, so that comes first.
#
# ### Two different non-analyticities, and only one of them is resolved
#
# Before going further it is worth separating two claims that sound alike and
# are not, because conflating them is the easiest mistake on this page.
#
# The linearisation section showed that ``\eta`` is not analytic in the
# **oscillation amplitude** ``\epsilon``: with ``\dot\gamma=O(\epsilon)`` and a
# non-integer ``a``, the correction scales as ``\epsilon^a``, larger than any
# linear term as ``\epsilon\to0``. **That claim stands, and nothing below
# repairs it.** It is
# why no small-amplitude expansion of this problem exists, and it is the reason
# the whole chain is built around evaluating Reid's relations exactly rather
# than perturbing them.
#
# What follows concerns a different question: whether ``\eta`` is a polynomial
# in the **strain-rate components**, which controls how its Legendre series
# behaves and therefore whether the coupling matrix is finite. The answer there
# is favourable, and it is favourable *regardless* of the amplitude
# non-analyticity, because the two involve different variables. A reader who
# takes the good news below as an answer to the amplitude non-analyticity has
# read it wrongly.
#
# ### Which exponents make ``\eta`` a polynomial
#
# We can now say precisely what is special about particular exponents, and
# it is a statement about **structure**, not about accuracy.
#
# The coupling matrix is exactly banded, and exactly truncatable, if and
# only if ``\eta`` is a **polynomial** in the components of ``\bm e``. Track
# the two nestings:
#
# ```math
# \eta = \eta_\infty + \Delta\eta\,(1+X)^{p},
# \qquad X = (\lambda_c\dot\gamma)^a,
# \qquad p = \frac{n-1}{a},
# \qquad \Delta\eta \equiv \eta_0-\eta_\infty .
# ```
#
# * ``X`` is a polynomial in ``\bm e`` **iff** ``a`` is a non-negative even
#   integer, because ``\dot\gamma^a=(\dot\gamma^2)^{a/2}`` and
#   ``\dot\gamma^2=2\bm e\!:\!\bm e`` is quadratic.
# * ``(1+X)^p`` is a polynomial in ``X`` **iff** ``p`` is a non-negative
#   integer.
#
# For a genuinely shear-thinning fluid ``n<1``, so ``p<0`` and the second
# condition never holds: ``\eta`` is not a polynomial in the strain rate, and
# its Legendre series does not terminate.
#
# That sounds like a barrier and is not one, because the coupling does not see
# the whole series. The Gaunt rule requires ``l'\le l+l''``, and the shape
# expansion stops at ``l=M``, so **every viscosity harmonic above
# ``l'=2M`` is orthogonal to every product ``P_lP_{l''}`` and cannot couple
# any pair of modes.** It is not small; it is absent.
#
# > The mode-coupling matrix is determined **exactly** by the first ``2M+1``
# > Legendre harmonics of the viscosity. Truncating there is not an
# > approximation -- everything above is annihilated by the angular integral.
#
# So the infinite series is a property of ``\eta``, not a limitation of the
# model. For ``M=50`` the coupling is fixed by 101 numbers per radius; for
# ``M=90``, by 181. Whether the matrix can be further *banded* -- kept to
# ``l'\le L_\eta`` with ``L_\eta\ll 2M``, for speed -- is a separate and
# genuinely empirical question, and *Truncating the viscosity spectrum*
# measures it.
#
# It is worth being precise about why the classical ``a=2`` theory closes.
# Not because ``a=2`` makes ``\eta`` polynomial -- it does not -- but because
# that theory additionally *truncates* ``(1+X)^p\approx 1+pX`` at first
# order, and ``X`` alone is polynomial when ``a=2``. The finiteness comes
# from the perturbative truncation, and ``a=2`` merely keeps the truncated
# object polynomial.
#
# ### The exception that matters: ``p=-1``
#
# There is one value of ``p`` that is special without being a non-negative
# integer, and it is the one the validation fluid has. If ``p=-1``, i.e.
# ``n=1-a``, then
#
# ```math
# \eta = \eta_\infty + \frac{\Delta\eta}{1+X}
# \qquad\Longleftrightarrow\qquad
# (1+X)\,\eta = (1+X)\,\eta_\infty + \Delta\eta .
# ```
#
# The constitutive law becomes an **algebraic constraint** that is *linear*
# in ``\eta`` and *linear* in ``X``, with no fractional power of anything
# left on the outside. This is the Cross model, and it is what makes this step
# is a genuine simplification rather than a relabelling.

let  #src
    println("\nSTEP 4b: which exponents make eta polynomial (hence EXACTLY banded)")  #src
    @variables Xs etainf Deta  #src
    lhs = (1 + Xs)*(etainf + Deta/(1 + Xs))  #src
    rhs = (1 + Xs)*etainf + Deta  #src
    ok, how = symbolic_zero(Symbolics.simplify(lhs - rhs), [Xs, etainf, Deta])  #src
    @assert ok "the Cross law failed to reduce to a polynomial constraint"  #src
    println("  ASSERTION 6 OK ($how): for p = (n-1)/a = -1 the constitutive law is")  #src
    println("    equivalent to the POLYNOMIAL constraint (1+X)*eta = (1+X)*eta_inf + Deta.")  #src
    println("    No fractional power remains on the outside. This is the Cross model.")  #src

    println("\n  Polynomiality bookkeeping (X = (lambda_c*gammadot)^a):")  #src
    for a in (0.5, 0.743, 1.0, 2.0, 4.0)  #src
        xpoly = (a >= 0) && isinteger(a/2)  #src
        @printf("    a=%-6.3f  X polynomial in e? %-5s\n", a, xpoly)  #src
    end  #src
    for (nn, aa) in ((0.257, 0.743), (0.5, 2.0), (0.0, 1.0))  #src
        pp = (nn - 1)/aa  #src
        @printf("    n=%-6.3f a=%-6.3f -> p=%-7.3f  (1+X)^p polynomial in X? %s\n",  #src
                nn, aa, pp, (pp >= 0 && isinteger(pp)))  #src
    end  #src
    for nn in (0.0, 0.257, 0.5, 0.9), aa in (0.5, 0.743, 1.0, 2.0, 4.0)  #src
        @assert (nn - 1)/aa < 0 "p = (n-1)/a must be negative for every shear-thinning fluid"  #src
    end  #src
    println("  ASSERTION 7 OK: p<0 for every shear-thinning fluid, so eta's Legendre")  #src
    println("    series does not terminate -- but the coupling is still exact, see below.")  #src
    println("    Physical meaning: no shear-thinning drop model can claim an EXACT")  #src
    println("    finite mode-coupling bandwidth; the truncation must carry an error bar.")  #src
end  #src


#
# Cross is not a different model; it is Carreau-Yasuda restricted to
# ``p=(n-1)/a=-1``, i.e. ``n=1-a``:
#
# ```math
# \eta = \eta_\infty + \frac{\eta_0-\eta_\infty}{1+(K\dot\gamma)^m},
# \qquad K=\lambda_c,\quad m=a .
# ```
#
# Two reasons this rung is worth taking, both structural:
#
# 1. By the argument above the constitutive law becomes an **algebraic constraint**
#    ``(1+X)\eta=(1+X)\eta_\infty+\Delta\eta`` -- linear in ``\eta``, linear
#    in ``X``, no outer fractional power. Symbolic manipulation becomes
#    tractable in a way it is not for general ``p``.
# 2. It is how the validation fluid is *actually characterised*: the fitting
#    script obtains ``\eta_0,\eta_\infty,K,m`` and only then converts to
#    Carreau-Yasuda parameters. Deriving Cross directly removes that
#    conversion, and with it a genuine trap -- under the Cross mapping the
#    Carreau-Yasuda exponent ``(1-n)/a`` is *exactly* 1, which is easy to
#    confuse with the thinning amplitude
#    ``\Delta=(\eta_0-\eta_\infty)/\eta_0``. The two coincide only when
#    ``\Delta\approx1``, as it happens to for this fluid.
#
# What Cross does **not** buy: ``X=(K\dot\gamma)^m`` still carries a
# fractional power of the *field* when ``m`` is not an even integer, so the
# viscosity spectrum is still formally infinite and ``L_\eta`` is still an
# empirical truncation. Cross simplifies the outer algebra, not the angular
# problem.

let  #src
    println("\nSTEP 9: Cross as the p = -1 slice of Carreau-Yasuda")  #src
    @variables gds Ks ms eta0s etainfs  #src
    a_sym = ms  #src
    n_sym = 1 - ms  #src
    cy = etainfs + (eta0s - etainfs)*(1 + (Ks*gds)^a_sym)^((n_sym - 1)/a_sym)  #src
    cross = etainfs + (eta0s - etainfs)/(1 + (Ks*gds)^ms)  #src
    ok, how = symbolic_zero(Symbolics.simplify(cy - cross), [gds, Ks, ms, eta0s, etainfs])  #src
    @assert ok "Carreau-Yasuda at n = 1-a did not reduce to Cross"  #src
    println("  ASSERTION 10 OK ($how): CY with n = 1 - a is IDENTICALLY the Cross model.")  #src
    println("    Cross is a slice of the general law, so the chain loses nothing by")  #src
    println("    descending to it, and the general case is recovered by freeing n.")  #src

    ETA_INF = 0.0037320997942061666  #src
    ETA_0   = 8.433817577956766  #src
    M_CROSS = 0.7430524574330837  #src
    K_CROSS = 18.48081673111359        ## s; the Cross time constant, from the same fit  #src
    Δ = (ETA_0 - ETA_INF)/ETA_0  #src
    eps_inf = ETA_INF/ETA_0  #src
    ## The claim is that the Cross law and the Carreau-Yasuda form the solver  #src
    ## implements are the SAME curve under lambda_c=K, a=m, eps_ST=1 -- not that  #src
    ## m/m equals one. Compare the two curves directly, across the shear rates  #src
    ## an impact actually samples.  #src
    cross(g)     = eps_inf + (1 - eps_inf)/(1 + (K_CROSS*g)^M_CROSS)  #src
    cy(g, eps)   = eps_inf + (1 - eps_inf)*(1 + (K_CROSS*g)^M_CROSS)^(-eps)  #src
    gs = exp10.(range(-3, 5; length=40))  #src
    worst_ok  = maximum(abs(cy(g, 1.0) - cross(g))/cross(g) for g in gs)  #src
    @assert worst_ok < 1e-12 "Cross is not the eps_ST = 1 slice of Carreau-Yasuda ($worst_ok)"  #src
    ## And the slot must be load-bearing: putting Delta there instead of 1 has  #src
    ## to change the curve materially, or the distinction would not matter.  #src
    worst_bad = maximum(abs(cy(g, Δ) - cross(g))/cross(g) for g in gs)  #src
    @assert worst_bad > 1e-3 "using Delta as the exponent should visibly change the curve"  #src
    println("  ASSERTION 11 OK: the Cross law is exactly the eps_ST = 1 slice of the")  #src
    println("    Carreau-Yasuda form the solver implements (agreement $(round(worst_ok, sigdigits=2)) over")  #src
    println("    eight decades of shear rate). Substituting the thinning amplitude")  #src
    println("    Delta = $(round(Δ, digits=6)) into the exponent slot instead shifts the curve by")  #src
    println("    up to $(round(100*worst_bad, sigdigits=3))%, so the two are not interchangeable even though")  #src
    println("    Delta is close to 1 for this fluid.")  #src
end  #src

# ## The Newtonian floor
#
# The bottom of the chain must reproduce the Newtonian result the solver
# already implements. Two checks against the *running solver*, rather than
# against this page's own algebra: Reid's exact finite-Ohnesorge coefficients
# must emerge when the viscosity stops depending on the flow, and they must
# reduce to Lamb as ``\mathrm{Oh}\to0``.
#
# Both are checked here. The first is the statement that switching the thinning
# off -- ``\lambda_c=0``, a fluid whose viscosity does not depend on the flow --
# must return ``\mathrm{Oh}_{\mathrm{eff}}=\mathrm{Oh}`` for every mode and
# recover Reid's coefficients through the shear-thinning code path rather than
# the Newtonian one. It does, to within the interpolation error of the
# tabulated lookup that path uses. That is the check this page most needs,
# because it is the one that exercises the generalisation itself; the Lamb
# limit below tests only the Newtonian floor underneath it.
#
# | ``\mathrm{Oh}`` | ``|\lambda_{\rm Reid}-\lambda_{\rm Lamb}|/\lambda_{\rm Lamb}`` |
# |:--|--:|
# | ``10^{-2}`` | ``5.17\times10^{-2}`` |
# | ``10^{-3}`` | ``1.60\times10^{-2}`` |
# | ``10^{-4}`` | ``5.05\times10^{-3}`` |
#
# Lamb is the ``\mathrm{Oh}\to0`` *asymptote*, not an approximation valid at
# fixed ``\mathrm{Oh}``, so the meaningful statement is that the discrepancy
# falls monotonically and reaches ``0.5\%``. At a working
# ``\mathrm{Oh}=0.05`` the live solver returns
# ``\lambda_2=0.2187`` against Lamb's ``0.2500`` -- a ``13\%`` difference,
# which is what the `:reid` viscous model exists to remove.

let  #src
    println("\nSTEP 10: the Newtonian floor, checked against the running solver")  #src
    errs = Float64[]  #src
    for Oh in (1e-2, 1e-3, 1e-4)  #src
        lam, om2, resid = reid_lambda_omega2(Oh, 2)  #src
        @assert resid < 1e-8 "Reid characteristic residual too large at Oh=$Oh"  #src
        lamb = Oh*(2 - 1)*(2*2 + 1)  #src
        push!(errs, abs(lam - lamb)/lamb)  #src
        @printf("  Oh=%-8.0e |lambda_Reid - lambda_Lamb|/lambda_Lamb = %.3e\n", Oh, errs[end])  #src
    end  #src
    ## Lamb is the Oh->0 ASYMPTOTE, not an approximation valid at fixed Oh, so  #src
    ## the meaningful test is that the discrepancy shrinks monotonically as Oh  #src
    ## falls and is small at the smallest Oh -- the same test test_reid.jl makes.  #src
    @assert all(errs[k+1] < errs[k] for k in 1:length(errs)-1) "Reid is not converging to Lamb"  #src
    @assert errs[end] < 0.01 "Reid-vs-Lamb discrepancy still >1% at Oh=1e-4"  #src
    lam_v, om2_v = drop_viscous_coeffs(6, 0.05, :reid)  #src
    lam_l, om2_l = drop_viscous_coeffs(6, 0.05, :lamb)  #src
    @assert all(isfinite, lam_v) && all(>(0), lam_v) "solver returned unusable Reid coefficients"  #src
    @printf("  live drop_viscous_coeffs(M=6, Oh=0.05, :reid) lambda_2 = %.6f (Lamb: %.6f)\n",  #src
            lam_v[1], lam_l[1])  #src
    println("  ASSERTION 12 OK: the bottom of the chain is the Reid/Lamb Newtonian")  #src
    println("    limit that julia/src/reid.jl already implements and validates.")  #src
    println("    Physical meaning of a failure: the hierarchy would not contain the")  #src
    println("    Newtonian case it claims to generalise.")  #src
end  #src

let  #src
    M = 5  #src
    Oh0 = 0.05  #src
    stx = STExactParams(M, Oh0, 0.0, 2.0, 0.5; viscous=:reid)  #src
    Adot = [0.05, -0.02, 0.031, 0.004]  #src
    oh_eff = oh_eff_all_coupled(stx, Oh0, Adot)  #src
    @assert all(o -> isapprox(o, Oh0; rtol=1e-12), oh_eff) "Oh_eff must equal Oh0 when lambda_c = 0"  #src
    lam_st, om2_st = lambda_omega2_from_oh_eff(stx, oh_eff)  #src
    lam_n, om2_n = drop_viscous_coeffs(M, Oh0, :reid)  #src
    worst = maximum(max(abs(lam_st[k]-lam_n[k])/lam_n[k], abs(om2_st[k]-om2_n[k])/om2_n[k])  #src
                    for k in eachindex(lam_n))  #src
    @assert worst < 5e-3 "the shear-thinning path does not reduce to Reid at lambda_c=0 ($worst)"  #src
    println("  ASSERTION 12b OK: with the thinning switched off (lambda_c = 0) the")  #src
    println("    shear-thinning path returns Oh_eff = Oh0 exactly, and its lambda_l,")  #src
    println("    omega_l^2 agree with the Newtonian Reid coefficients to")  #src
    println("    $(round(100*worst, sigdigits=2))% (tabulated-vs-exact interpolation error).")  #src
end  #src

# ## Numerical realisation
#
# The sections above are modelling concessions: each replaces the physics with
# different physics. What follows is of a different kind -- the model, as
# reduced, is still stated in continuous terms, and this is how those
# continuous objects are evaluated. Nothing here changes the model; it only
# changes how accurately the model is represented.
#
# ### Quadrature
#
# None of the integrals in the model has a closed form, because
# ``\eta(\dot\gamma)`` is not polynomial. They are evaluated on a **product
# Gauss-Legendre rule**: ``n_r`` nodes in ``x\in[0,1]`` and ``n_x`` nodes in
# ``\mu\in[-1,1]`` (`STExactParams`). At each node ``\dot\gamma`` is computed
# from the full superposition of active modes, ``\eta`` follows pointwise from
# the constitutive law, and a Legendre projection in ``\theta`` at each radius
# produces the ``\eta_{l'}(x)`` the matrices need.
#
# The **angular** rule has to scale with the truncation, and this is easy to
# get wrong. The integrand is a product of two strain bases, hence a polynomial
# of degree ``\le2M`` in ``\mu``, so an ``n_x``-node Gauss rule is exact only
# for ``n_x\ge M+\tfrac12``. A rule fixed at ``n_x=30`` is therefore exact at
# ``M=16`` and badly inexact at ``M=90``: the measured error in
# ``\mathrm{Oh}_{\mathrm{eff}}`` grows from ``1.6\times10^{-3}`` to
# ``8.4\times10^{-2}``. With the default ``n_x=\max(30,\,2M+2)`` the accuracy
# is independent of ``M``. `julia/test/test_oh_eff_quadrature.jl` asserts that
# independence rather than a fitted constant, so it still fails if the scaling
# is removed.
#
# The **radial** rule does not need to scale. Gauss nodes cluster towards the
# endpoints, which is where the ``x^{l-2}`` weight concentrates; ``n_r=20`` and
# ``n_r=400`` agree to six figures.
#
# ### Time stepping
#
# The reduced system is still quasi-linear -- the coefficients depend on the
# state through ``\eta`` -- so an implicit step would need the Jacobian of the
# coefficients as well as of the unknowns. It is closed instead by
# **extrapolation**: at the step from ``t_n`` to ``t_{n+1}`` the coefficients
# are evaluated not at the unknown state but at a second-order estimate of it,
#
# ```math
# \bm{\dot A}^{\,\ast} = (1+r)\,\bm{\dot A}_n - r\,\bm{\dot A}_{n-1},
# \qquad r=\frac{\Delta t_{n+1}}{\Delta t_n},
# ```
#
# exact whenever ``\bm{\dot A}`` is linear in ``t``, for non-uniform steps as
# well as uniform. Then ``\mathcal D(\bm{\dot A}^{\,\ast})`` is constant within
# the step, the system is linear in the unknowns, and the Jacobian is exact.
# This is an IMEX splitting: implicit in the modal amplitudes, explicit in the
# coefficients.
#
# The extrapolation order is not free. Holding the coefficients at the previous
# step instead -- ``\bm{\dot A}^{\,\ast}=\bm{\dot A}_n`` -- is a first-order
# splitting error inside an otherwise second-order BDF2 scheme, and it drags the
# observed order of the whole integration down to ``\approx1.3`` against a
# Newtonian control at ``2.0``. The linear extrapolation above restores it.
# `julia/test/test_convergence_order.jl` measures the order and holds it above
# ``1.5``.

# ## Summary of the chain
#
# The left column names the section; the model page supplies the first four
# rows, this page the rest.
#
# | model | assumption to get here | coupling structure |
# |:--|:--|:--|
# | exact free-surface generalized Newtonian | -- | -- |
# | linear in amplitude | ``\epsilon\ll1``; drops advection | ``\eta`` still fully nonlinear |
# | axisymmetric | axisymmetric forcing -- **exact** | ``Y_l^m\to P_l``, no ``m``-coupling ever |
# | poloidal + modal | change of variables | state is ``\{A_l\}``, ``l=2\ldots M`` |
# | **full coupled system** | none beyond the three above | banded interior BVP + dense surface matrices |
# | quasi-static interior | interior relaxed; parabolic ``\to`` eigenproblem | recovers ``\lambda_l``, ``\omega_l^2`` as instantaneous numbers |
# | temporal closure | instantaneous / period-averaged / Floquet | picks which ``\eta(t)`` channel survives |
# | truncate at ``L_\eta`` | discarded coupling small -- **measured false** | banded, but needs ``L_\eta\gtrsim M``: no saving |
# | ``\eta=\eta(x)`` | viscosity spherically symmetric -- **leading-order error** | **diagonal**; numerical radial BVP per mode |
# | ``\eta\to\eta_{\rm eff}`` | radial variation small | closed form -- **what the solver runs** |
# | ``\eta=\eta_0`` | no thinning | Reid |
# | ``\mathrm{Oh}\to0`` | inviscid | Lamb |
#
# A specific constitutive law is a slice through this table, not a rung of it,
# and can be taken at any row. The numerical realisation is not in the table at
# all: it is not a rung, because it does not change the model.
#
# ## What this file establishes
#
# **Free, and permanent.** A generalized Newtonian fluid cannot break
# axisymmetry, so Legendre polynomials suffice forever. BC1 carries no
# viscosity at all, and BC2 reduces to ``\mathcal L_2[U]=0`` for every fluid in
# the admissible class, because ``\eta\ge\eta_\infty>0`` cannot vanish. The
# normal-stress condition is where the rheology lands: once ``\eta`` varies
# with ``\theta`` its surface value is a field, and projecting it is what
# produces the coupling matrices. Under the spherically symmetric rung below
# that field collapses to the single number ``\eta_s``, which is the sense in
# which "BC3 changes by one coefficient" is true -- it is a property of that
# rung, not of the model.
#
# **Buys nothing.** Banding the coupling matrix at ``L_\eta``: the measured
# bandwidth is ``\gtrsim M`` for algebraically decaying spectra and ``\gg M``
# for real solver states, so the band is the whole matrix. (An exponentially
# decaying spectrum would band happily; this solver does not produce one.)
#
# **Priced, and not cheap.** Restricting to ``\eta=\eta(x)`` restores exact
# diagonality, but discards angular structure that is *comparable to or larger
# than* the mean viscosity. Parity decides where that error lands: odd
# harmonics vanish identically on the diagonal, so ``\eta_1`` -- the largest --
# does not corrupt each mode's own damping at all. It carries the entire
# mode-to-mode coupling instead, and setting it to zero does not shrink that
# coupling but removes it. The leading contamination of the diagonal is
# ``\eta_2``, at about a fifth of the mean.
#
# **The conclusion.** There is no cheap-and-accurate rung between the dense
# coupled system and a spherically symmetric viscosity. Everything implemented
# here so far sits below the latter and inherits at least its error. The
# options are to carry the coupled system and pay for it, or to quote the error
# alongside the result.
#
# **Still open.** Whether the coupled system is affordable at ``M\sim50`` (the
# Gaunt coefficients are geometry and precompute once; only the radial
# integrals of ``\eta_{l'}`` change per step); the eigenvalue problem for the
# variable-``\eta`` radial operator, which is a two-point boundary-value
# problem rather than a research question; and the ``m=2`` parametric channel,
# which the period-``\pi`` lemma places at exactly the principal resonance
# condition and which nothing here evaluates.
