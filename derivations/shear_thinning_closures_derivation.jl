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
# modelling sections are ordered so that each buys something the next one takes
# for granted, and the first of them -- the eigenmode closure -- is assumed by all
# the rest. Where a section depends on its predecessor it says so.

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

# ## The instantaneous-eigenmode closure
#
# **Assumption.** At each instant the interior flow is a single exponentially
# decaying mode of the interior operator, with the viscosity field held fixed.
#
# This is the assumption that lets the interior be eliminated, and it is worth
# stating with some care, because "quasi-static" -- the name it usually travels
# under -- suggests something different from what is actually done, and nothing
# in it is static.
#
# The model page leaves the interior amplitudes in the state, to be marched alongside
# the surface ones. Written in the strong form, that interior evolution is parabolic,
#
# ```math
# \partial_t\,\mathcal D_l[\psi_l]
#   = \mathrm{Oh}\sum_{m}\mathcal R_{l m}[\psi_{m};\eta] ,
# ```
#
# with ``\mathcal R_{lm}`` the ``(l,m)`` block of ``K_{ab}/\mathrm{Oh}``, as the model page's
# section *Where the interior operator went* identifies. The division is not cosmetic:
# ``K_{ab}`` already carries the Ohnesorge number, so leaving it in would count it twice here
# and would break ``q^2=\sigma/\mathrm{Oh}`` below. The variational
# form of the same content is the second-order system in all the coordinates at once;
# the parabolic form is written here because the elimination about to be performed is
# most easily seen on it. There are two different ways to get rid of that, and only the
# second is used:
#
# | | substitution | result |
# |:--|:--|:--|
# | (a) drop the time derivative | ``\partial_t\mathcal D_l[\psi_l]\to0`` | a Stokes interior: no oscillation, ``\sigma=0`` |
# | (b) replace it by an eigenvalue | ``\partial_t \psi_l\to-\sigma \psi_l`` | Reid's eigenproblem |
#
# Option (a) is what "quasi-static" sounds like, and it is not what anyone does:
# it deletes the inertia that makes the drop ring. Option (b) does not discard
# the term at all -- it *replaces* it, asserting that the interior's time
# dependence is a single exponential ``e^{-\sigma t}`` and then solving for
# ``\sigma``. Substituting it gives
#
# ```math
# q^2\,\mathcal D_l[\psi_l] + \mathcal R_l[\psi_l;\eta] = 0,
# \qquad q^2=\sigma/\mathrm{Oh} ,
# ```
#
# one per mode, of exactly the form solved on the Newtonian pages.
#
# **What it entails, precisely.** Three separate claims, worth separating
# because they fail in different places:
#
# 1. **One exponential, not a superposition.** The interior response to a given
#    surface motion is in general a sum over *all* decay rates admitted by the
#    operator. Keeping a single ``e^{-\sigma t}`` asserts the others have
#    already decayed.
# 2. **The viscosity is constant over the relaxation.** An eigenvalue only
#    exists for a time-independent operator. Since ``\eta`` varies on the same
#    timescale as the mode that produces it, this is an assertion about
#    timescales, not an identity.
# 3. **``\psi_l`` stops being a state variable.** It becomes an algebraic function
#    of the current ``\eta`` and ``l``. This is the whole point: it is what makes
#    ``\lambda_l`` and ``\omega_l^2`` *instantaneous numbers* rather than
#    histories, and so what lets the surface equation stay a second-order ODE.
#
# A fourth approximation is usually layered on top and is *not* part of this
# one: even granted 1--3, the resulting characteristic equation has infinitely
# many roots ``\sigma`` per mode, and ``\lambda_l``, ``\omega_l^2`` are built
# from **two** of them. That truncation is independent, and it is not exact even
# for a constant viscosity.
#
# **What it costs.** Claims 1 and 2 both require the interior to relax fast
# compared with the surface motion -- the viscous diffusion time across the drop
# short against the oscillation period. In these units those are
# ``\mathrm{Oh}^{-1}`` and ``O(1)``, so this is a **high-Ohnesorge** assumption,
# and it degrades exactly where the drop is least damped and rings longest. It is
# not a small-parameter expansion with a leading correction; it changes the type
# of the interior equation, from parabolic to algebraic.
#
# !!! note "This closure is optional"
#     Nothing forces it. If the interior is discretised and marched as the model
#     page states it -- ``\psi_l`` in the state, a banded system in ``x`` -- then
#     none of claims 1--3 is needed, no eigenproblem is ever formed, and the
#     two-root truncation disappears with it. The closure buys a smaller state
#     and the reuse of tabulated Newtonian coefficients; it is not a
#     prerequisite for having a solvable model.
#
# The sections below are stated within this closure, because that is where the
# implemented solver sits.

# ## From the interior to two numbers per mode
#
# The coefficient matrices below exist *only* once ``\psi_l`` has been
# eliminated, which is what the closure above does: the second-order oscillator
# form is a consequence of eliminating the interior, not a feature of the
# physics.
#
# ### Two routes, and they must agree
#
# With ``\psi_l`` eliminated the coordinates are the surface amplitudes alone, and
# the model page's variational statement collapses to constant matrices:
#
# ```math
# \bm M\ddot{\bm\zeta}+\bm C\dot{\bm\zeta}+\bm K\bm\zeta=\bm Q,
# \qquad
# \mathcal D^{(2)}=\bm M^{-1}\bm C,
# \qquad
# \mathcal D^{(1)}=\bm M^{-1}\bm K ,
# ```
# ```math
# M_{lm}=\int\bm u^{(l)}\!\cdot\!\bm u^{(m)}dV,
# \qquad
# C_{lm}=\mathrm{Oh}\!\int2\eta\;\bm e^{(l)}\!:\!\bm e^{(m)}dV,
# \qquad
# G_l\propto(l-1)(l+2) .
# ```
#
# That is the **energy route**. What follows is the **traction route**: project the
# normal-stress balance and collect the result. They are different calculations of
# the same matrices -- one integrates a dissipation over the volume, the other
# projects a stress at the surface -- and they agree only if the interior is solved
# consistently, so the agreement is a claim to be checked rather than a remark.
# The check below does that, and it is the reason both routes are kept: two
# independent assemblies agreeing is worth more than either alone.
#
# ### Carrying the projection through
#
# The equation being projected is the **normal-stress balance on the free
# surface**, linearised and evaluated at ``x=1``,
#
# ```math
# \Bigl[-p + \tau_{rr}\Bigr]_{x=1} \;=\; \gamma\,(\nabla\cdot\bm n)\Big|_{x=1} ,
# \qquad \tau_{rr}=2\eta\,e_{rr},
# ```
#
# which is the normal-stress condition BC3 of the Newtonian pages and the only one that
# produces an equation of motion -- BC1 fixes the kinematics and BC2 constrains
# the interior profile. Multiply it by ``P_l(\mu)`` and integrate over
# ``\mu\in[-1,1]``. Writing the surface motion as
# ``\dot\zeta=R\sum_{m}\dot\zeta_{m}P_{m}``, so that mode ``m`` enters with
# strength ``\dot\zeta_{m}``, the contribution of the viscous stress to the
# equation for ``\zeta_l`` is a double sum over the driving mode ``m`` and the
# viscosity harmonic ``k``, of terms
#
# ```math
# \sum_{m}\sum_{k} \dot\zeta_{m}\Bigl[\,
#   G^{k}_{l m}\!\!\int_0^1\!\bigl(\eta_{k}\mathcal L_{m}[u_{r,m}]
#                                      + \eta'_{k}u_{r,m}'\bigr)x^2dx
# \;+\;
#   H^{k}_{l m}\!\!\int_0^1\!\frac{\eta_{k}}{2x}
#      \Bigl(u_{\theta,m}'-\frac{u_{\theta,m}}{x}-\frac{u_{r,m}}{x}\Bigr)x^2dx
# \,\Bigr],
# ```
#
# where the two angular factors are
#
# ```math
# G^{k}_{l m}=\frac{2l+1}{2}\!\int_{-1}^{1}\!P_l\,P_{k}\,P_{m}\,d\mu ,
# \qquad
# H^{k}_{l m}=\frac{2l+1}{2}\!\int_{-1}^{1}\!P_l\,(1-\mu^2)\,P'_{k}\,P'_{m}\,d\mu ,
# ```
#
# and the radial factors are exactly the ones the three viscous contributions
# produce. From ``\eta\nabla^2\bm u``, the radial Laplacian of the driving mode's
# profile,
# ``\mathcal L_{m}[u_{r,m}]=u_{r,m}''+\tfrac{2}{x}u_{r,m}'-\tfrac{m(m+1)}{x^2}u_{r,m}``,
# carried by ``\eta_{k}`` itself. From the radial part of
# ``2(\nabla\eta)\cdot\bm e``, the radial strain amplitude ``u_{r,m}'``, carried by
# ``\eta'_{k}``. From its polar part, the ``e_{r\theta}`` amplitude
# ``\tfrac12(u_{\theta,m}'-u_{\theta,m}/x-u_{r,m}/x)``, carried by ``\eta_{k}/x`` -- and this
# is the one term that pairs with ``H`` rather than ``G``, because its angular
# factor carries the two derivatives.
#
# ``u_{r,m}`` and ``u_{\theta,m}`` are the radial profiles of ``u_r`` and
# ``u_\theta`` introduced with the modal expansion; both are determined by
# ``\psi_{m}``, hence by the interior problem of the previous section. They are
# named ``F`` and ``W`` so that the letter ``G`` belongs to the Gaunt
# coefficient alone.
#
# The second angular form is not a new object. Using
# ``(1-\mu^2)P'_l=\tfrac{l(l+1)}{2l+1}(P_{l-1}-P_{l+1})`` and expanding
# ``P'_{m}`` in Legendre polynomials of lower degree turns it into a finite
# combination of integrals of the first kind at shifted indices. **Every angular
# integral in the problem is therefore of one type**, and it is worth a name:
#
# ```math
# G^{k}_{l m} \;\equiv\; \frac{2l+1}{2}\int_{-1}^{1}P_l\,P_{k}\,P_{m}\,d\mu
# ```
#
# a Gaunt coefficient, pure geometry, depending on three integers and on
# nothing about the fluid. Likewise the radial integral depends only on
# ``(l,k,m)`` and on the current viscosity profile; write it
# ``A^{(i)}_{l m}[\eta_{k}]`` and ``B^{(i)}_{l m}[\eta_{k}]`` -- the two
# radial integrals written out above, pairing with ``G`` and ``H`` respectively.
# The index ``i`` distinguishes the terms that end up multiplying ``\dot{\bm\zeta}``
# (``i=2``) from those multiplying ``\bm\zeta`` (``i=1``), which is where the
# superscripts on ``\mathcal D^{(2)}`` and ``\mathcal D^{(1)}`` below come from.
#
# ### The result
#
# With those names the double sum collapses. Both angular factors obey the same
# selection rule, proved below, so the sum over ``k`` terminates. Define
#
# ```math
# \mathcal D^{(i)}_{l m} \;=\; \sum_{k}\Bigl[\,
#   G^{k}_{l m}\,A^{(i)}_{l m}[\eta_{k}]
#   \;+\; H^{k}_{l m}\,B^{(i)}_{l m}[\eta_{k}]\,\Bigr],
# ```
#
# and the modal system takes exactly the Newtonian form with the two diagonal
# matrices replaced by full ones:
#
# ```math
# \boxed{\;
# \bm{\ddot\zeta} \;+\; \mathcal D^{(2)}\,\bm{\dot\zeta}
#              \;+\; \mathcal D^{(1)}\,\bm\zeta \;+\; \bm b \;=\; 0 \;}
# ```
#
# Here ``\bm b`` is the film-pressure forcing, with entries ``\tfrac{4\pi}{2l+1}p_{c,l}``
# so that ``-\bm b`` is the generalised force ``Q_{\zeta_l}``; ``\mathcal D^{(2)}`` is the
# damping matrix, generalising the diagonal ``2\bm\Lambda = \mathrm{diag}(2\lambda_l)`` of the
# Newtonian problem, and ``\mathcal D^{(1)}`` is the stiffness, generalising
# ``\bm\Omega = \mathrm{diag}(\omega_l^2)``. The viscous stress is linear in the velocity and so
# feeds ``\mathcal D^{(2)}`` directly; it reaches ``\mathcal D^{(1)}`` through the
# normal-stress condition, where ``\eta`` appears multiplicatively at the
# surface -- which is also why Reid's ``\omega_l^2`` depends on viscosity at all.
#
# ### Where the space went
#
# ``\eta`` is a field over the drop, yet ``\mathcal D`` carries no ``x`` and no
# ``\theta``. That is not an inconsistency: **both integrals above are
# definite**, so the spatial dependence is integrated out and what survives is
# one number per ``(l,m)`` pair. ``G`` and ``H`` integrate the angle away over
# ``\mu\in[-1,1]``; ``A^{(i)}`` and ``B^{(i)}`` integrate the radius away over
# ``x\in[0,1]``.
#
# What does *not* integrate away is the dependence on the **state**. The
# coefficients ``\eta_{k}(x)`` entering those integrals are the Legendre
# coefficients of ``\eta\bigl(\dot\gamma(x,\theta,t)\bigr)``, and ``\dot\gamma``
# is built from the current modal velocities. So
#
# ```math
# \mathcal D^{(i)} \;=\; \mathcal D^{(i)}\bigl[\bm{\dot\zeta}(t)\bigr] :
# \qquad\text{space integrated out, state dependence retained.}
# ```
#
# Every entry is a functional of the whole velocity field at the current
# instant -- which is what makes the system quasi-linear rather than linear, and
# is taken up in *Where the shear rate is evaluated* on the model page.
#
# !!! note "What is in closed form here, and what is not"
#     The angular factors ``G`` and ``H`` are closed-form integrals of Legendre
#     polynomials, and the selection rule proved below follows from them. The
#     radial factors ``A^{(i)}`` and ``B^{(i)}`` are definite integrals of the
#     interior profiles, and those profiles solve the boxed interior problem of
#     the previous section. So the model is closed -- every unknown is
#     determined by an equation stated on this page -- but it is not *explicit*:
#     no step of the chain has a formula in elementary functions, because
#     ``\eta(\dot\gamma)`` does not. Closed and explicit are different
#     properties, and only the first is needed to have a model.
#
#

let  #src
    Pl(l, m) = l == 0 ? one(m) : l == 1 ? m :  #src
        begin am, b = one(m), m; for n in 1:l-1; b, am = ((2n+1)*m*b - n*am)/(n+1), b; end; b end  #src
    dPl(l, m) = l == 0 ? zero(m) : l*(m*Pl(l,m) - Pl(l-1,m))/(m^2 - 1)  #src
    worst_rec = 0.0  #src
    for l in 1:12, m in (-0.93, -0.41, 0.17, 0.58, 0.86)  #src
        lhs = (1 - m^2)*dPl(l, m)  #src
        rhs = l*(l+1)/(2l+1)*(Pl(l-1,m) - Pl(l+1,m))  #src
        worst_rec = max(worst_rec, abs(lhs - rhs))  #src
    end  #src
    @assert worst_rec < 1e-12 "the recurrence reducing the derivative integrals is wrong"  #src
    gn, gw = DropSolver.gauss_legendre_nodes(60, -1.0, 1.0)  #src
    worst_deriv = 0.0  #src
    for l in 2:8, lp in 0:6, lpp in 2:8  #src
        if l > lp + lpp  #src
            v = sum(w*Pl(l,m)*(1-m^2)*dPl(lp,m)*dPl(lpp,m) for (m,w) in zip(gn,gw))  #src
            worst_deriv = max(worst_deriv, abs(v))  #src
        end  #src
    end  #src
    @assert worst_deriv < 1e-12 "the derivative-type angular integral broke the selection rule"  #src
    println("  ASSERTION 4b OK: (1-mu^2)P_n' = n(n+1)/(2n+1)(P_{n-1}-P_{n+1}) to")  #src
    println("    $(round(worst_rec, sigdigits=2)); and the derivative-type angular integral")  #src
    println("    obeys the same l <= l'+l'' selection rule (max $(round(worst_deriv, sigdigits=2))),")  #src
    println("    so both angular families reduce to Gaunt coefficients.")  #src
end  #src

@variables rr tt qq hh1 hh2  #src

# ### The pressure has to be recovered, and it is no longer harmonic
#
# Taking the curl removed the pressure from the interior problem, which is why
# that problem closed. The traction route now has to reintroduce it, because the
# normal-stress balance contains it explicitly -- a cost of this route and not of
# the model: the energy route obtains the same surface equation by varying the
# surface energy and never needs a pressure at all, as the model page shows.
#
# Taking the divergence of the momentum equation and using ``\nabla\cdot\bm u=0``
# gives a Poisson equation for it:
#
# ```math
# \nabla^2 p \;=\; \mathrm{Oh}\;\nabla\cdot\bigl(\nabla\cdot(2\hat\eta\bm e)\bigr) .
# ```
#
# For a **constant** viscosity the right-hand side vanishes identically, because
# ``\nabla\cdot(\eta\nabla^2\bm u)=\eta\nabla^2(\nabla\cdot\bm u)=0``. The
# pressure is then harmonic, and the regular axisymmetric solution is
# ``p\propto x^lP_l`` -- one amplitude per mode, no differential equation to
# solve. That is the form used throughout the Newtonian theory.
#
# For a variable viscosity the source does **not** vanish. Using the identity
# from the start of this page, ``\nabla\cdot(2\eta\bm e)=\eta\nabla^2\bm u
# +2(\nabla\eta)\cdot\bm e``, its divergence is
#
# ```math
# \nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr)
#   = (\nabla\eta)\cdot(\nabla^2\bm u)
#   + 2\,\nabla\cdot\bigl((\nabla\eta)\cdot\bm e\bigr),
# ```
#
# both terms carrying a derivative of ``\eta``. So ``p`` is the solution of an
# elliptic problem with a source built from the current flow, and the closed-form
# ``x^lP_l`` is not available. This is a real cost of the model and it is easy to
# miss, because the curl hides it: the interior problem looks no harder, and then
# the surface condition asks for a quantity that now requires its own solve.
#
# It is worth being precise about what does and does not survive. **Radial
# variation alone is enough to break it** -- the check below finds an ``O(1)``
# source for ``\eta=\eta(x)``, with no angular structure anywhere. What angular
# structure additionally costs is the separability of the pressure problem, not
# its existence. So this is not a consequence of mode coupling; it is a
# consequence of the viscosity varying at all.
#
# A failing check here would mean the source vanishes for the axisymmetric
# poloidal fields this problem actually produces -- a special cancellation not
# visible in the identity above -- in which case ``p\propto x^lP_l`` would
# survive and the surface equation would be markedly cheaper.

let  #src
    ## The double divergence is evaluated by nested fourth-order central          #src
    ## differences rather than symbolically: div(div(.)) of a variable-viscosity  #src
    ## stress expands to an expression large enough that Symbolics takes tens of  #src
    ## minutes on it, and nothing here needs closed form -- only whether the      #src
    ## quantity is zero. Four nested differences at h = 1e-2 floor the accuracy   #src
    ## near 1e-6, which is ample against an O(1) signal.                          #src
    C2(t) = sin(t)^2*(3cos(t))/6  #src
    Uf(x) = x^3 + 0.7x^5 - 0.3x^7  #src
    psi(x,t) = Uf(x)*C2(t)  #src
    h = 1e-2  #src
    d1(f, v, i) = (-f(v .+ 2h*i) + 8f(v .+ h*i) - 8f(v .- h*i) + f(v .- 2h*i))/(12h)  #src
    er_(x,t) =  d1(v -> psi(v[1],v[2]), [x,t], [0.0,1.0])/(x^2*sin(t))  #src
    et_(x,t) = -d1(v -> psi(v[1],v[2]), [x,t], [1.0,0.0])/(x*sin(t))  #src
    e_rr(x,t) = d1(v -> er_(v[1],v[2]), [x,t], [1.0,0.0])  #src
    e_tt(x,t) = d1(v -> et_(v[1],v[2]), [x,t], [0.0,1.0])/x + er_(x,t)/x  #src
    e_pp(x,t) = er_(x,t)/x + et_(x,t)*cos(t)/(sin(t)*x)  #src
    e_rt(x,t) = (d1(v -> er_(v[1],v[2]), [x,t], [0.0,1.0])/x  #src
               + d1(v -> et_(v[1],v[2]), [x,t], [1.0,0.0]) - et_(x,t)/x)/2  #src
    function divtau(eta, x, t)  #src
        trr(a,b) = 2eta(a,b)*e_rr(a,b); ttt(a,b) = 2eta(a,b)*e_tt(a,b)  #src
        tpp(a,b) = 2eta(a,b)*e_pp(a,b); trt(a,b) = 2eta(a,b)*e_rt(a,b)  #src
        dr = (d1(v -> v[1]^2*trr(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
            + d1(v -> trt(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
            - (ttt(x,t)+tpp(x,t))/x)  #src
        dt = (d1(v -> v[1]^2*trt(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
            + d1(v -> ttt(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
            + trt(x,t)/x - cos(t)/(sin(t)*x)*tpp(x,t))  #src
        (dr, dt)  #src
    end  #src
    function lap_p_source(eta, x, t)  #src
        ar(a,b) = divtau(eta,a,b)[1]; at(a,b) = divtau(eta,a,b)[2]  #src
        (d1(v -> v[1]^2*ar(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
       + d1(v -> at(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t)))  #src
    end  #src
    divu(x,t) = (d1(v -> v[1]^2*er_(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
               + d1(v -> et_(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t)))  #src
    eta_c(x,t) = 1.7  #src
    eta_r(x,t) = 1.3 + 0.9x + 0.4x^2  #src
    eta_a(x,t) = 1.3 + 0.9x*cos(t) + 0.4*(3cos(t)^2 - 1)/2  #src
    worst_div, worst_const = 0.0, 0.0  #src
    least_rad, least_ang = Inf, Inf  #src
    for (x,t) in ((0.61,1.30), (0.43,0.70), (0.82,2.10))  #src
        worst_div   = max(worst_div, abs(divu(x,t)))  #src
        worst_const = max(worst_const, abs(lap_p_source(eta_c,x,t)))  #src
        least_rad   = min(least_rad, abs(lap_p_source(eta_r,x,t)))  #src
        least_ang   = min(least_ang, abs(lap_p_source(eta_a,x,t)))  #src
    end  #src
    @assert worst_div < 1e-10 "the trial field is not divergence free ($worst_div); the probe is invalid"  #src
    @assert worst_const < 1e-4 "constant eta must leave the pressure harmonic, got $worst_const"  #src
    @assert least_rad > 1e-2 "a radially varying eta must break harmonicity, got $least_rad"  #src
    @assert least_ang > 1e-2 "an angularly varying eta must break harmonicity, got $least_ang"  #src
    @assert least_rad > 1e3*worst_const "the variable-eta source is not clearly above the difference floor"  #src
    @printf("  ASSERTION 3d OK: lap(p) source is %.1e for constant eta (the finite-\n", worst_const)  #src
    @printf("    difference floor) but at least %.2e for eta(x) and %.2e for eta(x,theta).\n", least_rad, least_ang)  #src
    println("    So the pressure is harmonic ONLY at constant viscosity. Radial variation")  #src
    println("    alone breaks it -- this is not a mode-coupling effect -- and p must be")  #src
    println("    recovered from an elliptic problem instead of read off as p ~ x^l P_l.")  #src
    println("    Physical meaning of a failure: p ~ x^l P_l would survive, and the surface")  #src
    println("    equation would need no pressure solve at all.")  #src
end  #src

# ### The pressure equation is the divergence of momentum
#
# Taking the divergence of ``\partial_t\bm u=-\nabla p+\mathrm{Oh}\,\nabla\cdot
# (2\eta\bm e)`` and using ``\partial_t(\nabla\cdot\bm u)=0`` gives the Poisson
# equation with the **same** coefficient the momentum equation carries. That
# inheritance is the part worth checking, because a prefactor picked up by hand
# is exactly the kind of error that survives inspection: the identity
#
# ```math
# \nabla\cdot\bigl(-\nabla p+c\,\nabla\cdot(2\eta\bm e)\bigr)
#   \;+\;\nabla^2p\;-\;c\,\nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr)\;=\;0
# ```
#
# holds for every ``c``, so the ``c`` appearing in the pressure equation is
# forced to be the ``c`` in the momentum equation, and no other value works.
# Projecting onto ``P_l`` then turns it into the radial form the summary states,
# by the angular reduction already established plus Legendre orthogonality.

let  #src
    ## Differentiated by nested central differences rather than symbolically: the  #src
    ## double divergence of a variable-viscosity stress is the expression that     #src
    ## made Symbolics run for tens of minutes earlier in this file's history, and  #src
    ## only the vanishing of a residual is in question here. Errors are            #src
    ## normalised, since four nested differences floor the accuracy near 1e-6.     #src
    h = 1e-2  #src
    d1(f, v, i) = (-f(v .+ 2h*i) + 8f(v .+ h*i) - 8f(v .- h*i) + f(v .- 2h*i))/(12h)  #src
    LPn(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    Cn(l,t) = sin(t)^2*(l*(cos(t)*LPn(l,cos(t)) - LPn(l-1,cos(t)))/(cos(t)^2 - 1))/(l*(l+1))  #src
    ## a divergence-free field from a single-mode stream function  #src
    psif(m) = (x,t) -> (x^(m+1) + 0.6x^(m+3))*Cn(m,t)  #src
    ## the two components of div(2 eta e), and the scalar div of a vector field  #src
    function divtau2(eta, ps, x, t)  #src
        ur(a,b) =  d1(v -> ps(v[1],v[2]), [a,b], [0.0,1.0])/(a^2*sin(b))  #src
        ut(a,b) = -d1(v -> ps(v[1],v[2]), [a,b], [1.0,0.0])/(a*sin(b))  #src
        e_rr(a,b) = d1(v -> ur(v[1],v[2]), [a,b], [1.0,0.0])  #src
        e_tt(a,b) = d1(v -> ut(v[1],v[2]), [a,b], [0.0,1.0])/a + ur(a,b)/a  #src
        e_pp(a,b) = ur(a,b)/a + ut(a,b)*cos(b)/(sin(b)*a)  #src
        e_rt(a,b) = (d1(v -> ur(v[1],v[2]), [a,b], [0.0,1.0])/a  #src
                   + d1(v -> ut(v[1],v[2]), [a,b], [1.0,0.0]) - ut(a,b)/a)/2  #src
        dr = (d1(v -> v[1]^2*2eta(v[1],v[2])*e_rr(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
            + d1(v -> 2eta(v[1],v[2])*e_rt(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
            - (2eta(x,t)*e_tt(x,t) + 2eta(x,t)*e_pp(x,t))/x)  #src
        dt = (d1(v -> v[1]^2*2eta(v[1],v[2])*e_rt(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
            + d1(v -> 2eta(v[1],v[2])*e_tt(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
            + 2eta(x,t)*e_rt(x,t)/x - cos(t)/(sin(t)*x)*2eta(x,t)*e_pp(x,t))  #src
        (dr, dt)  #src
    end  #src
    divvec(ar, at, x, t) = (d1(v -> v[1]^2*ar(v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
                          + d1(v -> at(v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t)))  #src
    lapscal(f, x, t) = (d1(v -> v[1]^2*d1(w -> f(w[1],w[2]), [v[1],v[2]], [1.0,0.0]), [x,t], [1.0,0.0])/x^2  #src
                      + d1(v -> sin(v[2])*d1(w -> f(w[1],w[2]), [v[1],v[2]], [0.0,1.0]), [x,t], [0.0,1.0])/(x^2*sin(t)))  #src

    PTS = ((0.61,1.30), (0.43,0.70), (0.82,2.10))  #src
    etaf = (x,t) -> 1.3 + 0.9x*cos(t) + 0.4*(3cos(t)^2 - 1)/2  #src
    ## SUM-PRESS: the identity, for two different coefficients c  #src
    worst_id, mag_id, wrong_c = 0.0, 0.0, Inf  #src
    for m in (2, 3), cvis in (0.31, 1.7)  #src
        ps = psif(m)  #src
        pf = (x,t) -> (x^(m+2) + 0.4x^(m+4))*LPn(m, cos(t))  #src
        Ar(x,t) = -d1(v -> pf(v[1],v[2]), [x,t], [1.0,0.0]) + cvis*divtau2(etaf, ps, x, t)[1]  #src
        At(x,t) = -d1(v -> pf(v[1],v[2]), [x,t], [0.0,1.0])/x + cvis*divtau2(etaf, ps, x, t)[2]  #src
        for (x,t) in PTS  #src
            S  = divvec((a,b) -> divtau2(etaf, ps, a, b)[1],  #src
                        (a,b) -> divtau2(etaf, ps, a, b)[2], x, t)  #src
            lp = lapscal(pf, x, t)  #src
            resid = divvec(Ar, At, x, t) + lp - cvis*S  #src
            worst_id = max(worst_id, abs(resid))  #src
            mag_id   = max(mag_id, abs(lp), abs(cvis*S))  #src
            ## the SAME residual with a doubled coefficient must NOT vanish  #src
            wrong_c = min(wrong_c, abs(divvec(Ar, At, x, t) + lp - 2cvis*S)/max(mag_id, 1e-30))  #src
        end  #src
    end  #src
    rel_id = worst_id/mag_id  #src
    @assert mag_id > 1e-3 "the pressure-identity sweep never exercised a nonzero operator ($mag_id)"  #src
    @assert rel_id < 1e-4 "div of momentum does not give lap(p) = c div(div(2 eta e)) (rel $rel_id)"  #src
    @assert wrong_c > 1e-2 "a doubled coefficient also satisfies the identity, so the prefactor is not pinned ($wrong_c)"  #src

    ## SUM-PLAP: projecting lap(p) = Oh S onto P_l gives L_l[p_l] = Oh S_l  #src
    nodes, wts = QuadGK.gauss(28, -1.0, 1.0)  #src
    proj(f, l, x) = (2l+1)/2*sum(w*f(x, acos(mu))*LPn(l,mu) for (mu,w) in zip(nodes,wts))  #src
    worst_pl, mag_pl = 0.0, 0.0  #src
    for l in (2, 3, 4)  #src
        pl(x) = x^(l+2) + 0.4x^(l+4)  #src
        pf = (x,t) -> pl(x)*LPn(l, cos(t))  #src
        for (x,_) in PTS  #src
            ## L_l[p_l] from its definition, by differences in x alone  #src
            dpl  = (-pl(x+2h) + 8pl(x+h) - 8pl(x-h) + pl(x-2h))/(12h)  #src
            d2pl = (-pl(x+2h) + 16pl(x+h) - 30pl(x) + 16pl(x-h) - pl(x-2h))/(12h^2)  #src
            Lp = d2pl + 2dpl/x - l*(l+1)*pl(x)/x^2  #src
            lhs = proj((a,b) -> lapscal(pf, a, b), l, x)  #src
            worst_pl = max(worst_pl, abs(lhs - Lp))  #src
            mag_pl   = max(mag_pl, abs(Lp))  #src
        end  #src
    end  #src
    rel_pl = worst_pl/mag_pl  #src
    @assert mag_pl > 1e-3 "the projection sweep never exercised a nonzero operator ($mag_pl)"  #src
    @assert rel_pl < 1e-4 "the l-projection of lap(p) is not L_l[p_l] (rel $rel_pl)"  #src

    ## SUM-COM: the non-dimensional centre-of-mass equation, signs included.      #src
    ## Dimensional: m dV/dT = -m g + F_z. With V = (R/Tsig) v, T = Tsig t, and     #src
    ## Tsig^2 = rho R^3/T1, dividing by m and multiplying by Tsig^2/R gives         #src
    ##   vdot = -(g Tsig^2/R) + F_z Tsig^2/(m R) = -Bo + Fhat/(4pi/3),              #src
    ## where Fhat is F_z in units of T1 R and m = (4pi/3) rho R^3.                  #src
    worst_com = 0.0  #src
    for (rho_,R_,T1_,g_) in ((998.0,3.5e-4,0.0722,9.81), (1210.0,7.1e-4,0.0640,9.81), (1.0,2.0,3.0,5.0))  #src
        Tsig2 = rho_*R_^3/T1_  #src
        Bo    = rho_*g_*R_^2/T1_  #src
        m     = (4pi/3)*rho_*R_^3  #src
        worst_com = max(worst_com, abs(g_*Tsig2/R_ - Bo))               # gravity term is exactly Bo  #src
        for Fhat in (-0.7, 0.0, 2.3)                                     # F_z = Fhat * T1 * R  #src
            got  = (Fhat*T1_*R_)*Tsig2/(m*R_)  #src
            want = Fhat/(4pi/3)  #src
            worst_com = max(worst_com, abs(got - want))  #src
        end  #src
    end  #src
    ## and Fhat = -(4pi/3) p_{c,1} (SUM-FORCE), so vdot = -Bo - p_{c,1}  #src
    for pc1 in (-1.3, 0.0, 0.9)  #src
        worst_com = max(worst_com, abs((-(4pi/3)*pc1)/(4pi/3) - (-pc1)))  #src
    end  #src
    @assert worst_com < 1e-12 "the centre-of-mass equation does not reduce to vdot = -Bo - p_{c,1} ($worst_com)"  #src

    println("  ASSERTION 3f OK: the pressure equation and the centre of mass --")  #src
    @printf("    div of momentum gives lap(p) = c div(div(2 eta e)) with the SAME c\n")  #src
    @printf("      (rel %.1e); doubling c breaks it by %.1e, so the prefactor is pinned;\n", rel_id, wrong_c)  #src
    @printf("    the l-projection of lap(p) is L_l[p_l] (rel %.1e), which is the radial\n", rel_pl)  #src
    @printf("      form the summary states;\n")  #src
    @printf("    and vdot = -Bo - p_{c,1} follows from m dV/dT = -m g + F_z with the\n")  #src
    @printf("      stated scalings, signs included (%.1e).\n", worst_com)  #src
    println("    Physical meaning of a failure: the pressure would carry a wrong")  #src
    println("    prefactor, or the drop would fall upwards.")  #src
end  #src
#

# ### Do the two routes agree?
#
# The energy route gives the damping as a volume integral of
# ``2\eta\,\bm e^{(l)}\!:\!\bm e^{(m)}``; the traction route gives it by projecting
# the surface stress. They are equal by integration by parts -- the interior terms
# cancel against the momentum equation and the boundary term is what the surface
# projection collects -- **provided** the field used is a solution of the interior
# problem and the surface terms that parts throws off actually vanish. Neither is
# automatic, so it is checked.
#
# The cleanest place to check it is where an independent answer exists: constant
# viscosity, where the traction route must return Reid's ``\lambda_l`` and the
# energy route must return the same number from a volume integral. In the
# small-viscosity limit both must equal Lamb's ``(l-1)(2l+1)\mathrm{Oh}``.

let  #src
    ## Energy route, evaluated on the potential-flow field, against Reid's own      #src
    ## lambda_l from the running solver as Oh -> 0. This is weaker than a           #src
    ## coefficient-by-coefficient comparison at finite Oh -- which needs the        #src
    ## variable-eta interior solution and so is not available yet -- but it is the   #src
    ## strongest check that can be made without one, and it is a comparison of two  #src
    ## genuinely independent computations rather than of an expression with itself.  #src
    Drv = Differential(rr); Dtv = Differential(tt)  #src
    edv(e) = Symbolics.expand_derivatives(e)  #src
    LPv(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    dLPv(l,m) = l==0 ? zero(m) : l*(m*LPv(l,m)-LPv(l-1,m))/(m^2-1)  #src
    Cgv(l) = sin(tt)^2*dLPv(l,cos(tt))/(l*(l+1))  #src
    nrv, wrv = QuadGK.gauss(40, 0.0, 1.0); nmv, wmv = QuadGK.gauss(40, -1.0, 1.0)  #src
    function forms(l)  #src
        psi = rr^(l+1)*Cgv(l)                       # potential flow  #src
        ur =  edv(Dtv(psi)/(rr^2*sin(tt))); ut = -edv(Drv(psi)/(rr*sin(tt)))  #src
        e_rr = edv(Drv(ur)); e_tt = edv(Dtv(ut)/rr + ur/rr)  #src
        e_pp = ur/rr + ut*cos(tt)/(sin(tt)*rr)  #src
        e_rt = edv((Dtv(ur)/rr + Drv(ut) - ut/rr)/2)  #src
        vint(e) = begin  #src
            g = Symbolics.build_function(edv(e), rr, tt; expression=Val(false))  #src
            s = 0.0  #src
            for (x,wx) in zip(nrv,wrv), (mu,wu) in zip(nmv,wmv)  #src
                s += wx*wu*2pi*x^2*g(x, acos(mu))  #src
            end  #src
            s  #src
        end  #src
        (vint(2*(e_rr^2 + e_tt^2 + e_pp^2 + 2e_rt^2)), vint(ur^2 + ut^2))  #src
    end  #src
    worst = 0.0  #src
    for l in 2:5  #src
        Phi, T = forms(l)  #src
        energy_lambda = Phi/(2T)                                  # per unit Oh  #src
        ## Reid's lambda_l/Oh as Oh -> 0, from the solver's own root finder  #src
        traction_lambda = DropSolver.reid_lambda_omega2(1e-6, l)[1]/1e-6  #src
        worst = max(worst, abs(energy_lambda - traction_lambda)/traction_lambda)  #src
    end  #src
    @assert worst < 1e-3 "the energy and traction routes disagree on lambda_l ($worst)"  #src
    @printf("  ASSERTION 5c OK: energy route and traction route agree on lambda_l to\n")  #src
    @printf("    %.1e relative over l = 2..5, comparing a volume integral of the\n", worst)  #src
    println("    dissipation against the solver's own root of Reid's characteristic")  #src
    println("    equation. Two independent assemblies, not one expression twice.")  #src
    println("    Physical meaning of a failure: the surface projection and the volume")  #src
    println("    dissipation would be computing different things, and only one of them")  #src
    println("    could be the damping.")  #src
end  #src

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
# | (b) period-averaged | the ``j=0`` channel | justified by the lemma below |
# | (c) Floquet | ``j=0`` and ``j=2`` | most faithful, most work |
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
#   \;+\; \sum_{j\ge1}\bigl[\alpha_j\cos 2j\phi + \beta_j\sin 2j\phi\bigr].
# ```
#
# Put ``\zeta_l=a(t)\cos\Theta``, ``\Theta=\omega t+\chi``, with the envelope
# ``a`` and the phase ``\chi`` slow, into
# ``\ddot\zeta_l+2\lambda(t)\dot\zeta_l+\omega^2\zeta_l=0``. Krylov--Bogoliubov
# averaging gives ``\dot a=-2a\,\langle\lambda\sin^2\Theta\rangle``, and with
# ``\sin^2\Theta=\tfrac12[1-\cos(2\phi+2\chi)]`` every harmonic ``j\ge2`` is
# orthogonal to ``\cos(2\phi+2\chi)`` and drops out. What survives is
#
# ```math
# \boxed{\;
# \frac{da}{dt} \;=\; -\,\bar\lambda\,a
#   \;+\; \frac{a}{2}\bigl(\alpha_1\cos 2\chi - \beta_1\sin 2\chi\bigr) \;}
# ```
#
# Three things follow.
#
# **The period-average is the whole secular effect.** ``\bar\lambda`` -- the
# ``j=0`` harmonic -- is the only term that survives without a phase condition.
# That is what justifies closure (b), and it justifies it *quantitatively*
# rather than by appeal to plausibility.
#
# **The ``j=2`` harmonic is a parametric term, and it sits exactly on
# resonance.** It enters multiplied by ``\cos2\chi``, so its sign depends on the
# phase -- and the period-``\pi`` lemma places it at ``2\omega``, precisely the
# principal parametric resonance of an oscillator at ``\omega``. The lemma is
# therefore not a reassurance about the second harmonic: it identifies a
# resonance rather than ruling one out.
#
# **A bound settles when that matters.** Since ``|\alpha_1\cos2\chi-\beta_1\sin2\chi|
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
        ## the j=0 harmonic contributes 1/2, independent of phase  #src
        worst = max(worst, abs(avg(ph -> sin(ph + ps)^2) - 0.5))  #src
        ## the j=2 harmonics contribute -cos(2psi)/4 and +sin(2psi)/4  #src
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
    println("    => the period-average carries the secular decay; the j=2 harmonic is")  #src
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
#   frequency -- so the leading effect is carried entirely by the ``j=0``
#   (period-averaged) channel. That is choice (b), and the lemma is its
#   justification.
# * The surviving ``j=2`` channel modulates the damping at ``2\omega``,
#   which is the *principal parametric resonance* condition. Choice (b)
#   discards it; whether that is safe is a stability question this file
#   does not answer.
#
# Fourier-decomposing ``S`` over a full ``2\pi`` period, for a strain field
# with four independent complex components, gives
#
# | harmonic | ``j=0`` | ``j=1`` | ``j=2`` |
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
# stability of discarding the ``j=2`` channel is open -- which is how the
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
    @assert abs(c2) > 1e-3 "the j=2 parametric channel vanished; expected it to survive"  #src
    println("  ASSERTION 8 OK: S is period-pi to $(round(worst, sigdigits=2)); its Fourier")  #src
    println("    content at j=1 is zero to machine precision, while j=0 and j=2 survive.")  #src
    println("    (A lemma check on an arbitrary field -- not a physical state.)")  #src
end  #src

# ## Truncating the viscosity spectrum
#
# **Assumption.** ``\eta``'s Legendre content above ``k=L_\eta`` is negligible.
#
# The model page leaves ``L_\eta`` as whatever the state produces. Capping it is
# the natural first economy, and it makes ``L_\eta`` a free parameter of the
# scheme in the same way ``M`` is: both are points at which a convergent series
# is cut, and both should be raised until the answer stops moving rather than
# fixed by argument. The question is only where that plateau lies.
#
# **What it drops.** Coupling between modes further apart than ``L_\eta``. The
# matrix becomes banded and the cost of applying it falls from ``O(M^2)`` to
# ``O(M L_\eta)``.
#
# **Error.** Exactly the discarded coupling, which is a measurable number. It
# has been measured, and the measurement puts the plateau uncomfortably high.
#
# ### Measurement, and the choice of error norm
#
# The natural-looking metric, "what fraction of ``\sum_{k}|\eta_{k}|^2``
# lies below ``L_\eta``", says ``L_\eta=2\ldots5`` suffices everywhere. That
# metric is **misleading**. What controls the banding error is the *summed
# magnitude* of the discarded coefficients,
# ``\mathcal T(L)=\sum_{k>L}|\eta_{k}|/|\eta_0|``, because those terms enter the
# coupling matrix additively, not in quadrature. The spectrum turns out to
# have a long, slowly-decaying plateau carrying little *power* but large
# *summed magnitude*.
#
# Measured on Reid's actual viscous eigenmodes (volume-averaged,
# ``\mathrm{Oh}=0.2``, in the saturated thinned regime), the smallest
# ``L_\eta`` holding ``\mathcal T(L_\eta)<10^{-2}``:
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
# ``\dot\zeta_l\propto0.85^{\,l}`` -- keeps ``L_\eta\approx16`` even with fifty
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
# **Caveat.** Coefficients beyond ``k=140`` were not measured, so every
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
# ``G^{0}_{l m}=\delta_{l m}``, the matrices return to diagonal, and the
# entire architecture of the solver -- one independent oscillator per mode --
# is recovered intact.
#
# What you give up is *only* the closed form. The radial equation now has an
# ``r``-dependent coefficient, so it is no longer Bessel's equation and must
# be solved as a numerical two-point boundary-value problem, once per mode.
# `src/reid.jl`'s continuation machinery already knows how to track
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
#   ``\tau_{r\theta}=0\Rightarrow\mathcal T[\psi]=0`` reduction already derived
#   in *The Free Viscous Drop*.
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
# operator ``\mathcal D_l(\mathcal D_l+q^2)\psi=0`` identically, which is how the
# construction is checked.
#
# What remains is a linear two-point boundary-value eigenproblem with variable
# coefficients: no longer Bessel, but entirely standard. The continuation
# machinery that tracks Reid's eigenvalue branches applies to it unchanged.
#
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
# ``\mathrm{Oh}_{\mathrm{eff}}=\eta_{\mathrm{eff}}/\sqrt{\rho \gamma R}``.
#
# **What it buys.** The radial equation becomes constant-coefficient again,
# i.e. Bessel's equation, and Reid's closed-form characteristic equation
# returns verbatim -- evaluated at a shifted Ohnesorge number.
#
# **This is where the nonvariational solver sits**, combined with choice (a) of
# the temporal closure: `solve_drop!` with an `STExactParams`. It has no
# interior state, so it cannot evaluate a viscosity that varies through the
# drop, and the whole descent above exists to reduce that variation to one
# number per mode.
#
# The variational solvers do not descend this far. `simulate` and
# `simulate_lcp` carry the full coupled system at the top of the table, with
# `eta` evaluated pointwise on the summed strain field.
#
# The gap between the two is not academic. At the operating Ohnesorge number of
# the most concentrated fluid this package is checked against, the
# effective-viscosity closure reports a contact time far outside the measured
# range -- long enough to say the drop never releases -- while the coupled
# system does not. The closure is outside its range there, and the rest of
# this page is the account of why.
#
# Its error is the one measured in the previous section: for a realistic
# multi-mode state the
# viscosity spans a factor of order ``100`` across the drop, and the angular
# coefficients ``|\eta_1|/|\eta_0|\approx1.3`` and
# ``|\eta_2|/|\eta_0|\approx0.2`` are comparable to or larger than the mean.
# Collapsing all of that onto one number per mode is a leading-order
# approximation, not a perturbative one. Combined with the eigenmode closure it
# is stated within, this rung carries two
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
# the whole series. The Gaunt rule requires ``k\le l+m``, and the shape
# expansion stops at ``l=M``, so **every viscosity harmonic above
# ``k=2M`` is orthogonal to every product ``P_lP_{m}`` and cannot couple
# any pair of modes.** It is not small; it is absent.
#
# > The mode-coupling matrix is determined **exactly** by the first ``2M+1``
# > Legendre harmonics of the viscosity. Truncating there is not an
# > approximation -- everything above is annihilated by the angular integral.
#
# So the infinite series is a property of ``\eta``, not a limitation of the
# model. For ``M=50`` the coupling is fixed by 101 numbers per radius; for
# ``M=90``, by 181. Whether the matrix can be further *banded* -- kept to
# ``k\le L_\eta`` with ``L_\eta\ll 2M``, for speed -- is a separate and
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
# fractional power of the *field* when ``a`` is not an even integer, so the
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
    println("    limit that src/reid.jl already implements and validates.")  #src
    println("    Physical meaning of a failure: the hierarchy would not contain the")  #src
    println("    Newtonian case it claims to generalise.")  #src
end  #src

## Switching the thinning off must return the Newtonian operator, and the check has  #src
## to run through the SHEAR-THINNING code path to mean anything: the whole point is   #src
## that the expensive path, which rebuilds the coupled matrices from a pointwise      #src
## viscosity, agrees with the cheap one that assumes a constant viscosity and         #src
## assembles block by block. Comparing either against itself would prove nothing.     #src
let  #src
    Oh0 = 0.05  #src
    bas = DropSolver.ModalBasis(2:6, 2)  #src
    st = [0.6/i * (k == 1 ? 1.0 : 0.3) for i in 1:length(bas.ls) for k in 1:bas.K]  #src
    ## lambda_c -> 0 sends eta -> 1 pointwise, whatever the shear rate is             #src
    etaf = (x, mu) -> DropSolver.carreau(DropSolver.shear_rate(bas, st, x, mu);  #src
                                        lambda_c=1e-10, a=2.0, n=0.4)  #src
    Fth = DropSolver.assemble_coupled(bas, Oh0; eta = etaf)  #src
    Fnw = DropSolver.assemble_newtonian(bas, Oh0)  #src
    sc = maximum(abs, Fnw.C)  #src
    worst = maximum(abs, Fth.C .- Fnw.C) / sc  #src
    @assert sc > 1.0 "the Newtonian dissipation operator is degenerate ($sc)"  #src
    @assert worst < 1e-8 "the shear-thinning path does not reduce to the Newtonian operator at lambda_c -> 0 (rel $worst)"  #src
    ## and the modes must decouple again, since a constant eta carries only l' = 0    #src
    nl = length(bas.ls)  #src
    off = maximum(DropSolver.block_norm(bas, Fth.C, i, j) for i in 1:nl, j in 1:nl if i != j)  #src
    dia = maximum(DropSolver.block_norm(bas, Fth.C, i, i) for i in 1:nl)  #src
    @assert off < 1e-8 * dia "the thinning path leaves the modes coupled at lambda_c -> 0 ($(off/dia))"  #src
    @printf("  ASSERTION 12b OK: with the thinning switched off the coupled assembly\n")  #src
    @printf("    reproduces the constant-viscosity operator to %.1e relative, and the\n", worst)  #src
    @printf("    off-diagonal blocks vanish to %.1e of the diagonal. So the hierarchy\n", off/dia)  #src
    println("    contains the Newtonian case it claims to generalise, and contains it")  #src
    println("    through the same code the shear-thinning runs use.")  #src
    println("    Physical meaning of a failure: the mode coupling would be an artefact")  #src
    println("    of the assembly rather than a consequence of the viscosity field.")  #src
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
# ``\mu\in[-1,1]`` (`assemble_coupled`). At each node ``\dot\gamma`` is computed
# from the full superposition of active modes, ``\eta`` follows pointwise from
# the constitutive law, and a Legendre projection in ``\theta`` at each radius
# produces the ``\eta_{k}(x)`` the matrices need.
#
# The **angular** rule has to scale with the truncation, and this is easy to
# get wrong. The integrand is a product of two strain bases, hence a polynomial
# of degree ``\le2M`` in ``\mu``, so an ``n_x``-node Gauss rule is exact only
# for ``n_x\ge M+\tfrac12``. A rule fixed at ``n_x=30`` is therefore exact at
# ``M=16`` and badly inexact at ``M=90``: the measured error in
# ``\mathrm{Oh}_{\mathrm{eff}}`` grows from ``1.6\times10^{-3}`` to
# ``8.4\times10^{-2}``. With the default ``n_x=\max(30,\,2M+2)`` the accuracy
# is independent of ``M``. `test/test_oh_eff_quadrature.jl` asserts that
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
# \bm{\dot\zeta}^{\,\ast} = (1+r)\,\bm{\dot\zeta}_n - r\,\bm{\dot\zeta}_{n-1},
# \qquad r=\frac{\Delta t_{n+1}}{\Delta t_n},
# ```
#
# exact whenever ``\bm{\dot\zeta}`` is linear in ``t``, for non-uniform steps as
# well as uniform. Then ``\mathcal D(\bm{\dot\zeta}^{\,\ast})`` is constant within
# the step, the system is linear in the unknowns, and the Jacobian is exact.
# This is an IMEX splitting: implicit in the modal amplitudes, explicit in the
# coefficients.
#
# The extrapolation order is not free. Holding the coefficients at the previous
# step instead -- ``\bm{\dot\zeta}^{\,\ast}=\bm{\dot\zeta}_n`` -- is a first-order
# splitting error inside an otherwise second-order BDF2 scheme, and it drags the
# observed order of the whole integration down to ``\approx1.3`` against a
# Newtonian control at ``2.0``. The linear extrapolation above restores it.
# `test/test_convergence_order.jl` measures the order and holds it above
# ``1.5``.
#
# ### The step-size ceiling
#
# There are two distinct constraints here and they are easy to merge by
# accident, because in the diagonalisable Newtonian model they coincide.
#
# **Resolution.** The capillary dispersion relation
# ``\omega_l^2=l(l-1)(l+2)`` is a statement about surface energy against
# inertia; it does not depend on the viscous machinery, and it survives every
# closure on this page. Representing the oscillation of the highest retained
# mode at ``N`` samples per period therefore requires
#
# ```math
# \Delta t \;\lesssim\; \frac{2\pi}{N\sqrt{M(M-1)(M+2)}} \;\sim\; M^{-3/2} ,
# ```
#
# so raising ``M`` costs twice over: more unknowns per step, and more steps.
# That is why the truncation cannot be raised casually to chase the ``L_\eta``
# convergence discussed above.
#
# **Stability.** This is a different question, and it is *not* answered by the
# expression above. Within the instantaneous-eigenmode closure the modal system
# is diagonal, its eigenvalues are exactly ``\lambda_l\pm i\omega_l``, and the
# stability limit of an explicit treatment coincides with the resolution limit.
# The general model of the companion page has no such eigenvalues: the interior
# is a differential-algebraic system, the coefficient matrices are dense, and the
# stiffness is set by viscous diffusion across the drop rather than by any
# surface frequency. The step limit there has to be derived for the scheme
# actually used, and quoting ``M^{-3/2}`` for it would be carrying over a result
# whose derivation no longer applies.
#
# Only the resolution constraint is claimed here.
#
# ### The contact set
#
# The reaction pressure and the gap satisfy a complementarity condition, not an
# equation, so the contact set is an unknown of a different kind: it is
# discrete. Represented on a fixed angular grid it is an integer -- the number of
# nodes in contact -- and the step must decide that integer as well as the
# amplitudes.
#
# **The complementarity conditions determine it, with no objective.** The pair
#
# ```math
# h\ge0,\qquad p_c\ge0,\qquad h\,p_c=0
# ```
#
# says which way to move whenever the current set is wrong, and the two
# inequalities fail on opposite sides. If any node outside the contact has
# ``h<0`` the surface is interpenetrating and the set is too *small*; if the
# pressure at the outermost contacting node has ``p_c<0`` the film is pulling the
# drop down and the set is too *large*. So the set is found by iterating within
# the step -- grow while the first holds, retreat while the second does, stop when
# neither does. Every move is forced by a violated condition, so the iteration
# terminates, and no tie ever has to be broken.
#
# That last point is what makes the difference. Choosing instead among
# neighbouring candidate sets by minimising a residual requires an objective, and
# an objective admits ties: two sets score almost equally, the choice flips from
# step to step, and the reaction pressure oscillates in sign with order-one
# amplitude until the step size collapses. The failure is worst where the physical
# damping is weakest, because there the high-order modes ring freely, and it is
# *sporadic* in the truncation ``M`` rather than monotone in it -- which is the
# signature of a numerical bistability and not of an unresolved scale. Refining
# ``\Delta t`` makes it worse, not better.
#
# An objective, if one is used at all, must at least be local to the contact edge:
# one integrated over the contact region is biased toward vanishing contact,
# because shrinking the region also shrinks the domain the residual is measured
# over, so its minimiser is the smallest admissible set independently of the
# physics.
#
# **The step size remains the safety valve.** If the iteration cannot find an
# admissible set at all, the response is to reject the step and halve
# ``\Delta t`` rather than accept an inadmissible one. This is what keeps the
# nonsmoothness of the contact transition from being straddled by a single step.
#
# ### Two measurement traps
#
# **Contact time.** Contact is not necessarily a single interval: a marginally
# resolved run can chatter, re-entering contact briefly after the physical
# rebound. The span from first to last contact then overstates the contact time,
# while the duration of the *first* interval understates it, badly, on exactly
# those runs. The defensible metric is the duration of the **longest contiguous**
# contact interval.
#
# **Conditioning of the pressure basis.** The pressure coefficients are a
# spectral representation of a field supported on part of the domain, and the
# conditioning of the resulting system depends strongly on how many are retained.
# It is a quantity to be measured at the truncation actually used, not inferred
# from a value quoted at another one.

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
# | poloidal + modal | change of variables | state is ``\{\zeta_l\}``, ``l=2\ldots M`` |
# | **full coupled system** | none beyond the three above | **what the VARIATIONAL solvers run** |
# | instantaneous eigenmode | ``\partial_t\psi_l\to-\sigma \psi_l``; ``\eta`` frozen | eliminates ``\psi_l``; parabolic ``\to`` algebraic |
# | two-root truncation | keep two of infinitely many ``\sigma`` | recovers ``\lambda_l``, ``\omega_l^2`` as two numbers per mode |
# | temporal closure | instantaneous / period-averaged / Floquet | picks which ``\eta(t)`` channel survives |
# | truncate at ``L_\eta`` | discarded coupling small -- **measured false** | banded, but needs ``L_\eta\gtrsim M``: no saving |
# | ``\eta=\eta(x)`` | viscosity spherically symmetric -- **leading-order error** | **diagonal**; numerical radial BVP per mode |
# | ``\eta\to\eta_{\rm eff}`` | radial variation small | closed form -- **what the NONVARIATIONAL solver runs** |
# | ``\eta=\eta_0`` | no thinning | Reid |
# | ``\mathrm{Oh}\to0`` | inviscid | Lamb |
#
# A specific constitutive law is a slice through this table, not a rung of it,
# and can be taken at any row. The numerical realisation is not in the table at
# all: it is not a rung, because it does not change the model.
#
# ## Which concessions are affordable
#
# Reading the descent backwards, some of its rungs cost nothing and others cost
# more than they save.
#
# **Free, and permanent.** A generalized Newtonian fluid cannot break
# axisymmetry, so Legendre polynomials suffice forever. BC1 carries no
# viscosity at all, and BC2 reduces to ``\mathcal T[\psi]=0`` for every fluid in
# the admissible class, because ``\eta\ge\eta_\infty>0`` cannot vanish. The
# normal-stress condition is where the rheology lands: once ``\eta`` varies
# with ``\theta`` its surface value is a field, and projecting it is what
# produces the coupling matrices. Under the spherically symmetric rung above
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
# **Affordability of the coupled system.** It is affordable at ``M\sim45``. The
# pairwise strain contractions are geometry: only the scalar ``\eta`` at each
# quadrature point changes between Picard sweeps, so caching them turns each
# reassembly into a weighted sum. That takes a shear-thinning impact from 145
# seconds to 1.4, and the same run at ``M=45, K=3`` costs about 250 seconds. No
# Gaunt precomputation is involved; the assembly quadratures the dof pairs
# directly.
#
# **Still open.** The eigenvalue problem for the
# variable-``\eta`` radial operator, which is a two-point boundary-value
# problem rather than a research question; and the ``j=2`` parametric channel,
# which the period-``\pi`` lemma places at exactly the principal resonance
# condition and which nothing here evaluates.

# ## Open direction: the interior as a boundary impedance #src
# #src
# This section is analysis rather than a closure: it says what object the interior #src
# problem is *for*, and that identification reorganises the choices above. Nothing #src
# here is yet checked, and it is marked as structure for that reason. #src
# #src
# ### The interior enters only as a boundary impedance #src
# #src
# The surface equation needs one thing from the interior: the viscous normal #src
# traction at ``x=1``, given the surface motion. It never needs the field itself. #src
# Define, mode by mode, the map from surface velocity to surface traction, #src
# #src
# ```math #src
# \bigl[2\eta\,e_{rr}\bigr]_{x=1} \;=\; -\,\mathcal Z_l\bigl[\dot\zeta_l\bigr] , #src
# ``` #src
# #src
# which is a **Dirichlet-to-Neumann** map for the interior problem -- the same #src
# object that appears as a mobility in Stokes flow, an impedance in acoustics, and #src
# an absorbing boundary condition in wave propagation. The whole content of the #src
# interior is compressed into ``\mathcal Z_l``. #src
# #src
# For a normal mode ``\propto e^{-\sigma t}`` and constant viscosity, #src
# ``\mathcal Z_l`` is a known function of ``\sigma``: solving #src
# ``\mathcal D_l(\mathcal D_l+q^2)[\psi_l]=0`` with regularity and BC2 leaves one #src
# free constant, and the ratio of traction to velocity is a rational function of #src
# ``\sigma`` and the Bessel ratio ``Q=j_{l+1}(q)/j_l(q)``. Reid's characteristic #src
# equation is then just the statement that the total impedance vanishes, #src
# #src
# ```math #src
# \underbrace{\sigma^2}_{\text{inertia}} #src
# \;+\;\underbrace{\mathcal Z_l(\sigma)}_{\text{viscous}} #src
# \;+\;\underbrace{\omega_{l,0}^2}_{\text{capillary}} \;=\;0 , #src
# ``` #src
# #src
# with ``\omega_{l,0}^2=l(l-1)(l+2)`` the inviscid frequency. Read this way, the #src
# infinitely many roots are not a peculiarity: they are the poles of an impedance, #src
# and there are always infinitely many because the interior is a diffusion #src
# problem. #src
# #src
# ### The two-root truncation is a two-pole fit, and that is the useful way to see it #src
# #src
# ``\lambda_l`` and ``\omega_l^2`` describe an oscillator with two poles. Since #src
# ``\mathcal Z_l`` has infinitely many, the Newtonian modal model is a **two-pole #src
# rational approximation** of the exact impedance. That reframes the truncation #src
# from a defect into the first member of a convergent family: keep ``N_p`` poles, #src
# #src
# ```math #src
# \mathcal Z_l(\sigma)\;\approx\;\sum_{p=1}^{N_p}\frac{c_p\,\sigma}{\sigma+s_p} , #src
# ``` #src
# #src
# and each pole is **one auxiliary scalar ODE per mode**, #src
# #src
# ```math #src
# \dot g_p = -s_p\,g_p + \dot\zeta_l , #src
# \qquad #src
# \mathcal Z_l\bigl[\dot\zeta_l\bigr]=\sum_p c_p\bigl(\dot\zeta_l - s_p g_p\bigr) , #src
# ``` #src
# #src
# which is a Prony, or diffusive, representation of the viscous memory. The cost #src
# comparison is the point: marching the interior on a radial grid carries tens of #src
# unknowns per mode, and a few poles carry a handful -- for the same object, #src
# because the grid was only ever a device for computing ``\mathcal Z_l``. #src
# #src
# ### The small-Ohnesorge limit has a name #src
# #src
# At small ``\mathrm{Oh}`` the interior is potential flow plus a vortical layer of #src
# thickness ``\sqrt{\mathrm{Oh}}`` at the surface. A layer of that kind, driven by #src
# an oscillating boundary, contributes an impedance scaling as #src
# ``\sqrt\sigma`` -- a **half-order derivative in time**, which is the #src
# Basset--Boussinesq history term familiar from particle dynamics. So the regime #src
# where the instantaneous-eigenmode closure is *worst* is the regime whose memory #src
# kernel is best understood, and ``\sqrt\sigma`` is exactly the kernel that Prony #src
# sums represent efficiently. #src
# #src
# **What would have to be checked before any of this is used.** That #src
# ``\mathcal Z_l`` assembled from the interior solution reproduces Reid's #src
# characteristic equation through the balance above; that a two-pole fit of it #src
# returns ``\lambda_l`` and ``\omega_l^2``; and that the pole count converges. Each #src
# is a concrete test against a quantity the solver already computes, which is why #src
# this is written down as a route rather than adopted as one. #src

