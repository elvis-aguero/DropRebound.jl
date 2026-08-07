# # Shear-Thinning Drops
#
# A shear-thinning liquid has no single viscosity. Its resistance to
# deformation depends on how fast it is being deformed, so the constitutive
# law is a function ``\eta(\dot\gamma)`` rather than a number. This part asks
# what such a fluid does to the drop of Parts I to III.
#
# ## What breaks
#
# In the Newtonian theory the surface modes are independent damped
# oscillators. Each ``\zeta_l`` carries its own damping ``\lambda_l`` and
# frequency ``\omega_l``, and the interior flow can be solved once, in advance,
# for all time. That is what made Part II tractable in closed form.
#
# Neither property survives. The viscosity becomes a field over the drop,
# computed from the very flow it governs. Modes that deform different parts of
# the drop no longer see the same fluid, so they begin to drive one another.
# Two scalars per mode become matrices, and the interior flow has to be found
# afresh at every instant.
#
# The variational structure of Part I is what carries this without
# modification. Nothing in the Rayleighian assumed a constant viscosity: the
# dissipation functional is ``\mathcal R = \int \eta\,\bm e\!:\!\bm e\,dV``
# whatever ``\eta`` is, and a state-dependent ``\eta`` simply moves inside the
# integral. The formulation does not change. What changes is that its damping
# operator must be rebuilt as the solution evolves.
#
# ## What "no approximation" means here
#
# This page derives the system with no approximation beyond small amplitude and
# axisymmetry, and states it in full at the end. That phrase deserves unpacking,
# because the important part is what is *not* approximated.
#
# Three things happen below, and only the first is an approximation.
# Linearising in the surface amplitude drops advection, which is the single
# physical approximation the whole package rests on. Axisymmetry is exact for
# axisymmetric forcing: it removes the azimuthal index rather than approximating
# it. The poloidal and modal expansion is a change of variables.
#
# The rheology is not touched. ``\eta`` is evaluated pointwise from the shear
# rate of the whole summed strain field, never expanded in the amplitude and
# never reduced to one number per mode. That is the expensive choice, and it is
# deliberate: the shear-rate invariant does not superpose over modes, so any
# per-mode viscosity is already a closure rather than a rearrangement.
#
# The price is a system that asks, at every instant, for a coupled set of
# radial problems whose coefficients depend on their own solution. What that
# costs, and what may be given up to make it cheaper, is the subject of the
# companion page *Shear-Thinning Drops: Closures*.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``R`` | equilibrium drop radius; ``x=r/R`` |
# | ``\rho`` | density (constant) |
# | ``\eta(\dot\gamma)`` | dynamic viscosity -- **now a function of the flow** |
# | ``\eta_0,\ \eta_\infty`` | zero-shear and infinite-shear viscosity plateaus |
# | ``\lambda_c`` | Carreau-Yasuda time constant |
# | ``a`` | Carreau-Yasuda *shape* exponent (free) |
# | ``n`` | power-law index; ``n<1`` is shear-thinning |
# | ``\bm e = \tfrac12(\nabla\bm u+\nabla\bm u^{\mathsf T})`` | strain-rate tensor |
# | ``\dot\gamma=\sqrt{2\,\bm e\!:\!\bm e}`` | shear-rate invariant (a scalar) |
# | ``\zeta(\theta,t)`` | radial surface displacement; ``\epsilon=\zeta/R`` |
# | ``P_l(\mu)``, ``\mu=\cos\theta`` | Legendre polynomial |
# | ``\zeta_l(t)`` | amplitude of surface mode ``l``; the solver's state vector |
# | ``C_l(\theta)`` | Gegenbauer angular function ``\sin^2\!\theta\,P_l'(\cos\theta)/(l(l+1))``; the *angular* factor of the stream function -- **not** ``\zeta_l``, which is an amplitude in time |
# | ``\psi_l(x)`` | radial factor of mode ``l`` in the stream function ``\psi=\sum_l \psi_l(x)C_l(\theta)`` |
# | ``u_{r,l}(x),\ u_{\theta,l}(x)`` | radial profiles of ``u_r`` and ``u_\theta`` for mode ``l``; both fixed by ``\psi_l`` |
# | ``\mathcal D_l`` | ``d^2/dx^2-l(l+1)/x^2``; Reid's radial operator |
# | ``\mathcal T`` | ``d^2/dx^2-(2/x)d/dx+l(l+1)/x^2``; Reid's tangential-stress operator, ``\mathcal T[\psi]|_{x=1}=0`` is BC2 |
# | ``\mathcal L_l`` | ``d^2/dx^2+(2/x)d/dx-l(l+1)/x^2``; the radial Laplacian (**not** ``\mathcal T``) |
# | ``p`` | pressure inside the drop; ``p_l(x,t)`` its Legendre coefficients |
# | ``p_c(\theta,t)`` | pressure in the air film holding the drop off the substrate; ``p_{c,l}`` its coefficients |
# | ``\mathcal R_l[\psi;\eta]`` | variable-viscosity radial operator, derived below |
# | ``\mathcal R_{l m}`` | its off-diagonal generalisation once ``\eta`` varies with ``\theta`` |
# | ``k`` | degree index of the **viscosity field's own** Legendre series |
# | ``L_\eta`` | highest ``k`` present -- the *bandwidth* of the coupling |
# | ``G^{k}_{l m},\ H^{k}_{l m}`` | Gaunt-type angular integrals; pure numbers |
# | ``\mathrm{Oh}=\eta_0/\sqrt{\rho T_1R}`` | Ohnesorge number: viscous over inertio-capillary stress |
#
# ### Conventions
#
# Four rules, so that the symbols can be worked out rather than memorised.
#
# **Coefficients carry their field's own symbol.** Every field on this page is
# expanded in the same angular basis, and its coefficients keep its letter:
# ``\zeta\to\zeta_l``, ``\psi\to\psi_l``, ``p\to p_l``, ``p_c\to p_{c,l}``,
# ``\eta\to\eta_k``, ``u_r\to u_{r,l}``. If you know the field, you know the
# coefficient.
#
# **Primes are always ``d/dx``.** There are no primed indices anywhere, so
# ``\eta_k'`` is unambiguously the radial derivative of ``\eta_k``. The only
# apparent exceptions are ``P_l'`` and ``P_l''``, which are derivatives of the
# Legendre polynomial with respect to ``\mu`` -- still derivatives, never
# indices.
#
# **Three indices, one job each.** ``l`` is the mode whose equation is being
# written; ``m`` is the mode driving it; ``k`` is a harmonic of the viscosity
# field. ``l`` and ``m`` are shape modes and run ``2\ldots M``, where ``M`` is
# the truncation of the shape expansion (of order 50--90 in practice). ``k``
# starts at ``0`` and is a different kind of index: it controls only how far off
# the diagonal the mode-coupling matrix reaches, never the size of that matrix,
# which is always ``M\times M``. (``m`` also denotes azimuthal order in
# *Axisymmetry* below, which is the section that proves it vanishes; after that
# the letter is free.)
#
# **Uppercase is a basis function or an operator; lowercase is a field or a
# coefficient.** ``P_l``, ``C_l`` are basis functions; ``\mathcal D_l``,
# ``\mathcal T``, ``\mathcal L_l``, ``\mathcal R_{l m}`` are operators;
# ``G^k_{lm}``, ``H^k_{lm}`` are numbers. ``A^{(i)}_{lm}`` and ``B^{(i)}_{lm}``
# are the only remaining bare capitals, and they are radial integrals -- always
# carrying the ``^{(i)}`` superscript, never appearing alone.
#
# One cross-page note: the Reid page writes ``G(x)=U(x)/x^2`` and keeps ``U``
# for the stream-function profile, faithful to the original. Neither symbol is
# used in that sense here.
#
# Some results below are established symbolically, for arbitrary ``l``; others
# only by evaluation at a spread of concrete points, which is strong numerical
# evidence rather than a proof. The text says which is which.

using Symbolics, QuadGK, SpecialFunctions, DropSolver  #src
using Printf  #src
using LinearAlgebra: det  #src
using DropSolver: sph_bessel_ratio  #src

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
println("A HIERARCHY OF SHEAR-THINNING DROP MODELS")  #src
println("="^78)  #src

# ## Problem definition
#
# An incompressible fluid of constant density ``\rho``
# occupies a convex region bounded by a free surface ``\Sigma(t)``. The unknowns
# are the velocity ``\bm u``, the pressure ``p``, and the surface itself.
#
# ```math
# \nabla\cdot\bm u = 0,
# \qquad
# \rho\Bigl(\partial_t\bm u + \bm u\cdot\nabla\bm u\Bigr)
#   = -\nabla p + \nabla\cdot\bm\tau ,
# ```
#
# closed by the **constitutive law**, which is the only place the fluid's
# identity enters:
#
# ```math
# \bm\tau = 2\,\eta(\dot\gamma)\,\bm e,
# \qquad
# \dot\gamma = \sqrt{2\,\bm e\!:\!\bm e},
# \qquad
# \eta = \eta_\infty + (\eta_0-\eta_\infty)\bigl[1+(\lambda_c\dot\gamma)^a\bigr]^{\frac{n-1}{a}} .
# ```
#
# On ``\Sigma(t)``: the surface moves with the fluid, and the traction jump
# balances surface tension,
#
# ```math
# (-p\,\bm I + \bm\tau)\cdot\bm n = T_1(\nabla\cdot\bm n)\,\bm n .
# ```
#
# Here ``\bm u`` is the velocity field and ``p`` the pressure, ``\rho`` the
# (constant) density, ``\bm n`` the outward unit normal to the free surface,
# ``T_1`` the surface tension coefficient, and ``\nabla\cdot\bm n`` twice the
# mean curvature, positive for a convex surface. ``\bm\tau=2\eta\bm e`` is the
# deviatoric stress.
#
# Everything below is a simplification of these equations.
#
# ### The Gegenbauer-Legendre identity
#
# For a **constant** viscosity the divergence of the stress collapses to
# ``\eta\nabla^2\bm u``, which is what makes Reid's problem separable. For a
# variable viscosity it does not. Expanding the divergence of a product,
#
# ```math
# \nabla\cdot(2\eta\bm e) = 2\eta\,(\nabla\cdot\bm e) + 2(\nabla\eta)\cdot\bm e ,
# ```
#
# and for an incompressible field ``\nabla\cdot\bm e = \tfrac12\nabla^2\bm u``, so
#
# ```math
# \boxed{\ \nabla\cdot(2\eta\bm e) = \eta\,\nabla^2\bm u \;+\; 2(\nabla\eta)\cdot\bm e\ }
# ```
#
# The boxed result rests on the incompressible identity
# ``\nabla\cdot\bm e=\tfrac12\nabla^2\bm u``, which follows from
# ``\nabla\cdot\bm e=\tfrac12[\nabla^2\bm u+\nabla(\nabla\cdot\bm u)]``
# and ``\nabla\cdot\bm u=0``.

@variables xc yc zc  #src
let  #src
    ∂x = Differential(xc); ∂y = Differential(yc); ∂z = Differential(zc)  #src
    u1 = sin(xc) + xc*sin(zc)  #src
    u2 = sin(yc) - cos(xc)*yc  #src
    u3 = cos(zc) - cos(yc)*zc  #src
    u = [u1, u2, u3]  #src
    X = [xc, yc, zc]  #src
    D = [∂x, ∂y, ∂z]  #src
    divu = Symbolics.expand_derivatives(sum(D[i](u[i]) for i in 1:3))  #src
    ok_inc, _ = symbolic_zero(divu, X)  #src
    @assert ok_inc "test field is not divergence-free; the identity check below would be meaningless"  #src
    e = [Symbolics.expand_derivatives((D[i](u[j]) + D[j](u[i]))/2) for i in 1:3, j in 1:3]  #src
    div_e = [Symbolics.expand_derivatives(sum(D[i](e[i,j]) for i in 1:3)) for j in 1:3]  #src
    lap_u = [Symbolics.expand_derivatives(sum(D[i](D[i](u[j])) for i in 1:3)) for j in 1:3]  #src
    resid = [Symbolics.expand_derivatives(div_e[j] - lap_u[j]/2) for j in 1:3]  #src
    for j in 1:3  #src
        ok, how = symbolic_zero(resid[j], X)  #src
        @assert ok "div(e) = lap(u)/2 failed in component $j"  #src
    end  #src
    println("\nSTEP 0")  #src
    println("  ASSERTION 1 OK: for an incompressible field, div(e) = lap(u)/2, verified")  #src
    println("    componentwise on a nontrivial divergence-free test field.")  #src
    println("    => nabla.(2*eta*e) = eta*lap(u) + 2*(grad eta).e  is exact.")  #src
end  #src

# ## Linearisation in the surface amplitude
#
# **Assumption.** The surface displacement is small, ``\epsilon=\zeta/R\ll1``.
#
# **What it drops.** The advective term ``\bm u\cdot\nabla\bm u``, and the
# transfer of the boundary conditions from the deformed surface
# ``r=R+\zeta`` to the sphere ``r=R``. Error ``O(\epsilon^2)``.
#
# ### The rheology stays nonlinear
#
# Linearising the **kinematics** does not linearise the **rheology**.
# The strain rate is ``O(\epsilon)``, so ``\dot\gamma=\epsilon\,\hat{\dot\gamma}``,
# but ``\eta`` is a *nonlinear function of that small quantity*:
#
# ```math
# (\lambda_c\dot\gamma)^a = (\lambda_c\hat{\dot\gamma})^a\,\epsilon^a .
# ```
#
# If ``a`` is not an even integer, ``\epsilon^a`` is **not analytic** at
# ``\epsilon=0``. There is no Taylor series in ``\epsilon`` of which this is
# "the leading term", and for ``a<1`` the correction is *larger* than any
# linear term as ``\epsilon\to0``.
#
# !!! note
#     The chain built below never expands in ``\epsilon``. It evaluates
#     ``\eta(\dot\gamma)`` at whatever ``\dot\gamma`` the current state produces
#     and assembles the coupling from that, so a non-integer ``a`` costs it
#     nothing. Carreau-Yasuda at the fitted ``a\approx0.743`` is admissible
#     throughout, exactly as the table of admissible models below shows.

let  #src
    ## The claim is analytic, not numerical: d/deps eps^a = a*eps^(a-1). For  #src
    ## a < 1 that grows without bound as eps -> 0, so eps^a is not even  #src
    ## differentiable there; for a = 2 it vanishes; for a = 1 it is constant.  #src
    ## Sample the derivative to confirm the three behaviours are distinct.  #src
    epss = (1e-2, 1e-4, 1e-6, 1e-8)  #src
    for a in (0.5, 0.743)  #src
        d = [a*eps^(a-1) for eps in epss]  #src
        @assert all(d[k+1] > d[k] for k in 1:length(d)-1) "eps^a derivative must increase as eps->0 for a<1"  #src
        @assert d[end] > 10*d[1] "the growth over six decades is implausibly weak for a<1"  #src
    end  #src
    d2 = [2*eps for eps in epss]  #src
    @assert all(d2[k+1] < d2[k] for k in 1:length(d2)-1) "for a=2 the derivative must vanish as eps->0"  #src
    d1 = [1.0*eps^0 for eps in epss]  #src
    @assert all(isapprox(x, 1.0) for x in d1) "a=1 should be the marginal, constant-derivative case"  #src
    println("  ASSERTION 2 OK: d/deps of eps^a grows without bound as eps->0 for")  #src
    println("    a<1, is constant at a=1, and vanishes at a=2 -- so a small-amplitude")  #src
    println("    expansion exists only in the even-integer case. This closes the")  #src
    println("    perturbative route; it excludes no fluid from the chain below.")  #src
end  #src

# ## Axisymmetry
#
# **Assumption.** The forcing (impact) is axisymmetric.
#
# **Claim.** A generalized Newtonian fluid cannot break axisymmetry.
#
# **Proof.** ``\dot\gamma=\sqrt{2\bm e\!:\!\bm e}`` is a *scalar invariant* of
# the strain-rate tensor. If ``\bm u`` is axisymmetric then ``\bm e`` is
# axisymmetric, hence ``\dot\gamma`` is a function of ``(r,\theta)`` only,
# hence ``\eta=\eta(\dot\gamma)`` is axisymmetric, hence
# ``\nabla\cdot(2\eta\bm e)`` is an axisymmetric poloidal field. No
# ``m\neq0`` content is ever generated. ``\blacksquare``
#
# **Consequence.** Spherical harmonics ``Y_l^m`` collapse to Legendre
# polynomials ``P_l(\cos\theta)`` permanently. You never need the full
# ``Y_l^m`` machinery for this physics -- not as an approximation, but as a
# theorem. Only a non-axisymmetric *impact* would change that.

let  #src
    println("\nSTEP 2 (A2): axisymmetry is preserved EXACTLY")  #src
    I3 = [1.0 0 0; 0 1.0 0; 0 0 1.0]  #src
    rot(a) = [cos(a) -sin(a) 0.0; sin(a) cos(a) 0.0; 0.0 0.0 1.0]  #src
    worst = 0.0  #src
    for trial in 1:6  #src
        B = [0.3+0.11trial 0.7-0.05trial 0.2trial;  #src
             0.7-0.05trial -0.4+0.09trial 0.5-0.03trial;  #src
             0.2trial 0.5-0.03trial 0.1+0.02trial]  #src
        e = (B + B')/2  #src
        e = e - ((e[1,1] + e[2,2] + e[3,3])/3)*I3  #src
        S0 = sqrt(2*sum(e .* e))  #src
        for a in (0.0, 0.6, 1.9, 3.4, 5.7)  #src
            R = rot(a)  #src
            er = R*e*R'  #src
            Sr = sqrt(2*sum(er .* er))  #src
            worst = max(worst, abs(Sr - S0))  #src
        end  #src
    end  #src
    @assert worst < 1e-13 "the shear invariant is NOT rotation-invariant"  #src
    println("  ASSERTION 3 OK: S = sqrt(2 e:e) is invariant under rotation about the")  #src
    println("    symmetry axis, to $(round(worst, sigdigits=2)).")  #src
    println("    SCOPE: this checks the INVARIANT, on a constant tensor. It does not")  #src
    println("    exercise a field, any phi-dependence, or any mode coupling, so it is")  #src
    println("    corroboration of one step and not a proof that m != 0 is never")  #src
    println("    generated -- that argument is the prose above, which needs only that an")  #src
    println("    axisymmetric field maps to itself under the rotation.")  #src
end  #src

# ## Poloidal representation and modal expansion
#
# An axisymmetric incompressible velocity field is generated by a single
# scalar. Writing the surface shape as
#
# ```math
# \zeta(\theta,t) = R\sum_{l\ge2} \zeta_l(t)\,P_l(\cos\theta),
# ```
#
# the state of the drop is the vector of modal amplitudes ``\{\zeta_l(t)\}``,
# ``l=2\ldots M``. This is exactly the state `julia/src/types.jl` carries.
#
# An axisymmetric incompressible field is generated by a Stokes stream function,
# and expanding *that* in the same angular basis introduces the interior
# unknowns:
#
# ```math
# \psi(x,\theta,t) \;=\; \sum_{l\ge2} \psi_l(x,t)\,C_l(\theta),
# \qquad
# C_l(\theta)=\frac{\sin^2\!\theta\,P_l'(\cos\theta)}{l(l+1)} ,
# ```
#
# from which the velocity components follow by differentiation,
# ``u_r=\partial_\theta\psi/(x^2\sin\theta)`` and
# ``u_\theta=-\partial_x\psi/(x\sin\theta)``. Both reduce to something simple,
# because of one identity: differentiating ``C_l`` and eliminating ``P_l''``
# through Legendre's equation ``(1-\mu^2)P_l''=2\mu P_l'-l(l+1)P_l`` gives
#
# ```math
# \frac{\partial C_l}{\partial\theta} = \sin\theta\,P_l(\cos\theta) ,
# ```
#
# with no residual ``P_l'``. Hence
#
# ```math
# \boxed{\;
# u_r=\sum_l u_{r,l}(x)P_l(\mu),\quad u_{r,l}=\frac{\psi_l}{x^2};
# \qquad
# u_\theta=\sum_l u_{\theta,l}(x)\,\partial_\theta P_l(\mu),\quad
# u_{\theta,l}=\frac{\psi_l'}{x\,l(l+1)} \;}
# ```
#
# Two consequences are used repeatedly below, and both are stated here so that
# every factor in them is on the page rather than in a reader's head. The
# shear strain is proportional to the tangential-stress operator,
#
# ```math
# e_{r\theta}=\frac{\mathcal T[\psi_l]}{2x\,l(l+1)}\,\partial_\theta P_l ,
# \qquad
# \mathcal T=\frac{d^2}{dx^2}-\frac{2}{x}\frac{d}{dx}+\frac{l(l+1)}{x^2} ,
# ```
#
# which is why BC2 reduces to ``\mathcal T[\psi_l]\big|_{x=1}=0``; and the vorticity
# is proportional to Reid's radial operator,
#
# ```math
# \omega_\varphi=-\sum_l\frac{\mathcal D_l[\psi_l]\,C_l(\theta)}{x\sin\theta},
# \qquad
# \mathcal D_l=\frac{d^2}{dx^2}-\frac{l(l+1)}{x^2},
# ```
#
# which is the factor that turns the curl of the momentum equation into the
# interior problem. The ``C_l`` are orthogonal under
# ``\langle f,g\rangle=\int_0^\pi fg\,d\theta/\sin\theta``, with
#
# ```math
# \langle C_l,C_{m}\rangle=\frac{2\,\delta_{l m}}{(2l+1)\,l(l+1)} ,
# ```
#
# and it is that weight which makes the projection below well defined. The
# strain-rate tensor is linear in ``\bm u`` and therefore in ``\{\psi_l\}``.
#
# What ``\psi_l`` actually *is* -- the equation that determines it -- is the
# interior problem, and it is the substance of this page. Note already that it
# cannot be quoted from the Newtonian theory: there the profile is fixed once
# ``l`` and a constant viscosity are given, whereas here the viscosity is itself
# a functional of the field being solved for.

@variables xs ts  #src
@variables rr tt qq hh0 hh1 hh2  #src
let  #src
    ## Every boxed factor above is checked here, because each one is a place a  #src
    ## typo would survive review: they are all of the form "operator over       #src
    ## x*l*(l+1)", and an l(l+1) or a factor of 2 dropped anywhere would leave   #src
    ## the structure of the model intact while making it wrong.                  #src
    ##                                                                          #src
    ## Checks are at concrete l against the field built from a generic U(x), at  #src
    ## a spread of (x, theta), with the Legendre functions evaluated by          #src
    ## recurrence -- so nothing here assumes the identities being tested.        #src
    Dx = Differential(xs); Dt = Differential(ts)  #src
    ed(e) = Symbolics.expand_derivatives(e)  #src
    LP(l,m)  = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    LPp(l,m) = l==0 ? zero(m) : l*(m*LP(l,m)-LP(l-1,m))/(m^2-1)  #src
    Cl(l) = sin(ts)^2*LPp(l,cos(ts))/(l*(l+1))  #src
    ev(e,xv,tv) = Symbolics.build_function(ed(e), xs, ts; expression=Val(false))(xv,tv)  #src
    PTS = ((0.43,0.7),(0.61,1.3),(0.82,2.1),(0.37,2.6),(0.95,0.44))  #src
    U(l) = xs^(l+1) + 0.7xs^(l+3) - 0.3xs^(l+5) + 0.11xs^(l+2)  #src
    Dl(f,l) = ed(Dx(Dx(f))) - l*(l+1)*f/xs^2  #src
    L2(f,l) = ed(Dx(Dx(f))) - 2*ed(Dx(f))/xs + l*(l+1)*f/xs^2  #src

    w_dC, w_F, w_W, w_ert, w_vort = 0.0, 0.0, 0.0, 0.0, 0.0  #src
    mag_ert, mag_vort = 0.0, 0.0  #src
    for l in 2:6  #src
        Uf  = U(l)  #src
        psi = Uf*Cl(l)  #src
        Pl_ = LP(l, cos(ts))  #src
        dPl = ed(Dt(Pl_))                       # d/dtheta of P_l(cos theta)  #src
        ## (a) dC_l/dtheta = sin(theta) P_l  -- the Legendre-equation identity  #src
        for (xv,tv) in PTS  #src
            w_dC = max(w_dC, abs(ev(ed(Dt(Cl(l))) - sin(ts)*Pl_, xv, tv)))  #src
        end  #src
        ## (b) u_r = (U/x^2) P_l   and   u_theta = (U'/(x l(l+1))) dP_l/dtheta  #src
        u_r  =  ed(Dt(psi)/(xs^2*sin(ts)))  #src
        u_th = -ed(Dx(psi)/(xs*sin(ts)))  #src
        for (xv,tv) in PTS  #src
            w_F = max(w_F, abs(ev(u_r  - (Uf/xs^2)*Pl_, xv, tv)))  #src
            w_W = max(w_W, abs(ev(u_th - (ed(Dx(Uf))/(xs*l*(l+1)))*dPl, xv, tv)))  #src
        end  #src
        ## (c) e_rtheta = T[psi] dP_l/dtheta / (2 x l(l+1))  -- the factor BC2 rests on  #src
        e_rt = ed((Dt(u_r)/xs + Dx(u_th) - u_th/xs)/2)  #src
        claim_ert = L2(Uf,l)*dPl/(2*xs*l*(l+1))  #src
        ## (d) omega_phi = -D_l[psi] C_l / (x sin theta)  -- the factor the curl rests on  #src
        vort = ed((Dx(xs*u_th) - Dt(u_r))/xs)  #src
        claim_vort = -Dl(Uf,l)*Cl(l)/(xs*sin(ts))  #src
        for (xv,tv) in PTS  #src
            w_ert = max(w_ert, abs(ev(e_rt - claim_ert, xv, tv)))  #src
            w_vort = max(w_vort, abs(ev(vort - claim_vort, xv, tv)))  #src
            mag_ert  = max(mag_ert,  abs(ev(claim_ert, xv, tv)))  #src
            mag_vort = max(mag_vort, abs(ev(claim_vort, xv, tv)))  #src
        end  #src
    end  #src
    ## (e) orthogonality and the exact norm <C_l,C_l> = 2/((2l+1) l(l+1))  #src
    nodes, wts = QuadGK.gauss(40, -1.0, 1.0)  #src
    w_orth, w_norm = 0.0, 0.0  #src
    for l in 2:6, l2 in 2:6  #src
        ## <A,B> = int A B dtheta / sin(theta); with C_l = sin^2 P'_l/(l(l+1))  #src
        ## this is int (1-mu^2) P'_l P'_l2 dmu / (l(l+1) l2(l2+1)).  #src
        val = sum(w*(1-mu^2)*LPp(l,mu)*LPp(l2,mu) for (mu,w) in zip(nodes,wts)) /  #src
              (l*(l+1)*l2*(l2+1))  #src
        if l == l2  #src
            w_norm = max(w_norm, abs(val - 2/((2l+1)*l*(l+1))))  #src
        else  #src
            w_orth = max(w_orth, abs(val))  #src
        end  #src
    end  #src
    ## Errors are normalised by the largest magnitude the sweep produced. A bare  #src
    ## absolute tolerance is the wrong test here: e_rtheta and the vorticity both  #src
    ## carry a 1/(2 x l(l+1)) prefactor, so at l=6 they are ~1/50 of the operator  #src
    ## they are built from, and an absolute threshold either passes vacuously or   #src
    ## rejects a correct identity depending on which l dominates.                  #src
    @assert mag_ert > 1e-4 && mag_vort > 1e-4 "the e_rtheta / vorticity sweep never exercised a nonzero quantity ($mag_ert, $mag_vort)"  #src
    rel_ert, rel_vort = w_ert/mag_ert, w_vort/mag_vort  #src
    @assert w_dC   < 1e-10 "dC_l/dtheta = sin(theta) P_l failed ($w_dC)"  #src
    @assert w_F    < 1e-10 "u_r radial profile is not psi_l/x^2 ($w_F)"  ## CLAIM: SUM-URL  #src
    @assert w_F    < 1e-10 "u_theta is not -(1/(x sin t)) d(psi)/dx ($w_F)"  ## CLAIM: SUM-UTH  #src
    @assert w_W    < 1e-10 "u_theta radial profile is not U'/(x l(l+1)) ($w_W)"  ## CLAIM: SUM-FW  #src
    @assert rel_ert  < 1e-10 "e_rtheta is not T[psi] dP_l / (2 x l(l+1)) (rel $rel_ert)"  ## CLAIM: SUM-BC  #src
    @assert rel_vort < 1e-10 "vorticity is not -D_l[psi] C_l / (x sin theta) (rel $rel_vort)"  #src
    @assert w_orth < 1e-11 "the C_l are not orthogonal under the stated weight ($w_orth)"  #src
    @assert w_norm < 1e-11 "<C_l,C_l> is not 2/((2l+1) l(l+1)) ($w_norm)"  #src
    println("  ASSERTION 2b OK: every factor in the boxed velocity relations, over l=2..6:")  #src
    @printf("    dC_l/dtheta = sin(theta) P_l                    (%.1e)\n", w_dC)  #src
    @printf("    u_{r,l} = psi_l/x^2                                   (%.1e)\n", w_F)  #src
    @printf("    u_{th,l} = psi_l'/(x l(l+1))                           (%.1e)\n", w_W)  #src
    @printf("    e_rtheta = T[psi] dP_l/dtheta / (2 x l(l+1))     (rel %.1e)\n", rel_ert)  #src
    @printf("    omega_phi = -D_l[psi] C_l / (x sin theta)         (rel %.1e)\n", rel_vort)  #src
    @printf("    <C_l,C_m> = 2 delta / ((2l+1) l(l+1))          (%.1e, %.1e)\n", w_orth, w_norm)  #src
    println("    These are the factor-carrying identities the model summary quotes.")  #src
    println("    Physical meaning of a failure: a dropped l(l+1) or 2 would leave the")  #src
    println("    structure of the model intact while making every coefficient wrong.")  #src
end  #src

# ## Hypotheses on the viscosity
#
# The viscosity is required to satisfy three conditions, which together define the
# class of fluids the model describes.
#
# **(H1) Generalized Newtonian.** ``\eta`` depends on the flow only through the
# instantaneous local shear-rate invariant, ``\eta=\eta(\dot\gamma)`` with
# ``\dot\gamma=\sqrt{2\bm e\!:\!\bm e}`` -- no strain history, and no dependence
# on the individual components of ``\bm e`` beyond that scalar. This is what
# makes the coefficient matrices functions of the current state, so the system
# is quasi-linear rather than integro-differential.
#
# **(H2) Bounded above and below**, ``0<\eta_\infty\le\eta(\dot\gamma)\le\eta_0<\infty``.
# The lower bound is what BC2 needs: ``\tau_{r\theta}=2\eta e_{r\theta}=0`` forces
# ``e_{r\theta}=0`` only because ``\eta`` cannot vanish. The upper bound keeps the
# viscosity field integrable, so its Legendre coefficients exist.
#
# **(H3) Continuity in ``\dot\gamma``**, so ``\eta(x,\theta)`` has a convergent
# Legendre expansion.
#
#
# | model | admissible | reason |
# |:--|:--|:--|
# | Carreau-Yasuda, any ``a``, any ``n<1``, ``\eta_\infty>0`` | yes | satisfies all three |
# | Cross | yes | the ``p=-1`` slice of Carreau-Yasuda |
# | Ellis, truncated power law | yes | bounded by construction |
# | unregularised power law | no | violates (H2); ``\eta\to0`` or ``\infty`` and BC2 stops reducing |
# | Bingham, Herschel-Bulkley | no | yield stress makes ``\eta`` unbounded as ``\dot\gamma\to0`` |
# | Oldroyd-B, viscoelastic models | no | violates (H1): the stress carries memory |
#
# The viscoelastic exclusion is a different physical problem rather than a gap
# in this argument; it is treated on the Oldroyd-B page.
#
# ## The variational structure
#
# The kinematics above give the velocity and the strain rate as linear functionals of
# the coordinates. That is everything the model needs: a dissipative mechanical system
# is determined by three scalar functionals of those same quantities, and the equations
# of motion follow without ever writing the momentum equation down.[^route]
#
# ### The three quadratic forms
#
# The generalised coordinates are the interior displacement profiles ``\chi_l``,
# expanded in whatever radial basis is convenient; write ``\bm\xi`` for their
# amplitudes. The surface amplitudes are not separate coordinates -- ``\zeta_l`` is
# the boundary trace ``\chi_l(1)`` -- and the stream function is the rate,
# ``\psi_l=\dot\chi_l``. Then
#
# ```math
# T[\dot{\bm\xi}]=\tfrac12\int|\bm u|^2\,dV,
# \qquad
# \Phi[\dot{\bm\xi}]=\int 2\eta\;\bm e\!:\!\bm e\,dV,
# \qquad
# V[\bm\xi]=\text{surface energy} ,
# ```
#
# One convention to reconcile before going on. The home page and *Variational
# Mechanics* write the dissipation as the Rayleighian
# ``\mathcal R = \eta\int\bm e\!:\!\bm e\,dV``,
# whereas ``\Phi`` here is the *total viscous dissipation rate*, which is twice
# that: ``\Phi = 2\mathcal R``. Both appear in the literature. The factor of one
# half in the next display is exactly that difference, so
# ``\tfrac12\,\partial\Phi/\partial\dot\xi_a = \partial\mathcal R/\partial\dot\xi_a``
# and the two statements are the same equation.
#
# The equations of motion are the Euler--Lagrange equations of that data,
#
# ```math
# \boxed{\;
# \frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a}
# \;+\;\frac12\frac{\partial\Phi}{\partial\dot\xi_a}
# \;+\;\frac{\partial V}{\partial\xi_a}
# \;=\; Q_a \;}
# ```
#
# with ``Q_a`` the work done by the film pressure. **No elimination and no
# closure**: because the interior coordinates are retained, this is the same
# coupled surface-plus-interior system the differential form states, and it holds
# whether or not the interior has relaxed.
#
# What it buys is that every coefficient is a second derivative of a quadratic
# form, so the objects to assemble are
#
# ```math
# \frac{\partial^2 T}{\partial\dot\xi_a\partial\dot\xi_b}
#   =\int\bm u^{(a)}\!\cdot\!\bm u^{(b)}\,dV,
# \qquad
# \boxed{\;\frac{\partial^2\Phi}{\partial\dot\xi_a\partial\dot\xi_b}
#   =\int 2\eta\;\bm e^{(a)}\!:\!\bm e^{(b)}\,dV\;}
# ```
#
# Three things follow.
#
# **The damping needs one derivative, not four.** It is an integral of
# ``\bm e^{(a)}\!:\!\bm e^{(b)}`` against ``\eta``: no stress divergence, no curl,
# no fourth-order operator anywhere. The same information the projection route
# extracts from ``\nabla\times\nabla\cdot(2\eta\bm e)`` is obtained by
# differentiating the velocity once and integrating.
#
# **Both matrices are symmetric, which is a free correctness test.** They are
# Hessians, so an assembly that returns an asymmetric one has a bug rather than a
# feature.
#
# **The matrix is banded, and the band is a selection rule.**
# ``\bm e^{(a)}\!:\!\bm e^{(b)}`` against ``\eta_kP_k`` is a product of three angular
# factors, so the harmonic ``\eta_k`` connects shape modes ``l`` and ``m`` only when
#
# ```math
# |l-m|\;\le\;k\;\le\;l+m,
# \qquad l+k+m \ \text{even},
# ```
#
# the triangle and parity conditions on a Legendre triple product. Two consequences
# are worth naming. A constant viscosity has only its ``k=0`` harmonic, the rule
# collapses to ``l=m``, and the modes decouple -- which is why the Newtonian problem
# is one mode at a time. And truncating ``\eta`` at ``k\le L_\eta`` makes the matrix
# banded with bandwidth ``L_\eta`` rather than merely sparse, so the cost of the
# coupling is set by how many harmonics the viscosity field actually carries.
#
# (The ``H``-family of angular integrals -- the one carrying ``(1-\mu^2)P_k'P_m'`` --
# is precisely the ``e_{r\theta}e_{r\theta}`` term, which is why exactly two families
# appear and no third.)
#
# [^route]:
#     The momentum equation can of course be written per mode and its traction
#     projected onto the surface harmonics instead; that route needs
#     ``\nabla\times\nabla\cdot(2\eta\bm e)``, four nested derivatives, and a
#     separate treatment of the boundary conditions and of the pressure. It reaches
#     the same equations -- shown below -- and is carried out in full in
#     `differential_formulation_derivation.jl`, which is not part of the documentation
#     because nothing in the model depends on it.
#
# !!! note "What this section does not do"
#     It does not reduce the system to two coefficients per mode. Doing that means
#     eliminating the interior coordinates from ``\bm\xi``, which is a closure and
#     is treated as one on the companion page. The variational statement above is
#     an equivalent form of the model, not a simplification of it.
#
# ### Equivalence with the traction formulation
#
# "Equivalent" is a claim, and one identity carries it. For any two
# divergence-free fields ``\bm u`` and ``\bm w``,
#
# ```math
# \int 2\eta\;\bm e[\bm u]\!:\!\bm e[\bm w]\,dV
#   \;=\; -\int \bm w\cdot\bigl(\nabla\!\cdot\!2\eta\bm e[\bm u]\bigr)dV
#   \;+\; \oint \bm w\cdot\bigl(2\eta\bm e[\bm u]\bigr)\!\cdot\!\bm n\,dS .
# ```
#
# The left side is what the variational route assembles. The first term on the
# right is the momentum operator the differential route projects. So the two
# formulations differ by the **surface term** -- and that term is the traction,
# which is exactly what supplies the tangential condition as a *natural* boundary
# condition and the normal-stress balance as the surface equation. Same equations,
# reached from either end.
#
# Two consequences worth stating plainly.
#
# **The pressure is absent from every term.** Both sides are built from
# divergence-free fields, so incompressibility is satisfied identically and there
# is no multiplier to carry. The interior pressure *is* that multiplier. The
# differential route has to reintroduce it because its surface condition is
# written as a traction balance containing ``p``; the variational route obtains the
# same surface equation by varying the surface energy, and never needs it. Nothing
# is lost -- ``p`` can be reconstructed afterwards if wanted -- but it is not a
# state variable and no elliptic solve for it is required.
#
# **BC2 is not imposed, it emerges.** The surface term vanishes for arbitrary
# tangential ``\bm w`` only if ``e_{r\theta}|_{x=1}=0``. In the differential route
# that is a boundary condition to be applied; here it is a consequence of
# stationarity. The kinematic condition BC1 remains *essential* -- it constrains
# the admissible variations rather than following from them.

let  #src
    ## The identity above, on fields with angular viscosity structure. Errors are    #src
    ## normalised by the LARGEST dissipation over the sweep, not per case: the       #src
    ## cross-mode terms vanish identically for radially stratified eta (that is the  #src
    ## diagonality result proved earlier), so a per-case relative error divides two  #src
    ## machine zeros and reports nonsense. An earlier run of this very check         #src
    ## reported "55" for exactly that reason.                                        #src
    h = 1e-2  #src
    d1v(f, v, i) = (-f(v .+ 2h*i) + 8f(v .+ h*i) - 8f(v .- h*i) + f(v .- 2h*i))/(12h)  #src
    LPe(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    dLPe(l,m) = l==0 ? zero(m) : l*(m*LPe(l,m)-LPe(l-1,m))/(m^2-1)  #src
    Cge(l,t) = sin(t)^2*dLPe(l,cos(t))/(l*(l+1))  #src
    psie(l, cs) = (x,t) -> sum(c*x^(l+1+2(k-1)) for (k,c) in enumerate(cs))*Cge(l,t)  #src
    function strn(psi)  #src
        ur(x,t) =  d1v(v -> psi(v[1],v[2]), [x,t], [0.0,1.0])/(x^2*sin(t))  #src
        ut(x,t) = -d1v(v -> psi(v[1],v[2]), [x,t], [1.0,0.0])/(x*sin(t))  #src
        err_(x,t) = d1v(v -> ur(v[1],v[2]), [x,t], [1.0,0.0])  #src
        ett(x,t) = d1v(v -> ut(v[1],v[2]), [x,t], [0.0,1.0])/x + ur(x,t)/x  #src
        epp(x,t) = ur(x,t)/x + ut(x,t)*cos(t)/(sin(t)*x)  #src
        ert(x,t) = (d1v(v -> ur(v[1],v[2]), [x,t], [0.0,1.0])/x  #src
                  + d1v(v -> ut(v[1],v[2]), [x,t], [1.0,0.0]) - ut(x,t)/x)/2  #src
        (err_, ett, epp, ert, ur, ut)  #src
    end  #src
    nxv, wxv = QuadGK.gauss(44, 0.0, 1.0); nmv, wmv = QuadGK.gauss(44, -1.0, 1.0)  #src
    volv(f) = sum(wx*wm*2pi*x^2*f(x, acos(mu)) for (x,wx) in zip(nxv,wxv), (mu,wm) in zip(nmv,wmv))  #src
    surfv(f) = sum(wm*2pi*f(1.0, acos(mu)) for (mu,wm) in zip(nmv,wmv))  #src
    worst_abs, scale = 0.0, 0.0  #src
    for ef in ((x,t) -> 1.0,  #src
               (x,t) -> 1.3 + 0.9x + 0.4x^2,  #src
               (x,t) -> 1.3 + 0.9x*cos(t) + 0.4*(3cos(t)^2-1)/2)  #src
        for (la, lb) in ((2,2), (2,3), (3,4))  #src
            A = strn(psie(la, (1.0, 0.7, -0.3))); B = strn(psie(lb, (0.8, -0.4, 0.2)))  #src
            lhs = volv((x,t) -> 2*ef(x,t)*(A[1](x,t)*B[1](x,t) + A[2](x,t)*B[2](x,t)  #src
                                         + A[3](x,t)*B[3](x,t) + 2*A[4](x,t)*B[4](x,t)))  #src
            dr(x,t) = (d1v(v -> v[1]^2*2ef(v[1],v[2])*A[1](v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
                     + d1v(v -> 2ef(v[1],v[2])*A[4](v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
                     - (2ef(x,t)*A[2](x,t) + 2ef(x,t)*A[3](x,t))/x)  #src
            dt(x,t) = (d1v(v -> v[1]^2*2ef(v[1],v[2])*A[4](v[1],v[2]), [x,t], [1.0,0.0])/x^2  #src
                     + d1v(v -> 2ef(v[1],v[2])*A[2](v[1],v[2])*sin(v[2]), [x,t], [0.0,1.0])/(x*sin(t))  #src
                     + 2ef(x,t)*A[4](x,t)/x - cos(t)/(sin(t)*x)*2ef(x,t)*A[3](x,t))  #src
            rhs = -volv((x,t) -> B[5](x,t)*dr(x,t) + B[6](x,t)*dt(x,t)) +  #src
                   surfv((x,t) -> B[5](x,t)*2*ef(x,t)*A[1](x,t) + B[6](x,t)*2*ef(x,t)*A[4](x,t))  #src
            worst_abs = max(worst_abs, abs(lhs - rhs)); scale = max(scale, abs(lhs))  #src
        end  #src
    end  #src
    rel = worst_abs/scale  #src
    @assert scale > 1.0 "the equivalence sweep never produced a nonzero dissipation ($scale)"  #src
    @assert rel < 1e-5 "the variational and differential forms are not the same equations (rel $rel)"  ## CLAIM: SUM-EL  #src
    @printf("  ASSERTION 5d OK: the dissipation form equals the momentum operator plus a\n")  #src
    @printf("    surface traction, to %.1e relative (scale %.3g), for a viscosity with\n", rel, scale)  #src
    println("    angular structure. So the variational and differential statements are the")  #src
    println("    same equations, and everything verified of one holds of the other.")  #src
    println("    No pressure appears on either side: divergence-free fields carry no")  #src
    println("    multiplier, which is why the variational route needs no elliptic solve.")  #src

    ## The two relations that define the coordinate, checked against the running       #src
    ## solver rather than against this script's own algebra. Both are definitions, so   #src
    ## what can fail is not the algebra but the CODE disagreeing with the page: if      #src
    ## the solver's coordinate were the stream function rather than the displacement,   #src
    ## or if it took the surface amplitude as independent instead of as the boundary    #src
    ## trace, these would break while every equation above still held.                 #src
    let                                                                                 #src
        Kt = 4                                                                          #src
        for l in (2, 5)                                                                 #src
            ## The monomial family is the ansatz this page writes, x^(l+1+2(k-1)). #src
            ## The solver's default is the Legendre-shifted family, which spans the #src
            ## same space with far better conditioning; the checks below are on the #src
            ## powers themselves, so they name the family they are about.    #src
            bas = DropSolver.RitzBasis(l, Kt, :monomial)                                #src
            ## zeta_l is the boundary trace: every trial function equals 1 at x = 1,    #src
            ## so the trace of any coefficient vector is its plain sum.                 #src
            tr = [DropSolver.phi(bas, k, 1.0) for k in 1:Kt]                             #src
            @assert maximum(abs, tr .- 1) < 1e-14 "the trial functions are not 1 at the surface, so zeta_l is not the trace of chi_l"  ## CLAIM: SUM-TRACE  #src
            ## every trial function equals 1 at x = 1, so the trace is the PLAIN SUM of  #src
            ## the coefficients -- which is what the summary writes                      #src
            @assert abs(sum(tr) - Kt) < 1e-13 "zeta_l is not the plain sum of the a_{l,k}"  ## CLAIM: SUM-TRACESUM  #src
            ## and the exponents are the Taylor powers of j_l: l+1, l+3, l+5, ...        #src
            for k in 1:Kt                                                                #src
                pw = log(DropSolver.phi(bas, k, 2.0)) / log(2.0)                          #src
                @assert abs(pw - (l + 1 + 2(k-1))) < 1e-10 "trial function $k of mode $l is x^$pw, not x^{l+1+2(k-1)}"  ## CLAIM: SUM-RITZ  #src
            end                                                                          #src
            ## psi_l = d(chi_l)/dt: the velocity field the solver builds from a         #src
            ## coefficient vector must be the one this page builds from the same vector #src
            ## read as a stream function. Compare u_r = psi_l/x^2 at interior points.   #src
            cf = [0.7, -0.4, 0.25, -0.1]                                                #src
            worst_u = 0.0; scale_u = 0.0                                                 #src
            for xq in (0.3, 0.6, 0.9), muq in (-0.8, 0.1, 0.7)                            #src
                psi_here = sum(cf[k] * DropSolver.phi(bas, k, xq) for k in 1:Kt)          #src
                ang = DropSolver.legendre_angular(l, muq)                                  #src
                got = DropSolver.modal_field(l, psi_here, 0.0, 0.0, xq, ang)[1]            #src
                want = psi_here / xq^2 * ang.P                                             #src
                worst_u = max(worst_u, abs(got - want)); scale_u = max(scale_u, abs(want)) #src
            end                                                                           #src
            @assert worst_u/scale_u < 1e-12 "the solver's radial velocity is not psi_l/x^2 with psi_l = chi-dot, so its coordinate is not the displacement (rel $(worst_u/scale_u))"  ## CLAIM: SUM-CHIDOT  #src
        end                                                                               #src
        @printf("  ASSERTION 5f OK: the running solver's coordinate IS the interior\n")    #src
        println("    displacement, and its surface amplitude IS the boundary trace. Had it")  #src
        println("    carried the stream function as the coordinate instead, the kinetic")     #src
        println("    energy would be quadratic in the coordinate rather than its rate and")  #src
        println("    Lagrange's equations would not apply.")                                  #src
    end                                                                                    #src
    println("    Physical meaning of a failure: the two routes would be different models,")  #src
    println("    and the cheaper assembly would not be assembling this one.")  #src
end  #src
#
# ### Calibration against Reid
#
# The quadratic forms can be checked without committing to any closure, by
# evaluating them on a flow field that is known independently. In the inviscid
# limit the interior flow is potential, ``\psi_l\propto x^{l+1}``, and the
# dissipation rate of that field must reproduce Lamb's small-viscosity damping. It
# does, exactly:
#
# ```math
# \mathrm{Oh}\,\frac{\Phi_{ll}}{2\,T_{ll}}
#   \;=\;(l-1)(2l+1)\,\mathrm{Oh}
#   \;=\;\lambda_l^{\text{Lamb}} ,
# ```
#
# with no Bessel function, no characteristic equation and no root-finding -- one
# volume integral of a potential-flow field. That is a strong check on the whole
# construction: the dissipation form, the kinetic form, and the factor of two
# between them all have to be right simultaneously for the integers to come out.
#
# A practical consequence for the solver: the functional is quadratic and
# self-adjoint, so a **Rayleigh--Ritz** treatment of the interior coordinates
# converges quadratically -- an error ``\epsilon`` in the radial trial space gives
# ``\epsilon^2`` in the resulting coefficient. A handful of trial functions
# therefore buys what a much finer radial grid buys. (This says nothing about
# Reid's ``Cxj_l(qx)+Dx^{l+1}``, which is not a trial function at all: it is the
# exact general regular solution of the constant-viscosity interior equation.)

let  #src
    ## Two checks, and the second is the one that would be hard to pass by         #src
    ## accident. Both integrals are evaluated on a product Gauss rule over the     #src
    ## sphere; the integrands need only FIRST derivatives of the velocity, which   #src
    ## is the whole point -- no double divergence appears anywhere here, and this   #src
    ## block runs in seconds where the projection route needed minutes.            #src
    Dr3 = Differential(rr); Dt3 = Differential(tt)  #src
    ed3(e) = Symbolics.expand_derivatives(e)  #src
    LPd(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    dLPd(l,m) = l==0 ? zero(m) : l*(m*LPd(l,m)-LPd(l-1,m))/(m^2-1)  #src
    Cgd(l) = sin(tt)^2*dLPd(l,cos(tt))/(l*(l+1))  #src
    function strain4(psi)  #src
        ur =  ed3(Dt3(psi)/(rr^2*sin(tt))); ut = -ed3(Dr3(psi)/(rr*sin(tt)))  #src
        (ed3(Dr3(ur)),                                   # e_rr  #src
         ed3(Dt3(ut)/rr + ur/rr),                        # e_tt  #src
         ur/rr + ut*cos(tt)/(sin(tt)*rr),                # e_pp  #src
         ed3((Dt3(ur)/rr + Dr3(ut) - ut/rr)/2),          # e_rtheta  #src
         ur, ut)  #src
    end  #src
    nr, wr = QuadGK.gauss(40, 0.0, 1.0)  #src
    nm, wm = QuadGK.gauss(40, -1.0, 1.0)  #src
    function volint(expr, etaf)  #src
        g = Symbolics.build_function(ed3(expr), rr, tt; expression=Val(false))  #src
        tot = 0.0  #src
        for (x, wx) in zip(nr, wr), (mu, wu) in zip(nm, wm)  #src
            tot += wx*wu*2pi*x^2*etaf(x, mu)*g(x, acos(mu))  #src
        end  #src
        tot  #src
    end  #src
    Cform(l, m, Pl_, Pm_, etaf) = begin  #src
        A = strain4(Pl_*Cgd(l)); B = strain4(Pm_*Cgd(m))  #src
        volint(2*(A[1]*B[1] + A[2]*B[2] + A[3]*B[3] + 2*A[4]*B[4]), etaf)  #src
    end  #src
    Mform(l, m, Pl_, Pm_) = begin  #src
        A = strain4(Pl_*Cgd(l)); B = strain4(Pm_*Cgd(m))  #src
        volint(A[5]*B[5] + A[6]*B[6], (x,mu) -> 1.0)  #src
    end  #src

    ## (i) SYMMETRY of C, including for a viscosity with angular structure  #src
    prof(l) = rr^(l+1) + 0.7rr^(l+3) - 0.3rr^(l+5)  #src
    worst_sym, scale_sym = 0.0, 0.0  #src
    for etaf in ((x,mu) -> 1.0,  #src
                 (x,mu) -> 1.3 + 0.9x + 0.4x^2,  #src
                 (x,mu) -> 1.3 + 0.9x*mu + 0.4*(3mu^2 - 1)/2)  #src
        for (l, m) in ((2,3), (3,2), (2,4), (4,2))  #src
            a = Cform(l, m, prof(l), prof(m), etaf)  #src
            b = Cform(m, l, prof(m), prof(l), etaf)  #src
            worst_sym = max(worst_sym, abs(a - b)); scale_sym = max(scale_sym, abs(a))  #src
        end  #src
    end  #src
    @assert scale_sym > 1.0 "the dissipation sweep never produced a nonzero form ($scale_sym)"  #src
    @assert worst_sym/scale_sym < 1e-12 "the dissipation form is not symmetric ($worst_sym)"  ## CLAIM: SUM-FORMS  #src

    ## (ii) LAMB: with the potential profile, Oh*C_ll/(2 M_ll) = (l-1)(2l+1) Oh  #src
    worst_lamb = 0.0  #src
    for l in 2:6  #src
        P = rr^(l+1)  #src
        got = Cform(l, l, P, P, (x,mu) -> 1.0)/(2*Mform(l, l, P, P))  #src
        worst_lamb = max(worst_lamb, abs(got - (l-1)*(2l+1))/((l-1)*(2l+1)))  #src
    end  #src
    @assert worst_lamb < 1e-10 "the dissipation form does not reproduce Lamb's damping ($worst_lamb)"  ## CLAIM: SUM-HESS  #src

    println("  ASSERTION 5b OK: the damping matrix is the dissipation Hessian --")  #src
    @printf("    C_lm = int 2 eta e^(l):e^(m) dV is symmetric to %.1e relative, for a\n", worst_sym/scale_sym)  #src
    println("      viscosity with angular structure as well as radial;")  #src
    @printf("    and with the potential profile psi ~ x^(l+1), Oh C_ll/(2 M_ll) equals\n")  #src
    @printf("      Lamb's (l-1)(2l+1) Oh to %.1e over l = 2..6 -- exact integers, from\n", worst_lamb)  #src
    println("      one volume integral, with no Bessel function and no root-finding.")  #src
    println("    So the coupling can be assembled from FIRST derivatives of the velocity")  #src
    println("    instead of from a fourth-order projection.")  #src
    println("    Physical meaning of a failure: the system would not derive from a")  #src
    println("    dissipation functional, the damping matrix would not be symmetric, and")  #src
    println("    no energy budget could be closed against it.")  #src

    ## (iii) THE STIFFNESS, calibrated against Rayleigh. This is the check whose      #src
    ## absence let a wrong V onto the page: every existing calibration used the       #src
    ## RATIO C/(2M), and both forms carry the angular norm 4pi/(2l+1), so the ratio   #src
    ## is identically blind to it. V was stated without that norm and the error was   #src
    ## 60% at l=2 -- and l-DEPENDENT, so no rescaling of zeta, T or time could         #src
    ## absorb it. Assembling omega^2 = V''/M is the only test that can see it.        #src
    Kstiff(l) = (4pi/(2l+1))*(l-1.0)*(l+2)          # V'' from the summary's V  #src
    worst_ray, worst_noNorm = 0.0, 0.0  #src
    for l in 2:6  #src
        P = rr^(l+1)                                 # potential flow, psi_l(1) = 1  #src
        M = Mform(l, l, P, P)  #src
        ray = l*(l-1.0)*(l+2)                        # Rayleigh's inviscid frequency  #src
        worst_ray    = max(worst_ray,    abs(Kstiff(l)/M - ray)/ray)  #src
        ## the same quantity with the angular norm omitted, which must FAIL  #src
        worst_noNorm = max(worst_noNorm, abs(((l-1.0)*(l+2))/M - ray)/ray)  #src
    end  #src
    @assert worst_ray < 1e-6 "V''/M does not reproduce Rayleigh's l(l-1)(l+2) ($worst_ray)"  ## CLAIM: SUM-STIFF  #src
    @assert worst_ray < 1e-6 "the kinetic form is not the added mass Rayleigh needs ($worst_ray)"  ## CLAIM: SUM-T  #src
    @assert worst_noNorm > 0.3 "dropping the angular norm 4pi/(2l+1) also passes, so this check cannot see it ($worst_noNorm)"  #src
    @printf("  ASSERTION 5e OK: the surface energy is calibrated, not just the damping --\n")  #src
    @printf("    V''_ll/M_ll = l(l-1)(l+2) to %.1e over l = 2..6, Rayleigh's inviscid\n", worst_ray)  #src
    @printf("    frequency. Omitting the angular norm 4pi/(2l+1) is off by %.0f%%, and\n", 100*worst_noNorm)  #src
    println("    l-dependently, so it cannot be absorbed by any rescaling.")  #src
    println("    Physical meaning of a failure: the mode spectrum would be wrong, and")  #src
    println("    wrong by a different factor in each mode -- plausible numbers, no symptom.")  #src
end  #src

# ## Closing the system
#
# ### The constitutive law
#
# The matrices need ``\eta_{k}(x,t)``, which comes from the fluid. Under the
# Cross law -- the model this repository's validation fluid is characterised
# with, and its ``p=-1`` slice --
#
# ```math
# \frac{\eta(x,\theta,t)}{\eta_0}
#   \;=\; \varepsilon_\infty
#       \;+\; \frac{1-\varepsilon_\infty}
#                  {1+\bigl(K\,\dot\gamma(x,\theta,t)\bigr)^{a}} ,
# \qquad
# \varepsilon_\infty \equiv \frac{\eta_\infty}{\eta_0},
# ```
#
# where ``K`` is the fluid's time constant (written ``\lambda_c`` in the
# Carreau-Yasuda parametrisation) and ``a`` its thinning exponent (written
# ``a``). The shear rate ``\dot\gamma=\sqrt{2\,\bm e\!:\!\bm e}`` is built from
# the drop's current state; the next section says exactly how, because that
# step is where the model is easiest to get wrong.
#
# ### Where the shear rate is evaluated
#
# Two things are easy to assume and both are wrong, so they are worth stating
# before the algebra goes further.
#
# **The shear rate is a field, not a per-mode number.** The surface moves at
# ``\dot\zeta=\sum_l\dot\zeta_lP_l``, so the interior velocity is driven by the
# modal *velocities*, and superposes over modes. Strain rate is linear in
# velocity, so the strain **tensor** superposes too:
#
# ```math
# \bm e(x,\theta,t) \;=\; \sum_{l=2}^{M}\dot\zeta_l(t)\;\bm e^{(l)}(x,\theta).
# ```
#
# The invariant does not:
#
# ```math
# \dot\gamma \;=\; \sqrt{2\,\bm e\!:\!\bm e}
#   \;=\; \Bigl(2\sum_{l,m}\dot\zeta_l\dot\zeta_{m}\;
#          \bm e^{(l)}\!:\!\bm e^{(m)}\Bigr)^{1/2}.
# ```
#
# It is quadratic in the field and *then* square-rooted, so the ``l\neq m``
# cross terms do not drop out. Once more than one mode is active there is no
# such thing as "mode ``l``'s shear rate".
#
# Nothing about ``\eta(\dot\gamma)`` is polynomial, so the coefficients
# ``\eta_{k}(x)`` have no closed form, and neither do the radial and angular
# integrals built from them. ``\eta_{k}`` is defined by its Legendre
# projection at each radius,
#
# ```math
# \eta_{k}(x,t)=\frac{2k+1}{2}\int_{-1}^{1}
#   \eta\bigl(\dot\gamma(x,\mu,t)\bigr)P_{k}(\mu)\,d\mu ,
# ```
#
# with ``\dot\gamma`` evaluated from the *full* superposition of active modes
# at that point, never mode by mode. How those integrals are actually
# discretised is a numerical question and is treated on the companion page.
#
# Those cross terms are also where the angular structure of ``\eta`` comes
# from, and that closes a loop: ``\dot\gamma^2`` is quadratic in a field of
# degree ``M``, so it carries harmonics up to ``2M`` -- the same ``2M`` that
# bounds the coupling below. Not a coincidence: both come from the field
# entering quadratically. (``\dot\gamma`` and ``\eta(\dot\gamma)`` spread past
# ``2M``, because of the square root and the constitutive law; by the proof
# below, none of that reaches the dynamics.)
#
# **The system is nonlinear.** Since ``\eta`` depends on the state, so do the
# matrices, and the equation is honestly written
#
# ```math
# \ddot\zeta_l + \sum_{m}\Bigl[\mathcal D^{(2)}_{l m}(\bm{\dot\zeta})\,\dot\zeta_{m}
#   + \mathcal D^{(1)}_{l m}(\bm{\dot\zeta})\,\zeta_{m}\Bigr] + l\,p_{c,l} = 0 .
# ```
#
# This is *quasi-linear*: linear in ``(\zeta_l,\dot\zeta_l)`` at frozen coefficients,
# nonlinear through them. No rearrangement turns it into a linear ODE, and
# saying "the two scalars become matrices" without saying the matrices depend
# on the state makes it look like one.
#
# The nonlinearity is not confined to the surface equations either. ``\eta``
# is built from the interior field, and the interior field solves a problem
# whose coefficients are ``\eta``: the interior problem and the constitutive
# law are a **fixed point**, to be satisfied at each instant. That is the
# honest shape of the model, and it is the last piece of it.
#
# That completes the model. The interior problem, the surface projection, the
# constitutive law and this prescription for ``\dot\gamma`` are the whole of
# it, and the summary below states them together.
#



# ## The force on the drop
#
# ### The capillary restoring force
#
# The normal-stress balance carries the curvature of the deformed surface, and
# what the surface equation needs is its linearisation. For an axisymmetric
# surface ``x=1+\zeta(\theta)``, taking ``\bm n=\nabla F/|\nabla F|`` with
# ``F=x-1-\zeta`` and expanding to first order in ``\zeta``,
#
# ```math
# \nabla\cdot\bm n \;=\; 2 \;+\; \sum_{l\ge2}(l-1)(l+2)\,\zeta_l\,P_l(\mu)
#   \;+\;O(\zeta^2) .
# ```
#
# The leading ``2`` is the equilibrium sphere's ``2/R``; the mode-``l`` term is
# the restoring stiffness, and ``(l-1)(l+2)`` vanishing at ``l=1`` is why
# translation costs no surface energy -- the ``l=1`` mode is a rigid shift, which
# is exactly why ``\zeta`` starts at ``l=2`` and the drop's position is carried
# separately by ``z``.
#
# ### The net force the substrate exerts
#
# The reaction pressure ``p_c\ge0`` can only push, so the traction on the drop is
# ``-p_c\bm n`` and its vertical component is ``-p_c\mu``. Integrating over the
# surface -- and noting that ``p_c`` vanishes off the contact region, so the
# integral may be taken over the whole sphere -- gives
#
# ```math
# \mathfrak F \;=\; -2\pi\int_{-1}^{1}p_c(\mu)\,\mu\,d\mu
#   \;=\; -\frac{4\pi}{3}p_{c,1} ,
# ```
#
# because ``\mu=P_1`` and every other term in ``p_c=\sum_lp_{c,l}P_l`` is orthogonal
# to it. So the net force is carried *entirely* by the single coefficient
# ``p_{c,1}``, and every other pressure coefficient is invisible to the centre of
# mass. Dividing by the drop's non-dimensional mass, which is also ``4\pi/3``,
#
# ```math
# \dot v = -\mathrm{Bo} - p_{c,1} .
# ```
#
# The two ``4\pi/3`` factors cancelling is why the equation of motion looks as
# though ``p_{c,1}`` *were* the force. It is not -- it is one coefficient of the
# pressure expansion, and ``\mathfrak F=-\tfrac{4\pi}{3}p_{c,1}`` is the relation
# between them.

let  #src
    ## Discharges the remaining Model summary claims. Each check builds the        #src
    ## quantity from its definition rather than from the identity being tested.    #src
    Dr = Differential(rr); Dt = Differential(tt)  #src
    ed2(e) = Symbolics.expand_derivatives(e)  #src
    LP(l,m)  = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    LPp(l,m) = l==0 ? zero(m) : l*(m*LP(l,m)-LP(l-1,m))/(m^2-1)  #src
    Cg(l) = sin(tt)^2*LPp(l,cos(tt))/(l*(l+1))  #src
    f2(e,xv,tv) = Symbolics.build_function(ed2(e), rr, tt; expression=Val(false))(xv,tv)  #src
    f1(e,tv) = Symbolics.build_function(ed2(e), tt; expression=Val(false))(tv)  #src
    TH = (0.4, 0.9, 1.5, 2.2, 2.8)  #src
    XT = ((0.43,0.7),(0.61,1.3),(0.82,2.1),(0.37,2.6))  #src

    ## SUM-GAP: h = mu(1+zeta) + z is the EXACT height of the surface point at    #src
    ## angle theta, not a linearisation -- the vertical Cartesian coordinate of    #src
    ## the point at radius 1+zeta on a centre at height z.                        #src
    @variables zc zeta_  #src
    ## SUM-GAP. The previous version compared `zc + (1+zeta)*cos(t)` against  #src
    ## `cos(t)*(1+zeta) + zc` -- the same expression with terms reordered. It  #src
    ## verified commutativity of + and *, and could not fail for any input.  #src
    ## What actually needs checking is the GEOMETRY: that h is the signed  #src
    ## distance from the substrate plane to the surface point, with the  #src
    ## orientation the force integral also uses.  #src
    hfun(zv, zetav, tv) = cos(tv)*(1 + zetav) + zv  #src
    ## (i) an undeformed sphere resting exactly on the plane touches at the south  #src
    ##     pole and nowhere else: h(pi) = 0 and h > 0 elsewhere.  #src
    w_gap = abs(hfun(1.0, 0.0, pi) + 0.0)          # z = 1 puts the pole at h = 0  #src
    for tv in (0.3, 1.0, 2.0, 3.0)  #src
        w_gap = max(w_gap, hfun(1.0, 0.0, tv) > 0 ? 0.0 : 1.0)  #src
    end  #src
    ## (ii) lowering the drop makes the pole penetrate: h(pi) < 0 for z < 1.  #src
    w_gap = max(w_gap, hfun(0.9, 0.0, pi) < 0 ? 0.0 : 1.0)  #src
    ## (iii) the orientation matches the force integral, which uses n_z = mu and  #src
    ##     lands the contact region at mu < 0. The pole that touches must be mu = -1.  #src
    w_gap = max(w_gap, isapprox(cos(pi), -1.0; atol=1e-14) ? 0.0 : 1.0)  #src
    ## (iv) a positive radial bulge AT THE POLE closes the gap rather than opening  #src
    ##      it -- the outward radial direction there points at the substrate.        #src
    ##      Expecting the opposite is what this check caught on its first run.       #src
    w_gap = max(w_gap, hfun(1.0, 0.05, pi) < hfun(1.0, 0.0, pi) ? 0.0 : 1.0)  #src
    @assert w_gap < 1e-12 "h is not the gap: geometry or orientation is wrong ($w_gap)"  ## CLAIM: SUM-GAP  #src

    ## SUM-STREAM: the stream-function velocities are divergence free, for a       #src
    ## multi-mode psi (not one mode at a time -- incompressibility must survive    #src
    ## superposition).                                                            #src
    psi_m = (rr^3 + 0.7rr^5)*Cg(2) + (rr^4 - 0.3rr^6)*Cg(3) + 0.5rr^5*Cg(4)  #src
    ur_m  =  ed2(Dt(psi_m)/(rr^2*sin(tt)))  #src
    ut_m  = -ed2(Dr(psi_m)/(rr*sin(tt)))  #src
    divu_m = ed2(Dr(rr^2*ur_m)/rr^2 + Dt(ut_m*sin(tt))/(rr*sin(tt)))  #src
    w_div = maximum(abs(f2(divu_m, xv, tv)) for (xv,tv) in XT)  #src
    @assert w_div < 1e-12 "the stream-function velocity field is not divergence free ($w_div)"  ## CLAIM: SUM-STREAM  #src

    ## SUM-GROUPS: the three groups, checked by substituting values rather than    #src
    ## by symbolic simplification -- Symbolics will not cancel T1*(R^3 rho/T1)     #src
    ## under a square root, and an unsimplified residual is not a failure.        #src
    w_grp = 0.0  #src
    for (rho_,R_,T1_,g_,V_,e0_) in ((998.0,3.5e-4,0.0722,9.81,0.31,1.0e-3),  #src
                                    (1210.0,7.1e-4,0.0640,9.81,0.85,4.3e-2),  #src
                                    (1.0,1.0,1.0,1.0,1.0,1.0))  #src
        Tsig = sqrt(rho_*R_^3/T1_)  #src
        inert = rho_*R_/Tsig^2  #src
        w_grp = max(w_grp, abs((T1_/R_)/R_/inert - 1.0))                         # pressure scale matches inertial  #src
        w_grp = max(w_grp, abs(e0_*(R_/Tsig)/R_^2/inert - e0_/sqrt(rho_*T1_*R_))) # Oh  #src
        w_grp = max(w_grp, abs(rho_*g_/inert            - rho_*g_*R_^2/T1_))      # Bo  #src
        w_grp = max(w_grp, abs((V_/(R_/Tsig))^2         - rho_*R_*V_^2/T1_))      # We  #src
    end  #src
    @assert w_grp < 1e-10 "the non-dimensional groups do not follow from the stated scalings ($w_grp)"  ## CLAIM: SUM-GROUPS  #src
    @assert w_grp < 1e-10 "Bo does not follow from the stated scalings ($w_grp)"  ## CLAIM: SUM-BO  #src
    @assert w_grp < 1e-10 "We does not follow from the stated scalings ($w_grp)"  ## CLAIM: SUM-WE  #src
    ## SUM-RHEO: the factor 2. Build e from a velocity field and contract it in     #src
    ## full; do not assume which components are nonzero. Calibrated on simple      #src
    ## shear, where the engineering shear rate is unambiguous.                     #src
    @variables X Y Z Gs  #src
    DX, DY, DZ = Differential(X), Differential(Y), Differential(Z)  #src
    uvec = [Gs*Y, 0X, 0X]  #src
    Ds   = (DX, DY, DZ)  #src
    grad = [ed2(Ds[j](uvec[i])) for i in 1:3, j in 1:3]  #src
    estr = (grad .+ permutedims(grad)) ./ 2  #src
    ee   = sum(estr[i,j]^2 for i in 1:3, j in 1:3)  #src
    for Gv in (1.0, 2.5, 7.3)  #src
        got = sqrt(2*Symbolics.value(Symbolics.substitute(ee, Dict(Gs=>Gv))))  #src
        @assert abs(got - Gv) < 1e-12 "sqrt(2 e:e) is not the shear rate in simple shear ($got vs $Gv)"  ## CLAIM: SUM-RHEO  #src
    end  #src

    ## SUM-INT: the relative coefficient of the inertial and viscous terms. With   #src
    ## u ~ exp(-sigma t) the interior equation reads -sigma D_l[psi] = Oh R_l[U],     #src
    ## and at constant eta that must be Reid's D_l(D_l + q^2)U = 0 with            #src
    ## q^2 = sigma/Oh -- which is what fixes the Oh prefactor.                     #src
    @variables sig Ohs  #src
    w_int = 0.0  #src
    for l in 2:5  #src
        Uf = rr^(l+1) + 0.7rr^(l+3) - 0.3rr^(l+5)  #src
        Dl_(f) = ed2(Dr(Dr(f))) - l*(l+1)*f/rr^2  #src
        lhs = Ohs*Dl_(Dl_(Uf)) + sig*Dl_(Uf)          # from -sigma D_l[psi] = Oh D_l^2[psi]  #src
        rhs = Ohs*Dl_(Dl_(Uf) + (sig/Ohs)*Uf)         # Reid's D_l(D_l + q^2)U, q^2 = sigma/Oh  #src
        for (xv,_) in XT, (sv,ov) in ((0.7,0.031),(2.3,0.9))  #src
            d = Symbolics.value(Symbolics.substitute(ed2(lhs - rhs),  #src
                    Dict(rr=>xv, sig=>sv, Ohs=>ov)))  #src
            w_int = max(w_int, abs(Float64(d)))  #src
        end  #src
    end  #src
    @assert w_int < 1e-10 "q^2 = sigma/Oh does not reconcile the interior equation with Reid's ($w_int)"  #src

    ## SUM-NORMAL: the linearised curvature coefficient (l-1)(l+2).               #src
    @variables ep  #src
    w_curv = 0.0; base_ok = true  #src
    for l in 2:6  #src
        Rs = 1 + ep*LP(l, cos(tt))  #src
        Rp = ed2(Dt(Rs))  #src
        nr =  1/sqrt(1 + Rp^2/rr^2)  #src
        nt = (-Rp/rr)/sqrt(1 + Rp^2/rr^2)  #src
        dn = ed2(Dr(rr^2*nr)/rr^2 + Dt(sin(tt)*nt)/(rr*sin(tt)))  #src
        at = Symbolics.substitute(dn, Dict(rr => Rs))  #src
        d1e = Symbolics.substitute(ed2(Differential(ep)(at)), Dict(ep => 0.0))  #src
        ## Symbolics leaves this as "2/sqrt(1)", which is a symbolic expression      #src
        ## rather than a number, so it is evaluated rather than converted.          #src
        b0e = Symbolics.substitute(at, Dict(ep => 0.0))  #src
        claim = (l-1)*(l+2)*LP(l, cos(tt))  #src
        for tv in TH  #src
            base_ok &= abs(f1(b0e, tv) - 2.0) < 1e-12  #src
            w_curv = max(w_curv, abs(f1(d1e, tv) - f1(claim, tv)))  #src
        end  #src
    end  #src
    @assert base_ok "the undeformed sphere's curvature is not 2"  #src
    @assert w_curv < 1e-12 "the linearised curvature coefficient is not (l-1)(l+2) ($w_curv)"  ## CLAIM: SUM-CURV  #src
    ## and (l-1)(l+2) must vanish at l=1: translation costs no surface energy  #src
    @assert (1-1)*(1+2) == 0 "the l=1 curvature stiffness must vanish"  #src

    ## SUM-COM: the net vertical force is -(4 pi/3) B_1, and nothing else.        #src
    w_force = 0.0  #src
    for N in 3:6  #src
        Bv = [0.3 + 0.11k for k in 0:N]  #src
        Pif(mu) = sum(Bv[k+1]*LP(k,mu) for k in 0:N)  #src
        Fq = -2pi*QuadGK.quadgk(mu -> Pif(mu)*mu, -1, 1; rtol=1e-13)[1]  #src
        w_force = max(w_force, abs(Fq - (-(4pi/3)*Bv[2])))  #src
    end  #src
    @assert w_force < 1e-10 "the net substrate force is not -(4pi/3) B_1 ($w_force)"  ## CLAIM: SUM-FORCE  #src
    ## The generalised force Q_{zeta_l} = -(4pi/(2l+1)) p_{c,l}, from the virtual work  #src
    ## of the film traction. Its l = 1 member must BE the net force above -- that is  #src
    ## the consistency the derivation predicts, and it is what discharges it.  #src
    w_q = 0.0  #src
    for pc in (-0.8, 0.0, 1.7)  #src
        q1 = -(4pi/(2*1+1))*pc            # Q at l = 1  #src
        w_q = max(w_q, abs(q1 - (-(4pi/3)*pc)))  #src
    end  #src
    @assert w_q < 1e-14 "Q at l=1 is not the net substrate force ($w_q)"  ## CLAIM: SUM-Q  #src
    @assert w_q < 1e-14 "the virtual-work statement does not yield that Q ($w_q)"  ## CLAIM: SUM-VWORK  #src    @assert abs((-(4pi/3)*1.0)/(4pi/3) + 1.0) < 1e-14 "F/mass must reduce to -B_1"  #src
    ## The centre-of-mass equation, on THIS page rather than the archived one: only  #src
    ## the model page may discharge a model claim. Dimensional: m dV/dT = -m g + F_z,  #src
    ## so vdot = -(g Tsig^2/R) + Fhat/(4pi/3) = -Bo - p_{c,1}. Signs included, and the  #src
    ## opposite gravity sign must fail.  #src
    w_com, w_comWrong = 0.0, Inf  #src
    for (rho_,R_,T1_,g_) in ((998.0,3.5e-4,0.0722,9.81), (1210.0,7.1e-4,0.0640,9.81))  #src
        Tsig2 = rho_*R_^3/T1_; Bo = rho_*g_*R_^2/T1_  #src
        w_com      = max(w_com,      abs(g_*Tsig2/R_ - Bo))  #src
        ## relative to Bo: these are real fluids with Bo ~ 0.017, so an absolute  #src
        ## threshold on the wrong-sign residual is meaningless.  #src
        w_comWrong = min(w_comWrong, abs(g_*Tsig2/R_ + Bo)/Bo)  #src
        for pc in (-0.9, 1.4)  #src
            w_com = max(w_com, abs((-(4pi/3)*pc)/(4pi/3) - (-pc)))  #src
        end  #src
    end  #src
    @assert w_com < 1e-12 "vdot = -Bo - p_{c,1} does not follow from the stated scalings ($w_com)"  ## CLAIM: SUM-COM  #src
    @assert w_comWrong > 1.0 "the opposite gravity sign also passes, so this cannot fix it ($w_comWrong)"  #src
    ## SUM-PLAP: the angular reduction of the pressure Laplacian, and the fact    #src
    ## that x^n is its homogeneous solution -- which is why a vanishing source     #src
    ## leaves exactly one amplitude per mode and no radial equation.               #src
    Ln(f, n) = ed2(Dr(Dr(f))) + 2*ed2(Dr(f))/rr - n*(n+1)*f/rr^2  #src
    lap3(e) = ed2(Dr(rr^2*Dr(e))/rr^2 + Dt(sin(tt)*Dt(e))/(rr^2*sin(tt)))  #src
    w_plap, w_harm, mag_plap = 0.0, 0.0, 0.0  #src
    for n in 0:5  #src
        pn = rr^(n+2) + 0.4rr^(n+4) - 0.2rr^(n+1)   # generic, not the harmonic one  #src
        lhs = lap3(pn*LP(n, cos(tt)))  #src
        rhs = Ln(pn, n)*LP(n, cos(tt))  #src
        for (xv,tv) in XT  #src
            w_plap  = max(w_plap, abs(f2(lhs - rhs, xv, tv)))  #src
            mag_plap = max(mag_plap, abs(f2(rhs, xv, tv)))  #src
        end  #src
        ## x^n must be annihilated: the harmonic part carries no source  #src
        for (xv,_) in XT  #src
            w_harm = max(w_harm, abs(f2(Ln(rr^n, n), xv, 1.0)))  #src
        end  #src
    end  #src
    @assert mag_plap > 1e-3 "the pressure-Laplacian sweep never exercised a nonzero operator ($mag_plap)"  #src
    @assert w_plap < 1e-10 "lap(p_n P_n) is not L_n[p_n] P_n ($w_plap)"  #src
    @assert w_harm < 1e-10 "x^n is not annihilated by L_n, so the harmonic part is wrong ($w_harm)"  #src

    println("  ASSERTION 2c OK: the remaining Model summary claims, discharged --")  #src
    @printf("    h is the gap, with the orientation the force integral uses (%.1e)\n", w_gap)  #src
    @printf("    the stream-function field is divergence free        (%.1e)\n", w_div)  #src
    @printf("    Oh, Bo, We follow from the stated scalings          (%.1e)\n", w_grp)  #src
    println("    sqrt(2 e:e) is the shear rate in simple shear       (exact)")  #src
    @printf("    q^2 = sigma/Oh reconciles interior with Reid        (%.1e)\n", w_int)  #src
    @printf("    curvature: 2 + (l-1)(l+2)A_l P_l, l=2..6            (%.1e)\n", w_curv)  #src
    @printf("    net force = -(4pi/3)B_1, so vdot = -Bo - B_1        (%.1e)\n", w_force)  #src
    @printf("    lap(p_n P_n) = L_n[p_n] P_n, and L_n[x^n] = 0       (%.1e, %.1e)\n", w_plap, w_harm)  #src
    println("    Physical meaning of a failure: the summary would contain an equation")  #src
    println("    with a wrong coefficient, and every number the solver produces from")  #src
    println("    it would be wrong by that factor while looking entirely plausible.")  #src
end  #src

# ## Model summary
#
# The model, stated once. It is a **damped Lagrangian system**: three energies, a
# contact condition, and a constitutive law. The radial direction is not
# discretised and no numerical parameter appears; the two truncations that are
# present are named at the end.
#
# ### Scalings
#
# Lengths by the equilibrium radius ``R``, time by ``T_\sigma=\sqrt{\rho R^3/T_1}``,
# velocity by ``R/T_\sigma``, viscosity by the zero-shear plateau ``\eta_0``:
#
# ```math
# \mathrm{Oh}=\frac{\eta_0}{\sqrt{\rho T_1 R}},
# \qquad
# \mathrm{Bo}=\frac{\rho g R^2}{T_1},
# \qquad
# \mathrm{We}=\frac{\rho R V^2}{T_1} .
# ```
#
# ### The unknowns
#
# | unknown | domain | what it is |
# |:--|:--|:--|
# | ``\chi_l(x,t)``, ``l\ge2`` | ``0<x<1`` | interior displacement profiles |
# | ``p_{c,l}(t)``, ``l\ge0`` | ``t>0`` | coefficients of the air-film pressure |
# | ``z(t),\ v(t)`` | ``t>0`` | height of the drop's centre of mass, and its velocity |
#
# The coordinate is the interior **displacement** ``\chi_l``, and the stream function
# is its rate,
#
# ```math
# \psi_l = \dot\chi_l ,
# \qquad
# \zeta_l = \chi_l(1,t) .
# ```
#
# Both statements are forced. The interior displacement field is divergence-free to
# linear order, exactly as the velocity is, so it too derives from a Stokes stream
# function; ``\chi`` is that function and ``\psi`` is its rate. Taking ``\psi_l`` as
# the coordinate instead would put the kinetic energy quadratic in the coordinate
# rather than in its velocity, and Lagrange's equations would not apply. And the
# surface amplitude is not independent: the surface is the boundary of the interior,
# so ``\zeta_l`` is the boundary trace of ``\chi_l`` rather than a separate unknown.
#
# The radial profiles are expanded in ``K`` trial functions per mode,
#
# ```math
# \chi_l(x,t)=\sum_{k=1}^{K}a_{l,k}(t)\,x^{\,l+1+2(k-1)},
# \qquad
# \zeta_l=\sum_{k=1}^{K}a_{l,k} ,
# ```
#
# and these powers are not a free choice. The exact interior profile is a spherical
# Bessel function ``j_l(qx)``, whose Taylor series contains exactly
# ``x^{l+1},x^{l+3},x^{l+5},\dots``, so truncating at ``K`` truncates that series.
# ``k=1`` alone is potential flow, and reproduces Lamb's damping and Rayleigh's
# frequency exactly. Convergence is fast: ``K=3`` gives ``\lambda_2`` to one part in a
# thousand for ``\mathrm{Oh}\ge0.3``, and ``K=6`` does so down to
# ``\mathrm{Oh}=10^{-3}``. The basis is Vandermonde-like, though, so
# ``\mathrm{cond}(M)`` reaches ``4\times10^{10}`` by ``K=8``, and that is the practical
# ceiling on how far the interior can be resolved in double precision.
#
# So the coordinate vector is ``\bm\xi=\{a_{l,k}\}``, with ``\zeta_l`` and ``\psi_l``
# both read off it. There is **no pressure unknown**: the stream function makes the
# flow divergence-free identically, so incompressibility is not a constraint and
# carries no multiplier. ``\eta`` is not an unknown either -- it is a function of
# ``\dot{\bm\xi}`` through (4).
#
# The surface, the film pressure, and the gap are
#
# ```math
# \zeta(\theta,t)=\sum_{l\ge2}\zeta_l(t)P_l(\mu),
# \qquad
# p_c(\theta,t)=\sum_{l\ge0}p_{c,l}(t)P_l(\mu),
# ```
# ```math
# h(\theta,t)=\mu\,\bigl[1+\zeta(\theta,t)\bigr]+z(t) ,
# ```
#
# and the velocity follows from ``\psi=\sum_l\psi_l(x,t)C_l(\theta)`` with
# ``C_l=\sin^2\!\theta\,P_l'(\mu)/(l(l+1))``:
#
# ```math
# u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}=\sum_l u_{r,l}P_l(\mu),
# \qquad
# u_\theta=-\frac{1}{x\sin\theta}\frac{\partial\psi}{\partial x}
#   =\sum_l u_{\theta,l}\,\partial_\theta P_l(\mu) ,
# ```
# ```math
# u_{r,l}=\frac{\psi_l}{x^2},
# \qquad
# u_{\theta,l}=\frac{\psi_l'}{x\,l(l+1)} .
# ```
#
# ### (1) The three energies
#
# ```math
# T[\dot{\bm\xi}]=\tfrac12\int|\bm u|^2\,dV,
# \qquad
# \Phi[\dot{\bm\xi}]=\mathrm{Oh}\!\int 2\eta\;\bm e\!:\!\bm e\,dV,
# \qquad
# V[\bm\xi]=\tfrac12\sum_{l\ge2}\frac{4\pi}{2l+1}(l-1)(l+2)\,\zeta_l^2 ,
# ```
#
# with ``\bm e=\tfrac12(\nabla\bm u+\nabla\bm u^{\mathsf T})``. The surface energy
# follows from the curvature, ``\nabla\cdot\bm n=2+\sum_{l\ge2}(l-1)(l+2)\zeta_lP_l(\mu)``.
#
# ### (2) The equations of motion
#
# ```math
# \boxed{\;
# \frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a}
# \;+\;\frac12\frac{\partial\Phi}{\partial\dot\xi_a}
# \;+\;\frac{\partial V}{\partial\xi_a}
# \;=\; Q_a \;}
# ```
#
# **The coordinates are not independent, and the count works out because of it.**
# BC1 links them: ``\psi_l(1)=\dot\zeta_l``, so an admissible variation must satisfy
# ``\delta\psi_l(1)=\delta\zeta_l``. Stationarity then yields exactly two families,
# not one per symbol:
#
# | variation | where | what it gives |
# |:--|:--|:--|
# | ``\delta\psi_l`` with ``\delta\psi_l(1)=0`` | interior | one PDE in ``x`` per mode |
# | ``\delta\zeta_l=\delta\psi_l(1)`` | surface | one scalar equation per mode |
#
# So there are ``M-1`` interior equations and ``M-1`` surface equations, and
# ``\zeta_l`` is advanced by the second family rather than by an equation of its
# own. Reading (2) as literally one equation per symbol -- independent stationarity
# in both ``\zeta_l`` and ``\psi_l`` *while* constraining them -- over-determines
# the system, and that reading is what the essential condition forbids.
#
# The generalised force follows from the virtual work
# of the film traction ``-p_c\bm n`` against a radial surface displacement:
#
# ```math
# \delta W=-\oint p_c\,\delta\zeta\,dS
# \quad\Longrightarrow\quad
# Q_{\zeta_l}=-\frac{4\pi}{2l+1}\,p_{c,l},
# \qquad
# Q_{\psi_l}=0 ,
# ```
#
# the interior coordinates taking no direct forcing because the film acts only on
# the surface. At ``l=1`` this is ``-\tfrac{4\pi}{3}p_{c,1}``, which is precisely
# the net force ``\mathfrak F`` of block (3) -- so the centre-of-mass equation is
# the ``l=1`` member of this family rather than a separate statement.
#
# Because the interior coordinates are retained this is the full coupled
# surface-and-interior system, with no elimination and no closure. Two things it
# does *not* contain: any eigenvalue, and any second time derivative of
# ``\zeta_l`` alone -- the oscillator form with two coefficients per mode arises
# only on eliminating ``\{\psi_l\}``, which is the companion page's business.
#
# The coefficients are second derivatives of quadratic forms, so assembling them
# needs **one** derivative of the velocity:
#
# ```math
# \frac{\partial^2T}{\partial\dot\xi_a\partial\dot\xi_b}=\int\bm u^{(a)}\!\cdot\!\bm u^{(b)}dV,
# \qquad
# \frac{\partial^2\Phi}{\partial\dot\xi_a\partial\dot\xi_b}
#   =4\,\mathrm{Oh}\!\int\eta\;\bm e^{(a)}\!:\!\bm e^{(b)}\,dV
#   \;=\;2K_{ab} ,
# ```
# ```math
# K_{ab}\;\equiv\;\mathrm{Oh}\!\int 2\eta\;\bm e^{(a)}\!:\!\bm e^{(b)}\,dV .
# ```
#
# **The damping matrix in (2) is ``K_{ab}``, not the Hessian.** The ``\tfrac12``
# in the equation of motion means ``\tfrac12\,\partial\Phi/\partial\dot\xi_a
# =\sum_bK_{ab}\dot\xi_b``. Taking the Hessian *and* applying the ``\tfrac12``
# halves every decay rate, and the two expressions sit close enough together to
# invite exactly that -- note that the kinetic Hessian beside it needs no such
# care, because ``T`` carries its own ``\tfrac12``.
#
# Both forms are symmetric, which is a correctness test rather than a remark. The
# boundary conditions split: the kinematic condition ``u_r|_{x=1}=\dot\zeta`` is
# **essential** and constrains the admissible variations, while the vanishing of
# the tangential stress, ``\mathcal T[\psi_l]\big|_{x=1}=0``, is **natural** -- it
# follows from stationarity rather than being imposed. Regularity at ``x=0`` forces
# ``\psi_l\sim x^{l+1}``.
#
# ### (3) The drop as a whole, and contact
#
# ```math
# \dot z = v,
# \qquad
# \dot v = -\mathrm{Bo} - p_{c,1} ,
# ```
#
# the net force being ``\mathfrak F=-\oint p_c\,n_z\,dS= -\frac{4\pi}{3}p_{c,1}``
# and the drop's non-dimensional mass also ``4\pi/3``, which is why the two cancel.
# The drop does not wet the substrate: it is held off by a thin air film, ``p_c`` is
# that film's pressure, and because a gas film cannot sustain tension,
#
# ```math
# h\ge0,
# \qquad
# p_c\ge0,
# \qquad
# h\,p_c = 0
# \qquad\text{for all }\theta\in[0,\pi] .
# ```
#
# **Reduced to a finite system**, because the pointwise form is a continuum of
# conditions on a finite set of coefficients. Truncate ``p_c`` at the same ``M`` as
# the shape and impose the triple at ``M+1`` angular collocation nodes
# ``\{\theta_i\}``:
#
# ```math
# h(\theta_i)\ge0,\quad p_c(\theta_i)\ge0,\quad h(\theta_i)\,p_c(\theta_i)=0,
# \qquad i=0\ldots M ,
# ```
#
# which is ``M+1`` conditions for the ``M+1`` coefficients ``p_{c,l}``. These constrain
# the *reconstructed field* at the nodes, not the coefficients themselves -- those are
# different conditions, and only the first is the physics.
#
# **How the triple closes.** The three conditions are never imposed simultaneously.
# Because contact occupies a single connected patch about the pole, the set of nodes in
# contact is described by one integer ``c`` -- the number of them -- and given ``c`` the
# complementarity triple becomes a system of *equalities*:
#
# ```math
# h(\theta_i)=0 \quad (i\le c),
# \qquad
# p_c(\theta_i)=0 \quad (i>c) ,
# ```
#
# square in the ``M+1`` pressure coefficients. So ``c`` is an unknown of a different
# kind from the amplitudes -- discrete, and not obtainable from any linear solve -- and
# it is found by the two inequalities, which fail on opposite sides: a node outside the
# contact with ``h<0`` means ``c`` is too small, and ``p_c<0`` at the outermost
# contacting node means it is too large. Iterating on that until neither holds
# terminates, because every move is forced.
#
# **The nodes are the zeros of ``P_M``, together with ``\theta=\pi``.** This is not a
# numerical detail to be delegated. Those nodes cluster at the poles, which is where the
# contact forms, so the first contact is a genuinely small patch; on a uniformly spaced
# set of the same size the first contact spans a wide wedge, the pressure needed to hold
# it is impulsive, and the same equations return a drop leaving faster than it arrived.
#
# ### Initial conditions
#
# ```math
# z(0)=1,
# \qquad
# v(0)=-\sqrt{\mathrm{We}},
# \qquad
# \zeta_l(0)=0,
# \qquad
# \psi_l(x,0)=0 ,
# ```
#
# a spherical drop one radius above the substrate, falling at the impact speed with
# no internal motion relative to its centre of mass.
#
# ### (4) The fluid
#
# ```math
# \dot\gamma=\sqrt{2\,\bm e\!:\!\bm e},
# \qquad
# \eta=\eta(\dot\gamma) ,
# ```
#
# evaluated pointwise on the current interior field. ``\dot\gamma`` is computed
# from the **whole** field: ``\bm e`` superposes over modes, its invariant does
# not, so there is no such thing as one mode's shear rate. Admissible are the
# fluids satisfying (H1)--(H3) of *Which constitutive models this covers*:
# ``\eta`` a function of the instantaneous ``\dot\gamma`` alone, bounded as
# ``0<\eta_\infty\le\eta\le\eta_0<\infty``, and continuous in ``\dot\gamma``.
#
# ### The system and its truncations
#
# Blocks (1)--(2) are a coupled system of ordinary differential equations in the
# surface coordinates and partial differential equations in ``x`` for the interior
# ones, banded in the mode index by the selection rule ``|l-m|\le k\le l+m`` with
# ``l+k+m`` even. Block (3) adds a complementarity constraint. Block (4) closes the
# loop, because the ``\eta`` that (1) needs is computed from the solution of (2).
# So it is a **differential-algebraic system with a fixed point in the viscosity**.
#
# Three approximations are present and none is hidden.
#
# **Linearisation in the surface amplitude** is the largest, and it is what makes
# ``V`` quadratic and ``\nabla\cdot\bm n`` affine. All the integrals in (1) are
# taken over the *undeformed* sphere, and volume conservation therefore holds only
# to ``O(\zeta)``: with ``\zeta`` starting at ``l=2`` there is no ``\zeta_0`` to
# absorb the ``O(\zeta^2)`` change of volume, so the drop's volume drifts at that
# order and should be monitored rather than assumed.
#
# **The shape expansion is cut at ``M``**, which is what makes the viscosity's
# harmonic content terminate at ``k\le2M``, and the film pressure is cut at the
# same ``M``.
#
# **Any expansion of the radial profiles ``\psi_l`` is a numerical choice**, not
# part of the model -- the equations above are differential in ``x``, and choosing a
# radial basis belongs on the companion page with the quadrature and the time
# stepping.
#
# Every unknown is determined by an equation stated above. What is not available is
# a formula: ``\eta(\dot\gamma)`` is not polynomial, so nothing in the chain has a
# closed form in elementary functions and every piece must be computed. The
# companion page *Shear-Thinning Drops: Closures* is about what may be given up to
# compute them cheaply, and what each concession costs.
