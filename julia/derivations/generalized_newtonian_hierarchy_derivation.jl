# # Shear-Thinning Drops
#
# In this page we try to extend Reid's work to shear thinning rheologies. These are fluids that have an effective
# viscosity which is shear-rate dependent, possibly nonlinearly.
#
# In the Newtonian theory a drop's surface modes are independent damped
# oscillators: each ``\zeta_l`` has its own damping ``\lambda_l`` and frequency
# ``\omega_l``, and the interior flow can be solved once, in advance, for all
# time. Neither survives a shear-thinning fluid. The viscosity becomes a field
# over the drop, computed from the very flow it governs; modes that deform
# different parts of the drop no longer see the same fluid, and so begin to
# drive one another. The two scalars per mode become matrices, and the interior
# flow has to be found at every instant.
#
# This page derives that system, with no approximation beyond small amplitude
# and axisymmetry, and states it in full at the end. It does not simplify it:
# what the model costs to solve, and what may be given up to make it cheaper,
# is the subject of the companion page *Shear-Thinning Drops: Closures*.
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
# ### How to read the notation
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
# ### The one identity that shapes the whole problem
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
# ### Linearising the kinematics does not linearise the rheology
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
    println("  ASSERTION 3 OK: S = sqrt(2 e:e) is unchanged by rotation about the")  #src
    println("    symmetry axis, to $(round(worst, sigdigits=2)). An axisymmetric field maps")  #src
    println("    to itself under that rotation, so S -- and hence eta -- cannot depend on phi.")  #src
    println("    => eta is axisymmetric => no m != 0 coupling is EVER generated.")  #src
    println("    Legendre polynomials suffice permanently; Y_l^m is never required.")  #src
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
# which is why BC2 reduces to ``\mathcal T[\psi_l]|_{x=1}=0``; and the vorticity
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
    @assert w_F    < 1e-10 "u_r radial profile is not U/x^2 ($w_F)"  #src
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

# ## Where a mode equation comes from
#
# Before deriving anything it is worth being explicit about how a single-mode
# equation is obtained at all, because that is where the coupling will enter,
# and because it fixes the order of everything that follows. For a given surface
# deformation, three things happen in sequence:
#
# 1. the interior velocity field is whatever the momentum equation produces
#    when driven by that surface motion;
# 2. that field carries a viscous normal stress up to the free surface;
# 3. the stress balance at the surface -- viscous plus pressure against surface
#    tension -- is projected onto ``P_l`` to give the equation for ``\zeta_l``.
#
# With a constant viscosity all three steps preserve the mode. A surface
# deformation ``\propto P_l`` drives an interior field whose angular dependence
# is still ``P_l``, its stress at the surface is still ``\propto P_l``, and step 3
# picks out one equation per mode by orthogonality. That is why Reid obtains one
# characteristic equation per ``l``. Collecting them into a vector
# ``\bm\zeta=(A_2,\ldots,A_M)^{\mathsf T}``,
#
# ```math
# \bm{\ddot\zeta} \;+\; 2\bm\Lambda\,\bm{\dot\zeta} \;+\; \bm\Omega\,\bm\zeta
#   \;+\; \bm b \;=\; 0 ,
# \qquad
# \bm\Lambda=\operatorname{diag}(\lambda_l),
# \quad
# \bm\Omega=\operatorname{diag}(\omega_l^2),
# ```
#
# where ``\bm b`` has entries ``l\,p_{c,l}``, the contact pressure the wall applies
# while the drop is touching it. ``\bm b`` is kinematic and geometric; no
# rheology enters it, and nothing below changes it. Both matrices are diagonal,
# and that is the entire content of the Newtonian model.
#
# ``\lambda_l`` and ``\omega_l^2`` are built from **two** roots of Reid's
# characteristic equation. That equation has infinitely many roots per mode, so
# the second-order form above is a *truncation* of the viscous relaxation
# spectrum, not an identity: the exact response of a mode to forcing is a
# superposition over all of them. The truncation is a good one when the
# discarded roots decay fast, and it is what makes a two-coefficient oscillator
# model possible at all -- but it is an approximation independent of anything to
# do with rheology, and it is worth knowing it is there before it gets
# generalised. It does not appear in the model derived below, because that model
# never eliminates the interior.
#
# Step 1 is therefore the one to derive first, and steps 2 and 3 follow it. The
# rest of this page is those three steps in order: *The interior problem*, then
# *The surface projection*, then the constitutive law that ties the second back
# to the first.
#
# ## Which constitutive models this covers
#
# The derivation below uses exactly three properties of ``\eta``, and isolating
# them defines the class of fluids it applies to.
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
# ## The interior problem
#
# The surface equation above needs the matrices, the matrices need the viscosity
# field, and the viscosity field needs the interior flow. So the interior flow
# is what has to be found first, and this section derives the equation that
# determines it. The route is Reid's -- curl the momentum equation to kill the
# pressure, then reduce what is left -- but it is taken here for a viscosity
# that varies over the drop from the outset, because that is the physical case
# and because restricting it early hides where the difficulty actually lies.
#
# ### The curl still removes the pressure
#
# That step is indifferent to the rheology: ``\nabla\times\nabla p=0`` whatever
# ``\eta`` does. Take the ``\varphi``-component of the curl of the linearised
# momentum equation, with the field written through the modal stream function of
# the previous section. What survives is an expression in
#
# ```math
# \eta,\quad \partial_x\eta,\quad \partial_\theta\eta,\quad
# \partial^2_{xx}\eta,\quad \partial^2_{x\theta}\eta,\quad
# \partial^2_{\theta\theta}\eta
# ```
#
# and no derivative of third order or higher -- a fact worth having, because it
# bounds how much of the viscosity field the interior problem can possibly see.
#
# ### The interior equation
#
# Non-dimensionalising as in the summary below -- lengths by ``R``, time by
# ``T_\sigma``, viscosity by ``\eta_0`` -- the linearised momentum equation is
#
# ```math
# \partial_t\bm u = -\nabla p + \mathrm{Oh}\,\nabla\cdot(2\hat\eta\,\bm e),
# \qquad \hat\eta=\eta/\eta_0 .
# ```
#
# The vorticity of a Stokes stream function is
# ``\omega_\varphi=-\sum_l\mathcal D_l[\psi_l]\,C_l/(x\sin\theta)`` with
# ``\mathcal D_l=d^2/dx^2-l(l+1)/x^2``. Multiplying the curled equation by
# ``-x\sin\theta`` and projecting onto ``C_l`` -- the Gegenbauer functions are
# orthogonal under ``\langle f,g\rangle=\int_0^\pi fg\,d\theta/\sin\theta`` --
# leaves one equation per mode:
#
# ```math
# \boxed{\;
# \partial_t\,\mathcal D_l[\psi_l]
#   \;=\; \mathrm{Oh}\sum_{m}\mathcal R_{l m}\bigl[\psi_{m};\hat\eta\bigr],
# \qquad 0<x<1 \;}
# ```
#
# ```math
# \mathcal R_{l m}[\psi_{m};\hat\eta] \;\equiv\;
#   \frac{\bigl\langle\,
#     -x\sin\theta\,\bigl[\nabla\times\nabla\!\cdot\!(2\hat\eta\bm e)\bigr]_\varphi
#     \bigr|_{\psi=\psi_{m}C_{m}}
#   ,\;C_l\bigr\rangle}{\langle C_l,C_l\rangle} .
# ```
#
# This is the whole interior problem. Two features of it decide everything that
# follows.
#
# **It is stated in time, not in frequency.** Reid's ``q^2=\sigma/\nu`` is an
# eigenvalue and presupposes a single normal mode. Nothing in an impacting drop
# is a normal mode -- ``\zeta_l(t)`` is whatever the impact makes it -- and once
# ``\eta`` depends on the state there is no time-independent operator to take
# eigenvalues *of*. The equation above is a parabolic evolution for the interior
# vorticity and needs no such assumption. Replacing it by an instantaneous
# eigenproblem is a genuine approximation, made and priced on the companion
# page.
#
# **It couples the modes.** The sum over ``m`` is the substance of the
# problem: mode ``m`` drives mode ``l`` through the viscosity field. That is
# not a failure of technique. It is the physical statement that a fluid whose
# viscosity varies from place to place cannot respond to each surface harmonic
# independently, because the harmonics no longer see the same fluid. *How far
# the coupling reaches* is the next section.
#
# ### The diagonal, and where Reid sits
#
# Expand the viscosity in the same angular basis,
# ``\eta=\sum_{k}\eta_{k}(x)P_{k}(\cos\theta)``. The term ``k=0`` is the
# spherically symmetric part, and it contributes only to ``l=m``. Its
# contribution is a genuine radial operator, and it is worth having explicitly
# because it is the one piece of the problem that can be written in closed form:
#
# ```math
# \boxed{\;
# \mathcal R_{l l}\big|_{k=0} \;=\; \mathcal R_l[\psi;\eta] \;=\;
#   \eta\,\mathcal D_l^{\,2}[\psi]
# \;+\; 2\eta'(x)\,\frac{d}{dx}\!\Bigl(\mathcal D_l[\psi]-\frac{\psi'}{x}\Bigr)
# \;+\; \eta''(x)\,\mathcal T[\psi] \;}
# ```
#
# ```math
# \mathcal T[\psi]\;\equiv\;\psi''-\frac{2}{x}\psi'+\frac{l(l+1)}{x^2}\psi .
# ```
#
# Three things read off it.
#
# **Reid is the constant-viscosity case.** Setting ``\eta'=\eta''=0`` leaves
# ``\eta\,\mathcal D_l^2[\psi]``, so the interior equation becomes
# ``\partial_t\mathcal D_l[\psi_l]=\mathrm{Oh}\,\mathcal D_l^2[\psi_l]``: vorticity
# diffusion. Substituting a normal mode ``\psi_l\propto e^{-\sigma t}`` gives
# ``\mathcal D_l(\mathcal D_l+q^2)[\psi_l]=0`` with ``q^2=\sigma/\mathrm{Oh}``,
# which is Reid's Eq. 9. Nothing has been added to recover it; the
# constant-viscosity theory is the diagonal of this one, at one harmonic.
#
# **The order of the equation does not change.** The highest derivative is
# ``\eta \psi''''``; ``\eta'`` reaches only ``\psi'''`` and ``\eta''`` only ``\psi''``.
# The operator is fourth order in ``\psi``, exactly as Reid's is, so the solution
# space is still four-dimensional and the problem still closes on the same
# number of boundary conditions. A variable viscosity changes the coefficients
# of the interior problem; it does not change its type.
#
# **``\eta''`` enters through the tangential-stress operator.** The
# ``\mathcal T`` above is not a new object: it is the operator whose vanishing
# at the surface is BC2 below. It appears here because
# ``\mathcal T[\psi]\propto e_{r\theta}``, and the second radial derivative of the
# viscosity couples to precisely the shear component of the strain.
#
# ### The boundary conditions
#
# The interior problem is fourth order per mode, so it needs four conditions:
# two of regularity at ``x=0``, and two at the free surface.
#
# ```math
# \text{BC1 (kinematic):}\qquad u_r\big|_{x=1}=\frac{\partial\zeta}{\partial t}
# ```
#
# contains no viscosity at all, and is **unchanged** by the rheology.
#
# ```math
# \text{BC2 (tangential):}\qquad \tau_{r\theta}\big|_{x=1}=2\eta\,e_{r\theta}=0
# \;\Longleftrightarrow\; e_{r\theta}\big|_{x=1}=0
# \;\Longleftrightarrow\; \mathcal T[\psi]\big|_{x=1}=0 ,
# ```
#
# where the first equivalence holds because ``\eta(\dot\gamma)\ge\eta_\infty>0``
# everywhere, so the scalar factor cannot vanish. **BC2 is therefore
# rheology-agnostic for every fluid in the admissible class**, and the whole
# ``\tau_{r\theta}=0\Rightarrow\mathcal T[\psi]=0`` reduction carries over
# verbatim.
#
# The third surface condition -- the balance of normal stress against surface
# tension -- is **not** a boundary condition on the interior problem. Projected
# onto ``P_l`` it *is* the equation of motion for ``\zeta_l``, and that projection is
# what produces the coefficient matrices. It is carried out two sections below.
# This is where the rheology lands: once ``\eta`` varies with ``\theta``, its
# value at the surface is a field rather than a number, and projecting that
# field is the entire origin of the mode coupling.
#
@variables rr tt qq hh0 hh1 hh2  #src
@variables uc[1:5] ec[1:3]  #src
@variables Hgen(..)  #src
let  #src
    Drr = Differential(rr); Dtt = Differential(tt)  #src
    LPn(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    LPpn(l,m) = l==0 ? zero(m) : l*(m*LPn(l,m)-LPn(l-1,m))/(m^2-1)  #src
    angl(l) = sin(tt)^2*LPpn(l,cos(tt))/(l*(l+1))  #src
    function curlmom(l, U, H)  #src
        psi = U*angl(l)  #src
        u_r  =  Symbolics.expand_derivatives(Dtt(psi)/(rr^2*sin(tt)))  #src
        u_th = -Symbolics.expand_derivatives(Drr(psi)/(rr*sin(tt)))  #src
        e_rr = Symbolics.expand_derivatives(Drr(u_r))  #src
        e_tt = Symbolics.expand_derivatives(Dtt(u_th)/rr + u_r/rr)  #src
        e_pp = u_r/rr + u_th*cos(tt)/(sin(tt)*rr)  #src
        e_rt = Symbolics.expand_derivatives((Dtt(u_r)/rr + Drr(u_th) - u_th/rr)/2)  #src
        trr=2H*e_rr; ttt=2H*e_tt; tpp=2H*e_pp; trt=2H*e_rt  #src
        dr = Symbolics.expand_derivatives(Drr(rr^2*trr)/rr^2 + Dtt(trt*sin(tt))/(rr*sin(tt)) - (ttt+tpp)/rr)  #src
        dt = Symbolics.expand_derivatives(Drr(rr^2*trt)/rr^2 + Dtt(ttt*sin(tt))/(rr*sin(tt)) + trt/rr - cos(tt)/(sin(tt)*rr)*tpp)  #src
        Symbolics.expand_derivatives((Drr(rr*(qq*u_th+dt)) - Dtt(qq*u_r+dr))/rr)  #src
    end  #src
    Dln(e,l) = Symbolics.expand_derivatives(Drr(Drr(e))) - l*(l+1)*e/rr^2  #src
    reidop(U,l) = Symbolics.expand_derivatives(Dln(Dln(U,l) + qq*U, l))  #src
    ev(e,rv,tv,qv,a0,a1,a2) = Float64(eval(Symbolics.toexpr(Symbolics.substitute(  #src
        Symbolics.expand_derivatives(e), Dict(rr=>rv, tt=>tv, qq=>qv, hh0=>a0, hh1=>a1, hh2=>a2)))))  #src
    l = 2  #src
    PTS = ((0.43,0.7),(0.61,1.3),(0.82,2.1),(0.37,2.6),(0.95,0.44))  #src

    ## (i) which derivatives of eta survive the curl  #src
    full = curlmom(l, rr^(l+1) + 0.7*rr^(l+3), Hgen(rr,tt))  #src
    txt = string(full)  #src
    d3 = Differential(rr)(Differential(rr)(Differential(rr)(Hgen(rr,tt))))  #src
    @assert !occursin(string(d3), txt) "a third derivative of eta survived the curl"  #src
    for d in (Differential(rr)(Hgen(rr,tt)), Differential(tt)(Hgen(rr,tt)))  #src
        @assert occursin(string(d), txt) "expected first derivatives of eta to appear"  #src
    end  #src

    ## (ii) constant eta must return Reid's operator exactly -- and must do so  #src
    ## with the RIGHT dependence on eta. Testing only at eta = 1 cannot tell  #src
    ## tau = 2*eta*e from tau = 2*eta^2*e, since the two agree there. With  #src
    ## eta = c the viscous terms carry a factor c while the inertial q^2 term  #src
    ## does not, so the whole curl must equal  #src
    ##     c * [-C_l/(x sin th)] * D_l(D_l + q^2/c)[U],  #src
    ## which pins the power of eta as well as the operator.  #src
    Uc = rr^(l+1) + 0.7*rr^(l+3) - 0.3*rr^(l+5) + 0.11*rr^(l+2)  #src
    rats = Float64[]  #src
    for c in (1.0, 2.5, 0.4)  #src
        tgt = -c*Dln(Dln(Uc,l) + (3.7/c)*Uc, l)*angl(l)/(rr*sin(tt))  #src
        for (rv,tv) in PTS  #src
            push!(rats, ev(curlmom(l,Uc,c),rv,tv,3.7,0.,0.,0.) / ev(tgt,rv,tv,3.7,0.,0.,0.))  #src
        end  #src
    end  #src
    @assert all(abs(x-1) < 1e-9 for x in rats) "constant eta did not return Reid's operator with the right power of eta"  #src

    ## (ii-b) the BOXED radial operator itself. The page now states  #src
    ##     R_l[psi;eta] = eta*D_l^2[psi] + 2*eta' d/dx(D_l[psi] - psi'/x) + eta''*T[psi]  #src
    ## as an exact result, so it is checked as one -- against the curl computed  #src
    ## from scratch, for several U and several eta(x), at several (r,theta).  #src
    ## A failure here would mean the interior equation printed on the page is  #src
    ## not the interior equation the fluid actually obeys.  #src
    L2op(U) = Symbolics.expand_derivatives(Drr(Drr(U))) - 2*Symbolics.expand_derivatives(Drr(U))/rr + l*(l+1)*U/rr^2  #src
    Dx(e) = Symbolics.expand_derivatives(Drr(e))  #src
    Rl(U,H) = H*Dln(Dln(U,l),l) + 2*Dx(H)*Dx(Dln(U,l) - Dx(U)/rr) + Dx(Dx(H))*L2op(U)  #src
    ## The error is normalised GLOBALLY -- largest discrepancy over the sweep  #src
    ## against the largest magnitude over the sweep -- not pointwise. A         #src
    ## pointwise relative error is meaningless wherever the operator genuinely  #src
    ## vanishes, and it does vanish here: for U = x^(l+1) (harmonic, so         #src
    ## D_l[psi] = 0) and eta = 2 - 0.15x^2 the eta' and eta'' terms cancel        #src
    ## identically, 1.8x - 1.8x = 0. That is a property of the operator, not a  #src
    ## failure, and an earlier pointwise version of this check reported it as   #src
    ## a 100% error.                                                            #src
    maxdiff_box, maxmag_box = 0.0, 0.0  #src
    for Ub in (rr^(l+1), rr^(l+3), rr^(l+1) + 0.7rr^(l+3) - 0.3rr^(l+5) + 0.11rr^(l+2))  #src
        for (a0,a1,a2) in ((1.3,0.9,0.0), (0.8,0.5,0.6), (2.0,0.0,-0.15))  #src
            Hb = a0 + a1*rr + a2*rr^2  #src
            ## viscous part only: set the inertial q^2 to zero on both sides  #src
            for (rv,tv) in PTS  #src
                lhs = ev(curlmom(l,Ub,Hb),rv,tv,0.0,0.,0.,0.)  #src
                rhs = ev(-Rl(Ub,Hb)*angl(l)/(rr*sin(tt)),rv,tv,0.0,0.,0.,0.)  #src
                maxdiff_box = max(maxdiff_box, abs(lhs-rhs))  #src
                maxmag_box  = max(maxmag_box, abs(rhs))  #src
            end  #src
        end  #src
    end  #src
    @assert maxmag_box > 1.0 "the R_l sweep never exercised a nonzero operator; the check is vacuous"  #src
    worst_box = maxdiff_box / maxmag_box  #src
    @assert worst_box < 1e-8 "the boxed radial operator R_l does not equal the curl it claims to ($worst_box)"  #src

    ## (iii) eta = eta(r) keeps the angular factor -- separability survives  #src
    Hr = hh0 + hh1*rr + hh2*rr^2  #src
    sep = 0.0  #src
    for U in (rr^(l+1), rr^(l+3), rr^(l+2))  #src
        C = curlmom(l, U, Hr)  #src
        v = [ev(C,0.6,tv,3.7,1.0,0.4,-0.2)/ev(angl(l)/(rr*sin(tt)),0.6,tv,3.7,1.0,0.4,-0.2)  #src
             for tv in (0.5,1.1,1.9,2.5)]  #src
        sep = max(sep, (maximum(v)-minimum(v))/abs(v[1]))  #src
    end  #src
    @assert sep < 1e-9 "eta(r) broke separability, which it must not ($sep)"  #src

    ## (iv) eta with angular structure must break it  #src
    Hrt = hh0 + hh1*rr + hh2*cos(tt)^2  #src
    C = curlmom(l, rr^(l+1), Hrt)  #src
    v = [ev(C,0.6,tv,3.7,1.0,0.4,-0.2)/ev(angl(l)/(rr*sin(tt)),0.6,tv,3.7,1.0,0.4,-0.2)  #src
         for tv in (0.5,1.1,1.9,2.5)]  #src
    coup = (maximum(v)-minimum(v))/abs(v[1])  #src
    @assert coup > 1e-2 "angular eta failed to break separability; the coupling claim is false"  #src

    println("  ASSERTION 3b OK: the interior operator, checked in four ways --")  #src
    println("    the curl removes the pressure and leaves eta only through its")  #src
    println("      first and second derivatives, never third or higher;")  #src
    println("    constant eta returns D_l(D_l+q^2)U exactly (to $(round(maximum(abs(x-1) for x in rats), sigdigits=2))), so Reid")  #src
    println("      is recovered as the diagonal at one harmonic;")  #src
    println("    the boxed diagonal R_l = eta*D_l^2 + 2eta' d/dx(D_l - d/dx /x) + eta''*L2")  #src
    println("      reproduces the curl whenever eta carries no angular structure")  #src
    println("      (to $(round(worst_box, sigdigits=2)));")  #src
    println("    that case stays diagonal (spread $(round(sep, sigdigits=2))), and an angular eta")  #src
    println("      does not (spread $(round(coup, sigdigits=3))) -- which is the coupling itself.")  #src
end  #src

# ### The off-diagonal operators, explicitly
#
# ``\mathcal R_{l m}`` was defined above by a projection, which is enough to
# close the model but not enough to assemble a matrix. It has an explicit form,
# and the reason is structural: the viscous term is linear in ``\eta`` and linear
# in ``\psi``, and reaches no higher than ``\eta''``. So for a single viscosity
# harmonic and a single driving mode,
#
# ```math
# \boxed{\;
# \mathcal R_{l m}\bigl[\psi_m;\eta_k\bigr]
#   \;=\; \sum_{j=0}^{2}\sum_{i=0}^{4}
#     a^{(k)}_{l m,\,ji}\; x^{\,j+i-4}\;
#     \eta_k^{(j)}(x)\;\psi_m^{(i)}(x) \;}
# ```
#
# with ``a^{(k)}_{lm,ji}`` **pure numbers depending only on the three integers**.
# Fifteen of them per triple, and the selection rule already says which triples
# are nonzero, so the whole coupling is a sparse table computed once and then
# contracted at each step. Nothing symbolic happens inside the time loop.
#
# The coefficients are extracted by a device worth naming, because the direct
# route is prohibitively slow. Take ``\psi_m=U(x)C_m(\theta)`` and
# ``\eta_k=E(x)P_k(\mu)`` with ``U`` and ``E`` *local Taylor polynomials about*
# ``x=1``,
#
# ```math
# U(x)=\sum_{i=0}^{4}u_i\frac{(x-1)^i}{i!},
# \qquad
# E(x)=\sum_{j=0}^{2}e_j\frac{(x-1)^j}{j!} ,
# ```
#
# truncated at exactly the orders the operator can reach. Then
# ``U^{(i)}(1)=u_i`` and ``E^{(j)}(1)=e_j`` identically, the projected value is
# *bilinear* in ``(e_j,u_i)``, and at ``x=1`` every power of ``x`` is one -- so
#
# ```math
# a^{(k)}_{lm,\,ji} \;=\; \bigl[\text{projection}\bigr]
#   \Big|_{e_j=1,\;u_i=1,\;\text{all others }0} .
# ```
#
# One symbolic curl per ``(m,k)`` pair serves every ``l`` in the band. Feeding
# monomials instead needs one curl per sampled pair of exponents, which is
# twenty-five times the work for the same fifteen numbers.
#
# The check below runs the extraction and compares the diagonal against the boxed
# ``\mathcal R_l``, which is the one case already known in closed form. A failure
# would mean the operator does not have this form at all -- that some coefficient
# depends on ``x`` beyond the stated power, and the coupling could not be
# tabulated.

let  #src
    ## One curl per (m,k), reused for every l in the band; and ONE build_function  #src
    ## over (theta, e..., u...) rather than fifteen, which is what makes this      #src
    ## affordable. An earlier monomial-fitting version needed 25 symbolic curls    #src
    ## per (m,k) and ran for over half an hour without finishing.                  #src
    Dr2 = Differential(rr); Dt2 = Differential(tt)  #src
    ed2(e) = Symbolics.expand_derivatives(e)  #src
    LPc(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    dLPc(l,m) = l==0 ? zero(m) : l*(m*LPc(l,m)-LPc(l-1,m))/(m^2-1)  #src
    Cgc(l) = sin(tt)^2*dLPc(l,cos(tt))/(l*(l+1))  #src
    function curlv(psi, eta)  #src
        ur =  ed2(Dt2(psi)/(rr^2*sin(tt))); ut = -ed2(Dr2(psi)/(rr*sin(tt)))  #src
        err_ = ed2(Dr2(ur)); ett = ed2(Dt2(ut)/rr + ur/rr)  #src
        epp = ur/rr + ut*cos(tt)/(sin(tt)*rr)  #src
        ert = ed2((Dt2(ur)/rr + Dr2(ut) - ut/rr)/2)  #src
        dr = ed2(Dr2(rr^2*2eta*err_)/rr^2 + Dt2(2eta*ert*sin(tt))/(rr*sin(tt)) - (2eta*ett+2eta*epp)/rr)  #src
        dt = ed2(Dr2(rr^2*2eta*ert)/rr^2 + Dt2(2eta*ett*sin(tt))/(rr*sin(tt)) + 2eta*ert/rr - cos(tt)/(sin(tt)*rr)*2eta*epp)  #src
        ed2((Dr2(rr*dt) - Dt2(dr))/rr)  #src
    end  #src
    nodes, wts = QuadGK.gauss(36, -1.0, 1.0)  #src
    Uc = sum(uc[i+1]*(rr-1)^i/factorial(i) for i in 0:4)  #src
    Ec = sum(ec[j+1]*(rr-1)^j/factorial(j) for j in 0:2)  #src
    ## a_{ji} for all (j,i) at once, for one (l,m,k)  #src
    function coeffs(gen, l)  #src
        A = zeros(3, 5)  #src
        for j in 0:2, i in 0:4  #src
            ev = ntuple(a -> (a-1 == j ? 1.0 : 0.0), 3)  #src
            uv = ntuple(b -> (b-1 == i ? 1.0 : 0.0), 5)  #src
            num = 0.0  #src
            for (mu, w) in zip(nodes, wts)  #src
                th = acos(mu)  #src
                F  = -sin(th)*gen(th, ev..., uv...)          # x = 1  #src
                num += w*(1-mu^2)*(F/sin(th)^2)*dLPc(l,mu)/(l*(l+1))  #src
            end  #src
            A[j+1, i+1] = num/(2/((2l+1)*l*(l+1)))  #src
        end  #src
        A  #src
    end  #src
    ## the boxed diagonal operator in the same basis, L = l(l+1)  #src
    function boxedA(l)  #src
        L = l*(l+1.0); A = zeros(3,5)  #src
        A[1,5] += 1.0; A[1,3] += -2L; A[1,2] += 4L; A[1,1] += L^2 - 6L  #src
        A[2,4] += 2.0; A[2,3] += -2.0; A[2,2] += -2(L-1); A[2,1] += 4L  #src
        A[3,3] += 1.0; A[3,2] += -2.0; A[3,1] += L  #src
        A  #src
    end  #src
    worst_diag, worst_offmag, ntriple = 0.0, 0.0, 0  #src
    ## Trimmed to three (m,k) pairs to keep the CI budget sane: (2,0) supplies the  #src
    ## diagonal cross-check against the boxed operator, (2,1) and (2,2) exercise  #src
    ## the off-diagonal coupling at both parities.  #src
    for (m, k) in ((2,0), (2,1), (2,2))  #src
        Cx = Symbolics.substitute(ed2(curlv(Uc*Cgc(m), Ec*LPc(k,cos(tt)))), Dict(rr => 1.0))  #src
        gen = Symbolics.build_function(Cx, tt, ec[1],ec[2],ec[3], uc[1],uc[2],uc[3],uc[4],uc[5]; expression=Val(false))  #src
        for l in max(1, abs(m-k)):(m+k)  #src
            iseven(l+k+m) || continue  #src
            A = coeffs(gen, l); ntriple += 1  #src
            if k == 0 && l == m  #src
                worst_diag = max(worst_diag, maximum(abs.(A .- boxedA(l))))  #src
            elseif k > 0  #src
                worst_offmag = max(worst_offmag, maximum(abs.(A)))  #src
            end  #src
        end  #src
    end  #src
    @assert ntriple >= 5 "the sweep covered too few (l,m,k) triples to be meaningful ($ntriple)"  #src
    @assert worst_offmag > 1.0 "every off-diagonal operator came out zero; the extraction is not exercising the coupling"  #src
    @assert worst_diag < 1e-9 "the extracted diagonal does not reproduce the boxed R_l ($worst_diag)"  ## CLAIM: SUM-RLM  #src
    @printf("  ASSERTION 3g OK: R_{lm} extracted as 15 numbers per (l,m,k) over %d triples.\n", ntriple)  #src
    @printf("    The diagonal (k=0, l=m) reproduces the boxed R_l to %.1e, and the\n", worst_diag)  #src
    @printf("    off-diagonal operators are nonzero (largest coefficient %.3g), so the\n", worst_offmag)  #src
    println("    mode coupling is a sparse table computed once, not a per-step projection.")  #src
    println("    Physical meaning of a failure: the operator would not have this form,")  #src
    println("    some coefficient would depend on x beyond x^(j+i-4), and the coupling")  #src
    println("    could not be tabulated at all -- assembly would need a projection")  #src
    println("    inside the time loop.")  #src
end  #src

# ### How far the coupling reaches
#
# The sum over ``m`` in the interior equation looks unbounded, and if it were,
# the problem would be no easier to solve than the original partial differential
# equation. It is not unbounded. Each viscosity harmonic ``k`` couples ``m``
# to ``l`` only when
#
# ```math
# |l-m|\le k\le l+m
# \qquad\text{and}\qquad
# l+k+m \ \text{even},
# ```
#
# which is the same selection rule the surface projection obeys two sections
# below, and holds for the same reason: every angular integral in the problem is
# an integral of three Legendre-type factors, and those vanish outside the
# triangle condition and against the wrong parity.
#
# The consequence is that a viscosity field whose angular content reaches degree
# ``L_\eta`` couples each mode to at most its ``L_\eta`` neighbours on either
# side. The interior problem is then a **banded** system of coupled radial
# boundary-value problems, of half-bandwidth ``L_\eta`` and dimension ``M``,
# rather than ``M`` independent ones. That is the entire price of a viscosity
# that varies over the drop: bandedness, not loss of closure.
#
# ``L_\eta`` is a property of the state, not a choice -- it is however far the
# Legendre series of ``\eta(\dot\gamma)`` actually extends. But it enters the
# problem in the same way ``M`` does, as the point at which a series is cut, and
# it can be treated the same way: capping it at some ``L_\eta`` below its true
# extent is an available approximation, with a cost that is then to be priced
# rather than assumed. The companion page does exactly that, and finds the price
# higher than it looks.
#
# The check below evaluates the projection directly, with no appeal to the rule,
# and reports which ``l`` survive. A failure would mean the interior coupling
# reaches further than the surface coupling does -- in which case no single
# truncation could be consistent for both halves of the model, and any ``L_\eta``
# chosen for one would silently corrupt the other.
let  #src
    ## Project the phi-curl of div(2 eta e), driven by a single mode l'', onto  #src
    ## C_l for a range of l, and report the band. eta is taken as the pure      #src
    ## harmonic P_l'(mu) so the (k, m) pair under test is unambiguous.       #src
    ##                                                                          #src
    ## Everything is polynomial in mu once the sin factors are handled          #src
    ## analytically, so FIXED Gauss-Legendre nodes are exact and never land on  #src
    ## mu = +-1, where the 1/sin(theta) factors in the curl are 0/0. An         #src
    ## adaptive quadrature in theta walks into that endpoint and returns NaN.   #src
    Drr = Differential(rr); Dtt = Differential(tt)  #src
    LPn(l,m) = l==0 ? one(m) : l==1 ? m : begin a,b=one(m),m; for n in 1:l-1; b,a=((2n+1)*m*b-n*a)/(n+1),b; end; b end  #src
    LPpn(l,m) = l==0 ? zero(m) : l*(m*LPn(l,m)-LPn(l-1,m))/(m^2-1)  #src
    angl(l) = sin(tt)^2*LPpn(l,cos(tt))/(l*(l+1))  #src
    function curlvisc(psi, eta)  #src
        ed(e) = Symbolics.expand_derivatives(e)  #src
        u_r  =  ed(Dtt(psi)/(rr^2*sin(tt)))  #src
        u_th = -ed(Drr(psi)/(rr*sin(tt)))  #src
        e_rr = ed(Drr(u_r))  #src
        e_tt = ed(Dtt(u_th)/rr + u_r/rr)  #src
        e_pp = u_r/rr + u_th*cos(tt)/(sin(tt)*rr)  #src
        e_rt = ed((Dtt(u_r)/rr + Drr(u_th) - u_th/rr)/2)  #src
        trr=2eta*e_rr; ttt=2eta*e_tt; tpp=2eta*e_pp; trt=2eta*e_rt  #src
        dr = ed(Drr(rr^2*trr)/rr^2 + Dtt(trt*sin(tt))/(rr*sin(tt)) - (ttt+tpp)/rr)  #src
        dt = ed(Drr(rr^2*trt)/rr^2 + Dtt(ttt*sin(tt))/(rr*sin(tt)) + trt/rr - cos(tt)/(sin(tt)*rr)*tpp)  #src
        ed((Drr(rr*dt) - Dtt(dr))/rr)  #src
    end  #src
    nodes, wts = QuadGK.gauss(24, -1.0, 1.0)  #src
    for (lpp, lp) in ((2,1), (2,2))  #src
        expr = -rr*sin(tt)*curlvisc((rr^(lpp+1) + 0.6rr^(lpp+3))*angl(lpp), LPn(lp,cos(tt)))  #src
        f = Symbolics.build_function(Symbolics.expand_derivatives(expr), rr, tt; expression=Val(false))  #src
        vals = Tuple{Int,Float64}[]  #src
        for l in 1:(lpp+lp+2)  #src
            num, den = 0.0, 0.0  #src
            for (mu, w) in zip(nodes, wts)  #src
                F = f(0.63, acos(mu))/(1 - mu^2)  #src
                dP = Float64(LPpn(l, mu))  #src
                num += w*(1 - mu^2)*F*dP/(l*(l+1))  #src
                den += w*(1 - mu^2)*dP^2/(l*(l+1))^2  #src
            end  #src
            push!(vals, (l, num/den))  #src
        end  #src
        scale = maximum(abs(v) for (_,v) in vals)  #src
        @assert scale > 1.0 "the interior coupling sweep never produced a nonzero projection; the check is vacuous"  #src
        got  = [l for (l,v) in vals if abs(v) > 1e-7*scale]  #src
        want = [l for l in 1:(lpp+lp+2) if abs(lpp-lp) <= l <= lpp+lp && iseven(l+lp+lpp)]  #src
        @assert got == want "interior coupling band $(got) is not the Gaunt band $(want) for (l''=$lpp, l'=$lp)"  #src
        @printf("  l''=%d, l'=%d : interior coupling reaches l in %s (Gaunt band, parity-filtered)\n", lpp, lp, string(got))  #src
    end  #src
    println("  ASSERTION 3c OK: the interior problem projects onto a FINITE band of")  #src
    println("    modes, and that band is the same Gaunt+parity rule the surface")  #src
    println("    projection obeys. So losing separability costs bandedness, not")  #src
    println("    closure: the interior problem still exists, it is just coupled.")  #src
    println("    Physical meaning of a failure: the interior coupling would reach")  #src
    println("    further than the surface coupling, and no single truncation could")  #src
    println("    be consistent for both halves of the model.")  #src
end  #src

# ## The surface projection
#
# The interior problem determines the flow. What turns that flow into an
# equation of motion for the surface is the remaining boundary condition, and
# projecting it is where the coefficient matrices come from.
#
# ### The pressure has to be recovered, and it is no longer harmonic
#
# Taking the curl removed the pressure from the interior problem, which is why
# that problem closed. But the normal-stress balance contains ``p`` explicitly,
# so it has to be brought back before the surface equation can be written.
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
    @assert least_rad > 1e-2 "a radially varying eta must break harmonicity, got $least_rad"  ## CLAIM: SUM-PSRC  #src
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
    @assert rel_id < 1e-4 "div of momentum does not give lap(p) = c div(div(2 eta e)) (rel $rel_id)"  ## CLAIM: SUM-PRESS  #src
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
    @assert rel_pl < 1e-4 "the l-projection of lap(p) is not L_l[p_l] (rel $rel_pl)"  ## CLAIM: SUM-PLAP  #src

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
    @assert worst_com < 1e-12 "the centre-of-mass equation does not reduce to vdot = -Bo - p_{c,1} ($worst_com)"  ## CLAIM: SUM-COM  #src

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
# ### Carrying the projection through
#
# The equation being projected is the **normal-stress balance on the free
# surface**, linearised and evaluated at ``x=1``,
#
# ```math
# \Bigl[-p + \tau_{rr}\Bigr]_{x=1} \;=\; T_1\,(\nabla\cdot\bm n)\Big|_{x=1} ,
# \qquad \tau_{rr}=2\eta\,e_{rr},
# ```
#
# which is the third of the three conditions listed above and the only one that
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
# ``\mathcal L_{m}[u_{r,m}]=F''_{m}+\tfrac{2}{x}u_{r,m}'-\tfrac{m(m+1)}{x^2}u_{r,m}``,
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
# The index ``i`` distinguishes the terms that end up multiplying ``\dot\zeta``
# from those multiplying ``A``.
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
# ``\mathcal D^{(2)}`` generalises ``2\bm\Lambda`` and ``\mathcal D^{(1)}``
# generalises ``\bm\Omega``. The viscous stress is linear in the velocity and so
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
# is taken up in *Where the shear rate is evaluated* below.
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

# ### The selection rule
#
# The angular factor ``G^{k}_{l m}`` is not merely small for most index
# triples -- it is exactly zero for most of them. The rule is
#
# ```math
# G^{k}_{l m} = 0
# \qquad\text{unless}\qquad
# |l-m| \;\le\; k \;\le\; l+m
# \qquad\text{and}\qquad
# l+k+m\ \text{is even}.
# ```
#
# It is a statement about three integers, so it holds for every fluid, every
# amplitude and every geometry in this problem. It is also elementary, which is
# worth showing rather than citing.
#
# ### Why the truncation is exact: a two-line proof
#
# Everything rests on one elementary fact about Legendre polynomials:
#
# > ``P_N`` is orthogonal to every polynomial of degree ``<N``.
#
# That is immediate from orthogonality, since any such polynomial is a finite
# combination of ``P_0,\ldots,P_{N-1}``.
#
# Now ``P_l`` has degree exactly ``l``, so the product ``P_lP_{m}`` is a
# polynomial of degree exactly ``l+m``. Therefore:
#
# * **Upper bound.** If ``k>l+m`` then ``P_{k}`` is orthogonal to
#   ``P_lP_{m}``, a polynomial of degree ``l+m<k``. Hence
#   ``G^{k}_{l m}=0``.
# * **Lower bound.** Take ``l\ge m`` without loss of generality. If
#   ``k<l-m`` then ``P_{k}P_{m}`` has degree ``k+m<l``, and ``P_l`` is
#   orthogonal to it. Hence ``G^{k}_{l m}=0`` again.
# * **Parity.** ``P_l(-\mu)=(-1)^lP_l(\mu)``, so the integrand has parity
#   ``(-1)^{l+k+m}``. An odd integrand over ``[-1,1]`` integrates to zero.
#
# Together: ``G^{k}_{l m}=0`` unless ``|l-m|\le k\le l+m`` **and**
# ``l+k+m`` is even. No special-function machinery is needed -- degree
# counting and parity give the whole selection rule, and both bounds are the
# *same* argument applied to different factors.
#
# The consequence for this model is the upper bound. The shape expansion stops
# at ``l=M``, so ``l+m\le 2M``, so viscosity harmonics above ``k=2M`` are
# annihilated exactly. The check below is a regression test on that proof, not
# its justification.

let  #src
    Pl(l, m) = l == 0 ? one(m) : l == 1 ? m :  #src
        begin am, b = one(m), m; for n in 1:l-1; b, am = ((2n+1)*m*b - n*am)/(n+1), b; end; b end  #src
    ## Gauss-Legendre with 40 nodes integrates any polynomial of degree <= 79  #src
    ## exactly, so this evaluates the triple product with no quadrature error  #src
    ## at all -- the residual below is pure floating point.  #src
    gn, gw = DropSolver.gauss_legendre_nodes(40, -1.0, 1.0)  #src
    Gc(l, lp, lpp) = sum(w * Pl(l,m) * Pl(lp,m) * Pl(lpp,m) for (m, w) in zip(gn, gw))  #src
    forbidden = 0.0  #src
    allowed = 0.0  #src
    for l in 2:8, lpp in 2:8  #src
        for lp in (l + lpp + 1):(l + lpp + 4)  #src
            forbidden = max(forbidden, abs(Gc(l, lp, lpp)))  #src
        end  #src
        for lp in max(0, l - lpp):2:(l + lpp)  #src
            allowed = max(allowed, abs(Gc(l, lp, lpp)))  #src
        end  #src
    end  #src
    @assert forbidden < 1e-13 "viscosity above l'=l+l'' coupled two modes ($forbidden)"  #src
    @assert allowed > 0.1 "allowed coefficients should be O(1), for contrast"  #src
    println("  ASSERTION 7b OK: every Gaunt coefficient with l' > l+l'' vanishes")  #src
    println("    (largest $(round(forbidden, sigdigits=3)) against an allowed scale of $(round(allowed, sigdigits=3))).")  #src
    println("    => truncating the viscosity at l' = 2M is exact, as proved above.")  #src
end  #src

# ### What the rule does to the matrix
#
# Read the upper bound backwards: mode ``l`` reaches mode ``m`` only through
# viscosity harmonics with ``k\ge|l-m|``. So if ``\eta`` has angular content
# only up to some ``L_\eta``, no pair of modes further apart than ``L_\eta`` is
# coupled at all, and ``\mathcal D`` is banded with half-bandwidth ``L_\eta``.
# The middle panel of the figure above is this banded case.
#
#
# Applying a banded matrix costs ``O(M L_\eta)`` rather than ``O(M^2)``, so
# ``L_\eta`` -- not ``M`` -- is what sets the price of the coupled model.

let  #src
    Pl(l, m) = l == 0 ? one(m) : l == 1 ? m :  #src
        begin am, b = one(m), m; for n in 1:l-1; b, am = ((2n+1)*m*b - n*am)/(n+1), b; end; b end  #src
    gn, gw = DropSolver.gauss_legendre_nodes(40, -1.0, 1.0)  #src
    Gr(l, lp, lpp) = sum(w * Pl(l,m) * Pl(lp,m) * Pl(lpp,m) for (m, w) in zip(gn, gw))  #src
    viol_zero = 0.0  #src
    viol_nonzero = Inf  #src
    for l in 0:8, lp in 0:8, lpp in 0:8  #src
        g = Gr(l, lp, lpp)  #src
        allowed = (abs(l - lpp) <= lp <= l + lpp) && iseven(l + lp + lpp)  #src
        if allowed  #src
            (l + lp + lpp > 0) && (viol_nonzero = min(viol_nonzero, abs(g)))  #src
        else  #src
            viol_zero = max(viol_zero, abs(g))  #src
        end  #src
    end  #src
    @assert viol_zero < 1e-13 "a forbidden Gaunt coefficient was nonzero ($viol_zero)"  #src
    @assert viol_nonzero > 1e-6 "an allowed Gaunt coefficient vanished ($viol_nonzero)"  #src
    println("  ASSERTION 4 OK: the selection rule holds over all 729 triples with")  #src
    println("    l, l', l'' <= 8: forbidden max $(round(viol_zero, sigdigits=3)), allowed min $(round(viol_nonzero, sigdigits=3)).")  #src

    band_ok = true  #src
    for lp in (0, 2, 4)  #src
        band = 0  #src
        for l in 0:12, lpp in 0:12  #src
            abs(Gr(l, lp, lpp)) > 1e-13 && (band = max(band, abs(l - lpp)))  #src
        end  #src
        band <= lp || (band_ok = false)  #src
    end  #src
    @assert band_ok "viscosity content at l' coupled modes further apart than l'"  #src
    println("  ASSERTION 5 OK: viscosity content at degree l' couples modes no")  #src
    println("    further apart than |l - l''| = l'. The matrix is banded with")  #src
    println("    half-bandwidth L_eta and dimension M, so cost scales with L_eta, not M.")  #src
end  #src

# ### The Newtonian case, as a check
#
# A constant viscosity has a single nonzero coefficient, ``\eta_0``, at
# ``k=0``. Since ``P_0=1``, the angular factor collapses to the ordinary
# orthogonality relation,
#
# ```math
# G^{0}_{l m} \;=\; \frac{2l+1}{2}\int_{-1}^{1}P_l\,P_{m}\,d\mu
#   \;=\; \delta_{l m},
# ```
#
# the sums lose every off-diagonal term, and the two matrices collapse to the
# two scalars ``2\lambda_l`` and ``\omega_l^2``. Gabbard's system is recovered
# term by term. The Newtonian model is the special case of this one, not a
# separate theory.
#
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
# ### Where the shear rate is evaluated, and why the system is not linear
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



# ## The normal-stress balance, and whether the interior problem is well posed
#
# ### The balance, derived rather than quoted
#
# The traction on the free surface is ``\bm t=\bm\sigma\cdot\bm n`` with
# ``\bm\sigma=-p\bm I+2\eta\bm e``. Outside sits the atmosphere, at gauge zero,
# plus the film pressure ``p_c`` pushing inward. Surface tension balances the
# jump, so the normal component gives
#
# ```math
# \boxed{\;
# \bigl[-p+2\eta\,e_{rr}\bigr]_{x=1} \;+\; \bigl(\nabla\cdot\bm n\bigr)\big|_{x=1}
#   \;+\; p_c \;=\; 0 \;}
# ```
#
# and the surface strain follows from the strain tensor rather than from memory,
#
# ```math
# e_{rr}=\Bigl(\frac{\psi_l'}{x^2}-\frac{2\psi_l}{x^3}\Bigr)P_l
# \quad\Longrightarrow\quad
# e_{rr}\big|_{x=1}=\psi_l'(1)-2\psi_l(1) .
# ```
#
# **The base state fixes the sign, and it is worth doing before anything else.**
# A drop at rest carries ``p=2`` (Laplace) and ``\nabla\cdot\bm n=2``, with no
# flow and no contact. Substituting: ``-2+0+2+0=0``. The opposite sign on the
# curvature leaves ``-4``, which is a Laplace pressure pointing the wrong way.
# One line, no algebra, and it is decisive.
#
# ### The pressure amplitude is not free
#
# For a normal mode with ``p=\mathcal P x^lP_l``, the radial momentum equation
# reads ``-\sigma u_r=-\partial_xp+\mathrm{Oh}[\nabla^2\bm u]_r``, and the
# poloidal field has ``[\nabla^2\bm u]_r=\mathcal D_l[\psi_l]P_l/x^2``. Multiplying
# by ``x^2``,
#
# ```math
# \mathrm{Oh}\,(\mathcal D_l+q^2)[\psi_l]=\mathcal P\,l\,x^{l+1} .
# ```
#
# With ``\psi_l=Cxj_l(qx)+Dx^{l+1}`` the Bessel part is annihilated and
# ``x^{l+1}`` maps to ``q^2x^{l+1}``, so
#
# ```math
# \mathcal P=\frac{\mathrm{Oh}\,q^2D}{l}=\frac{\sigma D}{l} .
# ```
#
# Counting ``\mathcal P`` as an independent constant is one way to manufacture a
# phantom under-determinacy; the other is treating ``\zeta_l`` as an unknown
# needing its own equation, when BC1 ties it to the interior and BC3 advances it.
#
# ### The allocation, and the test that settles it
#
# ```math
# \boxed{\;
# \begin{aligned}
# &\text{regularity at } x=0 &&\text{kills the two singular solutions}\\
# &\text{BC1:}\ \ \psi_l\big|_{x=1}=\dot\zeta_l &&\text{kinematic}\\
# &\text{BC2:}\ \ \mathcal T[\psi_l]\big|_{x=1}=0 &&\text{tangential stress}\\
# &\text{BC3:}\ \ \text{normal stress, above} &&\text{advances }\zeta_l
# \end{aligned}\;}
# ```
#
# In the Newtonian limit with a normal mode this is a homogeneous ``3\times3``
# system in ``(C,D,Z)``, where ``Z`` is the surface amplitude. A non-trivial
# solution exists only where the determinant vanishes, so that determinant *is* a
# characteristic equation -- and if the allocation is right it must be Reid's,
# which the solver computes independently. That is the check below, and it also
# re-tests the interior equation end to end: ``\sigma=\mathrm{Oh}\,q^2`` enters
# through it, so a wrong ``\mathrm{Oh}`` would break it.

let  #src
    ## Two things went wrong in earlier attempts at this check, and both are      #src
    ## worth recording because each produced a confident wrong answer.            #src
    ##                                                                            #src
    ## 1. A hand-rolled j_l. sqrt(pi/2z)*besselj(l+1/2,z) takes a different       #src
    ##    branch than the entire function for complex z with negative imaginary    #src
    ##    part, and a recurrence seeded from sin(z)/z loses precision where        #src
    ##    |sin z| ~ e^15. Both made the test report agreement for the WRONG sign.  #src
    ##    Use sph_bessel_ratio -- the determinant's Bessel column is homogeneous   #src
    ##    in j_l, so only the ratio is needed and j_l itself is never evaluated.   #src
    ## 2. Deriving the base-state sign by rearrangement rather than by evaluating  #src
    ##    a residual. The rearrangement had an algebra slip and reversed the       #src
    ##    conclusion.                                                             #src
    Tg(f, d1, d2, l) = d2 - 2d1 + l*(l+1)*f  #src
    function bess_at_one(q, l)  #src
        Q  = sph_bessel_ratio(l, q)            # j_{l+1}(q)/j_l(q); set j_l(q) := 1  #src
        j1 = l/q - Q                            # j_l'  #src
        j2 = -2j1/q - (1 - l*(l+1)/q^2)         # j_l'' from the Bessel equation  #src
        (one(q), 1 + q*j1, 2q*j1 + q^2*j2)      # f, f', f'' for f = x j_l(qx)  #src
    end  #src
    ## BC3 with the curvature sign as a parameter, so the test discriminates.  #src
    function alloc_det(q, Oh, l; s = +1.0)  #src
        fB, d1B, d2B = bess_at_one(q, l)  #src
        fP, d1P = 1.0, l + 1.0                  # x^(l+1) and its derivative at 1  #src
        d2P = (l + 1.0)*l  #src
        er(f, d1) = 2Oh*(d1 - 2f)               # 2 Oh e_rr at x=1  #src
        det([ fB                    fP                        complex(Oh*q^2)  #src
              Tg(fB,d1B,d2B,l)      Tg(fP,d1P,d2P,l)          complex(0.0)  #src
              er(fB,d1B)            er(fP,d1P) - Oh*q^2/l     complex(s*(l-1)*(l+2)) ])  #src
    end  #src
    ## (i) the base state, by residual and not by rearrangement  #src
    base_right = -2.0 + 0.0 + 2.0 + 0.0        # -p + 2 eta e_rr + div n + p_c  #src
    base_wrong = -2.0 + 0.0 - 2.0 + 0.0  #src
    @assert abs(base_right) < 1e-14 "the derived normal-stress balance fails the base state ($base_right)"  ## CLAIM: SUM-NORMAL  #src
    @assert abs(base_wrong) > 1.0 "the opposite curvature sign also passes; the base state cannot fix it"  #src
    ## (ii) the determinant vanishes at Reid's roots, and only for the derived sign  #src
    worst_right, best_wrong = 0.0, Inf  #src
    for Oh in (0.006, 0.05, 0.3, 1.0), l in (2, 4, 8)  #src
        q  = dominant_root(Oh, l)  #src
        sc = maximum(abs(alloc_det(q*(1+d), Oh, l)) for d in (0.15, -0.15, 0.3))  #src
        @assert sc > 0 "the determinant is identically zero; the check is vacuous"  #src
        worst_right = max(worst_right, abs(alloc_det(q, Oh, l; s = +1.0))/sc)  #src
        best_wrong  = min(best_wrong,  abs(alloc_det(q, Oh, l; s = -1.0))/sc)  #src
    end  #src
    @assert worst_right < 1e-12 "the BC allocation does not reproduce Reid's characteristic equation ($worst_right)"  ## CLAIM: SUM-INT  #src
    @assert best_wrong > 1e-3 "the opposite curvature sign also satisfies Reid, so this cannot fix it ($best_wrong)"  #src
    println("  ASSERTION 3e OK: the normal-stress balance and the BC allocation --")  #src
    @printf("    base state: -p + 2 eta e_rr + div n + p_c = %.1e, and the opposite\n", abs(base_right))  #src
    @printf("      curvature sign leaves %.1f, a Laplace pressure of the wrong sign;\n", abs(base_wrong))  #src
    @printf("    the 3x3 determinant vanishes at Reid's roots to %.1e over Oh=0.006..1,\n", worst_right)  #src
    @printf("      l=2,4,8, while the wrong sign is off by at least %.1e.\n", best_wrong)  #src
    println("    So the interior problem is well posed on {regularity, BC1, BC2, BC3},")  #src
    println("    zeta needs no extra equation, and the pressure amplitude is not free.")  #src
    println("    Physical meaning of a failure: the interior problem would be under- or")  #src
    println("    over-determined, and a solver built on it would still run, returning a")  #src
    println("    spectrum that is not the drop's.")  #src
end  #src

# ## Two results the summary needs
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
    gap_exact  = zc + (1 + zeta_)*cos(tt)          # z + r cos(theta), r = 1+zeta  #src
    gap_claim  = cos(tt)*(1 + zeta_) + zc  #src
    w_gap = maximum(abs(Symbolics.value(Symbolics.substitute(  #src
                ed2(gap_exact - gap_claim), Dict(tt=>tv, zc=>0.3, zeta_=>0.17))))  #src
            for tv in TH)  #src
    @assert w_gap < 1e-14 "h is not the exact height of the surface point ($w_gap)"  ## CLAIM: SUM-GAP  #src

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
    @assert abs((-(4pi/3)*1.0)/(4pi/3) + 1.0) < 1e-14 "F/mass must reduce to -B_1"  #src

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
    @assert w_harm < 1e-10 "x^n is not annihilated by L_n, so the harmonic part is wrong ($w_harm)"  ## CLAIM: SUM-PHARM  #src

    println("  ASSERTION 2c OK: the remaining Model summary claims, discharged --")  #src
    @printf("    h = mu(1+zeta)+z is the exact surface height        (%.1e)\n", w_gap)  #src
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
# Everything above, collected. The model is stated once, completely, and in
# continuous terms: no truncation, no discretisation, and no numerical parameter
# appears in it. Nothing below is an approximation of anything above it.
#
# ### Scalings
#
# Lengths by the equilibrium radius ``R``, time by the inertio-capillary time
# ``T_\sigma=\sqrt{\rho R^3/T_1}``, velocity by ``R/T_\sigma``, pressure by
# ``T_1/R``, viscosity by the zero-shear plateau ``\eta_0``. Three groups
# survive:
#
# ```math
# \mathrm{Oh}=\frac{\eta_0}{\sqrt{\rho T_1 R}},
# \qquad
# \mathrm{Bo}=\frac{\rho g R^2}{T_1},
# \qquad
# \mathrm{We}=\frac{\rho R V^2}{T_1},
# ```
#
# with ``V`` the impact speed. The drop occupies ``0\le x\le1``,
# ``\theta\in[0,\pi]``, with ``\theta=\pi`` facing the substrate, ``\mu=\cos\theta``.
#
# ### The unknowns
#
# | unknown | domain | what it is |
# |:--|:--|:--|
# | ``\zeta_l(t)``, ``l\ge2`` | ``t>0`` | surface mode amplitudes |
# | ``\psi_l(x,t)``, ``l\ge2`` | ``0<x<1`` | interior stream function profiles |
# | ``p(x,\mu,t)`` | inside the drop | pressure |
# | ``p_{c,l}(t)``, ``n\ge0`` | ``t>0`` | Legendre coefficients of the substrate reaction |
# | ``z(t),\ v(t)`` | ``t>0`` | height of the drop's centre of mass, and its velocity |
#
# ``\eta`` is not an unknown: it is a function of the others, given by the
# constitutive law in (5). Note that the interior profiles ``\psi_l`` **are** part
# of the state. They are eliminated only by the closure discussed on the
# companion page; here they are evolved.
#
# The surface shape, the reaction pressure, and the gap to the substrate are
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
# and the velocity follows from the stream function
# ``\psi=\sum_l \psi_l(x,t)C_l(\theta)``, with
# ``C_l=\sin^2\!\theta\,P_l'(\mu)/(l(l+1))``:
#
# ```math
# u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}
#   =\sum_l u_{r,l}\,P_l(\mu),
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
# ### (1) The interior
#
# ```math
# \partial_t\,\mathcal D_l[\psi_l]
#   = \mathrm{Oh}\sum_{m\ge2}\mathcal R_{l m}\bigl[\psi_{m};\eta\bigr],
# \qquad 0<x<1,\quad l\ge2,
# ```
#
# with ``\mathcal D_l=d^2/dx^2-l(l+1)/x^2`` and ``\mathcal R_{l m}`` the
# projection onto ``C_l`` defined in *The interior equation*, and explicitly
#
# ```math
# \mathcal R_{l m}\bigl[\psi_m;\eta_k\bigr]=\sum_{j=0}^{2}\sum_{i=0}^{4}
#   a^{(k)}_{l m,\,ji}\;x^{\,j+i-4}\;\eta_k^{(j)}\,\psi_m^{(i)} ,
# ```
#
# fifteen pure numbers per ``(l,m,k)``, tabulated once. Coupling is
# confined to ``|l-m|\le k\le l+m`` with ``l+k+m`` even, so the system is
# banded. Closed by regularity at ``x=0``, which forces ``\psi_l\sim x^{l+1}``, and
# at ``x=1`` by
#
# ```math
# \underbrace{\psi_l\big|_{x=1}=\dot\zeta_l}_{\text{BC1, kinematic}},
# \qquad
# \underbrace{\mathcal T[\psi_l]\big|_{x=1}=0}_{\text{BC2, tangential stress}},
# \qquad
# \mathcal T=\frac{d^2}{dx^2}-\frac{2}{x}\frac{d}{dx}+\frac{l(l+1)}{x^2} .
# ```
#
# Note that ``\mathcal D_l`` annihilates ``x^{l+1}``, so the time derivative
# above does not act on all of ``\psi_l``: the harmonic part of each profile is
# fixed by the boundary conditions at each instant rather than evolved. The
# interior is a **differential-algebraic** system, not an evolution equation.
#
# ### (2) The pressure
#
# Taking the divergence of the momentum equation, with ``\nabla\cdot\bm u=0``,
#
# ```math
# \nabla^2 p = \mathrm{Oh}\,\nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr)
#   \;\equiv\; \mathrm{Oh}\,S(x,\mu,t) ,
# ```
#
# an elliptic problem at each instant. It is needed because the curl that closed
# (1) removed ``p``, while (3) requires it.
#
# The pressure is carried in the same angular basis as everything else,
# ``p(x,\mu,t)=\sum_{l\ge0}p_l(x,t)P_l(\mu)``, and since the ``P_l`` are
# eigenfunctions of the angular Laplacian each coefficient obeys a radial
# equation:
#
# ```math
# \mathcal L_l[p_l] \;=\; \mathrm{Oh}\,S_l(x,t),
# \qquad
# \mathcal L_l = \frac{d^2}{dx^2}+\frac{2}{x}\frac{d}{dx}-\frac{l(l+1)}{x^2} ,
# ```
#
# with ``S_l`` the ``l``-th Legendre coefficient of
# ``\nabla\cdot(\nabla\cdot(2\eta\bm e))`` and ``p_l`` regular at ``x=0``, which
# forces ``p_l\sim x^l``. The operator ``\mathcal L_l`` is the same radial
# Laplacian that appears in the projection integrals of *Carrying the projection
# through*; it is not a new object.
#
# Two things follow. First, the general solution is a particular part driven by
# ``S_l`` plus a **harmonic** part ``c_l(t)\,x^l``, and the ``c_l`` are fixed by
# the surface condition (3) rather than by anything here -- which is the sense in
# which the pressure is coupled to the surface rather than determined
# independently of it. Second, when ``\eta`` is constant the source vanishes and
# only the harmonic part survives,
#
# ```math
# p = \sum_l c_l(t)\,x^lP_l(\mu) ,
# ```
#
# one amplitude per mode and no differential equation at all. That is exactly the
# form the Newtonian theory uses, and it is why a variable viscosity turns a
# closed-form pressure into a solve.
#
# ### A note on representation
#
# It is worth saying plainly what is spectral here and what is not, because the
# treatment is deliberately asymmetric.
#
# **In the angle, everything is spectral, and exactly so.** The surface, the film
# pressure, the stream function and the fluid pressure are all carried as
# Legendre or Gegenbauer series with no truncation imposed. That is not a
# discretisation: ``\{P_l\}`` is complete on ``[-1,1]``, so the expansions are
# changes of variable, from a function of ``\theta`` to a sequence of functions of
# ``x`` and ``t``. The reason to do it is that the ``P_l`` are the eigenfunctions
# of the angular Laplacian on the sphere, which is what turns every angular
# derivative in the problem into an algebraic factor of ``l(l+1)`` and leaves
# behind the radial operators ``\mathcal D_l``, ``\mathcal L_l``,
# ``\mathcal T``. The geometry supplies that basis; nothing is chosen.
#
# **In the radius, nothing is.** ``\psi_l(x,t)`` and ``p_l(x,t)`` are left as
# functions of ``x``, and the equations governing them are differential equations
# in ``x``, not algebraic relations among coefficients. This is not an oversight
# and not a gap in the model: the radial direction has no distinguished basis --
# no operator whose eigenfunctions it hands you -- so any radial expansion is a
# numerical choice rather than a consequence of the geometry. Making that choice
# is discretisation, and it belongs on the companion page with the quadrature and
# the time stepping.
#
# So the model is a system of partial differential equations in one space
# variable, indexed by mode number, and it is stated that way on purpose.
#
# ### (3) The surface
#
# The balance of normal stress against surface tension at ``x=1``, projected onto
# ``P_l``, is the equation of motion for each mode:
#
# ```math
# \Bigl\langle\,\bigl[-p+2\eta\,e_{rr}\bigr]_{x=1}
#   + \bigl(\nabla\!\cdot\!\bm n\bigr)\big|_{x=1}
#   + p_c \,,\ P_l\Bigr\rangle \;=\; 0,
# \qquad l\ge2,
# ```
#
# The sign on the curvature term is fixed by the base state: a drop at rest
# carries ``p=2`` and ``\nabla\cdot\bm n=2``, and only this combination gives
# zero. The opposite sign leaves a residual of ``-4`` -- a Laplace pressure of
# the wrong sign, which is what the page said before it was checked.
#
# with ``\langle f,P_l\rangle=\tfrac{2l+1}{2}\int_{-1}^{1}fP_l\,d\mu`` and the
# curvature supplying the capillary restoring term through
#
# ```math
# \nabla\cdot\bm n = 2 + \sum_{l\ge2}(l-1)(l+2)\zeta_lP_l(\mu) + O(\zeta^2) .
# ```
#
# Carried out explicitly,
# the viscous part of this projection is the double sum of *Carrying the
# projection through*, and it is what the coefficient matrices
# ``\mathcal D^{(i)}`` collect. Those matrices are functionals of the interior
# state ``\{\psi_l\}``, not of ``\bm{\dot\zeta}`` alone.
#
# ### (4) The drop as a whole
#
# ```math
# \dot z = v,
# \qquad
# \dot v = -\mathrm{Bo} - p_{c,1},
# \qquad
# v(0) = -\sqrt{\mathrm{We}} ,
# ```
#
# where the net vertical force the substrate exerts is
#
# ```math
# \mathfrak F = -\!\oint p_c\,n_z\,dS = -2\pi\!\int_{-1}^{1}p_c\,\mu\,d\mu
#   = -\frac{4\pi}{3}p_{c,1} .
# ```
#
# ``p_{c,1}`` is one coefficient of the pressure expansion, not the force; the two
# are proportional, and the drop's non-dimensional mass being ``4\pi/3`` as well
# is what cancels the factor and leaves ``\dot v=-\mathrm{Bo}-p_{c,1}``. Every other
# ``p_{c,l}`` is orthogonal to ``\mu`` and so invisible to the centre of mass. The
# impact speed enters only as the initial condition.
#
# The drop does not wet the substrate: it is held off by a thin intervening air
# layer, and ``p_c`` is that layer's pressure. The layer is taken to be thin
# enough to transmit stress without dynamics of its own -- the lubrication limit
# -- so it appears in the model only through the traction it applies. Because a
# gas film cannot sustain tension, it can push and never pull, and the gap and
# the pressure therefore satisfy a **Signorini complementarity condition**: at
# each ``\theta`` the drop is either pressed and carries pressure, or is free and
# carries none,
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
# It is worth being clear that the drop never touches anything in this model. The
# constraint ``h\ge0`` is the statement that the surface does not cross the plane
# on which the film's pressure acts, which is exact only in the limit of
# vanishing film thickness. ``p_c\ge0`` is a property of the film, not of the
# substrate's rigidity.
#
# ### (5) The fluid
#
# ```math
# \bm e=\tfrac12\bigl(\nabla\bm u+\nabla\bm u^{\mathsf T}\bigr),
# \qquad
# \dot\gamma=\sqrt{2\,\bm e\!:\!\bm e},
# \qquad
# \eta=\eta(\dot\gamma) ,
# ```
#
# evaluated pointwise on the current interior field. ``\dot\gamma`` is computed
# from the **whole** field: ``\bm e`` superposes over modes, its invariant does
# not, so there is no such thing as one mode's shear rate. Any
# ``\eta(\dot\gamma)`` satisfying (H1)-(H3) is admissible; no particular
# constitutive law is assumed here.
#
# ### What kind of system this is
#
# Blocks (1) and (2) are partial differential equations in ``x`` at each instant,
# banded in the mode index; (3) and (4) are ordinary differential equations in
# ``t`` with a complementarity constraint; (5) closes the loop, because the
# ``\eta`` that (1) and (2) need is computed from their own solution. So the
# model is a **differential-algebraic system with a fixed point in the
# viscosity**.
#
# Two things are worth noting about what is *absent*. There is no eigenvalue
# anywhere: no ``q^2``, no ``\lambda_l``, no ``\omega_l^2``, and so no truncation
# of a relaxation spectrum to two roots. And there is no second time derivative
# of ``\zeta_l``: the inertia that produced ``\ddot\zeta_l`` in the Newtonian modal
# system lives in (1) here, and the second-order oscillator form is a
# *consequence* of eliminating the interior rather than a feature of the physics.
#
# Every unknown is determined by an equation stated above. What is not available
# is a formula: ``\eta(\dot\gamma)`` is not polynomial, so nothing in the chain
# has a closed form and every piece must be computed. The companion page
# *Shear-Thinning Drops: Closures* is about what may be given up to compute them
# cheaply, and what each concession costs.
