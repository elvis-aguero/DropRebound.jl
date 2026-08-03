# # Shear-Thinning Drops
#
# In this page we try to extend Reid's work to shear thinning rheologies. These are fluids that have an effective
# viscosity which is shear-rate dependent, possibly nonlinearly.
#
# In the Newtonian theory a drop's surface modes are independent damped
# oscillators: each ``A_l`` has its own damping ``\lambda_l`` and frequency
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
# | ``A_l(t)`` | amplitude of surface mode ``l``; the solver's state vector |
# | ``C_l(\theta)`` | Gegenbauer angular function ``\sin^2\!\theta\,P_l'(\cos\theta)/(l(l+1))``; the *angular* factor of the stream function -- **not** ``A_l``, which is an amplitude in time |
# | ``U_l(x)`` | radial factor of mode ``l`` in the stream function ``\psi=\sum_l U_l(x)C_l(\theta)`` |
# | ``F_l(x),\ W_l(x)`` | radial profiles of ``u_r`` and ``u_\theta`` for mode ``l``; both fixed by ``U_l`` |
# | ``\mathcal D_l`` | ``d^2/dx^2-l(l+1)/x^2``; Reid's radial operator |
# | ``\mathcal L_2`` | ``d^2/dx^2-(2/x)d/dx+l(l+1)/x^2``; Reid's tangential-stress operator, ``\mathcal L_2[U]|_{x=1}=0`` is BC2 |
# | ``\mathcal R_l[U;\eta]`` | variable-viscosity radial operator, derived below |
# | ``\mathcal R_{l l''}`` | its off-diagonal generalisation once ``\eta`` varies with ``\theta`` |
# | ``l'`` | degree index of the **viscosity field's own** Legendre series |
# | ``L_\eta`` | highest ``l'`` present -- the *bandwidth* of the coupling |
# | ``G^{l'}_{l l''},\ H^{l'}_{l l''}`` | Gaunt-type angular integrals; pure numbers |
# | ``\mathrm{Oh}=\eta_0/\sqrt{\rho T_1R}`` | Ohnesorge number: viscous over inertio-capillary stress |
#
# Two collisions are worth pointing out because the surrounding literature
# invites them. ``G^{l'}_{l l''}`` below is an *angular* integral of three
# Legendre polynomials -- a number, not a function; the radial velocity
# profiles are written ``F_l`` and ``W_l`` here precisely so that nothing
# called ``G`` carries an ``x``. And ``U_l`` is a radial profile, while ``A_l``
# is an amplitude in time; Reid's page writes ``G(x)=U(x)/x^2`` for a third
# object again, which does not appear on this page at all.
#
# One symbol is easy to misread. ``l`` and ``l''`` label modes of the drop's
# *shape*, and run ``2\ldots M``, where ``M`` is the truncation of the shape
# expansion (of order 50--90 in practice). ``l'`` labels harmonics of the
# *viscosity field* and is a different index entirely: it starts at ``0``, and
# *The full coupled system* below shows it controls only how far off the
# diagonal the mode-coupling
# matrix reaches -- never the size of that matrix, which is always
# ``M\times M``.
#
# Some results below are established symbolically, for arbitrary ``l``; others
# only by evaluation at a spread of concrete points, which is strong numerical
# evidence rather than a proof. The text says which is which.

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
# \zeta(\theta,t) = R\sum_{l\ge2} A_l(t)\,P_l(\cos\theta),
# ```
#
# the state of the drop is the vector of modal amplitudes ``\{A_l(t)\}``,
# ``l=2\ldots M``. This is exactly the state `julia/src/types.jl` carries.
#
# An axisymmetric incompressible field is generated by a Stokes stream function,
# and expanding *that* in the same angular basis introduces the interior
# unknowns:
#
# ```math
# \psi(x,\theta,t) \;=\; \sum_{l\ge2} U_l(x,t)\,C_l(\theta),
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
# u_r=\sum_l F_l(x)P_l(\mu),\quad F_l=\frac{U_l}{x^2};
# \qquad
# u_\theta=\sum_l W_l(x)\,\partial_\theta P_l(\mu),\quad
# W_l=\frac{U_l'}{x\,l(l+1)} \;}
# ```
#
# Two consequences are used repeatedly below, and both are stated here so that
# every factor in them is on the page rather than in a reader's head. The
# shear strain is proportional to the tangential-stress operator,
#
# ```math
# e_{r\theta}=\frac{\mathcal L_2[U_l]}{2x\,l(l+1)}\,\partial_\theta P_l ,
# \qquad
# \mathcal L_2=\frac{d^2}{dx^2}-\frac{2}{x}\frac{d}{dx}+\frac{l(l+1)}{x^2} ,
# ```
#
# which is why BC2 reduces to ``\mathcal L_2[U_l]|_{x=1}=0``; and the vorticity
# is proportional to Reid's radial operator,
#
# ```math
# \omega_\varphi=-\sum_l\frac{\mathcal D_l[U_l]\,C_l(\theta)}{x\sin\theta},
# \qquad
# \mathcal D_l=\frac{d^2}{dx^2}-\frac{l(l+1)}{x^2},
# ```
#
# which is the factor that turns the curl of the momentum equation into the
# interior problem. The ``C_l`` are orthogonal under
# ``\langle A,B\rangle=\int_0^\pi AB\,d\theta/\sin\theta``, with
#
# ```math
# \langle C_l,C_{l''}\rangle=\frac{2\,\delta_{l l''}}{(2l+1)\,l(l+1)} ,
# ```
#
# and it is that weight which makes the projection below well defined. The
# strain-rate tensor is linear in ``\bm u`` and therefore in ``\{U_l\}``.
#
# What ``U_l`` actually *is* -- the equation that determines it -- is the
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
        ## (c) e_rtheta = L2[U] dP_l/dtheta / (2 x l(l+1))  -- the factor BC2 rests on  #src
        e_rt = ed((Dt(u_r)/xs + Dx(u_th) - u_th/xs)/2)  #src
        claim_ert = L2(Uf,l)*dPl/(2*xs*l*(l+1))  #src
        ## (d) omega_phi = -D_l[U] C_l / (x sin theta)  -- the factor the curl rests on  #src
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
    @assert w_W    < 1e-10 "u_theta radial profile is not U'/(x l(l+1)) ($w_W)"  #src
    @assert rel_ert  < 1e-10 "e_rtheta is not L2[U] dP_l / (2 x l(l+1)) (rel $rel_ert)"  #src
    @assert rel_vort < 1e-10 "vorticity is not -D_l[U] C_l / (x sin theta) (rel $rel_vort)"  #src
    @assert w_orth < 1e-11 "the C_l are not orthogonal under the stated weight ($w_orth)"  #src
    @assert w_norm < 1e-11 "<C_l,C_l> is not 2/((2l+1) l(l+1)) ($w_norm)"  #src
    println("  ASSERTION 2b OK: every factor in the boxed velocity relations, over l=2..6:")  #src
    @printf("    dC_l/dtheta = sin(theta) P_l                    (%.1e)\n", w_dC)  #src
    @printf("    F_l = U_l/x^2                                   (%.1e)\n", w_F)  #src
    @printf("    W_l = U_l'/(x l(l+1))                           (%.1e)\n", w_W)  #src
    @printf("    e_rtheta = L2[U] dP_l/dtheta / (2 x l(l+1))     (rel %.1e)\n", rel_ert)  #src
    @printf("    omega_phi = -D_l[U] C_l / (x sin theta)         (rel %.1e)\n", rel_vort)  #src
    @printf("    <C_l,C_l''> = 2 delta / ((2l+1) l(l+1))          (%.1e, %.1e)\n", w_orth, w_norm)  #src
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
#    tension -- is projected onto ``P_l`` to give the equation for ``A_l``.
#
# With a constant viscosity all three steps preserve the mode. A surface
# deformation ``\propto P_l`` drives an interior field whose angular dependence
# is still ``P_l``, its stress at the surface is still ``\propto P_l``, and step 3
# picks out one equation per mode by orthogonality. That is why Reid obtains one
# characteristic equation per ``l``. Collecting them into a vector
# ``\bm A=(A_2,\ldots,A_M)^{\mathsf T}``,
#
# ```math
# \bm{\ddot A} \;+\; 2\bm\Lambda\,\bm{\dot A} \;+\; \bm\Omega\,\bm A
#   \;+\; \bm b \;=\; 0 ,
# \qquad
# \bm\Lambda=\operatorname{diag}(\lambda_l),
# \quad
# \bm\Omega=\operatorname{diag}(\omega_l^2),
# ```
#
# where ``\bm b`` has entries ``l\,B_l``, the contact pressure the wall applies
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
# ``\omega_\varphi=-\sum_l\mathcal D_l[U_l]\,C_l/(x\sin\theta)`` with
# ``\mathcal D_l=d^2/dx^2-l(l+1)/x^2``. Multiplying the curled equation by
# ``-x\sin\theta`` and projecting onto ``C_l`` -- the Gegenbauer functions are
# orthogonal under ``\langle A,B\rangle=\int_0^\pi AB\,d\theta/\sin\theta`` --
# leaves one equation per mode:
#
# ```math
# \boxed{\;
# \partial_t\,\mathcal D_l[U_l]
#   \;=\; \mathrm{Oh}\sum_{l''}\mathcal R_{l l''}\bigl[U_{l''};\hat\eta\bigr],
# \qquad 0<x<1 \;}
# ```
#
# ```math
# \mathcal R_{l l''}[U_{l''};\hat\eta] \;\equiv\;
#   \frac{\bigl\langle\,
#     -x\sin\theta\,\bigl[\nabla\times\nabla\!\cdot\!(2\hat\eta\bm e)\bigr]_\varphi
#     \bigr|_{\psi=U_{l''}C_{l''}}
#   ,\;C_l\bigr\rangle}{\langle C_l,C_l\rangle} .
# ```
#
# This is the whole interior problem. Two features of it decide everything that
# follows.
#
# **It is stated in time, not in frequency.** Reid's ``q^2=\sigma/\nu`` is an
# eigenvalue and presupposes a single normal mode. Nothing in an impacting drop
# is a normal mode -- ``A_l(t)`` is whatever the impact makes it -- and once
# ``\eta`` depends on the state there is no time-independent operator to take
# eigenvalues *of*. The equation above is a parabolic evolution for the interior
# vorticity and needs no such assumption. Replacing it by an instantaneous
# eigenproblem is a genuine approximation, made and priced on the companion
# page.
#
# **It couples the modes.** The sum over ``l''`` is the substance of the
# problem: mode ``l''`` drives mode ``l`` through the viscosity field. That is
# not a failure of technique. It is the physical statement that a fluid whose
# viscosity varies from place to place cannot respond to each surface harmonic
# independently, because the harmonics no longer see the same fluid. *How far
# the coupling reaches* is the next section.
#
# ### The diagonal, and where Reid sits
#
# Expand the viscosity in the same angular basis,
# ``\eta=\sum_{l'}\eta_{l'}(x)P_{l'}(\cos\theta)``. The term ``l'=0`` is the
# spherically symmetric part, and it contributes only to ``l=l''``. Its
# contribution is a genuine radial operator, and it is worth having explicitly
# because it is the one piece of the problem that can be written in closed form:
#
# ```math
# \boxed{\;
# \mathcal R_{l l}\big|_{l'=0} \;=\; \mathcal R_l[U;\eta] \;=\;
#   \eta\,\mathcal D_l^{\,2}[U]
# \;+\; 2\eta'(x)\,\frac{d}{dx}\!\Bigl(\mathcal D_l[U]-\frac{U'}{x}\Bigr)
# \;+\; \eta''(x)\,\mathcal L_2[U] \;}
# ```
#
# ```math
# \mathcal L_2[U]\;\equiv\;U''-\frac{2}{x}U'+\frac{l(l+1)}{x^2}U .
# ```
#
# Three things read off it.
#
# **Reid is the constant-viscosity case.** Setting ``\eta'=\eta''=0`` leaves
# ``\eta\,\mathcal D_l^2[U]``, so the interior equation becomes
# ``\partial_t\mathcal D_l[U_l]=\mathrm{Oh}\,\mathcal D_l^2[U_l]``: vorticity
# diffusion. Substituting a normal mode ``U_l\propto e^{-\sigma t}`` gives
# ``\mathcal D_l(\mathcal D_l+q^2)[U_l]=0`` with ``q^2=\sigma/\mathrm{Oh}``,
# which is Reid's Eq. 9. Nothing has been added to recover it; the
# constant-viscosity theory is the diagonal of this one, at one harmonic.
#
# **The order of the equation does not change.** The highest derivative is
# ``\eta U''''``; ``\eta'`` reaches only ``U'''`` and ``\eta''`` only ``U''``.
# The operator is fourth order in ``U``, exactly as Reid's is, so the solution
# space is still four-dimensional and the problem still closes on the same
# number of boundary conditions. A variable viscosity changes the coefficients
# of the interior problem; it does not change its type.
#
# **``\eta''`` enters through the tangential-stress operator.** The
# ``\mathcal L_2`` above is not a new object: it is the operator whose vanishing
# at the surface is BC2 below. It appears here because
# ``\mathcal L_2[U]\propto e_{r\theta}``, and the second radial derivative of the
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
# \;\Longleftrightarrow\; \mathcal L_2[U]\big|_{x=1}=0 ,
# ```
#
# where the first equivalence holds because ``\eta(\dot\gamma)\ge\eta_\infty>0``
# everywhere, so the scalar factor cannot vanish. **BC2 is therefore
# rheology-agnostic for every fluid in the admissible class**, and the whole
# ``\tau_{r\theta}=0\Rightarrow\mathcal L_2[U]=0`` reduction carries over
# verbatim.
#
# The third surface condition -- the balance of normal stress against surface
# tension -- is **not** a boundary condition on the interior problem. Projected
# onto ``P_l`` it *is* the equation of motion for ``A_l``, and that projection is
# what produces the coefficient matrices. It is carried out two sections below.
# This is where the rheology lands: once ``\eta`` varies with ``\theta``, its
# value at the surface is a field rather than a number, and projecting that
# field is the entire origin of the mode coupling.
#
@variables rr tt qq hh0 hh1 hh2  #src
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
    ##     R_l[U;eta] = eta*D_l^2[U] + 2*eta' d/dx(D_l[U] - U'/x) + eta''*L2[U]  #src
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
    ## D_l[U] = 0) and eta = 2 - 0.15x^2 the eta' and eta'' terms cancel        #src
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

# ### How far the coupling reaches
#
# The sum over ``l''`` in the interior equation looks unbounded, and if it were,
# the problem would be no easier to solve than the original partial differential
# equation. It is not unbounded. Each viscosity harmonic ``l'`` couples ``l''``
# to ``l`` only when
#
# ```math
# |l-l''|\le l'\le l+l''
# \qquad\text{and}\qquad
# l+l'+l'' \ \text{even},
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
    ## harmonic P_l'(mu) so the (l', l'') pair under test is unambiguous.       #src
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
# ``\dot\zeta=R\sum_{l''}\dot A_{l''}P_{l''}``, so that mode ``l''`` enters with
# strength ``\dot A_{l''}``, the contribution of the viscous stress to the
# equation for ``A_l`` is a double sum over the driving mode ``l''`` and the
# viscosity harmonic ``l'``, of terms
#
# ```math
# \sum_{l''}\sum_{l'} \dot A_{l''}\Bigl[\,
#   G^{l'}_{l l''}\!\!\int_0^1\!\bigl(\eta_{l'}\mathcal L_{l''}[F_{l''}]
#                                      + \eta'_{l'}F'_{l''}\bigr)x^2dx
# \;+\;
#   H^{l'}_{l l''}\!\!\int_0^1\!\frac{\eta_{l'}}{2x}
#      \Bigl(W'_{l''}-\frac{W_{l''}}{x}-\frac{F_{l''}}{x}\Bigr)x^2dx
# \,\Bigr],
# ```
#
# where the two angular factors are
#
# ```math
# G^{l'}_{l l''}=\frac{2l+1}{2}\!\int_{-1}^{1}\!P_l\,P_{l'}\,P_{l''}\,d\mu ,
# \qquad
# H^{l'}_{l l''}=\frac{2l+1}{2}\!\int_{-1}^{1}\!P_l\,(1-\mu^2)\,P'_{l'}\,P'_{l''}\,d\mu ,
# ```
#
# and the radial factors are exactly the ones the three viscous contributions
# produce. From ``\eta\nabla^2\bm u``, the radial Laplacian of the driving mode's
# profile,
# ``\mathcal L_{l''}[F_{l''}]=F''_{l''}+\tfrac{2}{x}F'_{l''}-\tfrac{l''(l''+1)}{x^2}F_{l''}``,
# carried by ``\eta_{l'}`` itself. From the radial part of
# ``2(\nabla\eta)\cdot\bm e``, the radial strain amplitude ``F'_{l''}``, carried by
# ``\eta'_{l'}``. From its polar part, the ``e_{r\theta}`` amplitude
# ``\tfrac12(W'_{l''}-W_{l''}/x-F_{l''}/x)``, carried by ``\eta_{l'}/x`` -- and this
# is the one term that pairs with ``H`` rather than ``G``, because its angular
# factor carries the two derivatives.
#
# ``F_{l''}`` and ``W_{l''}`` are the radial profiles of ``u_r`` and
# ``u_\theta`` introduced with the modal expansion; both are determined by
# ``U_{l''}``, hence by the interior problem of the previous section. They are
# named ``F`` and ``W`` so that the letter ``G`` belongs to the Gaunt
# coefficient alone.
#
# The second angular form is not a new object. Using
# ``(1-\mu^2)P'_n=\tfrac{n(n+1)}{2n+1}(P_{n-1}-P_{n+1})`` and expanding
# ``P'_{l''}`` in Legendre polynomials of lower degree turns it into a finite
# combination of integrals of the first kind at shifted indices. **Every angular
# integral in the problem is therefore of one type**, and it is worth a name:
#
# ```math
# G^{l'}_{l l''} \;\equiv\; \frac{2l+1}{2}\int_{-1}^{1}P_l\,P_{l'}\,P_{l''}\,d\mu
# ```
#
# a Gaunt coefficient, pure geometry, depending on three integers and on
# nothing about the fluid. Likewise the radial integral depends only on
# ``(l,l',l'')`` and on the current viscosity profile; write it
# ``A^{(i)}_{l l''}[\eta_{l'}]`` and ``B^{(i)}_{l l''}[\eta_{l'}]`` -- the two
# radial integrals written out above, pairing with ``G`` and ``H`` respectively.
# The index ``i`` distinguishes the terms that end up multiplying ``\dot A``
# from those multiplying ``A``.
#
# ### The result
#
# With those names the double sum collapses. Both angular factors obey the same
# selection rule, proved below, so the sum over ``l'`` terminates. Define
#
# ```math
# \mathcal D^{(i)}_{l l''} \;=\; \sum_{l'}\Bigl[\,
#   G^{l'}_{l l''}\,A^{(i)}_{l l''}[\eta_{l'}]
#   \;+\; H^{l'}_{l l''}\,B^{(i)}_{l l''}[\eta_{l'}]\,\Bigr],
# ```
#
# and the modal system takes exactly the Newtonian form with the two diagonal
# matrices replaced by full ones:
#
# ```math
# \boxed{\;
# \bm{\ddot A} \;+\; \mathcal D^{(2)}\,\bm{\dot A}
#              \;+\; \mathcal D^{(1)}\,\bm A \;+\; \bm b \;=\; 0 \;}
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
# one number per ``(l,l'')`` pair. ``G`` and ``H`` integrate the angle away over
# ``\mu\in[-1,1]``; ``A^{(i)}`` and ``B^{(i)}`` integrate the radius away over
# ``x\in[0,1]``.
#
# What does *not* integrate away is the dependence on the **state**. The
# coefficients ``\eta_{l'}(x)`` entering those integrals are the Legendre
# coefficients of ``\eta\bigl(\dot\gamma(x,\theta,t)\bigr)``, and ``\dot\gamma``
# is built from the current modal velocities. So
#
# ```math
# \mathcal D^{(i)} \;=\; \mathcal D^{(i)}\bigl[\bm{\dot A}(t)\bigr] :
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
# The angular factor ``G^{l'}_{l l''}`` is not merely small for most index
# triples -- it is exactly zero for most of them. The rule is
#
# ```math
# G^{l'}_{l l''} = 0
# \qquad\text{unless}\qquad
# |l-l''| \;\le\; l' \;\le\; l+l''
# \qquad\text{and}\qquad
# l+l'+l''\ \text{is even}.
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
# Now ``P_l`` has degree exactly ``l``, so the product ``P_lP_{l''}`` is a
# polynomial of degree exactly ``l+l''``. Therefore:
#
# * **Upper bound.** If ``l'>l+l''`` then ``P_{l'}`` is orthogonal to
#   ``P_lP_{l''}``, a polynomial of degree ``l+l''<l'``. Hence
#   ``G^{l'}_{l l''}=0``.
# * **Lower bound.** Take ``l\ge l''`` without loss of generality. If
#   ``l'<l-l''`` then ``P_{l'}P_{l''}`` has degree ``l'+l''<l``, and ``P_l`` is
#   orthogonal to it. Hence ``G^{l'}_{l l''}=0`` again.
# * **Parity.** ``P_l(-\mu)=(-1)^lP_l(\mu)``, so the integrand has parity
#   ``(-1)^{l+l'+l''}``. An odd integrand over ``[-1,1]`` integrates to zero.
#
# Together: ``G^{l'}_{l l''}=0`` unless ``|l-l''|\le l'\le l+l''`` **and**
# ``l+l'+l''`` is even. No special-function machinery is needed -- degree
# counting and parity give the whole selection rule, and both bounds are the
# *same* argument applied to different factors.
#
# The consequence for this model is the upper bound. The shape expansion stops
# at ``l=M``, so ``l+l''\le 2M``, so viscosity harmonics above ``l'=2M`` are
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
# Read the upper bound backwards: mode ``l`` reaches mode ``l''`` only through
# viscosity harmonics with ``l'\ge|l-l''|``. So if ``\eta`` has angular content
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
# ``l'=0``. Since ``P_0=1``, the angular factor collapses to the ordinary
# orthogonality relation,
#
# ```math
# G^{0}_{l l''} \;=\; \frac{2l+1}{2}\int_{-1}^{1}P_l\,P_{l''}\,d\mu
#   \;=\; \delta_{l l''},
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
# The matrices need ``\eta_{l'}(x,t)``, which comes from the fluid. Under the
# Cross law -- the model this repository's validation fluid is characterised
# with, and its ``p=-1`` slice --
#
# ```math
# \frac{\eta(x,\theta,t)}{\eta_0}
#   \;=\; \varepsilon_\infty
#       \;+\; \frac{1-\varepsilon_\infty}
#                  {1+\bigl(K\,\dot\gamma(x,\theta,t)\bigr)^{m}} ,
# \qquad
# \varepsilon_\infty \equiv \frac{\eta_\infty}{\eta_0},
# ```
#
# where ``K`` is the fluid's time constant (written ``\lambda_c`` in the
# Carreau-Yasuda parametrisation) and ``m`` its thinning exponent (written
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
# ``\dot\zeta=\sum_l\dot A_lP_l``, so the interior velocity is driven by the
# modal *velocities*, and superposes over modes. Strain rate is linear in
# velocity, so the strain **tensor** superposes too:
#
# ```math
# \bm e(x,\theta,t) \;=\; \sum_{l=2}^{M}\dot A_l(t)\;\bm e^{(l)}(x,\theta).
# ```
#
# The invariant does not:
#
# ```math
# \dot\gamma \;=\; \sqrt{2\,\bm e\!:\!\bm e}
#   \;=\; \Bigl(2\sum_{l,l''}\dot A_l\dot A_{l''}\;
#          \bm e^{(l)}\!:\!\bm e^{(l'')}\Bigr)^{1/2}.
# ```
#
# It is quadratic in the field and *then* square-rooted, so the ``l\neq l''``
# cross terms do not drop out. Once more than one mode is active there is no
# such thing as "mode ``l``'s shear rate".
#
# Nothing about ``\eta(\dot\gamma)`` is polynomial, so the coefficients
# ``\eta_{l'}(x)`` have no closed form, and neither do the radial and angular
# integrals built from them. ``\eta_{l'}`` is defined by its Legendre
# projection at each radius,
#
# ```math
# \eta_{l'}(x,t)=\frac{2l'+1}{2}\int_{-1}^{1}
#   \eta\bigl(\dot\gamma(x,\mu,t)\bigr)P_{l'}(\mu)\,d\mu ,
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
# \ddot A_l + \sum_{l''}\Bigl[\mathcal D^{(2)}_{l l''}(\bm{\dot A})\,\dot A_{l''}
#   + \mathcal D^{(1)}_{l l''}(\bm{\dot A})\,A_{l''}\Bigr] + l\,B_l = 0 .
# ```
#
# This is *quasi-linear*: linear in ``(A_l,\dot A_l)`` at frozen coefficients,
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
# | ``A_l(t)``, ``l\ge2`` | ``t>0`` | surface mode amplitudes |
# | ``U_l(x,t)``, ``l\ge2`` | ``0<x<1`` | interior stream function profiles |
# | ``p(x,\mu,t)`` | inside the drop | pressure |
# | ``B_n(t)``, ``n\ge0`` | ``t>0`` | Legendre coefficients of the substrate reaction |
# | ``z(t),\ v(t)`` | ``t>0`` | height of the drop's centre of mass, and its velocity |
#
# ``\eta`` is not an unknown: it is a function of the others, given by the
# constitutive law in (5). Note that the interior profiles ``U_l`` **are** part
# of the state. They are eliminated only by the closure discussed on the
# companion page; here they are evolved.
#
# The surface shape, the reaction pressure, and the gap to the substrate are
#
# ```math
# \zeta(\theta,t)=\sum_{l\ge2}A_l(t)P_l(\mu),
# \qquad
# \Pi(\theta,t)=\sum_{n\ge0}B_n(t)P_n(\mu),
# ```
# ```math
# h(\theta,t)=\mu\,\bigl[1+\zeta(\theta,t)\bigr]+z(t) ,
# ```
#
# and the velocity follows from the stream function
# ``\psi=\sum_l U_l(x,t)C_l(\theta)``, with
# ``C_l=\sin^2\!\theta\,P_l'(\mu)/(l(l+1))``:
#
# ```math
# u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}
#   =\sum_l F_l\,P_l(\mu),
# \qquad
# u_\theta=-\frac{1}{x\sin\theta}\frac{\partial\psi}{\partial x}
#   =\sum_l W_l\,\partial_\theta P_l(\mu) ,
# ```
# ```math
# F_l=\frac{U_l}{x^2},
# \qquad
# W_l=\frac{U_l'}{x\,l(l+1)} .
# ```
#
# ### (1) The interior
#
# ```math
# \partial_t\,\mathcal D_l[U_l]
#   = \mathrm{Oh}\sum_{l''\ge2}\mathcal R_{l l''}\bigl[U_{l''};\eta\bigr],
# \qquad 0<x<1,\quad l\ge2,
# ```
#
# with ``\mathcal D_l=d^2/dx^2-l(l+1)/x^2`` and ``\mathcal R_{l l''}`` the
# projection onto ``C_l`` defined in *The interior equation*. Coupling is
# confined to ``|l-l''|\le l'\le l+l''`` with ``l+l'+l''`` even, so the system is
# banded. Closed by regularity at ``x=0``, which forces ``U_l\sim x^{l+1}``, and
# at ``x=1`` by
#
# ```math
# \underbrace{U_l\big|_{x=1}=\dot A_l}_{\text{BC1, kinematic}},
# \qquad
# \underbrace{\mathcal L_2[U_l]\big|_{x=1}=0}_{\text{BC2, tangential stress}},
# \qquad
# \mathcal L_2=\frac{d^2}{dx^2}-\frac{2}{x}\frac{d}{dx}+\frac{l(l+1)}{x^2} .
# ```
#
# Note that ``\mathcal D_l`` annihilates ``x^{l+1}``, so the time derivative
# above does not act on all of ``U_l``: the harmonic part of each profile is
# fixed by the boundary conditions at each instant rather than evolved. The
# interior is a **differential-algebraic** system, not an evolution equation.
#
# ### (2) The pressure
#
# ```math
# \nabla^2 p = \mathrm{Oh}\,\nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr),
# \qquad\text{regular at } x=0 ,
# ```
#
# an elliptic problem at each instant, whose source vanishes only for constant
# ``\eta``. It is needed because the curl that closed (1) removed ``p``, while
# (3) requires it.
#
# ### (3) The surface
#
# The balance of normal stress against surface tension at ``x=1``, projected onto
# ``P_l``, is the equation of motion for each mode:
#
# ```math
# \Bigl\langle\,\bigl[-p+2\eta\,e_{rr}\bigr]_{x=1}
#   - \bigl(\nabla\!\cdot\!\bm n\bigr)\big|_{x=1}
#   + \Pi \,,\ P_l\Bigr\rangle \;=\; 0,
# \qquad l\ge2,
# ```
#
# with ``\langle f,P_l\rangle=\tfrac{2l+1}{2}\int_{-1}^{1}fP_l\,d\mu`` and
# ``\nabla\cdot\bm n`` twice the mean curvature of the deformed surface, whose
# linearisation supplies the capillary restoring term. Carried out explicitly,
# the viscous part of this projection is the double sum of *Carrying the
# projection through*, and it is what the coefficient matrices
# ``\mathcal D^{(i)}`` collect. Those matrices are functionals of the interior
# state ``\{U_l\}``, not of ``\bm{\dot A}`` alone.
#
# ### (4) The drop as a whole
#
# ```math
# \dot z = v,
# \qquad
# \dot v = -\mathrm{Bo} - \mathfrak F,
# \qquad
# \mathfrak F(t) = \oint \Pi\,n_z\,dS ,
# \qquad
# v(0) = -\sqrt{\mathrm{We}} ,
# ```
#
# where ``\mathfrak F`` is the net vertical force the substrate exerts, obtained
# by integrating the reaction pressure over the contact region. It is a
# functional of ``\Pi``, and therefore of the ``B_n``; it is **not** ``B_1``,
# which is one coefficient of that expansion. The impact speed enters only as
# the initial condition.
#
# The substrate is rigid and cannot pull, so the gap and the reaction satisfy a
# **Signorini complementarity condition**: at each ``\theta`` the drop either
# touches and carries pressure, or does not touch and carries none,
#
# ```math
# h\ge0,
# \qquad
# \Pi\ge0,
# \qquad
# h\,\Pi = 0
# \qquad\text{for all }\theta\in[0,\pi] .
# ```
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
# of ``A_l``: the inertia that produced ``\ddot A_l`` in the Newtonian modal
# system lives in (1) here, and the second-order oscillator form is a
# *consequence* of eliminating the interior rather than a feature of the physics.
#
# Every unknown is determined by an equation stated above. What is not available
# is a formula: ``\eta(\dot\gamma)`` is not polynomial, so nothing in the chain
# has a closed form and every piece must be computed. The companion page
# *Shear-Thinning Drops: Closures* is about what may be given up to compute them
# cheaply, and what each concession costs.
