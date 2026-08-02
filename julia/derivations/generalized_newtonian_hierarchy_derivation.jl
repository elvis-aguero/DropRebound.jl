# # Shear-Thinning Drops
#
# In this page we try to extend Reid's work to shear thinning rheologies. These are fluids that have an effective
# viscosity which is shear-rate dependent, possibly nonlinearly.
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
# | ``U(x)`` | radial factor of the stream function ``\psi=U(x)C_l(\theta)`` |
# | ``\mathcal D_l`` | ``d^2/dx^2-l(l+1)/x^2``; Reid's radial operator |
# | ``\mathcal L_2`` | ``d^2/dx^2-(2/x)d/dx+l(l+1)/x^2``; Reid's tangential-stress operator, ``\mathcal L_2[U]|_{x=1}=0`` is BC2 |
# | ``l'`` | degree index of the **viscosity field's own** Legendre series |
# | ``L_\eta`` | highest ``l'`` retained -- the *bandwidth* of the coupling |
# | ``\mathrm{Oh}=\eta_0/\sqrt{\rho T_1R}`` | Ohnesorge number: viscous over inertio-capillary stress |
# | ``\mathrm{Oh}_{\mathrm{eff}}`` | the single Ohnesorge number a scalar closure substitutes for the field |
#
# One symbol is easy to misread. ``l`` and ``l''`` label modes of the drop's
# *shape*, and run ``2\ldots M``, where ``M`` is the truncation of the shape
# expansion (of order 50--90 in practice). ``l'`` labels harmonics of the
# *viscosity field* and is a different index entirely: it starts at ``0``, and
# Step 4 shows it controls only how far off the diagonal the mode-coupling
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
# For each mode the radial velocity profile, and the strain components that
# follow from it, are Reid's; they are derived on the page *The Viscous Drop:
# Reid (1960)*. Only two of their properties are needed below: the profile is
# fixed once ``l`` and the viscosity are given, and the strain tensor is linear
# in the modal velocity ``\dot A_l``.

# ## The full coupled system
#
# This is the rung to remember. Everything below it is a simplification *of*
# this, and anything you later decide you need, you recover by climbing back
# here.
#
# ### Where a mode equation comes from
#
# It is worth being explicit about how a single-mode equation is obtained at
# all, because that is where the coupling will enter. For a given surface
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
# where ``\lambda_l`` and ``\omega_l^2`` are the two roots of Reid's
# characteristic equation and ``\bm b`` has entries ``l\,B_l``, the contact
# pressure the wall applies while the drop is touching it. ``\bm b`` is
# kinematic and geometric; no rheology enters it, and nothing below changes it.
# Both matrices are diagonal, and that is the entire content of the Newtonian
# model.
#
# ### Which constitutive models this covers
#
# The chain below uses exactly three properties of ``\eta``, and isolating them
# defines the class of fluids it applies to.
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
# ### Following Reid's route, one step at a time
#
# Reid's derivation has a fixed shape: linearise; take the curl of the momentum
# equation to eliminate the pressure; reduce what is left to a radial ODE; solve
# it; impose the three boundary conditions; read off the characteristic
# equation. The question is not whether a shear-thinning fluid has an analogous
# route -- it is *where along that route* the extra term
# ``2(\nabla\eta)\cdot\bm e`` first does damage. Taking the steps in order
# answers that, and the answer is later than one might expect.
#
# **(i) The curl still removes the pressure.** That step is indifferent to the
# rheology: ``\nabla\times\nabla p=0`` whatever ``\eta`` does. Represent the
# axisymmetric field by the Stokes stream function ``\psi=U(x)\,C_l(\theta)``,
# where
#
# ```math
# C_l(\theta)\;=\;\frac{\sin^2\!\theta\;P_l'(\cos\theta)}{l(l+1)},
# \qquad
# u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta},
# \qquad
# u_\theta=-\frac{1}{x\sin\theta}\frac{\partial\psi}{\partial x},
# ```
#
# is the Gegenbauer angular function -- the same field Reid writes as
# ``u_r=f(x)P_l(\cos\theta)``, since ``\partial_\theta C_l\propto\sin\theta P_l``.
# Taking the ``\varphi``-component of the curl of the momentum equation kills the
# pressure and leaves, for ``\eta=\eta(x)``, exactly
#
# ```math
# \bigl[\nabla\times\bm M\bigr]_\varphi
#   \;=\; -\frac{C_l(\theta)}{x\sin\theta}\;
#     \Bigl\{\,q^2\,\mathcal D_l[U] \;+\; \mathcal R_l[U;\eta]\,\Bigr\},
# \qquad
# \mathcal D_l \equiv \frac{d^2}{dx^2}-\frac{l(l+1)}{x^2},
# ```
#
# with the **variable-viscosity radial operator**
#
# ```math
# \boxed{\;
# \mathcal R_l[U;\eta] \;=\;
#   \eta\,\mathcal D_l^{\,2}[U]
# \;+\; 2\eta'(x)\,\frac{d}{dx}\!\Bigl(\mathcal D_l[U]-\frac{U'}{x}\Bigr)
# \;+\; \eta''(x)\,\mathcal L_2[U] \;}
# ```
#
# ```math
# \mathcal L_2[U]\;\equiv\;U''-\frac{2}{x}U'+\frac{l(l+1)}{x^2}U .
# ```
#
# This is the interior equation, and it is worth reading three things off it.
#
# **Reid is the constant-viscosity case.** Setting ``\eta'=\eta''=0`` leaves
# ``\eta\mathcal D_l^2[U]+q^2\mathcal D_l[U]``, i.e. Reid's
# ``\mathcal D_l(\mathcal D_l+q^2)[U]=0`` after dividing by ``\eta`` and
# absorbing it into ``q^2=\sigma/\nu``. Nothing was assumed to get there; it
# is recovered.
#
# **The order of the equation does not change.** The highest derivative is
# ``\eta U''''`` in the first term; ``\eta'`` reaches only ``U'''`` and
# ``\eta''`` only ``U''``. So ``\mathcal R_l`` is fourth order in ``U``, exactly
# as Reid's operator is, the solution space is still four-dimensional, and the
# problem still closes on the *same number* of boundary conditions. A variable
# viscosity changes the coefficients of the interior problem; it does not change
# its type.
#
# **``\eta''`` enters through Reid's own tangential-stress operator.** The
# ``\mathcal L_2`` appearing above is not a new object: it is the operator whose
# vanishing at the surface *is* Reid's BC2, ``\mathcal L_2[U]|_{x=1}=0``, derived
# on the page *The Viscous Drop: Reid (1960)*. It appears here because
# ``\mathcal L_2[U]\propto e_{r\theta}``, and the second radial derivative of the
# viscosity couples to precisely the shear component of the strain.
#
# **(ii) The boundary conditions, one at a time.** Reid closes the problem with
# three conditions. Under ``\eta=\eta(\dot\gamma)`` they fare as follows.
#
# ```math
# \text{BC1 (kinematic):}\qquad u_r\big|_{x=1}=\frac{\partial\zeta}{\partial t}
# ```
#
# contains no viscosity at all, and is **unchanged**.
#
# ```math
# \text{BC2 (tangential):}\qquad \tau_{r\theta}\big|_{x=1}=2\eta\,e_{r\theta}=0
# \;\Longleftrightarrow\; e_{r\theta}\big|_{x=1}=0
# \;\Longleftrightarrow\; \mathcal L_2[U]\big|_{x=1}=0 ,
# ```
#
# where the first equivalence holds because ``\eta(\dot\gamma)\ge\eta_\infty>0``
# everywhere, so the scalar factor cannot vanish. **BC2 is therefore
# rheology-agnostic for every generalized Newtonian fluid**, and Reid's whole
# ``\tau_{r\theta}=0\Rightarrow\mathcal L_2[U]=0`` chain carries over verbatim.
#
# ```math
# \text{BC3 (normal):}\qquad -p_{rr}\big|_{x=1}
#   = p+\delta p-2\eta_s\frac{\partial u_r}{\partial x}\Big|_{x=1},
# \qquad \eta_s\equiv\eta\bigl(\dot\gamma|_{x=1}\bigr),
# ```
#
# in which ``\eta`` is genuinely multiplicative, so it survives as the *surface
# value* ``\eta_s``. The form of the condition is unchanged; one coefficient
# becomes state-dependent.
#
#
# **(iii) A radially varying viscosity keeps one ODE per mode.** For
# ``\eta=\eta(x)`` the display above already says it: every ``\theta``-dependence
# sits in the single factor ``C_l(\theta)``, so dividing it out leaves a pure
# radial equation,
#
# ```math
# q^2\,\mathcal D_l[U] \;+\; \mathcal R_l[U;\eta] \;=\; 0 ,
# ```
#
# one such equation for each ``l``, with no reference to any other mode.
# Separability is intact and the modes do not talk to each other. **This is the
# rung Step 7 is about.**
#
# **(iv) An angularly varying viscosity destroys separability.** Let
# ``\eta=\eta(x,\theta)``. Now ``\partial_\theta\eta\neq0``, the curl acquires
# terms in ``\partial_\theta\eta`` and ``\partial^2_{\theta\theta}\eta``, and it
# no longer factors: the ratio
# ``[\nabla\times\bm M]_\varphi\big/C_l(\theta)`` is a function of ``\theta``,
# not a constant. There is no radial operator to write down, because there is no
# ``\theta``-independent thing left after dividing out the angular factor.
#
# It is worth saying exactly what is lost, because "one mode can no longer be
# isolated" is doing real work. Let ``\mathcal A`` be the linearised operator of
# the whole problem and ``\Lambda_\theta`` the angular Laplacian, whose
# eigenfunctions are the ``P_l``. For constant -- or purely radial -- viscosity,
#
# ```math
# [\mathcal A,\Lambda_\theta]=0 ,
# ```
#
# so the two commute, share an eigenbasis, and ``\mathcal A`` is *block diagonal*
# in ``\{P_l\}``:
#
# ```math
# \langle P_l|\mathcal A|P_{l''}\rangle \;=\; 0 \quad\text{for } l\neq l'' .
# ```
#
# That is what "isolates one mode" means, and it is exactly the simultaneous
# diagonalisation Reid's construction rests on. Once ``\eta`` depends on
# ``\theta``, ``[\mathcal A,\Lambda_\theta]\neq0``, no common eigenbasis exists,
# and those off-diagonal projections are nonzero -- which is not a failure of
# technique but the statement that mode ``l''`` genuinely drives mode ``l``. The
# ``\mathcal D`` matrices of the previous section are those projections.
#
# The mechanism is therefore located precisely. The extra term
# ``2(\nabla\eta)\cdot\bm e`` is carried through every one of Reid's steps: it is
# inert at the curl, it merely re-decorates the radial operator at (iii), and
# only when ``\eta`` acquires *angular* structure does it break the commutation
# that the whole construction rests on.
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

    println("  ASSERTION 3b OK: following Reid's route with a variable viscosity --")  #src
    println("    (i) the curl removes the pressure and leaves eta only through its")  #src
    println("        first and second derivatives, never third or higher;")  #src
    println("    (ii) constant eta returns D_l(D_l+q^2)U exactly (to $(round(maximum(abs(x-1) for x in rats), sigdigits=2)));")  #src
    println("    (ii-b) the boxed operator R_l = eta*D_l^2 + 2eta' d/dx(D_l - d/dx /x) + eta''*L2")  #src
    println("         reproduces the curl for variable eta(x) (to $(round(worst_box, sigdigits=2)));")  #src
    println("    (iii) eta(r) preserves separability (spread $(round(sep, sigdigits=2))), so one")  #src
    println("         radial ODE per mode survives -- this is the A7 rung;")  #src
    println("    (iv) eta(r,theta) destroys it (spread $(round(coup, sigdigits=3))), which is the coupling.")  #src
end  #src

# ### Carrying the projection through
#
# Project as before: multiply by ``P_l(\mu)`` and integrate over
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
#      \Bigl(G'_{l''}-\frac{G_{l''}}{x}-\frac{F_{l''}}{x}\Bigr)x^2dx
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
# ``\tfrac12(G'_{l''}-G_{l''}/x-F_{l''}/x)``, carried by ``\eta_{l'}/x`` -- and this
# is the one term that pairs with ``H`` rather than ``G``, because its angular
# factor carries the two derivatives.
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
# -- a Gaunt coefficient, pure geometry, depending on three integers and on
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
# !!! note "What is derived here and what is not"
#     The angular structure above is complete, and the selection rule proved
#     below follows from it. The radial integrals ``R^{(i)}`` are *defined* by
#     the expression above but are not evaluated in closed form on this page:
#     doing so means solving the interior problem of the previous section with
#     a variable ``\eta``, which is the open work Step 7 identifies. Step 7 does
#     construct that radial operator for ``l'=0``, where it can be checked
#     against Reid's.
#
# ### The Newtonian case, as a check
#
# A constant viscosity has one nonzero harmonic, ``\eta_0``, at ``l'=0``. Since
# ``P_0=1`` the angular factor is the ordinary orthogonality relation,
#
# ```math
# G^{0}_{l l''}=\frac{2l+1}{2}\int_{-1}^{1}P_l\,P_{l''}\,d\mu=\delta_{l l''},
# ```
#
# so both sums lose every off-diagonal term and ``\mathcal D^{(2)}``,
# ``\mathcal D^{(1)}`` collapse to ``2\bm\Lambda`` and ``\bm\Omega``. Gabbard's
# system returns term by term. The Newtonian model is the special case of this
# one, not a separate theory.
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
# ### Closing the system: the constitutive law
#
# The matrices need ``\eta_{l'}(x,t)``, which comes from the fluid. Under the
# Cross law -- the model this repository's validation fluid is characterised
# with, and the ``p=-1`` slice of Carreau-Yasuda derived in Step 9 --
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
# such thing as "mode ``l``'s shear rate", and treating each mode as carrying
# its own is a modelling error, not a simplification.
#
# In practice ``\dot\gamma`` is evaluated **pointwise on the ``(x,\theta)``
# quadrature grid** from the full superposition; ``\eta(\dot\gamma)`` follows
# pointwise from the constitutive law; and a Legendre projection in ``\theta``
# at each radius gives the coefficients ``\eta_{l'}(x)`` that the matrices need.
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
# The solver closes it by **lagging**: ``\mathcal D`` is evaluated at the
# previous step's ``\bm{\dot A}``, so within a single Newton step the
# coefficients are constants and the Jacobian is exact for that step. The
# nonlinearity is carried across steps rather than inside them, which costs one
# order in ``\Delta t`` and constrains the step size.
#
# That completes the model. The modal system, the matrix assembly, the
# constitutive law and this prescription for ``\dot\gamma`` are the whole of
# it; every step below either justifies one of those pieces or simplifies the
# double sum.
#


# ## Step 5 (A5) -- the temporal closure
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

# ## Step 6 (A6) -- truncate the viscosity spectrum at ``L_\eta``
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
# solver produces, the practical move is to skip this rung and either take L4
# dense or drop to A7.

# ## Step 7 (A7) -- keep only ``l'=0``: a spherically symmetric viscosity
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
# !!! warning "A7 is cheap, but it is not a small correction"
#     The same measurement that killed A6 also prices A7. At the physical
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
#     with a dozen modes beating, the nodal collapse described in Step 6
#     makes ``\eta`` span ``100\times``. Any
#     scalar-``\mathrm{Oh}_{\mathrm{eff}}`` model inherits this error.
#
# **A7 is the cheapest defensible rung, not an accurate one.** It keeps the
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
# diagonal in ``l`` exactly as Step 4 predicts, and each mode again reduces to a
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
# ## Step 8 (A8) -- freeze the radial profile: ``\eta(r)\to\eta_{\rm eff}``
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
# of Step 5.
#
# Its error is the one measured in Step 7: for a realistic multi-mode state the
# viscosity spans a factor of order ``100`` across the drop, and the angular
# coefficients ``|\eta_1|/|\eta_0|\approx1.3`` and
# ``|\eta_2|/|\eta_0|\approx0.2`` are comparable to or larger than the mean.
# Collapsing all of that onto one number per mode is a leading-order
# approximation, not a perturbative one. Combined with the temporal closure of
# Step 5(a), which is unlicensed at any amplitude, this rung carries two
# uncontrolled approximations at once -- which is worth knowing when reading its
# predictions.
#
# **How to undo it.** Step 7.

# ## Step 9 -- the Cross model as a rung, not an alternative
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
# genuinely empirical question, and Step 6 measures it.
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

# ## Step 10 -- the Newtonian floor, and a live cross-check
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

# ## Summary of the chain
#
# | rung | model | assumption to get here | coupling structure |
# |:--|:--|:--|:--|
# | **0** | exact free-surface generalized Newtonian | -- | -- |
# | **A1** | linear in amplitude | ``\epsilon\ll1``; drops advection | ``\eta`` still fully nonlinear |
# | **A2** | axisymmetric | axisymmetric forcing -- **exact** | ``Y_l^m\to P_l``, no ``m``-coupling ever |
# | **A3** | poloidal + modal | change of variables | state is ``\{a_l\}``, ``l=2\ldots M`` |
# | **L4** | **full coupled system** | none beyond A1--A3 | dense; bandwidth ``=`` spectral width of ``\eta`` |
# | **A5** | temporal closure | (a) instant / (b) period-avg / (c) Floquet | picks which ``\eta(t)`` channel survives |
# | **A6** | truncate at ``L_\eta`` | discarded coupling small -- **measured false** | banded, but needs ``L_\eta\gtrsim M``: no saving |
# | **A7** | ``\eta=\eta(r)`` | viscosity spherically symmetric -- **leading-order error** | **diagonal**; numerical radial BVP per mode |
# | **A8** | ``\eta\to\eta_{\rm eff}`` | radial variation small | Bessel closed form -- **current code** |
# | **A9** | ``\eta=\eta_0`` | no thinning | Reid |
# | **A10** | ``Oh\to0`` | inviscid | Lamb |
#
# Cross sits alongside the chain as the ``p=-1`` slice (Step 9), and can be
# taken at any rung.
#
# ## What this file establishes
#
# **Free, and permanent.** A2: a generalized Newtonian fluid cannot break
# axisymmetry, so Legendre polynomials suffice forever. BC1 and BC2 carry
# over from Reid untouched; BC3 changes by one coefficient. The scope of
# "adapt Reid for shear thinning" is one extra momentum term, ``2(\nabla\eta)\cdot\bm e``.
#
# **Ruled out.** A5(a), evaluating ``\eta`` instantaneously -- what the
# production code does -- is not licensed at any amplitude, because the
# period-``\pi`` lemma puts the modulation at exactly twice the mode
# frequency. A6, banding the coupling matrix, buys nothing: the measured
# ``L_\eta`` is ``\gtrsim M`` for algebraically decaying spectra and ``\gg M``
# for real solver states. (An exponentially decaying spectrum would band
# happily; this solver does not produce one.)
#
# **Priced, and not cheap.** A7 restores exact diagonality, but discards
# angular structure that is *comparable to or larger than* the mean
# viscosity. It is the cheapest defensible rung, not an accurate one.
#
# **The conclusion.** Between L4 (dense) and A7 (leading-order error) there
# is no cheap-and-accurate rung. Any model built on a single scalar
# ``\mathrm{Oh}_{\mathrm{eff}}`` per mode -- which is every shear-thinning
# model implemented here so far -- sits below A7 and therefore inherits at
# least its error. The available options are to take L4 dense and pay for
# it, or to quote A7's error alongside the result.
#
# **Still open.** Whether L4 is affordable at ``M\sim50`` (the Gaunt
# coefficients are geometry and precompute once; only the radial integrals
# of ``\eta_{l'}`` change per step); the stability cost of discarding the
# ``m=2`` parametric channel; and the eigenvalue problem for the
# variable-``\eta`` radial operator of A7.
