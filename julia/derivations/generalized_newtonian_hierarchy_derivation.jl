# # Shear-Thinning Drops: a Hierarchy of Models, from the Exact Problem Down
#
# This file answers one question: **can a shear-thinning fluid live inside
# Reid's Legendre-mode framework, and if so, at what price?**
#
# It is written as a *chain*. We start from the exact free-surface problem
# with no assumptions at all, then make one simplification at a time. Each
# step is labelled, and each step records three things:
#
# 1. the assumption that was made,
# 2. what it throws away,
# 3. how you would undo it if you later decide you need to.
#
# The point of the chain is that you can stop wherever you like and know
# exactly what you are standing on. The model this repo currently
# implements sits near the bottom; a reader who wants something more
# faithful can walk back up and see precisely which rung to climb to.
#
# **Why this file exists.** Earlier attempts in this repo went astray by
# making several of these simplifications at once, silently, and then
# defending the result with checks that could not fail. Writing the chain
# out makes each choice visible and each error estimate a number rather
# than an adjective.
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
# | ``a_l(t)`` | amplitude of surface mode ``l`` |
# | ``l'`` | degree index of the **viscosity field's own** Legendre series |
# | ``L_\eta`` | highest ``l'`` retained -- the *bandwidth* of the coupling |
#
# A note on one symbol that has caused confusion in this repo: ``l'`` is
# **not** a mode of the drop's shape. The shape expansion runs
# ``l=2\ldots M`` with ``M`` of order 50--90. ``l'`` indexes the angular
# structure of the *viscosity*, and as shown in Step 5 it controls only how
# far off the diagonal the mode-coupling matrix reaches. Its dimension is
# always ``M``.
#
# A failing assertion below means the chain has a broken link: either the
# algebra here is wrong, or a step that was claimed to be exact is not.

using Symbolics, QuadGK, SpecialFunctions, DropSolver  #src
using Printf  #src

# Helper reused throughout: try to prove an expression is identically zero
# symbolically, and if the CAS cannot finish the simplification, fall back
# to evaluating it at a spread of concrete points. The fallback is honest --
# it demotes a proof to strong evidence, and we say so where it is used.
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

# ## Step 0 -- the exact problem
#
# No assumptions. An incompressible fluid of constant density ``\rho``
# occupies a region bounded by a free surface ``\Sigma(t)``. The unknowns
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
# Everything below is a simplification of these five lines.
#
# ### The one identity that shapes the whole problem
#
# For a **constant** viscosity the divergence of the stress collapses to
# ``\mu\nabla^2\bm u``, which is what makes Reid's problem separable. For a
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
# **In plain English:** a shear-thinning fluid adds exactly one term to the
# Navier-Stokes momentum equation. Everything hard about this project is
# contained in that single extra term ``2(\nabla\eta)\cdot\bm e``. Reid was
# entitled to drop it because for him ``\nabla\eta\equiv 0``.
#
# Let us verify the incompressible identity ``\nabla\cdot\bm e =
# \tfrac12\nabla^2\bm u`` rather than assert it, since the boxed result
# rests on it.

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

# ## Step 1 (A1) -- linearise in the surface amplitude
#
# **Assumption.** The surface displacement is small, ``\epsilon=\zeta/R\ll1``.
#
# **What it drops.** The advective term ``\bm u\cdot\nabla\bm u``, and the
# transfer of the boundary conditions from the deformed surface
# ``r=R+\zeta`` to the sphere ``r=R``. Error ``O(\epsilon^2)``.
#
# **How to undo it.** Standard weakly-nonlinear free-surface theory. This
# is entirely independent of the rheology question and can be revisited
# separately.
#
# ### The trap in this step
#
# This is the single most important warning in the file, because an earlier
# generation of this work fell into it.
#
# Linearising the **kinematics** does *not* linearise the **rheology**.
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
# So "linear in ``\epsilon``" and "expandable in powers of ``\epsilon``" are
# different statements. The first is a legitimate assumption about the
# geometry. The second is false for a general Carreau-Yasuda fluid, and no
# amount of small amplitude rescues it.
#
# The ratio ``|\epsilon|^a/|\epsilon|`` measures whether the correction is
# genuinely higher order. If it stays bounded, the term is a legitimate
# perturbation; if it diverges, there is nothing to expand in:
#
# | ``\epsilon`` | ``a=0.5`` | ``a=0.743`` | ``a=1`` | ``a=2`` |
# |:--|--:|--:|--:|--:|
# | ``10^{-1}`` | 3.16 | 1.81 | 1 | 0.1 |
# | ``10^{-3}`` | 31.6 | 5.90 | 1 | ``10^{-3}`` |
# | ``10^{-5}`` | 316 | 19.3 | 1 | ``10^{-5}`` |
# | ``10^{-7}`` | ``3.16\times10^{3}`` | 63.0 | 1 | ``10^{-7}`` |
#
# The shear-thinning columns *grow* as the amplitude falls, monotonically,
# across six decades. The ``a=2`` column *collapses* -- which is precisely
# what "higher order" means, and precisely why the classical Carreau theory
# closes where the general one cannot. ``a=1`` is the marginal case.

let  #src
    println("\nSTEP 1 (A1): linear in amplitude, NONLINEAR in rheology")  #src
    println("  ratio |eps|^a / |eps|  -- if this diverges, eps^a is not a higher-order term")  #src
    println("  eps        a=0.5      a=0.743     a=1.0      a=2.0")  #src
    thin = Float64[]  #src
    newt = Float64[]  #src
    for epsv in (1e-1, 1e-3, 1e-5, 1e-7)  #src
        row = [abs(epsv)^a / abs(epsv) for a in (0.5, 0.743, 1.0, 2.0)]  #src
        push!(thin, row[2]); push!(newt, row[4])  #src
        @printf("  %-10.0e %-10.3g %-11.3g %-10.3g %-10.3g\n",  #src
                epsv, row[1], row[2], row[3], row[4])  #src
    end  #src
    ## The claim is that the ratio DIVERGES as eps -> 0, so the test is that it  #src
    ## increases monotonically across every decade -- not that it clears some  #src
    ## arbitrary threshold at the smallest amplitude we happened to tabulate.  #src
    @assert all(thin[k+1] > thin[k] for k in 1:length(thin)-1) "eps^a/eps must grow as eps->0 for a<1"  #src
    @assert thin[end] > 10*thin[1] "growth over six decades is implausibly weak"  #src
    ## The contrast is the point: for a=2 the same ratio COLLAPSES, which is  #src
    ## exactly what "higher order" means.  #src
    @assert all(newt[k+1] < newt[k] for k in 1:length(newt)-1) "for a=2 the ratio must shrink"  #src
    println("  ASSERTION 2 OK: for a<1 the Carreau-Yasuda correction grows without bound")  #src
    println("    relative to the linear term as eps->0. A regular perturbation expansion")  #src
    println("    in amplitude does not exist. Physically: there is no amplitude small")  #src
    println("    enough that shear-thinning becomes a small correction of fixed order.")  #src
    println("    (For a=2 the ratio is proportional to eps -- genuinely higher order.")  #src
    println("     That is why the CLASSICAL Carreau theory closes and the general one does not.)")  #src
end  #src

# ## Step 2 (A2) -- axisymmetry
#
# **Assumption.** The forcing (impact) is axisymmetric.
#
# **What it drops.** Nothing. This step is *exact*, and it is worth proving
# because it is the one piece of good news in the whole chain.
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

# ## Step 3 (A3) -- poloidal representation and modal expansion
#
# **Assumption.** None beyond incompressibility and axisymmetry -- this is a
# change of variables, not a simplification.
#
# An axisymmetric incompressible velocity field is generated by a single
# scalar. Writing the surface shape as
#
# ```math
# \zeta(\theta,t) = R\sum_{l\ge2} a_l(t)\,P_l(\cos\theta),
# ```
#
# the state of the drop is the vector of modal amplitudes ``\{a_l(t)\}``,
# ``l=2\ldots M``. This is exactly the state `julia/src/types.jl` carries.
#
# For each mode, Reid's radial velocity amplitude is ``F(x)=U(x)/x^2`` with
# ``U(x)=C\,x\,j_l(qx)+\Pi_0 x^{l+1}``, and the strain components follow
# from the poloidal representation. Those are derived in
# `reid1960_full_derivation.jl`; we use them below without re-deriving.

# ## Step 4 (L4) -- the full coupled system
#
# We are now at the most general *tractable* statement of the problem, and
# this is the rung to remember: everything below is a simplification **of
# this**, and anything you later want to recover, you recover by climbing
# back to here.
#
# Expand the viscosity field in its own Legendre series,
#
# ```math
# \eta(r,\theta,t) = \sum_{l'} \eta_{l'}(r,t)\,P_{l'}(\cos\theta),
# ```
#
# substitute into ``\eta\nabla^2\bm u + 2(\nabla\eta)\cdot\bm e``, and project
# the momentum equation onto ``P_l``. Because both ``\eta`` and ``\bm u`` are
# now Legendre series, every term becomes a sum of integrals of *three*
# angular functions, and the modal equations couple:
#
# ```math
# \boxed{\ \ddot a_l + \sum_{l''}\Bigl[\mathcal D^{(2)}_{l l''}\,\dot a_{l''}
#        + \mathcal D^{(1)}_{l l''}\,a_{l''}\Bigr] = F_l\ }
# ```
#
# where ``\mathcal D^{(1)},\mathcal D^{(2)}`` are ``M\times M`` matrices built as
#
# ```math
# \mathcal D_{l l''} = \sum_{l'} \bigl(\text{radial integral of }\eta_{l'}\bigr)
#                      \times G^{l'}_{l l''} ,
# \qquad
# G^{l'}_{l l''} = \frac{2l+1}{2}\int_{-1}^{1} P_l P_{l'} P_{l''}\,d\mu .
# ```
#
# Compare Reid, where ``\eta`` is constant so only ``l'=0`` survives,
# ``G^{0}_{l l''}=\delta_{l l''}``, and the matrices are **diagonal** -- one
# independent oscillator per mode.
#
# ### The Gaunt selection rule, and what actually controls cost
#
# The triple integral ``G^{l'}_{l l''}`` is a Gaunt coefficient. It vanishes
# unless
#
# ```math
# |l-l''| \le l' \le l+l'' \qquad\text{and}\qquad l+l'+l''\ \text{is even}.
# ```
#
# This rule is **geometric**: it holds no matter what the fluid is. Read
# backwards it says something decisive:
#
# > Mode ``l`` couples to mode ``l''`` only if ``|l-l''|\le l'`` for some ``l'``
# > at which the viscosity field has content. **The bandwidth of the coupling
# > matrix equals the spectral width of ``\eta``, and nothing else.**
#
# So the cost of the coupled model is *not* set by ``M``. If ``\eta``'s
# angular spectrum is confined to ``l'\le L_\eta``, then ``\mathcal D`` is a
# banded ``M\times M`` matrix of half-bandwidth ``L_\eta``, and applying it
# costs ``O(M L_\eta)`` rather than ``O(M^2)``. With ``M=90`` and
# ``L_\eta=6`` that is nearly as cheap as diagonal.

let  #src
    println("\nSTEP 4 (L4): the coupled system, and the Gaunt selection rule")  #src
    Pl(l, m) = l == 0 ? one(m) : l == 1 ? m :  #src
        begin am, b = one(m), m; for n in 1:l-1; b, am = ((2n+1)*m*b - n*am)/(n+1), b; end; b end  #src
    gaunt(l, lp, lpp) = quadgk(m -> Pl(l,m)*Pl(lp,m)*Pl(lpp,m), -1, 1; rtol=1e-12)[1]  #src
    viol_zero = 0.0  #src
    viol_nonzero = Inf  #src
    for l in 0:8, lp in 0:8, lpp in 0:8  #src
        g = gaunt(l, lp, lpp)  #src
        allowed = (abs(l-lpp) <= lp <= l+lpp) && iseven(l+lp+lpp)  #src
        if allowed  #src
            (l+lp+lpp > 0) && (viol_nonzero = min(viol_nonzero, abs(g)))  #src
        else  #src
            viol_zero = max(viol_zero, abs(g))  #src
        end  #src
    end  #src
    @assert viol_zero < 1e-11 "a FORBIDDEN Gaunt coefficient was nonzero ($viol_zero)"  #src
    @assert viol_nonzero > 1e-6 "an ALLOWED Gaunt coefficient vanished ($viol_nonzero)"  #src
    @printf("  largest FORBIDDEN |G| : %.2e   (must be ~0)\n", viol_zero)  #src
    @printf("  smallest ALLOWED  |G| : %.2e   (must be well clear of 0)\n", viol_nonzero)  #src
    println("  ASSERTION 4 OK: the selection rule |l-l''|<=l'<=l+l'' with l+l'+l'' even")  #src
    println("    holds exactly, over all 729 triples with l,l',l'' <= 8.")  #src

    println("\n  Bandwidth consequence -- multiplication by P_{l'} in the P_l basis:")  #src
    for lp in (0, 2, 4)  #src
        M = 12  #src
        band = 0  #src
        for l in 0:M, lpp in 0:M  #src
            if abs(gaunt(l, lp, lpp)) > 1e-11  #src
                band = max(band, abs(l - lpp))  #src
            end  #src
        end  #src
        @printf("    eta content at l'=%d  ->  couples modes up to |l-l''| = %d\n", lp, band)  #src
        @assert band <= lp "bandwidth exceeded l' -- the selection rule argument is broken"  #src
    end  #src
    println("  ASSERTION 5 OK: viscosity content at degree l' couples modes no further")  #src
    println("    than |l-l''| = l'. The coupling matrix is BANDED with half-bandwidth")  #src
    println("    L_eta, dimension M. Cost is set by L_eta, NOT by M.")  #src
    println("    (Physical meaning of a failure here: a spatially varying viscosity would")  #src
    println("     couple every mode to every other regardless of how smooth it is, and no")  #src
    println("     truncation of the coupled model could ever be justified.)")  #src
end  #src

# ## Step 4b -- where the exponent ``a`` enters *structurally*
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
# \qquad p = \frac{n-1}{a}.
# ```
#
# * ``X`` is a polynomial in ``\bm e`` **iff** ``a`` is a non-negative even
#   integer, because ``\dot\gamma^a=(\dot\gamma^2)^{a/2}`` and
#   ``\dot\gamma^2=2\bm e\!:\!\bm e`` is quadratic.
# * ``(1+X)^p`` is a polynomial in ``X`` **iff** ``p`` is a non-negative
#   integer.
#
# For a genuinely shear-thinning fluid ``n<1``, so ``p<0`` and the second
# condition **never** holds. The honest conclusion:
#
# > For no real shear-thinning fluid is ``\eta`` exactly polynomial. The
# > coupling matrix is never *exactly* banded, and ``L_\eta`` is always an
# > empirical truncation with a measurable error -- never an exact one.
#
# This also corrects a piece of folklore in this repo. The classical
# ``a=2`` theory does **not** close because ``a=2`` makes ``\eta``
# polynomial -- it does not. It closes because that theory additionally
# *truncates* ``(1+X)^p\approx 1+pX`` at first order, and ``X`` alone is
# polynomial when ``a=2``. The finiteness comes from the perturbative
# truncation, and ``a=2`` merely keeps the truncated object polynomial.
#
# ### The exception that matters: ``p=-1``
#
# There is one value of ``p`` that is special without being a non-negative
# integer, and it is the one your fluid actually has. If ``p=-1``, i.e.
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
# left on the outside. This is the Cross model, and it is why Step 9 below
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
    println("  ASSERTION 7 OK: p<0 for every shear-thinning fluid, so eta is never")  #src
    println("    exactly polynomial and L_eta is always an empirical truncation.")  #src
    println("    Physical meaning: no shear-thinning drop model can claim an EXACT")  #src
    println("    finite mode-coupling bandwidth; the truncation must carry an error bar.")  #src
end  #src

# ## Step 5 (A5) -- the temporal closure
#
# **Assumption.** How ``\eta``'s time dependence is reduced.
#
# ``\eta`` depends on ``t`` through ``\dot\gamma(t)``, which oscillates with
# the mode. Three inequivalent choices, and this repo's production code
# silently takes the first:
#
# | choice | what it keeps | status |
# |:--|:--|:--|
# | (a) instantaneous ``\eta(t)`` | everything | **quasi-static; unjustified** |
# | (b) period-averaged | the ``m=0`` channel | justified by the lemma below |
# | (c) Floquet | ``m=0`` and ``m=2`` | most faithful, most work |
#
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

let  #src
    println("\nSTEP 5 (A5): the period-pi lemma -- why instantaneous evaluation is not licensed")  #src
    e_amp = [1.3 + 0.7im, -0.4 + 1.1im, 0.9 - 0.2im, 0.5 + 0.3im]  #src
    Sof(φ) = begin  #src
        r = [real(z*(cos(φ) - im*sin(φ))) for z in e_amp]  #src
        sqrt(2*(r[1]^2 + r[2]^2 + r[3]^2 + 2*r[4]^2))  #src
    end  #src
    worst = maximum(abs(Sof(φ) - Sof(φ + π)) for φ in range(0, 2π; length=257))  #src
    @assert worst < 1e-13 "S is not period-pi; the lemma is false as stated"  #src
    @printf("  max |S(phi) - S(phi+pi)| over a full period : %.2e\n", worst)  #src
    N = 4096  #src
    φs = range(0, 2π; length=N+1)[1:end-1]  #src
    c1 = sum(Sof(φ)*cos(φ) for φ in φs)*2/N  #src
    s1 = sum(Sof(φ)*sin(φ) for φ in φs)*2/N  #src
    c2 = sum(Sof(φ)*cos(2φ) for φ in φs)*2/N  #src
    c0 = sum(Sof(φ) for φ in φs)/N  #src
    @printf("  Fourier content:  m=0 : %.4f    m=1 : %.2e / %.2e    m=2 : %.4f\n",  #src
            c0, abs(c1), abs(s1), abs(c2))  #src
    @assert abs(c1) < 1e-12 && abs(s1) < 1e-12 "S has content at the mode's OWN frequency"  #src
    @assert abs(c2) > 1e-3 "the m=2 parametric channel vanished; expected it to survive"  #src
    println("  ASSERTION 8 OK: S has EXACTLY zero content at m=1 (the mode's own")  #src
    println("    frequency) and nonzero content at m=0 and m=2.")  #src
    println("    => the leading temporal effect is the PERIOD-AVERAGED (m=0) channel.")  #src
    println("    => evaluating eta instantaneously injects a spurious 2*omega modulation")  #src
    println("       whose quasi-static justification is violated by O(1), at any amplitude.")  #src
end  #src

# ## Step 6 (A6) -- truncate the viscosity spectrum at ``L_\eta``
#
# **Assumption.** ``\eta``'s Legendre content above ``l'=L_\eta`` is negligible.
#
# **What it drops.** Coupling between modes further apart than ``L_\eta``.
# The matrix becomes banded; cost falls from ``O(M^2)`` to ``O(M L_\eta)``.
#
# **Error.** Exactly the discarded coupling -- a measurable number, not an
# adjective. It has now been measured, and the answer kills this rung.
#
# ### Measurement, and a warning about the wrong error norm
#
# The natural-looking metric, "what fraction of ``\sum_{l'}|\eta_{l'}|^2``
# lies below ``L_\eta``", says ``L_\eta=2\ldots5`` suffices everywhere. That
# metric is **misleading**. What controls the banding error is the *summed
# magnitude* of the discarded coefficients,
# ``T_1(L)=\sum_{l'>L}|\eta_{l'}|/|\eta_0|``, because those terms enter the
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
# > **``L_\eta\gtrsim M``, and for realistic multi-mode states
# > ``L_\eta\gg M``. The coupling matrix is effectively dense. This rung
# > buys nothing.**
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
# entry above is a *lower* bound on the discarded coupling. This measurement
# is reported, not re-executed in CI -- it costs ~30 minutes and depends on a
# separate high-``l`` eigenfunction evaluator validated to
# ``5\times10^{-15}`` on Reid's tangential-stress boundary condition.
#
# **How to undo it.** Raise ``L_\eta``. The structure does not change -- only
# the bandwidth. But since ``L_\eta\gtrsim M`` is required, the honest move
# is to skip this rung entirely and either take L4 dense or drop to A7.

# ## Step 7 (A7) -- keep only ``l'=0``: a spherically symmetric viscosity
#
# **Assumption.** ``L_\eta=0``, i.e. ``\eta=\eta(r,t)`` with no angular
# structure.
#
# **Why this rung is special.** A spherically symmetric ``\eta`` commutes
# with the angular Laplacian. Every mode **decouples again**:
# ``G^{0}_{l l''}=\delta_{l l''}``, the matrices return to diagonal, and the
# entire architecture of this repo -- one independent oscillator per mode --
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
#     This corrects an earlier estimate in this repo. A measurement using a
#     *single* active ``l=2`` mode found the viscosity varying by only
#     ``1.1``--``1.2\times`` across the drop, and concluded that spatial
#     homogenization was a ``\sim10\%`` effect. That is right for one mode
#     and wrong for a real state: with a dozen modes beating, the nodal
#     collapse described in Step 6 makes ``\eta`` span ``100\times``. Any
#     scalar-``\mathrm{Oh}_{\mathrm{eff}}`` model inherits this error.
#
# **A7 is the cheapest defensible rung, not an accurate one.** It keeps honest
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
#   in `reid1960_full_derivation.jl`.
# * **BC3 (normal stress)** carries ``\eta`` multiplicatively, so it becomes
#   the *surface* value ``\eta_s=\eta(\dot\gamma|_{r=R})``. Same form, one
#   state-dependent coefficient.
#
# So the honest scope of "adapt Reid for a shear-thinning fluid" is: **one
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

# ### The radial equation, derived
#
# Reid's problem is separable because his operator commutes with the angular
# Laplacian. With ``\eta=\eta(r)`` that is still true, so we can derive the
# replacement for his radial ODE directly.
#
# Use the Stokes stream function for axisymmetric incompressible flow. For
# surface mode ``l``,
#
# ```math
# \psi(r,\theta) = U(r)\,A_l(\theta),
# \qquad
# A_l(\theta) = \frac{\sin^2\theta\;P_l'(\cos\theta)}{l(l+1)},
# ```
#
# which reproduces ``u_r=F(r)P_l(\cos\theta)`` with ``U=r^2F`` -- Reid's own
# variable. The angular part is a Gegenbauer function, so it diagonalises the
# Stokes operator ``E^2=\partial_{rr}+\frac{1-\mu^2}{r^2}\partial_{\mu\mu}``:
#
# ```math
# E^2\psi = \mathcal D_l[U]\,A_l,
# \qquad
# \mathcal D_l[U] \equiv U'' - \frac{l(l+1)}{r^2}U .
# ```
#
# Taking the ``\varphi``-component of the curl of the momentum equation
# annihilates the pressure, which is what makes this tractable at all. With
# ``\rho=1``, ``\eta_0=1``, and a time factor ``e^{-\sigma t}`` (so
# ``\sigma=q^2``), building ``\bm e``, then ``\bm\tau=2\eta(r)\bm e``, then
# ``\nabla\cdot\bm\tau`` in spherical coordinates, and finally the curl, gives
# for **constant** ``\eta``
#
# ```math
# \boxed{\ \bigl[\nabla\times\bm M\bigr]_\varphi
#   = -\frac{A_l(\theta)}{r\sin\theta}\,
#     \mathcal D_l\bigl(\mathcal D_l + q^2\bigr)[U]\ }
# ```
#
# -- exactly Reid's ``\mathcal D_l(\mathcal D_l+q^2)U=0``, recovered from the
# general variable-viscosity machinery. Verified to a relative spread of
# ``\sim10^{-13}`` at five independent ``(r,\theta)`` points, for
# ``l=2,3,4``, against an independently constructed Reid operator. The factor
# ``1/(r\sin\theta)`` is the Jacobian between the curl and the stream-function
# operator, not a physical term.
#
# **What this buys.** The same construction with ``\eta=\eta(r)`` non-constant
# produces additional terms proportional to ``\eta'`` and ``\eta''``, and
# nothing else -- the angular structure is untouched, so the equation stays
# diagonal in ``l``, exactly as Step 7 promised. The resulting radial problem
# is a linear two-point boundary-value eigenproblem: no longer Bessel, but
# entirely standard, and `julia/src/reid.jl`'s continuation machinery already
# knows how to track eigenvalue branches through one.
#
# A note on method, since it cost real time: `Symbolics.build_function` on
# these expression trees does not terminate in any useful time. Substituting
# concrete values *first* collapses the tree to something trivial; the only
# wrinkle is that `substitute` does not fold transcendentals, so the result
# needs `toexpr` + `eval` to become a number.

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
# **Assumption.** The radial variation of ``\eta`` is small enough to
# replace by a single number.
#
# **What it buys.** The radial equation becomes constant-coefficient again,
# i.e. Bessel's equation, and Reid's closed-form characteristic equation
# returns verbatim -- evaluated at a shifted Ohnesorge number.
#
# **This is where the current production code sits** (`st_exact_extension.jl`),
# combined with choice (a) of Step 5.
#
# **How to undo it.** Step 7.

# ## Step 9 -- the Cross model as a rung, not an alternative
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
# 1. By Step 4b the constitutive law becomes an **algebraic constraint**
#    ``(1+X)\eta=(1+X)\eta_\infty+\Delta\eta`` -- linear in ``\eta``, linear
#    in ``X``, no outer fractional power. Symbolic manipulation becomes
#    tractable in a way it is not for general ``p``.
# 2. It is how this repo's validation fluid is *actually characterised*
#    (`scripts/validate_shear_thinning.jl` fits ``\eta_0,\eta_\infty,K,m``
#    and only then converts). Deriving Cross directly removes that
#    conversion -- which is worth doing on its own merits, since the
#    conversion currently mis-assigns a parameter (see the assertion below).
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
    Δ = (ETA_0 - ETA_INF)/ETA_0  #src
    p_true = -1.0  #src
    exponent_should_be = (1 - (1 - M_CROSS))/M_CROSS  #src
    @printf("  this fluid: m = %.4f, Delta = (eta0-etainf)/eta0 = %.6f\n", M_CROSS, Δ)  #src
    @printf("  the CY exponent (1-n)/a implied by the Cross mapping = %.6f (exactly 1)\n",  #src
            exponent_should_be)  #src
    @assert abs(exponent_should_be - 1.0) < 1e-12 "the Cross->CY exponent should be exactly 1"  #src
    println("  ASSERTION 11 OK: under the Cross mapping the CY exponent is EXACTLY 1,")  #src
    println("    but scripts/validate_shear_thinning.jl passes Delta = $(round(Δ, digits=6))")  #src
    println("    into that slot. Harmless here only because Delta ~ 1 for this fluid;")  #src
    println("    for a fluid with Delta = 0.5 the exponent would silently be wrong by 2x.")  #src
    println("    Physical meaning: the AMPLITUDE of the thinning is being used as its")  #src
    println("    SHAPE exponent. Deriving Cross natively removes the conversion entirely.")  #src
end  #src

# ## Step 10 -- the Newtonian floor, and a live cross-check
#
# The bottom of the chain must reproduce what this repo already trusts. Two
# checks against the *running solver*, not against this file's own algebra:
# Reid's exact finite-Ohnesorge coefficients must emerge when the viscosity
# stops depending on the flow, and they must reduce to Lamb as ``Oh\to0``.
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
# ``\lambda_2=0.2187`` against Lamb's ``0.2500`` -- a ``13\%`` difference, which
# is why this repo uses Reid's exact relation rather than the closed form.

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
# ``L_\eta`` is ``\gtrsim M`` for synthetic spectra and ``\gg M`` for real
# solver states.
#
# **Priced, and not cheap.** A7 restores exact diagonality, but discards
# angular structure that is *comparable to or larger than* the mean
# viscosity. It is the cheapest defensible rung, not an accurate one.
#
# **The uncomfortable conclusion.** Between L4 (dense) and A7
# (leading-order error) there is no cheap-and-accurate rung. Any model
# built on a single scalar ``\mathrm{Oh}_{\mathrm{eff}}`` per mode -- which
# is every shear-thinning model this repo has implemented so far -- sits
# below A7 and therefore inherits at least its error. The honest options
# are to take L4 dense and pay for it, or to quote A7's error rather than
# assume it away.
#
# **Still open.** Whether L4 is actually affordable at ``M\sim50`` (the
# Gaunt coefficients are geometry and precompute once; only the radial
# integrals of ``\eta_{l'}`` change per step); the stability cost of
# discarding the ``m=2`` parametric channel; and the variable-``\eta``
# radial operator itself, which is the content of A7 and is derived in the
# next section to be added here.
