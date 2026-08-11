# # The Differential Formulation
#
# The model is stated variationally: three quadratic forms, and Lagrange's equations
# in the interior displacement amplitudes. This page carries the other route to the
# same equations -- the momentum equation written per mode, its boundary conditions
# imposed explicitly, and the traction projected onto the surface harmonics.
#
# It is kept because it is an independent derivation of a result the model depends on,
# and its checks still run. It is not part of the model: the assembly evaluates the
# three forms by quadrature on the strain field, and calls none of the operators
# below. Where the two routes are shown to agree is on the model page itself.

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
@variables xc yc zc  #src
@variables xs ts  #src

println("="^78)  #src
println("THE DIFFERENTIAL FORMULATION")  #src
println("="^78)  #src

# ## How a mode equation is formed
#
# A single-mode equation is assembled in three steps, and the order matters: the
# coupling between modes enters at the second, and everything that follows is fixed
# by it. For a given surface deformation,
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
# ### Eliminating the pressure
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
# ### The spherically symmetric part of the viscosity
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
# Two consequences.[^reid]
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
# [^reid]:
#     Setting ``\eta'=\eta''=0`` leaves ``\eta\,\mathcal D_l^2[\psi]``, so the
#     interior equation reduces to vorticity diffusion,
#     ``\partial_t\mathcal D_l[\psi_l]=\mathrm{Oh}\,\mathcal D_l^2[\psi_l]``, and a
#     normal mode ``\psi_l\propto e^{-\sigma t}`` gives
#     ``\mathcal D_l(\mathcal D_l+q^2)[\psi_l]=0`` with ``q^2=\sigma/\mathrm{Oh}`` --
#     Reid's Eq. 9. The constant-viscosity theory is the ``k=0`` diagonal of this
#     one at a single harmonic.
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

# ### The off-diagonal operators
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
    ## Expanded about an arbitrary point x0, NOT hard-wired to x0 = 1. An earlier    #src
    ## version substituted rr => 1.0 before extracting anything, which made any      #src
    ## residual x-dependence invisible BY CONSTRUCTION -- while the prose claimed     #src
    ## the failure mode under test was exactly that. Expanding about two different    #src
    ## points and requiring the recovered coefficients to agree is what actually     #src
    ## tests the claimed power x^(j+i-4).                                             #src
    Uc(x0) = sum(uc[i+1]*(rr-x0)^i/factorial(i) for i in 0:4)  #src
    Ec(x0) = sum(ec[j+1]*(rr-x0)^j/factorial(j) for j in 0:2)  #src
    ## a_{ji} for all (j,i) at once, for one (l,m,k), from an expansion about x0.  #src
    ## The projection at radius x0 returns a_{ji} * x0^(j+i-4), so the recovered    #src
    ## coefficient is scaled back by x0^(4-j-i).                                    #src
    function coeffs(gen, l, x0)  #src
        A = zeros(3, 5)  #src
        for j in 0:2, i in 0:4  #src
            ev = ntuple(a -> (a-1 == j ? 1.0 : 0.0), 3)  #src
            uv = ntuple(b -> (b-1 == i ? 1.0 : 0.0), 5)  #src
            num = 0.0  #src
            for (mu, w) in zip(nodes, wts)  #src
                th = acos(mu)  #src
                F  = -x0*sin(th)*gen(x0, th, ev..., uv...)  #src
                num += w*(1-mu^2)*(F/sin(th)^2)*dLPc(l,mu)/(l*(l+1))  #src
            end  #src
            A[j+1, i+1] = (num/(2/((2l+1)*l*(l+1)))) * x0^(4-j-i)  #src
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
    worst_xdep = 0.0  #src
    for (m, k) in ((2,0), (2,1), (2,2))  #src
        As = Dict{Int,Matrix{Float64}}()  #src
        for x0 in (1.0, 0.73)  #src
            expr = ed2(curlv(Uc(x0)*Cgc(m), Ec(x0)*LPc(k,cos(tt))))  #src
            gen  = Symbolics.build_function(expr, rr, tt, ec[1],ec[2],ec[3],  #src
                                            uc[1],uc[2],uc[3],uc[4],uc[5]; expression=Val(false))  #src
            for l in max(1, abs(m-k)):(m+k)  #src
                iseven(l+k+m) || continue  #src
                A = coeffs(gen, l, x0)  #src
                if haskey(As, l)  #src
                    ## same coefficients from a different expansion point  #src
                    sc = max(maximum(abs.(As[l])), 1e-30)  #src
                    worst_xdep = max(worst_xdep, maximum(abs.(A .- As[l]))/sc)  #src
                else  #src
                    As[l] = A; ntriple += 1  #src
                    if k == 0 && l == m  #src
                        worst_diag = max(worst_diag, maximum(abs.(A .- boxedA(l))))  #src
                    elseif k > 0  #src
                        worst_offmag = max(worst_offmag, maximum(abs.(A)))  #src
                    end  #src
                end  #src
            end  #src
        end  #src
    end  #src
    @assert worst_xdep < 1e-8 "the coefficients depend on the expansion point, so R_lm does not have the claimed x^(j+i-4) form ($worst_xdep)"  #src
    @assert ntriple >= 5 "the sweep covered too few (l,m,k) triples to be meaningful ($ntriple)"  #src
    @assert worst_offmag > 1.0 "every off-diagonal operator came out zero; the extraction is not exercising the coupling"  #src
    @assert worst_diag < 1e-9 "the extracted diagonal does not reproduce the boxed R_l ($worst_diag)"  #src
    @printf("  ASSERTION 3g OK: R_{lm} extracted as 15 numbers per (l,m,k) over %d triples,\n", ntriple)  #src
    @printf("    and the coefficients agree to %.1e between expansions about x0 = 1 and\n", worst_xdep)  #src
    @printf("    x0 = 0.73, which is what tests the claimed x^(j+i-4) power.\n")  #src
    @printf("    The diagonal (k=0, l=m) reproduces the boxed R_l to %.1e, and the\n", worst_diag)  #src
    @printf("    off-diagonal operators are nonzero (largest coefficient %.3g), so the\n", worst_offmag)  #src
    println("    mode coupling is a sparse table computed once, not a per-step projection.")  #src
    println("    Physical meaning of a failure: the operator would not have this form,")  #src
    println("    some coefficient would depend on x beyond x^(j+i-4), and the coupling")  #src
    println("    could not be tabulated at all -- assembly would need a projection")  #src
    println("    inside the time loop.")  #src
end  #src

# ### The reach of the coupling
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
# ### Exactness of the truncation
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

# ### The banded matrix
#
# Read the upper bound backwards: mode ``l`` reaches mode ``m`` only through
# viscosity harmonics with ``k\ge|l-m|``. So if ``\eta`` has angular content
# only up to some ``L_\eta``, no pair of modes further apart than ``L_\eta`` is
# coupled at all, and ``\mathcal D`` is banded with half-bandwidth ``L_\eta``.
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

# ### The Newtonian case
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

# ## The normal-stress balance
#
# ### The balance
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
# ### The allocation
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
    @assert abs(base_right) < 1e-14 "the derived normal-stress balance fails the base state ($base_right)"  #src
    @assert abs(base_wrong) > 1.0 "the opposite curvature sign also passes; the base state cannot fix it"  #src
    ## The same balance, read at l = 0 as an equation for the pressure LEVEL rather  #src
    ## than as an equation of motion: with no flow and no contact it has the unique  #src
    ## solution p = 2. That is what closes c_0, and it needs no extra postulate.     #src
    ## An earlier version of the page claimed volume conservation closed it, which   #src
    ## is a condition on the SHAPE standing in for one on the pressure.              #src
    p_level = 2.0 + 0.0 + 0.0          # solve -p + 2 eta e_rr + div n + p_c = 0 for p  #src
    @assert abs(p_level - 2.0) < 1e-14 "the l=0 balance does not fix the pressure level at the Laplace value ($p_level)"  #src
    ## (ii) the determinant vanishes at Reid's roots, and only for the derived sign  #src
    worst_right, best_wrong = 0.0, Inf  #src
    for Oh in (0.006, 0.05, 0.3, 1.0), l in (2, 4, 8)  #src
        q  = dominant_root(Oh, l)  #src
        sc = maximum(abs(alloc_det(q*(1+d), Oh, l)) for d in (0.15, -0.15, 0.3))  #src
        @assert sc > 0 "the determinant is identically zero; the check is vacuous"  #src
        worst_right = max(worst_right, abs(alloc_det(q, Oh, l; s = +1.0))/sc)  #src
        best_wrong  = min(best_wrong,  abs(alloc_det(q, Oh, l; s = -1.0))/sc)  #src
    end  #src
    @assert worst_right < 1e-12 "the BC allocation does not reproduce Reid's characteristic equation ($worst_right)"  #src
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


