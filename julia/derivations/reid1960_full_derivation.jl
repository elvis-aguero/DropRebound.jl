# # Oscillations of a Viscous Liquid Drop
#
# A from-scratch, CAS-verified derivation of W. H. Reid's (1960) exact
# characteristic equation for the free oscillations of a viscous liquid
# drop -- the physics underneath every rheology model in DropSolver. The
# other derivations here (Oldroyd-B, Carreau-Yasuda) are corrections on top
# of what follows, so this page is the one to read first.
#
# **What problem are we solving?** Take a small liquid drop (water in air,
# say) that has been slightly deformed from its equilibrium spherical
# shape. Surface tension pulls it back; viscosity dissipates energy as it
# does. Two questions: *how fast does it oscillate, and how fast does that
# oscillation decay?* Answering both simultaneously, for **any** viscosity
# (not just the small- or large-viscosity limits), is what this derivation
# does.
#
# Three distinct pieces of physics have to be resolved at once, and this is
# why the problem is harder than it first looks:
#
# 1. the inviscid capillary oscillation modes, characterized by a single
#    frequency ``\sigma_{l;0}`` for each spherical-harmonic order ``l``;
# 2. viscous dissipation throughout the bulk of the drop;
# 3. the coupling between surface deformation, internal flow, and pressure,
#    all of which must be mutually consistent *at the boundary*.
#
# **Why does this matter for a numerical solver?** `julia/src/timestepper.jl`
# needs, for every surface mode ``l``, a damping rate and an oscillation
# frequency to put into that mode's equation of motion. Those two numbers
# are exactly the two roots of the equation derived here.
#
# ## Where this sits historically
#
# Three analytical milestones precede Reid's result:
#
# - **Rayleigh (1879) / Lamb (1932)** -- the inviscid frequencies
#   ``\sigma_{l;0}`` (Section 1 below).
# - **Lamb (1881)** -- the *small*-viscosity correction,
#   ``\sigma_{l;\nu} = (l-1)(2l+1)\nu/R^2 \pm i\,\sigma_{l;0}``. This is the
#   formula the production timestepper implements (Section 8), and the one
#   the companion page *Finite-Ohnesorge Coefficients* generalizes.
# - **Chandrasekhar (1959)** -- the full, arbitrary-viscosity solution, but
#   for a *self-gravitating* liquid globe rather than a surface-tension-held
#   drop.
#
# Reid (1960)'s contribution is to show that the surface-tension problem is
# *the same problem* as Chandrasekhar's. The insight is a rescaling: surface
# tension enters the entire calculation only through the single number
# ``\sigma_{l;0}``, so once the problem is non-dimensionalized the
# characteristic equation no longer remembers which restoring force produced
# it. A problem about surface tension turns out to be, mathematically,
# a problem about gravity. Section 7 is where that cancellation happens
# explicitly, and it is worth watching for.
#
# Downstream, Reid's equation is what Molaček & Bush (2012) compress into
# their per-mode coefficients ``A_m`` and ``D_m`` -- see Section 10 for that
# mapping, and *Finite-Ohnesorge Coefficients* for the parametrization
# DropSolver uses.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``R`` | equilibrium (undeformed) drop radius |
# | ``\rho`` | fluid density |
# | ``\mu = \rho\nu`` | dynamic / kinematic viscosity |
# | ``T_1`` | surface tension |
# | ``x = r/R`` | dimensionless radial coordinate |
# | ``\theta,\varphi`` | polar, azimuthal angle |
# | ``l`` | spherical-harmonic degree of the surface deformation |
# | ``Y_l^m(\theta,\varphi)`` | spherical harmonic (angular shape of the deformation) |
# | ``\epsilon(t) = \epsilon_0 e^{-\sigma t}`` | dimensionless amplitude of the deformation |
# | ``\sigma`` | complex decay rate; ``\mathrm{Re}(\sigma)>0`` is decay, ``\mathrm{Im}(\sigma)`` is the oscillation frequency |
# | ``\sigma_{l;0}`` | inviscid oscillation frequency of mode ``l`` (real, no viscosity) |
# | ``u_r(r,\theta,\varphi,t)`` | radial fluid velocity |
# | ``U(x)`` | dimensionless radial-velocity eigenfunction (all of the unknown physics lives in this one function) |
# | ``q^2 = \sigma R^2/\nu`` | viscous wavenumber -- ``\sigma`` measured in units of the viscous diffusion rate ``\nu/R^2`` |
# | ``\alpha^2 = \sigma_{l;0}R^2/\nu`` | inviscid frequency, in the same units |
# | ``j_l(z)`` | spherical Bessel function of the first kind, order ``l`` |
# | ``Q_{l+1/2}(q) = j_{l+1}(q)/j_l(q) = J_{l+3/2}(q)/J_{l+1/2}(q)`` | the one Bessel-function combination the whole problem reduces to |
#
using Symbolics            #src
using SpecialFunctions     #src
using DropSolver           #src

function symbolic_zero(expr)                                                     #src
    simplified = simplify(expr; expand=true)                                     #src
    is_symbolic_zero = isequal(simplified, 0) || isequal(simplified, 0.0)         #src
    vars = Symbolics.get_variables(expr)                                         #src
    is_numeric_zero = if isempty(vars)                                           #src
        ## A variable-free expression must be checked, not waved through: an     #src
        ## unconditional `true` here accepts any constant, so an assertion like  #src
        ## `symbolic_zero(bc1_solution - (-1))` would pass even for the wrong    #src
        ## sign. Evaluate it and compare against zero.                           #src
        v = Symbolics.value(simplified)                                          #src
        v isa Number && abs(v) < 1e-8                                            #src
    else                                                                         #src
        f = Symbolics.build_function(expr, vars...; expression=false)             #src
        test_vals = (0.37, 1.21, 2.03, 0.68, 1.59, 3.14, 0.91, 2.77)              #src
        all(abs(f((test_vals[mod1(i + k, length(test_vals))] for k in 1:length(vars))...)) < 1e-8  #src
            for i in 1:length(test_vals))                                        #src
    end                                                                          #src
    is_symbolic_zero || is_numeric_zero                                          #src
end                                                                              #src

# ## Section 1: The physical setup
#
# A drop of density ``\rho``, viscosity ``\mu=\rho\nu``, and surface tension
# ``T_1`` sits in equilibrium as a sphere of radius ``R``. With the external
# pressure set to zero, the Young-Laplace equation (the pressure jump across
# an interface with principal curvatures ``1/R_1+1/R_2`` is
# ``T_1(1/R_1+1/R_2)``, here ``R_1=R_2=R``) fixes the internal equilibrium
# pressure:
# ```math
# p_0 = \frac{2T_1}{R}.
# ```
#
# We perturb the surface with a single spherical-harmonic mode:
# ```math
# r = R\left[1 + \epsilon(t)\, Y_l^m(\theta,\varphi)\right], \qquad
# \epsilon(t) = \epsilon_0 e^{-\sigma t}, \qquad \epsilon_0 \ll 1,
# ```
# and look for the complex decay rate ``\sigma`` (the sign convention is
# ``\mathrm{Re}(\sigma)>0`` for decay -- Reid follows Chandrasekhar's
# convention here, not the more common ``e^{+i\omega t}``, specifically so
# that a damped oscillation is ``\sigma=\gamma+i\omega_d`` with both
# ``\gamma,\omega_d>0`` real and positive; in the inviscid limit
# ``\gamma\to0`` and ``\sigma\to\pm i\omega_d``).
#
# Here ``Y_l^m(\theta,\varphi)=P_l^m(\cos\theta)e^{-im\varphi}`` in the
# unnormalized convention. The normalization never matters: every equation
# below is linear and homogeneous in ``Y_l^m``, so the constant divides out,
# and the characteristic equation this page builds toward is independent of
# the choice.
#
# In the complete absence of viscosity, Rayleigh and Lamb's classical result
# gives the oscillation frequency of a spherical-harmonic mode of order
# ``l`` as
# ```math
# \sigma_{l;0}^2 = l(l-1)(l+2)\,\frac{T_1}{\rho R^3}.
# ```
# This is real (undamped) and independent of the azimuthal order ``m``. The
# ``l=2`` case is the oblate-prolate wobble familiar from a falling raindrop
# or a drop shaken loose from a tap. We will not re-derive this particular
# formula (it requires the inviscid energy functional, a separate
# calculation); we take it as given, exactly as Reid (1960) does, and use it
# only as a bookkeeping device: viscosity will enter the final answer only
# through this one number, ``\sigma_{l;0}``.
#
# Two sanity checks on the formula itself, before building anything on it.
# First, dimensions: it should be a rate squared, and
# ``[T_1/(\rho R^3)] = [\mathrm{N/m}]\,/\,[\mathrm{kg\,m^{-3}}\cdot\mathrm{m^3}]
# = [\mathrm{kg\,s^{-2}}]/[\mathrm{kg}] = [\mathrm{s^{-2}}]`` -- consistent.
# Second, the ``l=1`` mode is a rigid translation of the whole drop, which
# has no restoring force (translating a sphere changes neither its shape nor
# its surface energy), so ``\sigma_{1;0}^2`` must vanish identically. It
# does.

@variables l_sym T1 rho R                                                #src
sigma_l0_sq = l_sym * (l_sym - 1) * (l_sym + 2) * T1 / (rho * R^3)       #src
@assert symbolic_zero(substitute(sigma_l0_sq, Dict(l_sym => 1)))         #src

# ## Section 2: Mathematical preliminaries
#
# Two facts from mathematical physics do essentially all of the geometric
# work in what follows, and a third -- a structural statement about
# divergence-free vector fields -- is what lets the whole velocity field be
# described by one scalar function. None of the three is specific to viscous
# flow; they would show up in any spherically-symmetric wave or diffusion
# problem. They are isolated here, with the parts that are easy to get wrong
# derived in full and the standard, citable parts checked numerically.
#
# ### 2.1 Spherical harmonics are eigenfunctions of the angular Laplacian
#
# For axisymmetric problems (which is all we need: the deformed drop has no
# preferred azimuthal direction), the spherical harmonic reduces to a
# Legendre polynomial, ``Y_l^0(\theta) \propto P_l(\cos\theta)``. The
# **angular part** of the Laplacian (the piece that acts only on ``\theta``,
# for a function with no ``\varphi``-dependence) is
# ```math
# \nabla^2_{\text{angular}} f(\theta) = \frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\!\left(\sin\theta\,\frac{\partial f}{\partial\theta}\right).
# ```
# The claim we need is that ``P_l(\cos\theta)`` is an eigenfunction of this
# operator:
# ```math
# \nabla^2_{\text{angular}}\,P_l(\cos\theta) = -l(l+1)\,P_l(\cos\theta).
# ```
# The practical consequence, used constantly below, is that for any function
# of the separated form ``f(r)\,Y_l^m(\theta,\varphi)`` the full
# three-dimensional scalar Laplacian collapses to a purely radial operator:
# ```math
# \nabla^2\!\left[f(r)\,Y_l^m\right] = \left[\,f'' + \frac{2}{r}f' - \frac{l(l+1)}{r^2}\,f\,\right] Y_l^m.
# ```
#
# The eigenvalue claim splits cleanly into two independent facts, each
# checked below: (i) a change of variables from ``\theta`` to
# ``x=\cos\theta`` turns the angular Laplacian, acting on *any* function,
# into the "Legendre operator" ``\frac{d}{dx}\!\left[(1-x^2)\frac{dQ}{dx}\right]``;
# (ii) ``P_l(x)`` *specifically* satisfies Legendre's equation, so the
# Legendre operator acting on ``P_l`` gives exactly ``-l(l+1)P_l``.
#
# Fact (i) is checked for a generic cubic. That is not a restriction to
# cubics: both operators are linear, so an operator identity verified on
# ``1, x, x^2, x^3`` with independent coefficients holds generally.

@variables theta_sym c0 c1 c2 c3 x_leg                                                                               #src
Dth = Differential(theta_sym)                                                                                        #src
Dxleg = Differential(x_leg)                                                                                          #src
Q_theta = c0 + c1 * cos(theta_sym) + c2 * cos(theta_sym)^2 + c3 * cos(theta_sym)^3                                    #src
angular_lap_generic = expand_derivatives(Dth(Dth(Q_theta)) + (cos(theta_sym) / sin(theta_sym)) * Dth(Q_theta))        #src
Q_x = c0 + c1 * x_leg + c2 * x_leg^2 + c3 * x_leg^3                                                                  #src
legendre_op_generic = substitute(expand_derivatives(Dxleg((1 - x_leg^2) * Dxleg(Q_x))), Dict(x_leg => cos(theta_sym)))  #src
@assert symbolic_zero(angular_lap_generic - legendre_op_generic)                                                      #src

# Fact (ii) needs the Legendre polynomials themselves. They are built from
# Bonnet's three-term recursion ``(n+1)P_{n+1}(x)=(2n+1)xP_n(x)-nP_{n-1}(x)``
# -- the same construction used in the other derivations here -- and confirmed
# to satisfy Legendre's equation
# ``(1-x^2)P_l'' - 2xP_l' + l(l+1)P_l = 0`` at several concrete ``l``.
# (Concrete ``l``, because a CAS differentiates a polynomial of concrete
# degree readily but has no notion of a symbolic-degree polynomial.) This one
# definition is worth seeing, since the recursion is quoted by name several
# more times below:

function legendre_P(l::Int, x)
    l == 0 && return one(x)
    l == 1 && return x
    Pm1, P = one(x), x
    for n in 1:(l-1)
        P, Pm1 = ((2n + 1) * x * P - n * Pm1) / (n + 1), P
    end
    P
end
nothing # hide

for l_val in (2, 3, 4, 5)                                                                                                            #src
    Pl = legendre_P(l_val, x_leg)                                                                                                    #src
    legendre_eq_lhs = expand_derivatives((1 - x_leg^2) * Dxleg(Dxleg(Pl)) - 2 * x_leg * Dxleg(Pl) + l_val * (l_val + 1) * Pl)         #src
    @assert symbolic_zero(legendre_eq_lhs)                                                                                           #src
end                                                                                                                                  #src

# Both halves check out (for ``l=2,3,4,5``), so ``P_l(\cos\theta)`` is an
# eigenfunction of the angular Laplacian with eigenvalue ``-l(l+1)``.
#
# ### 2.2 The spherical Bessel substitution
#
# The ordinary Bessel equation of order ``\nu`` is
# ```math
# w'' + \frac{1}{z}w' + \left(1-\frac{\nu^2}{z^2}\right)w = 0,
# ```
# with solutions ``J_\nu(z)`` (first kind) and ``Y_\nu(z)`` (second kind,
# singular at the origin). The *spherical* Bessel equation of order ``l``,
# ```math
# v'' + \frac{2}{z}v' + \left(1-\frac{l(l+1)}{z^2}\right)v = 0,
# ```
# is what you get separating the Helmholtz equation ``(\nabla^2+k^2)F=0`` in
# spherical coordinates. Its regular solution is the spherical Bessel
# function of the first kind, ``j_l(z) = \sqrt{\pi/2z}\,J_{l+1/2}(z)``.
#
# The fact we need is that the ODE
# ```math
# U'' - \frac{l(l+1)}{x^2}U + q^2 U = 0
# ```
# (which we will meet again, forced, as Reid's Eq. 9) is solved by
# ``U(x) = x\,j_l(qx)``. The clean way to see this is a substitution:
# writing ``U=xv(x)`` should turn this ODE into the spherical Bessel
# equation for ``v`` (at ``q=1``; the general-``q`` case follows by
# ``x\to qx``), whose regular solution is ``v=j_l`` by definition. Working
# the chain rule, ``U'=v+xv'`` and ``U''=2v'+xv''``, so substituting gives
# ``x\left[v''+\frac{2}{x}v'+v-\frac{l(l+1)}{x^2}v\right]``, which is
# ``x`` times the spherical Bessel equation exactly.
#
# This identity is also *why* the velocity scaling ``u_r \propto U(x)/x^2``
# is chosen the way it is in Section 4: it is what makes ``x\,j_l(qx)`` the
# natural building block for everything downstream.

@variables x_sub                                                                                                             #src
@variables vfun(..)                                                                                                          #src
Dxs = Differential(x_sub)                                                                                                    #src
v = vfun(x_sub)                                                                                                              #src
U_from_v = x_sub * v                                                                                                         #src
U_ode_lhs = expand_derivatives(Dxs(Dxs(U_from_v)) - l_sym * (l_sym + 1) / x_sub^2 * U_from_v + U_from_v)                      #src
sph_bessel_eq = expand_derivatives(x_sub * (Dxs(Dxs(v)) + 2 / x_sub * Dxs(v) + (1 - l_sym * (l_sym + 1) / x_sub^2) * v))      #src
@assert symbolic_zero(U_ode_lhs - sph_bessel_eq)                                                                             #src

# We will also need the derivative of ``j_l`` at the surface. The standard
# recurrence ``j_l'(z) = \frac{l}{z}j_l(z) - j_{l+1}(z)``, divided through
# by ``j_l``, gives
# ```math
# \frac{q\,j_l'(q)}{j_l(q)} = l - q\,\frac{j_{l+1}(q)}{j_l(q)} = l - q\,Q_{l+1/2}(q),
# ```
# where ``Q_{l+1/2}(q) \equiv j_{l+1}(q)/j_l(q) = J_{l+3/2}(q)/J_{l+1/2}(q)``
# is Reid's ratio. This single substitution is what makes the final answer
# collapse into one Bessel combination rather than two; it is used in
# Section 7.
#
# ### 2.3 The poloidal decomposition
#
# Any divergence-free vector field ``\bm u`` decomposes uniquely into a
# **toroidal** part (of the form ``\nabla\times(\Psi\hat r)``, which has no
# radial component at all) and a **poloidal** part (of the form
# ``\nabla\times\nabla\times(\Phi\hat r)``, which does). For axisymmetric
# flow with no swirl -- our case, since the deformed drop has no preferred
# azimuthal direction and nothing sets it spinning -- the toroidal part
# vanishes identically and ``\bm u`` is purely poloidal.
#
# The consequence is the structural fact that makes this problem tractable:
# for a poloidal field with angular dependence ``Y_l^m``, the *entire*
# velocity field is determined by its radial component ``u_r``.
# Incompressibility then supplies ``u_\theta``, and there are no other free
# components. This is why a single scalar ODE for one function ``U(x)``,
# with ``u_r \propto U(x)/x^2``, can carry the whole flow -- and it is the
# statement invoked in Section 3 for the pressure gradient and in Section 4
# for the velocity itself.

# ## Section 3: Linearized governing equations and the pressure field
#
# Departures from equilibrium obey the linearized incompressible
# Navier-Stokes equations,
# ```math
# \frac{\partial \bm u}{\partial t} = -\nabla\frac{\delta p}{\rho} - \nu\,\nabla\times\nabla\times\bm u,
# \qquad \nabla\cdot\bm u = 0.
# ```
# (For incompressible flow ``\nabla\times\nabla\times\bm u = -\nabla^2\bm u``,
# so this is the same as the more familiar ``\partial_t\bm u = -\nabla(\delta p/\rho)+\nu\nabla^2\bm u``;
# both forms appear in the literature, and we use the second from here on.)
# With every field varying as ``e^{-\sigma t}`` (so ``\partial_t \to -\sigma``),
# ```math
# -\sigma\bm u = -\nabla\frac{\delta p}{\rho} + \nu\nabla^2\bm u.
# ```
# This is the equation everything downstream is built from.
#
# **Pressure satisfies Laplace's equation.** Taking the divergence of this
# equation: ``\nabla\cdot\bm u=0`` kills the left side, and
# ``\nabla\cdot(\nabla^2\bm u) = \nabla^2(\nabla\cdot\bm u) = 0`` kills the
# viscous term (divergence and Laplacian commute), leaving
# ```math
# \nabla^2\!\left(\frac{\delta p}{\rho}\right) = 0.
# ```
# We need the solution with angular dependence ``Y_l^m`` that is regular at
# the origin. Writing ``\delta p/\rho = \epsilon\,f(r)\,Y_l^m``, the scalar
# Laplacian formula from Section 2.1 turns this into the radial ODE
# ``f'' + (2/r)f' - l(l+1)f/r^2 = 0``, whose general solution is
# ``f = Ar^l + Br^{-(l+1)}`` for symbolic ``l``. Regularity at ``r=0``
# requires ``B=0``, so with ``x=r/R``,
# ```math
# \frac{\delta p}{\rho} = \epsilon\, P_0\, x^l\, Y_l^m,
# ```
# for some constant ``P_0`` (dimensions of velocity squared) fixed later by
# the boundary conditions.

@variables x l A B                                                                                                    #src
Dx_ = Differential(x)                                                                                                 #src
f_pressure = A * x^l + B * x^(-(l + 1))                                                                               #src
radial_laplace_residual = expand_derivatives(Dx_(Dx_(f_pressure)) + 2 / x * Dx_(f_pressure) - l * (l + 1) / x^2 * f_pressure)  #src
@assert symbolic_zero(radial_laplace_residual)                                                                        #src

# The pressure *gradient* is a poloidal, divergence-free field with angular
# structure ``Y_l^m``, so by Section 2.3 it is entirely determined by a
# single scalar function ``\Pi(x)`` through its radial component. Reid's
# defining relation for that scalar is
# ```math
# \left(\nabla\frac{\delta p}{\rho}\right)_r
#   = \epsilon_0\,\sigma^2 R\,\frac{\Pi(x)}{x^2}\,Y_l^m\,e^{-\sigma t},
# ```
# and computing the same radial derivative explicitly from
# ``\delta p/\rho = \epsilon P_0 x^l Y_l^m`` gives
# ```math
# \left(\nabla\frac{\delta p}{\rho}\right)_r
#   = \frac{\partial}{\partial r}\!\left(\epsilon P_0 x^l Y_l^m\right)
#   = \epsilon_0 e^{-\sigma t}\,P_0\,\frac{l}{R}\,x^{l-1}\,Y_l^m.
# ```
# Matching the two forces ``\Pi(x) = \Pi_0 x^{l+1}`` with
# ```math
# \Pi_0 = \frac{l}{\sigma^2 R^2}\,P_0.
# ```
# ``P_0`` (equivalently ``\Pi_0``) is the second of the three constants the
# boundary conditions will pin down; the third, ``C``, appears once the
# velocity field enters next.

# ## Section 4: The velocity field and its governing ODE (Reid's Eq. 9)
#
# This is the technical heart of the whole problem: one ODE that packages
# up the entire viscous, incompressible flow field consistent with the
# linearized momentum equation.
#
# Since the pressure gradient is purely poloidal, the momentum equation
# forces the velocity to be purely poloidal too (nothing drives a toroidal
# part -- Section 2.3). We write the radial velocity as
# ```math
# u_r = \epsilon_0\,\sigma R\,\frac{U(x)}{x^2}\,Y_l^m\,e^{-\sigma t}
# ```
# (the ``1/x^2`` scaling is conventional -- it is exactly what makes the
# final ODE for ``U`` take the clean Bessel-equation form derived in
# Section 2.2, which is the entire point of choosing it). Write
# ``G(x)=U(x)/x^2``, so ``u_r \propto G(x)``.
#
# !!! warning "Where the Newtonian assumption enters -- read this before adapting the page"
#     The momentum equation above carries a single viscous term,
#     ``\nu\nabla^2\bm u``: a *scalar, constant* ``\nu`` multiplying the
#     Laplacian of the velocity. That is the Newtonian constitutive law --
#     stress proportional to strain rate, with a viscosity that does not
#     depend on the flow itself. A shear-thinning fluid (Carreau-Yasuda)
#     replaces that scalar with a viscosity that depends on the local strain
#     rate, which in general no longer factors out of the Laplacian this
#     cleanly. Everything in this section and the next -- the ODE for
#     ``U``, both stress boundary conditions, and the characteristic
#     equation itself -- rests on this one substitution and does **not**
#     carry over unmodified. Section 10 spells out exactly which parts
#     survive a change of rheology and which do not.
#
# For an axisymmetric poloidal field, the incompressibility constraint
# lets you eliminate the angular (``u_\theta``) part of the vector
# Laplacian entirely in favor of radial derivatives of ``u_r``, giving
# ```math
# \left[\nabla^2\bm u\right]_r = \nabla^2 u_r + \frac{2}{r^2}u_r + \frac{2}{r}\frac{\partial u_r}{\partial r}.
# ```
# We derive this from two standard, textbook spherical-coordinate operator
# formulas (these are coordinate-system facts, not specific to this
# problem, so we cite them the way we'd cite "the Laplacian in spherical
# coordinates has this form" -- but everything built from them below is
# checked). For an axisymmetric field with no ``\varphi``-dependence and no
# ``u_\varphi``, the raw r-component of the vector Laplacian is
# ```math
# \left[\nabla^2\bm u\right]_r = \nabla^2 u_r - \frac{2}{r^2}u_r - \frac{2}{r^2}A,
# \qquad
# A \equiv \frac{1}{\sin\theta}\frac{\partial(\sin\theta\,u_\theta)}{\partial\theta},
# ```
# and incompressibility (``\nabla\cdot\bm u=0``) in the same coordinates is
# ```math
# \frac{1}{r^2}\frac{\partial(r^2 u_r)}{\partial r} + \frac{1}{r}A = 0.
# ```
# Solving the second for ``A`` gives ``A = -(2u_r + r\,\partial_r u_r)`` --
# purely in terms of ``u_r``'s radial dependence, no explicit ``u_\theta``
# needed. Substituting this into the raw formula reproduces the target
# formula above exactly, treating ``A``, ``u_r``, and ``\partial_r u_r`` as
# independent symbols related only by that one incompressibility identity.

@variables r_var A_ang ur_r ur_r_prime Lap_ur                                                                            #src
raw_vec_lap_r = Lap_ur - 2 / r_var^2 * ur_r - 2 / r_var^2 * A_ang                                                        #src
incompressibility_A = -(2 * ur_r + r_var * ur_r_prime)   # from (1/r^2)(r^2 u_r)' = -(1/r)*A, solved for A                #src
target_vec_lap_r = Lap_ur + 2 / r_var^2 * ur_r + 2 / r_var * ur_r_prime                                                  #src
vec_lap_residual = simplify(substitute(raw_vec_lap_r, Dict(A_ang => incompressibility_A)) - target_vec_lap_r; expand=true)  #src
@assert symbolic_zero(vec_lap_residual)                                                                                  #src

# Substituting ``u_r=\sigma R\,G(x)\,Y_l^m`` (the ``\epsilon_0 e^{-\sigma t}``
# prefactor cancels throughout, exactly as in the source) and using the
# scalar Laplacian formula (Section 2.1) for ``\nabla^2 u_r``, all three
# terms -- the scalar Laplacian, ``2u_r/r^2``, and ``(2/r)\partial_r u_r``
# -- combine into a single operator on ``G``:
# ```math
# \left[\nabla^2\bm u\right]_r \;\propto\; G'' + \frac{4}{x}G' + \frac{2-l(l+1)}{x^2}G.
# ```

@variables Gfun(..)                                                                                                   #src
G = Gfun(x)                                                                                                           #src
term_scalar_lap = Dx_(Dx_(G)) + 2 / x * Dx_(G) - l * (l + 1) / x^2 * G   # nabla^2 u_r, via Section 2.1                #src
term_2_over_r2 = 2 / x^2 * G                                             # (2/r^2) u_r                                 #src
term_2_over_r_deriv = 2 / x * Dx_(G)                                     # (2/r) du_r/dr                               #src
combined_operator = expand_derivatives(term_scalar_lap + term_2_over_r2 + term_2_over_r_deriv)                         #src
target_operator = expand_derivatives(Dx_(Dx_(G)) + 4 / x * Dx_(G) + (2 - l * (l + 1)) / x^2 * G)                       #src
@assert symbolic_zero(combined_operator - target_operator)                                                             #src

# The r-component of the momentum equation (after canceling ``Y_l^m
# e^{-\sigma t}`` and using ``x=r/R``, ``q^2=\sigma R^2/\nu``) is then
# ```math
# -q^2 G = -\frac{P_0\,l\,x^{l-1}}{\sigma\nu} + G'' + \frac{4}{x}G' + \frac{2-l(l+1)}{x^2}G.
# ```
# Now substitute ``G = U/x^2`` -- this is where the ``1/x^2`` scaling earns
# its keep. It turns the bracket above into ``[U'' - l(l+1)U/x^2]/x^2``
# exactly, verified below for symbolic ``l``.

@variables Ufun(..)                                                                                                   #src
Usym = Ufun(x)                                                                                                        #src
Gsub = Usym / x^2                                                                                                     #src
lhs_full = expand_derivatives(Dx_(Dx_(Gsub)) + 4 / x * Dx_(Gsub) + (2 - l * (l + 1)) / x^2 * Gsub)                     #src
rhs_full = expand_derivatives((Dx_(Dx_(Usym)) - l * (l + 1) / x^2 * Usym) / x^2)                                       #src
@assert symbolic_zero(lhs_full - rhs_full)                                                                             #src

# Multiplying the momentum equation through by ``x^2`` and using the
# pressure solution ``P_0 = \sigma^2 R^2 \Pi_0/l`` (Section 3) to write
# ``P_0 l/(\sigma\nu) = \sigma R^2\Pi_0/\nu = q^2\Pi_0``, the momentum
# equation becomes exactly **Reid's Eq. 9** -- the single equation that
# packages up the entire linearized, viscous, incompressible flow problem:

#md # ```@eval
#md # using Symbolics, Markdown
#md # @variables x l q Pi_0
#md # rhs = q^2 * Pi_0 * x^(l + 1)
#md # Markdown.parse("```math\n\\left[\\frac{d^2}{dx^2} - \\frac{l(l+1)}{x^2} + q^2\\right] U(x) = " *
#md #     Main.pretty_latex(rhs) * "\n```")
#md # ```

# ## Section 5: Solving the ODE
#
# Reid's Eq. 9 is a linear, second-order, inhomogeneous ODE. Its general
# solution is a particular solution (matching the ``x^{l+1}`` forcing) plus
# the general solution of the homogeneous equation.
#
# **Particular solution.** Try ``U_p = \Pi_0 x^{l+1}`` directly. The
# ``l(l+1)/x^2`` and ``q^2`` pieces of ``U_p''`` cancel the corresponding
# pieces on the left, leaving exactly the ``q^2\Pi_0 x^{l+1}`` forcing on
# the right -- so it is an exact solution, for symbolic ``l`` and ``q``.

@variables x l q Pi_0                                                                                                        #src
Dx_ = Differential(x)                                                                                                        #src
Up = Pi_0 * x^(l + 1)                                                                                                        #src
particular_residual = expand_derivatives(Dx_(Dx_(Up)) - l * (l + 1) / x^2 * Up + q^2 * Up - q^2 * Pi_0 * x^(l + 1))           #src
@assert symbolic_zero(particular_residual)                                                                                   #src

# **Homogeneous solution.** The homogeneous equation
# ``U''-l(l+1)U/x^2+q^2U=0`` is Section 2.2's ODE with ``x\to qx``: writing
# ``U_h = x\,v(qx)``, the chain rule gives ``U_h'=v+xqv'``,
# ``U_h''=2qv'+xq^2v''``, and substituting ``v`` satisfying the spherical
# Bessel equation at argument ``z=qx`` (``v''=-\tfrac{2}{z}v'-(1-\tfrac{l(l+1)}{z^2})v``)
# collapses the whole thing to zero. The check below treats ``v, v', v''``
# as related by exactly that one substitution rather than evaluating an
# actual Bessel function -- this is the general, function-independent
# statement, exactly analogous to Section 2.2. So ``U_h = C\,x\,j_l(qx)``
# for any constant ``C``.

@variables v vp vpp                                                                                                                     #src
z = q * x                                                                                                                               #src
Uh = x * v                                                                                                                              #src
Uh_pp = 2q * vp + x * q^2 * vpp   # U_h'' = 2q*vp + x*q^2*vpp, by the chain rule shown in prose above                                    #src
vpp_relation = -(2 / z) * vp - (1 - l * (l + 1) / z^2) * v                                                                               #src
homogeneous_residual = simplify(substitute(Uh_pp, Dict(vpp => vpp_relation)) - l * (l + 1) / x^2 * Uh + q^2 * Uh; expand=true)           #src
@assert isequal(homogeneous_residual, 0)                                                                                                #src

# The other homogeneous solution -- built from the second-kind spherical
# Bessel function ``n_l``, which diverges as ``z\to0`` -- is discarded on
# regularity grounds, the same argument as the pressure field's ``B``-term
# in Section 3. The general solution is therefore

#md # ```@eval
#md # using Symbolics, Markdown
#md # @variables x l q C Pi_0
#md # Markdown.parse("```math\n" * Main.pretty_latex(C*x) * "\\,j_l(qx) + " *
#md #     Main.pretty_latex(Pi_0*x^(l+1)) * "\n```")
#md # ```
#
# where ``C`` and ``\Pi_0`` are fixed by the three boundary conditions we
# derive next.

# ## Section 6: The three boundary conditions
#
# Three physical conditions at the drop surface (``x=1``) fix the two
# constants ``C, \Pi_0`` and, together, produce the characteristic
# equation.
#
# ### BC1: Kinematic condition
#
# The fluid at the surface must move *with* the surface -- the surface isn't
# a permeable membrane. The surface itself moves at
# ``\partial r/\partial t|_{\text{surface}} = -\sigma R\epsilon_0 e^{-\sigma t}Y_l^m``
# (differentiating the surface ansatz directly), while the fluid's radial
# velocity there is ``u_r|_{x=1} = \epsilon_0\sigma R\,U(1)\,Y_l^m e^{-\sigma t}``.
# Equating the two and canceling the common prefactor ``\epsilon_0\sigma``
# gives, exactly,
# ```math
# U(1) = -1.
# ```
# The condition is imposed at ``x=1``, the undeformed surface: evaluating
# it at the true deformed surface ``r=R[1+\epsilon Y_l^m]`` would only add
# ``O(\epsilon^2)`` corrections, negligible at this linear order.

@variables sigma_sym epsilon0 U1                                                                     #src
surface_velocity = -sigma_sym * epsilon0            # d/dt of the surface position                   #src
fluid_velocity_at_1 = epsilon0 * sigma_sym * U1      # u_r at x=1, with U(1) still unknown           #src
bc1_solution = Symbolics.solve_for(surface_velocity - fluid_velocity_at_1 ~ 0, U1)                   #src
@assert symbolic_zero(bc1_solution - (-1))                                                           #src

# ### BC2: Tangential stress condition
#
# At a free surface (no exterior viscous fluid), the tangential viscous
# stress must vanish:
# ```math
# \tau_{r\theta} = \mu\left[r\frac{\partial}{\partial r}\!\left(\frac{u_\theta}{r}\right) + \frac{1}{r}\frac{\partial u_r}{\partial\theta}\right] = 0 \qquad \text{at } x=1.
# ```
# This needs ``u_\theta``, which hasn't shown up until now. For an
# axisymmetric poloidal field, the standard Stokes stream function ``\psi(r,\theta)``
# gives ``u_r=(1/(r^2\sin\theta))\,\partial\psi/\partial\theta`` and
# ``u_\theta=-(1/(r\sin\theta))\,\partial\psi/\partial r``. Writing
# ``u_r=f(r)P_l(\cos\theta)`` and integrating in ``\theta`` (using the
# standard Legendre recurrence ``(2l+1)P_l=P_{l+1}'-P_{l-1}'``) gives
# ```math
# u_\theta = \frac{g(r)}{\sin\theta}\Big[P_{l+1}(\cos\theta)-P_{l-1}(\cos\theta)\Big], \qquad g(r)\equiv\frac{2f(r)+rf'(r)}{2l+1}.
# ```
# Substituting both into ``\tau_{r\theta}`` and using a second standard
# recurrence, ``(2l+1)(1-x^2)P_l'(x)=l(l+1)[P_{l-1}(x)-P_{l+1}(x)]``, every
# term collapses onto the single common angular factor
# ``\sin\theta\,P_l'(\cos\theta)`` -- exactly the shape a genuine
# tangential-stress condition should have -- times a purely radial
# coefficient. Setting that coefficient to zero at ``r=R`` turns out to be
# **exactly** ``\mathcal{L}_2[U]=0``:
# ```math
# \mathcal{L}_2[U] \equiv \left[\frac{d^2}{dx^2} - \frac{2}{x}\frac{d}{dx} + \frac{l(l+1)}{x^2}\right]U = 0 \qquad\text{at } x=1,
# ```
# so BC2, usually quoted in this operator form, follows from the
# stream-function representation directly.
#
# Both Legendre recurrences hold for ``l=2,3,4,5`` against the Bonnet-built
# polynomials of Section 2.1, so the stream-function ansatz does produce a
# pure tangential-stress condition proportional to ``\mathcal{L}_2[U]``.

@variables x_leg2                                                                                                  #src
Dxleg2 = Differential(x_leg2)                                                                                      #src
for l_val in (2, 3, 4, 5)                                                                                          #src
    Pl_a, Plp1_a, Plm1_a = legendre_P(l_val, x_leg2), legendre_P(l_val + 1, x_leg2), legendre_P(l_val - 1, x_leg2)  #src
    @assert symbolic_zero((2l_val + 1) * Pl_a - expand_derivatives(Dxleg2(Plp1_a) - Dxleg2(Plm1_a)))                #src
end                                                                                                                #src

for l_val in (2, 3, 4, 5)                                                                                          #src
    Pl_b, Plp1_b, Plm1_b = legendre_P(l_val, x_leg2), legendre_P(l_val + 1, x_leg2), legendre_P(l_val - 1, x_leg2)  #src
    lhs_b = expand_derivatives((2l_val + 1) * (1 - x_leg2^2) * Dxleg2(Pl_b))                                        #src
    @assert symbolic_zero(lhs_b - l_val * (l_val + 1) * (Plm1_b - Plp1_b))                                          #src
end                                                                                                                #src

# Building ``\tau_{r\theta}/\mu`` directly from the stream function (concrete
# ``l``, abstract radial function ``f(r)``) makes the collapse explicit:
# ``\tau_{r\theta}=0`` at ``r=R`` is equivalent to
# ```math
# R^2 f''(R) + 2R f'(R) + \big[l(l+1)-2\big] f(R) = 0,
# ```
# and imposing that single relation among ``f(R)``, ``f'(R)``, ``f''(R)``
# makes ``\tau_{r\theta}`` vanish identically in ``\theta``, not merely at one
# angle -- for ``l=2,\dots,6``.

@variables r_var theta_sym2 R_sym F0 F1 F2 Ffun(..)                                                                                                        #src
Dr_ = Differential(r_var)                                                                                                                                  #src
Dth_ = Differential(theta_sym2)                                                                                                                            #src
for l_val in (2, 3, 4, 5, 6)                                                                                                                               #src
    f = Ffun(r_var)                                                                                                                                        #src
    g = (2 * f + r_var * Dr_(f)) / (2 * l_val + 1)                                                                                                         #src
    Pl_th, Plp1_th, Plm1_th = legendre_P(l_val, cos(theta_sym2)), legendre_P(l_val + 1, cos(theta_sym2)), legendre_P(l_val - 1, cos(theta_sym2))            #src
    u_r_ = f * Pl_th                                                                                                                                       #src
    u_theta_ = (g / sin(theta_sym2)) * (Plp1_th - Plm1_th)                                                                                                 #src
    tau_over_mu = expand_derivatives(r_var * Dr_(u_theta_ / r_var) + (1 / r_var) * Dth_(u_r_))                                                              #src
    at_R = substitute(tau_over_mu, Dict(Dr_(Dr_(f)) => F2, Dr_(f) => F1, f => F0))                                                                          #src
    at_R = substitute(at_R, Dict(r_var => R_sym))                                                                                                          #src
    F2_conjectured = -(2 * R_sym * F1 + (l_val * (l_val + 1) - 2) * F0) / R_sym^2                                                                           #src
    @assert symbolic_zero(substitute(at_R, Dict(F2 => F2_conjectured)))                                                                                     #src
end                                                                                                                                                        #src

# Finally, translate ``f(r)``'s condition into ``U(x)``. Write
# ``f(r)=\kappa\,G(r/R)`` for some overall constant ``\kappa`` (which
# cancels, since the condition above is linear and homogeneous in ``f, f',
# f''``), with ``G(x)=U(x)/x^2`` from Section 4. The chain rule gives
# ``f(R)=\kappa G(1)``, ``f'(R)=\kappa G'(1)/R``, ``f''(R)=\kappa G''(1)/R^2``,
# so the condition becomes ``G''(1)+2G'(1)+[l(l+1)-2]G(1)=0``, which in turn
# is exactly ``\mathcal{L}_2[U]|_{x=1}=U''(1)-2U'(1)+l(l+1)U(1)=0`` for
# symbolic ``l``. BC2's operator is thus derived end to end: stream function
# ``\to`` ``\tau_{r\theta}=0`` ``\to`` the ``f``-ODE ``\to`` the ``G``-ODE
# ``\to`` ``\mathcal{L}_2[U]=0``.

@variables xx l q Pi_0 C U0 U1 U2 Ufun2(..)                                                                        #src
Dxx = Differential(xx)                                                                                             #src
Gexpr = Ufun2(xx) / xx^2                                                                                           #src
stepU(e) = substitute(e, Dict(Dxx(Dxx(Ufun2(xx))) => U2, Dxx(Ufun2(xx)) => U1, Ufun2(xx) => U0))                    #src
G0 = simplify(stepU(substitute(Gexpr, Dict(xx => 1))); expand=true)                                                 #src
G1 = simplify(substitute(stepU(expand_derivatives(Dxx(Gexpr))), Dict(xx => 1)); expand=true)                        #src
G2 = simplify(substitute(stepU(expand_derivatives(Dxx(Dxx(Gexpr)))), Dict(xx => 1)); expand=true)                   #src
G_condition = G2 + 2 * G1 + (l * (l + 1) - 2) * G0                                                                  #src
L2_target = U2 - 2 * U1 + l * (l + 1) * U0                                                                          #src
@assert symbolic_zero(G_condition - L2_target)                                                                      #src

# Evaluating ``\mathcal{L}_2`` on the two pieces of the general solution
# gives, exactly and for symbolic ``l``,
# ```math
# \mathcal{L}_2[U_p] = 2(l-1)(l+1)\,\Pi_0\,x^{l-1},
# \qquad
# \mathcal{L}_2[U_h]\big|_{x=1} = C\!\left[-q^2 j_l(q) + 2(l^2+l-1)j_l(q) - 2q\,j_l'(q)\right],
# ```
# the second after eliminating ``v''`` via the spherical Bessel equation at
# ``z=q``, and matching Reid's own stated coefficient. Setting the sum
# ``\mathcal{L}_2[U_h]|_{x=1} + \mathcal{L}_2[U_p]|_{x=1} = 0`` (BC2) gives
# one linear relation between ``C`` and ``\Pi_0``, used in Section 7.

x = xx                                                                                                             #src
Dx_ = Differential(x)                                                                                              #src
L2(U) = expand_derivatives(Dx_(Dx_(U)) - 2 / x * Dx_(U) + l * (l + 1) / x^2 * U)                                    #src
Up = Pi_0 * x^(l + 1)                                                                                              #src
L2_Up_target = 2 * (l - 1) * (l + 1) * Pi_0 * x^(l - 1)                                                            #src
@assert symbolic_zero(L2(Up) - L2_Up_target)                                                                       #src

@variables v vp vpp                                                                                                #src
z = q * x                                                                                                          #src
Uh_at_1, Uh_p_at_1, Uh_pp_at_1 = C * v, C * (v + q * vp), C * (2q * vp + q^2 * vpp)                                 #src
L2_Uh_at_1 = Uh_pp_at_1 - 2 * Uh_p_at_1 + l * (l + 1) * Uh_at_1                                                    #src
vpp_relation = -(2 / q) * vp - (1 - l * (l + 1) / q^2) * v                                                          #src
L2_Uh_target = C * (-q^2 * v + 2 * (l^2 + l - 1) * v - 2 * q * vp)                                                  #src
@assert isequal(simplify(substitute(L2_Uh_at_1, Dict(vpp => vpp_relation)) - L2_Uh_target; expand=true), 0)         #src

# ### BC3: Normal stress condition
#
# The most involved of the three: it couples pressure, viscous normal
# stress, and the curvature of the deformed surface through surface
# tension. The radial normal stress inside the drop is
# ``-p_{rr} = p + \delta p - 2\mu\,\partial u_r/\partial r``, and the free-surface
# condition is ``-p_{rr} = T_1(1/R_1+1/R_2)``.
#
# To first order in ``\epsilon``, the mean curvature of the perturbed
# surface is
# ```math
# \frac{1}{R_1}+\frac{1}{R_2} = \frac{1}{R}\Big[2 + (l-1)(l+2)\,\epsilon\,Y_l^m\Big].
# ```
# The ``2/R`` piece is just the equilibrium curvature (Young-Laplace,
# Section 1); the ``(l-1)(l+2)`` factor is the same combination that
# appears in the inviscid frequency ``\sigma_{l;0}^2`` -- both come from the
# curvature response to a harmonic deformation, so this is a consistency
# check worth noting, not a coincidence.
#
# **Where this comes from.** The starting point is the standard linearized
# mean-curvature formula of differential geometry, for a nearly-spherical
# surface ``r=R+\zeta(\theta,\varphi)`` with ``\zeta`` small:
# ```math
# \frac{1}{R_1}+\frac{1}{R_2} = \frac{2}{R} - \frac{1}{R^2}\Big[2\zeta + \Delta_\Omega\zeta\Big].
# ```
#
# !!! note "A dimensional subtlety in that formula"
#     ``\Delta_\Omega`` here is the Laplace-Beltrami operator on the *unit*
#     sphere, which is dimensionless. It is not the same object as Section
#     2.1's ``\nabla^2_{\text{angular}}``, which carries a ``1/r^2`` because
#     it is the angular piece of the full three-dimensional Laplacian. The
#     two agree on the only property used here -- both have eigenvalue
#     ``-l(l+1)`` on ``Y_l^m`` -- but confusing them costs a factor of
#     ``R^2``, and the ``1/R^2`` prefactor above is where that factor has
#     already been accounted for.
#
# The algebraic collapse from there is immediate: with
# ``\zeta=\epsilon R\,Y_l^m`` and the angular eigenvalue property of Section
# 2.1, this reproduces the boxed formula above exactly, for symbolic ``l``.

@variables l_sym eps_sym R_sym Yl                                                                                     #src
zeta = eps_sym * R_sym * Yl                                                                                           #src
angular_lap_zeta = eps_sym * R_sym * (-l_sym * (l_sym + 1) * Yl)   # Section 2.1's eigenvalue property applied to zeta #src
curvature_from_formula = 2 / R_sym - (1 / R_sym^2) * (2 * zeta + angular_lap_zeta)                                     #src
curvature_target = (1 / R_sym) * (2 + (l_sym - 1) * (l_sym + 2) * eps_sym * Yl)                                        #src
@assert symbolic_zero(simplify(curvature_from_formula - curvature_target; expand=true))                                #src

# At ``O(\epsilon^0)`` the equilibrium Young-Laplace balance is automatically
# satisfied. At ``O(\epsilon^1)``, using the pressure solution
# (``\delta p|_{x=1}=\rho\epsilon P_0 Y_l^m``) and
# ``u_r=\epsilon_0\sigma R\,G(x)Y_l^m e^{-\sigma t}`` for the viscous term,
# the one genuinely computational piece is ``\partial u_r/\partial r`` at
# the surface, which needs ``d(U/x^2)/dx`` at ``x=1``. That equals
# ``U'(1)-2U(1)`` exactly.

@variables x                                                                                        #src
@variables U(x)                                                                                     #src
Dx__ = Differential(x)                                                                              #src
dGdx_at_1 = substitute(expand_derivatives(Dx__(U / x^2)), Dict(x => 1))                             #src
target_at_1 = substitute(Dx__(U) - 2 * U, Dict(x => 1))                                             #src
@assert isequal(simplify(dGdx_at_1 - target_at_1; expand=true), 0)                                   #src

# Dividing through by ``\rho\epsilon Y_l^m`` and using ``\mu=\rho\nu``, BC3
# becomes

# ```math
# \frac{(l-1)(l+2)\,T_1}{\rho R}
#   \;=\; P_0 - 2\nu\sigma\bigl[\,U'(1) - 2U(1)\,\bigr],
# ```
#
# the only equation containing the surface tension ``T_1`` explicitly.
# Reid's key move, next, is to rescale so that ``T_1`` (and ``\rho``) drop
# out entirely, leaving a problem in terms of two purely dimensionless
# numbers.

# ## Section 7: The characteristic equation
#
# ### Eliminating the surface tension
#
# Define ``\alpha^2 \equiv \sigma_{l;0}R^2/\nu`` (the inviscid frequency,
# measured in units of the viscous diffusion rate ``\nu/R^2``) alongside
# ``q^2=\sigma R^2/\nu`` from Section 4. Using ``\sigma_{l;0}^2 =
# l(l-1)(l+2)T_1/(\rho R^3)`` (Section 1) to rewrite BC3's left side, and
# ``P_0=\sigma^2R^2\Pi_0/l`` (Section 3) to rewrite its right side, then
# multiplying through by ``lR^2/\nu^2``, the left side becomes
# ``\alpha^4 = \sigma_{l;0}^2R^4/\nu^2`` and the right side becomes purely a
# function of ``q`` and the boundary-condition unknowns -- both verified
# below for symbolic ``l``:

@variables l q alpha nu sigma R T1 rho sigma_l0 Pi_0 Uterm P0                                                                                            #src
bc3_lhs = (l - 1) * (l + 2) * T1 / (rho * R)                                                                                                             #src
bc3_rhs = P0 - 2 * nu * sigma * Uterm                                                                                                                    #src
lhs_rescaled = simplify(substitute(bc3_lhs, Dict(T1 => sigma_l0^2 * rho * R^3 / (l * (l - 1) * (l + 2)))) * l * R^2 / nu^2; expand=true)                  #src
rhs_rescaled = simplify(substitute(bc3_rhs, Dict(P0 => sigma^2 * R^2 * Pi_0 / l)) * l * R^2 / nu^2; expand=true)                                          #src
@assert symbolic_zero(lhs_rescaled - sigma_l0^2 * R^4 / nu^2)                                                                                             #src
rhs_in_q = simplify(substitute(rhs_rescaled, Dict(sigma => q^2 * nu / R^2)); expand=true)   # sigma = q^2*nu/R^2                                          #src
@assert symbolic_zero(rhs_in_q - q^2 * (q^2 * Pi_0 - 2 * l * Uterm))                                                                                      #src

#md # ```@eval
#md # using Symbolics, Markdown
#md # @variables l q Pi_0 alpha Uterm
#md # rhs = q^2*(q^2*Pi_0 - 2*l*Uterm)
#md # Markdown.parse("```math\n\\alpha^4 = " * Main.pretty_latex(rhs) * "\n```")
#md # ```
#
# ``T_1`` and ``\rho`` have cancelled entirely. This is Reid's key
# rescaling, and it is the reason the result carries over unchanged to any
# other spherically-symmetric restoring force: the physical mechanism enters
# only through ``\sigma_{l;0}``, hence only through ``\alpha``. Surface
# tension here; self-gravity in Chandrasekhar's version of the same problem.
#
# ### Solving BC1 + BC2 for ``C, \Pi_0``, and the characteristic equation itself
#
# BC1 (``U(1)=C j_l(q)+\Pi_0=-1``) and BC2 (Section 6, evaluated on the
# general solution) are two linear equations in the two unknowns
# ``C,\Pi_0``. Solving them and then using the Bessel recurrence
# ``qj_l'/j_l = l-qQ_{l+1/2}(q)`` from Section 2.2 to eliminate ``j_l'`` in
# favor of ``Q_{l+1/2}(q)`` gives
# ```math
# C = \frac{2(l-1)(l+1)}{j_l(q)\,q\,\bigl(2Q_{l+1/2}(q)-q\bigr)},
# ```
# which is Reid's Eq. 20.

@variables jl jlp Q C                                                                             #src
bc1_eq = C * jl + Pi_0 ~ -1                                                                       #src
bc2_eq = C * (-q^2 * jl + 2 * (l^2 + l - 1) * jl - 2 * q * jlp) + 2 * (l^2 - 1) * Pi_0 ~ 0         #src
bc_solution = Symbolics.solve_for([bc1_eq, bc2_eq], [C, Pi_0])                                    #src
Csol = simplify(substitute(bc_solution[1], Dict(jlp => jl * (l - q * Q) / q)); expand=true)        #src
Pi0sol = simplify(substitute(bc_solution[2], Dict(jlp => jl * (l - q * Q) / q)); expand=true)      #src
@assert symbolic_zero(Csol - 2 * (l^2 - 1) / (jl * q * (2Q - q)))                                  #src

# Finally: substitute ``C, \Pi_0`` into ``U'(1)-2U(1)`` (using ``U(1)=-1``
# from BC1), then into the ``\alpha^4`` relation above. This is the largest
# algebraic reduction in the derivation, and it is carried out symbolically
# here rather than by hand.
#
# In plain terms: every physically meaningful outcome of this whole
# derivation -- the oscillation frequency, the damping rate, whether a drop
# of a given size and viscosity oscillates or just squashes back down --
# follows from exactly this one equation:

#md # ```@eval
#md # using Symbolics, Markdown
#md # @variables l q Q
#md # rhs = (2*(l-1)/q^2) * (l + (l+1)*(q-2*l*Q)/(q-2*Q))
#md # Markdown.parse("```math\n\\frac{\\alpha^4}{q^4} + 1 = " * Main.pretty_latex(rhs) *
#md #     "\n```\nwhere \$q^2=\\sigma R^2/\\nu\$, \$\\alpha^2=\\sigma_{l;0}R^2/\\nu\$, and " *
#md #     "\$Q_{l+1/2}(q)=j_{l+1}(q)/j_l(q)\$.")
#md # ```

## NOTE: two SEPARATE substitute passes, not one combined dict -- C/Pi_0's                                                                        #src
## own raw solutions still contain jlp, and a single-pass substitute does                                                                         #src
## not recursively re-process values it has just substituted in.                                                                                  #src
Uprime1_minus_2U1 = substitute(C * (jl + q * jlp) + Pi_0 * (l + 1) + 2, Dict(C => bc_solution[1], Pi_0 => bc_solution[2]))                         #src
Uprime1_minus_2U1 = simplify(substitute(Uprime1_minus_2U1, Dict(jlp => jl * (l - q * Q) / q)); expand=true)                                        #src
alpha4_derived = simplify(q^2 * (q^2 * Pi0sol - 2 * l * Uprime1_minus_2U1); expand=true)                                                           #src
characteristic_eq_rhs = q^4 * ((2 * (l - 1) / q^2) * (l + (l + 1) * (q - 2 * l * Q) / (q - 2 * Q)) - 1)                                            #src
@assert symbolic_zero(alpha4_derived - characteristic_eq_rhs)                                                                                      #src

# This is Reid (1960)'s Eq. 19, and it is the derivation's main result: every
# damping rate and frequency DropSolver computes traces back to it.
#
# !!! note "Universality"
#     This equation is identical to Chandrasekhar's characteristic equation
#     for a viscous, self-gravitating liquid globe, with the
#     self-gravitational parameter identified with ``\alpha^2``. The
#     physical restoring force -- surface tension here, self-gravity there
#     -- enters only through ``\sigma_{l;0}``, which defines ``\alpha``.
#     The equation holds for any spherically symmetric restoring force.

# ## Section 8: The structure of the solutions
#
# The characteristic equation is transcendental in ``q^2=\sigma R^2/\nu``,
# so before extracting any number from it we need to know what its solution
# set looks like: how many roots there are, which of them are physical, and
# how their character changes with viscosity.
#
# ### Why there are infinitely many roots
#
# The Bessel function ``J_{l+1/2}(q)`` sitting in the denominator of
# ``Q_{l+1/2}(q)`` has infinitely many real zeros ``q_1<q_2<\cdots`` on the
# positive real axis. Near each of those zeros the characteristic equation
# has a pole, and each pole "captures" one pair of roots. Physically these
# are increasingly high radial overtones: modes with more and more nodal
# surfaces in the radial direction *inside* the drop.
#
# ### Which roots matter
#
# Roots with smaller ``\mathrm{Re}(\sigma)`` decay more slowly, so they are
# the last to die out and are the ones that govern the observable
# oscillation. The two roots with the smallest ``\mathrm{Re}(\sigma)``
# correspond to the fundamental surface oscillation at harmonic order ``l``;
# every higher root decays exponentially faster and is invisible on the
# timescale of interest.
#
# !!! warning "The sign convention inverts the usual instinct"
#     Here ``\epsilon\sim e^{-\sigma t}`` with ``\mathrm{Re}(\sigma)>0``
#     meaning decay, so the *smallest* positive real part is the dominant,
#     slowest-decaying mode. This is the opposite of ordinary stability
#     analysis, where one hunts for the most-positive growth rate. Reversing
#     it selects the wrong root; *Finite-Ohnesorge Coefficients* shows what
#     that looks like numerically.
#
# ### Two regimes: oscillatory and aperiodic
#
# For fixed ``l``, the character of the two dominant roots depends on
# ``\alpha^2``:
#
# - **Large ``\alpha^2``** (low viscosity, large drop): the two dominant
#   roots are complex conjugates, ``\sigma=\gamma\pm i\omega_d`` -- a damped
#   oscillation.
# - **Small ``\alpha^2``** (high viscosity, small drop): both dominant roots
#   are real and positive -- two aperiodic decaying modes, overdamped, no
#   oscillation at all.
#
# ![Root locus of the two dominant roots as viscosity increases: a complex
# conjugate pair sweeps in from the imaginary axis, meets on the real axis at
# a critical Ohnesorge number, and splits into two real
# roots.](../assets/reid_root_locus.png)
#
# *The two dominant roots in the complex ``\sigma`` plane, for ``l=2``, as
# ``\mathrm{Oh}`` increases along each branch. At low viscosity they are a
# conjugate pair close to ``\pm i\sigma_0`` -- almost pure oscillation. As
# viscosity rises the pair migrates right and inward, meets the real axis at
# the critical point, and separates into a fast-decaying root and a slow
# creeping one.*
#
# The transition happens at a critical ``\alpha^2``, determined numerically
# by Chandrasekhar. For ``l=2`` it sits at ``\sigma_{2;0}R^2/\nu = 3.69``
# with ``\sigma_{2;\nu}/\sigma_{2;0} = 0.968`` at the transition. Since
# ``\alpha^2 = \sqrt{l(l-1)(l+2)}/\mathrm{Oh}``, that is a critical
# Ohnesorge number ``\mathrm{Oh}_c = \sqrt{8}/3.69 \approx 0.766`` for
# ``l=2``. The companion page *Finite-Ohnesorge Coefficients* recovers both of
# Chandrasekhar's numbers to three digits from DropSolver's own root finder.
#
# !!! note "What the critical point means for water"
#     Converting the critical point into a critical *radius* is a useful
#     sanity check on scale. For water in air (``T_1=74\;\mathrm{dyn/cm}``,
#     ``\nu=0.01\;\mathrm{cm^2/s}``), substituting into
#     ``\alpha^2=\sigma_{2;0}R^2/\nu=3.69`` gives
#     ``R_c \approx 23\,\mathrm{nm}``, far below the continuum limit. Water
#     drops at any observable size are therefore deeply *underdamped*
#     (``\mathrm{Oh}\sim10^{-3}`` at millimetric scale), and the aperiodic
#     regime is reachable only with a genuinely viscous fluid -- such as the
#     ``\mathrm{Oh}_0\sim57`` shear-thinning fluid that motivates
#     *Finite-Ohnesorge Coefficients*.

# ## Section 9: The small-viscosity (Lamb) limit
#
# One limit is available in closed form: low viscosity, ``q\to\infty``.
#
# Standard large-argument Bessel asymptotics, taken here as given, put
# ``Q_{l+1/2}(q)/q\to 0`` between the poles of ``Q_{l+1/2}`` as
# ``q\to\infty``. Setting ``Q\to0`` collapses the bracket on the right of
# the characteristic equation to exactly ``l+(l+1)=2l+1``, for symbolic
# ``l``:
# ```math
# \frac{\alpha^4}{q^4} + 1 = \frac{2(l-1)(2l+1)}{q^2}.
# ```

@variables l q alpha Q                                                                         #src
characteristic_lhs = alpha^4 / q^4 + 1                                                          #src
characteristic_rhs_full = (2 * (l - 1) / q^2) * (l + (l + 1) * (q - 2 * l * Q) / (q - 2 * Q))    #src
characteristic_rhs_Q0 = substitute(characteristic_rhs_full, Dict(Q => 0))                       #src
@assert symbolic_zero(characteristic_rhs_Q0 - 2 * (l - 1) * (2l + 1) / q^2)                      #src

# Multiplying through by ``q^4`` turns this into a quadratic in
# ``q^2``:
# ```math
# q^4 - 2(l-1)(2l+1)\,q^2 + \alpha^4 = 0.
# ```
# Solved by the quadratic formula,
# ``q^2 = (l-1)(2l+1) \pm \sqrt{(l-1)^2(2l+1)^2-\alpha^4}``. For
# ``\alpha \gg 1`` (the low-viscosity regime this limit describes), the
# ``\alpha^4`` term dominates under the square root, giving
# ``q^2 \approx (l-1)(2l+1) \pm i\alpha^2`` to leading order -- and since
# ``\sigma = q^2\nu/R^2`` and ``\alpha^2\nu/R^2=\sigma_{l;0}`` by
# definition, this is exactly **Lamb's classical result**:
# ```math
# \sigma_{l;\nu} = (l-1)(2l+1)\,\frac{\nu}{R^2} \pm i\,\sigma_{l;0}.
# ```
# This is checked quantitatively, not merely asymptotically: for
# ``l=2,3,5,10`` the relative gap between the exact quadratic root and
# Lamb's leading-order formula shrinks monotonically as ``\alpha`` grows
# from ``10`` to ``10^4``. Lamb's formula is genuinely the
# ``\alpha\to\infty`` limit of the exact equation, not an unrelated
# approximation that happens to resemble it.

for l_val in (2, 3, 5, 10)                                                                                                       #src
    prev_err = Inf                                                                                                               #src
    for alpha_val in (10.0, 100.0, 1000.0, 10000.0)                                                                              #src
        exact_q2 = (l_val - 1) * (2l_val + 1) + im * sqrt(Complex(alpha_val^4 - (l_val - 1)^2 * (2l_val + 1)^2))                  #src
        lamb_q2 = (l_val - 1) * (2l_val + 1) + im * alpha_val^2                                                                  #src
        err = abs(exact_q2 - lamb_q2) / abs(lamb_q2)                                                                             #src
        @assert err < prev_err || err < 1e-6                                                                                     #src
        prev_err = err                                                                                                           #src
    end                                                                                                                          #src
end                                                                                                                              #src

# ### Cross-check against the running solver
#
# Everything above is algebra; this part compares it against numbers. In the
# small-Oh regime the solver takes its per-mode restoring and damping
# coefficients straight from Lamb -- ``\omega_l^2=l(l-1)(l+2)`` and
# ``2\lambda_l=2\,\mathrm{Oh}\,(l-1)(2l+1)`` -- so a small-Oh run should decay
# at the rate ``\mathrm{Oh}\,(l-1)(2l+1)``.
#
# The measurement itself is just the slope of the log-amplitude between the
# first and last saved frame:

extract_decay_rate(times, A_l) = -log(abs(A_l[end]) / abs(A_l[1])) / (times[end] - times[1])
nothing # hide

# Running the solver on a single excited mode and comparing:
#
# | ``\mathrm{Oh}`` | ``l`` | Lamb ``\gamma`` | measured ``\gamma`` | relative error |
# |--:|--:|--:|--:|--:|
# | 0.02 | 2 | 0.100 | 0.1031 | 3.1% |
# | 0.02 | 3 | 0.280 | 0.2853 | 1.9% |
# | 0.05 | 2 | 0.250 | 0.2533 | 1.3% |
#
# Agreement is better than 5% in every case, so the theory on this page and
# the code that runs agree in the regime where both apply.

for (Oh_val, l_val) in ((0.02, 2), (0.02, 3), (0.05, 2))                                       #src
    gamma_lamb = (l_val - 1) * (2l_val + 1) * Oh_val                                           #src
    M = l_val                                                                                  #src
    theta_vec = make_theta_vec(M)                                                              #src
    precomp = precompute_integrals(NaN, M)[1]                                                  #src
    sigma0 = sqrt(l_val * (l_val - 1) * (l_val + 2))                                           #src
    dt_osc = 2 * pi / (sigma0 * 40)                                                            #src
    cfg = SimConstants(M, M + 1, Oh_val, 1e-6, theta_vec, precomp, dt_osc)                     #src
    init = DropState(M)                                                                        #src
    init.A[l_val] = 0.05                                                                       #src
    init.z = 2.0                                                                               #src
    init.dt = dt_osc                                                                           #src
    T_period = 2 * pi / sigma0                                                                 #src
    times, states = solve_drop!(cfg, OBParams(), init;                                         #src
        t_end=6 * T_period, save_every=T_period / 50, dt_init=dt_osc)                          #src
    gamma_sim = extract_decay_rate(times, [s.A[l_val] for s in states])                         #src
    err = abs(gamma_sim - gamma_lamb) / gamma_lamb                                             #src
    @assert err < 0.05                                                                         #src
end                                                                                            #src

# ## Section 10: How this is used downstream
#
# ### Connection to Molaček & Bush (2012)
#
# In their quasi-static model of drop impact, Molaček & Bush parametrize the
# kinetic energy and viscous dissipation of each surface harmonic ``m``
# (their notation for our ``l``) with two coefficients ``A_m`` and ``D_m``,
# appearing in their Eq. 31:
# ```math
# \mathrm{K.E.} = \pi\rho R_0^5 \dot{B}^2 \sum_m A_m \frac{2b_m^2}{m(2m+1)},
# \qquad
# D = 8\pi\mu R_0^3 \dot{B}^2 \sum_m D_m \frac{m}{2m+1}b_m^2.
# ```
# Those coefficients are *defined* so that the two roots of the quadratic
# ```math
# A_m b^2 - 2a D_m b + 1 = 0, \qquad a \equiv \mathrm{Oh}\sqrt{m(m-1)(m+2)},
# ```
# reproduce the two dominant roots of the characteristic equation derived
# above, at harmonic order ``l=m``. The variable identification is
# ``q^2 = (b/a)\,\alpha^2``, under which ``W(b/a) = Q_{l+1/2}(q)/q``.
#
# So ``A_m`` and ``D_m`` are not a separate theory: they are the *result* of
# solving Reid's problem at each order, compressed into two numbers per mode
# by fitting a quadratic to the dominant eigenvalue pair. The recipe is
#
# 1. for each ``m`` and ``\mathrm{Oh}``, evaluate
#    ``\alpha^2=\mathrm{Oh}^{-1}\sqrt{m(m-1)(m+2)}``;
# 2. solve the characteristic equation numerically for the two roots with
#    the smallest ``\mathrm{Re}(\sigma)``;
# 3. read off ``A_m`` and ``D_m`` by matching the quadratic -- sum of roots
#    ``=2aD_m/A_m``, product of roots ``=1/A_m`` (Vieta);
# 4. tabulate or fit ``A_m(\mathrm{Oh})``, ``D_m(\mathrm{Oh})`` for use in
#    the Lagrangian equation of motion.
#
# **Why discard the higher roots?** The quasi-static assumption restricts
# the drop shape to a one-parameter family of equilibrium shapes. Only the
# fundamental surface mode at each harmonic order is representable in that
# family -- the higher radial overtones of Section 8 have internal nodal
# surfaces that simply do not exist in it, and in any case they decay on
# timescales far shorter than ``R^2/\nu``, well separated from the impact
# timescale.
#
# Steps 2 and 3 are exactly what *Finite-Ohnesorge Coefficients* carries out,
# in a slightly different (and, for the solver, more convenient) gauge.
#
# ### Summary of key equations
#
# | Equation | Physical content |
# |:--|:--|
# | ``\sigma_{l;0}^2 = l(l-1)(l+2)\,T_1/(\rho R^3)`` | inviscid capillary frequency |
# | ``q^2 = \sigma R^2/\nu`` | complex decay rate, scaled |
# | ``\alpha^2 = \sigma_{l;0}R^2/\nu`` | inviscid frequency, scaled |
# | ``Q_{l+1/2}(q) = J_{l+3/2}(q)/J_{l+1/2}(q)`` | the one Bessel ratio |
# | ``\alpha^4/q^4 + 1 = (2(l-1)/q^2)\left[l + (l+1)\frac{q-2lQ}{q-2Q}\right]`` | Reid's characteristic equation |
# | ``\sigma_{l;\nu} = (l-1)(2l+1)\nu/R^2 \pm i\,\sigma_{l;0}`` | small-``\nu`` (Lamb) limit |
#
# ### What survives a change of rheology, and what does not
#
# This matters because every other derivation here starts by
# modifying this one, and it is easy to modify too little.
#
# **Rheology-agnostic** -- depends only on incompressibility and geometry,
# carries over unchanged: the poloidal decomposition (Section 2.3), the
# pressure field satisfying Laplace's equation (Section 3), and BC1, the
# kinematic condition (Section 6).
#
# **Newtonian-specific** -- must be redone for any other constitutive law:
# the ``\nu\nabla^2\bm u`` term in the momentum equation is the obvious one,
# flagged in Section 4. But it is *not the only one*. BC2 and BC3 both also
# assume a constant ``\mu`` multiplying a linear strain rate
# (``\tau_{r\theta}=\mu[\cdots]`` and ``-2\mu\,\partial u_r/\partial r``
# respectively). A shear-thinning or viscoelastic model changes the momentum
# equation and both stress boundary conditions -- Oldroyd-B and
# Carreau-Yasuda each have to redo that work, not merely swap a coefficient.
#
# ### Where to go next
#
# Starting from the linearized Navier-Stokes equations for a perturbed
# spherical drop, this page derived -- symbolically, at every step where
# that was feasible -- the exact transcendental characteristic equation
# governing damped drop oscillations, and confirmed it reduces to Lamb's
# classical small-viscosity formula in the appropriate limit. That exact
# equation, not Lamb's asymptotic approximation to it, is the physically
# correct starting point for any drop whose Ohnesorge number isn't small --
# which includes every shear-thinning fluid the Carreau-Yasuda
# extension was built to handle (see *Carreau-Yasuda: Multi-Mode*),
# since shear-thinning can swing the effective Ohnesorge number across
# orders of magnitude within a single impact.
#
# See *Finite-Ohnesorge Coefficients* for what comes next: the numerical
# machinery that solves this transcendental equation robustly at finite Oh,
# the ``(\lambda_l(\mathrm{Oh}), \omega_l^2(\mathrm{Oh}))`` parametrization
# wired into `julia/src/reid.jl`, and a live cross-check against
# the running solver.
