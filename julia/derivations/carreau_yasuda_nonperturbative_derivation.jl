# # Non-Perturbative Carreau-Yasuda: Exact Effective-Oh from Instantaneous Shear
#
# ## Where the perturbative correction runs out
#
# `julia/src/st_extension.jl` builds its shear-thinning damping correction by
# Taylor-expanding the Carreau-Yasuda viscosity law to **first order in**
# ``\varepsilon_{ST}``:
#
# ```math
# \frac{\mu_{\mathrm{eff}}}{\mu_0}
#  = \bigl[1+(\lambda_c\dot\gamma)^{a}\bigr]^{-\varepsilon_{ST}}
#  \;\approx\; 1-\varepsilon_{ST}\,(\lambda_c\dot\gamma)^{a}.
# ```
#
# That requires ``\varepsilon_{ST}\ll1`` *and* ``(\lambda_c\dot\gamma)^a``
# not too large. Neither holds for the shear-thinning fluid used for
# validation here: ``\varepsilon_{ST}\approx0.9996`` -- not small at all --
# and the dimensionless ``\lambda_c\approx3.05\times10^{4}``, so
# ``(\lambda_c\dot\gamma)^a`` exceeds 1 for essentially any nonzero
# deformation rate the solver will see.
#
# The breakdown is qualitative rather than a loss of accuracy. The factor
# ``1-\varepsilon_{ST}(\lambda_c\dot\gamma)^{a}`` has no floor once
# ``\varepsilon_{ST}\approx1``: it goes **negative**, which means the damping
# term injects energy. §3 shows it is already negative at the smallest shear
# rate worth testing.
#
# ## What replaces it
#
# The effective Ohnesorge number for mode ``l`` becomes a function of the
# current shear rate rather than a fixed constant. Crucially, evaluating it
# at the instantaneous state is **exact** for this constitutive law, not a
# convenient approximation: Carreau-Yasuda has no relaxation time, so
# viscosity is an algebraic function of the instantaneous local shear rate
# rather than a memory kernel. There is nothing to average over -- which is
# what the perturbative derivation's phase-averaging step was doing.
#
# Both the damping ``\lambda_l`` and the restoring frequency
# ``\omega_l^2`` must then be re-evaluated at
# ``\mathrm{Oh}_{\mathrm{eff},l}`` through Reid's exact relations -- not just
# the damping term the way the perturbative correction does. The reason is in
# §4: ``\omega_l^2`` equals the inviscid ``l(l-1)(l+2)`` only in the
# small-``\mathrm{Oh}`` limit, and this fluid's rest-state
# ``\mathrm{Oh}_0\approx57`` is nowhere near small.
#
# ## Scope
#
# Two limitations carry through this page. True cross-mode coupling is not
# included: the characteristic shear rate uses only mode ``l``'s own velocity
# field, the same simplification the perturbative code makes through its
# per-mode ``\Gamma_l^{(a)}``, and the one
# `carreau_yasuda_multimode_derivation.jl` removes. And the coefficients
# derived here are free-decay coefficients, so substituting them into the
# forced, contact-coupled ODE is an additional approximation -- one shared
# with the `:reid` viscous model generally.

using DropSolver

# ## 1. The single-mode strain-rate field
#
# For a single active mode ``l`` with surface velocity ``\dot A_l
# P_l(\cos\theta)``, the interior velocity field is exactly the potential
# flow ``\phi=r^{l}P_l(x)\Phi_l`` with ``x=\cos\theta`` and
# ``\Phi_l=\dot A_l/l``, which is what reproduces
# ``u_r(1,\theta)=\dot A_l P_l(\cos\theta)``. This is the same inviscid
# mode-shape simplification the perturbative ``\Gamma_l^{(a)}`` machinery
# uses, not a new approximation.
#
# The standard axisymmetric spherical strain-rate formulas, after using
# Legendre's equation ``(1-x^2)P_l''-2xP_l'+l(l+1)P_l=0`` to eliminate
# ``P_l''``, collapse to (with ``X=P_l(x)``, ``X'=P_l'(x)``, and an overall
# prefactor ``r^{l-2}\Phi_l``)
#
# ```math
# e_{rr}=l(l-1)X,\quad
# e_{\theta\theta}=xX'-l^2X,\quad
# e_{\varphi\varphi}=lX-xX',\quad
# e_{r\theta}=-(l-1)\sin\theta\,X'.
# ```
#
# The first three sum to zero identically -- incompressibility -- confirmed
# numerically at ``l=2,3,5,8`` and several ``x``, to within ``10^{-10}``.
# Were that to fail, the field would not be divergence-free, and no shear
# rate derived from it would describe an incompressible flow.

function legendre_P_dP(l::Int, x::Float64)                                    #src
    l == 0 && return 1.0, 0.0                                                 #src
    Pm1, P = 1.0, x                                                           #src
    for n in 1:l-1                                                            #src
        P, Pm1 = ((2n + 1) * x * P - n * Pm1) / (n + 1), P                     #src
    end                                                                       #src
    Plm1 = l == 1 ? 1.0 : begin                                               #src
        Qm1, Q = 1.0, x                                                       #src
        for n in 1:l-2                                                        #src
            Q, Qm1 = ((2n + 1) * x * Q - n * Qm1) / (n + 1), Q                 #src
        end                                                                   #src
        Q                                                                     #src
    end                                                                       #src
    dP = l * (Plm1 - x * P) / (1 - x^2)                                       #src
    P, dP                                                                     #src
end                                                                           #src

function strain_components(l::Int, x::Float64)                                #src
    X, Xp = legendre_P_dP(l, x)                                               #src
    e_rr = l * (l - 1) * X                                                    #src
    e_thth = x * Xp - l^2 * X                                                 #src
    e_phph = l * X - x * Xp                                                   #src
    e_rth = -(l - 1) * sqrt(1 - x^2) * Xp                                     #src
    e_rr, e_thth, e_phph, e_rth                                               #src
end                                                                           #src

for l in (2, 3, 5, 8), x in (-0.7, -0.1, 0.3, 0.9)                            #src
    e_rr, e_thth, e_phph, _ = strain_components(l, x)                         #src
    @assert abs(e_rr + e_thth + e_phph) < 1e-10                               #src
end                                                                           #src
println("ASSERTION 1 OK: incompressibility holds for l=2,3,5,8 to <1e-10")     #src

# ## 2. The characteristic shear rate ``K_l``, exactly
#
# Squaring and contracting,
#
# ```math
# \dot\gamma^2(r,x)
#  = 2\bigl(e_{rr}^2+e_{\theta\theta}^2+e_{\varphi\varphi}^2+2e_{r\theta}^2\bigr)
#  = r^{2l-4}\,\Phi_l^2\,g_l(x),
# ```
#
# and ``g_l(x)`` is an *exact polynomial* in ``x``: Legendre polynomials are
# polynomials, and the ``\sin^2\theta=1-x^2`` factor carried by
# ``e_{r\theta}^2`` keeps it that way. So the volume-RMS characteristic shear
# rate
#
# ```math
# \dot\gamma_{\mathrm{char},l}(t)=K_l\,|\dot A_l(t)|,
# \qquad
# K_l^2=\frac{\langle\dot\gamma^2\rangle_{\mathrm{volume}}}{\dot A_l^2},
# ```
#
# can be evaluated by term-by-term integration in exact rational arithmetic:
# the radial part ``\int_0^1 r^{2l-2}\,dr=1/(2l-1)`` and the angular part
# ``\int_{-1}^{1}g_l(x)\,dx`` are both closed-form. No quadrature is
# involved, and no floating-point error enters:
#
# | ``l`` | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
# |:--|:--|:--|:--|:--|:--|:--|:--|:--|:--|
# | ``K_l^2`` | ``3`` | ``4`` | ``9/2`` | ``24/5`` | ``5`` | ``36/7`` | ``21/4`` | ``16/3`` | ``27/5`` |
#
# The whole table is one closed form,
#
# ```math
# \boxed{\;K_l^2 \;=\; \frac{6(l-1)}{l}\;}
# ```
#
# which matches the polynomial integration for every ``l`` from 2 to 40. The
# shape is the expected one: ``K_l^2`` rises monotonically and saturates at
# 6, so higher modes shear the fluid somewhat harder per unit surface
# velocity but not without bound, and ``K_1=0`` -- the ``l=1`` translation
# mode, which deforms nothing, produces no shear. Independent numerical
# quadrature reproduces the same values to machine precision.

## Polynomial arithmetic on Vector{Rational{BigInt}}: coeffs[i] multiplies x^(i-1). #src
_padd(a, b) = begin                                                           #src
    n = max(length(a), length(b))                                             #src
    [(i <= length(a) ? a[i] : 0 // 1) + (i <= length(b) ? b[i] : 0 // 1) for i in 1:n] #src
end                                                                           #src
_pscale(a, c) = [c * ai for ai in a]                                          #src
_pmulx(a) = vcat(0 // 1, a)                                                   #src
_pmul(a, b) = begin                                                           #src
    r = zeros(Rational{BigInt}, length(a) + length(b) - 1)                    #src
    for i in eachindex(a), j in eachindex(b)                                  #src
        r[i+j-1] += a[i] * b[j]                                               #src
    end                                                                       #src
    r                                                                         #src
end                                                                           #src
_ptrim(a) = begin                                                             #src
    n = length(a)                                                             #src
    while n > 1 && a[n] == 0 // 1                                             #src
        n -= 1                                                                #src
    end                                                                       #src
    a[1:n]                                                                    #src
end                                                                           #src

## BigInt rationals throughout: the Bonnet recursion's coefficients overflow  #src
## Int64 well before l=40, and this table is asserted exactly, not to a tol.  #src
function legendre_poly(l::Int)                                                #src
    P0 = [one(Rational{BigInt})]                                              #src
    l == 0 && return P0                                                       #src
    P1 = [zero(Rational{BigInt}), one(Rational{BigInt})]                      #src
    l == 1 && return P1                                                       #src
    Pm1, P = P0, P1                                                           #src
    for n in 1:l-1                                                            #src
        term1 = _pscale(_pmulx(P), (2n + 1) // (n + 1))                        #src
        term2 = _pscale(Pm1, n // (n + 1))                                    #src
        P, Pm1 = _ptrim(_padd(term1, _pscale(term2, -1 // 1))), P              #src
    end                                                                       #src
    P                                                                         #src
end                                                                           #src
polyderiv(a) = length(a) == 1 ? [0 // 1] : [i * a[i+1] for i in 1:length(a)-1] #src

function poly_int_m1_1(a)                                                     #src
    s = 0 // 1                                                                #src
    for (idx, c) in enumerate(a)                                              #src
        n = idx - 1                                                           #src
        iseven(n) && (s += c * 2 // (n + 1))                                  #src
    end                                                                       #src
    s                                                                         #src
end                                                                           #src

function K_l_squared_exact(l::Int)                                            #src
    X = legendre_poly(l)                                                      #src
    Xp = polyderiv(X)                                                         #src
    e_rr = _pscale(X, Rational{BigInt}(l * (l - 1)))                          #src
    e_thth = _padd(_pmulx(Xp), _pscale(X, Rational{BigInt}(-l^2)))            #src
    e_phph = _padd(_pscale(X, Rational{BigInt}(l)), _pscale(_pmulx(Xp), -1 // 1)) #src
    Xp2 = _pmul(Xp, Xp)                                                       #src
    one_minus_x2_Xp2 = _padd(Xp2, _pscale(_pmul([0 // 1, 0 // 1, 1 // 1], Xp2), -1 // 1)) #src
    g = _padd(_pmul(e_rr, e_rr), _pmul(e_thth, e_thth))                       #src
    g = _padd(g, _pmul(e_phph, e_phph))                                       #src
    g = _padd(g, _pscale(one_minus_x2_Xp2, Rational{BigInt}(2 * (l - 1)^2)))   #src
    g = _pscale(g, 2 // 1)                                                    #src
    (3 // (2 * (2l - 1))) * poly_int_m1_1(g) // (l^2)                         #src
end                                                                           #src

_K_TABLE_EXPECTED = Dict(2 => 3 // 1, 3 => 4 // 1, 4 => 9 // 2, 5 => 24 // 5,  #src
    6 => 5 // 1, 7 => 36 // 7, 8 => 21 // 4, 9 => 16 // 3, 10 => 27 // 5)      #src
for (l, expected) in _K_TABLE_EXPECTED                                        #src
    got = K_l_squared_exact(l)                                                #src
    @assert got == expected "l=$l: got $got, expected $expected"              #src
end                                                                           #src
@assert issorted([K_l_squared_exact(l) for l in 2:10])                        #src
## The closed form K_l^2 = 6(l-1)/l, checked well beyond the tabulated range. #src
for l in 2:40                                                                 #src
    @assert K_l_squared_exact(l) == 6 * (l - 1) // l "closed form fails at l=$l" #src
end                                                                           #src
println("ASSERTION 2 OK: K_l^2 = 6(l-1)/l exactly, for l=2..40")               #src

K_l(l::Int) = sqrt(Float64(K_l_squared_exact(l)))                             #src

# ## 3. ``\mathrm{Oh}_{\mathrm{eff},l}(t)``, and the sign of the perturbative factor
#
# The replacement is the exact law, untruncated:

Oh_eff(Oh0, lambda_c, a, eps_ST, gammadot_char) =
    Oh0 * (1 + (lambda_c * gammadot_char)^a)^(-eps_ST)

# with ``\dot\gamma_{\mathrm{char},l}(t)=K_l|\dot A_l(t)|`` from §2, and
# ``\varepsilon_{ST}=(1-n)/a`` per `julia/src/types.jl`'s `STParams`
# docstring.
#
# **The perturbative factor is a special case, not a different model.** Taylor-expanding
# ``[1+X]^{-\varepsilon_{ST}}`` to first order in ``\varepsilon_{ST}`` and
# then linearizing again in ``X=(\lambda_c\dot\gamma)^a`` gives exactly the
# ``1-\varepsilon_{ST}(\lambda_c\dot\gamma_{\mathrm{char},l})^a`` factor
# already multiplying ``D_2`` in `st_extension.jl`, with that file's per-mode
# geometric factor `Gamma_eff` playing the role of ``K_l^a``. Driving
# ``\varepsilon_{ST}`` through ``10^{-2},10^{-3},10^{-4}`` at fixed modest
# shear, the gap shrinks monotonically and falls below ``10^{-6}``.
#
# **For this fluid, the two are not close.** At
# ``\mathrm{Oh}_0=57.4``, ``\lambda_c=30507``, ``a=0.7431``,
# ``\varepsilon_{ST}=0.99956``:
#
# | ``\dot\gamma_{\mathrm{char}}`` | ``10^{-4}`` | ``10^{-3}`` | ``10^{-2}`` | ``0.1`` | ``1`` |
# |:--|:--|:--|:--|:--|:--|
# | perturbative ``1-\varepsilon_{ST}X`` | ``-1.29`` | ``-11.7`` | ``-69.1`` | ``-387`` | ``-2148`` |
# | exact ``\mathrm{Oh}_{\mathrm{eff}}/\mathrm{Oh}_0`` | 0.304 | 0.0732 | 0.0141 | 0.00258 | 0.000467 |
#
# The perturbative multiplier is already negative at the *smallest* shear
# rate in the table, and grows more negative from there -- a damping
# coefficient with the wrong sign, i.e. energy injection. The exact ratio
# stays in ``(0,1]`` for every shear rate tested out to ``10^3``, and does so
# by construction: ``[1+(\lambda_c\dot\gamma)^a]^{-\varepsilon_{ST}}`` is
# strictly between 0 and 1 for any positive ``\varepsilon_{ST},\lambda_c,a,
# \dot\gamma``. The difference between the two rows is one of sign, not of
# accuracy.

let Oh0 = 1.0, lambda_c = 0.3, a = 2.0                                        #src
    prev_err = Inf                                                            #src
    for eps_ST in (0.01, 0.001, 0.0001)                                       #src
        gammadot = 0.05   # small enough that (lambda_c*gammadot)^a << 1 too   #src
        exact = Oh_eff(Oh0, lambda_c, a, eps_ST, gammadot)                    #src
        perturbative = Oh0 * (1 - eps_ST * (lambda_c * gammadot)^a)           #src
        err = abs(exact - perturbative) / Oh0                                 #src
        @assert err < prev_err                                                #src
        prev_err = err                                                        #src
    end                                                                       #src
    @assert prev_err < 1e-6                                                   #src
end                                                                           #src
println("ASSERTION 3 OK: the exact law reduces to st_extension.jl's perturbative factor as eps_ST -> 0") #src

let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956             #src
    @assert (1 - eps_ST * (lambda_c * 1e-4)^a) < 0   # already negative at the smallest rate tabulated #src
    for gammadot_char in (1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 1000.0)           #src
        r = Oh_eff(Oh0, lambda_c, a, eps_ST, gammadot_char) / Oh0             #src
        @assert 0 < r <= 1                                                    #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 4 OK: perturbative multiplier goes negative; exact ratio stays in (0,1]") #src

# ## 4. Feeding ``\mathrm{Oh}_{\mathrm{eff},l}`` through Reid's exact relations
#
# Both `julia/src/timestepper.jl`'s ``D_1[l]`` (restoring) and ``D_2[l]``
# (damping) must be replaced by
# ``\lambda_l(\mathrm{Oh}_{\mathrm{eff},l})`` and
# ``\omega_l^2(\mathrm{Oh}_{\mathrm{eff},l})`` from `reid_lambda_omega2` --
# not ``D_2`` alone, the way a damping-only correction would.
#
# The size of that difference at this fluid's rest ``\mathrm{Oh}_0=57.4``:
#
# | ``l`` | ``\omega_l^2`` exact | inviscid ``l(l-1)(l+2)`` | deviation |
# |:--|:--|:--|:--|
# | 2 | 7.48 | 8 | 6.5% |
# | 3 | 25.48 | 30 | 15.1% |
# | 5 | 102.89 | 140 | 26.5% |
# | 10 | 665.4 | 1080 | 38.4% |
#
# The deviation grows monotonically with ``l`` and is already 6.5% at the
# fundamental. A damping-only correction silently assumes every entry in that
# last column is zero.

let Oh0 = 57.4, devs = Float64[]                                              #src
    for l in (2, 3, 5, 10)                                                    #src
        lam, om2, resid = reid_lambda_omega2(Oh0, l)                          #src
        @assert resid < 1e-8                                                  #src
        om0sq = Float64(l * (l - 1) * (l + 2))                                #src
        push!(devs, abs(om2 - om0sq) / om0sq)                                 #src
    end                                                                       #src
    @assert issorted(devs)          # deviation grows monotonically with l     #src
    @assert devs[1] > 0.03          # even at l=2, not negligible              #src
    @assert devs[end] > 0.3         # by l=10, over 30% -- not a minor correction #src
end                                                                           #src
println("ASSERTION 5 OK: omega_l^2 deviates from the inviscid value by 6.5% (l=2) to 38% (l=10)") #src

# ## 5. It integrates
#
# The final check is dynamic rather than algebraic: run a free ``l=2``
# oscillation at the real fluid's parameters -- ``\mathrm{Oh}_0=57.4``,
# ``\lambda_c=30507``, ``a=0.7431``, ``\varepsilon_{ST}=0.99956``, exactly
# where the perturbative correction has already gone negative -- recomputing
# ``\mathrm{Oh}_{\mathrm{eff},l}`` from the instantaneous ``\dot A_2`` at
# every step, and integrating semi-implicitly over six nominal periods.
#
# Starting from ``A_2=0.05``, the amplitude never exceeds its initial value
# and ends at ``3.9\times10^{-4}``: bounded, finite throughout, and decaying
# overall -- which is what the negative perturbative multiplier at these
# parameters cannot deliver.

function run_exact_st_oscillation(Oh0, lambda_c, a, eps_ST, l; M=l, A_init=0.05, #src
    t_end_periods=6.0, Bo=1e-6, viscous=:lamb)                                #src
    Kl = K_l(l)                                                               #src
    theta_vec = make_theta_vec(M)                                             #src
    precomp = precompute_integrals(NaN, M)[1]                                 #src
    sigma0 = sqrt(Float64(l * (l - 1) * (l + 2)))                             #src
    dt = 2 * pi / (sigma0 * 40)                                               #src
    cfg = SimConstants(M, M + 1, Oh0, Bo, theta_vec, precomp, dt; viscous=:lamb) #src
    A = zeros(M)                                                              #src
    A[l] = A_init                                                             #src
    Adot = zeros(M)                                                           #src
    t = 0.0                                                                   #src
    T_period = 2 * pi / sigma0                                                #src
    t_end = t_end_periods * T_period                                          #src
    A_hist = Float64[A[l]]                                                    #src
    t_hist = Float64[t]                                                       #src
    while t < t_end                                                           #src
        gammadot_char = Kl * abs(Adot[l])                                     #src
        Oh_eff_l = Oh_eff(Oh0, lambda_c, a, eps_ST, gammadot_char)            #src
        lam, om2, resid = viscous === :reid ? reid_lambda_omega2(Oh_eff_l, l) : #src
                          ((Oh_eff_l * (l - 1) * (2l + 1)), Float64(l * (l - 1) * (l + 2)), 0.0) #src
        @assert isfinite(lam) && isfinite(om2) && om2 > 0                     #src
        Adot_new = Adot[l] + dt * (-2 * lam * Adot[l] - om2 * A[l])  # semi-implicit Euler #src
        A_new = A[l] + dt * Adot_new                                          #src
        A[l], Adot[l] = A_new, Adot_new                                       #src
        t += dt                                                               #src
        push!(A_hist, A[l])                                                   #src
        push!(t_hist, t)                                                      #src
    end                                                                       #src
    t_hist, A_hist                                                            #src
end                                                                           #src

let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, l = 2       #src
    t_hist, A_hist = run_exact_st_oscillation(Oh0, lambda_c, a, eps_ST, l)     #src
    @assert all(isfinite, A_hist)                                             #src
    @assert maximum(abs.(A_hist)) < 1.0            # bounded, no blow-up       #src
    @assert abs(A_hist[end]) < abs(A_hist[1])      # net decay over the run    #src
    @assert abs(A_hist[end]) < 1e-3                # pins the reported endpoint #src
end                                                                           #src
println("ASSERTION 6 OK: free oscillation at the real fluid's parameters stays bounded and decays") #src

# ## Summary and scope
#
# Mode ``l``'s viscosity becomes ``\mathrm{Oh}_{\mathrm{eff},l}(t)``,
# evaluated from that mode's own instantaneous shear rate ``K_l|\dot A_l(t)|``
# -- with ``K_l`` an exact rational geometric constant -- through the exact,
# untruncated Carreau-Yasuda law. Evaluating at the instantaneous state is
# exact rather than approximate because this constitutive law has no memory.
# That ``\mathrm{Oh}_{\mathrm{eff},l}(t)`` then replaces *both* the damping
# and the restoring coefficient through Reid's exact relations. For the
# validation fluid the perturbative correction is already unphysical at the
# smallest shear rates tested; the untruncated law stays bounded across the
# whole range, and reproduces the perturbative factor exactly in the
# double-small-parameter limit where that factor is valid.
#
# The corresponding residual and Jacobian live in
# `julia/src/st_exact_extension.jl` -- a separate code path from
# `st_extension.jl`'s, because the ODE structure now needs ``D_1`` to vary
# and not only ``D_2``. Since ``\mathrm{Oh}_{\mathrm{eff},l}`` changes every
# step, `reid_lambda_omega2`'s continuation solve is too slow to call
# directly; the implementation interpolates a per-mode table built once
# (`build_reid_table`), with `:lamb` at the exact
# ``\mathrm{Oh}_{\mathrm{eff},l}`` as the cheaper alternative.
#
# The single-mode restriction is itself the subject of
# `carreau_yasuda_multimode_derivation.jl`, which supersedes this page's
# closure while keeping everything above it intact.
