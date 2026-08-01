# # Multi-Mode-Coupled Carreau-Yasuda: Shear Rate Is Kinematic
#
# This is the model the validation pipeline runs. It supersedes the
# single-mode closure of `carreau_yasuda_nonperturbative_derivation.jl`.
# Two further ingredients are needed before it behaves on a real impact -- a
# dealiasing rule for truncation ringing at contact onset (§6) and a finite
# infinite-shear viscosity plateau (§7) -- and both are derived below. §8
# reports where the model stands against the sampled experiments.
#
# ## What the single-mode closure misses
#
# The single-mode derivation computes mode ``l``'s effective viscosity from
# mode ``l``'s own velocity alone: ``\dot\gamma_{\mathrm{char},l}=K_l|\dot
# A_l|``. That is inadequate, because shear rate is a *field* quantity.
# At any point inside the drop, the local strain rate is the superposition of
# every active mode's velocity field at that point, regardless of what
# excited those modes. During contact several modes are driven at once -- not
# just the dominant ``l=2`` -- and their combined local shear rate can be far
# larger than any single mode's own characteristic value. A model that
# evaluates each mode's viscosity in isolation cannot see this at all.
#
# ## The replacement, in one line
#
# Mode ``l`` still needs a scalar effective viscosity, but it is now obtained
# from a genuine multi-mode Rayleigh-dissipation argument. With the total
# strain field written as a linear superposition,
#
# ```math
# e_{ij}(r,\theta,t) \;=\; \sum_k \dot A_k(t)\, e^{(k)}_{ij}(r,\theta),
# \qquad
# S(r,\theta,t) \;=\; \sqrt{2\,e_{ij}e_{ij}},
# ```
#
# the total dissipation is ``D=\int \mu_{\mathrm{eff}}(S)\,S^2\,dV`` and the
# generalized force on mode ``l`` is ``Q_l=\partial D/\partial\dot A_l``.
# Rather than use ``Q_l`` as an absolute damping coefficient -- which would
# require reconciling it against Reid's ``l``-dependent damping
# normalization, and §2 below shows that reconciliation does not exist for
# this basis -- only its *structure* is used, to define an effective
# viscosity **ratio**:
#
# ```math
# \boxed{\;
# \mathrm{Oh}_{\mathrm{eff},l}(t)
#   \;=\; \mathrm{Oh}_0\,
#   \frac{\displaystyle\int |w_l|\;\frac{\mu_{\mathrm{eff}}(S)}{\mu_0}\,dV}
#        {\displaystyle\int |w_l|\,dV},
# \qquad
# w_l \;\equiv\; \frac{\partial S^2}{\partial \dot A_l}
#   \;=\; 2\,\bigl(e_{ij} : e^{(l)}_{ij}\bigr).
# \;}
# ```
#
# The weight ``w_l`` is the sensitivity of the local dissipation to mode
# ``l``'s own motion; it is nonzero for ``k\neq l`` precisely when several
# modes are active at once, which is the coupling the single-mode closure
# omits. Taking ``|w_l|`` keeps the average well defined regardless of sign.
#
# Two properties make this construction safe, and both are checked below.
# It reduces exactly to the single-mode ``K_l`` result when only mode ``l``
# is active. And it returns ``\mathrm{Oh}_{\mathrm{eff},l}=\mathrm{Oh}_0``
# identically in the Newtonian limit, for *any* weight ``w_l`` whatsoever --
# because a ratio of ``\mu_{\mathrm{eff}}/\mu_0`` is 1 everywhere when there
# is no thinning. That second property is what lets §2's normalization defect
# cancel instead of propagating.

using DropSolver

# ## 1. Strain fields superpose; shear rate does not
#
# For the potential-flow representation used here, mode ``l``'s strain
# components are exactly (with ``x=\cos\theta``, ``X=P_l(x)``, and the
# trivial radial factor ``r^{l-2}`` divided out)
#
# ```math
# e_{rr} = (l-1)X, \quad
# e_{\theta\theta} = \frac{xX'-l^2X}{l}, \quad
# e_{\varphi\varphi} = \frac{lX-xX'}{l}, \quad
# e_{r\theta} = -\frac{l-1}{l}\sqrt{1-x^2}\,X' .
# ```
#
# Velocity fields add linearly across modes, so these components add too --
# each carrying its own ``\dot A_k`` and its own ``r^{k-2}``. Confirmed
# directly: for a two-mode state ``\{\dot A_2=0.3,\ \dot A_4=0.5\}``, the
# assembled total field agrees with term-by-term superposition to machine
# precision at every component.
#
# The invariant ``S=\sqrt{2e_{ij}e_{ij}}`` emphatically does **not**
# superpose -- it is quadratic in the field, and that is the entire source of
# the cross-mode coupling below. Were `total_strain` to stop being linear,
# every ``w_l`` and every ``\mathrm{Oh}_{\mathrm{eff},l}`` downstream would
# be built on a field that is not the physical one.

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

## Mode l's strain-rate tensor, per unit Adot_l and per unit r^(l-2).         #src
function strain_basis(l::Int, x::Float64)                                     #src
    X, Xp = legendre_P_dP(l, x)                                               #src
    e_rr = (l - 1) * X                                                        #src
    e_thth = (x * Xp - l^2 * X) / l                                           #src
    e_phph = (l * X - x * Xp) / l                                             #src
    e_rth = -(l - 1) / l * sqrt(1 - x^2) * Xp                                 #src
    e_rr, e_thth, e_phph, e_rth                                               #src
end                                                                           #src

## Superposition over all currently active modes, each with its own r^(l-2).  #src
function total_strain(active_modes, x::Float64, r::Float64)                   #src
    e_rr = e_thth = e_phph = e_rth = 0.0                                      #src
    for (l, Adot_l) in active_modes                                           #src
        Adot_l == 0.0 && continue                                             #src
        b = strain_basis(l, x)                                                #src
        rp = r^(l - 2)                                                        #src
        e_rr += Adot_l * rp * b[1]                                            #src
        e_thth += Adot_l * rp * b[2]                                          #src
        e_phph += Adot_l * rp * b[3]                                          #src
        e_rth += Adot_l * rp * b[4]                                           #src
    end                                                                       #src
    e_rr, e_thth, e_phph, e_rth                                               #src
end                                                                           #src

shear_rate_S(e_rr, e_thth, e_phph, e_rth) =                                   #src
    sqrt(2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2))                    #src

let x = 0.4, r = 0.7                                                          #src
    modes = Dict(2 => 0.3, 4 => 0.5)                                          #src
    e1 = total_strain(modes, x, r)                                            #src
    e2 = strain_basis(2, x) .* (0.3 * r^0) .+ strain_basis(4, x) .* (0.5 * r^2) #src
    @assert all(abs.(collect(e1) .- collect(e2)) .< 1e-12)                    #src
end                                                                           #src
println("ASSERTION 1 OK: two-mode total strain == term-by-term superposition") #src

# ## 2. Raw dissipation of this basis is not Lamb's damping
#
# It is tempting to use the bulk dissipation integral of the potential-flow
# field as an absolute damping coefficient. It does not work. Define
# ``H(k,l)=\int e^{(k)}_{ij}e^{(l)}_{ij}\,dV``. Lamb's classical result says
# mode ``l``'s damping goes as ``(l-1)(2l+1)``. The two disagree badly:
#
# | ``l`` | 2 | 3 | 4 | 5 |
# |:--|:--|:--|:--|:--|
# | ``H(l,l)`` | 6.283 | 8.378 | 9.425 | 10.053 |
# | ``(l-1)(2l+1)`` | 5 | 14 | 27 | 44 |
# | ratio | 1.257 | 0.598 | 0.349 | 0.229 |
#
# ``H(l,l)`` does grow with ``l``, but nowhere near in proportion: the ratio
# falls by more than a factor of five across four modes.
#
# The physical reading is that damping of a freely oscillating drop is not
# bulk dissipation of a bare potential flow -- it is dominated by the viscous
# correction to that flow, which this basis does not contain. So the raw
# integral is never used on its own here. It appears only inside the ratio of
# the boxed formula above, where the same normalization sits in the numerator
# and the denominator and cancels identically. §4's Newtonian-limit check
# tests that the cancellation is exact rather than approximate.

function H_kl(k::Int, l::Int; n_x::Int=200)                                   #src
    x_nodes, x_wts = DropSolver.gauss_legendre_nodes(n_x, -1.0, 1.0)          #src
    s = 0.0                                                                   #src
    for (x, w) in zip(x_nodes, x_wts)                                         #src
        a = strain_basis(k, x)                                                #src
        b = strain_basis(l, x)                                                #src
        s += (a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + 2 * a[4] * b[4]) * w   #src
    end                                                                       #src
    2 * pi * (1 / (k + l - 1)) * s                                            #src
end                                                                           #src

let                                                                           #src
    Hs = [H_kl(l, l) for l in (2, 3, 4, 5)]                                   #src
    lamb_scaling = [(l - 1) * (2l + 1) for l in (2, 3, 4, 5)]                 #src
    ratios = Hs ./ lamb_scaling                                               #src
    @assert issorted(Hs)                                # H does grow with l...  #src
    @assert maximum(ratios) / minimum(ratios) > 2.0     # ...but not as Lamb's   #src
end                                                                           #src
println("ASSERTION 2 OK: raw dissipation grows with l but NOT as (l-1)(2l+1)") #src

# ## 3. The viscosity law, with a finite infinite-shear plateau
#
# The local viscosity ratio entering the boxed average is the complete
# Cross-model form,
#
# ```math
# \frac{\mu_{\mathrm{eff}}(S)}{\mu_0}
#  \;=\; \frac{\eta_\infty}{\eta_0}
#      + \left(1-\frac{\eta_\infty}{\eta_0}\right)
#        \bigl[1+(\lambda_c S)^{a}\bigr]^{-\varepsilon_{ST}} .
# ```
#
# Setting ``\eta_\infty/\eta_0=0`` (the default) recovers the plain
# Carreau-Yasuda expression, so no existing caller changes behaviour. §7
# explains why the plateau term had to be added at all.

mu_ratio(S, lambda_c, a, eps_ST, eta_inf_ratio=0.0) =
    eta_inf_ratio + (1 - eta_inf_ratio) * (1 + (lambda_c * S)^a)^(-eps_ST)

## The two properties the formula is relied on for: it reduces to the plain   #src
## law at eta_inf_ratio=0, and it can never fall below the plateau.           #src
let lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, floor_ratio = 4.4e-4    #src
    for S in (1e-6, 1e-3, 1.0, 1e3, 1e9)                                      #src
        @assert mu_ratio(S, lambda_c, a, eps_ST) ≈ (1 + (lambda_c * S)^a)^(-eps_ST) #src
        @assert mu_ratio(S, lambda_c, a, eps_ST, floor_ratio) >= floor_ratio  #src
        @assert mu_ratio(S, lambda_c, a, eps_ST, floor_ratio) <= 1.0          #src
    end                                                                       #src
    @assert mu_ratio(1e12, lambda_c, a, eps_ST) < floor_ratio  # the plain law has no floor #src
end                                                                           #src
println("ASSERTION 3 OK: mu_ratio reduces to the plain law and respects the eta_inf floor") #src

# ## 4. The coupled effective Ohnesorge number
#
# The boxed average is evaluated by Gauss-Legendre quadrature in both
# ``r\in[0,1]`` and ``x=\cos\theta\in[-1,1]``, with the ``r^2`` volume
# element included.
#
# **Reduction to the single-mode result.** With only ``l=2`` active at
# ``\dot A_2=0.05``, and the fitted validation fluid
# (``\mathrm{Oh}_0=57.4``, ``\lambda_c=30507``, ``a=0.7431``,
# ``\varepsilon_{ST}=0.99956``), the coupled formula gives
# ``\mathrm{Oh}_{\mathrm{eff}}=0.1644``, against ``0.1644`` from the
# already-validated single-mode ``K_l`` formula -- agreement to better than
# the 1% the check demands, and in fact to the digits shown.
#
# **Newtonian limit.** With ``\lambda_c=0`` there is no thinning, and
# ``\mathrm{Oh}_{\mathrm{eff},l}=\mathrm{Oh}_0`` to a relative ``10^{-10}``
# -- for a single active mode and for three simultaneously active modes
# alike. This is the check that §2's normalization defect really does cancel:
# if the weight ``w_l`` leaked into the answer anywhere other than as a
# ratio, this equality would fail.
#
# **The coupling itself.** Hold mode 2 fixed at ``\dot A_2=0.02`` -- a modest
# velocity, typical early in a low-``\mathrm{We}`` impact -- and switch on a
# higher mode alongside it:
#
# | ``\dot A_5`` | 0 | 0.02 | 0.05 | 0.10 |
# |:--|:--|:--|:--|:--|
# | ``\mathrm{Oh}_{\mathrm{eff},2}`` (coupled) | 0.3238 | 0.2346 | 0.1437 | 0.0886 |
# | ``\mathrm{Oh}_{\mathrm{eff},2}`` (self-only) | 0.3238 | 0.3238 | 0.3238 | 0.3238 |
#
# Mode 2's own effective viscosity falls by a factor of 3.7 while mode 2's
# own velocity never changes. That is the entire point: shear is kinematic,
# so mode 2 *feels* the thinning that mode 5's velocity field produces at the
# same points in the fluid. The self-only model is blind to the whole column.

function w_l_field(l::Int, active_modes, x::Float64, r::Float64)              #src
    e_tot = total_strain(active_modes, x, r)                                  #src
    e_l = strain_basis(l, x) .* r^(l - 2)                                     #src
    2 * (e_tot[1] * e_l[1] + e_tot[2] * e_l[2] + e_tot[3] * e_l[3] + 2 * e_tot[4] * e_l[4]) #src
end                                                                           #src

function oh_eff_coupled(l::Int, active_modes, Oh0::Float64, lambda_c::Float64, #src
    a::Float64, eps_ST::Float64; n_r::Int=24, n_x::Int=40)                     #src
    r_nodes, r_wts = DropSolver.gauss_legendre_nodes(n_r, 0.0, 1.0)           #src
    x_nodes, x_wts = DropSolver.gauss_legendre_nodes(n_x, -1.0, 1.0)          #src
    num = 0.0                                                                 #src
    den = 0.0                                                                 #src
    for (ri, rw) in zip(r_nodes, r_wts), (xi, xw) in zip(x_nodes, x_wts)      #src
        wl = w_l_field(l, active_modes, xi, ri)                               #src
        wl == 0.0 && continue                                                 #src
        S = shear_rate_S(total_strain(active_modes, xi, ri)...)               #src
        weight = abs(wl) * rw * xw * ri^2   # r^2 from the volume element     #src
        num += weight * mu_ratio(S, lambda_c, a, eps_ST)                      #src
        den += weight                                                         #src
    end                                                                       #src
    den == 0.0 ? Oh0 : Oh0 * num / den                                        #src
end                                                                           #src

let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, l = 2, Adot_l = 0.05 #src
    gammadot_char = characteristic_shear_K(l) * abs(Adot_l)                   #src
    single_mode = Oh0 * (1 + (lambda_c * gammadot_char)^a)^(-eps_ST)          #src
    coupled = oh_eff_coupled(l, Dict(l => Adot_l), Oh0, lambda_c, a, eps_ST)  #src
    @assert abs(coupled - single_mode) / single_mode < 0.01                   #src
end                                                                           #src
println("ASSERTION 4 OK: one active mode reproduces the single-mode Oh_eff to <1%") #src

let Oh0 = 3.0, l = 3                                                          #src
    for active_modes in (Dict(2 => 0.3), Dict(2 => 0.3, 4 => 0.6, 5 => -0.2))  #src
        Oh_eff = oh_eff_coupled(l, active_modes, Oh0, 0.0, 2.0, 0.5)   # lambda_c=0 -> no thinning #src
        @assert isapprox(Oh_eff, Oh0; rtol=1e-10)                             #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 5 OK: Newtonian limit gives Oh_eff_l = Oh_0 exactly, coupled or not") #src

## The coupling demonstration tabulated above.                                #src
let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, l = 2, Adot2 = 0.02 #src
    self_only = oh_eff_coupled(l, Dict(2 => Adot2), Oh0, lambda_c, a, eps_ST)  #src
    vals = [oh_eff_coupled(l, Dict(2 => Adot2, 5 => Ah), Oh0, lambda_c, a, eps_ST) #src
            for Ah in (0.0, 0.02, 0.05, 0.1)]                                 #src
    @assert isapprox(vals[1], self_only; rtol=1e-12)      # Adot_5=0 must be the self-only case #src
    @assert all(i -> vals[i] > vals[i+1], 1:length(vals)-1)  # more high-mode shear -> more thinning #src
    @assert vals[end] < 0.3 * self_only                   # a factor >3, not a rounding effect #src
end                                                                           #src
println("ASSERTION 6 OK: a co-excited mode 5 cuts mode 2's Oh_eff by >3x") #src

# ## 5. Truncation ringing at contact onset
#
# The coupling of §1-§4, used on its own, interacts badly with the spectral
# truncation at the instant contact begins. Against the 20 sampled
# experiments it puts the contact-time median relative error at **62%**,
# against **16.8%** for the self-only closure it replaces, with predicted
# contact times collapsing to near zero on most samples.
#
# The cause is visible frame by frame at contact onset. The moment
# contact begins, the shape boundary condition pins a single collocation
# point to the plane -- a spatial kink. A finite Legendre truncation cannot
# represent a kink smoothly, and the un-representable part is not spread
# evenly across modes: it concentrates almost entirely in the single highest
# *retained* mode ``l=M``, because in a finite linear expansion that is the
# only place left for everything the lower, smoother modes cannot capture.
# This is the Gibbs phenomenon, and it has nothing to do with fluid
# dynamics -- it is what truncating any basis at finite order does when asked
# to represent a non-smooth boundary condition.
#
# Why this hits contact time specifically: this fluid's dimensionless
# ``\lambda_c\approx3\times10^4`` is astronomically large. Feeding even a
# tiny, purely numerical ``\dot A_M`` into §1's superposed shear field
# manufactures an enormous ``(\lambda_c S)^a``, which crashes
# ``\mu_{\mathrm{eff}}/\mu_0`` -- and hence ``\mathrm{Oh}_{\mathrm{eff}}``
# for *every* mode, not only mode ``M`` -- toward zero. The result is a
# spurious, near-total loss of viscosity at the exact instant contact
# begins, before any physically resolved mode has moved at all.
#
# On a live `solve_drop!` run at ``\mathrm{We}=0.7649`` the signature is
# unmistakable, at the first frame contact registers:
#
# | truncation | ``|\dot A_M|`` | ``\max_{l<M}|\dot A_l|`` | ratio |
# |:--|:--|:--|:--|
# | ``M=12`` | 1.15 | ``6.1\times10^{-15}`` | ``1.9\times10^{14}`` |
# | ``M=16`` | 0.286 | ``1.9\times10^{-15}`` | ``1.5\times10^{14}`` |
#
# Fourteen orders of magnitude, at both truncations. This amplitude
# signature is what the detector of §6 keys on; the check below requires
# only a factor of 100.

let                                                                           #src
    OH0 = 57.371648873370795                                                  #src
    LAMBDA_C = 30507.34501244818                                              #src
    A_SHAPE = 0.7430524574330837                                              #src
    EPS_ST = 0.9995574839318364                                               #src
    BO = 0.012                                                                #src
    We = 0.7649                                                               #src
    for M in (12, 16)                                                         #src
        stx = STExactParams(M, OH0, LAMBDA_C, A_SHAPE, EPS_ST; viscous=:reid)  #src
        dt_max = make_dt_max(M)                                               #src
        theta_vec = make_theta_vec(M)                                         #src
        precomp = precompute_integrals(NaN, M)[1]                             #src
        cfg = SimConstants(M, M + 1, OH0, BO, theta_vec, precomp, dt_max)     #src
        init = DropState(M)                                                   #src
        init.z = 1.05                                                         #src
        init.v = -sqrt(We)                                                    #src
        init.dt = dt_max                                                      #src
        init.cp = 0                                                           #src
        times, states = solve_drop!(cfg, OBParams(), init; stx=stx,           #src
            t_end=0.15, save_every=dt_max / 4)                                #src
        first_c = findfirst(s -> s.cp > 0, states)                            #src
        Adot = states[first_c].Adot[2:end]                                    #src
        @assert abs(Adot[end]) / maximum(abs.(Adot[1:end-1])) > 100           #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 7 OK: at contact onset the top mode carries >100x every other mode (M=12 and M=16)") #src

# ## 6. Dealiasing: neither index nor amplitude alone is enough
#
# The obvious fix -- drop roughly the top 10% of modes *by index* from the
# coupling field -- removes the ringing at the cost of discarding real
# physics. On a low-``\mathrm{We}`` run (``\mathrm{We}=0.0158``: a gentle
# impact with a long, slow contact) mode ``M`` itself reaches an amplitude
# comparable to mode 2's. Excluding it purely because its index is near
# ``M`` roughly halves the predicted coefficient of restitution for that
# sample: **0.21** with the index-based filter, against **0.74** with
# dealiasing disabled entirely and a measured **0.82**.
#
# Amplitude ratio alone does not work either. Compare a candidate against
# every *other* mode with no index restriction, and a live mid-mode case --
# mode 2 sitting at ``10^{-15}`` floating-point noise, not exactly zero,
# alongside a genuinely large mode 8 -- gets mode 8 flagged as an outlier. A
# plain "other modes" comparison has no way to distinguish noise sitting next
# to real signal from everything being genuinely silent.
#
# The rule that works uses both signals, each for what it is actually good
# at. **Index** decides which modes may even be candidates: only the
# truncation boundary, the same ``\sim`` 10% margin, because that is still
# the only place ringing concentrates. **Amplitude** then decides among the
# candidates, compared against the largest amplitude among the *trusted*
# (non-candidate) modes, floored at `RINGING_NOISE_FLOOR` ``=10^{-9}`` so a
# fully quiescent trusted set does not make every nonzero candidate look
# infinitely large. A candidate is excluded only if it exceeds
# `OUTLIER_FACTOR` ``=20`` times that trusted amplitude.
#
# Three cases pin the rule down, and no simpler rule passes all three:
#
# | case | state | required verdict |
# |:--|:--|:--|
# | ringing (§5, ``M=12``) | trusted modes at ``10^{-15}``, mode ``M`` at 0.673 | exclude mode ``M`` |
# | real low-``\mathrm{We}`` | live trace, mode ``M`` at 0.0163 next to mode 2 at ``-0.0159`` | keep mode ``M`` |
# | mid-mode | mode 2 at ``10^{-15}``, mode 8 (a trusted index) at 0.286 | keep mode 8 |
#
# Every mode, including a masked-out one, still receives its own
# ``\mathrm{Oh}_{\mathrm{eff}}`` for its own dynamics. The mask controls only
# what counts as *real shear* when building the field everyone's thinning is
# averaged against.

## Case A: the §5 ringing signature.                                          #src
let                                                                           #src
    Adot_ringing = fill(1e-15, 11)                                            #src
    Adot_ringing[end] = 0.673                                                 #src
    mask_a = DropSolver._ringing_outlier_mask(Adot_ringing)                   #src
    @assert mask_a[end] == true                                               #src
    @assert all(mask_a[1:end-1] .== false)                                    #src
end                                                                           #src

## Case B: live low-We data (We=0.0158, M=12, frame i=60 of the trace) --     #src
## mode M genuinely active, comparable to mode 2.                             #src
let                                                                           #src
    Adot_real = [-0.0159, 0.00719, -0.0137, 0.0135, -0.00666, 0.00287,        #src
        -0.00241, 0.002, -0.00129, 0.000523, 0.0163]                          #src
    @assert all(DropSolver._ringing_outlier_mask(Adot_real) .== false)        #src
end                                                                           #src

## Case C: the mid-mode failure of a pure amplitude-ratio rule.               #src
let                                                                           #src
    Adot_midmode = zeros(11)                                                  #src
    Adot_midmode[1] = 1e-15                                                   #src
    Adot_midmode[7] = 0.286   # l=8, a trusted (non-candidate) index at M=12   #src
    @assert DropSolver._ringing_outlier_mask(Adot_midmode)[7] == false        #src
end                                                                           #src
println("ASSERTION 8 OK: the index+trusted-amplitude mask gets all three cases right") #src

# ## 7. The missing physical floor on viscosity
#
# With dealiasing in place the contact-time median error falls from 62% to
# **29%** -- better, but still above the 16.8% of the closure it replaces.
# Part of the remaining gap has an identifiable cause: on a live run,
# several ``\mathrm{Oh}_{\mathrm{eff}}`` values sit *below the fluid's own
# physical minimum*.
#
# Plain Carreau-Yasuda as coded,
# ``\mu_{\mathrm{eff}}/\mu_0=[1+(\lambda_c\dot\gamma)^a]^{-\varepsilon_{ST}}``,
# sends ``\mu_{\mathrm{eff}}/\mu_0\to0`` exactly as ``\dot\gamma\to\infty``.
# The real fluid's infinite-shear viscosity ``\eta_\infty=0.00373`` Pa·s is
# small next to ``\eta_0=8.43`` Pa·s, but it is not zero. No amount of shear
# can push this fluid below
# ``\mathrm{Oh}_0\,\eta_\infty/\eta_0\approx0.0254``. The plain formula has no
# such floor, and on a live ``\mathrm{We}=0.7649`` run it goes straight
# through it:
#
# | quantity | value |
# |:--|:--|
# | physical floor ``\mathrm{Oh}_0\,\eta_\infty/\eta_0`` | 0.0254 |
# | min ``\mathrm{Oh}_{\mathrm{eff}}`` observed, no floor term | 0.0110 |
# | min ``\mathrm{Oh}_{\mathrm{eff}}`` observed, with floor term | 0.0364 |
#
# Both numbers are computed from the *same* recorded ``\dot A`` trajectory,
# so the comparison isolates the formula change from any difference in the
# resulting dynamics. Under-predicting viscosity by a factor of two at the
# contact region is not a cosmetic error: it lets the contact-region
# deformation swing to amplitudes the real fluid could never reach,
# eventually large enough to violate the small-deformation assumption the
# entire linearized shape-mode model rests on.
#
# The fix is §3's plateau term, with
# ``\eta_\infty/\eta_0=1-\varepsilon_{ST}\approx4.4\times10^{-4}`` for this
# fluid.

let                                                                           #src
    OH0 = 57.371648873370795                                                  #src
    LAMBDA_C = 30507.34501244818                                              #src
    A_SHAPE = 0.7430524574330837                                              #src
    EPS_ST = 0.9995574839318364                                               #src
    ETA_INF_OVER_ETA_0 = 1 - EPS_ST   # eta_inf/eta_0 under this fluid's Cross-model mapping #src
    BO = 0.012                                                                #src
    We = 0.7649                                                               #src
    M = 12                                                                    #src
    physical_floor = OH0 * ETA_INF_OVER_ETA_0                                 #src
    stx_no_floor = STExactParams(M, OH0, LAMBDA_C, A_SHAPE, EPS_ST; viscous=:reid) #src
    dt_max = make_dt_max(M)                                                   #src
    theta_vec = make_theta_vec(M)                                             #src
    precomp = precompute_integrals(NaN, M)[1]                                 #src
    cfg = SimConstants(M, M + 1, OH0, BO, theta_vec, precomp, dt_max)         #src
    init = DropState(M)                                                       #src
    init.z = 1.05                                                             #src
    init.v = -sqrt(We)                                                        #src
    init.dt = dt_max                                                          #src
    init.cp = 0                                                               #src
    times, states = solve_drop!(cfg, OBParams(), init; stx=stx_no_floor,      #src
        t_end=0.25, save_every=dt_max / 4)                                    #src
    min_no_floor = Inf                                                        #src
    for s in states                                                           #src
        Adot_vec = s.Adot[2:end]                                              #src
        any(x -> x != 0.0, Adot_vec) || continue                              #src
        min_no_floor = min(min_no_floor, minimum(oh_eff_all_coupled(stx_no_floor, OH0, Adot_vec))) #src
    end                                                                       #src
    @assert min_no_floor < physical_floor  # failing means the floor-less formula never violates the floor here #src
    stx_floor = STExactParams(M, OH0, LAMBDA_C, A_SHAPE, EPS_ST;              #src
        viscous=:reid, eta_inf_ratio=ETA_INF_OVER_ETA_0)                      #src
    min_floor = Inf                                                           #src
    for s in states   # same states -> same Adot trajectory, only the formula changes #src
        Adot_vec = s.Adot[2:end]                                              #src
        any(x -> x != 0.0, Adot_vec) || continue                              #src
        min_floor = min(min_floor, minimum(oh_eff_all_coupled(stx_floor, OH0, Adot_vec))) #src
    end                                                                       #src
    @assert min_floor >= physical_floor - 1e-9                                #src
end                                                                           #src
println("ASSERTION 9 OK: without the plateau term a live run drives Oh_eff below the physical floor; with it, never") #src

# ## 8. Where this leaves the model
#
# **What is settled.** Mode ``l``'s effective viscosity is no longer a
# function of ``\dot A_l`` alone. It is a dissipation-weighted average of the
# true local viscosity ratio ``\mu_{\mathrm{eff}}(S)/\mu_0``, where ``S`` is
# the shear-rate invariant of the full superposed strain field of every
# active mode, weighted by mode ``l``'s own sensitivity to that field. The
# construction is exact kinematics, uses the exact (non-Taylor-expanded)
# Carreau-Yasuda law, reduces to the validated single-mode result when only
# one mode is active, and returns ``\mathrm{Oh}_0`` identically in the
# Newtonian limit despite §2's normalization defect. It required two further
# ingredients that are physics, not tuning: a dealiasing rule for
# contact-onset truncation ringing (§6), and a finite infinite-shear plateau
# (§7).
#
# **What is not settled.** Against the 20 sampled experiments, the
# contact-time median relative error reads:
#
# | model | contact-time median error |
# |:--|:--|
# | self-only single-mode closure (superseded) | **16.8%** |
# | multi-mode coupling, without §6 and §7 | **62%** |
# | + dealiasing + ``\eta_\infty`` plateau | **29%** |
#
# The multi-mode model has *not* recovered the single-mode closure's
# contact-time accuracy. It is better founded -- the self-only closure
# misrepresents the physics in the way §4's coupling table makes concrete,
# and each of the two additions above corrects a specific defect rather than
# fitting a number -- but on this metric it remains nearly twice as far off
# as the model it replaces, and that gap is unexplained.
#
# Two candidate explanations have not been separated. Contact time is
# dominated by the contact-region dynamics, which is exactly where the
# truncation kink lives and where the mask of §6 is most active; a dealiasing
# rule that is right on all three test cases may still be wrong on the
# details of which shear is real during contact. And the collapse of a
# spatially varying ``\eta(r,\theta,t)`` onto one scalar
# ``\mathrm{Oh}_{\mathrm{eff},l}(t)`` remains an uncontrolled effective-medium
# step -- quantified at roughly a factor of two in local viscosity across the
# drop in `carreauYasuda_firstprinciples_derivation.jl`, which is not small
# enough to dismiss. Neither of these is closed here.
