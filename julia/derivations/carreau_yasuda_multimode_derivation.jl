#!/usr/bin/env julia
# ==============================================================================
# Multi-Mode-Coupled Non-Perturbative Carreau-Yasuda: Shear Rate Is Kinematic
#
# What was wrong with julia/derivations/carreau_yasuda_nonperturbative_derivation.jl:
# it computed mode l's effective viscosity from mode l's OWN Adot_l alone
# (gammadot_char_l = K_l*|Adot_l|). But shear rate is a purely kinematic field
# quantity -- at any point in the drop, the LOCAL strain rate is the
# superposition of EVERY active mode's velocity field at that point, whether
# those modes were excited by contact, a body force, or anything else. During
# contact, several modes are excited simultaneously (not just the dominant
# l=2), and their combined LOCAL shear rate at a given point can be much
# larger than any single mode's own characteristic shear rate would suggest --
# this is invisible to a model that treats each mode's viscosity in isolation.
#
# What this derives: mode l's damping still needs an effective viscosity, but
# now obtained from a genuine multi-mode Rayleigh-dissipation argument. Define
#     D_total(t) = integral over the drop of mu_eff(S(r,theta,t)) * S(r,theta,t)^2 dV
# where S(r,theta,t) = sqrt(2 e_ij e_ij) is the TRUE local shear-rate invariant
# of the SUPERPOSED strain field e_ij(r,theta,t) = sum_k Adot_k(t) * e_ij^(k)(r,theta)
# (linear superposition -- exact for this potential-flow representation, same
# basis already used in the single-mode derivation). The generalized force on
# mode l is Q_l = dD_total/dAdot_l. Rather than use Q_l as an absolute damping
# coefficient (which would require reconciling it with Reid's L-DEPENDENT
# damping normalization -- a known subtlety, see the note in Section 2), we
# use its STRUCTURE to define a properly multi-mode-coupled effective
# viscosity RATIO:
#     Oh_eff_l(t) = Oh_0 * <mu_eff(S)/mu_0>_{w_l-weighted},
#     w_l(r,theta,t) = dS^2/dAdot_l = 2 * sum_k Adot_k(t) * g_kl(r,theta)
# (the "sensitivity" of mode l's own dissipation contribution to the LOCAL
# shear field -- nonzero even for k != l when multiple modes are active
# simultaneously), weighted by |w_l| so the average is always well-defined.
# This reduces EXACTLY to the single-mode K_l result when only mode l is
# active, and to Oh_eff_l=Oh_0 identically in the Newtonian limit (mu_eff=mu_0
# everywhere) REGARDLESS of any normalization quirk in w_l -- a ratio of
# mu_eff/mu_0 is 1 everywhere when there's no thinning, independent of the
# weight used to average it. This sidesteps a genuine, separately-important
# subtlety (Section 2): the RAW dissipation integral of this potential-flow
# basis does NOT reproduce Reid/Lamb's l-dependent damping coefficient (it is
# then discovered EXACTLY here, again, having been found once already earlier
# in this session and not fully reconciled at the time) -- but since we only
# ever use it as a RATIO, that defect cancels rather than propagating into the
# physics.
# ==============================================================================

using DropSolver

# ------------------------------------------------------------------------------
# Section 1: strain-rate basis functions per mode (reused from the single-mode
# derivation), and the TOTAL (superposed) shear-rate field for an arbitrary
# active-mode vector.
# ------------------------------------------------------------------------------

function legendre_P_dP(l::Int, x::Float64)
    l == 0 && return 1.0, 0.0
    Pm1, P = 1.0, x
    for n in 1:l-1
        P, Pm1 = ((2n + 1) * x * P - n * Pm1) / (n + 1), P
    end
    Plm1 = l == 1 ? 1.0 : begin
        Qm1, Q = 1.0, x
        for n in 1:l-2
            Q, Qm1 = ((2n + 1) * x * Q - n * Qm1) / (n + 1), Q
        end
        Q
    end
    dP = l * (Plm1 - x * P) / (1 - x^2)
    P, dP
end

"""
    strain_basis(l, x) -> (e_rr, e_thth, e_phph, e_rth)

Strain-rate tensor components for mode l's own potential-flow field, PER UNIT
`Adot_l` and PER UNIT `r^(l-2)` (the r-dependence is trivial -- a single power
per mode -- and factored out here so multi-mode sums can apply it per term).
Matches julia/derivations/carreau_yasuda_nonperturbative_derivation.jl Section 1,
divided through by the r^(l-2) prefactor there.
"""
function strain_basis(l::Int, x::Float64)
    X, Xp = legendre_P_dP(l, x)
    e_rr = (l - 1) * X
    e_thth = (x * Xp - l^2 * X) / l
    e_phph = (l * X - x * Xp) / l
    e_rth = -(l - 1) / l * sqrt(1 - x^2) * Xp
    e_rr, e_thth, e_phph, e_rth
end

"""
    total_strain(active_modes, x, r) -> (e_rr, e_thth, e_phph, e_rth)

`active_modes` is a `Dict{Int,Float64}` (or any iterable of `(l, Adot_l)`
pairs) of currently-active modes. Superposes each mode's contribution
(linear in `Adot_l`, each with its own `r^(l-2)` radial scaling).
"""
function total_strain(active_modes, x::Float64, r::Float64)
    e_rr = e_thth = e_phph = e_rth = 0.0
    for (l, Adot_l) in active_modes
        Adot_l == 0.0 && continue
        b = strain_basis(l, x)
        rp = r^(l - 2)
        e_rr += Adot_l * rp * b[1]
        e_thth += Adot_l * rp * b[2]
        e_phph += Adot_l * rp * b[3]
        e_rth += Adot_l * rp * b[4]
    end
    e_rr, e_thth, e_phph, e_rth
end

shear_rate_S(e_rr, e_thth, e_phph, e_rth) =
    sqrt(2 * (e_rr^2 + e_thth^2 + e_phph^2 + 2 * e_rth^2))

println("="^78)
println("Section 1: total (superposed) shear rate from multiple active modes")
println("="^78)
println("""
Velocity (hence strain-rate) fields superpose linearly across modes for this
potential-flow representation -- confirmed directly: adding contributions
mode-by-mode and evaluating the full 2-mode field agree to machine precision.
""")

let x = 0.4, r = 0.7
    modes = Dict(2 => 0.3, 4 => 0.5)
    e1 = total_strain(modes, x, r)
    e2 = strain_basis(2, x) .* (0.3 * r^0) .+ strain_basis(4, x) .* (0.5 * r^2)
    @assert all(abs.(collect(e1) .- collect(e2)) .< 1e-12)
end
println("ASSERTION 1 OK: total_strain for {mode 2, mode 4} matches direct")
println("term-by-term superposition to machine precision.")

# ------------------------------------------------------------------------------
# Section 2: the raw dissipation integral does NOT reproduce Reid/Lamb's
# l-dependent damping -- confirmed again here (found once already earlier this
# session and not reconciled at the time), and why the RATIO-based Oh_eff_l
# formula in Section 3 is immune to this defect regardless of the details.
# ------------------------------------------------------------------------------

println()
println("="^78)
println("Section 2: raw bulk dissipation of this basis is NOT Lamb's damping")
println("="^78)

function H_kl(k::Int, l::Int; n_x::Int=200)
    x_nodes, x_wts = DropSolver.gauss_legendre_nodes(n_x, -1.0, 1.0)
    s = 0.0
    for (x, w) in zip(x_nodes, x_wts)
        a = strain_basis(k, x)
        b = strain_basis(l, x)
        s += (a[1] * b[1] + a[2] * b[2] + a[3] * b[3] + 2 * a[4] * b[4]) * w
    end
    radial = 1 / (k + l - 1)
    2 * pi * radial * s
end

let
    Hs = [H_kl(l, l) for l in (2, 3, 4, 5)]
    lamb_scaling = [(l - 1) * (2l + 1) for l in (2, 3, 4, 5)]
    ratios = Hs ./ lamb_scaling
    println("H(l,l) for l=2,3,4,5: ", round.(Hs, digits=4))
    println("Lamb's (l-1)(2l+1):    ", lamb_scaling)
    println("ratio H(l,l)/Lamb:     ", round.(ratios, digits=4))
    # H(l,l) DOES grow with l (an earlier, less careful check with a coarse
    # uniform-grid quadrature wrongly suggested it was l-independent -- caught
    # by testing l=3,4,5 explicitly here, not just l=2), but NOT in Lamb's
    # (l-1)(2l+1) proportion: the ratio H(l,l)/Lamb is not constant.
    @assert issorted(Hs)                                  # H grows with l...
    @assert maximum(ratios) / minimum(ratios) > 2.0        # ...but not proportionally to Lamb's scaling
end
println("ASSERTION 2 OK: this basis's raw dissipation integral grows with l,")
println("but not in Reid/Lamb's (l-1)(2l+1) proportion -- confirming (again)")
println("that bulk dissipation of a bare potential-flow field is not a stand-in")
println("for the viscous boundary-layer physics behind the actual damping, and")
println("must only be used as a RATIO (Section 3), where this defect cancels")
println("between numerator and denominator.")

# ------------------------------------------------------------------------------
# Section 3: Oh_eff_l(t), properly multi-mode-coupled.
# ------------------------------------------------------------------------------

println()
println("="^78)
println("Section 3: Oh_eff_l(t) via a dissipation-weighted, multi-mode average")
println("="^78)
println("""
w_l(r,x,t) = dS^2/dAdot_l = 2 * sum_k Adot_k(t) * [e^(k)(r,x) : e^(l)(r,x)]
           = 2 * (e_total(r,x,t) : e^(l)(r,x))   [directly, no need to expand the sum]
Oh_eff_l(t) = Oh_0 * [integral |w_l| * mu_eff(S)/mu_0 dV] / [integral |w_l| dV]
""")

function w_l_field(l::Int, active_modes, x::Float64, r::Float64)
    e_tot = total_strain(active_modes, x, r)
    e_l = strain_basis(l, x) .* r^(l - 2)
    2 * (e_tot[1] * e_l[1] + e_tot[2] * e_l[2] + e_tot[3] * e_l[3] + 2 * e_tot[4] * e_l[4])
end

"""
    oh_eff_coupled(l, active_modes, Oh0, lambda_c, a, eps_ST; n_r=24, n_x=40)

Multi-mode-coupled effective Oh for mode `l`, given the CURRENT active-mode
Adot vector. Numerically integrates over `r in [0,1]`, `x=cos(theta) in
[-1,1]` (Gauss-Legendre in both).
"""
function oh_eff_coupled(l::Int, active_modes, Oh0::Float64, lambda_c::Float64,
    a::Float64, eps_ST::Float64; n_r::Int=24, n_x::Int=40)
    r_nodes, r_wts = DropSolver.gauss_legendre_nodes(n_r, 0.0, 1.0)
    x_nodes, x_wts = DropSolver.gauss_legendre_nodes(n_x, -1.0, 1.0)
    num = 0.0
    den = 0.0
    for (ri, rw) in zip(r_nodes, r_wts), (xi, xw) in zip(x_nodes, x_wts)
        wl = w_l_field(l, active_modes, xi, ri)
        wl == 0.0 && continue
        e_tot = total_strain(active_modes, xi, ri)
        S = shear_rate_S(e_tot...)
        mu_ratio = (1 + (lambda_c * S)^a)^(-eps_ST)
        weight = abs(wl) * rw * xw * ri^2   # r^2 from the volume element
        num += weight * mu_ratio
        den += weight
    end
    den == 0.0 ? Oh0 : Oh0 * num / den
end

# --- Reduction check: single active mode matches the ORIGINAL (validated) K_l result ---
let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, l = 2, Adot_l = 0.05
    Kl = characteristic_shear_K(l)
    gammadot_char = Kl * abs(Adot_l)
    Oh_eff_single_mode_formula = Oh0 * (1 + (lambda_c * gammadot_char)^a)^(-eps_ST)
    Oh_eff_coupled_result = oh_eff_coupled(l, Dict(l => Adot_l), Oh0, lambda_c, a, eps_ST)
    err = abs(Oh_eff_coupled_result - Oh_eff_single_mode_formula) / Oh_eff_single_mode_formula
    println("  single mode l=2: Oh_eff(original K_l formula)=$(round(Oh_eff_single_mode_formula,digits=4))" *
            "  Oh_eff(coupled, 1 active mode)=$(round(Oh_eff_coupled_result,digits=4))  rel_err=$(round(err,digits=4))")
    @assert err < 0.01
end
println("ASSERTION 3 OK: with only one mode active, the coupled formula")
println("reproduces the original (already-validated) single-mode Oh_eff to <1%.")

# --- Newtonian limit: Oh_eff_l = Oh0 identically, regardless of coupling ---
let Oh0 = 3.0, l = 3
    for active_modes in (Dict(2 => 0.3), Dict(2 => 0.3, 4 => 0.6, 5 => -0.2))
        Oh_eff = oh_eff_coupled(l, active_modes, Oh0, 0.0, 2.0, 0.5)   # lambda_c=0 -> no thinning
        @assert isapprox(Oh_eff, Oh0; rtol=1e-10)
    end
end
println("ASSERTION 4 OK: Oh_eff_l = Oh_0 exactly (to 1e-10) in the Newtonian")
println("limit (lambda_c=0), for both a single active mode and several")
println("simultaneously active modes -- confirming the ratio construction is")
println("immune to Section 2's raw-dissipation defect.")

# --- The actual point: cross-mode coupling ENHANCES thinning beyond what
#     mode 2's own shear rate alone would suggest, once other modes are
#     excited (e.g. by contact) at comparable or even smaller amplitude. ---
println()
let Oh0 = 57.4, lambda_c = 30507.0, a = 0.7431, eps_ST = 0.99956, l = 2
    Adot2 = 0.02   # a modest mode-2 velocity, e.g. early in a low-We impact
    Oh_eff_self_only = oh_eff_coupled(l, Dict(2 => Adot2), Oh0, lambda_c, a, eps_ST)
    for Adot_high in (0.0, 0.02, 0.05, 0.1)
        active = Dict(2 => Adot2, 5 => Adot_high)
        Oh_eff = oh_eff_coupled(l, active, Oh0, lambda_c, a, eps_ST)
        println("  Adot_2=$Adot2, Adot_5=$Adot_high: Oh_eff_2(coupled)=$(round(Oh_eff,digits=4))" *
                "  (self-only would give $(round(Oh_eff_self_only,digits=4)))")
    end
end
println("""
As mode 5's amplitude grows (contact exciting a higher mode alongside the
dominant mode 2), mode 2's OWN effective Oh drops further than its
self-only shear rate would predict -- exactly the mechanism the self-only
model was missing: shear is kinematic, so mode 2 "feels" the extra thinning
that mode 5's velocity field creates at the same points in the drop.
""")

println()
println("="^78)
println("Summary")
println("="^78)
println("""
Fix, precisely: mode l's effective Oh is no longer a function of Adot_l
alone. It is a dissipation-weighted average of the TRUE local viscosity
ratio mu_eff(S_total)/mu_0, where S_total is the shear-rate invariant of the
FULL superposed strain field from every currently active mode, weighted by
mode l's own sensitivity (dS^2/dAdot_l) to that same field. This is exact
kinematics (linear superposition of potential-flow velocity fields), uses
only the exact (non-Taylor-expanded) Carreau-Yasuda law, reduces to the
already-validated single-mode result when only one mode is active, and to
Oh_eff_l=Oh_0 identically in the Newtonian limit regardless of the known
defect in this basis's raw (un-ratioed) dissipation normalization.

Not yet done: wiring this into julia/src/st_exact_extension.jl (replacing
the self-only oh_eff_lambda_omega2 with a version that takes the full
current mode-velocity vector, and reworking the residual/Jacobian since
Oh_eff_l for EVERY mode now depends on EVERY mode's Adot -- the Jacobian's
D2 block is no longer diagonal); revalidating against the 20 sampled
experiments.
""")
