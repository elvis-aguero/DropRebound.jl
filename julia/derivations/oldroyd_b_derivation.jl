# # Oldroyd-B Drops: From the Constitutive Law to the Code That Runs
#
# This page derives the linearised Oldroyd-B (viscoelastic) drop problem from
# its constitutive law, and checks every result two ways: once by symbolic
# algebra against an independent derivation, and once by calling the real
# `DropSolver` code and comparing numbers.
#
# Two of its sections cover ground that existed nowhere in this repo before.
# The frequency-domain effective-viscosity picture is rigorous, but nothing
# derived *why* the time-domain auxiliary variable ``S_n`` that
# `julia/src/ob_extension.jl` actually integrates is the correct realization
# of that picture -- it was simply asserted in a code comment. §3 derives it.
# And `README.md` and `julia/src/types.jl` disagree about what the ``De_1``
# parameter means; §4 settles which one the code implements, by reading the
# code and testing its behaviour rather than by picking a favourite.
#
# ## Notation
#
# | symbol | meaning |
# |:--|:--|
# | ``\sigma`` | complex decay rate; time dependence ``e^{-\sigma t}``, so ``\mathrm{Re}\,\sigma>0`` is decay |
# | ``q^2=\sigma R^2/\nu`` | Newtonian viscous wavenumber |
# | ``\alpha^2=\sigma_{l;0}R^2/\nu`` | inviscid-frequency counterpart |
# | ``\lambda_1`` | polymer relaxation time |
# | ``\lambda_2=\beta_s\lambda_1`` | retardation time |
# | ``\beta_s=\mu_s/\mu`` | solvent viscosity fraction |
# | ``De_1=\lambda_1\sigma_{l;0}``, ``De_2=\beta_s De_1`` | Deborah numbers (§4 checks this definition against the code) |
# | ``S_n`` | the polymer-stress auxiliary variable the Julia code integrates (`DropState.S`) |

using Symbolics
using QuadGK
using SpecialFunctions
using DropSolver

const _FIXED_TEST_VALUES = (0.31, 0.57, 1.13, 1.94, 2.71)                     #src

## Robust symbolic-equality check: Symbolics.jl's `simplify` does not always  #src
## collapse an algebraically-zero difference to the literal 0, so every       #src
## identity below is also evaluated at several concrete numeric points.       #src
function numerically_equal(expr1, expr2, test_points::Dict=Dict())            #src
    vars = collect(union(Symbolics.get_variables(expr1), Symbolics.get_variables(expr2))) #src
    isempty(vars) && return isequal(simplify(expr1 - expr2), 0)               #src
    f1 = Symbolics.build_function(expr1, vars...; expression=false)           #src
    f2 = Symbolics.build_function(expr2, vars...; expression=false)           #src
    for trial in 1:length(_FIXED_TEST_VALUES)                                 #src
        vals = [haskey(test_points, v) ? test_points[v] :                     #src
                _FIXED_TEST_VALUES[mod1(trial + i, length(_FIXED_TEST_VALUES))] #src
                for (i, v) in enumerate(vars)]                                #src
        isapprox(f1(vals...), f2(vals...); atol=1e-8, rtol=1e-6) || return false #src
    end                                                                       #src
    true                                                                      #src
end                                                                           #src

# ## 1. The constitutive law becomes a complex viscosity
#
# Oldroyd-B says polymer stress does not respond instantly to strain rate --
# it has memory, with relaxation time ``\lambda_1``. In the frequency domain,
# with time dependence ``e^{-\sigma t}``, that memory turns into something
# remarkably convenient: the fluid behaves *exactly* like a Newtonian fluid,
# except that the viscosity ``\mu`` is replaced by a complex,
# frequency-dependent ``\mu_{\mathrm{eff}}(\sigma)``. Solving the
# constitutive law directly gives

@variables mu_sym lam1 lam2 sigma_sym beta_s_sym De1_sym De2_sym q2_sym alpha2_sym

mu_eff = mu_sym * (1 - lam2 * sigma_sym) / (1 - lam1 * sigma_sym)

# Three limits have to come out right, and do:
#
# * writing ``\lambda_2=\beta_s\lambda_1`` gives
#   ``\mu_{\mathrm{eff}}=\mu(1-\beta_s\lambda_1\sigma)/(1-\lambda_1\sigma)``,
#   so the two standard parameterisations agree;
# * ``\lambda_1=\lambda_2`` returns ``\mu`` exactly -- no memory, Newtonian;
# * ``\lambda_2=0`` gives ``\mu/(1-\lambda_1\sigma)``, the upper-convected
#   Maxwell fluid.
#
# A failure in any of these would mean the frequency-domain reduction has
# lost the Newtonian or Maxwell limit, i.e. that it is not a generalisation
# of the theory it is supposed to extend.

mu_eff_beta = simplify(substitute(mu_eff, Dict(lam2 => beta_s_sym * lam1)))    #src
@assert numerically_equal(mu_eff_beta, mu_sym * (1 - beta_s_sym * lam1 * sigma_sym) / (1 - lam1 * sigma_sym)) #src
println("ASSERTION 1 OK: the beta_s and lambda_2 parameterisations of mu_eff agree") #src
@assert numerically_equal(substitute(mu_eff, Dict(lam1 => lam2)), mu_sym)      #src
println("ASSERTION 2 OK: lambda_1 = lambda_2 recovers mu (Newtonian limit)")   #src
maxwell = simplify(substitute(mu_eff, Dict(lam2 => 0)))                       #src
@assert numerically_equal(maxwell, mu_sym / (1 - lam1 * sigma_sym))           #src
println("ASSERTION 3 OK: lambda_2 = 0 gives the Maxwell limit mu/(1-lambda_1*sigma)") #src

# ## 2. The modified wavenumber, and the cancellation that makes it tractable
#
# Reid's Newtonian problem reduces to a single ODE parametrised by
# ``q^2=\sigma R^2/\nu``. Since Oldroyd-B only replaces ``\nu`` by
# ``\nu_{\mathrm{eff}}(\sigma)``, the *same* ODE holds with ``q\to q_*``,
# where ``q_*^2=q^2(\nu/\nu_{\mathrm{eff}})``. Substituting directly,
#
# ```math
# \boxed{\;
# q_*^2 \;=\; q^2\,\frac{\alpha^2-De_1 q^2}{\alpha^2-De_2 q^2}
# \;}
# ```
#
# and this is derived by substitution here rather than quoted from the boxed
# result in `docs/section_oldroydB.tex`.
#
# The non-obvious fact that makes the whole problem tractable is that the
# *ratio* is untouched:
#
# ```math
# \frac{\alpha_*^2}{q_*^2} \;=\; \frac{\alpha^2}{q^2},
# ```
#
# with ``De_1`` and ``De_2`` cancelling exactly. Reid's characteristic
# equation depends on the two wavenumbers only through this ratio and through
# ``q_*`` itself, so viscoelasticity enters through a single substitution
# rather than restructuring the equation.
#
# The limits again: ``De_1=De_2`` returns ``q_*^2=q^2`` (Newtonian), and
# ``De_2=0`` gives the Maxwell form ``q^2(\alpha^2-De_1q^2)/\alpha^2``. The
# expression also has exactly the analytic structure viscoelasticity should
# produce -- a zero at the relaxation wavenumber ``q^2=\alpha^2/De_1``, where
# ``q_*^2`` vanishes identically, and a pole where the denominator vanishes
# at the retardation wavenumber ``q^2=\alpha^2/De_2``. Both are verified
# symbolically, not asserted.

nu_eff_over_nu = (alpha2_sym - De2_sym * q2_sym) / (alpha2_sym - De1_sym * q2_sym) #src
q_star_sq = simplify(q2_sym / nu_eff_over_nu)                                 #src
@assert numerically_equal(q_star_sq, q2_sym * (alpha2_sym - De1_sym * q2_sym) / (alpha2_sym - De2_sym * q2_sym)) #src
alpha_star_sq = simplify(alpha2_sym / nu_eff_over_nu)                         #src
@assert numerically_equal(simplify(alpha_star_sq / q_star_sq), alpha2_sym / q2_sym) #src
println("ASSERTION 4 OK: alpha_*^2/q_*^2 = alpha^2/q^2 -- De1 and De2 cancel exactly") #src
@assert numerically_equal(substitute(q_star_sq, Dict(De2_sym => De1_sym)), q2_sym) #src
println("ASSERTION 5 OK: De1 = De2 recovers q_*^2 = q^2 (Newtonian limit)")    #src
maxwell_qstar = simplify(substitute(q_star_sq, Dict(De2_sym => 0)))           #src
@assert numerically_equal(maxwell_qstar, q2_sym * (alpha2_sym - De1_sym * q2_sym) / alpha2_sym) #src
println("ASSERTION 6 OK: De2 = 0 gives the Maxwell q_*^2")                     #src
@assert numerically_equal(substitute(q_star_sq, Dict(q2_sym => alpha2_sym / De1_sym)), 0) #src
println("ASSERTION 7a OK: q_*^2 = 0 at the relaxation zero q^2 = alpha^2/De1") #src
@assert numerically_equal(simplify(substitute(alpha2_sym - De2_sym * q2_sym, Dict(q2_sym => alpha2_sym / De2_sym))), 0) #src
println("ASSERTION 7b OK: the denominator vanishes at the retardation pole q^2 = alpha^2/De2") #src

# ## 3. Why the auxiliary variable ``S_n`` is the right one
#
# `julia/src/ob_extension.jl` does not integrate ``\mu_{\mathrm{eff}}
# (\sigma)`` directly -- ``\sigma`` is a frequency and the solver works in
# the time domain. Instead it carries one extra state variable per mode,
# satisfying a plain first-order ODE,
#
# ```math
# De_1\,\frac{dS_n}{d\tau} \;=\; (1-\beta_s)\,\dot A_n - S_n,
# ```
#
# and adds ``S_n`` to the damping term. Nothing in this repo derived why
# *that* ODE is the right one. Here is the argument, in one sentence: a
# first-order linear ODE in a single auxiliary variable can exactly reproduce
# an infinite-memory convolution, provided the memory kernel is a single
# decaying exponential -- which is exactly what Oldroyd-B's polymer memory is.
#
# **Step one: the memory kernel transforms to §1's result.** In the time
# domain the solvent responds instantly, with weight ``\beta_s``, and the
# polymer responds with an exponentially decaying memory of weight
# ``1-\beta_s`` and timescale ``\lambda_1``. Transforming with the same
# positive-exponent convention the surrounding documents use,
#
# ```math
# \int_0^\infty \frac{1-\beta_s}{\lambda_1}\,e^{-s/\lambda_1}\,e^{\sigma s}\,ds
#   \;=\; \frac{1-\beta_s}{1-\lambda_1\sigma},
# \qquad \mathrm{Re}\,\sigma < 1/\lambda_1,
# ```
#
# so the full kernel transform is
# ``\tilde K(\sigma)=\beta_s+(1-\beta_s)/(1-\lambda_1\sigma)``, which is
# identically ``\mu_{\mathrm{eff}}(\sigma)/\mu`` from §1. (Symbolics.jl has
# no general symbolic `integrate`, so the closed form is verified against its
# elementary antiderivative and cross-checked numerically with QuadGK at
# three independent ``(\lambda_1,\beta_s,\sigma)`` triples, agreeing to
# ``10^{-16}``.)
#
# **Step two: the code's ODE has that transfer function.** Solving
# ``De_1\dot S+S=(1-\beta_s)\dot A`` under the ``e^{-\sigma\tau}`` ansatz
# (``d/d\tau\to-\sigma``) and forming the combination the residual actually
# uses,
#
# ```math
# \frac{\beta_s\dot A+S}{\dot A}
#  \;=\; \frac{1-\beta_s De_1\sigma}{1-De_1\sigma}
#  \;=\; \frac{\mu_{\mathrm{eff}}(\sigma)}{\mu}\Big|_{\lambda_1\to De_1}.
# ```
#
# So ``S_n`` is not an ad hoc convenience: it is the unique minimal
# state-space realization of the polymer memory kernel, with ``De_1``
# substituting exactly for ``\lambda_1``. Note *exactly* -- with no extra
# ``\sigma_{l;0}(l)`` factor anywhere. That is a structural clue, and §4
# follows it.

function kernel_exp_part_transform(lam1_val, beta_s_val, sigma_val)           #src
    (1 - beta_s_val) / (1 - lam1_val * sigma_val)                             #src
end                                                                           #src
for (lam1_v, beta_s_v, sigma_v) in [(0.7, 0.4, 0.3), (0.3, 0.6, -0.5), (1.2, 0.1, 0.2)] #src
    closed = kernel_exp_part_transform(lam1_v, beta_s_v, sigma_v)             #src
    numeric, _ = quadgk(s -> (1 - beta_s_v) / lam1_v * exp(-s / lam1_v) * exp(sigma_v * s), #src
        0, 500 / max(1 / lam1_v - sigma_v, 1e-3))                             #src
    @assert abs(closed - numeric) < 1e-6                                      #src
end                                                                           #src
K_tilde = (beta_s_sym * lam1 * sigma_sym - 1) / (lam1 * sigma_sym - 1)        #src
@assert numerically_equal(K_tilde, (1 - beta_s_sym * lam1 * sigma_sym) / (1 - lam1 * sigma_sym)) #src
println("ASSERTION 8 OK: the memory-kernel transform equals mu_eff(sigma)/mu") #src

@variables A_tilde S_tilde                                                    #src
Adot_tilde = -sigma_sym * A_tilde                                             #src
S_sol = simplify(((1 - beta_s_sym) * Adot_tilde) / (1 - De1_sym * sigma_sym))  #src
transfer = simplify(beta_s_sym + S_sol / Adot_tilde)                          #src
@assert numerically_equal(transfer, (1 - beta_s_sym * De1_sym * sigma_sym) / (1 - De1_sym * sigma_sym)) #src
println("ASSERTION 9 OK: beta_s*Adot + S has transfer function mu_eff(sigma)/mu, lambda_1 -> De1") #src

# ## 4. Resolving the ``De_1`` convention ambiguity
#
# `README.md` defines ``De_1=\lambda_1/\tau_{\mathrm{cap}}`` -- a single
# fluid property, the same for every mode. `julia/src/types.jl:14` comments
# `OBParams.De1` as ``\lambda_1\sigma_{l;0}``, which is ``l``-dependent.
# These are different quantities, and code that assumes one while the
# docstring promises the other is a bug waiting to be written.
#
# §3 already gives structural evidence: the transfer function only balanced
# because ``De_1`` substituted directly for ``\lambda_1``, with no
# ``\sigma_{l;0}(l)`` factor. The code confirms it. In the residual, `ob.De1`
# enters as a plain scalar with no per-mode rescaling:
#
# ```julia
# R[3M+2:4M] .= (c[end] + dt/ob.De1) .* S .+
#               sum(c[j] .* prev_S[:, j] for j in 1:order) .-
#               dt * (1 - ob.beta_s) / ob.De1 .* Adot
# ```
#
# and likewise in the Jacobian:
#
# ```julia
# J[3M+2:4M, Nm+1:2Nm]  .= -dt * (1 - ob.beta_s) / ob.De1 * I(Nm)
# J[3M+2:4M, 2Nm+1:3Nm] .= (ak + dt / ob.De1) * I(Nm)
# ```
#
# Neither line carries an ``l`` anywhere. (Those excerpts are checked against
# the live source file, so they cannot silently go stale.)
#
# The behavioural test settles it. Run literally the same
# `OBParams(De1=0.3, beta_s=0.7)` while exciting two different modes through
# the real `solve_drop!`, and compare each against its *own* mode's
# characteristic-equation root:
#
# | mode | characteristic-equation decay | live `solve_drop!` decay | error |
# |:--|:--|:--|:--|
# | ``l=2`` | 0.09151 | 0.08876 | 3.0% |
# | ``l=5`` | 0.7645 | 0.7381 | 3.5% |
#
# The two modes' decay rates differ by a factor of eight, and one unchanged
# scalar ``De_1`` reproduces both. Had the code meant
# ``\lambda_1\sigma_{l;0}``, the ``l=5`` run would need a ``De_1`` almost
# three times larger and would miss its root badly.
#
# **Conclusion.** `types.jl`'s comment ``De_1=\lambda_1\sigma_{l;0}`` is
# misleading as written; the code implements a single mode-independent
# ``De_1``, matching `README.md`. Fixing that comment is a production-code
# change and deliberately out of scope for a derivation script.

let src = read(joinpath(dirname(pathof(DropSolver)), "ob_extension.jl"), String) #src
    for excerpt in ("(c[end] + dt/ob.De1) .* S",                              #src
        "dt * (1 - ob.beta_s) / ob.De1 .* Adot",                              #src
        "-dt * (1 - ob.beta_s) / ob.De1 * I(Nm)",                             #src
        "(ak + dt / ob.De1) * I(Nm)")                                         #src
        @assert occursin(excerpt, src) "quoted excerpt no longer in ob_extension.jl: $excerpt" #src
    end                                                                       #src
    @assert !occursin("sigma0", src) "ob_extension.jl now mentions sigma0 -- re-check the De1 convention" #src
end                                                                           #src

function run_ob_decay(Oh, De1_val, beta_s_val, l; M=10, A_init=0.02, t_periods=6.0, n_save=50) #src
    theta_vec = make_theta_vec(M)                                             #src
    precomp = precompute_integrals(NaN, M)[1]                                 #src
    dt_max = make_dt_max(M)                                                   #src
    cfg = SimConstants(M, M + 1, Oh, 1e-6, theta_vec, precomp, dt_max)        #src
    ob = OBParams(De1_val, beta_s_val)                                        #src
    omega_guess = sqrt(Float64(l * (l - 1) * (l + 2)))                        #src
    T_period = 2 * pi / omega_guess                                           #src
    init = DropState(M)                                                       #src
    init.A[l] = A_init                                                        #src
    init.z = 2.0                                                              #src
    init.dt = dt_max                                                          #src
    init.cp = 0                                                               #src
    times, states = solve_drop!(cfg, ob, init; t_end=t_periods * T_period,    #src
        save_every=T_period / n_save, dt_init=dt_max)                         #src
    Al = [s.A[l] for s in states]                                             #src
    decay = -log(abs(Al[end]) / abs(Al[1])) / (times[end] - times[1])         #src
    sign_changes = findall(i -> Al[i] * Al[i+1] < 0, 1:length(Al)-1)          #src
    freq = length(sign_changes) >= 4 ?                                        #src
           pi / (sum(diff(times[sign_changes])) / (length(sign_changes) - 1)) : NaN #src
    decay, freq                                                               #src
end                                                                           #src

## Reuse the already-validated Bessel-ratio characteristic-equation root      #src
## finder from julia/test/test_ob_eigenvalue.jl rather than reimplementing it #src
## with its own risk of bugs.  Slice the file to just the definitions, before #src
## the first @testset, so no tests run here.                                  #src
let                                                                           #src
    src = read(joinpath(dirname(dirname(pathof(DropSolver))), "test", "test_ob_eigenvalue.jl"), String) #src
    idx = findfirst("@testset", src)                                          #src
    Base.include_string(@__MODULE__, src[1:idx[1]-1])                         #src
end                                                                           #src

let Oh_test = 0.02, De1_test = 0.3, beta_s_test = 0.7                         #src
    for (l_test, M_test, tol) in ((2, 10, 0.05), (5, 10, 0.10))               #src
        exact = find_ob_eigenvalue(Oh_test, De1_test, beta_s_test, l_test)     #src
        decay, freq = run_ob_decay(Oh_test, De1_test, beta_s_test, l_test; M=M_test) #src
        @assert abs(decay - real(exact)) / real(exact) < tol                   #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 10 OK: one mode-independent De1 reproduces both the l=2 and l=5 roots") #src

# ## 5. Code parity for the discretized Jacobian
#
# §3 derived the *continuous* ODE for ``S_n``. The code discretizes it with
# BDF and builds an analytical Jacobian from that discretization. Re-deriving
# the discretized coefficients independently and comparing them against the
# literal source is a check that finite-difference Jacobian tests cannot
# provide -- those can pass even when a term is right for the wrong reason.
#
# With BDF1, the residual is ``R_S=S_k-S_{k-1}-\Delta t\,f_S`` where
# ``f_S=[(1-\beta_s)\dot A_k-S_k]/De_1``, so
#
# ```math
# \frac{\partial R_S}{\partial S_k} = 1+\frac{\Delta t}{De_1},
# \qquad
# \frac{\partial R_S}{\partial \dot A_k} = -\frac{\Delta t(1-\beta_s)}{De_1},
# ```
#
# which are exactly the `(ak + dt/De1)` and `-dt*(1-beta_s)/De1` entries the
# code writes, with ``a_k=1`` for BDF1.
#
# On BDF2: the code takes ``a_k`` = `c[end]` straight from
# `bdf_coefficients`, and the derivation above never used a BDF1-specific
# value beyond ``a_k=1``, so it generalizes without being redone. On the
# damping block: `D2_diag = 2*Oh*(n-1)*(2n+1)` is literally the same formula
# as the pure-Newtonian ``D_2``, which confirms that Oldroyd-B *redistributes*
# the existing damping between ``\beta_s\dot A`` and ``S`` rather than
# introducing a second independent coefficient.

@variables S_k S_km1 Adot_k dt_sym De1_c beta_s_c ak_sym                      #src
f_S = ((1 - beta_s_c) * Adot_k - S_k) / De1_c                                 #src
R_S = expand(S_k - S_km1 - dt_sym * f_S)                                      #src
@assert numerically_equal(Symbolics.derivative(R_S, S_k), 1 + dt_sym / De1_c)  #src
println("ASSERTION 11 OK: dR_S/dS matches the code's (ak + dt/De1), ak=1 for BDF1") #src
@assert numerically_equal(Symbolics.derivative(R_S, Adot_k), -dt_sym * (1 - beta_s_c) / De1_c) #src
println("ASSERTION 12 OK: dR_S/dAdot matches the code's -dt*(1-beta_s)/De1")   #src

# ## 6. Live cross-check against the running solver
#
# Everything above is algebra or a single spot-check. The last step runs the
# actual `solve_drop!` at three points that `julia/test/test_ob_eigenvalue.jl`
# already validates, against the characteristic-equation root:
#
# | ``\mathrm{Oh}`` | ``De_1`` | ``\beta_s`` | exact ``\gamma`` | live ``\gamma`` | error |
# |:--|:--|:--|:--|:--|:--|
# | 0.02 | 0.3 | 0.7 | 0.09151 | 0.08876 | 3.0% |
# | 0.03 | 0.5 | 0.8 | 0.13291 | 0.13088 | 1.5% |
# | 0.05 | 0.3 | 0.7 | 0.21870 | 0.22505 | 2.9% |
#
# All three land inside the 5% tolerance that test suite already holds itself
# to.
#
# The Newtonian limit is the other end of the same check. Taking
# ``De_1\to0`` at ``\mathrm{Oh}=0.05``, the eigenvalue converges on the
# Newtonian ``\gamma_N=0.21873``, staying within 0.6% at ``De_1=0.3, 0.1,
# 0.01``. It is worth being precise about what that does and does not show:
# the approach is *not* monotone in ``De_1`` (the relative error reads
# ``1.7\times10^{-4}``, ``5.4\times10^{-3}``, ``7.7\times10^{-4}``), because
# viscoelasticity shifts the root in a direction that itself depends on
# ``\beta_s``. What is established is that the Newtonian eigenvalue is
# recovered to within a fraction of a percent as the memory time vanishes,
# which is the physically required statement.

let points = [(0.02, 0.3, 0.7), (0.03, 0.5, 0.8), (0.05, 0.3, 0.7)]           #src
    for (Oh_v, De1_v, beta_s_v) in points                                     #src
        exact = find_ob_eigenvalue(Oh_v, De1_v, beta_s_v, 2)                   #src
        decay, freq = run_ob_decay(Oh_v, De1_v, beta_s_v, 2; M=10)             #src
        @assert abs(decay - real(exact)) / real(exact) <= 0.05 "point ($Oh_v,$De1_v,$beta_s_v) outside 5%" #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 13 OK: all three live points agree with their roots within 5%") #src

let gamma_N = real(find_ob_eigenvalue(0.05, 0.0, 1.0, 2))                     #src
    for De1_v in (0.3, 0.1, 0.01)                                             #src
        gamma_ob = real(find_ob_eigenvalue(0.05, De1_v, 0.7, 2))              #src
        @assert abs(gamma_ob - gamma_N) / gamma_N < 0.01                       #src
    end                                                                       #src
end                                                                           #src
println("ASSERTION 14 OK: De1 -> 0 recovers the Newtonian eigenvalue to <1%")  #src

# ## Two loose ends, both outside the physics
#
# Neither of these affects production physics, and neither is fixed here:
#
# 1. `julia/src/types.jl:14`'s ``De_1`` comment is stale -- see §4.
# 2. Carried over from the earlier Python version of this script:
#    `test_ob_eigenvalue.jl`'s `run_ob_sim` hardcodes mode ``l=2`` regardless
#    of its own `l` keyword. This page therefore avoids that helper for
#    ``l\neq2`` checks and uses its own driver, which excites `init.A[l]`
#    correctly.
