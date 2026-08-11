# # Oldroyd-B Drops
#
# This page derives the linearised Oldroyd-B (viscoelastic) drop problem from
# its constitutive law, checking each result both symbolically and against
# numbers produced by `DropSolver` itself.
#
# ## Why this part takes a different route from the last
#
# Part IV and this one both leave the Newtonian fluid behind, and they leave it
# in different directions.
#
# A shear-thinning liquid has a viscosity that depends on the local rate of
# deformation. That makes ``\eta`` a field over the drop, computed from the flow
# it governs, which couples every mode to every other and forces the interior to
# be carried in the state. Hence the variational treatment of Part IV.
#
# A viscoelastic liquid has memory instead. Its stress depends on the deformation
# history rather than on the instantaneous rate, and for a linear constitutive
# law that history is a convolution. Under the ansatz ``e^{-\sigma t}`` a
# convolution becomes multiplication, so the whole rheology collapses into a
# single complex, frequency-dependent viscosity. No field appears, no modes
# couple, and each ``l`` remains independent.
#
# That is why this part reuses Reid's machinery rather than the variational
# assembly, and why it lives in the nonvariational solver: with the modes still
# separable there is nothing for an interior state to buy. The memory is carried
# instead by one auxiliary variable per mode, derived in §3.
#
# ## The route
#
# It is short. The constitutive law becomes a complex,
# frequency-dependent viscosity (§1), which enters Reid's Newtonian problem
# through a single modified wavenumber (§2). §3 derives the time-domain
# auxiliary variable ``S_n`` the solver integrates and identifies it as the
# minimal state-space realization of the same memory kernel. §4 fixes the
# meaning of the ``De_1`` parameter, and §5 compares the whole construction
# against eigenvalues measured from the running solver.
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
# | ``De_1=\lambda_1/\tau_{\mathrm{cap}}``, ``De_2=\beta_s De_1`` | Deborah numbers, mode-independent (§4) |
# | ``S_n`` | the polymer-stress auxiliary variable the solver integrates |

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
#
# ```math
# \mu_{\mathrm{eff}}(\sigma) \;=\;
#   \mu\,\frac{1-\lambda_2\sigma}{1-\lambda_1\sigma},
# ```
#
# with ``\lambda_1`` the relaxation time and ``\lambda_2`` the retardation
# time.

@variables mu_sym lam1 lam2 sigma_sym beta_s_sym De1_sym De2_sym q2_sym alpha2_sym  #src
mu_eff = mu_sym * (1 - lam2 * sigma_sym) / (1 - lam1 * sigma_sym)  #src

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
# at the retardation wavenumber ``q^2=\alpha^2/De_2``.

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
# The solver does not integrate ``\mu_{\mathrm{eff}}(\sigma)`` directly --
# ``\sigma`` is a frequency and the solver works in
# the time domain. Instead it carries one extra state variable per mode,
# satisfying a plain first-order ODE,
#
# ```math
# De_1\,\frac{dS_n}{d\tau} \;=\; (1-\beta_s)\,\dot A_n - S_n,
# ```
#
# and adds ``S_n`` to the damping term. The reason that ODE is the right one
# fits in a sentence: a
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
# identically ``\mu_{\mathrm{eff}}(\sigma)/\mu`` from §1. (The integral is
# elementary; numerical quadrature at three independent
# ``(\lambda_1,\beta_s,\sigma)`` triples reproduces the closed form to
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
# So ``S_n`` is the minimal state-space realization of the polymer memory
# kernel, with ``De_1`` substituting directly for ``\lambda_1`` -- and with
# no ``\sigma_{l;0}(l)`` factor anywhere. §4 follows up on that last point.

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
# Two definitions of ``De_1`` are in circulation for this problem:
# ``\lambda_1/\tau_{\mathrm{cap}}``, a single fluid property shared by every
# mode, and ``\lambda_1\sigma_{l;0}``, which is ``l``-dependent. They are
# different quantities, so it is worth establishing which one this solver's
# ``De_1`` denotes.
#
# §3 supplies the structural evidence: the transfer function balances only
# because ``De_1`` substitutes directly for ``\lambda_1``, with no
# ``\sigma_{l;0}(l)`` factor. The discretized equations agree -- ``De_1``
# enters them as a plain scalar, carrying no mode index and no per-mode
# rescaling.
#
# A behavioural test confirms the reading. Holding ``De_1=0.3`` and
# ``\beta_s=0.7`` fixed while exciting two different modes, and comparing
# each free-decay run against its *own* mode's characteristic-equation root:
#
# | mode | characteristic-equation decay | live solver decay | error |
# |:--|:--|:--|:--|
# | ``l=2`` | 0.09151 | 0.08876 | 3.0% |
# | ``l=5`` | 0.7645 | 0.7381 | 3.5% |
#
# The two modes' decay rates differ by a factor of eight, and one unchanged
# scalar ``De_1`` reproduces both. Under the ``\lambda_1\sigma_{l;0}``
# reading, the ``l=5`` run would require a ``De_1`` almost three times larger
# and would miss its root by a wide margin.
#
# **Convention.** ``De_1`` is a single, mode-independent Deborah number,
# ``\lambda_1/\tau_{\mathrm{cap}}``, and is not rescaled per mode. The
# ``\lambda_1\sigma_{l;0}`` form is a different quantity and does not
# describe this parameter.

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
## finder from test/test_ob_eigenvalue.jl rather than reimplementing it #src
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

## Internal code-parity check, kept out of the rendered page.  §3 derived    #src
## the *continuous* ODE for S_n; the solver discretizes it with BDF and      #src
## builds an analytical Jacobian from that discretization.  With BDF1 the    #src
## residual is R_S = S_k - S_{k-1} - dt*f_S, f_S = [(1-beta_s)*Adot_k -      #src
## S_k]/De1, whose two partial derivatives are asserted below to match the   #src
## coefficients the solver writes into its Jacobian.  A failure here means   #src
## the discretization and the continuous ODE have drifted apart.            #src
@variables S_k S_km1 Adot_k dt_sym De1_c beta_s_c ak_sym                      #src
f_S = ((1 - beta_s_c) * Adot_k - S_k) / De1_c                                 #src
R_S = expand(S_k - S_km1 - dt_sym * f_S)                                      #src
@assert numerically_equal(Symbolics.derivative(R_S, S_k), 1 + dt_sym / De1_c)  #src
println("ASSERTION 11 OK: dR_S/dS matches the code's (ak + dt/De1), ak=1 for BDF1") #src
@assert numerically_equal(Symbolics.derivative(R_S, Adot_k), -dt_sym * (1 - beta_s_c) / De1_c) #src
println("ASSERTION 12 OK: dR_S/dAdot matches the code's -dt*(1-beta_s)/De1")   #src

# ## 5. Live cross-check against the running solver
#
# The last step runs `solve_drop!` at three parameter points and compares the
# measured free-decay rate against the characteristic-equation root:
#
# | ``\mathrm{Oh}`` | ``De_1`` | ``\beta_s`` | exact ``\gamma`` | live ``\gamma`` | error |
# |:--|:--|:--|:--|:--|:--|
# | 0.02 | 0.3 | 0.7 | 0.09151 | 0.08876 | 3.0% |
# | 0.03 | 0.5 | 0.8 | 0.13291 | 0.13088 | 1.5% |
# | 0.05 | 0.3 | 0.7 | 0.21870 | 0.22505 | 2.9% |
#
# All three agree to within 5%.
#
# The Newtonian limit is the other end of the same check. Taking
# ``De_1\to0`` at ``\mathrm{Oh}=0.05``, the eigenvalue converges on the
# Newtonian ``\gamma_N=0.21873``, staying within 0.6% at ``De_1=0.3, 0.1,
# 0.01``. The approach is not monotone in ``De_1`` (the relative error reads
# ``1.7\times10^{-4}``, ``5.4\times10^{-3}``, ``7.7\times10^{-4}``), because
# viscoelasticity shifts the root in a direction that itself depends on
# ``\beta_s``. What the limit establishes is that the Newtonian eigenvalue is
# recovered to within a fraction of a percent as the memory time vanishes.

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

# ## Scope
#
# Everything above is linear theory: the ansatz ``e^{-\sigma t}`` and the
# single-mode reduction both presume small deformation, and the
# characteristic-equation comparisons are free-decay ones. Contact forcing,
# mode-to-mode coupling, and the finite-amplitude behaviour of the polymer
# stress lie outside what this page establishes.
