using Test
using DropSolver

# ---------------------------------------------------------------------------
# Spherical Bessel function j_l(z) via upward recurrence (works for complex z)
# ---------------------------------------------------------------------------
function sph_bessel_j(l::Int, z::Number)
    if abs(z) < 1e-300
        return l == 0 ? complex(1.0) : complex(0.0)
    end
    # Small |z|: upward recurrence suffers catastrophic cancellation for l≥2
    # because the two large terms (each O(1/z^l)) nearly cancel.
    # Use leading-order power series j_l(z) = z^l / (2l+1)!!.
    if abs(z) < 0.1
        fac = one(float(real(z)))
        for k in 1:2:2l+1
            fac *= k
        end
        return z^l / fac
    end
    j0 = sin(z) / z
    l == 0 && return j0
    j1 = sin(z) / z^2 - cos(z) / z
    l == 1 && return j1
    jp, jc = j0, j1
    for n in 1:l-1
        jn = (2n + 1) / z * jc - jp
        jp, jc = jc, jn
    end
    return jc
end

# ---------------------------------------------------------------------------
# Stably compute j_{l+1}(z) / j_l(z) (the Bessel ratio Q used in the char eq).
#
# Two failure modes in the naive j_{l+1}/j_l computation:
#   1. Large Im(z) (> ~709): sin(z) overflows Float64 → NaN.
#   2. Small |z|: sph_bessel_j loses precision for l≥2 → wrong ratio.
#
# Fix: upward ratio recurrence r_{k+1} = (2k+1)/z - 1/r_k, starting from
#   r_1 = j_1/j_0 = 1/z - cot(z).
# After l steps the recurrence yields r_{l+1} = j_{l+1}/j_l.
# Guards handle the two pathological limits analytically.
# ---------------------------------------------------------------------------
function bessel_ratio(l::Int, z::Number)
    # Small z: j_l ~ z^l/(2l+1)!!, so j_{l+1}/j_l ~ z/(2l+3)
    if abs(z) < 0.1
        return z / (2l + 3)
    end
    # Large Im(z): sin(z) overflows for |Im(z)| > ~709.
    # Asymptotic: j_{l+1}/j_l → i (or -i) as Im(z) → ±∞.
    y = imag(z)
    y >  300 && return  im
    y < -300 && return -im
    # Miller backward ratio recurrence — stable for ALL l and complex z.
    # Upward ratio recurrence fails when l ≳ |z| (upward recurrence captures
    # the dominant Neumann solution instead of j_l in that regime).
    return _bessel_ratio_miller(l, z)
end

# Miller backward ratio recurrence for j_{l+1}(z)/j_l(z).
#
# Three-term recurrence: j_{n-1} = (2n+1)/z * j_n - j_{n+1}.
# Divide by j_n and define rho_n = j_n/j_{n+1}:
#   rho_{n-1} = (2n+1)/z - 1/rho_n
# Start from 1/rho_L ≈ 0 (j_{L+1} ≈ 0 for L >> l), run down to n = l+1.
# On exit, inv_rho = 1/rho_l = j_{l+1}/j_l.
# Never overflows: only ratio updates; convergence is geometric in (L-l).
function _bessel_ratio_miller(l::Int, z::Number)
    L = l + max(100, l + 1)    # L >> l guarantees convergence
    inv_rho = complex(0.0)     # 1/rho_L ≈ 0
    for n in L:-1:l+1
        rho     = (2n + 1)/z - inv_rho
        inv_rho = 1 / rho
    end
    return inv_rho    # = j_{l+1}/j_l
end

# ---------------------------------------------------------------------------
# Oldroyd-B characteristic equation residual.
#
# Convention: w = σ is the complex eigenvalue in capillary-inertial time units,
# with Re(σ) > 0 for a decaying mode (e^{-σt} convention, Re(σ) > 0 ↔ decay).
# The char eq is eq:char_OB from section_oldroydB.tex:
#   α⁴/q⁴ + 1 = (2(l-1)/q*²) * [l + (l+1)*(q* - 2l*Q)/(q* - 2*Q)]
# with α² = σ_{l;0}/Oh, q² = σ/Oh, q*² = q²*(α²-De1*q²)/(α²-β_s*De1*q²).
# De1 = λ₁·σ_{l;0} is in units of (capillary time)^{-1}.
# ---------------------------------------------------------------------------
function ob_char_eq_residual(w, Oh, De1, beta_s, l)
    sigma_l0 = sqrt(Float64(l * (l - 1) * (l + 2)))
    alpha2   = sigma_l0 / Oh
    q2       = w / Oh
    q2star   = q2 * (alpha2 - De1 * q2) / (alpha2 - beta_s * De1 * q2)
    qstar    = sqrt(q2star)   # complex square root

    Q = bessel_ratio(l, qstar)   # j_{l+1}(q*)/j_l(q*), stable for all q*

    lhs = (sigma_l0 / w)^2 + 1
    rhs = (2 * (l - 1) / q2star) * (l + (l + 1) * (qstar - 2 * l * Q) / (qstar - 2 * Q))
    return lhs - rhs
end

# ---------------------------------------------------------------------------
# Newton root-finder in the complex plane.
# Returns σ with Re(σ) = γ (positive decay rate) and Im(σ) = ω (frequency).
# Initial guess: Lamb limit in simulation units.
# ---------------------------------------------------------------------------
function find_ob_eigenvalue(Oh, De1, beta_s, l)
    sigma_l0 = sqrt(Float64(l * (l - 1) * (l + 2)))
    gamma_lamb = Float64((l - 1) * (2l + 1)) * Oh   # Lamb decay rate (positive)
    w = complex(gamma_lamb, sigma_l0)                # initial guess
    for _ in 1:500
        F  = ob_char_eq_residual(w, Oh, De1, beta_s, l)
        # Scale eps_c to |w| so the finite-difference remains accurate for
        # large eigenvalues (l >> 1 gives |w| ~ l², so a fixed eps_c causes
        # cancellation at O(eps_c/|w|) ~ machine epsilon).
        eps_c = max(1e-7, abs(w) * 1e-7)
        dF = (ob_char_eq_residual(w + eps_c * im, Oh, De1, beta_s, l) - F) / (eps_c * im)
        step = F / dF
        w -= step
        abs(F) < 1e-12 && break
    end
    return w   # Re(w) = γ (decay rate), Im(w) = ω (oscillation frequency)
end

# ---------------------------------------------------------------------------
# Helper: extract γ and ω from A₂ time series.
# γ > 0 (amplitude decays), ω > 0 (oscillation angular frequency).
# ---------------------------------------------------------------------------
function extract_decay_freq(times, A2)
    gamma = -log(abs(A2[end]) / abs(A2[1])) / (times[end] - times[1])

    sign_changes = findall(i -> A2[i] * A2[i+1] < 0, 1:length(A2)-1)
    omega = NaN
    if length(sign_changes) >= 4
        hp = diff(times[sign_changes])
        omega = π / (sum(hp) / length(hp))
    end
    return gamma, omega
end

# ---------------------------------------------------------------------------
# Shared simulation runner
# ---------------------------------------------------------------------------
function run_ob_sim(Oh, De1, beta_s; l=2, M=6, A2_init=0.05, t_end_periods=6)
    Bo = 1e-6
    sigma_exact = find_ob_eigenvalue(Oh, De1, beta_s, l)   # σ in simulation units
    omega_exact = imag(sigma_exact)

    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = make_dt_max(M)
    cfg       = SimConstants(M, M + 1, Oh, Bo, theta_vec, precomp, dt_max)
    ob        = OBParams(De1, beta_s)

    init      = DropState(M)
    init.A[2] = A2_init
    init.z    = 2.0
    init.v    = 0.0
    init.dt   = dt_max
    init.cp   = 0

    T_period = 2π / omega_exact
    t_end    = t_end_periods * T_period

    times, states = solve_drop!(cfg, ob, init;
                                t_end      = t_end,
                                save_every = T_period / 50,
                                dt_init    = dt_max)

    A2 = [s.A[2] for s in states]
    return times, A2, sigma_exact
end

# ---------------------------------------------------------------------------
# Tests: OB eigenvalue vs simulation.
#
# Parameters chosen so that the char eq and numerical simulation agree within
# 5% on γ (both use the same physical model in the linear regime).
# Key: use small Oh and moderate De1*(1-beta_s) so nonlinear/truncation errors
# are small and the linear eigenvalue dominates the A₂ time series.
# ---------------------------------------------------------------------------

@testset "OB eigenvalue: Oh=0.02 De1=0.3 beta_s=0.7" begin
    # Low Oh (Newtonian part well-resolved), moderate viscoelasticity
    Oh = 0.02; De1 = 0.3; beta_s = 0.7

    times, A2, sigma = run_ob_sim(Oh, De1, beta_s)

    gamma_exact = real(sigma)   # decay rate (positive)
    omega_exact = imag(sigma)   # oscillation frequency

    gamma_sim, omega_sim = extract_decay_freq(times, A2)

    @test abs(gamma_sim - gamma_exact) / gamma_exact < 0.05   # 5% tolerance
    if !isnan(omega_sim)
        @test abs(omega_sim - omega_exact) / omega_exact < 0.02   # 2% tolerance
    end
end

@testset "OB eigenvalue: Oh=0.03 De1=0.5 beta_s=0.8" begin
    # Medium Oh, moderate De1, high solvent fraction
    Oh = 0.03; De1 = 0.5; beta_s = 0.8

    times, A2, sigma = run_ob_sim(Oh, De1, beta_s)

    gamma_exact = real(sigma)
    omega_exact = imag(sigma)

    gamma_sim, omega_sim = extract_decay_freq(times, A2)

    @test abs(gamma_sim - gamma_exact) / gamma_exact < 0.05
    if !isnan(omega_sim)
        @test abs(omega_sim - omega_exact) / omega_exact < 0.02
    end
end

@testset "OB eigenvalue: Oh=0.05 De1=0.3 beta_s=0.7" begin
    # Higher Oh, same De1/beta_s structure as case 1
    Oh = 0.05; De1 = 0.3; beta_s = 0.7

    times, A2, sigma = run_ob_sim(Oh, De1, beta_s)

    gamma_exact = real(sigma)
    omega_exact = imag(sigma)

    gamma_sim, omega_sim = extract_decay_freq(times, A2)

    @test abs(gamma_sim - gamma_exact) / gamma_exact < 0.05
    if !isnan(omega_sim)
        @test abs(omega_sim - omega_exact) / omega_exact < 0.02
    end
end

# ---------------------------------------------------------------------------
# Convergence test: OB decay rate converges to the exact Newtonian char-eq
# value as De1 → 0 (not the Lamb approximation, the exact char eq root).
# ---------------------------------------------------------------------------
@testset "OB eigenvalue → Newtonian as De1→0" begin
    Oh = 0.05; beta_s = 0.7; l = 2

    # Exact Newtonian char-eq eigenvalue (De1=0 limit)
    sigma_N = find_ob_eigenvalue(Oh, 0.0, 1.0, l)
    gamma_N = real(sigma_N)

    # OB eigenvalue should approach gamma_N as De1 → 0
    prev_err = Inf
    for De1 in [0.3, 0.1, 0.01]
        sigma = find_ob_eigenvalue(Oh, De1, beta_s, l)
        gamma_ob = real(sigma)
        err = abs(gamma_ob - gamma_N) / gamma_N
        # Errors should decrease monotonically as De1 → 0
        @test err < prev_err || err < 0.01
        prev_err = err
    end
end

# ---------------------------------------------------------------------------
# Spherical Bessel function correctness and stability tests.
# ---------------------------------------------------------------------------

@testset "sph_bessel_j: small-z power series (no catastrophic cancellation)" begin
    # For l≥2, the upward recurrence cancels terms of magnitude ~1/z^l while
    # the result is O(z^l). The power series branch must be active for |z| < 0.1.
    z = 1e-4 + 0im
    @test sph_bessel_j(0, z) ≈ 1.0      atol=1e-8       # j_0 → 1
    @test sph_bessel_j(1, z) ≈ z/3      atol=1e-12      # j_1 ~ z/3
    @test sph_bessel_j(2, z) ≈ z^2/15   atol=1e-24      # j_2 ~ z²/15
    @test sph_bessel_j(3, z) ≈ z^3/105  atol=1e-36      # j_3 ~ z³/105

    # Verify the recurrence j_{l+1} = (2l+1)/z·j_l - j_{l-1} holds for moderate z
    z = 2.5 + 0.5im
    j0, j1, j2, j3 = [sph_bessel_j(l, z) for l in 0:3]
    @test j2 ≈ 3/z * j1 - j0   rtol=1e-12   # n=1 recurrence
    @test j3 ≈ 5/z * j2 - j1   rtol=1e-12   # n=2 recurrence
end

@testset "sph_bessel_j: known exact values" begin
    # j_0(nπ) = sin(nπ)/(nπ) = 0
    @test abs(sph_bessel_j(0, π))   < 1e-14
    @test abs(sph_bessel_j(0, 2π))  < 1e-14
    # j_1(π) = sin(π)/π² - cos(π)/π = 0 + 1/π = 1/π
    @test sph_bessel_j(1, π) ≈ 1/π   rtol=1e-12
end

@testset "bessel_ratio: large imaginary argument (overflow protection)" begin
    # sin(z) = sin(x)cosh(y) + i·cos(x)sinh(y) overflows for Im(z) > ~709.
    # bessel_ratio must return a finite result instead of NaN for all Im(z).
    for t in [100.0, 350.0, 710.0, 800.0, 1000.0]
        Q = bessel_ratio(2, complex(0.0, t))
        @test !isnan(real(Q)) && !isnan(imag(Q))
        @test !isinf(real(Q)) && !isinf(imag(Q))
        # Q → i as Im(z) → +∞; the ratio recurrence gives O(1/t) error so
        # use a tolerance commensurate with 1/t (the asymptotic guard at t≥300 gives exact i).
        @test abs(Q - im) < max(3/t, 1e-10)
    end
    # Negative imaginary argument → -i
    Q_neg = bessel_ratio(2, complex(0.0, -400.0))
    @test abs(Q_neg + im) < 1e-10
end

@testset "bessel_ratio: small z" begin
    # j_{l+1}/j_l ~ z/(2l+3) for small z
    z = 1e-5 + 0im
    @test bessel_ratio(0, z) ≈ z/3   rtol=1e-8
    @test bessel_ratio(1, z) ≈ z/5   rtol=1e-8
    @test bessel_ratio(2, z) ≈ z/7   rtol=1e-8
end

@testset "bessel_ratio: consistent with sph_bessel_j at moderate z" begin
    # For moderate z (no overflow, no cancellation), both methods must agree.
    for z in [1.0+0im, 2.0+1.0im, 0.5+2.0im, 3.14+0.1im]
        for l in [0, 1, 2, 3]
            Q_ratio  = bessel_ratio(l, z)
            Q_direct = sph_bessel_j(l+1, z) / sph_bessel_j(l, z)
            @test Q_ratio ≈ Q_direct   rtol=1e-8
        end
    end
end

# ---------------------------------------------------------------------------
# Newton vs grid search: verify Newton converges to the correct root.
#
# For each parameter set we run two independent solvers:
#   1. Newton from the Lamb initial guess (find_ob_eigenvalue).
#   2. Coarse grid search over (γ, ω) → Newton refinement from the
#      grid minimum (no prior knowledge of the root location).
# Both must reach the same root to machine precision.  Agreement confirms
# that Newton does not converge to a spurious local minimum and that the
# physical root is unique in the search region for all tested regimes.
# ---------------------------------------------------------------------------

function _grid_then_newton(Oh, De1, beta_s, l)
    sigma_l0   = sqrt(Float64(l * (l - 1) * (l + 2)))
    gamma_lamb = Float64((l - 1) * (2l + 1)) * Oh

    # Coarse grid: γ ∈ [0.2, 3.0] × γ_lamb, ω ∈ [0.6, 1.3] × σ_{l;0}
    best_F = Inf; best_w = complex(gamma_lamb, sigma_l0)
    for g_fac in range(0.2, 3.0; length=12)
        for o_fac in range(0.6, 1.3; length=12)
            w0 = complex(g_fac * gamma_lamb, o_fac * sigma_l0)
            f  = abs(ob_char_eq_residual(w0, Oh, De1, beta_s, l))
            if f < best_F
                best_F = f; best_w = w0
            end
        end
    end

    # Newton from the grid minimum
    w = best_w
    for _ in 1:500
        F    = ob_char_eq_residual(w, Oh, De1, beta_s, l)
        eps_c = max(1e-7, abs(w) * 1e-7)
        dF   = (ob_char_eq_residual(w + eps_c * im, Oh, De1, beta_s, l) - F) / (eps_c * im)
        step = F / dF
        w   -= step
        abs(F) < 1e-12 && break
    end
    return w
end

@testset "Newton converges to same root as grid search" begin
    # Parameters span Newtonian (De1=0) and Oldroyd-B (De1>0), low and moderate
    # Oh, and mode numbers l=2..10.  For each case both solvers must agree on
    # γ and ω to relative tolerance 1e-4 (they both apply Newton at the end,
    # so exact agreement implies they found the same basin).
    cases = [
        # (Oh,   De1,  beta_s, l)   description
        (0.02,  0.0,  1.0,  2),   # Newtonian baseline
        (0.05,  0.0,  1.0,  2),   # Newtonian, higher Oh
        (0.02,  0.0,  1.0,  3),   # Newtonian, l=3
        (0.02,  0.0,  1.0,  5),   # Newtonian, l=5
        (0.02,  0.0,  1.0,  10),  # Newtonian, l=10
        (0.02,  0.3,  0.7,  2),   # Oldroyd-B, moderate De1
        (0.03,  0.5,  0.8,  2),   # Oldroyd-B, higher De1
        (0.05,  0.3,  0.7,  2),   # Oldroyd-B, higher Oh
        (0.02,  0.1,  0.5,  3),   # Oldroyd-B, l=3
        (0.02,  0.3,  0.7,  5),   # Oldroyd-B, l=5
        (0.02,  0.3,  0.7,  10),  # Oldroyd-B, l=10
        (0.01,  0.2,  0.6,  10),  # Oldroyd-B, l=10, lower Oh
    ]
    for (Oh, De1, beta_s, l) in cases
        w_newton = find_ob_eigenvalue(Oh, De1, beta_s, l)
        w_grid   = _grid_then_newton(Oh, De1, beta_s, l)
        @test abs(w_newton - w_grid) / abs(w_newton) < 1e-4
    end
end

@testset "char eq residual ≈ 0 at the computed eigenvalue" begin
    # The Newton root-finder must satisfy |F(w*)| < 1e-8 at its own output.
    for (Oh, De1, beta_s) in [(0.02, 0.3, 0.7), (0.03, 0.5, 0.8), (0.05, 0.3, 0.7)]
        w = find_ob_eigenvalue(Oh, De1, beta_s, 2)
        res = ob_char_eq_residual(w, Oh, De1, beta_s, 2)
        @test abs(res) < 1e-8
        @test isfinite(real(w)) && isfinite(imag(w))
        @test real(w) > 0   # positive decay rate (damped mode)
        @test imag(w) > 0   # positive frequency (oscillatory)
    end
end

# ---------------------------------------------------------------------------
# High-mode-number coverage: Bessel ratio stability and char eq correctness.
#
# Two gaps that the l=2-only tests above leave open:
#   1. bessel_ratio(l, qstar) is called with l up to M in the solver.  For
#      large l the ratio recurrence runs l steps starting from r₁ = j₁/j₀.
#      We must confirm it stays finite and the char eq root exists for l≤200.
#   2. The residual damping D2[k] = 2Oh(k-1)(2k+1) and capillary D1[k] =
#      k(k-1)(k+2) are polynomial in k, but the simulation has only been
#      compared against theory at l=2.  Eigenvalue comparison at l=3,5,10
#      confirms the code is correct for higher modes.
# ---------------------------------------------------------------------------

@testset "char eq residual ≈ 0 for Newtonian l up to 200" begin
    # For each l, Newton must converge to a root with |F(w*)| < 1e-8.
    # This simultaneously tests:
    #   • bessel_ratio(l, qstar) is finite/stable for l up to 200 iterations
    #   • find_ob_eigenvalue converges (Lamb initial guess is good enough)
    #   • ob_char_eq_residual formula is correct for all these mode numbers
    #
    # Oh=0.001: small enough that Oh·(2l+1) ≤ 0.4 for l≤200, so the Lamb
    # initial guess is within Newton's basin of attraction for all l.
    # For l ≥ 50 the argument Im(qstar) > 300, so bessel_ratio returns im
    # analytically and the char eq reduces exactly to the Lamb quadratic.
    Oh = 0.001
    for l in [2, 3, 5, 10, 20, 50, 100, 200]
        w = find_ob_eigenvalue(Oh, 0.0, 1.0, l)
        res = ob_char_eq_residual(w, Oh, 0.0, 1.0, l)
        @test abs(res) < 1e-8
        @test isfinite(real(w)) && isfinite(imag(w))
        @test real(w) > 0
        @test imag(w) > 0
        gamma_lamb = Float64((l - 1) * (2l + 1)) * Oh
        @test abs(real(w) - gamma_lamb) / gamma_lamb < 0.20
    end
end

@testset "Newtonian eigenvalue from simulation: l=3, l=5, l=10" begin
    # Run a single-mode free oscillation for each l (M=l, only A[l] excited).
    # The simulation implements D1[k]=k(k-1)(k+2) and D2[k]=2Oh(k-1)(2k+1),
    # which are the Lamb approximation coefficients.  The correct reference is
    # therefore gamma_lamb=(l-1)(2l+1)Oh and omega_lamb=sqrt(l(l-1)(l+2)),
    # not the exact Newtonian char-eq eigenvalue (which includes corrections
    # from the full Bessel-function structure beyond the Lamb formula).
    Oh = 0.02; Bo = 1e-6
    for l in [3, 5, 10]
        gamma_lamb = Float64((l - 1) * (2l + 1)) * Oh
        sigma_l0   = sqrt(Float64(l * (l - 1) * (l + 2)))

        M         = l
        theta_vec = make_theta_vec(M)
        precomp   = precompute_integrals(NaN, M)[1]
        # BDF2 at 8 steps/period (make_dt_max default) gives ωh ≈ 0.79 rad,
        # which produces O(1) numerical damping — unusable for eigenvalue tests.
        # 40 steps/period (ωh ≈ 0.16) reduces the BDF2 truncation error in γ
        # to ~1%, well inside the 5% tolerance.
        dt_osc    = 2π / (sigma_l0 * 40)
        cfg       = SimConstants(M, M + 1, Oh, Bo, theta_vec, precomp, dt_osc)

        init      = DropState(M)
        init.A[l] = 0.05        # excite only mode l; all other modes start at 0
        init.z    = 2.0
        init.dt   = dt_osc

        T_period = 2π / sigma_l0
        times, states = solve_drop!(cfg, OBParams(), init;
                                    t_end      = 6 * T_period,
                                    save_every = T_period / 50,
                                    dt_init    = dt_osc)

        Al = [s.A[l] for s in states]
        gamma_sim, omega_sim = extract_decay_freq(times, Al)

        @test abs(gamma_sim - gamma_lamb) / gamma_lamb < 0.05
        if !isnan(omega_sim)
            @test abs(omega_sim - sigma_l0) / sigma_l0 < 0.02
        end
    end
end

# `run_impact` is the documented entry point, and the solver-choice page routes Oldroyd-B
# users to it -- but until the `ob` keyword existed it hardwired `OBParams()`, so every
# "viscoelastic" impact launched through it was silently Newtonian. These lock the route.
@testset "run_impact carries Oldroyd-B parameters through" begin
    b  = Backend(formulation = :nonvariational, contact = :tangency)
    kw = (We = 0.5, Bo = 0.019, Oh = 0.3038, M = 20, K = 2, t_max = 12.0)
    newt = run_impact(b; kw...)
    obr  = run_impact(b; kw..., ob = OBParams(0.6, 0.4))
    @test newt.ok && obr.ok
    ## a polymer with De1 = 0.6 and only 40% solvent must not reproduce the Newtonian bounce
    @test abs(obr.cor - newt.cor) > 1e-6
    ## and the Newtonian limits of the parameter must reproduce it.
    ##
    ## The tolerance is 1e-8, not machine epsilon: the assembly is threaded, so reduction
    ## order varies between runs and the same case reproduces to about 1e-11 rather than
    ## bitwise. That still separates the two cases by three orders of magnitude, since a
    ## genuine polymer moves the restitution by more than 1e-6 (asserted above).
    for p in (OBParams(0.0, 1.0), OBParams(0.0, 0.4), OBParams(0.6, 1.0))
        @test isapprox(run_impact(b; kw..., ob = p).cor, newt.cor; rtol = 1e-8)
    end
    ## the variational formulation has no polymer state and must refuse rather than ignore
    @test_throws ErrorException run_impact(Backend(); kw..., ob = OBParams(0.6, 0.4))
end
