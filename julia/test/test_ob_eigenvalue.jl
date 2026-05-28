using Test
using DropSolver

# ---------------------------------------------------------------------------
# Spherical Bessel function j_l(z) via upward recurrence (works for complex z)
# ---------------------------------------------------------------------------
function sph_bessel_j(l::Int, z::Number)
    if abs(z) < 1e-300
        return l == 0 ? complex(1.0) : complex(0.0)
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

    jl_q  = sph_bessel_j(l,     qstar)
    jl1_q = sph_bessel_j(l + 1, qstar)
    Q     = jl1_q / jl_q   # Q_{l+1/2}(qstar)

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
    for _ in 1:200
        F  = ob_char_eq_residual(w, Oh, De1, beta_s, l)
        eps_c = 1e-7
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
