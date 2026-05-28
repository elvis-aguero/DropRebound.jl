#!/usr/bin/env julia
# Sweep OB eigenvalue vs simulation: paper-quality table.
# Parameters chosen so that linear eigenvalue theory (char eq) agrees with
# the numerical simulation to within a few percent.
# Usage: julia --project=.. scripts/run_eigenvalue_sweep.jl

using Pkg; Pkg.activate(joinpath(@__DIR__, ".."))
using DropSolver
using Printf

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
# w = σ (complex eigenvalue in capillary-inertial units, Re(σ) > 0 for decay).
# De1 = λ₁·σ_{l;0} (Deborah number referenced to inviscid oscillation frequency).
# ---------------------------------------------------------------------------
function ob_char_eq_residual(w, Oh, De1, beta_s, l)
    sigma_l0 = sqrt(Float64(l * (l - 1) * (l + 2)))
    alpha2   = sigma_l0 / Oh
    q2       = w / Oh
    q2star   = q2 * (alpha2 - De1 * q2) / (alpha2 - beta_s * De1 * q2)
    qstar    = sqrt(q2star)

    jl_q  = sph_bessel_j(l,     qstar)
    jl1_q = sph_bessel_j(l + 1, qstar)
    Q     = jl1_q / jl_q

    lhs = (sigma_l0 / w)^2 + 1
    rhs = (2 * (l - 1) / q2star) * (l + (l + 1) * (qstar - 2 * l * Q) / (qstar - 2 * Q))
    return lhs - rhs
end

# ---------------------------------------------------------------------------
# Newton root-finder. Returns σ: Re(σ) = γ (decay rate), Im(σ) = ω (freq).
# ---------------------------------------------------------------------------
function find_ob_eigenvalue(Oh, De1, beta_s, l)
    sigma_l0 = sqrt(Float64(l * (l - 1) * (l + 2)))
    gamma_lamb = Float64((l - 1) * (2l + 1)) * Oh
    w = complex(gamma_lamb, sigma_l0)
    for _ in 1:200
        F  = ob_char_eq_residual(w, Oh, De1, beta_s, l)
        eps_c = 1e-7
        dF = (ob_char_eq_residual(w + eps_c * im, Oh, De1, beta_s, l) - F) / (eps_c * im)
        step = F / dF
        w -= step
        abs(F) < 1e-12 && break
    end
    return w
end

# ---------------------------------------------------------------------------
# Run one simulation and extract γ and ω
# ---------------------------------------------------------------------------
function run_one(Oh, De1, beta_s; l=2, M=6, t_end_periods=6)
    Fr = 1e6
    sigma_exact = find_ob_eigenvalue(Oh, De1, beta_s, l)
    omega_ex    = imag(sigma_exact)
    gamma_ex    = real(sigma_exact)

    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    dt_max    = make_dt_max(M)
    cfg       = SimConstants(M, M + 1, Oh, Fr, theta_vec, precomp, dt_max)
    ob        = OBParams(De1, beta_s)

    init      = DropState(M)
    init.A[2] = 0.05
    init.z    = 2.0
    init.v    = 0.0
    init.dt   = dt_max
    init.cp   = 0

    T_period = 2π / omega_ex
    t_end    = t_end_periods * T_period

    times, states = solve_drop!(cfg, ob, init;
                                t_end      = t_end,
                                save_every = T_period / 50,
                                dt_init    = dt_max)

    A2 = [s.A[2] for s in states]

    gamma_sim = -log(abs(A2[end]) / abs(A2[1])) / (times[end] - times[1])

    sign_changes = findall(i -> A2[i] * A2[i+1] < 0, 1:length(A2)-1)
    omega_sim = NaN
    if length(sign_changes) >= 4
        hp = diff(times[sign_changes])
        omega_sim = π / (sum(hp) / length(hp))
    end

    return gamma_ex, omega_ex, gamma_sim, omega_sim
end

# ---------------------------------------------------------------------------
# Sweep over a parameter grid.
# Oh values: 0.02, 0.03, 0.05 (underdamped regime, char eq valid)
# De1 values: 0.3, 0.5, 0.8 (moderate elasticity, De1*(1-beta_s) small)
# beta_s values: 0.7, 0.8 (high solvent fraction for convergence)
# ---------------------------------------------------------------------------
@printf("%-6s  %-5s  %-6s  %-11s  %-11s  %-11s  %-11s  %-10s  %-10s\n",
        "Oh", "De1", "beta_s", "γ_exact", "ω_exact", "γ_sim", "ω_sim",
        "err_γ%", "err_ω%")
println("-"^100)

for Oh in [0.02, 0.03, 0.05]
    for De1 in [0.3, 0.5, 0.8]
        for beta_s in [0.7, 0.8]
            gamma_ex, omega_ex, gamma_sim, omega_sim = run_one(Oh, De1, beta_s)

            err_gamma = isnan(gamma_sim) ? NaN : 100 * abs(gamma_sim - gamma_ex) / gamma_ex
            err_omega = isnan(omega_sim) ? NaN : 100 * abs(omega_sim - omega_ex) / omega_ex

            if isnan(omega_sim)
                @printf("%-6.3f  %-5.1f  %-6.3f  %-11.5f  %-11.5f  %-11.5f  %-11s  %-10.2f  %-10s\n",
                        Oh, De1, beta_s, gamma_ex, omega_ex, gamma_sim, "N/A", err_gamma, "N/A")
            else
                @printf("%-6.3f  %-5.1f  %-6.3f  %-11.5f  %-11.5f  %-11.5f  %-11.5f  %-10.2f  %-10.2f\n",
                        Oh, De1, beta_s, gamma_ex, omega_ex, gamma_sim, omega_sim, err_gamma, err_omega)
            end
        end
    end
end
