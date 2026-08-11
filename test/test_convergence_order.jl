# A posteriori convergence-order tests for the time integrator.
#
# Nothing else in this suite measures order. That mattered: the shear-thinning
# path evaluated its viscosity coefficients at the previous step -- a constant
# extrapolation -- which is a first-order splitting error sitting inside an
# otherwise second-order BDF2 scheme. Every existing test passed throughout,
# because they check values and invariants, not rates.
#
# Both tests below integrate a free oscillation with the drop held far from the
# wall, so no contact events occur and the time integrator is isolated. States
# are compared at a fixed time reached by interpolation, NOT at the final saved
# step: `solve_drop!` overshoots `t_end` by up to one step, so comparing final
# states measures the time offset rather than the truncation error, and reports
# a meaningless order.

const _CO_ETA_INF = 0.0037320997942061666
const _CO_ETA_0   = 8.433817577956766
const _CO_K       = 18.48081673111359
const _CO_M       = 0.7430524574330837
const _CO_BO      = 0.012
const _CO_R       = 0.0003
const _CO_SIGMA   = 0.0728
const _CO_G       = 9.81
const _CO_RHO     = _CO_BO * _CO_SIGMA / (_CO_G * _CO_R^2)
const _CO_TSIG    = sqrt(_CO_RHO * _CO_R^3 / _CO_SIGMA)
const _CO_OH0     = _CO_ETA_0 / sqrt(_CO_RHO * _CO_SIGMA * _CO_R)
const _CO_LAMC    = _CO_K / _CO_TSIG
const _CO_ETAR    = _CO_ETA_INF / _CO_ETA_0

"""
Mode amplitudes at `t_cmp`, linearly interpolated from a trajectory saved every
step. `M = 6`, drop parked at `z = 50` so it never contacts the wall.
"""
function _co_state_at(dt::Float64; thinning::Bool, t_cmp::Float64 = 0.8,
                      n_table::Int = 60)
    stx = thinning ? STExactParams(6, _CO_OH0, _CO_LAMC, _CO_M, 1.0;
                                   viscous = :reid, eta_inf_ratio = _CO_ETAR,
                                   n_table = n_table) : nothing
    Oh  = thinning ? _CO_OH0 : 0.05
    cfg = SimConstants(6, 7, Oh, _CO_BO, make_theta_vec(6),
                       precompute_integrals(NaN, 6)[1], dt; viscous = :reid)
    init = DropState(6)
    init.A[2] = 0.05
    init.z    = 50.0
    init.dt   = dt
    init.cp   = 0
    t, st = solve_drop!(cfg, OBParams(), init; stx = stx, t_end = 1.0,
                        save_every = dt, dt_init = dt)
    i = findlast(x -> x <= t_cmp, t)
    w = (t_cmp - t[i]) / (t[i+1] - t[i])
    (1 - w) .* st[i].A[2:end] .+ w .* st[i+1].A[2:end]
end

"Observed order between the two finest resolutions of a refinement sweep."
function _co_observed_order(; thinning::Bool)
    ref = _co_state_at(1e-2 / 256; thinning = thinning)
    errs = [maximum(abs.(_co_state_at(1e-2 * 2.0^(-k); thinning = thinning) .- ref))
            for k in 0:3]
    (log2(errs[end-1] / errs[end]), errs)
end

@testset "Time integrator: observed convergence order" begin

    @testset "Newtonian control is second order" begin
        # Constant coefficients: BDF2 with nothing split off. If this is not ~2
        # the harness is broken and the shear-thinning result below means
        # nothing, so it is checked first.
        ord, errs = _co_observed_order(thinning = false)
        @test all(errs[k+1] < errs[k] for k in 1:length(errs)-1)
        @test ord > 1.8
        @test ord < 2.2
    end

    @testset "Shear-thinning path does not collapse to first order" begin
        # The coefficients are extrapolated to second order rather than held at
        # the previous step (`_extrapolated_Adot`). Without that the observed
        # order here is ~1.3; with it the sequence rises through ~1.6 toward 2.
        # The bound is deliberately loose -- it guards against a REGRESSION to
        # first-order splitting, which is what a constant extrapolation gives.
        ord, errs = _co_observed_order(thinning = true)
        @test all(errs[k+1] < errs[k] for k in 1:length(errs)-1)
        @test ord >= 1.5
    end

    @testset "Extrapolation is exact for a linear history" begin
        # (1 + r) y_n - r y_{n-1} must reproduce a linearly-varying velocity
        # exactly, for a non-uniform step ratio as well as a uniform one.
        for (dt, dt_prev) in ((0.01, 0.01), (0.01, 0.004), (0.002, 0.005))
            h1 = DropState(4); h2 = DropState(4)
            h1.dt = dt_prev; h2.dt = dt
            slope = [0.3, -0.7, 1.1]
            h1.Adot[2:end] .= 0.0
            h2.Adot[2:end] .= slope .* dt_prev            # linear in t
            got  = DropSolver._extrapolated_Adot([h1, h2], dt)
            want = slope .* (dt_prev + dt)
            @test got ≈ want rtol = 1e-12
        end
    end

    @testset "Extrapolation degrades safely on a one-state history" begin
        # BDF1 startup: nothing to extrapolate from, so the previous value is
        # used. Must not error and must not invent a slope.
        h = DropState(4)
        h.dt = 0.01
        h.Adot[2:end] .= [0.2, -0.5, 0.9]
        @test DropSolver._extrapolated_Adot([h], 0.01) ≈ [0.2, -0.5, 0.9]
    end
end
