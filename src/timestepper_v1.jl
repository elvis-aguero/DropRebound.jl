"""
    solve_drop_v1!(cfg, init; t_end, save_every, dt_init) → (times, states)

v1 linearized timestepper with continuous contact angle θ*(t).

At each step:
  1. Build the 3M+1 linear system (v1 physics: full contact-force integral,
     continuous θ* from previous state, z_prev frozen in R7 coefficient).
  2. Solve directly (J \\ (-R(prev)) — one shot, no Newton iterations).
  3. Update θ* via scalar bisection + Newton polish on h(θ*) = 0.
  4. If free-surface penetration detected, halve dt and retry.
"""
function solve_drop_v1!(cfg::SimConstants, init::DropState;
                        t_end::Float64      = 10.0,
                        save_every::Float64 = 0.1,
                        dt_init::Float64    = cfg.dt_max)

    history    = DropState[deepcopy(init)]
    dt         = dt_init
    t          = init.t

    saved_times  = Float64[t]
    saved_states = DropState[deepcopy(init)]
    next_save    = t + save_every

    while t < t_end
        order    = min(length(history), 2)
        hist_slice = order == 2 ? history[end-1:end] : history[end:end]

        # Build and solve linear system in one shot
        candidate = deepcopy(history[end])
        candidate.t  = t + dt
        candidate.dt = dt

        J = build_jacobian_v1(history[end], hist_slice, dt, cfg)

        R = zeros(3cfg.M + 1)
        build_residual_v1!(R, history[end], hist_slice, dt, cfg)

        # Solve J·ΔX = -R(x_prev), then x_new = x_prev + ΔX.
        # For a linear system this is equivalent to solving J·x = rhs where
        # rhs absorbs the history terms (already in R via BDF coefficients).
        ΔX = J \ (-R)
        X_new = pack_X(history[end], cfg.M) .+ ΔX
        unpack_X!(candidate, X_new, cfg.M)

        # Update continuous contact angle
        update_theta_star!(candidate)

        # Validate: no free-surface penetration
        valid = true
        θs = candidate.theta_star
        if θs < π - 1e-10
            # Check that the drop surface is above the substrate on the free zone
            u_check, _ = gauss_legendre_nodes(10, cos(θs), 1.0)
            θ_check = acos.(clamp.(u_check, -1.0, 1.0))
            for θc in θ_check
                if drop_height(candidate, θc) < -1e-6
                    valid = false
                    break
                end
            end
        end

        if !valid
            dt /= 2
            if dt < 1e-14
                error("Time step underflow at t=$t")
            end
            continue
        end

        # Accept step
        t = candidate.t
        push!(history, candidate)
        if length(history) > 2
            popfirst!(history)
        end

        dt = min(dt * 1.1, cfg.dt_max)

        if t >= next_save
            push!(saved_times, t)
            push!(saved_states, deepcopy(candidate))
            next_save = t + save_every
        end
    end

    if saved_times[end] < t
        push!(saved_times, t)
        push!(saved_states, deepcopy(history[end]))
    end

    return saved_times, saved_states
end
