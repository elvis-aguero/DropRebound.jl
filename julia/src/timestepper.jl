"""
    solve_drop!(cfg, ob, init; t_end, save_every, dt_init) → (times, states)

Main adaptive BDF1/BDF2 time-integration loop with contact search.
At each step:
  1. Try cp in {cp_prev-1, cp_prev, cp_prev+1}
  2. For each cp, Newton-solve for the next state
  3. Accept the cp with smallest |contact_error|
  4. If all fail: halve dt and retry
  5. Ramp dt back toward dt_max at +10% per step
"""
function solve_drop!(cfg::SimConstants, ob::OBParams, init::DropState;
                     t_end::Float64      = 10.0,
                     save_every::Float64 = 0.1,
                     dt_init::Float64    = cfg.dt_max)

    history    = DropState[deepcopy(init)]
    dt         = dt_init
    t          = init.t

    saved_times  = Float64[t]
    saved_states = DropState[deepcopy(init)]
    next_save    = t + save_every
    step_count   = 0

    while t < t_end
        order   = min(length(history), 2)
        cp_prev = history[end].cp

        best_state = nothing
        best_err   = Inf

        for cp in max(0, cp_prev - 1) : cp_prev + 1
            # Build history slice for BDF
            hist_slice = order == 2 ? history[end-1:end] : history[end:end]

            X0 = pack_X(history[end], cfg.M)

            R! = (buf, Xv) -> begin
                s = deepcopy(history[end])
                unpack_X!(s, Xv, cfg.M)
                s.cp = cp
                fill!(buf, 0.0)
                build_residual!(buf, s, hist_slice, dt, cp, cfg, ob)
            end

            J_fn = Xv -> begin
                s = deepcopy(history[end])
                unpack_X!(s, Xv, cfg.M)
                s.cp = cp
                build_jacobian(s, hist_slice, dt, cp, cfg, ob)
            end

            X   = copy(X0)
            key = (cp, round(dt; sigdigits=6), order)
            converged = newton_solve!(X, R!, J_fn; cache_key=key)

            if converged
                candidate     = deepcopy(history[end])
                unpack_X!(candidate, X, cfg.M)
                candidate.t   = t + dt
                candidate.dt  = dt
                candidate.cp  = cp
                err = abs(contact_error(candidate, cfg.theta_vec, cp))
                if err < best_err
                    best_err   = err
                    best_state = candidate
                end
            end
        end

        if best_state === nothing || best_err == Inf
            dt /= 2
            if dt < 1e-14
                error("Time step underflow at t=$t")
            end
            clear_jac_cache!()
            continue
        end

        # Accept step
        t = best_state.t
        push!(history, best_state)
        step_count += 1
        if length(history) > 2
            popfirst!(history)
        end

        # Ramp dt back toward dt_max
        dt = min(dt * 1.1, cfg.dt_max)

        if t >= next_save
            push!(saved_times, t)
            push!(saved_states, deepcopy(best_state))
            next_save = t + save_every
        end
    end

    # Ensure the final state is always saved
    if saved_times[end] < t
        push!(saved_times, t)
        push!(saved_states, deepcopy(history[end]))
    end

    return saved_times, saved_states
end
