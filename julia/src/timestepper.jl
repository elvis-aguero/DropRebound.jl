"""
    solve_drop!(cfg, ob, init; t_end, save_every, dt_init) → (times, states)

Main adaptive BDF1/BDF2 time-integration loop with contact search.
At each step:
  1. Solve with cp=cp_prev (current contact count)
  2. If cp=0 result shows penetration (south pole below z=0), retry with cp=1
     using the penetrated state as the initial guess (warm start)
  3. Also try cp_prev±1 for smooth contact transitions
  4. Accept the cp with smallest |contact_error|; reject penetrating cp=0 states
  5. If all fail: halve dt and retry
  6. Ramp dt back toward dt_max at +10% per step
"""
function solve_drop!(cfg::SimConstants, ob::OBParams, init::DropState;
                     t_end::Float64      = 10.0,
                     save_every::Float64 = 0.1,
                     dt_init::Float64    = cfg.dt_max,
                     dt_min::Float64     = 1e-6)

    is_ob = ob.De1 > 0.0 && ob.beta_s < 1.0
    pack_fn      = is_ob ? pack_X_ob     : (s, M) -> pack_X(s, M)
    unpack_fn!   = is_ob ? unpack_X_ob!  : (s, X, M) -> unpack_X!(s, X, M)
    residual_fn! = is_ob ? build_residual_ob! : build_residual!
    jacobian_fn  = is_ob ? build_jacobian_ob  : build_jacobian

    rheol = is_ob ? "Oldroyd-B" : "Newtonian"
    @info "solve_drop! starting" M=cfg.M Oh=cfg.Oh Bo=cfg.Bo t_end=t_end rheology=rheol dt_max=cfg.dt_max

    history    = DropState[deepcopy(init)]
    dt         = dt_init
    t          = init.t
    prev_cp    = init.cp

    saved_times  = Float64[t]
    saved_states = DropState[deepcopy(init)]
    next_save    = t + save_every
    step_count   = 0

    while t < t_end
        order   = min(length(history), 2)
        cp_prev = history[end].cp

        best_state = nothing
        best_err   = Inf

        # Collect candidate initial guesses: history[end] for each cp, plus
        # a warm-start from cp=0 solution when penetration is detected
        warm_start_X = nothing   # warm-start X for cp=1 if cp=0 penetrates

        for cp in max(0, cp_prev - 1) : cp_prev + 1
            # Build history slice for BDF
            hist_slice = order == 2 ? history[end-1:end] : history[end:end]

            # Use warm-start if available (cp=0 penetrated → use as guess for cp=1)
            if cp == 1 && warm_start_X !== nothing
                X0 = warm_start_X
            else
                X0 = pack_fn(history[end], cfg.M)
            end

            R! = (buf, Xv) -> begin
                s = deepcopy(history[end])
                unpack_fn!(s, Xv, cfg.M)
                s.cp = cp
                fill!(buf, 0.0)
                residual_fn!(buf, s, hist_slice, dt, cp, cfg, ob)
            end

            J_fn = Xv -> begin
                s = deepcopy(history[end])
                unpack_fn!(s, Xv, cfg.M)
                s.cp = cp
                jacobian_fn(s, hist_slice, dt, cp, cfg, ob)
            end

            X   = copy(X0)
            key = (cfg.M, cp, round(dt; sigdigits=6), order, is_ob)
            converged = newton_solve!(X, R!, J_fn; cache_key=key)

            if converged
                candidate     = deepcopy(history[end])
                unpack_fn!(candidate, X, cfg.M)
                candidate.t   = t + dt
                candidate.dt  = dt
                candidate.cp  = cp
                err = abs(contact_error(candidate, cfg.theta_vec, cp))
                # Penetration guard: if cp=0 but south pole is below substrate, reject
                # and save the penetrated state as warm-start for cp=1
                if cp == 0 && drop_height(candidate, cfg.theta_vec[1]) < 0.0
                    warm_start_X = copy(X)   # use penetrated solution as warm start
                    err = Inf
                end
                if err < best_err
                    best_err   = err
                    best_state = candidate
                end
            end
        end

        if best_state === nothing || best_err == Inf
            dt /= 2
            @debug "Step rejected: halving dt" t=t new_dt=dt order=order cp_prev=cp_prev
            if dt < 10 * dt_min
                @warn "dt approaching minimum — possible numerical difficulty" t=t dt=dt dt_min=dt_min
            end
            if dt < dt_min
                error("Time step below dt_min=$dt_min at t=$t (dt=$dt); results would be numerically invalid")
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

        # Log contact events (cp changes)
        new_cp = best_state.cp
        if new_cp != prev_cp
            if new_cp > 0 && prev_cp == 0
                @info "Contact onset" t=t cp=new_cp
            elseif new_cp == 0 && prev_cp > 0
                @info "Lift-off" t=t cp_was=prev_cp
            else
                @info "Contact count changed" t=t cp_prev=prev_cp cp_new=new_cp
            end
            prev_cp = new_cp
        end

        @debug "Step accepted" t=t dt=dt cp=new_cp contact_err=best_err step=step_count

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

    @info "solve_drop! complete" steps=step_count t_final=t saved_frames=length(saved_times)
    return saved_times, saved_states
end
