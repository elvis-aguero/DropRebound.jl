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

Once in contact, the search window is limited to cp_prev-1 : cp_prev+1. Per
Gabbard et al. (2025, JFM, "Drop rebound at low Weber number") §C.1-C.2: the
contact point count must never change by more than one collocation point in
a single accepted step *while already in contact*. Testing only these three
candidates (rather than a wider window) is what makes the tangency-error
scoring in `contact_error` reliable — a wider window can accept a candidate
that merely scores lowest among a set that includes points where the
linearized height function's scoring is no longer trustworthy, rather than a
genuine local minimum. If no candidate in this window converges without
penetration, the step is rejected and dt is halved, exactly as the reference
algorithm requires — never widened.

Contact *onset* (cp_prev == 0) is the one genuine exception: a smooth convex
body touching a flat plane nucleates a contact patch whose radius initially
grows like √(t − t_onset), i.e. arbitrarily fast in continuum time, so no
achievable dt can resolve the very first instants of contact one collocation
point at a time — this is a real feature of rigid-contact kinematics, not a
discretization artifact. At cp_prev == 0 the full candidate range 0:N_angles-1
is searched and the true error-minimizing cp accepted directly (no incremental
widening, so there is no risk of settling on a merely-locally-best candidate).
  6. Ramp dt back toward dt_max at +10% per step
"""
function solve_drop!(cfg::SimConstants, ob::OBParams, init::DropState;
                     st::STParams        = STParams(),
                     t_end::Float64      = 10.0,
                     save_every::Float64 = 0.1,
                     dt_init::Float64    = cfg.dt_max,
                     dt_min::Float64     = cfg.dt_max * 1e-4)

    is_ob = ob.De1 > 0.0 && ob.beta_s < 1.0
    is_st = !is_ob && st.eps_ST > 0.0 && !isempty(st.Gamma)
    pack_fn      = is_ob ? pack_X_ob     : (s, M) -> pack_X(s, M)
    unpack_fn!   = is_ob ? unpack_X_ob!  : (s, X, M) -> unpack_X!(s, X, M)
    residual_fn! = is_ob ? build_residual_ob! : build_residual!
    jacobian_fn  = is_ob ? build_jacobian_ob  : build_jacobian

    rheol = is_ob ? "Oldroyd-B" : (is_st ? "Carreau" : "Newtonian")
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

        # cp search window: cp_prev-1 : cp_prev+1 while in contact (never
        # widened); the full range at contact onset (cp_prev == 0) — see
        # docstring above.
        cp_upper = cp_prev == 0 ? cfg.N_angles - 1 : min(cfg.N_angles - 1, cp_prev + 1)

        # Collect candidate initial guesses: history[end] for each cp, plus
        # a warm-start from cp=0 solution when penetration is detected
        warm_start_X = nothing   # warm-start X for cp=1 if cp=0 penetrates

        for cp in max(0, cp_prev - 1) : cp_upper
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
                if is_st
                    build_residual_st!(buf, s, hist_slice, dt, cp, cfg, ob, st)
                else
                    residual_fn!(buf, s, hist_slice, dt, cp, cfg, ob)
                end
            end

            J_fn = Xv -> begin
                s = deepcopy(history[end])
                unpack_fn!(s, Xv, cfg.M)
                s.cp = cp
                if is_st
                    build_jacobian_st(s, hist_slice, dt, cp, cfg, ob, st)
                else
                    jacobian_fn(s, hist_slice, dt, cp, cfg, ob)
                end
            end

            X   = copy(X0)
            # Carreau: shear_sq_lag changes each step, so skip Jacobian cache
            key = is_st ? nothing : (cfg.M, cp, round(dt; sigdigits=6), order, is_ob)
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
            @debug "Step rejected" t=t dt=dt cp_prev=cp_prev

            # No candidate in cp_prev-1:cp_prev+1 converged without penetration —
            # the contact patch needs to change by more than one collocation
            # point to remain valid. Per the reference algorithm, this is always
            # a time-step-size problem, never a reason to widen the cp search:
            # halve dt and retry so cp can only ever advance one point per step.
            dt /= 2
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
