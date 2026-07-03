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
  6. After ≥4 halvings with no valid step, restart at the original step dt with
     cp search expanded to cp_prev+2 (one extra level). This handles cases where
     the contact patch spans more collocation points than cp_prev+1 predicts — a
     physical need, not a time-step-size issue.
  7. Ramp dt back toward dt_max at +10% per step
"""
function solve_drop!(cfg::SimConstants, ob::OBParams, init::DropState;
                     st::STParams        = STParams(),
                     cl::CLParams        = CLParams(),
                     t_end::Float64      = 10.0,
                     save_every::Float64 = 0.1,
                     dt_init::Float64    = cfg.dt_max,
                     dt_min::Float64     = cfg.dt_max * 1e-4)

    is_ob = ob.De1 > 0.0 && ob.beta_s < 1.0
    is_st = !is_ob && st.eps_ST > 0.0 && !isempty(st.Gamma)
    # Contact-line physics changes ONLY the contact-set selection criterion (not M, b,
    # or the COM equation). See docs/DropRebound_ContactLine.tex §Discrete.
    #   is_cl_mob (ξ>0): Milestone-2 operator-split mobility residual.
    #   is_cl     (θ_e<π): Milestone-1 quasi-static contact-angle residual.
    # (θ_e=π, ξ=0) ⇒ neither ⇒ identical to GA.
    is_cl_mob = cl.xi > 0.0
    is_cl     = cl.theta_e < π - 1e-12 || is_cl_mob
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

    # Per-step adaptivity state (reset after each accepted step)
    dt_step_start  = dt  # dt at the start of the current time step attempt
    n_halvings     = 0   # halvings within the current step (or since last restart)
    cp_upper_extra = 0   # additional cp levels granted by restart cycles (0 = normal)

    # Milestone 2: continuous contact-line radius tracker (mobility law).
    r_c_track = contact_edge_radius(init, cfg.theta_vec)

    while t < t_end
        order   = min(length(history), 2)
        cp_prev = history[end].cp

        best_state = nothing
        best_err   = Inf
        cand_by_cp = Dict{Int,DropState}()   # valid candidates this step (mobility tracker)

        # cp search upper bound: normally cp_prev+1 (continuity);
        # each restart cycle grants one extra level up to N_angles-1
        cp_upper = min(cfg.N_angles - 1, cp_prev + 1 + cp_upper_extra)

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
                # Quasi-static selection residual (Milestone 1); reduces to GA tangency
                # at θ_e=π. The mobility case (ξ>0) reuses this to locate the
                # quasi-static target, then rate-limits it below.
                err = is_cl ?
                    contact_angle_error(candidate, cfg.theta_vec, cp, cl.theta_e) :
                    abs(contact_error(candidate, cfg.theta_vec, cp))
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
                # Milestone 2: keep every valid candidate for the mobility tracker
                if is_cl_mob && isfinite(err)
                    cand_by_cp[cp] = candidate
                end
            end
        end

        # Milestone 2: rate-limited contact-line mobility (continuous tracker).
        # The quasi-static pick (best_state) gives the target radius r_qs; the physical
        # contact line moves toward it at the finite rate set by ξ, decoupled from the
        # integer mesh. Limits: ξ→0 reaches r_qs (=Milestone 1); ξ→∞ stays put (pinned).
        if is_cl_mob && best_state !== nothing && !isempty(cand_by_cp)
            r_qs = contact_edge_radius(best_state, cfg.theta_vec)
            # apparent angle at the currently tracked line (candidate nearest r_c_track)
            cur_cp = argmin(cp -> abs(contact_edge_radius(cand_by_cp[cp], cfg.theta_vec)
                                      - r_c_track), collect(keys(cand_by_cp)))
            θ_d_now = cur_cp == 0 ? cl.theta_e :
                contact_angle(cand_by_cp[cur_cp],
                    (cfg.theta_vec[cur_cp] +
                     cfg.theta_vec[min(cur_cp+1, length(cfg.theta_vec))]) / 2)
            mob     = (dt / cl.xi) * (cos(cl.theta_e) - cos(θ_d_now))  # signed rate·dt
            desired = r_qs - r_c_track
            step    = sign(mob) == sign(desired) ?
                      copysign(min(abs(mob), abs(desired)), desired) : 0.0
            r_max   = maximum(contact_edge_radius(cand_by_cp[cp], cfg.theta_vec)
                              for cp in keys(cand_by_cp))
            r_c_new = clamp(r_c_track + step, 0.0, r_max)
            sel_cp  = argmin(cp -> abs(contact_edge_radius(cand_by_cp[cp], cfg.theta_vec)
                                       - r_c_new), collect(keys(cand_by_cp)))
            best_state = cand_by_cp[sel_cp]
            r_c_track  = r_c_new
        end

        if best_state === nothing || best_err == Inf
            n_halvings += 1
            @debug "Step rejected" t=t dt=dt n_halvings=n_halvings cp_prev=cp_prev cp_upper_extra=cp_upper_extra

            # Principled cp restart: after ≥2 halvings without a valid step,
            # reset to the original step dt and allow one more cp level. This
            # handles vigorous impacts where the contact patch jumps more than
            # one collocation point at once — a physical gap, not a time-step
            # problem. Repeat up to 4 times (cp_prev+5 max) before halving.
            # Threshold = 2 so the restart fires before dt can reach dt_min.
            if n_halvings >= 2 && cp_upper_extra < 4 &&
               cp_prev + 1 + cp_upper_extra < cfg.N_angles - 1
                dt             = dt_step_start
                cp_upper_extra += 1
                n_halvings     = 0
                @debug "cp restart: dt reset, cp_upper = cp_prev+$(1 + cp_upper_extra)" t=t dt=dt
                clear_jac_cache!()
                continue
            end

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

        # Ramp dt back toward dt_max and reset per-step adaptivity state
        dt             = min(dt * 1.1, cfg.dt_max)
        dt_step_start  = dt
        n_halvings     = 0
        cp_upper_extra = 0

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
