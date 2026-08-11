using Test
using DropSolver

@testset "Timestepper: contact-point continuity" begin
    # Gabbard et al. (2025, JFM, "Drop rebound at low Weber number") §C.1:
    # "if the model predicts that the number of contact mesh points m^k
    # changes by more than one in a single step, the step is rejected,
    # halved, and candidate solutions are recalculated." The discrete
    # contact search only ever tests cp_prev-1, cp_prev, cp_prev+1 per step
    # while already in contact, for exactly this reason: those three
    # candidates are the only ones for which the tangency-error scoring
    # (contact_error) is reliable, since the surrounding two points
    # (cp_prev-2, cp_prev+2) are what would confirm the accepted candidate
    # is a genuine local minimum rather than an edge artifact.
    #
    # The one exception is contact *onset* (cp_prev == 0): a smooth convex
    # body touching a flat plane nucleates a contact patch whose radius
    # initially grows like √(t − t_onset) — arbitrarily fast in continuum
    # time — so no achievable dt resolves the first instants of contact one
    # collocation point at a time. This is a genuine feature of rigid-contact
    # kinematics, not a discretization artifact, so cp is allowed to jump by
    # more than one *only* when the previous state was airborne (cp == 0).
    #
    # A restart heuristic previously allowed cp to jump by more than one
    # collocation point within a single accepted step even *during* contact
    # ("vigorous impacts"), which is exactly the mechanism behind spurious
    # energy injection: the solver accepted whichever widened-window
    # candidate scored lowest, with no guarantee that candidate was a true
    # local minimum rather than an artifact of the wrong branch of the
    # (by-then-unreliable) linearized height function. This case (M=15,
    # Oh=0.05, Bo=0.02, We=0.5, OB) is a concrete, fast-running instance
    # where that mechanism used to fire mid-contact.
    M = 15; Oh = 0.05; Bo = 0.02; We = 0.5
    dt_max    = make_dt_max(M)
    theta_vec = make_theta_vec(M)
    precomp   = precompute_integrals(NaN, M)[1]
    cfg       = SimConstants(M, M+1, Oh, Bo, theta_vec, precomp, dt_max)

    init = DropState(M)
    init.z = 1.05; init.v = -sqrt(We); init.dt = dt_max; init.cp = 0

    times, states = solve_drop!(cfg, OBParams(0.5, 0.5), deepcopy(init);
                                 t_end=2.0, save_every=1e-9)

    @test all(1:length(states)-1) do i
        cp_from, cp_to = states[i].cp, states[i+1].cp
        abs(cp_to - cp_from) <= 1 || cp_from == 0
    end
end
