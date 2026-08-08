# Conservation laws, tested on the trajectory the solver actually produces.
#
# Everything else in the suite is structural -- symmetry, rank, limits, agreement with
# Reid's coefficients. Those check that the operators are assembled correctly. None of
# them checks that INTEGRATING those operators conserves what the physics conserves,
# and that is a different failure mode: an assembly can be exactly right and a stepper
# still leak energy, gain momentum, or grow the drop.
#
# The statements below are all consequences of the model rather than of any
# discretisation, so each has a value it must approach as the step shrinks, and each
# fails loudly rather than drifting quietly.

using Test
using DropSolver
using LinearAlgebra

## the solver keeps these internal; the audits need them
basis_of(p) = DropSolver.basis(p)
force_column_of(p, j) = DropSolver.force_column(p, j)

const REF = (We = 1.0, Bo = 0.0189, Oh = 0.303767)

@testset "conservation" begin

    @testset "the energy budget closes: dE/dt = -xi-dot' C xi-dot" begin
        # THE test for a variational solver, because it checks M, C, G and the
        # integrator jointly and none of them in isolation. Dotting the equations of
        # motion with the velocity,
        #
        #     M xi-ddot + C xi-dot + G xi = 0
        #  => d/dt (1/2 xi-dot' M xi-dot + 1/2 xi' G xi) = - xi-dot' C xi-dot ,
        #
        # so the mechanical energy of a FREE drop falls at exactly the rate the
        # dissipation form gives. Written this way the identity is free of the factor
        # of two that has bitten this model twice -- it is derived from the equations
        # the solver integrates, not from the Lagrangian it was derived from.
        b = ModalBasis(2:12, 2)
        F = assemble_newtonian(b, 0.05)
        N = ndof(b)
        a = zeros(N); a[1] = 0.04; a[3] = -0.02      # a two-mode initial shape
        ad = zeros(N)

        E(a, ad) = 0.5 * dot(ad, F.M, ad) + 0.5 * dot(a, F.G, a)
        diss(ad) = dot(ad, F.C, ad)

        # BDF2 with a fixed step, no contact
        dt = 0.004
        prev_a, prev_ad = copy(a), copy(ad)
        curr_a, curr_ad = copy(a), copy(ad)
        E0 = E(curr_a, curr_ad)
        integ = 0.0                                  # running integral of dissipation
        first = true
        nsteps = 500
        for _ in 1:nsteps
            c0, c1, c2 = first ? (1.0, -1.0, 0.0) : (1.5, -2.0, 0.5)
            β = c0 / dt
            hv_a  = (c1 * curr_a  + c2 * prev_a)  / dt
            hv_ad = (c1 * curr_ad + c2 * prev_ad) / dt
            A = β^2 * F.M + β * F.C + F.G
            rhs = -F.M * (β * hv_a + hv_ad) - F.C * hv_a
            a_new = A \ rhs
            ad_new = β * a_new + hv_a
            ## trapezoid on the dissipation rate, matching the integrator's own order
            integ += 0.5 * dt * (diss(curr_ad) + diss(ad_new))
            prev_a, prev_ad = curr_a, curr_ad
            curr_a, curr_ad = a_new, ad_new
            first = false
        end
        E1 = E(curr_a, curr_ad)

        @test E1 < E0                                   # it really did dissipate
        @test integ > 0.1 * (E0 - E1)                   # and the budget is not vacuous
        ## the energy lost equals the dissipation integrated, to the integrator's order
        @test isapprox(E0 - E1, integ; rtol = 2e-3)

        ## and it must CLOSE BETTER as the step shrinks -- a fixed offset would mean a
        ## conservation error rather than a truncation error
        function budget_gap(dt, nsteps)
            pa, pad = copy(a), copy(ad); ca, cad = copy(a), copy(ad)
            e0 = E(ca, cad); ig = 0.0; fst = true
            for _ in 1:nsteps
                c0, c1, c2 = fst ? (1.0, -1.0, 0.0) : (1.5, -2.0, 0.5)
                β = c0 / dt
                hva  = (c1 * ca  + c2 * pa)  / dt
                hvad = (c1 * cad + c2 * pad) / dt
                A = β^2 * F.M + β * F.C + F.G
                an = A \ (-F.M * (β * hva + hvad) - F.C * hva)
                adn = β * an + hva
                ig += 0.5 * dt * (diss(cad) + diss(adn))
                pa, pad = ca, cad; ca, cad = an, adn; fst = false
            end
            abs((e0 - E(ca, cad)) - ig) / e0
        end
        g_coarse = budget_gap(0.008, 250)
        g_fine   = budget_gap(0.002, 1000)
        @test g_fine < g_coarse / 2                     # converging, not offset
    end

    @testset "the drop never returns more than it arrived with" begin
        # CoR <= 1 is the second law applied to the whole event: the substrate is
        # rigid and does no net positive work on the drop, so the rebound speed cannot
        # exceed the impact speed. This is the cheapest possible guard and the one this
        # solver most needed -- an earlier contact closure returned CoR = 12.3, and
        # nothing in the suite objected.
        for (We, Oh) in ((0.05, 0.05), (0.3, 0.05), (1.0, 0.1),
                         (0.3, 0.3), (1.0, 0.3), (2.0, 0.5))
            r = simulate(ImpactParams(We = We, Bo = 0.0189, Oh = Oh,
                                      M = 30, K = 3, t_max = 25.0))
            @test isfinite(r.cor)
            @test 0 < r.cor <= 1.0
            ## and the drop must actually have hit something
            @test maximum(r.cp) > 0
        end
    end

    @testset "vertical impulse balances the momentum change" begin
        # The centre-of-mass equation integrated once. Over the whole event,
        #
        #     m (v_out - v_in) = integral of F dt  -  m Bo (t_out - t_in),
        #
        # with m = 4 pi/3 and F = -(4 pi/3) p_c,1. Dividing by m, the statement is
        # that the velocity change equals the integrated film acceleration less the
        # free-fall contribution -- so it checks the film force's MAGNITUDE, which
        # every other test here only constrains in sign.
        p = ImpactParams(; REF..., M = 30, K = 3, t_max = 25.0)
        r = simulate(p)
        i0, i1 = 1, length(r.t)
        ## trapezoid on -p_c,1, which is the film's contribution to v-dot
        imp = 0.0
        for i in i0:(i1-1)
            dt = r.t[i+1] - r.t[i]
            imp += 0.5 * dt * ((-r.pc1[i]) + (-r.pc1[i+1]))
        end
        grav = -p.Bo * (r.t[i1] - r.t[i0])
        predicted = imp + grav
        actual = r.v[i1] - r.v[i0]
        @test abs(imp) > 0.5                            # the film did real work
        @test isapprox(actual, predicted; rtol = 5e-3)
        ## the film's net impulse must be upward: it can push, not pull, on balance
        @test imp > 0
    end

    @testset "the drop does not change volume" begin
        # Volume is conserved exactly by the physics, and the model conserves it at
        # first order by construction: the l = 0 mode is absent, so no combination of
        # shape amplitudes changes the volume linearly. What is NOT guaranteed is that
        # the second-order term stays small, and it is the quantity to watch, because a
        # slow volume drift is exactly how a shape expansion goes wrong without ever
        # producing an obviously bad number.
        #
        # For r = 1 + sum zeta_l P_l the volume is
        #     V/V_0 = 1 + 3 sum zeta_l^2/(2l+1) + O(zeta^3),
        # so the drift is bounded by the amplitudes themselves and must vanish with them.
        ## The volume is INTEGRATED from the represented surface rather than predicted
        ## from the amplitudes. Predicting it would only restate the second-order formula
        ## and compare it with itself -- which an earlier version of this test did, and
        ## which measures nothing. What has to be checked is that the integrated volume
        ## contains NO FIRST-ORDER error, because a first-order error is an l = 0
        ## component, and the only way one can appear is by leaking in from the contact
        ## constraint or the assembly.
        gn, gw = DropSolver.gauss_legendre_nodes(200, -1.0, 1.0)
        function volume_error(We)
            p = ImpactParams(We = We, Bo = REF.Bo, Oh = REF.Oh, M = 30, K = 3,
                             t_max = 25.0)
            r = simulate(p); ls = p.ls
            wres = 0.0; wn = 0.0
            for a in r.a
                z = surface_amplitudes(p, a)
                ## V = (2 pi/3) integral (1 + zeta)^3 dmu, exactly
                vol = (2pi / 3) * sum(gw[q] *
                        (1 + sum(z[i] * DropSolver.legendre_angular(ls[i], gn[q]).P
                                 for i in eachindex(ls)))^3 for q in eachindex(gn))
                rel = vol / (4pi / 3) - 1
                pred = 3 * sum(z[i]^2 / (2 * ls[i] + 1) for i in eachindex(ls))
                nrm = sqrt(sum(abs2, z))
                ## POINTWISE, both quantities at the same instant. Comparing the maximum
                ## residual with the maximum amplitude is meaningless when they occur at
                ## different times, which is what an earlier version of this test did --
                ## it made a clean cubic scaling look like an exponent of 2.2.
                ##
                ## Only where the drop is actually deformed. At t = 0 the amplitude is
                ## exactly zero, and a ratio there divides machine epsilon by 1e-27 and
                ## reports 4e11. This trap -- a near-zero denominator producing a
                ## meaningless relative error -- has now appeared five times in this
                ## project, and it is always the same fix: exclude the degenerate states
                ## rather than floor the denominator.
                if nrm > 0.05
                    wres = max(wres, abs(rel - pred) / nrm^3)
                end
                wn = max(wn, nrm)
            end
            (wres, wn)
        end

        res1, n1 = volume_error(1.0)
        res2, n2 = volume_error(0.25)
        @test n1 > 0.02 && n2 < n1                      # the drop deformed, and less so
        ## Everything the second-order formula does not explain is third order in the
        ## amplitude, with a coefficient of order one at BOTH impact energies. So the
        ## volume of the represented shape is exactly what a linear shape expansion
        ## implies, with no first-order component -- no l = 0 leak from the contact
        ## constraint or the assembly.
        @test res1 < 1.0
        @test res2 < 1.0
        ## Recorded rather than asserted away: the SECOND-order volume error is itself
        ## large at this Weber number -- the amplitude norm reaches 0.41, so the linear
        ## shape expansion misrepresents the volume by about ten per cent. That is a
        ## property of linearising in the amplitude rather than of the solver.
        ##
        ## It is NOT a hard ceiling on We, and an earlier commit message of mine claimed
        ## it was -- it attributed two failing high-Weber runs to "linearity breaking
        ## down". That was wrong, and the correction belongs here because the commit
        ## message cannot be edited. The reference implementation completes the same
        ## cases, Oh = 0.03 at We = 5 and 10, releasing cleanly with |zeta| = 0.99 and
        ## 1.22 -- inside an energy bound computed independently from the surface
        ## stiffness. Those runs failed here for a solver reason, not a modelling one.
        ## What the volume error bounds is ACCURACY at large amplitude, not admissibility.
        @test 0.05 < 0.6 * n1^2 < 0.20
    end

    @testset "the dimensionless groups are formed correctly" begin
        # The solver takes We, Bo and Oh directly, so its invariance under a change of
        # physical units is structural -- there is no dimensional input to vary, and a
        # test of it would be vacuous. What is NOT vacuous is the arithmetic one level
        # up, where physical properties become groups and a measured time in seconds
        # becomes a time in capillary units. That is where a real error lived: the
        # measured contact times are in seconds and comparing them raw is meaningless.
        #
        # Two different fluids, same groups, must therefore give the same answer.
        groups(R, sigma, rho, eta0, g, U) =
            (We = rho * R * U^2 / sigma,
             Bo = rho * g * R^2 / sigma,
             Oh = eta0 / sqrt(rho * sigma * R))

        g1 = groups(0.0003, 0.0728, 1000.0, 0.05, 9.81, 0.3)
        ## scale the length up by four and pick the other properties to hold all three
        ## groups fixed: Oh ~ eta/sqrt(rho sigma R), Bo ~ rho g R^2/sigma, We ~ rho R U^2/sigma
        R2 = 4 * 0.0003
        sigma2 = 0.0728 * 16.0            # Bo fixed at fixed rho, g: sigma ~ R^2
        eta2 = 0.05 * sqrt(16.0 * 4.0)    # Oh fixed: eta ~ sqrt(sigma R)
        U2 = 0.3 * sqrt(16.0 / 4.0)       # We fixed: U^2 ~ sigma/R
        g2 = groups(R2, sigma2, 1000.0, eta2, 9.81, U2)

        for k in (:We, :Bo, :Oh)
            @test isapprox(getfield(g1, k), getfield(g2, k); rtol = 1e-12)
        end
        ## and the capillary time really does scale as sqrt(rho R^3/sigma)
        tc1 = sqrt(1000.0 * 0.0003^3 / 0.0728)
        tc2 = sqrt(1000.0 * R2^3 / sigma2)
        @test isapprox(tc2 / tc1, sqrt(4.0^3 / 16.0); rtol = 1e-12)
    end
end

# Accuracy against measurement, as a hard ceiling rather than a report.
#
# The tests above check that the model is self-consistent. This one checks that it is
# RIGHT, against experiment, on both KPIs at once, with a bound that fails rather than
# prints. Fifteen per cent relative error on each.
#
# The points are named with their sources so the choice can be audited, and one known
# EXCEEDANCE is recorded explicitly rather than avoided by choosing only points that
# pass -- at high Ohnesorge the contact time is about twenty per cent low, which is the
# film-tension limitation the model has on record, and hiding it by omission would make
# this test a decoration.
@testset "within 15 per cent of experiment, both KPIs" begin
    TOL = 0.15

    ## Gabbard et al. (2025), water-glycerol, from julia/data/gabbard2025_*.csv,
    ## binned at the stated Oh and Bo. Newtonian, so K = 2 carries Reid's damping.
    newtonian_points = [
        (name = "Oh=0.037 We=0.296", We = 0.296,  Bo = 0.01875, Oh = 0.0373,
         cor = 0.7843, tc = 2.747),
        (name = "Oh=0.037 We=0.866", We = 0.8658, Bo = 0.01875, Oh = 0.0373,
         cor = 0.6851, tc = 2.486),
    ]
    for pt in newtonian_points
        r = simulate(ImpactParams(We = pt.We, Bo = pt.Bo, Oh = pt.Oh,
                                  M = 45, K = 2, t_max = 25.0))
        @test isfinite(r.cor)
        @test abs(r.cor - pt.cor) / pt.cor < TOL
        @test abs(r.tc  - pt.tc)  / pt.tc  < TOL
    end

    ## The 3000 ppm shear-thinning fluid, from julia/derivations/data. Its parameters
    ## come from its own Cross fit -- nothing is fitted to the impact data.
    @testset "3000 ppm shear-thinning" begin
        eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
        K_cross, m_cross = 18.48081673111359, 0.7430524574330837
        R, sigma, rho = 0.0003, 0.0728, 1000.0
        t_cap = sqrt(rho * R^3 / sigma)
        Oh_0 = eta_0 / sqrt(rho * sigma * R)
        etaf = gd -> carreau(gd; lambda_c = K_cross / t_cap, a = m_cross,
                             n = 1 - m_cross, eta_inf_ratio = eta_inf / eta_0)
        ## measured averages near these Weber numbers
        for (We, cor_e, tc_e) in ((0.0769, 0.8168, 3.099),
                                  (0.1912, 0.8037, 2.693),
                                  (0.4754, 0.7816, 2.288))
            r = simulate(ImpactParams(We = We, Bo = 0.012, Oh = Oh_0, M = 30, K = 3,
                                      eta = etaf, t_max = 25.0))
            @test isfinite(r.cor)
            @test abs(r.cor - cor_e) / cor_e < TOL
            @test abs(r.tc  - tc_e)  / tc_e  < TOL
        end
    end

    @testset "known exceedance: restitution at intermediate Ohnesorge" begin
        # Oh = 0.077, We = 0.336: CoR comes out 0.696 against 0.595 measured, seventeen
        # per cent HIGH -- the opposite sign to the fleet-wide bias, which is low. This
        # point is asserted at its measured error rather than left out of the list above,
        # because a fifteen-per-cent ceiling enforced only on the points that clear it is
        # not a ceiling. Contact time at the same point is inside the bound.
        r = simulate(ImpactParams(We = 0.3362, Bo = 0.01776, Oh = 0.0767,
                                  M = 45, K = 2, t_max = 25.0))
        err = abs(r.cor - 0.5952) / 0.5952
        @test err > TOL                                     # it does not clear the bound
        @test err < 0.22                                    # and it is this bad, not worse
        @test r.cor > 0.5952                                # high, not low
        @test abs(r.tc - 2.816) / 2.816 < TOL               # contact time is fine here
    end

    @testset "known exceedance: contact time at high Ohnesorge" begin
        # At high Ohnesorge the contact time is well outside fifteen per cent -- about
        # twenty low at Oh = 0.685 -- while the restitution is inside it. That pattern
        # is the signature of premature release: the measured contact time flattens to a
        # floor near the l = 2 period, 2 pi/sqrt(8) = 2.22, and the model goes below it.
        # Nothing in the solver holds the film in tension, so it lets go as soon as the
        # edge pressure turns. This test asserts the failure at its measured size, so
        # that FIXING it breaks the test and forces this comment to be revisited.
        r = simulate(ImpactParams(We = 1.383, Bo = 0.02714, Oh = 0.6849,
                                  M = 45, K = 2, t_max = 25.0))
        @test abs(r.cor - 0.1878) / 0.1878 < TOL          # restitution is fine
        err_tc = abs(r.tc - 2.752) / 2.752
        @test err_tc > TOL                                 # contact time is not
        @test err_tc < 0.30                                 # and it is this bad, not worse
        @test r.tc < 2.22                                   # below the l = 2 period floor
    end
end

# The two guards that came out of diagnosing a run which manufactured forty times the
# energy it arrived with. Both are cheap, both are physical, and either would have caught
# that failure the moment it appeared -- where the guards already in this file did not,
# because `CoR <= 1` only sweeps to We = 2 and the budget audit above switches contact off.
@testset "energy accounting through contact" begin
    REF2 = (We = 1.0, Bo = 0.0189, Oh = 0.303767)

    @testset "stored surface energy cannot exceed the energy that arrived" begin
        # A CEILING, and the cheapest possible statement of the second law for this
        # problem: the drop cannot deform beyond what its own kinetic energy pays for.
        # The softest retained mode is the cheapest place to put energy, so
        #
        #     max|zeta_l|  <=  sqrt( E_kin / min_l (V per unit zeta_l^2) )
        #
        # This is not a tolerance to be tuned -- exceeding it is impossible, so any
        # exceedance is a defect by construction. A diagnosis that took most of a day
        # would have been a one-line failure with this in place: the run in question
        # reached |zeta| = 9.7 against a ceiling of 1.44.
        for (We, Oh, M, K) in ((1.0, 0.303767, 45, 2), (0.3, 0.05, 30, 2), (2.0, 0.5, 30, 2))
            p = ImpactParams(We = We, Bo = REF2.Bo, Oh = Oh, M = M, K = K, t_max = 25.0)
            r = simulate(p)
            F = assemble_newtonian(basis_of(p), Oh)
            Ekin = 0.5 * (4pi/3) * We
            ## surface energy of a unit amplitude in each mode, taken one at a time
            N = length(r.a[1])
            cheapest = minimum(begin
                e = zeros(N); e[i] = 1.0
                0.5 * dot(e, F.G, e)
            end for i in 1:N if begin e = zeros(N); e[i] = 1.0; dot(e, F.G, e) end > 0)
            ceiling = sqrt(Ekin / cheapest)
            observed = maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a)
            @test observed > 0.01                 # the drop really deformed
            @test observed <= ceiling             # and not beyond what it could pay for
        end
    end

    @testset "the film's work accounts for the energy change" begin
        # The audit that distinguishes "a pulling film did work" from "energy appeared
        # from outside the equations" -- and it must be run WITH contact, because that is
        # the only regime where the film contributes at all. Validated at 2.4 per cent on
        # this case; a failing run missed by a factor of 156, which is what told me the
        # accounting, not the film, was the thing to look at.
        p = ImpactParams(; REF2..., M = 30, K = 3, t_max = 25.0)
        r = simulate(p)
        b = basis_of(p); F = assemble_newtonian(b, p.Oh); mass = 4pi/3
        E(a, ad, z, v) = 0.5*dot(ad, F.M, ad) + 0.5*dot(a, F.G, a) +
                         0.5*mass*v^2 + mass*p.Bo*z
        function film_power(i)
            Q = zeros(length(r.a[i]))
            for j in eachindex(r.pc[i])
                Q .+= r.pc[i][j] .* force_column_of(p, j)
            end
            dot(Q, r.adot[i]) - mass * r.pc[i][2] * r.v[i]
        end
        Wf = 0.0; Dis = 0.0
        for i in 1:(length(r.t)-1)
            h = r.t[i+1] - r.t[i]
            Wf  += 0.5*h*(film_power(i) + film_power(i+1))
            Dis += 0.5*h*(dot(r.adot[i], F.C, r.adot[i]) +
                          dot(r.adot[i+1], F.C, r.adot[i+1]))
        end
        dE = E(r.a[end], r.adot[end], r.z[end], r.v[end]) -
             E(r.a[1],   r.adot[1],   r.z[1],   r.v[1])
        @test Dis > 0                                    # it dissipated
        @test abs(dE) > 0.1                              # and the audit is not vacuous
        @test isapprox(dE, Wf - Dis; rtol = 0.10)        # measured: 2.4 per cent
    end
end
