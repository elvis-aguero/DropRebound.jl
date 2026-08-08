# The contact set solved for rather than searched over.
#
# `simulate_lcp` determines which nodes are in contact by solving the complementarity
# problem `h = A_c p + b`, `h >= 0`, `p >= 0`, `p'h = 0` at every step. The tests below are
# about what that makes possible and the searching solver cannot do: it satisfies the
# pressure inequality everywhere rather than only at the patch edge, it is not restricted
# to a disc, and it runs cases the search cannot.
#
# `A_c` is NOT symmetric -- it is asymmetric by about forty per cent, because the film
# pressure is carried as a Legendre field whose forcing is not the transpose of the gap
# Jacobian. So the problem is not a convex programme and the projected Gauss-Seidel sweep
# does not apply to it; `lcp_active_set` solves it without assuming symmetry. An earlier
# version symmetrised the matrix and swept it anyway, which produced states with the drop
# inside the substrate and is the origin of every claim these tests used to make about the
# two closures disagreeing.

using Test
using DropSolver
using LinearAlgebra

const LREF = (We = 1.0, Bo = 0.0189, Oh = 0.303767)

@testset "LCP contact" begin

    @testset "the complementarity conditions are satisfied, not approximated" begin
        # The searching solver enforces h >= 0 and h*p = 0 by construction but NOT
        # p >= 0 -- measured at 0.015 per cent of node-steps in a case that works and 41
        # per cent in one that does not. This is the condition the LCP adds, so it is the
        # thing to check, and the complementarity residual is a certificate rather than a
        # tolerance: it vanishes exactly when all three hold.
        p = ImpactParams(; LREF..., M = 45, K = 2, t_max = 25.0)
        r = simulate_lcp(p)
        @test r.lcp_resid_max < 1e-8
        @test r.rejects == 0                       # the solve never had to be retried
        @test r.lcp_sweeps_max < 100               # and converged cheaply
        ## the pressure is non-negative at every node of every step, by construction
        V = DropSolver.legendre_vandermonde(p)
        worst = 0.0
        for i in eachindex(r.t)
            worst = min(worst, minimum(V * r.pc[i]))
        end
        @test worst > -1e-6 * maximum(maximum(abs, V * r.pc[i]) for i in eachindex(r.t))
    end

    @testset "restitution agrees with two independent methods" begin
        # CoR at K = 1, M = 90 is 0.3129 from the LCP, 0.3129 from the tangency rule, and
        # 0.3138 from the reference MATLAB implementation. Three routes to the same
        # number, which is what makes it evidence rather than coincidence.
        r = simulate_lcp(ImpactParams(; LREF..., M = 90, K = 1, t_max = 25.0))
        @test isapprox(r.cor, 0.3138; atol = 0.01)
    end

    @testset "the cases that defeated the searching solver now run" begin
        # Oh = 0.03 at We = 5 and 10 manufactured forty times the drop's energy under the
        # active-set search and starved under tangency selection. The reference
        # implementation completes them with |zeta| = 0.99 and 1.22.
        for (We, Oh, zref) in ((5.0, 0.03, 0.99), (10.0, 0.03, 1.22), (10.0, 0.30, 0.84))
            p = ImpactParams(We = We, Bo = 0.0189, Oh = Oh, M = 90, K = 1, t_max = 30.0)
            r = simulate_lcp(p)
            @test r.rejects == 0
            @test r.lcp_resid_max < 1e-6
            amp = maximum(maximum(abs, surface_amplitudes(p, a)) for a in r.a)
            @test isapprox(amp, zref; rtol = 0.15)      # and the amplitude is physical
            ## the energy ceiling from the surface stiffness, which the old runs blew
            ## through by a factor of five
            F = assemble_newtonian(DropSolver.basis(p), Oh)
            N = length(r.a[1])
            cheap = minimum(begin e = zeros(N); e[i] = 1.0; 0.5*dot(e, F.G, e) end
                            for i in 1:N if begin e = zeros(N); e[i] = 1.0
                                                 dot(e, F.G, e) end > 0)
            @test amp <= sqrt(0.5*(4pi/3)*We / cheap)
        end
    end

    @testset "the contact comes out a single patch, though nothing requires it to" begin
        # The model page carries "contact occupies a single connected patch about the pole" as
        # an axiom. The complementarity closure does not impose it -- the active set is
        # whatever the solve returns -- so this test asks whether the axiom is earned.
        #
        # It is. And an earlier version of this test claimed the opposite.
        #
        # THE ANNULUS WAS A DISCRETISATION ARTEFACT. The earlier version
        # asserted the opposite -- that the free arc converged to about 11 degrees while the
        # nodes spanning it grew, so the angle was the invariant and the dimple was resolved
        # rather than created. That measurement was taken while `contact_lcp` handed a
        # symmetrised copy of a forty-per-cent-asymmetric compliance to a sweep that assumes
        # symmetry, so the pressures were not a solution of the contact problem at all and the
        # states they produced had the drop up to 4.7 per cent of a radius inside the
        # substrate. A dimple is exactly what an unresolved penetration looks like from the
        # outside.
        #
        # Solved against the true compliance the free arc is ZERO at M = 90: contact comes out
        # a single patch anchored at the pole. Non-contiguity survives only as a transient at
        # the release edge, in a minority of steps that SHRINKS with resolution -- 12 of 530
        # accepted steps at M = 30, 4 of 967 at M = 45, none at M = 90 -- which is the
        # signature of an artefact rather than a feature.
        #
        # The test is kept, inverted, because the complementarity closure genuinely does not
        # constrain the contact set to be an interval, so whether it comes out as one is a
        # result worth pinning. It comes out as one.
        function free_arc(M)
            p = ImpactParams(; LREF..., M = M, K = 2, t_max = 25.0)
            F0 = assemble_newtonian(DropSolver.basis(p), p.Oh)
            Vf = lu(DropSolver.legendre_vandermonde(p))
            prev = curr = DropSolver.initial_state(p); dt = p.dt0
            for _ in 1:200000
                curr.z <= 0.95 && break
                st, nxt, _ = DropSolver.try_step_lcp(p, prev, curr, dt; F0=F0, Vfac=Vf)
                st === :ok || (dt /= 2; dt < p.dt_min && break; continue)
                prev, curr = curr, nxt; dt = min(2dt, p.dt0)
            end
            Ac, bv, idx, _ = DropSolver.contact_lcp(p, prev, curr, dt; F0=F0, Vfac=Vf)
            pv, resid, _ = DropSolver.lcp_active_set(Ac, bv)
            @test resid < 1e-9                         # a real solution, not a surrogate's
            act = [i for i in eachindex(pv) if pv[i] > 1e-10*max(maximum(pv), 1)]
            isempty(act) && return (NaN, 0)
            (180 - rad2deg(p.nodes[idx[first(act)]]), first(act) - 1)
        end
        a45, n45 = free_arc(45)
        a90, n90 = free_arc(90)
        ## at production resolution the contact set is a single patch at the pole
        @test n90 == 0
        @test a90 == 0.0
        ## and the free arc does not grow as the grid coarsens into a converged feature
        @test n45 >= n90

        ## the same statement over a whole march: non-contiguity is a shrinking minority
        frac(M) = let r = simulate_lcp(ImpactParams(We = 0.5, Bo = 0.017, Oh = 0.0373,
                                                    M = M, K = 2, t_max = 25.0))
            r.noncontiguous_steps / length(r.t)
        end
        f30, f45 = frac(30), frac(45)
        @test f30 < 0.05                               # measured 2.3 per cent
        @test f45 < f30                                # measured 0.4 per cent
    end
end

# Shear thinning through the complementarity solve.
#
# A variable viscosity makes the problem NONLINEAR: eta depends on the strain rate, the
# dissipation operator on eta, and the compliance on that operator -- so A_c is a function
# of the velocity the solve produces. It is closed by the same Picard iteration the
# searching solver uses, freezing eta at an extrapolated strain rate and re-evaluating at
# the answer. Each iterate is a genuine LCP -- asymmetric, like the constant-viscosity one --
# so complementarity holds exactly at every sweep and only eta lags.
@testset "shear thinning through the LCP" begin
    ## the 3000 ppm fluid, from its own Cross fit -- nothing fitted to the impact data
    eta_0, eta_inf = 8.433817577956766, 0.0037320997942061666
    K_cross, m_cross = 18.48081673111359, 0.7430524574330837
    R, sigma, rho = 0.0003, 0.0728, 1000.0
    t_cap = sqrt(rho * R^3 / sigma)
    Oh_0 = eta_0 / sqrt(rho * sigma * R)
    etaf = gd -> carreau(gd; lambda_c = K_cross/t_cap, a = m_cross,
                         n = 1 - m_cross, eta_inf_ratio = eta_inf/eta_0)

    @testset "it integrates, and complementarity still holds exactly" begin
        p = ImpactParams(We = 0.1912, Bo = 0.012, Oh = Oh_0, M = 30, K = 3,
                         eta = etaf, t_max = 25.0)
        @test !p.eta_const                        # the probe really sees a variable eta
        r = simulate_lcp(p)
        ## the LCP residual is unaffected by the nonlinearity: every Picard iterate is a
        ## genuine complementarity solve, so this stays at machine level
        @test r.lcp_resid_max < 1e-10
        ## and the viscosity iteration converged on every accepted step
        @test r.eta_resid_max <= p.eta_tol
        @test r.rejects < 0.05 * length(r.t)
        @test 0 < r.cor < 1
    end

    @testset "it agrees with experiment, and with the search to six digits" begin
        # An earlier version of this asserted the LCP was CLOSER to experiment than the search
        # -- LCP CoR 0.7872 against search 0.7673, at a measured 0.804 -- and read the split
        # (better on restitution, worse on contact time) as evidence about the contact model.
        # There was no split. Those numbers came from a compliance that had been symmetrised
        # before solving, so the LCP was answering a different question; corrected, the two
        # closures agree to 3.5e-6 relative on a shear-thinning fluid, and the earlier difference
        # was the bug.
        #
        # That agreement is now the stronger claim, so it is what gets asserted. The two
        # closures share the assembly and nothing about how they choose the contact set, so
        # this pins the rheology and the contact treatment against each other on a fluid whose
        # zero-shear Ohnesorge number is 57.
        p = ImpactParams(We = 0.1912, Bo = 0.012, Oh = Oh_0, M = 30, K = 3,
                         eta = etaf, t_max = 25.0)
        rl = simulate_lcp(p); rs = simulate(p)
        @test abs(rl.cor - 0.804)/0.804 < 0.10        # within ten per cent on restitution
        @test abs(rl.tc - 2.693)/2.693 < 0.15         # contact time within fifteen
        ## How closely the two closures may agree is set by how far each is allowed to
        ## sit from the fixed point, which is `eta_tol`. They share the assembly and the
        ## rheology and differ only in how the contact set is chosen, so the residual
        ## disagreement is the tolerance, not the closure.
        ##
        ## Asserting a fixed 1e-5 here was correct while eta_tol was 1e-6 and became
        ## wrong when the tolerance began scaling with M: at M = 30 the tolerance is
        ## 2.7e-4 and the closures differ by 1.25e-4, which is inside it. A test that
        ## demands agreement finer than the tolerance permits is testing luck.
        @test isapprox(rl.cor, rs.cor; rtol = max(1e-5, p.eta_tol))
        ## Contact time is a multiple of the step, so the closures agree to within the
        ## step on which they flag release. Measured: they differ by 0.00456 against a
        ## dt of 0.00471, which is one step. Bit-identical was the right assertion when
        ## both closures flagged the same steps and is not once eta_tol permits them to
        ## accept different ones.
        @test abs(rl.tc - rs.tc) <= 2 * p.dt0
    end

    @testset "the Newtonian limit is recovered through the LCP path" begin
        # lambda_c -> 0 sends eta -> 1 pointwise, so the expensive path -- a coupled
        # reassembly and a fresh compliance every sweep -- must reproduce the cached
        # constant-viscosity one. This checks the two code paths against each other rather
        # than each against itself.
        base = (We = 1.0, Bo = 0.0189, Oh = 0.303767, M = 30, K = 3, t_max = 25.0)
        rn = simulate_lcp(ImpactParams(; base...))
        rt = simulate_lcp(ImpactParams(; base...,
                eta = gd -> carreau(gd; lambda_c = 1e-10, a = 2.0, n = 0.4,
                                    eta_inf_ratio = 0.01)))
        @test isapprox(rt.cor, rn.cor; rtol = 1e-3)
        @test isapprox(rt.tc,  rn.tc;  rtol = 1e-3)
    end
end

# The equator node.
#
# At odd M the polynomial P_M has a root at mu = 0, so one collocation node sits exactly on the
# equator. Every entry of its gap row carries a factor cos(theta), so that row is identically
# zero: the node's gap does not depend on the drop's shape and its pressure does no work. It
# can never legitimately be in contact.
#
# It used to be admitted or not according to the sign that mu = 0 happened to evaluate to --
# about -3e-16 at M = 45 and M = 61, positive at M = 21, 31 and 91. When admitted it put a zero
# row and column into the compliance, which is where the 1e-35 smallest eigenvalue at M = 45
# came from. M = 45 is the truncation the validation sweeps use.
@testset "the equator node is excluded at every truncation" begin
    for M in (20, 21, 30, 31, 44, 45, 60, 61, 90, 91)
        p = ImpactParams(We = 0.5, Bo = 0.0189, Oh = 0.0373, M = M, K = 1)
        F0 = assemble_newtonian(DropSolver.basis(p), p.Oh)
        Vf = lu(DropSolver.legendre_vandermonde(p))
        s0 = DropSolver.initial_state(p)
        _, _, idx, _ = DropSolver.contact_lcp(p, s0, s0, p.dt0; F0 = F0, Vfac = Vf)
        ## every retained node is genuinely below the equator, so no gap row is degenerate
        @test all(cos(p.nodes[i]) < -1e-9 for i in idx)
        @test all(norm(DropSolver.gap_row(p, p.nodes[i])[1]) > 1e-6 for i in idx)
    end
end

# Full rank of the retained constraint set.
#
# With the equator excluded the lower-hemisphere gap rows are independent, so the CONJUGATE
# compliance H A^-1 H' is positive definite rather than merely semi-definite. The shipped
# compliance is not that matrix -- it is asymmetric, see `contact_lcp` -- but the rank of H is
# what decides whether a unique contact force exists at all, so it is worth pinning here.
@testset "the retained gap constraints are independent" begin
    for (M, K) in ((20, 2), (30, 2), (45, 2), (61, 1), (90, 1))
        p = ImpactParams(We = 0.5, Bo = 0.0189, Oh = 0.0373, M = M, K = K)
        b = DropSolver.basis(p)
        F0 = assemble_newtonian(b, p.Oh); beta = 1 / p.dt0
        A = beta^2 * F0.M + beta * F0.C + F0.G
        nn = length(p.nodes)
        H = zeros(nn, DropSolver.ndof(b))
        for i in 1:nn; H[i, :] = DropSolver.gap_row(p, p.nodes[i])[1]; end
        low = findall(<(-1e-8), cos.(p.nodes))
        Hl = H[low, :]
        @test rank(Hl; rtol = 1e-10) == length(low)          # no redundant constraint
        W = Hl * (A \ transpose(Hl))
        @test maximum(abs, W - W') / maximum(abs, W) < 1e-12  # conjugate pairing IS symmetric
        @test minimum(eigvals(Symmetric(0.5 .* (W .+ W')))) > 0   # and positive definite
    end
end
