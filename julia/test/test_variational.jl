# The variational assembly, against the Newtonian theory it must reproduce.
#
# This is the first thing the solver has to get right: with a constant viscosity,
# the variational statement of the model -- three quadratic forms and the
# Euler-Lagrange equations, with the interior retained as part of the state -- must
# return Reid's exact coefficients. `reid_lambda_omega2` computes them by an
# entirely separate route (Newton on the characteristic equation, via a Bessel
# ratio), so agreement is evidence about the assembly and not about itself.
#
# It also ties the one loose end the derivation page could not: that the surface
# equation obtained by varying the surface energy agrees with the traction balance.
# If it did not, these numbers would not match.

using Test
using DropSolver
using LinearAlgebra

@testset "variational assembly" begin

    @testset "the forms are Hessians, so symmetric" begin
        for (l, K, Oh) in ((2, 4, 0.05), (5, 3, 0.3), (8, 5, 0.006))
            F = assemble(RitzBasis(l, K), Oh)
            @test F.M ≈ F.M'
            @test F.C ≈ F.C'
            @test isposdef(Symmetric(F.M))          # a kinetic energy
            @test all(eigvals(Symmetric(F.C)) .> -1e-12)  # a dissipation rate
        end
    end

    @testset "the surface energy is rank one" begin
        # It reaches the interior only through the boundary trace, so it cannot
        # resolve interior structure. A higher rank would mean the stiffness had
        # acquired a dependence on the profile shape, which surface tension has no
        # way to supply.
        for K in 2:6
            @test rank(assemble(RitzBasis(3, K), 0.1).G) == 1
        end
    end

    @testset "one trial function is the inviscid limit" begin
        # phi_1 = x^(l+1) is potential flow, so K = 1 must give Rayleigh's
        # frequency exactly and Lamb's damping exactly -- the two classical results
        # the variational route reproduces without a Bessel function in sight.
        for l in 2:6
            F = assemble(RitzBasis(l, 1), 1.0)
            ω² = F.G[1, 1] / F.M[1, 1]
            @test ω² ≈ l * (l - 1) * (l + 2)  rtol=1e-10
            λ_per_Oh = F.C[1, 1] / (2 * F.M[1, 1])
            @test λ_per_Oh ≈ (l - 1) * (2l + 1)  rtol=1e-10
        end
    end

    @testset "it reproduces Reid at finite Ohnesorge" begin
        # Across three decades of Oh and both dynamical regimes. At Oh = 1, l = 2 the
        # root pair has gone REAL (above the critical Ohnesorge), which is why the
        # coefficients are Vieta on the two slowest roots rather than on a conjugate
        # pair -- reading them as a conjugate pair is wrong by 64% there.
        for (Oh, l, K, tol) in ((0.006, 2, 8, 1e-5),
                                (0.05,  2, 8, 1e-10),
                                (0.3,   2, 8, 1e-10),
                                (1.0,   2, 8, 1e-10),   # overdamped: real roots
                                (0.05,  4, 8, 1e-6),
                                (0.05,  8, 8, 1e-5),
                                (0.3,   6, 8, 1e-7))
            λ, ω² = dominant_pair(RitzBasis(l, K), Oh)
            rλ, rω², _ = reid_lambda_omega2(Oh, l)
            @test λ  ≈ rλ   rtol=tol
            @test ω² ≈ rω²  rtol=tol
        end
    end

    @testset "refining the basis converges" begin
        # Quadratic convergence is the property that makes a handful of trial
        # functions worth as much as a fine radial grid. Check the error actually
        # falls rather than just being small at one K.
        Oh, l = 0.05, 2
        rλ, _, _ = reid_lambda_omega2(Oh, l)
        errs = [abs(dominant_pair(RitzBasis(l, K), Oh)[1] - rλ) / rλ for K in 2:6]
        @test errs[end] < errs[1] / 1e4
        @test issorted(errs; rev=true) || errs[end] < 1e-9
    end

    @testset "the basis has a conditioning limit, and it is visible" begin
        # The monomial trial functions are Vandermonde-like: cond(M) grows by roughly
        # three decades per two functions added. Whitening by M's Cholesky factor
        # buys usable room, but the limit is real and a caller must be able to see
        # it rather than discover it as a silently wrong answer -- which is what
        # happened before the whitening, with omega^2 collapsing to 1e-9 at K = 8.
        cs = [decay_rates(RitzBasis(2, K, :monomial), 0.05).cond_M for K in 2:8]
        @test issorted(cs)
        @test cs[1] < 1e3
        @test cs[end] > 1e9      # so an unbounded K is not safe, by construction
        ## THE DEFAULT BASIS IS NO LONGER THIS ONE. `:legendre` spans the same space with the
        ## same K=1 potential-flow limit and stays usable an order of magnitude further in K,
        ## which is why it is the default; the monomial ceiling above is what it replaced.
        cl = [decay_rates(RitzBasis(2, K, :legendre), 0.05).cond_M for K in 2:8]
        @test issorted(cl)
        @test cl[end] < cs[end] / 1e4
    end
end

@testset "shear thinning: the coupled assembly" begin

    @testset "constant viscosity leaves the modes uncoupled" begin
        # The diagonality result: with only its l = 0 harmonic, eta cannot move
        # energy between surface modes. Every off-diagonal block must be zero, and
        # this is what makes the single-mode assembly sufficient for Newtonian.
        b = ModalBasis(2:5, 3)
        F = assemble_coupled(b, 0.05)
        nl = length(b.ls)
        diag_scale = maximum(block_norm(b, F.C, i, i) for i in 1:nl)
        off = maximum(block_norm(b, F.C, i, j) for i in 1:nl, j in 1:nl if i != j)
        @test diag_scale > 1.0                 # not vacuous
        @test off < 1e-12 * diag_scale
        # the added mass is likewise block diagonal, by angular orthogonality
        offM = maximum(block_norm(b, F.M, i, j) for i in 1:nl, j in 1:nl if i != j)
        @test offM < 1e-12 * maximum(block_norm(b, F.M, i, i) for i in 1:nl)
    end

    @testset "an angular viscosity couples exactly the Gaunt band" begin
        # eta ~ P_k must couple l and m only when |l-m| <= k <= l+m with l+k+m even.
        # Checked against the assembled matrix, not against the rule's derivation.
        b = ModalBasis(2:6, 2)
        nl = length(b.ls)
        for k in 1:3
            F = assemble_coupled(b, 0.05;
                                 eta = (x, mu) -> DropSolver.legendre_angular(k, mu).P)
            sc = maximum(block_norm(b, F.C, i, j) for i in 1:nl, j in 1:nl)
            got = Set((b.ls[i], b.ls[j]) for i in 1:nl, j in 1:nl
                      if block_norm(b, F.C, i, j) > 1e-8 * sc)
            want = Set((l, m) for l in b.ls, m in b.ls
                       if abs(l - m) <= k <= l + m && iseven(l + k + m))
            @test got == want
        end
    end

    @testset "the shear rate does not superpose" begin
        # e superposes over modes; its invariant does not. So gammadot of a two-mode
        # state is NOT the sum of the single-mode shear rates -- if it were, "mode
        # l's shear rate" would be meaningful and the whole coupling would vanish.
        b = ModalBasis([2, 3], 1)
        x, mu = 0.7, -0.4
        g12 = shear_rate(b, [1.0, 1.0], x, mu)
        g1 = shear_rate(b, [1.0, 0.0], x, mu)
        g2 = shear_rate(b, [0.0, 1.0], x, mu)
        @test g1 > 0 && g2 > 0
        @test !isapprox(g12, g1 + g2; rtol=1e-6)
        @test !isapprox(g12, hypot(g1, g2); rtol=1e-6)   # nor is it Pythagorean
    end

    @testset "Carreau-Yasuda on a real state" begin
        b = ModalBasis(2:6, 2)
        nl = length(b.ls)
        a = [0.6 / i * (k == 1 ? 1.0 : 0.3) for i in 1:nl for k in 1:b.K]
        pts = [(x, mu) for x in (0.3, 0.6, 0.9), mu in (-0.95, -0.5, 0.0, 0.5, 0.95)]

        # the viscosity really does vary over the drop
        gd = [shear_rate(b, a, x, mu) for (x, mu) in pts]
        @test maximum(gd) / minimum(gd) > 3

        etaf = (x, mu) -> carreau(shear_rate(b, a, x, mu); lambda_c=8.0, a=2.0, n=0.4)
        ev = [etaf(x, mu) for (x, mu) in pts]
        @test all(0 .< ev .<= 1)                 # bounded by the plateaus, (H2)
        @test maximum(ev) / minimum(ev) > 2

        # and the coupling it induces is a leading-order effect, not a correction
        F = assemble_coupled(b, 0.05; eta = etaf)
        dmax = maximum(block_norm(b, F.C, i, i) for i in 1:nl)
        omax = maximum(block_norm(b, F.C, i, j) for i in 1:nl, j in 1:nl if i != j)
        @test omax / dmax > 0.05
        @test F.C ≈ F.C'                          # still a Hessian
    end

    @testset "the Newtonian limit is recovered" begin
        # lambda_c -> 0 sends eta -> 1 pointwise, so the coupled assembly must
        # collapse onto the constant-viscosity one and the modes must decouple again.
        b = ModalBasis(2:5, 2)
        nl = length(b.ls)
        a = fill(0.5, ndof(b))
        etaf = (x, mu) -> carreau(shear_rate(b, a, x, mu); lambda_c=1e-8, a=2.0, n=0.4)
        Fth = assemble_coupled(b, 0.05; eta = etaf)
        Fnw = assemble_coupled(b, 0.05)
        @test maximum(abs, Fth.C .- Fnw.C) < 1e-6 * maximum(abs, Fnw.C)
        off = maximum(block_norm(b, Fth.C, i, j) for i in 1:nl, j in 1:nl if i != j)
        @test off < 1e-6 * maximum(block_norm(b, Fth.C, i, i) for i in 1:nl)
    end

    @testset "single-mode and multi-mode assemblies agree" begin
        # Two independent code paths for the same object at constant viscosity.
        for l in (2, 4), K in (2, 3)
            Fs = assemble(RitzBasis(l, K), 0.1)
            Fm = assemble_coupled(ModalBasis([l], K), 0.1)
            @test Fs.M ≈ Fm.M  rtol=1e-10
            @test Fs.C ≈ Fm.C  rtol=1e-10
            @test Fs.G ≈ Fm.G  rtol=1e-10
        end
    end
end
