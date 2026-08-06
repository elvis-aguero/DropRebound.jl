# Every equation printed on a hand-authored documentation page, checked against the
# code that runs.
#
# The derivation chapters execute in CI, so their algebra cannot drift from the solver
# without something going red. The hand-authored pages (`index.md`, `variational.md`,
# `contact.md`, `resolution.md`) are prose, and prose does not execute. This file is
# what stops them stating a formula the solver does not implement.
#
# Each testset names the page and the claim. A failure means the page is wrong, or the
# code changed and the page was not updated. Either way the site is lying to a reader,
# which is the failure mode this file exists to catch.

using Test
using LinearAlgebra
using Random
using DropSolver
using DropSolver: RitzBasis, ModalBasis, phi, dphi, d2phi, legendre_angular,
                  modal_field, strain_at, stiffness_matrix, assemble_newtonian,
                  ndof, dofindex, radial_window, ImpactParams, gap, gap_row,
                  force_column, trace_vec, basis, pc_len

# A handful of representative modes and sample points. Nothing here depends on the
# specific values; they are spread out so that a formula that happens to hold at one
# point does not pass by luck.
const LS  = (2, 3, 5, 12)
const XS  = (0.17, 0.43, 0.88, 1.0)
const MUS = (-0.93, -0.31, 0.22, 0.77)

@testset "docs math: variational.md" begin

    # PAGE CLAIM: "a field written this way satisfies div u = 0 identically, for every
    # Psi". The page argues the trace of the strain-rate tensor cancels because
    # Legendre's equation supplies d2P/dth2 + cot(th) dP/dth = -l(l+1) P.
    #
    # If this fails, the trial space is not incompressible, and the entire justification
    # for the pressure having left the problem collapses: the weak form dropped the
    # pressure term by invoking div v = 0 on the test fields.
    @testset "incompressibility is identical, not approximate" begin
        for l in LS, x in XS, mu in MUS
            A = legendre_angular(l, mu)
            rb = RitzBasis(l, 4, :legendre)
            for k in 1:4
                f, df, d2f = phi(rb, k, x), dphi(rb, k, x), d2phi(rb, k, x)
                (_, _, e_rr, e_tt, e_pp, _) = modal_field(l, f, df, d2f, x, A)
                @test isapprox(e_rr + e_tt + e_pp, 0.0; atol = 1e-10)
            end
        end
    end

    # PAGE CLAIM: the Legendre identity the cancellation above rests on. Checked
    # separately so that a failure distinguishes "the identity is wrong" from "the
    # strain components are wrong".
    @testset "Legendre's equation in theta" begin
        for l in LS, mu in MUS
            A = legendre_angular(l, mu)
            cot_th = A.mu / A.sinth
            @test isapprox(A.d2Pdth + cot_th * A.dPdth, -l * (l + 1) * A.P; atol = 1e-10)
        end
    end

    # PAGE CLAIM: the published velocity and strain components,
    #   u_r     = f/x^2 P_l
    #   u_th    = f'/(x l(l+1)) dP/dth
    #   e_rr    = (f'/x^2 - 2f/x^3) P_l
    #   e_rth   = (1/(2 x l(l+1))) (f'' - 2f'/x + l(l+1) f/x^2) dP/dth
    # A reader who differentiates these by hand must land on what the solver assembles.
    @testset "published field components match modal_field" begin
        for l in LS, x in XS, mu in MUS
            A = legendre_angular(l, mu); L = l * (l + 1)
            rb = RitzBasis(l, 3, :legendre)
            for k in 1:3
                f, df, d2f = phi(rb, k, x), dphi(rb, k, x), d2phi(rb, k, x)
                (u_r, u_th, e_rr, _, _, e_rt) = modal_field(l, f, df, d2f, x, A)
                @test isapprox(u_r,  f / x^2 * A.P;                          atol = 1e-12)
                @test isapprox(u_th, df / (x * L) * A.dPdth;                 atol = 1e-12)
                @test isapprox(e_rr, (df / x^2 - 2f / x^3) * A.P;            atol = 1e-12)
                T_op = d2f - 2df / x + L * f / x^2
                @test isapprox(e_rt, T_op / (2 * x * L) * A.dPdth;           atol = 1e-12)
            end
        end
    end

    # PAGE CLAIM: "G is rank one for a single mode", with
    #   G_ab = (4 pi/(2l+1)) (l-1)(l+2) phi_a(1) phi_b(1).
    #
    # Physically: surface energy sees the interior only through the boundary trace, so
    # the stiffness cannot distinguish two interior fields with the same surface
    # amplitude. A failure of the rank test would mean capillarity is somehow reaching
    # into the interior, which would make the page's account of why refining K does not
    # refine the restoring force wrong.
    @testset "stiffness is rank one per mode, with the published coefficient" begin
        K = 4
        for l in LS
            b = ModalBasis([l], K, :legendre)
            G = stiffness_matrix(b)
            rb = RitzBasis(l, K, :legendre)
            tr = [phi(rb, k, 1.0) for k in 1:K]
            @test isapprox(G, (4pi / (2l + 1)) * (l - 1) * (l + 2) * (tr * tr'); atol = 1e-10)
            @test rank(G; atol = 1e-8 * maximum(abs, G)) == 1
        end
    end

    # PAGE CLAIM: "the matrices are Hessians of scalar functionals ... symmetric by
    # construction rather than by cancellation".
    @testset "M and C are symmetric" begin
        b = ModalBasis(2:8, 3, :legendre)
        M, C, G = assemble_newtonian(b, 0.05)
        for (name, A) in (("M", M), ("C", C), ("G", G))
            @test norm(A - A') <= 1e-10 * norm(A)
        end
    end
end

@testset "docs math: resolution.md" begin

    # PAGE CLAIM: "Setting K = 1 keeps only the first [term]. That is a potential flow."
    #
    # For a Stokes stream function Psi = f(x) C_l(theta), the flow is irrotational when
    # E^2 Psi = 0, that is f'' - l(l+1) f / x^2 = 0. The first trial function is
    # x^(l+1) in both bases, and it must satisfy this exactly.
    #
    # If it did not, the page's explanation of why K = 1 reproduces Lamb's over-damped
    # answer would be wrong: the whole argument is that one trial function admits no
    # vorticity, so the tangential-stress condition cannot be met.
    @testset "K = 1 is irrotational in both bases" begin
        for kind in (:monomial, :legendre), l in LS, x in XS
            rb = RitzBasis(l, 1, kind)
            f, d2f = phi(rb, 1, x), d2phi(rb, 1, x)
            @test isapprox(d2f - l * (l + 1) * f / x^2, 0.0; atol = 1e-9)
        end
    end

    # PAGE CLAIM: the exact radial profile is A x^(l+1) + B x j_l(q x), an irrotational
    # part plus a vortical one. The vortical part must satisfy the viscous radial
    # equation f'' - l(l+1) f/x^2 + q^2 f = 0.
    #
    # This is what makes the trial functions the Taylor terms of a Bessel function, and
    # so what makes the requirement K >~ |q|/1.5 the right shape of statement.
    @testset "x j_l(q x) solves the vortical equation" begin
        # Spherical Bessel by its power series, which converges quickly for the
        # arguments used here and avoids the instability of upward recurrence at small
        # argument. No new dependency: the package carries only a Bessel *ratio*.
        function sphj(l::Int, z::Real)
            dfact = 1.0
            for m in 3:2:(2l + 1)
                dfact *= m
            end
            s, t = 0.0, 1.0
            for k in 0:80
                s += t
                t *= -(z^2 / 2) / ((k + 1) * (2l + 2(k + 1) + 1))
            end
            z^l / dfact * s
        end
        @test isapprox(sphj(0, 1.7), sin(1.7) / 1.7; atol = 1e-12)   # series sanity

        h = 1e-4
        fv(l, q, x) = x * sphj(l, q * x)
        for l in LS, q in (1.3, 4.0, 9.5), x in (0.3, 0.6, 0.9)
            d2f = (fv(l, q, x + h) - 2 * fv(l, q, x) + fv(l, q, x - h)) / h^2
            resid = d2f - l * (l + 1) * fv(l, q, x) / x^2 + q^2 * fv(l, q, x)
            @test abs(resid) <= 1e-5 * max(1.0, abs(fv(l, q, x)))
        end
    end

    # PAGE CLAIM: "K(l) = 47.2 l^(-0.642)".
    @testset "the published radial window is the one the code uses" begin
        for l in (2, 5, 10, 20, 45, 90)
            @test radial_window(l; margin = 1.0) ==
                  clamp(floor(Int, 47.2 * float(l)^(-0.642)), 1, 40)
        end
    end

    # PAGE CLAIM: "Both families span the same space ... where both are numerically
    # sound they agree to eight decimal places."
    @testset "the two radial bases agree where both are conditioned" begin
        b_mono = ModalBasis(2:6, 3, :monomial)
        b_legd = ModalBasis(2:6, 3, :legendre)
        Mm, Cm, Gm = assemble_newtonian(b_mono, 0.05)
        Ml, Cl, Gl = assemble_newtonian(b_legd, 0.05)
        # Same span means the same generalised eigenvalues, not the same matrices.
        em = sort(real(eigvals(Gm, Mm)))
        el = sort(real(eigvals(Gl, Ml)))
        @test maximum(abs, em .- el) <= 1e-8 * maximum(abs, em)
    end
end

@testset "docs math: contact.md" begin

    p = ImpactParams(We = 0.5, Bo = 0.0189, Oh = 0.05, M = 12, K = 2)
    b = basis(p)

    # PAGE CLAIM: "h(theta) = z + cos(theta) (1 + sum_l zeta_l P_l(cos theta))".
    #
    # The page leans on two features of this expression: that it is affine in the state,
    # and that it carries a factor cos(theta) because the constraint is on vertical
    # clearance while zeta is a radial displacement. Both are used later in the
    # conjugacy argument, so an error here would propagate to that conclusion.
    @testset "the published clearance is the one imposed" begin
        Random.seed!(20260806)
        a = 0.01 .* randn(ndof(b)); z = 1.4
        for th in (pi, 2.7, 2.1, 1.6)
            mu = cos(th)
            zeta = [dot(trace_vec(p, l),
                        view(a, dofindex(b, i, 1):dofindex(b, i, p.K)))
                    for (i, l) in enumerate(b.ls)]
            published = z + mu * (1 + sum(zeta[i] * legendre_angular(l, mu).P
                                          for (i, l) in enumerate(b.ls)))
            @test isapprox(gap(p, a, z, th), published; atol = 1e-12)
        end
    end

    # PAGE CLAIM: "H_il = cos(theta_i) P_l(cos(theta_i))", the constraint Jacobian read
    # off the clearance. This is one of the two objects the symmetry argument compares.
    @testset "constraint Jacobian has the published entries" begin
        for th in (pi, 2.6, 1.9)
            row, cst = gap_row(p, th)
            mu = cos(th)
            @test isapprox(cst, mu; atol = 1e-14)
            for (i, l) in enumerate(b.ls)
                Pl = legendre_angular(l, mu).P
                tv = trace_vec(p, l)
                for k in 1:p.K
                    @test isapprox(row[dofindex(b, i, k)], mu * Pl * tv[k]; atol = 1e-12)
                end
            end
        end
    end

    # PAGE CLAIM: "Q_l = -(4 pi/(2l+1)) p_{c,l}", and "the l = 0 and l = 1 harmonics do
    # no work on the shape". The second is why the page can say the l = 1 harmonic acts
    # only through the centre-of-mass equation.
    @testset "generalised force has the published coefficient" begin
        for j in 1:pc_len(p)
            col = force_column(p, j)
            l_pc = j - 1                       # film harmonics run l = 0..M
            if l_pc in b.ls
                i = findfirst(==(l_pc), b.ls)
                tv = trace_vec(p, l_pc)
                for k in 1:p.K
                    @test isapprox(col[dofindex(b, i, k)],
                                   -(4pi / (2 * l_pc + 1)) * tv[k]; atol = 1e-12)
                end
            else
                # l = 0 changes volume and l = 1 translates: neither does work on a
                # retained shape mode, so the column is zero.
                @test norm(col) <= 1e-14
            end
        end
    end

    # PAGE CLAIM: "A_c = H A^-1 Q_n = -H A^-1 H^T is symmetric for any H, because A is."
    #
    # This is the pivot of the whole convexity section: it says asymmetry of the
    # compliance is a statement about conjugacy of the forcing, not about the
    # discretisation being sloppy. Stated as pure linear algebra, so it is checkable
    # without running a step.
    @testset "conjugate forcing would give a symmetric compliance" begin
        Random.seed!(20260806)
        n, m = 9, 4
        S = randn(n, n); A = S * S' + n * I          # symmetric positive definite
        H = randn(m, n)
        W = Diagonal(rand(m) .+ 0.5)                 # any positive diagonal
        Ac_conj = H * (A \ (-H' * W))
        # Symmetric up to the weighting: W^(1/2) scales it to an exactly symmetric form.
        Wh = sqrt.(W)
        sym = Wh * Ac_conj * inv(Wh)
        @test norm(sym - sym') <= 1e-8 * norm(sym)

        # And a NON-conjugate forcing generally does not. This is the control: without
        # it the test above would pass for a reason unrelated to conjugacy.
        Qn = randn(n, m)
        Ac_free = H * (A \ Qn)
        @test norm(Ac_free - Ac_free') > 1e-6 * norm(Ac_free)
    end
end
