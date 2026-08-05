# What an experimentalist expects to see, tested on our own solver.
#
# The rest of the suite checks the model against itself: operators against their symmetries,
# limits against Reid and Lamb, the stepper against the laws it should conserve. All of that
# can pass while the drop behaves in a way nobody has ever filmed. These five tests are the
# other direction -- each states something a person watching high-speed video of a drop
# bouncing off a non-wetting surface would report, and then asks whether we reproduce it.
#
# They are deliberately about SHAPES OF DEPENDENCE rather than single numbers. An experiment
# has scatter, so "restitution is 0.73" is weak evidence about a model; "contact time barely
# moves when the impact speed changes by a factor of fifteen, while restitution falls by a
# third over the same range" is strong evidence, because a model that has the physics wrong
# has no reason to get that contrast right.
#
# Every tolerance below was measured before it was asserted.

using Test
using DropSolver
using LinearAlgebra

basis_of_exp(p) = DropSolver.basis(p)

## Statistics is not a bundled stdlib in Julia 1.12, and only the mean is needed here
mean(v) = sum(v) / length(v)

"""Least-squares slope of log y against log x -- the exponent of a power law."""
function logslope(x, y)
    lx, ly = log.(x), log.(y)
    n = length(lx)
    (n*sum(lx.*ly) - sum(lx)*sum(ly)) / (n*sum(lx.^2) - sum(lx)^2)
end

spread(v) = (maximum(v) - minimum(v)) / mean(v)

@testset "experimental expectations" begin

    @testset "contact time barely depends on how fast the drop arrives" begin
        # THE CANONICAL MEASUREMENT on bouncing drops (Richard, Clanet & Quere, Nature 2002):
        # over roughly a decade of impact speed the contact time does not move, and its value
        # is about 2.6 times the capillary time sqrt(rho R^3 / sigma) -- which is our unit of
        # time, so the number to look for is 2.6 itself.
        #
        # This is a strong test because it is a NON-result: the drop spends longer on the
        # surface neither when it arrives gently nor when it arrives hard. A model that got
        # the stiffness or the added mass wrong would show contact time tracking the speed.
        # Measured here: 2.15 to 2.73 over a fifteenfold change in Weber number, a slope of
        # -0.085 where a naive contact-mechanics picture would give something of order 1/2.
        for Oh in (0.0373, 0.0767)
            Wes = exp.(range(log(0.2), log(3.0); length = 6))
            tcs = Float64[]; cors = Float64[]; ok = Float64[]
            for We in Wes
                p = ImpactParams(We = We, Bo = 0.02, Oh = Oh, M = 45, K = 2, t_max = 25.0)
                m = proximity_metrics(p, simulate(p))
                (isfinite(m.tc) && 0 < m.tc < 20) || continue
                push!(tcs, m.tc); push!(cors, m.cor); push!(ok, We)
            end
            @test length(tcs) == length(Wes)              # every case must run
            ## flat in We, and near the measured 2.6
            @test abs(logslope(ok, tcs)) < 0.15
            @test spread(tcs) < 0.30
            @test 2.0 < mean(tcs) < 3.0
            ## and the flatness is not because nothing is happening: restitution DOES move
            @test spread(cors) > 0.30
        end
    end

    @testset "contact time is set by capillarity, not by weight or by viscosity" begin
        # The same non-result from two other directions, and this is the sharpest contrast in
        # the suite. Surface tension provides the spring, so the time on the surface should be
        # a capillary time and nothing else. Make the drop twenty times heavier: contact time
        # must hardly move. Make it three hundred times more viscous: contact time must
        # hardly move THERE EITHER, even though restitution collapses.
        #
        # A model in which the film, gravity or the damping were setting the contact duration
        # would fail this while still fitting a restitution curve, which is why it is worth
        # asserting the insensitivity and the sensitivity together.
        We = 0.5
        Bos = (0.005, 0.02, 0.05, 0.1)                    # a factor of twenty in weight
        tb = Float64[]; cb = Float64[]
        for Bo in Bos
            p = ImpactParams(We = We, Bo = Bo, Oh = 0.0373, M = 45, K = 2, t_max = 25.0)
            m = proximity_metrics(p, simulate(p))
            push!(tb, m.tc); push!(cb, m.cor)
        end
        @test spread(tb) < 0.05                           # measured 1.6 per cent
        @test spread(cb) > 0.05                           # while restitution notices: 7 per cent

        Ohs = (0.001, 0.01, 0.03, 0.1, 0.3)               # a factor of three hundred
        to = Float64[]; co = Float64[]
        for Oh in Ohs
            p = ImpactParams(We = 0.2, Bo = 0.02, Oh = Oh, M = 45, K = 2, t_max = 30.0)
            m = proximity_metrics(p, simulate(p))
            push!(to, m.tc); push!(co, m.cor)
        end
        @test spread(to) < 0.05                           # measured 2.4 per cent
        @test spread(co) > 0.50                           # while restitution halves: 0.95 -> 0.50
    end

    @testset "where the lost energy goes crosses over as the drop thins" begin
        # A drop never returns all of its downward motion, and the reason is not the same at
        # both ends of the viscosity range -- which is exactly the sort of thing that is easy
        # to get wrong and hard to notice, because both ends produce a restitution below one.
        #
        #   viscous drop   : the deficit is HEAT. It leaves nearly unshaken.
        #   inviscid drop  : the deficit is MOTION STILL IN THE DROP. Almost nothing is
        #                    dissipated; the drop leaves oscillating, and an experimentalist
        #                    sees exactly that on the video.
        #
        # So restitution below one is not evidence of dissipation, and this test is what
        # separates the two. Measured: dissipation falls from 91 per cent of the arriving
        # kinetic energy at Oh = 0.3 to 2.4 per cent at Oh = 0.001, while the share of the
        # deficit still sitting in the shape modes rises from 1 per cent to 81 per cent.
        function partition(We, Oh)
            p = ImpactParams(We = We, Bo = 0.02, Oh = Oh, M = 20, K = 2, t_max = 30.0)
            r = simulate(p)
            F = assemble_newtonian(basis_of_exp(p), Oh); mass = 4pi/3
            KEin  = 0.5*mass*r.v[1]^2
            defic = KEin - 0.5*mass*r.v[end]^2
            Esh   = 0.5*dot(r.adot[end], F.M, r.adot[end]) + 0.5*dot(r.a[end], F.G, r.a[end])
            Dis = 0.0
            for i in 1:length(r.t)-1
                h = r.t[i+1] - r.t[i]
                Dis += 0.5*h*(dot(r.adot[i],   F.C, r.adot[i]) +
                              dot(r.adot[i+1], F.C, r.adot[i+1]))
            end
            (heat = Dis/KEin, shaken = Esh/defic, cor = 1 - defic/KEin)
        end
        for We in (0.2, 0.5)
            lo = partition(We, 0.001)      # as close to inviscid as the stepper is happy at
            hi = partition(We, 0.3)        # genuinely viscous
            ## the inviscid drop barely dissipates, and keeps the deficit as motion
            @test lo.heat   < 0.10
            @test lo.shaken > 0.70
            ## the viscous drop does the opposite
            @test hi.heat   > 0.80
            @test hi.shaken < 0.10
            ## and the crossover is monotone through the middle, not a coincidence at the ends
            hs = [partition(We, Oh).heat for Oh in (0.001, 0.01, 0.1, 0.3)]
            @test issorted(hs)
        end
    end

    @testset "gentle impacts deform in proportion to their speed" begin
        # In the linear regime the drop is a driven oscillator, so its deformation should be
        # proportional to the impact SPEED, not the impact energy:
        #
        #     max |zeta_2|  ~  sqrt(We)
        #
        # An exponent of one half, and nothing else will do -- 1 would mean the amplitude
        # follows the energy, 1/4 would be the flattening scaling of a violent impact. This
        # is measurable off a single frame per run and it pins the response to be LINEAR,
        # which is the assumption the whole small-amplitude model rests on. Measured 0.518 at
        # low Ohnesorge and 0.481 at high, over a twentyfold range of Weber number.
        for Oh in (0.0233, 0.3038)
            Wes = exp.(range(log(0.02), log(0.4); length = 5))
            amps = Float64[]
            for We in Wes
                p = ImpactParams(We = We, Bo = 0.02, Oh = Oh, M = 30, K = 2, t_max = 25.0)
                r = simulate(p)
                push!(amps, maximum(abs(surface_amplitudes(p, a)[1]) for a in r.a))
            end
            @test all(amps[i] < amps[i+1] for i in 1:length(amps)-1)   # monotone, obviously
            @test 0.44 < logslope(Wes, amps) < 0.57                    # and the exponent is 1/2
            @test maximum(amps) < 0.5                                  # still small-amplitude
        end
    end

    @testset "a drop left alone rings at the frequency drops are seen to ring at" begin
        # The oldest measurement in the subject: a free drop oscillates at Rayleigh's
        # frequency, omega_l^2 = l(l-1)(l+2) in capillary units -- 2.8284 for the ellipsoidal
        # l = 2 mode, 5.4772 for l = 3. Every other frequency check here is an EIGENVALUE of
        # the assembled operators; this one integrates the equations forward and counts zero
        # crossings, which is what a camera does, and it also checks that viscosity lowers the
        # frequency rather than raising it.
        #
        # A failing l = 2 would mean the stiffness or the added mass is wrong. An l = 2 that
        # passed while l = 3 failed would mean the modes are right individually but mislabelled
        # -- the surface would ring at the wrong shape, which no eigenvalue test detects.
        function ring(Oh, l)
            b = ModalBasis(2:6, 2)
            F = assemble_newtonian(b, Oh)
            n = DropSolver.ndof(b); j = (l-2)*2 + 1
            a = zeros(n); a[j] = 0.05; ad = zeros(n)      # excite mode l alone
            Mf = factorize(Matrix(F.M))
            dt = 1e-3; ts = Float64[]; zs = Float64[]
            for k in 1:60_000
                ad .+= dt .* (-(Mf \ (F.C*ad + F.G*a)))
                a  .+= dt .* ad
                push!(ts, k*dt); push!(zs, a[j] + a[j+1]) # zeta_l is the sum of its radial terms
            end
            zc = [k for k in 1:length(zs)-1 if zs[k]*zs[k+1] < 0]
            @assert length(zc) > 10
            2pi / (2*mean(diff(ts[zc])))
        end
        for l in (2, 3)
            ray = sqrt(l*(l-1)*(l+2))
            @test isapprox(ring(1e-4, l), ray; rtol = 1e-4)   # measured to 1e-5 relative
        end
        ## damping lowers the frequency -- omega_d = sqrt(omega_0^2 - gamma^2) < omega_0
        @test ring(0.01, 2) < ring(1e-4, 2)
        @test isapprox(ring(0.01, 2), sqrt(8.0); rtol = 3e-3)  # and only slightly: 0.08 per cent
    end
end
