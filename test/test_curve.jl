# The adaptive sampler the figure scripts share, and the one place it can make a real
# discontinuity look identical to an artefact of stopping at a declared boundary.
#
# `scripts/_curve.jl` is not part of the package -- it is tooling every restitution-vs-We
# figure calls -- but it decides which points a figure draws, so a fault in it reaches
# every one of those figures at once. Tested with synthetic curves rather than the solver,
# because the logic under test is about samples and gaps, not physics.

using Test

include(joinpath(@__DIR__, "..", "scripts", "_curve.jl"))

@testset "curve sampling and the ceiling-jump guard" begin

    @testset "adaptive_curve resolves a smooth roll-off" begin
        ## A logistic step, standing in for the geometric roll-off every restitution
        ## curve has: zero below threshold, a plateau above it.
        f(x) = 1 / (1 + exp(-40 * (log10(x) - log10(0.01))))
        xs, ys = adaptive_curve(f, 1e-4, 1.0; th = 0.02, x_max = 0.05, x_th = 0.005,
                                n0 = 6, maxpts = 200)
        @test issorted(xs)
        gaps = [abs(ys[i+1] - ys[i]) for i in 1:length(ys)-1
                if log10(xs[i+1] / xs[i]) > 0.005]
        @test all(<=(0.02), gaps)          # the y-tolerance held everywhere refinable
    end

    @testset "jumps finds a genuine step and nothing else" begin
        ## A curve with one deliberate step, narrower in x than x_th, and otherwise a
        ## smooth decline -- the shape a real jump has among ordinary bends.
        xs = [1.0, 1.5, 2.0, 2.5, 2.501, 3.0]
        ys = [0.90, 0.85, 0.80, 0.75, 0.60, 0.55]     # the step is between x=2.5 and 2.501
        js = jumps(xs, ys, 0.05; x_th = 0.01)
        @test length(js) == 1
        @test js[1][1] == 2.5 && js[1][2] == 2.501
        @test isapprox(js[1][3], 0.15)

        ## The ordinary bends (0.05 apart, well separated in x) must not be flagged.
        js_wide = jumps(xs, ys, 0.02; x_th = 0.01)
        @test all(j -> j[1] == 2.5, js_wide)     # only the narrow step qualifies
    end

    @testset "trim_ceiling_jump drops an artefact at the boundary, and only that" begin
        ## Smooth decline reaching the ceiling with no jump: left alone.
        xs = [1.0, 2.0, 3.0, 4.0, 5.0]
        ys = [0.90, 0.85, 0.80, 0.75, 0.70]
        xs2, ys2 = trim_ceiling_jump(xs, ys, 5.0, 0.05; x_th = 0.02)
        @test xs2 == xs && ys2 == ys

        ## The BSI symptom this guards against: smooth out to the point before the
        ## ceiling, then a step into the ceiling itself, narrower in x than x_th.
        xs = [1.0, 2.0, 3.0, 4.0, 4.99, 5.0]
        ys = [0.90, 0.85, 0.80, 0.75, 0.70, 0.55]
        xs2, ys2 = trim_ceiling_jump(xs, ys, 5.0, 0.05; x_th = 0.02)
        @test xs2 == xs[1:end-1] && ys2 == ys[1:end-1]

        ## A step of the same size, but NOT at the declared ceiling: an interior
        ## feature, not a boundary artefact, and must survive untouched.
        xs = [1.0, 2.0, 2.99, 3.0, 4.0]
        ys = [0.90, 0.85, 0.80, 0.65, 0.60]
        xs2, ys2 = trim_ceiling_jump(xs, ys, 5.0, 0.05; x_th = 0.02)
        @test xs2 == xs && ys2 == ys

        ## A step at the ceiling that is too WIDE in x to count as a jump (a genuine
        ## bend the sampler already resolved): left alone.
        xs = [1.0, 2.0, 3.0, 4.0, 5.0]
        ys = [0.90, 0.85, 0.80, 0.75, 0.55]
        xs2, ys2 = trim_ceiling_jump(xs, ys, 5.0, 0.05; x_th = 0.02)
        @test xs2 == xs && ys2 == ys
    end
end
