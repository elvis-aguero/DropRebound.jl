# The disk cache of finished impacts, and the one way it is easy to get wrong.
#
# `scripts/_runcache.jl` is not part of the package -- it is tooling the figure scripts
# share -- but it decides what every published figure is drawn from, so a fault in it is
# a fault in the figures. It is tested here rather than left to be discovered by a sweep
# that dies an hour in.
#
# The behaviour under test is concurrent storage of ONE key. That is not a contrived
# case: a sensitivity figure sweeps one property at a time about a reference liquid, so
# the reference appears in every panel, and `parallel_curves` puts those panels on
# different threads. Two threads then compute an identical parameter set at the same
# moment and race to write the same cache entry.

using Test
using DropSolver

include(joinpath(@__DIR__, "..", "scripts", "_runcache.jl"))

@testset "run cache" begin

    ## A key that no real sweep would produce, so a stray hit cannot make this pass.
    probe() = ImpactParams(We = 0.123456789, Bo = 0.00987654321, Oh = 0.0246813579,
                           M = 8, K = 1, t_max = 1.0)

    @testset "every trajectory-changing field is in the key" begin
        # A key built from a hand-picked subset is how a cache serves a stale run
        # silently. Changing any one field must change the key; if one is ever dropped
        # from `cache_key`, the run for a different physical problem is handed back and
        # the figure is wrong with nothing to show for it.
        ## `M` is a constructor keyword rather than a stored field -- it sets `ls` -- so
        ## the base case is kept as the keywords it was built from, not read back off
        ## the struct.
        base_kw = (We = 0.123456789, Bo = 0.00987654321, Oh = 0.0246813579,
                   M = 8, K = 1, t_max = 1.0)
        k0 = cache_key(ImpactParams(; base_kw...))
        for (f, v) in ((:We, 0.2), (:Bo, 0.02), (:Oh, 0.05), (:M, 10), (:K, 2),
                       (:t_max, 2.0))
            @test cache_key(ImpactParams(; merge(base_kw, (; f => v))...)) != k0
        end
        ## Same keywords, same key -- otherwise nothing would ever hit.
        @test cache_key(ImpactParams(; base_kw...)) == k0
        ## and the solver that produced it is part of the identity too
        @test cache_key(ImpactParams(; base_kw...); solver = :lcp) !=
              cache_key(ImpactParams(; base_kw...); solver = :active_set)
    end

    @testset "a shear-dependent viscosity is refused, not keyed" begin
        # `eta` is a closure, and a closure has no stable hash, so two different
        # rheologies would collide and the cache would answer with the wrong fluid.
        @test cacheable(probe())                       # constant viscosity: fine
        @test !cacheable(ImpactParams(We = 0.5, Bo = 0.02, Oh = 0.03, M = 8, K = 1,
                                      eta = gd -> 1 / (1 + gd)))
    end

    @testset "a truncated file is a miss, not a crash" begin
        # A job killed mid-write must not poison the cache for every later run.
        p = probe()
        path = cache_path(p)
        rm(path; force = true)
        n = 40
        t = collect(1.0:n)
        cache_store(path, t, t, t, t)
        @test cache_load(path) !== nothing

        full = read(path)
        write(path, full[1:(length(full) ÷ 3)])        # chop it mid-series
        @test cache_load(path) === nothing             # a miss
        rm(path; force = true)
        @test cache_load(path) === nothing             # and so is absent
    end

    @testset "concurrent writers of one key do not destroy each other" begin
        # THE REGRESSION. `cache_store` used to write to `path * ".tmp"`, a name unique
        # per KEY but not per WRITER, and then rename it into place. Two threads on the
        # same key interleave as: A writes tmp, B reopens and truncates that same tmp, A
        # renames it away, B renames and finds nothing there. The result is an ENOENT
        # that takes down a sweep after an hour of correct work, and because it depends
        # on thread timing it does not reproduce on demand.
        #
        # A failure here means any figure whose curves share a reference liquid can die
        # partway through, or -- worse, if the rename were made non-fatal without making
        # the scratch name unique -- can leave a half-written file that a later run reads
        # back as a trajectory.
        p = probe()
        path = cache_path(p)
        rm(path; force = true)
        n = 64
        t = collect(1.0:n); z = 2t; v = 3t; m = 4t

        Threads.@threads for _ in 1:64
            cache_store(path, t, z, v, m)
        end

        s = cache_load(path)
        @test s !== nothing                            # somebody won, and wrote it whole
        @test s.t == t && s.z == z && s.v == v && s.minh == m

        ## and nobody left scratch behind
        @test isempty(filter(f -> occursin("jl_", f) || endswith(f, ".tmp"),
                             readdir(dirname(path))))
        rm(path; force = true)
    end
end
