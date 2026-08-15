# A disk cache of finished impacts, so that changing the postprocessing never costs a
# simulation.
#
# WHY THIS IS CHEAP, AND WHY IT WAS NOT OBVIOUS. Restitution and contact time are read
# off four series: time, centre-of-mass height, centre-of-mass velocity, and the smallest
# gap between the surface and the substrate. Only the last touches the shape, and it is
# the expensive one -- but it does NOT depend on the contact threshold, which enters only
# as `minh .< h_thresh`. So one stored run can be re-scored at any threshold, under any
# later change of restitution convention, forever.
#
# What is deliberately NOT stored is the modal state at every step. It is two orders of
# magnitude larger and nothing downstream reads it. Storing the four series instead turns
# a run from tens of megabytes into a few hundred kilobytes.
#
# The key is every parameter that changes the trajectory. A key that omitted one would
# serve a stale run silently, which is the failure mode a cache has to be built against:
# the resolution belongs in the key for the same reason the basis kind belongs in the
# geometry cache's key one file over.
#
# Format is a small binary file rather than CSV: these are dense float series, and a text
# round trip costs more disk and more time than the simulation saves.

using Printf
using DropSolver

const CACHE_DIR = joinpath(@__DIR__, "..", "outputs", "jld", "runs")

"""
Everything that changes the trajectory.

Every field of `ImpactParams` goes into the key, rather than a hand-picked subset. A
subset is how a cache serves a stale run silently: omit `basis_kind` and a monomial run is
handed to a Legendre caller, which is the exact mistake the geometry cache one file over
records having made. Listing them all costs nothing and cannot go quietly wrong.

A run with a shear-dependent viscosity is refused rather than keyed. `eta` is always a
closure here -- even the Newtonian case carries one -- and a closure has no stable hash, so
two different rheologies would collide and the cache would answer with the wrong fluid. The
`eta_const` flag marks the case where the closure is a constant and `Oh` alone determines
it, which IS in the key. Checked before anything is written or read.
"""
cacheable(p::ImpactParams) = p.eta_const

function cache_key(p::ImpactParams; solver::Symbol = :active_set)
    parts = Any[solver]
    for f in fieldnames(typeof(p))
        f === :eta && continue                      # refused above, never keyed
        push!(parts, getfield(p, f))
    end
    string(hash(Tuple(parts)); base = 16)
end

function cache_path(p::ImpactParams; solver::Symbol = :active_set)
    mkpath(CACHE_DIR)
    joinpath(CACHE_DIR, "run_" * cache_key(p; solver = solver) * ".bin")
end

function cache_store(path, t, z, v, minh)
    ## The scratch name must be unique per writer, not just per key. Two curves on one
    ## figure can be the SAME liquid -- clean water is both the reference surface tension
    ## and the reference viscosity -- so two threads reach an identical key at the same
    ## moment. With a shared `path * ".tmp"` they interleave as: A writes tmp, B reopens
    ## and truncates the same tmp, A renames it into place, B renames and finds nothing
    ## there. That is an ENOENT that kills the run after an hour of correct work, and it
    ## depends on thread timing, so it does not reproduce on demand.
    ##
    ## `tempname` in the cache's own directory keeps the rename on one filesystem, which
    ## is what makes it atomic and what makes write-then-rename worth doing at all.
    tmp = tempname(dirname(path); cleanup = false)
    open(tmp, "w") do io                # write-then-rename: a killed job leaves no half file
        write(io, Int64(length(t)))
        for a in (t, z, v, minh); write(io, Float64.(a)); end
    end
    try
        mv(tmp, path; force = true)
    catch
        ## Another writer got there first with byte-identical content -- same key means
        ## same run. Drop ours rather than fail the sweep.
        rm(tmp; force = true)
    end
end

function cache_load(path)
    isfile(path) || return nothing
    try
        open(path, "r") do io
            n = read(io, Int64)
            (0 < n < 100_000_000) || return nothing
            rd() = (a = Vector{Float64}(undef, n); read!(io, a); a)
            (t = rd(), z = rd(), v = rd(), minh = rd())
        end
    catch
        nothing                          # a truncated file is a miss, not a crash
    end
end

"""
    cached_series(p; solver = :active_set, fallback = true) -> NamedTuple or nothing

The four series for this impact, simulated only if they are not already on disk.

`fallback` retries with the complementarity closure when the active set finishes still in
contact, which is its documented failure mode at high Weber number; the retry is part of
the trajectory and so part of what is cached.
"""
function cached_series(p::ImpactParams; solver::Symbol = :active_set, fallback::Bool = true)
    use = cacheable(p)
    path = use ? cache_path(p; solver = solver) : ""
    if use
        hit = cache_load(path)
        hit === nothing || return hit
    end

    r = solver === :lcp ? simulate_lcp(p) : simulate(p)
    if fallback && solver === :active_set &&
       !isfinite(r.cor) && isfinite(r.tc) && r.tc >= 0.9 * p.t_max
        r = simulate_lcp(p)
    end
    minh = min_gap_series(p, r)
    use && cache_store(path, r.t, r.z, r.v, minh)
    (t = r.t, z = r.z, v = r.v, minh = minh)
end

"""
    score(p, s; h_thresh = 0.02) -> NamedTuple

Restitution and contact time from cached series, through the same code path a live run
uses -- `proximity_metrics` is handed the stored `minh` rather than recomputing it.
"""
score(p::ImpactParams, s; h_thresh::Real = 0.02) =
    proximity_metrics(p, (t = s.t, z = s.z, v = s.v, a = nothing);
                      h_thresh = h_thresh, minh = s.minh)

"""Number of runs on disk, and their total size in megabytes."""
function cache_stats()
    isdir(CACHE_DIR) || return (n = 0, mb = 0.0)
    fs = filter(f -> endswith(f, ".bin"), readdir(CACHE_DIR; join = true))
    (n = length(fs), mb = sum(filesize, fs; init = 0) / 1e6)
end
