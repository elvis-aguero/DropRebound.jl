#!/usr/bin/env julia
# A provenance-first sweeper.
#
# Simulations here cost minutes, and the same points were being recomputed every time a
# comparison was rerun. This stores every completed run keyed by its full parameter set, so
# a repeated sweep does no work, and a sweep that is interrupted resumes where it stopped.
#
# THREE THINGS IT IS CAREFUL ABOUT.
#
# Provenance. Each row carries everything needed to reproduce it: the solver, every
# dimensionless group, the truncation, the rheology, the metric threshold, the git commit,
# and the timestamp. A KPI without the commit that produced it is not a result, because this
# solver's numbers have moved several times for reasons that were only visible in the diff.
#
# Memory. A run at M = 90, K = 2 holds roughly nine megabytes of trajectory -- the interior
# amplitudes, their rates and the pressure at every step. Multiplied by the thread count
# that is real, so each worker reduces its trajectory to scalars and drops it before the next
# case starts, and the concurrency is capped by an explicit budget rather than by the core
# count alone.
#
# Idempotence. Rows are appended and the manifest is read back before each sweep, so
# interrupting this is safe and rerunning it is free.
#
# Usage:
#   julia --project=julia -t auto julia/scripts/sweep.jl [--budget-mb N] [--force]

using Printf, Dates, Base.Threads
using DropSolver
include(joinpath(@__DIR__, "_stats.jl"))

const RESULTS = joinpath(@__DIR__, "..", "..", "results")
const STORE   = joinpath(RESULTS, "runs.csv")
const COLS = ["solver","basis","We","Bo","Oh","M","K","rheology","t_max","h_thresh",
              "cor","tc","cor_internal","tc_internal","min_gap","n_detected",
              "maxcp","rejects","gap_fraction","lcp_resid","eta_sweeps","wall_s",
              "commit","stamp"]

git_commit() = try
    strip(read(`git -C $(joinpath(@__DIR__,"..","..")) rev-parse --short HEAD`, String))
catch; "unknown" end

"""
The identity of a run: everything that changes its numbers, and nothing that does not.

A KNOWN HOLE, worth stating because it has already cost a set of results. The key covers the
parameters but not the SOLVER'S OWN VERSION, so a row stays valid-looking after the algorithm
behind it changes. When the complementarity closure was corrected -- it had been solving a
symmetrised surrogate of an asymmetric compliance -- all 38 stored `lcp` rows became wrong
while still matching their keys exactly, and a rerun would have reused every one of them. They
were deleted by hand.

Adding the commit to the key would fix it and also invalidate the whole store on every commit,
which defeats the point. The workable version is a per-solver algorithm tag bumped only when
the numbers can move; until that exists, delete the affected rows by hand after a solver
change and say so in the commit.
"""
key(row) = join((row.solver, row.basis, @sprintf("%.10g",row.We), @sprintf("%.10g",row.Bo),
                 @sprintf("%.10g",row.Oh), row.M, row.K, row.rheology,
                 @sprintf("%.10g",row.t_max), @sprintf("%.10g",row.h_thresh)), "|")

function load_manifest()
    seen = Set{String}()
    isfile(STORE) || return seen
    hdr = nothing
    for ln in eachline(STORE)
        f = split(chomp(ln), ',')
        if hdr === nothing; hdr = f; continue; end
        length(f) < 10 && continue
        push!(seen, join(f[1:10], "|"))
    end
    seen
end

"""Rough trajectory footprint of one run, in bytes: three vectors of length ndof-ish per step."""
function est_bytes(M, K, t_max, dt0)
    nsteps = ceil(Int, t_max/dt0)
    ndof = (M-1)*K
    nsteps * (2*ndof + (M+1) + 6) * 8
end

"""
One case, through the backend API rather than through a solver-specific call.

The `solver` column now holds `label(Backend)`, so a row says which formulation, which contact
closure and which forcing produced it. That is also what keys the store, which closes half of
the hole described above: switching backend can no longer be mistaken for the same run.
"""
function run_one(c)
    r = run_impact(c.backend; We=c.We, Bo=c.Bo, Oh=c.Oh, M=c.M, K=c.K, t_max=c.t_max,
                   h_thresh=c.h_thresh,
                   eta = c.backend.formulation === :variational ? c.eta : nothing,
                   eta_nonvar = c.backend.formulation === :variational ? nothing : c.eta_nonvar)
    ## reduce to scalars HERE, so the trajectory is collectable before the next case
    (solver=c.solver, basis=c.basis, We=c.We, Bo=c.Bo, Oh=c.Oh, M=c.M, K=c.K, rheology=c.rheology,
     t_max=c.t_max, h_thresh=c.h_thresh,
     cor=r.cor, tc=r.tc, cor_internal=r.diag.cor_internal, tc_internal=r.diag.tc_internal,
     min_gap=r.diag.min_gap, n_detected=r.diag.n_detected,
     maxcp=isempty(r.cp) ? 0 : maximum(r.cp), rejects=r.diag.rejects,
     gap_fraction=NaN, lcp_resid=r.diag.lcp_resid, eta_sweeps=r.diag.eta_sweeps,
     wall_s=r.wall, commit=git_commit(), stamp=string(now()))
end

fmt(v) = v isa AbstractFloat ? (isfinite(v) ? @sprintf("%.10g", v) : "NaN") : string(v)

"""
    sweep(cases; budget_mb, force) -> Int

Run every case not already in the store. Returns the number actually computed.
"""
function sweep(cases; budget_mb::Int = 4000, force::Bool = false)
    mkpath(RESULTS)
    if !isfile(STORE)
        open(STORE, "w") do io; println(io, join(COLS, ",")); end
    end
    seen = force ? Set{String}() : load_manifest()
    todo = [c for c in cases if !(key(c) in seen)]
    @printf("%d cases, %d already stored, %d to run\n", length(cases),
            length(cases)-length(todo), length(todo))
    isempty(todo) && return 0
    ## cap concurrency by the memory budget as well as by the threads available
    worst = maximum(est_bytes(c.M, c.K, c.t_max,
                              ImpactParams(We=c.We,Bo=c.Bo,Oh=c.Oh,M=c.M,K=c.K).dt0)
                    for c in todo)
    permit = max(1, min(nthreads(), floor(Int, budget_mb*1e6 / max(worst,1))))
    @printf("worst-case trajectory %.0f MB; %d threads available; running %d at a time\n",
            worst/1e6, nthreads(), permit)
    lk = ReentrantLock(); sem = Base.Semaphore(permit); done = Threads.Atomic{Int}(0)
    @threads for c in todo
        Base.acquire(sem)
        try
            row = run_one(c)
            lock(lk) do
                open(STORE, "a") do io
                    println(io, join((fmt(getfield(row, Symbol(k))) for k in COLS), ","))
                end
            end
            n = Threads.atomic_add!(done, 1) + 1
            @printf("[%3d/%3d] %-7s We=%-8.4g Oh=%-7.4g M=%-3d K=%d  CoR=%-8.4f tc=%-8.4f %.0fs\n",
                    n, length(todo), c.solver, c.We, c.Oh, c.M, c.K, row.cor, row.tc, row.wall_s)
        catch e
            @printf("FAILED  %-7s We=%.4g Oh=%.4g M=%d K=%d : %s\n",
                    c.solver, c.We, c.Oh, c.M, c.K, sprint(showerror, e))
        finally
            Base.release(sem)
        end
    end
    done[]
end

"""A case, with the rheology carried as both a label (for provenance) and a function."""
function case(; backend::Backend = Backend(), We, Bo, Oh, M, K, t_max = 25.0,
              h_thresh = 0.02, rheology = "newtonian", eta = gd -> 1.0, eta_nonvar = nothing,
              basis_kind::Symbol = :legendre)
    (solver=label(backend), basis=String(basis_kind), backend=backend, We=We, Bo=Bo, Oh=Oh, M=M, K=K, t_max=t_max,
     h_thresh=h_thresh, rheology=rheology, eta=eta, eta_nonvar=eta_nonvar,
     basis_kind=basis_kind)
end

## Carreau-Yasuda for the 3000 ppm fluid, from its own Cross fit
const ST = let e0=8.433817577956766, ei=0.0037320997942061666,
               Kc=18.48081673111359, mc=0.7430524574330837,
               R=0.0003, sg=0.0728, rh=1000.0
    tcap = sqrt(rh*R^3/sg)
    (Oh = e0/sqrt(rh*sg*R), Bo = 0.012,
     eta = gd -> carreau(gd; lambda_c=Kc/tcap, a=mc, n=1-mc, eta_inf_ratio=ei/e0),
     label = "carreau_3000ppm")
end

if abspath(PROGRAM_FILE) == @__FILE__
    budget = 4000; force = false
    for (i,a) in enumerate(ARGS)
        a == "--force" && (force = true)
        a == "--budget-mb" && i < length(ARGS) && (budget = parse(Int, ARGS[i+1]))
    end
    cases = Any[]
    ## Newtonian: the Gabbard bands
    for (Oh,Bo) in ((0.0233,0.0526),(0.0373,0.0188),(0.0767,0.0178),(0.2889,0.0158),(0.6849,0.0271))
        for We in exp.(range(log(0.02), log(3.0); length=6)),
            bk in (Backend(contact=:active_set), Backend(), Backend(forcing=:nodal))
            push!(cases, case(backend=bk, We=We, Bo=Bo, Oh=Oh, M=45, K=2))
        end
    end
    ## shear thinning: the 3000 ppm fluid
    for We in exp.(range(log(0.02), log(2.0); length=6)),
        bk in (Backend(contact=:active_set), Backend(), Backend(forcing=:nodal))
        push!(cases, case(backend=bk, We=We, Bo=ST.Bo, Oh=ST.Oh, M=14, K=2,
                          rheology=ST.label, eta=ST.eta))
    end
    n = sweep(cases; budget_mb=budget, force=force)
    @printf("\ncomputed %d new runs; store is %s\n", n, STORE)
end
