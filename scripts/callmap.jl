# A static call map of the package: which file uses which, and what is dead.
#
# The question this answers is "in what order do I read this?", which a directory
# listing cannot. It parses every source file with Julia's own parser rather than
# grepping, so a name inside a string or a comment is never mistaken for a call, and
# short-form definitions (`f(x) = ...`) are found as reliably as long-form ones.
#
# WHAT IT PRODUCES
#   * a reading order, by dependency depth: level 0 depends on nothing else in the
#     package, and each level only needs the ones below it
#   * a file-level graph, drawn as boxes listing their own definitions with arrows for
#     "uses a name defined in"
#   * a dead-code list: defined, exported nowhere, and called nowhere in src, test or
#     scripts
#
# CAVEATS, because a static map cannot see everything. Dispatch through a variable, a
# name reached only via `getproperty`, and anything called from a string are invisible.
# The dead list is therefore a list of CANDIDATES to look at, not a delete list.

using Printf

const ROOT = joinpath(@__DIR__, "..")
const SRC  = joinpath(ROOT, "src")
const USERS = [joinpath(ROOT, "test"), joinpath(ROOT, "scripts")]

"""Name of whatever is being defined by `sig`, or `nothing`."""
function defname(sig)
    sig isa Symbol && return sig
    sig isa Expr || return nothing
    isempty(sig.args)     && return nothing
    sig.head === :where   && return defname(sig.args[1])
    sig.head === :(::)    && return defname(sig.args[1])
    sig.head === :call    && return defname(sig.args[1])
    sig.head === :curly   && return defname(sig.args[1])
    sig.head === :(.)     && return length(sig.args) >= 2 ? defname(sig.args[2]) : nothing
    sig isa QuoteNode     && return defname(sig.value)
    nothing
end

"""Every name defined at any depth of `ex`, and every name referenced."""
function walk!(ex, defs::Set{Symbol}, calls::Dict{Symbol,Int})
    ## Every bare identifier counts as a use, not only the head of a call. A constant
    ## like DEFAULT_M is referenced as a value and never "called", and counting only
    ## call heads reported every such name as dead.
    ex isa Symbol && (calls[ex] = get(calls, ex, 0) + 1; return)
    ex isa Expr || return
    if ex.head === :function || ex.head === :macro
        n = defname(ex.args[1]); n !== nothing && push!(defs, n)
    elseif ex.head === :(=) && ex.args[1] isa Expr &&
           (ex.args[1].head === :call || ex.args[1].head === :where)
        n = defname(ex.args[1]); n !== nothing && push!(defs, n)   # short form
    elseif ex.head === :struct && length(ex.args) >= 2
        n = defname(ex.args[2]); n !== nothing && push!(defs, n)
    elseif ex.head === :const && !isempty(ex.args) && ex.args[1] isa Expr &&
           ex.args[1].head === :(=) && length(ex.args[1].args) >= 1
        n = defname(ex.args[1].args[1]); n !== nothing && push!(defs, n)
    elseif ex.head === :(.) && length(ex.args) >= 2 && ex.args[2] isa QuoteNode &&
           ex.args[2].value isa Symbol
        ## qualified reference: DropSolver.foo
        n = ex.args[2].value; calls[n] = get(calls, n, 0) + 1
    end
    for a in ex.args
        walk!(a, defs, calls)
    end
end

function scan(path)
    src = read(path, String)
    defs, calls = Set{Symbol}(), Dict{Symbol,Int}()
    try
        walk!(Meta.parseall(src), defs, calls)
    catch e
        @warn "could not parse" path exception = e
    end
    (defs, calls)
end

files = sort(filter(f -> endswith(f, ".jl"), readdir(SRC)))
DEF  = Dict{String,Set{Symbol}}()
CALL = Dict{String,Dict{Symbol,Int}}()
for f in files
    DEF[f], CALL[f] = scan(joinpath(SRC, f))
end

## exported names, which are public API and never "dead"
exported = Set{Symbol}()
let ex = Meta.parseall(read(joinpath(SRC, "DropSolver.jl"), String))
    function ex!(e)
        e isa Expr || return
        e.head === :export && (for a in e.args; a isa Symbol && push!(exported, a); end)
        foreach(ex!, e.args)
    end
    ex!(ex)
end

## who owns each name (first definition wins; collisions are reported)
owner = Dict{Symbol,String}()
for f in files, d in DEF[f]
    haskey(owner, d) || (owner[d] = f)
end

## file -> file edges
edges = Dict{String,Set{String}}(f => Set{String}() for f in files)
for f in files, (c, _) in CALL[f]
    o = get(owner, c, nothing)
    o === nothing && continue
    o == f && continue
    push!(edges[f], o)
end

# ---------------------------------------------------------------------------
# Reading order: dependency depth. Level 0 needs nothing else in the package.
# Cycles are broken by reporting them rather than by looping forever.
# ---------------------------------------------------------------------------
## In a function, not a top-level loop: assigning a counter inside a top-level
## `while` puts it in soft scope and Julia treats it as a fresh local.
function layer(files, edges)
    level = Dict{String,Int}()
    remaining = Set(files)
    lv = 0
    while !isempty(remaining) && lv < 50
        ready = [f for f in remaining if all(d -> !(d in remaining), edges[f])]
        isempty(ready) && break                  # a cycle: stop and report
        for f in ready; level[f] = lv; delete!(remaining, f); end
        lv += 1
    end
    for f in remaining; level[f] = lv; end       # cyclic group, all on one level
    (level, remaining)
end
level, remaining = layer(files, edges)

flushln(a...) = (println(a...); flush(stdout))
flushln("== reading order (level 0 depends on nothing else in the package) ==")
for l in 0:maximum(values(level))
    fs = sort([f for f in files if level[f] == l])
    isempty(fs) && continue
    flushln("\nlevel $l")
    for f in fs
        deps = sort(collect(edges[f]))
        @printf("  %-26s %3d defs   uses: %s\n", f, length(DEF[f]),
                isempty(deps) ? "-" : join(deps, ", "))
    end
end
!isempty(remaining) && flushln("\nCYCLE among: ", join(sort(collect(remaining)), ", "))

# ---------------------------------------------------------------------------
# Dead-code candidates: defined, not exported, not called anywhere.
# ---------------------------------------------------------------------------
## A name's own definition contributes one occurrence of it, so "referenced
## somewhere" means the total count exceeds the number of places defining it.
total = Dict{Symbol,Int}()
for f in files, (n, k) in CALL[f]; total[n] = get(total, n, 0) + k; end
for dir in USERS, f in readdir(dir)
    endswith(f, ".jl") || continue
    _, c = scan(joinpath(dir, f))
    for (n, k) in c; total[n] = get(total, n, 0) + k; end
end
ndef = Dict{Symbol,Int}()
for f in files, d in DEF[f]; ndef[d] = get(ndef, d, 0) + 1; end
used = Set(n for (n, k) in total if k > get(ndef, n, 0))

function report_dead(files, DEF, used, exported)
    flushln("\n== dead-code candidates (defined, not exported, never called) ==")
    ndead = 0
    for f in files
        d = sort([string(n) for n in DEF[f] if !(n in used) && !(n in exported)])
        isempty(d) && continue
        ndead += length(d)
        @printf("  %-26s %s\n", f, join(d, ", "))
    end
    ndead == 0 && flushln("  none")
    ndead
end
ndead = report_dead(files, DEF, used, exported)
@printf("\n%d files, %d definitions, %d dead candidates\n",
        length(files), sum(length(DEF[f]) for f in files), ndead)

## Files nothing uses. DropSolver.jl is excluded as a consumer: it includes and
## re-exports everything, so counting it would mask a file that nothing actually calls.
flushln("\n== used by nothing except the module itself ==")
others = filter(!=("DropSolver.jl"), files)
for f in others
    any(g -> f in edges[g], others) && continue
    @printf("  %-26s <- entry point, or dead\n", f)
end

## the machine-readable form, for the drawing step
open(joinpath(ROOT, "outputs", "csv", "callmap.csv"), "w") do io
    println(io, "from,to")
    for f in files, t in sort(collect(edges[f])); println(io, "$f,$t"); end
end
open(joinpath(ROOT, "outputs", "csv", "callmap_nodes.csv"), "w") do io
    println(io, "file,level,ndefs,defs")
    for f in files
        println(io, join((f, level[f], length(DEF[f]),
                          join(sort(string.(collect(DEF[f]))), " ")), ","))
    end
end
flushln("\nwrote outputs/csv/callmap.csv and outputs/csv/callmap_nodes.csv")
