# The claim ledger.
#
# The Model summary of "Shear-Thinning Drops" is the statement of the model, and
# the standard it is held to is that every equation in it is discharged by a
# machine-checked derivation rather than by prose. This test is what makes that
# enforceable: it is not itself a proof of anything, it is the bookkeeping that
# prevents an unproved equation from entering the summary unnoticed.
#
# It enforces three things:
#
#   1. COVERAGE.   Every display equation in the Model summary is claimed by a
#                  ledger entry. Adding an equation without registering it fails
#                  here, so the summary cannot grow silently.
#   2. ANCHORING.  Every ledger entry's anchor text is still present in the
#                  summary. Editing an equation breaks its anchor, so a changed
#                  equation cannot keep an old proof's endorsement.
#   3. DISCHARGE.  Every entry marked `:proved` points at a check that carries
#                  that claim's OWN id as a `## CLAIM: <id>` tag on the assertion
#                  line. An earlier version searched only for the assertion's
#                  NAME, which any neighbouring check satisfies -- and that hole
#                  passed five claims that were never verified.
#
# Entries carry one of three statuses, and the distinction is the point:
#
#   :proved  -- a symbolic or exact-numerical check exists and runs in CI.
#   :axiom   -- a modelling postulate, not a consequence of anything above it.
#               These can never become :proved and are listed so that the set of
#               things taken on faith is explicit and finite.
#   :open    -- believed true, not yet checked. This count is allowed to fall and
#               is not allowed to rise (see the ratchet at the end).

using Test

const DERIV_DIR = joinpath(@__DIR__, "..", "derivations")
const MODEL_PAGE = joinpath(DERIV_DIR, "generalized_newtonian_hierarchy_derivation.jl")

"""Return the Model summary section of the model page, as one string with the
Literate comment prefix stripped."""
function summary_text(path::AbstractString)
    lines = readlines(path)
    i = findfirst(l -> startswith(l, "# ## Model summary"), lines)
    i === nothing && error("the Model summary section has been renamed or removed")
    body = lines[i:end]
    join([replace(l, r"^#\s?" => "") for l in body], "\n")
end

"""Every ```math ... ``` display block in `text`."""
function math_blocks(text::AbstractString)
    blocks = String[]
    inblock = false
    buf = String[]
    for l in split(text, "\n")
        if !inblock && startswith(strip(l), "```math")
            inblock = true; empty!(buf)
        elseif inblock && strip(l) == "```"
            inblock = false; push!(blocks, join(buf, "\n"))
        elseif inblock
            push!(buf, l)
        end
    end
    inblock && error("unterminated ```math block in the Model summary")
    blocks
end

# id, status, anchor (must appear verbatim in the summary), discharged-by
const LEDGER = [
    (id = "SUM-GROUPS", status = :proved, anchor = raw"\mathrm{Oh}=\frac{\eta_0}",
     by = "ASSERTION 2c"),
    (id = "SUM-ZETA",   status = :axiom,  anchor = raw"\zeta(\theta,t)=\sum_{l\ge2}\zeta_l(t)P_l(\mu)",
     by = ""),  # choice of basis and of truncation floor l >= 2
    (id = "SUM-PI",     status = :axiom,  anchor = raw"p_c(\theta,t)=\sum_{l\ge0}p_{c,l}(t)P_l(\mu)",
     by = ""),
    # The angular reduction of the operator is checked, but the right-hand side
    # -- that the source is exactly Oh div(div(2 eta e)) -- is not.
    (id = "SUM-PLAP",   status = :open,   anchor = raw"\mathcal L_l[p_l] \;=\; \mathrm{Oh}\,S_l(x,t)",
     by = ""),
    (id = "SUM-PHARM",  status = :proved, anchor = raw"p = \sum_l c_l(t)\,x^lP_l(\mu)",
     by = "ASSERTION 2c"),
    (id = "SUM-GAP",    status = :proved, anchor = raw"h(\theta,t)=\mu\,\bigl[1+\zeta(\theta,t)\bigr]+z(t)",
     by = "ASSERTION 2c"),
    (id = "SUM-STREAM", status = :proved, anchor = raw"u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}",
     by = "ASSERTION 2c"),
    (id = "SUM-FW",     status = :proved, anchor = raw"u_{\theta,l}=\frac{\psi_l'}{x\,l(l+1)}",
     by = "ASSERTION 2b"),
    # Discharged end to end: the 3x3 determinant built from this equation plus
    # BC1-BC3 must vanish at roots the solver finds from Reid's characteristic
    # equation independently, and sigma = Oh q^2 enters through it, so a wrong
    # Oh prefactor breaks it.
    (id = "SUM-INT",    status = :proved, anchor = raw"\partial_t\,\mathcal D_l[\psi_l]",
     by = "ASSERTION 3e"),
    (id = "SUM-BC",     status = :proved, anchor = raw"\mathcal T[\psi_l]\big|_{x=1}=0",
     by = "ASSERTION 2b"),
    # ASSERTION 3d proves the SOURCE vanishes iff eta is constant (that is
    # SUM-PSRC below). It does not prove the equation: neither the Oh prefactor
    # nor the form of the source is checked.
    (id = "SUM-PRESS",  status = :open,   anchor = raw"\nabla^2 p = \mathrm{Oh}\,\nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr)",
     by = ""),
    (id = "SUM-PSRC",   status = :proved, anchor = raw"\;\equiv\; \mathrm{Oh}\,S(x,\mu,t)",
     by = "ASSERTION 3d"),
    (id = "SUM-FILM",   status = :axiom,  anchor = raw"h\,p_c = 0",
     by = ""),  # lubrication limit: the film transmits stress with no dynamics of its own
    # Derived from the traction jump, and the sign on the curvature term is fixed
    # by the base state: p = 2, div n = 2 at rest, and only one combination gives
    # zero. The form previously on the page failed that test by -4.
    (id = "SUM-NORMAL", status = :proved, anchor = raw"\bigl[-p+2\eta\,e_{rr}\bigr]_{x=1}",
     by = "ASSERTION 3e"),
    (id = "SUM-CURV",   status = :proved, anchor = raw"(l-1)(l+2)\zeta_lP_l(\mu)",
     by = "ASSERTION 2c"),
    # The force relation and the mass are checked (SUM-FORCE); the assembled
    # equation of motion is not, so a wrong sign on the gravity term would pass.
    (id = "SUM-COM",    status = :open,   anchor = raw"\dot v = -\mathrm{Bo} - p_{c,1}",
     by = ""),
    (id = "SUM-FORCE",  status = :proved, anchor = raw"= -\frac{4\pi}{3}p_{c,1}",
     by = "ASSERTION 2c"),
    (id = "SUM-NOPULL", status = :axiom,  anchor = raw"p_c\ge0",
     by = ""),  # a gas film cannot sustain tension
    (id = "SUM-RHEO",   status = :proved, anchor = raw"\dot\gamma=\sqrt{2\,\bm e\!:\!\bm e}",
     by = "ASSERTION 2c"),
]

# The number of :open entries may fall but never rise. Lower this when a claim is
# discharged; raising it requires deleting this comment and explaining why.
#
# It was briefly recorded as zero. That was wrong: the discharge test then only
# searched for an assertion's NAME, so five entries counted as proved on the
# strength of checks that verified something else nearby. With the tag-based test
# the count is honest, and it is not zero.
#
# Reduce this only by tagging a check that would falsify the equation in
# question. Do not reduce it by re-pointing an anchor or by citing a neighbouring
# assertion -- that is precisely what produced the false zero.
const OPEN_BUDGET = 3

@testset "claim ledger" begin
    text = summary_text(MODEL_PAGE)
    blocks = math_blocks(text)

    @testset "the summary has equations to check" begin
        @test length(blocks) >= 10
    end

    @testset "anchoring: every ledger anchor is still in the summary" begin
        for e in LEDGER
            @test occursin(e.anchor, text) ||
                  error("ledger entry $(e.id) is orphaned: its anchor no longer " *
                        "appears in the Model summary. If the equation changed, " *
                        "its proof must be revisited, not its anchor patched.")
        end
    end

    @testset "coverage: every summary equation is claimed" begin
        for (k, b) in enumerate(blocks)
            covering = [e.id for e in LEDGER if occursin(e.anchor, b)]
            @test !isempty(covering) ||
                  error("display equation $k of the Model summary is not claimed " *
                        "by any ledger entry:\n$(first(split(b, '\n')))\n" *
                        "Register it in LEDGER with a status.")
        end
    end

    @testset "discharge: every :proved claim is tagged at an actual check" begin
        # A `:proved` entry must point at a check that carries its OWN id as a
        # `## CLAIM: <id>` tag on the assertion line. An earlier version of this
        # test merely searched for the assertion's NAME anywhere in the sources,
        # which is satisfied by any check in the same file -- and that hole was
        # exploited: SUM-NORMAL was marked proved by an assertion that verified
        # the curvature term beside it and never touched the normal-stress
        # balance. Proof by adjacency passes a name search; it cannot pass this.
        tagged = Dict{String,Vector{String}}()
        for f in readdir(DERIV_DIR)
            endswith(f, ".jl") || continue
            for (i, line) in enumerate(readlines(joinpath(DERIV_DIR, f)))
                m = match(r"##\s*CLAIM:\s*([A-Z0-9\-]+)", line)
                m === nothing && continue
                occursin("@assert", line) ||
                    error("$f:$i tags CLAIM: $(m.captures[1]) but is not an " *
                          "assertion line. A tag must sit on the check itself.")
                push!(get!(tagged, m.captures[1], String[]), "$f:$i")
            end
        end
        for e in LEDGER
            if e.status === :proved
                @test haskey(tagged, e.id) ||
                      error("ledger entry $(e.id) is marked :proved but no " *
                            "assertion carries `## CLAIM: $(e.id)`. Either tag " *
                            "the check that falsifies it, or set it :open.")
            else
                @test !haskey(tagged, e.id) ||
                      error("$(e.id) is tagged at a check but is not marked " *
                            ":proved -- the ledger understates what is verified.")
            end
        end
        @info "claim tags found" tags = sort(collect(keys(tagged)))
    end

    @testset "the open set is not growing" begin
        n_open = count(e -> e.status === :open, LEDGER)
        @test n_open <= OPEN_BUDGET
        n_axiom = count(e -> e.status === :axiom, LEDGER)
        @info "claim ledger: $(count(e -> e.status === :proved, LEDGER)) proved, " *
              "$n_axiom axioms, $n_open open (budget $OPEN_BUDGET)"
    end
end
