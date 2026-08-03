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
#   3. DISCHARGE.  Every entry marked `:proved` names an assertion label that
#                  actually appears in the derivation sources.
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
    (id = "SUM-ZETA",   status = :axiom,  anchor = raw"\zeta(\theta,t)=\sum_{l\ge2}A_l(t)P_l(\mu)",
     by = ""),  # choice of basis and of truncation floor l >= 2
    (id = "SUM-PI",     status = :axiom,  anchor = raw"\Pi(\theta,t)=\sum_{n\ge0}B_n(t)P_n(\mu)",
     by = ""),
    (id = "SUM-GAP",    status = :proved, anchor = raw"h(\theta,t)=\mu\,\bigl[1+\zeta(\theta,t)\bigr]+z(t)",
     by = "ASSERTION 2c"),
    (id = "SUM-STREAM", status = :proved, anchor = raw"u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}",
     by = "ASSERTION 2c"),
    (id = "SUM-FW",     status = :proved, anchor = raw"W_l=\frac{U_l'}{x\,l(l+1)}",
     by = "ASSERTION 2b"),
    (id = "SUM-INT",    status = :proved, anchor = raw"\partial_t\,\mathcal D_l[U_l]",
     by = "ASSERTION 2c"),
    (id = "SUM-BC",     status = :proved, anchor = raw"\mathcal L_2[U_l]\big|_{x=1}=0",
     by = "ASSERTION 2b"),
    (id = "SUM-PRESS",  status = :proved, anchor = raw"\nabla^2 p = \mathrm{Oh}\,\nabla\cdot\bigl(\nabla\cdot(2\eta\bm e)\bigr)",
     by = "ASSERTION 3d"),
    (id = "SUM-NORMAL", status = :proved, anchor = raw"\bigl[-p+2\eta\,e_{rr}\bigr]_{x=1}",
     by = "ASSERTION 2c"),
    (id = "SUM-CURV",   status = :proved, anchor = raw"(l-1)(l+2)A_lP_l(\mu)",
     by = "ASSERTION 2c"),
    (id = "SUM-COM",    status = :proved, anchor = raw"\dot v = -\mathrm{Bo} - B_1",
     by = "ASSERTION 2c"),
    (id = "SUM-FORCE",  status = :proved, anchor = raw"= -\frac{4\pi}{3}B_1",
     by = "ASSERTION 2c"),
    (id = "SUM-SIGN",   status = :axiom,  anchor = raw"h\,\Pi = 0",
     by = ""),  # the substrate is rigid and cannot pull: a postulate
    (id = "SUM-RHEO",   status = :proved, anchor = raw"\dot\gamma=\sqrt{2\,\bm e\!:\!\bm e}",
     by = "ASSERTION 2c"),
]

# The number of :open entries may fall but never rise. Lower this when a claim is
# discharged; raising it requires deleting this comment and explaining why.
#
# It is now zero: every equation in the Model summary is discharged by a check
# that runs in CI, and the only things taken on faith are the three :axiom
# entries -- the choice of Legendre basis for the surface and the pressure, and
# the postulate that a rigid substrate cannot pull.
const OPEN_BUDGET = 0

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

    @testset "discharge: :proved entries name a real assertion" begin
        sources = join([read(joinpath(DERIV_DIR, f), String)
                        for f in readdir(DERIV_DIR) if endswith(f, ".jl")], "\n")
        for e in LEDGER
            if e.status === :proved
                @test !isempty(e.by)
                @test occursin(e.by, sources) ||
                      error("ledger entry $(e.id) claims to be discharged by " *
                            "\"$(e.by)\", which appears in no derivation source.")
            else
                @test isempty(e.by)  # only :proved entries cite a check
            end
        end
    end

    @testset "the open set is not growing" begin
        n_open = count(e -> e.status === :open, LEDGER)
        @test n_open <= OPEN_BUDGET
        n_axiom = count(e -> e.status === :axiom, LEDGER)
        @info "claim ledger: $(count(e -> e.status === :proved, LEDGER)) proved, " *
              "$n_axiom axioms, $n_open open (budget $OPEN_BUDGET)"
    end
end
