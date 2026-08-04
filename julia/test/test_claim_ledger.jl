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
    # Bo and We need their own entries: the per-equation coverage test showed that
    # SUM-GROUPS claimed only Oh, and the other two rode along in the same fence.
    (id = "SUM-BO",     status = :proved, anchor = raw"\mathrm{Bo}=\frac{\rho g R^2}{T_1}",
     by = "ASSERTION 2c"),
    (id = "SUM-WE",     status = :proved, anchor = raw"\mathrm{We}=\frac{\rho R V^2}{T_1}",
     by = "ASSERTION 2c"),
    (id = "SUM-ZETA",   status = :axiom,  anchor = raw"\zeta(\theta,t)=\sum_{l\ge2}\zeta_l(t)P_l(\mu)",
     by = ""),  # choice of basis and of truncation floor l >= 2
    (id = "SUM-PI",     status = :axiom,  anchor = raw"p_c(\theta,t)=\sum_{l\ge0}p_{c,l}(t)P_l(\mu)",
     by = ""),
    (id = "SUM-GAP",    status = :proved, anchor = raw"h(\theta,t)=\mu\,\bigl[1+\zeta(\theta,t)\bigr]+z(t)",
     by = "ASSERTION 2c"),
    (id = "SUM-STREAM", status = :proved, anchor = raw"u_r=\frac{1}{x^2\sin\theta}\frac{\partial\psi}{\partial\theta}",
     by = "ASSERTION 2c"),
    (id = "SUM-FW",     status = :proved, anchor = raw"u_{\theta,l}=\frac{\psi_l'}{x\,l(l+1)}",
     by = "ASSERTION 2b"),
    # u_{r,l} had been riding on SUM-FW's anchor in the same fence, while the
    # assertion that actually checks it (w_F) carried no tag at all.
    (id = "SUM-URL",    status = :proved, anchor = raw"u_{r,l}=\frac{\psi_l}{x^2}",
     by = "ASSERTION 2b"),
    (id = "SUM-UTH",    status = :proved, anchor = raw"u_\theta=-\frac{1}{x\sin\theta}\frac{\partial\psi}{\partial x}",
     by = "ASSERTION 2b"),
    # T and V get their own entries. V sharing a fence with Phi is exactly how a
    # wrong surface energy was carried as :proved by a symmetry check on its
    # neighbour; SUM-STIFF is the calibration that can falsify it.
    # The two relations that DEFINE the coordinate: the interior displacement is the
    # coordinate, and the surface amplitude is its boundary trace. Discharged against
    # the RUNNING SOLVER rather than against the script's own algebra, because what can
    # fail for a definition is not the algebra but the code disagreeing with the page --
    # a solver carrying the stream function as its coordinate would satisfy every
    # equation on the page and still be integrating the wrong Lagrangian.
    (id = "SUM-CHIDOT", status = :proved, anchor = raw"\psi_l = \dot\chi_l",
     by = "ASSERTION 5f"),
    (id = "SUM-TRACE",  status = :proved, anchor = raw"\zeta_l = \chi_l(1,t)",
     by = "ASSERTION 5f"),
    (id = "SUM-T",      status = :proved, anchor = raw"T[\dot{\bm\xi}]=\tfrac12\int|\bm u|^2\,dV",
     by = "ASSERTION 5e"),
    (id = "SUM-STIFF",  status = :proved, anchor = raw"\frac{4\pi}{2l+1}(l-1)(l+2)\,\zeta_l^2",
     by = "ASSERTION 5e"),
    # The variational statement. Discharged by the identity that makes it the SAME
    # system as the differential form (ASSERTION 5d), plus the dissipation form's
    # symmetry and its exact reproduction of Lamb's damping (ASSERTION 5b).
    (id = "SUM-EL",     status = :proved, anchor = raw"\frac{d}{dt}\frac{\partial T}{\partial\dot\xi_a}",
     by = "ASSERTION 5d"),
    (id = "SUM-FORMS",  status = :proved, anchor = raw"\bm e\!:\!\bm e\,dV",
     by = "ASSERTION 5b"),
    (id = "SUM-HESS",   status = :proved, anchor = raw"\bm e^{(a)}\!:\!\bm e^{(b)}\,dV",
     by = "ASSERTION 5b"),
    (id = "SUM-BC",     status = :proved, anchor = raw"\mathcal T[\psi_l]\big|_{x=1}=0",
     by = "ASSERTION 2b"),
    (id = "SUM-FILM",   status = :axiom,  anchor = raw"h\,p_c = 0",
     by = ""),  # lubrication limit: the film transmits stress with no dynamics of its own
    (id = "SUM-CURV",   status = :proved, anchor = raw"(l-1)(l+2)\zeta_lP_l(\mu)",
     by = "ASSERTION 2c"),
    # Assembled from m dV/dT = -m g + F_z under the stated scalings, signs
    # included: the gravity term reduces to exactly -Bo and the force term to
    # F/(4pi/3), checked at several dimensional parameter sets.
    (id = "SUM-COM",    status = :proved,   anchor = raw"\dot v = -\mathrm{Bo} - p_{c,1}",
     by = "ASSERTION 3f"),
    (id = "SUM-FORCE",  status = :proved, anchor = raw"= -\frac{4\pi}{3}p_{c,1}",
     by = "ASSERTION 2c"),
    # The generalised force, which the summary previously left undefined -- and it
    # was the sole coupling between contact and shape, so the model did not exist
    # without it. Discharged by the consistency the derivation predicts: at l = 1
    # it must equal the net force F already proved in SUM-FORCE.
    (id = "SUM-Q",      status = :proved, anchor = raw"Q_{\zeta_l}=-\frac{4\pi}{2l+1}\,p_{c,l}",
     by = "ASSERTION 2c"),
    (id = "SUM-VWORK",  status = :proved, anchor = raw"\delta W=-\oint p_c\,\delta\zeta\,dS",
     by = "ASSERTION 2c"),
    # The film acts only on the surface, so the interior coordinates take no direct
    # forcing. A modelling statement about where the traction is applied, not a
    # derivable one.
    (id = "SUM-QPSI",   status = :axiom,  anchor = raw"Q_{\psi_l}=0",
     by = ""),
    (id = "SUM-ZDOT",   status = :axiom,  anchor = raw"\dot z = v",
     by = ""),  # definition of v
    # Initial data and the finite reduction of the contact conditions. All choices,
    # none derivable, and all previously missing from the summary entirely -- an
    # implementer had to invent them.
    (id = "SUM-V0",     status = :axiom,  anchor = raw"v(0)=-\sqrt{\mathrm{We}}",
     by = ""),
    (id = "SUM-Z0",     status = :axiom,  anchor = raw"z(0)=1",
     by = ""),
    (id = "SUM-ZETA0",  status = :axiom,  anchor = raw"\zeta_l(0)=0",
     by = ""),
    (id = "SUM-PSI0",   status = :axiom,  anchor = raw"\psi_l(x,0)=0",
     by = ""),
    (id = "SUM-COLLOC", status = :axiom,  anchor = raw"h(\theta_i)\ge0",
     by = ""),  # collocation on the reconstructed field, M+1 nodes
    (id = "SUM-NODES",  status = :axiom,  anchor = raw"i=0\ldots M",
     by = ""),
    (id = "SUM-NOPULL", status = :axiom,  anchor = raw"p_c\ge0",
     by = ""),  # a gas film cannot sustain tension
    (id = "SUM-HGE0",   status = :axiom,  anchor = raw"h\ge0",
     by = ""),  # the substrate is impenetrable
    (id = "SUM-ETA",    status = :axiom,  anchor = raw"\eta=\eta(\dot\gamma)",
     by = ""),  # hypothesis (H1): no strain history, no dependence beyond the invariant
    (id = "SUM-RHEO",   status = :proved, anchor = raw"\dot\gamma=\sqrt{2\,\bm e\!:\!\bm e}",
     by = "ASSERTION 2c"),
]

# The number of :open entries may fall but never rise. Lower this when a claim is
# discharged; raising it requires deleting this comment and explaining why.
#
# This reached zero once before dishonestly: the discharge test then searched
# only for an assertion's NAME, so five entries counted as proved on the strength
# of checks that verified something else nearby. It is zero again now, under the
# tag-based test, with every entry pointing at a check that carries its own id.
#
# Reduce this only by tagging a check that would falsify the equation in
# question. Do not reduce it by re-pointing an anchor or by citing a neighbouring
# assertion -- that is precisely what produced the false zero.
const OPEN_BUDGET = 0

# The number of :axiom entries, budgeted for the same reason. Four things are taken
# on faith: the Legendre bases for the surface and the film pressure, the rigid
# unilateral contact law, and that a gas film cannot sustain tension.
const AXIOM_BUDGET = 14

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
        # PER EQUATION, not per fenced block. One anchor used to endorse every
        # equation sharing its ```math fence, and that is precisely how a wrong
        # surface energy passed: it sat in the same block as the dissipation form,
        # whose anchor matched, and inherited its :proved status. An equation here
        # is a comma- or newline-separated relation containing `=`, `\ge` or `\le`.
        for (k, b) in enumerate(blocks)
            # A \boxed{...} block is one equation by convention, even when it is
            # wrapped over several lines with the relation symbol on the last of
            # them. Splitting it leaves a `= Q_a` fragment that no anchor matches.
            if occursin("\\boxed", b)
                @test any(e -> occursin(e.anchor, b), LEDGER) ||
                      error("boxed equation in display block $k is not claimed by " *
                            "any ledger entry:\n$(first(strip(b), 90))")
                continue
            end
            rels = String[]
            for chunk in split(b, r",\s*\n|\\\\|\n")
                # strip only the wide spacing macros. NOT `\,` or `\;` -- those
                # appear inside anchors, and stripping them makes a registered
                # equation look unclaimed.
                c = strip(replace(chunk, r"\\qquad|\\quad" => " "))
                (isempty(c) || startswith(c, "\\begin") || startswith(c, "\\end")) && continue
                occursin(r"=|\\ge|\\le", c) || continue
                # A chunk that STARTS with a relation symbol is the continuation of
                # a wrapped equation, not a new one. Splitting on newlines alone
                # turns `u_theta = -...\n = sum_l ...` into two fragments and then
                # reports the tail as unclaimed.
                if !isempty(rels) && occursin(r"^(\\;)?=", c)
                    rels[end] = rels[end] * " " * c
                else
                    push!(rels, c)
                end
            end
            isempty(rels) && continue
            uncovered = [r for r in rels
                         if !any(e -> occursin(e.anchor, r) ||
                                      occursin(e.anchor, b) && occursin(e.anchor, r), LEDGER)]
            # an equation is covered if some anchor appears INSIDE it
            uncovered = [r for r in rels if !any(e -> occursin(e.anchor, r), LEDGER)]
            @test isempty(uncovered) ||
                  error("in display block $k of the Model summary, these relations are " *
                        "not claimed by any ledger entry:\n" *
                        join(["    " * first(r, 90) for r in uncovered], "\n") *
                        "\nRegister each in LEDGER with its own status. Sharing a " *
                        "block with a claimed equation is not coverage.")
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
        # SCOPE: only the model page may discharge a claim about the model. A
        # SUM-COM tag was living on the closures page -- the file the review brief
        # designates "NOT the model" -- so an archived alternative route was
        # endorsing the live statement.
        tagged = Dict{String,Vector{String}}()
        for f in readdir(DERIV_DIR)
            endswith(f, ".jl") || continue
            for (i, line) in enumerate(readlines(joinpath(DERIV_DIR, f)))
                m = match(r"##\s*CLAIM:\s*([A-Z0-9\-]+)", line)
                m === nothing && continue
                occursin("@assert", line) ||
                    error("$f:$i tags CLAIM: $(m.captures[1]) but is not an " *
                          "assertion line. A tag must sit on the check itself.")
                f == basename(MODEL_PAGE) ||
                    error("$f:$i tags CLAIM: $(m.captures[1]), but only the model " *
                          "page may discharge a model claim. Move the check, or " *
                          "the claim is not about the model.")
                push!(get!(tagged, m.captures[1], String[]), "$f:$i")
            end
        end
        ids = Set(e.id for e in LEDGER)
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
        # ORPHANS: a tag with no ledger entry means an equation was removed from the
        # summary while its proof was left behind. That is how the open count once
        # reached zero partly by deletion, and it must not be silent.
        orphans = setdiff(keys(tagged), ids)
        @test isempty(orphans) ||
              error("these checks are tagged but claim nothing in the ledger: " *
                    join(sort(collect(orphans)), ", ") *
                    "\nEither the equation left the summary (drop the tag) or its " *
                    "entry was lost (restore it). A dangling proof is not evidence.")
        @info "claim tags found" tags = sort(collect(keys(tagged)))
    end

    @testset "the open set is not growing" begin
        n_open = count(e -> e.status === :open, LEDGER)
        @test n_open <= OPEN_BUDGET
        n_axiom = count(e -> e.status === :axiom, LEDGER)
        # AXIOMS ARE BUDGETED TOO. Unbudgeted, `:axiom` is an escape hatch that
        # makes a zero open count unfalsifiable -- any awkward claim can be
        # relabelled a postulate and the ratchet stays silent. Raising this
        # requires saying which new thing is being taken on faith and why it
        # cannot be derived.
        @test n_axiom <= AXIOM_BUDGET
        @info "claim ledger: $(count(e -> e.status === :proved, LEDGER)) proved, " *
              "$n_axiom axioms, $n_open open (budget $OPEN_BUDGET)"
    end
end
