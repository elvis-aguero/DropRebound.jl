using Latexify

"""
    _apply(s, pat, transform) -> String

Like `replace(s, pat => transform)`, but `transform` receives a real
`RegexMatch` (capture groups accessible via `m[1]`, `m[2]`, ...) rather than
the plain matched substring `Base.replace` itself passes to a function
replacement -- and the returned string is spliced in verbatim, with no
further escape-sequence reinterpretation (unlike `SubstitutionString`,
which re-parses backslashes in both the pattern AND in any already-escaped
LaTeX command being substituted in, corrupting it).
"""
function _apply(s::AbstractString, pat::Regex, transform)
    io = IOBuffer()
    lastpos = firstindex(s)
    for m in eachmatch(pat, s)
        write(io, s[lastpos:prevind(s, m.offset)])
        write(io, transform(m))
        lastpos = m.offset + ncodeunits(m.match)
    end
    write(io, s[lastpos:end])
    String(take!(io))
end

"""
    pretty_latex(expr) -> String

Latexify `expr` (a Symbolics.jl expression), then clean up the common case
Latexify itself doesn't handle: a variable named `greekname_subscript` (e.g.
`mu_0`, `lambda_c`, `eps_ST`) renders as a literal, unstyled `\\mathtt{...}`
block instead of `\\mu_{0}`, `\\lambda_{c}`, `\\varepsilon_{ST}`. This is
purely a DISPLAY transform on the string Latexify already produced from the
real Symbolics object -- it never touches the underlying expression, so the
math itself still can't drift from what's being asserted against; only how
its variable names are typeset changes.
"""
function pretty_latex(expr)
    s = latexify(expr)
    for (ascii, cmd) in GREEK_NAMES
        sub_pat = Regex(raw"\\mathtt\{" * ascii * raw"\\_([a-zA-Z0-9]+)\}")
        s = _apply(s, sub_pat, m -> cmd * "_{" * m[1] * "}")
        bare_pat = Regex(raw"\\mathtt\{" * ascii * raw"\}")
        s = _apply(s, bare_pat, m -> cmd)
    end
    # any other name_sub -> \mathrm{name}_{sub} (upright roman instead of
    # unstyled typewriter, still readable, no information lost)
    fallback_pat = Regex(raw"\\mathtt\{([a-zA-Z]+)\\_([a-zA-Z0-9]+)\}")
    s = _apply(s, fallback_pat, m -> "\\mathrm{" * m[1] * "}_{" * m[2] * "}")
    s
end

const GREEK_NAMES = [
    "alpha" => "\\alpha", "beta" => "\\beta", "gamma" => "\\gamma", "delta" => "\\delta",
    "eps" => "\\varepsilon", "epsilon" => "\\varepsilon", "zeta" => "\\zeta", "eta" => "\\eta",
    "theta" => "\\theta", "lambda" => "\\lambda", "mu" => "\\mu", "nu" => "\\nu",
    "xi" => "\\xi", "pi" => "\\pi", "rho" => "\\rho", "sigma" => "\\sigma", "tau" => "\\tau",
    "phi" => "\\phi", "chi" => "\\chi", "psi" => "\\psi", "omega" => "\\omega",
    "Oh" => "\\mathrm{Oh}", "gammadot" => "\\dot{\\gamma}",
]
