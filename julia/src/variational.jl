# Variational assembly of the shear-thinning drop model.
#
# This implements the model as stated in the "Model summary" of the derivation page
# *Shear-Thinning Drops*: a damped Lagrangian system in the interior displacement
# coordinates, with the surface amplitude determined by the kinematic condition.
#
# WHY THE COORDINATES ARE DISPLACEMENTS, not the stream function itself. The
# summary's state is (zeta_l, psi_l), but psi_l is a VELOCITY-like quantity -- the
# kinetic energy is quadratic in psi, and the kinematic condition reads
# psi_l(1) = d(zeta_l)/dt. A Lagrangian needs (q, qdot), so the coordinate is the
# interior DISPLACEMENT chi_l with
#
#     psi_l = d(chi_l)/dt ,      zeta_l = chi_l(1) ,
#
# the second following from integrating the kinematic condition. Everything is then
# in chi alone, and the surface amplitude is its boundary trace -- which is also the
# resolution of the coordinate-count question: there is no independent equation for
# zeta_l, because zeta_l is not an independent coordinate.
#
# With chi_l(x,t) = sum_k a_k(t) phi_k(x) on a Ritz basis satisfying regularity,
#
#     M ddot(a) + K dot(a) + G a = Q ,
#
#     M_ab = int u^(a) . u^(b) dV                     (kinetic / added mass)
#     K_ab = Oh int 2 eta e^(a):e^(b) dV              (dissipation)
#     G_ab = (4pi/(2l+1)) (l-1)(l+2) phi_a(1) phi_b(1)  (surface energy)
#
# G is RANK ONE for a single mode, because the surface energy sees the interior only
# through the boundary trace. That is a structural feature worth knowing: it means
# the stiffness cannot resolve interior structure, and only the trace phi_a(1)
# matters to it.
#
# All three are assembled by ONE derivative of the velocity -- no stress divergence,
# no curl, no fourth-order operator -- which is the point of the variational form.
# Strain components are evaluated from their analytic modal expressions, so there is
# no finite-difference error anywhere in the assembly.

using LinearAlgebra
# Gauss-Legendre nodes come from the package's own `integrals.jl` rather than a new
# dependency -- it is already used by the contact quadrature.

"""
    RitzBasis(l, K, kind = :legendre)

Radial trial functions for mode `l`, `k = 1..K`.

    :monomial   phi_k(x) = x^(l+1+2(k-1))
    :legendre   phi_k(x) = x^(l+1) * P_{k-1}(2x^2 - 1)

THE TWO SPAN THE SAME SPACE. Both are `x^(l+1)` times a polynomial of degree `k-1` in `x^2`,
so every statement about what the model can represent is identical, and both give the same
answer wherever both are numerically sound. What differs is conditioning, and by a lot.

The factor `x^(l+1)` builds in regularity at the origin, which is what the model requires of
`psi_l`, so no boundary condition has to be imposed there. `phi_1` is the potential-flow
profile in either kind -- `P_0 = 1` -- so `K = 1` reproduces Lamb exactly and larger `K` adds
interior structure. Every member equals one at the surface, since `P_j(1) = 1`, so the trace
vector is all ones either way and nothing downstream changes.

WHY THE DEFAULT IS `:legendre`. The monomials are a Vandermonde family: the per-mode mass
matrix loses about one and a half digits per added function at `l = 2` and two and a half at
`l = 90`, so double precision is exhausted at `K` between 3 and 8 depending on the mode. Past
that the solver does not degrade gracefully -- at `M = 45, K = 6` it returns a restitution of
102, and on a shear-thinning fluid `K >= 4` gives NaN or a plausible-looking 0.998. The
radial resolution needed to resolve the vortical field grows like `|q| ~ l^{3/4}/sqrt(Oh)`,
so the requirement and the monomial ceiling collide well inside the range of a generic
shear-thinning fluid.
"""
struct RitzBasis
    l::Int
    K::Int
    kind::Symbol
end
RitzBasis(l::Int, K::Int) = RitzBasis(l, K, :legendre)

npow(b::RitzBasis, k::Int) = b.l + 1 + 2 * (k - 1)

"""
    radial_window(l; margin = 0.9) -> Int

How many radial functions mode `l` can carry before its block stops being solvable.

The `x^(l+1)` factor confines every trial function to a layer of thickness `~1/l` at the
surface, so at high `l` they crowd together no matter what polynomial multiplies them and the
mode's mass matrix loses conditioning. That ceiling is a property of the mode, and it is what
makes a single global `K` wasteful: a `K` small enough for `l = 90` starves `l = 2`, and a `K`
large enough for `l = 2` makes `l = 90` singular.

MEASURED, on the `:legendre` basis, as the largest `K` keeping the normalised condition number
of the mode's mass matrix below 1e10 -- six digits of headroom in double precision:

    l       2   3   4   5   6   8  10  12  16  20  25  30  40  50  60  75  90
    K_max  40  29  21  16  14  11   9   8   7   6   5   5   4   4   4   4   3

A power law fits that to within one function everywhere:

    K(l) = 47.2 * l^(-0.642)

and the `margin` shaves it so the law never exceeds the measured ceiling -- at `margin = 0.9`
it overshoots at none of the seventeen truncations above. The exponent is the useful part: the
window closes like `l^(-0.64)`, much more slowly than the `1/l^2` one might guess, which is why
a staircase has real room to redistribute at low `l`.

The requirement runs the other way -- resolving the vortical field needs `K ~ |q| ~
l^{3/4}/sqrt(Oh)`, which GROWS with `l`. So the two curves cross, and past the crossing the
high modes are under-resolved by construction. That is tolerable only because those modes carry
almost none of the energy: the top mode holds about 3e-5 of the surface energy in a
representative run.

WHAT THIS DOES NOT PROMISE. It bounds the conditioning of the mode's MASS matrix and nothing
else. It is not a guarantee that a run will complete: on the 3000 ppm shear-thinning fluid,
`M = 60, K = 3` fails outright even though `radial_window(60) = 3` says that truncation is
within the window, while `M = 45, K = 3` runs. Whatever ends that run is not the conditioning
this law measures -- a variable viscosity makes the dissipation operator state-dependent and
the Picard closure nonlinear, neither of which enters here. Treat the window as a ceiling that
must not be exceeded, not as a licence to go up to it.
"""
radial_window(l::Integer; margin::Real = 0.9, kmin::Int = 1, kmax::Int = 40) =
    clamp(floor(Int, margin * 47.2 * float(l)^(-0.642)), kmin, kmax)

"""Legendre `P_j(u)` and its first two derivatives, by the standard recurrences."""
function legendre_uderivs(j::Int, u::Real)
    j == 0 && return (1.0, 0.0, 0.0)
    p0, p1 = 1.0, float(u)
    d0, d1 = 0.0, 1.0
    for n in 1:(j-1)
        p0, p1 = p1, ((2n + 1) * u * p1 - n * p0) / (n + 1)
        d0, d1 = d1, ((2n + 1) * (u * d1 + p0) - n * d0) / (n + 1)
    end
    ## P_j'' from the Legendre equation: (1-u^2)P'' - 2u P' + j(j+1) P = 0
    dd = abs(1 - u^2) < 1e-12 ? j*(j+1)*(j+1)*j/8 * one(p1) :
         (2u * d1 - j*(j+1) * p1) / (1 - u^2)
    (p1, d1, dd)
end

function phi(b::RitzBasis, k::Int, x::Real)
    b.kind === :monomial && return x^npow(b, k)
    m = b.l + 1
    P, _, _ = legendre_uderivs(k - 1, 2x^2 - 1)
    x^m * P
end

function dphi(b::RitzBasis, k::Int, x::Real)
    if b.kind === :monomial
        n = npow(b, k); return n * x^(n - 1)
    end
    m = b.l + 1
    P, dP, _ = legendre_uderivs(k - 1, 2x^2 - 1)
    m * x^(m - 1) * P + 4 * x^(m + 1) * dP
end

function d2phi(b::RitzBasis, k::Int, x::Real)
    if b.kind === :monomial
        n = npow(b, k); return n * (n - 1) * x^(n - 2)
    end
    m = b.l + 1
    P, dP, ddP = legendre_uderivs(k - 1, 2x^2 - 1)
    m * (m - 1) * x^(m - 2) * P + 4 * (2m + 1) * x^m * dP + 16 * x^(m + 2) * ddP
end

"""Legendre `P_l(mu)`, and the theta-derivatives the strain components need."""
function legendre_angular(l::Int, mu::Float64)
    p0, p1 = 1.0, mu
    if l == 0
        P, dPdmu = 1.0, 0.0
    elseif l == 1
        P, dPdmu = mu, 1.0
    else
        for n in 1:(l - 1)
            p0, p1 = p1, ((2n + 1) * mu * p1 - n * p0) / (n + 1)
        end
        P = p1
        dPdmu = l * (mu * p1 - p0) / (mu^2 - 1)
    end
    # Legendre's equation supplies the second derivative without a second recurrence
    d2Pdmu = (2mu * dPdmu - l * (l + 1) * P) / (1 - mu^2)
    s = sqrt(max(1 - mu^2, 0.0))
    dPdth = -s * dPdmu                        # d/dtheta
    d2Pdth = -mu * dPdmu + s^2 * d2Pdmu       # d2/dtheta2
    (P = P, dPdth = dPdth, d2Pdth = d2Pdth, sinth = s, mu = mu)
end

"""
Strain and velocity components of the field generated by `psi = f(x) C_l(theta)`,
from their analytic modal forms. Returns `(u_r, u_th, e_rr, e_tt, e_pp, e_rt)`.
"""
function modal_field(l::Int, f::Real, df::Real, d2f::Real, x::Real, A)
    L = l * (l + 1)
    u_r = f / x^2 * A.P
    u_th = df / (x * L) * A.dPdth
    e_rr = (df / x^2 - 2f / x^3) * A.P
    e_tt = df / (x^2 * L) * A.d2Pdth + f / x^3 * A.P
    cot_th = A.sinth > 1e-12 ? A.mu / A.sinth : 0.0
    e_pp = f / x^3 * A.P + df / (x^2 * L) * cot_th * A.dPdth
    T_op = d2f - 2df / x + L * f / x^2          # the tangential-stress operator
    e_rt = T_op / (2 * x * L) * A.dPdth
    (u_r, u_th, e_rr, e_tt, e_pp, e_rt)
end

"""
    assemble(b::RitzBasis, Oh; eta = (x, mu) -> 1.0, nx = 40, nmu = 40)

The three matrices of the variational statement, for one surface mode.

`eta` is the viscosity field in units of the zero-shear plateau, evaluated at
`(x, mu)`. For a Newtonian fluid it is `1`; for a generalized Newtonian one it is
`eta(gammadot)` on the current state, which makes `K` state-dependent.
"""
function assemble(b::RitzBasis, Oh::Real; eta = (x, mu) -> 1.0,
                  nx::Int = 40, nmu::Int = 40)
    l, K = b.l, b.K
    xs, wxs = gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = gauss_legendre_nodes(nmu, -1.0, 1.0)
    M = zeros(K, K)
    Kd = zeros(K, K)
    for (x, wx) in zip(xs, wxs), (mu, wmu) in zip(mus, wmus)
        A = legendre_angular(l, mu)
        w = wx * wmu * 2pi * x^2                     # dV = 2 pi x^2 dx dmu
        ev = eta(x, mu)
        F = ntuple(k -> modal_field(l, phi(b, k, x), dphi(b, k, x),
                                    d2phi(b, k, x), x, A), K)
        for a in 1:K, c in a:K
            Fa, Fc = F[a], F[c]
            M[a, c] += w * (Fa[1] * Fc[1] + Fa[2] * Fc[2])
            Kd[a, c] += w * 2 * ev * (Fa[3] * Fc[3] + Fa[4] * Fc[4] +
                                      Fa[5] * Fc[5] + 2 * Fa[6] * Fc[6])
        end
    end
    for a in 1:K, c in 1:a-1        # both forms are Hessians, hence symmetric
        M[a, c] = M[c, a]; Kd[a, c] = Kd[c, a]
    end
    Kd .*= Oh
    # surface energy: rank one, seeing the interior only through the trace at x = 1
    tr = [phi(b, k, 1.0) for k in 1:K]
    G = (4pi / (2l + 1)) * (l - 1) * (l + 2) * (tr * tr')
    (M = M, C = Kd, G = G)
end

"""
    decay_rates(b::RitzBasis, Oh; kwargs...) -> Vector{ComplexF64}

Returns `(sigma, cond_M)`: roots of the quadratic eigenvalue problem `(sigma^2 M - sigma C + G) a = 0`,
i.e. the free-decay rates of mode `l` for `chi ~ exp(-sigma t)`.

Linearised to a standard generalized eigenproblem of twice the size. The physical
root is the one with the smallest real part -- the least damped.
"""
function decay_rates(b::RitzBasis, Oh::Real; kwargs...)
    F = assemble(b, Oh; kwargs...)
    K = b.K
    # WHITEN by the Cholesky factor of M before linearising. The monomial trial
    # functions are Vandermonde-like, so M's condition number grows fast with K and
    # the raw pencil returns spurious near-zero roots: at K = 8 omega^2 collapses to
    # 1e-9 and lambda halves, silently, while K = 5 is accurate to 1e-11. Whitening
    # moves the usable range out; it does not make the basis good, and `cond_M` is
    # returned so a caller can see when it has run out of room.
    Lc = cholesky(Symmetric(F.M)).L
    Ci = Lc \ (F.C / Lc')
    Gi = Lc \ (F.G / Lc')
    # [sigma^2 I - sigma Ci + Gi] y = 0  linearised:
    Abig = [Ci  -Gi; Matrix(I, K, K)  zeros(K, K)]
    Bbig = Matrix(I, 2K, 2K)
    (sigma = eigvals(Abig, Bbig), cond_M = cond(F.M))
end

"""
    dominant_pair(b::RitzBasis, Oh; kwargs...) -> (lambda, omega2)

The least-damped root pair, in the convention the Newtonian literature uses:
`zeta.. + 2 lambda zeta. + omega^2 zeta = 0`, returning `(lambda, omega^2)`.
Comparable directly against `reid_lambda_omega2`.

Note which `omega^2` this is. For that oscillator the roots are
`sigma = lambda +- i sqrt(omega^2 - lambda^2)`, so by Vieta `omega^2` is their
PRODUCT, `|sigma|^2` -- not `Im(sigma)^2`, which is the damped frequency
`omega^2 - lambda^2`. Reporting the latter leaves a discrepancy against Reid of
exactly `lambda^2/omega^2` that does not shrink as the basis is refined: a plateau
at `6.1e-3` for `Oh = 0.05` and `1.6e-1` for `Oh = 0.3`, while `lambda` itself was
converging to machine precision. A convergence study is what exposes a convention
error of this kind, because a wrong convention converges to the wrong number.
"""
function dominant_pair(b::RitzBasis, Oh::Real; kwargs...)
    R = decay_rates(b, Oh; kwargs...)
    σ = R.sigma
    finite = filter(s -> isfinite(s) && real(s) > 1e-12, σ)
    length(finite) < 2 && return (NaN, NaN)
    # VIETA ON THE TWO SLOWEST ROOTS, which is what lambda and omega^2 are. For an
    # underdamped mode those two are a conjugate pair and this reduces to
    # (Re sigma, |sigma|^2); above the critical Ohnesorge the pair goes REAL and the
    # conjugate reading breaks -- at Oh = 1, l = 2 it is wrong by 64%, and no amount
    # of basis refinement helps, because the two roots are then different numbers.
    # Taking the two slowest and applying Vieta covers both regimes.
    idx = partialsortperm(real.(finite), 1:2)
    s1, s2 = finite[idx[1]], finite[idx[2]]
    (real(s1 + s2) / 2, real(s1 * s2))
end

# ---------------------------------------------------------------------------
# Multi-mode assembly, which is where shear thinning actually bites.
#
# For a constant viscosity the modes do not talk to each other: eta has only its
# l = 0 harmonic, the Gaunt coefficient is a delta, and every block off the
# diagonal vanishes. That is why the single-mode assembly above is enough for the
# Newtonian case. A shear-thinning fluid destroys it: eta = eta(gammadot) inherits
# the angular structure of the strain field, its harmonics k >= 1 are non-zero, and
# the dissipation form acquires off-diagonal blocks. Those blocks ARE the physics
# of the shear-thinning problem, so the assembly has to carry them.
#
# gammadot is computed from the WHOLE field, never mode by mode: the strain tensor
# superposes over modes but its invariant does not, so there is no such thing as one
# mode's shear rate once more than one is active.

"""
    ModalBasis(ls, K)

Radial trial functions `K` per surface mode, over the modes `ls` (e.g. `2:M`).
Degrees of freedom are ordered mode-major: `(l, k) -> (i-1)*K + k` for the `i`-th
mode in `ls`.
"""
struct ModalBasis
    ls::Vector{Int}
    K::Int
    kind::Symbol
end
ModalBasis(ls, K::Int, kind::Symbol = :legendre) = ModalBasis(collect(ls), K, kind)
ndof(b::ModalBasis) = length(b.ls) * b.K
dofindex(b::ModalBasis, i::Int, k::Int) = (i - 1) * b.K + k

"""
    strain_at(b::ModalBasis, x, mu) -> Matrix

The six field components of every degree of freedom at one quadrature point, as a
`ndof x 6` array in the order `(u_r, u_th, e_rr, e_tt, e_pp, e_rt)`. Assembling any
of the quadratic forms is then a contraction over rows.
"""
function strain_at(b::ModalBasis, x::Float64, mu::Float64)
    F = zeros(ndof(b), 6)
    for (i, l) in enumerate(b.ls)
        A = legendre_angular(l, mu)
        rb = RitzBasis(l, b.K, b.kind)
        for k in 1:b.K
            f, df, d2f = phi(rb, k, x), dphi(rb, k, x), d2phi(rb, k, x)
            F[dofindex(b, i, k), :] .= modal_field(l, f, df, d2f, x, A)
        end
    end
    F
end

"""Full double contraction `e:e` of a strain given as `(e_rr, e_tt, e_pp, e_rt)`."""
ddot_strain(a, c) = a[3]*c[3] + a[4]*c[4] + a[5]*c[5] + 2*a[6]*c[6]

"""
    shear_rate(b::ModalBasis, a, x, mu) -> Float64

`gammadot = sqrt(2 e:e)` of the superposed field with velocity amplitudes `a`.
Built from the summed tensor, not from a per-mode sum of invariants.
"""
function shear_rate(b::ModalBasis, a::AbstractVector, x::Float64, mu::Float64)
    F = strain_at(b, x, mu)
    e = zeros(6)
    for d in 1:ndof(b), c in 1:6
        e[c] += a[d] * F[d, c]
    end
    sqrt(max(2 * ddot_strain(e, e), 0.0))
end

"""
Cached geometry for the coupled assembly, keyed by basis and quadrature.

WHY THIS EXISTS. With a variable viscosity the operator is rebuilt on every Picard sweep, and
profiling put essentially the entire cost of a shear-thinning run there: 101 ms per sweep
against 0.1 ms for the whole contact solve, so a single impact spent 145 of its 145 seconds
assembling. But almost nothing in the assembly depends on the state. `strain_at` is geometry,
the mass matrix has no `eta` in it at all, and the stiffness is a closed form. The only thing
that changes between sweeps is the scalar `eta` at each quadrature point.

So the pairwise contractions `e^p : e^q` are computed once per basis and stored, and a sweep
becomes a weighted sum over quadrature points -- `npairs * nq` multiply-adds, no tensor
algebra and no allocation.
"""
struct CoupledGeometry
    M::Matrix{Float64}          # mass, state-independent
    G::Matrix{Float64}          # stiffness, state-independent
    D::Matrix{Float64}          # (npairs x nq) contractions e^p:e^q at each quadrature point
    F::Array{Float64,3}         # (N x 6 x nq) the strain basis itself
    w::Vector{Float64}          # (nq) quadrature weights, including the 2 pi x^2 Jacobian
    xs::Vector{Float64}         # (nq) radial coordinate of each point
    mus::Vector{Float64}        # (nq) angular coordinate
end

## The key MUST include the basis kind. It did not, and a monomial-populated cache was served
## to a :legendre run -- which made two different bases return identical numbers and looked
## exactly like a confirmation that they span the same space. Same failure as keying a results
## store on parameters but not on what produced them.
const COUPLED_CACHE = Dict{Tuple{Vector{Int},Int,Symbol,Int,Int},CoupledGeometry}()

"""
Rough footprint of the cached geometry, in bytes.

Two tables: the pair contractions, `npairs x nq`, which dominate, and the strain basis itself,
`N x 6 x nq`, which is smaller but is what lets the viscosity be evaluated without rebuilding
anything.
"""
coupled_cache_bytes(N::Int, nx::Int, nmu::Int) =
    ((N*(N+1) ÷ 2) + 6N) * nx * nmu * 8

"""Above this the geometry is not cached and the direct loop is used instead."""
const COUPLED_CACHE_BUDGET = 500_000_000

"""Footprints above this are reported once, because a silent 300 MB allocation is a surprise."""
const COUPLED_CACHE_WARN = 100_000_000

const COUPLED_CACHE_ANNOUNCED = Set{Tuple{Vector{Int},Int,Symbol,Int,Int}}()

function coupled_geometry(b::ModalBasis, nx::Int, nmu::Int)
    key = (b.ls, b.K, b.kind, nx, nmu)
    haskey(COUPLED_CACHE, key) && return COUPLED_CACHE[key]
    N = ndof(b)
    xs, wxs = gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = gauss_legendre_nodes(nmu, -1.0, 1.0)
    nq = nx * nmu; npairs = N*(N+1) ÷ 2
    bytes = coupled_cache_bytes(N, nx, nmu)
    if bytes >= COUPLED_CACHE_WARN && !(key in COUPLED_CACHE_ANNOUNCED)
        push!(COUPLED_CACHE_ANNOUNCED, key)
        @warn "caching the coupled geometry: this is a large allocation" *
              " (raise the truncation further and it will be dropped for the slow path)" ndof=N nx nmu megabytes=round(bytes/1e6, digits=1) budget_MB=COUPLED_CACHE_BUDGET÷1_000_000
    end
    D = Matrix{Float64}(undef, npairs, nq)
    Fs = Array{Float64,3}(undef, N, 6, nq)
    w = Vector{Float64}(undef, nq)
    xv = Vector{Float64}(undef, nq); muv = Vector{Float64}(undef, nq)
    Mm = zeros(N, N)
    q = 0
    for (ix, x) in enumerate(xs), (im, mu) in enumerate(mus)
        q += 1
        F = strain_at(b, x, mu)
        ww = wxs[ix] * wmus[im] * 2pi * x^2
        w[q] = ww; xv[q] = x; muv[q] = mu
        @inbounds for d in 1:N, c in 1:6
            Fs[d, c, q] = F[d, c]
        end
        k = 0
        for p in 1:N, r in p:N
            k += 1
            D[k, q] = ddot_strain(view(F, p, :), view(F, r, :))
            Mm[p, r] += ww * (F[p,1]*F[r,1] + F[p,2]*F[r,2])
        end
    end
    for p in 1:N, r in 1:p-1
        Mm[p, r] = Mm[r, p]
    end
    g = CoupledGeometry(Mm, stiffness_matrix(b), D, Fs, w, xv, muv)
    COUPLED_CACHE[key] = g
    g
end

"""
    shear_rate_at(g::CoupledGeometry, a, q) -> Float64

Shear-rate invariant at cached quadrature point `q`, from the stored strain basis.

This is the half of the assembly the first cache missed. `shear_rate` rebuilds `strain_at` --
the whole `N x 6` tensor basis -- on every call, and the viscosity closure calls it once per
quadrature point per Picard sweep, so the geometry that was cached out of the pair loop was
still being rebuilt on the way in. Here it is a mat-vec against a table that already exists.
"""
function shear_rate_at(g::CoupledGeometry, a::AbstractVector, q::Int)
    N = length(a)
    e1 = e2 = e3 = e4 = e5 = e6 = 0.0
    @inbounds for d in 1:N
        ad = a[d]
        ad == 0 && continue
        e1 += ad*g.F[d,1,q]; e2 += ad*g.F[d,2,q]; e3 += ad*g.F[d,3,q]
        e4 += ad*g.F[d,4,q]; e5 += ad*g.F[d,5,q]; e6 += ad*g.F[d,6,q]
    end
    sqrt(max(2 * ddot_strain((e1,e2,e3,e4,e5,e6), (e1,e2,e3,e4,e5,e6)), 0.0))
end

"""The surface-energy Hessian, which has a closed form and no quadrature."""
function stiffness_matrix(b::ModalBasis)
    Gm = zeros(ndof(b), ndof(b))
    for (i, l) in enumerate(b.ls)
        rb = RitzBasis(l, b.K, b.kind)
        tr = [phi(rb, k, 1.0) for k in 1:b.K]
        blk = (4pi / (2l + 1)) * (l - 1) * (l + 2) * (tr * tr')
        rng = dofindex(b, i, 1):dofindex(b, i, b.K)
        Gm[rng, rng] .= blk
    end
    Gm
end

"""
    assemble_coupled(b::ModalBasis, Oh; eta = (x, mu) -> 1.0, nx = 40, nmu = 48)

The three forms over all retained modes. `eta` is the viscosity in units of the
zero-shear plateau, as a function of position -- for a shear-thinning fluid it is
`eta(shear_rate(...))` on the current state, which is what makes `C` couple modes.

The pairwise strain contractions are cached per basis (see [`CoupledGeometry`](@ref)), so a
repeated call with a different `eta` costs a weighted sum rather than a fresh quadrature.
"""
function assemble_coupled(b::ModalBasis, Oh::Real; eta = (x, mu) -> 1.0,
                          eta_rate = nothing, state = nothing,
                          nx::Int = 40, nmu::Int = 48)
    N = ndof(b)
    ## Cache the geometry when it fits the budget; fall back to the direct loop when it does
    ## not, so a very large truncation degrades in speed rather than in memory.
    if coupled_cache_bytes(N, nx, nmu) <= COUPLED_CACHE_BUDGET
        g = coupled_geometry(b, nx, nmu)
        npairs = N*(N+1) ÷ 2
        acc = zeros(npairs)
        ## `eta_rate` + `state` is the fast path: the shear rate comes from the cached strain
        ## basis instead of rebuilding it. Passing `eta` as a function of position still works
        ## and is what a purely spatial viscosity wants, but for a shear-thinning fluid it
        ## makes the cache pointless, because the closure rebuilds what the cache stores.
        fast = eta_rate !== nothing && state !== nothing
        @inbounds for q in eachindex(g.w)
            ev = fast ? eta_rate(shear_rate_at(g, state, q)) : eta(g.xs[q], g.mus[q])
            c = g.w[q] * 2 * ev
            c == 0 && continue
            @simd for k in 1:npairs
                acc[k] += c * g.D[k, q]
            end
        end
        Cm = Matrix{Float64}(undef, N, N)
        k = 0
        @inbounds for p in 1:N, r in p:N
            k += 1
            Cm[p, r] = acc[k]; Cm[r, p] = acc[k]
        end
        Cm .*= Oh
        return (M = g.M, C = Cm, G = g.G)
    end
    xs, wxs = gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = gauss_legendre_nodes(nmu, -1.0, 1.0)
    Mm, Cm = zeros(N, N), zeros(N, N)
    for (x, wx) in zip(xs, wxs), (mu, wmu) in zip(mus, wmus)
        F = strain_at(b, x, mu)
        w = wx * wmu * 2pi * x^2
        ev = eta(x, mu)
        for p in 1:N, q in p:N
            Mm[p, q] += w * (F[p,1]*F[q,1] + F[p,2]*F[q,2])
            Cm[p, q] += w * 2 * ev * ddot_strain(view(F, p, :), view(F, q, :))
        end
    end
    for p in 1:N, q in 1:p-1
        Mm[p, q] = Mm[q, p]; Cm[p, q] = Cm[q, p]
    end
    Cm .*= Oh
    (M = Mm, C = Cm, G = stiffness_matrix(b))
end

"""
    block_norm(b::ModalBasis, A, i, j) -> Float64

Largest entry of the `(i, j)` mode block of `A`. The measure of whether modes
`b.ls[i]` and `b.ls[j]` are coupled.
"""
function block_norm(b::ModalBasis, A::AbstractMatrix, i::Int, j::Int)
    ri = dofindex(b, i, 1):dofindex(b, i, b.K)
    rj = dofindex(b, j, 1):dofindex(b, j, b.K)
    maximum(abs, view(A, ri, rj))
end

"""
    carreau(gammadot; eta_inf_ratio, lambda_c, a, n) -> Float64

Carreau-Yasuda viscosity in units of the zero-shear plateau. `n < 1` is
shear-thinning; `eta_inf_ratio` is the infinite-shear plateau, which (H2) requires
to be strictly positive.
"""
function carreau(gd::Real; eta_inf_ratio::Real = 0.01, lambda_c::Real = 1.0,
                 a::Real = 2.0, n::Real = 0.5)
    eta_inf_ratio + (1 - eta_inf_ratio) * (1 + (lambda_c * gd)^a)^((n - 1) / a)
end

"""
    assemble_newtonian(b::ModalBasis, Oh) -> (M, C, G)

The constant-viscosity operator, assembled block by block.

At constant `eta` the three forms are block diagonal in `l` -- `eta` carries only its
`l = 0` harmonic, so it cannot move energy between surface modes -- and they are also
independent of the state and therefore of time. Both facts are already established:
the block diagonality is `test_variational.jl`'s "constant viscosity leaves the modes
uncoupled", and this routine is checked against `assemble_coupled` there too.

Building the blocks separately is what makes a realistic mode count affordable. The
coupled routine costs `ndof^2` quadratures because it cannot assume the structure;
at `l` up to 40 that is 78 x 78 pairs on a 40 x 48 grid, per assembly, and the time
stepper wants an assembly per contact candidate per step. Block by block it is
`length(ls)` small dense problems, built once before the march begins.
"""
function assemble_newtonian(b::ModalBasis, Oh::Real; nx::Int = 40)
    N = ndof(b)
    M = zeros(N, N); C = zeros(N, N); G = zeros(N, N)
    for (i, l) in enumerate(b.ls)
        F = assemble(RitzBasis(l, b.K, b.kind), Oh; nx = nx)
        for ka in 1:b.K, kb in 1:b.K
            ia = dofindex(b, i, ka); ib = dofindex(b, i, kb)
            M[ia, ib] = F.M[ka, kb]
            C[ia, ib] = F.C[ka, kb]
            G[ia, ib] = F.G[ka, kb]
        end
    end
    (M = M, C = C, G = G)
end
