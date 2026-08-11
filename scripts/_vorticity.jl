# The vorticity of the interior flow, and the blue palettes used to draw it.
#
# Shared by `figure_vorticity.jl` (space-time map) and `animate_vorticity.jl` (the
# drop with its interior coloured), so the two cannot drift apart.
#
# THE FORMULA
#
# The interior velocity is a Stokes stream function. For one surface mode `l` with
# radial trial function `f(x)`,
#
#     Psi = f(x) C_l(theta),        C_l = -sin(theta) P_l'(theta) / L,   L = l(l+1)
#
# which is what makes the solver's field read u_r = f/x^2 P_l and u_th = f'/(x L) P_l'.
# Axisymmetric flow has one vorticity component, omega_phi = -(1/(x sin th)) E^2 Psi,
# with E^2 = d_xx + (sin th / x^2) d_th( (1/sin th) d_th ).
#
# The angular operator acting on C_l returns -L C_l -- that is the Gegenbauer equation,
# and it is why this basis is the natural one -- so E^2 Psi = (f'' - L f/x^2) C_l and
#
#     omega_phi = (1/(x L)) (f'' - L f / x^2) dP_l/dtheta .
#
# The sin(theta) cancels, so this is finite on the axis. Note the bracket: a trial
# function with f'' = L f/x^2 carries no vorticity at all. The first radial function is
# f = x^(l+1), for which f'' = l(l+1) x^(l-1) = L f/x^2 exactly, so K = 1 is
# irrotational and every bit of rotation drawn comes from k >= 2. `assert_irrotational`
# checks that before anything is plotted, because if it failed the picture would be
# showing a basis artefact rather than the flow.
#
# The amplitudes are `adot`, not `a`: the mass matrix is built from the field components
# as a kinetic energy, so it is the RATE that multiplies the field to give the velocity.
# Using `a` would draw the vorticity of the displacement, which is not a physical
# quantity.

using Printf
using DropSolver

"""
    vorticity(p, adot, x, mu) -> Float64

Azimuthal vorticity of the interior flow at radius `x` and `mu = cos(theta)`.
"""
function vorticity(p::ImpactParams, adot::AbstractVector, x::Real, mu::Real)
    b = DropSolver.basis(p)
    w = 0.0
    for (i, l) in enumerate(b.ls)
        L  = l * (l + 1)
        A  = DropSolver.legendre_angular(l, mu)
        rb = DropSolver.RitzBasis(l, p.K, p.basis_kind)
        s  = 0.0
        for k in 1:p.K
            f   = DropSolver.phi(rb, k, x)
            d2f = DropSolver.d2phi(rb, k, x)
            s  += adot[DropSolver.dofindex(b, i, k)] * (d2f - L * f / x^2)
        end
        w += s * A.dPdth / (x * L)
    end
    w
end

"""Enstrophy, the volume integral of `omega^2`, by Gauss-Legendre in both directions."""
function enstrophy(p::ImpactParams, adot::AbstractVector; nx = 24, nmu = 32)
    xs, wxs   = DropSolver.gauss_legendre_nodes(nx, 0.0, 1.0)
    mus, wmus = DropSolver.gauss_legendre_nodes(nmu, -1.0, 1.0)
    s = 0.0
    for (x, wx) in zip(xs, wxs), (mu, wmu) in zip(mus, wmus)
        s += wx * wmu * 2pi * x^2 * vorticity(p, adot, x, mu)^2   # dV = 2 pi x^2 dx dmu
    end
    s
end

"""
    radial_kernel(p) -> (xs, Phi)

`Phi[i, k, n] = (phi_k''(x_n) - L phi_k(x_n)/x_n^2) / (x_n L)` for surface mode `i`,
on the radial grid `xs`. Everything in the vorticity except the amplitudes and the
angle, so a frame costs one small contraction instead of a fresh basis evaluation per
point. Used by the animation, which needs the field on a grid at every frame.
"""
function radial_kernel(p::ImpactParams; nx::Int = 90)
    b   = DropSolver.basis(p)
    xs  = collect(range(1e-3, 1.0; length = nx))     # x = 0 is a removable singularity
    Phi = zeros(length(b.ls), p.K, nx)
    for (i, l) in enumerate(b.ls)
        L  = l * (l + 1)
        rb = DropSolver.RitzBasis(l, p.K, p.basis_kind)
        for k in 1:p.K, (n, x) in enumerate(xs)
            Phi[i, k, n] = (DropSolver.d2phi(rb, k, x) -
                            L * DropSolver.phi(rb, k, x) / x^2) / (x * L)
        end
    end
    (xs, Phi)
end

"""
    angular_kernel(p, mus) -> (P, dPdth)

`P[i, j]` and `dPdth[i, j]` for surface mode `i` at `mus[j]`. Independent of the state,
so it is built once for a whole animation. `P` gives the drop's outline, `dPdth` the
vorticity.
"""
function angular_kernel(p::ImpactParams, mus::AbstractVector)
    ls = DropSolver.basis(p).ls
    P  = zeros(length(ls), length(mus)); D = similar(P)
    for (i, l) in enumerate(ls), (j, mu) in enumerate(mus)
        A = DropSolver.legendre_angular(l, mu)
        P[i, j] = A.P; D[i, j] = A.dPdth
    end
    (P, D)
end

"""
    vorticity_grid(p, adot, Phi, D) -> Matrix

`|omega|` on the (radius, angle) grid the two kernels were built on. The sum over modes
and radial functions is a contraction, so a frame is a couple of small matrix products
rather than a basis evaluation per grid point.
"""
function vorticity_grid(p::ImpactParams, adot::AbstractVector, Phi, D)
    b  = DropSolver.basis(p)
    nx = size(Phi, 3)
    G  = zeros(length(b.ls), nx)                     # g_l(x_n)
    for (i, _) in enumerate(b.ls), k in 1:p.K
        a = adot[DropSolver.dofindex(b, i, k)]
        a == 0 && continue
        @views G[i, :] .+= a .* Phi[i, k, :]
    end
    abs.(transpose(G) * D)                           # (nx x nmu)
end

"""The irrotational check: `K = 1` must carry no vorticity at all."""
function assert_irrotational(; We = 0.5, Bo = 0.019, Oh = 0.0373, M = 8, tol = 1e-10)
    p1 = ImpactParams(We = We, Bo = Bo, Oh = Oh, M = M, K = 1, t_max = 1.0)
    v  = [vorticity(p1, ones(DropSolver.ndof(DropSolver.basis(p1))), x, mu)
          for x in (0.3, 0.7, 0.95), mu in (-0.9, -0.2, 0.5)]
    worst = maximum(abs, v)
    @assert worst < tol "K = 1 is not irrotational: max |omega| = $worst"
    println("self-check: K = 1 carries no vorticity (max |omega| = ",
            @sprintf("%.1e", worst), ")")
    worst
end

# ---------------------------------------------------------------------------
# Palettes in the manner of MATLAB's `abyss` and `sky`.
#
# Neither ships with Julia and neither can be copied out of MATLAB here, so these are
# hand-built gradients with the same character: `abyss` runs from near-black navy up
# through blue to a pale sky, `sky` is the light half of that range alone. Both are
# sequential, which is what |omega| wants.
#
# Requires `Plots` to be loaded by the including script, for `cgrad` and `RGB`.
# ---------------------------------------------------------------------------
const ABYSS = cgrad([RGB(0.000, 0.012, 0.086), RGB(0.008, 0.106, 0.278),
                     RGB(0.016, 0.239, 0.478), RGB(0.055, 0.400, 0.667),
                     RGB(0.302, 0.612, 0.816), RGB(0.749, 0.882, 0.945)])
const SKY   = cgrad([RGB(0.055, 0.220, 0.400), RGB(0.145, 0.420, 0.667),
                     RGB(0.353, 0.639, 0.831), RGB(0.643, 0.831, 0.925),
                     RGB(0.878, 0.949, 0.980)])
const DEEP  = RGB(0.016, 0.239, 0.478)      # a line colour from the same family

# ---------------------------------------------------------------------------
# THE LOCAL VISCOSITY FIELD
#
# A generalized Newtonian fluid has no single viscosity: eta depends on the local shear
# rate, so during an impact the drop is a map of stiff and soft regions that moves. That
# map is what an "effective viscosity" is being asked to stand in for, and drawing it is
# the most direct statement of how much is being averaged over.
#
# `shear_rate` evaluates gammadot = sqrt(2 e:e) one point at a time, rebuilding the
# strain basis at every call. On an animation grid that is around 10^10 operations and
# is not worth waiting for. But the strain basis does not depend on the state: only the
# amplitudes do. So the basis is built once as a matrix and every frame becomes a single
# GEMV, which is a few milliseconds.
# ---------------------------------------------------------------------------

"""
    strain_kernel(p, xs, mus) -> Matrix

The six strain components of every degree of freedom at every grid point, flattened to a
`(npoints*6) x ndof` matrix. State-independent, so it is built once per animation.
"""
function strain_kernel(p::ImpactParams, xs, mus)
    b = DropSolver.basis(p); N = DropSolver.ndof(b)
    K = zeros(length(xs) * length(mus) * 6, N)
    row = 0
    for x in xs, mu in mus
        F = DropSolver.strain_at(b, x, mu)
        for c in 1:6
            row += 1
            @views K[row, :] .= F[:, c]
        end
    end
    K
end

"""
    shear_rate_grid(Kern, adot, nx, nmu) -> Matrix

`gammadot` on the grid `strain_kernel` was built for, in solver units. One matrix-vector
product, then the second invariant of the summed tensor -- not a per-mode sum of
invariants, which would be a different and wrong quantity.
"""
function shear_rate_grid(Kern, adot::AbstractVector, nx::Int, nmu::Int)
    e = Kern * adot                      # all six components at all points
    out = zeros(nx, nmu)
    idx = 0
    for i in 1:nx, j in 1:nmu
        a = e[idx+1]; bb = e[idx+2]; c = e[idx+3]
        d = e[idx+4]; f = e[idx+5]; g = e[idx+6]
        idx += 6
        ## ddot_strain uses components 3..6 as (e_rr, e_tt, e_pp, e_rt); 1..2 are the
        ## velocities and take no part in the rate of strain.
        out[i, j] = sqrt(max(2 * (c*c + d*d + f*f + 2*g*g), 0.0))
    end
    out
end
