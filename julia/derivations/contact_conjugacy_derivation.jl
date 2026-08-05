# # What the film pressure is conjugate to
#
# This settles, symbolically, why the contact compliance in `simulate_lcp` comes out
# asymmetric, and what would have to change for it to be symmetric. It is not a style
# question: symmetry decides whether the contact problem is a convex programme with a unique
# solution or a general complementarity problem, and therefore which algorithm is legitimate.
#
# THE SHORT VERSION. A constraint and its Lagrange multiplier are conjugate by construction:
# the generalised force of a constraint `h(xi) >= 0` is the transpose of its own Jacobian.
# The solver imposes a constraint on the VERTICAL gap at COLLOCATION NODES, but forces the
# drop with a RADIAL pressure projected in the LEGENDRE basis. Those are two independent
# departures from conjugacy -- direction, and projection-versus-evaluation -- and either alone
# destroys symmetry.
#
# Everything below is derived rather than asserted, and the pieces that are claims about the
# shipped code are checked against `DropSolver` directly at the end.

using Symbolics
using LinearAlgebra
using QuadGK
using DropSolver

@variables mu z lam                                     #src

# ## 1. The gap, and its Jacobian
#
# The surface sits at radius `r = 1 + zeta(theta)` with `zeta = sum_l zeta_l P_l(mu)`, so a
# surface point is at height
#
# ```math
# h(\theta) = z + \big(1 + \zeta(\theta)\big)\cos\theta ,
# ```
#
# measured from the substrate. This is a VERTICAL distance: it is what "the drop must not
# enter the wall" means, and it is what the solver constrains.

Pl(l, x) = l == 0 ? one(x) : l == 1 ? x :
           ((2l - 1) * x * Pl(l - 1, x) - (l - 1) * Pl(l - 2, x)) / l   #src

const LS = 2:5                                                          #src
zetas = [Symbolics.variable(:z, l) for l in LS]                         #src
zeta_field = sum(zetas[i] * Pl(l, mu) for (i, l) in enumerate(LS))      #src
h_sym = z + (1 + zeta_field) * mu                                       #src

# The Jacobian of the constraint with respect to the shape amplitudes is then

dh = [Symbolics.derivative(h_sym, zetas[i]) for i in eachindex(LS)]      #src
for (i, l) in enumerate(LS)                                             #src
    @assert iszero(Symbolics.simplify(dh[i] - mu * Pl(l, mu)))           #src
end                                                                      #src

# ```math
# \frac{\partial h(\theta)}{\partial \zeta_l} = \cos\theta\, P_l(\cos\theta) ,
# ```
#
# and that factor of ``\cos\theta`` is the whole story of this page. It is present because the
# gap is vertical while ``\zeta_l`` displaces the surface radially. It is exactly what
# `gap_row` computes, which is checked against the shipped code in §5.

# ## 2. The forcing the model actually uses
#
# The model treats the film pressure as a smooth field on the surface, ``p_c(\theta) = \sum_k
# p_{c,k}P_k(\cos\theta)``, acting along the surface normal — radially, to first order in the
# deformation. Its virtual work against a radial surface displacement is
#
# ```math
# \delta W = -\!\int p_c\,\delta\zeta \,dA = -2\pi\!\int_{-1}^{1} p_c(\mu)\,\delta\zeta(\mu)\,d\mu ,
# ```
#
# and orthogonality ``\int_{-1}^1 P_kP_l\,d\mu = 2\delta_{kl}/(2l+1)`` collapses it to one term
# per mode.

pcs = [Symbolics.variable(:p, k) for k in LS]                            #src
pc_field = sum(pcs[i] * Pl(l, mu) for (i, l) in enumerate(LS))           #src
## virtual work, integrated exactly in the Legendre basis                #src
W_radial = -2pi * sum(                                                   #src
    pcs[i] * zetas[j] * quadgk(m -> Pl(LS[i], m) * Pl(LS[j], m), -1, 1)[1] #src
    for i in eachindex(LS), j in eachindex(LS))                          #src
Q_radial = [Symbolics.derivative(W_radial, zetas[j]) for j in eachindex(LS)] #src
for (j, l) in enumerate(LS)                                              #src
    @assert abs(Symbolics.value(Symbolics.substitute(                    #src
        Q_radial[j] + (4pi / (2l + 1)) * pcs[j],                          #src
        Dict(p => 1.0 for p in pcs)))) < 1e-10                            #src
end                                                                       #src

# ```math
# \boxed{\;Q_{\zeta_l} = -\frac{4\pi}{2l+1}\,p_{c,l}\;}
# ```
#
# In words: each shape mode is pushed only by the pressure harmonic of the same order, with a
# coefficient set by the norm of that harmonic. This is the forcing in `force_column`, and it
# is *correct* — it is the exact virtual work of a smooth radial pressure field. Nothing here
# is a bug. The problem is what happens when it is paired with a nodal gap constraint.

# ## 3. What conjugacy requires
#
# Now impose the constraint properly. Contact is `h_i >= 0` at the `n` collocation nodes, with
# multipliers ``\lambda_i``. Writing the gap as ``\bm h = \bm H\bm\xi + \bm c`` and the step
# operator as ``\bm A = \beta^2\bm M + \beta\bm C + \bm G`` (symmetric, because each of the
# three is a Hessian of a quadratic form), the stationarity condition of
#
# ```math
# \tfrac12\bm\xi^{\mathsf T}\bm A\bm\xi - \bm f^{\mathsf T}\bm\xi - \sum_i \lambda_i h_i(\bm\xi)
# ```
#
# is ``\bm A\bm\xi = \bm f + \bm H^{\mathsf T}\bm\lambda``. The multiplier's generalised force
# is the TRANSPOSE OF THE CONSTRAINT JACOBIAN — not a modelling choice, a consequence of
# differentiating the constraint. Eliminating ``\bm\xi``,
#
# ```math
# \boxed{\;\bm h = \bm W\bm\lambda + \bm b,\qquad \bm W = \bm H\bm A^{-1}\bm H^{\mathsf T}\;}
# ```
#
# `W` is symmetric for any `H` whenever `A` is symmetric, and positive semi-definite whenever
# `A` is positive definite:

let                                                                       #src
    ## symbolic proof of both properties on a small dense instance        #src
    @variables a11 a12 a22 h11 h12 h21 h22                               #src
    A = [a11 a12; a12 a22]; H = [h11 h12; h21 h22]                        #src
    Wm = H * inv(A) * transpose(H)                                        #src
    ## symmetry, identically in every entry                              #src
    @assert iszero(Symbolics.simplify(Wm[1, 2] - Wm[2, 1]; expand = true)) #src
    ## and the quadratic form is A^-1 evaluated on y = H'lam, hence non-negative for A PD    #src
    @variables l1 l2                                                      #src
    lamv = [l1; l2]; y = transpose(H) * lamv                              #src
    q1 = Symbolics.simplify(transpose(lamv) * Wm * lamv; expand = true)    #src
    q2 = Symbolics.simplify(transpose(y) * inv(A) * y; expand = true)      #src
    @assert iszero(Symbolics.simplify(q1 - q2; expand = true))             #src
end                                                                        #src

# ```math
# \bm\lambda^{\mathsf T}\bm W\bm\lambda
#   = \big(\bm H^{\mathsf T}\bm\lambda\big)^{\mathsf T}\bm A^{-1}\big(\bm H^{\mathsf T}\bm\lambda\big) \ge 0 .
# ```
#
# So with the conjugate forcing the contact problem is the KKT system of
#
# ```math
# \min_{\bm\lambda\ge0}\ \tfrac12\bm\lambda^{\mathsf T}\bm W\bm\lambda + \bm b^{\mathsf T}\bm\lambda ,
# ```
#
# a convex quadratic programme: a solution exists, it is unique when `W` is definite, and a
# projected Gauss-Seidel sweep converges to it. None of that is available otherwise.

# ## 4. Where the shipped pairing departs, and by exactly how much
#
# The shipped force is §2's, re-expressed in nodal pressures through the Vandermonde inverse:
# ``\bm Q_n = \bm Q_{\rm modal}\bm V^{-1}``. Ask when that could equal ``\bm H^{\mathsf T}\bm D``
# for some positive diagonal `D` — the most general rescaling that preserves symmetry after the
# change of variables ``\tilde{\bm\lambda} = \bm D^{1/2}\bm p``.
#
# Do the same virtual work by NODAL QUADRATURE instead of exact projection. With Gauss weights
# `w_i` and ``\delta\zeta(\mu_i) = \sum_l P_l(\mu_i)\delta\zeta_l``,
#
# ```math
# \delta W = -2\pi\sum_i w_i\,p_i\,\delta\zeta(\mu_i)
# \quad\Longrightarrow\quad
# Q_{\zeta_l} = -2\pi\sum_i w_i\,p_i\,P_l(\mu_i) ,
# ```
#
# whose `i`-th column is ``-2\pi w_i P_l(\mu_i)``. Compare the constraint Jacobian's `i`-th
# row, ``\mu_i P_l(\mu_i)``: the two differ by the single factor ``-2\pi w_i/\mu_i``, which is
# diagonal. So

let                                                                        #src
    nq = 6                                                                 #src
    mus, ws = DropSolver.gauss_legendre_nodes(nq, -1.0, 1.0)                #src
    lset = 2:5                                                             #src
    H = [m * Pl(l, m) for m in mus, l in lset]                              #src
    Qquad = [-2pi * w * Pl(l, m) for l in lset, (m, w) in zip(mus, ws)]      #src
    D = Diagonal([-2pi * ws[i] / mus[i] for i in 1:nq])                      #src
    @assert maximum(abs, Qquad - transpose(H) * D) < 1e-12                   #src
    ## and the resulting compliance is symmetric after the diagonal is split evenly           #src
    Asym = Symmetric(Matrix(I, length(lset), length(lset)) .+ 0.3)           #src
    Wq = H * (Asym \ Qquad)                                                  #src
    Dh = Diagonal(sqrt.(abs.(diag(D))))                                      #src
    Wsc = Dh * (H * (Asym \ transpose(H))) * Dh                              #src
    @assert maximum(abs, Wsc - transpose(Wsc)) < 1e-12                        #src
    @assert maximum(abs, Wq - (H * (Asym \ transpose(H))) * D) < 1e-12         #src
end                                                                           #src

# ```math
# \boxed{\;\bm Q_{\rm quad} = \bm H^{\mathsf T}\bm D,\qquad D_i = -\frac{2\pi w_i}{\cos\theta_i}\;>0
# \ \text{ on the lower hemisphere}\;}
# ```
#
# and ``\bm D^{1/2}\bm H\bm A^{-1}\bm H^{\mathsf T}\bm D^{1/2}`` is symmetric positive
# semi-definite. **Evaluating the pressure work by quadrature at the same nodes where the gap
# is collocated restores conjugacy; projecting it exactly in the Legendre basis does not.**
#
# The reason is that exact projection replaces a load at one node by the Galerkin force of the
# polynomial that interpolates a spike there — a function spread over the whole sphere. The
# gap constraint is local; that force is not; a local constraint and a global force cannot be
# transposes of one another.

# ## 5. Against the shipped solver
#
# Three claims about `DropSolver`, each checked rather than described: `gap_row` is §1's
# Jacobian, `force_column` is §2's forcing, and the assembled compliance is NOT of the form
# ``\bm H^{\mathsf T}\bm D`` for any diagonal `D` — measured as the best-fit residual, so it is
# a statement about how far from conjugate the pairing is, not merely that it differs.

let                                                                          #src
    p = ImpactParams(We = 0.5, Bo = 0.02, Oh = 0.0373, M = 12, K = 1)         #src
    b = DropSolver.basis(p); npc = DropSolver.pc_len(p)                        #src
    ## (a) gap_row reproduces mu*P_l(mu)                                      #src
    th = p.nodes[3]; row, mu_c = DropSolver.gap_row(p, th)                     #src
    @assert isapprox(mu_c, cos(th); atol = 1e-14)                              #src
    for (i, l) in enumerate(b.ls)                                              #src
        @assert isapprox(row[i], cos(th) * Pl(l, cos(th)); atol = 1e-12)         #src
    end                                                                        #src
    ## (b) force_column reproduces -(4pi/(2l+1)) on its own mode and nothing else            #src
    for j in 1:npc                                                             #src
        l = DropSolver.pc_l(j); col = DropSolver.force_column(p, j)             #src
        i = findfirst(==(l), b.ls)                                              #src
        if i === nothing                                                       #src
            @assert all(iszero, col)                                            #src
        else                                                                   #src
            @assert isapprox(col[i], -(4pi / (2l + 1)); atol = 1e-12)            #src
            @assert count(!iszero, col) == 1                                     #src
        end                                                                     #src
    end                                                                        #src
    ## (c) how far the shipped nodal force is from conjugate                   #src
    Vfac = lu(DropSolver.legendre_vandermonde(p))                               #src
    Vinv = Vfac \ Matrix{Float64}(I, npc, npc)                                  #src
    Qn = hcat((DropSolver.force_column(p, j) for j in 1:npc)...) * Vinv          #src
    nn = length(p.nodes)                                                        #src
    H = zeros(nn, DropSolver.ndof(b))                                           #src
    for i in 1:nn; H[i, :] = DropSolver.gap_row(p, p.nodes[i])[1]; end            #src
    low = findall(<(0.0), cos.(p.nodes))                                        #src
    ## best diagonal scaling per node, then the residual it leaves            #src
    worst = 0.0                                                                 #src
    for i in low                                                                #src
        hrow = H[i, :]; qcol = Qn[:, i]                                          #src
        d = dot(hrow, qcol) / dot(hrow, hrow)          # least-squares D_i        #src
        worst = max(worst, norm(qcol - d * hrow) / norm(qcol))                    #src
    end                                                                         #src
    @assert worst > 0.5                                                          #src
    println("shipped force vs best conjugate fit: worst relative residual ",     #src
            round(worst; digits = 3))                                             #src
    ## and therefore the assembled compliance is not symmetric                  #src
    F0 = assemble_newtonian(b, p.Oh)                                             #src
    beta = 1 / p.dt0                                                             #src
    A = beta^2 * F0.M + beta * F0.C + F0.G                                        #src
    S = (H*(A \ Qn))[low, low]                                                    #src
    asym = maximum(abs, S - transpose(S)) / maximum(abs, S)                        #src
    @assert asym > 0.1                                                            #src
    println("assembled compliance relative asymmetry: ", round(asym; digits = 3))  #src
    ## while the conjugate one is symmetric to rounding                          #src
    Wc = (H*(A \ transpose(H)))[low, low]                                          #src
    @assert maximum(abs, Wc - transpose(Wc)) / maximum(abs, Wc) < 1e-12             #src
end                                                                                #src

# ## What this does and does not settle
#
# It settles the mathematics: conjugate forcing gives a symmetric positive semi-definite
# compliance and a convex programme, the shipped pairing is not conjugate, and the repair is to
# evaluate the pressure work at the collocation nodes rather than project it exactly.
#
# It does not settle the physics, which is a separate question with a real answer either way.
# The film pressure is a normal stress on the drop's surface, so radial is the honest direction
# for it; the constraint is a vertical distance to a flat plate, so vertical is the honest
# direction for the multiplier. They agree to ``O((\theta-\pi)^2)`` near the pole and differ at
# ``O(1)`` for a wide patch. Adopting the conjugate pairing therefore changes the model at
# finite truncation — the two converge as ``M\to\infty`` — and the identity that currently ties
# the forcing to Lamb's added mass, ``Q/M = l`` exactly, is a property of the radial-Legendre
# pairing and would have to be re-derived for the nodal one.
#
# So this is a choice, and the cost of each side is now explicit rather than implicit.
