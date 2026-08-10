# Appendix: Identities and Standard Results

The chapters use a small number of standard results repeatedly. Each is proved here so that
no page has to interrupt itself, and so that a reader who wants to check one step can find it
in one place rather than in whichever chapter happened to need it first.

Nothing here is original. What is here is complete: every identity the main line of argument
leans on, with enough of a proof that a reader with basic fluid mechanics, Lagrangian mechanics
and linear algebra can verify it without a reference book.

Throughout, ``\mu = \cos\theta``, ``P_l`` is the Legendre polynomial of degree ``l``, and
``L = l(l+1)``.

## A. Legendre polynomials

### A.1 Legendre's equation, in ``\mu`` and in ``\theta``

``P_l`` satisfies

```math
\frac{d}{d\mu}\left[(1-\mu^2)\frac{dP_l}{d\mu}\right] + L\,P_l = 0 .
```

Substituting ``\mu = \cos\theta``, so that ``d/d\mu = -(1/\sin\theta)\,d/d\theta`` and
``1-\mu^2 = \sin^2\theta``, gives the form used whenever an angular Laplacian appears:

```math
\boxed{\;\frac{1}{\sin\theta}\frac{d}{d\theta}\!\left(\sin\theta\,\frac{dP_l}{d\theta}\right)
= -L\,P_l\;}
```

Equivalently ``\partial_\theta^2 P_l + \cot\theta\,\partial_\theta P_l = -L\,P_l``. This is the
identity that makes the trial fields exactly incompressible in *Variational Mechanics*, and the
one that collapses the Gegenbauer function in A.4.

### A.2 Orthogonality, and the factor ``4\pi/(2l+1)``

```math
\int_{-1}^{1} P_l(\mu)P_m(\mu)\,\mathrm{d}\mu = \frac{2}{2l+1}\,\delta_{lm} .
```

On an axisymmetric spherical surface of radius ``R``, ``\mathrm{d}S = R^2\,\mathrm{d}\varphi\,
\mathrm{d}\mu`` and the azimuthal integral contributes ``2\pi``, so

```math
\boxed{\;\oint P_l\,P_m\,\mathrm{d}S = \frac{4\pi R^2}{2l+1}\,\delta_{lm}\;}
```

Every ``4\pi/(2l+1)`` in the corpus is this factor: in the surface energy, in the generalised
force ``Q_l = -\tfrac{4\pi}{2l+1}p_{c,l}``, and in the stiffness matrix.

### A.3 The derivative identity

```math
\boxed{\;(1-\mu^2)\,P_l'(\mu) \;=\; \frac{L}{2l+1}\bigl[P_{l-1}(\mu)-P_{l+1}(\mu)\bigr]\;}
```

This is what turns a product of *derivatives* of Legendre polynomials into a sum of ordinary
Legendre products, and it is the step that lets the selection rule of A.5 cover the
``e_{r\theta}`` contraction as well as the diagonal one. It follows from the standard recurrences
``(2l+1)\mu P_l = (l+1)P_{l+1} + lP_{l-1}`` and ``(1-\mu^2)P_l' = l(P_{l-1}-\mu P_l)`` by
eliminating ``\mu P_l``.

### A.4 The Gegenbauer function

Define

```math
\Gamma_l(\theta) \;=\; -\frac{\sin\theta}{L}\,\frac{dP_l}{d\theta} .
```

Then

```math
\boxed{\;\frac{d\Gamma_l}{d\theta} = \sin\theta\,P_l\;}
```

*Proof.* Differentiate and apply A.1:
``\Gamma_l' = -\tfrac1L\,\tfrac{d}{d\theta}\!\left(\sin\theta\,\tfrac{dP_l}{d\theta}\right)
= -\tfrac1L\,(-L\sin\theta P_l) = \sin\theta P_l``. ∎

This is why the Stokes stream function ``\Psi = f(x)\Gamma_l(\theta)`` produces a radial velocity
proportional to ``P_l``, so that an interior field drives exactly the surface harmonic it is
paired with.

### A.5 Triple products: the triangle and parity rules

```math
\int_{-1}^{1} P_l P_m P_k \,\mathrm{d}\mu \;=\; 0
\qquad\text{unless}\qquad
|l-m| \le k \le l+m
\quad\text{and}\quad
l+m+k \ \text{is even} .
```

The triangle condition is the statement that ``P_lP_m`` is a polynomial of degree ``l+m``
containing only degrees ``\ge |l-m|``, so it is orthogonal to ``P_k`` outside that range. The
parity condition follows from ``P_l(-\mu) = (-1)^lP_l(\mu)``: if ``l+m+k`` is odd the integrand
is odd and the integral over a symmetric interval vanishes. These are the ``m = 0`` case of the
Gaunt coefficients.

## B. Geometry of a nearly spherical surface

Write ``r(\theta) = R\bigl(1 + \zeta_0 + \sum_{l\ge2}\zeta_lP_l(\mu)\bigr)`` with all
``\zeta \ll 1``.

### B.1 Volume conservation forces a second-order shift

Incompressibility fixes the volume, and that cannot be arranged with the ``l\ge2`` amplitudes
alone: they are second-order sources of volume. The mean radius must absorb it.

*Proof.* With ``u = \zeta_0 + \sum_{l\ge2}\zeta_lP_l``,

```math
\mathcal V = \frac13\oint r^3\,\mathrm{d}\Omega
 = \frac{R^3}{3}\oint\bigl(1+3u+3u^2+O(u^3)\bigr)\mathrm{d}\Omega .
```

Now ``\oint\mathrm{d}\Omega = 4\pi``, ``\oint u\,\mathrm{d}\Omega = 4\pi\zeta_0`` since every
``P_{l\ge1}`` integrates to zero, and ``\oint u^2\mathrm{d}\Omega = 4\pi\sum_{l\ge2}
\zeta_l^2/(2l+1)`` to leading order by A.2. Setting ``\mathcal V = \tfrac43\pi R^3`` gives
``3\zeta_0 + 3\sum\zeta_l^2/(2l+1) = 0``, that is

```math
\boxed{\;\zeta_0 = -\sum_{l\ge2}\frac{\zeta_l^2}{2l+1}\;}
```
∎

Omitting this term is not a small error: it changes the surface-energy coefficient from
``(l-1)(l+2)`` to ``l(l+1)+2``, and with it every restoring force and every frequency.

### B.2 First variation of area

For a normal displacement ``\delta\zeta`` of a surface with mean-curvature sum
``\kappa = \nabla_s\!\cdot\bm n``,

```math
\boxed{\;\delta|\partial\Omega| = \oint \kappa\,\delta\zeta\,\mathrm{d}S\;}
```

This is the standard first-variation-of-area formula, and it is the reason capillarity enters
the weak form as the derivative of an energy rather than as a force that has to be modelled. For
a sphere ``\kappa = 2/R``, so a volume-preserving displacement, which satisfies
``\oint\delta\zeta\,\mathrm{d}S = 0`` at first order, costs no area at first order: the leading
capillary energy is quadratic, as B.1 and the surface energy require.

## C. Stokes stream function

For axisymmetric incompressible flow, ``\Psi`` defined by

```math
u_r = \frac{1}{r^2\sin\theta}\frac{\partial\Psi}{\partial\theta} ,
\qquad
u_\theta = -\frac{1}{r\sin\theta}\frac{\partial\Psi}{\partial r}
```

satisfies ``\nabla\cdot\bm u = 0`` identically. With ``\Psi = f(r)\Gamma_l(\theta)`` and A.4,

```math
u_r = \frac{f}{r^2}P_l ,
\qquad
u_\theta = \frac{f'}{rL}\frac{dP_l}{d\theta} ,
```

which are the components used throughout.

### C.1 The ``E^2`` operator, and vorticity

With ``E^2 = \partial_r^2 + \frac{\sin\theta}{r^2}\partial_\theta\!\left(\frac{1}{\sin\theta}
\partial_\theta\right)``, the angular part acting on ``\Gamma_l`` returns ``-L\,\Gamma_l``, by
A.4 followed by A.1. Hence

```math
E^2\Psi = \left(f'' - \frac{L}{r^2}f\right)\Gamma_l ,
\qquad
\boxed{\;\omega_\varphi = -\frac{E^2\Psi}{r\sin\theta}
 = \frac{1}{rL}\left(f'' - \frac{L}{r^2}f\right)\frac{dP_l}{d\theta}\;}
```

The ``\sin\theta`` cancels, so the vorticity is finite on the axis. A radial function with
``f'' = Lf/r^2`` carries none: ``f = r^{l+1}`` is such a function, which is why the first radial
trial function is irrotational and all vorticity comes from the rest.

## D. The reduced oscillator

### D.1 Vieta, and what ``\lambda_l`` and ``\omega_l^2`` mean

Under the convention ``A_l \propto e^{-\sigma t}``, the equation
``\ddot A_l + 2\lambda_l\dot A_l + \omega_l^2A_l = 0`` has characteristic equation

```math
\sigma^2 - 2\lambda_l\sigma + \omega_l^2 = 0 ,
\qquad\text{so}\qquad
\lambda_l = \tfrac12(\sigma_1+\sigma_2) , \qquad \omega_l^2 = \sigma_1\sigma_2 .
```

The **minus** sign in the middle term is a consequence of the ``e^{-\sigma t}`` convention, not a
typographical slip; with ``e^{+\sigma t}`` it would be a plus and the decay rates would come out
negative. Given any root pair, these two formulas define a unit-mass oscillator reproducing it.

## E. Complementarity

### E.1 The quadratic programme and the LCP are the same statement

For symmetric ``\bm A``, the problem ``\min_{\bm p\ge0}\tfrac12\bm p^{\mathsf T}\bm A\bm p +
\bm b^{\mathsf T}\bm p`` has Karush-Kuhn-Tucker conditions: attach a multiplier ``\bm s\ge0`` to
``\bm p\ge0``, then

```math
\underbrace{\bm A\bm p + \bm b - \bm s = \bm 0}_{\text{stationarity}},
\qquad \bm p\ge0,\ \bm s\ge0, \qquad s_ip_i = 0 .
```

Stationarity identifies ``\bm s = \bm A\bm p + \bm b``, which is the clearance ``\bm h``.
Substituting gives ``\bm h\ge0``, ``\bm p\ge0``, ``h_ip_i=0``: the complementarity problem
unchanged. **The multiplier of the pressure's sign constraint is the clearance.**

### E.2 P-matrices

``\bm A`` is a **P-matrix** if every principal minor is positive. For a P-matrix, the linear
complementarity problem has **exactly one** solution for every ``\bm b`` (Samelson, Thrall and
Wesler; see Cottle, Pang and Stone, *The Linear Complementarity Problem*, Thm 3.3.7), and
principal pivoting terminates finitely under a least-index rule.

A cheap sufficient condition, which does not require symmetry: if the symmetric part
``\tfrac12(\bm A + \bm A^{\mathsf T})`` is positive definite, then ``\bm A`` is a P-matrix.

*Proof, in two steps.* First, ``\bm x^{\mathsf T}\bm A\bm x =
\bm x^{\mathsf T}\tfrac12(\bm A+\bm A^{\mathsf T})\bm x > 0`` for every ``\bm x \ne \bm 0``,
since the antisymmetric part contributes nothing to a quadratic form. Second, any matrix with
that property has positive determinant: its eigenvalues are either real and positive (take
``\bm x`` the eigenvector) or occur in complex-conjugate pairs, whose product
``|\lambda|^2`` is positive, so the determinant, being the product of all of them, is
positive. Finally, every principal submatrix inherits the property, because restricting
``\bm x`` to the corresponding coordinate subspace is a special case of the same inequality.
Hence every principal minor is positive. ∎

This is the certificate *Contact* uses, because it is computable at any truncation
while checking ``2^n-1`` minors is not.

## F. Generalized Newtonian dissipation

### F.1 The dissipation potential

For ``\eta = \eta(\dot\gamma)`` with ``\dot\gamma = \sqrt{2\,\bm e\!:\!\bm e}``, the generalised
force whose divergence is ``\bm\tau = 2\eta(\dot\gamma)\bm e`` is the gradient of

```math
\boxed{\;\mathcal R = \int_\Omega W(\dot\gamma)\,\mathrm{d}V ,
\qquad W(\dot\gamma) = \int_0^{\dot\gamma}\!\eta(s)\,s\,\mathrm{d}s \;}
```

*Proof.* ``\partial\dot\gamma/\partial\dot\xi_a = 2\,\bm e\!:\!\bm e^{(a)}/\dot\gamma`` from
``\dot\gamma^2 = 2\bm e\!:\!\bm e``, so
``\partial W/\partial\dot\xi_a = \eta\dot\gamma\cdot 2\bm e\!:\!\bm e^{(a)}/\dot\gamma
= 2\eta\,\bm e\!:\!\bm e^{(a)}``. ∎

It is **not** the gradient of ``\int\eta(\dot\gamma)\,\bm e\!:\!\bm e\,\mathrm{d}V``, which
carries an extra ``\dot\gamma\eta'`` and is wrong by about thirty per cent on a Carreau fluid
with ``\lambda_c = 3``, ``a = 0.75``, ``n = 0.25``, ``r_\infty = 10^{-3}``.
For constant ``\eta`` the two coincide, since ``W = \tfrac12\eta\dot\gamma^2 = \eta\,\bm e\!:\!\bm e``.
