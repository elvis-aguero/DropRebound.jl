# Glossary of Symbols

Every symbol used in the chapters, with its meaning, its size where it is an array, and where
it is introduced. Sizes are given in terms of the truncation ``M`` and the number of radial
trial functions per mode ``K``.

A reader who meets an unfamiliar symbol mid-chapter should be able to resolve it here without
hunting backwards through the text. Where the same letter carries two meanings, that is recorded
rather than quietly tolerated: see [Collisions](@ref) at the end.

## Fields and coordinates

| symbol | size | meaning |
|---|---|---|
| ``\bm u`` | field | velocity |
| ``\bm v`` | field | test velocity field: incompressible, compatible with the kinematic condition |
| ``p`` | field | pressure inside the drop |
| ``\bm e`` | tensor | rate of strain, ``\tfrac12(\nabla\bm u + \nabla\bm u^{\mathsf T})`` |
| ``\bm\Sigma`` | tensor | Cauchy stress, ``-p\bm I + 2\eta\bm e`` |
| ``\dot\gamma`` | scalar | shear rate, the invariant ``\sqrt{2\,\bm e\!:\!\bm e}`` |
| ``\Psi`` | field | Stokes stream function |
| ``\zeta(\theta,t)`` | field | radial surface displacement from the sphere |
| ``\zeta_l`` | scalar | its Legendre amplitude; the boundary trace ``\chi_l(1)`` |
| ``\chi_l(x,t)`` | field | radial profile of the interior displacement for mode ``l`` |
| ``\psi_l`` | field | its rate, ``\dot\chi_l`` |
| ``\bm\xi`` | ``(M-1)K`` | **the state**: interior amplitudes ``a_{l,k}``, one per mode ``l\in\{2..M\}`` and radial function ``k``. Not the surface |
| ``z``, ``v`` | scalars | centre-of-mass height and velocity |
| ``\theta``, ``\mu`` | | polar angle and ``\cos\theta``; ``\theta=\pi`` is the contact pole |
| ``x`` | | radius scaled by the drop radius, ``r/R`` |

## Material and dimensionless groups

| symbol | meaning |
|---|---|
| ``\rho``, ``\gamma``, ``\eta`` | density, **surface tension**, dynamic viscosity |
| ``R``, ``v_0`` | undeformed radius, impact speed |
| ``\tau_c`` | inertio-capillary time ``\sqrt{\rho R^3/\gamma}``. A time *scale*, not a period: the ``l=2`` period is ``2.22\,\tau_c`` |
| ``\mathrm{We}`` | Weber, ``\rho v_0^2R/\gamma`` |
| ``\mathrm{Oh}`` | Ohnesorge, ``\eta/\sqrt{\rho\gamma R}``, built from the **zero-shear** viscosity |
| ``\mathrm{Bo}`` | Bond, ``\rho gR^2/\gamma`` |
| ``\eta_0``, ``\eta_\infty`` | zero-shear and infinite-shear viscosity plateaus |
| ``r_\infty`` | ``\eta_\infty/\eta_0``, the solver's `eta_inf_ratio` |
| ``\varepsilon_{ST}`` | ``1-r_\infty``, the fraction of viscosity the fluid can shed. **Not** `eta_inf_ratio` |
| ``\lambda_c``, ``a``, ``n`` | Carreau-Yasuda time constant (non-dimensional), shape exponent, power-law index |
| ``K``, ``m`` | Cross time constant (seconds) and exponent, before conversion |
| ``\mathrm{De}_1``, ``\beta_s`` | Oldroyd-B relaxation time ``\lambda_1/\tau_c`` (mode-independent) and solvent fraction |

## The reduced system

| symbol | size | meaning |
|---|---|---|
| ``T``, ``V``, ``\mathcal R`` | scalars | kinetic energy, potential (surface + gravity), dissipation potential |
| ``W(\dot\gamma)`` | scalar | ``\int_0^{\dot\gamma}\eta(s)s\,\mathrm{d}s``, the integrand of ``\mathcal R`` |
| ``\Phi`` | scalar | total viscous dissipation rate, ``2\mathcal R`` |
| ``\bm M`` | ``((M-1)K)^2`` | mass matrix. **The boldface is this matrix, not the truncation** ``M`` |
| ``\bm C``, ``K_{ab}`` | ``((M-1)K)^2`` | damping matrix, ``\mathrm{Oh}\!\int2\eta\,\bm e^{(a)}\!:\!\bm e^{(b)}\mathrm{d}V`` |
| ``\bm G`` | ``((M-1)K)^2`` | stiffness; rank one per mode |
| ``\bm Q``, ``Q_a`` | ``(M-1)K`` | generalised force of the film pressure |
| ``\bm u^{(a)}``, ``\bm e^{(a)}`` | fields | velocity and strain of trial function ``a`` |
| ``\phi_k(x)`` | | ``k``-th radial trial function; ``\phi_k(1)`` its surface trace |
| ``\lambda_l``, ``\omega_l^2`` | scalars | damping and squared frequency of the reduced mode ``l`` |
| ``\sigma`` | scalar | complex decay rate, under ``e^{-\sigma t}`` |
| ``q``, ``\alpha`` | scalars | ``q^2=\sigma R^2/\nu``, ``\alpha^2=\sigma_{l;0}R^2/\nu`` |
| ``\mathcal A_l``, ``D_l`` | scalars | Molaček and Bush's inertia and dissipation coefficients |
| ``\mathcal R_{lm}`` | block | interior operator; the ``(l,m)`` block of ``K_{ab}/\mathrm{Oh}`` |
| ``\eta_k`` | scalar | ``k``-th Legendre harmonic of the viscosity field. **Not** a viscosity plateau |
| ``L_\eta`` | integer | angular bandwidth of the viscosity field |
| ``G^k_{lm}``, ``H^k_{lm}`` | scalars | the two families of angular integrals in the traction-route assembly |
| ``A^{(i)}_{lm}``, ``B^{(i)}_{lm}`` | scalars | the radial integrals pairing with ``G`` and ``H`` |
| ``\mathcal D^{(2)}``, ``\mathcal D^{(1)}`` | ``(M-1)^2`` | damping and stiffness of the eliminated surface system; generalise ``2\bm\Lambda`` and ``\bm\Omega`` |
| ``\bm\Lambda``, ``\bm\Omega`` | diagonal | ``\mathrm{diag}(\lambda_l)`` and ``\mathrm{diag}(\omega_l^2)`` of the Newtonian problem |
| ``\mathcal T(L)`` | scalar | the truncation error metric of the bandwidth measurement |

## Contact

| symbol | size | meaning |
|---|---|---|
| ``p_c(\theta,t)`` | field | film pressure on the lower surface |
| ``p_{c,l}`` | scalar | its Legendre coefficient |
| ``\bm p`` | ``n`` | film pressure at the retained nodes; the LCP unknown |
| ``\bm\lambda`` | ``n`` | multipliers conjugate to the nodal clearances; the unknown under `force_mode = :nodal`. Load-like, but not the vertical force carried by the pressure |
| ``h``, ``\bm h``, ``\bm g`` | ``n`` | clearance between surface and substrate |
| ``n`` | integer | number of retained contact nodes, ``\lfloor M/2\rfloor+1`` |
| ``\bm V`` | ``(M{+}1)^2`` | Legendre Vandermonde, ``V_{ij}=P_{j-1}(\mu_i)`` |
| ``\bm H`` | ``n\times(M-1)K`` | clearance Jacobian, ``H_{i,(l,k)}=\cos\theta_iP_l(\mu_i)\phi_k(1)``, restricted to the retained nodes |
| ``\bm b_0`` | ``n`` | undeformed clearance along each ray, ``\cos\theta_i``, restricted |
| ``\bm A`` | ``((M-1)K)^2`` | step operator ``\beta^2\bm M+\beta\bm C+\bm G`` |
| ``\bm Q_n`` | ``(M-1)K\times n`` | nodal pressures to generalised forces, ``\bm Q_{\text{modal}}\bm V^{-1}\bm E`` |
| ``\bm A_c`` | ``n^2`` | compliance, ``\partial\bm h/\partial\bm p`` at frozen history |
| ``\bm W`` | ``n^2`` | its symmetric counterpart under nodal forcing |
| ``\bm b`` | ``n`` | clearance the step would produce with no contact force |
| ``\beta`` | scalar | ``c_0/\Delta t``; BDF2 has ``c_0=(1+2r)/(1+r)``, BDF1 has ``c_0=1`` |
| ``\bm 1`` | ``n`` | column of ones, ``\partial\bm h/\partial z`` |
| ``D_i`` | scalar | ``-2\pi w_i/\cos\theta_i``, the diagonal relating multiplier to pressure |
| ``\bm E`` | ``(M{+}1)\times n`` | selects the retained contact nodes |
| ``w_i`` | scalar | quadrature weight at node ``i`` |

## Collisions

The same letter carries more than one meaning in places. Each is unavoidable without breaking
with the literature, so each is flagged where it occurs.

| letter | meanings | how to tell |
|---|---|---|
| ``M`` | truncation degree; mass matrix | the matrix is boldface ``\bm M`` |
| ``K`` | radial functions per mode; Cross time constant; damping ``K_{ab}`` | context: an integer count, a time in seconds, or a subscripted matrix |
| ``\eta_0`` | zero-shear plateau; ``k=0`` harmonic of the viscosity field | the harmonic is written ``\eta_k`` with a running index |
| ``a`` | trial-function index; Carreau-Yasuda shape exponent | index in ``\bm u^{(a)}``, exponent in ``(\lambda_c\dot\gamma)^a`` |
| ``\lambda`` | decay rate ``\lambda_l``; Carreau time ``\lambda_c``; multiplier ``\bm\lambda`` | subscript ``l``, subscript ``c``, or boldface |
| ``\nu`` | kinematic viscosity; Bessel order in Appendix D.2 of *The Free Viscous Drop* | the Bessel order appears only inside that subsection |
| ``\bm b`` | film-pressure forcing (shear-thinning pages); free-flight clearance (*Contact*) | context: a forcing term in an oscillator, or a clearance |
| ``p`` | interior pressure; the exponent ``(n-1)/a`` on the closures page | the exponent is always written in an exponent position |
| ``\sigma`` | decay rate | surface tension is ``\gamma`` everywhere in this corpus |
| ``\gamma``, ``\dot\gamma`` | surface tension; shear rate | the dot |
| ``A_l`` | modal amplitude (reduced model); Molaček-Bush inertia | the latter is script, ``\mathcal A_l`` |
