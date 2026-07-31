# Mathematical Preliminaries

Two facts from mathematical physics do essentially all of the geometric work
in the rest of this derivation. Neither is specific to viscous flow -- they
are facts about spherical harmonics and Bessel functions that would show up
in any spherically-symmetric wave or diffusion problem. Both are checked
below against the actual computer-algebra objects verified in
`julia/derivations/reid1960_full_derivation.jl`, not merely asserted.

## Spherical harmonics

A *spherical harmonic* of degree ``l`` and order ``m`` is the function

```math
Y_l^m(\theta,\varphi) = P_l^m(\cos\theta)\,e^{-im\varphi},
```

where ``P_l^m`` is the associated Legendre polynomial, ``\theta`` is the
polar angle from the ``z``-axis, and ``\varphi`` is the azimuthal angle. (We
use the unnormalized convention; the characteristic equation this chapter
builds toward is independent of the normalization choice.)

The key property we use is that spherical harmonics are *eigenfunctions of
the Laplace-Beltrami operator* (the angular part of the Laplacian):

```math
\nabla^2_{\text{angular}}\,Y_l^m = -\frac{l(l+1)}{r^2}\,Y_l^m.
```

This means that for any function of the form ``f(r)\,Y_l^m(\theta,\varphi)``,
the full scalar Laplacian is

```math
\nabla^2\!\left[f(r)\,Y_l^m\right] = \left[\,f'' + \frac{2}{r}f' - \frac{l(l+1)}{r^2}\,f\,\right] Y_l^m,
```

where primes denote ``d/dr``. For axisymmetric problems we take ``m=0``, so
``Y_l^0 \propto P_l(\cos\theta)`` is just the Legendre polynomial.

!!! note "Live check, not just a citation"
    The eigenvalue property above splits into two independent facts: a
    change of variables from ``\theta`` to ``x=\cos\theta`` turns the
    angular Laplacian into the operator ``\frac{d}{dx}\!\left[(1-x^2)\frac{dQ}{dx}\right]``
    for *any* function; and the Legendre polynomials specifically satisfy
    ``(1-x^2)P_l'' - 2xP_l' + l(l+1)P_l = 0``. Both are verified in the
    companion script for symbolic ``l``. Concretely, for ``l=3``, built from
    the same Bonnet-recursion code this repo's other derivations use:

```@eval
using Symbolics, Markdown
@variables x
function legendre_P(l::Int, xv)
    l == 0 && return one(xv)
    l == 1 && return xv
    Pm1, P = one(xv), xv
    for n in 1:(l - 1)
        P, Pm1 = ((2n + 1) * xv * P - n * Pm1) / (n + 1), P
    end
    P
end
P3 = simplify(legendre_P(3, x); expand=true)
Dx = Differential(x)
residual = simplify(expand_derivatives((1 - x^2) * Dx(Dx(P3)) - 2 * x * Dx(P3) + 3 * 4 * P3); expand=true)
Markdown.parse("```math\nP_3(x) = " * Main.pretty_latex(P3) *
    "\n```\nand substituting into Legendre's equation gives, literally, `" * string(residual) * "` -- exactly zero, not approximately.")
```

## Spherical Bessel functions

The ordinary Bessel equation of order ``\nu`` is

```math
w'' + \frac{1}{z}w' + \left(1-\frac{\nu^2}{z^2}\right)w = 0,
```

with solutions ``J_\nu(z)`` (Bessel functions of the first kind) and
``Y_\nu(z)`` (second kind, singular at ``z=0``). The *spherical* Bessel
equation of order ``l``,

```math
v'' + \frac{2}{z}v' + \left(1-\frac{l(l+1)}{z^2}\right)v = 0,
```

arises when separating the Helmholtz equation ``(\nabla^2+k^2)F=0`` in
spherical coordinates. Its regular solution is the spherical Bessel function
of the first kind, ``j_l(z) = \sqrt{\pi/2z}\,J_{l+1/2}(z)``.

!!! tip "A substitution we will need repeatedly"
    Any function ``U(x) = x\,j_l(qx)`` satisfies

    ```math
    U'' - \frac{l(l+1)}{x^2}\,U + q^2 U = 0.
    ```

    This follows by writing ``U=xv(x)`` with ``x=qr/R``: then
    ``U'=v+xv'`` and ``U''=2v'+xv''``, and substituting gives
    ``x\left[v''+\frac{2}{x}v'+q^2v-\frac{l(l+1)}{x^2}v\right]=0``, which is
    exactly the spherical Bessel equation above (at argument ``qx``). This
    identity -- verified symbolically for symbolic ``l`` in the companion
    script, not just at concrete values -- is what makes ``U(x)=x\,j_l(qx)``
    the natural building block for everything that follows, and is why the
    scaling ``u_r \propto U(x)/x^2`` is chosen the way it is once we reach
    the velocity-field ODE.

We will also need the derivative identity. Using the standard recurrence
``j_l'(z) = \frac{l}{z}j_l(z) - j_{l+1}(z)``,

```math
\frac{q\,j_l'(q)}{j_l(q)} = l - q\,\frac{j_{l+1}(q)}{j_l(q)}.
```

Reid defines the ratio ``Q_{l+1/2}(q) \equiv j_{l+1}(q)/j_l(q) = J_{l+3/2}(q)/J_{l+1/2}(q)``,
so this becomes ``qj_l'(q)/j_l(q) = l - qQ_{l+1/2}(q)`` -- the one Bessel
combination the entire characteristic equation reduces to, later in this
chapter.

## The poloidal decomposition

Any divergence-free vector field ``\bm u`` can be uniquely decomposed into a
*toroidal* part (involving only ``\nabla\times(\Psi\hat r)``, with no radial
component) and a *poloidal* part (involving ``\nabla\times\nabla\times(\Phi\hat r)``,
with non-zero radial component). For axisymmetric flow with no swirl, the
toroidal part vanishes, and ``\bm u`` is purely poloidal.

For a poloidal field with angular dependence ``Y_l^m``, the entire velocity
field is determined by its radial component ``u_r``. Incompressibility then
gives ``u_\theta``, and there are no other free components. This is why it
suffices to write down a single scalar ODE for ``U(x)``, where
``u_r \propto U(x)/x^2`` -- which is exactly what the next chapter does.
