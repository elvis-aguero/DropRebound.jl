royalsocietypublishing.org/journal/rspa

## Research


**Cite this article:** Zrnić D, Brenn G. 2024
Weakly nonlinear shape oscillations of
viscoelastic drops. _Proc. R. Soc. A_ **480** :
20230887.
[https://doi.org/10.1098/rspa.2023.0887](https://doi.org/10.1098/rspa.2023.0887)


Received: 27 January 2024
Accepted: 3 July 2024


**Subject Category:**
Physics


**Subject Areas:**
fluid mechanics


**Keywords:**
drop-shape oscillations, weakly nonlinear
analysis, viscoelastic liquid, quasi-periodic
motion, onset of aperiodic motion


**Author for correspondence:**
Günter Brenn
e-mail: guenter.brenn@tugraz.at


Electronic supplementary material is available
[online at https://doi.org/10.6084/](https://doi.org/10.6084/m9.figshare.c.7399609)
[m9.figshare.c.7399609.](https://doi.org/10.6084/m9.figshare.c.7399609)


# Weakly nonlinear shape oscillations of viscoelastic drops

#### Dino Zrnić and Günter Brenn

Institute of Fluid Mechanics and Heat Transfer, Inffeldgasse 25/F, Graz 8010,
Austria


[GB, 0000-0001-7576-0790](http://orcid.org/0000-0001-7576-0790)


Axisymmetric shape oscillations of a viscoelastic drop
in a vacuum are studied by a weakly nonlinear
analysis. The two-lobed initial drop deformation
mode is studied. The Oldroyd-B model is used for
characterizing the drop liquid rheological behaviour.
The equations of motion and the solutions up to
second order are developed, where elastic effects
appear. Solutions of the characteristic equation are
validated against the decay rate and frequency
of damped shape oscillations of polymer solution
drops in an acoustic levitator. The theory shows
enhancing or dampening of the nonlinear behaviour
and enhanced mode coupling, as compared with
the Newtonian case. The study reveals an excess
time in the prolate shape and a frequency change,
together with a quasi-periodicity of the oscillations.
The Fourier power spectra of traces of the drop
north pole position in time show mode coupling
as a nonlinear effect. At moderate stress-relaxation
Deborah number, the resultant drop motion for
large Ohnesorge number, suggesting aperiodic drop
behaviour, may nonetheless be oscillatory, where the
drop elasticity is owing to the liquid bulk elasticity
instead of surface tension. A method is developed
to predict the oscillatory or aperiodic behaviour
of a drop of a given size from a rheologically
characterized viscoelastic Oldroyd-B liquid.

### 1. Introduction

For more than one and a half centuries, shape oscillations of drops have been scientifically investigated. The
motivations are fundamental interest and the influence
on mass, momentum and energy transport across the


© 2024 The Authors. Published by the Royal Society under the terms of the
Creative Commons Attribution License [http://creativecommons.org/licenses/](http://creativecommons.org/licenses/by/4.0/)
[by/4.0/, which permits unrestricted use, provided the original author and](http://creativecommons.org/licenses/by/4.0/)
source are credited.


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

drop surface. Rayleigh presented in the appendix to his paper on the capillary phenomena of
jets an analysis of linear shape oscillations of an inviscid drop in a vacuum around a spherical
equilibrium state [1].
Rayleigh obtained an equation for the angular oscillation frequency of the drop deformed for
a given deformation mode assuming natural values greater than unity corresponding to the
number of lobes on the drop surface. A generalization for viscous drop liquid and an ambientm
medium characterized by its density is owing to Lamb [2,3]. The analysis predicts the Ohne
sorge number 0 = /( ) [1/2] of the drop (dynamic viscosity, surface tension, drop radius
and density Oℎ) at the onset of aperiodic motion. Chandrasekhar analysed small-amplitudeμ σaρ μ σ
ashape oscillations of a viscous, self-gravitating globe in a vacuum, developing the characteristicρ
equation of the globe to determine the complex angular oscillation frequency [4]. Reid showed
that the results of this analysis are equivalent for the oscillations of a drop if the restoring
gravity force is replaced by the force owing to surface tension [5]. Another generalizing step
in the analysis is owing to Miller & Scriven [6] and Basaran [7], who accounted for both the
viscous and the inertial influences from the ambient medium hosting the viscous oscillating
drop. In the work [6], the interfacial dilatational viscosity and elasticity of the drop were
addressed for the first time. The important aspect of the initiation of the oscillations was
studied by Prosperetti [8]. He analysed the drop-shape oscillations, formulating an initial-value
problem. The solutions show that the normal-mode approach may predict aperiodic drop
motion in a range of Ohnesorge numbers, where the solutions of the initial-value problem
show the motion to turn into periodic in time [8]. The most important results from analyses of
drop-shape oscillations are the angular frequency and decay rate of the oscillations, together
with the time-dependent shapes of the drops.
One of the first theoretical analyses of small-amplitude axisymmetric shape oscillations
of a viscoelastic liquid drop in microgravity is owing to Khismatullin & Nadim [9]. The
authors derived the characteristic equation for a viscoelastic drop and found that the equation
has an infinite number of roots, depending critically on the values of the stress relaxation
and deformation retardation times, as well as the surface tension. Asymptotic analyses of
the characteristic equation in the low- and high-viscosity limits and for low, moderate and
large elasticities are performed, showing the periodic and damped or aperiodic behaviours
of the drops. A more recent study numerically analysed freely decaying shape oscillations of
viscoelastic liquid drops in another viscoelastic host liquid, using a Navier–Stokes/Oldroyd-B
immersed boundary technique [10]. For low viscosity, the fluid viscoelasticity was found to
modulate the drop oscillatory shape relaxation, raising the oscillation frequency and reducing
the decay rate, when the fluid relaxation time is above a critical value. The oscillation behaviour
of small drops of a viscoelastic, shear-thinning aqueous solution applied in drop-on-demand
inkjet printing was experimentally investigated by Hoath _et al_ . [11]. The effects of surface
charge and the electrical properties of the viscoelastic liquids on drop oscillation and instability
characteristics were investigated by Li _et al_ . [12].
Properties of drop-shape oscillations and their dependency on the physical properties of
the drop liquid may be used for liquid material characterization. Using this oscillating drop
method, the liquid surface tension against the ambient gas [13–15] and the interfacial tension
between immiscible liquids [16], as well as the liquid dynamic viscosity [17–19], are measured.
For drops from viscoelastic liquids, including surfactant solutions, surface rheological properties were measured [20–24]. The shape oscillations of a drop in the presence of surfactants
were studied by Lu & Apfel [25]. The authors showed that, with interfacial Gibbs elasticity and
viscosity, shape oscillations are dampened more strongly than with just an interfacial tension
as the interfacial property. Investigations of the dynamic surface tension of a drop against
its ambient gaseous host medium, using a growing-drop technique, are owing to Zhang _et_
_al_ . [26]. Shape oscillations of the pendant drop prevent the measurements in the first 20 ms
after the formation of the interface. Kovalchuk _et al_ . [27] provided a comprehensive review of
oscillating drop and bubble techniques as published before the year 2001. The paper Ravera

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```



_et al_ . [28] updated this with a focus on interfacial dilational rheology by the oscillating bubble
and drop methods. The conclusion was that the accuracy achieved with these methods allows
for a reliable measurement of the dilational viscoelasticity of interfaces as a function of the
frequency. For Newtonian drops, the dynamic surface tension, and for shear-thinning drops,
both the transient shear viscosity and the dynamic surface tension were measured, using the
oscillating drop method in Yang _et al_ . [29]. Measurements were performed on free-falling,
oscillating drops. The authors also showed that the dynamic viscosity greatly affects the liquid
behaviour at high shear.
The study by Brenn & Plohl [30] used the damped drop-shape oscillation method for
measuring the deformation retardation time of polymeric liquids, thus introducing a drop-oscillation-based experimental method for determining this polymeric timescale. In their studies,
they relied on Khismatullin & Nadim [9] and Brenn & Teichtmeister [31] in theoretically
analysing small-amplitude axisymmetric shape oscillations of viscoelastic drops in a gas using
the Jeffreys model as the rheological constitutive equation (RCE) of the liquid. The study
showed that the measured values of the deformation retardation time corresponding to the
damped shape oscillations investigated deviate strongly from the values predicted by the
viscoelastic stress splitting method, which is often used in simulations of viscoelastic liquid
flow. Ponce-Torres _et al_ . [32] studied the break-up of surfactant-laden drops, concluding from
shape oscillations of satellite droplets on their surfactant content. Drop-shape oscillations
induced by acoustic drop levitation were used for measuring the rheological material behaviour of blood in the work [33]. The technique makes it possible to assess blood viscosity,
including its changes, and can be used for monitoring blood rheology in sickle cell and other
haematological diseases. Most recently, the frequency-response diagram of a soft viscoelastic
drop driven to shape oscillations was proposed as the basis for a rheometric experiment, using
either an ultrasonically levitated drop or to be conducted in a microgravity environment [34].
The drop rheological behaviour was characterized using the Kelvin–Voigt and Maxwell models
for the soft gels and polymeric fluids, respectively. The frequency-response diagram from the
theoretical analysis is the basis for the use of the drop-oscillation experimental set-up as a
‘drop vibration rheometer’, which is proposed for future experiments. These findings show
that the literature on nonlinear shape oscillations for viscoelastic drops is sparse. For further
understanding of the measured nonlinear drop-shape oscillations described above, a reliable
theoretical method is needed to describe the time evolution of the drop shapes for viscoelastic
liquids.
In the present work, we analyse in detail the time behaviour of nonlinear shape oscillations
of viscoelastic drops. Steered by our earlier work [35,36], we study nonlinear effects in shape
oscillations of viscoelastic drops in a vacuum, which are initially deformed by the two-lobed
mode = 2. The basis is a weakly nonlinear analysis of the oscillations, which we restrict to
axisymmetric drop shape. The analysis is carried to the second order of approximation sincem
elastic effects appear at that order. This kind of analysis was not done in the literature before.
The viscoelasticity of the drop liquid is characterized by the Oldroyd-B RCE. The analysis is
based on series expansions of the flow field variables and the drop shape with respect to a
deformation parameter. The nonlinear phenomena, which consist of an asymmetry of the times
spent in different states of drop deformation and an oscillation frequency change for varying
initial drop deformation, are revealed, covering a wide range of drop Ohnesorge numbers as
found in many technical applications. The frequencies of the oscillatory motions are determined
by the characteristic equation of the drop. The quasi-periodicity of the motion, which results
from mode coupling, is unveiled. The time dependency of the solutions for the flow field
variables and the drop shape, at a given approximation order, are governed by frequencies in
products of solutions of the respective lower-order approximations. The results are relevant
for mass, momentum and energy transport across the oscillating drop surface. Applications in
production technologies, such as inkjet printing [37], containerless materials processing based
on single drops [38], spray cooling [39,40], in-air microfluidics [41,42] and food production, will

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```



benefit from the results, since they allow for a more accurate prediction of evaporation rates,
drag coefficients and heating or cooling rates of the drops.
In the following section, the problem is formulated, and the equations of motion with their
boundary and initial conditions are derived up to the second order of approximation. In §3, the
solutions of the governing equations for the respective orders of approximation are presented.
Results from the analysis are presented and discussed in §4. In §5, the conclusions from the
work are summarized.

### 2. Formulation of the problem

We study the weakly nonlinear shape oscillations of a viscoelastic liquid drop as sketched in
figure 1. The drop is axisymmetric with respect to the azimuthal coordinate _ϕ_ of the spherical
coordinate system. The liquid is treated as incompressible, and its viscoelastic rheological
behaviour is represented by the Oldroyd-B model. The dynamic influence from an ambient
medium is neglected, i.e. we analyse the drop motion in a vacuum. Body forces are not





medium is neglected, i.e. we analyse the drop motion in a vacuum. Body forces are not

accounted for, since in many processes with droplets the Froude number = ( / ) [1/2] is large
enough to allow for this neglect. The problem is formulated in spherical coordinates (Fr σ ρga [2],, _ϕ_ ) to
account for its geometry. r θ
The equations of motion with their initial and boundary conditions are non-dimensionalized



The equations of motion with their initial and boundary conditions are non-dimensionalized

with the undeformed drop radius, the capillary time scale ( / ) [1/2], the capillary pressure /

and the extra stress ( / [3] ) [1/2] for length, time, pressure and extra stress, respectively. Here, a ρa [3] σ σ a



with the undeformed drop radius and the extra stress 0( / ) [1/2] for length, time, pressure and extra stress, respectively. Here,, the capillary time scale ( [3] / ) [1/2], the capillary pressure /

is the liquid density, μ the air–liquid interfacial tension and σ ρa [3] 0 the zero-shear viscosity of theρ
drop liquid. The drop surface is described as the place where σ μ (, ) = 1 + (, ) (cf. figure 1),
with the non-dimensional deformation against the undisturbed spherical shape.rs θ t η θ t
For the problem at hand, the equation of continuity and the momentum balances in theη
radial and polar angular directions, and, read




[= 0,]



∂
∂ur [+][ u][r]
~~t~~



∂



Oℎ






[+]



∂ 2




[+]



2

[−] [∂][p]

[2]





∂
∂uθ [+][ u][r]
~~t~~



∂
∂u ~~r~~ θ [+] [u] ~~r~~ [θ]



Oℎ



∂
=     - [1]

[2]



∂





(2.2)




[+]



(2.3)



where 0 = 0/( ) [1/2] is the Ohnesorge number, the characteristic dimensionless parameter

by the non-dimensional RCE of the Oldroyd-B fluid



+ 1▽ = 2 + 2▽, (2.4)
**τ** De **τ** **D** De **D**



,



where and are the extra-stress and rate-of-deformation tensors, respectively, the latter
**τ** **D**



where and

defined with the velocity gradient tensor as
→



as

[∇]


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

_z_





![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-4-1.png)







**Figure 1.** Sketch of the geometry of a liquid drop under deformation at mode 2 (adapted from Smuda _et al_ . [43]).



The Deborah numbers 1 = 1( / ) [1/2] and 2 = 2( / ) [1/2] represent the non-dimensional
stress relaxation and deformation retardation times of the liquid, respectively. The upper-con-De λ σ ρa [3] De λ σ ρa [3]



vected derivative ▽
**A**



of a tensor is given as
**A**




- = - T, (2.6)
**A** [d] **[A]** [−] →v→⋅ [∇] **A** **A** ⋅∇→v→
~~dt~~



where / is the material derivative. The alternative of modelling the drop liquid as a
second-order fluid would represent the time behaviour of the extra stress differently.d **A** dt
The equations of change are solved by applying initial conditions and boundary conditions.
The material rate of deformation of the drop surface is equal to the radial velocity component of
the deformed drop surface. The kinematic boundary condition therefore reads



∂

ur = [dη] [=] [∂] ∂ [η] [+] [u][θ] ∂η at r = 1 + η . (2.7)
~~dt~~ ~~t~~ ~~r~~ ~~θ~~



∂
at = 1 + .
∂η
r η
~~θ~~



The shear stress at the surface of a drop in a vacuum vanishes since momentum transfer across
the drop surface is not possible. The related dynamic boundary condition therefore reads



× = [→] 0 at = 1 + . (2.8)

[→] n ⋅ **τ** [→] n r η



The unit normal vector on the drop surface, pointing outward, is determined by

[→] n



1

with = −1 − (, ) = 0.
| | →F∇ F r η t θ
~~→F∇~~



1
= with = −1 − (, ) = 0. (2.9)

[→] n | | →F∇ F r η t θ



The tensor of deformation-induced stress in equation (2.8) is formulated for an incompressible Oldroyd-B fluid. Furthermore, the normal stress on the drop surface, which consists of
pressure, viscoelastic normal stress and capillary stress, vanishes for the same reason as shear
stress. The normal-stress boundary condition therefore is obtained as




- + 0 + = 0 at = 1 + . (2.10)
p Oℎ [→] n ⋅ **τ** ⋅n [→] →⋅n→∇ r η



For the deformed drop, the divergence of the normal unit vector on the drop surface in this
equation is obtained as [35]


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-5-1.png)











![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-5-2.png)



at = 1 + .







In a weakly nonlinear analysis of the oscillating drop problem, the flow field variables and the
drop surface shape are formulated as series expansions with respect to a small deformation
parameter, denoted 0. For details, the reader is referred to our earlier work [36]. The formulation, e.g. for the extra stress tensor, readsη



= 1 0 + 2 02 + … . (2.12)
**τ** **τ** η **τ** η



1 0 2 0

The series expansions converge for sufficiently small deformation parameter 0. In the weakly
nonlinear analysis, the boundary conditions are formulated on the deformed drop surface. Theη
values of the field variables on the deformed drop surface are obtained from Taylor expansions,
such as, for example, for
**τ**





Applying this framework to the equations of motion, the RCE and the boundary conditions and
representing the flow properties and their derivatives as given for the example of in (2.13),
sets of first- and second-order equations of motion and boundary conditions are obtained, **τ**
which consist of terms with the deformation parameter 0 to the first and second powers,
respectively. All the terms together, which are linear in the parameter, represent the first-orderη
equations, and the terms with quadratic 0 are the second-order equations.
The analysis assumes that the non-equilibrium starting the drop motion consists of an initialη
drop surface that differs from the spherical shape. This drop shape at the beginning of the
motion is described by a Legendre polynomial of degree with the amplitude 0. Furthermore,
it is assumed that the surface is initially at rest. The mode number m equals the number of lobesη
of the drop surface along the polar angle . When calculating the volume of the deformed drop,m
the deformed drop shape is obtained as θ



(, 0) = 1 + (, 0) = 1 + 3 02

2 η

rs θ η θ

~~m~~



3 02 3

[1]
2 η + 1 [+] ~~2~~ [η][0]
~~m~~



1
(cos ) [3] (cos )
Pm θ d θ



−1/3
+ 0 (cos )
η Pm θ



−1

η Pm θ η
~~m~~



−1



1
(cos ) [3] (cos ) … .
Pm θ d θ ∓





(2.14)



Each term with 0 to a given power represents the factor ensuring volume conservation at

frequencies for the oscillation modes zero and unity are needed, which are not readily found in
the characteristic equation of the drop. This problem can be solved using the fact that the centre
of drop mass does not move away from its initial position



1



1



1 + 4 0 1 + 02 6 12 + 4 2 + 03 4 13 + 12 1 2 + 4 3 + …,
η η η η η η η η η η xdx



where = cos . The physical reason is that there is no resultant force acting on the drop that


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

##### (a) First-order equations





The equations for the two different approximation orders are presented here without details
of their derivation, in order not to repeat details given in previous publications [35,36]. The
first-order equations of change read




[+]




[= 0,]







Oℎ





∂ 1 = - [1]
∂uθ
~~t~~ ~~r~~



∂ 1




[+]



(2.17)


(2.18)



Oℎ




[+]




[τ][rθ][1] τ



The RCE of the first order becomes



∂ 1
∂ **τ** [= 2]

**[D]** [1][ +][ De][2]
~~t~~





**τ** De **[D]** [1][ +][ De][2]



The first-order boundary conditions at = 1 read
r



1 = [∂][η][1] 1 = 0, (2.20)
∂ [,]
r rθ
u τ
~~t~~



p Oℎ τ η



and the first-order initial conditions are



1(, 0) = (cos ), ∂∂η1 [(][θ][, 0) = 0.] (2.22)
η θ Pm θ
~~t~~


##### (b) Second-order equations



The second-order equations of change read




[= 0,]



∂ 2






[+]





2
1
θ
u




[+]




[+]



2
1

~~r~~



∂ 1
~~θ~~



(2.24)


(2.25)



∂ 2



∂ 1
~~r~~ ~~r~~



∂ 2
∂p [=]
~~θ~~




[τ][rθ][2] τ





∂ 1 .
~~θ~~ ~~r~~




- 1
r
u



The RCE of the second order is


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

2 + 1
**τ** De



∂ 2
∂ **τ** [−2]
~~t~~









(2.26)





The second-order boundary conditions at = 1 are
r



∂ 1
∂ur [,]
~~r~~



∂ 1
∂η [−] [η][1]
~~θ~~



r2 − [∂] ∂ [η][2] [=] ∂∂η1 [−] [η][1] ∂∂ur1 [,] (2.27)
u [u][θ][1]
~~t~~ ~~r~~ ~~θ~~ ~~r~~



τ τ



τ η τ τ







∂ 1
∂η [cot][ θ][ + 2][η][1]
~~θ~~







∂ 1
~~r~~




- 0 1
Oℎ η



∂ 1 −2 [1]
∂τrr
~~r~~ ~~r~~



(2.29)



and the initial conditions of the second order are



1 ∂ 2
2(, 0) = − 2 + 1 [,] ∂η [(][θ][, 0) = 0.] (2.30)
η θ
~~m~~ ~~t~~



In the next step, the solutions of the first-order equations will be presented.


### 3. Solutions of the governing equations
##### (a) First-order solutions



The linear problem is governed by the first-order equations. For viscous, Newtonian drop
liquid, well-known linear solutions are owing to Lamb [2,3] and Chandrasekhar [4]. The
two-dimensional flow field allows the method of the Stokesian stream function to be applied
to determine the velocity and pressure fields. The stream function (,, ) is related to the two
velocity components 1 and 1 by [44] ψ r θ t
r θ
u u



∂ 1
and 1 =
∂ψ θ sin
u
~~θ~~ ~~r~~ ~~θ~~



∂
∂ψ [.]
~~r~~



1 ∂ 1 ∂
1 = − and 1 = (3.1)
r sin ∂ψ θ sin ∂ψ [.]
u u
~~r~~ [2] ~~θ~~ ~~θ~~ ~~r~~ ~~θ~~ ~~r~~



These formulations ensure that the first-order velocity field satisfies the continuity equation.
The first-order drop surface deformation is governed by the Legendre polynomial of the
initial deformation. The solution is therefore sought in the form



1(, ) = 1 (cos ) −αmt, (3.2)
η θ t η [^] Pm θ e



with the first-order initial surface amplitude 1 and the complex angular frequency for the

deformation mode . η [^] αm



with the first-order initial surface amplitude 1 and the complex angular frequency for the

deformation mode . [^] m
To find the first-order solutions of the equations of motion, we first solve the RCE of the firstm
order. Given the exponential dependency of all the flow field variables on time via exp( − ),
we find for the extra-stress tensor of first order αmt



1 = 2 [1 −] [De][2][α][m] 1 =: 2 1 1 . (3.3)
1 − 1
**τ** ~~De~~ ~~αm~~ **D** β **D**



This means that, at first order, the extra-stress tensor of the viscoelastic fluid differs from the
Newtonian material just by a frequency-dependent factor 1 in front of the rate-of-deformation
β


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

tensor. This fluid is therefore formally identical to a Newtonian one, so that all the first-order
results obtained in our previous paper [36] apply here, just with the Ohnesorge number 0 in
those equations replaced by the modified Ohnesorge number =: 1 0. We therefore applyOℎ
those previous results as follows. Oℎv β Oℎ
The curl of the vectorial momentum equation consisting of equations (2.17) and (2.18) yields
the partial differential equation









for the stream function, where the operator is given as [44]




[+]











The solution is the stream function = 1 + 2 [45], where



1 2

1 = 1 sin [2] ′ (cos ) −αmt and 2 = 2 ( ) sin [2] ′ (cos ) −αmt. (3.6)
ψ C mr [m][ + 1] θPm θ e ψ C mqrjm qr θPm θ e



1 1 2 2

In these solutions, ′ (cos ) is the first derivative of the Legendre polynomial with respect to
its argument. The symbol Pm θ denotes a spherical Bessel function of the first kind and order Pm .
jm m



its argument. The symbol denotes a spherical Bessel function of the first kind and order .

The parameter in its argument equals / . The velocity components in the radial and the
polar angular directions are obtained asq αm Oℎv



and



1 = 1 ( + 1) + 2 ( + 1) [j][m][(][qr][)]  - + 1( ) sin ′ (cos ) −αmt, (3.8)
uθ C m m r [m][ −1] C mq [2] m jm qr θPm θ e
~~qr~~





respectively. The two integration constants 1 and 2 are determined by the first-order
kinematic and zero-shear stress boundary conditions and readC m C m







The first-order pressure field is obtained by integration of the momentum equations as



1 = - 1 ( + 1) (cos ) −αmt . (3.10)
p C m m αmr [m] Pm θ e



The solutions for the velocity and pressure fields allow the boundary condition (2.21) for the
normal stress to be formulated, yielding the characteristic equation of the drop



2



2



+ 1

This equation determines the complex angular oscillation frequency, where we have denoted



This equation determines the complex angular oscillation frequency, where we have denoted

the non-dimensional Rayleigh frequency, 0 = [ ( −1)( + 2)] [1/2] . The equation is formally
m
identical to the results of Lamb, Chandrasekhar [α m 2,4m ] and Khismatullin & Nadim [m 9]. The
spherical Bessel functions are taken at the value of their arguments. We can express the square

of the argument, [2], as q



of the argument,, as
q [2]



= m, 0 Ωm0 1 −Ω1 −ΩmDe12αm, 0, 0 [,] (3.12)
q [2] α ~~Oℎ~~ ~~mDe~~ ~~αm~~



Ω
m0
~~Oℎ~~



1 −Ω 1, 0
1 −Ωm ~~m~~ De ~~De~~ 2α ~~αm~~ m, 0 [,]


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

**Table 1.** Measured properties of the aqueous Praestol 2500 solutions and measured drop radius. Density = 10 [3] / .







![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-9-1.png)



|0.3|0.043|0.073|0.078|0.925 x 10−3|0.16|
|---|---|---|---|---|---|
|0.5|0.35|0.074|0.14|1.08 x 10−3|1.238|


where Ω = /, 0. This form is identical to the results by Brenn & Teichtmeister [31] and by
Khismatullin & Nadim [m αm αm 9]. For 1 = 0 and 2 = 0, equation (3.12) reduces to the Newtonian
case. De De
In the following section, the solutions of the characteristic equation will be shown and their
importance will be addressed. A method of finding the appropriate solutions with the use of
experiments will be proposed.

##### (b) Solutions of the characteristic equation and first-order flow fields


For identifying solutions of the characteristic equation representing correctly the time behaviour of oscillating drops, both in terms of frequency and decay rate, shape oscillations of
individual drops from aqueous solutions of the polyacrylamide Praestol 2500, levitated in an
ultrasonic resonator, were investigated experimentally [30]. The measured zero-shear viscosity
0, surface tension and stress-relaxation time 1, characterizing the drop liquid for different
solute mass fractions (SMF), together with the drop radius μ σ λ, are given in table 1. Owing to
the small solute contents, the density for all SMF equals 1000 kg/ma [3] . These data are used
to compute the Deborah numbers 1. The deformation retardation time ρ 2 will be computed
from the Deborah number 2 selected for solving the characteristic equation of the drop. TheDe λ
procedure will be explained in detail below.De
For the given drop radii and liquids with SMF of 0.3 and 0.5 wt%, the values of the Deborah
number 1 are 23.71 and 33.93, respectively. The Ohnesorge numbers for the former and latter
cases are De 0 = 0.16 and 0 = 1.238, respectively. The solutions of the characteristic equation
are non-dimensional complex angular frequencies ΩOℎ Oℎ, where is the deformation mode. The
imaginary and real parts of the complex angular frequency are the frequency and the decaym m
rate of the drop-shape oscillations, respectively. They are shown in figure 2 as functions of 0
for the deformation modes = 2 and = 4. In order to show the behaviour of the drops up toOℎ
aperiodic states, the Ohnesorge number ranges extend to 8 and 30 in m m figure 2 _a,b,c,d_, respectively.
For fixed Deborah and Ohnesorge numbers, the characteristic equation yields an infinite
number of roots [9]. We show only one pair of solutions since the others assume large
real values, which are not relevant since the corresponding motions disappear rapidly. The
remaining problem is to select a pair of complex conjugate values appropriate for representing
the drop-shape oscillations observed in the experiments. For doing this, 2 is varied such that
the decay rate of the oscillations agrees with the value from the experiment. The procedure ofDe
selection will be explained in the next paragraph. Figure 3 shows traces of the drop north pole
obtained from the experiments with the two viscoelastic liquid drops in table 1. For comparison,
in the same figure, we show the first-order solutions of the drop north pole position (0, )
as a function of time. Frequency and decay rate of the linear solutions are obtained from thers t
complex conjugate angular frequency for the corresponding liquid shown in figure 2. The
surface deformation amplitudes and phase angles at are set according to the experimental data.
Length and time in the experimental data are non-dimensionalized with the droplet radius and
capillary timescale, respectively. The theoretical and experimental data agree very well within
the uncertainty of the measured drop radius of, which is shown by the error bars.

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( _a_ )



1.05


0.75


0.45


0.15


−0.15


−0.45


−0.75



( _b_ )
4


3.5


3





![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-10-1.png)



−1.05

|Col1|Col2|
|---|---|
|||
|||
|||
|||
|||
|||
|||
|||
|||
|||
|||
|||
|||


0 2 4 6 8


_Oh_
0



2.5


2


1.5


1


0.5


3.5


3


2.5


2


1.5


1


0.5


|m = 2<br>m = 4<br>Exp. data|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|
|---|---|---|---|---|---|---|---|---|
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data||||||
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data||||||
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||



( _c_ ) ( _d_ )
1.05



0 2 4 6 8

_Oh_
0



0.75


0.45


0.15


−0.15


−0.45


−0.75


−1.05


|Col1|Col2|Col3|
|---|---|---|
||||
||||
||||
||||
||||
||||
||||
||||
||||
||||
||||
||||
||||


|m = 2<br>m = 4<br>Exp. data|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|
|---|---|---|---|---|---|---|---|---|---|
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|||||||
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|||||||
|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|m = 2<br>m = 4<br>Exp. data|||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||



0 5 10 15 20 25 30



![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-10-2.png)



0 5 10 15 20 25 30 0 5 10 15


_Oh_ _Oh_
0 0



_Oh_ _Oh_
0 0

**Figure 2.** Solutions of the characteristic equation for polymer mass factions of ( _a_, _b_ ) 0.3 wt% ( 1 = 23.71) and ( _c_, _d_ ) 0.5
wt% ( 1 = 33.93), as functions of the Ohnesorge number. ( _a_, _c_ ) Imaginary and ( _b_, _d_ ) real parts of the non-dimensionalDe
frequency ΩDe .
m



_Oh_
0



1.06


1.02


1


0.98


0.96

|0.3 wt% − First order|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|Col11|Col12|
|---|---|---|---|---|---|---|---|---|---|---|---|
|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|0.3 wt% − First order<br>− Exp. data<br>0.5 wt% − First order<br>− Exp. data|
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||



0 1 2 3 4 5 6 7 8



t

**Figure 3.** The trace of the drop north pole in time obtained from theory for = 2 and experimental data.
m



t



In figure 2, we show for each mode only the pair of complex conjugate values for the
given range of Ohnesorge numbers. For calculation of the presented solutions, the value of them
Deborah number 2 is selected such that, for the fixed Ohnesorge number 0 and Deborah
number 1, the real part of the solution ΩDe 2 matches the experimental non-dimensional decayOℎ
De


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

rate for the deformation mode = 2. This applies the oscillating drop method of Brenn &
Plohl [30] in a simplified manner, laying the focus on the decay rate. The selected Deborahm
numbers are 2 = 3.5 and 2 = 1.26 for the two drops in table 1 with SMF of 0.3 and 0.5 wt%,
respectively. The former and latter numbers, with the measured properties in De De table 1, yield
the deformation retardation times 2 of 0.0115 and 0.0052 s, respectively. For the same drops,
the measured oscillation frequencies are λ = 133.28 Hz and = 113.12 Hz, respectively, and the
decay rates are 2, = 33.13 s [−1] and 2, = 51.02 sf [−1] . Therefore, the complex angular frequenciesf

 - = 2π are known. The measured non-dimensional complex angular frequency Ωα r α r one





decay rates are 2, = 33.13 s [−1] and 2, = 51.02 s [−1] . Therefore, the complex angular frequencies
r r

2* = 2, + 2π are known. The measured non-dimensional complex angular frequency Ω2 one
αformulates as α r i f */ -, where - is the dimensional Rayleigh frequency. The results with the



2* = 2, + 2π are known. The measured non-dimensional complex angular frequency Ω2 one
r

formulates as 2*/ 2,0*, where 2,0* is the dimensional Rayleigh frequency. The results with the
corresponding Ohnesorge numbers are shown by the black circles in α α α figure 2. The error bars
—barely visible in the diagrams—represent uncertainties in the frequency and decay rate of
±2 and ±11%, respectively, following the work by Brenn & Plohl [30]. The comparison yields
good agreement between the measured and the calculated (solid black line) non-dimensional
complex angular frequencies for the mode of initial deformation = 2. For drops with different
0, but the same = 2, as well as 1 and 2, therefore, the values of Ωm 2 can be taken from
those lines. For the second-order approximation, the complex conjugate solution for Oℎ m De De = 4 in
figure 2 will be needed. This will be detailed later. m
Since there exist pairs of solutions of the characteristic equation, the first-order solutions are



Since there exist pairs of solutions of the characteristic equation, the first-order solutions are

rewritten, denoting the frequency with the positive imaginary part as ( ), and the one with the
p

( ) αm



rewritten, denoting the frequency with the positive imaginary part as ( ), and the one with the

negative imaginary part as ( ). The resulting first-order deformation of the drop surface reads
n
αm



η1(θ, t) = η [^] (1p)e−αm(p)t + η^(1n)e−αm(n)t Pm(cos θ) . (3.13)



This formulation allows the two first-order initial conditions to be used for determining the



amplitudes η [^] (1p) and η^(1n). We obtain



( )

( ) m ( ) [.]

  ~~p~~ ~~n~~
~~αm~~ ~~αm~~



( ) ( )
( ) ( )
1p =  - ( α) −mn ( ) and 1n = ( α) −mp ( ) [.] (3.14)
η [^] ~~p~~ ~~n~~ η [^] ~~p~~ ~~n~~

~~αm~~ ~~αm~~ ~~αm~~ ~~αm~~



( ) ( )

( )
( α) −mn ( ) and 1n = ( α) −mp
~~p~~ ~~n~~ η [^] ~~p~~
~~αm~~ ~~αm~~ ~~αm~~



The first-order stream function, velocity components and pressure are given in the electronic
supplementary materials. With these findings, the first-order solutions, which are formally
identical to the Newtonian ones obtained by Zrnić _et al_ . [36], are complete. One important
difference is the modified Ohnesorge number entering the solutions, which depends on the
frequency and on the two viscoelastic timescales.


##### (c) Second-order solutions



The second-order solutions are developed starting from pressure. We search for the general
solution for the second-order pressure as a sum of two contributions



2(,, ) = 21(,, ) + 22(,, ), (3.15)
p r θ t p r θ t p r θ t



where the subscript ‘21’ denotes the solution of the inhomogeneous system of second-order
equations, and the subscript ‘22’ is the solution of the homogeneous system [36].
The solution with subscript ‘21’ is first determined. We analyse the second-order equations
of motion (2.23)–(2.25) and form the divergence of the vectorial momentum equation, which,
with the use of the continuity equation (2.23), eliminates the second-order velocities from the
momentum equation. The differential equation for the second-order pressure 21 is obtained as
p



21 is obtained as

Δ 21 = - (( 1 1 [) +][ Oℎ] 21 [],] (3.16)
p ∇⋅ [→] u [⋅∇][)][u→] - [∇⋅] [[][∇⋅] **[τ]**



which we rewrite, using the Lamé identity, into the form


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

2
Δ 21 = - 1 /2 − 1 [×] 1 [−] [Oℎ] 0 [∇⋅] **[τ]** 21 [.] (3.17)
p ∇⋅∇u [→] u→ [∇] [×] [u→]





The present flow field allows equation (3.17) to be rewritten as



p [→] u ▽⋅







2
This is a Poisson equation for the modified pressure P21 = 21 + 1/2 [36,46]. We determine
p [→] v



the structure of the dependencies on the polar angular coordinate and time by inspecting

2
the square of the first-order velocity vector, 1, which appears in the modified pressure. The
inspected structure of the dependencies follows also from the products of first-order terms on [→] v
the right-hand side of equation (3.18). Both groups of terms contain functions exp( −2 ) for
only positive ( _p_ ) or negative ( _n_ ) time dependency. The RCE (2.26) corresponding to the timeαmt



the structure of the dependencies on the polar angular coordinate and time by inspecting



only positive ( _p_ ) or negative ( _n_ ) time dependency. The RCE (2.26) corresponding to the time

dependencies according to exp(−2 ( ) ) or exp(−2 ( ) ), with the first-order extra stress tensor 1
p n
known, reads αm t αm t τ



21 = 2 2 21 + 2 1 1 1 − 1 1 [−] **[D]** 1 1T . (3.19)
→
**τ** β **D** γ [→] v [⋅∇] **D** →v→∇ [⋅] **[D]** [⋅∇] →v→



The coefficients in the equation are



2 = [1 −2][α][m][De][2] and 1 = (3.20)
1 −2 1 1 −2 1 [.]
β ~~αmDe~~ γ [De][2][ −] ~~α~~ [β] ~~m~~ [1] ~~De~~ [De][1]



For the positive ( _p_ ) or negative ( _n_ ) time dependency, we substitute the second-order extra stress
tensor from equation (3.19) into the right-hand side of equation (3.18) and obtain








[→] v [⋅∇] **D** →v→∇ [⋅] **[D]** [⋅∇] →v→



since 21 = 0. We substitute the stream function and its constituent 2 into the







2
+





2

 


+





2 +



(3.22)



+









~~m~~

→



where = cos . In equation (3.22), the coefficients 1 and 2 and the quantity represent
x( ) θ( ) C m C m ( ) q ( )



where = cos . In equation (3.22), the coefficients 1 and 2 and the quantity represent

either ( ) or ( ), which correspond to the time dependencies with exp(−2 ( ) ) or exp(−2 ( ) ),
p n p n
αm αm αm t ( ) αm t



either ( ) or ( ), which correspond to the time dependencies with exp(−2 ( ) ) or exp(−2 ( ) ),

respectively. One further time dependency comes from products of terms with exp(− ( ) ) and
p

( ) ( ) + ( ) αm t



respectively. One further time dependency comes from products of terms with exp(− ( ) ) and

exp(− ( ) ), which infer the time dependency with exp[−( ( ) + ( )) ]. The RCE for this time
n p n
dependency is obtained asαm t αm αm t







+



,



(3.23)




```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

and



( ) ( )





Substituting the second-order extra stress tensor from (3.23)

(3.18), the equation for the time dependency with exp[−( ( ) + ( )) ] reads
p n
αm αm t



Substituting the second-order extra stress tensor from (3.23) into the right-hand side of equation





∂ψ∂2(n)



( )
2n

[+]
∂



![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-13-4.png)





1 ~~n~~ 0





1





.



We substitute the stream function and its constituent 2 into the right-hand side of equation



( ( ))
n
d q [(][n][)] rjm
~~dr~~



ΔP(21pn) = - C2(pm)C2(nm) q [(][p][)2] + q [(][n][)2]

~~r~~ [2]



( ( ))
p
d q [(][p][)] rjm
~~dr~~



( (p))( (n) (n)) + 1(n) 2(p) + 1 (p) [2] m −2 [d][(][q][(][p][)][rj] m(p)) +



−2 2(p) 2(n)
C mC m q [(][p][)2] q [(][n][)2]

~~r~~ [2]



~~r~~ + [2]



C mC m m q r





+ 2p 2n



2
(1 − ) +
x [2]







+







.



(3.27)



The further development of the right-hand sides of the Poisson equations for the modified
pressure for all time dependencies is detailed in the electronic supplementary material owing to
lengthy expressions.
The solution P21, of the Laplace equation for the modified pressure is first determined,
H



P21,

which, for the time dependency according to exp (−2 ( ) ), reads
p
αm t


###### ∑=L 0 D21(p)lrlPl(cos θ)e−2αm(p)t .


###### P21, H = − ∑= 0 21( )l l l(cos ) −2αm( )t . (3.28)

l



In this equation, the 21(p) are expansion coefficients. The sums of these solutions, arising
from the application of the Legendre polynomial orthogonality, are taken up to a maximumD l
summation index (upper-case letter). Particular solutions P21, (,, ) for the right-hand sides
I r θ t

in equations (3.22) and (3.27) in the forms of known functions seem to be out of reach.
We therefore approximate the right-hand sides, replacing the spherical Bessel functions and
products of Legendre functions by their series expansions for a given value of the mode of
initial drop deformation. Their forms are presented in the electronic supplementary material.m


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

The particular solutions are searched in the forms of power series of the radial coordinate and
of series of Legendre polynomials as functions of the polar angle. The approximation accuracy
increases with the order of approximation by the series expansions.
The solution for the modified pressure is composed of the sum of the general solution of the
Laplace equation and a particular solution of the Poisson equation as





21(,, ) = 21, (,, ) + 21, 1(,, ) + 21, 2(,, ), (3.29)
P P P P
r θ t H r θ t I r θ t I r θ t



1 2

where the particular solution 21, 1(,, ) satisfies the Poisson equation (3.22) for the modified
P
I r θ t



pressure with only the right-hand side term in curled brackets, while the particular solution
21, 2(,, ) satisfies the Poisson equation with only the right-hand side term with the diverP
I r θ t



gence. Analogous solutions are sought for the Poisson equation (3.27).
The second-order pressure contribution 21 accounts for all three different time dependencies and reads p



p21(r, θ, t) = p21(p)(r, θ, t) + p21(n)(r, θ, t) + p21(pn)(r, θ, t). (3.30)



21 21 21 21

We determine the contribution to the pressure solution for the time dependency ( ) with the

2 p



2
definition of the modified pressure P21 = 21 + 1/2 and obtain
p [→] v



2k + 2m + 2 ( ) +
Pl x


###### 21(p),, = − −2αm(p)t ∑=L 0

p r θ t e

2 −1 l



= 0 = 0


( )
N p(2 )

= 1 δk m q [(][p][)] r



2

N m
###### ∑ ∑

= 0 = 0

l



2
m ( ) ( )
###### ∑ 2p p

= 0 G kl q r



2 + 2
###### k m ( ) + ∑N
Pl x = 1

[2] k



2

N m
###### ∑ ∑

= 0 = 0

l



( )
###### L 21p l ( ) − ∑N

= 0 D lr Pl x = 0

k



2k + 2m (2 )( ) +
m
P x


###### + ∑

= 0
k



2 −1
m ( ) ( )
###### ∑ p p

= 0 δkl q r
l



(2 ) (2 )

= 0 = 0 = 1



(3.31)





to large expressions, the contribution 21, 2 for all time dependencies is also detailed in
P
I

the electronic supplementary material. The coefficients in these pressure solutions for each
time dependency are determined as presented in the electronic supplementary material. The
second-order solution 21 consists of two parts: the first two terms—the one in curled brackets

( ) p



second-order solution 21 consists of two parts: the first two terms—the one in curled brackets

and the term P21,(p) 2(,, )—are obtained by solving the Poisson equation for the modified
I r θ t

pressure. The second part represents the dynamic pressure caused by the first-order velocity
field.
The corresponding second-order velocities are determined from the pressure field, the RCE
equation (3.19) and the equations of motion (2.23)–(2.25). We obtain the radial momentum
equation, with the modified pressure P21 introduced for the time dependency ( ), as
p





( ) 2

∂ [2] ( )

∂ ~~r~~ [2] [(][r][2][u][r] p [21] ) + ~~Oℎ~~ [2][α] m0p ~~β~~ r2( ~~p~~



∂P21,(p) 1
I



∂
~~r~~



∂ [2]



( ) 2









∂




[)] +





1

sin [2]
~~r~~ [2] ~~θ~~



1
( )
~~β~~ 1 ~~p~~ ~~Oℎ~~ 0







=



∂



p





(
2

∂
~~t~~



1 0



(3.32)



where the tensor









.


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

The equation for the time dependency ( ) is analogous. For the time dependency according to

exp [−( ( ) + ( )) ], the radial momentum equation readsn



The equation for the time dependency (

exp [−( ( ) + ( )) ], the radial momentum equation reads
p n
αm αm t









∂ [2]





21,

∂
~~r~~

















=











I







+



∂P21,(pn)1
I



( )
I     








where







(3.34)


(3.35)





+


.



The further development of the second-order radial momentum equation for all the time
dependencies is detailed in the electronic supplementary material. The equation is solved for all
the time dependencies of the second-order radial velocity 21 to obtain
r
u

21(,, ) = (p21)(,, ) + (n21)(,, ) + (pn21 )(,, ). (3.36)
ur r θ t ur r θ t ur r θ t ur r θ t

The solutions 21 for all the time dependencies are given as the sum of the general solutions
r
of the homogeneous differential equations and one particular solution of the inhomogeneousu
equations. As an example, we present the solution for the positive time dependency ( _p_ ), which
reads



21(,, ) = (p21)(,, ) + (n21)(,, ) + (pn21 )(,, ). (3.36)
ur r θ t ur r θ t ur r θ t ur r θ t



( )
( + 1) ( ) p t +
l l Pl x e [−2][α][m]



( )
###### p21,, = − ∑L
ur r θ t = 0
l





L



Pl x





= 0


( ) ( )
1p p
K kl q r



2

N m
###### ∑ ∑

= 0 = 0

2l −1


###### + ∑

= 0
k 2



2
m



2 + 2 + 1 ( ) ( )



2k + 2m + 1 ( ) +
Pl x



= 0 = 0

2 −1

N m
###### ∑

= 0 = 0

l


###### + ∑

= 0
k
###### + ∑N

= 1
k



= 0

2 −1
m ( ) ( )
###### ∑ 3p p

= 0 K kl q r
l



L



= 0 = 0


( )
N
3p(2 )

= 1 K k m q [(][p][)] r



2k + 2m + 1P(2m)(x) e−2αm(p)t + Bpr(p21)(r, θ, t)



(3.37)


###### =: − ∑

= 0
l



( ) ( ) ( ) ( ) −2 ( ) ( ) ( )
A21plfarp + B21plfbrp Pl(cos θ)e αmp t + Aprp21(r, θ, t) + Bprp21(r, θ, t),



( ) ( ) ( )
where we have introduced the new argument 2p = p / 2p 0 of the spherical Bessel

functions. The radial velocity components, together with the solutions q αm β Oℎ (,, ) and



where we have introduced the new argument 2( ) ( )/ 2( ) 0 of the spherical Bessel
m

functions. The radial velocity components, together with the solutions 21(,, ) and

(,, ) for every time dependency, are given and explained in the electronic supplementaryApr r θ t



functions. The radial velocity components, together with the solutions 21(,, ) and
r

21(,, ) for every time dependency, are given and explained in the electronic supplementary
material.Bpr r θ t
The known radial velocity component 21 allows the second-order velocity component in
r
the direction of the polar angle to be determined, using the continuity u equation (2.23). The
solution is represented as



21(,, ) = (p21) (,, ) + (n21) (,, ) + (pn21)(,, ) . (3.38)
uθ r θ t uθ r θ t uθ r θ t uθ r θ t



The solution for the positive time dependency ( _p_ ) reads


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( )
###### p21,, = − ∑L
uθ r θ t = 0
l



L



21(p)( + 1) l −1 +



= 0







2 + 2 + 1 + 1( ) − −1( )







( ) ( )
1p p



2

N m
###### ∑

= 0 = 0

2l




###### + ∑

= 0
k
###### + ∑N

= 0
k
###### + ∑N

= 1
k



2
m



= 0

2 −1
m ( ) ( )
###### ∑ 3p p

= 0 L kl q r
l



= 0 = 0

2

N m
###### ∑

= 0 = 0

l



2


2 + 2 + 1 ( ) ( )
k m + 2p p
L kl q r



= 0 = 0


( )
N
3p(2 )

= 1 L k m q [(][p][)] r





L



(3.39)



( )
###### + p21(,, ) =: − ∑L
Bpθ r θ t = 0

( ) (l )



( ) ( ) ( ) ( ) 1 −2 ( )





= 0

+ (p21) (,, ) + (p21) (,, ) .
Apθ r θ t Bpθ r θ t



21 21

The functions,, and, defined for each oscillation frequency, are introduced to
facilitate the expressions of the flow variables. The solutions for the polar angular velocityfar faθ fbr fbθ
component for the two other time dependencies are given in the electronic supplementary



component for the two other time dependencies are given in the electronic supplementary

material. The solutions (p21) (,, ) and (p21) (,, ), together with (p21)(,, ) and (p21)(,, ),

respectively, satisfy the second-order continuity equation.Apθ r θ t Bpθ r θ t Apr r θ t Bpr r θ t



respectively, satisfy the second-order continuity equation.
The structure of the above-given solutions determines the second-order deformation of the
drop surface 21 as
η



L





l



The deformed shape is governed by all the Legendre polynomials in the sum. The reason is
that, in the zero normal-stress boundary condition (2.29), pressure and surface deformation
are coupled. The set of unknown coefficients 21, 21 and 21 are determined for each time
dependency and summation index from the second-order boundary conditions A l B l H l (2.27)–(2.29).
The electronic supplementary material details the calculation of these coefficients.l
The second contributions to the second-order solutions, with subscript ‘22’, are determined
from the homogeneous forms of equations (2.23), (2.24) and (2.25). This set of equations has
the same structure as at first order, where the complex angular frequency is replaced by 2
(with the number ‘2’ in the subscript indicating the second-order solution, and αm a deformationα k

mode number), and the deformation amplitude of the drop surface [^] by [^] k. Both replaced



(with the number ‘2’ in the subscript indicating the second-order solution, and a deformation

mode number), and the deformation amplitude of the drop surface 1 by 22 . Both replaced

variables are under a sum with the summation index, implying that all possible deformationη [^] η [^] k



mode number), and the deformation amplitude of the drop surface 1 by 22 . Both replaced

variables are under a sum with the summation index, implying that all possible deformation [^] [^] k
modes affect the results. k
For the exponential dependency on time via exp (−k 2 ), positive ( _p_ ) or negative ( _n_ ), we find
the second-order solutions of the extra-stress tensor, with subscript ‘22’α kt



1 − 2 2

= 0K 1 − De ~~De~~ 1α ~~α~~ 2k ~~k~~ **D** 22k =: 2 ∑k = 0K



K 2 22 .

= 0 β k **D** k



1 − 2 2

**τ** 22 = 2 ∑k = 0K 1 − De ~~De~~ 1α ~~α~~ 2k ~~k~~ **D** 22k =: 2 ∑k = 0K β2k **D** 22k . (3.41)



The second-order pressure, velocity and drop surface deformation ‘22’ are linear combinations
of eigensolutions of the differential equations with the summation index . They exhibit forms
close to the solutions of the first order. As an example, the velocity component k 22 in the radial
r
direction reads u


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

22(,, ) =  - ∑
ur r θ t = 0
k









×



× (cos ) − ∑
Pk θ = 0
k







(cos ) =



K



(3.42)



=: - ∑
= 0
k





= 0

with the argument 2 = 2 / 2 0 and the definition of the function = + 1( 2 )/ ( 2 ). The

function for each oscillation frequency is introduced in the interest of compact formulations.q k α k β kOℎ Qk jm q k jm q k



with the argument 2 = 2 / 2 0 and the definition of the function = + 1( 2 )/ ( 2 ). The
k k k k m k m k

function 2 for each oscillation frequency is introduced in the interest of compact formulations.

The solutions fr and are given in the electronic supplementary material. The drop surface



function 2 for each oscillation frequency is introduced in the interest of compact formulations.
r

The solutions 22 and 22 are given in the electronic supplementary material. The drop surface
θ
deformation ‘22’ readsu p



K



η [^] (22p)k e−α2(pk)t + η^(22n)k e−α2(nk)t Pk(cos θ).


###### η22(θ, t) = ∑ = 0K η [^] (22p)k e−α2(pk)t + η^(22n)k e−α2(nk)t Pk(cos θ). (3.43)

k



The characteristic equation obtained from the boundary condition (2.29) determines the
complex conjugate angular frequencies 2 . This equation is formally identical to the first-order
equation (3.11), but it is solved by replacing the mode α k of initial drop deformation by the
m [^] ( ) ^( )



equation (3.11), but it is solved by replacing the mode of initial drop deformation by the

summation index, as required by the solution. The amplitudes (22p) and ^(22n) are calculated
using the second-order initial conditions k (2.30) [46]. The sum of the two contributions toη [^] k η k
the second-order drop deformation must satisfy the second-order initial conditions for every
summation index . These contributions therefore depend on each other. Applying the initial
k [^] ( ) ^( )



summation index . These contributions therefore depend on each other. Applying the initial

conditions, one obtains a set of two equations with the two unknown coefficients (22p) and ^(22n) .
The method applied for determining these coefficients is the same as on page 10 of our paperη [^] k η k
Zrnić _et al_ . [36]. We therefore do not repeat its description here.
This analysis shows that the cases of a volumetric drop pulsation ( = 0) and a translation
motion of the drop ( = 1) drop out from the solutions ‘22’. This is plausible, since the case k = 0
is impossible for incompressible liquid, and the case k = 1 would be possible only owing to ak
resultant force acting on the drop, which does not exist here.k



resultant force acting on the drop, which does not exist here.

The two amplitudes (22p) and ^(22n) for values ≥2 are determined using the second-order
initial condition (2.30). The coefficient η [^] k η k 21 is known from the boundary conditions, and thek

( ) H k

[^]



( )
equation for the amplitude 22p is
η [^] k





( )

liquid drop-shape oscillations.


### 4. Results and discussion

In this section, we evaluate the solutions derived above. First, we verify the conservation of
the drop volume for the two drops of table 1, varying the deformation parameter 0. Then
we present traces of the drop aspect ratio and deformed drop shapes for different viscoelasticη
liquid material properties and compare the results on the nonlinear behaviour with data for
a corresponding viscous case with the mode of initial deformation = 2. We then show two
well-known nonlinear effects in drop-shape oscillations, which are the excess time in the dropm
prolate form and the frequency change with varying initial deformation. For selected traces of

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

the drop north pole motion, we analyse the quasi-periodicity of motion and the Fourier power
spectra. We conclude by a study of the onset of aperiodic behaviour.

##### (a) Volume conservation


First, we verify the conservation of the drop volume, because the weakly nonlinear theory
represents all the field variables, together with the drop surface shape, as truncated power
series of a small deformation parameter 0. Consequently, the non-dimensional drop volume
may deviate from its exact value of = 4π/3. The deviation is determined analytically as aη
function of time and reads Vs







1



3(, ) −1 cos =: ( 0, ),
rs θ t d θ R η t



where the non-dimensional drop volume ( ) is determined as







1



1 + 1(, ) 0 + 2(, ) 02 [3] cos .
η θ t η η θ t η d θ



The values of are all positive, oscillate in time and increase with 0. Figure 4 presents
deformation = 2 for the two viscoelastic drops in t tmax table 1. In both cases, for deformation
parameters m0 ≤0.35, the volume deviates by less than 2% from the exact value. This is
considered acceptable for the present analysis.η


##### (b) Deformed drop shapes

Figure 5 presents the drop aspect ratio (0, )/ (π/2, ) as a function of time for the mode of
initial deformation = 2 for three different drops: the two viscoelastic liquid drops from rs t rs t table
1, with 0 = 0.16 and m 0 = 1.238, and a viscous drop with = 0.048. The value of the viscous
drop Ohnesorge number is selected such that the linear value of its decay rate is equal to thatOℎ Oℎ Oℎ
value for the viscoelastic drop with 0 = 1.238. For this comparison, all the results presented
account for the solutions up to the second order of approximation. For the viscoelastic dropOℎ
with 0 = 1.238, the Ohnesorge number is well beyond the critical value of 0.7665 of the
Newtonian case. In that case, the result of the characteristic equation suggests that such dropsOℎ
have an aperiodic response to the initial deformation. This is not the case, however, for the far
higher Ohnesorge number 0 = 1.238 in the viscoelastic case, which figure 5 confirms, showing
shape oscillations of the drop. The driver of this viscoelastic drop oscillation is the elasticity, inOℎ
contrast to the viscous case, where it is surface tension. Another interesting comparison below
the critical value is between the viscous drop with = 0.048 and the viscoelastic drop with
0 = 1.238. In these two cases, the linear values of the decay rate are the same, but initiallyOℎ
figure 5Oℎ shows significantly different peaks of the drop aspect ratio between the two cases.
The difference can be explained by the involved elasticity, which plays a significant role in
large deformations, raises the maximum amplitudes and changes the resultant decay rate and
frequency of the oscillation in comparison with the viscous case. The difference in the minimum
aspect ratio in the first oblate state between these two drops is around 12%. For the two
compared drops, the minimum and maximum peaks of the aspect ratio are more pronounced
for the viscoelastic drop than for the viscous drop throughout the oscillations. This behaviour
is explained by the elasticity as an extra driving force. As a further consequence, these peaks
always appear earlier in time than for the viscous drop, i.e. the oscillation frequency of the
viscoelastic drop is higher.


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

0.025


0.02


0.015


0.01


0.005






|Col1|Oh = 0.16<br>0<br>Oh = 1.238<br>0|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|
|---|---|---|---|---|---|---|---|---|---|
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||
|||||||||||



0 0.1 0.2 0.3 0.4

_η_                    

**Figure 4.** Maximum relative volume deviation from the exact value for the viscoelastic liquid drops from table 1 with
0 = 0.16 and 0 = 1.238 at = 2 as a function of the deformation parameter 0.
Oℎ Oℎ m η


1.7



1.5


1.3


1.1


0.9


0.7


0.5



|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|Viscous − Oh = 0.048<br>Viscoelastic − Oh = 0.16<br>0<br>− Oh = 1.238<br>0|Col12|Col13|Col14|Col15|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||
||||||||||||||||


0 3 6 9 12



t

**Figure 5.** The drop aspect ratio over time for three liquid drops with = 2, 0 = 0.35, from the theory up to the second
order. m η


The present results from the analysis allow the deformed drop shapes up to second order to
be plotted, representing nonlinear effects. Figure 6 _a,b_ show the meridional sections of deformed
drops for the mode of initial drop deformation = 2, i.e. for prolate-to-oblate oscillations, with
the relatively large initial aspect ratio of the drop of m / = 1.65 of the three liquid drops in
figure 5. The degrees of Legendre polynomials contributing to the solutions of the second orderL W
are 0, 2 and 4. The drop shapes are represented for the first oscillation period at the states
of maximum deformation, omitting the initial state. The time instants of the maximum states
shown in figure 6 _a,b_ are related to the trace of the north pole data in figure 5 for 0 = 1.238. The
corresponding case assumes the maximum deformed states at the earliest times in comparisonOℎ
with the drops with = 0.048 and 0 = 0.16. Owing to the dampening, the first oblate and
second prolate shapes are, of course, less deformed than the initial one. With ongoing time, asOℎ Oℎ
oscillation amplitudes decrease, the shapes converge to the linear ones.


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( _a_ ) ( _b_ )


0 0







![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-20-1.png)







_π_ /2









**Figure 6.** The ( _a_ ) first oblate and ( _b_ ) second prolate shapes for the initial deformation mode = 2 with 0 = 0.35 of three

the first minimum and second maximum drop aspect ratio, seen in figure 5, for the viscoelastic drop with 0 = 1.238. The
results include solutions up to the second order. Oℎ


##### (c) Timescales of the drop-shape oscillations

Two well-known nonlinear effects in drop-shape oscillations with = 2 are the asymmetry of
the times the drop spends in the oblate and prolate deformed states and a frequency changem
with varying drop deformation. One expects the former effect to be less pronounced in the
present viscoelastic than in the inviscid case, but similar to the viscous case. The second-order
solutions can represent these phenomena. In figure 7 _a_, the second-order solutions for = 2
show the time the drop spends in the elongated form, determined from the first oscillation,m
except for the inviscid analysis that is based on the solutions up to third order and data are
based on the first 10 oscillations [35]. The shaded area presents the s.d. from the mean value
denoted with the solid black line. The relatively high Ohnesorge number of 0.1 yields only
a few oscillations, whereas the second and following periods exhibit very small amplitudes
and thus do not show nonlinearities. The relatively small Ohnesorge number of 0.048 shows
more pronounced excess time and in comparison with other viscous and viscoelastic liquid
drop cases presents a case with the strongest nonlinear effects. In the comparison of the two
viscoelastic cases with 0 = 0.16 and 0 = 1.238, we see a difference, but not as significant as
in the comparison between the two viscous cases, although the Ohnesorge number differenceOℎ Oℎ
is smaller. Also, both viscoelastic cases yield significantly smaller excess times for the present
range of initial drop aspect ratios. Comparing the viscoelastic drop with 0 = 1.238 and the
viscous drop with = 0.048, which have the same linear decay rate, the latter drop has moreOℎ
than twice as large excess time than the former for the whole range of initial aspect ratiosOℎ
studied. This means that the elasticity reduces the excess time in favour of an increase in the
maximum amplitudes of the drop north pole motion.

Figure 7 _b_ shows the frequency change of the drop north pole motion for the initial deformation mode = 2 as a function of the initial drop aspect ratio. The results are presented in the
same way as in m figure 7 _a_ . The frequency change is normalized by the imaginary part of the
first-order oscillation frequency 2, for the viscous and viscoelastic cases. In the inviscid case
also shown, the frequency change is normalized by the Rayleigh frequency. The results by Zrnićα i
& Brenn [35] for the inviscid case are shown to quantify the difference between the inviscid on
the one hand and the viscous and viscoelastic drops on the other. We do not see a significant
change in frequency for the second-order solutions, which suggests that an analysis up to the
third order of approximation could add to these results.
The next subject of our study is the oscillatory motion of the drop north pole for the mode
of initial deformation = 2. Figure 8 _a_ depicts the drop north pole position for the deformation
m


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( _a_ ) ( _b_ )





0


−0.025


−0.05


−0.075


−0.1



62


60


58


56


54


52



Third ord. − _Oh_ = 0



![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-21-2.png)




|0.048|Col2|Col3|Col4|
|---|---|---|---|
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||
|0.048<br>  0.16<br>  1.238<br>  0.048<br> = 0.1<br> = 0.1||||



![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-21-1.png)

−0.125
1 1.1 1.2 1.3 1.4 1.5 1.6 1.7



50



1.2 1.3 1.4 1.5 1.6 1.7 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7

Initial drop aspect ratio ( _L_ / _W_ ) Initial drop aspect ratio ( _L_ / _W_ )



**Figure 7.** ( _a_ ) Percentage of the oscillation period spent in prolate form and ( _b_ ) change in frequency as functions of the initial
drop aspect ratio / for = 2. The maximum deformation parameter 0 = 0.4 corresponds to the maximum value on
the _x_ -axis / = 1.78.L W m η
L W


parameter 0 = 0.35 as a function of time for the same three liquid drops: the two viscoelastic
liquid drops from η table 1, and the viscous drop with = 0.048. The first four periods are
denoted by, 1–, 4 in figure 8 _a_ . These oscillation periods are evaluated for the whole of the twoOℎ
p p
oscillations in t figure 8t _a_ . The viscoelastic results from this analysis show that, with ongoing time,
the period length changes, leading to a slight decrease in the frequency in time. Accounting for
the third-order solution in the viscous case leads to an increase in the frequency in time [36].
These behaviours indicate a time dependency of the frequency for nonlinear drop oscillations,
as stated by Prosperetti [8] for the Newtonian case. The data in figure 8 _b_ show a non-monotonic
convergence of the oscillation periods towards the linear values. For the Ohnesorge numbers
= 0.16 and = 1.238, the linear theory predicts the non-dimensional oscillation periods of
2.222 and 2.202, respectively. These values are depicted by the red dashed and blue dottedOℎ Oℎ
horizontal lines in the figure. At the start of the motion, for 0 = 0.16 and 0 = 1.238 with
the larger deformation parameter 0 = 0.35, the period is shorter than the linear value by −0.59Oℎ Oℎ
and −0.68%, respectively. Another interesting finding is that the beat between different modesη
of oscillation may make the period length suddenly increase, even during cycles with small
deformations. The quasi-periodicity of the drop-shape oscillations is caused by a beat between
various oscillation modes excited in the course of the motion [35]. This beat phenomenon,
however, is dampened as the viscosity increases, as shown in figure 8 _b_, where the data for the
higher Ohnesorge number fluctuate much less than those for the lower one.
Next, we determine the damping rate from the results for the oscillations depicted in
figure 8 _a_ . The damping rate is given by pairs of data points representing states of maximumαr
positive or negative deformation against the spherical shape. These points define an oscillation
period. Examples of such oscillations are marked by the starting and ending times of the
periods, 1 and, 3 in figure 8 _a_ . For every such pair of drop north pole positions at an earlier
p p
time 1 and a later time t t 2, the decay rate is deduced as = ln [ (0, 2)/ (0, 1)]/( 1 − 2). The decay
rate is plotted in t figure 8t _c_ against the starting time of an oscillation for the three differentαr η t η t t t
liquid drops. The viscoelastic drop with 0 = 0.16 converges to the value of 0.109, which is
the real part of the complex angular frequency Oℎ 2 predicted by the linear theory, as depicted
by the red dashed horizontal line in figure 8 _c_ . The damping rate for the viscous and viscoelas-α
tic drops with = 0.048 and 0 = 1.238, respectively, converges in time to the blue dotted
horizontal line at the value of 0.21, which is the real part of the complex angular frequency Oℎ Oℎ 2
α


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( _a_ )


( _b_ )


( _c_ )


( _d_ )



1.4


1.2


1


0.8


0.6


2.43


2.38


2.33


2.28


2.23


2.18


0.25


0.2


0.15


0.1


0.05


0.1


0.075


0.05


0.025










|t|Col2|Col3|Col4|Col5|Col6|Col7|Second ord. viscous − Oh = 0.048|Col9|Col10|Col11|Col12|Col13|Col14|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|~~_t_~~|~~_t_~~|_p,_1<br>|_p,_1<br>||||Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|Third ord.<br>Second ord. viscous −_Oh_0 = 0.048<br>Second ord. viscoelastic −_Oh_0 = 0.16<br> −_Oh_0 = 1.238|
|||~~_t_~~_p,_2|~~_t_~~_p,_2|||||||||||
|||||||||||||||
|||||||||||||||
|||||||||||||||
|||||||||||||||
|||||~~_t_~~||||||||||
|||||_p,_4<br>||||||||||
|||||~~_t_~~_p,_3||||||||||


|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Second ord. viscous − Oh = 0.048<br>0<br>Third ord.<br>Second ord. viscoelastic − Oh = 0.16<br>0<br>− Oh = 1.238<br>0|Col10|Col11|Col12|Col13|Col14|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|||||||||||||||
|||||||||||||||
|||||||||||||||
|||||||||||||||
|||||||||||||||



![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-22-1.png)

![](/Users/eaguerov/.cache/pymupdf-mcp/images/rspa_2023_0887/rspa.2023.0887.pdf-22-2.png)

0 3 6 9 12



_t_

**Figure 8.** ( _a_ ) Position of the drop north pole for three different liquid drops with the deformation parameter 0 = 0.35.
Double arrows mark periods starting from different time instants. ( _b_ ) Periods measured from different starting timesη
during the oscillation are shown in (tp _a_ ). ( _c_ ) Damping rate is determined for different pairs of displacement maxima or minima.tp
( _d_ ) Loss angle _ϕ_ ( ). Mode of initial deformation = 2.
t m


for 0 = 1.238 predicted by the linear theory. Comparing the phenomenon between the two
viscoelastic liquid drops, we conclude that, in nonlinear drop oscillations, viscosity affects theOℎ
damping rate differently. The larger viscosity does not necessarily dampen the system more
rapidly, so that, at 0 = 0.35, the drop with 0 = 1.238 approaches the linear value similarly fast
as the drop with η 0 = 0.16. For the former and latter viscoelastic liquid drops, the dampingOℎ
rate from the first period Oℎ, 1 deviates around −22.4 and −29.3%, respectively, from the linear
p
value. The viscoelastic liquid drop with t 0 = 1.238 is dampened faster than the viscous drop
with the same decay rate at first order. The values of the decay rate extracted from the firstOℎ
oscillation period, 1, using up to the second-order solution, are similar, assuming 0.163 and
p
t



_t_


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

0.17 for a viscoelastic drop with 0 = 1.238 and the viscous drop with 0 = 0.048, respectively.
In all cases, the return of a nonlinear damped to a linear motion with ongoing time is seen,Oℎ Oℎ
since the deformations decrease in time. The trend of the decay rate to increase in time agrees
with the effect found for Newtonian drop-shape oscillations by Zrnić _et al_ . [36], which is seen
in the Newtonian data in figure 8 _c_ again. With ongoing time, the data converge to the linear
respective value for the given drop liquid, where the values for = 0.048 and for 0 = 1.238
are the same, since the Newtonian drop with = 0.048 was so designed. It is interesting to seeOℎ Oℎ
that the decay rates of the viscoelastic drop oscillations converge faster to the linear values thanOℎ
the Newtonian. The effect of elasticity is the explanation for this behaviour.
The viscous influence on the viscoelastic drop-shape oscillations is reflected by the loss angle
_ϕ_ ( ) = arctan ( / ) of the motion also, as depicted in figure 8 _d_ . For the viscoelastic drops, the
angle is higher for the higher Ohnesorge number t αr αi 0, and, analogously to the decay rate, it
converges faster to the linear values for the viscoelastic than for the Newtonian drops. TheOℎ
linear values are 0.0385 and 0.073 for 0 = 0.16 and 1.238, respectively. The loss angles in the
first oscillation period, 1 for the former and the latter Oℎ 0 deviate from the respective linear
p
values by 29.3 and 21.9%, respectively.t Oℎ

##### (d) Power spectra


For the drop-shape oscillations studied, we are interested in the oscillation frequencies
contributing to the motions. For quantifying this, we compute the Fourier power spectra of the
frequency involved in the data traces in time. The finite-time frequency spectrum of a function
of time ( ) is given in an analytical form as
g t





2
t



( ) .
g t e [−][iαt] dt



2

( ) = ( ) . (4.3)
g [^] α 1 g t e [−][iαt] dt



1
t



1

The integral is defined for the time interval between 1 and 2. The Fourier transform becomes
more important for finding the frequencies, the more the quasi-periodicity influences the timet t
series in motions with large deformations. The quasi-periodicity depends on the liquid drop
material properties also. The frequencies are then not easily found from a ‘manual’ inspection
of the time series. Figure 9 shows the Fourier power spectra of the frequency for the motions
of the drop north pole in figure 8 _a_ . The data present power densities for the three different
liquid drops with the deformation parameter 0 = 0.35. The spectra exhibit zero values at zero
frequency because the underlying data represent the deviations of the north pole position fromη
its time average. The spectra are shown as obtained from an analysis of the first oscillation
period only and from analysing the period of time 0 ≤ ≤12 in order to see the evolution of the
frequency spectra with ongoing time during the oscillations.t
For the viscous and viscoelastic oscillating drops with = 0.048 and 0 = 1.238 studied,
the method yields similar power spectra of the oscillation frequency, as represented in theOℎ Oℎ
diagrams by the solid black and dotted blue lines, respectively. The spectra of the first oscillations exhibit peaks at frequencies deviating from the solutions of the characteristic equation
of the drop, i.e. from the linear oscillation frequency. The complex angular frequencies, which
represent the solutions of the characteristic equation, are 0.109 ± 2.827 and 0.21 ± 2.853 for
the viscoelastic drops with 0 = 0.16 and 0 = 1.238, respectively. The deviations of the peaki i
frequency for the first oscillation period from the linear oscillation frequency are 2.8 and 1.7%Oℎ Oℎ
for the viscoelastic drops with 0 = 0.16 and 0 = 1.238, respectively. Surprisingly, the peak
frequencies are not significantly affected by the change of the viscosity and the two viscoelasticOℎ Oℎ
timescales. They rather impact the decay rate, as seen in the previous paragraph. Even the
higher frequency corresponding to the mode = 4 is around 7.8 in the diagram for the first
oscillation and is even less represented in the power spectrum.k


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

( _a_ )


0.15


0.12


0.09



( _b_ )











0.06


0.03

|Col1|Col2|2.9|91|Viscous  − Oh = 0.048|Col6|Col7|Col8|Col9|Col10|Col11|Col12|
|---|---|---|---|---|---|---|---|---|---|---|---|
|||||Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||
|||||||||||||



0 3 6 9 12



1.6


1.2


0.8


0.4


|Col1|Col2|2.82|21 Viscous    − Oh = 0.048|Col5|Col6|Col7|Col8|Col9|
|---|---|---|---|---|---|---|---|---|
||||Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|Viscoelastic −_Oh_0 = 0.16<br>− _Oh_0 = 1.238|
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||
||||||||||



0 3 6 9



_α_ _α_



**Figure 9.** Squared real parts of the Fourier power spectra of the oscillation frequency for the drop north pole motion
corresponding to figure 5 for = 2 and 0 = 0.35: results in ( _a_ ) for the first oscillation and ( _b_ ) for the interval 0 ≤ ≤12.
The results include solutions up to the second order.m η t


The spectra evaluated for the time interval 0 ≤ ≤12 exhibit secondary peaks below the
dominant peak frequency. These frequencies may be interpreted as beat frequencies, whicht
develop with ongoing time and are therefore not visible in the first oscillation period. The
frequency of 1.6 appears in all three spectra, and the value around 0.4 is found in the spectra for
the two viscoelastic cases. The oscillations from the first period interact with those developing
later in time.

##### (e) State of onset of aperiodic motion


The state of onset of aperiodic behaviour of viscoelastic drops is studied. For Newtonian
liquids, in the mode of initial drop deformation = 2, the cut-off Ohnesorge number, where
the oscillation frequency vanishes,, = 0.76647 (notwithstanding the result by Prosperettim

[8] for the initial-value problem). Figure 2Oℎco N _b,d_ shows that viscoelastic drops may oscillate at far
higher 0 than the cut-off value of, for Newtonian drops. Viscoelastic drops with, e.g.
1 = 33.93 and the ratio Oℎ 2/ 1 = 0.03714 of the two Deborah numbers, oscillate at OhnesorgeOℎco N
numbers up to De 0 = 27. In the present section, first, the Ohnesorge number De De 0 of viscoelastic drops with varying Oℎ 1 and 2 at the cut-off state is studied. The data are obtainedOℎ
from solutions of the characteristic De Deequation (3.11) of the drop. Figure 10 _a_ depicts the critical
Ohnesorge number 0, of the viscoelastic drop as a function of the Deborah number 1 for
the value of 2/ 1 = 0.1 as an example. The data are represented empirically by an equation ofOℎ crit De
the form De De



0, = 0 + 1 1a2, (4.4)
Oℎ crit a a De



0, 0 1 1

where the values of the depend on the ratio 2/ 1. This dependency may be represented by
empirical functions given in the electronic supplementary materials. The data cover the rangeai De De
0.037 ≤ 2/ 1 ≤0.2. For a given Deborah number 1, the critical 0 decreases with increasing
2/ 1. The increasing deformation retardation time hinders the oscillatory motion, thereforeDe De De Oℎ
requiring a drop with lower De De 0.
Secondly, the cut-off state is characterized by the corresponding decay rate of theOℎ
drop surface deformation. The difference between the viscoelastic cut-off decay rate and
the Newtonian value 2,, = 2.72338 is depicted in figure 10 _b_ as a function of the
αr co N


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

`(` _a_ `)`


`(` _b_ `)`


```
15

13

11

 9

 7

2.7

1.8

0.9

```




|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|Col11|Col12|Col13|Col14|Col15|Col16|Col17|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
||||||||||||||||||
||||||||||||||||||
||||||||||||||||||
||||||||||||||||||
||||||||||||||||||
||||||||||||||||||
||||||||||||||||||

```
15 30 45 60 75 90

```

_De_
```
                  1

```

|Col1|Col2|Col3|Col4|Col5|Col6|Col7|Col8|Col9|Col10|Col11|Col12|Col13|Col14|Col15|Col16|Col17|Col18|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||
|||||||||||||||||||


```
0 0.03 0.06 0.09 0.12 0.15 0.18

```


_Ohv_ `–` _Ohco,N_

**Figure 10.** The state of viscoelastic drops at the cut-off of aperiodic motion: ( _a_ ) the critical Ohnesorge number 0, as a
function of the drop Deborah number 1, and ( _b_ ) the decay rate 2, of the deformation as a function of the modifiedOℎ crit
drop Ohnesorge number, both for the example that De 2/ 1 = 0.1 and the mode of initial drop deformation αr crit = 2.
v
Oℎ De De m


difference between the modified Ohnesorge number and the Newtonian value at cut-off
, = 0.76647 for the example 2/ 1 = 0.1. The values apply to the mode of initial deforma-Oℎv
tion Oℎco N = 2. These critical states to aperiodic drop behaviour are represented empirically by aDe De
power law, readingm



_Ohv_ `–` _Ohco,N_



2,  - 2,, = 1  -, b2. (4.5)
αr crit αr co N b Oℎv Oℎco N



b2.



2, 2,, 1,

The decay rate increases with the modified Ohnesorge number . The latter quantity, on the
other hand, contains the decay rate, so that the depicted relation is implicit. The coefficientsOℎv
1 and 2 again depend on the ratio of the Deborah numbers 2/ 1. The dependencies are
brepresented by empirical functions given in the electronic supplementary materials.b De De
With these findings, the oscillatory or aperiodic behaviour of a drop with a given size from
a rheologically characterized viscoelastic liquid, behaving according to the Oldroyd-B model,
may be predicted. Furthermore, the rate of aperiodic decay of drop surface deformations in
the cut-off state is obtained. Knowing all the material data 0,,, 1 and 2 of the drop liquid,
the viscoelastic timescale (or Deborah number) ratio, 2/ 1 ≡μ σ2/ ρ λ1, is known. Using this ratio,λ
the three coefficients in the empirical equation (4.4)λ for the critical Ohnesorge number mayλ De De
be determined, so that the relationship between the critical Ohnesorge number ai 0, and the
Deborah number 1 is known. This relationship allows the size to be determined, at whichOℎ crit
a drop from the given liquid assumes the critical state. Comparing this critical size with the sizeDe dcrit
of a given drop, allows the drop behaviour to be predicted: drops with - will perform
damped shape oscillations, while deformations of drops with d ≤ will decay aperiodically.d dcrit
The critical drop size now allows all the characteristic numbers d d0crit, 1 and 2 to be calculated. This allows the decay rate of the drop with the critical size to be determined from Oℎ De De equation


```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

(4.5). Knowing all the characteristic non-dimensional numbers governing the shape oscillations,
the actual oscillation frequency and decay rate of dampened shape oscillations for a drop with
oscillatory behaviour (i.e. with a size        - ), or the rate of aperiodic decay for a drop with
≤, are obtained as a solution of the characteristic d dcrit equation (3.11).
d dcrit

### 5. Conclusions

A weakly nonlinear analysis of axisymmetric shape oscillations of viscoelastic Oldroyd-B drops
in a dynamically inert environment was performed. The aim was to reveal the detailed time
behaviour of the nonlinear shape oscillations for the two-lobed mode = 2 of initial drop
deformation, carrying the derivations to the second order of approximation. The characteristicm
equation yields an infinite number of solutions for the complex angular oscillation frequency.
The appropriate pair of complex conjugate solutions is selected according to experimental
data for damped linear shape oscillations of aqueous poly(acrylamide) solution drops in an
acoustic levitator. The linear solutions from the theory agree well with the experimental
data. The theory reveals two nonlinear effects, which are a time asymmetry and a frequency
change for viscoelastic drops with increasing initial drop deformation. The excess time the
drop spends in the prolate mode is reduced by more than one-half for the whole range of
initial drop deformations in comparison with a viscous drop with the same first-order decay
rate. Comparing the trace of the drop north pole motion in time against the corresponding
viscous case, elastic effects are found to reduce the minimum drop aspect ratio in the oblate
state by around 12% against the viscous drop. The oscillation frequency shows the same trend
to vary with increasing initial drop deformation as known from both the inviscid and the
viscous cases. The present second-order theory predicts a slight increase in the frequency
with the initial drop deformation. A third-order approximation might add to this finding. The
variation of the oscillation frequency with the drop deformation makes it vary in time during
dampened oscillations, rendering the oscillations quasi-periodic. Likewise, the damping rate of
the oscillations varies with the drop deformation, converging to the linear value within the first
1.5 oscillation periods. The Fourier power spectra of the north pole motion show an increase in
the peak frequency within the first oscillation period against the solution of the characteristic
equation. Secondary peaks appear later in time owing to mode coupling and the consequent
beat between the different oscillation mode frequencies. The state of onset of aperiodic motion
is characterized by solutions of the characteristic equation, using empirical equations for the
cut-off Ohnesorge number 0, ( 1, 2) and the decay rate 2, . These findings allow the
oscillatory or aperiodic motion of a drop with a given size from a rheologically characterizedOℎ crit De De αr crit
Oldroyd-B liquid to be predicted.


Data accessibility. This article has no additional data.
Supplementary material is available online [47].
Declaration of AI use. We have not used AI-assisted technologies in creating this article.
Authors’ contributions. D.Z.: data curation, investigation, methodology, writing—original draft; G.B.: conceptuali
zation, funding acquisition, supervision, writing—review and editing.
Both authors gave final approval for publication and agreed to be held accountable for the work
performed therein.
Conflict of interest. We declare we have no competing interests.
Funding. The authors acknowledge Graz University of Technology for funding the open-access publication of
this article.
Acknowledgements. Financial support of this research project by the Austrian Science Fund (FWF) through
project number I3326-N32 in the DACH framework is gratefully acknowledged. The authors acknowledge
the excellent cooperation with Prof. Dr. M. Oberlack and his group at Darmstadt University of Technology
in Darmstadt, Germany.

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```



## References

1. Rayleigh JWS. 1879 On the capillary phenomena of jets. _Proc. R. Soc. Lond._ **29** [, 71–97. (doi:10.](http://dx.doi.org/10.1098/rspl.1879.0015)
[1098/rspl.1879.0015)](http://dx.doi.org/10.1098/rspl.1879.0015)
2. Lamb H. 1881 On the oscillations of a viscous spheroid. _Proc. Lond. Math. Soc._ **13**, 51–66. (doi:

[10.1112/plms/s1-13.1.51)](http://dx.doi.org/10.1112/plms/s1-13.1.51)
3. Lamb H. 1932 Hydrodynamics, 6th edn. Cambridge, UK: Cambridge University Press.
4. Chandrasekhar S. 1959 The oscillations of a viscous liquid globe. _Proc. Lond. Math. Soc._ **9**,
[141–149. (doi:10.1112/plms/s3-9.1.141)](http://dx.doi.org/10.1112/plms/s3-9.1.141)
5. Reid WH. 1960 The oscillations of a viscous liquid drop. _Quart. Appl. Math._ **18**, 86–89. (doi:

[10.1090/qam/114449)](http://dx.doi.org/10.1090/qam/114449)
6. Miller CA, Scriven LE. 1968 The oscillations of a fluid drop immersed in another fluid. _J._
_Fluid Mech._ **32** [, 417–435. (doi:10.1017/S0022112068000832)](http://dx.doi.org/10.1017/S0022112068000832)
7. Basaran OA, Scott TC, Byers CH. 1989 Drop oscillations in liquid-liquid systems. _AIChE J._
**35** [, 1263–1270. (doi:10.1002/aic.690350805)](http://dx.doi.org/10.1002/aic.690350805)
8. Prosperetti A. 1980 Free oscillations of drops and bubbles: the initial-value problem. _J. Fluid_
_Mech._ **100** [, 333–347. (doi:10.1017/S0022112080001188)](http://dx.doi.org/10.1017/S0022112080001188)
9. Khismatullin DB, Nadim A. 2001 Shape oscillations of a viscoelastic drop. _Phys. Rev. E_ **63**,
[061508. (doi:10.1103/PhysRevE.63.061508)](http://dx.doi.org/10.1103/PhysRevE.63.061508)
10. Chrispell JC, Cortez R, Khismatullin DB, Fauci LJ. 2011 Shape oscillations of a droplet in an
Oldroyd-B fluid. _Physica D. Nonlinear Phenomena._ **240** [, 1593–1601. (doi:10.1016/j.physd.2011.](http://dx.doi.org/10.1016/j.physd.2011.03.004)
[03.004)](http://dx.doi.org/10.1016/j.physd.2011.03.004)
11. Hoath SD etal. 2015 Oscillations of aqueous PEDOT:PSS fluid droplets and the properties of
complex fluids in drop-on-demand Inkjet printing. _J. Nonnewton. Fluid Mech._ **223**, 28–36.
[(doi:10.1016/j.jnnfm.2015.05.006)](http://dx.doi.org/10.1016/j.jnnfm.2015.05.006)
12. Li F, Yin XY, Yin XZ. 2019 Small-amplitude shape oscillation and linear instability of an
electrically charged viscoelastic liquid droplet. _J. Nonnewton. Fluid Mech._ **264** [, 85–97. (doi:10.](http://dx.doi.org/10.1016/j.jnnfm.2018.10.001)
[1016/j.jnnfm.2018.10.001)](http://dx.doi.org/10.1016/j.jnnfm.2018.10.001)
13. Trinh EH, Marston PL, Robey JL. 1988 Acoustic measurement of the surface tension of
levitated drops. _J. Colloid Interface Sci._ **124** [, 95–103. (doi:10.1016/0021-9797(88)90329-3)](http://dx.doi.org/10.1016/0021-9797(88)90329-3)
14. Hiller WJ, Kowalewski TA. 1989 Surface tension measurements by the oscillating droplet
method. _PCH. Physico. Chem. Hyd._ **11**, 103–112.
15. Stückrad B, Hiller WJ, Kowalewski TA. 1993 Measurement of dynamic surface tension by
the oscillating drop method. _Exp. Fluids._ **15** [, 332–340. (doi:10.1007/BF00223411)](http://dx.doi.org/10.1007/BF00223411)
16. Hsu CJ, Apfel RE. 1985 A technique for measuring interfacial tension by quadrupole
oscillation of drops. _J. Colloid Interface Sci._ **107** [, 467–476. (doi:10.1016/0021-9797(85)90199-7)](http://dx.doi.org/10.1016/0021-9797(85)90199-7)
17. Egry I, Lohöfer G, Seyhan I, Schneider S, Feuerbacher B. 1998 Viscosity of eutectic
Pd78Cu6Si16 measured by the oscillating drop technique in microgravity. _Appl. Phys. Lett._
**73** [, 462–463. (doi:10.1063/1.121900)](http://dx.doi.org/10.1063/1.121900)
18. Perez M, Salvo L, Suéry M, Bréchet Y, Papoular M. 2000 Contactless viscosity measurement
by oscillations of gas-levitated drops. _Phys. Rev. E_ **61** [, 2669–2675. (doi:10.1103/PhysRevE.61.](http://dx.doi.org/10.1103/PhysRevE.61.2669)
[2669)](http://dx.doi.org/10.1103/PhysRevE.61.2669)
19. Matsumoto T, Nakano T, Fujii H, Kamai M, Nogi K. 2002 Precise measurement of liquid
viscosity and surface tension with an improved oscillating drop method. _Phys. Rev. E_ **65**,
[031201. (doi:10.1103/PhysRevE.65.031201)](http://dx.doi.org/10.1103/PhysRevE.65.031201)
20. Gottier GN, Amundson NR, Flumerfelt RW. 1986 Transient dilation of bubbles and drops:
theoretical basis for dynamic interfacial measurements. _J. Colloid Interface Sci._ **114**, 106–130.
[(doi:10.1016/0021-9797(86)90244-4)](http://dx.doi.org/10.1016/0021-9797(86)90244-4)
21. Tian Y, Holt RG, Apfel RE. 1995 Investigations of liquid surface rheology of surfactant
solutions by droplet shape oscillations: theory. _Phys. Fluids (1994)._ **7** [, 2938–2949. (doi:10.](http://dx.doi.org/10.1063/1.868671)
[1063/1.868671)](http://dx.doi.org/10.1063/1.868671)
22. Apfel RE etal. 1997 Free oscillations and surfactant studies of superdeformed drops in
microgravity. _Phys. Rev. Lett._ **78** [, 1912–1915. (doi:10.1103/PhysRevLett.78.1912)](http://dx.doi.org/10.1103/PhysRevLett.78.1912)
23. Tian YR, Holt RG, Apfel RE. 1997 Investigation of liquid surface rheology of surfactant
solutions by droplet shape oscillations: experiments. _J. Colloid Interface Sci._ **187** [, 1–10. (doi:10.](http://dx.doi.org/10.1006/jcis.1996.4698)
[1006/jcis.1996.4698)](http://dx.doi.org/10.1006/jcis.1996.4698)

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```



24. Aske N, Orr R, Sjöblom J. 2002 Dilatational elasticity moduli of water-crude oil interfaces
using the oscillating pendant drop. _J. Dispers. Sci. Technol._ **23** [, 809–825. (doi:10.1081/DIS-](http://dx.doi.org/10.1081/DIS-120015978)
[120015978)](http://dx.doi.org/10.1081/DIS-120015978)
25. Lu HL, Apfel RE. 1991 Shape oscillations of drops in the presence of surfactants. _J. Fluid_
_Mech._ **222** [, 351. (doi:10.1017/S0022112091001131)](http://dx.doi.org/10.1017/S0022112091001131)
26. Zhang XG, Harris MT, Basaran OA. 1994 Measurement of dynamic surface tension by a
growing drop technique. _J. Colloid Interface Sci._ **168** [, 47–60. (doi:10.1006/jcis.1994.1392)](http://dx.doi.org/10.1006/jcis.1994.1392)
27. Kovalchuk VI, Krägel J, Aksenenko EV, Loglio G, Liggieri L. 2001 Oscillating bubble and
drop techniques. In Novel methods to study interfacial layers, pp. 485–516. Amsterdam,
Netherlands: Elsevier.
28. Ravera F, Loglio G, Kovalchuk VI. 2010 Interfacial dilational rheology by oscillating bubble/
drop methods. _Curr. Opin. Colloid Interface Sci._ **15** [, 217–228. (doi:10.1016/j.cocis.2010.04.001)](http://dx.doi.org/10.1016/j.cocis.2010.04.001)
29. Yang L et al. 2014 Determination of dynamic surface tension and viscosity of nonNewtonian fluids from drop oscillations. _Phys. Fluids (1994)._ **26** [, 113103. (doi:10.1063/1.](http://dx.doi.org/10.1063/1.4901823)
[4901823)](http://dx.doi.org/10.1063/1.4901823)
30. Brenn G, Plohl G. 2015 The oscillating drop method for measuring the deformation
retardation time of viscoelastic liquids. _J. Nonnewton. Fluid Mech._ **223** [, 88–97. (doi:10.1016/j.](http://dx.doi.org/10.1016/j.jnnfm.2015.05.011)
[jnnfm.2015.05.011)](http://dx.doi.org/10.1016/j.jnnfm.2015.05.011)
31. Brenn G, Teichtmeister S. 2013 Linear shape oscillations and polymeric time scales of
viscoelastic drops. _J. Fluid Mech._ **733** [, 504–527. (doi:10.1017/jfm.2013.452)](http://dx.doi.org/10.1017/jfm.2013.452)
32. Ponce-Torres A, Montanero JM, Herrada MA, Vega EJ, Vega JM. 2017 Influence of the
surface viscosity on the breakup of a surfactant-laden drop. _Phys. Rev. Lett._ **118**, 024501. (doi:
[10.1103/PhysRevLett.118.024501)](http://dx.doi.org/10.1103/PhysRevLett.118.024501)
33. Hosseinzadeh VA, Brugnara C, Holt RG. 2018 Shape oscillations of single blood drops:
applications to human blood and sickle cell disease. _Sci. Rep._ **8** [, 16794. (doi:10.1038/s41598-](http://dx.doi.org/10.1038/s41598-018-34600-7)
[018-34600-7)](http://dx.doi.org/10.1038/s41598-018-34600-7)
34. Tamim SI, Bostwick JB. 2021 Oscillations of a soft viscoelastic drop. _NPJ. Microgravity._ **7**, 42.
[(doi:10.1038/s41526-021-00169-1)](http://dx.doi.org/10.1038/s41526-021-00169-1)
35. Zrnić D, Brenn G. 2021 Weakly nonlinear shape oscillations of inviscid drops. _J. Fluid Mech._
**923** [, A9. (doi:10.1017/jfm.2021.568)](http://dx.doi.org/10.1017/jfm.2021.568)
36. Zrnić D, Berglez P, Brenn G. 2022 Weakly nonlinear shape oscillations of a Newtonian drop.
_Phys. Fluids (1994)._ **34** [, 043103. (doi:10.1063/5.0085070)](http://dx.doi.org/10.1063/5.0085070)
37. Calvert P. 2001 Inkjet printing for materials and devices. _Chem. Mater._ **13** [, 3299–3305. (doi:10.](http://dx.doi.org/10.1021/cm0101632)
[1021/cm0101632)](http://dx.doi.org/10.1021/cm0101632)
38. Pandey K, Prabhakaran D, Basu S. 2019 Review of transport processes and particle selfassembly in acoustically levitated nanofluid droplets. _Phys. Fluids (1994)._ **31** [, 112102. (doi:10.](http://dx.doi.org/10.1063/1.5125059)
[1063/1.5125059)](http://dx.doi.org/10.1063/1.5125059)
39. Kim J. 2007 Spray cooling heat transfer: the state of the art. _Int. J. Heat Fluid Flow._ **28**, 753–
[767. (doi:10.1016/j.ijheatfluidflow.2006.09.003)](http://dx.doi.org/10.1016/j.ijheatfluidflow.2006.09.003)
40. Breitenbach J, Roisman IV, Tropea C. 2018 From drop impact physics to spray cooling
models: a critical review. _Exp. Fluids_ **59** [, 55. (doi:10.1007/s00348-018-2514-3)](http://dx.doi.org/10.1007/s00348-018-2514-3)
41. Kamperman T, Trikalitis VD, Karperien M, Visser CW, Leijten J. 2018 Ultrahigh-throughput
production of monodisperse and multifunctional janus microparticles using in-air
microfluidics. _ACS Appl. Mater. Interfaces_ **10** [, 23433–23438. (doi:10.1021/acsami.8b05227)](http://dx.doi.org/10.1021/acsami.8b05227)
42. Baumgartner D, Brenn G, Planchette C. 2020 Effects of viscosity on liquid structures
produced by in-air microfluidics. _Phys. Rev. Fluids_ **5** [, 103602. (doi:10.1103/PhysRevFluids.5.](http://dx.doi.org/10.1103/PhysRevFluids.5.103602)
[103602)](http://dx.doi.org/10.1103/PhysRevFluids.5.103602)
43. Smuda M, Kummer F, Oberlack M, Zrnić D, Brenn G. 2024 From weakly to strongly
nonlinear viscous drop shape oscillations: an analytical and numerical study. _Phys. Rev._
_Fluids._ **9** [, 1–30. (doi:10.1103/PhysRevFluids.9.063601)](http://dx.doi.org/10.1103/PhysRevFluids.9.063601)
44. Bird RB, Stewart WE, Lightfoot EN. 1962 Transport phenomena. New York, NY: John Wiley
& Sons.
45. Tomotika S. 1935 On the instability of a cylindrical thread of a viscous liquid surrounded by
another viscous fluid. _Proc. R. Soc. Lond. A._ **150** [, 322–337. (doi:10.1098/rspa.1935.0104)](http://dx.doi.org/10.1098/rspa.1935.0104)
46. Renoult MC, Brenn G, Plohl G, Mutabazi I. 2018 Weakly nonlinear instability of a
Newtonian liquid jet. _J. Fluid Mech._ **856** [, 169–201. (doi:10.1017/jfm.2018.677)](http://dx.doi.org/10.1017/jfm.2018.677)

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```



47. Zrnić D, Brenn G. 2024 Data from: weakly nonlinear shape oscillations of viscoelastic drops.
[Figshare. (doi:10.6084/m9.figshare.c.7399609)](http://dx.doi.org/10.6084/m9.figshare.c.7399609)

```
Downloaded from http://royalsocietypublishing.org/rspa/article-pdf/doi/10.1098/rspa.2023.0887/512897/rspa.2023.0887.pdf
by guest
on 26 May 2026

```

