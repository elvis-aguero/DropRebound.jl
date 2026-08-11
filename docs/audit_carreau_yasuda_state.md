# Audit: state of the Carreau-Yasuda extension, and what to adapt from Reid

Date: 2026-08-01. Scope: every artifact that claims to extend Reid (1960)
to shear-thinning fluids — `src/st_extension.jl`,
`src/st_exact_extension.jl`, the four `derivations/carreau*.jl`
scripts, `derivations/cross_fluid_derivation.jl`,
`docs/section_carreau.tex`, `docs/section_carreauYasuda_firstprinciples.tex`,
the `docs/src/carreau_yasuda_fp/` chapters (since consolidated away), and the CY
test files.

Purpose: decide what survives into a first-principles model and what has to
be rebuilt.

---

## 0. Verdict

The Newtonian foundation is sound. The shear-thinning layer on top of it is
not a derivation — it is a chain of modeling choices, each individually
defensible, none derived, with self-checks that cannot fail. It currently
passes its own empirical target by a median metric that conceals a large
failure tail.

Three findings drive everything below.

1. **The strain-rate field is built from the wrong velocity.** Production
   uses inviscid potential flow `φ ∝ r^l P_l`, not Reid's actual viscous
   eigenfunction `U(x) = C·x·j_l(qx) + Π₀·x^{l+1}`. There is no Bessel
   function anywhere in the shear-thinning path.
2. **A time-varying `Oh_eff(t)` is substituted into eigenvalues derived
   for constant viscosity, with no justification stated anywhere.** The
   substitution is not merely unargued; the assumption it needs is violated
   by two to three orders of magnitude.
3. **Nothing tests the physics.** Zero of 49 assertions in the CY test
   files constrain the model against anything independent of its own
   algebra, and the only artifact that touches experimental data is run by
   neither `runtests.jl` nor CI.

A fourth finding concerns priority between the two spatial/temporal defects.
The superseded Carreau-Yasuda chapters named spatial homogenization "the
genuinely open question" and said direct evaluation at `Oh_eff(t)` "was never
the weak point". The second half of that is wrong — the temporal closure is a
measured ~80% error. But an earlier draft of this audit then over-corrected,
calling spatial homogenization a mere ~10% effect on the strength of a
single-mode measurement. Measured on realistic multi-mode states it is a
**leading-order** error. Both defects are real and neither dominates — see §3,
which records the correction.

---

## 1. What is solid, and carries over unchanged

### 1.1 Reid's Newtonian theory

`src/reid.jl` is the strongest part of the repo. It is the only
component with genuine, independent validation, in `test/test_reid.jl`:

- `test_reid.jl:27-40` — convergence to Lamb's `λ_l = Oh(l-1)(2l+1)` as
  `Oh → 0`, with the error required to *decrease monotonically* across
  three decades and finish below 1%.
- `test_reid.jl:42-53` — agreement with Molaček & Bush's published
  high-`Oh` limit `D_l = (l-1)(2l²+4l+3)/(l²(2l+1))` to 0.1%, at five
  modes.

Two independent analytic limits, both external to the code. Keep all of it.
The continuation strategy in `dominant_root` and `build_reid_table` (branch
tracking by small steps in `Oh`) is load-bearing and well documented; do not
simplify it.

### 1.2 What Reid gives a *generalized Newtonian* fluid for free

This is the part worth stating precisely, because the repo currently
under-claims it. the Reid chapter (since consolidated into
`derivations/reid1960_full_derivation.jl`) said
that a shear-thinning model "changes the momentum equation AND both stress
boundary conditions". That is too pessimistic. Working through it with
`η = η(γ̇)` in place of constant `μ`:

- **BC1 (kinematic), `u_r = ∂ζ/∂t` at the surface.** Contains no viscosity
  at all. Unchanged.
- **BC2 (tangential stress).** Written in
  in Reid's boundary-condition section as
  `τ_rθ = μ[ r ∂_r(u_θ/r) + (1/r)∂_θ u_r ] = 0`. The bracket is exactly
  `2e_rθ` — purely kinematic. For a free surface with no exterior fluid the
  condition is `τ_rθ = 2η e_rθ = 0`, and since `η(γ̇) ≥ η_∞ > 0` everywhere,
  this gives `e_rθ = 0` **regardless of whether η is constant**. BC2 is
  rheology-agnostic for *any* generalized Newtonian fluid, and therefore the
  whole `τ_rθ = 0 → L₂[U] = 0` chain already derived in
  `derivations/reid1960_full_derivation.jl` (Assertions 12-15) carries
  over verbatim. This is a real saving and should be claimed.
- **BC3 (normal stress).** `-p_rr = p + δp - 2μ ∂u_r/∂r`. Here η is genuinely
  multiplicative, so it becomes the *surface* value `η_s ≡ η(γ̇|_{r=R})`.
  The form is unchanged; one coefficient becomes state-dependent.
- **Momentum equation.** This is where the real change lives. Cauchy plus
  incompressibility gives
  `∇·(2η e) = 2(∇η)·e + η∇²u`,
  so Reid's `ν∇²u` acquires exactly one extra term, `2(∇η)·e`. Reid's
  derivation drops it, correctly, because for him `∇η ≡ 0`.

So the honest scope of "adapt Reid for Carreau-Yasuda" is: **one new term in
the momentum equation, one coefficient in BC3, and nothing else.** That is a
far more tractable problem than the current docs imply, and it makes the
"spatial homogenization" question concrete — collapsing `η(x,θ)` to a scalar
*is* the act of dropping `2(∇η)·e`. They are the same approximation.

`2(∇η)·e` appears nowhere in this repo. That is the gap.

### 1.3 General results worth keeping

From `docs/section_carreauYasuda_firstprinciples.tex` and its backing
script — all three are genuine, general, and independent of any fluid's
parameters:

- The **adjoint sensitivity shortcut** for Reid's operator,
  `Y'(1) = (1/j_l(q₀))∫₀¹ x j_l(q₀x) RHS(x) dx`, verified against direct RK4
  shooting.
- The **viscous strain-rate tensor** built from Reid's actual `U(x)`, with
  incompressibility verified for `l = 2..8`.
- The **period-π lemma**: `S = √(2 e_ij e_ij)` from any single-frequency
  axisymmetric poloidal field is exactly period-π, so it has identically
  zero content at the mode's own frequency. This is correct and it is the
  key structural fact for §3.

From the legacy scripts, these survive independent of the `a=2` and
small-`ε_ST` assumptions:

- `C(a) = (2/√π)·Γ((a+3)/2)/Γ((a+4)/2)` — the Wallis ratio
  `⟨|sin|^{a+2}⟩/⟨sin²⟩`, exact for any real `a > -2`
  (`carreau_yasuda_derivation.jl:397-415`).
- `norm_p(p) = Γ(p+½)/(√π Γ(p+1)) = ⟨|sin|^{2p}⟩` (`:708`).
- The force identity `d/dḃ[-|ḃ|^{a+2}] = -(a+2)|ḃ|^a ḃ` (`:350-365`) — note
  the `(a+2)`, see §2.4.
- The Reid base-flow machinery in `carreau_yasuda_derivation.jl:170-227`:
  the ODE verification for `U`, the 2×2 BC solve reproducing Reid eq. (38),
  and `g = (2f + x f')/(l(l+1))`.
- The strain-component *expressions* at `:239-243` / `:660-663` — correct
  given a radial amplitude `f`; only the `f` fed into them is wrong (§2.1).
- The Gabbard energy result `A₂ = √(5We/12)` plus its genuine `solve_drop!`
  cross-check (`cross_fluid_derivation.jl:174-211`) — rheology-independent,
  and the only place in the legacy pair that meets the CLAUDE.md contract.
- Infrastructure: `legendre_arrays`, the dependency-free Gauss-Legendre
  rule, `reid_char`/`find_eigenvalue`.

---

## 2. What is broken

### 2.1 Wrong velocity field (two separate instances)

**Production.** `st_exact_extension.jl:39-52` builds the strain basis from
"mode `l`'s own **potential-flow** field", with radial dependence `r^(l-2)`
applied at `:163-167`. Reid's rotational (Bessel) part is absent.

The repo's own first-principles script already says this is wrong —
`carreauYasuda_firstprinciples_derivation.jl:319-324`:

> the existing heuristic scripts estimate shear rate from the INVISCID
> potential-flow mode shape […] Reid's damping normalization comes from the
> homogeneous (Bessel) part of `U(x)`, which IS the viscous correction to
> potential flow.

Size of the error, measured as `‖U(x)/U(1) − x^{l+1}‖`:

| Oh | l=2 | l=4 |
|---|---|---|
| 57.4 | 0.147 | 0.338 |
| 0.05 | 0.057 | 0.095 |

15–34% at the operating point, and it does **not** shrink at large `Oh`
(Reid's `q` saturates near 2.67 for `l=2`). Since the fluid sits deep in the
power-law branch, `Oh_eff ∝ S^{-0.743}`, so a multiplicative error `c` in
`S` propagates as `c^{-0.743}` — a factor-2 error in `S` moves `Oh_eff` by
~37%. Nothing cancels it.

**Legacy, separately.** `carreau_yasuda_derivation.jl` feeds the *poloidal
scalar* `U` into the strain builder where the *radial amplitude* `f = U/x²`
belongs (`:282`, `:718-722`), contradicting this repo's own
`reid1960_full_derivation.jl:312,317` (`u_r ∝ U(x)/x²`). The check is
physical: the inviscid `l=2` mode is uniform straining, so `γ̇²` must be
constant. Under `f = U` it is not; under `f = U/x²` it is exactly 3.
Consequence: `H(θ) = 3cos⁴θ + 11cos²θ + 13` is an artifact, and **every**
numerical `Γ_l^(a)` in the repo is wrong — including `Γ₂ = 1783566/385`,
which is hardcoded into `test_carreau.jl:5` and `test_carreau_yasuda.jl:9`.

These are two distinct bugs. Fixing one does not fix the other.

### 2.2 Instantaneous `Oh_eff` fed into constant-viscosity eigenvalues

`carreau_yasuda_nonperturbative_derivation.jl:20-25` argues that because
Carreau-Yasuda has no relaxation time, evaluating at the instantaneous state
"is exact for this constitutive law, not merely a convenient approximation."

The premise is true; the conclusion does not follow. The constitutive law is
memoryless — but `λ_l` and `ω_l²` are eigenvalues of a boundary-value
problem solved over the mode's entire spatial and temporal structure. They
are not pointwise functions of the state. Substituting `Oh_eff(t)` into them
is a quasi-static approximation, and no quasi-static, adiabatic, WKB, or
separation-of-timescales argument appears anywhere in the repo (verified by
grep across all `.jl`, `.tex`, `.md`).

The assumption is violated as badly as it can be. Measured directly on the
production code (real fluid, `l=2`, modal velocity amplitude 0.05):

| phase/π | Ȧ₂ | Oh_eff | λ₂ |
|---|---|---|---|
| 0.0 | 0.05 | 0.190 | 0.72 |
| 0.25 | 0.0354 | 0.238 | 0.88 |
| 0.5 | ~0 | **57.37** | **203.9** |
| 0.75 | −0.0354 | 0.238 | 0.88 |

`λ₂` swings by a factor of **284**, twice per period, reaching exactly `Oh₀`
at every velocity zero-crossing. `ω₂²` moves only ~3%, so the entire effect
is in the damping — the quantity whose behaviour the repo's own chapter 1
already shows is the badly-conditioned one.

This is not a narrow spike. At `Ȧ₂ = 0.001` (2% of peak) `λ₂` is already
10.3, still 14× its peak-velocity value.

Two independent reasons the substitution is not licensed:

- **Timescale.** Quasi-static substitution needs `Oh_eff` to vary slowly
  compared with the eigenvalue's own timescale. By the period-π lemma
  (already proved in this repo), `Oh_eff` varies with period π — *exactly
  twice* the mode's own frequency, the fastest modulation possible. The
  ratio is O(1) by construction, at every amplitude. There is no limit in
  which this becomes small.
- **Boundary layer.** `λ_l(Oh)` encodes the structure of the vorticity
  layer, thickness `δ ~ √(ν/ω)`, whose equilibration time `δ²/ν ~ 1/ω` is by
  construction the same order as the modulation period. The layer physically
  cannot re-form twice per period, so the flow field cannot adiabatically
  track a viscosity swinging over three orders of magnitude.

**Consequence, measured.** Comparing the cycle-averaged damping the code
produces against an energy-consistent RMS-shear closure:

| amplitude | ⟨λ(Oh_eff(t))⟩ | λ(Oh_eff(rms)) | ratio |
|---|---|---|---|
| 0.2 | 0.681 | 0.412 | 1.65 |
| 0.05 | 1.583 | 0.883 | 1.79 |
| 0.01 | 4.542 | 2.578 | 1.76 |

The code over-damps by ~75–80%, robustly across amplitudes. Worse, the
integral is dominated by the region where it is least defensible: **the 4%
of the cycle nearest a velocity zero-crossing carries 40% of the
cycle-averaged damping** (20% of the cycle carries 70%). The answer is
controlled by a modeling artifact, not by the bulk of the motion.

### 2.3 The `Oh_eff` closure is postulated

`st_exact_extension.jl:274-281` computes
`Oh_eff_l = Oh₀ · ⟨μ_eff/μ₀⟩` weighted by `|w_l|`, `w_l = ∂S²/∂Ȧ_l`.

- `w_l` is **not** `∂D_total/∂Ȧ_l`, though the derivation header sets up
  `D_total = ∫ μ_eff(S) S² dV` and calls `w_l` a dissipation sensitivity.
  The chain-rule term through `μ_eff`'s own `S`-dependence is dropped. Since
  `S μ'(S)/μ → −a·ε_ST = −0.743`, the omitted term is ~37% of the retained
  one and opposite in sign. `D_total` is never differentiated anywhere.
- The `|·|` is a modeling choice made for numerical convenience
  (`multimode:31`), not physics: regions where mode `l` opposes the total
  field have their contribution flipped rather than subtracted.
- The header's claim that using `Oh_eff` "only ever as a RATIO" makes the
  known normalization defect "cancel rather than propagating into the
  physics" is **wrong as stated**. The ratio cancels the *scale* of the
  weight. It cannot cancel an error in `S`, because `S` enters the
  *argument* of a nonlinear function. And the defect it is meant to excuse
  is not a scalar: the basis's raw dissipation integral departs from
  Lamb's `(l-1)(2l+1)` by a factor spanning **13.6× across l=2..8**, i.e. a
  shape error mode by mode.
- The assertion offered as proof that the ratio is safe (`multimode:236-245`)
  sets `λ_c = 0`, making `μ_ratio ≡ 1` pointwise, then checks that the
  weighted mean of the constant 1 is 1. Any weight function passes.
- The single-mode reduction check (`multimode:222-233`, "<1%") is run at
  `l = 2` only — the one mode where `S` is spatially uniform, so the two
  formulas coincide identically by construction. Sweeping `l` gives 4.2%
  (l=3), 9.4% (l=4), 14.1% (l=5), 24.5% (l=8).

### 2.4 The legacy weakly-nonlinear result does not generalize to non-integer `a`

`carreau_yasuda_derivation.jl` proves its expansion is the correct Taylor
term by repeated differentiation in `ε` (`:112-128`), tested at `a = 2, 3, 4`
— a technique only defined for integer `a`. For the fitted `a = 0.743`,
`(1+ε^a)^p` is not analytic at `ε = 0` and there is no Taylor series to be
the correct term of.

`st_extension.jl:28-31` nevertheless applies the formula with `a` symbolic.
Four specific gaps:

1. **A known `a`-dependent factor is computed and discarded.** ASSERTION 16
   (`:350-365`) proves the generalized force carries `(a+2)`; the Newtonian
   one carries `2`; the ratio `(a+2)/2` is 2 at `a=2` but 1.37 at `a=0.743`.
   Four lines later `:369` writes the coefficient as exactly 1.
2. Every validating assertion is a one-point match at `a = 2` against a
   prior Python notebook. Any `F(a)` with `F(2)=1` passes all of them.
3. **The multimode sum has no derivation.** `Γ_l^(a)` is defined for a
   single mode; `st_extension.jl:31` sums `Σ_k Γ_eff[k]|Ȧ_k|^a`. But `γ̇²`
   of a superposition carries cross-terms, and raising to `(a+2)/2` is not
   additive over modes for any `a`.
4. The extra `Λ^a = (λ_c σ_{l;0})^a` factor is asserted
   (`section_carreau.tex:118`) and never carried through a dimensional
   argument; `γ̇` is already normalized by `ḃ_l/R`.

Validity conditions, all violated by the fitted fluid: `ε_ST Λ² a² ≪ 1`
(`section_carreau.tex:171`) fails by ~6 orders of magnitude, and
`st_extension.jl:34` has no floor, so the correction flips sign and becomes
*anti-damping*; `Oh ≪ 1` for the averaging (actual `Oh₀ = 57.4`, grossly
overdamped, no oscillation to average); no contact (production applies it
during impact anyway).

`section_carreau.tex:158` is also simply wrong: it attributes the factor 3/4
to "the time average of `cos²ωt sin²ωt`", which is 1/8. The correct origin is
`⟨sin⁴⟩/⟨sin²⟩ = 3/4`, which is what the script computes at `:421-423`.

### 2.5 Parameter conflation: `ε_ST` used in two incompatible roles

`st_exact_extension.jl:87` documents `eps_ST = (1-n)/a` — a pure exponent,
with the amplitude `Δ` factored out separately as `(1 - eta_inf_ratio)`. But
`cross_fluid_derivation.jl:345` defines the conversion recipe as
`eps_ST = [(μ₀-μ_∞)/μ₀]·(1-n)/a = Δ·(1-n)/a`, and
`scripts/validate_shear_thinning.jl:39` passes exactly that (`Δ = 0.99956`)
into the exponent slot.

Under the Cross→CY mapping `(1-n)/a = 1` identically, so the exponent
*should* be 1.0. It is 0.99956. For this fluid the induced error is ~0.2% —
harmless. But the two definitions are not the same quantity, and for a fluid
with `Δ = 0.5` the exponent would silently become 0.5, changing `μ_eff/μ₀`
by a factor of ~15 at operating shear. This is a live landmine, not a
present error.

Related: `cross_fluid_derivation.jl:344-346` recommends `n = 1-m`, which
gives `n ≤ 0` for every `m ≥ 1` — outside the `n ∈ (0,1]` domain that
`carreau_yasuda_derivation.jl:86` requires, for three of the four `m` values
its own prototype runs.

### 2.6 Numerical hacks presented as physics

`_ringing_outlier_mask` (`st_exact_extension.jl:211-229`) deletes modal
amplitude from the shear field using three tuned constants
(`OUTLIER_FACTOR = 20.0`, `RINGING_NOISE_FLOOR = 1e-9`,
`candidate_fraction = 0.1`), fitted against three cases, two of them
synthetic. Its six tests are regressions on those same cases. The derivation
itself documents that an earlier version of this filter misclassified real
physics and "cut predicted CoR roughly in half".

### 2.7 The "improvement" chain made predictions worse, and does not say so

From `carreau_yasuda_multimode_derivation.jl:295-297` and `:450-451`:
contact-time median relative error went **16.8% (self-only) → 62%
(multi-mode) → 29% (multi-mode + dealiasing)**. The chain never returns to
the baseline it replaced. After the `η_∞` floor was added, no number is
reported at all.

`:373-375` also reports 0.74 (no dealiasing) vs 0.21 (index filter) against
a measured 0.82 — i.e. the cited evidence shows *no filter* is closest to
the measurement, which argues against the mechanism the section introduces.
The script does not engage with this.

### 2.8 No physics tests

Across `test_carreau.jl`, `test_carreau_yasuda.jl`,
`test_carreau_yasuda_exact.jl`: **49 assertions, 0 physics.** 15 tautologies,
25 smoke checks, 9 regressions against the code's own prior output.

- No test compares against experimental data or a published shear-thinning
  result.
- No test would catch the wrong velocity field — and
  `test_carreau_yasuda_exact.jl:12` would **fail if it were corrected**,
  because it pins the potential-flow `K_l` table in place.
- The Newtonian-limit tests in `test_carreau.jl:14-19` and
  `test_carreau_yasuda.jl:51-57` verify the early-return guard
  `st.eps_ST == 0.0 && return` (`st_extension.jl:23`), not a limit.
  `test_carreau_yasuda_exact.jl:16-40` is genuine but uses `viscous=:lamb`;
  **the production `:reid` path has no Newtonian-limit test and no
  finite-difference Jacobian check anywhere.**
- `test_carreau_yasuda.jl:7-11` hardcodes `439.1737`, `1416.7890`,
  `51233.2220` attributed to a Python notebook (which CLAUDE.md forbids),
  reproduced nowhere in the repo, and never asserted against — any could be
  off by 10× undetected. They are also invalid regardless, per §2.1.
- `test_carreau_yasuda_exact.jl:115-117` references a `dealiasing_cutoff`
  function and a "dedicated dealiasing test below" that do not exist.

`scripts/validate_shear_thinning.jl` — the only artifact touching the
experimental data in `derivations/data/metrics_3000ppm.csv` — contains
zero assertions, is absent from `runtests.jl`, and is referenced by no CI
job.

### 2.9 CLAUDE.md rheology-contract violations

Item 2 requires a live cross-check against the running solver.

- `carreau_yasuda_derivation.jl` has `using DropSolver` at `:42` and never
  calls it. All 30 assertions are against its own algebra or a prior Python
  notebook. **Fails.**
- `carreau_yasuda_nonperturbative_derivation.jl` labels Section 5 a "live
  cross-check"; it builds `SimConstants` at `:358`, never uses it, and
  hand-rolls a two-line Euler step. `solve_drop!` is never called. It also
  runs `viscous=:lamb`, so `ω²` is the *constant* inviscid value — the run
  "validating" the model is one in which the thing Section 4 argues must
  vary does not vary. **Fails, and is mislabeled.**
- `carreau_yasuda_multimode_derivation.jl` and
  `cross_fluid_derivation.jl` do call `solve_drop!`. **Pass.**

---

## 3. Spatial homogenization: an earlier estimate here was wrong

**This section previously argued that spatial homogenization is a ~10% effect
and that the temporal closure (~80%) dominates it. The first half of that is
wrong, and it is corrected here rather than quietly deleted.**

The ~10% figure came from evaluating the Carreau-Yasuda law on Reid's viscous
mode with a **single active mode** `l=2`, where the period-averaged viscosity
varies by only 1.1-1.2x across the drop. That measurement is correct, and it
is irrelevant, because a drop after impact has dozens of modes live.

Measured properly, on Reid's actual eigenmodes with realistic multi-mode
spectra, the smallest bandwidth `L_eta` holding the discarded coupling below
1% is:

| modal spectrum | `L_eta` |
|:--|--:|
| single mode `l=2` | 8 |
| `~ l^-2`, M=30 | 32 |
| `~ l^-2`, M=50 | 50 |
| `~ l^-1`, M=30 | 126 |
| real solver state, M=30 | 133 |
| real solver state, M=20 | 139 |

`L_eta >= M` always, and `L_eta >> M` for real states. Two consequences:

1. **Banding the mode-coupling matrix buys nothing.** At M=30 you would need a
   half-bandwidth of ~133 to keep 1% -- a full matrix.
2. **Spatial homogenization is a leading-order error, not a 10% one.** At the
   physical operating point, `|eta_1|/|eta_0| ~ 1.3` and `|eta_2|/|eta_0| ~ 0.2`:
   the angular structure of the viscosity is *comparable to or larger than its
   mean*. Every model this repo has implemented collapses that onto one scalar
   per mode.

**Why the spectrum decays so slowly**, and this is physical rather than
numerical: `L_eta` tracks the *contrast* of the viscosity field. The contrast
comes from near-nodal points of the superposed strain field, where the shear
rate collapses and the fluid snaps back toward unthinned. That is a
near-discontinuity in `cos(theta)`, and a near-discontinuity has a Legendre
spectrum decaying only algebraically. One mode has few nodes; a dozen modes
beating against each other have a dense nodal set.

A methodological note worth carrying forward: the natural metric -- fraction of
`sum |eta_l'|^2` below `L_eta` -- says `L_eta = 2..5` suffices everywhere, and is
**misleading**. Coupling coefficients enter the matrix additively, not in
quadrature, so the controlling quantity is the summed magnitude
`sum_{l' > L} |eta_l'| / |eta_0|`. The spectrum has a long plateau carrying
little power but large summed magnitude.

**Caveat.** Coefficients beyond `l' = 140` were not measured, so every entry
above is a *lower* bound on the discarded coupling.

So the two defects are not 80% vs 10%. The temporal closure is a measured ~80%
error with a proof-backed fix; the spatial collapse is a leading-order error
with no cheap fix at all. Both are real, and neither dominates the other.

## 4. Empirical state

`scripts/validate_shear_thinning.jl 20 1`, current `main`:

```
median relative error, CoR: 0.1845  (20/20 solved)
median relative error, tc:  0.1325  (20/20 solved)
```

Both under the stated 20% target. But the median conceals a bimodal
distribution:

- **5/20 contact-time predictions are wrong by ~10×** — at We = 0.32, 0.91,
  1.90, 2.25, 2.84, predicted `t_c` collapses to 6e-5…1.4e-4 s against
  ~1e-3 s measured (`tc_err` 0.89–0.96). The failures are concentrated at
  high We.
- **6/20 CoR predictions are wrong by 33–94%** — worst at low We
  (We = 0.0064 predicts CoR 0.054 against 0.825 measured).

The model is roughly right in the middle of the We range and badly wrong in
both wings — the signature of a mis-specified closure that has been tuned to
a median. Any future change must be judged on the tail, not the median,
and CoR and `t_c` must be reported per-We-decade.

---

## 5. Recommended path

An earlier draft of this section ranked the fixes by proximity to the existing
code, which put a cheap patch first. That is anchor bias, and the ranking below
is by sharpness instead.

### Already closed

- **Falsifiable tests on the production path.** A Newtonian-limit test and a
  finite-difference Jacobian check now exist for `viscous=:reid` (previously
  both existed only for `:lamb`, i.e. not for the path anything actually runs).
  Suite: 1366 -> 1378. The FD check fails by 2.5e-2 if the Jacobian overwrite
  is dropped, so it can genuinely fail.
- **The `eps_ST` conflation** (§2.5) is fixed and derived rather than
  hardcoded. Effect on validation: CoR 0.1845 -> 0.1843. A correctness fix, not
  a performance one, exactly as predicted.
- **A hierarchy derivation** now states every rung from the exact problem down
  to Lamb, with the assumption, the discarded term, and the undo path at each
  step: `derivations/generalized_newtonian_hierarchy_derivation.jl`.

### The real decision

The measurement in §3 forces it. There is **no cheap-and-accurate rung**:

| option | cost | accuracy |
|:--|:--|:--|
| patch the current A8 (envelope closure) | ~free | fixes ~80% temporal error; leaves wrong basis + leading-order spatial collapse |
| **A7**, `eta = eta(r)` | numerical radial BVP per mode; still diagonal | derived, uses Reid's real eigenfunction; leading-order angular error remains |
| **L4** dense | full `MxM` coupling | the only accurate option |

**Recommendation: price L4 before choosing.** The Gaunt coefficients are pure
geometry and precompute once; only the radial integrals of `eta_l'` change per
step. Assembling `D = sum_l' eta_l' G^l'` at `M=50, L=140` is order 350k
operations per residual -- plausibly affordable, but this is an estimate and
has not been measured. Measuring it is a day's work and it decides the
architecture. If L4 is affordable, build L4. If not, build A7 **and quote its
error** rather than assuming it away.

The envelope closure is still worth shipping as an interim -- it is nearly
free and strictly better than what runs today -- but it should not be mistaken
for progress toward a first-principles model. It is a better patch on the
wrong rung.

### Independently worth doing, whichever rung wins

- **Replace the potential-flow basis with Reid's `U(x)`** (§2.1). 15-34% error
  in the strain field, ~37% in `Oh_eff` per factor-2 in `S`. A validated
  scale-free high-`l` evaluator now exists (BC2 residual 5.4e-15 up to `l=50`).
  Note this breaks `test_carreau_yasuda_exact.jl:12` *for being right*; that
  test pins the wrong field and must be retired, not accommodated.
- **Derive the closure or label it.** Either compute `Q_l = dD_total/dAdot_l`
  properly -- keeping the `mu'(S)` chain-rule term and resolving the
  normalization against Reid's `(l-1)(2l+1)` instead of side-stepping it with a
  ratio -- or keep the current weighting and state it as a postulate with a
  measured error bar. The present middle position is the worst of both.
- **Promote `validate_shear_thinning.jl`** to an asserted, CI-run test with
  per-We-decade tolerances, so §4's failure tail is visible rather than hidden
  behind a median.
- **Housekeeping.** Every `Gamma_l^(a)` number is invalid (§2.1) and is still
  hardcoded into two test files; `section_carreau.tex:158` misattributes the
  3/4 factor; `reid1960_expanded-3.tex:809` states a critical radius of 0.23 mm
  where its own constants give ~23 nm.

## 6. One-line summary

Reid's Newtonian core is sound and gives more to a generalized Newtonian
fluid than the docs currently claim — BC1 and BC2 carry over untouched, BC3
changes by one coefficient, and the momentum equation gains exactly one
term. Everything built on top of it for Carreau-Yasuda is a postulated
closure evaluated on the wrong velocity field at the wrong instant, tested
only against itself, and tuned to a median that hides a 25–30% catastrophic
failure rate.
