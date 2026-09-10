# Incompressible rungs — implementation ledger

Status: **both incompressible rungs implemented and gated, `Pi` and `Q` both
evolved in rung 4.** Two claims in an earlier version of this ledger are
corrected below and marked as such: that `Q` could be dropped (Section 1) and
that the timestep carries no acoustic limb (Section 2). **Rung 4 is
trustworthy for `K <= 0.5` only**; it diverges at `K >= 1` for a reason not
yet established (Section 4c). Rungs 1--3 are unaffected.

Branch: `dfmm_3d_metal`. Scope: `NDIM=3`, single rank, `levelmin == levelmax`,
periodic, uniform grid of `N = 2^L` cells per side. Any other configuration is
rejected at startup, not silently approximated.

---

## 0. Why a fourth rung, and what it is for

The blowup construction is **incompressible**. Running it as a compressible
gas tests several things at once, and one identity makes the confound
unavoidable for a one-mode initial condition (`namelist/dfmm_blowup3d.nml`):

```
Kn * Ma = 0.494 K ,      K = 8 mu / (p Delta) = 8 tau / Delta
```

independent of box size, mode number and `p_0`. At `K = 1` neither Mach nor
Knudsen can be brought below 0.70. So the compressible runs cannot separate
"the Navier--Stokes reduction fails" from "the flow became transonic".

Four rungs separate it, on the *same* grid and the *same* initial condition:

| rung | velocity constraint | stress closure | what it isolates |
|---|---|---|---|
| 1 | `div u = 0` | `mu = p_0 tau`, Newtonian | the reference; the thing predicted to blow up |
| 2 | compressible | `dfmm_closure='ns'` | the cost of compressibility alone |
| 3 | compressible | evolved `Pi`, `Q` | closure + compressibility (implemented, `doc/dfmm_3d.md`) |
| 4 | `div u = 0` | evolved `Pi`, `Q` | the closure alone -- **incompressible dfmm** |

Rung 4 needs `DFMM>=2` to be the incompressible limit of rung 3. At `DFMM=1`
it is the incompressible limit of the ten-moment system -- still a usable
rung, but not the one that closes the 2x2. `incomp_validate` prints which.

Rungs 2 and 3 already exist and differ *only* in the closure selector, which
was the point of implementing `dfmm_closure='ns'` inside the same kernel.
Rungs 1 and 4 differ from each other in exactly the same way, so the 2x2 is
closed and every cell of it shares a discretisation.

**The pressure-gauge objection, and why rung 4 answers it.** Incompressible
Navier--Stokes has no gas pressure: its `p` is a Lagrange multiplier whose
additive gauge is free, which is why the blowup note says its "additive gauge
cannot establish a negative absolute gas pressure". Check 7 -- the decisive
one -- therefore cannot even be *asked* of rung 1. It can be asked of rung 4,
because an incompressible dfmm carries two separate objects:

* `phi`, the multiplier that enforces `div u = 0`. Gauge-free part only.
* `p_0`, the **kinetic** gas pressure, `p_0 = rho_0 theta_0`, which sets
  `mu = p_0 tau` and against which `Pi` is measured.

So `lam_min(p_0 I + Pi) >= 0` is a well-posed, gauge-free statement in the
incompressible setting, and `Pi_NS = -2 p_0 tau S0` is a well-posed
extrapolation to audit. Rung 4 is the only rung in which the note's check 7
can be tested without compressibility in the way.

---

## 1. Governing equations as discretised

Rung 1 (incompressible Navier--Stokes):

```
d_t u = -(u.grad) u - grad phi + nu lap u ,    div u = 0 ,   nu = p_0 tau / rho_0
```

Rung 4 (incompressible dfmm), with `Pi` and `Q` carrying their own dynamics:

```
d_t u      = -(u.grad) u - grad phi - (1/rho_0) div Pi ,   div u = 0

D Pi_ij/Dt = -2 p_0 S_ij - [Pi_ik G_jk + Pi_jk G_ik]^dev
             - d_k Q_ijk + (2/3) delta_ij div q               - Pi_ij/tau_Pi

D Q_ijk/Dt = -theta_0 ( d_k Pi_ij + d_j Pi_ik + d_i Pi_jk )
             -(1/rho_0) sum_l ( Pi_kl d_l Pi_ij + Pi_jl d_l Pi_ik
                                                + Pi_il d_l Pi_jk )
             -[ Q_jkl G_il + Q_ikl G_jl + Q_ijl G_kl ]        - Q_ijk/tau_q
```

with `theta_0 = p_0/rho_0`, `tau_q = tau_Pi/Pr`, and `q_i = Q_ijj/2`. `S0 = S`
exactly, because the deviatoric projection subtracts `(1/3) delta_ij div u`.

**`Q_ijk` is carried, and an earlier version of this ledger was wrong to drop
it.** The argument for dropping it was that `theta = p_0/rho_0` is a constant
of the motion, so `grad theta = 0` and the first-order Chapman--Enskog heat
flux `q_CE = -(5/2) tau_q p grad theta` vanishes identically. That is true,
and it bounds the *equilibrium* value of one contraction of `Q` -- not the
evolved tensor. The `Q` production above is nonzero wherever `Pi` varies in
space, which in this rung is everywhere. Dropping `Q` therefore did not drop a
term that vanishes; it silently replaced the incompressible limit of the
*twenty*-moment system with the incompressible limit of the *ten*-moment one,
which is a different closure and breaks the premise of the 2x2 -- rung 4 is
supposed to differ from rung 3 in the velocity constraint alone.

The `Q` production is the analytically combined `-d_l R_ijkl + T_Q1` of
`doc/dfmm_3d.md` Section 4, with the `d_l rho` term dropped since `rho` is
constant. In the compressible code `R` must stay inside the flux, because it
carries the wave speeds, and `T_Q1` is a separate source, so the cancellation
is only discrete and Gate 5 there measures the residual. Here there are no
Riemann fluxes, so the cancellation is exact by construction.

Contracting the production gives the leading heat flux in this rung:

```
q_i -> -tau_q theta_0 d_j Pi_ij      (= +tau_q tau theta_0 p_0 lap u_i
                                        when Pi = -2 p_0 tau S, div u = 0)
```

a **second-order, Burnett-order** heat flux driven by the velocity Laplacian
rather than by a temperature gradient. So `q` here is not zero; it is
`O(tau^2)` instead of `O(tau)`. Because its first-order value is *exactly*
zero, the measured `|q|` in this rung is by itself the departure from the
first-order closure -- it plays the role that `||q - q_CE||` plays in the
compressible ledger, with `q_CE = 0`. Gate 8 checks it against the closed
form above.

**A ten-moment rung 4 is still selectable, at `DFMM=1`.** It is a legitimate
object, and `incomp_validate` prints which of the two is running, because a
`DFMM=1` incompressible run is *not* comparable with a `DFMM>=2` compressible
one.

`rho = rho_0` and `p = p_0` are constants of the motion by construction, not
by a numerical constraint, so the density and energy slots of `uold` carry
`rho_0` and `rho_0|u|^2/2 + 3 p_0/2` exactly. That is deliberate: condinit,
the snapshot writer, `rd_cell`, and every dfmm diagnostic keep working with no
special case, and a rung-1 snapshot is directly comparable with a rung-3 one.

Chapman--Enskog consistency check: eliminating `Pi` from rung 4 in the
small-`tau` limit gives `Pi -> -2 p_0 tau S0`, hence
`-(1/rho_0) div Pi -> nu lap u` for `div u = 0`, which is rung 1. So rung 4
reduces to rung 1 as `tau -> 0` at fixed `nu`, and that is a gate. Carrying
`Q` does not spoil it -- the `-d_k Q` term is `O(Kn^2)` relative to
`Pi/tau_Pi` -- and measurably *sharpens* it: Gate 6 below converges about
three times faster with `Q` than without.

---

## 2. Discretisation, and why it is spectral

**Spatial operators are spectral**: every derivative, the projection, and the
curl are computed by FFT on the `N = 2^L` grid. Three reasons, in order of
weight:

1. **The projection is exact.** The Leray projector in Fourier space is
   `u_hat <- (I - k k^T/|k|^2) u_hat`, after which `i k . u_hat = 0` to
   machine precision, identically, with no solver tolerance and no iteration.
   This is the property the rung exists to provide. A collocated
   *finite-difference* projection cannot be exact: `div` and `grad` over `2 dx`
   are adjoint to a Laplacian whose stencil decouples the grid into eight
   sub-lattices, so the discrete Laplacian has a checkerboard null space. The
   standard fixes are a staggered (MAC) velocity or an approximate projection;
   spectral avoids the choice entirely.
2. **No numerical viscosity to confound the comparison.** The whole point of
   the four-rung design is that the rungs differ only in the closure. A
   Godunov advection contributes a numerical viscosity of order `c dx / 2`,
   and comparing a run whose dissipation is numerical against one whose
   dissipation is `mu = p_0 tau` would confound exactly the quantity under
   study. This is the same argument that put `dfmm_closure='ns'` inside the
   dfmm kernel instead of writing a separate viscous solver.
3. **No external dependency.** RAMSES levels are always powers of two, so a
   self-contained radix-2 FFT suffices. `bin/Makefile` reaches FFTW only under
   `TURB=1` and at a hardcoded cluster path; a 60-line in-tree transform is
   more portable than that, and it is verified against a direct DFT (Gate 1).

**Nonlinear term in rotational form.** Since
`(u.grad) u = grad(|u|^2/2) - u x omega` with `omega = curl u`, and the
projector annihilates gradients,

```
P[ -(u.grad) u ] = P[ u x omega ] .
```

Computing `u x omega` rather than `u.grad u` means the discrete nonlinear term
conserves kinetic energy exactly in the inviscid, unaliased limit. For a
blowup study that matters: any energy change is then physical or from the
dealiasing truncation, never from the advection scheme.

**Dealiasing: the 2/3 rule.** Modes with `|k_d| > (2/3) k_max` in any
direction are zeroed after the nonlinear product. This is Orszag's standard
truncation and it is the one deliberate departure from "no dissipation" in the
scheme. It is a *projection*, not a dissipation: it is idempotent, it removes
no energy from the retained modes on a resolved field, and its effect is
reported as the fraction of energy in the truncated band so that an
under-resolved run is visible rather than silent.

**Time integration: explicit RK2 (Heun) with a projection per stage.**

```
a(u)   = P[ u x omega ] + P[ stress(u) ]
u_1    = P[ u^n + dt a(u^n) ]
u^n+1  = P[ u^n + (dt/2)( a(u^n) + a(u_1) ) ]
```

Chorin's first-order projection would be acceptable and is more common, but it
carries an O(dt) splitting error in the velocity, and the deliverable of this
study is the *time* at which an indicator crosses a threshold. RK2 costs one
extra projection per step and removes that leading error. Projecting inside
each stage rather than only at the end keeps every intermediate field
divergence-free, so the nonlinear term is never evaluated on a field with
spurious compression.

**Timestep.** Three limbs:

```
dt = courant * min( dx / (ndim (max|u| + c_mom)) ,  dx^2 / (2 ndim nu) )
```

The viscous limb applies only in rung 1: rung 4's stress relaxation is stiff
and handled by an exact exponential map, so it carries no parabolic
constraint, exactly as in `doc/dfmm_3d.md` Section 4.

`c_mom` is a **correction to an earlier claim in this ledger**, which said the
timestep is advective and viscous, "not acoustic -- the substantive gain over
filtering a compressible step". That is true of the ten-moment rung and false
of the twenty-moment one. Once `Q` is evolved, `Pi`'s flux *is* `Q` and `Q`'s
production is `theta_0 grad Pi`, so the pair is hyperbolic: a `z`-directed
wave in `(Pi_xx, Q_xxz)` obeys

```
d_t Pi_xx = -d_z Q_xxz ,      d_t Q_xxz = -theta_0 d_z Pi_xx
```

i.e. it propagates at the thermal speed. So

```
c_mom = sqrt( (3 + sqrt 6) lam_max(P) / rho_0 ) ,    lam_max(P) <= p_0 + 2 max|Pi|
```

using the same `CSCOEF = 3 + sqrt(6)` the compressible kernel uses, and the
cheap upper bound on `lam_max(P)` rather than a per-cell eigenvalue, which
would buy nothing. `c_mom = 0` for rung 1 and for a `DFMM=1` rung 4.

**Removing compressibility does not remove the acoustic timestep** once the
third moment is evolved; it removes it from the velocity equation only.
Measured on Taylor--Green at level 5, this costs a factor 3.4 in step count
(24 steps -> 81). That is a real cost of getting the closure right, and it is
the honest version of the advantage this paragraph used to claim.

---

## 3. Where it sits in the code

The solver operates on `uold` for one level through the oct hash, so it is a
*step*, not a fork of the mesh:

```
r_incomp_step(pst, ilevel)
  gather   uold(2:4)/rho_0            -> u(3, N, N, N)   float64
  [rung 4] gather uold(ipi..ipi+4)    -> Pi(5, N, N, N)
  RK2 with a projection per stage
  scatter  rho_0 u -> uold(2:4),  rho_0 -> uold(1),
           rho_0|u|^2/2 + 3 p_0/2 -> uold(5)
```

Gather/scatter is by `ckey`, which for `levelmin == levelmax` is a bijection
onto the uniform lattice, so the mapping is exact and checkable (Gate 2
verifies that scatter(gather(x)) is the identity to the bit).

Everything is float64 on the host. That is a second reason to want this rung:
it is the float64 reference that `doc/dfmm_3d.md` Section 6 asks for when it
records float32 adequacy for the Stage-4 phase-space sector as *assumed*.

---

## 4. Gates

Reproducers: `namelist/taylorgreen3d.nml` for Gates 4 and 6,
`namelist/incomp_blowup3d.nml` for Gate 7. Gates 1, 3 and 5 are
standalone driver tests of the operators, with no namelist.

Gates 4 and 6 were re-run from the committed namelist after it was written, to
check that the namelist alone reproduces them: Gate 4 gives `max|u|` ratio
**1.000000668**, `max|u - u_exact|` **5.65e-7**, `u_z` identically zero, `rho`
uniform to the bit, and `t = 0.100000000000002`; Gate 6 gives
**1.029621 / 1.000745 / 1.000166** at `max|Pi|/p_0 = ` 0.2435 / 0.0247 /
0.0062, with `div u <= 1.75e-15` and `n(Gamma<0) = 0` throughout.

| Gate | Setup | Criterion | Result |
|---|---|---|---|
| 1 FFT | random field, `N = 8, 16, 32` | round-trip and a direct DFT | round-trip **1.6e-16**; vs direct DFT **2.1e-14**; spectral `d/dx sin(2kx)` **3.3e-14** against a signal of 12.6 |
| 2 Gather | `levelmin = levelmax` | the `ckey` map must hit every lattice site exactly once | asserted at runtime: the gather counts the sites it filled and aborts unless it is `n^3` |
| 3 Projection | random field, `N = 32` | `k . u_hat = 0`, idempotent | spectral divergence **0.776 -> 2.6e-16**; idempotence **2.8e-16** |
| 4 Taylor--Green | `u = (sin kx cos ky, -cos kx sin ky, 0)`, exact NS solution decaying as `exp(-2 nu k^2 t)`; `nu = 0.02`, `t = 0.101` | amplitude *and* shape, since the nonlinear term is a pure gradient that the projector must absorb exactly | standalone `N = 8/16/32`: relative error **2.3e-5 / 5.3e-6 / 6.8e-7**, i.e. the RK2 rate at `dt ~ dx`. In RAMSES at `N = 32`: `max\|u\|` ratio **1.0000007**, `max\|u - u_exact\|` **6.0e-7**, `u_z` **identically zero**, `rho` uniform to the bit, `div u` **1.3e-15** |
| 5 Energy | smooth 3D field, `nu = 0`, 50 steps | inviscid rotational form must conserve energy | `E/E_0 - 1 = ` **2.3e-7**; energy in the truncated band **1.6e-29**, i.e. fully resolved |
| 6 CE limit | rung 4 on Taylor--Green at `tau = 0.02 / 0.002 / 0.0005`, against the rung-1 analytic decay at `nu = p_0 tau/rho_0` | rung 4 must reduce to rung 1 as `tau -> 0` | twenty-moment: ratio to the Navier--Stokes answer **1.03587 / 1.000428 / 1.000054** at `max\|Pi\|/p_0 = ` 0.2142 / 0.0247 / 0.0062. Ten-moment, same runs at `DFMM=1`: 1.029621 / 1.000745 / 1.000166. Carrying `Q` converges **~3x faster** at small `tau` (it supplies part of the missing dissipation) and departs **further** from Navier--Stokes at `tau = 0.02`, where Navier--Stokes is the wrong reference. Both limbs are `O(tau)` to `1`, as required |
| 7 Strain box | rungs 1 and 4 on the blowup initial condition, Family-B `K = 1`, level 5, one deformation time | the indicators must respond, with `div u` held | rung 4 (twenty-moment), minima over the run: `min lam(p_0 I + Pi)/p_0 = ` **0.5482**, `max\|Pi\|/p_0 = ` 0.3633, `min g(rank) = ` **0.6187**, `n(Gamma<0) = 0`, `div u <= ` **2.3e-15**, `E_trunc/E = ` 3.1e-7. Rung 1 at the same `tau` decays (`max\|u\|` 0.837 final) where rung 4 does not (**1.031** final, 1.407 peak) |
| 8 Burnett heat flux | rung 4 on Taylor--Green, `tau = 0.02 / 0.002 / 0.0005`, against the closed form `q_i = -tau_q theta_0 d_j Pi_ij` with `Pi = -2 p_0 tau S` | `q` is *exactly* zero at first order here, so this tests the `Q` production against an analytic second-order value with no leading-order term to hide behind | measured `max\|q\|` = 5.776e-2 / 4.675e-4 / 2.932e-5 against predicted 4.007e-2 / 4.618e-4 / 2.921e-5, i.e. ratios **1.4415 / 1.0122 / 1.0036** -- converging to 1 as it must, with the `tau = 0.02` point 44% high at `max\|Pi\|/p_0 = 0.21` where a second-order prediction should fail. `tau^2` scaling confirmed to **0.8%** between the two small-`tau` points. `max\|Pi\|` itself matches `-2 p_0 tau S` to 0.8 / 0.8 / 0.2% |

**The headline of Gate 7, corrected.** An earlier version of this ledger
compared `min lam(P)/p = 0.655` (compressible, `DFMM=4`) against **0.433**
from the incompressible rung and concluded that removing compressibility costs
0.22 of realizability margin "at almost the same anisotropy". That comparison
was invalid: the incompressible run was ten-moment and the compressible one
twenty-moment, so it was reading a difference in *closure order* as a
difference in *compressibility*. With `Q` carried on both sides, at `K = 1`:

| | compressible (`doc/dfmm_3d.md` Gate 9) | incompressible (Gate 7) |
|---|---|---|
| `min lam(P)/p_0` | 0.655 | **0.548** |
| `max \|Pi\|/p_0` | 0.584 | **0.363** |
| `min g(rank)` | 0.625 | **0.619** |
| `n(Gamma<0)` | 0 | 0 |

Three readings, in order of confidence:

1. **Removing compressibility still costs realizability margin, but half as
   much as reported** -- 0.107 rather than 0.222 -- and it now does so at a
   *lower* stress anisotropy (0.363 vs 0.584) rather than a comparable one. So
   per unit anisotropy the incompressible rung is worse off by more than the
   raw margins suggest, while in absolute terms the gap is smaller. The
   mechanism is unchanged and still physically sensible: a compressible gas
   can relieve strain by expanding and can raise `p` by viscous heating, both
   of which widen the cone; an incompressible one can do neither.
2. **The rank indicator agrees to 1%** between the two rungs, 0.619 vs 0.625.
   These are independent implementations of the dual-frame sector -- Metal
   float32 finite-volume against host float64 spectral -- on the same initial
   condition at the same `K`. It is *not* a strict cross-check, because the
   physics genuinely differs, but it is the strongest consistency signal the
   phase-space sector has, and it did not hold before `Q` was carried (0.717
   vs 0.625).
3. **The flow does not decay.** Rung 1 reaches `max|u| = 0.837` at one
   deformation time; rung 4 reaches **1.031**, having peaked at 1.407. The
   evolved stress lags the strain, so it extracts less energy than the
   Newtonian stress does, and the strain keeps working.

**And the result that only exists because `Q` is carried:** `max |q|/(p_0 c_0)`
reaches **0.732** at `K = 1`, in a flow with **no temperature gradient at
all**. The first-order Chapman--Enskog heat flux is identically zero here, so
this is entirely second-order-and-beyond, and at `K = 1` it is not a
correction -- it is `O(1)`. Any argument that an isothermal incompressible
flow needs no heat-flux sector fails quantitatively at the `K` this study
cares about. This is also the number that most sharply distinguishes rung 4
from rung 1, since rung 1 has no `q` to report.

**A defect Gate 6 exposed, now fixed.** Running `incomp_stress='moment'`
under a `DFMM=0` binary segfaulted on the first step: the moment closure reads
`Pi` out of the hydro state, those slots do not exist, and `incomp_validate`
had no opinion about it. Since the two incompressible rungs are *designed* to
differ in that one namelist entry, the mismatch is the single most likely user
error in the whole four-rung setup, and a SIGSEGV is a bad way to report it.
`incomp_validate` now rejects `incomp_stress='moment'` when `ndfmm < 5` and
names the build flag to change. Verified: the run stops at startup with
`incomp_stress='moment' needs the Pi fields; rebuild with DFMM>=1`.

**Two reporting caveats.**

* `econs` is meaningless in these rungs and should be ignored. The internal
  energy is pinned at `3 p_0/2` by construction, so viscous dissipation leaves
  the total-energy budget rather than heating the gas. `mcons` is exactly zero
  and the kinetic energy is the quantity to watch.
* Watch `E_trunc/E`. It is `2.7e-16` for rung 1 on the strain box but
  `3.1e-7` for rung 4, because advecting `Pi` and `Q` generates finer scales
  than the velocity carries -- an order of magnitude worse than the `2.7e-8`
  the ten-moment rung gave, since `Q`'s production differentiates `Pi`. That
  is still fully resolved, but it is the number that will announce
  under-resolution first as `K` rises, and it will do so sooner now.
* The twenty-moment rung costs a factor **3.4** in step count over the
  ten-moment one at the same `tau` (Taylor--Green level 5: 81 steps vs 24),
  from the `c_mom` limb in Section 2. Budget for it.

---

## 4a. The low-resolution vortex sweep

`namelist/vortex_sweep3d.nml`, level 4 (16^3), all four rungs at six values
of `K`, each run to `t_star = blowup_delta` with `tout_exact` on the quarter
ladder. `t_star` is the construction's singular time in its own clock -- the
note parameterises by the remaining time `Delta = t_star - t`, and the initial
condition freezes the field at `Delta = blowup_delta`. There is **no forcing**,
so the flow does not actually blow up at `t_star`; it is the reference clock,
not an event.

Minimum of `lam_min(P)/p` over the run -- the note's check 7:

| | K=0.1 | K=0.25 | K=0.5 | K=1 | K=2 | K=3 |
|---|---|---|---|---|---|---|
| Navier--Stokes prediction `1-K` | 0.900 | 0.750 | 0.500 | 0.000 | -1.000 | -2.000 |
| rung 1 incomp + NS | *n/a* | *n/a* | *n/a* | *n/a* | *n/a* | *n/a* |
| rung 2 compr + NS | 0.906 | 0.766 | 0.531 | **0.063** | **-0.875** | **-1.812** |
| rung 3 compr + moment | 0.926 | 0.845 | 0.748 | 0.620 | 0.480 | **0.401** |
| rung 4 incomp + moment | 0.920 | 0.829 | 0.714 | *diverges* | *diverges* | *diverges* |

`n(lam<0)` cells rises 0 / 0 / 0 / 0 / **336** / **1504** for rung 2 and is
**0 everywhere** for rung 3 *at this level and this closure order*. Both
qualifiers are essential: Section 4d refines the ladder and rung 3's zero does
**not** survive at twenty moments, while Section 4e shows it does survive, and
converges, at ten.

**The headline of the sweep, as it stands after Sections 4d and 4e.** Rung 2
tracks the Newtonian prediction `1-K` closely and crosses zero between `K = 1`
and `K = 2`, exactly where the note says the extrapolation must fail; that half
is resolution-converged and is the solid result.

Rung 3's zero on this grid went through two corrections. It is **not** a
property of level 4 alone, but neither is it a property of "the moment
closure":

1. At the **twenty**-moment order used for this table, rung 3's zero is a
   level-4 artifact -- refining to levels 5 and 6 makes it diverge, because the
   twenty-moment system is not hyperbolic on the states this flow visits
   (Section 4d).
2. At **ten**-moment order, where the system *is* hyperbolic wherever `P > 0`,
   rung 3 runs clean at every level tested and `min lam(P)/p` **converges** to
   0.39 at `K = 2` and 0.31 at `K = 3` (Section 4e). So the zero is real, once
   the closure is one whose initial-value problem is well posed.

The statement the data supports is therefore neither the original nor the
intermediate retraction, but a sharper one:

> An evolved moment closure removes the negative-variance state that
> Navier--Stokes produces -- but only in a **compressible** gas, and only at
> ten-moment order. Hold the density fixed (rung 4) and the same closure
> crosses at the same `K` as Navier--Stokes; go to twenty moments and the
> system stops being well posed before it can answer.

The mechanism first offered for the zero -- that the nonlinear `-[Pi G]^dev`
term limits the anisotropy the Newtonian extrapolation grows without bound --
is real and is why rung 3's `max|Pi|/p` sits a factor ~3 below rung 2's at
every `K`. But it is not sufficient on its own: Section 4e shows the cone has
to be *growing* too, which is what compression supplies.

Rung 1 is marked *n/a*, not 1.000. It carries no `Pi`, so
`min lam(p_0 I + Pi)/p_0` is identically 1 and the diagnostic is vacuous --
which is the pressure-gauge argument of Section 0 showing up in the output:
check 7 cannot be *asked* of rung 1. The 1.000 the code prints should be read
as "not applicable".

Other indicators at level 4, maxima over the run:

| | K=0.1 | K=0.25 | K=0.5 | K=1 | K=2 | K=3 |
|---|---|---|---|---|---|---|
| `max\|Pi\|/p` rung 2 | 0.115 | 0.288 | 0.578 | 1.160 | 2.322 | 3.488 |
| `max\|Pi\|/p` rung 3 | 0.101 | 0.229 | 0.401 | 0.654 | 0.965 | 1.150 |
| `max\|q\|/(p c_s)` rung 3 | 0.024 | 0.078 | 0.189 | 0.425 | 0.823 | 1.113 |
| `min g(rank)` rung 3 | 0.875 | 0.805 | 0.737 | 0.645 | 0.385 | **0.000** |

So the ordering of the failures at level 4 is: the **rank indicator** `g`
collapses first (rung 3 at `K = 3`), the Newtonian `lam_min` crosses zero next
(rung 2 between `K = 1` and 2), and the evolved `lam_min` does not cross in
this range. Only the first two orderings survive refinement, and the `g`
collapse turns out to be the *early warning* of the rung-3 failure rather than
a separate phenomenon -- at level 5, `K = 3`, `min g` reaches zero at step 41
while `min lam` is still 0.367, and `lam` crosses only at step 157. Level 4 is
adequate for ordering `g` against rung 2's crossing; it is **not** adequate
for the claim that rung 3 does not cross, nor for the
values: `doc/dfmm_3d.md` Section 7 shows `sigma_max(dL/dx)`
and `min g` are resolution-limited extrema, and `|rho/det J - 1|` is 0.31 at
level 4 against 0.043 at level 6, so `Kn_local` from this sweep is a lower
bound.

---

## 4b. Corrected during this work, recorded so the reasoning is not lost

* **`Q_ijk` was dropped on a valid argument about the wrong quantity.**
  `grad theta = 0` in this rung, so the *first-order* Chapman--Enskog heat
  flux is identically zero -- true, and irrelevant to whether the evolved
  tensor is zero. `Q`'s production is `-theta_0 (d_k Pi_ij + ...)`, nonzero
  wherever `Pi` varies. Dropping it turned rung 4 into the incompressible
  limit of the *ten*-moment system, which silently broke the premise of the
  2x2 and, via the invalid Gate-7 comparison above, produced a headline result
  that was half closure order and half compressibility. Lesson: an argument
  that a *closure's equilibrium value* vanishes says nothing about whether the
  *evolved* variable vanishes -- which is the entire premise of a moment
  method, and so exactly the mistake this project should not make.
* **The timestep does carry an acoustic limb.** "Advective and viscous, not
  acoustic -- the substantive gain over filtering a compressible step" was
  true of the ten-moment rung and false once `Q` is evolved: `Pi`'s flux is
  `Q`, `Q`'s production is `theta_0 grad Pi`, and the pair propagates at the
  thermal speed. Cost: 3.4x in step count. Removing compressibility removes
  the acoustic constraint from the *velocity* equation only.
* **`incomp_scatter` left the `Q` slots of `unew` unwritten.** `unew` is not
  initialised from `uold` on this path and `r_set_uold` copies `unew -> uold`
  immediately afterwards, so any dfmm slot the incompressible step does not
  write propagates whatever was in the buffer. It happened to be zero, so this
  was latent rather than active, but it would have become a live bug the
  moment `Q` was evolved. The scatter now copies every slot above 5 from
  `uold` first and then overwrites what it evolved.
* **`incomp_stress='moment'` under a `DFMM=0` binary segfaulted.** See the
  note under Gate 7 above; `incomp_validate` now rejects it at startup.
* **Three spectral operators were differentiating at full spectrum.**
  `incomp_pigrad` and `incomp_qdiv` (new with `Q`) and `incomp_divpi`
  (pre-existing) all fed products without truncating first, where
  `incomp_advect` had always truncated. The comment in `incomp_rhs` claiming
  the 2/3 rule was applied "here and nowhere else" was wrong and is corrected.
  The rule now stated in that comment is the one to keep: **any operator whose
  output feeds a product must truncate its input.**

---

## 4c. Rung 4 diverges at K >= 1, and it is not the aliasing

At level 4 and `K >= 1`, rung 4 grows `|Pi|/p_0` without bound -- 105 by step
~250 at `K = 1`, NaN a few hundred steps later.

**What that growth actually was.** The guard originally stopped the run at
`|Pi|/p_0 > 100`, which fires long after the state is meaningless, and the
`-inf` it then reported was post-mortem noise. With the guard moved to the
physically meaningful boundary (`lam_min(P) < 0`, plus a warning at the
Section 4d hyperbolicity threshold) the same run gives a **measurable**
crossing instead:

| event | `t` | `max\|Pi\|/p_0` | `min lam(P)/p_0` | `min g` |
|---|---|---|---|---|
| hyperbolicity warning | 0.00225 | 0.281 | 0.614 | 0.638 |
| `g` collapses / cone crossing | 0.22415 = `t_star/2` | 0.491 -> 0.523 | 0.018 -> **-0.046** | 0.000 |

So rung 4 at `K = 1` leaves the realizability cone at `|Pi|/p_0 = 0.52` with
`lam_min/p_0 = -0.046`, at exactly half the reference time, and
`E_trunc/E = 1.1e-13` there -- a **fully resolved** crossing, not aliasing.
The `|Pi|/p_0 = 105` was what happened afterwards, integrating a non-hyperbolic
system. The rank indicator `g` reaching zero in the same step as the crossing
(with `n(Gamma<0) = 16`) repeats the ordering seen in rung 3: `g` is the early
warning, `lam_min` is the event.

What it is not:

* **Not the aliasing.** Two real omissions were found and fixed on the way
  here -- `incomp_pigrad`/`incomp_qdiv` and then `incomp_divpi` were
  differentiating at full spectrum instead of truncating first. Those
  produced a genuine grid-scale mode growing 2.9x per step from round-off,
  reaching `E_trunc/E = 0.68`. With all four operators truncating,
  `E_trunc/E` is **1e-13** across the whole sweep and the divergence is in a
  **resolved** mode. The fixes were necessary and are not the cure.
* **Not a CFL limit.** The growth rate per *step* falls with `dt` (2.88 at
  `courant = 0.4`, 1.13 at 0.1) but the rate per unit *time* is
  ~170--310 either way. Quartering `dt` delays it and does not remove it.
* **Not specific to `Q`.** A `DFMM=1` ten-moment rung 4 at `K = 1` shows the
  same growth, ten orders of magnitude in 100 steps. It stays *bounded* there
  only because the flow decays -- which is why one deformation time of Gate 7
  never showed it. So this was latent in the ten-moment rung before `Q`
  landed; `Q` supplies enough extra high-`k` input to make it fatal.

The growth is `Q`-led (`|Q|` reaches 8e15 while `|Pi|` reaches 6e9) and the
feedback is closed: `Q` grows, `Pi` grows, `div Pi` drives `u`, `G` grows, `Q`
grows faster. Two candidate explanations, neither yet established:

1. **A genuine instability of the closure at large strain.** `T_Q2` amplifies
   `Q` at `~3|G|` against BGK relaxation at `1/tau_q`, and
   `3|G| tau_q = 3 (4/Delta)(tau/Pr)` is 1.13 at `K = 0.5` (stable) and 2.25
   at `K = 1` (unstable). The twenty-moment Gaussian closure is hyperbolic
   only in a neighbourhood of equilibrium, so losing it at large anisotropy
   would be a property of the closure, not of this discretisation.
2. **Masked in rung 3 by numerical diffusion.** Rung 3 has the same `T_Q2`
   and is stable to `K = 3`, but it carries `R_ijkl` inside an HLL flux, and a
   compressible gas can also relieve strain by expanding. This spectral rung
   has no numerical diffusion at all -- its design virtue and, here, its
   exposure.

**It gets worse with refinement**, which is the most consequential thing
measured about it. Fitting the exponential growth of `max|Q|` over the second
half of the run at `K = 1`:

| level | `dt` | growth rate | `max E_trunc/E` |
|---|---|---|---|
| 4 (16^3) | 2.0e-3 | **38.4** /time | 1.1e-13 |
| 5 (32^3) | 9.3e-4 | **55.6** /time | 2.8e-13 |

So this is not a low-resolution artefact that refinement will remove -- it
grows. That rules the divergence *out* as a k-independent amplification, which
is what candidate 1 above would predict if `T_Q2` were the mechanism: the
`T_Q2` rate is `3|G| - 1/tau_q`, and for this field
`3|G| tau_q = 18 tau/Delta = 2.25 K` exactly, so the predicted rate at `K = 1`
is `26.8 - 11.9 = 15` /time, resolution-independent. The measured 38--56 /time
is a factor 2.5--4 high and *k*-dependent. **So the first hypothesis is
quantitatively wrong and should be discarded.**

A rate that grows with `k` is the signature of a continuum problem that is
**ill-posed at high wavenumber** -- i.e. loss of hyperbolicity of the
twenty-moment closure, whose flux Jacobian has real eigenvalues only in a
neighbourhood of equilibrium. Imaginary characteristic speeds give a growth
rate proportional to `k`, and a spectral scheme with no numerical diffusion
resolves it as soon as the grid admits it. That would also explain rung 3:
it carries the same closure but damps high `k` with an HLL flux, and a
compressible gas can relieve strain by expanding.

**That last inference was wrong, and Section 4d retracts it.** Rung 3 does not
survive -- it fails at level 5 and 6 for the same reason, and the HLL damping
only postpones the failure to a finer grid. Read "explains why rung 3 survives"
throughout this section as "explains why rung 3 survives *at level 4*".

**The decisive diagnostic has now been run, and the answer is the first
branch: the closure has lost hyperbolicity on the states this flow visits.**
Scanning the *incompressible* twenty-moment principal symbol over strain
shapes puts the loss at

|  | `\|Pi\|/p_0` at which hyperbolicity is lost |
|---|---|
| worst shape found | **0.277** |
| best shape found | **0.408** |
| the sweep's strain shape | 0.340 |
| with `Q != 0` | 0.305 |

against the four rung-4 sweep points reaching `max|Pi|/p_0` = 0.113 / 0.268 /
0.497 / diverge. The first two sit below the threshold and are clean; the third
is above it and the fourth diverges. **So rung 4's divergence is
ill-posedness of the constrained system, not a bug, and no integration fix
will help it.** Section 4d does the same for the compressible member and finds
the corresponding statement for rung 3.

The construction of the symbol had to be corrected first, and the correction
matters enough to record: the original attempt finite-differenced the *flux*
Jacobian alone, exactly as proposed in the paragraph above. That is wrong --
`-2p S0_ij` contains derivatives of `u` and so belongs in the **principal
part**, as do `[Pi G]^dev`, `T_Q1` and `T_Q2`; only the BGK terms are
algebraic. The flux-only version gave the ten-moment fastest speed as the
sound speed 1.291 instead of `sqrt(3)` and declared the Gaussian closure
non-hyperbolic, contradicting Levermore. Rebuilt quasilinearly -- the system is
linear in the directional derivative, so `A[:,c] = -Op(V0, dV = e_c)` needs no
differencing at all -- it reproduces `sqrt(3) = 1.732051` and
`sqrt(3+sqrt(6)) = 2.334414` exactly. Two claims derived from the broken
version were withdrawn: a loss of hyperbolicity at `|Pi|/p_0 = 0.25`, and that
`CSCOEF = 3+sqrt(6)` was 15% too small.

Had the eigenvalues stayed real, the fault would have been in the
integration: the velocity is advanced by RK2 with `Pi` frozen and the moments
then by a single AP step, which is only first-order consistent overall and is
not a stable pairing for a hyperbolic system, so the fix would have been a
two-stage treatment of the coupled `(u, Pi, Q)` system. That work is now
**not** worth doing for its own sake -- it would not make an ill-posed system
well posed -- though it remains the right thing to do if rung 4 is rebuilt at
ten-moment order, where the system *is* hyperbolic and the first-order
splitting is then the leading error.

An earlier framing of this section proposed freezing the velocity as the
diagnostic. The eigenvalue test above is better: it is offline, it needs one
snapshot rather than a modified solver, and it answers the question directly
rather than by elimination.

**Superseded distinguishing test**, kept only to record why it was dropped: run rung 4 at `K = 1`
with the velocity frozen (`Pi` and `Q` advected and sourced but not fed back
into `u`). If `Q` still diverges, it is (1), a closure property, and worth
reporting as such. If it does not, the feedback loop is doing it and the
integration of the coupled `(u, Pi, Q)` system needs the two-stage treatment
the velocity already gets -- currently the velocity is advanced by RK2 with
`Pi` frozen and then the moments by a single AP step, which is only
first-order consistent overall.

Until that is settled, **rung 4 results are trustworthy for `K <= 0.5` only**.
Rung 1 is unaffected and rungs 2--3 do not use this solver at all -- but
rung 3 is **not** thereby sound: Section 4d shows it fails under refinement
for the same underlying reason, the twenty-moment closure's hyperbolic region
being too small for the states this problem visits. Read 4c and 4d together:
they are one obstruction seen in the incompressible and compressible members
of the same closure.

---

## 4d. The resolution ladder, and why rung 3 fails too

Section 4a's headline was a level-4 statement. Repeating rungs 2 and 3 at
`K = 2` and `K = 3` on levels 4, 5, 6 (16^3, 32^3, 64^3), everything else
fixed:

| run | min `lam(P)/p` | max `\|Pi\|/p` | max `n(lam<0)` | steps | first step with `lam<0` |
|---|---|---|---|---|---|
| rung 2, K=2, L4 | -0.875 | 2.322 | 336 | 205 | 1 |
| rung 2, K=2, L5 | -0.968 | 2.421 | 2912 | 822 | 1 |
| rung 2, K=2, L6 | -0.992 | 2.443 | 23264 | 3292 | 1 |
| rung 2, K=3, L4 | -1.812 | 3.488 | 1504 | 232 | 1 |
| rung 2, K=3, L5 | -1.952 | 3.632 | 11936 | 933 | 1 |
| rung 2, K=3, L6 | -1.988 | 3.665 | 95120 | 3742 | 1 |
| rung 3, K=2, L4 | **0.480** | 0.965 | 0 | 117 | never |
| rung 3, K=2, L5 | **0.448** | 1.107 | 0 | 238 | never |
| rung 3, K=2, L6 | *diverges* | 4.5e35 | 11904 | 1084 | 334 |
| rung 3, K=3, L4 | **0.401** | 1.150 | 0 | 103 | never |
| rung 3, K=3, L5 | *diverges* | 1.4e21 | 128 | 213 | 158 |
| rung 3, K=3, L6 | *diverges* | 3.5e35 | 1003 | 263 | 252 |

Two opposite behaviours, and the contrast is the whole point.

**Rung 2 converges.** `min lam(P)/p` tightens monotonically onto the Newtonian
prediction `1 - K`: -0.875 / -0.968 / -0.992 toward -1 at `K = 2`, and
-1.812 / -1.952 / -1.988 toward -2 at `K = 3`. `max|Pi|/p` converges too
(2.322 / 2.421 / 2.443). The violation appears at step 1 at every level.

`n(lam<0)` grows slightly *faster* than the cell count -- 336 / 2912 / 23264
is 8.7x then 8.0x for 8x the cells -- so the violating **volume fraction rises
with refinement** rather than staying fixed: at `K = 3`, 1.17% / 1.81% / 2.17%
of cells at L4 / L5 / L6, with decreasing increments (0.64, 0.36 percentage
points) consistent with converging to a couple of per cent. An earlier version
of this section called it "a fixed fraction of the volume"; that was wrong, and
the direction matters. A grid-scale or shock-capturing artifact confined to a
fixed physical surface would give `n(lam<0) ~ n_side^2`, hence a volume
fraction *shrinking* like `1/n_side`. The measured fraction grows, so the
violation is a genuine volume-filling property of the closure, not a surface
artifact -- which is the conclusion the original phrasing was reaching for, by
a route that did not support it.

Where those cells sit is worth recording too, since it bounds how much of this
is shock physics. Taking the largest relative pressure jump to the six face
neighbours as a shock indicator, the rung-2 violating cells sit at the
96th-100th percentile of that indicator at L5 and L6 (median 98.8) -- i.e.
concentrated in the strongest pressure-gradient region, which is where the
strain is largest and so exactly where Newtonian stress should fail first.
Their absolute jump is only 0.15-0.22 though, far from the order-unity jump of
a strong shock. So the Navier--Stokes closure's negative-variance state
is a converged, resolution-independent, quantitatively predicted property.
**This is the solid result of the study** and it does not depend on anything
in this section.

**Rung 3 fails the other way: refining makes it worse.** Realizable at L4 for
both `K`, and at L5 for `K = 2`; diverges at L5 for `K = 3` and at L6 for
both. A failure that *appears* under refinement at fixed physical setup is the
signature of an ill-posed or marginally-posed system, not of under-resolution
-- under-resolution is cured by refining. Level 4's zero was numerical
diffusion holding the state inside a region the closure cannot actually
sustain.

### The obstruction: positive-definite `P` is not sufficient at twenty moments

Levermore's theorem gives hyperbolicity of the **ten**-moment Gaussian closure
wherever `P > 0`. It says nothing about twenty. Testing the principal symbol
(the validated construction of Section 4c, which reproduces `sqrt(3)` and
`sqrt(3+sqrt(6))`) on the states the code actually visited, read back from the
snapshots:

| run | out | worst `\|Pi\|/p` | non-hyperbolic, of 24 most-strained cells |
|---|---|---|---|
| rung 3, K=3, L4 | 2 | 0.535 | 0/24 |
| rung 3, K=3, L4 | 3 | 0.356 | 0/24 |
| rung 3, K=3, L4 | 4 | 0.203 | 0/24 |
| rung 3, K=3, L5 | 2 | 0.619 | 0/24 |
| rung 3, K=3, L5 | 3 | 0.403 | 0/24 |
| rung 3, K=3, L5 | 4 | 0.613 | **24/24**, max `\|Im\|/\|lam\| = 0.77` |
| rung 3, K=2, L5 | 2..5 | 0.473 | 0/24 |

Cross-tabulating realizability against hyperbolicity per cell on the failing
snapshot (200 most-strained cells of 32768):

| | hyperbolic | non-hyperbolic |
|---|---|---|
| `lam_min(P) >= 0` | 152 | **32** |
| `lam_min(P) < 0` | 0 | 16 |

Every unrealizable cell is non-hyperbolic, as it must be. But **32 cells are
realizable and still non-hyperbolic**, at `lam_min(P)/p = 0.180` and
`|Pi|/p = 0.412` -- comfortably inside the cone. (All 32 report the identical
value: they are copies of one state under the initial condition's symmetry
group.) So the hyperbolic region is a *strict subset* of the realizability
cone.

It is driven by `Q`, not `Pi`. Taking that cell and scaling `Q -> alpha Q`:

| alpha | 0 | 0.1 | 0.2 | 0.5 | 1.0 |
|---|---|---|---|---|---|
| max `\|Im\|/\|lam\|` | 1.2e-16 | 0.364 | 0.498 | 0.625 | 0.689 |

At `alpha = 0` the state is a valid ten-moment state and is hyperbolic to
round-off, exactly as Levermore requires. The mirror test -- `Pi -> beta Pi`
at full `Q` -- stays non-hyperbolic all the way to `beta = 0` (0.398). So
`Q` alone destroys hyperbolicity here and `Pi` alone does not.

Mapping the boundary in the `(|Pi|_inf/p, |Q|)` plane with the failing cell's
shapes, `Q` normalised by its thermal scale `p^{3/2}/rho^{1/2}`:

| `\|Pi\|_inf/p` | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | >=0.6 |
|---|---|---|---|---|---|---|---|
| `lam_min(P)/p` | 1.000 | 0.801 | 0.602 | 0.403 | 0.204 | 0.005 | <0 |
| max hyperbolic `\|Q\|` | 0.896 | 0.653 | 0.425 | 0.233 | 0.084 | 0.000 | 0 |

The hyperbolic region is a bounded neighbourhood of equilibrium in `Q` that
shrinks to nothing exactly as `lam_min(P) -> 0`, and it touches the
realizability boundary only at `Q = 0`, where the ten-moment guarantee takes
over. That is the obstruction, and it is a property of the twenty-moment
closure, not of this discretisation.

**Consequences for the 2x2.** Rung 3 and rung 4 fail for the *same* reason at
different places: the twenty-moment system has no hyperbolic neighbourhood
large enough for the states this problem visits. Section 4c found the
incompressible twenty-moment system loses hyperbolicity above
`|Pi|/p_0 = 0.277-0.408`; rung 3 is the compressible member and loses it in a
`Q`-dependent region inside the cone. Both are salvageable at `DFMM=1`
(ten-moment), which is hyperbolic wherever `P > 0` -- so **the 2x2 should be
closed at ten-moment order**, and the twenty-moment runs reported as what they
are: a demonstration that the closure's well-posed region does not reach the
states of interest.

### What was checked and what is still open

* The exponent in the boundary table is **not** universal. `|Q|_max` looked
  like `0.91 (lam_min(P)/p)^{3/2}` to 2% for the failing cell's shape pair,
  but fitting six random `(Pi, Q)` shape pairs over `lam/p` in [0.15, 0.95]
  gives slopes 0.16, 0.90, 1.17, 0.95, 1.12, 0.16 -- nowhere near 1.5. The
  3/2 law was a coincidence of one shape pair and is retracted. Only the
  qualitative statements (bounded in `Q`, shrinking with `lam_min`, vanishing
  at the cone boundary) are shape-robust.
* The equilibrium twenty-moment symbol's repeated eigenvalues
  (`+-1.732` x2, `+-1.000` x3, `0` x6) are all **semisimple**, geometric
  multiplicity equal to algebraic. So the system is *not* merely weakly
  hyperbolic at equilibrium and there is a genuine hyperbolic neighbourhood --
  the hypothesis that an arbitrarily small `Q` splits a defective pair is
  wrong and is retracted.
* **The in-kernel diagnostic is verified against an independent
  reimplementation.** Reading the snapshots back and recomputing
  `lam_min(p I + Pi)` with a float64 eigensolver reproduces the kernel's
  `n(lam<0)` **exactly** at every time-matched output -- 48/48 (rung 2, K=3,
  L4), 5696/5696 (L6), 0/0 (rung 3, K=3, L4, all outputs), 16/16 (rung 3,
  K=3, L5, out 4) -- and `min lam/p` to 3-4 digits, the residual being the
  float32 kernel against a float64 eigensolve on a near-singular eigenvalue.
* **Still not a true float64 verification.** The Metal build is `NPRE=4` and
  `output_hydro.f90` writes `real(...,kind=4)` regardless of `NPRE`, so both
  the run and the snapshot are single precision; the recheck verifies the
  *diagnostic*, not the *precision*. The `K = 3` rung-2 crossing is far too
  large to be precision-sensitive, but the rung-3 L5 onset -- 16 cells at
  `lam/p = -0.29` -- deserves an `NPRE=8` rerun before publication. That is
  what the CUDA branch is for.

### Reader trap: `rd_cell` returns primitives

`utils/py/ramses.py`'s `rd_cell` returns **primitive** variables. Slot 5 is
`p`, not total energy. Treating it as energy and forming
`eint = E/rho - |u|^2/2` gives, on a `t = 0` snapshot whose initial condition
is exactly uniform `rho = p = 1`, `Pi = 0`, an `eint` ranging over
[-1.80, +1.00] where 1.5 is required -- and then a spurious 752 violating
cells in a run the kernel reports as clean throughout. Always run the `t = 0`
control first: it must give `rho = p = 1` exactly, `|Pi| = 0`, and
`min lam(P)/p = 1`.

The dfmm slots need care in the other direction. `output_hydro.f90:125-137`
deliberately does **not** divide the density-like block (`Pi_ij`, `Q_ijk`) by
`rho`, while it does divide the mass-like tower (`rho L_i`, `rho Sxx`,
`rho Sxv`). So snapshot slots 6..20 hold `Pi` and `Q` themselves, and slots
21..38 hold `L`, `Sxx`, `Sxv` per unit mass.

---

## 4e. Closing the 2x2 at ten-moment order, and what it actually shows

Section 4d ended with the recommendation to close the 2x2 at `DFMM=1`, since
the ten-moment system is hyperbolic wherever `P > 0` and the twenty-moment one
is not. Done: `DFMM=1 INIT=BLOWUP` (NVAR = 10, NDFMM = 5), same namelists,
same `K` calibration, levels 4-6.

### Rung 3, ten-moment: converges, and never leaves the cone

| run | min `lam(P)/p` | max `\|Pi\|/p` | `n(lam<0)` | completed |
|---|---|---|---|---|
| K=2, L4 | 0.4400 | 1.166 | 0 | yes |
| K=2, L5 | 0.4039 | 1.337 | 0 | yes |
| K=2, L6 | **0.3935** | 1.388 | 0 | yes |
| K=3, L4 | 0.3530 | 1.379 | 0 | yes |
| K=3, L5 | 0.3168 | 1.570 | 0 | yes |
| K=3, L6 | **0.3076** | 1.623 | 0 | yes |

Compare the twenty-moment column of Section 4d, which diverged at L5 for
`K = 3` and at L6 for both. At ten moments every run completes, `n(lam<0)` is
zero throughout, and `min lam(P)/p` **converges** -- increments of -0.036 then
-0.010 at `K = 2`, and -0.036 then -0.009 at `K = 3`. `max|Pi|/p` converges
too, and note it reaches **1.62**, further from equilibrium than the
twenty-moment runs ever got (1.15) before failing. That is Levermore's theorem
doing its job: with only ten moments, staying inside the cone is sufficient for
hyperbolicity, so the system can be driven far from equilibrium without the
integration losing meaning.

**This is the resolution-converged form of the claim Section 4a originally
made, and at ten-moment order it holds.** It also settles a question 4d could
not: the twenty-moment failure is specifically the `Q` sector, not a defect in
the shared transport or source machinery, because the identical setup with `Q`
removed runs clean at every level tested.

### Rung 4, ten-moment: crosses between K = 1 and K = 2

| run | min `lam(P)/p_0` | max `\|Pi\|/p_0` | outcome |
|---|---|---|---|
| K=1, L4 | 0.4648 | 1.310 | survives to `t_star` |
| K=1, L5 | 0.4323 | 1.390 | survives to `t_star` |
| K=2, L4 | **-0.0064** | 1.006 | **crosses** |
| K=2, L5 | **-0.0047** | 1.004 | **crosses** |
| K=3, L4 | **-0.0044** | 1.003 | **crosses** |

The crossings all occur at `max|Pi|/p_0 ~ 1.005`, which is structural rather
than coincidental: in the incompressible rung `P = p_0 I + Pi` with `p_0`
constant, so `lam_min(P)/p_0 = 1 + lam_min(Pi)/p_0` and the cone boundary sits
at `|Pi|/p_0 = 1` exactly. The new guard (Section 4c) catches these at
`lam/p_0 = -0.006` instead of letting them run to `|Pi|/p_0 = 105`, so the
crossing is now a measurement.

### The 2x2 finally pays off: both ingredients are load-bearing

With all four cells well posed:

| rung | closure | crosses check 7? |
|---|---|---|
| 1 incompressible + Newtonian | *not askable* -- carries no `Pi` (Section 0) |
| 2 compressible + Newtonian | **yes**, between `K = 1` and 2, converging to `1-K` |
| 3 compressible + moments (10) | **no**, converged, `min lam/p = 0.31` at `K = 3` |
| 4 incompressible + moments (10) | **yes**, between `K = 1` and 2 |

Read along the rows and columns:

* changing **only the closure** (rung 2 -> rung 3) removes the violation;
* changing **only the compressibility** (rung 3 -> rung 4) puts it back.

So neither ingredient alone is the answer, and the original framing -- "the
moment closure does not produce the negative-variance state" -- was picking out
one of two necessary conditions. The statement the data supports is:

> The evolved moment closure removes the negative-variance state **only in a
> compressible gas**. Compression raises `p`, which enlarges the realizability
> cone faster than the evolved `Pi` grows into it. Hold the density fixed and
> the cone stops growing, and the same closure crosses at the same `K` as
> Navier--Stokes.

**The rung-3/rung-4 contrast is Mach-matched, which is what makes it a
controlled comparison.** At `K = 2` the compressible rung attains a maximum
local Mach number of 1.4970 and the incompressible rung at `incomp_p0 = 1`
attains 1.502 -- 0.3% apart, both close to the construction's `Ma = 1.1 sqrt(K)
= 1.556`. Same `K`, same Mach number, same closure order, same initial
condition; the only difference is whether the gas may compress, and that alone
flips check 7.

### Caveat that must travel with any rung-4 number: it depends on `incomp_p0`

Rung 4 has a free parameter the compressible rungs do not. Scanning it at
`K = 2`, L4, with the velocity field unchanged:

| `incomp_p0` | `c_s` | max `\|u\|` | `Ma` | outcome |
|---|---|---|---|---|
| 0.25 | 0.6455 | 1.944 | 3.01 | crosses, `lam/p_0 = -0.058` |
| 1.0 | 1.2910 | 1.939 | 1.50 | crosses, `lam/p_0 = -0.006` |
| 4.0 | 2.5820 | 1.934 | 0.75 | **survives to `t_star`** |

`max|u|` is the same to 0.5% in all three, so this is not the velocity field
changing -- it is the Mach number. The mechanism: `Pi` relaxes toward
`-2 p_0 tau S0`, so `Pi ~ p_0`, and the momentum equation feels it through
`-(1/rho_0) d_j Pi_ij`, a force that therefore also scales with `p_0`. The
back-reaction of stress on velocity strengthens with `p_0`, damps the strain,
and limits `Pi`. So `|Pi|/p_0` is **not** invariant under rescaling `p_0`, and

> a rung-4 answer to check 7 is a statement at one Mach number, not a
> gauge-free statement.

This matters because rung 4 exists precisely to ask check 7 of an
incompressible flow, which Section 0's pressure-gauge argument says cannot be
asked of incompressible Navier--Stokes. Rung 4 does make the question
*well posed* -- `p_0` is a genuine thermodynamic pressure here, not a Lagrange
multiplier -- but it does not make it *parameter-free*. Any published rung-4
result must state `incomp_p0` and, better, be Mach-matched to the compressible
rung it is being compared with, as the `K = 2` comparison above is.

---

## 5. References

* `doc/dfmm_3d.md` — the compressible rungs, the field ledger the gather
  reads, and the source terms rung 4 reuses.
* `~/Downloads/before_blowup_ideal_gas_pedagogical.pdf` — the five checks; its
  Section on the pressure multiplier is the argument in Section 0 above.
* Chorin (1968), Temam (1969) — the projection method.
* Orszag (1971) — the 2/3 dealiasing rule.
* Almgren, Bell & Szymczak (1996) — approximate projections, and why a
  collocated exact projection is not available.
