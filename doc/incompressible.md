# Incompressible rungs — implementation ledger

Status: **both incompressible rungs implemented and gated, `Pi` and `Q` both
evolved in rung 4.** Two claims in an earlier version of this ledger are
corrected below and marked as such: that `Q` could be dropped (Section 1) and
that the timestep carries no acoustic limb (Section 2).

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
