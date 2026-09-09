# Incompressible rungs — implementation ledger

Status: **ledger frozen; both incompressible rungs implemented and gated.**

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

Rung 4 (incompressible dfmm), with `Pi` carrying its own dynamics:

```
d_t u   = -(u.grad) u - grad phi - (1/rho_0) div Pi ,   div u = 0
D Pi/Dt = -2 p_0 S0 - [Pi_ik G_jk + Pi_jk G_ik]^dev - Pi/tau_Pi
```

`rho = rho_0` and `p = p_0` are constants of the motion by construction, not
by a numerical constraint, so the density and energy slots of `uold` carry
`rho_0` and `rho_0|u|^2/2 + 3 p_0/2` exactly. That is deliberate: condinit,
the snapshot writer, `rd_cell`, and every dfmm diagnostic keep working with no
special case, and a rung-1 snapshot is directly comparable with a rung-3 one.

Chapman--Enskog consistency check: eliminating `Pi` from rung 4 in the
small-`tau` limit gives `Pi -> -2 p_0 tau S0`, hence
`-(1/rho_0) div Pi -> nu lap u` for `div u = 0`, which is rung 1. So rung 4
reduces to rung 1 as `tau -> 0` at fixed `nu`, and that is a gate.

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

**Timestep.** Advective and viscous, not acoustic -- the substantive gain over
filtering a compressible step:

```
dt = courant * min( dx / max|u| ,  dx^2 / (2 ndim nu) )
```

with the viscous limb applied only in rung 1. Rung 4's stress is hyperbolic
and relaxed by the same exponential map as rung 3, so it needs no parabolic
limit, exactly as in `doc/dfmm_3d.md` Section 4 -- and the timestep is where
that advantage shows up.

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

| Gate | Setup | Criterion | Result |
|---|---|---|---|
| 1 FFT | random field, `N = 8, 16, 32` | round-trip and a direct DFT | round-trip **1.6e-16**; vs direct DFT **2.1e-14**; spectral `d/dx sin(2kx)` **3.3e-14** against a signal of 12.6 |
| 2 Gather | `levelmin = levelmax` | the `ckey` map must hit every lattice site exactly once | asserted at runtime: the gather counts the sites it filled and aborts unless it is `n^3` |
| 3 Projection | random field, `N = 32` | `k . u_hat = 0`, idempotent | spectral divergence **0.776 -> 2.6e-16**; idempotence **2.8e-16** |
| 4 Taylor--Green | `u = (sin kx cos ky, -cos kx sin ky, 0)`, exact NS solution decaying as `exp(-2 nu k^2 t)`; `nu = 0.02`, `t = 0.101` | amplitude *and* shape, since the nonlinear term is a pure gradient that the projector must absorb exactly | standalone `N = 8/16/32`: relative error **2.3e-5 / 5.3e-6 / 6.8e-7**, i.e. the RK2 rate at `dt ~ dx`. In RAMSES at `N = 32`: `max\|u\|` ratio **1.0000007**, `max\|u - u_exact\|` **6.0e-7**, `u_z` **identically zero**, `rho` uniform to the bit, `div u` **1.3e-15** |
| 5 Energy | smooth 3D field, `nu = 0`, 50 steps | inviscid rotational form must conserve energy | `E/E_0 - 1 = ` **2.3e-7**; energy in the truncated band **1.6e-29**, i.e. fully resolved |
| 6 CE limit | rung 4 on Taylor--Green at `tau = 0.02 / 0.002 / 0.0005`, against the rung-1 analytic decay at `nu = p_0 tau/rho_0` | rung 4 must reduce to rung 1 as `tau -> 0` | ratio to the Navier--Stokes answer **1.0295 / 1.00075 / 1.000166** at `max\|Pi\|/p_0 = ` 0.243 / 0.0247 / 0.0062. The 3% deviation at `tau = 0.02` is *the physics*: the evolved stress lags the strain and therefore dissipates less than Newtonian |
| 7 Strain box | rungs 1 and 4 on the blowup initial condition, Family-B `K = 1`, level 5, one deformation time | the indicators must respond, with `div u` held | rung 4: `min lam(p_0 I + Pi)/p_0 = ` **0.433**, `max\|Pi\|/p_0 = ` 0.567, `min g(rank) = ` 0.717, `n(Gamma<0) = 0`, `div u = ` **1.4e-15**. Rung 1 at the same `tau` decays faster (`max\|u\|` 0.835 vs 0.999) |

**The headline of Gate 7.** At `K = 1` the *compressible* moment run
(`doc/dfmm_3d.md` Gate 9) gives `min lam(P)/p = 0.655` at
`max |Pi|/p = 0.584`; the incompressible moment run gives **0.433** at almost
the same anisotropy, `0.567`. So removing compressibility makes the
realizability margin *worse* at fixed stress. That is a physically sensible
reading -- a compressible gas can relieve strain by expanding and can raise
`p` by viscous heating, both of which widen the cone, and an incompressible
one can do neither -- and it is exactly the kind of statement the four-rung
design exists to isolate. It also means the compressible runs were, if
anything, *optimistic* about check 7.

**Two reporting caveats.**

* `econs` is meaningless in these rungs and should be ignored. The internal
  energy is pinned at `3 p_0/2` by construction, so viscous dissipation leaves
  the total-energy budget rather than heating the gas. `mcons` is exactly zero
  and the kinetic energy is the quantity to watch.
* Watch `E_trunc/E`. It is `1.7e-16` for rung 1 on the strain box but
  `2.7e-8` for rung 4, because advecting `Pi` generates finer scales than the
  velocity carries. That is still fully resolved, but it is the number that
  will announce under-resolution first as `K` rises.

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
