# dfmm in 3D on the Metal backend — implementation ledger

Status: **ledger frozen; Stages 1-4 implemented and gated -- the dual-frame
method is complete.** This document fixes the
equations, field layout, flux ledger, source ledger, realizability rule and
diagnostics *before* physics code is written, following the practice used for
the mini-RAMSES multimoment work. Any change to the equations below is a
change to this file first.

Branch: `dfmm_3d_metal`. Backend: Apple Metal (`COMPILER=METAL`), `NDIM=3` only.

---

## 0. Objective and what "dfmm in 3D" means here

Target science question: run the OpenAI finite-time-blowup construction for
incompressible Navier--Stokes as a *gas* problem and measure which of the
assumptions behind the Navier--Stokes reduction fails first, where, and by how
much. The pedagogical note `before_blowup_ideal_gas_pedagogical.pdf` lists five
separable checks; three of them are directly instrumentable by a moment scheme:

| Note's check | What must be evolved to see it | dfmm diagnostic |
|---|---|---|
| 1. Spatial locality (`Kn` not small) | local mean free path vs. gradient scale | `Kn_local` |
| 2. Stress cannot follow strain (`t_c/t_def` not small) | full pressure tensor with its own relaxation | `dev_NS = ‖Pi - Pi_NS‖ / p` |
| 4. Viscous heating needs an energy equation | compressible energy equation | already in RAMSES |
| 7. **Pressure cannot have negative variance** | full pressure tensor, positivity monitored | `lam_min(P) / p` |
| (paper's indicator 1) | evolved heat flux vs. Fourier prediction | `‖q - q_CE‖ / (p c_s)` |
| (paper's indicator 2) | phase-space rank collapse | per-axis `gamma` |

Check 7 is the decisive one in the note: the Newtonian extrapolation gives
`P_NS(0,t) = diag(p+4mu/Delta, p+4mu/Delta, p-8mu/Delta)`, which predicts a
*negative* axial pressure once `Delta < 8 mu / p`. That is kinematically
impossible for any non-negative distribution function. A ten-moment scheme
carries `P_ij` as a primary variable and therefore reports this failure as a
measured quantity instead of an extrapolation.

**Which dfmm.** The reference implementations are two different schemes with the
same name:

* `~/dfmm/py-1d` — 1D **Eulerian finite-volume** eight-field moment scheme
  (HLL fluxes, operator-split BGK). This is what the moment-scheme paper
  documents, and the paper explicitly states that the Lagrangian coordinate
  "can be added to any finite-volume scheme as an additional conserved field,
  with fluxes computed from the same Riemann solver used for the gas
  variables."
* `~/dfmm` top level (Julia) — a **variational Lagrangian mass-coordinate**
  discretisation on hierarchical grids, with a 3D extension (M3-7) and Berry
  connection. Its discretisation (implicit Hamilton--Pontryagin Newton solves
  on a Lagrangian mesh) has no correspondence to RAMSES's Eulerian AMR Godunov
  architecture.

This port is the **3D generalisation of the Eulerian finite-volume scheme**.
Where the Julia work is used, it is used as a cross-check on the 3D tensor
algebra (`src/berry.jl`, `src/cholesky_DD_3d.jl`), not as an architecture.

---

## 1. Kinetic starting point

BGK-relaxed kinetic equation for one species, with `c = v - u`:

```
d_t f + v.grad_x f + g.grad_v f = -(f - f_M)/tau
```

Central moments, all mass-weighted:

```
rho     = m Int f d3v
rho u_i = m Int v_i f d3v
P_ij    = m Int c_i c_j f d3v            (6 components, symmetric)
Q_ijk   = m Int c_i c_j c_k f d3v        (10 components, symmetric)
R_ijkl  = m Int c_i c_j c_k c_l f d3v    (closed, not evolved)
p       = tr P / 3 ,   Pi_ij = P_ij - p delta_ij   (5 independent, traceless)
q_i     = Q_ijj / 2                      (contracted heat flux vector)
theta   = p / rho                        (temperature in energy units)
```

Fourth-moment closure is **Wick/Gaussian**, the exact 3D generalisation of the
1D reference's `M4 = rho u^4 + 6 u^2 P + 4 u Q + 3 P^2/rho`:

```
R_ijkl = ( P_ij P_kl + P_ik P_jl + P_il P_jk ) / rho
```

The maximum-entropy (polynomial-exponent) closure of the paper's Section 6 is
**not** ported. It needs a per-cell nonlinear Newton solve on 24-point
Gauss--Hermite quadrature, the paper itself reports its practical gain over
Wick as modest, and its `kappa in [1,5]` / `|s| <= 1.4` box caps admit states
that no distribution realises (the cone requires `kappa >= 1 + s^2`, so the
corner `(1.4, 1)` is unrealisable). It is unsuitable for a GPU kernel and
would import a repair we do not want.

### Central-moment transport theorem

Every source term below follows from one identity. For any polynomial
`phi(c)`, using `Du_m/Dt = g_m - (1/rho) d_l P_lm`:

```
d_t(rho<phi>) + d_l(rho u_l <phi>) + d_l(rho <c_l phi>)
    = <d phi / d c_m> d_l P_lm  -  rho <c_l d phi / d c_m> d_l u_m  +  collisions
```

* `phi = c_i c_j` gives the `P_ij` equation.
* `phi = c_i c_j c_k` gives the `Q_ijk` equation.

Both are written out in Section 3.

---

## 2. Field ledger

RAMSES conserved index 1..5 keeps its **exact** standard meaning, so cooling,
gravity synchronisation, AMR restriction/prolongation, refinement flagging,
`upload_kernel`'s internal-energy path and the output writer all continue to
work unmodified. Everything dfmm adds sits at index 6 and above.

| ivar | symbol | kind | stage |
|---|---|---|---|
| 1 | `rho` | conserved | base |
| 2..4 | `rho u_i` | conserved | base |
| 5 | `E = rho|u|^2/2 + 3p/2` | conserved (`gamma = 5/3`) | base |
| 6..10 | `Pi_xx, Pi_yy, Pi_xy, Pi_xz, Pi_yz` | density-like, `Pi_zz = -(Pi_xx+Pi_yy)` | **1** |
| 11..20 | `Q_ijk` (10 comps) | density-like | **2** |
| 21..23 | `rho D_i`, `D_i = L_i - x_i` | mass-like | **3** |
| 24..29 | `rho Sxx_ij` (6 comps) | mass-like | **4** |
| 30..38 | `rho Sxv_ij` (9 comps, not symmetric) | mass-like | **4** |

`NVAR` = 10 / 20 / 23 / 38 at the four stages, selected by `DFMM=n` in
`bin/Makefile` (`NDFMM` = 5 / 15 / 18 / 33). Every stage stays selectable so
its gate results below remain reproducible, and because each sector is
expensive enough to be worth switching off when it is not needed.

**Why the displacement and not the label.** Stage 3 stores `D_i = L_i - x_i`,
not `L_i`. On a periodic box a periodic flow satisfies
`L(x + Lbox e) = L(x) + Lbox e`, so `L` carries a jump of one box length across
the wrap plane. A centred difference of `L` there reports a deformation tensor
too large by `Lbox/(2 dx)` -- the whole grid -- and an upwind advection of that
jump smears it, corrupting a band of cells permanently. `D` is periodic and
smooth, starts at zero, and stays small, so

```
d L_i / d x_j = delta_ij + d D_i / d x_j
```

is clean everywhere and float32 carries `D` at its own magnitude. The price is
one source term, `D D_i / Dt = -u_i`, which follows from `D L_i / Dt = 0`.

**Why `Sxv` is not symmetrised.** Its antisymmetric part is the phase-space
packet's angular momentum, which any flow with vorticity generates from an
isotropic start. Storing 6 instead of 9 would discard physics, not storage.
The Stage-4 gate below measures the antisymmetric part directly, because that
is the part the `-Sxv G^T` term generates.

`Q_ijk` component order is
`xxx, yyy, zzz, xxy, xxz, yyx, yyz, zzx, zzy, xyz`. The kernel carries the
symmetric-index map explicitly (`DF_QIJK`, `DF_QMAP`), so every tensor
contraction is written once as a loop over the ten packed slots rather than
component by component.

The contracted heat flux is `q_i = Q_ijj / 2`, i.e.

```
q_x = (Q_xxx + Q_yyx + Q_zzx)/2
q_y = (Q_xxy + Q_yyy + Q_zzy)/2
q_z = (Q_xxz + Q_yyz + Q_zzz)/2
```

"Density-like" means the transport flux is `u_k X` and AMR restriction is a
plain volume average (RAMSES's existing conservative average is already
correct). "Mass-like" means the stored variable is `rho` times a specific
quantity and the flux is `u_k (rho X)`, i.e. RAMSES's passive-scalar
convention.

### Why `Pi` and not the full second raw moment

The alternative is to store the six components of `E_ij = rho u_i u_j + P_ij`
and derive `E_tot = tr E / 2`. Rejected, on two counts:

1. `E_ij` is **not** a conserved quantity once collisions act — BGK relaxes
   `Pi` — so exact conservation of the anisotropic part buys nothing. Only
   `rho`, `rho u_i` and `E_tot` are physically conserved, and those stay in
   exact conservation form here.
2. Recovering `Pi_ij = E_ij - rho u_i u_j - p delta_ij` is a catastrophic
   cancellation at high Mach number. The blowup construction drives
   `|u| -> large` by design, and the Metal kernels are float32. Carrying `Pi`
   directly keeps the anisotropy at its own magnitude.

This matches the mini-RAMSES multimoment work, which evolves central variables
with explicit "conversion sources" for the same reason.

---

## 3. Flux ledger

All fluxes are the `k`-th component (flux in direction `x_k`).

**Mass** (unchanged):

```
F[rho]_k = rho u_k
```

**Momentum** (RAMSES flux plus one term):

```
F[rho u_i]_k = rho u_i u_k + p delta_ik + Pi_ik
```

**Total energy** (RAMSES flux plus two terms). This is exactly `tr/2` of the
second-moment flux, which is why index 5 stays compatible:

```
F[E]_k = u_k (E + p) + u_i Pi_ik + q_k
```

**Anisotropic pressure**:

```
F[Pi_ij]_k = u_k Pi_ij + Q_ijk - (2/3) delta_ij q_k
```

Contracting on `ij` gives `Q_iik - 2 q_k = 0`, so the flux of the traceless
block is itself traceless and storing only five components stays consistent:
the implicit `F[Pi_zz]_k = -(F[Pi_xx]_k + F[Pi_yy]_k)` is exact, and because
HLL is linear in the fluxes and states it preserves that relation.

**Third central moment**:

```
F[Q_ijk]_l = u_l Q_ijk + R_ijkl        (R by Wick, Section 1)
```

**Passive tower** (Stages 3--4), one flux form for all of them:

```
F[rho X]_k = rho X u_k          X in { L_i, Sxx_ij, Sxv_ij }
```

Checks:

* Setting `Pi = 0`, `Q = 0` reduces the mass/momentum/energy fluxes to
  RAMSES's Euler fluxes identically. This is a regression gate, not a claim.
* Contracting `F[Pi_ij]_k` on `ij` gives zero, as it must for a traceless
  field.
* Reducing to one dimension (`i=j=k=l=x`, no transverse structure) reproduces
  the reference `_common.py:hll_edge_flux` term by term, including
  `R_xxxx = 3 P_xx^2 / rho`.

---

## 4. Source ledger

Let `G_ij = d_j u_i` be the velocity gradient, `S_ij = (G_ij + G_ji)/2` the
strain rate and `S0_ij = S_ij - (1/3) delta_ij G_kk` its deviatoric part.

**Anisotropic pressure** (Stage 1):

```
S[Pi_ij] = -2 p S0_ij
           - [ Pi_ik G_jk + Pi_jk G_ik ]^dev
           - Pi_ij / tau_Pi
```

The first term is the production of anisotropy by strain, the second its
advective distortion, the third BGK relaxation. Quasi-steady balance of the
first and third terms gives `Pi -> -2 p tau_Pi S0`, i.e. Newtonian viscous
stress with `mu = p tau_Pi`. This is precisely the extrapolation the blowup
note audits, so `Pi_NS = -2 p tau_Pi S0` is computed alongside `Pi` as the
Stage-1 closure diagnostic.

**Third central moment**:

```
S[Q_ijk] = (1/rho) ( P_jk d_l P_li + P_ik d_l P_lj + P_ij d_l P_lk )    (= T_Q1)
           - [ Q_jkl G_il + Q_ikl G_jl + Q_ijl G_kl ]                   (= T_Q2)
           - Q_ijk / tau_q
```

Equilibrium check (`Pi = 0`, `Q = 0`): the first source term and `-d_l R_ijkl`
do **not** cancel; their sum is
`-(delta_ij d_k + delta_ik d_j + delta_jk d_i) (p^2/rho) d ln theta`, i.e. a
temperature gradient generates heat flux, as it must. Contracting and
balancing against `-q/tau_q` gives

```
q_i -> -(5/2) tau_q p d_i theta
```

which is the BGK Chapman--Enskog heat flux with `kappa = (5/2) tau p (k_B/m)`
and hence `Pr = 1`.

`T_Q1` and `-d_l R_ijkl` cancel *analytically* down to that temperature
gradient, and the cancellation is worth writing out because it determines where
each term may be evaluated. Since `d_l R_ijkl` contains
`(P_ij d_l P_kl + P_ik d_l P_jl + P_jk d_l P_il)/rho` — exactly `T_Q1` — the
divergence-of-`P` parts cancel identically, leaving

```
-d_l R_ijkl + T_Q1 =
    -(1/rho)   sum_l ( d_l P_ij P_kl + d_l P_ik P_jl + d_l P_jk P_il )
    +(1/rho^2) sum_l ( P_ij P_kl + P_ik P_jl + P_il P_jk ) d_l rho
```

which for `P = p I` collapses to
`-(delta_ij p d_k theta + delta_ik p d_j theta + delta_jk p d_i theta)`. This
combined form is what the MUSCL predictor uses, so no O(1) term is left
unbalanced in the reconstruction. In the update itself the two are *not*
combined: `R` stays in the flux, because it carries the characteristic
structure that sets the wave speed, and `T_Q1` is a source. Their discrete
difference is then a truncation error that converges under refinement, which
Gate 5 measures directly.

Two relaxation times are therefore exposed:

```
tau_Pi = tau            (mu = p tau)
tau_q  = tau / Pr_target
```

Setting `Pr_target = 2/3` reproduces the hard-sphere Prandtl number that the
blowup note's Eq. (16) uses, so the note's coefficients
`mu/p = (5/4) t_c`, `kappa/(rho c_p) = (3/2) nu` can be matched directly. This
means `tau = (5/4) t_c` when calibrating against the note.

**Phase-space sector** (Stage 4) -- the second frame. With the 6x6 phase-space
covariance

```
M = [[Sxx, Sxv], [Sxv^T, Svv]] ,    Svv = P/rho ,
```

the Liouville evolution of the local packet under the linearised flow
`Jac = [[0, I], [0, -G]]` is `D M/Dt = Jac M + M Jac^T`, giving

```
D Sxx / Dt = Sxv + Sxv^T
D Sxv / Dt = Svv - Sxv G^T
D Svv / Dt = -G Svv - Svv G^T      (already carried by p and Pi)
```

1D reduction: `Sxx' = 2 Sxv`, `Sxv' = Svv - Sxv du/dx`, `Svv' = -2 Svv du/dx`,
matching Eqs. (9)--(11) of the paper.

**Why the covariance and not the Cholesky factors.** The 1D reference
(`py-1d/dfmm/schemes/cholesky.py`) carries the factors `alpha, beta` of
`Sigma = L L^T` with `L = [[alpha, 0], [beta, gamma]]`, so that

```
Sxx = alpha^2 ,  Sxv = alpha beta ,  Svv = beta^2 + gamma^2 ,
D alpha / Dt = beta ,  D beta / Dt = gamma^2/alpha - (du/dx) beta .
```

Substituting shows the two formulations are *analytically identical*:

```
Sxx' = 2 alpha alpha' = 2 alpha beta                    = 2 Sxv
Sxv' = alpha' beta + alpha beta' = beta^2 + gamma^2 - G alpha beta
                                                        = Svv - G Sxv
```

so carrying the covariance directly gives up nothing, and gains three things:
the system is linear in the evolved variables, there is no `1/alpha` needing a
floor, and realizability needs no clip anywhere because `gamma` is no longer a
state variable. `gamma` becomes a pure diagnostic -- Section 5.

**Relaxation.** Collisions randomise velocity, so they destroy the
position-velocity correlation on the collision time; they do not move
particles, so `Sxx` is untouched. `Sxv` therefore relaxes with `tau_Pi` under
the same asymptotic-preserving map as `Pi` and `Q`, and the choice of map
matters: it gives the collisional equilibrium

```
Sxv -> tau theta   =>   D Sxx / Dt = 2 tau theta ,
```

i.e. Brownian spreading with diffusivity `tau theta = nu`. The reference
applies its exponential decay *after* the explicit source instead, which sends
`Sxv -> 0` and freezes `Sxx`, losing the diffusive limit entirely. The two
agree in the collisionless regime the reference's own test problems use.
Gate 8 measures the diffusive limit directly and recovers `D = 0.0098` against
`nu = tau theta = 0.01`.

**Time centring is a correctness requirement, not a refinement.** `Sxv` is
advanced first and `Sxx` then uses the *trapezoidal* average of `Sxv` across
the step. In a uniform flow the exact solution is `Sxv = theta t`,
`Sxx = sigma_x0^2 + theta t^2`, so the Schur complement of Section 5 has
`Gamma/Svv = sigma_x0^2/(sigma_x0^2 + theta t^2) > 0` for all time. Forward
Euler on `Sxx` gives `Sxx = sigma_x0^2 + theta(t^2 - n dt^2)` instead -- a
one-signed lag growing linearly in the step count -- so

```
Gamma/Svv = 1 - t^2 / (sigma_x0^2 + t^2 - n dt^2)
```

crosses zero once `n dt^2 = sigma_x0^2`. That was measured at step 11 of the
Stage-4 uniform-flow gate at level 4 before the trapezoid was used, and the
first six steps matched the formula above to five digits. It is a purely
numerical rank collapse and it is indistinguishable in a real run from the
physical one the diagnostic exists to detect. With the trapezoid the
telescoping sum gives `Sxx_N = sigma_x0^2 + theta t^2` exactly. This is the
Stormer--Verlet pairing appropriate to a position-velocity pair, and it is
what replaces the reference's structural guarantee.

**No Berry connection is needed.** The Julia work
(`src/berry.jl`, `src/cholesky_DD_3d.jl`) parameterises the Cholesky factor as
`R(theta) diag(alpha)` and carries a Berry 1-form for the gauge freedom of the
principal-axis frame. That gauge exists only because the factorisation is not
unique. The covariance `M` is the physical object and has no gauge freedom, so
in this formulation there is nothing to connect. The Julia files were used as
a cross-check on the 3D tensor algebra, as Section 0 says, not as an
architecture.

**Lagrangian displacement** (Stage 3): `D D_i / Dt = -u_i`, initialised to
zero (Section 2). Applied unconditionally rather than under `dfmm_source`:
`dfmm_source` selects whether the *closure* is driven -- the strain production
of `Pi`, the production of `Q` -- which is what the Euler-reduction gate
switches off. The second frame's evolution is not a closure model but the
definition of the frame; it has no modelling content to disable, and the tower
has no feedback whatsoever on the hyperbolic core, so leaving it on cannot
perturb that gate. The same reasoning applies to the Liouville terms above,
and it makes `dfmm_source=.false.` a clean manufactured solution: with `Pi`
and `Q` pinned at zero, `Svv = theta I` exactly and every component of `Sxx`
and `Sxv` has a closed form. That is Gate 7.

**Operator order per step**, fixed:

1. Unsplit Godunov transport of the full state (Section 3 fluxes).
2. Gradient sources from the old-time stencil: `S[Pi]` strain terms,
   `S[Q]` terms, phase-space Liouville terms.
3. Exact-exponential BGK relaxation of `Pi` and `Q` (see Section 5).
4. Realizability audit and diagnostics.

---

## 5. Realizability

The single physical constraint is `P_ij` positive semidefinite, since
`a^T P a = m Int (a.c)^2 f d3v >= 0` for every direction `a`. Equivalently
`lam_min(p I + Pi) >= 0`. This is check 7 of the blowup note, so it is
**instrumented and reported, never silently repaired**.

Three separate mechanisms, in order of preference:

1. **Asymptotic-preserving relaxation map.** The production terms and the
   BGK sink are integrated *together*. For `dX/dt = rate - X/tau` with `rate`
   frozen over the step the exact solution is

   ```
   X <- X_old d + tau (1 - d) rate,      d = exp(-dt/tau)
   ```

   and `rate` must be the *whole* non-stiff right-hand side: the production
   terms `T` **plus the transport rate the Godunov step already applied**,
   which is recoverable as `(unew - uold)/dt` because `set_unew` copies
   `uold -> unew` before transport.

   This is unconditionally stable and asymptotic-preserving in both terms:
   `dt << tau` gives `X_transported + dt T`, and `dt >> tau` gives
   `X -> tau (transport + T)`. Applying the source explicitly and *then*
   multiplying by `d` — the naive splitting — instead sends `X -> 0` whenever
   `dt >> tau`, which would destroy precisely the limit this study has to
   measure deviations from. `tau <= 0` means collisionless: `d = 1` and
   `tau(1-d) -> dt`, so the map degrades to `X_transported + dt T`, the
   correct explicit update.

   Including the transport rate is **not optional for `Q`**. Its flux carries
   the Wick fourth moment `R_ijkl`, which is O(1) and cancels all but a
   temperature gradient against `T_Q1` (Section 4); relaxing toward `T_Q1`
   alone whenever `dt >> tau_q` would replace the Fourier heat flux
   `-(5/2) tau_q p grad theta` by the unrelated quantity `theta grad p`. For
   `Pi` the same omission is only O(tau^2) — a Burnett-order correction, and
   indeed Gate 3 is unchanged to the digits printed by including it — but the
   two blocks are treated identically for uniformity.

   For the same reason the stiff production terms are **not** folded into the
   MUSCL predictor — `-2 p S0 - [Pi G]^dev` for `Pi`, and `T_Q2` for `Q`. An
   unrelaxed half-step copy of a stiff term overshoots the interface states
   when `dt >> tau`. The 1D reference makes the same choice: sources follow the
   flux update. Everything non-stiff *is* retained, so transport stays second
   order: `-(1/rho) d_k Pi_ik`, `-(2/3) Pi_kl d_l u_k`, `-(2/3) div q`,
   `-d_k Q_ijk + (2/3) delta_ij div q`, and the combined
   `-d_l R_ijkl + T_Q1` of Section 4.

2. **Realizability-preserving numerical flux.** The transport step is the one
   place that can leave the cone. Following the mini-RAMSES result, the HLL
   wave speed is not merely a signal-speed estimate: the split states
   `U +/- F/a` are tested for `lam_min(P) >= 0` and `a` is doubled until both
   are admissible. This changes numerical viscosity and the timestep, not the
   state, and it is the correct statement of transport positivity for a moment
   system (bound the flux, not the speed).

   **Implemented — with two corrections that are not optional, both found by
   measurement rather than by inspection.**

   *The enlargement must be capped.* `dt` comes from `dfmm_cmpdt_kernel` using
   the **un-enlarged** signal speed, so a face flux built with
   `a > a0/courant_factor` violates the very CFL condition `dt` was chosen to
   satisfy. An uncapped doubling loop is unconditionally unstable: eight
   doublings (256x) destroyed a previously exactly-conservative run in a single
   step. `DF_ABOOST_MAX = 1.25 <= 1/courant_factor`. Where the cap is not
   enough, the state is allowed to leave the cone and be *reported*, which is
   the whole point of the instrument.

   *The enlargement must be CONTINUOUS in the state, not a branch.* An
   oct-boundary face is reconstructed and solved independently by the two
   threadgroups that own the adjoining octs, and conservation depends on both
   arriving at the same flux. A discrete "enlarge if inadmissible" test turns a
   round-off difference in the cone margin into a *finite* difference in the
   wave speed, hence a non-telescoping flux. Measured: `mcons` went from 0 to
   -2.3e-2 in two steps on a case that is otherwise exact, at every `tau`, with
   the source term off, and at every resolution. The implementation therefore
   computes a smooth margin (`df_cone_margin`) and widens by a continuous ramp
   `w = 1 + (DF_ABOOST_MAX - 1) clamp(-10 m, 0, 1)`, applied symmetrically to
   `SL` and `SR`. Verified: Gates 3 and 5 are **bit-for-bit unchanged** — the
   widening never fires on smooth low-anisotropy flow — and conservation
   returns to round-off.

2b. **A permanent self-consistency guard on the diagnostic itself.**
   `P = p I + Pi` with `Pi` traceless, so Weyl's inequality forces
   `lam_min(P)/p >= 1 - ||Pi||_F/p` cell by cell. `dfmm_diag_kernel` records
   `max over cells of ((1 - ani) - lam)`, which must be `<= 0`: two flops, and
   it paid for itself immediately. It exposed a defect in which the
   realizability report divided by `max(p, 1e-30)` while building `P` from the
   unfloored `p`, so an empty cell (`rho = smallr`,
   `p = rho smallc2/gamma = 6e-31`) reported `1/gamma = 0.6` with
   `||Pi||/p = 0` — a state-independent constant that masked the true minimum
   in every run. Both now use one pressure. `DFMM_DIAG_RAW=1` dumps the raw
   accumulator slots and the offending cell.

   The minimum is also **not** carried as a single offset float. The only
   lock-free float atomics available are min/max on the raw bit pattern, which
   order IEEE-754 correctly for non-negative values only; minimising
   `lam/p + 2` silently clips anything below `-2` — invisible in the gate
   problems, catastrophic here, where `|Pi|/p` reaches order ten.

3. **Diagnostics, not clipping.** `lam_min(P)/p` is recorded per level as a
   minimum (not a maximum of slack — a max is structurally blind to cells at
   the cone boundary). A negative value is reported and, under
   `dfmm_fatal_realizability=.true.`, aborts the run with the cell count and
   violation magnitude. There is no clip, projection, rollback or accepted
   negative tolerance.

Note on what is *not* claimed: the 1D reference is not repair-free either — it
recovers `gamma` with a `max(., 0)` and clips `|beta| <= 0.999 sqrt(Svv)`. The
Stage-4 phase-space sector inherits that problem and will be handled by a
Schur-complement eigenvalue audit on `Svv - Sxv^T Sxx^{-1} Sxv`, reported the
same way.

### Phase-space realizability

`M` is realizable exactly when `Sxx > 0` and the Schur complement

```
Gamma = Svv - Sxv^T Sxx^-1 Sxv
```

is positive semidefinite. In 1D, `Gamma = Svv - Sxv^2/Sxx = Svv - beta^2 =
gamma^2`, so `Gamma` is the direct generalisation of the reference's
rank-collapse variable, and

```
g = sqrt( lam_min( Svv^-1 Gamma ) )
```

generalises its normalised indicator `gamma/sqrt(Svv)`. `g = 1` means position
and velocity are uncorrelated; `g -> 0` means the packet has collapsed onto a
phase-space filament, at which point no Gaussian closure can describe it. This
is the paper's second indicator. It is computed as `lam_min` of
`chol(Svv)^-1 Gamma chol(Svv)^-T`, a congruence transform with the same
spectrum as `Svv^-1 Gamma` but symmetric, so the existing closed-form
eigenvalue routine applies.

Reported as the **minimum** over cells, for the same reason as
`lam_min(P)/p`: a maximum of the slack is structurally blind to the cells that
have already collapsed. `lam_min(Gamma) < 0` is counted and reported, and --
unlike the pressure cone -- it is **never** repaired, because nothing in the
update requires `Gamma >= 0`. The exact Liouville flow does preserve `M >= 0`
(it is `M(t) = Phi M(0) Phi^T` with `Phi' = Jac Phi`), but only if all three
blocks obey the Liouville equations. At Stage 4 `Svv` does not: it additionally
carries `div Q`, the BGK relaxation of `Pi`, and its own transport. So
`Gamma < 0` is the statement that the pressure tensor the moment system
delivers is no longer consistent with any single Gaussian phase-space packet
-- which is exactly the closure failure the indicator exists to report. Gate 9
measures it at `K = 3` and confirms it in float64, not as a float32 artifact.

There is one numerical caveat worth stating. `Gamma` is a difference of two
quantities that approach each other as the packet spreads, so once
`g^2 < ~1e-6` float32 cannot resolve its sign. In the collisionless limit `g`
decays as `sigma_x0/(t sqrt(theta))` without bound -- free streaming genuinely
correlates position and velocity completely -- so a long collisionless run
will reach that floor. It is a floor on the diagnostic, not on the state.

### Wave speeds

Along a face normal `n`, with `P_nn = n_i P_ij n_j`:

```
c_n = sqrt( (3 + sqrt 6) P_nn / rho )
```

`3 + sqrt 6 ~= 5.449` is the reference scheme's `CSCOEF`, the largest
characteristic of the Wick-closed four-moment subsystem. It bounds the
ten-moment value `sqrt(3 P_nn / rho)`, so it was safe at Stage 1 (a deliberate
factor 1.35 of margin, chosen so the timestep would not change when `Q` landed)
and is the correct value at Stage 2. Anisotropy therefore enters the timestep directly, which is intended:
`Pi` growing along one axis shortens `dt`.

---

## 6. Metal kernel architecture

The constraint that shapes everything: **Apple GPUs allow 32 KiB of
threadgroup memory per threadgroup.** The existing hydro kernel holds a
6x6x6 primitive stencil plus six 3x2x2 interface arrays in threadgroup memory:

```
per field:  6^3 * 4 B (stencil) + 6 * 12 * 4 B (interfaces) = 1152 B
```

so the budget is `32768 / 1152 = 28` fields. Measured cost by stage:

| stage | NVAR | threadgroup bytes | fits |
|---|---|---|---|
| 1 | 10 | 11 584 | yes |
| 2 | 20 | 23 104 | yes |
| 3 | 23 | 26 560 | yes |
| 4 | 38 | 43 840 | **no** |

Stage 4 does not fit, so **Stages 3 and 4 both use a two-kernel split**,
reusing the pattern already present in this repository for MHD UCT face reuse
(`mtl_godunov` -> `uct_face_product` kernel writes a global face buffer ->
`uct_reuse` kernel consumes it). Stage 3 would fit in one kernel, but using
the same split at both stages means one code path, and it means the Stage-3
numerics are unchanged when the Stage-4 fields land:

* `dfmm_integrator_kernel` — the 20-field hyperbolic core, unchanged from
  Stage 2, plus a writeout of the face mass fluxes to a global buffer. That
  buffer is 36 floats per oct (12 faces per direction), indexed by the
  *relative* subgrid index within a dispatch, so it is `num_subgrids` long,
  not `ngridmax`. Written *after* the fine-face zeroing, so a face hidden
  behind refinement arrives as an exact zero and the tower inherits the AMR
  bookkeeping for free.
* `dfmm_passive_kernel` — advects the mass-like tower using those stored
  fluxes, so `rho X` transport is exactly consistent with `rho` transport by
  construction rather than by recomputation.

Measured threadgroup cost of the passive kernel, `(4 + DF_NMASS)` stencil
fields (rho, u_i, then the tower) plus `DF_NMASS` interface fields:
6976 B at Stage 3, 24256 B at Stage 4. Both fit.

**Upwind, not HLL, for the tower.** Running `rho X` through HLL as an extra
field gives `F = a_L X_L + a_R X_R` with `a_L + a_R = F[rho]` but `a_R < 0`,
so a static fluid (`F[rho] = 0`, `S_L = -c`, `S_R = +c`) still produces the
flux `(c/2)(X_L - X_R)`. That smears `X` at the sound speed where nothing
moves. `Pi` and `Q` have to accept that -- their fluxes are genuinely
non-advective, so they must ride the acoustic waves -- but for a field whose
exact evolution is `D X / Dt = 0` it would destroy `d L_i / d x_j` within a
sound-crossing time, which is the one thing Stage 3 exists to measure. The
tower therefore uses

```
F[rho X] = F[rho] * X_upwind ,   upwind side chosen by sign(F[rho]),
```

the standard finite-volume treatment of an advected scalar, exact for
`F[rho] = 0`. It needs no frame rotation, since it is a scalar multiple of a
lab-frame flux.

Field storage inside the kernels is `float f[N][6][6][6]` indexed by field
rather than named struct members, so the tensor loops stay compact.

### Precision

The Metal path is float32 (`NPRE=4`). The paper's `gamma` diagnostic spans six
decades, and the blowup problem is deliberately ill-conditioned. Stage 1--2
quantities (`Pi/p`, `q/(p c_s)`, `lam_min(P)/p`) are all O(1) ratios and are
safe in float32. The Stage-4 phase-space sector is not obviously safe and its
dynamic range will be measured against a float64 CPU reference before any
claim is made from it.

---

## 7. Staging and gates

Each stage is independently useful and has to pass its gate before the next
starts.

**Stage 1 — ten-moment Gaussian closure.** `NVAR=10`. **Implemented and
gated.** Adds `Pi_ij`, its flux contributions to momentum and energy, the
strain source, asymptotic-preserving BGK relaxation, and the `lam_min(P)` /
`dev_NS` diagnostics.

Gates, all run on an Apple M3 with `COMPILER=METAL NDIM=3 HYDRO=1 DFMM=1`:

| Gate | Setup | Criterion | Result |
|---|---|---|---|
| 1 Euler reduction | `sedov3d`, `gamma=5/3`, `riemann=hll`, `dfmm_source=.false.` | conservation matches baseline | `mcons=0`, `econs=4.97e-8` — **identical** to the `DFMM=0` baseline |
| 2 Relaxation | uniform box, `u=0`, `Pi_xx=0.1`, `tau=0.05`, 40 steps | `\|Pi\|` decays as `exp(-t/tau)` | max relative error **3.3e-4** over 2.5 decades of decay, at the 4-digit print resolution |
| 3 Navier--Stokes limit | shear `V0=0.1`, `tau=1e-4`, level 4/5/6 | `\|Pi\| -> 2 p tau \|S0\|` and `dev_NS -> 0` | `\|Pi\|`/analytic = 1.0102 / 1.0046 / 1.0015; `dev_NS/\|Pi\|` = 7.7e-3 / 4.1e-3 / 1.9e-3 — both converge under refinement |
| 4 Collisionless | shear `V0=0.5`, `tau=-1`, 60 steps | stays realizable | `\|Pi\|/p` grows to 0.54, `min lam(P)/p = 0.644 > 0`, `n(lam<0)=0`, `econs=3.8e-8` |

Reproduce with `namelist/dfmm_shear3d.nml` (see its header for the per-gate
parameter settings).

Two results worth stating explicitly:

* **The `dt` reduction relative to Euler is physics, not a defect.** The
  ten-moment system's fastest characteristic along a face normal is
  `sqrt(3 P_nn/rho)`, not the adiabatic `sqrt(gamma p/rho)`; with the pressure
  tensor free, compression along one axis does not share energy with the
  transverse directions on the fast timescale. Stage 1 additionally uses the
  Stage-2-correct `CSCOEF = 3 + sqrt(6)`, so `dt` is a further factor 1.35
  below the ten-moment requirement. That is deliberate — it avoids a
  wave-speed change when `Q_ijk` lands.
* **`dev_NS` behaves as an indicator should.** It is 0.8% of `\|Pi\|` at
  `tau = 1e-4` (Navier--Stokes adequate) and rises to **76%** at `tau = 0.1`
  (Navier--Stokes inadequate), on the same flow. This is the quantity the
  blowup study reads.

**Stage 2 — evolved third moment.** `NVAR=20`, `DFMM=2`. **Implemented and
gated.** Adds `Q_ijk`, the Wick fourth moment, `tau_q`, and `‖q - q_CE‖` as the
paper's first indicator.

Gates, all run on an Apple M3 with `COMPILER=METAL NDIM=3 HYDRO=1 DFMM=2`:

| Gate | Setup | Criterion | Result |
|---|---|---|---|
| 5 Fourier limit | isobaric `rho = 1 + 0.1 sin(2 pi z)`, `p = 1`, `tau = 1e-4`, `Pr = 2/3`, level 4/5/6 | `q -> -(5/2) tau_q p grad theta` | `\|q\|`/analytic = **1.093 / 1.022 / 1.002**; `dev_q/\|q\|` = 5.6e-2 / 2.5e-2 / 1.3e-2 |
| 3' NS limit, restated | shear `V0=0.1`, `tau=1e-4`, level 4/5/6 | unchanged from Stage 1 | `\|Pi\|`/analytic = 1.0100 / 1.0045 / 1.0015; `dev_NS/\|Pi\|` = 7.8e-3 / 4.0e-3 / 1.9e-3 — **identical to Stage 1** |
| 2' Relaxation | `tau=0.05`, `Pi_xx=0.1`, 40 steps | `exp(-t/tau)` decay | step ratio 0.8671 vs `exp(-dt/tau) = 0.86699`, rel. error 1.3e-4 |
| 4' Collisionless | shear `V0=0.5`, `tau=-1`, 60 steps | stays realizable | `min lam(P)/p = 0.9475 > 0`, `n(lam<0)=0` |

The Fourier gate is a clean isobaric conduction test: with `p` uniform there is
no pressure gradient, so the gas stays nearly static (the thermal diffusion
time `1/chi = 1/(1.5 tau) ~ 6700` is four orders of magnitude beyond the run)
while `q` relaxes to its Fourier value within the first step. The measured
amplitude converges to the analytic one at 0.2%, which is the statement that
`tau_q = tau_Pi / Pr` with `Pr = 2/3` reproduces the note's hard-sphere
conductivity `kappa/(rho c_p) = (3/2) nu`.

Note on the Stage-2 collisionless result: `\|Pi\|/p` saturates at 0.13 rather
than the Stage-1 value 0.54 on the same flow, with `\|q\|/(p c_s)` reaching
0.20. That is not a regression — with the third moment free, collisionless
anisotropy is partly carried away as heat flux instead of accumulating in
`Pi`. It is a physical difference between the ten- and twenty-moment
collisionless response, and the realizability margin is what the gate tests.

The 1D reduction against `~/dfmm/py-1d` Sod and cold-sinusoid is *not* yet
run; it is the one Stage-2 gate outstanding and is recorded as such in
Section 8.

**Stage 3 — Lagrangian displacement.** `NVAR=23`, `DFMM=3`. **Implemented and
gated.** Adds `rho D_i` with `D D_i/Dt = -u_i`, so
`d L_i/d x_j = delta_ij + d D_i/d x_j` gives the deformation tensor, and with
it the accumulated compression factor `sigma_max(dL/dx)` and a properly
defined local Knudsen number.

**Stage 4 — phase-space covariance.** `NVAR=38`, `DFMM=4`. **Implemented and
gated.** Adds `rho Sxx_ij` and `rho Sxv_ij`, the Liouville sources, the
asymptotic-preserving relaxation of `Sxv`, and the Schur-complement rank
indicator `g`.

Gates, all on an Apple M3 with `COMPILER=METAL NDIM=3 HYDRO=1 DFMM=4`:

| Gate | Setup | Criterion | Result |
|---|---|---|---|
| 6 Uniform advection | `dfmm_ic_uadv=0.3`, `tau=-1`, level 4, `t=0.2024` | closed form in a flow with no velocity gradient | `D_i` vs `-u t` to **2.0e-8** (3 ULP); `Sxv = theta t I` to 2.2e-7 with off-diagonals **exactly zero**; `Sxx = sigma_x0^2 + theta t^2` to 6e-8; `g` vs `sigma_x0/sqrt(sigma_x0^2+t^2)` to 2.6e-5; every field uniform to the bit |
| 7 Frozen shear | `dfmm_ic_shear=0.01`, `dfmm_source=.false.`, `tau=-1`, level 4/5/6 | isolates `-Sxv G^T`: `Sxv_zx = -g theta t^2/2`, `Sxx_xz = -g theta t^3/6` | `Sxv_zx`/exact-discrete = 0.9695 / 0.9957 / **1.00003**, spread 5.4e-2 / 1.6e-2 / **4.2e-3**; `Sxx_xz` = 0.9766 / 0.9973 / **1.00024**, spread 4.2e-2 / 1.1e-2 / **2.3e-3**; `Sxv_xz`, `Sxv_xy`, `Sxx_xy` **exactly zero**; `\|Pi\|/p <= 3e-8` throughout |
| 8 Diffusive limit | `dfmm_ic_uadv=0.3`, `tau=0.01`, level 4, `t=0.506` (`t/tau=51`) | AP map must give `Sxv -> tau theta` and `Sxx -> sigma_x0^2 + 2 tau theta t` | `Sxv = 0.0099999979` vs `0.01` (**2.1e-7**); `Sxx` to 6.4e-4; implied `D = 0.0098` vs `nu = tau theta = 0.01`; `min g = 0.9951`, `n(Gamma<0) = 0` |
| 9 Strain box | `INIT=BLOWUP`, Family-B `K = 0.1 / 1 / 3`, level 5, one deformation time | nontrivial 3D advection and deformation; indicators must respond | `mcons = 0`, `econs <= 1.1e-7`; `sigma_max(dL/dx)` = 2.24 / 2.93 / 3.34 so `Kn_local` = 0.051 / 0.212 / 0.418; `min g` = 0.987 / 0.625 / **0.000** with `n(Gamma<0)` = 0 / 0 / **448**; `min lam(P)/p` = 0.957 / 0.655 / 0.365 (K=1 value reproduces the Stage-2 ladder) |
| 3', 5' Regression | Gates 3 and 5 rerun under a `DFMM=4` binary | the hyperbolic core must be untouched | `max \|Pi\|/p` = 8.975 / 8.926 / 8.899e-5 and `max \|q\|/(p c_s)` = 1.960 / 1.872 / 1.847e-4 at level 4/5/6 -- **bit-for-bit** the Stage-2 values |

Reproduce Gates 6-8 and 3'/5' with `namelist/dfmm_shear3d.nml` and Gate 9
with `namelist/dfmm_blowup3d.nml`; the incompressible counterparts of
Gates 6-9 are `namelist/taylorgreen3d.nml` and `namelist/incomp_blowup3d.nml`
(`doc/incompressible.md` Section 4).

Three results worth stating explicitly:

* **Gate 6 found a real defect and the fix is in Section 4.** Before the
  trapezoidal pairing, the first six printed values of `g` matched
  `1 - n^2 dt^2/(sigma_x0^2 + n(n-1) dt^2)` to five digits -- confirming the
  Liouville sector -- and then crossed zero at step 11. The agreement is what
  made the diagnosis possible: the discrete solution was exactly right, so the
  scheme was wrong.
* **`Gamma < 0` at `K = 3` is real, not round-off.** Recomputed from the
  snapshot in float64, `lam_min(Svv^-1 Gamma)` reaches `-0.141` in 448 of
  32768 cells (1.4%), with Cauchy--Schwarz saturation 1.14. At `K = 1` the
  margin is comfortable (`min = 0.390`). So the second indicator fires one
  rung *after* the first, and before `lam_min(P)` does.
* **`sigma_max(dL/dx)` and `min g` are resolution-limited peaks, not
  converged numbers.** At `K = 1` they run 2.29 / 2.93 / 3.44 and
  0.716 / 0.625 / 0.529 at level 4/5/6: each is an extremum over cells of a
  quantity that keeps sharpening, so refining finds more of it. `Kn_local` is
  therefore a *lower bound* that grows with resolution -- which is itself the
  statement the note's check 1 makes. By contrast `min lam(P)/p` converges
  cleanly (0.6836 / 0.6554 / 0.6486, differences 2.8e-2 then 0.7e-2), and the
  Lagrangian residual `|rho/det(dL/dx) - 1|` falls 9.7e-2 / 7.5e-2 / 4.3e-2,
  so the deformation map itself carries a ~4% error at level 6.

**Two build traps worth recording.**

* `bin/Makefile` does not treat `INIT=` as a dependency of `condinit.o`, so
  changing `INIT=` without `make clean` silently relinks the *previous*
  problem's initial condition. Doing this produced a uniform static state for
  a blowup namelist -- `ekin = eint`, `|Pi|/p = 0` -- which looks exactly like
  a broken initial condition. Always `make clean` when changing `INIT=`.
* An unrecognised `DFMM=` value used to fall through to `NDFMM = 0` silently.
  This bit a build-verification sweep written in zsh: zsh does **not**
  word-split an unquoted `$var`, so `make ... $cfg` with
  `cfg="DFMM=4 INIT=BLOWUP"` reaches make as **one** argument, make reads it as
  `DFMM = "4 INIT=BLOWUP"`, no `ifeq` branch matches, and the build succeeds as
  a plain-hydro binary. Two sweeps reported "all stages build" while building
  `NDFMM = 0` every time. The selector now `$(error)`s on any value outside
  `0..4`, and the sweep greps the recorded `-DNDFMM=` out of the build log
  rather than trusting the exit code.

**Stage 5 — the blowup problem.** Initial and forcing conditions from the
construction, run as a compressible gas at prescribed `Kn`, with all three
indicators recorded. The deliverable is the *ordering* of the failures: which
of the note's checks fires first, at what `Delta/t_c`, and whether the
measured `lam_min(P)` crosses zero where the Newtonian extrapolation says it
should.

---

## 8. What is verified vs. assumed

Verified in this repository at the time of writing:
* Baseline `COMPILER=METAL NDIM=3 HYDRO=1` builds and `sedov3d` runs to
  completion with `mcons = econs = 0`. Required a build fix: the gfortran-driven
  link cannot resolve clang's `objc_msgSend$<selector>` stubs, so the bridge is
  now compiled with `-fno-objc-msgsend-selector-stubs` and linked with `-lobjc`.
  This was a pre-existing break on this toolchain, not caused by dfmm.
* Threadgroup-memory budget arithmetic in Section 6, from the existing
  `local_subgrid_t` / `interfaces_*_t` declarations.
* All four Stage-1 gates, all four Stage-2 gates, and all five Stage-3/4
  gates in Section 7, from fresh runs.
* That the hyperbolic core is untouched by Stages 3 and 4: Gates 3 and 5
  reproduce their Stage-2 numbers bit-for-bit under a `DFMM=4` binary.
* That the incompressible rungs are the incompressible limit of the *same*
  moment system, not a lower one: `doc/incompressible.md` Gates 6 and 8. The
  cross-rung consistency signal is `min g(rank)` at `K = 1`, 0.619 from the
  host float64 spectral incompressible rung against 0.625 from the Metal
  float32 finite-volume compressible one.
* That every stage selector still builds after the shared-file edits that
  Stages 3, 4 and the incompressible rungs made to `amr_commons.f90`,
  `read_params.f90`, `bin/Makefile`, `pm/newdt_fine.f90`,
  `hydro/courant_fine.f90`, `hydro/godunov_fine.f90`, `hydro/condinit.f90` and
  `hydro/hydro_parameters.f90`: clean builds of `DFMM=1/2/3/4 INIT=DFMMTEST`,
  `DFMM=0/4 INIT=BLOWUP` and `DFMM=0/4 INIT=TAYLORGREEN`, each confirmed to
  carry the intended `-DNDFMM=` rather than only to exit zero. See the second
  build trap in Section 7 for why that last clause is not redundant.
* That AMR prolongation (`refine.metal` / `interpol_hydro.f90`) and
  restriction (`upload_kernel`) already treat all `NVAR` fields with a
  conservative linear interpolation and a plain volume average, which is the
  correct handling for the density-like `Pi`. No dfmm-specific AMR code was
  needed. Note this is verified by inspection, not yet by a refined run.

Corrected during Stage 1, recorded so the reasoning is not lost:
* `cons_from_prim` / `prim_from_cons` multiply `ivar > 5+nener` by `rho`
  ("passive scalar density"). That is wrong for `Pi`, which is density-like.
  The loops now stop at `nvar-ndfmm`.
* RAMSES's `ekin` reporting slot carries the *total* energy, not the kinetic
  part: `update_time.f90` forms `econs` from `g%ekin_tot` alone. Reporting the
  true kinetic energy there showed a spurious `econs ~ 0.86`.
* The strain source must be integrated *together* with the BGK sink
  (Section 5). Naive splitting passed the relaxation gate but silently failed
  the Navier--Stokes limit whenever `dt >> tau`.

Corrected during Stage 2:
* The asymptotic-preserving map must relax toward the *whole* non-stiff
  right-hand side, transport rate included (Section 5). The Stage-1 form
  relaxed only toward the production term. For `Pi` that is an O(tau^2)
  omission and Gate 3 is unchanged by fixing it; for `Q` it would have
  replaced the Fourier heat flux by an unrelated quantity, because `R_ijkl`
  in the `Q` flux is O(1). This is the same class of error as the Stage-1 one
  and was found by asking what the stiff limit of *each* term is, not by a
  failing gate.
* `T_Q1` and `-d_l R_ijkl` must both be present in the MUSCL predictor or
  both absent: each is O(1) and they cancel to O(grad theta). The predictor
  uses their analytically combined form so the cancellation is exact there.

Corrected during Stages 3 and 4:
* `Sxx` must use the trapezoidal average of `Sxv` across the step, not its
  value at the start (Section 4). Found by Gate 6, where the forward-Euler
  form drove the rank indicator through zero at a predictable step number.
* The Lagrangian label must be stored as the displacement `L_i - x_i`, not as
  `L_i` (Section 2), or the periodic wrap poisons the deformation tensor.
* The tower must be advected by upwinding on the stored face mass flux, not by
  passing `rho X` through HLL as an extra field (Section 6).
* `output_hydro.f90` and `cons_from_prim` / `prim_from_cons` now distinguish
  the density-like block from the mass-like tower, which is why
  `ndfmm_dens` / `ndfmm_mass` exist in `hydro_parameters.f90`. The Stage-1 fix
  recorded above ("the loops now stop at `nvar-ndfmm`") was correct only while
  every dfmm field was density-like.

Outstanding:
* The 1D reduction against `~/dfmm/py-1d` Sod and cold-sinusoid. The 3D tensor
  algebra has been checked against the reference *analytically* (`R_xxxx =
  3 P_xx^2/rho`, `CSCOEF`, the `q_i = Q_ijj/2` contraction, and in Section 4
  the full `alpha`/`beta` <-> `Sxx`/`Sxv` equivalence) but not yet by running
  the two codes on the same problem. Note that a *quantitative* comparison of
  the moment sector is not available even in principle: `py-1d` carries
  `P_perp` with flux `u P_perp`, dropping `Q_yyx`, so it is a reduced model
  rather than the 1D restriction of the 20-moment system. The dual-frame
  sector, by contrast, is identical and Section 4 proves it.

Assumed, to be measured:
* That the Wick closure is adequate through the interesting part of the blowup
  problem. The paper's own honest finding in 1D is that the max-entropy
  upgrade "validates the framework more than it extends the working range", so
  the expectation is that Wick plus the realizability diagnostic is the right
  instrument, and that the diagnostic firing *is* the result.
* float32 adequacy for Stage 4. Partly measured now: the state variables are
  fine (Gates 6-8 agree with closed forms to 2e-7), but the *diagnostic*
  `Gamma` is a cancelling difference and loses its sign below `g^2 ~ 1e-6`
  (Section 5). Gate 9's violation was confirmed in float64 for that reason.

---

## 9. References

* `~/Downloads/moment scheme paper.pdf` (= `~/dfmm/specs/01_methods_paper.tex`
  is a *different*, later variational document; the moment-scheme paper is the
  1D Eulerian one) — Sections 2--7 give the 1D scheme this generalises.
* `~/Downloads/before_blowup_ideal_gas_pedagogical.pdf` — the five checks and
  the hard-sphere coefficients.
* `~/dfmm/py-1d/dfmm/schemes/cholesky.py`, `schemes/_common.py` — reference
  fluxes, sources, BGK map, realizability clip.
* `~/phrike` branch `dfmm`, `dfmm_plan.md` — a pseudospectral 1D port of the
  same eight-field system; useful for its operator-splitting discussion.
* `~/dfmm/reference/notes_M3_7_3d_extension.md`, `src/berry.jl`,
  `src/cholesky_DD_3d.jl` — 3D tensor algebra cross-check for Stage 4.
* Obsidian `20-projects/mini-ramses/mini-ramses — DFMM dust plan.md` and
  `— realizable multimoment dust findings.md` — prior art on moment transport
  positivity, the split-state wave-speed rule, and why max-entropy caps were
  rejected.
