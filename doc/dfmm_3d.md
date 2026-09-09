# dfmm in 3D on the Metal backend — implementation ledger

Status: **ledger frozen; Stages 1 and 2 implemented and gated.** This document fixes the
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
| 21..23 | `rho L_i` | mass-like | **3** |
| 24..29 | `rho Sxx_ij` (6 comps) | mass-like | **4** |
| 30..38 | `rho Sxv_ij` (9 comps, not symmetric) | mass-like | **4** |

`NVAR` = 10 / 20 / 23 / 38 at the four stages, selected by `DFMM=n` in
`bin/Makefile` (`NDFMM` = 5 / 15 / 18 / 33). Stage 1 stays selectable as
`DFMM=1` so its gate results below remain reproducible, and because the `Q`
sector is expensive enough to be worth switching off when it is not needed.

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

**Phase-space sector** (Stage 4). With the 6x6 phase-space covariance blocks
`A = Sxx`, `B = Sxv`, `C = Svv = P/rho`, the Liouville evolution of the local
phase-space packet under the linearised flow `Jac = [[0, I], [0, -G]]` is
`d M/dt = Jac M + M Jac^T`, giving

```
D Sxx / Dt = Sxv + Sxv^T
D Sxv / Dt = Svv - Sxv G^T
D Svv / Dt = -G Svv - Svv G^T      (already carried by the hydro sector)
```

1D reduction: `Sxx' = 2 Sxv`, `Sxv' = Svv - Sxv du/dx`, `Svv' = -2 Svv du/dx`,
matching Eqs. (9)--(11) of the paper.

**Lagrangian coordinate** (Stage 3): no source, `D L_i / Dt = 0`, initialised
to `L_i = x_i`.

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

Stages 1--3 therefore extend the single-kernel design. Stage 4 splits into two
kernels, reusing the pattern already present in this repository for MHD UCT
face reuse (`mtl_godunov` -> `uct_face_product` kernel writes a global face
buffer -> `uct_reuse` kernel consumes it):

* `dfmm_core_kernel` — the 20-field hyperbolic core, and it writes the face
  mass fluxes to a global buffer.
* `dfmm_passive_kernel` — advects the 18-field passive tower using those
  stored mass fluxes, so `rho X` transport is exactly consistent with `rho`
  transport by construction rather than by recomputation.

Field storage inside the kernels is `float f[NDFMM][6][6][6]` indexed by field
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

**Stage 3 — Lagrangian coordinate.** `NVAR=23`. Pure passive advection;
`d L_i / d x_j` gives the deformation tensor and hence the local compression
ratio and a properly defined local Knudsen number.

Gate: exact advection in uniform flow; `L_i` deviation from `x_i - u t` at
round-off.

**Stage 4 — phase-space covariance.** `NVAR=38`. Two-kernel split. Per-axis
`gamma` rank-collapse diagnostic.

Gate: 1D reduction against the reference `alpha`/`beta` trajectories;
Schur-complement positivity audit.

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
* All four Stage-1 gates and all four run Stage-2 gates in Section 7, from
  fresh runs.
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

Outstanding at Stage 2:
* The 1D reduction against `~/dfmm/py-1d` Sod and cold-sinusoid. The 3D tensor
  algebra has been checked against the reference *analytically* (`R_xxxx =
  3 P_xx^2/rho`, `CSCOEF`, the `q_i = Q_ijj/2` contraction) but not yet by
  running the two codes on the same problem.

Assumed, to be measured:
* That the Wick closure is adequate through the interesting part of the blowup
  problem. The paper's own honest finding in 1D is that the max-entropy
  upgrade "validates the framework more than it extends the working range", so
  the expectation is that Wick plus the realizability diagnostic is the right
  instrument, and that the diagnostic firing *is* the result.
* float32 adequacy for Stage 4 (Section 6).

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
