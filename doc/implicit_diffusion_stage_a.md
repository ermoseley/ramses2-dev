# Implicit non-ideal MHD: spatial prerequisite checkpoint

Status on 2026-10-05: **no implicit solver is implemented or qualified**.
This branch starts from the split non-ideal MHD framework at
`6e6fc3ec0dab40628ecdc3e2ff43411def042b4d` (`gpu-non-ideal-sts`),
the branch associated with the Claude conversation “RAMSES2 diffusion
framework port.” Existing solver modes and their defaults are unchanged.

The reference branch tested synchronized all-level explicit/STS stages at
`8cb067f6`, then reverted that coupling at `b540fe69`. The last inspected
head was `cb9a979e`, adding the STS Gegenbauer parameter `nimhd_alpha`.
This checkpoint has incorporated none of those commits. Their inspection
does not establish the implicit operator's energy form or macro scheduling.

The reviewed October 4 implicit-diffusion plan requires a consistent,
accretive composite spatial operator and a local energy partition before
solver integration. In particular, it says: “If those requirements cannot
be met, stop at the operator design; do not accept an inconsistent diagonal
substitute.” The following experiments test that prerequisite; they do not
validate an implicit RAMSES executable.

## Uniform operator

The actual frozen `cmpnimhd_mhd` stencil fails accretivity for a nonuniform,
divergence-free background and a 1:10000 AD coefficient jump. On a periodic
4-cubed mesh with unit spacing, its relative Frobenius symmetry defect is
0.1618173, and the minimum eigenvalue of its symmetric part is
-79.1930144656. Literal evaluation of the source expressions agrees with
the host matrix construction to 1.78e-15.

A CUDA microtest extracted the source routine and changed only the magnetic
samples that define the frozen tensor, separating them from the trial
field used to evaluate current. On an H200, a normalized trial vector gave
`x_H_x = -7.9193014465632473E+01`. This is a frozen-operator counterexample,
not an observed nonlinear production trajectory.

The proposed uniform repair retains native Ohmic conductances, gathers all
three current components to each edge orientation for AD, applies its
positive tensor, and deposits every component with the transpose of the
gather. Equal orientation weights of 1/3 pass the tested symmetry, energy,
nullspace and second-order smooth consistency checks. The short-wave
spectrum changes, so original STS and
STS with this revised stencil must remain separate comparison cases.

## Composite interface

A cell-local mimetic face mass is SPD and reproduces zero current for a
constant magnetic field, but fails affine-current consistency. Changing
its scalar stabilization is insufficient. A wider two-dimensional
reconstruction passes commuting-curl, face-flux conservation and
constant/affine/quadratic current checks at a finite refinement corner.
However, it uses a finite-domain moment fit and broad reconstruction
support; it is not a qualified local three-dimensional construction.

A three-dimensional refined-cube audit has 350 cells, 1206 faces, 1386
edges and 531 vertices. Its topological products `D C` and `C G` vanish
exactly. At the center of a coarse face, a hanging vertex has an incoming
normal edge on only the fine side. Consequently a positive diagonal edge
metric cannot give zero vertex divergence for every constant current.
Since `G^T C^T = 0`, changing only the face mass cannot repair this defect.
This observation does not exclude a coupled edge metric.

A fixed 78-edge block for the isolated refined cube does give a positive
edge metric that satisfies the necessary constant and affine-current
vertex constraints. Its eigenvalues are 0.801 to 3.027, and its moment
residual is 1.8e-14. This is a feasible component, not a complete operator:
the corresponding face metric and a bounded block partition for adjacent
refined parents remain to be constructed.

An independent local current reconstruction, with radius three fine cells
and at most 165 source edges per row, reproduces quadratic-field currents
to 1.53e-14 and gives positive energy through transpose deposition.
Nevertheless its complete diffusion operator has an affine-field error
0.24383 at the interface. Thus reproducing the current alone is insufficient:
the transpose electric-field deposition must satisfy its own moments.

A larger paired moment construction fixes both moment conditions on a
14-cell patch and reproduces the full Ohmic operator on quadratic fields
to 5.3e-14. It nevertheless creates sixteen extra Ohmic null directions.
For example, the curl-free field `(y,x,0)` has a nonzero adjoint-current
numerator that the fitted current map is forced to annihilate. A
moment-preserving completion cannot remove this null. The face metric
itself must preserve the polynomial gradient kernel before fitting the
current reconstruction.

The final joint face/edge fit exposes an exact obstruction for this
14-cell support. Take the harmonic potential
`phi=(x-1)(y-1)(z-1)`. Its gradient has nonzero fine-face fluxes inside the
patch, but zero normal-flux averages on every exterior coarse face. The
compact face vector has norm 1.5, exactly zero divergence, and is in the
range of the physical curl (reconstruction residual 5.6e-17). The exterior
native field contributes exactly zero adjoint current on the patch.
Reproducing zero current would therefore require this nonzero curl-range
vector to have zero magnetic energy, contradicting a positive definite
face metric. This explains the zero magnetic polynomial-energy eigenvalue
in the joint fit; it is not a least-squares convergence failure.

This certificate assumes the 66-face block and unchanged native metric
outside it. A wider construction resolves this particular obstruction.
A 3-by-3-by-3 parent neighborhood contains 34 leaf cells, 138 faces and
186 edges. Joint face/edge action moments admit positive patch shares
after subtracting the unchanged native exterior contributions. Their
minimum eigenvalues are 0.22153 and 0.56975. Constant-vector energy and
heat normalize to the patch volume, 216. The completed operator reproduces
quadratic magnetic-field currents and the full Ohmic diffusion action to
about 1.7e-13. Its positive face and edge metrics retain the physical
Ohmic nullspace. The earlier compact harmonic alias is absent on this
wider support.

Two adjacent refined parents also admit positive overlapping patch
shares, with each native cell contribution split by its patch coverage.
The two 41-cell patches have 165 faces and 220 edges each; their assembled
operator moment error is 1.2e-12. This finite joint fit establishes
feasibility for that geometry. It does not supply a bounded independent
assembly rule for arbitrary refinement patterns, nor a local inverse of
the assembled sparse edge metric. A production implementation cannot
substitute the isolated dense inverse for that missing construction.

Two independent patch-load assembly rules fail the adjacent-parent test
before positivity: their harmonic bilinear loads are not symmetric. The
joint fit's existence therefore does not justify either local rule.
A subsequent Opus proposal uses shared commuting maps to virtual fine
cells. However, an independent exact audit rejects its uncorrected splice
to native exterior metrics. Replacing one width-H native cube by eight
native half-cubes changes the magnetic energy of the curl-free affine
fields `(x,-y,0)` and `(y,x,0)` by `-H^5/4` and `H^5/8`, respectively
(before the physical factor 1/2). With unchanged exterior loads, exact
zero current would require both changes to vanish. A commuting, SPD host
example also has a constant-current error proportional to `1/H` and
nonconvergent affine-current errors. Higher interpolation order does not
repair that energy mismatch. A signed symmetric correction can cancel
the energy defects while retaining SPD in that example, but the tested
correction still fails full current consistency. This rejects that
particular splice, not all bounded local interface metrics.

Testing the actual six-face boundary correction with fixed affine-exact
traces also cancels the complete affine energy Gram and retains SPD
(minimum eigenvalue 0.80743). It nevertheless gives a constant-field
adjoint-current residual of exactly `5/64` on a seam edge. This simple
boundary form therefore does not supply the missing complete operator.

A sparse coarse/detail construction subsequently passes the isolated real
leaf-complex test without a moment fit. Commuting restriction/prolongation
maps `R`, `P` and detail projection `D=I-P R` give the positive metrics
`M=R^T M_coarse R+D^T M_native D` for faces and edges. Their current and
quadratic diffusion errors are 5.33e-15 and 4.45e-15; changed-block minimum
eigenvalues are 0.48672 and 0.73258. Root independently checked the saved
maps, mass identities, moments and eigenvalues. Changes remain confined to
the isolated 138-face/186-edge neighborhood.

Extending this formula unchanged through a large refined region changes
the native fine operator: `H P_face=P_face H_coarse`, so a prolonged fine
`(+,+,-,-)` mode damps at `1/h^2` instead of native `2/h^2`. Adjusting only
the positive detail masses cannot change that identity. The construction
therefore remains an isolated prototype under the native-interior contract.

A genuine scalar one-dimensional coarse/fine transition exposes a separate
symmetry obstruction. Integrated interval fields use native face mass
`1/length` and node mass `(left_length+right_length)/2`. Any compact symmetric
correction retaining those exterior masses fails simultaneous constant,
affine and quadratic current moments. The exact witness `(x^2,-2*x,1)`
annihilates every correction by symmetry but pairs with the required load
to -1; root reproduced this certificate. This result uses the physical
nonuniform mesh, without a coarse/detail projection or translated caps.
A separate mathematical review extends this obstruction to a periodic
planar three-dimensional interface with compact metric changes and native
diagonal exterior metrics. Even affine-current consistency plus quadratic
Ohmic-operator consistency is incompatible there. The adjoint identity
closes the possible curl-free current-error escape. This conditional result
does not exclude alternative bulk metric factorizations or prove an AD-only
impossibility. Opposing interface loads in an isolated refined island can
cancel, explaining why the finite cube fits remain feasible.

The pointwise interface-operator gate is stronger than solution convergence.
A separate four-grid steady scalar test with native metrics has bounded
`O(1)` layer truncation yet finest volume-weighted L1/L2 solution rates
2.0064/1.9997. Its integrated linear residual is 3.86e-14. Root inspected
the source and independently checked the saved solutions and residuals.
This is scalar one-dimensional evidence, not a three-dimensional AD or
time-dependent qualification. Alternatively, consistent nodal P1 edge
masses give exact quadratic currents on the nonuniform chain, but change
regular bulk metrics and evolution. Independent Astra and Opus 5.5 high
reviews confirm the conditional obstruction and recommend the
solution-convergence gate.

On October 5 the user adopted that amendment, retaining native bulk.
The interface must have at least first-order physical current/EMF
consistency, exact constant-field current, symmetric positive energy,
the physical nullspace and exact topology. Bounded conservative `O(1)`
interface-operator errors are allowed only with demonstrated second-order
global L1/L2 solution convergence for the actual three-dimensional Ohmic
and AD operators, including planar interfaces, refined cubes, edges,
corners and the declared coefficient variants. Complete-operator layer
errors remain reported. The scalar example and mathematical reviews do
not satisfy those tests. Stage A construction resumes under this contract;
solver integration remains conditional on its spatial, thermodynamic
and cost gates.

The subsequent Opus 5.5 high and Astra review found published elliptic
finite-volume precedent for second-order solution convergence with
zeroth-order local truncation. That precedent does not qualify CT/AD.
The fixed-domain 3D tests must measure the evolved magnetic field and its
actual reconstructed current, EMF and deposited heat, reporting their
orders separately. Exact curl structure and a conserved total do not
replace quantitative weak cancellation or solution convergence. Include
the full operator halo, deep native interiors, oblique/tangential pure AD
and the declared coefficient variants; shrinking an isolated fixture
cannot supply this evidence. If the required accuracy fails, consistent
compatible mass metrics with sparse mixed equations are a fallback, with
a new 3D construction and the existing thermal/cost gates still required.
See [Diskin and Thomas's primary examples](https://ntrs.nasa.gov/api/citations/20110016434/downloads/20110016434.pdf).

A finite positive AD reconstruction has also been tested on the wider
isolated patch. Its gather and full weighted transpose reproduce the
tested axis and oblique quadratic-field actions to about 7.4e-14 and
retain the resolved constant/affine parallel-current nulls. A smooth
extended stationary field has second-order current and first-order
interface-operator errors under homothetic refinement. This is local
consistency evidence, not general AMR assembly, variable nonlinear
coefficient qualification or fixed-domain PDE convergence.

Stage A remains incomplete: general patch assembly, extended AD
nullspace/variable-coefficient qualification, thermal admissibility, and
fixed-domain spatial convergence still require qualification before
solver integration.

## Thermal accounting

For the revised uniform stencil, a fixed local redistribution of edge work
gives nonnegative quadrature heat and a conservative Poynting flux. Host
checks reproduce the cell energy identity, the common-stage time defect,
and inexact-stage residual work to floating-point precision; the smooth
flux test approaches second order.

The compatible face energy differs from the magnetic energy subtracted by
RAMSES when recovering primitive pressure. The change in that staggering
difference can exceed positive Joule heat. Acceptance therefore still
requires raw internal-energy checks and atomic bounded retries, as the
plan specifies.

The wider isolated patch admits a positive volume-weighted cell partition
of its face and edge metric shares. Smooth harmonic and nonharmonic
quadratic tests give at least first-order interface heat consistency and
second-order Poynting-flux corrections under local homothetic refinement.
Native exterior Yee flux loads are essential to the exact cell balance.
An arbitrary divergence-free field, with a deliberately perturbed current,
closes the ledger including signed constitutive residual work to 2.5e-14.
That routing uses a fixed 124-cell patch-plus-native-halo support. These
are local truncation and algebraic tests, not fixed-domain PDE convergence
or pressure/recovery qualification.

The two adjacent patches also pass a finite heat/flux test using separate
bounded trees and a shared overlap anchor. Nonzero individual patch loads
cancel at that anchor. A changing quadratic field gives first-order heat
and second-order flux errors, while native exterior fluxes are unchanged.
The random-field residual ledger closes to 2.31e-14. The measured flux
error coefficient is large (39.23 at unit scale), so order alone does not
establish a useful production accuracy envelope. General edge-local
physical-reference/cycle consistency remains a separate gate; a growing
patch-component tree is neither implemented nor required for conservation.

The isolated volume-weighted partition fails a smooth low-beta raw-pressure
test. With background Bz=1, a compact smooth CT perturbation of amplitude
1e-3, pure Ohmic mobility and initial internal-energy density 1e-10,
exact Padé endpoints at steps .01, .005 and .0025 have minima
`-1.19277e-5`, `-5.99405e-6` and `-3.00455e-6`. The signed staggering
change is linear in perturbation amplitude, while Joule heat is quadratic;
stage and energy-ledger errors are around 1e-14. These are three advances
from the same state, not a successful retry sequence. Repartitioning this
same compatible metric cannot universally remove the leading defect
under the unchanged primitive-energy convention, although finite-pressure
acceptance envelopes remain possible.

A separate signed conservative redistribution through a fixed 34-cell
tree repairs this case's patch cells without changing CT, exterior fluxes
or incoming internal-energy reservoirs. It is not accepted as a production
change. Consecutive .0025 substeps reject the second substep because eight
unchanged exterior cells become negative. Thus it does not complete the
four-substep recovery of a .01 interval. Always applying that correction
also changes the energy generator at first order in time. Positive metric
heat and conservation alone do not qualify pressure recovery.

## Hardware evidence and remaining work

S3DF CPU build job `39865645` and Hopper runtime job `39865646` completed
with all steps exiting zero. Runtime used one H200 NVL on `sdfhopper003`,
UUID `GPU-e2ab22ff-8d0f-babe-6394-aee3b28854f5`, NVHPC 25.5, FP64,
`nsubgrid=2`, and `-cuda -gpu=cc90,nofma -O0`.

The CUDA audit source SHA256 is
`fa47b7950ca990fbfecd9001f2afb99ac6db15e857b398a96e6238c1109081fb`;
the binary SHA256 is
`ccbfbe323c2d03d1dc5617de052644efdedd59e5167c968e236fdb319b1fd174`.
Remote evidence is under
`/sdf/scratch/users/e/emoseley/ramses2-implicit-20261004`.
The local scripts, logs and manifests are preserved outside the repository
in `/Users/moseley/ramses-development/artifacts/2026-10-04-implicit-stage-a`.

The wider positive Ohmic candidate subsequently passed a separate CUDA
audit: CPU build `39870588` and one-H200 runtime `39870589`, with every
step completing at exit zero. The GPU was an H200 NVL on `sdfhopper003`,
UUID `GPU-6c476df3-8287-b58d-4d1e-c5a0fefae691`; the compiler and flags
were again NVHPC 25.5, FP64, `-cuda -gpu=cc90,nofma -O0`.
Across 33 cases, relative GPU/host current and operator errors were
8.33e-17 and 9.71e-17; quadratic operator error was 1.62e-13,
weighted energy/symmetry defects were 2.27e-13/1.14e-13, and divergence
and weighted-gradient null errors were 3.68e-15/1.21e-14.
This executes the actual face metric, transpose curl, fixed edge-block
inverse and curl on device, with native diagonal exterior actions.
It is an isolated operator audit, not an implicit solve or a performance
measurement. Its source SHA256 is
`720024c932dabbc1ed91d7094f3969d06fe8b628829109a0e779f128cce7c2a9`;
binary SHA256 is
`089d9a4767c8fbe198249872ed7ca57f46fad162d2c60b9b200d4bc7f49c3f28`.
Remote evidence is in
`/sdf/scratch/users/e/emoseley/ramses2-implicit-wide-20261004`,
with local copies in the artifact directory's `wide-cuda` subdirectory.

The next required result is a compact, geometry-local 3D interface
construction that passes complete-operator consistency, positive energy,
native-interior and thermodynamic-partition tests. Neither an implicit
namelist mode nor a level-local substitute is enabled before that result.
After it passes, the existing coefficient hook and solver selector provide
the integration points. The synchronized macro step must settle and lock
the all-level timestep before its first diffusion half step, suppress
level-local diffusion, restrict total energy conservatively, and restore
the entire macro-step state on retry. The present timestep routine mutates
state, and the present upload can restrict internal energy; neither can
be reused unchanged for that contract.
