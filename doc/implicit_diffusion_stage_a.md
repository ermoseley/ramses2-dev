# Implicit non-ideal MHD: spatial prerequisite checkpoint

Status on 2026-10-04: **no implicit solver is implemented or qualified**.
This branch starts from the split non-ideal MHD framework at
`6e6fc3ec0dab40628ecdc3e2ff43411def042b4d` (`gpu-non-ideal-sts`),
the branch associated with the Claude conversation “RAMSES2 diffusion
framework port.” Existing solver modes and their defaults are unchanged.

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
outside it. It does not rule out a larger reconstruction patch or another
compatible discretization. Such a construction, its assembly for adjacent
refined parents, and its energy partition remain unresolved. Stage A has
not passed, so solver integration stops at the boundary specified by the
plan.

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
plan specifies. An interface heat partition remains unqualified.

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
