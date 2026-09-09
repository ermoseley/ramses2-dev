/*
 * metal/kernels/dfmm.metal
 *
 * dfmm (dual-frame moment method) in 3D — Stage 1: ten-moment Gaussian
 * closure.  See doc/dfmm_3d.md for the frozen equation, flux, source and
 * realizability ledger; this file implements Sections 3, 4 and 5 of it.
 *
 * Scope: NDIM=3, float32, HLL (or LLF) numerical flux only.  HLLC is not
 * provided: its middle-state construction is defined for the Euler system and
 * has no standard contact reconstruction for the anisotropic pressure, and the
 * 1D reference implementation this generalises is HLL.
 *
 * State (NVAR = 5 + NDFMM, NDFMM = 5):
 *   ivar 1     rho
 *   ivar 2..4  rho u_i
 *   ivar 5     E = rho|u|^2/2 + 3p/2          (gamma = 5/3 enforced)
 *   ivar 6..10 Pi_xx, Pi_yy, Pi_xy, Pi_xz, Pi_yz   with Pi_zz = -(Pi_xx+Pi_yy)
 *
 * This is a self-contained translation unit: it deliberately does not include
 * hydro.metal, so the baseline hydro kernels cannot be perturbed by anything
 * here.  The few accessors it shares with hydro.metal are redeclared static.
 *
 * Global buffer layout:
 *   dfmm_cmpdt:       [0]=grid [1]=uold [2]=data_buf [3]=head_idx [4]=num_octs
 *                     [5]=dx [6]=smallr [7]=smallc2 [8]=courant_factor
 *                     [9]=constant_gravity(3)
 *   dfmm_integrator:  [0]=grid [1]=uold [2]=unew [3]=nbor [4]=head_idx
 *                     [5]=num_subgrids [6]=ngridmax [7]=ilevel [8]=levelmin
 *                     [9]=levelmax [10]=smallr [11]=smallc2 [12]=dt [13]=dx
 *                     [14]=slope [15]=riemann [16]=constant_gravity(3)
 *                     [17]=father
 *   dfmm_source:      [0]=grid [1]=uold [2]=unew [3]=nbor [4]=head_idx
 *                     [5]=num_octs [6]=smallr [7]=smallc2 [8]=dt [9]=dx
 *                     [10]=tau_pi [11]=source_on
 *   dfmm_diag:        [0]=grid [1]=uold [2]=nbor [3]=diag [4]=head_idx
 *                     [5]=num_octs [6]=smallr [7]=smallc2 [8]=dx [9]=tau_pi
 */

#include <metal_stdlib>
#include <metal_atomic>
using namespace metal;

#include "../metal_types.h"
#include "metal_utils.h"

#ifdef DFMM

/* ---------------------------------------------------------------------------
 * Compile-time constants
 * --------------------------------------------------------------------------*/
#ifndef NVAR
#define NVAR 10
#endif
#ifndef NDFMM
#define NDFMM 5
#endif

/* Number of fields carried through reconstruction and the Riemann solve. */
#define DF_NV (5 + NDFMM)

/* Primitive/conserved vector slots. */
#define DI_RHO 0
#define DI_UX  1
#define DI_UY  2
#define DI_UZ  3
#define DI_P   4
#define DI_PI  5          /* DI_PI+0..4 = Pi_xx, Pi_yy, Pi_xy, Pi_xz, Pi_yz */

/* Deviatoric pressure component slots within the 5-vector. */
#define PI_XX 0
#define PI_YY 1
#define PI_XY 2
#define PI_XZ 3
#define PI_YZ 4

constant int   TWOTONDIM = 8;
constant int   NSUBGRID  = 1;
constant int   NSUBGRIDP2 = 3;   /* NSUBGRID + 2 */

/* Riemann solver ids (mirror hydro_parameters.f90).  Only LLF (1) and HLL (2)
 * are meaningful for the moment system; read_params.f90 rejects anything else
 * in a DFMM build, and HLL is the default here. */
constant int   DF_SOLVER_LLF = 1;

/* Monatomic ideal gas: the ten-moment trace identity tr P = 3p = 2 rho e
 * fixes gamma.  read_params.f90 rejects any other value in a DFMM build. */
constant float DF_GAMMA   = 5.0f / 3.0f;
constant float DF_GM1     = 2.0f / 3.0f;

/* Largest characteristic of the Wick-closed moment subsystem, CSCOEF in
 * ~/dfmm/py-1d/dfmm/schemes/_common.py.  Bounds the ten-moment value 3, so it
 * is safe at Stage 1 and correct once Q_ijk is evolved at Stage 2. */
constant float DF_CSCOEF  = 3.0f + 2.449489742783178f;   /* 3 + sqrt(6) */

/* ---------------------------------------------------------------------------
 * Global-buffer accessors.  Fortran column-major uold(cell, ivar, oct).
 * --------------------------------------------------------------------------*/
static inline float df_u_get(device const float *u, int oct_1, int ivar_1, int cell_1) {
    return u[(oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)];
}
static inline void df_u_set(device float *u, int oct_1, int ivar_1, int cell_1, float v) {
    u[(oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1)] = v;
}
static inline int df_u_flat(int oct_1, int ivar_1, int cell_1) {
    return (oct_1-1)*(NVAR)*TWOTONDIM + (ivar_1-1)*TWOTONDIM + (cell_1-1);
}

/* ---------------------------------------------------------------------------
 * Threadgroup storage.  Budget (doc/dfmm_3d.md Section 6):
 *   subgrid    DF_NV * 6^3 * 4 B = 8640 B   (+ 64 B refined)
 *   interfaces 6 * DF_NV * 12 * 4 B = 2880 B
 *   total                            11584 B  of the 32768 B Apple limit.
 * --------------------------------------------------------------------------*/
struct df_subgrid_t {
    float v[DF_NV][6][6][6];
    bool  refined[4][4][4];
};
struct df_ix_t { float v[DF_NV][3][2][2]; };
struct df_iy_t { float v[DF_NV][2][3][2]; };
struct df_iz_t { float v[DF_NV][2][2][3]; };

/* ===========================================================================
 * Symmetric-tensor helpers
 *
 * sym6 packing order is (xx, yy, zz, xy, xz, yz).
 * ========================================================================= */

/* Expand the five stored deviatoric components to full sym6, using
 * Pi_zz = -(Pi_xx + Pi_yy). */
static inline void df_pi_to_sym6(thread const float *pi, thread float *s) {
    s[0] = pi[PI_XX];
    s[1] = pi[PI_YY];
    s[2] = -(pi[PI_XX] + pi[PI_YY]);
    s[3] = pi[PI_XY];
    s[4] = pi[PI_XZ];
    s[5] = pi[PI_YZ];
}

/* Full pressure tensor P = p I + Pi, as sym6. */
static inline void df_P_sym6(float p, thread const float *pi, thread float *P) {
    df_pi_to_sym6(pi, P);
    P[0] += p; P[1] += p; P[2] += p;
}

/* Normal component P_nn for n = x, y, z (idim = 0, 1, 2). */
static inline float df_P_nn(float p, thread const float *pi, int idim) {
    if (idim == 0) return p + pi[PI_XX];
    if (idim == 1) return p + pi[PI_YY];
    return p - (pi[PI_XX] + pi[PI_YY]);
}

/* Smallest eigenvalue of a symmetric 3x3 given as sym6 (trigonometric form,
 * branch-free apart from the degenerate-diagonal case). */
static inline float df_lam_min_sym6(thread const float *A) {
    float p1 = A[3]*A[3] + A[4]*A[4] + A[5]*A[5];
    float q  = (A[0] + A[1] + A[2]) / 3.0f;
    if (p1 <= 0.0f) return min(A[0], min(A[1], A[2]));
    float d0 = A[0] - q, d1 = A[1] - q, d2 = A[2] - q;
    float p2 = d0*d0 + d1*d1 + d2*d2 + 2.0f*p1;
    float pp = sqrt(p2 / 6.0f);
    float ipp = 1.0f / max(pp, 1e-30f);
    /* B = (A - q I)/pp ; r = det(B)/2 */
    float b0 = d0*ipp,   b1 = d1*ipp,   b2 = d2*ipp;
    float b3 = A[3]*ipp, b4 = A[4]*ipp, b5 = A[5]*ipp;
    float detB = b0*(b1*b2 - b5*b5) - b3*(b3*b2 - b5*b4) + b4*(b3*b5 - b1*b4);
    float r = clamp(0.5f*detB, -1.0f, 1.0f);
    float phi = acos(r) / 3.0f;
    /* eig1 is the largest, eig3 the smallest */
    float eig1 = q + 2.0f*pp*cos(phi);
    float eig3 = q + 2.0f*pp*cos(phi + 2.0943951023931953f);  /* +2pi/3 */
    return min(eig3, min(eig1, 3.0f*q - eig1 - eig3));
}

/* ---------------------------------------------------------------------------
 * Frame rotation for y/z faces.
 *
 * The Riemann solve is written once for an x-normal face; y and z faces are
 * rotated into that frame and the resulting flux rotated back.  The velocity
 * rotation matches the baseline hydro kernel exactly; the deviatoric pressure
 * rotates as the rank-2 tensor it is.
 *
 *   y-face: (x',y',z') = (y,z,x)
 *   z-face: (x',y',z') = (z,x,y)
 * --------------------------------------------------------------------------*/
static inline void df_rotate_fwd(thread float *w, int idim) {
    if (idim == 0) return;
    float ux = w[DI_UX], uy = w[DI_UY], uz = w[DI_UZ];
    float a = w[DI_PI+PI_XX], b = w[DI_PI+PI_YY];
    float c = w[DI_PI+PI_XY], d = w[DI_PI+PI_XZ], e = w[DI_PI+PI_YZ];
    float zz = -(a + b);
    if (idim == 1) {
        w[DI_UX] = uy; w[DI_UY] = uz; w[DI_UZ] = ux;
        w[DI_PI+PI_XX] = b;    /* Pi'_x'x' = Pi_yy */
        w[DI_PI+PI_YY] = zz;   /* Pi'_y'y' = Pi_zz */
        w[DI_PI+PI_XY] = e;    /* Pi'_x'y' = Pi_yz */
        w[DI_PI+PI_XZ] = c;    /* Pi'_x'z' = Pi_yx */
        w[DI_PI+PI_YZ] = d;    /* Pi'_y'z' = Pi_zx */
    } else {
        w[DI_UX] = uz; w[DI_UY] = ux; w[DI_UZ] = uy;
        w[DI_PI+PI_XX] = zz;   /* Pi'_x'x' = Pi_zz */
        w[DI_PI+PI_YY] = a;    /* Pi'_y'y' = Pi_xx */
        w[DI_PI+PI_XY] = d;    /* Pi'_x'y' = Pi_zx */
        w[DI_PI+PI_XZ] = e;    /* Pi'_x'z' = Pi_zy */
        w[DI_PI+PI_YZ] = c;    /* Pi'_y'z' = Pi_xy */
    }
}

/* Rotate a flux vector from the x-normal frame back to the lab frame.  The
 * momentum slots carry momentum-flux components and the Pi slots carry the
 * fluxes of the rotated Pi components, so both invert with the transpose. */
static inline void df_rotate_flux_back(thread float *w, int idim) {
    if (idim == 0) return;
    float fx = w[DI_UX], fy = w[DI_UY], fz = w[DI_UZ];
    float g0 = w[DI_PI+PI_XX], g1 = w[DI_PI+PI_YY];
    float g2 = w[DI_PI+PI_XY], g3 = w[DI_PI+PI_XZ], g4 = w[DI_PI+PI_YZ];
    if (idim == 1) {
        /* mx'->my, my'->mz, mz'->mx */
        w[DI_UY] = fx; w[DI_UZ] = fy; w[DI_UX] = fz;
        /* F(Pi_yy)=g0, F(Pi_zz)=g1 => F(Pi_xx) = -(g0+g1) */
        w[DI_PI+PI_YY] = g0;
        w[DI_PI+PI_XX] = -(g0 + g1);
        w[DI_PI+PI_YZ] = g2;
        w[DI_PI+PI_XY] = g3;
        w[DI_PI+PI_XZ] = g4;
    } else {
        /* mx'->mz, my'->mx, mz'->my */
        w[DI_UZ] = fx; w[DI_UX] = fy; w[DI_UY] = fz;
        /* F(Pi_zz)=g0, F(Pi_xx)=g1 => F(Pi_yy) = -(g0+g1) */
        w[DI_PI+PI_XX] = g1;
        w[DI_PI+PI_YY] = -(g0 + g1);
        w[DI_PI+PI_XZ] = g2;
        w[DI_PI+PI_YZ] = g3;
        w[DI_PI+PI_XY] = g4;
    }
}

/* ===========================================================================
 * Conserved <-> primitive
 * ========================================================================= */
static inline void df_cons_to_prim(thread const float *c, thread float *w,
                                    float smallr, float smallc2) {
    float rho = max(c[DI_RHO], smallr);
    w[DI_RHO] = rho;
    w[DI_UX]  = c[DI_UX] / rho;
    w[DI_UY]  = c[DI_UY] / rho;
    w[DI_UZ]  = c[DI_UZ] / rho;
    float ekin = 0.5f*(c[DI_UX]*c[DI_UX] + c[DI_UY]*c[DI_UY] + c[DI_UZ]*c[DI_UZ]) / rho;
    w[DI_P]   = max(DF_GM1*(c[DI_P] - ekin), rho*smallc2/DF_GAMMA);
    for (int m = 0; m < NDFMM; m++) w[DI_PI+m] = c[DI_PI+m];
}

static inline void df_prim_to_cons(thread const float *w, thread float *c) {
    float rho = w[DI_RHO];
    c[DI_RHO] = rho;
    c[DI_UX]  = rho*w[DI_UX];
    c[DI_UY]  = rho*w[DI_UY];
    c[DI_UZ]  = rho*w[DI_UZ];
    c[DI_P]   = w[DI_P]/DF_GM1
              + 0.5f*rho*(w[DI_UX]*w[DI_UX] + w[DI_UY]*w[DI_UY] + w[DI_UZ]*w[DI_UZ]);
    for (int m = 0; m < NDFMM; m++) c[DI_PI+m] = w[DI_PI+m];
}

/* ===========================================================================
 * Physical flux in the x-normal frame (doc/dfmm_3d.md Section 3)
 *
 *   F[rho]   = rho ux
 *   F[mx]    = rho ux ux + p + Pi_xx
 *   F[my]    = rho uy ux + Pi_xy
 *   F[mz]    = rho uz ux + Pi_xz
 *   F[E]     = ux (E + p) + (ux Pi_xx + uy Pi_xy + uz Pi_xz)      [+ q_x at Stage 2]
 *   F[Pi_ij] = ux Pi_ij                                          [+ Q_ijx - (2/3)d_ij q_x]
 * ========================================================================= */
static inline void df_phys_flux(thread const float *w, thread float *F) {
    float rho = w[DI_RHO], ux = w[DI_UX], uy = w[DI_UY], uz = w[DI_UZ], p = w[DI_P];
    float pxx = w[DI_PI+PI_XX], pxy = w[DI_PI+PI_XY], pxz = w[DI_PI+PI_XZ];
    float E = p/DF_GM1 + 0.5f*rho*(ux*ux + uy*uy + uz*uz);

    F[DI_RHO] = rho*ux;
    F[DI_UX]  = rho*ux*ux + p + pxx;
    F[DI_UY]  = rho*uy*ux + pxy;
    F[DI_UZ]  = rho*uz*ux + pxz;
    F[DI_P]   = ux*(E + p) + (ux*pxx + uy*pxy + uz*pxz);
    for (int m = 0; m < NDFMM; m++) F[DI_PI+m] = ux*w[DI_PI+m];
}

/* ===========================================================================
 * HLL / LLF flux for the DF_NV-field state.
 *
 * Wave speeds use the anisotropic normal sound speed
 *   c = sqrt(CSCOEF * P_xx / rho),  P_xx = p + Pi_xx,
 * so pressure anisotropy feeds the numerical viscosity and (via cmpdt) the
 * timestep, as intended.
 * ========================================================================= */
static inline void df_hll_flux(thread float *wl, thread float *wr,
                                float smallr, float smallc2, bool llf,
                                thread float *F) {
    wl[DI_RHO] = max(wl[DI_RHO], smallr);
    wr[DI_RHO] = max(wr[DI_RHO], smallr);
    float smallp = smallc2 / DF_GAMMA;
    wl[DI_P] = max(wl[DI_P], smallp*wl[DI_RHO]);
    wr[DI_P] = max(wr[DI_P], smallp*wr[DI_RHO]);

    float Pl = max(wl[DI_P] + wl[DI_PI+PI_XX], smallp*wl[DI_RHO]);
    float Pr = max(wr[DI_P] + wr[DI_PI+PI_XX], smallp*wr[DI_RHO]);
    float cl = sqrt(DF_CSCOEF*Pl/wl[DI_RHO]);
    float cr = sqrt(DF_CSCOEF*Pr/wr[DI_RHO]);

    float SL, SR;
    if (llf) {
        float a = max(abs(wl[DI_UX]) + cl, abs(wr[DI_UX]) + cr);
        SL = -a; SR = a;
    } else {
        SL = min(wl[DI_UX] - cl, wr[DI_UX] - cr);
        SR = max(wl[DI_UX] + cl, wr[DI_UX] + cr);
    }

    float FL[DF_NV], FR[DF_NV], UL[DF_NV], UR[DF_NV];
    df_phys_flux(wl, FL);
    df_phys_flux(wr, FR);

    if (SL >= 0.0f) { for (int k = 0; k < DF_NV; k++) F[k] = FL[k]; return; }
    if (SR <= 0.0f) { for (int k = 0; k < DF_NV; k++) F[k] = FR[k]; return; }

    df_prim_to_cons(wl, UL);
    df_prim_to_cons(wr, UR);
    float inv = 1.0f/(SR - SL);
    for (int k = 0; k < DF_NV; k++)
        F[k] = (SR*FL[k] - SL*FR[k] + SL*SR*(UR[k] - UL[k]))*inv;
}

/* ===========================================================================
 * Slopes and MUSCL-Hancock trace
 * ========================================================================= */
static inline float df_slope(float left, float middle, float right, int slope) {
    float sl = middle - left;
    float sr = right  - middle;
    float sc = 0.5f*(sl + sr);
    float f  = float(slope);
    if (sl*sr <= 0.0f) return 0.0f;
    if (sl > 0.0f) return min(f*min(sl, sr), sc);
    else           return max(f*max(sl, sr), sc);
}

/* Load the 6x6x6 primitive stencil from uold. */
static void df_load_subgrid(
    device const oct_t *grid, device const float *uold, device const int *nbor,
    constant float *constant_gravity, device const float *f,
    int head_idx, int block_idx, int thread_idx,
    float smallr, float smallc2, float dt, uint threads_per_tg,
    threadgroup df_subgrid_t &ls)
{
    const int work_size  = 2*NSUBGRID + 4;
    const int total_work = work_size*work_size*work_size;

    for (int wi = thread_idx; wi < total_work; wi += int(threads_per_tg)) {
        int i_sg, j_sg, k_sg;
        index_1Dto3D(wi/8, work_size/2, work_size/2, i_sg, j_sg, k_sg);
        int subgrid_idx = head_idx + block_idx;
        int ind_nbor    = wi/8 + 1;
        int source_idx  = nbor_get(nbor, subgrid_idx, ind_nbor);
        int cell_idx    = wi%8 + 1;

        int ib, jb, kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2*i_sg, j = jb + 2*j_sg, k = kb + 2*k_sg;

        float c[DF_NV], w[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++)
            c[iv] = df_u_get(uold, source_idx, iv + 1, cell_idx);
        df_cons_to_prim(c, w, smallr, smallc2);

#ifdef GRAV
        w[DI_UX] += f[(source_idx-1)*3*8 + 0*8 + (cell_idx-1)]*0.5f*dt;
        w[DI_UY] += f[(source_idx-1)*3*8 + 1*8 + (cell_idx-1)]*0.5f*dt;
        w[DI_UZ] += f[(source_idx-1)*3*8 + 2*8 + (cell_idx-1)]*0.5f*dt;
#else
        w[DI_UX] += constant_gravity[0]*0.5f*dt;
        w[DI_UY] += constant_gravity[1]*0.5f*dt;
        w[DI_UZ] += constant_gravity[2]*0.5f*dt;
#endif
        for (int iv = 0; iv < DF_NV; iv++) ls.v[iv][i][j][k] = w[iv];

        if (i >= 1 && i <= 2*NSUBGRID+2 && j >= 1 && j <= 2*NSUBGRID+2
                                        && k >= 1 && k <= 2*NSUBGRID+2)
            ls.refined[i-1][j-1][k-1] = (grid[source_idx-1].refined[cell_idx-1] != 0);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* MUSCL-Hancock predictor.
 *
 * Slope convention follows the baseline hydro kernel: sx[] holds the *half*
 * slope (dx/2) d_x w, and the caller passes dtdx = dt/dx, so `w + dtdx*src`
 * is the half-step predicted state.  Every dfmm source term below therefore
 * also uses the bare half-slopes, which keeps the scaling uniform.
 *
 * Primitive sources are the linearised Euler ones plus, from
 * doc/dfmm_3d.md Sections 3-4:
 *   velocity:  -(1/rho) d_k Pi_ik           (from the momentum flux term Pi_ik)
 *   pressure:  -(2/3) Pi_kl d_l u_k         (from the energy flux term u_i Pi_ik)
 *   Pi_ij:     -u.grad Pi_ij - Pi_ij div u  (transport, density-like)
 *              -2 p S0_ij - [Pi_ik G_jk + Pi_jk G_ik]^dev   (strain production)
 *
 * The first two are *flux*-derived and are handled conservatively by the
 * Riemann solve; they appear here only because the predictor works in
 * primitive variables.  The strain production is a genuine source and is
 * applied over the full step by dfmm_source_kernel; including it at half
 * weight here only sharpens the interface states.
 *
 * BGK relaxation is deliberately absent from the predictor: it is applied as
 * an exact exponential map, mirroring the 1D reference.
 */
static void df_trace(
    threadgroup df_subgrid_t &ls, int thread_idx, uint threads_per_tg,
    threadgroup df_ix_t &lx, threadgroup df_ix_t &rx,
    threadgroup df_iy_t &ly, threadgroup df_iy_t &ry,
    threadgroup df_iz_t &lz, threadgroup df_iz_t &rz,
    float smallr, float smallc2, float dtdx, int slope, int source_on)
{
    const int work_size  = 2*NSUBGRID + 2;
    const int total_work = work_size*work_size*work_size;
    float smallp = smallr*smallc2;

    for (int wi = thread_idx; wi < total_work; wi += int(threads_per_tg)) {
        int i, j, k;
        index_1Dto3D(wi, work_size, work_size, i, j, k);
        i += 1; j += 1; k += 1;

        float w[DF_NV], sx[DF_NV], sy[DF_NV], sz[DF_NV], src[DF_NV], q[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++) {
            w [iv] = ls.v[iv][i][j][k];
            sx[iv] = 0.5f*df_slope(ls.v[iv][i-1][j][k], w[iv], ls.v[iv][i+1][j][k], slope);
            sy[iv] = 0.5f*df_slope(ls.v[iv][i][j-1][k], w[iv], ls.v[iv][i][j+1][k], slope);
            sz[iv] = 0.5f*df_slope(ls.v[iv][i][j][k-1], w[iv], ls.v[iv][i][j][k+1], slope);
        }

        float rho = w[DI_RHO], ux = w[DI_UX], uy = w[DI_UY], uz = w[DI_UZ], p = w[DI_P];
        float irho = 1.0f/rho;
        float divs = sx[DI_UX] + sy[DI_UY] + sz[DI_UZ];

        src[DI_RHO] = -(ux*sx[DI_RHO] + uy*sy[DI_RHO] + uz*sz[DI_RHO]) - divs*rho;
        src[DI_UX]  = -(ux*sx[DI_UX] + uy*sy[DI_UX] + uz*sz[DI_UX]) - sx[DI_P]*irho;
        src[DI_UY]  = -(ux*sx[DI_UY] + uy*sy[DI_UY] + uz*sz[DI_UY]) - sy[DI_P]*irho;
        src[DI_UZ]  = -(ux*sx[DI_UZ] + uy*sy[DI_UZ] + uz*sz[DI_UZ]) - sz[DI_P]*irho;
        src[DI_P]   = -(ux*sx[DI_P] + uy*sy[DI_P] + uz*sz[DI_P]) - divs*DF_GAMMA*p;
        /* Pi is density-like: transport contributes u.grad Pi + Pi div u */
        for (int m = 0; m < NDFMM; m++)
            src[DI_PI+m] = -(ux*sx[DI_PI+m] + uy*sy[DI_PI+m] + uz*sz[DI_PI+m])
                           - divs*w[DI_PI+m];

        if (source_on != 0) {
            /* G[a][b] = half-slope of u_a in direction b, i.e. (dx/2) d_b u_a */
            float G[3][3] = {
                {sx[DI_UX], sy[DI_UX], sz[DI_UX]},
                {sx[DI_UY], sy[DI_UY], sz[DI_UY]},
                {sx[DI_UZ], sy[DI_UZ], sz[DI_UZ]}};
            float divu = G[0][0] + G[1][1] + G[2][2];

            float pis[6]; df_pi_to_sym6(&w[DI_PI], pis);
            float PI[3][3] = {{pis[0], pis[3], pis[4]},
                              {pis[3], pis[1], pis[5]},
                              {pis[4], pis[5], pis[2]}};

            /* Pi_kl d_l u_k  -> Pi_ab G_ab */
            float pig = 0.0f;
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) pig += PI[a][b]*G[a][b];

            /* d_k Pi_ik from the half-slopes of the five stored components */
            float dpi[3];
            dpi[0] = sx[DI_PI+PI_XX] + sy[DI_PI+PI_XY] + sz[DI_PI+PI_XZ];
            dpi[1] = sx[DI_PI+PI_XY] + sy[DI_PI+PI_YY] + sz[DI_PI+PI_YZ];
            dpi[2] = sx[DI_PI+PI_XZ] + sy[DI_PI+PI_YZ]
                   - (sz[DI_PI+PI_XX] + sz[DI_PI+PI_YY]);

            src[DI_UX] -= dpi[0]*irho;
            src[DI_UY] -= dpi[1]*irho;
            src[DI_UZ] -= dpi[2]*irho;
            src[DI_P]  -= (2.0f/3.0f)*pig;

            /* The strain production -2 p S0 - [Pi G]^dev is deliberately NOT
             * added here.  It is stiff: it balances against -Pi/tau_Pi, and
             * the two must be integrated together to stay
             * asymptotic-preserving (dfmm_source_kernel does that exactly).
             * Adding it here at an unrelaxed weight would overshoot the
             * interface states whenever dt >> tau_Pi.  The 1D reference makes
             * the same choice: sources follow the flux update, they are not
             * folded into the reconstruction.
             * `divu` is retained above only for the flux-derived terms. */
            (void)divu;
        }

        for (int iv = 0; iv < DF_NV; iv++) q[iv] = w[iv] + dtdx*src[iv];

#define DF_STORE(ARR, I0, I1, I2, S, SGN)                                      \
        {                                                                      \
            for (int iv = 0; iv < DF_NV; iv++)                                 \
                ARR.v[iv][I0][I1][I2] = q[iv] SGN S[iv];                       \
            if (ARR.v[DI_RHO][I0][I1][I2] < smallr)                            \
                ARR.v[DI_RHO][I0][I1][I2] = ls.v[DI_RHO][i][j][k];             \
            if (ARR.v[DI_P][I0][I1][I2] < smallp)                              \
                ARR.v[DI_P][I0][I1][I2] = ls.v[DI_P][i][j][k];                 \
        }

        if (i > 1 && (j > 1 && j < work_size) && (k > 1 && k < work_size))
            DF_STORE(rx, i-2, j-2, k-2, sx, -)
        if (i < work_size && (j > 1 && j < work_size) && (k > 1 && k < work_size))
            DF_STORE(lx, i-1, j-2, k-2, sx, +)
        if ((i > 1 && i < work_size) && j > 1 && (k > 1 && k < work_size))
            DF_STORE(ry, i-2, j-2, k-2, sy, -)
        if ((i > 1 && i < work_size) && j < work_size && (k > 1 && k < work_size))
            DF_STORE(ly, i-2, j-1, k-2, sy, +)
        if ((i > 1 && i < work_size) && (j > 1 && j < work_size) && k > 1)
            DF_STORE(rz, i-2, j-2, k-2, sz, -)
        if ((i > 1 && i < work_size) && (j > 1 && j < work_size) && k < work_size)
            DF_STORE(lz, i-2, j-2, k-1, sz, +)
#undef DF_STORE
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * Riemann driver — one x-normal solver, y/z rotated in and out
 * ========================================================================= */
static void df_riemann_driver(
    threadgroup df_ix_t &lx, threadgroup df_ix_t &rx,
    threadgroup df_iy_t &ly, threadgroup df_iy_t &ry,
    threadgroup df_iz_t &lz, threadgroup df_iz_t &rz,
    int thread_idx, uint threads_per_tg,
    float smallr, float smallc2, int riemann)
{
    const int ias = (2*NSUBGRID+1)*(2*NSUBGRID)*(2*NSUBGRID);
    bool llf = (riemann == DF_SOLVER_LLF);

    for (int wi = thread_idx; wi < 3*ias; wi += int(threads_per_tg)) {
        int i, j, k, idim;
        float L[DF_NV], R[DF_NV], F[DF_NV];

        if (wi < ias) {
            idim = 0;
            index_1Dto3D(wi, 2*NSUBGRID+1, 2*NSUBGRID, i, j, k);
            for (int iv = 0; iv < DF_NV; iv++) { L[iv] = lx.v[iv][i][j][k]; R[iv] = rx.v[iv][i][j][k]; }
        } else if (wi < 2*ias) {
            idim = 1;
            index_1Dto3D(wi - ias, 2*NSUBGRID, 2*NSUBGRID+1, i, j, k);
            for (int iv = 0; iv < DF_NV; iv++) { L[iv] = ly.v[iv][i][j][k]; R[iv] = ry.v[iv][i][j][k]; }
        } else {
            idim = 2;
            index_1Dto3D(wi - 2*ias, 2*NSUBGRID, 2*NSUBGRID, i, j, k);
            for (int iv = 0; iv < DF_NV; iv++) { L[iv] = lz.v[iv][i][j][k]; R[iv] = rz.v[iv][i][j][k]; }
        }

        df_rotate_fwd(L, idim);
        df_rotate_fwd(R, idim);
        df_hll_flux(L, R, smallr, smallc2, llf, F);
        df_rotate_flux_back(F, idim);

        if (idim == 0)      { for (int iv = 0; iv < DF_NV; iv++) lx.v[iv][i][j][k] = F[iv]; }
        else if (idim == 1) { for (int iv = 0; iv < DF_NV; iv++) ly.v[iv][i][j][k] = F[iv]; }
        else                { for (int iv = 0; iv < DF_NV; iv++) lz.v[iv][i][j][k] = F[iv]; }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* ===========================================================================
 * Fine-face flux zeroing, conservative update, coarse correction
 * ========================================================================= */
static void df_zero_fine_fluxes(
    threadgroup const df_subgrid_t &ls, int thread_idx, uint threads_per_tg,
    threadgroup df_ix_t &fx, threadgroup df_iy_t &fy, threadgroup df_iz_t &fz)
{
    const int ias = (2*NSUBGRID+1)*(2*NSUBGRID)*(2*NSUBGRID);
    for (int wi = thread_idx; wi < 3*ias; wi += int(threads_per_tg)) {
        int i, j, k;
        if (wi < ias) {
            index_1Dto3D(wi, 2*NSUBGRID+1, 2*NSUBGRID, i, j, k);
            if (ls.refined[i][j+1][k+1] || ls.refined[i+1][j+1][k+1])
                for (int iv = 0; iv < DF_NV; iv++) fx.v[iv][i][j][k] = 0.0f;
        } else if (wi < 2*ias) {
            index_1Dto3D(wi - ias, 2*NSUBGRID, 2*NSUBGRID+1, i, j, k);
            if (ls.refined[i+1][j][k+1] || ls.refined[i+1][j+1][k+1])
                for (int iv = 0; iv < DF_NV; iv++) fy.v[iv][i][j][k] = 0.0f;
        } else {
            index_1Dto3D(wi - 2*ias, 2*NSUBGRID, 2*NSUBGRID, i, j, k);
            if (ls.refined[i+1][j+1][k] || ls.refined[i+1][j+1][k+1])
                for (int iv = 0; iv < DF_NV; iv++) fz.v[iv][i][j][k] = 0.0f;
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
}

/* Conservative update — transcription of conservative_update in hydro.metal,
 * looped over DF_NV fields. */
static void df_conservative_update(
    device float *unew, device const int *nbor,
    threadgroup const df_ix_t &fx, threadgroup const df_iy_t &fy,
    threadgroup const df_iz_t &fz,
    int head_idx, int block_idx, int thread_idx, uint threads_per_tg, float dtdx)
{
    const int work_size  = 2*NSUBGRID;
    const int total_work = work_size*work_size*work_size;

    for (int work_idx = thread_idx; work_idx < total_work; work_idx += int(threads_per_tg)) {
        int subgrid_idx = head_idx + block_idx;

        int i_sg, j_sg, k_sg;
        index_1Dto3D(work_idx/8, work_size/2, work_size/2, i_sg, j_sg, k_sg);
        i_sg += 1; j_sg += 1; k_sg += 1;

        int ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
        int oct_idx  = nbor_get(nbor, subgrid_idx, ind_nbor);

        int cell_idx = work_idx%8 + 1;
        int i, j, k;
        index_1Dto3D(cell_idx - 1, 2, 2, i, j, k);
        i += 2*(i_sg - 1);
        j += 2*(j_sg - 1);
        k += 2*(k_sg - 1);

        for (int iv = 0; iv < DF_NV; iv++) {
            float upd = (fx.v[iv][i][j][k] - fx.v[iv][i+1][j  ][k  ])*dtdx
                      + (fy.v[iv][i][j][k] - fy.v[iv][i  ][j+1][k  ])*dtdx
                      + (fz.v[iv][i][j][k] - fz.v[iv][i  ][j  ][k+1])*dtdx;
            df_u_set(unew, oct_idx, iv + 1, cell_idx,
                     df_u_get(unew, oct_idx, iv + 1, cell_idx) + upd);
        }
    }
}

/* Coarse-level flux correction — transcription of coarse_cell_update in
 * hydro.metal, looped over DF_NV fields.  Only ghost-cache neighbours
 * (source_idx > ngridmax) carry a coarse father that needs correcting; the
 * coarse cell is located from the *source* oct's ckey relative to its father,
 * not from the current oct. */
static void df_coarse_cell_update(
    device float *unew, device const int *nbor, device const oct_t *grid,
    device const int *father,
    threadgroup const df_ix_t &fx, threadgroup const df_iy_t &fy,
    threadgroup const df_iz_t &fz,
    int head_idx, int block_idx, int ngridmax, int thread_idx, float cfs)
{
    const int ias = NSUBGRID*NSUBGRID;
    if (thread_idx >= 6*ias) return;

    int subgrid_idx = head_idx + block_idx;
    int face     = thread_idx/ias;
    int work_idx = thread_idx%ias;
    int j_raw = work_idx%NSUBGRID;
    int k_raw = work_idx/NSUBGRID;

    float acc[DF_NV];
    for (int iv = 0; iv < DF_NV; iv++) acc[iv] = 0.0f;
    int ind_nbor = 0;

    if (face == 0) {
        int j_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] -= fx.v[iv][0][j][k]*cfs;
        ind_nbor = 1 + 0 + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 1) {
        int j_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] += fx.v[iv][2*NSUBGRID][j][k]*cfs;
        ind_nbor = 1 + (NSUBGRID+1) + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 2) {
        int i_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] -= fy.v[iv][i][0][k]*cfs;
        ind_nbor = 1 + i_sg + NSUBGRIDP2*0 + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 3) {
        int i_sg = j_raw + 1, k_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int k = 2*k_sg-2; k <= 2*k_sg-1; k++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] += fy.v[iv][i][2*NSUBGRID][k]*cfs;
        ind_nbor = 1 + i_sg + NSUBGRIDP2*(NSUBGRID+1) + NSUBGRIDP2*NSUBGRIDP2*k_sg;
    } else if (face == 4) {
        int i_sg = j_raw + 1, j_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] -= fz.v[iv][i][j][0]*cfs;
        ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*0;
    } else {
        int i_sg = j_raw + 1, j_sg = k_raw + 1;
        for (int i = 2*i_sg-2; i <= 2*i_sg-1; i++)
            for (int j = 2*j_sg-2; j <= 2*j_sg-1; j++)
                for (int iv = 0; iv < DF_NV; iv++) acc[iv] += fz.v[iv][i][j][2*NSUBGRID]*cfs;
        ind_nbor = 1 + i_sg + NSUBGRIDP2*j_sg + NSUBGRIDP2*NSUBGRIDP2*(NSUBGRID+1);
    }

    int source_idx = nbor_get(nbor, subgrid_idx, ind_nbor);
    if (source_idx > ngridmax) {
        int father_idx = father[source_idx - 1];
        int ic = grid[source_idx-1].ckey[0] - 2*grid[father_idx-1].ckey[0];
        int jc = grid[source_idx-1].ckey[1] - 2*grid[father_idx-1].ckey[1];
        int kc = grid[source_idx-1].ckey[2] - 2*grid[father_idx-1].ckey[2];
        int cell_idx = 1 + ic + 2*jc + 4*kc;
        for (int iv = 0; iv < DF_NV; iv++)
            atomic_add_float((device atomic_uint *)&unew[df_u_flat(father_idx, iv + 1, cell_idx)],
                             acc[iv]);
    }
}

/* ===========================================================================
 * dfmm_integrator_kernel — transport only (Section 3 fluxes)
 * ========================================================================= */
kernel void dfmm_integrator_kernel(
    device const oct_t *grid             [[buffer(0)]],
    device const float *uold             [[buffer(1)]],
    device float       *unew             [[buffer(2)]],
    device const int   *nbor             [[buffer(3)]],
    constant int       &head_idx         [[buffer(4)]],
    constant int       &num_subgrids     [[buffer(5)]],
    constant int       &ngridmax         [[buffer(6)]],
    constant int       &ilevel           [[buffer(7)]],
    constant int       &levelmin         [[buffer(8)]],
    constant int       &levelmax         [[buffer(9)]],
    constant float     &smallr           [[buffer(10)]],
    constant float     &smallc2          [[buffer(11)]],
    constant float     &dt               [[buffer(12)]],
    constant float     &dx               [[buffer(13)]],
    constant int       &slope            [[buffer(14)]],
    constant int       &riemann          [[buffer(15)]],
    constant float     *constant_gravity [[buffer(16)]],
    device const int   *father           [[buffer(17)]],
    device const float *f                [[buffer(18)]],
    constant int       &source_on        [[buffer(19)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_subgrids) return;
    float dtdx = dt/dx;

    threadgroup df_subgrid_t ls;
    threadgroup df_ix_t lx, rx;
    threadgroup df_iy_t ly, ry;
    threadgroup df_iz_t lz, rz;

    df_load_subgrid(grid, uold, nbor, constant_gravity, f,
                    head_idx, int(block_idx), int(thread_idx),
                    smallr, smallc2, dt, threads_per_tg, ls);

    df_trace(ls, int(thread_idx), threads_per_tg, lx, rx, ly, ry, lz, rz,
             smallr, smallc2, dtdx, slope, source_on);

    df_riemann_driver(lx, rx, ly, ry, lz, rz, int(thread_idx), threads_per_tg,
                      smallr, smallc2, riemann);

    if (ilevel < levelmax)
        df_zero_fine_fluxes(ls, int(thread_idx), threads_per_tg, lx, ly, lz);

    df_conservative_update(unew, nbor, lx, ly, lz,
                           head_idx, int(block_idx), int(thread_idx),
                           threads_per_tg, dtdx);

    if (ilevel > levelmin)
        df_coarse_cell_update(unew, nbor, grid, father, lx, ly, lz,
                              head_idx, int(block_idx), ngridmax,
                              int(thread_idx), dtdx/float(TWOTONDIM));
}

/* ===========================================================================
 * dfmm_source_kernel — strain production and exact BGK relaxation
 *
 * Operator order (doc/dfmm_3d.md Section 4): this runs after transport and
 * before uold is overwritten.  Velocity gradients are taken from uold, not
 * unew: set_unew zeroes unew in virtual boundaries, so unew's neighbour data
 * is not valid for a centred difference.  This mirrors the 1D reference, which
 * evaluates du/dx on the pre-update state.
 * ========================================================================= */
kernel void dfmm_source_kernel(
    device const oct_t *grid     [[buffer(0)]],
    device const float *uold     [[buffer(1)]],
    device float       *unew     [[buffer(2)]],
    device const int   *nbor     [[buffer(3)]],
    constant int       &head_idx [[buffer(4)]],
    constant int       &num_octs [[buffer(5)]],
    constant float     &smallr   [[buffer(6)]],
    constant float     &smallc2  [[buffer(7)]],
    constant float     &dt       [[buffer(8)]],
    constant float     &dx       [[buffer(9)]],
    constant float     &tau_pi   [[buffer(10)]],
    constant int       &source_on[[buffer(11)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_octs) return;

    /* 6x6x6 velocity stencil (only three fields: 3*216*4 = 2592 B). */
    threadgroup float uu[3][6][6][6];

    const int work_size  = 6;
    const int total_work = work_size*work_size*work_size;
    for (int wi = int(thread_idx); wi < total_work; wi += int(threads_per_tg)) {
        int i_sg, j_sg, k_sg;
        index_1Dto3D(wi/8, 3, 3, i_sg, j_sg, k_sg);
        int subgrid_idx = head_idx + int(block_idx);
        int source_idx  = nbor_get(nbor, subgrid_idx, wi/8 + 1);
        int cell_idx    = wi%8 + 1;
        int ib, jb, kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2*i_sg, j = jb + 2*j_sg, k = kb + 2*k_sg;

        float rho = max(df_u_get(uold, source_idx, 1, cell_idx), smallr);
        uu[0][i][j][k] = df_u_get(uold, source_idx, 2, cell_idx)/rho;
        uu[1][i][j][k] = df_u_get(uold, source_idx, 3, cell_idx)/rho;
        uu[2][i][j][k] = df_u_get(uold, source_idx, 4, cell_idx)/rho;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int oct_idx = head_idx + int(block_idx);
    float inv2dx = 0.5f/dx;

    /* Exact solution of  dPi/dt = T - Pi/tau  over dt with T frozen:
     *     Pi <- Pi*d + tau*T*(1 - d),      d = exp(-dt/tau).
     * This is asymptotic-preserving.  For dt << tau it reduces to
     * Pi + dt*T; for dt >> tau it gives Pi -> tau*T, which is the
     * Navier-Stokes stress -2 p tau S0.  Naive "add dt*T then multiply by d"
     * instead sends Pi -> 0 for dt >> tau and would destroy exactly the
     * limit the blowup study has to measure against.
     * tau <= 0 means collisionless: no relaxation at all. */
    bool collisional = (tau_pi > 0.0f);
    float decay = 1.0f, tau_omd = dt;
    if (collisional) {
        float xr = dt/tau_pi;
        float omd;
        if (xr < 1.0e-4f) { omd = xr*(1.0f - 0.5f*xr); decay = 1.0f - omd; }
        else              { decay = exp(-xr);          omd = 1.0f - decay; }
        tau_omd = tau_pi*omd;     /* -> dt for small xr, -> tau for large xr */
    }

    for (int cell = int(thread_idx); cell < TWOTONDIM; cell += int(threads_per_tg)) {
        int cell_idx = cell + 1;
        if (grid[oct_idx-1].refined[cell] != 0) continue;

        int ib, jb, kb;
        index_1Dto3D(cell, 2, 2, ib, jb, kb);
        int i = ib + 2, j = jb + 2, k = kb + 2;   /* centre 2x2x2 of the stencil */

        /* Post-transport state of this cell */
        float c[DF_NV], w[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++) c[iv] = df_u_get(unew, oct_idx, iv + 1, cell_idx);
        df_cons_to_prim(c, w, smallr, smallc2);
        float p = w[DI_P];

        float pi[NDFMM];
        for (int m = 0; m < NDFMM; m++) pi[m] = w[DI_PI+m];

        if (source_on != 0) {
            /* G_ij = d_j u_i, centred, from uold */
            float G[3][3];
            for (int a = 0; a < 3; a++) {
                G[a][0] = (uu[a][i+1][j][k] - uu[a][i-1][j][k])*inv2dx;
                G[a][1] = (uu[a][i][j+1][k] - uu[a][i][j-1][k])*inv2dx;
                G[a][2] = (uu[a][i][j][k+1] - uu[a][i][j][k-1])*inv2dx;
            }

            float pis[6]; df_pi_to_sym6(pi, pis);
            float PI[3][3] = {{pis[0], pis[3], pis[4]},
                              {pis[3], pis[1], pis[5]},
                              {pis[4], pis[5], pis[2]}};

            float T[3][3];
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
                float S = 0.5f*(G[a][b] + G[b][a]);
                float PG = 0.0f;
                for (int cc = 0; cc < 3; cc++) PG += PI[a][cc]*G[b][cc] + PI[b][cc]*G[a][cc];
                T[a][b] = -2.0f*p*S - PG;
            }
            float trT = (T[0][0] + T[1][1] + T[2][2])/3.0f;
            T[0][0] -= trT; T[1][1] -= trT; T[2][2] -= trT;

            pi[PI_XX] = pi[PI_XX]*decay + tau_omd*T[0][0];
            pi[PI_YY] = pi[PI_YY]*decay + tau_omd*T[1][1];
            pi[PI_XY] = pi[PI_XY]*decay + tau_omd*T[0][1];
            pi[PI_XZ] = pi[PI_XZ]*decay + tau_omd*T[0][2];
            pi[PI_YZ] = pi[PI_YZ]*decay + tau_omd*T[1][2];
        } else {
            /* No strain production: pure relaxation toward isotropy. */
            for (int m = 0; m < NDFMM; m++) pi[m] *= decay;
        }

        for (int m = 0; m < NDFMM; m++) df_u_set(unew, oct_idx, 6 + m, cell_idx, pi[m]);
    }
}

/* ===========================================================================
 * dfmm_cmpdt_kernel — anisotropic timestep and conserved-quantity sums
 *
 *   ctot = sum_d ( |u_d| + sqrt(CSCOEF P_dd / rho) )
 *
 * so anisotropy shortens dt directly, which is the intended behaviour.
 * ========================================================================= */
kernel void dfmm_cmpdt_kernel(
    device const oct_t   *grid             [[buffer(0)]],
    device const float   *uold             [[buffer(1)]],
    device atomic_uint   *data_buf         [[buffer(2)]],
    constant int         &head_idx         [[buffer(3)]],
    constant int         &num_octs         [[buffer(4)]],
    constant float       &dx               [[buffer(5)]],
    constant float       &smallr           [[buffer(6)]],
    constant float       &smallc2          [[buffer(7)]],
    constant float       &courant_factor   [[buffer(8)]],
    constant float       *constant_gravity [[buffer(9)]],
    device const float   *f                [[buffer(10)]],
    uint tid  [[thread_position_in_threadgroup]],
    uint bid  [[threadgroup_position_in_grid]],
    uint lane [[thread_index_in_simdgroup]],
    uint sg_idx [[simdgroup_index_in_threadgroup]])
{
    float mass_loc = 0.0f, etot_loc = 0.0f, eint_loc = 0.0f, eani_loc = 0.0f;
    float dt_loc = HUGE_VALF;

    int oct_offset = int(bid*1024 + tid)/TWOTONDIM;
    if (oct_offset < num_octs) {
        int oct_idx  = head_idx + oct_offset;
        int cell_idx = int(tid)%TWOTONDIM + 1;
        if (grid[oct_idx-1].refined[cell_idx-1] == 0) {
            float c[DF_NV], w[DF_NV];
            for (int iv = 0; iv < DF_NV; iv++) c[iv] = df_u_get(uold, oct_idx, iv + 1, cell_idx);
            df_cons_to_prim(c, w, smallr, smallc2);

            float vol = dx*dx*dx;
            mass_loc = w[DI_RHO]*vol;
            etot_loc = c[DI_P]*vol;
            eint_loc = w[DI_P]/DF_GM1*vol;

            /* Anisotropy energy measure: |Pi|_F / p, volume weighted. */
            float pis[6]; df_pi_to_sym6(&w[DI_PI], pis);
            float fro = pis[0]*pis[0] + pis[1]*pis[1] + pis[2]*pis[2]
                      + 2.0f*(pis[3]*pis[3] + pis[4]*pis[4] + pis[5]*pis[5]);
            eani_loc = sqrt(fro)*vol;

            float ctot = 0.0f;
            for (int d = 0; d < 3; d++) {
                float Pnn = max(df_P_nn(w[DI_P], &w[DI_PI], d), w[DI_RHO]*smallc2/DF_GAMMA);
                ctot += abs(w[DI_UX + d]) + sqrt(DF_CSCOEF*Pnn/w[DI_RHO]);
            }

            float grav;
#ifdef GRAV
            grav = abs(f[(oct_idx-1)*3*8 + 0*8 + (cell_idx-1)])
                 + abs(f[(oct_idx-1)*3*8 + 1*8 + (cell_idx-1)])
                 + abs(f[(oct_idx-1)*3*8 + 2*8 + (cell_idx-1)]);
#else
            grav = abs(constant_gravity[0]) + abs(constant_gravity[1]) + abs(constant_gravity[2]);
#endif
            grav = max(grav*dx/(ctot*ctot), 0.0001f);
            dt_loc = dx/ctot*(sqrt(1.0f + 2.0f*courant_factor*grav) - 1.0f)/grav;
        }
    }

    threadgroup float tg_dt[32], tg_mass[32], tg_etot[32], tg_eint[32], tg_eani[32];
    float wdt = simd_min(dt_loc);
    float wm  = simd_sum(mass_loc);
    float we  = simd_sum(etot_loc);
    float wi  = simd_sum(eint_loc);
    float wa  = simd_sum(eani_loc);
    if (lane == 0) {
        tg_dt[sg_idx] = wdt; tg_mass[sg_idx] = wm;
        tg_etot[sg_idx] = we; tg_eint[sg_idx] = wi; tg_eani[sg_idx] = wa;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (sg_idx == 0) {
        float bdt = simd_min(tg_dt[lane]);
        float bm  = simd_sum(tg_mass[lane]);
        float be  = simd_sum(tg_etot[lane]);
        float bi  = simd_sum(tg_eint[lane]);
        float ba  = simd_sum(tg_eani[lane]);
        if (lane == 0) {
            atomic_add_float     (&data_buf[0], bm);
            atomic_add_float     (&data_buf[1], be);
            atomic_add_float     (&data_buf[2], bi);
            atomic_add_float     (&data_buf[3], ba);
            atomic_min_float_bits(&data_buf[4], bdt);
        }
    }
}

/* ===========================================================================
 * dfmm_diag_kernel — closure-quality and realizability diagnostics
 *
 * Reported per level (doc/dfmm_3d.md Section 0):
 *   diag[0]  min over cells of lam_min(P)/p          (realizability margin)
 *   diag[1]  max over cells of ||Pi - Pi_NS||_F / p  (Navier-Stokes deviation)
 *   diag[2]  max over cells of ||Pi||_F / p          (anisotropy amplitude)
 *   diag[3]  count of cells with lam_min(P) < 0      (as a float)
 *
 * A minimum is used for the realizability margin deliberately: a maximum of
 * the slack is structurally blind to cells sitting at the cone boundary.
 * ========================================================================= */
kernel void dfmm_diag_kernel(
    device const oct_t *grid     [[buffer(0)]],
    device const float *uold     [[buffer(1)]],
    device const int   *nbor     [[buffer(2)]],
    device atomic_uint *diag     [[buffer(3)]],
    constant int       &head_idx [[buffer(4)]],
    constant int       &num_octs [[buffer(5)]],
    constant float     &smallr   [[buffer(6)]],
    constant float     &smallc2  [[buffer(7)]],
    constant float     &dx       [[buffer(8)]],
    constant float     &tau_pi   [[buffer(9)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_octs) return;

    threadgroup float uu[3][6][6][6];
    for (int wi = int(thread_idx); wi < 216; wi += int(threads_per_tg)) {
        int i_sg, j_sg, k_sg;
        index_1Dto3D(wi/8, 3, 3, i_sg, j_sg, k_sg);
        int subgrid_idx = head_idx + int(block_idx);
        int source_idx  = nbor_get(nbor, subgrid_idx, wi/8 + 1);
        int cell_idx    = wi%8 + 1;
        int ib, jb, kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2*i_sg, j = jb + 2*j_sg, k = kb + 2*k_sg;
        float rho = max(df_u_get(uold, source_idx, 1, cell_idx), smallr);
        uu[0][i][j][k] = df_u_get(uold, source_idx, 2, cell_idx)/rho;
        uu[1][i][j][k] = df_u_get(uold, source_idx, 3, cell_idx)/rho;
        uu[2][i][j][k] = df_u_get(uold, source_idx, 4, cell_idx)/rho;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int oct_idx = head_idx + int(block_idx);
    float inv2dx = 0.5f/dx;

    float lam_min_loc = HUGE_VALF, dev_max = 0.0f, ani_max = 0.0f, nbad = 0.0f;

    for (int cell = int(thread_idx); cell < TWOTONDIM; cell += int(threads_per_tg)) {
        int cell_idx = cell + 1;
        if (grid[oct_idx-1].refined[cell] != 0) continue;

        int ib, jb, kb;
        index_1Dto3D(cell, 2, 2, ib, jb, kb);
        int i = ib + 2, j = jb + 2, k = kb + 2;

        float c[DF_NV], w[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++) c[iv] = df_u_get(uold, oct_idx, iv + 1, cell_idx);
        df_cons_to_prim(c, w, smallr, smallc2);
        float p = max(w[DI_P], 1e-30f);

        float P[6]; df_P_sym6(w[DI_P], &w[DI_PI], P);
        float lam = df_lam_min_sym6(P)/p;
        lam_min_loc = min(lam_min_loc, lam);
        if (lam < 0.0f) nbad += 1.0f;

        float pis[6]; df_pi_to_sym6(&w[DI_PI], pis);
        float ani = sqrt(pis[0]*pis[0] + pis[1]*pis[1] + pis[2]*pis[2]
                       + 2.0f*(pis[3]*pis[3] + pis[4]*pis[4] + pis[5]*pis[5]))/p;
        ani_max = max(ani_max, ani);

        /* Pi_NS = -2 p tau S0, the Newtonian extrapolation the blowup note
         * audits.  dev = ||Pi - Pi_NS||_F / p is the Stage-1 closure
         * indicator: small means Navier-Stokes is an adequate description. */
        if (tau_pi > 0.0f) {
            float G[3][3];
            for (int a = 0; a < 3; a++) {
                G[a][0] = (uu[a][i+1][j][k] - uu[a][i-1][j][k])*inv2dx;
                G[a][1] = (uu[a][i][j+1][k] - uu[a][i][j-1][k])*inv2dx;
                G[a][2] = (uu[a][i][j][k+1] - uu[a][i][j][k-1])*inv2dx;
            }
            float divu = G[0][0] + G[1][1] + G[2][2];
            float dsum = 0.0f;
            float PIm[3][3] = {{pis[0], pis[3], pis[4]},
                               {pis[3], pis[1], pis[5]},
                               {pis[4], pis[5], pis[2]}};
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
                float S0 = 0.5f*(G[a][b] + G[b][a]) - ((a == b) ? divu/3.0f : 0.0f);
                float d = PIm[a][b] - (-2.0f*p*tau_pi*S0);
                dsum += d*d;
            }
            dev_max = max(dev_max, sqrt(dsum)/p);
        }
    }

    threadgroup float tg[4][32];
    uint lane = thread_idx % 32u, sg = thread_idx / 32u;
    float a0 = simd_min(lam_min_loc);
    float a1 = simd_max(dev_max);
    float a2 = simd_max(ani_max);
    float a3 = simd_sum(nbad);
    if (lane == 0) { tg[0][sg] = a0; tg[1][sg] = a1; tg[2][sg] = a2; tg[3][sg] = a3; }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (sg == 0 && lane == 0) {
        uint ng = (threads_per_tg + 31u)/32u;
        float b0 = HUGE_VALF, b1 = 0.0f, b2 = 0.0f, b3 = 0.0f;
        for (uint g = 0; g < ng; g++) {
            b0 = min(b0, tg[0][g]); b1 = max(b1, tg[1][g]);
            b2 = max(b2, tg[2][g]); b3 += tg[3][g];
        }
        /* diag[0] is stored as an offset positive quantity so the bitwise
         * atomic min is valid: lam_min/p can be negative, and
         * atomic_min_float_bits only orders positive IEEE-754 floats. */
        atomic_min_float_bits(&diag[0], b0 + 2.0f);
        atomic_max_float     (&diag[1], b1);
        atomic_max_float     (&diag[2], b2);
        atomic_add_float     (&diag[3], b3);
    }
}

#endif /* DFMM */
