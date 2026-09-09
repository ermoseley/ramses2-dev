/*
 * metal/kernels/dfmm.metal
 *
 * dfmm (dual-frame moment method) in 3D.  See doc/dfmm_3d.md for the frozen
 * equation, flux, source and realizability ledger; this file implements
 * Sections 3, 4 and 5 of it.
 *
 * Two stages are selectable at compile time via NDFMM (bin/Makefile DFMM=n):
 *   NDFMM=5   Stage 1, ten-moment Gaussian closure: Pi_ij evolved, Q_ijk = 0.
 *   NDFMM=15  Stage 2, adds the third central moment Q_ijk with the Wick
 *             fourth moment R_ijkl = (P_ij P_kl + P_ik P_jl + P_il P_jk)/rho.
 *
 * Scope: NDIM=3, float32, HLL (or LLF) numerical flux only.  HLLC is not
 * provided: its middle-state construction is defined for the Euler system and
 * has no standard contact reconstruction for the anisotropic pressure, and the
 * 1D reference implementation this generalises is HLL.
 *
 * State (NVAR = 5 + NDFMM):
 *   ivar 1      rho
 *   ivar 2..4   rho u_i
 *   ivar 5      E = rho|u|^2/2 + 3p/2         (gamma = 5/3 enforced)
 *   ivar 6..10  Pi_xx, Pi_yy, Pi_xy, Pi_xz, Pi_yz  with Pi_zz = -(Pi_xx+Pi_yy)
 *   ivar 11..20 Q_xxx Q_yyy Q_zzz Q_xxy Q_xxz Q_yyx Q_yyz Q_zzx Q_zzy Q_xyz
 *               (Stage 2 only)
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
 *                     [10]=tau_pi [11]=source_on [12]=tau_q
 *   dfmm_diag:        [0]=grid [1]=uold [2]=nbor [3]=diag [4]=head_idx
 *                     [5]=num_octs [6]=smallr [7]=smallc2 [8]=dx [9]=tau_pi
 *                     [10]=tau_q
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

/* Third central moment Q_ijk (Stage 2).  Ten independent components in the
 * packing order fixed by doc/dfmm_3d.md Section 2:
 *   0 xxx  1 yyy  2 zzz  3 xxy  4 xxz  5 yyx  6 yyz  7 zzx  8 zzy  9 xyz  */
#if NDFMM >= 15
#define DF_HAVE_Q 1
#define DF_NQ     10
#define DI_Q      (DI_PI + 5)
#else
#define DF_HAVE_Q 0
#define DF_NQ     0
#endif

#if DF_HAVE_Q
/* (i,j,k) of each packed slot, and the packed slot of an arbitrary (i,j,k)
 * addressed as 9i + 3j + k.  Q is fully symmetric, so DF_QMAP is invariant
 * under every permutation of its three indices. */
constant int DF_QIJK[DF_NQ][3] = {
    {0,0,0}, {1,1,1}, {2,2,2},
    {0,0,1}, {0,0,2}, {1,1,0}, {1,1,2}, {2,2,0}, {2,2,1},
    {0,1,2}};
constant int DF_QMAP[27] = {
    0,3,4,  3,5,9,  4,9,7,
    3,5,9,  5,1,6,  9,6,8,
    4,9,7,  9,6,8,  7,8,2};
#define DF_QS(i,j,k) (DF_QMAP[9*(i) + 3*(j) + (k)])

/* sym6 slot of (a,b) for the (xx,yy,zz,xy,xz,yz) packing. */
constant int DF_S6[3][3] = {{0,3,4},{3,1,5},{4,5,2}};
#endif

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

/* Offset for the realizability atomic; see dfmm_diag_kernel. */
constant float DF_LAM_OFF = 16.0f;

/* Largest admissible enlargement of the HLL wave speed, as a multiple of the
 * signal-speed estimate.  Bounded by 1/courant_factor; see df_hll_flux. */
constant float DF_ABOOST_MAX = 1.25f;

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
 *   subgrid    DF_NV * 6^3 * 4 B        (+ 64 B refined)
 *   interfaces 6 * DF_NV * 12 * 4 B
 *   total      DF_NV * 1152 B + 64 B  =  11584 B at Stage 1 (DF_NV = 10)
 *                                        23104 B at Stage 2 (DF_NV = 20)
 * against the 32768 B Apple per-threadgroup limit.
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

/* Expand sym6 to a row-major 3x3 addressed as M[3*a + b]. */
static inline void df_mat_from_sym6(thread const float *s, thread float *M) {
    M[0] = s[0]; M[1] = s[3]; M[2] = s[4];
    M[3] = s[3]; M[4] = s[1]; M[5] = s[5];
    M[6] = s[4]; M[7] = s[5]; M[8] = s[2];
}

#if DF_HAVE_Q
/* Contracted heat flux q_i = Q_ijj / 2. */
static inline void df_heat_flux(thread const float *Q, thread float *qv) {
    qv[0] = 0.5f*(Q[0] + Q[5] + Q[7]);   /* Q_xxx + Q_xyy + Q_xzz */
    qv[1] = 0.5f*(Q[3] + Q[1] + Q[8]);   /* Q_yxx + Q_yyy + Q_yzz */
    qv[2] = 0.5f*(Q[4] + Q[6] + Q[2]);   /* Q_zxx + Q_zyy + Q_zzz */
}
#endif

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
#if DF_HAVE_Q
    /* Q is a fully symmetric rank-3 tensor, so the same axis map applies to
     * all three of its indices:  Q'_abc = Q_{s(a) s(b) s(c)}, with s the
     * primed -> unprimed map used above for the velocity. */
    {
        int sg[3];
        if (idim == 1) { sg[0] = 1; sg[1] = 2; sg[2] = 0; }
        else           { sg[0] = 2; sg[1] = 0; sg[2] = 1; }
        float Qo[DF_NQ];
        for (int m = 0; m < DF_NQ; m++) Qo[m] = w[DI_Q+m];
        for (int m = 0; m < DF_NQ; m++) {
            int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
            w[DI_Q+m] = Qo[DF_QS(sg[qi], sg[qj], sg[qk])];
        }
    }
#endif
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
#if DF_HAVE_Q
    /* Inverse of the forward Q permutation: the flux slot of the primed
     * component Q'_abc is the flux slot of Q_{s(a) s(b) s(c)}.  s is a
     * permutation and DF_QMAP is symmetric, so this hits each of the ten
     * packed slots exactly once. */
    {
        int sg[3];
        if (idim == 1) { sg[0] = 1; sg[1] = 2; sg[2] = 0; }
        else           { sg[0] = 2; sg[1] = 0; sg[2] = 1; }
        float Fo[DF_NQ];
        for (int m = 0; m < DF_NQ; m++) Fo[m] = w[DI_Q+m];
        for (int m = 0; m < DF_NQ; m++) {
            int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
            w[DI_Q + DF_QS(sg[qi], sg[qj], sg[qk])] = Fo[m];
        }
    }
#endif
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
 *   F[E]     = ux (E + p) + (ux Pi_xx + uy Pi_xy + uz Pi_xz) + q_x
 *   F[Pi_ij] = ux Pi_ij + Q_ijx - (2/3) d_ij q_x
 *   F[Q_ijk] = ux Q_ijk + R_ijkx                       (Stage 2, Wick R)
 *
 * The Q and q terms vanish identically at Stage 1.  Contracting F[Pi_ij] on
 * ij gives Q_iix - 2 q_x = 0, so the flux of the traceless Pi block stays
 * traceless and the two-component storage remains consistent.
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

#if DF_HAVE_Q
    thread const float *Q = &w[DI_Q];
    float qv[3];
    df_heat_flux(Q, qv);

    F[DI_P] += qv[0];

    float tq = (2.0f/3.0f)*qv[0];
    F[DI_PI+PI_XX] += Q[DF_QS(0,0,0)] - tq;
    F[DI_PI+PI_YY] += Q[DF_QS(1,1,0)] - tq;
    F[DI_PI+PI_XY] += Q[DF_QS(0,1,0)];
    F[DI_PI+PI_XZ] += Q[DF_QS(0,2,0)];
    F[DI_PI+PI_YZ] += Q[DF_QS(1,2,0)];

    float Ps[6]; df_P_sym6(p, &w[DI_PI], Ps);
    float PM[9]; df_mat_from_sym6(Ps, PM);
    float irho = 1.0f/rho;
    for (int m = 0; m < DF_NQ; m++) {
        int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
        F[DI_Q+m] += (PM[3*qi+qj]*PM[3*qk+0]
                    + PM[3*qi+qk]*PM[3*qj+0]
                    + PM[3*qi+0 ]*PM[3*qj+qk])*irho;
    }
#endif
}

/* Is a conserved state realizable?  Requires positive density, positive
 * pressure, and a positive-semidefinite pressure tensor -- the last being the
 * statement that a^T P a = m Int (a.c)^2 f d3v >= 0 for every direction a.
 * Pi is density-like, so its conserved and primitive slots coincide. */
static inline bool df_state_ok(thread const float *U, float smallr) {
    float rho = U[DI_RHO];
    if (!(rho > smallr)) return false;
    float ekin = 0.5f*(U[DI_UX]*U[DI_UX] + U[DI_UY]*U[DI_UY]
                     + U[DI_UZ]*U[DI_UZ])/rho;
    float p = DF_GM1*(U[DI_P] - ekin);
    if (!(p > 0.0f)) return false;
    float Pm[6]; df_P_sym6(p, &U[DI_PI], Pm);
    return df_lam_min_sym6(Pm) >= 0.0f;
}

/* Signed distance of a conserved state from the realizability cone, as
 * lam_min(P)/p.  Positive inside, negative outside; a smooth function of the
 * state apart from the density and pressure floors, which is what lets it be
 * used inside a numerical flux without breaking conservation. */
static inline float df_cone_margin(thread const float *U, float smallr) {
    float rho = max(U[DI_RHO], smallr);
    float ekin = 0.5f*(U[DI_UX]*U[DI_UX] + U[DI_UY]*U[DI_UY]
                     + U[DI_UZ]*U[DI_UZ])/rho;
    float p = DF_GM1*(U[DI_P] - ekin);
    if (!(p > 0.0f)) return -1.0f;
    float Pm[6]; df_P_sym6(p, &U[DI_PI], Pm);
    return df_lam_min_sym6(Pm)/p;
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
    df_prim_to_cons(wl, UL);
    df_prim_to_cons(wr, UR);

    /* ------------------------------------------------------------------
     * Realizability-preserving wave speed (doc/dfmm_3d.md Section 5 item 2).
     *
     * Transport is the only place the scheme can leave the realizability
     * cone, and the correct statement of transport positivity for a moment
     * system is to bound the FLUX, not the signal speed.  The split states
     *     V_L = U_L - F_L/a ,   V_R = U_R + F_R/a
     * are the states whose convex combinations the update forms, so it is
     * they -- not the cell states -- that must have a positive-semidefinite
     * pressure tensor.  Enlarge a until both do.  V -> U as a -> infinity, so
     * this terminates whenever the reconstructed states are themselves
     * admissible; the bound of eight doublings is a backstop for when they
     * are not.
     *
     * Enlarging a raises the numerical viscosity of that face and nothing
     * else -- the state is never modified, clipped or projected.
     *
     * THE ENLARGEMENT IS CAPPED, and the cap is not a tuning knob.  dt is set
     * by dfmm_cmpdt_kernel from the UN-enlarged signal speed, so a face flux
     * built with speed a > a0/courant_factor violates the CFL condition that
     * dt was chosen to satisfy.  An uncapped doubling loop is therefore
     * unconditionally unstable: measured on the smooth Delta=1.644, tau=1e-3
     * case, eight doublings (256x) destroyed a previously exactly-conservative
     * run in a single step.  DF_ABOOST_MAX must stay at or below
     * 1/courant_factor; read_params.f90 enforces courant_factor <= 0.8.
     *
     * When the cap is not enough, the flux is left at the plain HLL speeds and
     * the state is allowed to leave the cone, to be REPORTED by
     * dfmm_diag_kernel.  Silently enlarging further would trade a visible
     * physical diagnostic for an invisible numerical instability.
     * ------------------------------------------------------------------ */
    float a0 = max(max(-SL, SR), 1.0e-20f);
    {
        float ia = 1.0f/a0;
        float VL[DF_NV], VR[DF_NV];
        for (int k = 0; k < DF_NV; k++) {
            VL[k] = UL[k] - FL[k]*ia;
            VR[k] = UR[k] + FR[k]*ia;
        }
        /* Signed margin of the worst split state, in units of its own
         * pressure.  Negative means that state is outside the cone. */
        float m = min(df_cone_margin(VL, smallr), df_cone_margin(VR, smallr));
        /* CONTINUOUS widening.  It must be continuous in the state, not a
         * branch: an oct-boundary face is reconstructed and solved
         * independently by the two threadgroups that own the adjoining octs,
         * and conservation depends on both arriving at the same flux.  A
         * discrete "enlarge if inadmissible" test turns a round-off
         * difference in the margin into a finite difference in the wave
         * speed, and hence into a non-telescoping flux: measured, that
         * destroyed mass and energy conservation (mcons went from 0 to
         * -2.3e-2 in two steps) on a case that is otherwise exact. */
        float w = 1.0f + (DF_ABOOST_MAX - 1.0f)*clamp(-m*10.0f, 0.0f, 1.0f);
        SL -= (w - 1.0f)*a0;
        SR += (w - 1.0f)*a0;
    }

    if (SL >= 0.0f) { for (int k = 0; k < DF_NV; k++) F[k] = FL[k]; return; }
    if (SR <= 0.0f) { for (int k = 0; k < DF_NV; k++) F[k] = FR[k]; return; }

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
 *              -(2/3) div q                 (from the energy flux term q_k)
 *   Pi_ij:     -u.grad Pi_ij - Pi_ij div u  (transport, density-like)
 *              -d_k Q_ijk + (2/3) d_ij div q          (from the Pi_ij flux)
 *   Q_ijk:     -u.grad Q_ijk - Q_ijk div u  (transport, density-like)
 *              -d_l R_ijkl + T_Q1                     (see below)
 *
 * All of these are *flux*-derived or non-stiff, and are handled
 * conservatively by the Riemann solve; they appear here only because the
 * predictor works in primitive variables.
 *
 * What is deliberately absent, in both cases because it is stiff and must be
 * integrated together with its own relaxation sink (dfmm_source_kernel does
 * that exactly):
 *   Pi_ij:  -2 p S0_ij - [Pi_ik G_jk + Pi_jk G_ik]^dev  and  -Pi_ij/tau_Pi
 *   Q_ijk:  -[Q_jkl G_il + Q_ikl G_jl + Q_ijl G_kl]     and  -Q_ijk/tau_q
 * An unrelaxed half-step copy of a stiff production term overshoots the
 * interface states whenever dt >> tau.  The 1D reference makes the same
 * choice: sources follow the flux update.
 *
 * T_Q1 = (1/rho)(P_jk d_l P_li + P_ik d_l P_lj + P_ij d_l P_lk) is *not*
 * stiff, and it is kept here because it very nearly cancels d_l R_ijkl: the
 * divergence-of-P parts cancel analytically, leaving only a temperature
 * gradient (doc/dfmm_3d.md Section 4).  Dropping one of the pair would leave
 * an O(1) unbalanced term in the reconstruction.
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

#if DF_HAVE_Q
            thread const float * const sl[3] = { sx, sy, sz };

            /* div q from the half-slopes:  d_k q_k = (1/2) d_k Q_kll */
            float divq = 0.0f;
            for (int kk = 0; kk < 3; kk++)
                for (int ll = 0; ll < 3; ll++)
                    divq += 0.5f*sl[kk][DI_Q + DF_QS(kk,ll,ll)];

            /* d_k Q_ijk for the five stored (ij) */
            float dQ[5];
            dQ[PI_XX] = sx[DI_Q+DF_QS(0,0,0)] + sy[DI_Q+DF_QS(0,0,1)] + sz[DI_Q+DF_QS(0,0,2)];
            dQ[PI_YY] = sx[DI_Q+DF_QS(1,1,0)] + sy[DI_Q+DF_QS(1,1,1)] + sz[DI_Q+DF_QS(1,1,2)];
            dQ[PI_XY] = sx[DI_Q+DF_QS(0,1,0)] + sy[DI_Q+DF_QS(0,1,1)] + sz[DI_Q+DF_QS(0,1,2)];
            dQ[PI_XZ] = sx[DI_Q+DF_QS(0,2,0)] + sy[DI_Q+DF_QS(0,2,1)] + sz[DI_Q+DF_QS(0,2,2)];
            dQ[PI_YZ] = sx[DI_Q+DF_QS(1,2,0)] + sy[DI_Q+DF_QS(1,2,1)] + sz[DI_Q+DF_QS(1,2,2)];

            src[DI_PI+PI_XX] -= dQ[PI_XX] - (2.0f/3.0f)*divq;
            src[DI_PI+PI_YY] -= dQ[PI_YY] - (2.0f/3.0f)*divq;
            src[DI_PI+PI_XY] -= dQ[PI_XY];
            src[DI_PI+PI_XZ] -= dQ[PI_XZ];
            src[DI_PI+PI_YZ] -= dQ[PI_YZ];
            src[DI_P]        -= (2.0f/3.0f)*divq;

            /* -d_l R_ijkl + T_Q1, with the div P terms already cancelled:
             *   -(1/rho)  sum_l ( d_l P_ij P_kl + d_l P_ik P_jl + d_l P_jk P_il )
             *   +(1/rho^2)sum_l ( P_ij P_kl + P_ik P_jl + P_il P_jk ) d_l rho
             * For P = p I this collapses to
             *   -( d_ij p d_k theta + d_ik p d_j theta + d_jk p d_i theta ),
             * i.e. exactly the term that drives the Fourier heat flux. */
            float PM[9];
            {
                float Ps[6];
                df_P_sym6(p, &w[DI_PI], Ps);
                df_mat_from_sym6(Ps, PM);
            }
            float sP[3][6];
            for (int l = 0; l < 3; l++) {
                df_pi_to_sym6(&sl[l][DI_PI], sP[l]);
                sP[l][0] += sl[l][DI_P];
                sP[l][1] += sl[l][DI_P];
                sP[l][2] += sl[l][DI_P];
            }
            float dr[3]  = { sx[DI_RHO], sy[DI_RHO], sz[DI_RHO] };
            float irho2  = irho*irho;
            for (int m = 0; m < DF_NQ; m++) {
                int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
                float acc = 0.0f;
                for (int l = 0; l < 3; l++) {
                    float t1 = sP[l][DF_S6[qi][qj]]*PM[3*qk+l]
                             + sP[l][DF_S6[qi][qk]]*PM[3*qj+l]
                             + sP[l][DF_S6[qj][qk]]*PM[3*qi+l];
                    float t2 = PM[3*qi+qj]*PM[3*qk+l]
                             + PM[3*qi+qk]*PM[3*qj+l]
                             + PM[3*qi+l ]*PM[3*qj+qk];
                    acc += -t1*irho + t2*dr[l]*irho2;
                }
                src[DI_Q+m] += acc;
            }
#endif
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
 * dfmm_source_kernel — production terms and exact BGK relaxation
 *
 * Operator order (doc/dfmm_3d.md Section 4): this runs after transport and
 * before uold is overwritten.  Gradients are taken from uold, not unew:
 * set_unew only copies uold -> unew for octs at ilevel, so unew's neighbour
 * data is not valid for a centred difference.  This mirrors the 1D reference,
 * which evaluates du/dx on the pre-update state.
 *
 * Both moment blocks are advanced by the same asymptotic-preserving map.  For
 * dX/dt = rate - X/tau with `rate` frozen over the step the exact solution is
 *
 *     X <- X_old d + tau (1 - d) rate,        d = exp(-dt/tau),
 *
 * and `rate` must be the *whole* non-stiff right-hand side: the production
 * terms T plus the transport rate the Godunov step already applied, which is
 * recoverable as (unew - uold)/dt.
 *
 * Including the transport rate is not optional for Q.  Its flux carries the
 * Wick fourth moment R_ijkl, which is O(1) and cancels all but a temperature
 * gradient against T_Q1; dropping it whenever dt >> tau_q would replace the
 * Fourier heat flux -(5/2) tau_q p grad theta by an unrelated quantity.  For
 * Pi the same omission is only O(tau^2), a Burnett-order correction, but the
 * two blocks are treated identically here for uniformity.
 *
 * tau <= 0 means collisionless: d = 1 and tau(1-d) -> dt, so the map degrades
 * to X_new = X_transported + dt T, which is the correct explicit update.
 * ========================================================================= */

/* Stage 1 needs only the velocity stencil.  Stage 2 additionally needs p and
 * the five stored Pi components, to form div P for the Q production term.
 *   Stage 1:  3 * 216 * 4 =  2592 B of threadgroup memory
 *   Stage 2:  9 * 216 * 4 =  7776 B
 */
#if DF_HAVE_Q
#define DF_NSRC 10        /* u(3), p, Pi(5), theta */
#else
#define DF_NSRC 5         /* u(3), p, theta */
#endif
#define DF_ST_P  3
#define DF_ST_TH 4
#if DF_HAVE_Q
#define DF_ST_PI 5
#endif

/* Closure selector (mirrors read_params.f90).
 *   0  evolve  -- Pi and Q are evolved with the AP relaxation map below.
 *   1  ns      -- Pi and Q are OVERWRITTEN each step with their
 *                 Chapman-Enskog values, which turns the identical flux
 *                 machinery into compressible Navier-Stokes-Fourier with
 *                 mu = p tau_Pi and Pr = tau_Pi/tau_q.  See the comment on
 *                 the branch below for why this is the right comparison run.
 */
constant int DF_CLOSURE_NS = 1;


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
    constant float     &tau_q    [[buffer(12)]],
    constant int       &closure  [[buffer(13)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_octs) return;

    threadgroup float st[DF_NSRC][6][6][6];

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

        float cc[5];
        for (int iv = 0; iv < 5; iv++)
            cc[iv] = df_u_get(uold, source_idx, iv + 1, cell_idx);
        float rho = max(cc[DI_RHO], smallr);
        st[0][i][j][k] = cc[DI_UX]/rho;
        st[1][i][j][k] = cc[DI_UY]/rho;
        st[2][i][j][k] = cc[DI_UZ]/rho;
        float ekin = 0.5f*(cc[DI_UX]*cc[DI_UX] + cc[DI_UY]*cc[DI_UY]
                         + cc[DI_UZ]*cc[DI_UZ])/rho;
        float pl = max(DF_GM1*(cc[DI_P] - ekin), rho*smallc2/DF_GAMMA);
        st[DF_ST_P ][i][j][k] = pl;
        st[DF_ST_TH][i][j][k] = pl/rho;
#if DF_HAVE_Q
        for (int m = 0; m < 5; m++)
            st[DF_ST_PI+m][i][j][k] = df_u_get(uold, source_idx, 6 + m, cell_idx);
#endif
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int oct_idx = head_idx + int(block_idx);
    float inv2dx = 0.5f/dx;
    float invdt  = 1.0f/dt;

    /* Relaxation weights: `decay` = exp(-dt/tau), `tau_omd` = tau(1 - decay).
     * The series form below the cutoff keeps tau_omd -> dt to float32
     * accuracy, where exp(-x) - 1 loses all of its significant digits. */
    float decay_pi = 1.0f, tomd_pi = dt;
    if (tau_pi > 0.0f) {
        float xr = dt/tau_pi, omd;
        if (xr < 1.0e-4f) { omd = xr*(1.0f - 0.5f*xr); decay_pi = 1.0f - omd; }
        else              { decay_pi = exp(-xr);       omd = 1.0f - decay_pi; }
        tomd_pi = tau_pi*omd;
    }
#if DF_HAVE_Q
    float decay_q = 1.0f, tomd_q = dt;
    if (tau_q > 0.0f) {
        float xr = dt/tau_q, omd;
        if (xr < 1.0e-4f) { omd = xr*(1.0f - 0.5f*xr); decay_q = 1.0f - omd; }
        else              { decay_q = exp(-xr);        omd = 1.0f - decay_q; }
        tomd_q = tau_q*omd;
    }
#endif

    for (int cell = int(thread_idx); cell < TWOTONDIM; cell += int(threads_per_tg)) {
        int cell_idx = cell + 1;
        if (grid[oct_idx-1].refined[cell] != 0) continue;

        int ib, jb, kb;
        index_1Dto3D(cell, 2, 2, ib, jb, kb);
        int i = ib + 2, j = jb + 2, k = kb + 2;   /* centre 2x2x2 of the stencil */

        /* Post-transport state of this cell, and the pre-transport moments so
         * that the transport rate (xt - xo)/dt can be recovered. */
        float c[DF_NV], w[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++) c[iv] = df_u_get(unew, oct_idx, iv + 1, cell_idx);
        df_cons_to_prim(c, w, smallr, smallc2);
        float p = w[DI_P];

        /* ------------------------------------------------------------------
         * Navier-Stokes-Fourier closure.
         *
         * Instead of evolving Pi and Q, overwrite them with their
         * Chapman-Enskog values every step:
         *     Pi_ij  = -2 p tau_Pi S0_ij
         *     Q_ijk  = -tau_q p ( d_jk d_i theta + d_ik d_j theta
         *                       + d_ij d_k theta )   =>  q = -(5/2) tau_q p grad theta
         * The flux ledger already carries Pi_ik in the momentum flux and
         * u_i Pi_ik + q_k in the energy flux, so this turns the SAME kernel
         * into a compressible Navier-Stokes-Fourier solver with explicit
         * viscosity mu = p tau_Pi and Prandtl number tau_Pi/tau_q -- not an
         * inviscid Euler run with numerical viscosity standing in for it.
         *
         * This is deliberately the comparison run: identical initial
         * condition, identical grid, identical Riemann solver, identical
         * transport coefficients.  The only thing that differs from the
         * moment run is the closure, which is exactly the variable under
         * study.  A separately written viscous-flux implementation would
         * confound the closure with the discretisation.
         *
         * Note the price: this branch reintroduces a parabolic stability
         * constraint dt < dx^2 / (2 ndim D), which dfmm_cmpdt_kernel applies
         * only when this closure is selected.  The evolved moment system is
         * hyperbolic and has no such constraint -- one of the practical
         * reasons for preferring it.
         * ---------------------------------------------------------------- */
        if (closure == DF_CLOSURE_NS) {
            float G[3][3];
            for (int a = 0; a < 3; a++) {
                G[a][0] = (st[a][i+1][j][k] - st[a][i-1][j][k])*inv2dx;
                G[a][1] = (st[a][i][j+1][k] - st[a][i][j-1][k])*inv2dx;
                G[a][2] = (st[a][i][j][k+1] - st[a][i][j][k-1])*inv2dx;
            }
            float divu = G[0][0] + G[1][1] + G[2][2];
            float c2 = -2.0f*p*max(tau_pi, 0.0f);
            float S0[3][3];
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++)
                S0[a][b] = 0.5f*(G[a][b] + G[b][a]) - ((a == b) ? divu/3.0f : 0.0f);

            df_u_set(unew, oct_idx, 6 + PI_XX, cell_idx, c2*S0[0][0]);
            df_u_set(unew, oct_idx, 6 + PI_YY, cell_idx, c2*S0[1][1]);
            df_u_set(unew, oct_idx, 6 + PI_XY, cell_idx, c2*S0[0][1]);
            df_u_set(unew, oct_idx, 6 + PI_XZ, cell_idx, c2*S0[0][2]);
            df_u_set(unew, oct_idx, 6 + PI_YZ, cell_idx, c2*S0[1][2]);
#if DF_HAVE_Q
            float gth[3] = {
                (st[DF_ST_TH][i+1][j][k] - st[DF_ST_TH][i-1][j][k])*inv2dx,
                (st[DF_ST_TH][i][j+1][k] - st[DF_ST_TH][i][j-1][k])*inv2dx,
                (st[DF_ST_TH][i][j][k+1] - st[DF_ST_TH][i][j][k-1])*inv2dx};
            float cq = -max(tau_q, 0.0f)*p;
            for (int m = 0; m < DF_NQ; m++) {
                int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
                float v = ((qj == qk) ? gth[qi] : 0.0f)
                        + ((qi == qk) ? gth[qj] : 0.0f)
                        + ((qi == qj) ? gth[qk] : 0.0f);
                df_u_set(unew, oct_idx, 11 + m, cell_idx, cq*v);
            }
#endif
            continue;
        }

        float xo[NDFMM], xt[NDFMM], T[NDFMM];
        for (int m = 0; m < NDFMM; m++) {
            xo[m] = df_u_get(uold, oct_idx, 6 + m, cell_idx);
            xt[m] = w[DI_PI+m];
            T[m]  = 0.0f;
        }

        if (source_on != 0) {
            /* G_ij = d_j u_i, centred, from uold */
            float G[3][3];
            for (int a = 0; a < 3; a++) {
                G[a][0] = (st[a][i+1][j][k] - st[a][i-1][j][k])*inv2dx;
                G[a][1] = (st[a][i][j+1][k] - st[a][i][j-1][k])*inv2dx;
                G[a][2] = (st[a][i][j][k+1] - st[a][i][j][k-1])*inv2dx;
            }

            float pis[6]; df_pi_to_sym6(xt, pis);
            float PIm[9]; df_mat_from_sym6(pis, PIm);

            /* S[Pi_ij] = -2 p S0_ij - [Pi_ik G_jk + Pi_jk G_ik]^dev */
            float TP[3][3];
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
                float S  = 0.5f*(G[a][b] + G[b][a]);
                float PG = 0.0f;
                for (int cc = 0; cc < 3; cc++)
                    PG += PIm[3*a+cc]*G[b][cc] + PIm[3*b+cc]*G[a][cc];
                TP[a][b] = -2.0f*p*S - PG;
            }
            float trT = (TP[0][0] + TP[1][1] + TP[2][2])/3.0f;
            TP[0][0] -= trT; TP[1][1] -= trT; TP[2][2] -= trT;

            T[PI_XX] = TP[0][0];
            T[PI_YY] = TP[1][1];
            T[PI_XY] = TP[0][1];
            T[PI_XZ] = TP[0][2];
            T[PI_YZ] = TP[1][2];

#if DF_HAVE_Q
            /* div P from uold: d_l P_li = d_i p + d_l Pi_li */
            float gp[3] = {
                (st[DF_ST_P][i+1][j][k] - st[DF_ST_P][i-1][j][k])*inv2dx,
                (st[DF_ST_P][i][j+1][k] - st[DF_ST_P][i][j-1][k])*inv2dx,
                (st[DF_ST_P][i][j][k+1] - st[DF_ST_P][i][j][k-1])*inv2dx};
            float gpi[3][5];
            for (int m = 0; m < 5; m++) {
                gpi[0][m] = (st[DF_ST_PI+m][i+1][j][k] - st[DF_ST_PI+m][i-1][j][k])*inv2dx;
                gpi[1][m] = (st[DF_ST_PI+m][i][j+1][k] - st[DF_ST_PI+m][i][j-1][k])*inv2dx;
                gpi[2][m] = (st[DF_ST_PI+m][i][j][k+1] - st[DF_ST_PI+m][i][j][k-1])*inv2dx;
            }
            float divP[3];
            divP[0] = gp[0] + gpi[0][PI_XX] + gpi[1][PI_XY] + gpi[2][PI_XZ];
            divP[1] = gp[1] + gpi[0][PI_XY] + gpi[1][PI_YY] + gpi[2][PI_YZ];
            divP[2] = gp[2] + gpi[0][PI_XZ] + gpi[1][PI_YZ]
                    - (gpi[2][PI_XX] + gpi[2][PI_YY]);

            /* P and Q of the post-transport state */
            float Ps[6]; df_P_sym6(p, xt, Ps);
            float PM[9]; df_mat_from_sym6(Ps, PM);
            thread const float *Qc = &xt[5];

            /* S[Q_ijk] = (1/rho)( P_jk d_l P_li + P_ik d_l P_lj + P_ij d_l P_lk )
             *            - ( Q_jkl G_il + Q_ikl G_jl + Q_ijl G_kl )
             * The -d_l R_ijkl companion of the first term is carried by the
             * flux and enters here through the transport rate. */
            for (int m = 0; m < DF_NQ; m++) {
                int qi = DF_QIJK[m][0], qj = DF_QIJK[m][1], qk = DF_QIJK[m][2];
                float prod = (PM[3*qj+qk]*divP[qi]
                            + PM[3*qi+qk]*divP[qj]
                            + PM[3*qi+qj]*divP[qk])/w[DI_RHO];
                float dist = 0.0f;
                for (int l = 0; l < 3; l++)
                    dist += Qc[DF_QS(qj,qk,l)]*G[qi][l]
                          + Qc[DF_QS(qi,qk,l)]*G[qj][l]
                          + Qc[DF_QS(qi,qj,l)]*G[qk][l];
                T[5+m] = prod - dist;
            }
#endif
        }

        for (int m = 0; m < 5; m++) {
            float rate = (xt[m] - xo[m])*invdt + T[m];
            df_u_set(unew, oct_idx, 6 + m, cell_idx, xo[m]*decay_pi + tomd_pi*rate);
        }
#if DF_HAVE_Q
        for (int m = 5; m < NDFMM; m++) {
            float rate = (xt[m] - xo[m])*invdt + T[m];
            df_u_set(unew, oct_idx, 6 + m, cell_idx, xo[m]*decay_q + tomd_q*rate);
        }
#endif
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
    constant float       &tau_pi           [[buffer(11)]],
    constant float       &tau_q            [[buffer(12)]],
    constant int         &closure          [[buffer(13)]],
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

            /* Parabolic constraint, needed ONLY under the Navier-Stokes
             * closure, where Pi and q are algebraic in the gradients and the
             * momentum/energy fluxes are therefore genuinely diffusive:
             *     dt <= courant * dx^2 / (2 ndim D),   D = max(nu, chi)
             * with nu = mu/rho = tau_Pi p/rho and, since
             * kappa = (5/2) tau_q p (k_B/m) and c_p = (5/2)(k_B/m),
             * chi = kappa/(rho c_p) = tau_q p/rho.
             *
             * The evolved moment system needs no such limit: its transport is
             * hyperbolic with the finite signal speed used above, and the
             * relaxation is integrated by an exact exponential map.  That is a
             * real advantage of the moment formulation and it is worth being
             * able to see it in the timestep. */
            if (closure == DF_CLOSURE_NS) {
                float dif = max(max(tau_pi, tau_q), 0.0f)*w[DI_P]/w[DI_RHO];
                if (dif > 0.0f)
                    dt_loc = min(dt_loc, courant_factor*dx*dx/(2.0f*3.0f*dif));
            }
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
 *   diag[0]  min over cells of max(lam_min(P)/p, 0)   (see below)
 *   diag[6]  max over cells of max(-lam_min(P)/p, 0)  (see below)
 *   diag[1]  max over cells of ||Pi - Pi_NS||_F / p   (Navier-Stokes deviation)
 *   diag[2]  max over cells of ||Pi||_F / p           (anisotropy amplitude)
 *   diag[3]  count of cells with lam_min(P) < 0       (as a float)
 *   diag[4]  max over cells of ||q - q_CE|| / (p c_s) (Fourier deviation)
 *   diag[5]  max over cells of ||q|| / (p c_s)        (heat-flux amplitude)
 *
 * Slots 4 and 5 are written only at Stage 2, where q is a primary variable;
 * at Stage 1 they stay at their initialised zero and the host does not print
 * them.  Pi_NS = -2 p tau_Pi S0 and q_CE = -(5/2) tau_q p grad theta are the
 * Newtonian and Fourier extrapolations the blowup note audits, so the two
 * deviations are the quantities the study reads.
 *
 * A minimum is used for the realizability margin deliberately: a maximum of
 * the slack is structurally blind to cells sitting at the cone boundary.
 *
 * lam_min(P)/p is reported as a POSITIVE PAIR, not as one offset float.  The
 * only lock-free float atomics available here are min/max on the raw bit
 * pattern, which order IEEE-754 correctly for non-negative values only.  The
 * obvious dodge -- atomically minimising lam/p + 2 -- silently clips anything
 * below -2, because a negative sum sets the sign bit and its bit pattern then
 * compares as a huge unsigned.  That clip is invisible in the gate problems,
 * where |Pi|/p << 1, and catastrophic on the blowup problem, where |Pi|/p
 * reaches order ten: it would report a floor of -2 for the one quantity the
 * whole study is about.  Splitting into
 *     diag[0] = min over cells of max( lam/p, 0)
 *     diag[6] = max over cells of max(-lam/p, 0)
 * keeps both atomics on non-negative values and is exact at any magnitude.
 * The host reconstructs  min lam/p = (diag[6] > 0) ? -diag[6] : diag[0].
 * ========================================================================= */

/* u_x, u_y, u_z, and at Stage 2 also theta = p/rho for grad theta. */
#if DF_HAVE_Q
#define DF_NDIAG 4
#else
#define DF_NDIAG 3
#endif

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
    constant float     &tau_q    [[buffer(10)]],
    device float       *dbg      [[buffer(11)]],
    uint block_idx      [[threadgroup_position_in_grid]],
    uint thread_idx     [[thread_position_in_threadgroup]],
    uint threads_per_tg [[threads_per_threadgroup]])
{
    if (int(block_idx) >= num_octs) return;

    threadgroup float st[DF_NDIAG][6][6][6];
    for (int wi = int(thread_idx); wi < 216; wi += int(threads_per_tg)) {
        int i_sg, j_sg, k_sg;
        index_1Dto3D(wi/8, 3, 3, i_sg, j_sg, k_sg);
        int subgrid_idx = head_idx + int(block_idx);
        int source_idx  = nbor_get(nbor, subgrid_idx, wi/8 + 1);
        int cell_idx    = wi%8 + 1;
        int ib, jb, kb;
        index_1Dto3D(cell_idx - 1, 2, 2, ib, jb, kb);
        int i = ib + 2*i_sg, j = jb + 2*j_sg, k = kb + 2*k_sg;

        float cc[5];
        for (int iv = 0; iv < 5; iv++)
            cc[iv] = df_u_get(uold, source_idx, iv + 1, cell_idx);
        float rho = max(cc[DI_RHO], smallr);
        st[0][i][j][k] = cc[DI_UX]/rho;
        st[1][i][j][k] = cc[DI_UY]/rho;
        st[2][i][j][k] = cc[DI_UZ]/rho;
#if DF_HAVE_Q
        float ekin = 0.5f*(cc[DI_UX]*cc[DI_UX] + cc[DI_UY]*cc[DI_UY]
                         + cc[DI_UZ]*cc[DI_UZ])/rho;
        st[3][i][j][k] = max(DF_GM1*(cc[DI_P] - ekin), rho*smallc2/DF_GAMMA)/rho;
#endif
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    int oct_idx = head_idx + int(block_idx);
    float inv2dx = 0.5f/dx;

    float lam_min_loc = HUGE_VALF, dev_max = 0.0f, ani_max = 0.0f, nbad = 0.0f;
    float devq_max = 0.0f, qmax = 0.0f, viol_max = 0.0f;

    for (int cell = int(thread_idx); cell < TWOTONDIM; cell += int(threads_per_tg)) {
        int cell_idx = cell + 1;
        if (grid[oct_idx-1].refined[cell] != 0) continue;

        int ib, jb, kb;
        index_1Dto3D(cell, 2, 2, ib, jb, kb);
        int i = ib + 2, j = jb + 2, k = kb + 2;

        float c[DF_NV], w[DF_NV];
        for (int iv = 0; iv < DF_NV; iv++) c[iv] = df_u_get(uold, oct_idx, iv + 1, cell_idx);
        df_cons_to_prim(c, w, smallr, smallc2);
        /* ONE pressure, used for both the tensor and the normalisation.
         * df_cons_to_prim has already floored it at rho*smallc2/gamma > 0, so
         * no extra floor is needed here -- and an extra floor is exactly the
         * bug this replaces: dividing by max(p, 1e-30) while building P from
         * the unfloored p made an empty cell (rho = smallr, p = 6e-31) report
         * lam = (rho smallc2/gamma)/1e-30 = 1/gamma = 0.6 with ||Pi||/p = 0,
         * a state-independent constant that masked the real minimum. */
        float p = w[DI_P];

        float P[6]; df_P_sym6(p, &w[DI_PI], P);
        float lam = df_lam_min_sym6(P)/p;
        /* A NaN must count as a violation rather than be swallowed by min/max,
         * whose NaN handling would otherwise report a healthy state for a run
         * that has already failed. */
        if (!(lam == lam)) { nbad += 1.0f; lam = -HUGE_VALF; }
        lam_min_loc = min(lam_min_loc, lam);
        if (lam < 0.0f) nbad += 1.0f;

        float pis[6]; df_pi_to_sym6(&w[DI_PI], pis);
        float ani = sqrt(pis[0]*pis[0] + pis[1]*pis[1] + pis[2]*pis[2]
                       + 2.0f*(pis[3]*pis[3] + pis[4]*pis[4] + pis[5]*pis[5]))/p;
        ani_max = max(ani_max, ani);

        /* Self-consistency of the two realizability measures.  P = p I + Pi
         * with Pi traceless, so Weyl's inequality forces
         *     lam_min(P)/p >= 1 - ||Pi||_F/p
         * cell by cell.  viol > 0 is therefore impossible and its appearance
         * means one of the two is being computed wrongly, not that the state
         * is unusual.  Kept permanently: it costs two flops and it is the
         * check that would have caught the offset-atomic clip immediately. */
        viol_max = max(viol_max, (1.0f - ani) - lam);
        if (viol_max > 1.0e-3f) {
            dbg[0] = w[DI_RHO];  dbg[1] = w[DI_P];   dbg[2] = lam;
            dbg[3] = ani;        dbg[4] = pis[0];    dbg[5] = pis[1];
            dbg[6] = pis[2];     dbg[7] = pis[3];    dbg[8] = pis[4];
            dbg[9] = pis[5];     dbg[10] = float(oct_idx); dbg[11] = float(cell);
        }

        bool need_grad = (tau_pi > 0.0f);
#if DF_HAVE_Q
        need_grad = need_grad || (tau_q > 0.0f);
#endif
        float G[3][3];
        if (need_grad) {
            for (int a = 0; a < 3; a++) {
                G[a][0] = (st[a][i+1][j][k] - st[a][i-1][j][k])*inv2dx;
                G[a][1] = (st[a][i][j+1][k] - st[a][i][j-1][k])*inv2dx;
                G[a][2] = (st[a][i][j][k+1] - st[a][i][j][k-1])*inv2dx;
            }
        }

        /* Pi_NS = -2 p tau S0, the Newtonian extrapolation the blowup note
         * audits.  dev = ||Pi - Pi_NS||_F / p is the Stage-1 closure
         * indicator: small means Navier-Stokes is an adequate description. */
        if (tau_pi > 0.0f) {
            float divu = G[0][0] + G[1][1] + G[2][2];
            float dsum = 0.0f;
            float PIm[9]; df_mat_from_sym6(pis, PIm);
            for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) {
                float S0 = 0.5f*(G[a][b] + G[b][a]) - ((a == b) ? divu/3.0f : 0.0f);
                float d = PIm[3*a+b] - (-2.0f*p*tau_pi*S0);
                dsum += d*d;
            }
            dev_max = max(dev_max, sqrt(dsum)/p);
        }

#if DF_HAVE_Q
        /* q_CE = -(5/2) tau_q p grad theta is the Fourier heat flux of the
         * BGK Chapman-Enskog expansion (doc/dfmm_3d.md Section 4).  Both
         * measures are normalised by p c_s, the natural heat-flux scale. */
        {
            float qv[3]; df_heat_flux(&w[DI_Q], qv);
            float nrm = 1.0f/(p*sqrt(DF_GAMMA*p/w[DI_RHO]));
            qmax = max(qmax, sqrt(qv[0]*qv[0] + qv[1]*qv[1] + qv[2]*qv[2])*nrm);
            if (tau_q > 0.0f) {
                float gth[3] = {
                    (st[3][i+1][j][k] - st[3][i-1][j][k])*inv2dx,
                    (st[3][i][j+1][k] - st[3][i][j-1][k])*inv2dx,
                    (st[3][i][j][k+1] - st[3][i][j][k-1])*inv2dx};
                float dsum = 0.0f;
                for (int d = 0; d < 3; d++) {
                    float e = qv[d] + 2.5f*tau_q*p*gth[d];
                    dsum += e*e;
                }
                devq_max = max(devq_max, sqrt(dsum)*nrm);
            }
        }
#endif
    }

    threadgroup float tg[9][32];
    uint lane = thread_idx % 32u, sg = thread_idx / 32u;
    float a0 = simd_min(max(lam_min_loc, 0.0f));
    float a6 = simd_max(max(-lam_min_loc, 0.0f));
    float a8 = simd_max(viol_max);
    /* Primary realizability report.  DF_LAM_OFF - lam is positive for every
     * lam < DF_LAM_OFF, so a single positive-float atomic max carries the
     * minimum exactly, with no sign-bit hazard and no clip on the negative
     * side.  Only lam > DF_LAM_OFF saturates, which is the uninteresting
     * direction (a hugely positive minimum eigenvalue). */
    float a7 = simd_min(lam_min_loc);
    float a1 = simd_max(dev_max);
    float a2 = simd_max(ani_max);
    float a3 = simd_sum(nbad);
    float a4 = simd_max(devq_max);
    float a5 = simd_max(qmax);
    if (lane == 0) {
        tg[0][sg] = a0; tg[1][sg] = a1; tg[2][sg] = a2;
        tg[3][sg] = a3; tg[4][sg] = a4; tg[5][sg] = a5;
        tg[6][sg] = a6;
        tg[7][sg] = a7;
        tg[8][sg] = a8;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    if (sg == 0 && lane == 0) {
        uint ng = (threads_per_tg + 31u)/32u;
        float b0 = HUGE_VALF, b1 = 0.0f, b2 = 0.0f, b3 = 0.0f;
        float b4 = 0.0f, b5 = 0.0f, b6 = 0.0f, b0lam = HUGE_VALF, b8 = 0.0f;
        for (uint g = 0; g < ng; g++) {
            b0lam = min(b0lam, tg[7][g]);
            b8 = max(b8, tg[8][g]);
            b0 = min(b0, tg[0][g]); b1 = max(b1, tg[1][g]);
            b2 = max(b2, tg[2][g]); b3 += tg[3][g];
            b4 = max(b4, tg[4][g]); b5 = max(b5, tg[5][g]);
            b6 = max(b6, tg[6][g]);
        }
        atomic_min_float_bits(&diag[0], b0);
        atomic_max_float     (&diag[7], max(DF_LAM_OFF - b0lam, 0.0f));
        atomic_max_float     (&diag[1], b1);
        atomic_max_float     (&diag[2], b2);
        atomic_add_float     (&diag[3], b3);
        atomic_max_float     (&diag[4], b4);
        atomic_max_float     (&diag[5], b5);
        atomic_max_float     (&diag[6], b6);
        atomic_max_float     (&diag[8], b8);
    }
}

#endif /* DFMM */
