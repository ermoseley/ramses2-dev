! Launch configuration of hydro_integrator_kernel in pure-hydro builds, measured with HLLC on a
! uniform periodic 3D Sedov blast at 512^3 on an H200 NVL and checked on an A100.
! launch_bounds(T,B): T threads per block, B blocks per SM, so at most 65536/(T*B) registers.
#ifndef HYDRO_LAUNCH_H
#define HYDRO_LAUNCH_H

#ifndef NSUBGRID
#define NSUBGRID 1
#endif

#ifndef MHD
#if defined(NPRE) && NPRE == 4
! Single precision: no bound. 40 registers unbounded; any bound gives 56, -2 to -7%.
#define HYDRO_LB
#if NSUBGRID == 1
#define HYDRO_THREADS 64
#elif NSUBGRID == 2
! 128 threads: +10% over 256.
#define HYDRO_THREADS 128
#else
! 512 = cells per block (8^3).
#define HYDRO_THREADS 512
#endif
#else
! Double precision: register-limited (92 registers unbounded at nsubgrid=1).
#if NSUBGRID == 1
! 64 registers, 16 blocks per SM: +15%.
#define HYDRO_THREADS 64
#define HYDRO_LB launch_bounds(64,16)
#elif NSUBGRID == 2
! 64 registers, 4 blocks per SM: +32%.
#define HYDRO_THREADS 256
#define HYDRO_LB launch_bounds(256,4)
#else
! Shared memory (207 KB) allows 1 block per SM. 640 threads (96 registers): +6% over 512.
#define HYDRO_THREADS 640
#define HYDRO_LB launch_bounds(640,1)
#endif
#endif
#endif

#endif
