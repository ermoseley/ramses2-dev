! Launch configuration of hydro_integrator_kernel in pure-hydro builds, measured with HLLC on a
! uniform periodic 3D Sedov blast at 512^3 on an H200 NVL and checked on an A100. launch_bounds(T,B)
! limits registers to 65536/(T*B).
#ifndef HYDRO_LAUNCH_H
#define HYDRO_LAUNCH_H

#ifndef NSUBGRID
#define NSUBGRID 1
#endif

#ifndef MHD
#if defined(NPRE) && NPRE == 4
! Single precision is fastest without launch_bounds: a bound lets the compiler use 56 registers
! instead of 40, which costs 2-7%.
#define HYDRO_LB
#if NSUBGRID == 1
#define HYDRO_THREADS 64
#elif NSUBGRID == 2
#define HYDRO_THREADS 128
#else
#define HYDRO_THREADS 512
#endif
#else
! Double precision: a 64-register limit fits 16 blocks of 64 threads (nsubgrid=1) or 4 of 256
! (nsubgrid=2) per SM. At nsubgrid=4 shared memory (207 KB) allows one block per SM, so more
! threads hide latency: 640 threads with up to 96 registers.
#if NSUBGRID == 1
#define HYDRO_THREADS 64
#define HYDRO_LB launch_bounds(64,16)
#elif NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_LB launch_bounds(256,4)
#else
#define HYDRO_THREADS 640
#define HYDRO_LB launch_bounds(640,1)
#endif
#endif
#endif

#endif
