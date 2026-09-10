module hydro_parameters
  use amr_parameters, only: ndim

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=5+nener
#else
  integer,parameter::nvar=NVAR
#endif

#ifdef MHD
  integer,parameter::nprim=NVAR+3
  integer,parameter::ie=8
#else
  integer,parameter::nprim=NVAR
  integer,parameter::ie=5
#endif

  ! dfmm fields (doc/dfmm_3d.md Section 2).  They occupy ivar 6 .. 5+ndfmm,
  ! which the Metal kernel hardwires, and split into two kinds:
  !
  !   density-like (ndfmm_dens fields, ivar 6 .. 5+ndfmm_dens)
  !       transport flux u_k X, so the stored value is X itself and it must
  !       NOT be multiplied by rho when converting primitive <-> conserved.
  !         Pi_ij at ivar  6..10   (Stage 1, ndfmm =  5)
  !         Q_ijk at ivar 11..20   (Stage 2, ndfmm = 15)
  !
  !   mass-like (ndfmm_mass fields, ivar 6+ndfmm_dens .. 5+ndfmm)
  !       transport flux u_k rho X, i.e. RAMSES's passive-scalar convention,
  !       so the stored value is rho X and it IS multiplied by rho.
  !         rho L_i    at ivar 21..23   (Stage 3, ndfmm = 18)
  !         rho Sxx_ij at ivar 24..29   (Stage 4, ndfmm = 33)
  !         rho Sxv_ij at ivar 30..38
  !
  ! D L_i / Dt = 0 is what forces the mass-like form for the tower: a
  ! density-like field with no source obeys D X / Dt = -X div u instead.
#ifdef DFMM
  integer,parameter::ndfmm=NDFMM
  integer,parameter::ipi=6
#if NDFMM>15
  integer,parameter::ndfmm_dens=15
#else
  integer,parameter::ndfmm_dens=NDFMM
#endif
  integer,parameter::ndfmm_mass=ndfmm-ndfmm_dens
#if NDFMM>=15
  integer,parameter::iq=11
#else
  integer,parameter::iq=0
#endif
#if NDFMM>=18
  integer,parameter::il=21
#else
  integer,parameter::il=0
#endif
#if NDFMM>=33
  integer,parameter::isxx=24
  integer,parameter::isxv=30
#else
  integer,parameter::isxx=0
  integer,parameter::isxv=0
#endif
#else
  integer,parameter::ndfmm=0
  integer,parameter::ndfmm_dens=0
  integer,parameter::ndfmm_mass=0
  integer,parameter::ipi=0
  integer,parameter::iq=0
  integer,parameter::il=0
  integer,parameter::isxx=0
  integer,parameter::isxv=0
#endif

#ifdef NION
  integer,parameter::nion=NION  ! # of ionization fractions species
#else
  integer,parameter::nion=1
#endif

  integer,parameter::solver_llf=1
  integer,parameter::solver_hll=2
  integer,parameter::solver_hllc=3
  integer,parameter::solver_hlld=4
  integer,parameter::solver_roe=5
  integer,parameter::solver_upwind=6
#ifdef _METAL
  integer,parameter::solver_uct_hlld=7
#endif

  integer,parameter::solver2d_llf=1
  integer,parameter::solver2d_hllf=2
  integer,parameter::solver2d_hlla=3
  integer,parameter::solver2d_hlld=4
  integer,parameter::solver2d_roe=5
  integer,parameter::solver2d_upwind=6

end module hydro_parameters

module const

  ! Some useful constant
  real(kind=8),parameter ::bigreal = 1.0d+30
  real(kind=8),parameter ::zero = 0.0
  real(kind=8),parameter ::one = 1.0
  real(kind=8),parameter ::two = 2.0
  real(kind=8),parameter ::three = 3.0
  real(kind=8),parameter ::four = 4.0
  real(kind=8),parameter ::two3rd = 0.6666666666666667
  real(kind=8),parameter ::half = 0.5
  real(kind=8),parameter ::third = 0.33333333333333333
  real(kind=8),parameter ::forth = 0.25
  real(kind=8),parameter ::sixth = 0.16666666666666667

end module const
