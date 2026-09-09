module incomp_step_module
  !---------------------------------------------------------------------------
  ! RAMSES-facing incompressible step (doc/incompressible.md Section 3).
  !
  ! Operates on uold for one level through the Cartesian key, so it is a
  ! *step*, not a fork of the mesh.  Requires levelmin == levelmax, periodic
  ! boundaries, a power-of-two grid and no refinement; anything else is
  ! rejected at startup rather than silently approximated.
  !
  ! rho = rho_0 and p = p_0 are constants of the motion by construction, so
  ! the density and energy slots of uold carry rho_0 and
  ! rho_0|u|^2/2 + 3 p_0/2 exactly.  That is deliberate: condinit, the
  ! snapshot writer and every dfmm diagnostic keep working with no special
  ! case, and an incompressible snapshot is directly comparable with a
  ! compressible one.
  !---------------------------------------------------------------------------
  use incomp_ops_module
  use incomp_solver_module
  implicit none
  private
  public :: incomp_validate, incomp_cmpdt, incomp_step, incomp_nu

contains

  integer function incomp_nside(r)
    use amr_commons, only: run_t
    type(run_t)::r
    incomp_nside = 2**r%levelmin
  end function incomp_nside

  real(kind=8) function incomp_nu(r)
    ! nu = mu / rho_0 with mu = p_0 tau, the SAME calibration the compressible
    ! rungs use (doc/dfmm_3d.md Section 4), so that a given dfmm_tau means the
    ! same physical viscosity in every rung.  That is the whole point of the
    ! four-rung design.
    use amr_commons, only: run_t
    type(run_t)::r
    if(r%dfmm_tau>0.0d0)then
       incomp_nu = r%incomp_p0*r%dfmm_tau/r%incomp_rho0
    else
       incomp_nu = 0.0d0
    endif
  end function incomp_nu

  subroutine incomp_validate(r)
    use amr_commons, only: run_t
    use incomp_fft_module, only: is_pow2
    type(run_t)::r
    integer::n
    if(.not.r%incompressible) return
    if(r%levelmin/=r%nlevelmax)then
       write(*,*)'incompressible requires levelmin = levelmax; got ', &
            r%levelmin, r%nlevelmax
       stop 1
    endif
    n = incomp_nside(r)
    if(.not.is_pow2(n))then
       write(*,*)'incompressible requires a power-of-two grid; got n = ',n
       stop 1
    endif
    if(r%incomp_p0<=0.0d0 .or. r%incomp_rho0<=0.0d0)then
       write(*,*)'incompressible requires incomp_p0 > 0 and incomp_rho0 > 0'
       stop 1
    endif
    if(trim(r%incomp_stress)/='viscous' .and. trim(r%incomp_stress)/='moment')then
       write(*,*)"incomp_stress must be 'viscous' or 'moment'; got ", &
            trim(r%incomp_stress)
       stop 1
    endif
    write(*,'(" incompressible: n=",I5,"  p0=",1pe10.3,"  rho0=",1pe10.3, &
         & "  nu=",1pe10.3,"  stress=",A)') &
         n, r%incomp_p0, r%incomp_rho0, incomp_nu(r), trim(r%incomp_stress)
  end subroutine incomp_validate

  subroutine incomp_gather(r,m,ilevel,n,u,pi5,want_pi)
    ! uold -> uniform lattice, by Cartesian key.  For levelmin == levelmax the
    ! key is a bijection onto the lattice, so this is exact; the caller checks
    ! that every site was written.
    use amr_commons, only: run_t, mesh_t
    use amr_parameters, only: ndim, twotondim
    use hydro_parameters, only: ipi
    type(run_t)::r
    type(mesh_t)::m
    integer,intent(in)::ilevel,n
    real(kind=8),intent(out)::u(3,n,n,n)
    real(kind=8),intent(out)::pi5(5,n,n,n)
    logical,intent(in)::want_pi
    integer::igrid,ind,idim,ic(3),nstride,i,j,k,iv
    integer::nfilled
    real(kind=8)::irho

    u = 0.0d0; pi5 = 0.0d0
    nfilled = 0
    irho = 1.0d0/r%incomp_rho0
    do igrid=m%head(ilevel),m%tail(ilevel)
       do ind=1,twotondim
          do idim=1,ndim
             nstride = 2**(idim-1)
             ic(idim) = 2*m%grid(igrid)%ckey(idim) + MOD((ind-1)/nstride,2)
          end do
          i = ic(1)+1; j = ic(2)+1; k = ic(3)+1
          if(i<1.or.i>n.or.j<1.or.j>n.or.k<1.or.k>n)then
             write(*,*)'incomp_gather: cell outside the lattice',i,j,k,n
             stop 1
          endif
          u(1,i,j,k) = m%uold(ind,2,igrid)*irho
          u(2,i,j,k) = m%uold(ind,3,igrid)*irho
          u(3,i,j,k) = m%uold(ind,4,igrid)*irho
          if(want_pi)then
             do iv=1,5
                pi5(iv,i,j,k) = m%uold(ind,ipi+iv-1,igrid)
             end do
          endif
          nfilled = nfilled+1
       end do
    end do
    if(nfilled/=n**3)then
       write(*,*)'incomp_gather: filled ',nfilled,' of ',n**3, &
            ' lattice sites -- the level is not a complete uniform grid'
       stop 1
    endif
  end subroutine incomp_gather

  subroutine incomp_scatter(r,m,ilevel,n,u)
    ! Uniform lattice -> unew, restoring the compressible slots to their
    ! incompressible values so that every downstream tool keeps working.
    !
    ! Writes unew, not uold: the step runs where the hyperbolic solver would,
    ! and r_set_uold copies unew -> uold immediately afterwards, so a write to
    ! uold would be clobbered.
    use amr_commons, only: run_t, mesh_t
    use amr_parameters, only: ndim, twotondim
    type(run_t)::r
    type(mesh_t)::m
    integer,intent(in)::ilevel,n
    real(kind=8),intent(in)::u(3,n,n,n)
    integer::igrid,ind,idim,ic(3),nstride,i,j,k,iv
    real(kind=8)::rho0,p0,ek

    rho0 = r%incomp_rho0
    p0   = r%incomp_p0
    do igrid=m%head(ilevel),m%tail(ilevel)
       do ind=1,twotondim
          do idim=1,ndim
             nstride = 2**(idim-1)
             ic(idim) = 2*m%grid(igrid)%ckey(idim) + MOD((ind-1)/nstride,2)
          end do
          i = ic(1)+1; j = ic(2)+1; k = ic(3)+1
          ek = 0.5d0*rho0*(u(1,i,j,k)**2+u(2,i,j,k)**2+u(3,i,j,k)**2)
          m%unew(ind,1,igrid) = rho0
          m%unew(ind,2,igrid) = rho0*u(1,i,j,k)
          m%unew(ind,3,igrid) = rho0*u(2,i,j,k)
          m%unew(ind,4,igrid) = rho0*u(3,i,j,k)
          m%unew(ind,5,igrid) = ek + 1.5d0*p0
          ! Everything above ivar 5 is carried through unchanged; the moment
          ! sector is advanced separately.
          do iv=6,size(m%unew,2)
             m%unew(ind,iv,igrid) = m%uold(ind,iv,igrid)
          end do
       end do
    end do
  end subroutine incomp_scatter

  subroutine incomp_cmpdt(r,g,m,ilevel,mass,ekin,eint,eani,dt)
    ! Advective and viscous timestep, plus the conserved sums in the slots
    ! update_time.f90 expects.  The ekin slot carries the TOTAL energy, as in
    ! the compressible path.
    use amr_commons, only: run_t, global_t, mesh_t
    type(run_t)::r
    type(global_t)::g
    type(mesh_t)::m
    integer,intent(in)::ilevel
    real(kind=8),intent(out)::mass,ekin,eint,eani,dt
    integer::n
    real(kind=8),allocatable::u(:,:,:,:),pi5(:,:,:,:)
    real(kind=8)::vol,e

    n = incomp_nside(r)
    allocate(u(3,n,n,n),pi5(5,n,n,n))
    call incomp_gather(r,m,ilevel,n,u,pi5,.false.)
    ! Note: the caller is responsible for the host copy being current.  On the
    ! Metal path r_courant_fine runs after r_set_uold, which has just been
    ! followed by incomp_step's own device upload, and the timestep only needs
    ! max|u|, so a stale copy would show up immediately as a wrong dt.
    vol  = (r%boxlen/dble(n))**3
    e    = incomp_energy(u,n)*r%incomp_rho0
    mass = r%incomp_rho0*r%boxlen**3
    eint = 1.5d0*r%incomp_p0*r%boxlen**3
    ekin = eint + e*r%boxlen**3
    eani = 0.0d0
    dt   = incomp_dt(u,n,r%boxlen,incomp_nu(r),r%courant_factor)
    deallocate(u,pi5)
  end subroutine incomp_cmpdt

  subroutine incomp_step(sim,ilevel,dt)
    use ramses_commons, only: ramses_t
#ifdef _METAL
    use metal_runner, only: metal_uold_to_host, metal_unew_to_device
#endif
    type(ramses_t)::sim
    integer,intent(in)::ilevel
    real(kind=8),intent(in)::dt
    integer::n
    logical::use_pi
    real(kind=8),allocatable::u(:,:,:,:),pi5(:,:,:,:)
    real(kind=8)::dv,ef

    n = incomp_nside(sim%r)
    use_pi = (trim(sim%r%incomp_stress)=='moment')
    allocate(u(3,n,n,n),pi5(5,n,n,n))
#ifdef _METAL
    ! The device holds the state between steps; bring it home first.
    call metal_uold_to_host(sim)
#endif
    call incomp_gather(sim%r,sim%m,ilevel,n,u,pi5,use_pi)
    call incomp_step_rk2(u,n,sim%r%boxlen,dt,incomp_nu(sim%r),pi5, &
         sim%r%incomp_rho0,use_pi)
    call incomp_scatter(sim%r,sim%m,ilevel,n,u)
#ifdef _METAL
    call metal_unew_to_device(sim)
#endif
    if(sim%r%incomp_diag)then
       dv = incomp_div_max(u,n,sim%r%boxlen)
       ef = incomp_dealias_frac(u,n,sim%r%boxlen)
       write(*,'(" incomp level=",I2,"  spectral div ",1pe10.3, &
            & "  E_trunc/E ",1pe10.3,"  max|u| ",1pe10.3)') &
            ilevel, dv, ef, maxval(abs(u))
    endif
    deallocate(u,pi5)
  end subroutine incomp_step

end module incomp_step_module
