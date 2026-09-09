module incomp_solver_module
  !---------------------------------------------------------------------------
  ! Fractional-step incompressible solver core (doc/incompressible.md Sec. 2).
  !
  ! One explicit RK2 (Heun) step with a Leray projection inside each stage:
  !
  !     a(u)  = P[ u x omega + stress ]
  !     u_1   = P[ u^n + dt a(u^n) ]
  !     u^n+1 = P[ u^n + (dt/2)( a(u^n) + a(u_1) ) ]
  !
  ! The nonlinear term is in rotational form.  Since
  ! (u.grad)u = grad(|u|^2/2) - u x omega and the projector annihilates
  ! gradients, P[-(u.grad)u] = P[u x omega], and the discrete term then
  ! conserves kinetic energy exactly in the inviscid unaliased limit.  For a
  ! blowup study that matters: any energy change is physical or from the
  ! dealiasing truncation, never from the advection.
  !
  ! `stress` is supplied by the caller so the same core serves both rungs:
  !   rung 1 (incompressible Navier-Stokes)  stress = nu lap u
  !   rung 4 (incompressible dfmm)           stress = -(1/rho_0) div Pi
  !---------------------------------------------------------------------------
  use incomp_ops_module
  implicit none
  private
  public :: incomp_rhs, incomp_step_rk2, incomp_dt, incomp_divpi

contains

  subroutine incomp_rhs(u,a,n,boxlen,nu,pi5,rho0,use_pi)
    ! a = P[ u x omega + stress ], dealiased.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen,nu,rho0
    real(kind=8),intent(in)::u(3,n,n,n)
    real(kind=8),intent(in)::pi5(5,n,n,n)
    logical,intent(in)::use_pi
    real(kind=8),intent(out)::a(3,n,n,n)
    real(kind=8),allocatable::w(:,:,:,:),s(:,:,:,:)
    integer::i,j,k

    allocate(w(3,n,n,n))
    call incomp_curl(u,w,n,boxlen)

    ! u x omega
    do k=1,n
       do j=1,n
          do i=1,n
             a(1,i,j,k) = u(2,i,j,k)*w(3,i,j,k) - u(3,i,j,k)*w(2,i,j,k)
             a(2,i,j,k) = u(3,i,j,k)*w(1,i,j,k) - u(1,i,j,k)*w(3,i,j,k)
             a(3,i,j,k) = u(1,i,j,k)*w(2,i,j,k) - u(2,i,j,k)*w(1,i,j,k)
          end do
       end do
    end do
    deallocate(w)

    ! The nonlinear product is the only place aliasing is generated, so the
    ! 2/3 truncation is applied here and nowhere else.
    call incomp_dealias(a,n,boxlen)

    allocate(s(3,n,n,n))
    if(use_pi)then
       call incomp_divpi(pi5,s,n,boxlen)
       a = a - s/rho0
    else if(nu>0.0d0)then
       call incomp_lap(u,s,n,boxlen,3)
       a = a + nu*s
    endif
    deallocate(s)

    call incomp_project(a,n,boxlen)

  end subroutine incomp_rhs

  subroutine incomp_divpi(pi5,s,n,boxlen)
    ! s_i = d_j Pi_ij, spectrally, from the five stored deviatoric components
    ! in the dfmm packing (xx, yy, xy, xz, yz) with Pi_zz = -(Pi_xx + Pi_yy).
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::pi5(5,n,n,n)
    real(kind=8),intent(out)::s(3,n,n,n)
    complex(kind=8),allocatable::ph(:,:,:,:),sh(:,:,:,:)
    real(kind=8)::kv(n),kx,ky,kz
    complex(kind=8)::ii,pxx,pyy,pzz,pxy,pxz,pyz
    integer::i,j,k,d

    ii = cmplx(0.0d0,1.0d0,kind=8)
    call incomp_wavenumbers(n,boxlen,kv)
    allocate(ph(n,n,n,5),sh(n,n,n,3))
    do d=1,5
       ph(:,:,:,d) = cmplx(pi5(d,:,:,:),0.0d0,kind=8)
       call fft3d_wrap(ph(:,:,:,d),n,-1)
    end do
    do k=1,n
       kz=kv(k)
       do j=1,n
          ky=kv(j)
          do i=1,n
             kx=kv(i)
             pxx=ph(i,j,k,1); pyy=ph(i,j,k,2)
             pxy=ph(i,j,k,3); pxz=ph(i,j,k,4); pyz=ph(i,j,k,5)
             pzz=-(pxx+pyy)
             sh(i,j,k,1) = ii*(kx*pxx + ky*pxy + kz*pxz)
             sh(i,j,k,2) = ii*(kx*pxy + ky*pyy + kz*pyz)
             sh(i,j,k,3) = ii*(kx*pxz + ky*pyz + kz*pzz)
          end do
       end do
    end do
    do d=1,3
       call fft3d_wrap(sh(:,:,:,d),n,1)
       s(d,:,:,:) = dble(sh(:,:,:,d))
    end do
    deallocate(ph,sh)
  end subroutine incomp_divpi

  subroutine fft3d_wrap(a,n,isign)
    use incomp_fft_module, only: fft3d
    integer,intent(in)::n,isign
    complex(kind=8),intent(inout)::a(n,n,n)
    call fft3d(a,n,isign)
  end subroutine fft3d_wrap

  subroutine incomp_step_rk2(u,n,boxlen,dt,nu,pi5,rho0,use_pi)
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen,dt,nu,rho0
    real(kind=8),intent(inout)::u(3,n,n,n)
    real(kind=8),intent(in)::pi5(5,n,n,n)
    logical,intent(in)::use_pi
    real(kind=8),allocatable::a0(:,:,:,:),a1(:,:,:,:),u1(:,:,:,:)

    allocate(a0(3,n,n,n),a1(3,n,n,n),u1(3,n,n,n))
    call incomp_rhs(u,a0,n,boxlen,nu,pi5,rho0,use_pi)
    u1 = u + dt*a0
    call incomp_project(u1,n,boxlen)
    call incomp_rhs(u1,a1,n,boxlen,nu,pi5,rho0,use_pi)
    u = u + 0.5d0*dt*(a0+a1)
    call incomp_project(u,n,boxlen)
    deallocate(a0,a1,u1)
  end subroutine incomp_step_rk2

  real(kind=8) function incomp_dt(u,n,boxlen,nu,courant,cmom)
    ! Three limbs: advective always, viscous for rung 1, and -- when the third
    ! moment is evolved -- the hyperbolic speed of the Pi <-> Q pair.
    !
    ! That last limb corrects a claim this function used to make.  Rung 4's
    ! *stress relaxation* is stiff and handled by an exact exponential map, so
    ! it carries no parabolic constraint; but once Q is carried, Pi's flux is Q
    ! and Q's production is theta_0 grad Pi, and the pair propagates at the
    ! thermal speed.  Removing compressibility removes the acoustic constraint
    ! from the velocity equation only.  Pass cmom = 0 for rung 1 or for a
    ! ten-moment rung 4, where Pi genuinely only advects.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen,nu,courant,cmom
    real(kind=8),intent(in)::u(3,n,n,n)
    real(kind=8)::dx,umax,dta,dtv
    dx   = boxlen/dble(n)
    umax = max(maxval(abs(u)),1.0d-30)
    dta  = courant*dx/(3.0d0*(umax+max(cmom,0.0d0)))
    if(nu>0.0d0)then
       dtv = courant*dx*dx/(2.0d0*3.0d0*nu)
       incomp_dt = min(dta,dtv)
    else
       incomp_dt = dta
    endif
  end function incomp_dt

end module incomp_solver_module
