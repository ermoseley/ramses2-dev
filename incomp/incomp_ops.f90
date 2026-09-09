module incomp_ops_module
  !---------------------------------------------------------------------------
  ! Spectral operators for the incompressible rungs (doc/incompressible.md
  ! Section 2): the Leray projector, the curl, the Laplacian, the 2/3-rule
  ! dealiasing truncation, and the divergence diagnostic.
  !
  ! Fields are real(8) u(3, n, n, n) in physical space.  Every routine that
  ! needs spectral space transforms, acts, and transforms back, so the caller
  ! never sees Fourier coefficients.  That costs transforms but it keeps the
  ! solver readable and it is not the bottleneck: the nonlinear term dominates.
  !---------------------------------------------------------------------------
  use incomp_fft_module, only: fft3d, is_pow2
  implicit none
  private
  public :: incomp_wavenumbers, incomp_project, incomp_curl, incomp_lap
  public :: incomp_div_max, incomp_dealias_frac, incomp_energy, incomp_dealias

contains

  subroutine incomp_wavenumbers(n,boxlen,kv)
    ! kv(j) is the signed wavenumber of index j: 2 pi/L * (j-1) for the first
    ! half, 2 pi/L * (j-1-n) for the second.  The Nyquist mode j = n/2+1 is
    ! given its negative alias, and is zeroed by every operator below that
    ! would otherwise make a real field complex.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(out)::kv(n)
    integer::j
    real(kind=8)::twopi
    twopi = 8.0d0*atan(1.0d0)
    do j=1,n
       if(j-1 <= n/2)then
          kv(j) = twopi/boxlen*dble(j-1)
       else
          kv(j) = twopi/boxlen*dble(j-1-n)
       endif
    end do
    kv(n/2+1) = 0.0d0   ! Nyquist: no signed derivative exists for a real field
  end subroutine incomp_wavenumbers

  subroutine incomp_project(u,n,boxlen)
    ! Leray projection:  u_hat <- (I - k k^T / |k|^2) u_hat.
    ! Exact by construction, so k . u_hat = 0 to machine precision with no
    ! iteration and no tolerance.  The k = 0 mode (the mean flow) is left
    ! alone: it is already divergence free and removing it would not conserve
    ! momentum.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(inout)::u(3,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:,:)
    real(kind=8)::kv(n),kx,ky,kz,k2,ik2
    complex(kind=8)::kdotu
    integer::i,j,k,d

    call incomp_wavenumbers(n,boxlen,kv)
    allocate(uh(n,n,n,3))
    do d=1,3
       uh(:,:,:,d) = cmplx(u(d,:,:,:),0.0d0,kind=8)
       call fft3d(uh(:,:,:,d),n,-1)
    end do

    do k=1,n
       kz = kv(k)
       do j=1,n
          ky = kv(j)
          do i=1,n
             kx = kv(i)
             k2 = kx*kx+ky*ky+kz*kz
             if(k2<=0.0d0) cycle
             ik2 = 1.0d0/k2
             kdotu = kx*uh(i,j,k,1)+ky*uh(i,j,k,2)+kz*uh(i,j,k,3)
             uh(i,j,k,1) = uh(i,j,k,1) - kx*kdotu*ik2
             uh(i,j,k,2) = uh(i,j,k,2) - ky*kdotu*ik2
             uh(i,j,k,3) = uh(i,j,k,3) - kz*kdotu*ik2
          end do
       end do
    end do

    do d=1,3
       call fft3d(uh(:,:,:,d),n,1)
       u(d,:,:,:) = dble(uh(:,:,:,d))
    end do
    deallocate(uh)
  end subroutine incomp_project

  subroutine incomp_curl(u,w,n,boxlen)
    ! w = curl u, spectrally.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(3,n,n,n)
    real(kind=8),intent(out)::w(3,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:,:),wh(:,:,:,:)
    real(kind=8)::kv(n),kx,ky,kz
    complex(kind=8)::ii
    integer::i,j,k,d

    ii = cmplx(0.0d0,1.0d0,kind=8)
    call incomp_wavenumbers(n,boxlen,kv)
    allocate(uh(n,n,n,3),wh(n,n,n,3))
    do d=1,3
       uh(:,:,:,d) = cmplx(u(d,:,:,:),0.0d0,kind=8)
       call fft3d(uh(:,:,:,d),n,-1)
    end do
    do k=1,n
       kz=kv(k)
       do j=1,n
          ky=kv(j)
          do i=1,n
             kx=kv(i)
             wh(i,j,k,1) = ii*(ky*uh(i,j,k,3)-kz*uh(i,j,k,2))
             wh(i,j,k,2) = ii*(kz*uh(i,j,k,1)-kx*uh(i,j,k,3))
             wh(i,j,k,3) = ii*(kx*uh(i,j,k,2)-ky*uh(i,j,k,1))
          end do
       end do
    end do
    do d=1,3
       call fft3d(wh(:,:,:,d),n,1)
       w(d,:,:,:) = dble(wh(:,:,:,d))
    end do
    deallocate(uh,wh)
  end subroutine incomp_curl

  subroutine incomp_lap(u,l,n,boxlen,ncomp)
    ! l = laplacian(u), spectrally, for ncomp components.
    integer,intent(in)::n,ncomp
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(ncomp,n,n,n)
    real(kind=8),intent(out)::l(ncomp,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:)
    real(kind=8)::kv(n),k2
    integer::i,j,k,d
    call incomp_wavenumbers(n,boxlen,kv)
    allocate(uh(n,n,n))
    do d=1,ncomp
       uh = cmplx(u(d,:,:,:),0.0d0,kind=8)
       call fft3d(uh,n,-1)
       do k=1,n
          do j=1,n
             do i=1,n
                k2 = kv(i)**2+kv(j)**2+kv(k)**2
                uh(i,j,k) = -k2*uh(i,j,k)
             end do
          end do
       end do
       call fft3d(uh,n,1)
       l(d,:,:,:) = dble(uh)
    end do
    deallocate(uh)
  end subroutine incomp_lap

  subroutine incomp_dealias(a,n,boxlen)
    ! Orszag's 2/3 rule, applied in place to a physical-space field: zero
    ! every mode with |k_d| > (2/3) k_max in any direction.  A projection, not
    ! a dissipation -- idempotent, and it removes nothing from the retained
    ! band of a resolved field.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(inout)::a(3,n,n,n)
    complex(kind=8),allocatable::ah(:,:,:)
    real(kind=8)::kv(n),kcut,kmax
    integer::i,j,k,d
    call incomp_wavenumbers(n,boxlen,kv)
    kmax = maxval(abs(kv))
    kcut = 2.0d0/3.0d0*kmax
    allocate(ah(n,n,n))
    do d=1,3
       ah = cmplx(a(d,:,:,:),0.0d0,kind=8)
       call fft3d(ah,n,-1)
       do k=1,n
          do j=1,n
             do i=1,n
                if(abs(kv(i))>kcut .or. abs(kv(j))>kcut .or. abs(kv(k))>kcut) &
                     ah(i,j,k) = cmplx(0.0d0,0.0d0,kind=8)
             end do
          end do
       end do
       call fft3d(ah,n,1)
       a(d,:,:,:) = dble(ah)
    end do
    deallocate(ah)
  end subroutine incomp_dealias

  real(kind=8) function incomp_div_max(u,n,boxlen)
    ! max |k . u_hat| normalised by max(|k| |u_hat|), the dimensionless
    ! spectral divergence.  This is the quantity the rung claims to hold at
    ! machine precision.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(3,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:,:)
    real(kind=8)::kv(n),kx,ky,kz,kmag,umag,dmax,scal
    complex(kind=8)::kdotu
    integer::i,j,k,d
    call incomp_wavenumbers(n,boxlen,kv)
    allocate(uh(n,n,n,3))
    do d=1,3
       uh(:,:,:,d) = cmplx(u(d,:,:,:),0.0d0,kind=8)
       call fft3d(uh(:,:,:,d),n,-1)
    end do
    dmax=0.0d0; scal=0.0d0
    do k=1,n
       kz=kv(k)
       do j=1,n
          ky=kv(j)
          do i=1,n
             kx=kv(i)
             kmag = sqrt(kx*kx+ky*ky+kz*kz)
             kdotu = kx*uh(i,j,k,1)+ky*uh(i,j,k,2)+kz*uh(i,j,k,3)
             umag = sqrt(abs(uh(i,j,k,1))**2+abs(uh(i,j,k,2))**2+abs(uh(i,j,k,3))**2)
             dmax = max(dmax,abs(kdotu))
             scal = max(scal,kmag*umag)
          end do
       end do
    end do
    deallocate(uh)
    if(scal>0.0d0)then
       incomp_div_max = dmax/scal
    else
       incomp_div_max = 0.0d0
    endif
  end function incomp_div_max

  real(kind=8) function incomp_dealias_frac(u,n,boxlen)
    ! Fraction of the kinetic energy sitting in the band the 2/3 rule
    ! truncates.  Reported so that an under-resolved run is visible rather
    ! than silently damped.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(3,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:,:)
    real(kind=8)::kv(n),kcut,kmax,etot,ecut,e
    integer::i,j,k,d
    call incomp_wavenumbers(n,boxlen,kv)
    kmax = maxval(abs(kv)); kcut = 2.0d0/3.0d0*kmax
    allocate(uh(n,n,n,3))
    do d=1,3
       uh(:,:,:,d) = cmplx(u(d,:,:,:),0.0d0,kind=8)
       call fft3d(uh(:,:,:,d),n,-1)
    end do
    etot=0.0d0; ecut=0.0d0
    do k=1,n
       do j=1,n
          do i=1,n
             e = abs(uh(i,j,k,1))**2+abs(uh(i,j,k,2))**2+abs(uh(i,j,k,3))**2
             etot = etot+e
             if(abs(kv(i))>kcut .or. abs(kv(j))>kcut .or. abs(kv(k))>kcut) ecut=ecut+e
          end do
       end do
    end do
    deallocate(uh)
    if(etot>0.0d0)then
       incomp_dealias_frac = ecut/etot
    else
       incomp_dealias_frac = 0.0d0
    endif
  end function incomp_dealias_frac

  real(kind=8) function incomp_energy(u,n)
    ! Volume-averaged kinetic energy per unit mass, |u|^2/2.
    integer,intent(in)::n
    real(kind=8),intent(in)::u(3,n,n,n)
    incomp_energy = 0.5d0*sum(u*u)/dble(n)**3
  end function incomp_energy

end module incomp_ops_module
