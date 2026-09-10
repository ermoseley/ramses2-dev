module incomp_fft_module
  !---------------------------------------------------------------------------
  ! In-place iterative radix-2 complex FFT, and the 3D transform built from it.
  !
  ! Written in-tree rather than linked against FFTW deliberately.  RAMSES
  ! levels are always powers of two, so radix-2 is sufficient; bin/Makefile
  ! reaches FFTW only under TURB=1 and at a hardcoded cluster path, so an
  ! in-tree transform is the more portable of the two and it is verified
  ! against a direct DFT (doc/incompressible.md Gate 1).
  !
  ! Convention:  isign = -1  forward,  sum_x a(x) exp(-2 pi i k x / n)
  !              isign = +1  inverse,  UNNORMALISED
  ! fft3d divides by n^3 on the inverse, so ifft3(fft3(a)) = a.
  !---------------------------------------------------------------------------
  implicit none
  private
  public :: fft1d, fft3d, is_pow2

contains

  logical function is_pow2(n)
    integer,intent(in)::n
    is_pow2 = (n>0) .and. (iand(n,n-1)==0)
  end function is_pow2

  subroutine fft1d(a,n,isign)
    ! Cooley-Tukey, decimation in time, with an explicit bit-reversal
    ! permutation followed by log2(n) butterfly passes.
    integer,intent(in)::n,isign
    complex(kind=8),intent(inout)::a(n)
    integer::i,j,m,mmax,istep
    real(kind=8)::theta,twopi
    complex(kind=8)::w,wp,tmp

    twopi = 8.0d0*atan(1.0d0)

    ! Bit-reversal permutation
    j = 1
    do i=1,n-1
       if(j>i)then
          tmp=a(j); a(j)=a(i); a(i)=tmp
       endif
       m = n/2
       do while(m>=1 .and. j>m)
          j = j-m
          m = m/2
       end do
       j = j+m
    end do

    ! Butterflies
    mmax = 1
    do while(mmax<n)
       istep = 2*mmax
       theta = dble(isign)*twopi/dble(istep)
       do m=1,mmax
          w = cmplx(cos(theta*dble(m-1)),sin(theta*dble(m-1)),kind=8)
          do i=m,n,istep
             j = i+mmax
             tmp  = w*a(j)
             a(j) = a(i)-tmp
             a(i) = a(i)+tmp
          end do
       end do
       mmax = istep
    end do

  end subroutine fft1d

  subroutine fft3d(a,n,isign)
    ! Three passes of n^2 one-dimensional transforms.  a is indexed (i,j,k)
    ! with i fastest, matching Fortran storage, so the x pass is contiguous.
    integer,intent(in)::n,isign
    complex(kind=8),intent(inout)::a(n,n,n)
    integer::i,j,k
    complex(kind=8),allocatable::line(:)
    real(kind=8)::scal

    allocate(line(n))

    do k=1,n
       do j=1,n
          line = a(:,j,k)
          call fft1d(line,n,isign)
          a(:,j,k) = line
       end do
    end do
    do k=1,n
       do i=1,n
          line = a(i,:,k)
          call fft1d(line,n,isign)
          a(i,:,k) = line
       end do
    end do
    do j=1,n
       do i=1,n
          line = a(i,j,:)
          call fft1d(line,n,isign)
          a(i,j,:) = line
       end do
    end do

    if(isign>0)then
       scal = 1.0d0/(dble(n)**3)
       a = a*scal
    endif

    deallocate(line)

  end subroutine fft3d

end module incomp_fft_module
