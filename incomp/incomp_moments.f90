module incomp_moments_module
  !---------------------------------------------------------------------------
  ! Moment sector of the incompressible dfmm rung (doc/incompressible.md).
  !
  ! With div u = 0 and rho = rho_0, p = p_0 constant, the compressible
  ! ten-moment system of doc/dfmm_3d.md Section 4 collapses to
  !
  !     D Pi_ij / Dt = -2 p_0 S_ij - [Pi_ik G_jk + Pi_jk G_ik]^dev - Pi_ij/tau
  !
  ! with two exact simplifications that are worth naming, because they are
  ! what make this rung clean rather than merely different:
  !
  !   * S0 = S exactly.  The deviatoric projection of the strain subtracts
  !     (1/3) delta_ij div u, which vanishes identically, so there is no trace
  !     removal and no associated truncation error.
  !   * Pi's transport is pure advection.  A density-like field obeys
  !     D X/Dt = -X div u, so with div u = 0 the density-like and mass-like
  !     forms coincide and Pi advects like a scalar.  In the compressible code
  !     Pi's flux additionally carries Q_ijk - (2/3) delta_ij q_k; here q is
  !     identically zero (see below), so nothing is dropped.
  !
  ! **Q_ijk is deliberately absent, on physical grounds rather than as a scope
  ! cut.**  theta = p_0/rho_0 is a constant of the motion in this rung, so
  ! grad theta = 0 and the Chapman-Enskog heat flux q_CE = -(5/2) tau_q p
  ! grad theta vanishes identically.  The paper's first indicator,
  ! ||q - q_CE||, therefore has nothing to measure in an incompressible
  ! isothermal setting; the two indicators that *do* have content here are
  ! check 7, lam_min(p_0 I + Pi)/p_0, and the phase-space rank collapse below,
  ! and both are carried.
  !
  ! The second frame carries over unchanged from doc/dfmm_3d.md Section 4:
  !     D D_i / Dt   = -u_i
  !     D Sxx / Dt   = Sxv + Sxv^T
  !     D Sxv / Dt   = Svv - Sxv G^T - Sxv/tau ,     Svv = (p_0 I + Pi)/rho_0
  ! with the same trapezoidal pairing of Sxx to Sxv that the compressible
  ! code needs, for the same reason (a forward-Euler Sxx lags by n dt^2 and
  ! drives the Schur complement through zero).
  !
  ! Everything here is float64 on the host, which also makes this rung the
  ! float64 reference that doc/dfmm_3d.md Section 6 records as outstanding for
  ! the Stage-4 phase-space sector.
  !---------------------------------------------------------------------------
  use incomp_fft_module, only: fft3d
  use incomp_ops_module
  implicit none
  private
  public :: incomp_velgrad, incomp_advect, incomp_moment_step
  public :: incomp_cone_min, incomp_rank_min

contains

  subroutine incomp_velgrad(u,gr,n,boxlen)
    ! gr(a,b,...) = d u_a / d x_b, spectrally.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(3,n,n,n)
    real(kind=8),intent(out)::gr(3,3,n,n,n)
    complex(kind=8),allocatable::uh(:,:,:),th(:,:,:)
    real(kind=8)::kv(n),kd
    complex(kind=8)::ii
    integer::i,j,k,a,b
    ii = cmplx(0.0d0,1.0d0,kind=8)
    call incomp_wavenumbers(n,boxlen,kv)
    allocate(uh(n,n,n),th(n,n,n))
    do a=1,3
       uh = cmplx(u(a,:,:,:),0.0d0,kind=8)
       call fft3d(uh,n,-1)
       do b=1,3
          do k=1,n
             do j=1,n
                do i=1,n
                   if(b==1)then
                      kd=kv(i)
                   else if(b==2)then
                      kd=kv(j)
                   else
                      kd=kv(k)
                   endif
                   th(i,j,k) = ii*kd*uh(i,j,k)
                end do
             end do
          end do
          call fft3d(th,n,1)
          gr(a,b,:,:,:) = dble(th)
       end do
    end do
    deallocate(uh,th)
  end subroutine incomp_velgrad

  subroutine incomp_advect(u,f,adv,n,boxlen,nc)
    ! adv = -u_j d_j f, for nc components of f, spectrally, dealiased.
    integer,intent(in)::n,nc
    real(kind=8),intent(in)::boxlen
    real(kind=8),intent(in)::u(3,n,n,n),f(nc,n,n,n)
    real(kind=8),intent(out)::adv(nc,n,n,n)
    complex(kind=8),allocatable::fh(:,:,:),th(:,:,:)
    real(kind=8),allocatable::df(:,:,:)
    real(kind=8)::kv(n),kd,kcut,kmax
    complex(kind=8)::ii
    integer::i,j,k,c,b
    ii = cmplx(0.0d0,1.0d0,kind=8)
    call incomp_wavenumbers(n,boxlen,kv)
    kmax = maxval(abs(kv)); kcut = 2.0d0/3.0d0*kmax
    allocate(fh(n,n,n),th(n,n,n),df(n,n,n))
    adv = 0.0d0
    do c=1,nc
       fh = cmplx(f(c,:,:,:),0.0d0,kind=8)
       call fft3d(fh,n,-1)
       ! Truncate the field before differentiating it, so the product below
       ! generates no aliased content in the retained band.
       do k=1,n
          do j=1,n
             do i=1,n
                if(abs(kv(i))>kcut.or.abs(kv(j))>kcut.or.abs(kv(k))>kcut) &
                     fh(i,j,k)=cmplx(0.0d0,0.0d0,kind=8)
             end do
          end do
       end do
       do b=1,3
          do k=1,n
             do j=1,n
                do i=1,n
                   if(b==1)then
                      kd=kv(i)
                   else if(b==2)then
                      kd=kv(j)
                   else
                      kd=kv(k)
                   endif
                   th(i,j,k) = ii*kd*fh(i,j,k)
                end do
             end do
          end do
          call fft3d(th,n,1)
          df = dble(th)
          adv(c,:,:,:) = adv(c,:,:,:) - u(b,:,:,:)*df
       end do
    end do
    deallocate(fh,th,df)
  end subroutine incomp_advect

  subroutine incomp_moment_step(u,pi5,dsp,sxx,sxv,n,boxlen,dt,p0,rho0,tau, &
                                have_frame)
    ! One operator-split step of the moment sector, applied after the velocity
    ! has been advanced: advection is explicit and spectral, the stiff
    ! relaxation uses the same asymptotic-preserving exponential map as the
    ! compressible code, and Sxx is paired to Sxv by the trapezoid.
    integer,intent(in)::n
    real(kind=8),intent(in)::boxlen,dt,p0,rho0,tau
    real(kind=8),intent(in)::u(3,n,n,n)
    real(kind=8),intent(inout)::pi5(5,n,n,n)
    real(kind=8),intent(inout)::dsp(3,n,n,n),sxx(6,n,n,n),sxv(9,n,n,n)
    logical,intent(in)::have_frame
    real(kind=8),allocatable::gr(:,:,:,:,:)
    real(kind=8),allocatable::adv5(:,:,:,:),adv3(:,:,:,:)
    real(kind=8),allocatable::adv6(:,:,:,:),adv9(:,:,:,:)
    real(kind=8)::G(3,3),PI3(3,3),TP(3,3),S,PG,trT
    real(kind=8)::decay,tomd,xr,omd
    real(kind=8)::pio(5),rate(5)
    real(kind=8)::vo(9),vn(9),vr(9),svv(3,3),txx(6)
    integer::i,j,k,a,b,c,m

    ! Relaxation weights.  The series form below the cutoff keeps
    ! tau*(1-decay) -> dt accurately where exp(-x)-1 loses its digits.
    decay = 1.0d0; tomd = dt
    if(tau>0.0d0)then
       xr = dt/tau
       if(xr<1.0d-8)then
          omd = xr*(1.0d0-0.5d0*xr); decay = 1.0d0-omd
       else
          decay = exp(-xr); omd = 1.0d0-decay
       endif
       tomd = tau*omd
    endif

    allocate(gr(3,3,n,n,n))
    call incomp_velgrad(u,gr,n,boxlen)

    allocate(adv5(5,n,n,n))
    call incomp_advect(u,pi5,adv5,n,boxlen,5)
    if(have_frame)then
       allocate(adv3(3,n,n,n),adv6(6,n,n,n),adv9(9,n,n,n))
       call incomp_advect(u,dsp,adv3,n,boxlen,3)
       call incomp_advect(u,sxx,adv6,n,boxlen,6)
       call incomp_advect(u,sxv,adv9,n,boxlen,9)
    endif

    do k=1,n
       do j=1,n
          do i=1,n
             G = gr(:,:,i,j,k)

             ! ---- Pi ----
             pio = pi5(:,i,j,k)
             PI3(1,1)=pio(1); PI3(2,2)=pio(2); PI3(3,3)=-(pio(1)+pio(2))
             PI3(1,2)=pio(3); PI3(2,1)=pio(3)
             PI3(1,3)=pio(4); PI3(3,1)=pio(4)
             PI3(2,3)=pio(5); PI3(3,2)=pio(5)
             do a=1,3
                do b=1,3
                   ! div u = 0, so S0 = S exactly: no trace removal needed.
                   S  = 0.5d0*(G(a,b)+G(b,a))
                   PG = 0.0d0
                   do c=1,3
                      PG = PG + PI3(a,c)*G(b,c) + PI3(b,c)*G(a,c)
                   end do
                   TP(a,b) = -2.0d0*p0*S - PG
                end do
             end do
             trT = (TP(1,1)+TP(2,2)+TP(3,3))/3.0d0
             TP(1,1)=TP(1,1)-trT; TP(2,2)=TP(2,2)-trT; TP(3,3)=TP(3,3)-trT
             rate(1)=adv5(1,i,j,k)+TP(1,1)
             rate(2)=adv5(2,i,j,k)+TP(2,2)
             rate(3)=adv5(3,i,j,k)+TP(1,2)
             rate(4)=adv5(4,i,j,k)+TP(1,3)
             rate(5)=adv5(5,i,j,k)+TP(2,3)
             pi5(:,i,j,k) = pio*decay + tomd*rate

             if(.not.have_frame) cycle

             ! ---- Lagrangian displacement:  D D_i/Dt = -u_i ----
             do a=1,3
                dsp(a,i,j,k) = dsp(a,i,j,k) + dt*(adv3(a,i,j,k)-u(a,i,j,k))
             end do

             ! ---- Sxv:  Svv - Sxv G^T - Sxv/tau, Svv = (p0 I + Pi)/rho0 ----
             ! Svv is formed from the UPDATED Pi, so the two sectors are at the
             ! same time level.
             pio = pi5(:,i,j,k)
             svv(1,1)=(p0+pio(1))/rho0
             svv(2,2)=(p0+pio(2))/rho0
             svv(3,3)=(p0-pio(1)-pio(2))/rho0
             svv(1,2)=pio(3)/rho0; svv(2,1)=svv(1,2)
             svv(1,3)=pio(4)/rho0; svv(3,1)=svv(1,3)
             svv(2,3)=pio(5)/rho0; svv(3,2)=svv(2,3)
             vo = sxv(:,i,j,k)
             do a=1,3
                do b=1,3
                   m = 3*(a-1)+b
                   vr(m) = adv9(m,i,j,k) + svv(a,b)
                   do c=1,3
                      vr(m) = vr(m) - vo(3*(a-1)+c)*G(b,c)
                   end do
                end do
             end do
             vn = vo*decay + tomd*vr
             sxv(:,i,j,k) = vn

             ! ---- Sxx:  Sxv + Sxv^T, trapezoid in Sxv ----
             txx(1) = vo(1)+vn(1)
             txx(2) = vo(5)+vn(5)
             txx(3) = vo(9)+vn(9)
             txx(4) = 0.5d0*(vo(2)+vo(4)+vn(2)+vn(4))
             txx(5) = 0.5d0*(vo(3)+vo(7)+vn(3)+vn(7))
             txx(6) = 0.5d0*(vo(6)+vo(8)+vn(6)+vn(8))
             do m=1,6
                sxx(m,i,j,k) = sxx(m,i,j,k) + dt*(adv6(m,i,j,k)+txx(m))
             end do
          end do
       end do
    end do

    deallocate(gr,adv5)
    if(have_frame) deallocate(adv3,adv6,adv9)

  end subroutine incomp_moment_step

  real(kind=8) function incomp_cone_min(pi5,n,p0)
    ! min over cells of lam_min(p_0 I + Pi)/p_0 -- the blowup note's check 7,
    ! well posed here because p_0 is the kinetic gas pressure and not the
    ! Lagrange multiplier (doc/incompressible.md Section 0).
    integer,intent(in)::n
    real(kind=8),intent(in)::p0
    real(kind=8),intent(in)::pi5(5,n,n,n)
    real(kind=8)::P(3,3),lam,lmin
    integer::i,j,k
    lmin = huge(1.0d0)
    do k=1,n
       do j=1,n
          do i=1,n
             P(1,1)=p0+pi5(1,i,j,k)
             P(2,2)=p0+pi5(2,i,j,k)
             P(3,3)=p0-pi5(1,i,j,k)-pi5(2,i,j,k)
             P(1,2)=pi5(3,i,j,k); P(2,1)=P(1,2)
             P(1,3)=pi5(4,i,j,k); P(3,1)=P(1,3)
             P(2,3)=pi5(5,i,j,k); P(3,2)=P(2,3)
             lam = sym3_lam_min(P)
             lmin = min(lmin,lam/p0)
          end do
       end do
    end do
    incomp_cone_min = lmin
  end function incomp_cone_min

  real(kind=8) function incomp_rank_min(pi5,sxx,sxv,n,p0,rho0,nbad)
    ! min over cells of g = sqrt(max(lam_min(Svv^-1 Gamma),0)) with
    ! Gamma = Svv - Sxv^T Sxx^-1 Sxv, and the count of cells with Gamma
    ! indefinite.  In float64, so it is not subject to the cancellation floor
    ! the float32 Metal diagnostic hits (doc/dfmm_3d.md Section 5).
    integer,intent(in)::n
    real(kind=8),intent(in)::p0,rho0
    real(kind=8),intent(in)::pi5(5,n,n,n),sxx(6,n,n,n),sxv(9,n,n,n)
    integer,intent(out)::nbad
    real(kind=8)::AA(3,3),VV(3,3),SS(3,3),Ai(3,3),Gam(3,3),TT(3,3),MM(3,3)
    real(kind=8)::det,g2,gmin
    integer::i,j,k,a,b
    gmin = huge(1.0d0); nbad = 0
    do k=1,n
       do j=1,n
          do i=1,n
             AA(1,1)=sxx(1,i,j,k); AA(2,2)=sxx(2,i,j,k); AA(3,3)=sxx(3,i,j,k)
             AA(1,2)=sxx(4,i,j,k); AA(2,1)=AA(1,2)
             AA(1,3)=sxx(5,i,j,k); AA(3,1)=AA(1,3)
             AA(2,3)=sxx(6,i,j,k); AA(3,2)=AA(2,3)
             do a=1,3
                do b=1,3
                   VV(a,b)=sxv(3*(a-1)+b,i,j,k)
                end do
             end do
             SS(1,1)=(p0+pi5(1,i,j,k))/rho0
             SS(2,2)=(p0+pi5(2,i,j,k))/rho0
             SS(3,3)=(p0-pi5(1,i,j,k)-pi5(2,i,j,k))/rho0
             SS(1,2)=pi5(3,i,j,k)/rho0; SS(2,1)=SS(1,2)
             SS(1,3)=pi5(4,i,j,k)/rho0; SS(3,1)=SS(1,3)
             SS(2,3)=pi5(5,i,j,k)/rho0; SS(3,2)=SS(2,3)
             call sym3_inv(AA,Ai,det)
             if(det<=0.0d0)then
                nbad = nbad+1
                gmin = 0.0d0
                cycle
             endif
             TT = matmul(Ai,VV)
             Gam = SS - matmul(transpose(VV),TT)
             ! Generalized eigenvalue lam_min(S^-1 Gam) via a congruence with
             ! S^-1/2; S is SPD whenever the pressure cone holds.
             call sym3_inv(SS,MM,det)
             if(det<=0.0d0)then
                nbad = nbad+1
                gmin = 0.0d0
                cycle
             endif
             g2 = sym3_lam_min(matmul(MM,Gam))
             if(g2<0.0d0) nbad = nbad+1
             gmin = min(gmin,sqrt(max(g2,0.0d0)))
          end do
       end do
    end do
    incomp_rank_min = gmin
  end function incomp_rank_min

  subroutine sym3_inv(A,Ai,det)
    real(kind=8),intent(in)::A(3,3)
    real(kind=8),intent(out)::Ai(3,3),det
    real(kind=8)::id
    det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
        - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
        + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
    Ai = 0.0d0
    if(det<=0.0d0) return
    id = 1.0d0/det
    Ai(1,1)=(A(2,2)*A(3,3)-A(2,3)*A(3,2))*id
    Ai(2,2)=(A(1,1)*A(3,3)-A(1,3)*A(3,1))*id
    Ai(3,3)=(A(1,1)*A(2,2)-A(1,2)*A(2,1))*id
    Ai(1,2)=(A(1,3)*A(3,2)-A(1,2)*A(3,3))*id; Ai(2,1)=Ai(1,2)
    Ai(1,3)=(A(1,2)*A(2,3)-A(1,3)*A(2,2))*id; Ai(3,1)=Ai(1,3)
    Ai(2,3)=(A(1,3)*A(2,1)-A(1,1)*A(2,3))*id; Ai(3,2)=Ai(2,3)
  end subroutine sym3_inv

  real(kind=8) function sym3_lam_min(Ain)
    ! Smallest eigenvalue of a 3x3 matrix that is symmetric to round-off,
    ! by the closed-form trigonometric solution of its characteristic cubic.
    real(kind=8),intent(in)::Ain(3,3)
    real(kind=8)::B(3,3),p1,q,p2,pp,r,phi,e1,e3,detB
    integer::ia,ib
    do ia=1,3
       do ib=1,3
          B(ia,ib) = 0.5d0*(Ain(ia,ib)+Ain(ib,ia))
       end do
    end do
    p1 = B(1,2)**2+B(1,3)**2+B(2,3)**2
    q  = (B(1,1)+B(2,2)+B(3,3))/3.0d0
    if(p1<=0.0d0)then
       sym3_lam_min = min(B(1,1),min(B(2,2),B(3,3)))
       return
    endif
    p2 = (B(1,1)-q)**2+(B(2,2)-q)**2+(B(3,3)-q)**2+2.0d0*p1
    pp = sqrt(p2/6.0d0)
    B(1,1)=B(1,1)-q; B(2,2)=B(2,2)-q; B(3,3)=B(3,3)-q
    B = B/pp
    detB = B(1,1)*(B(2,2)*B(3,3)-B(2,3)*B(3,2)) &
         - B(1,2)*(B(2,1)*B(3,3)-B(2,3)*B(3,1)) &
         + B(1,3)*(B(2,1)*B(3,2)-B(2,2)*B(3,1))
    r = max(-1.0d0,min(1.0d0,0.5d0*detB))
    phi = acos(r)/3.0d0
    e1 = q+2.0d0*pp*cos(phi)
    e3 = q+2.0d0*pp*cos(phi+2.0943951023931953d0)
    sym3_lam_min = min(e3,min(e1,3.0d0*q-e1-e3))
  end function sym3_lam_min

end module incomp_moments_module
