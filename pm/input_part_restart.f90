module input_part_restart_module

  type :: out_input_star_t
     real(kind=8)::mass
  end type out_input_star_t

contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_input_part_restart(pst)
  use mdl_module
  use amr_parameters, only: ndim
  use ramses_commons, only: pst_t
  use output_amr_module, only: input_header
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE, TREE_TYPE
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to read and dispatch particles
  ! from a Ramses restart file.
  !--------------------------------------------------------------------
  integer::icpu,ilun,ncpu_file
  integer(kind=8)::npart_tot_file,npart_tot_check
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::file_head,file_part
  integer,allocatable,dimension(:)::npart_file
  type(out_input_star_t)::output

  if(pst%s%r%verbose)write(*,*)'Entering input_part_restart'

  if(pst%s%r%part)then

     ! Read particle files header
     call title(pst%s%r%nrestart,nchar)
     file_head='backup_'//TRIM(nchar)//'/part_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Restart snapshot has ",I12," DM particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=PART_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part='backup_'//TRIM(nchar)//'/part.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_restart(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in dark matter=',output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%star)then

     ! Read star particle files header
     call title(pst%s%r%nrestart,nchar)
     file_head='backup_'//TRIM(nchar)//'/star_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Restart snapshot has ",I12," star particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=STAR_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part='backup_'//TRIM(nchar)//'/star.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_restart(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in stars=',output%mass
     pst%s%g%mass_star_tot=output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%sink)then

     ! Read star particle files header
     call title(pst%s%r%nrestart,nchar)
     file_head='backup_'//TRIM(nchar)//'/sink_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Restart snapshot has ",I12," sink particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=SINK_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part='backup_'//TRIM(nchar)//'/sink.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_restart(pst,npart_file,ncpu_file+1,output,2)
     write(*,*)'Total mass in sinks=',output%mass
     pst%s%g%mass_sink_tot=output%mass

     ! Deallocate local array
     deallocate(npart_file)

  endif

  if(pst%s%r%tree)then

     ! Read star particle files header
     call title(pst%s%r%nrestart,nchar)
     file_head='backup_'//TRIM(nchar)//'/tree_header.txt'
     call input_header(pst%s%r,pst%s%g,file_head,npart_tot_file,ncpu_file)
     write(*,'(" Restart snapshot has ",I12," tree particles")')npart_tot_file

     ! Allocate local array
     allocate(npart_file(0:ncpu_file))

     ! Read number of particles in each file
     npart_tot_check=0
     npart_file(0)=TREE_TYPE
     do icpu=1,ncpu_file
        call title(icpu,ncharcpu)
        file_part='backup_'//TRIM(nchar)//'/tree.'//TRIM(ncharcpu)
        ilun=10
        open(unit=ilun,file=TRIM(file_part),access="stream",action="read",form='unformatted')
        read(ilun,POS=5)npart_file(icpu)
        npart_tot_check=npart_tot_check+npart_file(icpu)
        close(ilun)
     end do
     if(npart_tot_check.NE.npart_tot_file)then
        write(*,*)' Input file corrupted'
        call mdl_abort(pst%s%mdl)
     endif

     ! Call recursive slave routine
     call r_input_part_restart(pst,npart_file,ncpu_file+1,output,2)

     ! Deallocate local array
     deallocate(npart_file)

  endif

end subroutine m_input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_input_part_restart(pst,input_array,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE, TREE_TYPE
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  type(out_input_star_t)::output, next_output

  integer::rID
  integer::p_type

  !--------------------------------------------------------------------
  ! This routine is the recursive slave procedure to read and dispatch
  ! particles from a Ramses restart file.
  !--------------------------------------------------------------------

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_PART_RESTART,pst%iUpper+1,input_size,output_size,input_array)
     call r_input_part_restart(pst%pLower,input_array,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
     p_type=input_array(1)
     if(p_type==PART_TYPE)then
        call input_part_restart(pst%s%r,pst%s%g,pst%s%p,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==STAR_TYPE)then
        call input_part_restart(pst%s%r,pst%s%g,pst%s%star,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==SINK_TYPE)then
        call input_part_restart(pst%s%r,pst%s%g,pst%s%sink,input_size-1,input_array(2:input_size),output%mass)
     endif
     if(p_type==TREE_TYPE)then
        call input_part_restart(pst%s%r,pst%s%g,pst%s%tree,input_size-1,input_array(2:input_size),output%mass)
     endif
  endif

end subroutine r_input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_part_restart(r,g,p,ncpu_file,npart_file,mpart_loc)
  use amr_parameters, only: ndim, dp, i8b
  use amr_commons, only: run_t, global_t
  use pm_parameters, only: PART_TYPE, STAR_TYPE, SINK_TYPE, TREE_TYPE
  use pm_commons, only: part_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(part_t)::p
  real(kind=8)::mpart_loc
  integer::ncpu_file
  integer,dimension(1:ncpu_file)::npart_file
  !------------------------------------------------------------
  ! Read particles positions and velocities from a Ramses
  ! restart file and allocate particle-based arrays.
  !------------------------------------------------------------
  integer::ipart,ipart_old
  integer::i,idim,icpu,ileft,iright,nrest
  integer::nbyte=storage_size(1.0_dp)/8
  integer::istart,iend
  integer(kind=8)::nleft,nright,npart_tot,ipos,iskip
  integer(kind=8),dimension(0:ncpu_file)::ncum_file

  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8

  character(LEN=5)::nchar,ncharcpu
  character(LEN=10)::prefix
  character(LEN=80)::file_part

  !-------------------------------------
  ! Compute local particle number
  !-------------------------------------
  ncum_file(0)=0
  do icpu=1,ncpu_file
     ncum_file(icpu)=ncum_file(icpu-1)+npart_file(icpu)
  end do
  npart_tot=ncum_file(ncpu_file)
  p%npart_tot=npart_tot

  p%npart=npart_tot/g%ncpu
  nleft=(g%myid-1)*p%npart
  nright=g%myid*p%npart
  nrest=npart_tot-p%npart*g%ncpu
  if(g%myid.LE.nrest)then
     p%npart=p%npart+1
     nleft=nleft+(g%myid-1)
     nright=nright+g%myid
  else
     nleft=nleft+nrest
     nright=nright+nrest
  endif

  ! Compute interval of file to open for current process
  ileft=0
  iright=-1
  if(nright.GT.nleft)then
     do icpu=1,ncpu_file
        if(icpu>1)then
           if(ncum_file(icpu).GT.nleft.AND.ncum_file(icpu-1).LT.nright)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        else
           if(ncum_file(icpu).GT.nleft)then
              if(ileft==0)ileft=icpu
              iright=MAX(icpu,iright)
           endif
        endif
     end do
  endif

  ! Determine prefix based on particle type
  if(p%type==PART_TYPE)prefix='part'
  if(p%type==STAR_TYPE)prefix='star'
  if(p%type==SINK_TYPE)prefix='sink'
  if(p%type==TREE_TYPE)prefix='tree'

  ! Loop over relevant files
  ipart=0
  ipart_old=0
  mpart_loc=0.0d0
  call title(r%nrestart,nchar)
  do icpu=ileft,iright
     if(icpu>1)then
        istart=MAX(nleft-ncum_file(icpu-1),0)+1
        iend=MIN(nright-ncum_file(icpu-1),npart_file(icpu))
     else
        istart=nleft+1
        iend=MIN(nright,npart_file(icpu))
     endif

     ! Open the particle file
     call title(icpu,ncharcpu)
     file_part='backup_'//TRIM(nchar)//'/'//TRIM(prefix)//'.'//TRIM(ncharcpu)
     open(unit=10,file=TRIM(file_part),access="stream",action="read",form='unformatted')

     iskip=9

     allocate(xdp(istart:iend))

     ! Read positions
     do idim=1,ndim
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%xp(ipart,idim)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     end do

     ! Read velocities
     do idim=1,ndim
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%vp(ipart,idim)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     end do

     ! Read masses
     ipos=iskip+nbyte*int(istart-1,8)
     read(10,POS=ipos)xdp
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%mp(ipart)=xdp(i)
        mpart_loc=mpart_loc+xdp(i)
     end do
     iskip=iskip+nbyte*int(npart_file(icpu),8)

     ! Read metallicity
     if(allocated(p%zp))then
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%zp(ipart)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     endif

     ! Read accelerations
     if(allocated(p%fp))then
        do idim=1,ndim
           ipos=iskip+nbyte*int(istart-1,8)
           read(10,POS=ipos)xdp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%fp(ipart,idim)=xdp(i)
           end do
           iskip=iskip+nbyte*int(npart_file(icpu),8)
        end do
     endif

     ! Read Angular Momenta
     if(allocated(p%jp))then
        do idim=1,ndim
           ipos=iskip+nbyte*int(istart-1,8)
           read(10,POS=ipos)xdp
           ipart=ipart_old
           do i=istart,iend
              ipart=ipart+1
              p%jp(ipart,idim)=xdp(i)
           end do
           iskip=iskip+nbyte*int(npart_file(icpu),8)
        end do
     endif

     ! Read birth time
     if(allocated(p%tp))then
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%tp(ipart)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     endif

     ! Read merging time
     if(allocated(p%tm))then
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%tm(ipart)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     endif

     ! Read size
     if(allocated(p%size))then
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%size(ipart)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     endif

     ! Read charge
     if(allocated(p%charge))then
        ipos=iskip+nbyte*int(istart-1,8)
        read(10,POS=ipos)xdp
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%charge(ipart)=xdp(i)
        end do
        iskip=iskip+nbyte*int(npart_file(icpu),8)
     endif

     deallocate(xdp)

     allocate(isp(istart:iend))

     ! Read level
     ipos=iskip+4*int(istart-1,8)
     read(10,POS=ipos)isp
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%levelp(ipart)=isp(i)
     end do
     iskip=iskip+4*int(npart_file(icpu),8)

     deallocate(isp)

     allocate(isp8(istart:iend))

     ! Read identity
#ifndef LONGINT
     ipos=iskip+4*int(istart-1,8)
#else
     ipos=iskip+8*int(istart-1,8)
#endif
     read(10,POS=ipos)isp8
     ipart=ipart_old
     do i=istart,iend
        ipart=ipart+1
        p%idp(ipart)=isp8(i)
     end do
#ifndef LONGINT
     iskip=iskip+4*int(npart_file(icpu),8)
#else
     iskip=iskip+8*int(npart_file(icpu),8)
#endif

     ! Read merging identity
     if(allocated(p%idm))then
#ifndef LONGINT
        ipos=iskip+4*int(istart-1,8)
#else
        ipos=iskip+8*int(istart-1,8)
#endif
        read(10,POS=ipos)isp8
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%idm(ipart)=isp8(i)
        end do
#ifndef LONGINT
        iskip=iskip+4*int(npart_file(icpu),8)
#else
        iskip=iskip+8*int(npart_file(icpu),8)
#endif
     endif

     ! Read tracking identity
     if(allocated(p%idt))then
#ifndef LONGINT
        ipos=iskip+4*int(istart-1,8)
#else
        ipos=iskip+8*int(istart-1,8)
#endif
        read(10,POS=ipos)isp8
        ipart=ipart_old
        do i=istart,iend
           ipart=ipart+1
           p%idt(ipart)=isp8(i)
        end do
#ifndef LONGINT
        iskip=iskip+4*int(npart_file(icpu),8)
#else
        iskip=iskip+8*int(npart_file(icpu),8)
#endif
     endif

     deallocate(isp8)

     ! Close the particle file
     close(10)
     ipart_old=ipart

  end do
  ! End loop over files

  ! Put all particles in levelmin
  p%headp=p%npart+1
  p%tailp=p%npart
  p%headp(r%levelmin)=1
  p%tailp(r%levelmin)=p%npart
  if(ANY(.not.r%periodic(1:ndim)))then
     p%headp(r%levelmin-1)=1
     p%tailp(r%levelmin-1)=0
  endif

end subroutine input_part_restart
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module input_part_restart_module
