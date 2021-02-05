


subroutine pe_decomposition
  use main_module   
  implicit none
  integer :: n,tag=0,iloc(20),ierr,color,key
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status

! ----------------------------------
!      domain decomposition for each PE
! ----------------------------------
   if (n_pes>1) then
      if (n_pes_x*n_pes_y*n_pes_z*n_pes_m*n_pes_l*n_pes_k /= n_pes ) &
        call halt_stop(' n_pes_x x n_pes_y x n_pes_z x n_pes_k x n_pes_l x n_pes_m not equal number of PEs')

      !n_pes_k and nk always even or one
      if (nk>1 .and. mod(nk,2)/=0) call halt_stop(' nk must be even or one ')
      if (n_pes_k>1 .and. mod(n_pes_k,2)/=0) call halt_stop(' n_pes_k must be even or one')

      !n_pes_l and nl always even or one
      if (nl>1 .and. mod(nl,2)/=0) call halt_stop(' nl must be even or one ')
      if (n_pes_l>1 .and. mod(n_pes_l,2)/=0) call halt_stop(' n_pes_l must be even or one')
      
      !n_pes_m and nm always even or one
      if (nm>1.and.mod(nm,2)/=0) call halt_stop(' nm must be even or one')
      if (n_pes_m>1 .and. mod(n_pes_m,2)/=0) call halt_stop(' n_pes_m must be even or one')

      x_blk  = (nx-1)/n_pes_x + 1    ! extent of each block
      y_blk  = (ny-1)/n_pes_y + 1    ! extent of each block
      z_blk  = (nz-1)/n_pes_z + 1    ! extent of each block
      k_blk  = (nk-1)/n_pes_k + 1
      l_blk  = (nl-1)/n_pes_l + 1
      m_blk  = (nm-1)/n_pes_m + 1   
      
      my_blk_x = mod(my_pe,n_pes_x)+1                                     ! number of PE in x dir.
      my_blk_y = mod(my_pe/n_pes_x,n_pes_y)+1                             ! number of PE in y dir.
      my_blk_z = mod(my_pe/(n_pes_x*n_pes_y),n_pes_z) + 1                 ! number of PE in z dir.
      my_blk_k = mod(my_pe/(n_pes_x*n_pes_y*n_pes_z),n_pes_k) + 1         ! number of PE in k-dir.
      my_blk_l = mod(my_pe/(n_pes_x*n_pes_y*n_pes_z*n_pes_k),n_pes_l) + 1 ! number of PE in l-dir.
      my_blk_m =     my_pe/(n_pes_x*n_pes_y*n_pes_z*n_pes_k*n_pes_l)  + 1 ! number of PE in m-dir.
      
      xs_pe = (my_blk_x-1)*x_blk + 1 ! start index in y dir of this PE
      xe_pe = min(my_blk_x*x_blk,nx)
      ys_pe = (my_blk_y-1)*y_blk + 1 ! start index in y dir of this PE
      ye_pe = min(my_blk_y*y_blk,ny)
      zs_pe = (my_blk_z-1)*z_blk + 1 ! start index in z dir of this PE
      ze_pe = min(my_blk_z*z_blk,nz)
 
      if (n_pes_k>1) then
       if (my_blk_k >= n_pes_k/2+1) then
         ks_pe = nk/2+1+(my_blk_k-n_pes_k/2-1)*k_blk
         ke_pe = min(ks_pe+k_blk-1,nk)
       else
         ke_pe = nk/2-(n_pes_k/2-my_blk_k)*k_blk
         ks_pe = max(ke_pe-k_blk+1,1)
       endif
      else
       ks_pe = 1
       ke_pe = nk
      endif 
 
      if (n_pes_l>1) then
       if (my_blk_l >= n_pes_l/2+1) then
         ls_pe = nl/2+1+(my_blk_l-n_pes_l/2-1)*l_blk
         le_pe = min(ls_pe+l_blk-1,nl)
       else
         le_pe = nl/2-(n_pes_l/2-my_blk_l)*l_blk
         ls_pe = max(le_pe-l_blk+1,1)
       endif
      else
       ls_pe = 1
       le_pe = nl
      endif
     
      if (n_pes_m>1) then
       if (my_blk_m >= n_pes_m/2+1) then
         ms_pe = nm/2+1+(my_blk_m-n_pes_m/2-1)*m_blk
         me_pe = min(ms_pe+m_blk-1,nm)
       else
         me_pe = nm/2-(n_pes_m/2-my_blk_m)*m_blk
         ms_pe = max(me_pe-m_blk+1,1)
       endif
      else
       ms_pe = 1
       me_pe = nm
      endif
          
! ----------------------------------
!     last block might have been truncated
! ----------------------------------
      x_blk = xe_pe-xs_pe+1 
      y_blk = ye_pe-ys_pe+1 
      z_blk = ze_pe-zs_pe+1 
      k_blk = ke_pe-ks_pe+1 
      l_blk = le_pe-ls_pe+1 
      m_blk = me_pe-ms_pe+1      
 
! ---------------------------------- 
! communicator for all PEs of the same x,y,z block
! ---------------------------------- 
     color = my_blk_x + (my_blk_y + my_blk_z*n_pes_y)*n_pes_x
     key   = my_blk_k + (my_blk_l + my_blk_m*n_pes_l)*n_pes_k
     call mpi_comm_split(MPI_COMM_WORLD,color,key,MPI_comm_phys,ierr)
     
! ----------------------------------
!     check for incorrect domain decomposition
! ----------------------------------
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      if (my_blk_x==n_pes_x .and. xs_pe>xe_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in x-direction'
       print*,' choose other number of PEs in x-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_y==n_pes_y .and. ys_pe>ye_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in y-direction'
       print*,' choose other number of PEs in y-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_z==n_pes_z .and. zs_pe>ze_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in z-direction'
       print*,' choose other number of PEs in z-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_k==n_pes_k .and. ks_pe>ke_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in k-direction'
       print*,' choose other number of PEs in k-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_l==n_pes_l .and. ls_pe>le_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in l-direction'
       print*,' choose other number of PEs in l-direction'
       call halt_stop(' in pe_decomposition')
      endif
      if (my_blk_m==n_pes_m .and. ms_pe>me_pe) then
       print*,' ERROR:'
       print*,' domain decompositon impossible in m-direction'
       print*,' choose other number of PEs in m-direction'
       call halt_stop(' in pe_decomposition')
      endif
 
      
   else
       n_pes_x = 1;n_pes_y = 1;n_pes_z = 1; n_pes_k = 1; n_pes_l = 1; n_pes_m = 1
       x_blk = nx; y_blk = ny; z_blk = nz; k_blk = nk; l_blk = nl; m_blk = nm
       my_blk_x = 1;my_blk_y = 1;my_blk_z = 1 ; my_blk_k = 1 ; my_blk_l = 1; my_blk_m = 1
       xs_pe = 1; xe_pe = nx
       ys_pe = 1; ye_pe = ny
       zs_pe = 1; ze_pe = nz
       ks_pe = 1; ke_pe = nk
       ls_pe = 1; le_pe = nl
       ms_pe = 1; me_pe = nm 
   endif
! ----------------------------------
!      print out the PE decomposition, let PE 0 talk
! ----------------------------------
   if (my_pe==0) print'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
        ' domain size is nx=',nx,' x ny=',ny,' x nz=',nz,' x nk=',nk,' x nl=',nl,' x nm=',nm
   do n=0,n_pes-1
     if (n==0) then
       iloc(1:12) = (/xs_pe,xe_pe,ys_pe,ye_pe,zs_pe,ze_pe,ks_pe,ke_pe,ls_pe,le_pe,ms_pe,me_pe/)
     else
       if (my_pe==n) then
          iloc(1:12) = (/xs_pe,xe_pe,ys_pe,ye_pe,zs_pe,ze_pe,ks_pe,ke_pe,ls_pe,le_pe,ms_pe,me_pe/)
          call mpi_send(iloc,12,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
       endif
       if (my_pe==0) call mpi_recv(iloc,12,mpi_integer,n,tag,MPI_COMM_WORLD,status,ierr)
     endif
     if (my_pe==0) print'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
           ' domain of PE #',n,': x=',iloc(1),':',iloc(2),' y=',iloc(3),':',iloc(4),' z=',iloc(5),':',iloc(6), &
                               ' k=',iloc(7),':',iloc(8),' l=',iloc(9),':',iloc(10),' m=',iloc(11),':',iloc(12)
   enddo

   if (my_pe==0) print*,' '

   do n=0,n_pes-1
     if (n==0) then
       iloc(1:6) = (/my_blk_x,my_blk_y,my_blk_z,my_blk_k,my_blk_l,my_blk_m/)
     else
       if (my_pe==n) then
          iloc(1:6) = (/my_blk_x,my_blk_y,my_blk_z,my_blk_k,my_blk_l,my_blk_m/)
          call mpi_send(iloc,6,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
       endif
       if (my_pe==0) call mpi_recv(iloc,6,mpi_integer,n,tag,MPI_COMM_WORLD,status,ierr)
     endif
     if (my_pe==0) print'(a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a,i4)', &
          ' pe#',n,': my_blk_x=',iloc(1),' my_blk_y=',iloc(2),' my_blk_z=',iloc(3) ,&
                    ' my_blk_k=',iloc(4),' my_blk_l=',iloc(5),' my_blk_m=',iloc(6) 
   enddo

   if (my_pe==0) print*,' '
   call fortran_barrier
   
   
end subroutine pe_decomposition




 subroutine my_mpi_init
!--------------------------------------------------------------
!     intitialize mpi system for model
!--------------------------------------------------------------
  use main_module   
  implicit none
  integer :: nlen,ierr
  include "mpif.h"
  character (len=MPI_MAX_PROCESSOR_NAME) :: pname
      
  call MPI_Comm_rank(MPI_COMM_WORLD, my_pe, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, n_pes, ierr)
  call MPI_Get_processor_name(pname, nlen, ierr)
  call my_mpi_test
 end subroutine my_mpi_init


subroutine halt_stop(string)
!--------------------------------------------------------------
!     controlled stop, should not be called from python
!--------------------------------------------------------------
      implicit none
      character*(*) :: string
      integer :: ierr,code,my_pe
      include "mpif.h"
      call mpi_comm_rank(MPI_COMM_WORLD,my_pe,ierr)
      print*,' global pe #',my_pe,' : ',string
      print*,' global pe #',my_pe,' aborting '
      code=99
      call MPI_ABORT(mpi_comm_world, code, IERR)
end subroutine halt_stop



subroutine fortran_barrier
!--------------------------------------------------------------
!     A barrier for the local sub domain
!     for use in fortran part only
!--------------------------------------------------------------
  
  implicit none
  include "mpif.h"
  integer :: ierr
  call mpi_barrier(MPI_COMM_WORLD, ierr)
end subroutine fortran_barrier



subroutine my_mpi_test
!--------------------------------------------------------------
!     test some basic mpi routines
!--------------------------------------------------------------
      implicit none
      integer :: my_pe=-1,all_pes,xint,xint2,ierr
      real*8    :: xreal,xreal2
      include "mpif.h"
!   get some mpi infos first
      call mpi_comm_rank(MPI_COMM_WORLD  ,my_pe,ierr)
      if (my_pe==0) print*,' testing mpi routines'
      call mpi_comm_size(MPI_COMM_WORLD,all_pes,ierr)
!   try first global barrier
      call mpi_barrier(MPI_COMM_WORLD, ierr)
!   try broadcasting
      xreal = 1.0
      call mpi_bcast(xreal,1,mpi_real8,0,MPI_COMM_WORLD ,ierr)
      xint = 1
      call mpi_bcast(xint,1,mpi_integer,0,MPI_COMM_WORLD,ierr)
!   check results of broadcasting
      if (xreal /= 1.0 ) then
       print*,'fatal: MPI test failed on broadcasting reals for PE #',my_pe
       stop
      endif
      if (xint /= 1 ) then
       print*,'fatal: MPI test failed on broadcasting integer for PE #',my_pe
       stop
      endif
      call mpi_barrier(MPI_COMM_WORLD , ierr)
!   try global sum
      xreal = 2.0
      call mpi_allreduce(xreal,xreal2,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD ,ierr)
      xint = 2
      call mpi_allreduce(xint,xint2,1,mpi_integer,MPI_SUM,MPI_COMM_WORLD,ierr)
!   check results 
      xreal = xreal2/all_pes
      if (xreal /= 2.0 ) then
       print*,'fatal: MPI test failed on global sum (real) for PE #',my_pe
       stop
      endif
      xint = xint2/all_pes
      if (xint /= 2.0 ) then
       print*,'fatal: MPI test failed on global sum (int) for PE #',my_pe
       stop
      endif
      call mpi_barrier(MPI_COMM_WORLD , ierr)
end subroutine my_mpi_test


subroutine pe0_bcast_int(a,len)
!--------------------------------------------------------------
!     Broadcast an integer vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      integer, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_integer,0,MPI_COMM_WORLD,ierr)
end subroutine pe0_bcast_int


subroutine pe0_bcast(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      real*8, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_real8,0,MPI_COMM_WORLD,ierr)
end subroutine pe0_bcast



subroutine pe0_bcast_complex(a,len)
!--------------------------------------------------------------
!     Broadcast a vector from pe0 to all other pe
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer, intent(in) :: len
      complex*16, intent(inout) :: a(len)
      integer :: ierr
      include "mpif.h"
      call mpi_bcast(a,len,mpi_double_complex,0,MPI_COMM_WORLD,ierr)
end subroutine pe0_bcast_complex


subroutine bcast_real(x,len,pe)
!--------------------------------------------------------------
!     Broadcast a real vector from PE pe to others
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer :: len,ierr,pe
      real*8 :: x(len)
      include "mpif.h"
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(x,len,mpi_real8,pe,MPI_COMM_WORLD,ierr)
end subroutine bcast_real


subroutine bcast_integer(x,len,pe)
!--------------------------------------------------------------
!     Broadcast an integer vector from PE pe to others
!--------------------------------------------------------------
      use main_module
      implicit none
      integer :: len,ierr,pe
      integer :: x(len)
      include "mpif.h"
      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_bcast(x,len,mpi_integer,pe,MPI_COMM_WORLD,ierr)
end subroutine bcast_integer



subroutine global_max(x)
!--------------------------------------------------------------
!     Get the max of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MAX,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
 end subroutine global_max


subroutine global_min(x)
!--------------------------------------------------------------
!     Get the min of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_MIN,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
end subroutine global_min


subroutine global_sum(x)
!--------------------------------------------------------------
!     Do a sum of real x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8,intent(inout)    :: x
      real*8    :: x_sym,x_sym2
      integer :: ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_real8,MPI_SUM,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
end subroutine global_sum






subroutine global_max_int(x)
!--------------------------------------------------------------
!     Get the max of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MAX,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
 end subroutine global_max_int


subroutine global_min_int(x)
!--------------------------------------------------------------
!     Get the min of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_MIN,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
end subroutine global_min_int


subroutine global_sum_int(x)
!--------------------------------------------------------------
!     Do a sum of integer x over all PEs in sub domain
!--------------------------------------------------------------
      use main_module   
      implicit none
      integer,intent(inout)    :: x
      integer    :: x_sym,x_sym2,ierr
      include "mpif.h"
      x_sym = x
      call mpi_allreduce(x_sym,x_sym2,1,mpi_integer,MPI_SUM,MPI_COMM_WORLD       ,ierr)
      x = x_sym2
end subroutine global_sum_int


subroutine wavenumber_sum(x)
!--------------------------------------------------------------
! sum x over wavenumber dimensions
!--------------------------------------------------------------
 use main_module   
 implicit none
 real*8, intent(inout)  :: x(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8                 :: x2(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 integer :: len,ierr
 include "mpif.h"
 len = (x_blk+2*onx)*(y_blk+2*onx)*(z_blk+2*onx)
 call mpi_allreduce(x,x2,len,mpi_real8,MPI_SUM,MPI_COMM_phys,ierr)
 x=x2
end subroutine wavenumber_sum



subroutine border_exchg(a)
!--------------------------------------------------------------
! Exchange overlapping areas of array a in all PEs 
!--------------------------------------------------------------
  use main_module   
  implicit none
  real*8, intent(inout)  :: a(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx, &
                              ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx)
  integer               :: xs,xe,ys,ye,ms,me,zs,ze,ls,le,ks,ke                            
  integer  ::  tag=0, ierr,i,len,back,forth
  include "mpif.h"
  integer,dimension(MPI_STATUS_SIZE)  :: Status
     
  ks=ks_pe; ke=ke_pe; ls=ls_pe; le=le_pe; ms=ms_pe; me=me_pe
  xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe
   
  ! exchange across m borders
  if (n_pes_m > 1) then
     back   = my_pe - n_pes_x*n_pes_z*n_pes_y*n_pes_l*n_pes_k
     forth  = my_pe + n_pes_x*n_pes_z*n_pes_y*n_pes_l*n_pes_k
     len=(x_blk)*(y_blk)*(z_blk)*(k_blk)*(l_blk)
     do i=1,onx
       if (my_blk_m /=1 )       call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,ms+i-1),&
                                               len,mpi_real8,back ,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_m /= n_pes_m) call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,me+i  ) &
                                              ,len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_m /= n_pes_m) call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,me-i+1), &
                                               len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_m /=1 )       call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,ms-i  ), &
                                               len,mpi_real8,back ,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif

  ! exchange across l borders
   if (n_pes_l > 1) then
     back   = my_pe - n_pes_x*n_pes_z*n_pes_y*n_pes_k
     forth  = my_pe + n_pes_x*n_pes_z*n_pes_y*n_pes_k
     len=(x_blk)*(y_blk)*(z_blk)*(k_blk)*(m_blk)
     do i=1,onx
       if (my_blk_l /=1 )       call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,ls+i-1,ms:me),&
                                               len,mpi_real8,back ,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_l /= n_pes_l) call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,le+i  ,ms:me),&
                                               len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_l /= n_pes_l) call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,le-i+1,ms:me),&
                                               len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_l /=1 )       call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,ls-i  ,ms:me),&
                                               len,mpi_real8,back ,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif
  
  
  ! exchange across k borders
  if (n_pes_k > 1) then
     back   = my_pe - n_pes_x*n_pes_z*n_pes_y
     forth  = my_pe + n_pes_x*n_pes_z*n_pes_y
     len=(x_blk)*(y_blk)*(z_blk)*(l_blk)*(m_blk)
     do i=1,onx
       if (my_blk_k /=1 )       call mpi_send(a(xs:xe,ys:ye,zs:ze,ks+i-1,ls:le,ms:me),&
                                              len,mpi_real8,back ,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_k /= n_pes_k) call mpi_recv(a(xs:xe,ys:ye,zs:ze,ke+i  ,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_k /= n_pes_k) call mpi_send(a(xs:xe,ys:ye,zs:ze,ke-i+1,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_k /=1 )       call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks-i  ,ls:le,ms:me),&
                                              len,mpi_real8,back ,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif

  ! exchange across z borders
  if ( n_pes_z > 1) then 
     back  = my_pe - n_pes_x*n_pes_y 
     forth = my_pe + n_pes_x*n_pes_y 
     len=(x_blk)*(y_blk)*(k_blk)*(l_blk)*(m_blk)
     do i=1,onx
       if (my_blk_z /=1 )       call mpi_send(a(xs:xe,ys:ye,zs+i-1,ks:ke,ls:le,ms:me), &
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_z /= n_pes_z) call mpi_recv(a(xs:xe,ys:ye,ze+i  ,ks:ke,ls:le,ms:me), &
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_z /= n_pes_z) call mpi_send(a(xs:xe,ys:ye,ze-i+1,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_z /=1 )       call mpi_recv(a(xs:xe,ys:ye,zs-i  ,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif

  ! exchange across y borders
  if ( n_pes_y > 1) then
     back = my_pe - n_pes_x
     forth = my_pe + n_pes_x
     len=(x_blk)*(z_blk)*(k_blk)*(l_blk)*(m_blk)
     do i=1,onx
       if (my_blk_y /=1 )       call mpi_send(a(xs:xe,ys+i-1,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_y /= n_pes_y) call mpi_recv(a(xs:xe,ye+i  ,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_y /= n_pes_y) call mpi_send(a(xs:xe,ye-i+1,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_y /=1 )       call mpi_recv(a(xs:xe,ys-i  ,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif
 
  ! exchange across x borders
  if ( n_pes_x > 1) then
     back = my_pe-1
     forth = my_pe+1
     len=(y_blk)*(z_blk)*(k_blk)*(l_blk)*(m_blk)
     do i=1,onx
       if (my_blk_x /=1 )       call mpi_send(a(xs+i-1,ys:ye,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_x /= n_pes_x) call mpi_recv(a(xe+i  ,ys:ye,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,status,ierr)
     enddo
     do i=1,onx
       if (my_blk_x /= n_pes_x) call mpi_send(a(xe-i+1,ys:ye,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,forth,tag,MPI_COMM_WORLD,ierr)
       if (my_blk_x /=1 )       call mpi_recv(a(xs-i  ,ys:ye,zs:ze,ks:ke,ls:le,ms:me),&
                                              len,mpi_real8,back,tag,MPI_COMM_WORLD,status,ierr)
     enddo
  endif
 
end subroutine border_exchg



subroutine pe0_recv_arr_phys(a)
 !
 ! PEs with my_blk_k=my_blk_l=my_blk_m=1 send data of 3D array in phys space to PE 0
 !
 use main_module   
 implicit none
 real*8, intent(inout) :: a(nx,ny,nz)
 integer               :: xs,xe,ys,ye,iproc,zs,ze
 integer               :: tag=0, ierr,len,pe_phys,flag
 
 include "mpif.h"
 integer, dimension(MPI_STATUS_SIZE) :: Status
 xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe
 
 call MPI_Comm_rank(MPI_COMM_phys, pe_phys, ierr)
 
 do iproc=1,n_pes-1
     
  if ( my_pe == iproc ) then
   if (pe_phys==0) then 
      flag = 1 
      call mpi_send(flag,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(xs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(xe,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(ys,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(ye,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(zs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
      call mpi_send(ze,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr) 
      len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)
      call mpi_send(a(xs:xe,ys:ye,zs:ze),len,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else
      flag = 0
      call mpi_send(flag,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
   endif
  endif
  if ( my_pe == 0 ) then
    call mpi_recv(flag,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
    if (flag == 1) then 
        call mpi_recv(xs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(xe,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ys,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ye,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(zs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ze,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)
        call mpi_recv(a(xs:xe,ys:ye,zs:ze),len,mpi_real8,iproc,tag,MPI_COMM_WORLD,Status,ierr)
     endif   
  endif
      
 enddo     
end subroutine pe0_recv_arr_phys

 
 
 
 
 

subroutine pe0_recv_arr_mdim(a,i)
!--------------------------------------------------------------
!     all PEs send their data of a 5D array to PE0
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8, intent(inout) :: a(nx,ny,nz,nk,nl)
      integer               :: i,xs,xe,ys,ye,iproc,zs,ze,ls,le,ks,ke
      integer               :: tag=0, ierr,len,flag
      include "mpif.h"
      integer, dimension(MPI_STATUS_SIZE) :: Status
      ks=ks_pe; ke=ke_pe; ls=ls_pe; le=le_pe; 
      xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe
      
      do iproc=1,n_pes-1
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if ( my_pe == iproc) then
        if (i>=ms_pe.and.i<=me_pe ) then
         flag = 1
         call mpi_send(flag,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(xs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(xe,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ys,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ye,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(zs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ze,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ks,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ke,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(ls,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         call mpi_send(le,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
         len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)*(le-ls+1)*(ke-ks+1)
         call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le),len,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
        else
         flag = 0
         call mpi_send(flag,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        endif
       endif
       if ( my_pe == 0 ) then
        call mpi_recv(flag,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        if (flag==1) then
         call mpi_recv(xs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(xe,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ys,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ye,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(zs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ze,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ks,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ke,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(ls,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         call mpi_recv(le,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
         len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)*(le-ls+1)*(ke-ks+1)
         call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le),len,mpi_real8,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        endif
       endif
       call mpi_barrier(MPI_COMM_WORLD ,ierr)
      enddo     
end subroutine pe0_recv_arr_mdim



subroutine pe0_recv_arr(a)
!--------------------------------------------------------------
!     all PEs send their data of a 6D array to PE0
!--------------------------------------------------------------
      use main_module   
      implicit none
      real*8, intent(inout) :: a(nx,ny,nz,nk,nl,nm)
      integer               :: xs,xe,ys,ye,ms,me,iproc,zs,ze,ls,le,ks,ke
      integer               :: tag=0, ierr,len
      include "mpif.h"
      integer, dimension(MPI_STATUS_SIZE) :: Status
      ks=ks_pe; ke=ke_pe; ls=ls_pe; le=le_pe; ms=ms_pe; me=me_pe
      xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe
      
      do iproc=1,n_pes-1
       call mpi_barrier(MPI_COMM_WORLD,ierr)
       if ( my_pe == iproc ) then
        call mpi_send(xs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(xe,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ys,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ye,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(zs,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ze,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ks,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ke,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ls,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(le,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(ms,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        call mpi_send(me,1,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
        len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)*(me-ms+1)*(le-ls+1)*(ke-ks+1)
        call mpi_send(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,ms:me),len,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
       endif
       if ( my_pe == 0 ) then
        call mpi_recv(xs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(xe,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ys,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ye,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(zs,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ze,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ks,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ke,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ls,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(le,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(ms,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        call mpi_recv(me,1,mpi_integer,iproc,tag,MPI_COMM_WORLD,Status,ierr)
        len=(xe-xs+1)*(ye-ys+1)*(ze-zs+1)*(me-ms+1)*(le-ls+1)*(ke-ks+1)
        call mpi_recv(a(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,ms:me),len,mpi_real8,iproc,tag,MPI_COMM_WORLD,Status,ierr)
       endif
       call mpi_barrier(MPI_COMM_WORLD ,ierr)
      enddo     
end subroutine pe0_recv_arr





subroutine reflection_condition(a)
!--------------------------------------------------------------
! account for reflection boundary condition at bottom and surface
!--------------------------------------------------------------
 use main_module   
 implicit none
 real*8, intent(inout)  :: a(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                             ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx)
 integer  ::  tag=0, ierr,i,len,adress
 include "mpif.h"
 integer,dimension(MPI_STATUS_SIZE)  :: Status
 integer               :: xs,xe,ys,ye,ms,me,zs,ze,ls,le,ks,ke
    
  ks=ks_pe; ke=ke_pe; ls=ls_pe; le=le_pe; ms=ms_pe; me=me_pe
  xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe

if (enable_upper_reflection) then

  if (my_blk_z==n_pes_z.and.n_pes_m>1) then
  
   len = m_blk*(x_blk)*(y_blk)*(k_blk)*(l_blk)   
   ! buddy is my_blk_x, my_blk_y, my_blk_z, my_blk_k, my_blk_l, (n_pes_m-my_blk_m+1)
   adress = my_blk_x-1 + (my_blk_y-1)*n_pes_x + (my_blk_z-1)*n_pes_x*n_pes_y  &
          + (my_blk_k-1)*n_pes_x*n_pes_y*n_pes_z  + (my_blk_l-1)*n_pes_x*n_pes_y*n_pes_z*n_pes_k &
          +  (n_pes_m-my_blk_m)*n_pes_x*n_pes_y*n_pes_z*n_pes_k*n_pes_l 
   
   do i=1,onx
    ! from positive m to negative
    if (my_blk_m> n_pes_m/2)   call mpi_send(a(xs:xe,ys:ye,nz-i+1,ks:ke,ls:le,ms:me: 1), &
                                            len,mpi_real8,adress,tag,MPI_COMM_WORLD,ierr)   
    if (my_blk_m<=n_pes_m/2) call mpi_recv(a(xs:xe,ys:ye,nz+i  ,ks:ke,ls:le,me:ms:-1), &
                                            len,mpi_real8,adress,tag,MPI_COMM_WORLD,status,ierr) 
    ! from negative m to positive
    if (my_blk_m<=n_pes_m/2)  call mpi_send(a(xs:xe,ys:ye,nz-i+1,ks:ke,ls:le,ms:me: 1), &
                                            len,mpi_real8,adress,tag,MPI_COMM_WORLD,ierr)    
    if (my_blk_m> n_pes_m/2)  call mpi_recv(a(xs:xe,ys:ye,nz+i  ,ks:ke,ls:le,me:ms:-1), &
                                            len,mpi_real8,adress,tag,MPI_COMM_WORLD,status,ierr)    
   enddo 
     
  else if (my_blk_z==n_pes_z.and.n_pes_m==1) then  
   do i=1,onx
     a(xs:xe,ys:ye,nz+i,ks:ke,ls:le,ms_pe:me_pe) = a(xs:xe,ys:ye,nz-i+1,ks:ke,ls:le,me_pe:ms_pe:-1)
   enddo
  endif

endif

if (enable_lower_reflection) then

  if (my_blk_z==1.and.n_pes_m>1) then
  
   len = m_blk*(x_blk)*(y_blk)*(k_blk)*(l_blk)  
   adress = my_blk_x-1 + (my_blk_y-1)*n_pes_x + (my_blk_z-1)*n_pes_x*n_pes_y  &
          + (my_blk_k-1)*n_pes_x*n_pes_y*n_pes_z  + (my_blk_l-1)*n_pes_x*n_pes_y*n_pes_z*n_pes_k &
          +  (n_pes_m-my_blk_m)*n_pes_x*n_pes_y*n_pes_z*n_pes_k*n_pes_l 
   
   do i=1,onx
    ! from positive m to negative
    if (my_blk_m> n_pes_m/2) call mpi_send(a(xs:xe,ys:ye,1+i-1,ks:ke,ls:le,ms_pe:me_pe: 1),&
                                           len,mpi_real8,adress,tag,MPI_COMM_WORLD,ierr)
    if (my_blk_m<=n_pes_m/2) call mpi_recv(a(xs:xe,ys:ye,1-i  ,ks:ke,ls:le,me_pe:ms_pe:-1),&
                                           len,mpi_real8,adress,tag,MPI_COMM_WORLD,status,ierr) 
    ! from negative m to positive
    if (my_blk_m<=n_pes_m/2) call mpi_send(a(xs:xe,ys:ye,1+i-1,ks:ke,ls:le,ms_pe:me_pe: 1),&
                                           len,mpi_real8,adress,tag,MPI_COMM_WORLD,ierr)
    if (my_blk_m> n_pes_m/2) call mpi_recv(a(xs:xe,ys:ye,1-i  ,ks:ke,ls:le,me_pe:ms_pe:-1),&
                                           len,mpi_real8,adress,tag,MPI_COMM_WORLD,status,ierr)
   enddo
        
  else if (my_blk_z==1.and.n_pes_m==1) then  
   do i=1,onx
     a(xs:xe,ys:ye,1-i,ks:ke,ls:le,ms_pe:me_pe) = a(xs:xe,ys:ye,1+i-1,ks:ke,ls:le,me_pe:ms_pe:-1)
   enddo
  endif
  
endif


end subroutine reflection_condition



