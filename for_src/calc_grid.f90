


subroutine calc_grid
!---------------------------------------------------------------------------------
!    setup grid in x,y,z,k,l,m directions
!---------------------------------------------------------------------------------
 use main_module   
 implicit none
 integer :: i,n
 real*8, dimension(1-onx:nx+onx) :: dxt_gl,dxu_gl,xt_gl,xu_gl
 real*8, dimension(1-onx:ny+onx) :: dyt_gl,dyu_gl,yt_gl,yu_gl
 real*8, dimension(1-onx:nz+onx) :: dzt_gl,dzu_gl,zt_gl,zu_gl
 real*8, dimension(1-onx:nk+onx) :: dkt_gl,dku_gl,kt_gl,ku_gl
 real*8, dimension(1-onx:nl+onx) :: dlt_gl,dlu_gl,lt_gl,lu_gl
 real*8, dimension(1-onx:nm+onx) :: dmt_gl,dmu_gl,mt_gl,mu_gl
 include "mpif.h"
 integer :: tag = 99, ierr, sender, ind(3)
 integer,dimension(MPI_STATUS_SIZE)  :: Status
      
!--------------------------------------------------------------
! transfer from locally defined grid variables to global ones
!--------------------------------------------------------------
 
  dxt_gl(xs_pe:xe_pe) = dxt(xs_pe:xe_pe)
  do n=1,n_pes_x-1
   if (my_blk_x==n+1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/xs_pe,xe_pe,x_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dxt(xs_pe:xe_pe),x_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dxt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_x>1) call pe0_bcast(dxt_gl,nx+2*onx)
  
  dyt_gl(ys_pe:ye_pe) = dyt(ys_pe:ye_pe)
  do n=1,n_pes_y-1
   if (my_blk_x==1.and.my_blk_y==n+1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/ys_pe,ye_pe,y_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dyt(ys_pe:ye_pe),y_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dyt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_y>1) call pe0_bcast(dyt_gl,ny+2*onx)
  
  dzt_gl(zs_pe:ze_pe) = dzt(zs_pe:ze_pe)
  do n=1,n_pes_z-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==n+1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/zs_pe,ze_pe,z_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dzt(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dzt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_z>1) call pe0_bcast(dzt_gl,nz+2*onx)

  dkt_gl(ks_pe:ke_pe) = dkt(ks_pe:ke_pe)
  do n=1,n_pes_k-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==n+1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/ks_pe,ke_pe,k_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dkt(ks_pe:ke_pe),k_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dkt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_k>1) call pe0_bcast(dkt_gl,nk+2*onx)

  dlt_gl(ls_pe:le_pe) = dlt(ls_pe:le_pe)
  do n=1,n_pes_l-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==n+1.and.my_blk_m==1) then
    call mpi_send((/ls_pe,le_pe,l_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dlt(ls_pe:le_pe),l_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z*n_pes_k
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dlt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_l>1) call pe0_bcast(dlt_gl,nl+2*onx)

  dmt_gl(ms_pe:me_pe) = dmt(ms_pe:me_pe)
  do n=1,n_pes_m-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==n+1) then
    call mpi_send((/ms_pe,me_pe,m_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dmt(ms_pe:me_pe),m_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z*n_pes_k*n_pes_l
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dmt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
  enddo
  if (n_pes_m>1) call pe0_bcast(dmt_gl,nm+2*onx)

  ! fill boundary values
  do i=1,onx
    dxt_gl(nx+i) = dxt_gl(nx); dxt_gl(1-i) = dxt_gl(1) 
    dyt_gl(ny+i) = dyt_gl(ny); dyt_gl(1-i) = dyt_gl(1) 
    dzt_gl(nz+i) = dzt_gl(nz); dzt_gl(1-i) = dzt_gl(1) 
    dkt_gl(nk+i) = dkt_gl(nk); dkt_gl(1-i) = dkt_gl(1) 
    dlt_gl(nl+i) = dlt_gl(nl); dlt_gl(1-i) = dlt_gl(1) 
    dmt_gl(nm+i) = dmt_gl(nm); dmt_gl(1-i) = dmt_gl(1) 
  enddo

  ! construct x grid 
  call u_centered_grid(dxt_gl,dxu_gl,xt_gl,xu_gl,nx+2*onx)
  xt_gl = xt_gl - dxu_gl(1); xu_gl = xu_gl - dxu_gl(1)
  ! transfer to locally defined variables
  xt(xs_pe-onx:xe_pe+onx)  = xt_gl(xs_pe-onx:xe_pe+onx)
  xu(xs_pe-onx:xe_pe+onx)  = xu_gl(xs_pe-onx:xe_pe+onx)
  dxu(xs_pe-onx:xe_pe+onx) = dxu_gl(xs_pe-onx:xe_pe+onx)
  dxt(xs_pe-onx:xe_pe+onx) = dxt_gl(xs_pe-onx:xe_pe+onx)

  ! construct y grid 
  call u_centered_grid(dyt_gl,dyu_gl,yt_gl,yu_gl,ny+2*onx)
  yt_gl = yt_gl - dyu_gl(1); yu_gl = yu_gl - dyu_gl(1)
  ! transfer to locally defined variables
  yt(ys_pe-onx:ye_pe+onx)  = yt_gl(ys_pe-onx:ye_pe+onx)
  yu(ys_pe-onx:ye_pe+onx)  = yu_gl(ys_pe-onx:ye_pe+onx)
  dyu(ys_pe-onx:ye_pe+onx) = dyu_gl(ys_pe-onx:ye_pe+onx)
  dyt(ys_pe-onx:ye_pe+onx) = dyt_gl(ys_pe-onx:ye_pe+onx)

  ! construct z grid 
  call u_centered_grid(dzt_gl,dzu_gl,zt_gl,zu_gl,nz+2*onx)
  zt_gl = zt_gl - zu_gl(nz); zu_gl = zu_gl - zu_gl(nz) 
  ! transfer to locally defined variables
  zt(zs_pe-onx:ze_pe+onx)  = zt_gl(zs_pe-onx:ze_pe+onx)
  zu(zs_pe-onx:ze_pe+onx)  = zu_gl(zs_pe-onx:ze_pe+onx)
  dzu(zs_pe-onx:ze_pe+onx) = dzu_gl(zs_pe-onx:ze_pe+onx)
  dzt(zs_pe-onx:ze_pe+onx) = dzt_gl(zs_pe-onx:ze_pe+onx)

  ! construct k grid 
  call u_centered_grid(dkt_gl,dku_gl,kt_gl,ku_gl,nk+2*onx)
  if (nk>1) then 
    kt_gl = kt_gl - ku_gl(nk/2); ku_gl = ku_gl - ku_gl(nk/2)
  else
    kt_gl = k_fixed; ku_gl = k_fixed
  endif
  ! transfer to locally defined variables
  kt(ks_pe-onx:ke_pe+onx)  = kt_gl(ks_pe-onx:ke_pe+onx)
  ku(ks_pe-onx:ke_pe+onx)  = ku_gl(ks_pe-onx:ke_pe+onx)
  dku(ks_pe-onx:ke_pe+onx) = dku_gl(ks_pe-onx:ke_pe+onx)
  dkt(ks_pe-onx:ke_pe+onx) = dkt_gl(ks_pe-onx:ke_pe+onx)

  ! construct l grid 
  call u_centered_grid(dlt_gl,dlu_gl,lt_gl,lu_gl,nl+2*onx)
  if (nl>1) then 
    lt_gl = lt_gl - lu_gl(nl/2); lu_gl = lu_gl - lu_gl(nl/2)
  else
    lt_gl = l_fixed; lu_gl = l_fixed
  endif
  ! transfer to locally defined variables
  lt(ls_pe-onx:le_pe+onx)  = lt_gl(ls_pe-onx:le_pe+onx)
  lu(ls_pe-onx:le_pe+onx)  = lu_gl(ls_pe-onx:le_pe+onx)
  dlu(ls_pe-onx:le_pe+onx) = dlu_gl(ls_pe-onx:le_pe+onx)
  dlt(ls_pe-onx:le_pe+onx) = dlt_gl(ls_pe-onx:le_pe+onx)
  
  ! construct m grid 
  call u_centered_grid(dmt_gl,dmu_gl,mt_gl,mu_gl,nm+2*onx)
  if (nm>1) then 
   mt_gl = mt_gl - mu_gl(nm/2); mu_gl = mu_gl - mu_gl(nm/2)
  else
   mt_gl = m_fixed; mu_gl = m_fixed
  endif
  ! transfer to locally defined variables
  mt(ms_pe-onx:me_pe+onx)  = mt_gl(ms_pe-onx:me_pe+onx)
  mu(ms_pe-onx:me_pe+onx)  = mu_gl(ms_pe-onx:me_pe+onx)
  dmu(ms_pe-onx:me_pe+onx) = dmu_gl(ms_pe-onx:me_pe+onx)
  dmt(ms_pe-onx:me_pe+onx) = dmt_gl(ms_pe-onx:me_pe+onx)
  

  ! print some infos
  if (my_pe==0) then
   print*,''
   print'(a)','Grid setup: '
   print*,''
   print'(a,i4,a)',' with ',nx,' grid points in x-direction:'
   print'(a,f8.1,a,f8.1,a)',' from xt(1) = ',xt_gl(1) ,' m to xt(nx) = ',xt_gl(nx),' m'
   print'(a,f8.1,a,f8.1,a)',' with dxt(1)= ',dxt_gl(1),' m to dxt(nx)= ',dxt_gl(nx),' m'
   print*,''
   print'(a,i4,a)',' with ',ny,' grid points in y-direction:'
   print'(a,f8.1,a,f8.1,a)',' from yt(1) = ',yt_gl(1) ,' m to yt(ny) = ',yt_gl(ny),' m'
   print'(a,f8.1,a,f8.1,a)',' with dyt(1)= ',dyt_gl(1),' m to dyt(ny)= ',dyt_gl(ny),' m'
   print*,''
   print'(a,i4,a)',' with ',nz,' grid points in z-direction:'
   print'(a,f8.2,a,f8.2,a)',' from zt(1) = ',zt_gl(1) ,' m to zt(nz) = ',zt_gl(nz),' m'
   print'(a,f8.2,a,f8.2,a)',' with dzt(1)= ',dzt_gl(1),' m to dzt(nz)= ',dzt_gl(nz),' m'
   print*,''
   print'(a,f8.2,a,f8.2,a)',' from zu(1) = ',zu_gl(1) ,' m to zu(nz) = ',zu_gl(nz),' m'
   print'(a,f8.2,a,f8.2,a)',' with dzu(1)= ',dzu_gl(1),' m to dzu(nz)= ',dzu_gl(nz),' m'
   print*,''
   print'(a,i4,a)',' with ',nk,' grid points in k-direction:'
   print'(a,e12.4,a,e12.4,a)',' from kt(1) = ',kt_gl(1) ,' 1/m to kt(nk) = ',kt_gl(nk),' 1/m'
   print'(a,e12.4,a,e12.4,a)',' with dkt(1)= ',dkt_gl(1),' 1/m to dkt(nk)= ',dkt_gl(nk),' 1/m'
   print*,''
   print'(a,i4,a)',' with ',nl,' grid points in l-direction:'
   print'(a,e12.4,a,e12.4,a)',' from lt(1) = ',lt_gl(1), ' 1/m to lt(nl) = ',lt_gl(nl),' 1/m'
   print'(a,e12.4,a,e12.4,a)',' with dlt(1)= ',dlt_gl(1),' 1/m to dlt(nl)= ',dlt_gl(nl),' 1/m'
   print*,''
   print'(a,i4,a)',' with ',nm,' grid points in m-direction:'
   print'(a,e12.4,a,e12.4,a)',' from mt(1) = ',mt_gl(1) ,' 1/m to mt(nm) = ',mt_gl(nm),' 1/m'
   print'(a,e12.4,a,e12.4,a)',' with dmt(1)= ',dmt_gl(1),' 1/m to dmt(nm)= ',dmt_gl(nm),' 1/m'
   print*,'' 
  endif 
  
  if (my_pe==0 .and. enable_show_grid_details) then
   print*, ''
   print*,' x grid details: '
   do i=1,nx
     print'(a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2 )',&
                        ' xt(',i,')=',xt_gl(i),' xu(',i,')=',xu_gl(i), &
                        ' dxt(',i,')=',dxt_gl(i),' dxu(',i,')=',dxu_gl(i)
   enddo
   print*, ''
   print*,' y grid details: '
   do i=1,ny
     print'(a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2 )',&
                        ' yt(',i,')=',yt_gl(i),' yu(',i,')=',yu_gl(i), &
                        ' dyt(',i,')=',dyt_gl(i),' dyu(',i,')=',dyu_gl(i)
   enddo  
   print*, ''
   print*,' z grid details: '
   do i=1,nz
     print'(a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2,a,i3,a,f8.2 )',&
                        ' zt(',i,')=',zt_gl(i),' zu(',i,')=',zu_gl(i), &
                        ' dzt(',i,')=',dzt_gl(i),' dzu(',i,')=',dzu_gl(i)
   enddo

   print*,'' 
   print*,' k grid details: '
   do i=1,nk
     print'(a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4 )',&
                        ' kt(',i,')=',kt_gl(i),' ku(',i,')=',ku_gl(i), &
                        ' dkt(',i,')=',dkt_gl(i),' dku(',i,')=',dku_gl(i)
   enddo
   print*, ''
   print*,' l grid details: '
   do i=1,nl
     print'(a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4 )',&
                        ' lt(',i,')=',lt_gl(i),' lu(',i,')=',lu_gl(i), &
                        ' dlt(',i,')=',dlt_gl(i),' dlu(',i,')=',dlu_gl(i)
   enddo   
   print*, ''
   print*,' m grid details: '
   do i=1,nm
     print'(a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4,a,i3,a,e12.4 )',&
                        ' mt(',i,')=',mt_gl(i),' mu(',i,')=',mu_gl(i), &
                        ' dmt(',i,')=',dmt_gl(i),' dmu(',i,')=',dmu_gl(i)
   enddo
  endif 
  
  ! check the grid
  do i=2,nz
    if (dzu_gl(i)-dzu_gl(i-1)> 1d-12 .or. dzt_gl(i)-dzt_gl(i-1) > 1d-12) then
     if (my_pe==0) then
      print*,' Error in vertical grid setup, dzu/dzt do not decrease monotonically '
      print*,' problem is at i=',i,dzt_gl(i)-dzt_gl(i-1),dzu_gl(i)-dzu_gl(i-1)
     endif 
     call fortran_barrier
     call halt_stop(' in calc_grid ')
    endif
  enddo
  
  
  do i=2,nm/2
    if (dmu_gl(i)-dmu_gl(i-1)> 1d-12  .or. dmt_gl(i)-dmt_gl(i-1) > 1d-12 ) then
     if (my_pe==0) then
      print*,' Error in m grid setup, dmu/dmt do not decrease monotonically '
      print*,' problem is at i=',i,dmt_gl(i)-dmt_gl(i-1),dmu_gl(i)-dmu_gl(i-1)
     endif 
     call fortran_barrier
     call halt_stop(' in calc_grid ')
    endif
  enddo 
  
  do i=2,nl/2
    if (dlu_gl(i)-dlu_gl(i-1)> 1d-12   .or. dlt_gl(i)-dlt_gl(i-1) > 1d-12 ) then
     if (my_pe==0) then
      print*,' Error in m grid setup, dlu/dlt do not decrease monotonically '
      print*,' problem is at i=',i,dlt_gl(i)-dlt_gl(i-1),dlu_gl(i)-dlu_gl(i-1)
     endif 
     call fortran_barrier
     call halt_stop(' in calc_grid ')
    endif
  enddo  

end subroutine calc_grid




subroutine u_centered_grid(dyt,dyu,yt,yu,n)
!---------------------------------------------------------------------------------
! setup u-centered grid based in Delta yt and the relations
! dyt_i = yu_i - yu_i-1 , yu_i = 0.5(yt_i+yt_(i+1)) , dyu_i = yt_(i+1)-yt_i
!---------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: dyt(n)
  real*8, intent(out) :: yu(n),yt(n),dyu(n)
  integer :: i 
  yu(1)=0
  do i=2,n  
   yu(i)=yu(i-1)+dyt(i)
  enddo
  yt(1)=yu(1)-dyt(1)*0.5
  do i=2,n
   yt(i) = 2*yu(i-1) - yt(i-1)
  enddo
  do i=1,n-1 
   dyu(i)= yt(i+1)-yt(i)
  enddo
  dyu(n)=2*dyt(n)- dyu(max(1,n-1))
end subroutine u_centered_grid

