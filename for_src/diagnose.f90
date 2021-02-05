

subroutine init_diag
 use main_module
 implicit none
 if (my_pe==0) print*,'initialize snapshot file'
 call init_snap3D_cdf
 
 if (enable_write_6D_single) then 
   call init_snap6D_cdf_single
 else
   call init_snap6D_cdf_many
 endif  
 
end subroutine init_diag




subroutine diagnose
 use main_module
 use timing_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,itdimid,ilen,itimeid,id,i,j,k,l,n,m
 real*8 :: fxa,umax,vmax,wmax,kmax,lmax,mmax,eint,fint,mflint,mfl_max
 real*8 :: cfl_x,cfl_y,cfl_z,cfl_k,cfl_l,cfl_m
 real*8 :: bfl_xs,bfl_xe,bfl_ys,bfl_ye,bfl_zs,bfl_ze,bfl_ks,bfl_ke,bfl_ls,bfl_le,bfl_ms,bfl_me
 real*8 :: E_int(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: meanfl_int(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_ks(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_ke(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_ls(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_le(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_ms(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 real*8 :: diss_me(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx)
 
 real*8,allocatable :: aloc(:,:,:,:,:), bloc(:,:,:), cloc(:,:,:,:,:,:)
  
 if (my_pe==0) then
  print*,' '
  print'(a,f15.2,a)',' writing snapshot at ',itt*dt,' s'
 endif

 ! calculate integral values
 umax=0;vmax=0;wmax=0;kmax=0;lmax=0;mmax=0;eint=0;fint=0;mflint=0;mfl_max=0
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
       fxa = dxt(i)*dyt(j)*dzt(k)*dkt(l)*dlt(n)*dmt(m)
       eint   = eint   + dE(i,j,k,l,n,m)*fxa
       mflint = mflint + meanfl(i,j,k,l,n,m)*fxa  
       mfl_max = max( mfl_max,  abs(meanfl(i,j,k,l,n,m))/max(1e-12,E(i,j,k,l,n,m)) )
       umax = max(umax,abs(xdot(i,j,k,l,n,m))) 
       vmax = max(vmax,abs(ydot(i,j,k,l,n,m))) 
       wmax = max(wmax,abs(zdot(i,j,k,l,n,m))) 
       kmax = max(kmax,abs(kdot(i,j,k,l,n,m))) 
       lmax = max(lmax,abs(ldot(i,j,k,l,n,m))) 
       mmax = max(mmax,abs(mdot(i,j,k,l,n,m))) 
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 ! communicate over pes
 call global_sum(eint)
 call global_sum(mflint)
 call global_max(mfl_max)
 call global_max(umax);
 call global_max(vmax);
 call global_max(wmax)
 call global_max(kmax);
 call global_max(lmax);
 call global_max(mmax)
 
 bfl_xs=0;bfl_xe=0;bfl_ys=0;bfl_ye=0;bfl_zs=0;bfl_ze=0
 bfl_ks=0;bfl_ke=0;bfl_ls=0;bfl_le=0;bfl_ms=0;bfl_me=0
 
 ! integrate energy in wavenumber space
 E_int = 0.;meanfl_int = 0
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
       E_int(i,j,k) = E_int(i,j,k) + E(i,j,k,l,n,m)*dlt(n)*dmt(m)*dkt(l)
       meanfl_int(i,j,k) = meanfl_int(i,j,k) + meanfl(i,j,k,l,n,m)*dkt(l)*dlt(n)*dmt(m)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 call wavenumber_sum(E_int)
 call wavenumber_sum(meanfl_int)
 
 ! calculate integrated energy fluxes at x-boundaries
 if (my_blk_x==1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       bfl_xs = bfl_xs + flux_x(0,j,k,l,n,m)*dyt(j)*dzt(k)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif     
 if (my_blk_x==n_pes_x) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe 
       bfl_xe = bfl_xe + flux_x(nx,j,k,l,n,m)*dyt(j)*dzt(k)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif  
 
 ! calculate integrated energy fluxes at y-boundaries
 if (my_blk_y==1)  then 
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe 
      do i=xs_pe,xe_pe 
       bfl_ys = bfl_ys + flux_y(i,0,k,l,n,m)*dxt(i)*dzt(k)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif  
 if (my_blk_y==n_pes_y) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe 
      do i=xs_pe,xe_pe  
       bfl_ye = bfl_ye + flux_y(i,ny,k,l,n,m)*dxt(i)*dzt(k)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif  
 
 ! calculate integrated energy fluxes at z-boundaries
 if (my_blk_z==1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe
       bfl_zs = bfl_zs + flux_z(i,j,0,l,n,m)*dxt(i)*dyt(j)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif      
 if (my_blk_z==n_pes_z) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe
       bfl_ze = bfl_ze + flux_z(i,j,nz,l,n,m)*dxt(i)*dyt(j)*dkt(l)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif

 diss_ks=0; diss_ke=0.
 ! calculate integrated energy fluxes at k-boundaries
 if (my_blk_k==1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
       bfl_ks = bfl_ks + flux_k(i,j,k,0,n,m)*dxt(i)*dyt(j)*dzt(k)*dlt(n)*dmt(m) 
       diss_ks(i,j,k) = diss_ks(i,j,k) + flux_k(i,j,k,0,n,m)*dlt(n)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 if (my_blk_k==n_pes_k) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
       bfl_ke = bfl_ke + flux_k(i,j,k,nk,n,m)*dxt(i)*dyt(j)*dzt(k)*dlt(n)*dmt(m)  
       diss_ke(i,j,k) = diss_ke(i,j,k) + flux_k(i,j,k,nk,n,m)*dlt(n)*dmt(m)  
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 
 diss_ls=0; diss_le=0.
 ! calculate integrated energy fluxes at l-boundaries
 if (my_blk_l==1) then
  do m=ms_pe,me_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
       bfl_ls = bfl_ls + flux_l(i,j,k,l,0,m)*dxt(i)*dyt(j)*dzt(k)*dkt(l)*dmt(m) 
       diss_ls(i,j,k) = diss_ls(i,j,k) + flux_l(i,j,k,l,0,m)*dkt(l)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif 
 if (my_blk_l==n_pes_l) then
  do m=ms_pe,me_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
        bfl_le = bfl_le + flux_l(i,j,k,l,nl,m)*dxt(i)*dyt(j)*dzt(k)*dkt(l)*dmt(m) 
        diss_le(i,j,k) = diss_le(i,j,k) + flux_l(i,j,k,l,nl,m)*dkt(l)*dmt(m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 endif 
 
 diss_ms = 0.; diss_me=0.
 ! calculate integrated energy fluxes at m-boundaries
 if (my_blk_m==1) then 
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe
       bfl_ms = bfl_ms + flux_m(i,j,k,l,n,0)*dxt(i)*dyt(j)*dzt(k)*dkt(l)*dlt(n) 
       diss_ms(i,j,k) =  diss_ms(i,j,k) + flux_m(i,j,k,l,n,0)*dkt(l)*dlt(n)
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 if (my_blk_m==n_pes_m) then
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe
        bfl_me = bfl_me + flux_m(i,j,k,l,n,nm)*dxt(i)*dyt(j)*dzt(k)*dkt(l)*dlt(n)
        diss_me(i,j,k) =  diss_me(i,j,k) + flux_m(i,j,k,l,n,nm)*dkt(l)*dlt(n)
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 
 ! communicate over pes
 call global_sum(bfl_xs);call global_sum(bfl_xe)
 call global_sum(bfl_ys);call global_sum(bfl_ye)
 call global_sum(bfl_zs);call global_sum(bfl_ze)
 call global_sum(bfl_ks);call global_sum(bfl_ke)
 call global_sum(bfl_ls);call global_sum(bfl_le)
 call global_sum(bfl_ms);call global_sum(bfl_me)

 call wavenumber_sum(diss_ke)
 call wavenumber_sum(diss_ks)
 call wavenumber_sum(diss_le)
 call wavenumber_sum(diss_ls)
 call wavenumber_sum(diss_me)
 call wavenumber_sum(diss_ms)
 
 ! check maximal CFL numbers
 cfl_x=0;cfl_y=0.;cfl_z=0;cfl_k=0;cfl_l=0;cfl_m=0
 do i=xs_pe,xe_pe
  cfl_x = max(cfl_x,umax*dt/dxt(i))
 enddo 
 do i=ys_pe,ye_pe
  cfl_y = max(cfl_y,vmax*dt/dyt(i))
 enddo 
 do i=zs_pe,ze_pe
  cfl_z = max(cfl_z,wmax*dt/dzt(i))
 enddo 
 do i=ks_pe,ke_pe
  cfl_k = max(cfl_k,kmax*dt/dkt(i))
 enddo 
 do i=ls_pe,le_pe
  cfl_l = max(cfl_l,lmax*dt/dlt(i))
 enddo 
 do i=ms_pe,me_pe
  cfl_m = max(cfl_m,mmax*dt/dmt(i))
 enddo 
 
 ! print out
 if (my_pe==0) then
     print'(a,f8.2,f8.2,f8.2,f8.2,f8.2,f8.2)',' max CFL in x,y,z,k,l,m-dir = ', cfl_x,cfl_y,cfl_z,cfl_k,cfl_l,cfl_m 
     print'(a,f8.2)'       ,' mean flow interaction scale ',dt/mfl_max
     print'(a,e14.6)'      ,' total dE          : ',eint
     print'(a,e14.6)'      ,' dE by mean-flow   : ',mflint
     print'(a,e14.6,e14.6)',' by flux at min/max x  : ',bfl_xs,bfl_xe
     print'(a,e14.6,e14.6)',' by flux at min/max y  : ',bfl_ys,bfl_ye
     !print'(a,e14.6,e14.6)',' by flux at min/max z  : ',bfl_zs,bfl_ze
     print'(a,e14.6,e14.6)',' by flux at min/max k  : ',bfl_ks,bfl_ke
     print'(a,e14.6,e14.6)',' by flux at min/max l  : ',bfl_ls,bfl_le
     print'(a,e14.6,e14.6)',' by flux at min/max m  : ',bfl_ms,bfl_me    
     print'(a,e14.6)',' residual : ',eint+mflint+bfl_xe+bfl_ye+bfl_ze+bfl_ke+bfl_le+bfl_me &
                                                -bfl_xs-bfl_ys-bfl_zs-bfl_ks-bfl_ls-bfl_ms      
 endif

 call tic('diag_write')
 
 ! write 6D arrays to file

 if (enable_write_6D_single) then
 
  ! only first PE writes to a single output file
  if (my_pe==0) then
   iret=nf_open('snap6D.cdf',NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   iret=nf_inq_varid(ncid,'Time',itimeid)
   ilen=ilen+1
   print*,'time steps in file snap6D.cdf : ',ilen
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,itt*dt)
   
  endif
  !call pe0_bcast_int(ilen,1)

  if (enable_write_6D_single_himem) then
  
   ! this version needs a lot of memory
   allocate( cloc(nx,ny,nz,nk,nl,nm) ,stat=i); 
   if (i/=0) then
    call halt_stop('ERROR: allocate failed in diagnose.f90 (1) ')
   else 
    cloc=0.
    cloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) = &
       E(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe)
    call pe0_recv_arr(cloc)
    if (my_pe==0) then
     iret = nf_inq_varid(ncid,'E',id)
     iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1,ilen/),(/nx,ny,nz,nk,nl,nm,1/),cloc)
    endif
  
    cloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) = &
       meanfl(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe)
    call pe0_recv_arr(cloc)
    if (my_pe==0) then
     iret = nf_inq_varid(ncid,'meanfl',id)
     iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1,ilen/),(/nx,ny,nz,nk,nl,nm,1/),cloc)
    endif 
    deallocate(cloc)
   endif
   
  else
  
   ! this version needs less memory but is slower
   allocate( aloc(nx,ny,nz,nk,nl),stat = i ); 
   if (i/=0) then
    call halt_stop('ERROR: allocate failed in diagnose.f90 (2) ')
   else 
    aloc=0.
    do i=1,nm
     if (i>=ms_pe.and.i<=me_pe) &
      aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       E(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
     call pe0_recv_arr_mdim(aloc,i)
     if (my_pe==0) then
      iret = nf_inq_varid(ncid,'E',id)
      iret = nf_put_vara_double(ncid,id,(/1,1,1,1,1,i,ilen/),(/nx,ny,nz,nk,nl,1,1/),aloc)
     endif
    enddo 
 
    do i=1,nm
     if (i>=ms_pe.and.i<=me_pe) &
      aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
        meanfl(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
     call pe0_recv_arr_mdim(aloc,i)
     if (my_pe==0) then
      iret = nf_inq_varid(ncid,'meanfl',id)
      iret = nf_put_vara_double(ncid,id,(/1,1,1,1,1,i,ilen/),(/nx,ny,nz,nk,nl,1,1/),aloc)
     endif
    enddo 
    deallocate(aloc)
   endif 
  endif
  
  if (my_pe==0) call ncclos (ncid, iret)  
   
 else
  ! each PE writes to individual file
  iret = nf_open(netcdf_6D_file,NF_WRITE,ncid)
  iret = nf_set_fill(ncid, NF_NOFILL, iret)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  iret=nf_inq_dimid(ncid,'Time',itdimid)
  iret=nf_inq_dimlen(ncid, itdimid,ilen)
  iret=nf_inq_varid(ncid,'Time',itimeid)
  ilen=ilen+1
  if (my_pe==0) print*,'time steps in file ',netcdf_6D_file(1:len_trim(netcdf_6D_file)),' : ',ilen
  iret= nf_put_vara_double(ncid,itimeid,ilen,1,itt*dt)
  iret = nf_inq_varid(ncid,'E',id)
  iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1,ilen/),(/x_blk,y_blk,z_blk,k_blk,l_blk,m_blk,1/), &
         E(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
  iret = nf_inq_varid(ncid,'meanfl',id)       
  iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1,ilen/),(/x_blk,y_blk,z_blk,k_blk,l_blk,m_blk,1/), &
         meanfl(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
  call ncclos (ncid, iret)  
 endif



 ! write physical arrays to 3D netcdf file, only first PE writes to a single file
 if (my_pe==0) then
   iret=nf_open('snap3D.cdf',NF_WRITE,ncid)
   iret=nf_set_fill(ncid, NF_NOFILL, iret)
   iret=nf_inq_dimid(ncid,'Time',itdimid)
   iret=nf_inq_dimlen(ncid, itdimid,ilen)
   iret=nf_inq_varid(ncid,'Time',itimeid)
   ilen=ilen+1
   print*,'time steps in file snap3D.cdf : ',ilen
   iret= nf_put_vara_double(ncid,itimeid,ilen,1,itt*dt)
 endif

 allocate( bloc(nx,ny,nz) ); bloc=0.
 
 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  E_int(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'E',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  meanfl_int(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'meanfl',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 
 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_ks(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_ks',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_ke(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_ke',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_ls(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_ls',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_le(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_le',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif
 
 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_ms(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_ms',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif

 bloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe) =  diss_me(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe)
 call pe0_recv_arr_phys(bloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'diss_me',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1,ilen/),(/nx,ny,nz,1/),bloc)
 endif
 
 deallocate(bloc)
 
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'mflint',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,mflint)
   iret=nf_inq_varid(ncid,'bflux_xs',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_xs)
   iret=nf_inq_varid(ncid,'bflux_xe',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_xe)
   iret=nf_inq_varid(ncid,'bflux_ys',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ys)
   iret=nf_inq_varid(ncid,'bflux_ye',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ye)
   iret=nf_inq_varid(ncid,'bflux_ks',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ks)
   iret=nf_inq_varid(ncid,'bflux_ke',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ke)
   iret=nf_inq_varid(ncid,'bflux_ls',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ls)
   iret=nf_inq_varid(ncid,'bflux_le',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_le)
   iret=nf_inq_varid(ncid,'bflux_ms',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_ms)
   iret=nf_inq_varid(ncid,'bflux_me',id)
   iret= nf_put_vara_double(ncid,id,ilen,1,bfl_me)
 endif
 
 if (my_pe==0)  then
  call ncclos (ncid, iret)
  print*,'done with diagnostics for now'
 endif 
 call toc('diag_write')
end subroutine diagnose




