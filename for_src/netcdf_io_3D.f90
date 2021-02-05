

subroutine init_snap3D_cdf
! create output file for 3D fields only
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,xtdim,ytdim,ztdim
 integer :: xudim,yudim,zudim,itimeid,itimedim,id,i,j,k
 character :: name*24, unit*16
 real*8,allocatable :: aloc(:,:,:)
 
 call fortran_barrier()
 if (my_pe==0) then
  iret = NF_CREATE ('snap3D.cdf',OR(NF_SHARE,NF_64BIT_OFFSET),ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xtdim  = ncddef(ncid, 'xt', nx, iret)
  xudim  = ncddef(ncid, 'xu', nx, iret)
  ytdim  = ncddef(ncid, 'yt', ny, iret)
  yudim  = ncddef(ncid, 'yu', ny, iret)
  ztdim  = ncddef(ncid, 'zt', nz, iret)
  zudim  = ncddef(ncid, 'zu', nz, iret)
  
  iTimedim  = ncddef(ncid, 'Time', nf_unlimited, iret)
  itimeid  = ncvdef (ncid,'Time', NF_DOUBLE,1,itimedim,iret)
  name = 'Time '; unit = 'seconds'
  call ncaptc(ncid, itimeid, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, itimeid, 'units',     NCCHAR, 16, unit, iret) 
  call ncaptc(ncid, iTimeid,'time_origin',NCCHAR, 20,'01-JAN-1900 00:00:00', iret)

  id  = ncvdef (ncid,'xt', NF_DOUBLE,1,xtdim,iret)
  name = 'distance'; unit = 'm'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'xu', NF_DOUBLE,1,xudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'yt', NF_DOUBLE,1,ytdim,iret)
  name = 'distance'; unit = 'm'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'yu', NF_DOUBLE,1,yudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'zt', NF_DOUBLE,1,ztdim,iret)
  name = 'distance'; unit = 'm'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'zu', NF_DOUBLE,1,zudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  

  

 
  id  = ncvdef (ncid,'Nsqr',NF_DOUBLE,3,(/xtdim,ytdim,ztdim/),iret)
  name = 'buoyancy frequency'; unit = '1/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'U',NF_DOUBLE,3,(/xtdim,ytdim,ztdim/),iret)
  name = 'mean flow'; unit = 'm/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'V',NF_DOUBLE,3,(/xtdim,ytdim,ztdim/),iret)
  name = 'mean flow'; unit = 'm/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'E',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy density'; unit = 'm^2/s^2'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'meanfl',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_ks',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_ke',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_ls',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_le',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_ms',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'diss_me',NF_DOUBLE,4,(/xtdim,ytdim,ztdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^2/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_xs',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across min x'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_xe',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across max x'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_ys',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across min y'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_ye',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across max y'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_ks',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across min k'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_ke',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across max k'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'bflux_ls',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across min l'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_le',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across max l'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_ms',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across min m'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'bflux_me',NF_DOUBLE,1,itimedim,iret)
  name = 'Energy flux across max m'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'mflint',NF_DOUBLE,1,itimedim,iret)
  name = 'mean flow energy flux'; unit = ' '
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'f0',NF_DOUBLE,1,f0)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'k_fixed',NF_DOUBLE,1,k_fixed)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'l_fixed',NF_DOUBLE,1,l_fixed)
  iret = NF_PUT_ATT_DOUBLE(ncid,NF_GLOBAL,'m_fixed',NF_DOUBLE,1,m_fixed)

  call ncendf(ncid, iret)
 endif
 
 allocate( aloc(nx,ny,nz) ); aloc=0.
 
 do k=zs_pe,ze_pe 
  do j=ys_pe,ye_pe 
   do i=xs_pe,xe_pe 
    aloc(i,j,k) = Nsqr(i,j,k)
   enddo
  enddo 
 enddo 
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'Nsqr',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1/),(/nx,ny,nz/),aloc)
 endif 

 do k=zs_pe,ze_pe 
  do j=ys_pe,ye_pe 
   do i=xs_pe,xe_pe 
    aloc(i,j,k) = U(i,j,k)
   enddo
  enddo 
 enddo 
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'U',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1/),(/nx,ny,nz/),aloc)
 endif 
 
 do k=zs_pe,ze_pe 
  do j=ys_pe,ye_pe 
   do i=xs_pe,xe_pe 
    aloc(i,j,k) = V(i,j,k)
   enddo
  enddo 
 enddo 
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'V',id)
   iret= nf_put_vara_double(ncid,id,(/1,1,1/),(/nx,ny,nz/),aloc)
 endif 
 

 do i=xs_pe,xe_pe 
  aloc(i,1,:) = xt(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'xt',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/nx/),aloc(:,1,1))
 endif 
 
 do i=xs_pe,xe_pe 
  aloc(i,1,:) = xu(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'xu',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/nx/),aloc(:,1,1))
 endif 

 do i=ys_pe,ye_pe 
  aloc(1,i,:) = yt(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'yt',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/ny/),aloc(1,:,1))
 endif 
 
 do i=ys_pe,ye_pe 
  aloc(1,i,:) = yu(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'yu',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/ny/),aloc(1,:,1))
 endif 

 do i=zs_pe,ze_pe 
  aloc(1,:,i) = zt(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'zt',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/nz/),aloc(1,1,:))
 endif 
 
 do i=zs_pe,ze_pe
  aloc(1,:,i) = zu(i)
 enddo
 call pe0_recv_arr_phys(aloc)
 if (my_pe==0) then
   iret=nf_inq_varid(ncid,'zu',id)
   iret= nf_put_vara_double(ncid,id,(/1/),(/nz/),aloc(1,1,:))
 endif 

 if (my_pe==0) then
  call ncclos (ncid, iret)
 endif
 call fortran_barrier()

 deallocate( aloc)
end subroutine init_snap3D_cdf






subroutine dvcdf(ncid,ivarid,name,iname,unit,iunit,spval)
 implicit none
 include "netcdf.inc"
 integer ncid,ivarid,iname,iunit,iret
 character (len=*) :: name, unit
 real*8 :: spval
 real*8 :: vv
 vv=spval
 call ncaptc(ncid,ivarid, 'long_name', NCCHAR,iname , name, iret) 
     if (iret.ne.0) print*,nf_strerror(iret)
 call ncaptc(ncid,ivarid, 'units',     NCCHAR,iunit, unit, iret) 
        if (iret.ne.0) print*,nf_strerror(iret)
 call ncapt (ncid,ivarid, 'missing_value',NCDOUBLE,1,vv,iret)
        if (iret.ne.0) print*,nf_strerror(iret)
 call ncapt (ncid,ivarid, '_FillValue', NCDOUBLE, 1,vv, iret)
      if (iret.ne.0) print*,nf_strerror(iret)
end subroutine dvcdf

