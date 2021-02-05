
subroutine init_snap6D_cdf_many
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,xtdim,ytdim,ztdim,ktdim,ltdim,mtdim,id
 integer :: xudim,yudim,zudim,kudim,ludim,mudim,itimedim,itimeid
 character :: name*24, unit*16

 write(netcdf_6D_file,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a,i5,a)')  &
  'snap6D_i=',xs_pe,':',xe_pe,'_j=',ys_pe,':',ye_pe,'_k=',zs_pe,':',ze_pe, &
        '_l=',ks_pe,':',ke_pe,'_n=',ls_pe,':',le_pe,'_m=',ms_pe,':',me_pe,'.cdf'
 call replace_space_zero(netcdf_6D_file)
 
 if (my_pe==0) print*,'creating file ',netcdf_6D_file(1:len_trim(netcdf_6D_file))

 !iret = NF_CREATE (netcdf_6D_file,OR(NF_SHARE,NF_64BIT_OFFSET),ncid)
 iret = NF_CREATE (netcdf_6D_file,IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
 iret=nf_set_fill(ncid, NF_NOFILL, iret)
 xtdim  = ncddef(ncid, 'xt', x_blk, iret)
 xudim  = ncddef(ncid, 'xu', x_blk, iret)
 ytdim  = ncddef(ncid, 'yt', y_blk, iret)
 yudim  = ncddef(ncid, 'yu', y_blk, iret)
 ztdim  = ncddef(ncid, 'zt', z_blk, iret)
 zudim  = ncddef(ncid, 'zu', z_blk, iret)
 ktdim  = ncddef(ncid, 'kt', k_blk, iret)
 kudim  = ncddef(ncid, 'ku', k_blk, iret)
 ltdim  = ncddef(ncid, 'lt', l_blk, iret)
 ludim  = ncddef(ncid, 'lu', l_blk, iret)
 mtdim  = ncddef(ncid, 'mt', m_blk, iret)
 mudim  = ncddef(ncid, 'mu', m_blk, iret)
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
 id  = ncvdef (ncid,'dxt', NF_DOUBLE,1,xtdim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'xu', NF_DOUBLE,1,xudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

 id  = ncvdef (ncid,'yt', NF_DOUBLE,1,ytdim,iret)
 name = 'distance'; unit = 'm'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'dyt', NF_DOUBLE,1,ytdim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'yu', NF_DOUBLE,1,yudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

 id  = ncvdef (ncid,'zt', NF_DOUBLE,1,ztdim,iret)
 name = 'distance'; unit = 'm'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'dzt', NF_DOUBLE,1,zudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'zu', NF_DOUBLE,1,zudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  
 id  = ncvdef (ncid,'kt', NF_DOUBLE,1,ktdim,iret)
 name = 'wavenumber'; unit = '1/m'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'ku', NF_DOUBLE,1,kudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'dkt', NF_DOUBLE,1,kudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 
 id  = ncvdef (ncid,'lt', NF_DOUBLE,1,ltdim,iret)
 name = 'wavenumber'; unit = '1/m'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'lu', NF_DOUBLE,1,ludim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'dlt', NF_DOUBLE,1,ludim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 
 id  = ncvdef (ncid,'mt', NF_DOUBLE,1,mtdim,iret)
 name = 'wavenumber'; unit = '1/m'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'mu', NF_DOUBLE,1,mudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 id  = ncvdef (ncid,'dmt', NF_DOUBLE,1,mudim,iret)
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
 
 id  = ncvdef (ncid,'xdot',NF_DOUBLE,6,(/xudim,ytdim, ztdim, ktdim,ltdim, mtdim/),iret)
 name = 'hor. group velocity'; unit = 'm/s'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'ydot',NF_DOUBLE,6,(/xtdim,yudim, ztdim, ktdim,ltdim, mtdim/),iret)
 name = 'hor. group velocity'; unit = 'm/s'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'zdot',NF_DOUBLE,6,(/xtdim,ytdim, zudim, ktdim,ltdim, mtdim/),iret)
 name = 'vertical group velocity'; unit = 'm/s'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'kdot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, kudim,ltdim, mtdim/),iret)
 name = 'hor. wavenr. change'; unit = '1/ms'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'ldot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ludim, mtdim/),iret)
 name = 'hor. wavenr. change'; unit = '1/ms'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'mdot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim,ktdim, ltdim, mudim/),iret)
 name = 'vert. wavenr. change'; unit = '1/ms'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'omega_i',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim/),iret)
 name = 'intrinsic frequency'; unit = '1/s'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 id  = ncvdef (ncid,'omega_e',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim/),iret)
 name = 'extrinsic frequency'; unit = '1/s'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

 id  = ncvdef (ncid,'E',NF_DOUBLE,7,(/xtdim,ytdim,ztdim,ktdim,ltdim,mtdim,itimedim/),iret)
 name = 'energy density'; unit = 'm^5/s^2'  ! int E dm dk dl = m^2/^2 -> 
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

 id  = ncvdef (ncid,'meanfl',NF_DOUBLE,7,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim,itimedim/),iret)
 name = 'energy change'; unit = 'm^5/s^3'
 call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
 call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
 call ncendf(ncid, iret)

 iret = nf_inq_varid(ncid,'xt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/x_blk/),xt(xs_pe:xe_pe))
 iret = nf_inq_varid(ncid,'xu',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/x_blk/),xu(xs_pe:xe_pe))
 iret = nf_inq_varid(ncid,'dxt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/x_blk/),dxt(xs_pe:xe_pe))
 
 iret = nf_inq_varid(ncid,'yt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/y_blk/),yt(ys_pe:ye_pe))
 iret = nf_inq_varid(ncid,'yu',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/y_blk/),yu(ys_pe:ye_pe))
 iret = nf_inq_varid(ncid,'dyt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/y_blk/),dyt(ys_pe:ye_pe))
 
 iret = nf_inq_varid(ncid,'zt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/z_blk/),zt(zs_pe:ze_pe))
 iret = nf_inq_varid(ncid,'zu',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/z_blk/),zu(zs_pe:ze_pe))
 iret = nf_inq_varid(ncid,'dzt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/z_blk/),dzt(zs_pe:ze_pe))
 
 iret = nf_inq_varid(ncid,'kt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/k_blk/),kt(ks_pe:ke_pe))
 iret = nf_inq_varid(ncid,'ku',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/k_blk/),ku(ks_pe:ke_pe))
 iret = nf_inq_varid(ncid,'dkt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/k_blk/),dkt(ks_pe:ke_pe))

 iret = nf_inq_varid(ncid,'lt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/l_blk/),lt(ls_pe:le_pe))
 iret = nf_inq_varid(ncid,'lu',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/l_blk/),lu(ls_pe:le_pe))
 iret = nf_inq_varid(ncid,'dlt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/l_blk/),dlt(ls_pe:le_pe))
 
 iret = nf_inq_varid(ncid,'mt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/m_blk/),mt(ms_pe:me_pe))
 iret = nf_inq_varid(ncid,'mu',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/m_blk/),mu(ms_pe:me_pe))
 iret = nf_inq_varid(ncid,'dmt',id)
 iret = nf_put_vara_double(ncid,id,(/1/),(/m_blk/),dmt(ms_pe:me_pe))
 
 iret = nf_inq_varid(ncid,'xdot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       xdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'ydot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       ydot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'zdot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       zdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'kdot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       kdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'ldot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       ldot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'mdot',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       mdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 iret = nf_inq_varid(ncid,'omega_i',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       omega_i(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )

 iret = nf_inq_varid(ncid,'omega_e',id)
 iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,1/),(/x_blk,y_blk,z_blk,k_blk,l_blk/), &
       omega_e(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) )
 
 call ncclos (ncid, iret)  

end subroutine init_snap6D_cdf_many




subroutine init_snap6D_cdf_single
 use main_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,xtdim,ytdim,ztdim,ktdim,ltdim,mtdim,id,n,i
 integer :: xudim,yudim,zudim,kudim,ludim,mudim,itimedim,itimeid
 character :: name*24, unit*16
 real*8, dimension(1-onx:nx+onx) :: dxt_gl,xt_gl,xu_gl
 real*8, dimension(1-onx:ny+onx) :: dyt_gl,yt_gl,yu_gl
 real*8, dimension(1-onx:nz+onx) :: dzt_gl,zt_gl,zu_gl
 real*8, dimension(1-onx:nk+onx) :: dkt_gl,kt_gl,ku_gl
 real*8, dimension(1-onx:nl+onx) :: dlt_gl,lt_gl,lu_gl
 real*8, dimension(1-onx:nm+onx) :: dmt_gl,mt_gl,mu_gl 
 include "mpif.h"
 integer :: tag = 99, ierr, sender, ind(3)
 integer,dimension(MPI_STATUS_SIZE)  :: Status
 real*8,allocatable :: aloc(:,:,:,:,:)
  
 
 if (my_pe==0) then
  iret = NF_CREATE ('snap6D.cdf',IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
  !iret = NF_CREATE ('snap6D.cdf',OR(NF_SHARE,NF_64BIT_OFFSET),ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xtdim  = ncddef(ncid, 'xt', nx, iret)
  xudim  = ncddef(ncid, 'xu', nx, iret)
  ytdim  = ncddef(ncid, 'yt', ny, iret)
  yudim  = ncddef(ncid, 'yu', ny, iret)
  ztdim  = ncddef(ncid, 'zt', nz, iret)
  zudim  = ncddef(ncid, 'zu', nz, iret)
  ktdim  = ncddef(ncid, 'kt', nk, iret)
  kudim  = ncddef(ncid, 'ku', nk, iret)
  ltdim  = ncddef(ncid, 'lt', nl, iret)
  ludim  = ncddef(ncid, 'lu', nl, iret)
  mtdim  = ncddef(ncid, 'mt', nm, iret)
  mudim  = ncddef(ncid, 'mu', nm, iret)
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
  id  = ncvdef (ncid,'dxt', NF_DOUBLE,1,xudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'yt', NF_DOUBLE,1,ytdim,iret)
  name = 'distance'; unit = 'm'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'yu', NF_DOUBLE,1,yudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'dyt', NF_DOUBLE,1,yudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'zt', NF_DOUBLE,1,ztdim,iret)
  name = 'distance'; unit = 'm'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'zu', NF_DOUBLE,1,zudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'dzt', NF_DOUBLE,1,zudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'kt', NF_DOUBLE,1,ktdim,iret)
  name = 'wavenumber'; unit = '1/m'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'ku', NF_DOUBLE,1,kudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'dkt', NF_DOUBLE,1,kudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'lt', NF_DOUBLE,1,ltdim,iret)
  name = 'wavenumber'; unit = '1/m'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'lu', NF_DOUBLE,1,ludim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'dlt', NF_DOUBLE,1,ludim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'mt', NF_DOUBLE,1,mtdim,iret)
  name = 'wavenumber'; unit = '1/m'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'mu', NF_DOUBLE,1,mudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  id  = ncvdef (ncid,'dmt', NF_DOUBLE,1,mudim,iret)
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'xdot',NF_DOUBLE,6,(/xudim,ytdim, ztdim, ktdim,ltdim, mtdim/),iret)
  name = 'hor. group velocity'; unit = 'm/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'ydot',NF_DOUBLE,6,(/xtdim,yudim, ztdim, ktdim,ltdim, mtdim/),iret)
  name = 'hor. group velocity'; unit = 'm/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'zdot',NF_DOUBLE,6,(/xtdim,ytdim, zudim, ktdim,ltdim, mtdim/),iret)
  name = 'vertical group velocity'; unit = 'm/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'kdot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, kudim,ltdim, mtdim/),iret)
  name = 'hor. wavenr. change'; unit = '1/ms'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'ldot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ludim, mtdim/),iret)
  name = 'hor. wavenr. change'; unit = '1/ms'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'mdot',NF_DOUBLE,6,(/xtdim,ytdim,ztdim,ktdim, ltdim, mudim/),iret)
  name = 'vert. wavenr. change'; unit = '1/ms'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'omega_i',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim/),iret)
  name = 'intrinsic frequency'; unit = '1/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  id  = ncvdef (ncid,'omega_e',NF_DOUBLE,6,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim/),iret)
  name = 'extrinsic frequency'; unit = '1/s'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'E',NF_DOUBLE,7,(/xtdim,ytdim,ztdim,ktdim,ltdim,mtdim,itimedim/),iret)
  name = 'energy density'; unit = 'm^5/s^2'  ! int E dm dk dl = m^2/^2 -> 
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 

  id  = ncvdef (ncid,'meanfl',NF_DOUBLE,7,(/xtdim,ytdim,ztdim, ktdim,ltdim, mtdim,itimedim/),iret)
  name = 'energy change'; unit = 'm^5/s^3'
  call ncaptc(ncid, id, 'long_name', NCCHAR, 24, name, iret) 
  call ncaptc(ncid, id, 'units',     NCCHAR, 16, unit, iret) 
  
  call ncendf(ncid, iret)
 endif
 
 xt_gl(xs_pe:xe_pe) = xt(xs_pe:xe_pe)
 xu_gl(xs_pe:xe_pe) = xu(xs_pe:xe_pe)
 dxt_gl(xs_pe:xe_pe) = dxt(xs_pe:xe_pe)
 do n=1,n_pes_x-1
   if (my_blk_x==n+1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/xs_pe,xe_pe,x_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(xt(xs_pe:xe_pe),x_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(xu(xs_pe:xe_pe),x_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dxt(xs_pe:xe_pe),x_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(xt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(xu_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dxt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
 
 yt_gl(ys_pe:ye_pe) = yt(ys_pe:ye_pe)
 yu_gl(ys_pe:ye_pe) = yu(ys_pe:ye_pe)
 dyt_gl(ys_pe:ye_pe) = dyt(ys_pe:ye_pe)
 do n=1,n_pes_y-1
   if (my_blk_x==1.and.my_blk_y==n+1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/ys_pe,ye_pe,y_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(yt(ys_pe:ye_pe),y_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(yu(ys_pe:ye_pe),y_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dyt(ys_pe:ye_pe),y_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(yt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(yu_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dyt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo

 zt_gl(zs_pe:ze_pe) = zt(zs_pe:ze_pe)
 zu_gl(zs_pe:ze_pe) = zu(zs_pe:ze_pe)
 dzt_gl(zs_pe:ze_pe) = dzt(zs_pe:ze_pe)
 do n=1,n_pes_z-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==n+1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/zs_pe,ze_pe,z_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(zt(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(zu(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dzt(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(zt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(zu_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dzt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
  
 kt_gl(ks_pe:ke_pe) = kt(ks_pe:ke_pe)
 ku_gl(ks_pe:ke_pe) = ku(ks_pe:ke_pe)
 dkt_gl(ks_pe:ke_pe) = dkt(ks_pe:ke_pe)
 do n=1,n_pes_k-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==n+1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/ks_pe,ke_pe,k_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(kt(ks_pe:ke_pe),k_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(ku(ks_pe:ke_pe),k_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dkt(ks_pe:ke_pe),k_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(kt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(ku_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dkt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
 
 lt_gl(ls_pe:le_pe) = lt(ls_pe:le_pe)
 lu_gl(ls_pe:le_pe) = lu(ls_pe:le_pe)
 dlt_gl(ls_pe:le_pe) = dlt(ls_pe:le_pe)
 do n=1,n_pes_l-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==n+1.and.my_blk_m==1) then
    call mpi_send((/ls_pe,le_pe,l_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(lt(ls_pe:le_pe),l_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(lu(ls_pe:le_pe),l_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dlt(ls_pe:le_pe),l_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z*n_pes_k
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(lt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(lu_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dlt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
 
 mt_gl(ms_pe:me_pe) = mt(ms_pe:me_pe)
 mu_gl(ms_pe:me_pe) = mu(ms_pe:me_pe)
 dmt_gl(ms_pe:me_pe) = dmt(ms_pe:me_pe)
 do n=1,n_pes_m-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==n+1) then
    call mpi_send((/ms_pe,me_pe,m_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(mt(ms_pe:me_pe),m_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(mu(ms_pe:me_pe),m_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dmt(ms_pe:me_pe),m_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y*n_pes_z*n_pes_k*n_pes_l
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(mt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(mu_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dmt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
 
 if (my_pe==0) then 
  iret = nf_inq_varid(ncid,'xt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nx/),xt_gl(1:nx))
  iret = nf_inq_varid(ncid,'xu',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nx/),xu_gl(1:nx))
  iret = nf_inq_varid(ncid,'dxt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nx/),dxt_gl(1:nx))
  
  iret = nf_inq_varid(ncid,'yt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/ny/),yt_gl(1:ny))
  iret = nf_inq_varid(ncid,'yu',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/ny/),yu_gl(1:ny))
  iret = nf_inq_varid(ncid,'dyt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/ny/),dyt_gl(1:ny))
  
  iret = nf_inq_varid(ncid,'zt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nz/),zt_gl(1:nz))
  iret = nf_inq_varid(ncid,'zu',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nz/),zu_gl(1:nz))
  iret = nf_inq_varid(ncid,'dzt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nz/),dzt_gl(1:nz))
  
  iret = nf_inq_varid(ncid,'kt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nk/),kt_gl(1:nk))
  iret = nf_inq_varid(ncid,'ku',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nk/),ku_gl(1:nk))
  iret = nf_inq_varid(ncid,'dkt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nk/),dkt_gl(1:nk))
  
  iret = nf_inq_varid(ncid,'lt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nl/),lt_gl(1:nl))
  iret = nf_inq_varid(ncid,'lu',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nl/),lu_gl(1:nl))
  iret = nf_inq_varid(ncid,'dlt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nl/),dlt_gl(1:nl))
  
  iret = nf_inq_varid(ncid,'mt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nm/),mt_gl(1:nm))
  iret = nf_inq_varid(ncid,'mu',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nm/),mu_gl(1:nm))
  iret = nf_inq_varid(ncid,'dmt',id)
  iret = nf_put_vara_double(ncid,id,(/1/),(/nm/),dmt_gl(1:nm))
  
 
 endif


 allocate( aloc(nx,ny,nz,nk,nl) ,stat=i); 
 if (i/=0) then
  call halt_stop('ERROR: allocate failed in netcdf_io_6D.f90 (1) ')
 else 
  aloc=0.
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
        xdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'xdot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       ydot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'ydot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       zdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'zdot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       kdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'kdot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       ldot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'ldot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       mdot(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'mdot',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       omega_i(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'omega_i',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe) = &
       omega_e(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'omega_e',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,1,i/),(/nx,ny,nz,nk,nl,1/),aloc)
   endif
  enddo 
 
  deallocate(aloc)
 endif
 
 if (my_pe==0) call ncclos (ncid, iret) 
end subroutine init_snap6D_cdf_single




subroutine replace_space_zero(name)
      implicit none
      character (len=*) :: name
      integer  :: i
      do i=1,len_trim(name)
          if (name(i:i)==' ')name(i:i)='0'
      enddo
end subroutine replace_space_zero


