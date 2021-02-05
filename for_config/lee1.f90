
! model with z,m, and k dimension with lee wave flux at bottom
! mean flow and N as in lee wave paper
! constant grid spacing in all directions
! rayleigh damping is optional

module config_module
 implicit none
 real*8, allocatable :: flux_lee(:,:,:,:,:)
end module config_module



subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=2
 nx = 1
 ny = 1
 nz = 20!int(20*fac)
 nk = int(24*fac) ! has to be even or one
 nl = 1!int(22*fac) ! has to be even or one
 nm = int(26*fac) ! has to be even or one
 f0 = 1e-4
 dt = 300.0/fac
 Lz = 2000.0
 Ly = 240d3
 Lx = 240d3
 snapint = 3600.
 runlen = snapint*40000
 
 l_fixed = 0.
 
 enable_upper_reflection = .true.
 enable_lower_reflection = .true.
 
 enable_write_6D_single = .true.
 enable_write_6D_single_himem = .true.
end subroutine set_parameter


subroutine set_grid
 use main_module   
 implicit none

 dxt = Lx/nx
 dyt = Ly/ny 
 dzt = Lz/nz 
 dkt = 0.004/nk*2
 dlt = 0.02/nk*2
 dmt = 0.07/nm*2

end subroutine set_grid


subroutine set_initial_conditions
 use main_module   
 use config_module
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: kh,om_lee,om_lee2,m_lee,F_top,freq_fac
 real*8,parameter :: h_rms = 10., nu = 0.8, k_s = 1/10e3

 do k=zs_pe-onx,ze_pe+onx 
  do j=ys_pe-onx,ye_pe+onx   
   do i=xs_pe-onx,xe_pe+onx
    U(i,j,k) =  0.1+0.10*(1+tanh( (zt(k)+750)/150))/(2.)
    Nsqr(i,j,k) = ( 2e-4 + (30e-4)*exp(zt(k)/800.)  )**2 
    !Nsqr(i,j,k) = (30*f0)**2 
   enddo
  enddo
 enddo

 if (my_blk_z == 1       ) then
  do k=1-onx,0
   U(:,:,k) = U(:,:,1)
   V(:,:,k) = V(:,:,1)
  enddo
 endif 
 if (my_blk_z == n_pes_z ) then
  do k=nz+1,nz+onx
   U(:,:,k) = U(:,:,nz)
   V(:,:,k) = V(:,:,nz)
  enddo
 endif 
 

 do m=ms_pe-onx,me_pe+onx 
  do n=ls_pe-onx,le_pe+onx 
   do l=ks_pe-onx,ke_pe+onx 
    do k=zs_pe-onx,ze_pe+onx 
     do j=ys_pe-onx,ye_pe +onx  
      do i=xs_pe-onx,xe_pe+onx 
       E(i,j,k,l,n,m) = 0d0
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo    
 
 allocate(flux_lee(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,&
                   ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx)); 
 flux_lee = 0.                  
 
 if (my_blk_z==1) then
  k=1
  do m=ms_pe-onx+1,me_pe+onx 
   do n=ls_pe-onx,le_pe+onx 
    do l=ks_pe-onx,ke_pe+onx 
     kh = sqrt(kt(l)**2+lt(n)**2)
     do j=ys_pe-onx,ye_pe+onx  
      do i=xs_pe-onx,xe_pe+onx 
        om_lee = -(kt(l)*U(i,j,k) + lt(n)*V(i,j,k))
        om_lee2 = om_lee**2
        
        if ( om_lee2>f0**2 .and. om_lee2<Nsqr(i,j,k) ) then
         m_lee  = -kh*sqrt(Nsqr(i,j,k) - om_lee2)/sqrt(om_lee2 - f0**2)
         F_top = h_rms**2*nu/(pi*k_s**2)/(1.+ (kh/k_s)**2 )**(nu+1.0)
         freq_fac = 4*pi**2*sqrt(Nsqr(i,j,k) - om_lee2)*sqrt(om_lee2 - f0**2)*abs(om_lee)/kh                          
         flux_lee(i,j,l,n,m) = freq_fac*F_top*exp(-(mt(m)-m_lee)**2/(2*dmt(m))**2)/(2*dmt(m)*sqrt(pi))
         
        endif
      enddo
     enddo
    enddo
   enddo
  enddo    
  
  !if (nl==1) flux_lee = flux_lee/dlt(1)
  
 endif

 call write_flux_lee
end subroutine set_initial_conditions










subroutine write_flux_lee
 use main_module
 use config_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,xtdim,ytdim,ktdim,ltdim,mtdim,id,i
 real*8,allocatable :: aloc(:,:,:,:,:)
  
 
 if (my_pe==0) then
  iret = NF_CREATE ('flux_lee.cdf',IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xtdim  = ncddef(ncid, 'xt', nx, iret)
  ytdim  = ncddef(ncid, 'yt', ny, iret)
  ktdim  = ncddef(ncid, 'kt', nk, iret)
  ltdim  = ncddef(ncid, 'lt', nl, iret)  
  mtdim  = ncddef(ncid, 'mt', nm, iret)
  id  = ncvdef (ncid,'flux_lee',NF_DOUBLE,5,(/xtdim,ytdim, ktdim,ltdim, mtdim/),iret)
  call ncendf(ncid, iret)
 endif


 allocate( aloc(nx,ny,nz,nk,nl) ,stat=i); 
 if (i/=0) then
  call halt_stop('ERROR: allocate failed in write_flux_lee (1) ')
 else 
  aloc=0.
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe   ,ks_pe:ke_pe,ls_pe:le_pe) = &
        flux_lee(xs_pe:xe_pe,ys_pe:ye_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'flux_lee',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,i/),(/nx,ny,nk,nl,1/),aloc(:,:,zs_pe,:,:))
   endif
  enddo 
  
  deallocate(aloc)
 endif
 
 if (my_pe==0) call ncclos (ncid, iret) 
end subroutine write_flux_lee



subroutine set_forcing
 use main_module   
 use config_module
 implicit none
 integer :: i,j,k,l,n,m
 
 if (my_blk_z==1) then
  k=1
  do m=ms_pe,me_pe 
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe 
     do j=ys_pe,ye_pe 
      do i=xs_pe,xe_pe 
       E(i,j,k,l,n,m) = E(i,j,k,l,n,m) + dt*flux_lee(i,j,l,n,m)/dzt(k)       
      enddo
     enddo
    enddo
   enddo
  enddo    
 endif
 
 do m=ms_pe,me_pe 
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe 
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe 
      do i=xs_pe,xe_pe 
       !E(i,j,k,l,n,m) = E(i,j,k,l,n,m) - dt*0.01*f0*E(i,j,k,l,n,m) 
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo   
 
 
end subroutine set_forcing


