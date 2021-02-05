
! Model with z,m dimension and flux with GM-spectrum at bottom
! mean flow and N as in critical layer paper
! constant grid spacing in all directions

module config_module
 implicit none
 real*8, allocatable :: E_gm(:,:,:,:,:)
end module config_module

subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=1
 nx = 1
 ny = 1
 nz = int(160*fac)
 nk = 1 !60
 nl = 1
 nm = int(300*fac) ! has to be even
 f0 = 1e-4
 dt = 30.
 Ly = 100d3
 Lx = 100d3
 Lz = 2000.0
 snapint = 400*50.0
 runlen = snapint*500 

 k_fixed = 25.5e-3
 l_fixed = 0.
 enable_upper_reflection = .false.
 enable_lower_reflection = .false.
 
 enable_write_6D_single = .true.
 enable_write_6D_single_himem = .true.
 
end subroutine set_parameter


subroutine set_grid
 use main_module   
 implicit none
 dxt = Lx/nx
 dyt = Ly/ny 
 dzt = Lz/nz
 dmt = 2*pi/Lz/2.
 dkt = 0.5e-3 
 dlt = 0.5e-3 
end subroutine set_grid


subroutine set_initial_conditions
 use main_module 
 use config_module
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: om,Jac,mstar,funcA,funcB,na,cstar,c1,kh
 real*8, parameter :: E0=3d-3, nb = 2./pi
 integer, parameter :: jstar=10

 do k=zs_pe-onx,ze_pe+onx 
  do j=ys_pe-onx,ye_pe+onx   
   do i=xs_pe-onx,xe_pe+onx
    Nsqr(i,j,k) = (30e-4)**2
    U(i,j,k) =  0.2*exp( -(zt(k)+1000)**2/200**2)
   enddo
  enddo 
 enddo
 
 if (my_blk_z == 1       ) then
  do k=1-onx,0
   U(:,:,k) = U(:,:,1)
  enddo
 endif 
 if (my_blk_z == n_pes_z ) then
  do k=nz+1,nz+onx
   U(:,:,k) = U(:,:,nz)
  enddo
 endif 
 
 allocate(E_gm(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx)); 
 E_gm = 0. 
 
 na = 2*gamma(1.)/gamma(0.5)/gamma(0.5)
 c1 = 0.
 do k=zs_pe,ze_pe
    c1 = c1 +  sqrt(Nsqr(xs_pe,ys_pe,k))*dzt(k)/pi
 enddo
 call global_sum(c1)
 cstar = c1/jstar
 
 if (my_blk_z==1) then
 
  do m=ms_pe-onx,me_pe+onx 
   do n=ls_pe-onx,le_pe+onx 
    do l=ks_pe-onx,ke_pe+onx 
     k = 1
     do j=ys_pe-onx,ye_pe+onx  
      do i=xs_pe-onx,xe_pe+onx    
        kh = sqrt( kt(l)**2 +lt(n)**2)
        om = sqrt( (Nsqr(i,j,k)*kh**2+f0**2*mt(m)**2 )/( kh**2+mt(m)**2) )
        Jac = (Nsqr(i,j,k)-om**2)**2/( mt(m)**2*om*(Nsqr(i,j,k)-f0**2) )
        mstar = sqrt((Nsqr(i,j,k)-om**2))/cstar
        funcA = na/(1.+(mt(m)/mstar)**2)
        funcB = nb*abs(f0)/(om*sqrt(om**2-f0**2))
        E_gm(i,j,l,n,m) = E0*funcA/mstar*funcB*Jac!/(4*pi)
      enddo
     enddo
    enddo
   enddo
  enddo 
 endif
  call write_e_gm
 
end subroutine set_initial_conditions


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
       E(i,j,k,l,n,m) = E(i,j,k,l,n,m) - dt*max(0d0,zdot(i,j,k,l,n,m))* &
                               (E(i,j,k,l,n,m)-E_gm(i,j,l,n,m))/dzt(k)    
      enddo
     enddo
    enddo
   enddo
  enddo    
 endif

end subroutine set_forcing






subroutine write_E_gm
 use main_module
 use config_module
 implicit none
 include "netcdf.inc"
 integer :: ncid,iret,xtdim,ytdim,ktdim,ltdim,mtdim,id,i
 real*8,allocatable :: aloc(:,:,:,:,:)
  
 
 if (my_pe==0) then
  iret = NF_CREATE ('E_gm.cdf',IOR(NF_CLOBBER,NF_64BIT_OFFSET),ncid)
  iret=nf_set_fill(ncid, NF_NOFILL, iret)
  xtdim  = ncddef(ncid, 'xt', nx, iret)
  ytdim  = ncddef(ncid, 'yt', ny, iret)
  ktdim  = ncddef(ncid, 'kt', nk, iret)
  ltdim  = ncddef(ncid, 'lt', nl, iret)  
  mtdim  = ncddef(ncid, 'mt', nm, iret)
  id  = ncvdef (ncid,'E_gm',NF_DOUBLE,5,(/xtdim,ytdim, ktdim,ltdim, mtdim/),iret)
  call ncendf(ncid, iret)
 endif


 allocate( aloc(nx,ny,nz,nk,nl) ,stat=i); 
 if (i/=0) then
  call halt_stop('ERROR: allocate failed in write_E_gm (1) ')
 else 
  aloc=0.
 
  do i=1,nm
   if (i>=ms_pe.and.i<=me_pe) &
    aloc(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe   ,ks_pe:ke_pe,ls_pe:le_pe) = &
        E_gm(xs_pe:xe_pe,ys_pe:ye_pe,ks_pe:ke_pe,ls_pe:le_pe,i)
   call pe0_recv_arr_mdim(aloc,i)
   if (my_pe==0) then
    iret = nf_inq_varid(ncid,'E_gm',id)
    iret= nf_put_vara_double(ncid,id,(/1,1,1,1,i/),(/nx,ny,nk,nl,1/),aloc(:,:,zs_pe,:,:))
   endif
  enddo 
  
  deallocate(aloc)
 endif
 
 if (my_pe==0) call ncclos (ncid, iret) 
end subroutine write_E_gm

