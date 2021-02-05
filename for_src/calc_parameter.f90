
subroutine calc_parameter
 use main_module
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: kh2
 real*8 :: aloc(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx, &
                ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) 
 
 if (my_pe==0) print*,' calculating in/extrinsic frequencies '
 
 ! calculate ex- and intrinsic frequency from dispersion relation on T grid points
 do m=ms_pe-onx,me_pe+onx
  do n=ls_pe-onx,le_pe+onx
   do l=ks_pe-onx,ke_pe+onx
    do k=zs_pe-onx,ze_pe+onx
     do j=ys_pe-onx,ye_pe+onx  
      do i=xs_pe-onx,xe_pe+onx 
       kh2 = lt(n)**2 + kt(l)**2
       omega_i(i,j,k,l,n,m) = sqrt(  (Nsqr(i,j,k)*kh2+f0**2*mt(m)**2)/(mt(m)**2+kh2) )
       omega_e(i,j,k,l,n,m) = omega_i(i,j,k,l,n,m) + kt(l)*U(i,j,k) + lt(n)*V(i,j,k)
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 
 if (my_pe==0) print*,' calculating velocities in x/k direction'
 
 if (nx>1.and.nk>1) then
 
  ! average omega in x and k direction
  do i=xs_pe-onx,xe_pe+onx-1
   aloc(i,:,:,:,:,:) = 0.5*( omega_e(i,:,:,:,:,:) + omega_e(i+1,:,:,:,:,:) ) 
  enddo
  call border_exchg(aloc)
  do i=ks_pe-onx,ke_pe+onx-1
   aloc(:,:,:,i,:,:) = 0.5*( omega_e(:,:,:,i,:,:) + omega_e(:,:,:,i+1,:,:) ) 
  enddo
  call border_exchg(aloc)
 
  ! xdot = d/dk  omega
  do i=ks_pe-onx+1,ke_pe+onx
   xdot(:,:,:,i,:,:) = ( aloc(:,:,:,i,:,:) - aloc(:,:,:,i-1,:,:) )/dkt(i) 
  enddo
  call border_exchg(xdot)
 
  ! here we set the advection velocity at boundary  to zero
  !if (my_blk_x==n_pes_x)  xdot(nx,:,:,:,:,:) = 0
  !if (my_blk_x==1)        xdot(0, :,:,:,:,:) = 0.

  ! here we set the advection velocity at boundary to finite values to 
  ! transport energy out of the domain. But only allow for vel. out of domain
  if (my_blk_x==n_pes_x) then
     xdot(nx,:,:,:,:,:) = max(0d0,xdot(nx,:,:,:,:,:)) 
     do i=1,onx
      xdot(nx+i,:,:,:,:,:)  = xdot(nx,:,:,:,:,:) 
     enddo
  endif
  if (my_blk_x==1) then
     xdot(0,:,:,:,:,:) = -max(0d0,-xdot(0,:,:,:,:,:))
     do i=1,onx-1
      xdot(0-i,:,:,:,:,:)  = xdot(0,:,:,:,:,:)
     enddo
  endif
  
  ! kdot = -d/dx omega
  do i=xs_pe-onx+1,xe_pe+onx
   kdot(i,:,:,:,:,:) = -( aloc(i,:,:,:,:,:) - aloc(i-1,:,:,:,:,:) )/dxt(i) 
  enddo
  call border_exchg(kdot)

  ! here we set the advection velocity at boundary to finite values to 
  ! transport energy out of the domain. But only allow for vel. out of domain
  if (my_blk_k==n_pes_k) then
     kdot(:,:,:,nk,:,:) = max(0d0,kdot(:,:,:,nk,:,:)) 
     do i=1,onx
      kdot(:,:,:,nk+i,:,:)  = kdot(:,:,:,nk,:,:) 
     enddo
  endif
  if (my_blk_k==1) then
     kdot(:,:,:,0,:,:) = -max(0d0,-kdot(:,:,:,0,:,:))
     do i=1,onx-1
      kdot(:,:,:,0-i,:,:)  = kdot(:,:,:,0,:,:)
     enddo
  endif
  
 else !nx==1 and nk ==1
 
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe
        kh2 = lt(n)**2 + kt(l)**2
        xdot(i,j,k,l,n,m)    = U(i,j,k) + (Nsqr(i,j,k)-f0**2)/omega_i(i,j,k,l,n,m)*mt(m)**2/(mt(m)**2+kh2)**2*kt(l)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif

 if (my_pe==0) print*,' calculating velocities in y/l direction'

 if (ny>1.and.nl>1) then
 
  ! average omega in y and l direction
  do i=ys_pe-onx,ye_pe+onx-1
   aloc(:,i,:,:,:,:) = 0.5*( omega_e(:,i,:,:,:,:) + omega_e(:,i+1,:,:,:,:) ) 
  enddo
  call border_exchg(aloc)
  do i=ls_pe-onx,le_pe+onx-1
   aloc(:,:,:,:,i,:) = 0.5*( omega_e(:,:,:,:,i,:) + omega_e(:,:,:,:,i+1,:) ) 
  enddo
  call border_exchg(aloc)
 
  ! ydot = d/dl  omega
  do i=ls_pe-onx+1,le_pe+onx
   ydot(:,:,:,:,i,:) = ( aloc(:,:,:,:,i,:) - aloc(:,:,:,:,i-1,:) )/dlt(i) 
  enddo
  call border_exchg(ydot)
 
  ! here we set the advection velocity at boundary  to zero 
  !if (my_blk_y==n_pes_y)  ydot(:,ny,:,:,:,:)  = 0
  !if (my_blk_y==1)        ydot(:,0,:,:,:,:) = 0.

  ! here we set the advection velocity at boundary to finite values to 
  ! transport energy out of the domain. But only allow for vel. out of domain
  if (my_blk_y==n_pes_y) then
     ydot(:,ny,:,:,:,:) = max(0d0,ydot(:,ny,:,:,:,:)) 
     do i=1,onx
      ydot(:,ny+i,:,:,:,:)  = ydot(:,ny,:,:,:,:) 
     enddo
  endif
  if (my_blk_y==1) then
     ydot(:,0,:,:,:,:) = -max(0d0,-ydot(:,0,:,:,:,:))
     do i=1,onx-1
      ydot(:,0-i,:,:,:,:)  = ydot(:,0,:,:,:,:)
     enddo
  endif

  ! ldot = -d/dy omega
  do i=ys_pe-onx+1,ye_pe+onx
   ldot(:,i,:,:,:,:) = -( aloc(:,i,:,:,:,:) - aloc(:,i-1,:,:,:,:) )/dyt(i) 
  enddo
  call border_exchg(ldot)

  ! here we set the advection velocity at boundary to finite values to 
  ! transport energy out of the domain. But only allow for vel. out of domain
  if (my_blk_l==n_pes_l) then
     ldot(:,:,:,:,nl,:) = max(0d0,ldot(:,:,:,:,nl,:)) 
     do i=1,onx
      ldot(:,:,:,:,nl+i,:)  = ldot(:,:,:,:,nl,:) 
     enddo
  endif
  if (my_blk_l==1) then
     ldot(:,:,:,:,0,:) = -max(0d0,-ldot(:,:,:,:,0,:))
     do i=1,onx-1
      ldot(:,:,:,:,0-i,:)  = ldot(:,:,:,:,0,:)
     enddo
  endif

 else
 
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe
        kh2 = lt(n)**2 + kt(l)**2
        ydot(i,j,k,l,n,m)    = V(i,j,k) + (Nsqr(i,j,k)-f0**2)/omega_i(i,j,k,l,n,m)*mt(m)**2/(mt(m)**2+kh2)**2*lt(n)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif

 if (my_pe==0) print*,' calculating velocities in z/m direction'

 ! average omega in z and m direction
 do i=zs_pe-onx,ze_pe+onx-1
   aloc(:,:,i,:,:,:) = 0.5*( omega_e(:,:,i,:,:,:) + omega_e(:,:,i+1,:,:,:) ) 
 enddo
 call border_exchg(aloc)
 do i=ms_pe-onx,me_pe+onx-1
   aloc(:,:,:,:,:,i) = 0.5*( omega_e(:,:,:,:,:,i) + omega_e(:,:,:,:,:,i+1) ) 
 enddo
 call border_exchg(aloc)
 
 ! zdot = d/dm omega
 do i=ms_pe-onx+1,me_pe+onx
   zdot(:,:,:,:,:,i) = ( aloc(:,:,:,:,:,i) - aloc(:,:,:,:,:,i-1) )/dmt(i) 
 enddo
 call border_exchg(zdot)
 
 
 if (my_blk_z==n_pes_z) then
   if (.not.enable_upper_reflection) zdot(:,:,nz,:,:,:) = max(0d0,zdot(:,:,nz,:,:,:))
   do i=1,onx
      zdot(:,:,nz+i,:,:,:)  = zdot(:,:,nz,:,:,:)
   enddo
 endif 
  
 if (my_blk_z==1) then
  if (.not.enable_lower_reflection) zdot(:,:,0,:,:,:) = -max(0d0,-zdot(:,:,0,:,:,:))
  do i=1,onx-1
   zdot(:,:,0-i,:,:,:)  = zdot(:,:,0,:,:,:)
  enddo
 endif
 
 ! mdot = -d/dz omega
 do i=zs_pe-onx+1,ze_pe+onx
   mdot(:,:,i,:,:,:) = -( aloc(:,:,i,:,:,:) - aloc(:,:,i-1,:,:,:) )/dzt(i) 
 enddo
 call border_exchg(mdot)
 
 ! here we set the advection velocity at boundary to finite values to 
 ! transport energy out of the domain. But only allow for vel. out of domain
 if (my_blk_m==n_pes_m) then
     mdot(:,:,:,:,:,nm)  = max(0d0,mdot(:,:,:,:,:,nm))
     do i=1,onx
      mdot(:,:,:,:,:,nm+i)  = mdot(:,:,:,:,:,nm)
     enddo
 endif
 if (my_blk_m==1) then
     mdot(:,:,:,:,:,0) = -max(0d0,-mdot(:,:,:,:,:,0))
     do i=1,onx-1
      mdot(:,:,:,:,:,0-i)  = mdot(:,:,:,:,:,0)
     enddo
 endif
 
 ! pre-calculate growth rate for mean-flow interaction
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
         lambda_x(i,j,k,l,n,m) = 1./omega_i(i,j,k,l,n,m)*kt(l)*( &
                   (xdot(i,j,k,l,n,m)-U(i,j,k))*(U(i+1,j,k)-U(i-1,j,k))/(2*dxt(i))   &
                  +(ydot(i,j,k,l,n,m)-V(i,j,k))*(V(i+1,j,k)-V(i-1,j,k))/(2*dxt(i))   )
         lambda_y(i,j,k,l,n,m) = 1./omega_i(i,j,k,l,n,m)*lt(n)*( &
                   (xdot(i,j,k,l,n,m)-U(i,j,k))*(U(i,j+1,k)-U(i,j-1,k))/(2*dyt(j))  &
                  +(ydot(i,j,k,l,n,m)-V(i,j,k))*(V(i,j+1,k)-V(i,j-1,k))/(2*dyt(j))  ) 
         lambda_z(i,j,k,l,n,m) = 1./omega_i(i,j,k,l,n,m)*   &
                              zdot(i,j,k,l,n,m)*( kt(l)*(U(i,j,k+1)-U(i,j,k-1))/(2*dzt(k))   &
                                                 +lt(n)*(V(i,j,k+1)-V(i,j,k-1))/(2*dzt(k)) )                
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 
 
 if (my_pe==0) print*,' done '
end subroutine calc_parameter

