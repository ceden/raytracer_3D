

subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=1
 nx = 1
 ny = int(10*fac)
 nz = int(20*fac)
 nk = 1
 nl = int(60*fac)   ! has to be even or one
 nm = int(60*fac) ! has to be even or one
 f0 = 1e-4
 dt = 40.0/fac
 Lz = 2000.0
 Ly = 20000.0
 Lx = 20000.0
 snapint = 20*50.0
 runlen = snapint*500

 k_fixed = 1e-3
 
 enable_write_6D_single = .true.
 enable_write_6D_single_himem = .true.
end subroutine set_parameter


subroutine set_grid
 use main_module   
 implicit none
 dxt = Lx/nx
 dyt = Ly/ny
 dzt = Lz/nz
 dkt = pi/Lx
 dlt = 10*pi/Ly
 dmt = pi/Lz/4.
end subroutine set_grid


subroutine set_initial_conditions
 use main_module   
 implicit none
 integer :: i,j,k,l,n,m
 
 do k=zs_pe-onx,ze_pe+onx 
  do j=ys_pe-onx,ye_pe+onx   
   do i=xs_pe-onx,xe_pe+onx
       !U(i,j,k) =  0.50*exp( -(zt(k)+1000)**2/400**2)*exp( -(yt(j)-Ly/2. )**2/2000**2)
       !U(i,j,k) = -0.2*(1+tanh( (zt(k)+750)/150))/(2.)*exp( -(yt(j)-3*Ly/4. )**2/5000**2) &
       !         +0.2*(1+tanh( (zt(k)+750)/150))/(2.)*exp( -(yt(j)-1*Ly/4. )**2/5000**2)
       U(i,j,k) =  0.2*(1+tanh( (zt(k)+750)/150))/(2.)*exp( -(yt(j)-Ly/2. )**2/5000**2)
       Nsqr(i,j,k) =  (30*f0)**2
       !Nsqr(i,j,k) =  (30*f0)**2*max(0.,1.-exp((zt(k)+100)/50.) )  
       !Nsqr(i,j,k) =  ( 10*f0*(1+4*exp(zt(k)/200.)) )**2  *max(0.,1.-exp((zt(k)+100)/50.) )  ! 
   enddo
  enddo
 enddo

 
 

 do m=ms_pe-onx,me_pe+onx 
  do n=ls_pe-onx,le_pe+onx 
   do l=ks_pe-onx,ke_pe+onx 
    do k=zs_pe-onx,ze_pe+onx 
     do j=ys_pe-onx,ye_pe +onx  
      do i=xs_pe-onx,xe_pe+onx 

       
       !if (mt(m)<0) E(i,j,k,l,n,m)=1.+sin(mt(m)*2*pi/(nm*dmt(ms_pe)/2/4.) )*sin(zt(k)*2*pi/(Lz/3.)) &
       !                              *sin(lt(n)*2*pi/(nl*dlt(ls_pe)/2/3.) )*sin(yt(j)*2*pi/(Ly/4.))
       !if (mt(m)<0) E(i,j,k,l,n,m)=1.+sin(lt(n)*2*pi/(nl*dlt(ls_pe)/2/6.) )*sin(zt(k)*2*pi/(Lz/6.))
       !if (mt(m)<0) E(i,j,k,l,n,m)=1.+sin(zt(k)*2*pi/(Lz/6.))*sin(yt(j)*2*pi/(Ly/6.)) 
       if (mt(m)<0) E(i,j,k,l,n,m)=1.
                  
          
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end subroutine set_initial_conditions

