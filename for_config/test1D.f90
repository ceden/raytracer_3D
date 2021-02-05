

subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=2
 nx = 1
 ny = 1
 nz = int(20*fac)
 nk = 1   ! has to be even or one
 nl = 1   ! has to be even or one
 nm = int(60*fac) ! has to be even or one
 f0 = 1e-4
 dt = 5.0/fac
 Lx = 2000.0
 Lz = 2000.0
 Ly = 2000.0
 snapint = 50.0*10
 runlen = snapint*30
 
 k_fixed=1e-3

end subroutine set_parameter


subroutine set_grid
 use main_module   
 implicit none
 dxt = Lx/nx
 dyt = Ly/ny
 dzt = Lz/nz
 dkt = pi/Lx
 dlt = pi/Ly
 dmt = pi/Lz/4.
end subroutine set_grid



subroutine set_initial_conditions
 use main_module   
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: z0,dz,m0,dm

 do k = zs_pe-onx,ze_pe+onx
  do j = ys_pe-onx,ye_pe+onx
   do i = xs_pe-onx,xe_pe+onx
    !U(i,j) =  0.50*exp( -(zt(k)+1000)**2/200**2)
    U(i,j,k) =  0.50*(1+tanh( (zt(k)+750)/150))/(2.)
    Nsqr(i,j,k) =  (30*f0)**2
    !Nsqr(i,j) =  (30*f0)**2*max(0.,1.-exp((zt(k)+100)/50.) )  
    !Nsqr(i,j) =  ( 10*f0*(1+4*exp(zt(k)/200.)) )**2  *max(0.,1.-exp((zt(k)+100)/50.) )  ! 
   enddo
  enddo
 enddo
 z0 = -200.;dz = 50.0
 m0 = -0.02;dm = 4*pi/Lz
 
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe     
        !print*,'pe=',my_pe,' i=',i,' j=',j,' k=',k,l,n,m
       !E(i,j,k,l,tau) = exp(  - (mt(l)-m0)**2/dm**2 - (zt(j)-z0)**2/dz**2 )
       !if (mt(l)<=0) E(i,j,k,l,tau)=1.
       !if (mt(l)>0) E(i,j,k,l,tau)=1.
       if (mt(m)<0) &
        E(i,j,k,l,n,m) = 1.+sin(mt(m)*2*pi/(nm*dmt(ms_pe)/2/6.) ) * sin(zt(k)*2*pi/(Lz/6.))
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo

end subroutine set_initial_conditions



