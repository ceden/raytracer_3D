
! model in z,m with GM spectrum as initial condition
! mean flow and N as in lee wave paper
! constant grid in z, grid spacing in m is unequal


subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=4
 nx = 1
 ny = 1
 nz = int(60*fac)
 nk = 1!int(18*fac) ! has to be even or one
 nl = 1!int(20*fac) ! has to be even or one
 nm = int(80*fac) ! has to be even or one
 f0 = 1e-4
 dt = 20.0/fac
 Lz = 2000.0
 Ly = 240d3
 Lx = 240d3
 snapint = 5*50.0
 runlen = snapint*400
 
 k_fixed = 50e-3
 l_fixed = 0.
 
 enable_write_6D_single = .true.
 enable_write_6D_single_himem = .true.
end subroutine set_parameter


subroutine set_grid
 use main_module   
 implicit none
 real,external :: rtbis,mfunc
 real :: m(0:nm),dm(nm),dm0

 dxt = Lx/nx
 dyt = Ly/ny 
 dzt = Lz/nz 
 dkt = pi/Lx
 dlt = pi/Ly

 dm0 =  rtbis(mfunc,1e-6,1e-1,1e-8)
 call mgrid(dm0,m,dm) 
 dm(nm/2+1:nm) = dm(nm/2:1:-1)
 dmt(ms_pe:me_pe) = dm(ms_pe:me_pe)
end subroutine set_grid


subroutine set_initial_conditions
 use main_module   
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: om,Jac,kh2,mstar,funcA,funcB,na,Lm,c1,cstar
 real*8, parameter :: E0=3d-3, nb = 2./pi, jstar = 15.
 real*8 :: z, N_func
 N_func(z) = z*0+30*f0
 !N_func(z) = 10*1d-4*(1+4*exp(z/200.))
 
 do k=zs_pe-onx,ze_pe+onx 
  do j=ys_pe-onx,ye_pe+onx   
   do i=xs_pe-onx,xe_pe+onx
     !U(i,j,k) =  0.4*cos( zt(k)/Lz*2*pi)
     U(i,j,k) =  0.10*(1+tanh( (zt(k)+750)/150))/(2.)
     !U(i,j,k) =  0.5*exp( -(zt(k)+1000)**2/200**2)
     Nsqr(i,j,k) = N_func(zt(k))**2 
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
 
 Lm = 0.
 do i=ms_pe,me_pe
  if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==1.and.my_blk_k==1.and.my_blk_l==1.)  Lm = Lm + dmt(i)
 enddo
 call global_sum(Lm)
 
 c1 = 0.
 do k=zs_pe,ze_pe
    c1 = c1 +  N_func(zt(k))*dzt(k)/pi
 enddo
 call global_sum(c1)
 cstar = c1/jstar
 
 na = 2*gamma(1.)/gamma(0.5)/gamma(0.5)
 do m=ms_pe-onx,me_pe+onx 
  do n=ls_pe-onx,le_pe+onx 
   do l=ks_pe-onx,ke_pe+onx 
    do k=zs_pe-onx,ze_pe+onx 
     do j=ys_pe-onx,ye_pe +onx  
      do i=xs_pe-onx,xe_pe+onx 
        kh2 = kt(l)**2+lt(n)**2
        om = sqrt( (N_func(zt(k))**2*kh2+f0**2*mt(m)**2 )/( kh2+mt(m)**2) )
        Jac = (N_func(zt(k))**2-om**2)**2/( mt(m)**2*om*(N_func(zt(k))**2-f0**2) )
        mstar = sqrt((N_func(zt(k))**2-om**2))/cstar
        funcA = na/(1.+(mt(m)/mstar)**2)
        funcB = nb*abs(f0)/(om*sqrt(om**2-f0**2))
        E(i,j,k,l,n,m) = E0*funcA/mstar*funcB*Jac*N_func(zt(k))/N_func(0d0)/(4*pi)     
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo


end subroutine set_initial_conditions

subroutine set_forcing
end subroutine set_forcing


function zfunc(dz0)
 use main_module 
 implicit none
 real :: zfunc,dz0
 real :: z(0:nz),dz(nz)
 call zgrid(dz0,z,dz)
 zfunc = z(nz)+dz(nz)/2.
end function zfunc

subroutine zgrid(dz0,z,dz)
 use main_module 
 implicit none
 real :: z(0:nz),dz(nz),x,dz0
 integer :: i
 z(0) = -real(Lz)-dz0/2
 do i=1,nz
   x = (i-1.)/(nz-1.) 
   dz(i) = dz0*(1-x**1.5*0.8)
   z(i) = z(i-1) + dz(i)
 enddo
end subroutine zgrid 



function mfunc(dm0)
 use main_module 
 implicit none
 real :: mfunc,dm0
 real :: m(0:nm),dm(nm)
 call mgrid(dm0,m,dm)
 mfunc = m(nm/2)+dm(nm/2)
end function mfunc

subroutine mgrid(dm0,m,dm)
 use main_module 
 implicit none
 real :: m(0:nm),dm(nm),x,dm0
 integer :: i
 m(0) = -1.0
 do i=1,nm/2
   x = (i-1.)/(nm/2-1.) 
   dm(i) = dm0*(1-x**1.5*0.8)
   m(i) = m(i-1) + dm(i)
 enddo
end subroutine mgrid 

function lfunc(dl0)
 use main_module 
 implicit none
 real :: lfunc,dl0
 real :: l(0:nm),dl(nm)
 call lgrid(dl0,l,dl)
 lfunc = l(nl/2)+dl(nl/2)
end function lfunc

subroutine lgrid(dl0,l,dl)
 use main_module 
 implicit none
 real :: l(0:nl),dl(nl),x,dl0
 integer :: i
 l(0) = -0.02
 do i=1,nl/2
   x = (i-1.)/(nl/2-1.) 
   dl(i) = dl0*(1-x**1.5*0.8)
   l(i) = l(i-1) + dl(i)
 enddo
end subroutine lgrid 




   FUNCTION rtbis(func,x1,x2,xacc)
      INTEGER JMAX
      REAL rtbis,x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (JMAX=40)
      INTEGER j
      REAL dx,f,fmid,xmid
      fmid=func(x2)
      f=func(x1)
      if(f*fmid.ge.0.) call halt_stop( 'root must be bracketed in rtbis')
      if(f.lt.0.)then
        rtbis=x1
        dx=x2-x1
      else
        rtbis=x2
        dx=x1-x2
      endif
      do 11 j=1,JMAX
        dx=dx*.5
        xmid=rtbis+dx
        fmid=func(xmid)
        if(fmid.le.0.)rtbis=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.) return
11    continue
      call halt_stop( 'too many bisections in rtbis') 
      END

