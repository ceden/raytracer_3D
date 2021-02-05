
! model with y,z,k,l,m dimensions for Brava eddy
! U as for Brava eddy and first baroclinic mode 
! N is idealized
! z,k,l,m grid with unequal spacing, y constant grid spacing

subroutine set_parameter
 use main_module   
 implicit none
 real*8 :: fac=1
 nx = 1
 ny = int(10*fac)
 nz = int(20*fac)
 nk = int(20*fac) ! has to be even or one
 nl = int(20*fac) ! has to be even or one
 nm = int(20*fac) ! has to be even or one
 f0 = 1e-4
 dt = 20.0/fac
 Lz = 2000.0
 Ly = 240d3
 Lx = 240d3
 snapint = 20*50.0
 runlen = snapint*500
 
 enable_write_6D_single = .true.
 enable_write_6D_single_himem = .true.
end subroutine set_parameter




subroutine set_grid
 use main_module   
 implicit none
 real,external :: zfunc,rtbis,mfunc,lfunc
 real :: dz0,z(0:nz),dz(nz)
 real :: m(0:nm),dm(nm),dm0
 real :: l(0:nl),dl(nm),dl0

 dxt = Lx/nx
 dyt = Ly/ny 
 
 dz0 =  rtbis(zfunc,10.0,1000.0,1e-8)
 call zgrid(dz0,z,dz)
 dzt(zs_pe:ze_pe) = dz(zs_pe:ze_pe) 
 
 dl0 =  rtbis(lfunc,1e-6,1e-1,1e-8)
 call lgrid(dl0,l,dl) 
 dl(nl/2+1:nl) = dl(nl/2:1:-1)
 dlt(ls_pe:le_pe) = dl(ls_pe:le_pe)
 dkt(ks_pe:ke_pe) = dl(ks_pe:ke_pe)
 
 dm0 =  rtbis(mfunc,1e-6,1e-1,1e-8)
 call mgrid(dm0,m,dm) 
 dm(nm/2+1:nm) = dm(nm/2:1:-1)
 dmt(ms_pe:me_pe) = dm(ms_pe:me_pe)
 
end subroutine set_grid


subroutine set_initial_conditions
 use main_module   
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: c1
 real*8 :: r,om,Jac,kh2,mstar,funcA,funcB,na,cstar 
 real*8, parameter :: E0=3d-3, nb = 2./pi, jstar = 10.
 real*8, dimension(1-onx:nz+onx) :: dzt_gl,zt_gl,wn
 include "mpif.h"
 integer :: tag = 99, ierr, sender, ind(3)
 integer,dimension(MPI_STATUS_SIZE)  :: Status
 real*8 :: z, N_func
 N_func(z) = 10*1d-4*(1+4*exp(z/200.))
 
 
 
 zt_gl(zs_pe:ze_pe) = zt(zs_pe:ze_pe)
 dzt_gl(zs_pe:ze_pe) = dzt(zs_pe:ze_pe)
 do n=1,n_pes_z-1
   if (my_blk_x==1.and.my_blk_y==1.and.my_blk_z==n+1.and.my_blk_k==1.and.my_blk_l==1.and.my_blk_m==1) then
    call mpi_send((/zs_pe,ze_pe,z_blk/),3,mpi_integer,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(zt(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
    call mpi_send(dzt(zs_pe:ze_pe),z_blk,mpi_real8,0,tag,MPI_COMM_WORLD,ierr)
   else if (my_pe==0) then
    sender = n*n_pes_x*n_pes_y
    call mpi_recv(ind,3,mpi_integer,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(zt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
    call mpi_recv(dzt_gl(ind(1):ind(2)),ind(3),mpi_real8,sender,tag,MPI_COMM_WORLD,Status,ierr)
   endif
 enddo
 if (n_pes_z>1) call pe0_bcast(zt_gl,nz+2*onx)
 if (n_pes_z>1) call pe0_bcast(dzt_gl,nz+2*onx)
 
 c1 = 0.
 do k=1,nz
    c1 = c1 +  N_func(zt_gl(k))*dzt_gl(k)/pi
 enddo
 wn=0
 do k=1,nz
    wn(k) = wn(k-1) +  N_func(zt_gl(k))*dzt_gl(k)/c1
 enddo
 do k=1,nz
  wn(k) = cos(wn(k))*sqrt(2*N_func(zt_gl(k)) /(c1*pi) )
 enddo
 wn(1-onx:0) = wn(1)
 wn(nz+1:nz+onx) = wn(nz)
 !where(wn >0) wn=0
 
 do k=zs_pe-onx,ze_pe+onx 
  do j=ys_pe-onx,ye_pe+onx   
   do i=xs_pe-onx,xe_pe+onx
     r = abs(yt(j)-Ly/2) 
     if ( r < 45d3 ) then
       U(i,j,k) =  -0.2*r/45d3*wn(k)*sqrt(Lz/2)
       !U(i,j,k) =  0.4*r/45d3*exp(zt(k)/250.) 
     else
       U(i,j,k) =  -0.2*exp( -(r-45d3 )/33d3 )*wn(k)*sqrt(Lz/2)
       !U(i,j,k) =  0.4*exp( -(r-45d3 )/33d3 )*exp(zt(k)/250.) 
     endif
     if (yt(j)>Ly/2) U(i,j,k) = -U(i,j,k) 
     Nsqr(i,j,k) =   N_func(zt(k))**2
   enddo
  enddo
 enddo

 na = 2*gamma(1.)/gamma(0.5)/gamma(0.5)
 cstar = c1/jstar
 do m=ms_pe-onx,me_pe+onx 
  do n=ls_pe-onx,le_pe+onx 
   do l=ks_pe-onx,ke_pe+onx 
    do k=zs_pe-onx,ze_pe+onx 
     do j=ys_pe-onx,ye_pe +onx  
      do i=xs_pe-onx,xe_pe+onx 
        ! set GM spectrum
        kh2 = kt(l)**2+lt(n)**2
        om = sqrt( (Nsqr(i,j,k)*kh2+f0**2*mt(m)**2 )/( kh2+mt(m)**2) )
        Jac = (Nsqr(i,j,k)-om**2)**2/( mt(m)**2*om*(Nsqr(i,j,k)-f0**2) )
        mstar = sqrt((Nsqr(i,j,k)-om**2))/cstar
        funcA = na/(1.+(mt(m)/mstar)**2)
        funcB = nb*abs(f0)/(om*sqrt(om**2-f0**2))
        E(i,j,k,l,n,m) = E0*funcA/mstar*funcB*Jac*N_func(zt(k))/N_func(0d0)/(4*pi)
        
        ! zero out waves with om<f 
        !r = (om+kt(l)*U(i,j,k))/abs(f0) 
        !E(i,j,k,l,n,m) = E(i,j,k,l,n,m) *(1+tanh((r-1.01)*20.))/(2.)  
        !if ( r <=1.01 ) E(i,j,k,l,n,m) = 0d0
       
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
 m(0) = -0.1
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





