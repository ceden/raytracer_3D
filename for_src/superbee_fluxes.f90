


subroutine superbee_fluxes
 use main_module
 implicit none
 integer :: i,j,k,l,n,m
 real*8 :: ucfl,Rjp,Rjm,Rj,Cr
 real*8 :: Limiter
 Limiter(Cr)=max(0.D0,max(min(1.D0,2.D0*Cr), min(2.D0,Cr))) 
 
 if (nx>1) then
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe-1,xe_pe   
         uCFL = ABS(xdot(i,j,k,l,n,m)*dt/dxt(i))
         Rjp = E(i+2,j,k,l,n,m) - E(i+1,j,k,l,n,m)
         Rj  = E(i+1,j,k,l,n,m) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i-1,j,k,l,n,m)
         IF (Rj.NE.0.) THEN
          IF (xdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (xdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_x(i,j,k,l,n,m) = xdot(i,j,k,l,n,m)*(E(i+1,j,k,l,n,m)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(xdot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
 if (ny>1) then
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe-1,ye_pe  
      do i=xs_pe,xe_pe 
         uCFL = ABS(ydot(i,j,k,l,n,m)*dt/dyt(j))
         Rjp = E(i,j+2,k,l,n,m) - E(i,j+1,k,l,n,m)
         Rj  = E(i,j+1,k,l,n,m) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i,j-1,k,l,n,m)
         IF (Rj.NE.0.) THEN
          IF (ydot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (ydot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_y(i,j,k,l,n,m) = ydot(i,j,k,l,n,m)*(E(i,j+1,k,l,n,m)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(ydot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
 if (nz>1) then
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe-1,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
         uCFL = ABS(zdot(i,j,k,l,n,m)*dt/dzt(k))
         Rjp = E(i,j,k+2,l,n,m) - E(i,j,k+1,l,n,m)
         Rj  = E(i,j,k+1,l,n,m) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i,j,k-1,l,n,m)
         IF (Rj.NE.0.) THEN
          IF (zdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (zdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_z(i,j,k,l,n,m) = zdot(i,j,k,l,n,m)*(E(i,j,k+1,l,n,m)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(zdot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
       enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
 if (nk>1) then
 do m=ms_pe,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe-1,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
         uCFL = ABS(kdot(i,j,k,l,n,m)*dt/dkt(l))
         Rjp = E(i,j,k,l+2,n,m) - E(i,j,k,l+1,n,m)
         Rj  = E(i,j,k,l+1,n,m) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i,j,k,l-1,n,m)
         IF (Rj.NE.0.) THEN
          IF (kdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (kdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_k(i,j,k,l,n,m) = kdot(i,j,k,l,n,m)*(E(i,j,k,l+1,n,m)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(kdot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
 if (nl>1) then
 do m=ms_pe,me_pe
  do n=ls_pe-1,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe 
         uCFL = ABS(ldot(i,j,k,l,n,m)*dt/dlt(n))
         Rjp = E(i,j,k,l,n+2,m) - E(i,j,k,l,n+1,m)
         Rj  = E(i,j,k,l,n+1,m) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i,j,k,l,n-1,m)
         IF (Rj.NE.0.) THEN
          IF (ldot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (ldot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_l(i,j,k,l,n,m) = ldot(i,j,k,l,n,m)*(E(i,j,k,l,n+1,m)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(ldot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
       enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
 if (nm>1) then
 do m=ms_pe-1,me_pe
  do n=ls_pe,le_pe
   do l=ks_pe,ke_pe
    do k=zs_pe,ze_pe
     do j=ys_pe,ye_pe  
      do i=xs_pe,xe_pe  
         uCFL = ABS(mdot(i,j,k,l,n,m)*dt/dmt(m))
         Rjp = E(i,j,k,l,n,m+2) - E(i,j,k,l,n,m+1)
         Rj  = E(i,j,k,l,n,m+1) - E(i  ,j,k,l,n,m)
         Rjm = E(i  ,j,k,l,n,m) - E(i,j,k,l,n,m-1)
         IF (Rj.NE.0.) THEN
          IF (mdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm/Rj; ELSE; Cr=Rjp/Rj; ENDIF
         ELSE
          IF (mdot(i,j,k,l,n,m).GT.0) THEN; Cr=Rjm*1.E20; ELSE; Cr=Rjp*1.E20; ENDIF
         ENDIF
         Cr=Limiter(Cr)
         flux_m(i,j,k,l,n,m) = mdot(i,j,k,l,n,m)*(E(i,j,k,l,n,m+1)+E(i,j,k,l,n,m))*0.5d0   &
                          -ABS(mdot(i,j,k,l,n,m))*((1.-Cr)+uCFL*Cr)*Rj*0.5d0
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 endif
 
end subroutine superbee_fluxes


