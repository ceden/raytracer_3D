


subroutine integrate
 use main_module
 use timing_module
 implicit none
 integer :: i,j,k,l,n,m
 integer               :: xs,xe,ys,ye,ms,me,zs,ze,ls,le,ks,ke
 
 ks=ks_pe; ke=ke_pe; ls=ls_pe; le=le_pe; ms=ms_pe; me=me_pe
 xs=xs_pe; xe=xe_pe; ys=ys_pe; ye=ye_pe; zs=zs_pe; ze=ze_pe

 call tic('bc exchange')
 call border_exchg(E)  
 call toc('bc exchange')
 
 call tic('reflection')
 ! E boundary condition in z-dir: reflection
 call reflection_condition(E) 
 call toc('reflection')
 
 ! E boundary conditions in x dir.: free flux out of the domain
 if (my_blk_x==n_pes_x .and. nx>1) then
     do i=1,onx
      E(nx+i,ys:ye,zs:ze,ks:ke,ls:le,ms:me)  = E(nx,ys:ye,zs:ze,ks:ke,ls:le,ms:me)
     enddo
 endif 
 if (my_blk_x==1.and. nx>1) then
     do i=1,onx
      E(1-i,ys:ye,zs:ze,ks:ke,ls:le,ms:me)  = E(1,ys:ye,zs:ze,ks:ke,ls:le,ms:me)
     enddo
 endif
 
 ! E boundary conditions in y dir.:  free flux out of the domain
 if (my_blk_y==n_pes_y.and. ny>1) then
     do i=1,onx
      E(xs:xe,ny+i,zs:ze,ks:ke,ls:le,ms:me)  = E(xs:xe,ny,zs:ze,ks:ke,ls:le,ms:me)
     enddo
 endif 
 if (my_blk_y==1.and. ny>1) then
     do i=1,onx
      E(xs:xe,1-i,zs:ze,ks:ke,ls:le,ms:me)  = E(xs:xe,1,zs:ze,ks:ke,ls:le,ms:me)
     enddo
 endif

 ! E boundary conditions in z dir.:  reflection or free flux out of the domain
 if (my_blk_z==n_pes_z .and. .not. enable_upper_reflection) then
    do i=1,onx
      E(xs:xe,ys:ye,nz+i,ks:ke,ls:le,ms:me)  = E(xs:xe,ys:ye,nz,ks:ke,ls:le,ms:me)
    enddo
 endif 
 if (my_blk_z==1 .and. .not. enable_lower_reflection) then
     do i=1,onx
      E(xs:xe,ys:ye,1-i,ks:ke,ls:le,ms:me)  = E(xs:xe,ys:ye,1,ks:ke,ls:le,ms:me)
     enddo
 endif


 ! E boundary conditions in k dir.: free flux
 if (my_blk_k==n_pes_k.and. nk>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,nk+i,ls:le,ms:me)  = E(xs:xe,ys:ye,zs:ze,nk,ls:le,ms:me)
     enddo
 endif 
 if (my_blk_k==1.and. nk>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,1-i,ls:le,ms:me)  = E(xs:xe,ys:ye,zs:ze,1,ls:le,ms:me)
     enddo
 endif

 ! E boundary conditions in l dir.: free flux
 if (my_blk_l==n_pes_l.and. nl>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,ks:ke,nl+i,ms:me)  = E(xs:xe,ys:ye,zs:ze,ks:ke,nl,ms:me)
     enddo
 endif 
 if (my_blk_l==1.and. nl>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,ks:ke,1-i,ms:me)  = E(xs:xe,ys:ye,zs:ze,ks:ke,1,ms:me)
     enddo
 endif
 
 ! E boundary conditions in m dir.: free flux
 if (my_blk_m==n_pes_m.and. nm>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,nm+i)  = E(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,nm)
     enddo
 endif 
 if (my_blk_m==1.and. nm>1) then
     do i=1,onx
      E(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,1-i)  = E(xs:xe,ys:ye,zs:ze,ks:ke,ls:le,1)
     enddo
 endif


 call tic('fluxes') 
 ! different advection schemes for fluxes 
 call superbee_fluxes
 call toc('fluxes') 
 
 ! time tendency due to flux divergence and mean flow interaction
 call tic('time tendency') 

 dE(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe) = 0.
 if (nx>1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_x(i,j,k,l,n,m) - flux_x(i-1,j,k,l,n,m))/dxt(i)      
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
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_y(i,j,k,l,n,m) - flux_y(i,j-1,k,l,n,m))/dyt(j)      
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
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_z(i,j,k,l,n,m) - flux_z(i,j,k-1,l,n,m))/dzt(k)      
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
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_k(i,j,k,l,n,m) - flux_k(i,j,k,l-1,n,m))/dkt(l)      
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 
 if (nl>1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_l(i,j,k,l,n,m) - flux_l(i,j,k,l,n-1,m))/dlt(n)      
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 
 if (nm>1) then
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - (flux_m(i,j,k,l,n,m) - flux_m(i,j,k,l,n,m-1))/dmt(m)      
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
 endif
 
 meanfl(xs_pe:xe_pe,ys_pe:ye_pe,zs_pe:ze_pe,ks_pe:ke_pe,ls_pe:le_pe,ms_pe:me_pe)=0.
 
 if (nx>1) then
 
  do m=ms_pe,me_pe
   do n=ls_pe,le_pe
    do l=ks_pe,ke_pe
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe 
         meanfl(i,j,k,l,n,m) = meanfl(i,j,k,l,n,m)  + E(i,j,k,l,n,m)*lambda_x(i,j,k,l,n,m)         
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
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe 
         meanfl(i,j,k,l,n,m) = meanfl(i,j,k,l,n,m)  + E(i,j,k,l,n,m)*lambda_y(i,j,k,l,n,m)         
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
     do k=zs_pe,ze_pe
      do j=ys_pe,ye_pe  
       do i=xs_pe,xe_pe              
        meanfl(i,j,k,l,n,m) = meanfl(i,j,k,l,n,m)  + E(i,j,k,l,n,m)*lambda_z(i,j,k,l,n,m)                    
       enddo
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
       dE(i,j,k,l,n,m) = dE(i,j,k,l,n,m) - meanfl(i,j,k,l,n,m) 
        E(i,j,k,l,n,m) =  E(i,j,k,l,n,m) + dt*dE(i,j,k,l,n,m)                       
      enddo
     enddo
    enddo
   enddo
  enddo
 enddo
 
 call toc('time tendency')  
end subroutine integrate 


