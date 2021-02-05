
subroutine write_restart
 use main_module
 implicit none
 integer :: io=7
 character :: fname*80 

 write(fname,'(a,i9,a)')  'restart_pe=',my_pe,'.dta'
 call replace_space_zero(fname)
 
 if (my_pe==0) print'(a,a)',' writing to restart file ',fname(1:len_trim(fname))
 if (my_pe==0) print'(a,i12)',' at itt = ',itt
 call get_free_iounit(io)    
 open(io,file=fname,form='unformatted',status='unknown')
 write(io,err=10) nx,ny,nz,nk,nl,nm,itt
 write(io,err=10) xs_pe,xe_pe,ys_pe,ye_pe,zs_pe,ze_pe
 write(io,err=10) ks_pe,ke_pe,ls_pe,le_pe,ms_pe,me_pe
 write(io,err=10) E,dE,meanfl
 close(io)
 
 if (my_pe==0) then
  call get_free_iounit(io)
  open(io,file='ritt',form='formatted',status='unknown')
  write(io,*) itt
  close(io)
 endif
 
 
 return
 10 continue
 print'(a)',' Warning: error writing restart file'
end subroutine write_restart




subroutine read_restart
 use main_module
 implicit none
 include 'mpif.h'
 character :: fname*80 
 logical :: file_exists
 integer :: io=7,nx_,ny_,nz_,nk_,nl_,nm_,ierr
 integer :: xs_,xe_,ys_,ye_,zs_,ze_
 integer :: ks_,ke_,ls_,le_,ms_,me_
 
 write(fname,'(a,i9,a)')  'restart_pe=',my_pe,'.dta'
 call replace_space_zero(fname)
 
 inquire ( FILE=fname, EXIST=file_exists )
 if (.not. file_exists) then
   if (my_pe==0) print'(a,a)',' found no restart file ',fname(1:len_trim(fname))
   return
 endif

 if (my_pe==0) print'(a,a)',' reading restart file ',fname(1:len_trim(fname))

 call get_free_iounit(io)    
 open(io,file=fname,form='unformatted',status='old',err=10)
 read(io,err=10) nx_,ny_,nz_,nk_,nl_,nm_,itt

 if (my_pe==0) print*,' itt = ',itt
 !runlen = runlen+itt*dt
 if (my_pe==0) print*,' adjusting run length to = ',runlen

 if (nx/=nx_ .or. ny/=ny_ .or. nz/= nz_ .or. nk/= nk_ .or. nl/= nl_ .or. nm/= nm_) then 
       if (my_pe==0) then
        print*,' read from restart dimensions: ',nx_,ny_,nz_,nk_,nl_,nm_
        print*,' does not match dimensions   : ',nx,ny,nz,nk,nl,nm
       endif
       call MPI_ABORT(mpi_comm_world, 99, IERR)
 endif
  
 read(io,err=10) xs_,xe_,ys_,ye_,zs_,ze_
 read(io,err=10) ks_,ke_,ls_,le_,ms_,me_

 if (xs_/=xs_pe.or.xe_/=xe_pe.or.ys_/=ys_pe.or.ye_/=ye_pe.or.zs_/=zs_pe.or.ze_/=ze_pe.or. &
     ks_/=ks_pe.or.ke_/=ke_pe.or.ls_/=ls_pe.or.le_/=le_pe.or.ms_/=ms_pe.or.me_/=me_pe) then
       if (my_pe==0) then
        print*,' read from restart PE boundaries: ',xs_,xe_,ys_,ye_,zs_,ze_,ks_,ke_,ls_,le_,ms_,me_
        print*,' which does not match           : ',xs_pe,xe_pe,ys_pe,ye_pe,zs_pe,ze_pe,&
                                                    ks_pe,ke_pe,ls_pe,le_pe,ms_pe,me_pe
       endif
       call MPI_ABORT(mpi_comm_world, 99, IERR)
 endif

 read(io,err=10)  E,dE,meanfl
 close(io)
 return
 10 continue
 print'(a)',' Warning: error reading restart file'
 call MPI_ABORT(mpi_comm_world, 99, IERR)
end subroutine read_restart








subroutine get_free_iounit (nu)
!-----------------------------------------------------------------------
!     returns the first free IO unit number in nu
!-----------------------------------------------------------------------
      implicit none
      integer nu,n
      logical in_use
      character (len=80) :: name
      
      do n=7,99
        inquire (n, OPENED=in_use, NAME=name)
        if (.not. in_use) then
          nu = n
          return
         endif
      enddo
end subroutine get_free_iounit




