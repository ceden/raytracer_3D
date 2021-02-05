

program main
use main_module
use timing_module
implicit none
integer :: enditt,ierr,n
character*80 :: arg

 !----------------------------------------------------------------------
 ! initialize parallelisation with MPI
 !----------------------------------------------------------------------
 call mpi_init(ierr)
 call my_mpi_init

 !----------------------------------------------------------------------
 ! read in domain decomposition
 !----------------------------------------------------------------------
 n_pes_x = 1; n_pes_y = 1; n_pes_z = 1; n_pes_k = 1; n_pes_l = 1; n_pes_m = 1
 if (n_pes>1) then
   if (iargc() < 6) then
      call halt_stop(' not enough command line input')
   endif
   call getarg(1,arg); read(arg,*) n_pes_x
   call getarg(2,arg); read(arg,*) n_pes_y
   call getarg(3,arg); read(arg,*) n_pes_z
   call getarg(4,arg); read(arg,*) n_pes_k
   call getarg(5,arg); read(arg,*) n_pes_l
   call getarg(6,arg); read(arg,*) n_pes_m
 else
   print*,' using a single cpu'
 endif
 if (my_pe==0.and.n_pes>1) print'(/a,i4,a,i4,a,i4,a,i4,a,i4,a,i4,a)', &
                    'using ',n_pes_x,' x ',n_pes_y ,' x ',n_pes_z,' x ',n_pes_k,' x ',n_pes_l,' x ',n_pes_m,' PEs'

 call set_parameter
 call pe_decomposition()
 if (my_pe==0) print*,' allocating memory'
 call allocate_main_module()
 
 if (my_pe==0) print*,' setting up grid'
 call set_grid()
 call calc_grid
 
 if (my_pe==0) print*,' setting initial conditions'
 call set_initial_conditions
 call calc_parameter
 
 call read_restart
 
 enditt = itt+int(runlen/dt)
 if (my_pe==0) then
     print*,' '
     print'(a,e8.2,a)',' setup integration for ',runlen,' s'
     print'(a,i10,a,i10)',' from time step ',itt,' to ',enditt
     print'(a,e8.2,a)',' snapshot intervall is ',snapint,' s'
     print'(a,i8,a)',' this is any ',int(snapint/dt),' time steps'
     print*,' '
 endif

 call init_diag() 
 
 ! main loop 
 do while (itt < enditt) 
      call tic('main loop')
      
      
      call tic('integrate')
      call set_forcing
      call integrate
      call toc('integrate')
      
      call tic('diag')
      if ( mod(itt,int(snapint/dt))  == 0 .or. itt == 0)  call diagnose
      call toc('diag')
    
      itt=itt+1
      call toc('main loop')
 enddo
 
 call write_restart
 
 !--------------------------------------------------------------
 !     show timing results here
 !--------------------------------------------------------------
 do n = 0,0!n_pes
       call fortran_barrier
       if (my_pe == n) then
        print'(/,a,i4)','Timing summary for PE #',my_pe 
        print'(a,f12.1,a)',' costs for measuring      = ',timing_secs('tictoc'),' s'
        
        print'(a,f12.1,a)',' main loop time summary   = ',timing_secs('main loop') ,' s'
        print'(a,f12.1,a)','        diagnostics       = ',timing_secs('diag') ,' s'
        print'(a,f12.1,a)','            writing       = ',timing_secs('diag_write') ,' s'
        print'(a,f12.1,a)','        integrate         = ',timing_secs('integrate') ,' s'
        print'(a,f12.1,a)','            fluxes        = ',timing_secs('fluxes') ,' s'
        print'(a,f12.1,a)','            bc exchange   = ',timing_secs('bc exchange') ,' s'
        print'(a,f12.1,a)','            reflection    = ',timing_secs('reflection') ,' s'
        print'(a,f12.1,a)','            time tendency = ',timing_secs('time tendency') ,' s'
       endif
 enddo
     
 !----------------------------------------------------------------------
 ! cancel parallelisation and quit
 !----------------------------------------------------------------------
 if (my_pe==0) print'(/a/)','cancelling MPI service'
 call mpi_finalize(ierr)
 
end program main




