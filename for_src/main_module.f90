

module main_module
 implicit none

!---------------------------------------------------------------------------------
! some fixed parameter
!---------------------------------------------------------------------------------
  real*8, parameter :: version = 0.1
  real*8, parameter :: pi       = 3.14159265358979323846264338327950588
!---------------------------------------------------------------------------------
! model parameter
!---------------------------------------------------------------------------------
  integer :: nx ! number of grid points in y-direction
  integer :: ny ! number of grid points in y-direction
  integer :: nz ! number of grid points in z-direction
  integer :: nk ! number of grid points in k-direction (must be even or zero)
  integer :: nm ! number of grid points in m-direction (must be even or zero)
  integer :: nl ! number of grid points in l-direction (must be even or zero)
  real*8 :: dt ! time step in seconds
  real*8 :: Lx ! hor. extent in m
  real*8 :: Ly ! hor. extent in m
  real*8 :: Lz ! water depth in m
  real*8 :: f0 = 1d-4 ! Coriolis parameter in 1/s
  real*8 :: k_fixed=0d0,l_fixed=0d0,m_fixed=0d0 ! wavenumbers in case that nk/nl/nm=1
  real*8 :: sign_of_omega_i = +1d0
!---------------------------------------------------------------------------------
! switches for general model setup
!---------------------------------------------------------------------------------  
  logical :: enable_upper_reflection = .true.
  logical :: enable_lower_reflection = .true.
!---------------------------------------------------------------------------------
! time stepping parameter
!---------------------------------------------------------------------------------
  integer :: itt = 0                       ! current time step number
  real*8 :: runlen = 0.                    ! length of integration in s
  real*8 :: snapint = 0.                   ! snapshot interval in s
!---------------------------------------------------------------------------------
! model fields
!---------------------------------------------------------------------------------
  real*8, allocatable :: E(:,:,:,:,:,:)  ! wave energy
  real*8, allocatable :: dE(:,:,:,:,:,:) ! rate of change in wave energy
  real*8, allocatable :: xdot(:,:,:,:,:,:) ! zonal group velocity in m/s
  real*8, allocatable :: ydot(:,:,:,:,:,:) ! meridional group velocity in m/s
  real*8, allocatable :: zdot(:,:,:,:,:,:) ! vertical group velocity in m/s
  real*8, allocatable :: kdot(:,:,:,:,:,:) ! refraction parameter in 1/ms
  real*8, allocatable :: ldot(:,:,:,:,:,:) ! refraction parameter in 1/ms
  real*8, allocatable :: mdot(:,:,:,:,:,:) ! refraction parameter in 1/ms
  real*8, allocatable :: omega_i(:,:,:,:,:,:)   ! intrinsic frequency in 1/s 
  real*8, allocatable :: omega_e(:,:,:,:,:,:)   ! extrinsic frequency in 1/s
  real*8, allocatable :: flux_x(:,:,:,:,:,:)    ! energy flux in x dir
  real*8, allocatable :: flux_y(:,:,:,:,:,:)    ! energy flux in y dir
  real*8, allocatable :: flux_z(:,:,:,:,:,:)    ! energy flux in z dir
  real*8, allocatable :: flux_k(:,:,:,:,:,:)    ! energy flux in k dir
  real*8, allocatable :: flux_l(:,:,:,:,:,:)    ! energy flux in l dir
  real*8, allocatable :: flux_m(:,:,:,:,:,:)    ! energy flux in m dir
  real*8, allocatable :: lambda_x(:,:,:,:,:,:)  ! mean-flow interaction growth rate
  real*8, allocatable :: lambda_y(:,:,:,:,:,:)  ! mean-flow interaction growth rate
  real*8, allocatable :: lambda_z(:,:,:,:,:,:)  ! mean-flow interaction growth rate
  
  real*8, allocatable :: meanfl(:,:,:,:,:,:)  ! energy change due to mean-flow interaction
  real*8, allocatable :: Nsqr(:,:,:)          ! stability frequency in 1/s
  real*8, allocatable :: U(:,:,:)             ! mean flow  in m/s
  real*8, allocatable :: V(:,:,:)             ! mean flow  in m/s
  
  real*8, allocatable :: xt(:),xu(:)    ! x position of T and U points in m 
  real*8, allocatable :: dxt(:),dxu(:)  ! x grid thickness for T and U points in m 
  real*8, allocatable :: yt(:),yu(:)    ! y position of T and U points in m 
  real*8, allocatable :: dyt(:),dyu(:)  ! y grid thickness for T and U points in m 
  real*8, allocatable :: zt(:),zu(:)    ! vertical position of T and U points in m 
  real*8, allocatable :: dzt(:),dzu(:)  ! vertical grid thickness for T and U points in m 
  real*8, allocatable :: mt(:),mu(:)    ! m position of T and U points in 1/m 
  real*8, allocatable :: dmt(:),dmu(:)  ! m grid thickness for T and U points in 1/m  
  real*8, allocatable :: lt(:),lu(:)    ! same for wavenumber l 
  real*8, allocatable :: dlt(:),dlu(:)  
  real*8, allocatable :: kt(:),ku(:)    ! same for wavenumber k
  real*8, allocatable :: dkt(:),dku(:)  
!---------------------------------------------------------------------------------
!     Parallel domain setup
!---------------------------------------------------------------------------------
  integer :: n_pes     ! total number of processors
  integer :: my_pe     ! index of this processor from 0 to n_pes-1
  integer :: n_pes_x   ! total number of processors in x direction
  integer :: n_pes_y   ! total number of processors in y direction
  integer :: n_pes_z   ! total number of processors in z direction
  integer :: n_pes_k   ! total number of processors in k direction
  integer :: n_pes_l   ! total number of processors in l direction
  integer :: n_pes_m   ! total number of processors in m direction
  integer :: my_blk_x  ! index of this processor in x direction from 1 to n_pes_x
  integer :: my_blk_y  ! index of this processor in y direction from 1 to n_pes_y
  integer :: my_blk_z  ! index of this processor in x direction from 1 to n_pes_z
  integer :: my_blk_m  ! index of this processor in m direction from 1 to n_pes_m
  integer :: my_blk_l  ! same for l
  integer :: my_blk_k  ! same for k
  integer :: x_blk     ! grid points of domain decompostion in x direction 
  integer :: y_blk     ! grid points of domain decompostion in y direction 
  integer :: z_blk     ! grid points of domain decompostion in z direction 
  integer :: m_blk     ! grid points of domain decompostion in m direction
  integer :: l_blk     ! grid points of domain decompostion in l direction
  integer :: k_blk     ! grid points of domain decompostion in k direction
  
  integer :: xs_pe     ! start index of grid points in x direction of this processor
  integer :: xe_pe     ! end index of grid points in x direction of this processor
  integer :: ys_pe     ! start index of grid points in y direction of this processor
  integer :: ye_pe     ! end index of grid points in y direction of this processor
  integer :: zs_pe     ! start index of grid points in z direction of this processor
  integer :: ze_pe     ! end index of grid points in z direction of this processor
  integer :: ms_pe     ! start index of grid points in m direction of this processor
  integer :: me_pe     ! end index of grid points in m direction of this processor
  integer :: ls_pe     ! same for l
  integer :: le_pe
  integer :: ks_pe     ! same for k
  integer :: ke_pe
  integer :: onx=2     ! number of overlapping points in all directions
  integer :: mpi_comm_phys ! mpi communication for all PEs of same x,y,z block
  
  
  !---------------------------------------------------------------------------------
  ! switches for diagnostic output
  !---------------------------------------------------------------------------------
  logical :: enable_show_grid_details     = .false. ! print out grid details
  logical :: enable_write_6D_single       = .true. ! write 6D fields to a single file
                                                   ! alternative: each PE writes to its own file
  logical :: enable_write_6D_single_himem = .true. ! high-speed, but high memory option
  character*160 :: netcdf_6D_file                  ! netcdf output file name for 6D fields
  
end module main_module



subroutine  allocate_main_module
 use main_module
 implicit none 
 integer :: status
 allocate( E(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
             ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10            
 allocate(dE(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
             ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx),stat=status ); 
 if (status/=0) goto 10         
 allocate(xdot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10              
 allocate(ydot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(zdot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(kdot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(ldot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status);
 if (status/=0) goto 10
 allocate(mdot(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
               ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10     
 allocate(omega_i(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                  ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(omega_e(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                  ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status);
 if (status/=0) goto 10             
 allocate(flux_x(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(flux_y(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(flux_z(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(flux_k(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(flux_l(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(flux_m(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate(lambda_x(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                   ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status);   
 if (status/=0) goto 10
 allocate(lambda_y(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                   ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status);                     
 if (status/=0) goto 10
 allocate(lambda_z(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                   ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status);                
 if (status/=0) goto 10
 allocate(meanfl(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx,&
                 ks_pe-onx:ke_pe+onx,ls_pe-onx:le_pe+onx,ms_pe-onx:me_pe+onx) ,stat=status); 
 if (status/=0) goto 10             
 allocate( Nsqr(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx),stat=status ); 
 if (status/=0) goto 10
 allocate( U(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx) ,stat=status);
 if (status/=0) goto 10
 allocate( V(xs_pe-onx:xe_pe+onx,ys_pe-onx:ye_pe+onx,zs_pe-onx:ze_pe+onx) ,stat=status);
 if (status/=0) goto 10
 allocate( xt(xs_pe-onx:xe_pe+onx) , xu(xs_pe-onx:xe_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate( dxt(xs_pe-onx:xe_pe+onx) , dxu(xs_pe-onx:xe_pe+onx) ,stat=status);
 if (status/=0) goto 10
 allocate( yt(ys_pe-onx:ye_pe+onx) , yu(ys_pe-onx:ye_pe+onx) ,stat=status); 
 if (status/=0) goto 10
 allocate( dyt(ys_pe-onx:ye_pe+onx) , dyu(ys_pe-onx:ye_pe+onx) ,stat=status)
 if (status/=0) goto 10
 allocate( zt(zs_pe-onx:ze_pe+onx)  , zu(zs_pe-onx:ze_pe+onx) ,stat=status)
 if (status/=0) goto 10
 allocate( dzt(zs_pe-onx:ze_pe+onx) , dzu(zs_pe-onx:ze_pe+onx),stat=status )
 if (status/=0) goto 10
 allocate( kt(ks_pe-onx:ke_pe+onx) , ku(ks_pe-onx:ke_pe+onx) ,stat=status)
 if (status/=0) goto 10
 allocate( dkt(ks_pe-onx:ke_pe+onx) , dku(ks_pe-onx:ke_pe+onx) ,stat=status)
 if (status/=0) goto 10
 allocate( lt(ls_pe-onx:le_pe+onx) , lu(ls_pe-onx:le_pe+onx) ,stat=status)
 if (status/=0) goto 10
 allocate( dlt(ls_pe-onx:le_pe+onx) , dlu(ls_pe-onx:le_pe+onx),stat=status )
 if (status/=0) goto 10
 allocate( mt(ms_pe-onx:me_pe+onx) , mu(ms_pe-onx:me_pe+onx),stat=status )
 if (status/=0) goto 10
 allocate( dmt(ms_pe-onx:me_pe+onx) , dmu(ms_pe-onx:me_pe+onx),stat=status )
 if (status/=0) goto 10

 E=0;dE=0;xdot=0; ydot=0; zdot=0; kdot=0;  ldot=0; mdot=0;    omega_i=0; omega_e=0; 
 flux_x=0;flux_y=0;flux_z=0;flux_k=0;flux_l=0;flux_m=0;
 lambda_x=0; lambda_y=0; lambda_z=0; meanfl=0;Nsqr=0; U=0; V=0;
 xt = 0;xu = 0;dxt = 0; dxu = 0
 yt = 0;yu = 0; dyt = 0; dyu = 0; zt=0; zu=0; 
 dzt=0; dzu=0; kt = 0; ku = 0; dkt = 0; dku = 0
 lt = 0; lu = 0; dlt = 0; dlu = 0; mt=0; mu=0; dmt=0; dmu=0
  return
  
10 continue
  call halt_stop('ERROR: not enough memory to allocate model variables ')
  
end subroutine  allocate_main_module
