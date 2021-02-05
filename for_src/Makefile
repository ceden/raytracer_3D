
include ../make.inc

all_files = timing_module.o main_module.o main.o parallel_mpi.o restart.o \
            integrate.o superbee_fluxes.o calc_grid.o calc_parameter.o \
            diagnose.o netcdf_io_3D.o netcdf_io_6D.o config.o

all:  ${all_files} 
	$(F90)  ${all_files} $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -o ../bin/model.x


clean :
	 rm -f *.o *.mod 

timing_module.o: timing_module.f90  
	$(F90) $(F90FLAGS) -c timing_module.f90
main_module.o    : main_module.f90 timing_module.o
	 	$(F90) $(F90FLAGS) -c main_module.f90 
main.o         : main.f90 main_module.o  timing_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c main.f90
integrate.o    : integrate.f90 main_module.o   timing_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c integrate.f90    
superbee_fluxes.o : superbee_fluxes.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c superbee_fluxes.f90    
calc_grid.o      : calc_grid.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c calc_grid.f90 
calc_parameter.o : calc_parameter.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c calc_parameter.f90 
restart.o : restart.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c restart.f90         
config.o       : config.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(CDFFLAGS)  -c config.f90
parallel_mpi.o : parallel_mpi.f90 main_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c parallel_mpi.f90
diagnose.o  : diagnose.f90 main_module.o  timing_module.o
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -c diagnose.f90
netcdf_io_3D.o  : netcdf_io_3D.f90 main_module.o  timing_module.o
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -c netcdf_io_3D.f90       
netcdf_io_6D.o  : netcdf_io_6D.f90 main_module.o  timing_module.o
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -c netcdf_io_6D.f90               
