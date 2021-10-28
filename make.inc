
#########################################################################
# platform dependent options
#########################################################################
 
# for mistral at DKRZ
#cdf = /sw/rhel6-x64/netcdf/netcdf_fortran-4.4.2-gcc48/bin/nf-config
#F90 = mpif90
#CDFFLAGS = -I`$(cdf) --includedir`  `$(cdf) --flibs`
#MPIFLAGS =
#F90FLAGS =   -fconvert=big-endian -O3 -Wall

# for thunder maybe outdated
#cdf = /sw/wheezy-x64/netcdf_fortran-4.4.3-static-gccsys/bin/nf-config
#F90 = mpif90
#CDFFLAGS = -I`$(cdf) --includedir`  `$(cdf) --flibs`
#MPIFLAGS =
#F90FLAGS =   -fconvert=big-endian -O3 -Wall

# for ubuntu
#F90 = gfortran 
#CDFFLAGS = -I/usr/include -lnetcdf -lnetcdff 
#MPIFLAGS = -I/usr/include/mpi -lmpi -lmpi_f77 
#F90FLAGS = -Wall -O3  -fbounds-check -finit-real=nan -fallow-argument-mismatch

# for OSX
F90 =  mpif90
CDFFLAGS = -I/usr/local/include -lnetcdff -lnetcdf -lnetcdf
MPIFLAGS =  
F90FLAGS = -Wall -O3 -fbounds-check -finit-real=nan -fallow-argument-mismatch


