#!/bin/csh
#set exp = `pwd`
module purge
module load gcc  openmpi/2.0.0-gccsys
#setenv HOSTTYPE thunder
export HOSTTYPE="thunder"
