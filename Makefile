#FFLAGS       = -real_size 128
#FFLAGS      = -check_bounds -real_size 128
FFLAGS       = -O3 -mcmodel=large
#FFLAGS       = -O
#FFLAGS      = -check_bounds
#Used for debugging
#FFLAGS       = -O3 -mcmodel=large  -g -fbacktrace -Wall -fcheck=all

multiscat:		multiscat.o scatsub.o diagsub.o potsub.o
			gfortran ${FFLAGS} -o multiscat multiscat.o scatsub.o diagsub.o potsub.o

multiscat.o:		multiscat.f90 multiscat.inc
			gfortran -c ${FFLAGS} -o multiscat.o multiscat.f90

diagsub.o:		diagsub.f
			gfortran -c ${FFLAGS} -o diagsub.o diagsub.f

scatsub.o:		scatsub.f multiscat.inc
			gfortran -c ${FFLAGS} -o scatsub.o scatsub.f

potsub.o:		potsub.f90 multiscat.inc
			gfortran -c ${FFLAGS} -o potsub.o potsub.f90

