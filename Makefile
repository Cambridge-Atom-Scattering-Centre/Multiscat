#FFLAGS       = -real_size 128
#FFLAGS      = -check_bounds -real_size 128
FFLAGS       = -O3 -mcmodel=large
FORTRAN_DIR  = fortran
#FFLAGS       = -O
#FFLAGS      = -check_bounds
#Used for debugging
#FFLAGS       = -O3 -mcmodel=large  -g -fbacktrace -Wall -fcheck=all

multiscat:		io_loaders.o multiscat.o scatsub.o diagsub.o
			gfortran ${FFLAGS} -o multiscat io_loaders.o multiscat.o scatsub.o diagsub.o

multiscat.o:		${FORTRAN_DIR}/multiscat.f90 ${FORTRAN_DIR}/multiscat.inc io_loaders.o
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o multiscat.o ${FORTRAN_DIR}/multiscat.f90

io_loaders.o:		${FORTRAN_DIR}/io_loaders.f90
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o io_loaders.o ${FORTRAN_DIR}/io_loaders.f90

diagsub.o:		${FORTRAN_DIR}/diagsub.f
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o diagsub.o ${FORTRAN_DIR}/diagsub.f

scatsub.o:		${FORTRAN_DIR}/scatsub.f ${FORTRAN_DIR}/multiscat.inc
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o scatsub.o ${FORTRAN_DIR}/scatsub.f


