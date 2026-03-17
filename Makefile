#FFLAGS       = -real_size 128
#FFLAGS      = -check_bounds -real_size 128
FFLAGS       = -O3 -mcmodel=large
FORTRAN_DIR  = fortran
#FFLAGS       = -O
#FFLAGS      = -check_bounds
#Used for debugging
#FFLAGS       = -O3 -mcmodel=large  -g -fbacktrace -Wall -fcheck=all

multiscat:		multiscat.o scatsub.o diagsub.o potsub.o
			gfortran ${FFLAGS} -o multiscat multiscat.o scatsub.o diagsub.o potsub.o

pot2lobatto:		scatsub.o diagsub.o potsub.o ${FORTRAN_DIR}/pot2lobatto.f90 ${FORTRAN_DIR}/multiscat.inc
			gfortran ${FFLAGS} -I${FORTRAN_DIR} -o pot2lobatto ${FORTRAN_DIR}/pot2lobatto.f90 scatsub.o diagsub.o potsub.o

multiscat.o:		${FORTRAN_DIR}/multiscat.f90 ${FORTRAN_DIR}/multiscat.inc
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o multiscat.o ${FORTRAN_DIR}/multiscat.f90

diagsub.o:		${FORTRAN_DIR}/diagsub.f
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o diagsub.o ${FORTRAN_DIR}/diagsub.f

scatsub.o:		${FORTRAN_DIR}/scatsub.f ${FORTRAN_DIR}/multiscat.inc
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o scatsub.o ${FORTRAN_DIR}/scatsub.f

potsub.o:		${FORTRAN_DIR}/potsub.f90 ${FORTRAN_DIR}/multiscat.inc
			gfortran -c ${FFLAGS} -I${FORTRAN_DIR} -o potsub.o ${FORTRAN_DIR}/potsub.f90

