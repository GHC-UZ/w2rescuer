#FC=ifort
#FFLAGS= # -openmp -prefetch -unroll-aggressive -axS -ipo -ip -ipo-jobs4 -mtune=pentium4  -openmp-report2  -fast # -fp-model precise -mp -mp1 -warn all
FC=gfortran
FFLAGS= -O2 -Wno-unused-variable  
LFLAGS= -L./  -l blas
#LFLAGS=  
EXE=my_1dsolver


SRC= main.f90 euler.f90  bound.f90 qinit.f90  entropyl.f90 flux_dot.f90 flux_roe.f90 limiters.f90  grid_x.f90 initial_bed.f90 initial_cond.f90 postpro.f90  friction.f90  gauges.f90  flux_hll.f90 flux_lw.f90 matmul2.f90

#
#
# No need to edit below this line
#

.SUFFIXES:
.SUFFIXES: .f90 .o

OBJ=	$(SRC:.f90=.o)
MOD=	$(SRC:.f90=.mod)

.f90.o:
	$(FC) $(FFLAGS) $< -c

all:	$(EXE)

$(EXE):	$(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LFLAGS)

$(OBJ):	$(MF)

clean:
	#@echo $(OBJ) $(EXE) $(MOD) 
	rm -f $(OBJ) $(EXE) $(MOD) 
