
#-------------------------------------------------------------------------------

PROG = ./run_script
LIB = /home/ecaas/anaconda3/envs/nc-fortran/lib
INCS = /home/ecaas/anaconda3/envs/nc-fortran/include

OBJ = main.o mycmimMod.o fluxMod.o initMod.o \
	testMod.o writeMod.o paramMod.o dispmodule.o \
	shr_kind_mod.o
#-------------------------------------------------------------------------------
FC = gfortran
FFLAGS = -g -fcheck=all -fbacktrace -ffpe-trap=zero,invalid,overflow,underflow -O0 -Wall -ffree-line-length-'none'\
-Wl,-rpath=$(LIB)

EXTRAFFLAGS =  -Wno-unused-dummy-argument -Wno-unused-parameter
FINCLUDES = -I$(INCS)

$(PROG): $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) \
	-L$(LIB) -lnetcdff -lnetcdf -lmfhdf -ldf -lhdf5_hl -lhdf5 -lrt -lpthread -lz -ldl -lm -lcurl \
	-o $(PROG)

mycmimMod.o: mycmimMod.f90 paramMod.o dispmodule.o fluxMod.o initMod.o writeMod.o \
testMod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) mycmimMod.f90

fluxMod.o: fluxMod.f90 paramMod.o dispmodule.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) fluxMod.f90

initMod.o: initMod.f90 paramMod.o dispmodule.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) initMod.f90

testMod.o: testMod.f90 paramMod.o dispmodule.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) testMod.f90

writeMod.o: writeMod.f90 paramMod.o dispmodule.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) writeMod.f90

main.o: main.f90 mycmimMod.o
	 $(FC) $(FFLAGS) -c $(FINCLUDES) main.f90

paramMod.o: paramMod.f90 shr_kind_mod.o
	$(FC) $(FFLAGS) -c $(FINCLUDES) paramMod.f90

dispmodule.o: dispmodule.f90
	$(FC) $(FFLAGS) -c $(FINCLUDES) dispmodule.f90

shr_kind_mod.o: shr_kind_mod.f90
	$(FC) $(FFLAGS) -c $(FINCLUDES) shr_kind_mod.f90
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
strip:
	rm -f *.o

#-------------------------------------------------------------------------------
clean: strip
	rm -f $(PROG)
#-------------------------------------------------------------------------------
