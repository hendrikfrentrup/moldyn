FC=gfortran
EXE=../test/emd.exe

OBJ_SRC=md.o \
	getControl.o \
	getStructure.o \
	vdwTables.o    \
	setup.o        \
	computeForces.o\
	integrate_nvtGauss.o\
        integrate_nveLeapFrog.o\
	writeConfig.o \
	resetAccum.o  \
	sample.o      \
	evolutionProperties.o \

OBJ_MOD=global.o \
	random.o

$(EXE): $(OBJ_SRC) $(OBJ_MOD)
	$(FC) -o $@ $(OBJ_SRC) $(OBJ_MOD)

%.o: %.f90
	$(FC) -c -O3 $<  # -g -fbounds-check -fbacktrace $<

$(OBJ_SRC):$(OBJ_MOD)

clean:
	rm *.o *.mod  $(EXE)


	

