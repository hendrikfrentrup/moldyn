#FC=ifort
FC=gfortran 
EXE=../build/nemd.bin

OBJ_SRC=nemd.o \
	getControl.o \
	getStructure.o \
	vdwTables.o    \
	setup.o        \
	computeForces.o\
	integrate_nvtGauss.o\
	writeConfig.o \
	resetAccum.o  \
	sample.o      \
	evolutionProperties.o \
	gradients.o \
	profiles.o \

OBJ_MOD=global.o \
	random.o

$(EXE): $(OBJ_SRC) $(OBJ_MOD)
#	$(FC) -o $@ $(OBJ_SRC) $(OBJ_MOD) #-g -pg -fbounds-check -fbacktrace
	$(FC) -o $@ $(OBJ_SRC) $(OBJ_MOD) # -g -check bounds

%.o: %.f90
#	$(FC)  -c -O3 $<  
	$(FC)  -c -O3 $<

$(OBJ_SRC):$(OBJ_MOD)

clean:
	rm *.o *.mod 


	

