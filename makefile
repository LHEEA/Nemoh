# makefile written by Christophe Peyrard from EDF R&D

#COMPILATEUR  
FC=ifort
#OPTIONS  
FFLAGS= -c
#LIBRAIRIES  
#LFLAGS= -L./ -lVAWTBEMT -lstdc++

#SOURCES FORTRAN preProc(modules de preprocessing)
SRCP=./Common/Identification.f90\
./Common/Environment.f90\
./preProcessor/Mesh.f90\
./preProcessor/BodyConditions.f90\
./preProcessor/Integration.f90\
./preProcessor/Main.f90\

# LISTE DES .o preProc
#TRANSFORME f90 en o  
OBJP=$(SRCP:.f90=.o)

#Liste pour transformer ./*/*.o en .o dans le OBJP (cf Yoann pour automatisation)  
OBJP2=Identification.o\
Environment.o\
Mesh.o\
BodyConditions.o\
Integration.o\
Main.o\


#SOURCES FORTRAN Solver(modules de preprocessing)
SRCS=./Solver/Core/COM_VAR.f90\
./Common/Environment.f90\
./Common/Identification.f90\
./Common/Mesh.f90\
./Solver/Core/BODYCONDITIONS.f90\
./Solver/Core/PREPARE_MESH.f90\
./Solver/Core/INITIALIZATION.f90\
./Solver/Core/OUTPUT.f90\
./Solver/Core/ELEMENTARY_FNS.f90\
./Solver/Core/M_SOLVER.f90\
./Solver/Core/ALLOCATE_DATA.f90\
./Solver/Core/COMPUTE_GREEN_INFD.f90\
./Solver/Core/SOLVE_BEM_INFD_DIRECT.f90\
./Solver/Core/COMPUTE_GREEN_FD.f90\
./Solver/Core/SOLVE_BEM_FD_DIRECT.f90\
./Solver/Core/SOLVE_BEM.f90\
./Solver/Core/COMPUTE_KOCHIN.f90\
./Solver/Core/COMPUTE_GREEN_FREESURFACE.f90\
./Solver/Core/COMPUTE_POTENTIAL_DOMAIN.f90\
./Solver/NEMOH.f90\
./Solver/Core/DEALLOCATE_DATA.f90\
#./Solver/Core/SOLVE_BEM_FD_ITERATIVE.f90
#./Solver/Core/SOLVE_BEM_INFD_ITERATIVE.f90


# LISTE DES .o preProc
#TRANSFORME f90 en o  
OBJS=$(SRCS:.f90=.o)

#Liste pour transformer ./*/*.o en .o dans le OBJS (cf Yoann pour automatisation)  
OBJS2=COM_VAR.o\
Environment.o\
Identification.o\
Mesh.o\
BODYCONDITIONS.o\
PREPARE_MESH.o\
INITIALIZATION.o\
OUTPUT.o\
ELEMENTARY_FNS.o\
M_SOLVER.o\
ALLOCATE_DATA.o\
COMPUTE_GREEN_INFD.o\
SOLVE_BEM_INFD_DIRECT.o\
COMPUTE_GREEN_FD.o\
SOLVE_BEM_FD_DIRECT.o\
SOLVE_BEM.o\
COMPUTE_KOCHIN.o\
COMPUTE_GREEN_FREESURFACE.o\
COMPUTE_POTENTIAL_DOMAIN.o\
NEMOH.o\
DEALLOCATE_DATA.o\
#./Solver/Core/SOLVE_BEM_FD_ITERATIVE.o
#./Solver/Core/SOLVE_BEM_INFD_ITERATIVE.o


#SOURCES FORTRAN preProc(modules de preprocessing)
SRCO=./Common/Identification.f90\
./Common/Environment.f90\
./Common/Results.f90\
./Common/Mesh.f90\
./postProcessor/Compute_RAOs.f90\
./postProcessor/IRF.f90\
./postProcessor/Plot_WaveElevation.f90\
./postProcessor/Main.f90\

# LISTE DES .o preProc
#TRANSFORME f90 en o  
OBJO=$(SRCO:.f90=.o)

#Liste pour transformer ./*/*.o en .o dans le OBJP (cf Yoann pour automatisation)  
OBJO2=Identification.o\
Environment.o\
Results.o\
Mesh.o\
Compute_RAOs.o\
IRF.o\
Plot_WaveElevation.o\
Main.o\

build: pre solver post clean

#
#Build preProc executable
pre:	preProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
preProc:	$(OBJP) 
		$(FC) -o preProc $(OBJP2)


#
#Build solver executable
solver:	Nemoh
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
Nemoh:	$(OBJS) 
		$(FC) -o solver $(OBJS2)


#
#Build postProc executable
post:	postProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
postProc:	$(OBJO) 
		$(FC) -o postProc $(OBJO2)

# Rules for .f comiplation
.f.o:
	$(FC) $(FFLAGS) $<
%.o:	%.f90
	$(FC) $(FFLAGS) $<

# Remove *.o and main executable
clean:
	rm *.o *.mod   
