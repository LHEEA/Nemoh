# makefile written by Christophe Peyrard from EDF R&D
# extended to OS X by Yi-Hsiang Hu & Eliot Quin from NREL

#COMPILATEUR
gtest=$(shell which gfortran 2> /dev/null | grep -o gfortran)
itest=$(shell which ifort 2> /dev/null | grep -o ifort)

outputdir=./bin

ifeq ($(gtest), gfortran)
	FC=gfortran
	FFLAGS=-cpp -DGNUFORT -O2 -ffree-line-length-none -c
endif

ifeq ($(itest), ifort)
	FC=ifort
	FFLAGS=-c -cpp
endif

#SOURCES FORTRAN Mesh(modules de maillage)
SRCM=./Common/Identification.f90\
./Mesh/calCol.f90\
./Mesh/coque.f90\
./Mesh/ExMaillage.f90\
./Mesh/hydre.f90\
./Mesh/Mailleur.f90\
./Mesh/mesh.f90\

# LISTE DES .o preProc
#TRANSFORME f90 en o
OBJM=$(SRCM:.f90=.o)

#Liste pour transformer ./*/*.o en .o dans le OBJP
OBJM2:=$(subst .f90,.o,$(notdir ${SRCM}))

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

#Liste pour transformer ./*/*.o en .o dans le OBJP
OBJP2:=$(subst .f90,.o,$(notdir ${SRCP}))

#SOURCES FORTRAN Solver(modules de preprocessing)
SRCS=./Solver/Core/COM_VAR.f90\
./Common/Environment.f90\
./Common/Identification.f90\
./Common/Mesh.f90\
./Solver/Core/Bodyconditions.f90\
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

#Liste pour transformer ./*/*.o en .o dans le OBJS
OBJS2:=$(subst .f90,.o,$(notdir ${SRCS}))


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

#Liste pour transformer ./*/*.o en .o dans le OBJP
OBJO2:=$(subst .f90,.o,$(notdir ${SRCO}))

build: bin msh pre solver post clean

bin:
	mkdir -p $(outputdir)

#
#Build Mesh executable
msh:	mesh
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
mesh:	$(OBJM)
		$(FC) -o $(outputdir)/mesh $(OBJM2)
#
#Build preProc executable
pre:	preProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
preProc:	$(OBJP)
		$(FC) -o $(outputdir)/preProc $(OBJP2)


#
#Build solver executable
solver:	Nemoh
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
Nemoh:	$(OBJS)
		$(FC) -o $(outputdir)/solver $(OBJS2)


#
#Build postProc executable
post:	postProc
#Rules to Build MAIN EXECUTABLE  (dependances et regle d'execution)
postProc:	$(OBJO)
		$(FC) -o $(outputdir)/postProc $(OBJO2)

# Rules for .f comiplation
.f.o:
	$(FC) $(FFLAGS) $<
%.o:	%.f90
	$(FC) $(FFLAGS) $<


#Copy to local bin directory
install: build
	cp $(outputdir)/* ~/bin/

# Remove *.o and main executable
clean:
	@rm -f *.o *.mod


DOCKER_RUN:=\
	docker run --rm \
	-u $(shell id -u):$(shell id -g) \
	-v $(shell pwd):/opt/share \
	-w /opt/share \
	gcc:8.2.0 /bin/bash -c

# Build targets with docker for linux arm64 architecture
docker_build:
	${DOCKER_RUN} "make"

# Run demo with docker for linux arm64 architecture
docker_demo:
	${DOCKER_RUN} \
	"cd Verification/Cylinder && ../../bin/preProc && \
	cd .. && ../bin/solver && \
	../bin/postProc && \
	ls"
