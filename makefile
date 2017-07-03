# Makefile written by Christophe Peyrard from EDF R&D
# Extended to OS X by Yi-Hsiang Hu & Eliot Quin from NREL
# Clean up by Matthieu Ancellin (only tested with gfortran on Linux for the moment)

# Compiler
gtest=$(shell which gfortran 2> /dev/null | grep -o gfortran)
itest=$(shell which ifort 2> /dev/null | grep -o ifort)
# M.A.: this test does not work by me...

MOD_DIR=/tmp/

ifeq ($(gtest), gfortran)
	FC=gfortran
	FFLAGS=  -c                                     # No linker (yet)
	FFLAGS+= -g                                     # Add extra indormations for debugging
	FFLAGS+= -O2                                    # Optimization level
	FFLAGS+= -J$(MOD_DIR)                           # Where to put .mod files
	FFLAGS+= -cpp -DGNUFORT -ffree-line-length-none # Run preprocessor
endif

ifeq ($(itest), ifort)
	FC=ifort
	FFLAGS=-c -cpp
endif

# Output directory
outputdir=./bin

# Default rule: build all
.PHONY: all clean_all remake
all:		mesh preProc solver postProc

clean_all:	clean_mesh clean_preProc clean_solver clean_postProc
			@rm -f $(MOD_DIR)/*.mod

remake:		clean_all all

# Rule to compile f90 file
%.o:	%.f90
		@$(FC) $(FFLAGS) $< -o $@

##################
#  Meshing tool  #
##################

# Sources (relative to DIRM)
SRCM=./Common/Identification.f90\
./Mesh/calCol.f90\
./Mesh/coque.f90\
./Mesh/ExMaillage.f90\
./Mesh/hydre.f90\
./Mesh/Mailleur.f90\
./Mesh/mesh.f90

OBJM=$(SRCM:.f90=.o)

# Rules to build
mesh:		$(OBJM)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/mesh $(OBJM)
			@echo "Meshing tool compilation successful!"

clean_mesh:
			@rm -f $(OBJM)
			@rm -f $(outputdir)/mesh

###################
#  Pre-processor  #
###################

# Sources
SRCP=./Common/Constants.f90\
./Common/Elementary_functions.f90\
./Common/Identification.f90\
./Common/Environment.f90\
./preProcessor/Mesh.f90\
./preProcessor/BodyConditions.f90\
./preProcessor/Integration.f90\
./preProcessor/Main.f90

OBJP=$(SRCP:.f90=.o)

# Rules to build
preProc:	$(OBJP)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/preProc $(OBJP)
			@echo "Preprocessor compilation successful!"

clean_preProc:
			@rm -f $(OBJP)
			@rm -f $(outputdir)/preProc

############
#  Solver  #
############

# Sources
SRCS=./Common/Constants.f90\
./Common/Elementary_functions.f90\
./Common/Bodyconditions.f90\
./Common/Environment.f90\
./Common/Mesh.f90\
./Common/Face.f90\
./Solver/Core/OUTPUT.f90\
./Solver/Core/M_SOLVER.f90\
./Solver/Core/INITIALIZE_GREEN_2.f90\
./Solver/Core/GREEN_1.f90\
./Solver/Core/GREEN_2.f90\
./Solver/Core/SOLVE_BEM_DIRECT.f90\
./Solver/Core/KOCHIN.f90\
./Solver/Core/FREESURFACE.f90\
./Solver/Core/FORCES.f90\
./Solver/NEMOH.f90

OBJS=$(SRCS:.f90=.o)

# Rules to build
solver:		$(OBJS)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/solver $(OBJS)
			@echo "Solver compilation succesful!"

clean_solver:
			@rm -f $(OBJS)
			@rm -f $(outputdir)/solver

####################
#  Post-processor  #
####################

# Sources
SRCO=./Common/Constants.f90\
./Common/Elementary_functions.f90\
./Common/Identification.f90\
./Common/Environment.f90\
./Common/Results.f90\
./Common/Mesh.f90\
./postProcessor/Compute_RAOs.f90\
./postProcessor/IRF.f90\
./postProcessor/Plot_WaveElevation.f90\
./postProcessor/Main.f90

OBJO=$(SRCO:.f90=.o)

# Rules to build
postProc:	$(OBJO)
			@test -d $(outputdir) || mkdir $(outputdir)
			@$(FC) -o $(outputdir)/postProc $(OBJO)
			@echo "Postprocessor compilation succesful!"

clean_postProc:
			@rm -f $(OBJO)
			@rm -f $(outputdir)/postProc

################
#  Test cases  #
################

.PHONY: run_cylinder clean_cylinder
run_cylinder: preProc solver postProc
	$(MAKE) -C Verification/Cylinder/ run

clean_cylinder:
	$(MAKE) -C Verification/Cylinder/ clean

.PHONY: run_nonsymmetrical clean_nonsymmetrical
run_nonsymmetrical: preProc solver postProc
	$(MAKE) -C Verification/NonSymmetrical/ run

clean_nonsymmetrical:
	$(MAKE) -C Verification/NonSymmetrical/ clean

.PHONY: test clean_test
test: preProc solver postProc
	@echo ""
	@echo "Sphere"
	@$(MAKE) --silent -C Verification/QuickTests/1_Sphere/                     test
	@echo ""
	@echo "Sphere using y-symmetry"
	@$(MAKE) --silent -C Verification/QuickTests/2_SymmetricSphere/            test
	@echo ""
	@echo "Sphere in finite depth"
	@$(MAKE) --silent -C Verification/QuickTests/3_FiniteDepthSphere/          test
	@echo ""
	@echo "Sphere in finite depth using y-symmetry"
	@$(MAKE) --silent -C Verification/QuickTests/4_SymmetricFiniteDepthSphere/ test
	@echo ""
	@echo "Alien sphere"
	@$(MAKE) --silent -C Verification/QuickTests/5_AlienSphere/                test

clean_test:
	@$(MAKE) --silent -C Verification/QuickTests/1_Sphere/                     clean
	@$(MAKE) --silent -C Verification/QuickTests/2_SymmetricSphere/            clean
	@$(MAKE) --silent -C Verification/QuickTests/3_FiniteDepthSphere/          clean
	@$(MAKE) --silent -C Verification/QuickTests/4_SymmetricFiniteDepthSphere/ clean
	@$(MAKE) --silent -C Verification/QuickTests/5_AlienSphere/                clean
