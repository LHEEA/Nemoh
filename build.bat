
gfortran -g -c Common/Identification.f90
gfortran -g -c Common/Environment.f90
gfortran -g -ffree-line-length-none -cpp -DGNUFORT -c preProcessor/Mesh.f90
gfortran -g -ffree-line-length-none -c preProcessor/BodyConditions.f90
gfortran -g -c preProcessor/Integration.f90
gfortran -g -ffree-line-length-none -c preProcessor/Main.f90
gfortran -g Identification.o Environment.o Mesh.o BodyConditions.o Integration.o Main.o -o nemoh_preproc


gfortran -c Common/Identification.f90
gfortran -cpp -DGNUFORT -c Common/Mesh.f90
gfortran -ffree-line-length-none -c Solver/Core/Bodyconditions.f90
gfortran -c Solver/Core/COM_VAR.f90
gfortran -ffree-line-length-none -c Solver/Core/OUTPUT.f90
gfortran -c Solver/Core/ELEMENTARY_FNS.f90
gfortran -c Solver/Core/COMPUTE_GREEN_INFD.f90
gfortran -c Solver/Core/M_SOLVER.f90
gfortran -c Solver/Core/SOLVE_BEM_INFD_DIRECT.f90
gfortran -c Solver/Core/COMPUTE_GREEN_FD.f90
gfortran -c Solver/Core/SOLVE_BEM_FD_DIRECT.f90
gfortran -c Solver/Core/COMPUTE_KOCHIN.f90
gfortran -c Solver/Core/COMPUTE_GREEN_FREESURFACE.f90
gfortran -c Solver/Core/COMPUTE_POTENTIAL_DOMAIN.f90
gfortran -c Solver/Core/SOLVE_BEM.f90
gfortran -c Solver/Core/PREPARE_MESH.f90
gfortran -c Solver/Core/ALLOCATE_DATA.f90
gfortran -c Solver/Core/DEALLOCATE_DATA.f90
gfortran -c Solver/Core/INITIALIZATION.f90
gfortran -cpp -DGNUFORT -ffree-line-length-none -c Solver/NEMOH.f90

gfortran Mesh.o Identification.o Bodyconditions.o COM_VAR.o OUTPUT.o ELEMENTARY_FNS.o COMPUTE_GREEN_INFD.o M_SOLVER.o SOLVE_BEM_INFD_DIRECT.o COMPUTE_GREEN_FD.o SOLVE_BEM_FD_DIRECT.o SOLVE_BEM.o PREPARE_MESH.o INITIALIZATION.o ALLOCATE_DATA.o COMPUTE_KOCHIN.o COMPUTE_GREEN_FREESURFACE.o COMPUTE_POTENTIAL_DOMAIN.o DEALLOCATE_DATA.o NEMOH.o -o nemoh

