Written by Yi-Hsiang Yu (NREL), Feb 2016
--------------------------------------------------------------

For all the OS X machine you need to

	1. Install gfortran if you don't have one
	https://gcc.gnu.org/wiki/GFortranBinaries
	2. Use the new makefile -> make -f makefile
	3. Copy the exe file to the matlab routines folder under Mesh and Nemoh
	4. Fix a MATLAB Lib issue by creating a soft link for the DYLIB: 
	ln -s /usr/local/lib/libgfortran.3.dylib /Applications/MATLAB_R20XXx.app/sys/os/maci64/libgfortran.3.dylib

Note: Four files has been changed from the original package:
	1. makefile
	2. ./matlab routines/axiMesh.m
	2. ./matlab routines/Mesh.m
	3. ./matlab routines/Nemoh.m 