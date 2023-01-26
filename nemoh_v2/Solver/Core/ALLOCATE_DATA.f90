!--------------------------------------------------------------------------------------
!
!   Copyright 2014 Ecole Centrale de Nantes, 1 rue de la Noë, 44300 Nantes, France
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!       http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License. 
!
!   Contributors list:
!   - J. Singh
!   - G. Delhommeau 
!
!--------------------------------------------------------------------------------------
    SUBROUTINE ALLOCATE_DATA

      USE COM_VAR

	! Array of cordinates of points
	  ALLOCATE(X(NP),Y(NP),Z(NP))
	! Array of centre of gravity for each panel
	!  ALLOCATE(XG(2**NSYMY*NFA),YG(2**NSYMY*NFA),ZG(2**NSYMY*NFA))
	  ALLOCATE(XG(2*NFA),YG(2*NFA),ZG(2*NFA))
	! Array for surface area of the panels
	  ALLOCATE(AIRE(NFA))
	!vertices of the panel, size NFA(no of facettes)
	  ALLOCATE(M1(NFA),M2(NFA),M3(NFA),M4(NFA))
	! Normal vector
	  ALLOCATE(XN(NFA),YN(NFA),ZN(NFA))
	! DIST(I) is the maximum distance (projection on XY plane!) of a panel I from other points 
	  ALLOCATE(DIST(NFA),TDIS(NFA))
	!  
	  ALLOCATE(ZIGB(NFA),ZIGS(NFA)) !sources
      ALLOCATE(ZPB(NFA),ZPS(NFA))   !potential
	! Linear complex matrix to be solved by direct solver
        !For GMRES workspace is calculated inside the main subroutine
	  if(Indiq_solver .eq. 0) ALLOCATE(ZIJ(NFA,2*NFA),AInv(NFA,NFA,2**NSYMY))

    END SUBROUTINE
