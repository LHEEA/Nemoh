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
	  if(Indiq_solver .eq. 0)ALLOCATE(ZIJ(NFA,NFA+1))

    END SUBROUTINE
