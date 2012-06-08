  MODULE COM_VAR

    IMPLICIT NONE

    ! --- Environment --------------------------------------------------
    ! Volumique mass of the fluid in KG/M**3
    REAL :: RHO
    ! Gravity
    REAL :: G
    ! depth of the domain
    REAL :: Depth
    ! Coordinates where the waves are measured
    REAL :: XEFF,YEFF,ZEFF    

    ! --- Geometrical quantities ---------------------------------------
    ! Mesh file
    CHARACTER*80 :: MESHFILE
    INTEGER :: LFILE
    ! No. of points in the surface mesh
    INTEGER :: NP
    ! No of total panels in the surface mesh
    INTEGER :: NFA
    !
    INTEGER :: IMX
    INTEGER :: IXX
    !!!!!!!!!
    !Symmetry: 0 for no symmetry and1 for symmetry
    INTEGER:: NSYMY
    !
    REAL :: ZER
    ! DIST(I) is the maximum distance (projection on XY plane!) of a panel I from other points 
    REAL, DIMENSION(:), ALLOCATABLE :: DIST,TDIS
    !vertices of the panel, size NFA(no of facettes)
    INTEGER, DIMENSION(:), ALLOCATABLE :: M1,M2,M3,M4
    ! Array of cordinates of points of size NP (no of points)
    REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z
    ! Normal vector
    REAL, DIMENSION(:), ALLOCATABLE :: XN,YN,ZN
    ! Array of centre of gravity for each panel
    REAL, DIMENSION(:), ALLOCATABLE :: XG,YG,ZG
    ! Array for surface area of the panels
    REAL, DIMENSION(:), ALLOCATABLE :: AIRE
    
    ! --- Boundary value problem ---------------------------------------
    ! normal velocity array as input
!    REAL, DIMENSION(:), ALLOCATABLE :: NVEL
    ! period
    REAL :: T
    ! Computed and computed potential on original (B) and symmetric boundary(S)
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ZPB,ZPS
    ! Source ditribution
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB,ZIGS
    
    ! --- Solver ------------------------------------------------------
    ! Which solver: (0) direct solver, (1): GMRES
    INTEGER:: Indiq_solver
    !Some GMRES parameters
    INTEGER::IRES,MAXIT
    REAL::TOL_GMRES
   ! Linear complex matrix to be solved
    COMPLEX, DIMENSION(:, :), ALLOCATABLE :: ZIJ
    REAL :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
    ! Variable for storage of Greens function
    REAL :: SP1,SM1,SP2,SM2 
    ! Variable for storage of Greens function
    REAL :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2 
    REAL :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
    !Values used in interpolation of infinite part of the Greens function
    REAL :: XR(328),XZ(46),APD1X(328,46),APD1Z(328,46),APD2X(328,46),APD2Z(328,46)   
    INTEGER:: NQ
    REAL:: CQ(101),QQ(101),AMBDA(31),AR(31)
    
    ! --- Reading and writing units -----------------------------------
    ! File for visualization of potential
    INTEGER :: Sav_potential
       
  END MODULE COM_VAR