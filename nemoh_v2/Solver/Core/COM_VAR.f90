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
!   - G. Delhommeau
!   - J. Singh
!   - A. Babarit
!
!--------------------------------------------------------------------------------------
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
    ! Source distribution
    COMPLEX, DIMENSION(:), ALLOCATABLE :: ZIGB,ZIGS
    
    ! --- Solver ------------------------------------------------------
    ! Which solver: (0) direct solver, (1): GMRES
    INTEGER:: Indiq_solver
    ! Linear complex matrix to be solved
    COMPLEX, DIMENSION(:,:), ALLOCATABLE :: ZIJ
    ! Storage of inverter matrix to speed up computations (influence coefficients depending only on w)
    COMPLEX, DIMENSION(:,:,:), ALLOCATABLE :: AInv
    REAL w_previous
    !
    REAL :: FSP,FSM,VSXP,VSYP,VSZP,VSXM,VSYM,VSZM
    ! Variable for storage of Greens function
    REAL :: SP1,SM1,SP2,SM2 
    ! Variable for storage of Greens function
    REAL :: VSXP1,VSXP2,VSYP1,VSYP2,VSZP1,VSZP2 
    REAL :: VSXM1,VSXM2,VSYM1,VSYM2,VSZM1,VSZM2
    !Values used in interpolation of infinite part of the Greens function
    INTEGER,PARAMETER :: NPINTE=251
    REAL :: XR(328),XZ(46),APD1X(328,46),APD1Z(328,46),APD2X(328,46),APD2Z(328,46)   
    INTEGER:: NQ
    REAL:: CQ(NPINTE),QQ(NPINTE),AMBDA(31),AR(31)
    
    ! --- Reading and writing units -----------------------------------
    ! File for visualization of potential
    INTEGER :: Sav_potential
       
  END MODULE COM_VAR