!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - BVP solver - January 2014
!
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
!   - P. Guével
!   - J.C. Daubisse
!   - J. Singh
!   - A. Babarit  
!
!--------------------------------------------------------------------------------------
!
    PROGRAM Main
!   
    USE MIdentification
    USE MMesh
    USE MBodyConditions
#ifndef GNUFORT
    USE iflport
#endif
    USE SOLVE_BEM
    USE OUTPUT
    USE INITIALIZATION
!    
    IMPLICIT NONE
!   ID
    TYPE(TID)               :: ID
!   Array sizes
    INTEGER                 :: NFA,NSYMY    ! Number of panels and symmetry about xOz plane
!   Body conditions
    TYPE(TBodyConditions)   :: BodyConditions
!   Kochin function
    INTEGER                 :: Ntheta    
    REAL,DIMENSION(:),ALLOCATABLE :: Theta
!   Meshes
    TYPE(TMesh) :: Mesh 
    TYPE(TMesh) :: MeshFS
!   Nemoh
    REAL                    :: T
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL,PRESSURE
    COMPLEX,DIMENSION(:),ALLOCATABLE :: HKochin
!   Results
    REAL :: XEFF,YEFF
    INTEGER :: Nintegration
    REAL,DIMENSION(:,:),ALLOCATABLE :: NDS
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: Force
    REAL,DIMENSION(:),ALLOCATABLE :: line
!   Locals
    REAL                    :: PI
    INTEGER			        :: c,d,M,l,i,j,k
    REAL                    :: Discard
    COMPLEX,PARAMETER       :: II=CMPLX(0.,1.)

!
!   --- Initialisation -------------------------------------------------------------------------------------------------------------------------------------------------------------------
!     
    PI=4.*ATAN(1.)
    WRITE(*,*) ' '
    WRITE(*,'(A,$)') '  -> Initialisation ' 
!   Read case ID
    CALL ReadTID(ID,'ID.dat')
!   Read Mesh
    CALL ReadTMesh(Mesh,ID)
!   Read Body Conditions
    CALL ReadTBodyConditions(BodyConditions,Mesh%Npanels*2**Mesh%Isym,ID%ID(1:ID%lID)//'/Normalvelocities.dat')
!   Initialise Nemoh
    CALL INITIALIZE(ID,NFA,NSYMY,XEFF,YEFF,Mesh)
    ALLOCATE(NVEL(NFA*2**NSYMY),PRESSURE(NFA*2**NSYMY))
    WRITE(*,'(A,$)') '.'
!   Initialise Force matrix
    print *, ID%ID(1:ID%lID)//'/mesh/Integration.dat'
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Integration.dat')
    READ(10,*) Nintegration
    ALLOCATE(NDS(Nintegration,NFA*2**NSYMY))
    DO i=1,Nintegration
        READ(10,*) (NDS(i,j),j=1,NFA*2**NSYMY)
    END DO
    CLOSE(10)
!   Initialise Kochin function calculation
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Kochin.dat')
    READ(10,*) Ntheta 
    ALLOCATE(Theta(Ntheta))
    IF (Ntheta.GT.0) THEN        
        DO j=1,Ntheta
            READ(10,*) Theta(j)
        END DO        
    END IF
    ALLOCATE(HKochin(NTheta))
    CLOSE(10)
!   Initialise free surface calculation points
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/Freesurface.dat')
    READ(10,*) MeshFS%Npoints,MeshFS%Npanels 
    IF (MeshFS%Npoints.GT.0) THEN
        CALL CreateTMesh(MeshFS,MeshFS%Npoints,MeshFS%Npanels,1)
        DO j=1,MeshFS%Npoints
            READ(10,*) MeshFS%X(1,j),MeshFS%X(2,j)
        END DO
        DO j=1,MeshFS%Npanels
            READ(10,*) MeshFS%P(1,j),MeshFS%P(2,j),MeshFS%P(3,j),MeshFS%P(4,j)
        END DO
    END IF
    CLOSE(10)
!   Initialise results table
    ALLOCATE(Force(Nintegration,Bodyconditions%Nproblems))
    Force(:,:)=0.
    WRITE(*,*) '. Done !'
    WRITE(*,*) ' '
!
!   --- Solve BVPs and calculate forces -------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Solve BVPs and calculate forces ' 
    WRITE(*,*) ' '
    DO j=1,BodyConditions%Nproblems 
        WRITE(*,'(A,I5,A,I5,A,$)') ' Problem ',j,' / ',BodyConditions%Nproblems,' .'
        DO c=1,Mesh%Npanels*2**Mesh%Isym
            NVEL(c)=BodyConditions%NormalVelocity(c,j)
        END DO
!       Solve BVP
        CALL SOLVE_BVP(j,ID,2.*PI/BodyConditions%Omega(j),NVEL,PRESSURE,BodyConditions%Switch_Kochin(j),NTheta,Theta,HKochin,BodyConditions%Switch_Freesurface(j),MeshFS,BodyConditions%Switch_Potential(j))
!       Calculate force coefficients
        DO i=1,Nintegration
            DO c=1,Mesh%Npanels*2**Mesh%Isym
                Force(i,j)=Force(i,j)-PRESSURE(c)*NDS(i,c)
            END DO
        END DO
        WRITE(*,*) '. Done !'      
    END DO    
    WRITE(*,*) ' ' 
    CLOSE(10)
    CLOSE(12)
!
!   --- Save results -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Save results ' 
    WRITE(*,*) ' '
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/Forces.dat')
    ALLOCATE(line(BodyConditions%Nproblems*2))
    DO c=1,Nintegration
        DO j=1,BodyConditions%Nproblems
            IF (Bodyconditions%Switch_type(j).NE.-1.) THEN
                line(2*j-1)=ABS(Force(c,j))
                line(2*j)=ATAN2(IMAG(Force(c,j)),REAL(Force(c,j)))
            ELSE
                line(2*j-1)=IMAG(Force(c,j))/BodyConditions%Omega(j)
                line(2*j)=-REAL(Force(c,j))
            END IF
        END DO
        WRITE(10,*) (line(j),j=1,2*BodyConditions%Nproblems)      
    END DO
    CLOSE(10)    
!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    DEALLOCATE(NVEL,PRESSURE)    
    DEALLOCATE(Force)
    DEALLOCATE(Theta,HKochin)
    IF (MeshFS%Npoints.GT.0) CALL DeleteTMesh(MeshFS)
    CALL DEALLOCATE_DATA
!
    END PROGRAM Main
!----------------------------------------------------------------
