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
!   - A. Babarit  
!
!--------------------------------------------------------------------------------------
    MODULE MMesh
!    
    TYPE TMesh
        INTEGER :: Isym                             ! Symmetry about the xOz plane (1 for yes)
        INTEGER :: Npoints                          ! Total number of points in mesh
        INTEGER :: Npanels                          ! Total number of panels in mesh
        INTEGER :: Nbodies                          ! Total number of bodies in mesh
        REAL,ALLOCATABLE :: X(:,:)                  ! Nodes coordinates
        REAL,ALLOCATABLE :: N(:,:)                  ! Normal vectors
        REAL,ALLOCATABLE :: XM(:,:)                 ! Centre of panels
        INTEGER,ALLOCATABLE :: P(:,:)               ! Connectivities
        INTEGER,ALLOCATABLE :: cPanel(:)            ! To which body belongs the panel
        REAL,ALLOCATABLE :: A(:)                    ! Area of panel
    END TYPE TMesh
!    
    CONTAINS 
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTMesh(Mesh,Npoints,Npanels,Nbodies)
        IMPLICIT NONE
        TYPE(TMesh) :: Mesh
        INTEGER :: Npoints,Npanels,Nbodies
        Mesh%Npoints=Npoints
        Mesh%Npanels=Npanels
        Mesh%Nbodies=Nbodies
        ALLOCATE(Mesh%X(3,Mesh%Npoints),Mesh%N(3,Mesh%Npanels),Mesh%XM(3,Mesh%Npanels))
        ALLOCATE(Mesh%P(4,Mesh%Npanels),Mesh%cPanel(Mesh%Npanels),Mesh%A(Mesh%Npanels))
        Mesh%X(:,:)=0.
        Mesh%N(:,:)=0.
        Mesh%XM(:,:)=0.
        Mesh%P(:,:)=0 
        Mesh%cPanel(:)=0
        Mesh%A(:)=0.
        END SUBROUTINE CreateTMesh
!       --- 
        SUBROUTINE CopyTMesh(MeshTarget,MeshSource)
        IMPLICIT NONE
        INTEGER :: i,j
        TYPE(TMesh) :: MeshTarget,MeshSource      
        CALL CreateTMesh(MeshTarget,MeshSource%Npoints,MeshSource%Npanels,MeshSource%Nbodies)
        MeshTarget%Isym=MeshSource%Isym
        DO j=1,MeshTarget%Npoints
            DO i=1,3
                MeshTarget%X(i,j)=MeshSource%X(i,j)
            END DO
        END DO
        DO j=1,MeshTarget%Npanels
            DO i=1,3
                MeshTarget%N(i,j)=MeshSource%N(i,j)
                MeshTarget%XM(i,j)=MeshSource%XM(i,j)
            END DO
            DO i=1,4
                MeshTarget%P(i,j)=MeshSource%P(i,j)
            END DO
            MeshTarget%cPanel(j)=MeshSource%cPanel(j)
            MeshTarget%A(j)=MeshSource%A(j)
        END DO
        END SUBROUTINE CopyTMesh
!       ---
        SUBROUTINE ReadTMesh(Mesh,ID)
#ifndef GNUFORT
        USE iflport
#endif
        USE MIdentification
        IMPLICIT NONE
        TYPE(TMesh) :: Mesh
        TYPE(TID) :: ID
        REAL    :: Lad
        INTEGER :: i,j,c,d,k
        INTEGER :: Npoints,Npanels,Nbodies
!       Initialize the mesh structure
        Npoints=0
        Npanels=0
        Nbodies=0
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/L10.dat')
        READ(10,*)
        READ(10,*) i,Npoints,Npanels,Nbodies
         CLOSE(10) 
        CALL CreateTMesh(Mesh,Npoints,Npanels,Nbodies)
!       Read mesh
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/L12.dat')
        READ(10,*) i,Mesh%Isym
        IF (i.NE.2) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' The mesh file format is not correct. '
            WRITE(*,*) ' Stopping'
            STOP
        END IF
        DO i=1,Npoints
            READ(10,*) j,(Mesh%X(k,i),k=1,3)
        END DO
        READ(10,*)
        DO i=1,Npanels
            READ(10,*) (Mesh%P(k,i),k=1,4)
        END DO
        READ(10,*)
        CLOSE(10)
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/L10.dat')
        READ(10,*)
        READ(10,*) j
        IF (j.NE.Mesh%Isym) THEN
            WRITE(*,*) ' '
            WRITE(*,*) ' The mesh file format is not correct. '
            WRITE(*,*) ' Stopping'
            STOP
        END IF
        DO i=1,Npanels
            READ(10,*) Mesh%cPanel(i),(Mesh%XM(k,i),k=1,3),(Mesh%N(k,i),k=1,3),Mesh%A(i)
        END DO
        CLOSE(10)  
        END SUBROUTINE ReadTMesh
!       --- 
        SUBROUTINE DeleteTMesh(Mesh)
        IMPLICIT NONE
        TYPE(TMesh) :: Mesh
        DEALLOCATE(Mesh%X,Mesh%N,Mesh%P,Mesh%XM,Mesh%A,Mesh%cPanel)
        END SUBROUTINE DeleteTMesh  
!       ---
    END MODULE MMesh
        
