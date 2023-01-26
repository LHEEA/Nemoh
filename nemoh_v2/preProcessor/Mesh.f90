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
        INTEGER,ALLOCATABLE :: LastPanel(:)         ! Last panel of each body
        REAL,ALLOCATABLE :: CG(:,:)                 ! Gravity centre of each body
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
        ALLOCATE(Mesh%LastPanel(Mesh%Nbodies))
        ALLOCATE(Mesh%CG(3,Mesh%Nbodies))
        Mesh%X(:,:)=0.
        Mesh%N(:,:)=0.
        Mesh%XM(:,:)=0.
        Mesh%P(:,:)=0 
        Mesh%cPanel(:)=0
        Mesh%A(:)=0.
        Mesh%LastPanel(:)=0
        Mesh%CG(:,:)=0.
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
        DO j=1,MeshTarget%NBodies
            MeshTarget%LastPanel(j)=MeshSource%LastPanel(j)    
            DO i=1,3
                MeshTarget%CG(i,j)=MeshSource%CG(i,j)   
            END DO
        END DO
        END SUBROUTINE CopyTMesh
!       ---
        SUBROUTINE ReadTMesh(Mesh,ID)

#ifdef GNUFORT
#define DIRECTORY FILE
#endif

#ifndef GNUFORT
        USE iflport
#endif
        USE MIdentification
        IMPLICIT NONE
        TYPE(TMesh) :: Mesh
        TYPE(TID) :: ID
        REAL    :: Lad,Length
        INTEGER :: i,j,c,d
        INTEGER :: Npoints,Npanels,Nbodies
        CHARACTER*80 :: meshfile,line
        INTEGER :: lfile
        INTEGER :: M,N
        REAL,DIMENSION(3) :: U,V,W1,W2
        REAL :: A1,A2
        REAL :: calNorme
        REAL :: tX,tY
        LOGICAL :: ex
!       Initialize the mesh structure
        Npoints=0
        Npanels=0
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*) Nbodies
        DO c=1,Nbodies
            READ(10,*)
            READ(10,*)
            READ(10,*) M,N
            Npoints=Npoints+M
            Npanels=Npanels+N
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
        END DO
        CLOSE(10)    
        CALL CreateTMesh(Mesh,Npoints,Npanels,Nbodies)
!       Read surfacic mesh
        Npanels=0
        Npoints=0
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*)
        READ(10,*) 
        DO c=1,Nbodies
            READ(10,*) meshfile
            tX=0.
            tY=0.
            READ(10,*) meshfile
              lfile=LNBLNK(meshfile)
            OPEN(11,FILE=meshfile(1:lfile))
            READ(11,*) M,N
            IF ((c.GT.1).AND.(N.NE.Mesh%Isym)) THEN
                WRITE(*,*) ' Error: there is an inconsistency in the mesh files regarding the xOz symmetries'  
                STOP     
            ELSE
                Mesh%Isym=N
            END IF          
            READ(10,*) M,N
            DO i=1,M
                READ(11,*) d,(Mesh%X(j,Npoints+i),j=1,3)
                Mesh%X(1,Npoints+i)=Mesh%X(1,Npoints+i)+tX
                Mesh%X(2,Npoints+i)=Mesh%X(2,Npoints+i)+tX
            END DO            
            READ(11,*)
            DO i=1,N
                READ(11,*) (Mesh%P(j,Npanels+i),j=1,4)
                DO j=1,4
                    Mesh%P(j,Npanels+i)=Mesh%P(j,Npanels+i)+Npoints
                END DO
                Mesh%cPanel(Npanels+i)=c
            END DO 
            CLOSE(11)
            Npoints=Npoints+M
            Npanels=Npanels+N
            mesh%LastPanel(c)=Npanels
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
            READ(10,*) N
            DO i=1,N
                READ(10,*)
            END DO
        END DO
        CLOSE(10)
!       Perform additional calculations (normal vectors, areas, ...)
        DO i=1,Mesh%Npanels
!           Surface of panel, centre and normal vector
            DO j=1,3
                U(j)=Mesh%X(j,Mesh%P(2,i))-Mesh%X(j,Mesh%P(1,i))
                V(j)=Mesh%X(j,Mesh%P(4,i))-Mesh%X(j,Mesh%P(2,i))
            END DO
            CALL CrossProduct(U,V,W1)
            A1=0.5*calNorme(W1)
            DO j=1,3
                U(j)=Mesh%X(j,Mesh%P(4,i))-Mesh%X(j,Mesh%P(3,i))
                V(j)=Mesh%X(j,Mesh%P(2,i))-Mesh%X(j,Mesh%P(3,i))
            END DO
            CALL CrossProduct(U,V,W2)
            A2=0.5*calNorme(W2)
            Mesh%A(i)=A1+A2
            IF (Mesh%A(i).LT.1.0E-07) THEN
                WRITE(*,'(A,I6,A,E14.7)') 'Error: surface of panel ',i,' is too small (',Mesh%A(i),')'
                STOP
            END IF 
            DO j=1,3
                Mesh%XM(j,i)=1./3*(Mesh%X(j,Mesh%P(1,i))+Mesh%X(j,Mesh%P(2,i))+Mesh%X(j,Mesh%P(4,i)))*A1/Mesh%A(i)+1./3*(Mesh%X(j,Mesh%P(2,i))+Mesh%X(j,Mesh%P(3,i))+Mesh%X(j,Mesh%P(4,i)))*A2/Mesh%A(i)
            END DO
            U=W1+W2
            DO j=1,3
                Mesh%N(j,i)=U(j)/calNorme(U)
            END DO
        END DO
!       Export mesh
        INQUIRE (DIRECTORY=ID%ID(1:ID%lID)//'/mesh', EXIST=ex) 
        IF (.NOT.ex) M=SYSTEM('mkdir '//ID%ID(1:ID%lID)//'/mesh')
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/L12.dat')
        WRITE(10,'(I3,X,I3)') 2,Mesh%Isym
        DO i=1,Mesh%Npoints
            WRITE(10,'(I6,3(X,E14.7))') i,(Mesh%X(j,i),j=1,3)    
        END DO
        WRITE(10,'(I6,3(X,E14.7))') 0,0.,0.,0.
        DO i=1,Mesh%Npanels
            WRITE(10,'(4(X,I6))') (Mesh%P(j,i),j=1,4) 
        END DO
        WRITE(10,'(4(X,I6))') 0,0,0,0
        CLOSE(10)
        OPEN(10,FILE=ID%ID(1:ID%lID)//'/mesh/L10.dat')
        WRITE(10,*) 
        WRITE(10,'(4(X,I6))') Mesh%Isym,Mesh%Npoints,Mesh%Npanels,Mesh%Nbodies
        DO i=1,Mesh%Npanels
            WRITE(10,'(I5,7(X,E14.7))') Mesh%cPanel(i),(Mesh%XM(j,i),j=1,3),(Mesh%N(j,i),j=1,3),Mesh%A(i)
        END DO
        WRITE(10,*) 0
        CLOSE(10)
        OPEN(10,file=ID%ID(1:ID%lID)//'/mesh/Mesh.tec')
        WRITE(10,'(A)') 'VARIABLES="X" "Y" "Z" "NX" "NY" "NZ" "A"'
	    WRITE(10,*) 'ZONE N=',Mesh%Npoints,', E=',Mesh%Npanels,' , F=FEPOINT,ET=QUADRILATERAL'
	    DO i=1,Mesh%Npoints
		    WRITE(10,'(7(2X,E14.7))') (Mesh%X(j,i),j=1,3),0.0,0.0,0.0,0.
	    END DO
	    DO i=1,Mesh%Npanels
		    WRITE(10,'(I6,3(2X,I6))') (Mesh%P(j,i),j=1,4) 
	    END DO
	    WRITE(10,*) 'ZONE t="normales", F=POINT, I=',Mesh%Npanels
	    DO i=1,Mesh%Npanels
		    WRITE(10,'(7(2X,E14.7))') (Mesh%XM(j,i),j=1,3),(Mesh%N(j,i),j=1,3),Mesh%A(i) 
	    END DO
	    CLOSE(10)
	    INQUIRE (DIRECTORY=ID%ID(1:ID%lID)//'/results', EXIST=ex)
	    IF (.NOT.ex) M=SYSTEM('mkdir '//ID%ID(1:ID%lID)//'/results')	    
        END SUBROUTINE ReadTMesh
!       --- 
        SUBROUTINE DeleteTMesh(Mesh)
        IMPLICIT NONE
        TYPE(TMesh) :: Mesh
        DEALLOCATE(Mesh%X,Mesh%N,Mesh%P,Mesh%XM,Mesh%A,Mesh%cPanel,Mesh%LastPanel,Mesh%CG)
        END SUBROUTINE DeleteTMesh  
!       ---
    END MODULE MMesh
!   ---    
    SUBROUTINE CrossProduct(U,V,W)
    IMPLICIT NONE
    REAL,DIMENSION(3) :: U,V,W
    W(1)=U(2)*V(3)-U(3)*V(2)
    W(2)=U(3)*V(1)-U(1)*V(3)
    W(3)=U(1)*V(2)-U(2)*V(1)
    END SUBROUTINE         
!   ---
    REAL FUNCTION calNorme(W)
    IMPLICIT NONE
    REAL,DIMENSION(3) :: W
    calNorme=SQRT(W(1)*W(1)+W(2)*W(2)+W(3)*W(3))
    END FUNCTION             
