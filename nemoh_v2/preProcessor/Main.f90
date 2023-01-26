!--------------------------------------------------------------------------------------
!
!   NEMOH V1.0 - preProcessor - January 2014
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
!   - A. Babarit  
!
!--------------------------------------------------------------------------------------
!
    PROGRAM Main
!
    USE MEnvironment
    USE MIdentification
    USE MMesh
    USE BodyConditions
    USE Integration
!
    IMPLICIT NONE
!
    TYPE(TID) :: ID                     ! Calculation identification data
    TYPE(TMesh) :: Mesh                 ! Mesh data
    TYPE(TEnvironment) :: Environment   ! Environment data   
!   Wave frequencies
    INTEGER :: Nw
    REAL :: wmin,wmax
    REAL,DIMENSION(:),ALLOCATABLE :: w
!   TYPE TCase
    TYPE TCase
        INTEGER :: Body
        INTEGER :: ICase
        INTEGER :: Mode
        REAL,DIMENSION(3) :: Direction,Axis 
    END TYPE
!   Radiation cases
    INTEGER :: Nradiation
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: RadCase
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: NormalVelocity
!   Diffraction cases
    INTEGER :: Nbeta
    REAL :: betamin,betamax
    REAL,DIMENSION(:),ALLOCATABLE :: beta
    COMPLEX,DIMENSION(:),ALLOCATABLE :: Pressure
!   Force integration cases
    INTEGER :: Switch_Potential
    INTEGER :: Nintegration
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: IntCase
    REAL,DIMENSION(:),ALLOCATABLE :: NDS
    REAL,DIMENSION(:,:),ALLOCATABLE :: FNDS 
!   Froude Krylov forces
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: FKforce 
!   Free surface visualisation
    INTEGER :: Switch_FreeSurface
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
!   Kochin function
    INTEGER :: Switch_Kochin
    INTEGER :: NTheta
    REAL :: Thetamin,Thetamax
!   Other local variables
    INTEGER :: M,N
    INTEGER :: i,j,c,d,k
    INTEGER :: jrad,jint
    REAL :: PI
!
    PI=4.*ATAN(1.)
!
!   --- Initialize and read input datas ----------------------------------------------------------------------------------------
!
    CALL ReadTID(ID,'ID.dat')
    CALL ReadTMesh(Mesh,ID)   
    CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/Nemoh.cal')
!   Read number of radiation cases, number of force integration cases, range of considered periods and diffraction cases
!   an directions for Kochin function
    Nradiation=0
    Nintegration=0
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
    DO c=1,7
        READ(10,*)
    END DO
    DO c=1,Mesh%Nbodies
        DO i=1,3
            READ(10,*)
        END DO
        READ(10,*) M
        Nradiation=Nradiation+M
        DO i=1,M
            READ(10,*)
        END DO
        READ(10,*) M
        Nintegration=Nintegration+M
        DO i=1,M
            READ(10,*)
        END DO
        READ(10,*) M
        DO i=1,M
            READ(10,*)
        END DO
    END DO
    READ(10,*)
    READ(10,*) Nw,wmin,wmax
    ALLOCATE(w(Nw))
    IF (Nw.GT.1) THEN
        DO j=1,Nw
            w(j)=wmin+(wmax-wmin)*(j-1)/(Nw-1)
        END DO
    ELSE
        w(1)=wmin
    END IF
    READ(10,*) Nbeta,betamin,betamax
    ALLOCATE(beta(Nbeta))
    IF (Nbeta.GT.1) THEN
        DO j=1,Nbeta
            beta(j)=(betamin+(betamax-betamin)*(j-1)/(Nbeta-1))*PI/180.
        END DO
    ELSE
        beta(1)=betamin*PI/180.
    END IF
    READ(10,*)
    READ(10,*)
    READ(10,*) Switch_Potential
    READ(10,*) Ntheta,thetamin,thetamax
    IF (NTheta.GT.0) THEN
        Switch_Kochin=1
    ELSE
        Switch_Kochin=0
    END IF
    READ(10,*) Nx,Ny,Lx,Ly
    IF (Nx.GT.0) THEN
        Switch_FreeSurface=1
    ELSE
        Switch_FreeSurface=0
    END IF
    CLOSE(10)
!   re-read input file and store radiation and integration cases
    ALLOCATE(RadCase(Nradiation))
    ALLOCATE(IntCase(Nintegration))
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
    DO c=1,7
        READ(10,*)
    END DO
    jrad=0
    jint=0
    DO c=1,Mesh%Nbodies
        DO i=1,3
            READ(10,*)
        END DO
        READ(10,*) M
        DO i=1,M
            READ(10,*) RadCase(jrad+i)%ICase,(RadCase(jrad+i)%Direction(j),j=1,3),(RadCase(jrad+i)%Axis(j),j=1,3)
            RadCase(jrad+i)%Body=c
            RadCase(jrad+i)%Mode=i
        END DO
        jrad=jrad+M
        READ(10,*) M
        DO i=1,M
            READ(10,*) IntCase(jint+i)%ICase,(IntCase(jint+i)%Direction(j),j=1,3),(IntCase(jint+i)%Axis(j),j=1,3)
            IntCase(jint+i)%Body=c
            IntCase(jint+i)%Mode=i
        END DO
        jint=jint+M
        READ(10,*) M
        DO i=1,M
            READ(10,*)
        END DO
    END DO
    CLOSE(10)
!   Print summary of calculation case
    WRITE(*,*) ' '
    WRITE(*,*) ' Summary of calculation'
    WRITE(*,*) ' '
    IF (Environment%Depth.GT.0.) THEN
        WRITE(*,'(A,F7.2,A)') '  ->  Water depth = ',Environment%Depth,' m'
    ELSE
        WRITE(*,'(A)') '  ->  Infinite water depth'
    END IF
    WRITE(*,'(A,I5,A,F7.4,A,F7.4)') '  ->',Nw,' wave frequencies from ',w(1),' to ',w(Nw)
    WRITE(*,'(A,I5,A,F7.4,A,F7.4)') '  ->',Nbeta,' wave directions from  ',beta(1),' to ',beta(Nbeta)
    WRITE(*,'(A,I5,A)') '  ->',Nradiation,' radiation problems'
    WRITE(*,'(A,I5,A)') '  ->',Nintegration,' forces'
    WRITE(*,*) ' '
!
!   --- Generate force integration file ----------------------------------------------------------------------------------------
!
    ALLOCATE(FNDS(Nintegration,Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(NDS(Mesh%Npanels*2**Mesh%Isym))
    DO j=1,Nintegration
        CALL ComputeNDS(Mesh,IntCase(j)%Body,IntCase(j)%Icase,IntCase(j)%Direction,IntCase(j)%Axis,NDS)
        DO c=1,Mesh%Npanels*2**Mesh%Isym
            FNDS(j,c)=NDS(c)
        END DO
    END DO 
    DEALLOCATE(NDS)
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/mesh/Integration.dat')
    WRITE(11,*) Nintegration
    DO j=1,Nintegration
        WRITE(11,*) (FNDS(j,c),c=1,Mesh%Npanels*2**Mesh%Isym)
    END DO
    CLOSE(11)
!
!   --- Generate body conditions and calculate FK forces ----------------------------------------------------------------------------------------
!
    ALLOCATE(NVEL(Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(PRESSURE(Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(FKForce(Nw,Nbeta,Nintegration))
    ALLOCATE(NormalVelocity(Mesh%Npanels*2**Mesh%Isym,(Nbeta+Nradiation)*Nw))
    DO i=1,Nw
        DO j=1,Nbeta
            CALL ComputeDiffractionCondition(Mesh,w(i),Beta(j),Environment,PRESSURE,NVEL)
            DO c=1,Mesh%Npanels*2**Mesh%Isym
                NormalVelocity(c,j+(i-1)*(Nbeta+Nradiation))=NVEL(c)
            END DO 
!           Calculate the corresponding FK forces
            DO k=1,Nintegration
                FKForce(i,j,k)=0.
                DO c=1,Mesh%nPanels*2**Mesh%Isym
                    FKForce(i,j,k)=FKForce(i,j,k)-PRESSURE(c)*FNDS(k,c)    
                END DO
            END DO
        END DO
        DO j=1,Nradiation
            CALL ComputeRadiationCondition(Mesh,RadCase(j)%Body,RadCase(j)%Icase,RadCase(j)%Direction,RadCase(j)%Axis,NVEL)
            DO c=1,Mesh%Npanels*2**Mesh%Isym
                NormalVelocity(c,j+Nbeta+(i-1)*(Nbeta+Nradiation))=NVEL(c)
            END DO 
        END DO
    END DO
    CLOSE(11)        
    DEALLOCATE(PRESSURE,NVEL,FNDS)
!
!   --- Save body conditions ----------------------------------------------------------------------------------------
!
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/Normalvelocities.dat')
    WRITE(11,*) (Nbeta+Nradiation)*Nw
    WRITE(11,*) ((w(i),j=1,Nbeta+Nradiation),i=1,Nw)
    DO i=1,Nw
        WRITE(11,*) (/ (/(beta(j),j=1,Nbeta)/),(-1.,j=1,Nradiation) /)
    ENDDO
    WRITE(11,*) ((Switch_Potential,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_Freesurface,j=1,Nbeta+Nradiation),i=1,Nw)
    WRITE(11,*) ((Switch_Kochin,j=1,Nbeta+Nradiation),i=1,Nw)
    DO c=1,Mesh%Npanels*2**Mesh%Isym
        WRITE(11,*) (REAL(NormalVelocity(c,j)),IMAG(NormalVelocity(c,j)),j=1,(Nbeta+Nradiation)*Nw)
    END DO
    CLOSE(11)        
    DEALLOCATE(NormalVelocity)
!
!   --- Save FK forces ----------------------------------------------------------------------------------------
!
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/FKForce.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    DO k=1,Nintegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',IntCase(k)%Body,k,')" "angle(F',IntCase(k)%Body,k,')"'
    END DO
    DO c=1,Nbeta
        WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="FKforce - beta = ',beta(c)*180./PI,'",I=',Nw,',F=POINT'
        DO i=1,Nw
            WRITE(10,'(80(X,E14.7))') w(i),(ABS(FKForce(i,c,k)),ATAN2(IMAG(FKForce(i,c,k)),REAL(FKForce(i,c,k))),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/FKForce.dat')
    DO k=1,Nintegration
        WRITE(10,*) ((ABS(FKForce(i,c,k)),ATAN2(IMAG(FKForce(i,c,k)),REAL(FKForce(i,c,k))),c=1,Nbeta),(0.*c,0.*c,c=1,Nradiation),i=1,Nw) 
    END DO
    CLOSE(10)
    DEALLOCATE(FKForce)
!
!   --- Generate Free Surface visualisation file ----------------------------------------------------------------------
!
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/mesh/Freesurface.dat')
    WRITE(11,*) Nx*Ny,(Nx-1)*(Ny-1)
    DO i=1,Nx
        DO j=1,Ny
            WRITE(11,'(3(X,E14.7))') -0.5*Lx+Lx*(i-1)/(Nx-1),-0.5*Ly+Ly*(j-1)/(Ny-1),0.
        END DO
    END DO  
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(11,'(4(X,I7))') j+(i-1)*Ny,j+1+(i-1)*Ny,j+1+i*Ny,j+i*Ny
        END DO
    END DO 
    CLOSE(11)
!
!   --- Generate Kochin file ----------------------------------------------------------------------------------------
!
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/mesh/Kochin.dat')
    WRITE(11,*) NTheta    
    IF (Ntheta.GT.0) THEN
        IF (NTheta.GT.1) THEN
            DO j=1,NTheta
                WRITE(11,*) (Thetamin+(Thetamax-Thetamin)*(j-1)/(NTheta-1))*PI/180.
            END DO
        ELSE
            WRITE(11,*) Thetamin*PI/180.
        END IF
    END IF
    CLOSE(11)
!   
!   --- Save index of cases ----------------------------------------------------------------------------------------------
!
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/index.dat')
    WRITE(10,*) Nw,Nbeta,Nradiation,Nintegration,Ntheta
    WRITE(10,*) '--- Force ---'
    DO k=1,Nintegration
        WRITE(10,*) k,IntCase(k)%Body,IntCase(k)%Mode
    END DO
    WRITE(10,*) '--- Motion ---'
    DO k=1,Nradiation
        WRITE(10,*) k,RadCase(k)%Body,RadCase(k)%Mode
    END DO 
    WRITE(10,*) (Beta(k),k=1,Nbeta)
    WRITE(10,*) (w(k),k=1,Nw)
    WRITE(10,*) ((Thetamin+(Thetamax-Thetamin)*(k-1)/(NTheta-1))*PI/180.,k=1,Ntheta)
    CLOSE(10)
!
!   --- Finalize ----------------------------------------------------------------------------------------------------
! 
    DEALLOCATE(RadCase,IntCase,Beta)    
!
    END PROGRAM Main
