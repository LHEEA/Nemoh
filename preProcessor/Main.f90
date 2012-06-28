! ----------------------------------------------
!
!   Definition of calculation cases
!
! ----------------------------------------------
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
    END TYPE
!   Radiation cases
    INTEGER :: Nradiation
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: RadCase
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL
!   Diffraction cases
    INTEGER :: Nbeta
    REAL :: betamin,betamax
    REAL,DIMENSION(:),ALLOCATABLE :: beta
    COMPLEX,DIMENSION(:),ALLOCATABLE :: Pressure
!   Force integration cases
    INTEGER :: Nintegration
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: IntCase
    REAL,DIMENSION(:),ALLOCATABLE :: NDS
    REAL,DIMENSION(:,:),ALLOCATABLE :: FNDS 
!   Froude Krylov forces
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: FKforce 
!   Free surface visualisation
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
!   Kochin function
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
    CALL ReadTEnvironment(Environment,ID%ID(1:ID%lID)//'/aquaplus.cal')
!   Read number of radiation cases, number of force integration cases, range of considered periods and diffraction cases
!   an directions for Kochin function
    Nradiation=0
    Nintegration=0
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/aquaplus.cal')
    DO c=1,7
        READ(10,*)
    END DO
    DO c=1,Mesh%Nbodies
        DO i=1,8
            READ(10,*)
        END DO
        READ(10,*) M
        Nradiation=Nradiation+M
        READ(10,*)
        READ(10,*) M
        Nintegration=Nintegration+M
        READ(10,*)
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
    READ(10,*)
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
    READ(10,*) Ntheta,thetamin,thetamax
    READ(10,*)
    READ(10,*) Nx,Ny,Lx,Ly
    CLOSE(10)
!   re-read input file and store radiation and integration cases
    ALLOCATE(RadCase(Nradiation))
    ALLOCATE(IntCase(Nintegration))
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/aquaplus.cal')
    DO c=1,7
        READ(10,*)
    END DO
    jrad=0
    jint=0
    DO c=1,Mesh%Nbodies
        DO i=1,8
            READ(10,*)
        END DO
        READ(10,*) M
        READ(10,*) (RadCase(jrad+i)%ICase,i=1,M)
        DO i=1,M
            RadCase(jrad+i)%Body=c
        END DO
        jrad=jrad+M
        READ(10,*) M
        READ(10,*) (IntCase(jint+i)%ICase,i=1,M)
        DO i=1,M
            IntCase(jint+i)%Body=c
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
    WRITE(*,'(A,I5,A)') '  ->',Nintegration,' integration cases'
    WRITE(*,*) ' '
!
!   --- Generate force integration file ----------------------------------------------------------------------------------------
!
    ALLOCATE(FNDS(Nintegration,Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(NDS(Mesh%Npanels*2**Mesh%Isym))
    DO j=1,Nintegration
        CALL ComputeNDS(Mesh,IntCase(j)%Body,IntCase(j)%Icase,NDS)
        DO c=1,Mesh%Npanels*2**Mesh%Isym
            FNDS(j,c)=NDS(c)
        END DO
    END DO 
    DEALLOCATE(NDS)
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/Integration.dat')
    WRITE(11,*) 0,Nintegration
    DO j=1,Nintegration
        WRITE(11,*) 0.,IntCase(j)%Body,IntCase(j)%Icase,(FNDS(j,c),c=1,Mesh%Npanels*2**Mesh%Isym)
    END DO
    CLOSE(11)
!
!   --- Generate radiations cases file ----------------------------------------------------------------------------------------
!
    ALLOCATE(NVEL(Mesh%Npanels*2**Mesh%Isym))
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/RadiationCases.dat')
    WRITE(11,*) Nw,Nradiation
    DO i=1,Nw
        DO j=1,Nradiation
            CALL ComputeRadiationCondition(Mesh,RadCase(j)%Body,RadCase(j)%Icase,NVEL)
            WRITE(11,*) 2.*PI/w(i),RadCase(j)%Body,RadCase(j)%Icase,(REAL(NVEL(c)),IMAG(NVEL(c)),c=1,Mesh%Npanels*2**Mesh%Isym)
        END DO
    END DO
    CLOSE(11)
!
!   --- Generate diffraction cases file ----------------------------------------------------------------------------------------
!
    ALLOCATE(PRESSURE(Mesh%Npanels*2**Mesh%Isym))
    ALLOCATE(FKForce(Nw,Nbeta,Nintegration))
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/DiffractionCases.dat')
    WRITE(11,*) Nw,Nbeta
    DO j=1,Nw
        DO i=1,Nbeta
            CALL ComputeWave(Mesh,w(j),Beta(i),Environment,PRESSURE,NVEL)
            WRITE(11,*) 2.*PI/w(j),Beta(i),0,(-REAL(NVEL(c)),-IMAG(NVEL(c)),c=1,Mesh%Npanels*2**Mesh%Isym)
!           Calculate the corresponding FK forces
            DO k=1,Nintegration
                FKForce(j,i,k)=0.
                DO c=1,Mesh%nPanels*2**Mesh%Isym
                    FKForce(j,i,k)=FKForce(j,i,k)+PRESSURE(c)*FNDS(k,c)    
                END DO
            END DO
        END DO    
    END DO
    CLOSE(11)
    DEALLOCATE(PRESSURE,NVEL,FNDS)
!
!   --- Save FK forces ----------------------------------------------------------------------------------------
!
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/FKForce.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    DO k=1,Nintegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',IntCase(k)%Body,IntCase(k)%Icase,')" "angle(F',IntCase(k)%Body,IntCase(k)%Icase,')"'
    END DO
    DO c=1,Nbeta
        WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="FKforce - beta = ',beta(c)*180./PI,'",I=',Nw,',F=POINT'
        DO i=1,Nw
            WRITE(10,'(80(X,E14.7))') w(i),(ABS(FKForce(i,c,k)),ATAN2(IMAG(FKForce(i,c,k)),REAL(FKForce(i,c,k))),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
!
!   --- Generate Free Surface visualisation file ----------------------------------------------------------------------
!
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/FS.dat')
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
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/Kochin.dat')
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
!   --- Finalize ----------------------------------------------------------------------------------------
! 
    DEALLOCATE(FKForce)
!
    END PROGRAM Main