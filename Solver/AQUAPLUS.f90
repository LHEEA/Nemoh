! ----------------------------------------------
!
!   Solve the BVPs
!
! ----------------------------------------------
!
    PROGRAM Main
!   
    USE MIdentification
    USE MResults
    USE iflport
    USE SOLVE_BEM
    USE OUTPUT
    USE INITIALIZATION
!    
    IMPLICIT NONE
!   ID
    TYPE(TID)               :: ID
!   Array sizes
    INTEGER                 :: NFA,NSYMY    ! Number of panels and symmetry about xOz plane
!   Wave periods
    INTEGER :: N
    REAL,DIMENSION(:),ALLOCATABLE :: Period
!   TYPE TCase
    TYPE TCase
        INTEGER :: Body
        INTEGER :: ICase
    END TYPE
!   Kochin function
    INTEGER                 :: Ntheta    
    REAL,DIMENSION(:),ALLOCATABLE :: Theta
!   Free surface visualisation
    TYPE(TMeshFS) :: MeshFS
!   Radiation cases
    INTEGER :: Nradiation
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: RadCase
    TYPE(TResults) :: RadiationResults
!   Diffraction cases
    INTEGER :: Nbeta
    REAL,DIMENSION(:),ALLOCATABLE :: beta
    TYPE(TResults) :: DiffractionResults
    REAL,DIMENSION(:,:,:),ALLOCATABLE :: AFKForce,PFKForce
    COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: ExcitationForce
!   Force integration cases
    INTEGER :: Nintegration
    TYPE(TCase),DIMENSION(:),ALLOCATABLE :: IntCase
     REAL,DIMENSION(:,:),ALLOCATABLE :: FNDS 
!   Aquaplus
    REAL                    :: T
    REAL,DIMENSION(:),ALLOCATABLE    :: RVEL,IVEL
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL,PRESSURE
    COMPLEX,DIMENSION(:),ALLOCATABLE :: HKochin
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
!   Initialise Aquaplus
    CALL INITIALIZE(ID,NFA,NSYMY)
    WRITE(*,'(A,$)') '.'
!   Initialise radiation case studies
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/RadiationCases.dat')
    READ(10,*) N,Nradiation
    ALLOCATE(Period(N),RadCase(NRadiation))
    CLOSE(10)
    ALLOCATE(RVEL(2**NSYMY*NFA),IVEL(2**NSYMY*NFA))
    ALLOCATE(NVEL(2**NSYMY*NFA),PRESSURE(2**NSYMY*NFA))
    WRITE(*,'(A,$)') '.'
!   Initialise diffraction case studies
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/DiffractionCases.dat')
    READ(10,*) M,Nbeta
    IF (M.NE.N) THEN
        WRITE(*,*) 'Error: number of periods for diffraction is different from for radiation'
    END IF
    ALLOCATE(beta(Nbeta))
    CLOSE(10)
    WRITE(*,'(A,$)') '.'
!   Initialise integration
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Integration.dat')
    READ(10,*) M,Nintegration
    ALLOCATE(IntCase(Nintegration),FNDS(Nintegration,NFA*2**NSYMY))
    DO j=1,Nintegration
        READ(10,*) Discard,IntCase(j)%Body,IntCase(j)%Icase,(FNDS(j,c),c=1,NFA*2**NSYMY)
    END DO
    CLOSE(10)
!   Initialise Kochin function calculation
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Kochin.dat')
    READ(10,*) Ntheta 
    IF (Ntheta.GT.0) THEN
        ALLOCATE(Theta(Ntheta))
        DO j=1,Ntheta
            READ(10,*) Theta(j)
        END DO
        ALLOCATE(HKochin(NTheta))
    END IF
    CLOSE(10)
!   Initialise free surface calculation points
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/FS.dat')
    READ(10,*) MeshFS%Npoints,MeshFS%Npanels
    IF (MeshFS%Npoints.GT.0) THEN
        CALL CreateTMeshFS(MeshFS,MeshFS%Npoints,MeshFS%Npanels)
        DO j=1,MeshFS%Npoints
            READ(10,*) MeshFS%XC(j),MeshFS%YC(j),MeshFS%ZC(j)
        END DO
        DO j=1,MeshFS%Npanels
            READ(10,*) (MeshFS%P(k,j),k=1,4)
        END DO
    END IF
    CLOSE(10)
    CALL CreateTResults(RadiationResults,N,Nradiation,Nintegration,Ntheta)
    CALL CreateTResults(DiffractionResults,N,Nbeta,Nintegration,Ntheta)   
    WRITE(*,*) '. Done !'
    WRITE(*,*) ' '
!
!   --- Solve BVPs and calculate forces -------------------------------------------------------------------------------------------------------------------------------------------------
!
    WRITE(*,*) ' -> Solve BVPs and calculate forces ' 
    WRITE(*,*) ' '
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/RadiationCases.dat')
    READ(10,*) M,M
    OPEN(12,FILE=ID%ID(1:ID%lID)//'/DiffractionCases.dat')
    READ(12,*) M,M
    DO j=1,N 
!       --> Radiation cases
        DO i=1,Nradiation            
!           Read Normal velocities        
            READ(10,*) Period(j),RadCase(i)%Body,RadCase(i)%Icase,(RVEL(c),IVEL(c),c=1,NFA*2**NSYMY)
            RadiationResults%Period(j)=Period(j)  
            RadiationResults%Rcase(i)=10.*RadCase(i)%Body+1.*RadCase(i)%Icase
            IF (i.EQ.1) THEN
                WRITE(*,'(A,F8.4,A,$)') ' Calculation for period = ',Period(j),' s '
            ELSE
                WRITE(*,'(A,$)') '.'
            END IF
            DO c=1,NFA*2**NSYMY
                NVEL(c)=CMPLX(RVEL(c),IVEL(c))
            END DO
!           Solve BVP
            CALL SOLVE_BVP(ID,Period(j),NVEL,PRESSURE,NTheta,Theta,HKochin,MeshFS)
!           Calculate force coefficient
            DO c=1,Nintegration
                RadiationResults%Iintegration(c)=10.*IntCase(c)%Body+1.*IntCase(c)%Icase
                RadiationResults%Force(j,i,c)=0.
                DO d=1,2**NSYMY*NFA
                    RadiationResults%Force(j,i,c)=RadiationResults%Force(j,i,c)+PRESSURE(d)*FNDS(c,d)
                END DO
            END DO
!           Retain Kochin coefficients
            DO c=1,Ntheta
                RadiationResults%Theta(c)=Theta(c)
                RadiationResults%HKochin(j,i,c)=HKochin(c)
            END DO            
        END DO
 !       --> Diffraction cases
        DO i=1,Nbeta
!           Read Normal velocities        
            READ(12,*) Period(j),Beta(i),Discard,(RVEL(c),IVEL(c),c=1,NFA*2**NSYMY)
            DiffractionResults%Period(j)=Period(j)
            DiffractionResults%Rcase(i)=Beta(i)
            IF ((i.EQ.1).AND.(Nradiation.EQ.0)) THEN
                WRITE(*,'(A,F8.4,A,$)') ' Calculation for period = ',Period(j),' s '
            ELSE
                WRITE(*,'(A,$)') '.'
            END IF
            DO c=1,NFA*2**NSYMY
                NVEL(c)=CMPLX(RVEL(c),IVEL(c))
            END DO
!           Solve BVP
            CALL SOLVE_BVP(ID,Period(j),NVEL,PRESSURE,NTheta,Theta,HKochin,MeshFS)
!           Calculate force coefficient
            DO c=1,Nintegration
                DiffractionResults%Iintegration(c)=10.*IntCase(c)%Body+1.*IntCase(c)%Icase
                DiffractionResults%Force(j,i,c)=0.
                DO d=1,2**NSYMY*NFA
                    DiffractionResults%Force(j,i,c)=DiffractionResults%Force(j,i,c)+PRESSURE(d)*FNDS(c,d)
                END DO
            END DO
!           Retain Kochin coefficients
            DO c=1,Ntheta
                DiffractionResults%Theta(c)=Theta(c)
                DiffractionResults%HKochin(j,i,c)=HKochin(c)
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
    CALL SaveTResults(RadiationResults,ID%ID(1:ID%lID)//'/results/Radiation.dat')
    CALL SaveTResults(DiffractionResults,ID%ID(1:ID%lID)//'/results/Diffraction.dat')
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/RadiationCoefficients.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    DO k=1,Nintegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"A',IntCase(k)%Body,IntCase(k)%Icase,'" "B',IntCase(k)%Body,IntCase(k)%Icase,'"'
    END DO
    DO i=1,Nradiation
        WRITE(10,'(A,I4,A,I4,A,I6,A)') 'Zone t="Motion of body ',RadCase(i)%Body,' in DoF',RadCase(i)%Icase,'",I=',N,',F=POINT'
        DO j=1,N
            WRITE(10,'(80(X,E14.7))') 2.*PI/Period(j),(IMAG(RadiationResults%Force(j,i,k))*Period(j)/(2.*PI),-REAL(RadiationResults%Force(j,i,k)),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/RadiationKochin.tec')
    WRITE(10,'(A)') 'VARIABLES="Theta (rad)" "ABS(HKochin)" "ANGLE(HKochin)"'
    DO i=1,Nradiation
        DO j=1,N
            WRITE(10,'(A,I4,A,F7.4,A,I6,A)') 'Zone t="Motion of body ',RadCase(i)%Body,' at frequency ',2.*PI/Period(j),' rad/s",I=',Ntheta,',F=POINT'
            DO k=1,Ntheta
                WRITE(10,'(80(X,E14.7))') Theta(k),ABS(RadiationResults%HKochin(j,i,k)),ATAN2(IMAG(RadiationResults%HKochin(j,i,k)),REAL(RadiationResults%HKochin(j,i,k)))
            END DO
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/DiffractionForce.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    DO k=1,Nintegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',IntCase(k)%Body,IntCase(k)%Icase,')" "angle(F',IntCase(k)%Body,IntCase(k)%Icase,')"'
    END DO
    DO i=1,Nbeta
        WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="Diffraction force - beta = ',beta(i)*180./PI,'",I=',N,',F=POINT'
        DO j=1,N
            WRITE(10,'(80(X,E14.7))') 2.*PI/Period(j),(ABS(DiffractionResults%Force(j,i,k)),ATAN2(IMAG(DiffractionResults%Force(j,i,k)),REAL(DiffractionResults%Force(j,i,k))),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/DiffractionKochin.tec')
    WRITE(10,'(A)') 'VARIABLES="Theta (rad)" "ABS(HKochin)" "ANGLE(HKochin)"'
    DO i=1,Nbeta
        DO j=1,N
            WRITE(10,'(A,F7.3,A,F7.4,A,I6,A)') 'Zone t="Diffraction - beta = ',beta(i)*180./PI,' at frequency ',2.*PI/Period(j),' rad/s",I=',Ntheta,',F=POINT'
            DO k=1,Ntheta
                WRITE(10,'(80(X,E14.7))') Theta(k),ABS(DiffractionResults%HKochin(j,i,k)),ATAN2(IMAG(DiffractionResults%HKochin(j,i,k)),REAL(DiffractionResults%HKochin(j,i,k)))
            END DO
        END DO
    END DO
    CLOSE(10)
    ALLOCATE(ExcitationForce(N,Nbeta,Nintegration),AFKForce(N,Nbeta,Nintegration),PFKForce(N,Nbeta,Nintegration))
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/FKForce.tec')
    READ(10,*)
    DO k=1,Nintegration
        READ(10,*)
    END DO
    DO i=1,Nbeta
        READ(10,*) 
        DO j=1,N
            READ(10,*) Discard,(AFKForce(j,i,k),PFKForce(j,i,k),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/ExcitationForce.tec')
    WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
    DO k=1,Nintegration
        WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',IntCase(k)%Body,IntCase(k)%Icase,')" "angle(F',IntCase(k)%Body,IntCase(k)%Icase,')"'
    END DO
    DO i=1,Nbeta
        WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="Excitation force - beta = ',beta(i)*180./PI,'",I=',N,',F=POINT'
        DO j=1,N
            DO k=1,Nintegration
                ExcitationForce(j,i,k)=DiffractionResults%Force(j,i,k)+AFKForce(j,i,k)*(COS(PFKForce(j,i,k))+II*SIN(PFKForce(j,i,k)))
            END DO
            WRITE(10,'(80(X,E14.7))') 2.*PI/Period(j),(ABS(ExcitationForce(j,i,k)),ATAN2(IMAG(ExcitationForce(j,i,k)),REAL(ExcitationForce(j,i,k))),k=1,Nintegration)
        END DO
    END DO
    CLOSE(10)
    
!
!   --- Finalize -------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!
    DEALLOCATE(NVEL,PRESSURE,RadCase,IntCase,Period,FNDS,ExcitationForce,AFKForce,PFKForce)    
    CALL DeleteTResults(RadiationResults)
    CALL DeleteTResults(DiffractionResults)
    IF (Ntheta.GT.0) DEALLOCATE(Theta,HKochin)
    IF (MeshFS%Npoints.GT.0) CALL DeleteTMeshFS(MeshFS)
    CALL DEALLOCATE_DATA
!
    END PROGRAM Main
!----------------------------------------------------------------

    SUBROUTINE CREK(NPINTE)
    USE COM_VAR
    IMPLICIT NONE
    INTEGER:: I,J,JZ,IR,NPINTE
    REAL:: AKZ,AKR
    JZ=46
    IR=328
    DO 8006 J=1,JZ
      AKZ=AMIN1(10**(J/5.-6),10**(J/8.-4.5),16.)
      XZ(J)=-AKZ     
8006 CONTINUE
    XR(1)=0.
    DO 8007 I=2,IR
      IF(I.LT.40)THEN
	AKR=AMIN1(10**((I-1.)/5-6),4./3.+ABS((I-32.)/3.))
      ELSE
	AKR=4./3.+ABS((I-32.)/3.)
      ENDIF
      XR(I)=AKR
8007 CONTINUE

    DO 8009 J=1,JZ
      DO 8008 I=1,IR
	CALL VNS(NPINTE,XZ(J),XR(I),I,J)                                           
8008 CONTINUE
8009 CONTINUE

    RETURN

    END SUBROUTINE
!---------------------------------------------------------                                                                    

      SUBROUTINE VNS(NPINTE,AKZ,AKR,I,J) 
      USE COM_VAR
      USE ELEMENTARY_FNS
    IMPLICIT NONE         
      INTEGER:: I,J,NPINTE,IT
      REAL:: AKZ,AKR,PI,CT,ST,CQIT,TETA
      REAL:: QQT(NPINTE),CQT(NPINTE)
      REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
      COMPLEX:: IM,C1,C2,ZIK,GZ,CEX            
      IM=(0.,1.)                                                                
      PI=4.*ATAN(1.)
      CALL COFINT(NPINTE,CQT,QQT)         
      FD1JX=0.                                                              
      FD1JZ=0.                                                              
      FD2JX=0.                                                              
      FD2JZ=0.                                                              
      DO 30 IT=1,NPINTE                                           
	TETA=QQT(IT)                                                           
	CQIT=CQT(IT)                                                   
	CT=COS(TETA)                                                              
	ST=SIN(TETA)                                                              
	ZIK=AKZ+IM*AKR*CT    
	IF(REAL(ZIK)+30.)2,2,1
      2 CEX=(0.,0.)
	GOTO 3 
      1 CEX=CEXP(ZIK)                                               
      3 GZ=GG(ZIK,CEX)                      
	C1=CQIT*(GZ-1./ZIK)
	C2=CQIT*CEX
	FD1JX=FD1JX+CT*AIMAG(C1)                                 
	FD1JZ=FD1JZ+REAL(C1)                                     
	FD2JX=FD2JX+CT*AIMAG(C2)                                 
	FD2JZ=FD2JZ+REAL(C2)                                     
   30 CONTINUE                                                                
      APD1X(I,J)=FD1JX                                               
      APD1Z(I,J)=FD1JZ                                               
      APD2X(I,J)=FD2JX                                               
      APD2Z(I,J)=FD2JZ                                              
      RETURN 
                                                                   
      END SUBROUTINE  
!-------------------------------------------------------------------

      SUBROUTINE COFINT(NPINTE,CQT,QQT)  
      USE COM_VAR
    IMPLICIT NONE 
      INTEGER :: J,NPINTE
      REAL:: PI,QQT(NPINTE),CQT(NPINTE)

      PI=4.*ATAN(1.)
      DO 160 J=1,NPINTE   
      QQT(J)=-PI/2.+(J-1.)/(NPINTE-1.)*PI
      IF(J-1)161,161,162
  161 CQT(J)=PI/(3.*(NPINTE-1.))
      GOTO 160
  162 IF(J-NPINTE)163,161,161
  163 IF(MOD(J,2))164,165,164
  164 CQT(J)=2./(3.*(NPINTE-1.))*PI
      GOTO 160
  165 CQT(J)=4./(3.*(NPINTE-1.))*PI
  160 CONTINUE

      RETURN                                                                    
      END SUBROUTINE   
!------------------------------------------------------------------