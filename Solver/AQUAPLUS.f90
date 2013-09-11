! ----------------------------------------------
!
!   Solve the BVPs
!
! ----------------------------------------------
!
    PROGRAM Main
!   
    USE MIdentification
    USE MMesh
    USE MBodyConditions
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
!   Body conditions
    TYPE(TBodyConditions)   :: BodyConditions
!   Kochin function
    INTEGER                 :: Ntheta    
    REAL,DIMENSION(:),ALLOCATABLE :: Theta
!   Meshes
    TYPE(TMesh) :: Mesh 
    TYPE(TMesh) :: MeshFS
!   Aquaplus
    REAL                    :: T
    COMPLEX,DIMENSION(:),ALLOCATABLE :: NVEL,PRESSURE
    COMPLEX,DIMENSION(:),ALLOCATABLE :: HKochin
!   Results
    REAL :: XEFF,YEFF
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
!   Initialise Aquaplus
    CALL INITIALIZE(ID,NFA,NSYMY,XEFF,YEFF,Mesh)
    ALLOCATE(NVEL(NFA*2**NSYMY),PRESSURE(NFA*2**NSYMY))
    WRITE(*,'(A,$)') '.'
!   Initialise Kochin function calculation
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mesh/Kochin.dat')
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
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Mesh/FS.dat')
    READ(10,*) MeshFS%Npoints    
    IF (MeshFS%Npoints.GT.0) THEN
        CALL CreateTMesh(MeshFS,MeshFS%Npoints,1,1)
        DO j=1,MeshFS%Npoints
            READ(10,*) MeshFS%X(1,j),MeshFS%X(2,j)
        END DO
    END IF
    CLOSE(10)
!   Initialise results table
    ALLOCATE(Force(6*Mesh%Nbodies,Bodyconditions%Nproblems))
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
        DO c=1,Mesh%Npanels
            Force(1+6*(Mesh%cPanel(c)-1),j)=Force(1+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*Mesh%N(1,c)*Mesh%A(c)
            Force(2+6*(Mesh%cPanel(c)-1),j)=Force(2+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*Mesh%N(2,c)*Mesh%A(c)
            Force(3+6*(Mesh%cPanel(c)-1),j)=Force(3+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*Mesh%N(3,c)*Mesh%A(c)
            Force(4+6*(Mesh%cPanel(c)-1),j)=Force(4+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*(-(Mesh%XM(3,c)-0.)*Mesh%N(2,c)+(Mesh%XM(2,c)-YEFF)*Mesh%N(3,c))*Mesh%A(c)
            Force(5+6*(Mesh%cPanel(c)-1),j)=Force(5+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*(-(Mesh%XM(1,c)-XEFF)*Mesh%N(3,c)+(Mesh%XM(3,c)-YEFF)*Mesh%N(1,c))*Mesh%A(c)
            Force(6+6*(Mesh%cPanel(c)-1),j)=Force(6+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c)*(-(Mesh%XM(2,c)-YEFF)*Mesh%N(1,c)+(Mesh%XM(1,c)-YEFF)*Mesh%N(2,c))*Mesh%A(c)
        END DO
        IF (Mesh%Isym.EQ.1) THEN
            DO c=1,Mesh%Npanels
                Force(1+6*(Mesh%cPanel(c)-1),j)=Force(1+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c+Mesh%Npanels)*Mesh%N(1,c)*Mesh%A(c)
                Force(2+6*(Mesh%cPanel(c)-1),j)=Force(2+6*(Mesh%cPanel(c)-1),j)+PRESSURE(c+Mesh%Npanels)*Mesh%N(2,c)*Mesh%A(c)
                Force(3+6*(Mesh%cPanel(c)-1),j)=Force(3+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c+Mesh%Npanels)*Mesh%N(3,c)*Mesh%A(c)
                Force(4+6*(Mesh%cPanel(c)-1),j)=Force(4+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c+Mesh%Npanels)*((Mesh%XM(3,c)-0.)*Mesh%N(2,c)-(Mesh%XM(2,c)-YEFF)*Mesh%N(3,c))*Mesh%A(c)
                Force(5+6*(Mesh%cPanel(c)-1),j)=Force(5+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c+Mesh%Npanels)*(-(Mesh%XM(1,c)-XEFF)*Mesh%N(3,c)+(Mesh%XM(3,c)-YEFF)*Mesh%N(1,c))*Mesh%A(c)
                Force(6+6*(Mesh%cPanel(c)-1),j)=Force(6+6*(Mesh%cPanel(c)-1),j)-PRESSURE(c+Mesh%Npanels)*((Mesh%XM(2,c)-YEFF)*Mesh%N(1,c)-(Mesh%XM(1,c)-YEFF)*Mesh%N(2,c))*Mesh%A(c)
            END DO    
        END IF    
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
    DO c=1,6*Mesh%Nbodies
        DO j=1,BodyConditions%Nproblems
            IF (Bodyconditions%Switch_type(j).EQ.0) THEN
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