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
!   - J. Singh
!   - A. Babarit 
!
!--------------------------------------------------------------------------------------
MODULE SOLVE_BEM
!
IMPLICIT NONE
!
CONTAINS
!
    SUBROUTINE SOLVE_BVP(ProblemNumber,ID,Period,NVEL,PRESSURE,Switch_Kochin,NTheta,Theta,HKochin,Switch_FS,MeshFS,Switch_potential)
!
        USE MIDENTIFICATION
        USE COM_VAR
        USE MMesh
        USE OUTPUT
        USE SOLVE_BEM_INFD_DIRECT
        USE SOLVE_BEM_FD_DIRECT
!
        IMPLICIT NONE
!       Inputs/outputs
        INTEGER :: ProblemNumber
        TYPE(TID) :: ID
        REAL :: Period
        COMPLEX,DIMENSION(*) :: NVEL,PRESSURE        
        INTEGER :: Switch_Kochin
        INTEGER :: NTheta
        REAL,DIMENSION(*) :: THeta
        COMPLEX,DIMENSION(*) :: HKochin
        INTEGER :: Switch_FS
        TYPE(TMesh) :: MeshFS
        INTEGER :: Switch_potential
!       For solver (DIRECT or GMRES)
        INTEGER :: NEXP
!       Locals
        INTEGER :: i,j
        REAL :: PI,W
        COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
        REAL :: kwave
        REAL :: AKH,AMH
        COMPLEX,DIMENSION(1+MeshFS%Npoints) :: PHI,ETA
!
        PI=4.*ATAN(1.)
!       Define period
        T=Period
        W=2.*PI/T
        kwave=w*w/g
!       Select calculation in function of water depth
        IF ((Depth.EQ.0.).OR.(kwave*Depth.GE.20)) THEN
!           Calculate wave number
            kwave=w*w/g
            AMH=kwave*Depth
!           Solve with direct method ?
            IF (Indiq_solver .eq. 0) CALL SOLVE_POTENTIAL_INFD_DIRECT(NVEL)        
!           Solve using GMRES ?
            IF (Indiq_solver .eq. 1) THEN
                WRITE(*,*) ' Iterative solver is not available'
                STOP
             END IF
        ELSE
!           Calculate wave number
            AKH=w*w*Depth/G                                                 
	        AMH=X0(AKH)                                                                                                   
	        kwave=AMH/Depth 
!           Solve with direct method ?
            IF (Indiq_solver .eq. 0) CALL SOLVE_POTENTIAL_FD_DIRECT(NVEL,AMH,NEXP)        
!           Solve using GMRES ?
            IF (Indiq_solver .eq. 1) THEN
                WRITE(*,*) ' Iterative solver is not available'
                STOP
            END IF        
        END IF
!       Assemble pressure
        DO i=1,NFA
            PRESSURE(i)=RHO*II*W*ZPB(i) !*AIRE(i)
            IF (NSYMY.EQ.1) THEN
                PRESSURE(i+NFA)=RHO*II*W*ZPS(i) !*AIRE(i)    
            END IF
        END DO
!       Compute Kochin Functions
        IF (Switch_Kochin.EQ.1) THEN
            DO j=1,NTheta
                CALL COMPUTE_KOCHIN(kwave,Theta(j),HKochin(j))
            END DO
            CALL WRITE_KOCHIN(ID,ProblemNumber,HKochin,NTheta,Theta)
        END IF
!       Save free surface elevation
        IF (Switch_FS.EQ.1) THEN
            DO j=1,MeshFS%Npoints
                CALL COMPUTE_POTENTIAL_DOMAIN(ID,PHI(j),MeshFS%X(1,j),MeshFS%X(2,j),0.,kwave,AMH,NEXP)
                ETA(j)=II*W/G*PHI(j)
            END DO
            CALL WRITE_FS(ID,ProblemNumber,ETA,MeshFS)
        END IF
!       Save output
        IF (Switch_Potential.EQ.1) THEN
            CALL WRITE_POTENTIAL(ID,ProblemNumber,PRESSURE)
        END IF    
!    
    END SUBROUTINE SOLVE_BVP
!----------------------------------------------------------------
!    SUBROUTINE CREK(NPINTE)
!
!        USE COM_VAR
!        IMPLICIT NONE
!        INTEGER:: I,J,JZ,IR,NPINTE
!        REAL:: AKZ,AKR
!        JZ=46
!        IR=328
!        DO J=1,JZ
!            AKZ=AMIN1(10**(J/5.-6),10**(J/8.-4.5),16.)
!            XZ(J)=-AKZ     
!        END DO
!        XR(1)=0.
!        DO I=2,IR
!            IF(I.LT.40)THEN
!	            AKR=AMIN1(10**((I-1.)/5-6),4./3.+ABS((I-32.)/3.))
!            ELSE
!	            AKR=4./3.+ABS((I-32.)/3.)
!            ENDIF
!            XR(I)=AKR
!        END DO
!        DO J=1,JZ
!            DO I=1,IR
!	            CALL VNS(NPINTE,XZ(J),XR(I),I,J)                                           
!            END DO
!        END DO
!
!        RETURN
!
!    END SUBROUTINE
!---------------------------------------------------------                                                                    
!    SUBROUTINE VNS(NPINTE,AKZ,AKR,I,J)    
!
!        USE COM_VAR
!        USE ELEMENTARY_FNS       
!   
!        INTEGER:: I,J,NPINTE,IT
!        REAL:: AKZ,AKR,PI,CT,ST,CQIT,TETA
!        REAL:: QQT(NPINTE),CQT(NPINTE)
!        REAL:: FD1JX,FD1JZ,FD2JX,FD2JZ
!        COMPLEX:: IM,C1,C2,ZIK,GZ,CEX     
!        IM=(0.,1.)                                                                
!        PI=4.*ATAN(1.)
!        CALL COFINT(NPINTE,CQT,QQT)         
!        FD1JX=0.                                                              
!        FD1JZ=0.                                                              
!        FD2JX=0.                                                              
!        FD2JZ=0.                                                              
!        DO IT=1,NPINTE                                           
!	        TETA=QQT(IT)                                                           !
!	        CQIT=CQT(IT)                                                   
!	        CT=COS(TETA)                                                              
!	        ST=SIN(TETA)                                                              
!	        ZIK=AKZ+IM*AKR*CT    
!	        IF(REAL(ZIK)+30.)2,2,1
!              2 CEX=(0.,0.)
!	        GOTO 3 
!              1 CEX=CEXP(ZIK)                                               
!              3 GZ=GG(ZIK,CEX)                      
!	        C1=CQIT*(GZ-1./ZIK)
!	        C2=CQIT*CEX
!	        FD1JX=FD1JX+CT*AIMAG(C1)                                 
!	        FD1JZ=FD1JZ+REAL(C1)                                     
!	        FD2JX=FD2JX+CT*AIMAG(C2)                                 
!	        FD2JZ=FD2JZ+REAL(C2)                                     
!        END DO                                                                
!        APD1X(I,J)=FD1JX                                               
!        APD1Z(I,J)=FD1JZ                                               
!        APD2X(I,J)=FD2JX                                               
!        APD2Z(I,J)=FD2JZ                                              
!        RETURN 
!!                                                                   
!    END SUBROUTINE  
!-------------------------------------------------------------------
!    SUBROUTINE COFINT(NPINTE,CQT,QQT)   
!
!    USE COM_VAR
!
!    INTEGER :: J,NPINTE
!    REAL:: PI,QQT(NPINTE),CQT(NPINTE)
!    PI=4.*ATAN(1.)
!      DO 160 J=1,NPINTE   
!      QQT(J)=-PI/2.+(J-1.)/(NPINTE-1.)*PI
!      IF(J-1)161,161,162
!  161 CQT(J)=PI/(3.*(NPINTE-1.))
!      GOTO 160
!  162 IF(J-NPINTE)163,161,161
!  163 IF(MOD(J,2))164,165,164
!  164 CQT(J)=2./(3.*(NPINTE-1.))*PI
!      GOTO 160
!  165 CQT(J)=4./(3.*(NPINTE-1.))*PI
!  160 CONTINUE
!    RETURN                                                                    
!
!    END SUBROUTINE   
!------------------------------------------------------------------
END MODULE