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
MODULE BodyConditions
!
IMPLICIT NONE
!
CONTAINS
!-- SUBROUTINE ComputeRadiationCondition
    SUBROUTINE ComputeRadiationCondition(Mesh,c,iCase,Direction,Axis,NVEL)  
    USE MMesh
    IMPLICIT NONE
    TYPE(TMesh) :: Mesh
    INTEGER :: c,iCase
    REAL,DIMENSION(3) :: Direction,Axis
    COMPLEX,DIMENSION(:) :: NVEL
    REAL,DIMENSION(3) :: VEL
    INTEGER :: i
    SELECT CASE (iCase)    
    CASE (1)
!       Degree of freedom is a translation      
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                VEL(1)=Direction(1)
                VEL(2)=Direction(2)
                VEL(3)=Direction(3)
                NVEL(i)=CMPLX(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
            ELSE
                NVEL(i)=CMPLX(0.,0.)
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(i).EQ.c) THEN
                    VEL(1)=Direction(1)
                    VEL(2)=Direction(2)
                    VEL(3)=Direction(3)
                    NVEL(i+Mesh%Npanels)=CMPLX(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
                ELSE
                    NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
                END IF
            END IF
        END DO
    CASE (2)
!       Degree of freedom is a rotation
        DO i=1,Mesh%Npanels
            IF (Mesh%cPanel(i).EQ.c) THEN
                VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(Mesh%XM(2,i)-Axis(2))
                VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
                VEL(3)=Direction(1)*(Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))                
                NVEL(i)=CMPLX(Mesh%N(1,i)*VEL(1)+Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
            ELSE
                NVEL(i)=CMPLX(0.,0.)
            END IF
            IF (Mesh%iSym.EQ.1) THEN
                IF (Mesh%cPanel(i).EQ.c) THEN
                    VEL(1)=Direction(2)*(Mesh%XM(3,i)-Axis(3))-Direction(3)*(-Mesh%XM(2,i)-Axis(2))
                    VEL(2)=Direction(3)*(Mesh%XM(1,i)-Axis(1))-Direction(1)*(Mesh%XM(3,i)-Axis(3))
                    VEL(3)=Direction(1)*(-Mesh%XM(2,i)-Axis(2))-Direction(2)*(Mesh%XM(1,i)-Axis(1))                
                    NVEL(i+Mesh%Npanels)=CMPLX(Mesh%N(1,i)*VEL(1)-Mesh%N(2,i)*VEL(2)+Mesh%N(3,i)*VEL(3),0.)
                ELSE
                    NVEL(i+Mesh%Npanels)=CMPLX(0.,0.)
                END IF
            END IF
        END DO
    CASE (3)
        WRITE(*,*) 'Error: radiation case 3 not implemented yet'
        STOP
    CASE DEFAULT
        WRITE(*,*) 'Error: unknown radiation case'
        STOP
    END SELECT
    END SUBROUTINE
!-- SUBROUTINE ComputeDiffractionCondition
    SUBROUTINE ComputeDiffractionCondition(Mesh,w,beta,Environment,PRESSURE,NVEL)
    USE MEnvironment
    USE MMEsh
    IMPLICIT NONE
!   Inputs/outputs
    TYPE(TMesh)             :: Mesh
    REAL                    :: w,beta           ! Wave period, direction and wavenumber
    TYPE(TEnvironment)      :: Environment      ! Environment
    COMPLEX,DIMENSION(*)    :: PRESSURE,NVEL    ! Pressure and normal velocities on panels
!   Locals
    REAL :: PI,X0,CIH,SIH
    REAL :: kwaveh,kwave
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    REAL :: wbar
    INTEGER :: i,j
    COMPLEX :: tmp
    COMPLEX :: Phi,p,Vx,Vy,Vz
    PI=4.*ATAN(1.) 
!   Compute wavenumber
!    kwave=w**2/Environment%g
!    IF ((Environment%Depth.GT.0.).AND.(kwave*Environment%Depth.LE.20)) THEN    
!    IF (Environment%Depth.GT.0.) THEN                                              
!        kwaveh=X0(kwave*Environment%Depth) 
!        kwave=kwaveh/Environment%Depth
!    END IF 
    kwave=Wavenumber(w,Environment)
!   Compute potential and normal velocities
    DO i=1,2**Mesh%Isym*Mesh%Npanels
        IF (i.LE.Mesh%Npanels) THEN
!            wbar=(Mesh%XM(1,i)-Environment%xeff)*COS(Beta)+(Mesh%XM(2,i)-Environment%YEFF)*SIN(Beta)
!            PRESSURE(i)=Environment%g/w*CEXP(II*kwave*wbar)
!            NVEL(i)=PRESSURE(i)*(kwave*(COS(Beta)*Mesh%N(1,i)+SIN(Beta)*Mesh%N(2,i))*CIH(kwave,Mesh%XM(3,i),Environment%Depth)-II*kwave*Mesh%N(3,i)*SIH(kwave,Mesh%XM(3,i),Environment%Depth))            
!            PRESSURE(i)=PRESSURE(i)*CIH(kwave,Mesh%XM(3,i),Environment%Depth)
            CALL Compute_Wave(kwave,w,beta,Mesh%XM(1,i),Mesh%XM(2,i),Mesh%XM(3,i),Phi,p,Vx,Vy,Vz,Environment) 
            PRESSURE(i)=p
            NVEL(i)=-(Vx*Mesh%N(1,i)+Vy*Mesh%N(2,i)+Vz*Mesh%N(3,i))    
        ELSE
!            wbar=(Mesh%XM(1,i-Mesh%Npanels)-Environment%XEFF)*COS(Beta)+(-Mesh%XM(2,i-Mesh%Npanels)-Environment%YEFF)*SIN(Beta)
!            PRESSURE(i)=-Environment%g/w*CEXP(II*kwave*wbar)
!            NVEL(i)=PRESSURE(i)*(II*kwave*(COS(Beta)*Mesh%N(1,i-Mesh%Npanels)+SIN(Beta)*(-1.*Mesh%N(2,i-Mesh%Npanels)))*CIH(kwave,Mesh%XM(3,i-Mesh%Npanels),Environment%Depth)+kwave*Mesh%N(3,i-Mesh%Npanels)*SIH(kwave,Mesh%XM(3,i-Mesh%Npanels),Environment%Depth)) 
!            PRESSURE(i)=PRESSURE(i)*CIH(kwave,Mesh%XM(3,i-Mesh%Npanels),Environment%Depth)
            CALL Compute_Wave(kwave,w,beta,Mesh%XM(1,i-Mesh%Npanels),-Mesh%XM(2,i-Mesh%Npanels),Mesh%XM(3,i-Mesh%Npanels),Phi,p,Vx,Vy,Vz,Environment) 
            PRESSURE(i)=p
            NVEL(i)=-(Vx*Mesh%N(1,i-Mesh%Npanels)-Vy*Mesh%N(2,i-Mesh%Npanels)+Vz*Mesh%N(3,i-Mesh%Npanels))  
        END IF        
!        PRESSURE(i)=Environment%RHO*w*II*PRESSURE(i)
   END DO
   END SUBROUTINE ComputeDiffractionCondition
!-- 
END MODULE
!   PROGRAMMES AQUADYN : (15-04-85) P.GUEVEL,J-C.DAUBISSE,G.DELHOMMEAU  C
!-- FUNCTION X0   
    REAL FUNCTION X0(AK)
    REAL::AK
    INTEGER::ITOUR,IITER    
    REAL::XS,XI,XM,PAS,EPS,VAL1,VAL2,VALM                                      
    EPS=5.E-6                                                                 
    ITOUR=0                                                                   
    XI=0.                                                                     
    XS=AK
    VAL1= AK-XI*TANH(XI)
    VAL2= AK-XS*TANH(XS) 
    DO WHILE ((ITOUR.LT.1000).AND.(VAL1*VAL2.GT.0)) 
        ITOUR=ITOUR+1 
        XI=XS
        XS=XS*2
        VAL1= AK-XI*TANH(XI)
        VAL2= AK-XS*TANH(XS) 
    END DO  
    IF (ITOUR.GE.1000) THEN
        WRITE(*,*) 'Error: unable to find the wavenumber'
        STOP
    END IF
    IITER=0
    DO WHILE ((IITER.LT.1000).AND.(ABS((XS-XI)/XM).GT.EPS))     
        IITER=IITER+1
        XM=(XI+XS)*0.5   
        VAL1= AK-XI*TANH(XI)
        VALM= AK-XM*TANH(XM)
        VAL2= AK-XS*TANH(XS)
        IF (VAL1*VALM.GE.0.) THEN
            XI=XM
        ELSE
            XS=XM
        END IF
    END DO  
    IF (IITER.GE.1000) THEN
        WRITE(*,*) 'Error: unable to find the wavenumber'
        STOP
    END IF  
    X0=XM 
    RETURN                                                                    
    END FUNCTION