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
MODULE MEnvironment
! Definition of TYPE Environment
TYPE TEnvironment
    REAL :: RHO         ! Sea water density
    REAL :: G           ! Gravity constant
    REAL :: Depth       ! Water depth
    REAL :: XEFF,YEFF   ! Coordinates of point where the incident wave is measured
END TYPE TEnvironment
!
CONTAINS
!   Read Environment data from file 
    SUBROUTINE ReadTEnvironment(Environment,file) 
    IMPLICIT NONE  
    TYPE(TEnvironment) :: Environment
    CHARACTER*(*) :: file
    OPEN(10,FILE=file)
    READ(10,*)
    READ(10,*) Environment%RHO
    READ(10,*) Environment%G
    READ(10,*) Environment%Depth
    READ(10,*) Environment%XEFF,Environment%YEFF
    CLOSE(10)    
    END SUBROUTINE
!   Calculate wave number for frequency w and given depth
    REAL FUNCTION Wavenumber(w,Environment)
    IMPLICIT NONE
    REAL :: w
    TYPE(TEnvironment) :: Environment
    REAL :: k0,x0,xg,xd,xc
    INTEGER,PARAMETER :: Nitemx=10000
    INTEGER :: Nite
    Wavenumber=w*w/Environment%g
    x0=Wavenumber*Environment%Depth
    IF ((x0.LE.20.).AND.(x0.GT.0.)) THEN
        xg=0.
        xd=x0
        Nite=0
        DO WHILE ((Nite.LT.Nitemx).AND.((x0-xg*TANH(xg))*(x0-xd*TANH(xd)).GT.0.))
            xg=xd
            xd=2.*xd
            Nite=Nite+1
        END DO
        Nite=0
        IF (Nite.GE.Nitemx) THEN
            WRITE(*,*) 'Error: unable to find the wavenumber'
            STOP
        END IF
        xc=0.5*(xd+xg)
        DO WHILE ((Nite.LT.Nitemx).AND.(ABS(xd-xg)/ABS(xc).GE.1.0E-06))
            xc=0.5*(xd+xg)
            IF ((x0-xg*TANH(xg))*(x0-xc*TANH(xc)).GT.0.) THEN
                xg=xc
            ELSE
                xd=xc
            END IF
            Nite=Nite+1
        END DO
        IF (Nite.GE.Nitemx) THEN
            WRITE(*,*) 'Error: unable to find the wavenumber'
            STOP
        END IF  
        Wavenumber=xc/Environment%Depth    
    END IF
    END FUNCTION Wavenumber
!   Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)
    SUBROUTINE Compute_Wave(k,w,beta,x,y,z,Phi,p,Vx,Vy,Vz,Environment)
    IMPLICIT NONE
    REAL :: k,w,beta,x,y,z
    COMPLEX :: Phi,p,Vx,Vy,Vz
    TYPE(TEnvironment) :: Environment
    REAL :: wbar,CIH,SIH
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    wbar=(x-Environment%XEFF)*COS(Beta)+(y-Environment%YEFF)*SIN(Beta)
    Phi=-II*Environment%g/w*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    p=Environment%rho*Environment%g*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vx=Environment%g/w*k*COS(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vy=Environment%g/w*k*SIN(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vz=-II*Environment%g/w*k*SIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    END SUBROUTINE Compute_wave

END MODULE 

!   COSH(k*(z+h))/COSH(kh)
    REAL FUNCTION CIH(k,z,h)
    IMPLICIT NONE
    REAL k,z,h
    IF ((k*h.LE.20).AND.(k*h.GT.0.)) THEN
        CIH=COSH(k*(z+h))/COSH(k*h)
    ELSE
        CIH=EXP(k*z)
    ENDIF
    END FUNCTION CIH
!   SINH(k*(z+h))/COSH(kh)
    REAL FUNCTION SIH(k,z,h)
    IMPLICIT NONE
    REAL k,z,h
    IF ((k*h.LE.20).AND.(k*h.GT.0.)) THEN
        SIH=SINH(k*(z+h))/COSH(k*h)
    ELSE
        SIH=EXP(k*z)
    ENDIF
    END FUNCTION SIH
