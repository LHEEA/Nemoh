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

  USE Elementary_functions, ONLY: CIH, SIH

  IMPLICIT NONE  

  ! Definition of TYPE Environment
  TYPE TEnvironment
    REAL :: RHO         ! Sea water density
    REAL :: G           ! Gravity constant
    REAL :: Depth       ! Water depth
    REAL :: Xeff, Yeff  ! Coordinates of point where the incident wave is measured
  END TYPE TEnvironment

CONTAINS

  SUBROUTINE ReadTEnvironment(Environment, file) 
    ! Read Environment data from file 

    CHARACTER(LEN=*),   INTENT(IN)  :: file
    TYPE(TEnvironment), INTENT(OUT) :: Environment

    INTEGER :: u

    OPEN(NEWUNIT=u, FILE=file, FORM='FORMATTED', STATUS='OLD')
    READ(u,*)
    READ(u,*) Environment%RHO
    READ(u,*) Environment%G
    READ(u,*) Environment%Depth
    READ(u,*) Environment%Xeff, Environment%Yeff
    CLOSE(u)    

  END SUBROUTINE


  REAL FUNCTION Wavenumber(w, Environment)
    ! Calculate wave number for frequency w and given depth
    ! To be merge with X0 on Elementary_functions.

    REAL,               INTENT(IN) :: w
    TYPE(TEnvironment), INTENT(IN) :: Environment

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

  SUBROUTINE Compute_Wave(k,w,beta,x,y,z,Phi,p,Vx,Vy,Vz,Environment)
    ! Calculate the complex potential, pressure and fluid velocities for a regular wave eta=sin(k*wbar-wt)

    REAL :: k,w,beta,x,y,z
    COMPLEX :: Phi,p,Vx,Vy,Vz
    TYPE(TEnvironment) :: Environment
    REAL :: wbar
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)

    wbar=(x-Environment%XEFF)*COS(Beta)+(y-Environment%YEFF)*SIN(Beta)
    Phi=-II*Environment%g/w*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    p=Environment%rho*Environment%g*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vx=Environment%g/w*k*COS(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vy=Environment%g/w*k*SIN(beta)*CIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
    Vz=-II*Environment%g/w*k*SIH(k,z,Environment%Depth)*CEXP(II*k*wbar)
  END SUBROUTINE Compute_wave

END MODULE 
