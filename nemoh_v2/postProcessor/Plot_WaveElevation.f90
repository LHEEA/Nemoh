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
    SUBROUTINE Plot_WaveElevation(ID,Environment,iw,iBeta,RAOS,Results)
!    
    USE MIdentification
    USE MResults
    USE MEnvironment
!    
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TID) :: ID
    INTEGER :: iw,iBeta
    TYPE(TResults) :: Results
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,*) :: RAOs
    TYPE(TEnvironment) :: Environment
!   Locals
    CHARACTER*20 :: lookfor
    INTEGER :: Nx,Ny
    REAL :: Lx,Ly
    REAL,DIMENSION(:),ALLOCATABLE :: X,Y
    REAL :: r,theta
    COMPLEX,DIMENSION(:,:),ALLOCATABLE :: etaI,etaP,eta
    INTEGER :: j,i,k,M,N,d,l
    REAL :: PI
    REAL :: w,kwave,CIH
    COMPLEX,PARAMETER :: II=CMPLX(0.,1.)
    COMPLEX :: HKleft,HKright,HKochin,Potential,p,Vx,Vy,Vz
!
    PI=4.*ATAN(1.)
!   Read data
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/Nemoh.cal')
    READ(10,'(A20)') lookfor
    DO WHILE (lookfor.NE.'--- Post processing ')
        READ(10,'(A20)') lookfor
    END DO
    DO i=1,3
        READ(10,*)
    END DO
    READ(10,*) Nx,Ny,Lx,Ly
    CLOSE(10)
!   Calculate and save wave elevations
    ALLOCATE(X(Nx),Y(Ny),etaI(Nx,Ny),etaP(Nx,Ny),eta(Nx,Ny))
    DO i=1,Nx
        X(i)=-0.5*Lx+Lx*(i-1)/(Nx-1)
    END DO
    DO i=1,Ny
        Y(i)=-0.5*Ly+Ly*(i-1)/(Ny-1)
    END DO
    w=Results%w(iw)
    kwave=Wavenumber(w,Environment)
    DO i=1,Nx
        DO j=1,Ny
            r=SQRT((X(i)-Environment%XEFF)**2+(Y(j)-Environment%YEFF)**2)
            theta=ATAN2((Y(j)-Environment%YEFF),(X(i)-Environment%XEFF))
            k=1
            DO WHILE ((k.LT.Results%Ntheta-1).AND.(Results%theta(k+1).LT.theta))
                k=k+1
            END DO
            IF (k.EQ.Results%Ntheta) THEN
                WRITE(*,*) ' Error: range of theta in Kochin coefficients is too small'
                STOP
            END IF
            CALL Compute_Wave(kwave,w,Results%beta(iBeta),X(i),Y(j),0.,Potential,p,Vx,Vy,Vz,Environment)
            EtaI(i,j)=1./Environment%G*II*w*Potential
            HKleft=0.
            HKright=0.                
            DO l=1,Results%Nradiation
                HKleft=HKleft+RAOs(l,iw,iBeta)*Results%HKochinRadiation(iw,l,k)
                HKright=HKright+RAOs(l,iw,iBeta)*Results%HKochinRadiation(iw,l,k+1)
            END DO
            HKleft=HKleft+Results%HKochinDiffraction(iw,iBeta,k)
            HKright=HKright+Results%HKochinDiffraction(iw,iBeta,k+1)
            HKochin=HKleft+(HKright-HKleft)*(theta-Results%theta(k))/(Results%theta(k+1)-Results%theta(k))
            IF (r.GT.0) THEN
                Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
            ELSE
                Potential=0.
            END IF
            EtaP(i,j)=1./Environment%G*II*w*Potential 
            Eta(i,j)=EtaI(i,j)+EtaP(i,j)
        END DO
    END DO
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/WaveField.tec')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "etaI_C" "etaI_S" "etaP_C" "etaC_S" "etaI_C+etaP_C" "etaI_S+etaI_P" "|etaP|" "|etaI+etaP|"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - w =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Ny-1),' , F=FEPOINT,ET=QUADRILATERAL'
    DO i=1,Nx
        DO j=1,Ny
            WRITE(10,'(10(X,E14.7))') X(i),Y(j),REAL(etaI(i,j)),IMAG(etaI(i,j)),REAL(etaP(i,j)),IMAG(etaP(i,j)),REAL(etaI(i,j)+etaP(i,j)),IMAG(etaI(i,j)+etaP(i,j)),ABS(etaP(i,j)),ABS(eta(i,j))
        END DO
    END DO
    DO i=1,Nx-1
        DO j=1,Ny-1
            WRITE(10,'(I5,3(2X,I5))') j+(i-1)*Ny,j+i*Ny,j+1+i*Ny,j+1+(i-1)*Ny
		END DO  
    END DO     
    CLOSE(10)
    DEALLOCATE(X,Y,etaI,etaP,eta)
!
    END SUBROUTINE Plot_WaveElevation
!   
    SUBROUTINE Initialize_Plot_WaveElevation(Switch_Plot_WaveElevation,namefile)   
    IMPLICIT NONE	    
    CHARACTER*(*) :: namefile
    CHARACTER*20 :: lookfor
    CHARACTER*80 :: discard
    INTEGER :: i
    REAL :: Switch_Plot_WaveElevation
    OPEN(10,FILE=namefile)
    READ(10,'(A20)') lookfor
    DO WHILE (lookfor.NE.'--- Post processing ')
        READ(10,'(A20,A)') lookfor,discard
    END DO    
    DO i=1,3
        READ(10,*)
    END DO
    CLOSE(10)
    END SUBROUTINE  Initialize_Plot_WaveElevation
