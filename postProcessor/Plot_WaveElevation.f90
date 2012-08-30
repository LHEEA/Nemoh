    SUBROUTINE Plot_WaveElevation(ID,Environment,iPeriod,iBeta,RAOS,RadiationResults,DiffractionResults)
!    
    USE MIdentification
    USE MResults
    USE MEnvironment
!    
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TID) :: ID
    INTEGER :: iPeriod,iBeta
    TYPE(TResults) :: RadiationResults,DiffractionResults
    COMPLEX,DIMENSION(RadiationResults%Nintegration,DiffractionResults%Nperiod,*) :: RAOs
    TYPE(TEnvironment) :: Environment
!   Locals
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
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/aquaplus.cal')
    DO i=1,6
        READ(10,*)
    END DO
    READ(10,*) M
    DO i=1,M
        DO j=1,8
            READ(10,*)
        END DO
        READ(10,*) N
        DO j=1,N
            READ(10,*)
        END DO
        READ(10,*) N
        DO j=1,N
            READ(10,*)
        END DO
        READ(10,*) N
        DO j=1,N
            READ(10,*)
        END DO
    END DO
    DO i=1,7
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
    w=2.*PI/DiffractionResults%Period(iPeriod)
    kwave=Wavenumber(w,Environment)
    DO i=1,Nx
        DO j=1,Ny
            r=SQRT((X(i)-Environment%XEFF)**2+(Y(j)-Environment%YEFF)**2)
            theta=ATAN2((Y(j)-Environment%YEFF),(X(i)-Environment%XEFF))
            k=1
            DO WHILE ((k.LT.RadiationResults%Ntheta-1).AND.(RadiationResults%theta(k+1).LT.theta))
                k=k+1
            END DO
            IF (k.EQ.RadiationResults%Ntheta) THEN
                WRITE(*,*) ' Error: range of theta in Kochin coefficients is too small'
                STOP
            END IF
            CALL Compute_Wave(kwave,w,DiffractionResults%Rcase(iBeta),X(i),Y(j),0.,Potential,p,Vx,Vy,Vz,Environment)
            EtaI(i,j)=1./Environment%G*II*w*Potential
            HKleft=0.
            HKright=0.                
            DO l=1,RadiationResults%Ncase
                HKleft=HKleft+RAOs(l,iPeriod,iBeta)*RadiationResults%HKochin(iPeriod,l,k)
                HKright=HKright+RAOs(l,iPeriod,iBeta)*RadiationResults%HKochin(iPeriod,l,k+1)
            END DO
            HKleft=HKleft+DiffractionResults%HKochin(iPeriod,iBeta,k)
            HKright=HKright+DiffractionResults%HKochin(iPeriod,iBeta,k+1)
            HKochin=HKleft+(HKright-HKleft)*(theta-RadiationResults%theta(k))/(RadiationResults%theta(k+1)-RadiationResults%theta(k))
            IF (r.GT.0) THEN
                Potential=SQRT(kwave/(2.*PI*r))*CIH(kwave,0.,Environment%Depth)*CEXP(II*(kwave*r-0.25*PI))*HKochin
            ELSE
                Potential=0.
            END IF
            EtaP(i,j)=-II*1./Environment%G*II*w*Potential ! -II factor because of the free surface definition (Real(-II*A*exp(II*(kx-wt)))
            Eta(i,j)=EtaI(i,j)+EtaP(i,j)
        END DO
    END DO
    OPEN(10,FILE=ID%ID(1:ID%lID)//'/results/WaveField.tec')
    WRITE(10,'(A)') 'VARIABLES="X" "Y" "etaI_C" "etaI_S" "etaP_C" "etaC_S" "etaI_C+etaP_C" "etaI_S+etaI_P" "|etaP|" "|etaI+etaP|"'
    WRITE(10,'(A,E14.7,A,I6,A,I6,A)') 'ZONE t="Wave frequency - Period =',w,'",N=',Nx*Ny,', E=',(Nx-1)*(Nx-1),' , F=FEPOINT,ET=QUADRILATERAL'
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