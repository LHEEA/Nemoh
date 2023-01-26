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
!   - A. Babarit / Ecole Centrale de Nantes 
!   - C. Peyrard / EDF R&D
!
!--------------------------------------------------------------------------------------
    MODULE MResults
!
    TYPE TResults
        INTEGER :: Nw,Nbeta,Nradiation,Nintegration
        REAL,DIMENSION(:),ALLOCATABLE :: w
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: IndxForce
        INTEGER,DIMENSION(:,:),ALLOCATABLE :: IndxRadiation
        REAL,DIMENSION(:),ALLOCATABLE :: beta
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: DiffractionForce
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: FroudeKrylovForce
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: AddedMass
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: RadiationDamping
        INTEGER :: Ntheta
        REAL,DIMENSION(:),ALLOCATABLE :: Theta
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: HKochinDiffraction
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: HKochinRadiation   
    END TYPE TResults
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTResults(Results,Nw,Nbeta,Nradiation,Nintegration,Ntheta)
        IMPLICIT NONE
        TYPE(TResults) :: Results
        INTEGER :: Nw,Nbeta,Nradiation,Nintegration,Ntheta
        Results%Nw=Nw
        Results%Nbeta=Nbeta
        Results%Nradiation=Nradiation
        Results%Nintegration=Nintegration
        ALLOCATE(Results%w(Nw),Results%IndxForce(Nintegration,3),Results%beta(Nbeta),Results%IndxRadiation(Nintegration,3))
        ALLOCATE(Results%DiffractionForce(Nw,Nbeta,Nintegration),Results%FroudeKrylovForce(Nw,Nbeta,Nintegration),Results%AddedMass(Nw,Nradiation,Nintegration),Results%RadiationDamping(Nw,Nradiation,Nintegration))
        Results%Ntheta=Ntheta
        IF (Results%Ntheta.GT.0) THEN            
            ALLOCATE(Results%Theta(Ntheta),Results%HKochinDiffraction(Nw,Nbeta,Ntheta),Results%HKochinRadiation(Nw,Nradiation,Ntheta))
        END IF
        END SUBROUTINE CreateTResults
!       --- 
        SUBROUTINE CopyTResults(ResultsTarget,ResultsSource)
        IMPLICIT NONE
        INTEGER :: i,j,k
        TYPE(TResults) :: ResultsTarget,ResultsSource      
        CALL CreateTResults(ResultsTarget,ResultsSource%Nw,ResultsSource%Nbeta,ResultsSource%Nradiation,ResultsSource%Nintegration,ResultsSource%Ntheta)
        DO i=1,ResultsTarget%Nw
            ResultsTarget%w(i)=ResultsSource%w(i)
        END DO
        DO i=1,ResultsTarget%Nbeta
            ResultsTarget%beta(i)=ResultsSource%beta(i)
        END DO
        DO i=1,ResultsTarget%Nintegration
            DO j=1,3
                ResultsTarget%IndxForce(i,j)=ResultsSource%IndxForce(i,j)
                ResultsTarget%IndxRadiation(i,j)=ResultsSource%IndxRadiation(i,j)
            END DO
        END DO
        DO j=1,ResultsTarget%Nw
            DO i=1,ResultsTarget%Nbeta
                DO k=1,ResultsTarget%Nintegration
                    ResultsTarget%DiffractionForce(j,i,k)=ResultsSource%DiffractionForce(j,i,k)
                    ResultsTarget%FroudeKrylovForce(j,i,k)=ResultsSource%FroudeKrylovForce(j,i,k)
                END DO
            END DO
            DO i=1,ResultsTarget%Nbeta
                DO k=1,ResultsTarget%Nintegration
                    ResultsTarget%AddedMass(j,i,k)=ResultsSource%AddedMass(j,i,k)
                    ResultsTarget%RadiationDamping(j,i,k)=ResultsSource%RadiationDamping(j,i,k)
                END DO
            END DO
        END DO
        DO k=1,ResultsTarget%Ntheta
            ResultsTarget%Theta(k)=ResultsSource%Theta(k)
            DO j=1,ResultsTarget%Nw
                DO i=1,ResultsTarget%Nbeta               
                    ResultsTarget%HKochinDiffraction(j,i,k)=ResultsSource%HKochinDiffraction(j,i,k)
                END DO
                DO i=1,ResultsTarget%Nradiation              
                    ResultsTarget%HKochinRadiation(j,i,k)=ResultsSource%HKochinRadiation(j,i,k)
                END DO
            END DO
        END DO
        END SUBROUTINE CopyTResults
!       ---
        SUBROUTINE ReadTResults(Results,namefile,nameindex,namefileFK)
        IMPLICIT NONE	    
        TYPE(TResults) :: Results
        CHARACTER*(*) :: namefile,nameindex,namefileFK
        INTEGER :: Nw,Nbeta,Nradiation,Nintegration,Ntheta
        INTEGER :: i,j,k,c
        REAL :: discard,amplitude,phase
        REAL,DIMENSION(:),ALLOCATABLE :: line
        OPEN(10,FILE=nameindex)
        READ(10,*) Nw,Nbeta,Nradiation,Nintegration,Ntheta
        CALL CreateTResults(Results,Nw,Nbeta,Nradiation,Nintegration,Ntheta)
        READ(10,*)
        DO k=1,Nintegration
            READ(10,*) Results%IndxForce(k,1),Results%IndxForce(k,2),Results%IndxForce(k,3)
        END DO
        READ(10,*)
        DO k=1,Nradiation
            READ(10,*) Results%IndxRadiation(k,1),Results%IndxRadiation(k,2),Results%IndxRadiation(k,3)
        END DO
        READ(10,*) (Results%beta(k),k=1,Nbeta)
        READ(10,*) (Results%w(k),k=1,Nw)
        READ(10,*) (Results%theta(k),k=1,Ntheta)
        CLOSE(10)
        ALLOCATE(line(2*Nw*(Nbeta+Nradiation)))
        OPEN(10,FILE=namefile)        
        DO k=1,Nintegration
            READ(10,*) (line(j),j=1,2*Nw*(Nbeta+Nradiation))
            c=1
            DO i=1,Nw
                DO j=1,Nbeta
                    Results%DiffractionForce(i,j,k)=line(c)*CEXP(CMPLX(0.,1.)*line(c+1))
                    c=c+2
                END DO
                DO j=1,Nradiation
                    Results%AddedMass(i,j,k)=line(c)
                    Results%RadiationDamping(i,j,k)=line(c+1)
                    c=c+2
                 END DO
            END DO
        END DO
        CLOSE(10)
        DEALLOCATE(line)
        ALLOCATE(line(2*Nintegration+1))
        OPEN(10,FILE=namefileFK)
        READ(10,*)
        DO k=1,Nintegration
            READ(10,*)
        END DO
        DO j=1,Nbeta
            READ(10,*)
            DO i=1,Nw
                READ(10,*) (line(k),k=1,1+2*Nintegration)
                c=2
                DO k=1,Nintegration
                    Results%FroudeKrylovForce(i,j,k)=line(c)*CEXP(CMPLX(0.,1.)*line(c+1))
                    c=c+2
                END DO
            END DO
        END DO
        CLOSE(10)
        DEALLOCATE(line)
        END SUBROUTINE ReadTResults
!       ---
        SUBROUTINE SaveTResults(Results,namedir)
        IMPLICIT NONE	    
        TYPE(TResults) :: Results
        CHARACTER*(*) :: namedir
        INTEGER :: i,j,k,l
        REAL :: PI
        PI=4.*ATAN(1.0)
        OPEN(10,FILE=namedir//'/RadiationCoefficients.tec')
        WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
        DO k=1,Results%Nintegration
            WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"A',Results%IndxForce(k,2),Results%IndxForce(k,3),'" "B',Results%IndxForce(k,2),Results%IndxForce(k,3),'"'
        END DO
        DO j=1,Results%Nradiation
            WRITE(10,'(A,I4,A,I4,A,I6,A)') 'Zone t="Motion of body ',Results%IndxRadiation(j,2),' in DoF',Results%IndxRadiation(j,3),'",I=',Results%Nw,',F=POINT'
            DO i=1,Results%Nw
                WRITE(10,'(80(X,E14.7))') Results%w(i),(Results%AddedMass(i,j,k),Results%RadiationDamping(i,j,k),k=1,Results%Nintegration)
            END DO
        END DO
        CLOSE(10)
        OPEN(10,FILE=namedir//'/DiffractionForce.tec')
        WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
        DO k=1,Results%Nintegration
            WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',Results%IndxForce(k,2),Results%IndxForce(k,3),')" "angle(F',Results%IndxForce(k,2),Results%IndxForce(k,3),')"'
        END DO
        DO j=1,Results%Nbeta
            WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="Diffraction force - beta = ',Results%beta(j)*180./(4.*ATAN(1.0)),' deg",I=',Results%Nw,',F=POINT'
            DO i=1,Results%Nw
                WRITE(10,'(80(X,E14.7))') Results%w(i),(ABS(Results%DiffractionForce(i,j,k)),ATAN2(IMAG(Results%DiffractionForce(i,j,k)),REAL(Results%DiffractionForce(i,j,k))),k=1,Results%Nintegration)
            END DO
        END DO
        CLOSE(10)
        OPEN(10,FILE=namedir//'/ExcitationForce.tec')
        WRITE(10,'(A)') 'VARIABLES="w (rad/s)"'
        DO k=1,Results%Nintegration
            WRITE(10,'(A,I4,I4,A,I4,I4,A)') '"abs(F',Results%IndxForce(k,2),Results%IndxForce(k,3),')" "angle(F',Results%IndxForce(k,2),Results%IndxForce(k,3),')"'
        END DO
        DO j=1,Results%Nbeta
            WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="Diffraction force - beta = ',Results%beta(j)*180./(4.*ATAN(1.0)),' deg",I=',Results%Nw,',F=POINT'
            DO i=1,Results%Nw
!				Adrien Combourieu: in ExcitationForce.tec, the phase is given in radians.
                WRITE(10,'(80(X,E14.7))') Results%w(i),(ABS(Results%DiffractionForce(i,j,k)+Results%FroudeKrylovForce(i,j,k)),ATAN2(IMAG(Results%DiffractionForce(i,j,k)+Results%FroudeKrylovForce(i,j,k)),REAL(Results%DiffractionForce(i,j,k)+Results%FroudeKrylovForce(i,j,k))),k=1,Results%Nintegration)
            END DO
        END DO
        CLOSE(10)
!       Added by Christophe Peyrard
!       Save hydrodynamic database with Aquaplus format
        OPEN(10,FILE=namedir//'/CA.dat')
        WRITE(10,'(A,I5)') 'Nb de periode : ',Results%Nw
        DO l=1,Results%Nw
	        WRITE(10,'(F7.4)') Results%w(l)
	        DO j=1,Results%Nradiation
	            WRITE(10,'(6(X,E13.6))') (Results%RadiationDamping(l,j,k),k=1,Results%Nradiation)
	        END DO
        END DO
        CLOSE(10)
        OPEN(10,FILE=namedir//'/CM.dat')
        WRITE(10,'(A,I5)') 'Nb de periode : ',Results%Nw
        DO l=1,Results%Nw
	        WRITE(10,'(F7.4)') Results%w(l)
	        DO j=1,Results%Nradiation
	            WRITE(10,'(6(X,E13.6))') (Results%AddedMass(l,j,k),k=1,Results%Nradiation)
	        END DO
        END DO
        CLOSE(10)
        OPEN(10,FILE=namedir//'/Fe.dat')
        WRITE(10,'(A)') 'VARIABLES="Period (s)" "|Fx| (N/m)" "|Fy| (N/m)" "|Fz| (N/m)" "|Cx| (N)" "|Cy| (N)" "|Cz| (N)" "ang(Fx) (°)" "ang(Fy) (°)" "ang(Fz) (°)" "ang(Cx) (°)" "ang(Cy) (°)" "ang(Cz) (°)"'
        WRITE(10,'(A,I2,A)') 'Zone t="Corps ',1,'"'
        DO j=1,Results%Nbeta
            WRITE(10,'(A,F7.3,A,I6,A)') 'Zone t="Diffraction force - beta = ',Results%beta(j)*180./PI,' deg",I=',Results%Nw,',F=POINT'
            DO l=1,Results%Nw
!           Be careful of the phase referential when comparing with Aquaplus
!           In Nemoh, PHI = -g/w * CIH CEXP(k.x), so ETA= i w/g * PHI = -i CEXP(k.x)
!           Thus, in time domain ETA = sin(k.x -wt) 
!           Forces/Phase shift in Fe.dat are given with the following form : F= Force * sin(-wt-phase), so a phase shift of -PI/2 is necessary
                WRITE(10,'(F7.4,6(X,E13.6),6(X,F7.2))') Results%w(l),(ABS(Results%DiffractionForce(l,j,k)+Results%FroudeKrylovForce(l,j,k)),k=1,Results%Nintegration),  &
                                                            (180.D0/PI*( ATAN2(IMAG(Results%DiffractionForce(l,j,k)+Results%FroudeKrylovForce(l,j,k)),        &
!                                                            (-180.D0/PI*( ATAN2(IMAG(Results%DiffractionForce(l,j,k)+Results%FroudeKrylovForce(l,j,k)),        &
                                                                               REAL(Results%DiffractionForce(l,j,k)+Results%FroudeKrylovForce(l,j,k)))),        &
        !                                                                         -PI/2D0  ),                                                                         &
                                                               k=1,Results%Nintegration)
            END DO
        END DO
        CLOSE(10)
!       End of addition
        END SUBROUTINE SaveTResults
!       --- 
        SUBROUTINE DeleteTResults(Results)
        IMPLICIT NONE
        TYPE(TResults) :: Results
        DEALLOCATE(Results%w,Results%beta,Results%IndxForce,Results%IndxRadiation)
        DEALLOCATE(Results%DiffractionForce,Results%FroudeKrylovForce,Results%AddedMass,Results%RadiationDamping)
        IF (Results%Ntheta.GT.0) DEALLOCATE(Results%Theta,Results%HKochinRadiation,Results%HKochinDiffraction)
        END SUBROUTINE DeleteTResults  
!       ---
END MODULE