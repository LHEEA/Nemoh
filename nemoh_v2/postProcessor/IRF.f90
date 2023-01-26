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
    MODULE MIRF
!
    TYPE TIRF
        INTEGER :: Switch,Ntime,Nradiation,Nintegration
        REAL,DIMENSION(:),ALLOCATABLE :: Time
        REAL,DIMENSION(:,:,:),ALLOCATABLE :: K
        REAL,DIMENSION(:,:),ALLOCATABLE :: AddedMass 
    END TYPE TIRF
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTIRF(IRF,Ntime,Nradiation,Nintegration)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        INTEGER :: Ntime,Nradiation,Nintegration
        IRF%Ntime=Ntime
        IRF%Nradiation=Nradiation
        ALLOCATE(IRF%Time(0:Ntime-1),IRF%AddedMass(Nradiation,Nintegration),IRF%K(0:Ntime-1,Nradiation,Nintegration))
        END SUBROUTINE CreateTIRF
!       --- 
        SUBROUTINE CopyTIRF(IRFTarget,IRFSource)
        IMPLICIT NONE
        INTEGER :: i,j,k
        TYPE(TIRF) :: IRFTarget,IRFSource      
        CALL CreateTIRF(IRFTarget,IRFSource%Ntime,IRFSource%Nradiation,IRFSource%Nintegration)
        DO i=0,IRFTarget%Ntime-1
            IRFTarget%Time(i)=IRFSource%Time(i)
        END DO
        DO j=1,IRFTarget%Nradiation
            DO k=1,IRFTarget%Nintegration
                IRFTarget%Addedmass(j,k)=IRFSource%Addedmass(j,k)
            END DO
        END DO
        DO i=0,IRFTarget%Ntime-1
            DO j=1,IRFTarget%Nradiation
                DO k=1,IRFTarget%Nintegration
                    IRFTarget%K(i,j,k)=IRFSource%K(i,j,k)
                END DO
            END DO
        END DO
        END SUBROUTINE CopyTIRF
!       ---
        SUBROUTINE Initialize_IRF(IRF,Results,namefile)   
        USE MResults   
        IMPLICIT NONE	    
        TYPE(TIRF) :: IRF
        TYPE(TResults) :: Results
        CHARACTER*(*) :: namefile
        CHARACTER*20 :: lookfor
        CHARACTER*80 :: discard
        REAL :: dt,tf
        INTEGER :: i
        OPEN(10,FILE=namefile)
        READ(10,'(A20)') lookfor
        DO WHILE (lookfor.NE.'--- Post processing ')
            READ(10,'(A20,A)') lookfor,discard
        END DO    
        READ(10,*) IRF%Switch,dt,tf
        CLOSE(10)
        IF (IRF%Switch.EQ.1) THEN
            IRF%Ntime=INT(tf/dt)
            IRF%Nradiation=Results%Nradiation
            IRF%Nintegration=Results%Nintegration            
            CALL CreateTIRF(IRF,IRF%Ntime,IRF%Nradiation,IRF%Nintegration)
            DO i=0,IRF%Ntime-1
                IRF%Time(i)=i*dt
            END DO
        END IF        
        END SUBROUTINE  Initialize_IRF
!       ---
        SUBROUTINE Compute_IRF(IRF,Results)
        USE MResults
        IMPLICIT NONE	    
        TYPE(TIRF) :: IRF
        TYPE(TResults) :: Results
        REAL,DIMENSION(Results%Nw) :: CM
        INTEGER :: i,j,k,l
        REAL :: PI
        PI=4.*ATAN(1.0)
        DO i=0,IRF%Ntime-1
            DO j=1,Results%Nradiation
                DO k=1,Results%Nintegration
                    IRF%K(i,j,k)=0.
                    DO l=1,Results%Nw-1
!                        IRF%K(i,j,k)=IRF%K(i,j,k)-REAL(Results%Force(l,j,k))*COS(2.*PI/Results%Period(l)*IRF%Time(i))*2.*PI*ABS(1./Results%Period(l)-1./Results%Period(l+1))
                        IRF%K(i,j,k)=IRF%K(i,j,k)+0.5*(Results%RadiationDamping(l,j,k)*COS(Results%w(l)*IRF%Time(i))+Results%RadiationDamping(l+1,j,k)*COS(Results%w(l+1)*IRF%Time(i)))*(Results%w(l+1)-Results%w(l))                    
                    END DO 
                    IRF%K(i,j,k)=IRF%K(i,j,k)*2./PI         
                END DO
            END DO
        END DO
        DO j=1,Results%Nradiation
            DO k=1,Results%Nintegration
                IRF%AddedMass(j,k)=0.
                DO l=1,Results%Nw
                    CM(l)=0.
                    DO i=0,IRF%Ntime-2
 !                       CM(l)=CM(l)+IRF%K(i,j,k)*SIN(2.*PI/Results%Period(l)*IRF%Time(i))*(IRF%Time(i+1)-IRF%Time(i))
                       CM(l)=CM(l)+0.5*(IRF%K(i,j,k)*SIN(Results%w(l)*IRF%Time(i))+IRF%K(i+1,j,k)*SIN(Results%w(l)*IRF%Time(i+1)))*(IRF%Time(i+1)-IRF%Time(i))
                    END DO  
                    CM(l)=Results%AddedMass(l,j,k)+CM(l)/Results%w(l)
                    IRF%AddedMass(j,k)=IRF%AddedMass(j,k)+CM(l)   
                END DO         
                IRF%AddedMass(j,k)=IRF%AddedMass(j,k)/Results%Nw
            END DO
        END DO
        END SUBROUTINE Compute_IRF
!       ---
        SUBROUTINE Save_IRF(IRF,namefile)
        IMPLICIT NONE	    
        TYPE(TIRF) :: IRF
        CHARACTER*(*) :: namefile
        INTEGER :: i,j,k
        OPEN(10,FILE=namefile)
        WRITE(10,*) 'VARIABLES="Time (s)"'
        DO k=1,IRF%Nintegration
            WRITE(10,'(A,I4,A,I4,A)') '"AddedMass ',k,'" "IRF ',k,'"'
        END DO
        DO j=1,IRF%Nradiation
            WRITE(10,'(A,I4,A,I6,A)') 'Zone t="DoF ',j,'",I=',IRF%Ntime,',F=POINT'
            DO i=0,IRF%Ntime-1                
                WRITE(10,'(80(X,E14.7))') IRF%Time(i),(IRF%AddedMass(j,k),IRF%K(i,j,k),k=1,IRF%Nintegration)
            END DO
        END DO
        CLOSE(10)
        END SUBROUTINE Save_IRF
!       --- 
        SUBROUTINE DeleteTIRF(IRF)
        IMPLICIT NONE
        TYPE(TIRF) :: IRF
        DEALLOCATE(IRF%Time,IRF%K,IRF%AddedMass)
        END SUBROUTINE DeleteTIRF  
!       ---
END MODULE