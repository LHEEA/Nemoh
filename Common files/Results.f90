    MODULE MResults
!
    TYPE TResults
        INTEGER :: Nperiod,Ncase,Nintegration
        REAL,DIMENSION(:),ALLOCATABLE :: Period
        REAL,DIMENSION(:),ALLOCATABLE :: Rcase
        INTEGER,DIMENSION(:),ALLOCATABLE :: Iintegration
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: Force
        INTEGER :: Ntheta
        REAL,DIMENSION(:),ALLOCATABLE :: Theta
        COMPLEX,DIMENSION(:,:,:),ALLOCATABLE :: HKochin        
    END TYPE TResults
!
    CONTAINS
!
!       Operators for creation, copy, initialisation and destruction
!
        SUBROUTINE CreateTResults(Results,Nperiod,Ncase,Nintegration,Ntheta)
        IMPLICIT NONE
        TYPE(TResults) :: Results
        INTEGER :: Nperiod,Ncase,Nintegration,Ntheta
        Results%Nperiod=Nperiod
        Results%Ncase=Ncase
        Results%Nintegration=Nintegration
        ALLOCATE(Results%Period(Nperiod),Results%Rcase(Ncase),Results%Iintegration(Nintegration))
        ALLOCATE(Results%Force(Nperiod,Ncase,Nintegration))
        Results%Ntheta=Ntheta
        IF (Results%Ntheta.GT.0) THEN            
            ALLOCATE(Results%Theta(Ntheta),Results%HKochin(Nperiod,Ncase,Ntheta))
        END IF
        END SUBROUTINE CreateTResults
!       --- 
        SUBROUTINE CopyTResults(ResultsTarget,ResultsSource)
        IMPLICIT NONE
        INTEGER :: i,j,k
        TYPE(TResults) :: ResultsTarget,ResultsSource      
        CALL CreateTResults(ResultsTarget,ResultsSource%Nperiod,ResultsSource%Ncase,ResultsSource%Nintegration,ResultsSource%Ntheta)
        DO i=1,ResultsTarget%Nperiod
            ResultsTarget%Period(i)=ResultsSource%Period(i)
        END DO
        DO i=1,ResultsTarget%Ncase
            ResultsTarget%Rcase(i)=ResultsSource%Rcase(i)
        END DO
        DO i=1,ResultsTarget%Nintegration
            ResultsTarget%Iintegration(i)=ResultsSource%Iintegration(i)
        END DO
        DO j=1,ResultsTarget%Nperiod
            DO i=1,ResultsTarget%Ncase
                DO k=1,ResultsTarget%Nintegration
                    ResultsTarget%Force(j,i,k)=ResultsSource%Force(j,i,k)
                END DO
            END DO
        END DO
        DO k=1,ResultsTarget%Ntheta
            ResultsTarget%Theta(k)=ResultsSource%Theta(k)
            DO j=1,ResultsTarget%Nperiod
                DO i=1,ResultsTarget%Ncase                
                    ResultsTarget%HKochin(j,i,k)=ResultsSource%HKochin(j,i,k)
                END DO
            END DO
        END DO
        END SUBROUTINE CopyTResults
!       ---
        SUBROUTINE ReadTResults(Results,namefile)
        IMPLICIT NONE	    
        TYPE(TResults) :: Results
        CHARACTER*(*) :: namefile
        INTEGER :: Nperiod,Ncase,Nintegration,Ntheta
        INTEGER :: i,j,k
        OPEN(10,FILE=namefile)
        READ(10,*) Nperiod,Ncase,Nintegration,Ntheta
        CALL CreateTResults(Results,Nperiod,Ncase,Nintegration,Ntheta)
        READ(10,*) (Results%Period(i),i=1,Nperiod)
        DO k=1,Nintegration
            DO j=1,Ncase
                READ(10,*) Results%Iintegration(k),Results%Rcase(j),(Results%Force(i,j,k),i=1,Nperiod)
            END DO
        END DO
        DO k=1,Ntheta
            DO j=1,Ncase
                READ(10,*) Results%Theta(k),(Results%HKochin(i,j,k),i=1,Nperiod)
            END DO
        END DO
        CLOSE(10)
        END SUBROUTINE ReadTResults
!       ---
        SUBROUTINE SaveTResults(Results,namefile)
        IMPLICIT NONE	    
        TYPE(TResults) :: Results
        CHARACTER*(*) :: namefile
        INTEGER :: i,j,k
        OPEN(10,FILE=namefile)
        WRITE(10,*) Results%Nperiod,Results%Ncase,Results%Nintegration,Results%Ntheta
        WRITE(10,*) (Results%Period(i),i=1,Results%Nperiod)
        DO k=1,Results%Nintegration
            DO j=1,Results%Ncase
                WRITE(10,*) Results%Iintegration(k),Results%Rcase(j),(Results%Force(i,j,k),i=1,Results%Nperiod)
            END DO
        END DO
        DO k=1,Results%Ntheta
            DO j=1,Results%Ncase
                WRITE(10,*) Results%Theta(k),(Results%HKochin(i,j,k),i=1,Results%Nperiod)
            END DO
        END DO
        CLOSE(10)
        END SUBROUTINE SaveTResults
!       --- 
        SUBROUTINE DeleteTResults(Results)
        IMPLICIT NONE
        TYPE(TResults) :: Results
        DEALLOCATE(Results%Force,Results%Period,Results%Rcase,Results%Iintegration)
        IF (Results%Ntheta.GT.0) DEALLOCATE(Results%Theta,Results%HKochin)
        END SUBROUTINE DeleteTResults  
!       ---
END MODULE