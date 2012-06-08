    SUBROUTINE Compute_RAOs(RAOS,RadiationResults,DiffractionResults)
!    
    USE MResults
!    
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TResults) :: RadiationResults,DiffractionResults
    COMPLEX,DIMENSION(RadiationResults%Nintegration,DiffractionResults%Nperiod,*) :: RAOs
!   Locals
    INTEGER :: j,i,k
    REAL :: PI
!
    PI=4.*ATAN(1.)
    DO k=1,RadiationResults%Ncase
        DO j=1,DiffractionResults%Nperiod
            DO i=1,DiffractionResults%Ncase
                RAOs(k,j,i)=CMPLX(0.,0.)
            END DO
        END DO
    END DO  
!
    END SUBROUTINE Compute_RAOs