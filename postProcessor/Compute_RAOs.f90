    SUBROUTINE Compute_RAOs(RAOS,Results)
!    
    USE MResults
!    
    IMPLICIT NONE
!
!   Inputs/outputs
    TYPE(TResults) :: Results
    COMPLEX,DIMENSION(Results%Nintegration,Results%Nw,*) :: RAOs
!   Locals
    INTEGER :: j,i,k
    REAL :: PI
!
    PI=4.*ATAN(1.)
    DO k=1,Results%Nradiation
        DO j=1,Results%Nw
            DO i=1,Results%Nbeta
                RAOs(k,j,i)=CMPLX(0.,0.)
            END DO
        END DO
    END DO  
!
    END SUBROUTINE Compute_RAOs