MODULE MIdentification
! Definition of TYPE TID
TYPE TID
    CHARACTER*80 :: ID
    INTEGER :: lID
END TYPE TID
!
CONTAINS
!   Read ID data from file
    SUBROUTINE ReadTID(ID,file)
    IMPLICIT NONE
    CHARACTER*(*) :: file
    TYPE(TID) :: ID
    OPEN(10,FILE=file)
    READ(10,*) ID%lID
    READ(10,*) ID%ID
    CLOSE(10)  
    END SUBROUTINE
!    
END MODULE MIdentification

