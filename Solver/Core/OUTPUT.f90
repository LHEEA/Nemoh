MODULE OUTPUT

CONTAINS
!
    SUBROUTINE WRITE_KOCHIN(ID,Pbnumber,HKochin,Ntheta)
    USE MIDENTIFICATION
    IMPLICIT NONE
    TYPE(TID) :: ID
    INTEGER :: Pbnumber
    COMPLEX,DIMENSION(*) :: HKochin
    INTEGER :: Ntheta,i
    CHARACTER*5 :: str
    WRITE(str,'(I5)') Pbnumber
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/Kochin.'//str//'.dat')
    DO i=1,Ntheta
        WRITE(11,*) ABS(HKochin(i)),ATAN2(IMAG(HKochin(i)),REAL(HKochin(i)))
    END DO
    CLOSE(11)
    END SUBROUTINE WRITE_KOCHIN
!
    SUBROUTINE WRITE_FS(ID,Pbnumber,PHI,MeshFS)
    USE MIDENTIFICATION
    USE MMesh
    IMPLICIT NONE
    INTEGER :: Pbnumber
    TYPE(TID) :: ID
    TYPE(TMesh) :: MeshFS
    INTEGER:: i
    COMPLEX,DIMENSION(*) :: PHI
    CHARACTER*5 :: str
    WRITE(str,'(I5)') Pbnumber
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/freesurface.'//str//'.dat')
    DO i=1,MeshFS%Npoints
        WRITE(11,*) ABS(PHI(i)),ATAN2(IMAG(PHI(i)),REAL(PHI(i))) 
    END DO
    CLOSE(11)
    END SUBROUTINE WRITE_FS 
!
    SUBROUTINE WRITE_POTENTIAL(ID,Pbnumber)
    USE MIDENTIFICATION
    USE COM_VAR
    IMPLICIT NONE
    INTEGER :: Pbnumber
    TYPE(TID) :: ID
    INTEGER:: i
    CHARACTER*5 :: str
    WRITE(str,'(I5)') Pbnumber
    OPEN(11,FILE=ID%ID(1:ID%lID)//'/results/potential.'//str//'.dat')
    DO i=1,NFA
        WRITE(11,*) ABS(ZPB(i)),ATAN2(IMAG(ZPB(i)),REAL(ZPB(i))) 
    END DO
    IF (NSYMY.EQ.1) THEN
        DO i=1,NFA
            WRITE(11,*) ABS(ZPS(i)),ATAN2(IMAG(ZPS(i)),REAL(ZPS(i))) 
        END DO
    END IF
    CLOSE(11) 
    END SUBROUTINE
  
 END MODULE